#!/usr/bin/env python3
"""
Parse a QRS CIF (multi-block), select structures with energy below median,
relax them with MACE (no lattice optimization) in parallel, and plot results.

Usage: python scripts/relax_qrs_mace.py /path/to/QRS-openffall.cif
       python scripts/relax_qrs_mace.py /path/to/Tests-0607/XAFQON
"""
import os
import re
from io import StringIO, BytesIO
import multiprocessing as mp
import numpy as np
import argparse
import time

import matplotlib.pyplot as plt
from ase.io import read, write
import warnings

warnings.filterwarnings('ignore', category=DeprecationWarning, module='spglib')
warnings.filterwarnings('ignore', category=UserWarning, module=r'ase\.io')
warnings.filterwarnings('ignore', category=UserWarning, module='torch_dftd')

from pyxtal.interface.ase_opt import (
    ASE_relax,
    resolve_model,
    DEFAULT_MODELS,
    get_calculator,
    ensure_mace_model_cached,
)
from pyxtal.util import ase2pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
from pyxtal.db import database as PyDatabase


def parse_qrs_cif(path):
    """Yield tuples (label, block_text, energy)
    Each block runs from a line starting with 'data_' until a line with '#END'.
    Energy is parsed from a line like '#Energy: -36.790237 eV/cell'.
    """
    blocks = []
    with open(path, 'r') as f:
        lines = f.readlines()

    cur = []
    cur_label = None
    cur_energy = None
    in_block = False
    for ln in lines:
        if ln.startswith('data_'):
            # start a new block
            in_block = True
            cur = [ln]
            cur_label = ln.strip()
            cur_energy = None
            continue
        if in_block:
            cur.append(ln)
            if ln.startswith('#Energy:'):
                # parse float
                try:
                    tok = ln.strip().split()[1]
                    cur_energy = float(tok)
                except Exception:
                    # try alternative parsing
                    try:
                        cur_energy = float(ln.strip().split(':',1)[1].split()[0])
                    except Exception:
                        cur_energy = None
            if ln.strip() == '#END':
                blocks.append((cur_label, ''.join(cur), cur_energy))
                in_block = False
                cur = []
                cur_label = None
                cur_energy = None
    return blocks


def get_matched_labels_from_file(path):
    """Return a set of labels in a matched CIF or label list file."""
    labels = set()
    def normalize_label(lbl):
        if not lbl.startswith('data_'):
            return lbl
        # Strip extra matched-label suffixes like '-QRandom-0.17-0.39'.
        first = lbl.find('-QRandom-')
        if first != -1:
            second = lbl.find('-QRandom-', first + 1)
            if second != -1:
                lbl = lbl[:second]
        # Normalize energy suffixes such as '-e-46.615' -> '-e46.615'.
        lbl = re.sub(r'-e-([0-9])', r'-e\1', lbl)
        return lbl
    # Try parsing as a CIF file first
    try:
        blocks = parse_qrs_cif(path)
        labels.update([normalize_label(lbl) for lbl, txt, eng in blocks if lbl])
        if labels:
            return labels
    except Exception:
        pass

    # Fallback: treat as a plain label list, but only accept lines that start with 'data_'
    with open(path, 'r') as f:
        for ln in f:
            line = ln.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('data_'):
                labels.add(normalize_label(line.split()[0]))
    return labels


def guess_ref_code_from_cif(path):
    """Try to infer a reference code from the CIF file if none was provided."""
    ref_code = None
    with open(path, 'r') as f:
        for ln in f:
            text = ln.strip()
            if not text:
                continue
            low = text.lower()
            if low.startswith('#ref_code') or low.startswith('#ref-code') or low.startswith('#ref code'):
                if ':' in text:
                    ref_code = text.split(':', 1)[1].strip()
                else:
                    ref_code = text.split(None, 1)[1].strip() if len(text.split(None, 1)) > 1 else None
                break
            if low.startswith('#ref:') or low.startswith('#ref '):
                ref_code = text.split(':', 1)[1].strip() if ':' in text else text.split(None, 1)[1].strip()
                break
            if low.startswith('_ref_code') or low.startswith('_ref-code') or low.startswith('_reference'):
                parts = text.split()
                if len(parts) >= 2:
                    ref_code = parts[1].strip().strip("'\"")
                break
            if low.startswith('references:'):
                ref_path = text.split(':', 1)[1].strip().strip("'\"")
                if ref_path:
                    ref_code = os.path.basename(os.path.dirname(ref_path))
                    if ref_code:
                        break
            if low.startswith('data_'):
                label = text[5:]
                if '-g' in label:
                    ref_code = label.split('-g', 1)[0]
                    break
    if ref_code:
        ref_code = ref_code.strip().strip("'\"")
    if not ref_code:
        base_dir = os.path.basename(os.path.dirname(os.path.abspath(path)))
        if base_dir and not base_dir.lower().startswith('tests'):
            ref_code = base_dir
    return ref_code


def _get_pmg_no_h(entry, cache):
    idx = entry['idx']
    if idx not in cache:
        atoms = read(StringIO(entry['text']), format='cif')
        pmg = ase2pymatgen(atoms)
        if hasattr(pmg, 'remove_species'):
            pmg.remove_species(['H'])
        cache[idx] = pmg
    return cache[idx]


def _pmg_for_matching(pmg, max_sites=50):
    """Return a pymatgen structure truncated for fast StructureMatcher comparisons."""
    if max_sites is None or max_sites <= 0 or len(pmg) <= max_sites:
        return pmg
    from pymatgen.core import Structure
    return Structure.from_sites(pmg.sites[:max_sites], validate_proximity=False)


def _structures_likely_match(pmg1, pmg2, volume_tol=0.05, max_sites=50):
    m1 = _pmg_for_matching(pmg1, max_sites)
    m2 = _pmg_for_matching(pmg2, max_sites)
    if len(m1) != len(m2):
        return False
    if m1.composition.reduced_formula != m2.composition.reduced_formula:
        return False
    v1, v2 = pmg1.volume, pmg2.volume
    if v1 <= 0 or v2 <= 0:
        return False
    return abs(v1 - v2) / max(v1, v2) <= volume_tol


def _deduplicate_selected(selected, matcher, e_tol=1e-3, progress_every=25, match_max_sites=50):
    """Drop duplicate structures among selected entries (lowest energy kept)."""
    import bisect

    unique_selected = []
    rep_energies = []
    dup_map = {}
    pmg_cache = {}
    n_fit = 0
    dedup_start = time.perf_counter()
    if match_max_sites and match_max_sites > 0:
        print(f'Deduplication uses first {match_max_sites} sites per structure for matching')

    for i, ent in enumerate(sorted(selected, key=lambda x: x['energy']), start=1):
        if progress_every and i % progress_every == 0:
            elapsed = time.perf_counter() - dedup_start
            print(
                f'Deduplicating [{i}/{len(selected)}]: '
                f'{len(unique_selected)} unique so far, {n_fit} structure matches, '
                f'elapsed {elapsed:.1f}s',
                flush=True,
            )

        try:
            pmg1 = _get_pmg_no_h(ent, pmg_cache)
        except Exception:
            unique_selected.append(ent)
            rep_energies.append(ent['energy'])
            continue

        lo = bisect.bisect_left(rep_energies, ent['energy'] - e_tol)
        hi = bisect.bisect_right(rep_energies, ent['energy'] + e_tol)
        is_dup = False
        for rep in unique_selected[lo:hi]:
            try:
                pmg2 = _get_pmg_no_h(rep, pmg_cache)
            except Exception:
                continue
            if not _structures_likely_match(pmg1, pmg2, max_sites=match_max_sites):
                continue
            n_fit += 1
            if matcher.fit(
                _pmg_for_matching(pmg1, match_max_sites),
                _pmg_for_matching(pmg2, match_max_sites),
            ):
                dup_map[ent['idx']] = rep['idx']
                is_dup = True
                break

        if not is_dup:
            unique_selected.append(ent)
            rep_energies.append(ent['energy'])

    if len(selected) >= progress_every:
        elapsed = time.perf_counter() - dedup_start
        print(
            f'Deduplication done: {len(unique_selected)} unique from {len(selected)} selected '
            f'({n_fit} structure matches), elapsed {elapsed:.1f}s',
            flush=True,
        )
    return unique_selected, dup_map


def _configure_worker_threads():
    """Avoid CPU oversubscription when running multiple PyTorch workers."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    try:
        import torch
        torch.set_num_threads(1)
    except ImportError:
        pass


def _worker_ping(_):
    return os.getpid()


def _init_worker(calculator, model, quick):
    """Load the calculator once per worker process (safe with spawn)."""
    import warnings
    warnings.filterwarnings('ignore', category=DeprecationWarning, module='spglib')
    warnings.filterwarnings('ignore', category=UserWarning, module=r'ase\.io')
    warnings.filterwarnings('ignore', category=UserWarning, module='torch_dftd')

    _configure_worker_threads()
    t0 = time.perf_counter()
    get_calculator(calculator, model=model, quick=quick)
    print(
        f'Worker {os.getpid()} loaded {calculator} model={model} in {time.perf_counter() - t0:.2f}s',
        flush=True,
    )


def worker_relax(args):
    """Worker to relax a single CIF block.
    args: (label, block_text, orig_energy, step, fmax, label_idx, calculator, model, quick)
    Returns (label, label_idx, original_energy, model_energy or None, status, relaxed_cif, elapsed_s)
    """
    label, block_text, orig_energy, step, fmax, label_idx, calculator, model, quick = args
    start_time = time.time()
    try:
        atoms = read(StringIO(block_text), format='cif')
    except Exception as e:
        return (label, label_idx, orig_energy, None, f"read_error: {e}", None, time.time() - start_time)

    try:
        # Run chosen MLP relaxation with no lattice optimization
        relaxed = ASE_relax(
            atoms,
            calculator=calculator,
            model=model,
            quick=quick,
            opt_lat=False,
            step=step,
            fmax=fmax,
            max_time=30.0,
            label=label,
            logfile=os.devnull,
        )
        elapsed = time.time() - start_time
        if relaxed is None:
            return (label, label_idx, orig_energy, None, 'relax_failed', None, elapsed)
        model_energy = relaxed.get_potential_energy()
        print(f'{label[-45:]} idx={label_idx:5d}, {calculator} {model_energy:.3f}, time {elapsed:.1f}s')
        # export relaxed structure as CIF text for downstream matching
        sio = BytesIO()
        write(sio, relaxed, format='cif')
        relaxed_cif = sio.getvalue().decode('utf-8')
        return (label, label_idx, orig_energy, model_energy, 'ok', relaxed_cif, elapsed)
    except Exception as e:
        return (label, label_idx, orig_energy, None, f'relax_exception: {e}', None, time.time() - start_time)


def find_qrs_job_dirs(root):
    """Return (input_dir, output_dir) pairs under root that contain QRS-openffall.cif."""
    root = os.path.abspath(root)
    direct_cif = os.path.join(root, 'QRS-openffall.cif')
    if os.path.isfile(direct_cif):
        return [(root, root)]

    jobs = []
    for name in sorted(os.listdir(root)):
        sub = os.path.join(root, name)
        if os.path.isdir(sub) and os.path.isfile(os.path.join(sub, 'QRS-openffall.cif')):
            jobs.append((sub, sub))
    return jobs


def collect_jobs(cif_path, out_dir=None):
    """Resolve one or many QRS jobs from a CIF path and optional output directory."""
    cif_path = os.path.abspath(cif_path)

    if os.path.isfile(cif_path):
        parent = os.path.dirname(cif_path)
        job_out = os.path.abspath(out_dir) if out_dir else parent
        return [(parent, job_out)]

    if out_dir is not None:
        out_dir = os.path.abspath(out_dir)
        if os.path.isdir(out_dir) and not os.path.isfile(os.path.join(out_dir, 'QRS-openffall.cif')):
            jobs = []
            for sub, _ in find_qrs_job_dirs(out_dir):
                jobs.append((sub, sub))
            if jobs:
                return jobs

    if os.path.isdir(cif_path):
        if not os.path.isfile(os.path.join(cif_path, 'QRS-openffall.cif')):
            jobs = find_qrs_job_dirs(cif_path)
            if jobs:
                if out_dir is not None:
                    out_root = os.path.abspath(out_dir)
                    return [
                        (sub, os.path.join(out_root, os.path.basename(sub)))
                        for sub, _ in jobs
                    ]
                return jobs
        job_out = out_dir if out_dir is not None else cif_path
        return [(cif_path, os.path.abspath(job_out))]

    raise FileNotFoundError(f"No QRS-openffall.cif found in {cif_path}")


def main(cif_path, nproc=4, step=200, fmax=0.1, out_dir=None, db_file=None, ref_code=None, matched_cif=None, cutoff_pct=50.0, e_tol=1e-3, calculator='MACE', model=None, quick=False, max_unique=None, match_max_sites=50, pool=None):
    main_start = time.time()
    input_folder = None
    if os.path.isdir(cif_path):
        input_folder = os.path.abspath(cif_path)
        inferred_cif = os.path.join(input_folder, 'QRS-openffall.cif')
        if not os.path.exists(inferred_cif):
            raise FileNotFoundError(f"No QRS-openffall.cif found in {input_folder}")
        cif_path = inferred_cif
        if matched_cif is None:
            inferred_matched = os.path.join(input_folder, 'QRS-openff-matched.cif')
            if os.path.exists(inferred_matched):
                matched_cif = inferred_matched
                print(f'Auto-detected matched CIF file: {matched_cif}')
            else:
                print(f'No QRS-openff-matched.cif found in {input_folder}; continuing without matched file')
    if out_dir is None:
        if input_folder is not None:
            out_dir = input_folder
        else:
            out_dir = os.path.dirname(os.path.abspath(cif_path)) or '.'
    out_dir = os.path.abspath(out_dir)
    model = resolve_model(calculator, model=model, quick=quick)
    default_model = DEFAULT_MODELS.get(calculator)
    model_tag = f'{calculator}_{model}' if model != default_model else calculator
    print(f'Output directory: {out_dir}')
    print(f'Calculator: {calculator}, model: {model}' + (' (quick preset)' if quick else ''))

    blocks = parse_qrs_cif(cif_path)
    if not blocks:
        print('No blocks found in', cif_path)
        return

    # Extract energies and labels
    entries = []
    for i, (lbl, txt, eng) in enumerate(blocks):
        entries.append({'label': lbl, 'text': txt, 'energy': eng, 'idx': i})

    energies = [e['energy'] for e in entries if e['energy'] is not None]
    if not energies:
        print('No energies parsed')
        return
    # compute percentile cutoff (default 50% median)
    threshold = float(np.percentile(energies, cutoff_pct))
    print(f'Parsed {len(entries)} blocks, {cutoff_pct:.1f} percentile energy = {threshold:.6f} eV')

    # select those with energy < threshold
    selected = [e for e in entries if e['energy'] is not None and e['energy'] <= threshold]
    print(f'{len(selected)} selected (energy <= {cutoff_pct:.1f} percentile)')

    # Use a single StructureMatcher across deduplication and reference matching.
    matcher = StructureMatcher(ltol=0.3, stol=0.3, angle_tol=5.0)
    unique_selected, dup_map = _deduplicate_selected(
        selected, matcher, e_tol=e_tol, match_max_sites=match_max_sites,
    )

    print(f'{len(unique_selected)} unique structures after deduplication (e_tol={e_tol})')

    if max_unique is not None and len(unique_selected) > max_unique:
        print(f'Limiting to first {max_unique} lowest-energy unique structures (from {len(unique_selected)})')
        kept_ids = {e['idx'] for e in unique_selected[:max_unique]}
        unique_selected = unique_selected[:max_unique]
        dup_map = {dup: rep for dup, rep in dup_map.items() if rep in kept_ids}

    os.makedirs(out_dir, exist_ok=True)

    # Build duplicate groups mapping: rep_idx -> [dup_idx,...]
    dup_groups = {}
    for dup_idx, rep_idx in dup_map.items():
        dup_groups.setdefault(rep_idx, []).append(dup_idx)

    # Load reference structure and test original CIF matches first
    ref_pmg = None
    original_match_ids = set()
    relaxed_match_ids = set()
    if ref_code is None:
        ref_code = guess_ref_code_from_cif(cif_path)
        if ref_code is not None:
            print(f'Inferred ref_code={ref_code} from CIF path')
    if matched_cif is not None:
        try:
            matched_labels = get_matched_labels_from_file(matched_cif)
            label_to_idx = {entry['label']: entry['idx'] for entry in entries}
            for label in matched_labels:
                if label in label_to_idx:
                    original_match_ids.add(label_to_idx[label])
            print(f'Loaded {len(original_match_ids)} original matches from {matched_cif}')
        except Exception as e:
            print(f'Failed to load matched_cif {matched_cif}: {e}')
    if db_file is not None and ref_code is not None:
        db = PyDatabase(db_file)
        ref_xtal = db.get_pyxtal(ref_code)
        ref_pmg = ref_xtal.to_pymatgen()
        ref_pmg.remove_species(["H"]) if hasattr(ref_pmg, 'remove_species') else None
        if matched_cif is None:
            for entry in entries:
                atoms_orig = read(StringIO(entry['text']), format='cif')
                pmg_orig = ase2pymatgen(atoms_orig)
                pmg_orig.remove_species(["H"]) if hasattr(pmg_orig, 'remove_species') else None
                if matcher.fit(pmg_orig, ref_pmg):
                    original_match_ids.add(entry['idx'])
            print(f'Loaded reference {ref_code} and found {len(original_match_ids)} original CIF matches')
        else:
            print(f'Loaded reference {ref_code} for relaxed matching only; original matches come from matched_cif')

    # Run relaxations in parallel on the deduplicated list.
    # Distribute low-energy structures evenly across worker chunks, then use
    # a large chunksize so each worker receives a balanced block.
    sorted_tasks = sorted(
        [
            (e['label'], e['text'], e['energy'], step, fmax, e['idx'], calculator, model, quick)
            for e in unique_selected
        ],
        key=lambda x: x[2]
    )
    buckets = [[] for _ in range(nproc)]
    for i, task in enumerate(sorted_tasks):
        buckets[i % nproc].append(task)
    tasks = [task for bucket in buckets for task in bucket]
    chunksize = max(1, len(tasks) // max(nproc, 1))
    workers = min(nproc, len(tasks))
    print(f'Running {len(tasks)} relaxation tasks for {len(unique_selected)} unique structures with calculator={calculator}, model={model}')
    print(f'Using chunksize={chunksize} for {workers} workers')
    if pool is not None:
        results = pool.map(worker_relax, tasks, chunksize=chunksize)
    elif workers <= 1:
        _configure_worker_threads()
        get_calculator(calculator, model=model, quick=quick)
        results = [worker_relax(task) for task in tasks]
    else:
        # spawn avoids macOS fork crashes after MPS/PyTorch init in the parent.
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            processes=workers,
            initializer=_init_worker,
            initargs=(calculator, model, quick),
        ) as job_pool:
            results = job_pool.map(worker_relax, tasks, chunksize=chunksize)
    print(f'Completed {len(results)} relaxation results')
    total_relax_time = sum([r[6] for r in results if len(r) > 6 and r[6] is not None])
    print(f'Total relaxation wall time: {total_relax_time:.1f} s')

    # collect results
    id_all = [i for i in range(len(entries))]
    orig_energy_all = [entries[i]['energy'] for i in id_all]

    selected_ids = []
    model_energies = {}
    time_costs = {}
    # store relaxed cif text per idx for matching
    relaxed_cifs = {}
    for res in results:
        # res: (label, idx, orig, mace, status, relaxed_cif, elapsed_s)
        lbl, idx, orig, mace, status, relaxed_cif, elapsed_s = res
        if mace is not None:
            # store for representative
            selected_ids.append(idx)
            model_energies[idx] = mace
            time_costs[idx] = elapsed_s
            relaxed_cifs[idx] = relaxed_cif
            # propagate to duplicates
            for dup in dup_groups.get(idx, []):
                selected_ids.append(dup)
                model_energies[dup] = mace
                time_costs[dup] = elapsed_s
                relaxed_cifs[dup] = relaxed_cif
        else:
            print(f'Relax failed for {lbl} idx={idx}: {status}, time {elapsed_s:.1f}s')

    # Prepare relaxed match IDs set (may be empty)
    if ref_pmg is not None:
        try:
            for idx, relaxed_cif in relaxed_cifs.items():
                try:
                    atoms_rel = read(StringIO(relaxed_cif), format='cif')
                    pmg_rel = ase2pymatgen(atoms_rel)
                    pmg_rel.remove_species(["H"]) if hasattr(pmg_rel, 'remove_species') else None
                    if matcher.fit(pmg_rel, ref_pmg):
                        relaxed_match_ids.add(idx)
                except Exception:
                    continue
            if relaxed_match_ids:
                print(f'Matched {len(relaxed_match_ids)} relaxed structures to {ref_code}:', sorted(list(relaxed_match_ids)))
        except Exception as e:
            print('Relaxed reference matching failed:', e)

    # Create a 2x1 stacked plot
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Top: ID vs original energy (all entries)
    xs = id_all
    ys = orig_energy_all
    below_x = [i for i in xs if ys[i] is not None and ys[i] < threshold]
    above_x = [i for i in xs if ys[i] is not None and ys[i] >= threshold]
    if below_x:
        ax_top.plot(
            below_x,
            [ys[i] for i in below_x],
            'o-',
            markersize=6,
            label=f'OpenFF (<= {cutoff_pct:.1f}%)',
        )
    if above_x:
        ax_top.scatter(
            above_x,
            [ys[i] for i in above_x],
            marker='o',
            color='grey',
            s=36,
            label=f'OpenFF (> {cutoff_pct:.1f}%)',
        )
    # highlight original CIF matches on top if present
    if original_match_ids:
        match_x_top = sorted(list(original_match_ids))
        match_y_top = [orig_energy_all[i] for i in match_x_top]
        ax_top.scatter(
            match_x_top,
            match_y_top,
            s=160,
            marker='*',
            c='crimson',
            edgecolors='black',
            linewidths=0.8,
            zorder=5,
            label='OpenFF match',
        )
    ax_top.set_ylabel('Energy (eV)')
    ax_top.grid(True, alpha=0.3)
    ax_top.legend(loc=3)
    ys_min = min(ys)
    ys_median = np.median(ys)
    ax_top.set_ylim(ys_min - 1, ys_median + 2.0)

    # Bottom: ID vs selected model energy for selected structures
    sel_x = sorted(selected_ids)
    sel_y = [model_energies[i] for i in sel_x]
    ax_bot.plot(sel_x, sel_y, 's-', color='#008B8B', 
                label=f'{calculator} energy ({len(unique_selected)} unique)', 
                markersize=6, zorder=1)
    # Highlight relaxed structure matches on bottom if present
    if sel_x:
        matched_sel = [i for i in sel_x if i in relaxed_match_ids]  
        if matched_sel:
            ax_bot.scatter(matched_sel, [model_energies[i] for i in matched_sel], s=160, marker='*', 
                           color='crimson', edgecolors='black', 
                           lw=0.8, zorder=5, label=f'{calculator} match')

    ax_bot.set_xlabel('ID')
    ax_bot.set_ylabel(f'{calculator} Energy (eV)')
    ax_bot.grid(True, alpha=0.3)
    ax_bot.legend(loc=2)

    fig.tight_layout()
    out_png = os.path.join(out_dir, f'{model_tag}_energy_{cutoff_pct:.1f}.png')
    fig.savefig(out_png, dpi=200)
    print('Saved plot to', out_png) 

    # Also save CSV of results
    import csv
    csvp = os.path.join(out_dir, f'{model_tag}_energy_{cutoff_pct:.1f}.csv')
    with open(csvp, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['idx','label','orig_energy','mace_energy','time_cost_s','status','original_matched','relaxed_matched'])
        res_map = {r[1]: r for r in results}
        for e in entries:
            idx = e['idx']
            orig_matched = 'yes' if idx in original_match_ids else ''
            relaxed_matched = 'yes' if idx in relaxed_match_ids else ''
            if idx in res_map:
                lbl, ii, orig, mace, status, relaxed_cif, elapsed_s = res_map[idx]
                writer.writerow([ii, lbl, orig, mace if mace is not None else '', f'{elapsed_s:.3f}' if elapsed_s is not None else '', status, orig_matched, relaxed_matched])
            elif idx in model_energies:
                lbl = e['label']
                ii = idx
                orig = e['energy']
                mace = model_energies[idx]
                status = 'propagated'
                elapsed_s = time_costs.get(idx, None)
                writer.writerow([ii, lbl, orig, mace if mace is not None else '', f'{elapsed_s:.3f}' if elapsed_s is not None else '', status, orig_matched, relaxed_matched])
            else:
                writer.writerow([idx, e['label'], e['energy'], '', '', 'not_selected', orig_matched, relaxed_matched])
    print('Saved CSV to', csvp)
    print(f'Total workflow elapsed time: {time.time() - main_start:.1f} s')


def default_nproc():
    """Use SLURM cpus-per-task when available, else local CPU count."""
    for key in ("SLURM_CPUS_PER_TASK", "SLURM_JOB_CPUS_PER_NODE"):
        val = os.environ.get(key)
        if val and val.isdigit() and int(val) > 0:
            return int(val)
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, NotImplementedError):
        return mp.cpu_count() or 1


if __name__ == '__main__':
    # spawn is required on macOS when PyTorch/MPS calculators are used with multiprocessing.
    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass
    
    parser = argparse.ArgumentParser(description='Relax QRS CIF blocks with MACE and compare energies')
    parser.add_argument('cif', help='Path to QRS multi-block CIF, a case folder, or a parent folder of case folders')
    parser.add_argument(
        '--nproc',
        type=int,
        default=None,
        help='Number of parallel workers (default: SLURM_CPUS_PER_TASK or available CPUs)',
    )
    parser.add_argument('--pct', type=float, default=10.0, help='Energy percentile cutoff (e.g. 10 for 10%%)')
    parser.add_argument('--db', default="pyxtal/database/test.db", help='Path to test.db for reference matching')
    parser.add_argument('--matched-cif', default=None, help='Path to file listing matched CIF labels, one per line')
    parser.add_argument('--e-tol', type=float, default=1e-3, help='Energy tolerance (eV) for grouping duplicates')
    parser.add_argument('--max-unique', type=int, default=None, help='Maximum number of unique structures to relax (lowest energy first)')
    parser.add_argument(
        '--match-max-sites',
        type=int,
        default=50,
        help='Compare only the first N sites during deduplication (default: 50; 0 = use full structure)',
    )
    parser.add_argument('--step', type=int, default=200, help='FIRE steps for relaxation')
    parser.add_argument('--fmax', type=float, default=0.1, help='Fmax for FIRE relax')
    parser.add_argument('--calculator', default='MACE', choices=['MACE', 'MACEOFF', 'ORB-V3', 'ORBV3', 'ORB'], help='MLP calculator to use for relaxation')
    parser.add_argument(
        '--model',
        default=None,
        help='Model variant (MACE/MACEOFF: small|medium|large; ORB v3: conservative-inf-omat|direct-20-omat|direct-inf-omat; ORB v2: orb-v2|orb-d3-v2|orb-d3-sm-v2|orb-d3-xs-v2)',
    )
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Use a cheaper/faster model preset for quick testing (ORB: direct-20-omat, MACEOFF: small)',
    )
    parser.add_argument('--out-dir', default=None, help='Output directory (defaults to input folder; if set without QRS-openffall.cif, process all subdirs)')
    args = parser.parse_args()
    if args.nproc is None:
        args.nproc = default_nproc()
    print(f'Using {args.nproc} parallel workers')
    if args.quick and args.model is not None:
        parser.error('Use only one of --quick or --model')

    common_kwargs = dict(
        nproc=args.nproc,
        step=args.step,
        fmax=args.fmax,
        db_file=args.db,
        ref_code=None,
        matched_cif=args.matched_cif,
        cutoff_pct=args.pct,
        e_tol=args.e_tol,
        max_unique=args.max_unique,
        match_max_sites=args.match_max_sites,
        calculator=args.calculator,
        quick=args.quick,
    )

    jobs = collect_jobs(args.cif, args.out_dir)
    if not jobs:
        parser.error(f'No QRS-openffall.cif found under {args.cif}')

    model = resolve_model(args.calculator, model=args.model, quick=args.quick)
    _configure_worker_threads()
    t0 = time.perf_counter()
    cache_path = ensure_mace_model_cached(args.calculator, model=model, quick=args.quick)
    cache_elapsed = time.perf_counter() - t0
    if cache_path:
        print(f'Model cache ready in {cache_elapsed:.2f}s: {cache_path}')
    else:
        print(f'Model cache check finished in {cache_elapsed:.2f}s')

    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=args.nproc,
        initializer=_init_worker,
        initargs=(args.calculator, model, args.quick),
    ) as pool:
        t0 = time.perf_counter()
        pool.map(_worker_ping, range(args.nproc), chunksize=1)
        pool_elapsed = time.perf_counter() - t0
        print(
            f'Worker pool ready in {pool_elapsed:.1f}s '
            f'({args.nproc} workers; model loaded once per worker, reused across all jobs)'
        )
        for i, (job_dir, job_out_dir) in enumerate(jobs, start=1):
            if len(jobs) > 1:
                print(f'\n=== [{i}/{len(jobs)}] Processing {job_dir} ===')
            main(
                job_dir,
                out_dir=job_out_dir,
                pool=pool,
                model=model,
                **common_kwargs,
            )
