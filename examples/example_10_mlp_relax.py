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

# Suppress specific noisy warnings from ASE CIF parser and torch_dftd
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings(
    "ignore",
    message="crystal system 'monoclinic' is not interpreted",
    category=UserWarning,
)
warnings.filterwarnings(
    "ignore",
    message="Creating a tensor from a list of numpy.ndarrays is extremely slow",
    category=UserWarning,
)


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


def _structure_max_dist(matcher, pmg1, pmg2):
    """Return StructureMatcher max site distance, or None if unmatched."""
    try:
        rmsd = matcher.get_rms_dist(pmg1, pmg2)
    except Exception:
        return None
    return rmsd[1] if rmsd is not None else None


def _rms_match(matcher, pmg1, pmg2, max_rmsd=0.3):
    """True only if StructureMatcher max site distance is strictly below ``max_rmsd``.

    Uses ``get_rms_dist`` → ``(rms, max_dist)``. A pair with ``max_dist >= max_rmsd``
    (or no mapping) is treated as not a match — same criterion as QRS.
    """
    max_dist = _structure_max_dist(matcher, pmg1, pmg2)
    return max_dist is not None and max_dist < max_rmsd


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


def _deduplicate_selected(
    selected,
    matcher,
    e_tol=1e-3,
    progress_every=25,
    match_max_sites=50,
    max_rmsd=0.3,
):
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
    print(f'Deduplication max_rmsd cutoff: {max_rmsd} (not a match if max site dist >= {max_rmsd})')

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
            if _rms_match(
                matcher,
                _pmg_for_matching(pmg1, match_max_sites),
                _pmg_for_matching(pmg2, match_max_sites),
                max_rmsd=max_rmsd,
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


def _entry_is_structural_dup(ent, reps, matcher, pmg_cache, match_max_sites, max_rmsd):
    """Return rep idx if ``ent`` matches any representative, else None."""
    try:
        pmg1 = _get_pmg_no_h(ent, pmg_cache)
    except Exception:
        return None
    pmg1_m = _pmg_for_matching(pmg1, match_max_sites)
    for rep in reps:
        try:
            pmg2 = _get_pmg_no_h(rep, pmg_cache)
        except Exception:
            continue
        if not _structures_likely_match(pmg1, pmg2, max_sites=match_max_sites):
            continue
        if _rms_match(
            matcher,
            pmg1_m,
            _pmg_for_matching(pmg2, match_max_sites),
            max_rmsd=max_rmsd,
        ):
            return rep['idx']
    return None


def _expand_unique_from_higher_energy(
    unique_selected,
    all_entries,
    selection_threshold,
    matcher,
    max_unique,
    match_max_sites=50,
    max_rmsd=0.3,
    progress_every=25,
    max_higher_scan=2000,
):
    """Scan higher-energy structures when the low-energy pool dedups to too few.

    ``max_higher_scan`` caps how many above-threshold candidates are examined
    (energy-sorted, lowest first). Use ``0`` to disable expansion entirely.
    """
    if max_unique is None or len(unique_selected) >= max_unique:
        return unique_selected, {}

    if selection_threshold is None:
        return unique_selected, {}

    if max_higher_scan is not None and int(max_higher_scan) <= 0:
        print(
            f'Only {len(unique_selected)} unique in low-energy pool; '
            f'skipping higher-energy expansion (max_higher_scan={max_higher_scan})',
            flush=True,
        )
        return unique_selected, {}

    higher = sorted(
        (
            e for e in all_entries
            if e.get('energy') is not None and e['energy'] > selection_threshold
        ),
        key=lambda x: x['energy'],
    )
    if not higher:
        print(
            f'Only {len(unique_selected)} unique in low-energy pool; '
            f'no higher-energy structures to scan',
            flush=True,
        )
        return unique_selected, {}

    n_available = len(higher)
    if max_higher_scan is not None:
        higher = higher[: int(max_higher_scan)]

    print(
        f'Only {len(unique_selected)} unique in low-energy pool; '
        f'scanning {len(higher)}/{n_available} higher-energy structures '
        f'(>{selection_threshold:.6f} eV) for up to {max_unique} total unique'
        + (
            f' [capped by max_higher_scan={max_higher_scan}]'
            if len(higher) < n_available
            else ''
        ),
        flush=True,
    )

    kept = list(unique_selected)
    pmg_cache = {}
    extra_dup_map = {}
    scan_start = time.perf_counter()

    for i, ent in enumerate(higher, start=1):
        if len(kept) >= max_unique:
            break
        if progress_every and i % progress_every == 0:
            elapsed = time.perf_counter() - scan_start
            print(
                f'Higher-energy scan [{i}/{len(higher)}]: '
                f'{len(kept)} unique so far, elapsed {elapsed:.1f}s',
                flush=True,
            )

        rep_idx = _entry_is_structural_dup(
            ent, kept, matcher, pmg_cache, match_max_sites, max_rmsd,
        )
        if rep_idx is not None:
            extra_dup_map[ent['idx']] = rep_idx
        else:
            kept.append(ent)

    kept.sort(key=lambda x: x['energy'])
    if len(kept) > max_unique:
        kept = kept[:max_unique]

    n_added = len(kept) - len(unique_selected)
    elapsed = time.perf_counter() - scan_start
    print(
        f'Higher-energy expansion: +{n_added} unique -> {len(kept)} total '
        f'({len(higher)} scanned / {n_available} available, {elapsed:.1f}s)',
        flush=True,
    )
    if len(kept) < max_unique and len(higher) < n_available:
        print(
            f'Note: stopped at {len(kept)}/{max_unique} unique due to '
            f'max_higher_scan={max_higher_scan}; raise it to search further',
            flush=True,
        )
    return kept, extra_dup_map


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
    args: (label, block_text, orig_energy, step, fmax, label_idx, calculator, model, quick, relax_timeout_min)
    Returns (label, label_idx, original_energy, model_energy or None, status, relaxed_cif, elapsed_s)
    """
    (
        label,
        block_text,
        orig_energy,
        step,
        fmax,
        label_idx,
        calculator,
        model,
        quick,
        relax_timeout_min,
    ) = args
    start_time = time.time()
    short_label = label[-45:] if len(label) > 45 else label
    print(f'Relax start idx={label_idx:5d} pid={os.getpid()} {short_label}', flush=True)
    try:
        atoms = read(StringIO(block_text), format='cif')
    except Exception as e:
        return (label, label_idx, orig_energy, None, f"read_error: {e}", None, time.time() - start_time)

    try:
        relaxed = ASE_relax(
            atoms,
            calculator=calculator,
            model=model,
            quick=quick,
            opt_lat=False,
            step=step,
            fmax=fmax,
            max_time=relax_timeout_min,
            label=label,
            logfile=None,
        )
        elapsed = time.time() - start_time
        if relaxed is None:
            return (label, label_idx, orig_energy, None, 'relax_failed', None, elapsed)
        model_energy = relaxed.get_potential_energy()
        print(
            f'Relax done  idx={label_idx:5d}, {calculator} {model_energy:.3f}, time {elapsed:.1f}s {short_label}',
            flush=True,
        )
        sio = BytesIO()
        write(sio, relaxed, format='cif')
        relaxed_cif = sio.getvalue().decode('utf-8')
        return (label, label_idx, orig_energy, model_energy, 'ok', relaxed_cif, elapsed)
    except Exception as e:
        return (label, label_idx, orig_energy, None, f'relax_exception: {e}', None, time.time() - start_time)


def _log_relax_result(completed, total, result, batch_start):
    lbl, idx, _orig, mace, status, _cif, elapsed = result
    short = lbl[-45:] if len(lbl) > 45 else lbl
    detail = f'{mace:.3f} eV' if mace is not None else status
    elapsed_s = f'{elapsed:.1f}s' if elapsed is not None else 'n/a'
    print(
        f'Relax [{completed}/{total}] idx={idx:5d} {detail} '
        f'task={elapsed_s} batch={time.perf_counter() - batch_start:.1f}s {short}',
        flush=True,
    )


def _relax_tasks(
    tasks,
    *,
    pool=None,
    workers=1,
    chunksize=1,
    relax_timeout_min=15.0,
    stall_timeout_min=20.0,
    heartbeat_s=60.0,
):
    """Run relaxations with per-task logging and stall detection."""
    total = len(tasks)
    if total == 0:
        return [], False

    batch_start = time.perf_counter()
    last_done = batch_start
    last_heartbeat = batch_start
    completed = 0
    results = []
    stall_limit = stall_timeout_min * 60.0
    pool_needs_restart = False

    print(
        f'Relax batch: {total} tasks, relax_timeout={relax_timeout_min:.1f} min, '
        f'stall_timeout={stall_timeout_min:.1f} min',
        flush=True,
    )

    if pool is None or workers <= 1:
        for task in tasks:
            result = worker_relax(task)
            completed += 1
            last_done = time.perf_counter()
            _log_relax_result(completed, total, result, batch_start)
            results.append(result)
        return results, False

    async_map = {}
    for task in tasks:
        ar = pool.apply_async(worker_relax, (task,))
        async_map[ar] = task
    pending = set(async_map.keys())

    while pending:
        made_progress = False
        for ar in list(pending):
            if not ar.ready():
                continue
            result = ar.get()
            pending.remove(ar)
            completed += 1
            last_done = time.perf_counter()
            _log_relax_result(completed, total, result, batch_start)
            results.append(result)
            made_progress = True

        if not pending:
            break

        if made_progress:
            continue

        now = time.perf_counter()
        stall = now - last_done
        if stall >= heartbeat_s and now - last_heartbeat >= heartbeat_s:
            print(
                f'Relax [{completed}/{total}] waiting... no finish for {stall:.0f}s '
                f'(stall limit {stall_limit:.0f}s)',
                flush=True,
            )
            last_heartbeat = now

        if stall < stall_limit:
            time.sleep(min(5.0, max(1.0, heartbeat_s - (now - last_heartbeat))))
            continue

        ar = next(iter(pending))
        task = async_map[ar]
        label, _, orig_energy, _, _, label_idx, *_ = task
        pending.remove(ar)
        completed += 1
        last_done = time.perf_counter()
        pool_needs_restart = True
        skip_result = (label, label_idx, orig_energy, None, 'stall_skipped', None, None)
        results.append(skip_result)
        print(
            f'Relax stall: skipping stuck task idx={label_idx} '
            f'({len(pending)} still running in pool)',
            flush=True,
        )
        _log_relax_result(completed, total, skip_result, batch_start)

    if pool_needs_restart:
        print(
            f'WARNING: Skipped {sum(1 for r in results if r[4] == "stall_skipped")} stuck tasks; '
            f'worker pool will be recreated for the next job.',
            flush=True,
        )

    return results, pool_needs_restart


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


def _robust_energy_ylim(energies, low_pad=1.0, high_pct=99.0, high_pad=2.0):
    """
    Y-limits that keep the bulk of energies visible while clipping rare
    very-high-energy outliers (e.g. failed / high-strain MLP relaxations).
    """
    arr = np.asarray([e for e in energies if e is not None], dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    ymin = float(np.min(arr))
    yhi = float(np.percentile(arr, high_pct))
    if yhi <= ymin:
        yhi = ymin + 1.0
    return ymin - float(low_pad), yhi + float(high_pad)


def _result_paths(out_dir, model_tag, cutoff_pct):
    base = f'{model_tag}_energy_{cutoff_pct:.1f}'
    matched = f'{model_tag}_matched_{cutoff_pct:.1f}.cif'
    return (
        os.path.join(out_dir, f'{base}.png'),
        os.path.join(out_dir, f'{base}.csv'),
        os.path.join(out_dir, matched),
    )


def _matched_data_header(orig_label, calculator, rms, max_msd):
    """Build a QRS-style matched ``data_`` header (without the ``data_`` prefix)."""
    base = orig_label[5:] if orig_label.startswith('data_') else orig_label
    # Drop a previous match suffix (e.g. -QRandom-0.09-0.16) if present.
    first = base.find('-QRandom-')
    if first != -1:
        second = base.find('-QRandom-', first + 1)
        if second != -1:
            base = base[:second]
    return f'{base}-{calculator}-{rms:.2f}-{max_msd:.2f}'


def _format_matched_cif_block(ase_cif_text, header, energy):
    """
    Rewrite an ASE CIF block with a QRS-style header / energy / #END.

    Header lines look like::

        data_xafpay03-g7-p36-QRandom-d1.335-spg14-e-47.008-MACEOFF-0.09-0.16
        #Energy: -210836.124246 eV/cell
    """
    from pyxtal.constants import logo

    lines = []
    body_started = False
    for ln in ase_cif_text.splitlines(keepends=True):
        if not body_started:
            if ln.startswith('data_'):
                body_started = True
                lines.append(f'data_{header}\n')
                lines.append(f'#Energy: {energy} eV/cell\n')
                continue
            # Drop ASE preamble before the first data_ block.
            continue
        if ln.strip() == '#END':
            continue
        lines.append(ln)
    if not body_started:
        lines.append(f'data_{header}\n')
        lines.append(f'#Energy: {energy} eV/cell\n')
        lines.append(ase_cif_text)
    text = ''.join(lines)
    if not text.endswith('\n'):
        text += '\n'
    text += '#END\n\n'
    return logo + text


def write_mlp_matched_cif(
    path,
    *,
    entries,
    relaxed_cifs,
    model_energies,
    match_rms,
    calculator,
    dup_map=None,
):
    """
    Write all unique MLP-relaxed matches to ``path`` with QRS-style headers.

    ``match_rms`` maps idx -> (rms, max_msd). Duplicate geometries (propagated
    from a representative) are written once.
    """
    dup_map = dup_map or {}
    idxs = []
    for idx in sorted(match_rms):
        rep = dup_map.get(idx, idx)
        if rep != idx and rep in match_rms:
            continue
        if idx not in relaxed_cifs or idx not in model_energies:
            continue
        idxs.append(idx)

    if not idxs:
        return 0

    label_by_idx = {e['idx']: e['label'] for e in entries}
    with open(path, 'w') as fh:
        for idx in idxs:
            rms, max_msd = match_rms[idx]
            header = _matched_data_header(
                label_by_idx.get(idx, f'idx{idx}'),
                calculator,
                rms,
                max_msd,
            )
            block = _format_matched_cif_block(
                relaxed_cifs[idx], header, model_energies[idx]
            )
            fh.write(block)
    return len(idxs)


def _select_for_relaxation(entries, cutoff_pct, min_selected=1000):
    """Select lowest-energy entries, starting at cutoff_pct and raising if needed."""
    valid = [e for e in entries if e['energy'] is not None]
    if not valid:
        return [], None, cutoff_pct

    energies = [e['energy'] for e in valid]
    initial_threshold = float(np.percentile(energies, cutoff_pct))
    selected = [e for e in valid if e['energy'] <= initial_threshold]

    if len(selected) >= min_selected or len(valid) <= len(selected):
        return selected, initial_threshold, cutoff_pct

    target = min(min_selected, len(valid))
    selection_threshold = sorted(valid, key=lambda x: x['energy'])[target - 1]['energy']
    selected = [e for e in valid if e['energy'] <= selection_threshold]
    effective_pct = 100.0 * len(selected) / len(valid)
    print(
        f'Only {len([e for e in valid if e["energy"] <= initial_threshold])} structures '
        f'at or below {cutoff_pct:.1f} percentile; '
        f'raising cutoff to {effective_pct:.1f} percentile '
        f'({len(selected)} structures, threshold {selection_threshold:.6f} eV)',
        flush=True,
    )
    return selected, selection_threshold, effective_pct


def main(cif_path, nproc=4, step=200, fmax=0.1, out_dir=None, db_file=None, ref_code=None, matched_cif=None, cutoff_pct=50.0, e_tol=1e-3, calculator='MACE', model=None, quick=False, max_unique=None, max_higher_scan=2000, match_max_sites=50, dedup_max_rmsd=0.3, ref_max_rmsd=0.6, pool=None, relax_timeout_min=15.0, stall_timeout_min=20.0, force=False):
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
    out_png, out_csv, out_matched_cif = _result_paths(out_dir, model_tag, cutoff_pct)
    if not force and os.path.isfile(out_png) and os.path.isfile(out_csv):
        case = os.path.basename(out_dir.rstrip(os.sep))
        print(
            f'Skipping {case}: {os.path.basename(out_png)} and {os.path.basename(out_csv)} already exist',
            flush=True,
        )
        return False
    print(f'Output directory: {out_dir}')
    print(f'Calculator: {calculator}, model: {model}' + (' (quick preset)' if quick else ''))
    print(
        f'RMSD cutoffs: dedup={dedup_max_rmsd}, ref-match={ref_max_rmsd}',
        flush=True,
    )

    blocks = parse_qrs_cif(cif_path)
    if not blocks:
        print('No blocks found in', cif_path)
        return False

    # Extract energies and labels
    entries = []
    for i, (lbl, txt, eng) in enumerate(blocks):
        entries.append({'label': lbl, 'text': txt, 'energy': eng, 'idx': i})

    energies = [e['energy'] for e in entries if e['energy'] is not None]
    if not energies:
        print('No energies parsed')
        return False
    print(
        f'Parsed {len(entries)} blocks, {cutoff_pct:.1f} percentile energy = '
        f'{float(np.percentile(energies, cutoff_pct)):.6f} eV',
        flush=True,
    )
    selected, threshold, plot_pct = _select_for_relaxation(entries, cutoff_pct)
    if not selected:
        print('No energies parsed')
        return False
    print(f'{len(selected)} structures selected for relaxation', flush=True)

    # Use a single StructureMatcher across deduplication and reference matching.
    matcher = StructureMatcher(ltol=0.3, stol=0.3, angle_tol=5.0)
    unique_selected, dup_map = _deduplicate_selected(
        selected, matcher, e_tol=e_tol, match_max_sites=match_max_sites,
        max_rmsd=dedup_max_rmsd,
    )

    print(f'{len(unique_selected)} unique structures after deduplication (e_tol={e_tol})')

    if max_unique is not None and len(unique_selected) < max_unique:
        unique_selected, extra_dup = _expand_unique_from_higher_energy(
            unique_selected,
            entries,
            threshold,
            matcher,
            max_unique,
            match_max_sites=match_max_sites,
            max_rmsd=dedup_max_rmsd,
            max_higher_scan=max_higher_scan,
        )
        dup_map.update(extra_dup)

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
            # Only scan the energy-selected pool — never the full CIF (can be 10^4–10^5).
            scan_entries = selected
            print(
                f'Matching {len(scan_entries)} selected OpenFF structures to '
                f'{ref_code} (skipping {len(entries) - len(scan_entries)} above cutoff)',
                flush=True,
            )
            for entry in scan_entries:
                atoms_orig = read(StringIO(entry['text']), format='cif')
                pmg_orig = ase2pymatgen(atoms_orig)
                pmg_orig.remove_species(["H"]) if hasattr(pmg_orig, 'remove_species') else None
                if _rms_match(matcher, pmg_orig, ref_pmg, max_rmsd=ref_max_rmsd):
                    original_match_ids.add(entry['idx'])
            print(f'Loaded reference {ref_code} and found {len(original_match_ids)} original CIF matches')
        else:
            print(f'Loaded reference {ref_code} for relaxed matching only; original matches come from matched_cif')

    # Run relaxations in parallel on the deduplicated list.
    # Distribute low-energy structures evenly across worker chunks, then use
    # a large chunksize so each worker receives a balanced block.
    sorted_tasks = sorted(
        [
            (
                e['label'],
                e['text'],
                e['energy'],
                step,
                fmax,
                e['idx'],
                calculator,
                model,
                quick,
                relax_timeout_min,
            )
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
    pool_terminated = False
    if pool is not None:
        results, pool_terminated = _relax_tasks(
            tasks,
            pool=pool,
            workers=workers,
            chunksize=chunksize,
            relax_timeout_min=relax_timeout_min,
            stall_timeout_min=stall_timeout_min,
        )
    elif workers <= 1:
        _configure_worker_threads()
        get_calculator(calculator, model=model, quick=quick)
        results, pool_terminated = _relax_tasks(
            tasks,
            pool=None,
            workers=1,
            relax_timeout_min=relax_timeout_min,
            stall_timeout_min=stall_timeout_min,
        )
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            processes=workers,
            initializer=_init_worker,
            initargs=(calculator, model, quick),
        ) as job_pool:
            results, pool_terminated = _relax_tasks(
                tasks,
                pool=job_pool,
                workers=workers,
                chunksize=chunksize,
                relax_timeout_min=relax_timeout_min,
                stall_timeout_min=stall_timeout_min,
            )
    print(f'Completed {len(results)} relaxation results')
    total_relax_time = sum([r[6] for r in results if len(r) > 6 and r[6] is not None])
    print(f'Total relaxation wall time: {total_relax_time:.1f} s')

    task_idxs = {task[5] for task in tasks}
    completed_idxs = {r[1] for r in results}
    stall_skipped_idxs = task_idxs - completed_idxs
    if stall_skipped_idxs:
        print(
            f'Skipped {len(stall_skipped_idxs)} relaxation tasks due to stall timeout: '
            f'{sorted(stall_skipped_idxs)}',
            flush=True,
        )

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
            if status != 'stall_skipped':
                elapsed_str = f'{elapsed_s:.1f}s' if elapsed_s is not None else 'n/a'
                print(f'Relax failed for {lbl} idx={idx}: {status}, time {elapsed_str}')

    # Prepare relaxed match IDs set (may be empty)
    relaxed_match_rms = {}  # idx -> (rms, max_msd)
    if ref_pmg is not None:
        try:
            for idx, relaxed_cif in relaxed_cifs.items():
                try:
                    atoms_rel = read(StringIO(relaxed_cif), format='cif')
                    pmg_rel = ase2pymatgen(atoms_rel)
                    if hasattr(pmg_rel, 'remove_species'):
                        pmg_rel.remove_species(["H"])
                    try:
                        rmsd = matcher.get_rms_dist(pmg_rel, ref_pmg)
                    except Exception:
                        rmsd = None
                    if rmsd is not None and rmsd[1] < ref_max_rmsd:
                        relaxed_match_ids.add(idx)
                        relaxed_match_rms[idx] = (float(rmsd[0]), float(rmsd[1]))
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
    below_x = [i for i in xs if ys[i] is not None and ys[i] <= threshold]
    above_x = [i for i in xs if ys[i] is not None and ys[i] > threshold]
    if below_x:
        ax_top.plot(
            below_x,
            [ys[i] for i in below_x],
            'o-',
            markersize=6,
            label=f'OpenFF (<= {plot_pct:.1f}%)',
        )
    if above_x:
        ax_top.scatter(
            above_x,
            [ys[i] for i in above_x],
            marker='o',
            color='grey',
            s=36,
            label=f'OpenFF (> {plot_pct:.1f}%)',
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
    bot_ylim = _robust_energy_ylim(sel_y)
    if bot_ylim is not None:
        ax_bot.set_ylim(*bot_ylim)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    print('Saved plot to', out_png) 

    # Also save CSV of results
    import csv
    with open(out_csv, 'w', newline='') as f:
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
            elif idx in stall_skipped_idxs:
                writer.writerow([idx, e['label'], e['energy'], '', '', 'stall_skipped', orig_matched, relaxed_matched])
            else:
                writer.writerow([idx, e['label'], e['energy'], '', '', 'not_selected', orig_matched, relaxed_matched])
    print('Saved CSV to', out_csv)

    if relaxed_match_rms:
        n_written = write_mlp_matched_cif(
            out_matched_cif,
            entries=entries,
            relaxed_cifs=relaxed_cifs,
            model_energies=model_energies,
            match_rms=relaxed_match_rms,
            calculator=calculator,
            dup_map=dup_map,
        )
        print(f'Saved {n_written} matched structure(s) to {out_matched_cif}')
    elif os.path.isfile(out_matched_cif):
        # Avoid stale matches from a previous run with different settings.
        os.remove(out_matched_cif)
        print(f'Removed stale matched CIF: {out_matched_cif}')

    print(f'Total workflow elapsed time: {time.time() - main_start:.1f} s')
    return pool_terminated


def _cpu_limit():
    """Use SLURM cpus-per-task when available, else local CPU count."""
    for key in ("SLURM_CPUS_PER_TASK", "SLURM_JOB_CPUS_PER_NODE"):
        val = os.environ.get(key)
        if val and val.isdigit() and int(val) > 0:
            return int(val)
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, NotImplementedError):
        return mp.cpu_count() or 1


def _slurm_total_mem_mb():
    """Return allocated node memory in MB from SLURM, if available."""
    node_mem = os.environ.get("SLURM_MEM_PER_NODE")
    if node_mem and node_mem.isdigit() and int(node_mem) > 0:
        return int(node_mem)

    cpu_mem = os.environ.get("SLURM_MEM_PER_CPU")
    cpus = os.environ.get("SLURM_CPUS_PER_TASK") or os.environ.get("SLURM_JOB_CPUS_PER_NODE")
    if cpu_mem and cpu_mem.isdigit() and cpus and cpus.isdigit():
        total = int(cpu_mem) * int(cpus)
        if total > 0:
            return total
    return None


def resolve_nproc(requested=None, mem_per_worker_gb=4.0):
    """Resolve worker count from CPU limit and SLURM memory allocation."""
    cpu_limit = _cpu_limit()
    nproc = cpu_limit if requested is None else min(requested, cpu_limit)

    mem_mb = _slurm_total_mem_mb()
    if mem_mb is not None and mem_per_worker_gb > 0:
        mem_per_worker_mb = int(mem_per_worker_gb * 1024)
        mem_limit = max(1, mem_mb // mem_per_worker_mb)
        if nproc > mem_limit:
            print(
                f'Capping workers from {nproc} to {mem_limit} based on '
                f'SLURM memory ({mem_mb} MB, {mem_per_worker_gb:g} GB/worker)',
                flush=True,
            )
            nproc = mem_limit
    return nproc


def default_nproc(mem_per_worker_gb=4.0):
    """Use SLURM cpus and memory when available, else local CPU count."""
    return resolve_nproc(None, mem_per_worker_gb=mem_per_worker_gb)


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
        help='Number of parallel workers (default: SLURM_CPUS_PER_TASK, capped by SLURM memory)',
    )
    parser.add_argument(
        '--mem-per-worker-gb',
        type=float,
        default=4.0,
        help='Estimated peak RAM per worker in GB for SLURM memory capping (default: 4.0)',
    )
    parser.add_argument('--pct', type=float, default=5.0, help='Energy percentile cutoff (e.g. 10 for 10%%)')
    parser.add_argument('--db', default="pyxtal/database/test.db", help='Path to test.db for reference matching')
    parser.add_argument('--matched-cif', default=None, help='Path to file listing matched CIF labels, one per line')
    parser.add_argument('--e-tol', type=float, default=1e-3, help='Energy tolerance (eV) for grouping duplicates')
    parser.add_argument(
        '--max-unique',
        type=int,
        default=None,
        help=(
            'Target number of unique structures to relax (lowest energy first). '
            'If the low-energy pool dedups to fewer, scans higher-energy '
            'structures until this count is reached (see --max-higher-scan).'
        ),
    )
    parser.add_argument(
        '--max-higher-scan',
        type=int,
        default=2000,
        help=(
            'Max higher-energy candidates to scan when filling --max-unique '
            '(default: 2000; 0 disables higher-energy expansion)'
        ),
    )
    parser.add_argument(
        '--match-max-sites',
        type=int,
        default=50,
        help='Compare only the first N sites during deduplication (default: 50; 0 = use full structure)',
    )
    parser.add_argument(
        '--dedup-max-rmsd',
        type=float,
        default=0.3,
        help=(
            'Max site distance for deduplication / higher-energy uniqueness '
            '(default: 0.3; not a duplicate if max_dist >= this value)'
        ),
    )
    parser.add_argument(
        '--ref-max-rmsd',
        type=float,
        default=0.6,
        help=(
            'Max site distance for matching OpenFF/MLP structures to the '
            'experimental reference (default: 0.6; not a match if max_dist >= this)'
        ),
    )
    parser.add_argument(
        '--max-rmsd',
        type=float,
        default=None,
        help='Deprecated alias for --ref-max-rmsd',
    )
    parser.add_argument('--step', type=int, default=300, help='FIRE steps for relaxation')
    parser.add_argument('--fmax', type=float, default=0.05, help='Fmax for FIRE relax')
    parser.add_argument(
        '--relax-timeout',
        type=float,
        default=15.0,
        help='Max wall time per structure relaxation in minutes (default: 15)',
    )
    parser.add_argument(
        '--stall-timeout',
        type=float,
        default=20.0,
        help='Abort batch if no relaxation finishes within this many minutes (default: 20)',
    )
    parser.add_argument('--calculator', default='MACEOFF', choices=['MACE', 'MACEOFF', 'ORB-V3', 'ORBV3', 'ORB'], help='MLP calculator to use for relaxation')
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
    parser.add_argument(
        '--force',
        action='store_true',
        help='Recompute even if the output PNG and CSV already exist',
    )
    parser.add_argument('--out-dir', default=None, help='Output directory (defaults to input folder; if set without QRS-openffall.cif, process all subdirs)')
    args = parser.parse_args()
    args.nproc = resolve_nproc(args.nproc, mem_per_worker_gb=args.mem_per_worker_gb)
    print(f'Using {args.nproc} parallel workers')
    if args.quick and args.model is not None:
        parser.error('Use only one of --quick or --model')

    ref_max_rmsd = args.ref_max_rmsd
    if args.max_rmsd is not None:
        print(
            f'Warning: --max-rmsd is deprecated; using it as --ref-max-rmsd={args.max_rmsd}',
            flush=True,
        )
        ref_max_rmsd = args.max_rmsd

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
        max_higher_scan=args.max_higher_scan,
        match_max_sites=args.match_max_sites,
        dedup_max_rmsd=args.dedup_max_rmsd,
        ref_max_rmsd=ref_max_rmsd,
        calculator=args.calculator,
        quick=args.quick,
        relax_timeout_min=args.relax_timeout,
        stall_timeout_min=args.stall_timeout,
        force=args.force,
    )

    def make_worker_pool():
        worker_pool = ctx.Pool(
            processes=args.nproc,
            initializer=_init_worker,
            initargs=(args.calculator, model, args.quick),
        )
        t0 = time.perf_counter()
        worker_pool.map(_worker_ping, range(args.nproc), chunksize=1)
        pool_elapsed = time.perf_counter() - t0
        print(
            f'Worker pool ready in {pool_elapsed:.1f}s '
            f'({args.nproc} workers; model loaded once per worker, reused across all jobs)'
        )
        return worker_pool

    ctx = mp.get_context("spawn")
    pool = make_worker_pool()
    try:
        for i, (job_dir, job_out_dir) in enumerate(jobs, start=1):
            if len(jobs) > 1:
                print(f'\n=== [{i}/{len(jobs)}] Processing {job_dir} ===')
            try:
                pool_terminated = main(
                    job_dir,
                    out_dir=job_out_dir,
                    pool=pool,
                    model=model,
                    **common_kwargs,
                )
            except RuntimeError as exc:
                if 'stalled' not in str(exc).lower():
                    raise
                print(f'WARNING: {exc}', flush=True)
                pool_terminated = True
            if pool_terminated:
                print('Recreating worker pool after stall timeout.', flush=True)
                try:
                    pool.terminate()
                    pool.join()
                except Exception:
                    pass
                pool = make_worker_pool()
    finally:
        try:
            pool.terminate()
            pool.join()
        except Exception:
            pass
