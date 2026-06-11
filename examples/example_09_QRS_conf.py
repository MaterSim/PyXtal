"""
Example: Run QRS with pregenerated molecular conformers.

This script iterates over entries in pyxtal/database/test.db, precomputes a pool
of pyxtal_molecule conformers via generate_molecules, and passes them to QRS.
With the lattice fixed and conformers pregenerated, QRS samples only Wyckoff-site
fractional coordinates and molecular orientations for each trial crystal.

Results are appended to Tests/qrs_results.csv.
"""

import csv
import os
from time import perf_counter

import matplotlib.pyplot as plt

from pyxtal.constants import single_smiles
from pyxtal.db import database
from pyxtal.molecule import generate_molecules, pyxtal_molecule
from pyxtal.optimize import QRS


def get_vector_dimension(xtal):
    """Return the total dimension of sampling vectors for this crystal."""
    return sum(len(site.get_bounds()) for site in xtal.mol_sites)


def build_sites_from_reference(ref_xtal):
    """Build QRS site labels grouped by molecule type from a reference crystal."""
    sites = [[] for _ in range(len(ref_xtal.numMols))]
    for site in ref_xtal.mol_sites:
        sites[site.type].append(site.wp.get_label())
    return sites


def filter_similar_molecules(mols, rmsd_tol=0.5):
    """Remove near-duplicate pyxtal_molecule conformers using pairwise RMSD."""
    unique_mols = []
    for mol in mols:
        is_duplicate = False
        for ref_mol in unique_mols:
            rmsd, _ = mol.get_rmsd2(mol.mol.cart_coords, ref_mol.mol.cart_coords)
            if rmsd < rmsd_tol:
                is_duplicate = True
                break
        if not is_duplicate:
            unique_mols.append(mol)
    return unique_mols


def get_match_points(qrs):
    """Map QRS match (gen, pop) records to visited-structure IDs and energies."""
    if not hasattr(qrs, "matches") or not hasattr(qrs, "stats"):
        return [], []

    match_keys = set()
    for match in qrs.matches:
        if len(match) >= 2:
            match_keys.add((int(match[0]), int(match[1])))

    if not match_keys:
        return [], []

    match_ids = []
    match_energies = []
    visited_id = 0
    for gen in range(qrs.N_gen):
        for pop in range(qrs.N_pop):
            energy = float(qrs.stats[gen, pop, 0])
            if energy < qrs.E_max:
                visited_id += 1
                if (gen, pop) in match_keys:
                    match_ids.append(visited_id)
                    match_energies.append(energy)

    return match_ids, match_energies


def get_visited_energies(qrs):
    """Return visited energies in the same (gen, pop) order used for match IDs."""
    if not hasattr(qrs, "stats"):
        return []

    energies = []
    for gen in range(qrs.N_gen):
        for pop in range(qrs.N_pop):
            energy = float(qrs.stats[gen, pop, 0])
            if energy < qrs.E_max:
                energies.append(energy)
    return energies


def select_delta_angle(molecules, composition=None):
    """Return delta_angle(s) based on molecule sizes.

    If `composition` is provided and has multiple components, return a list
    of delta angles (one per component). Otherwise return a single float.

    Rules (per component):
      - max_atoms <= 3 -> 90.0
      - max_atoms < 10 -> 60.0
      - max_atoms < 20 -> 45.0
      - max_atoms < 30 -> 30.0
      - max_atoms < 40 -> 20.0
      - else           -> 15.0
    """
    # Default fallback
    DEFAULT = 15.0
    if not molecules:
        if composition is None or len(composition) <= 1:
            return DEFAULT
        return [DEFAULT for _ in composition]

    # If composition provided and multi-component, compute per-component
    if composition is not None and len(composition) > 1:
        delta_list = []
        for idx, count in enumerate(composition):
            # guard when molecules list does not align with composition
            pool = molecules[idx] if idx < len(molecules) else None
            if not pool:
                delta_list.append(DEFAULT)
                continue
            mol = pool[0]
            try:
                atom_count = len(mol.mol)
            except Exception:
                atom_count = len(getattr(mol, "atoms", []))
            if atom_count <= 3:
                delta_list.append(90.0)
            elif atom_count < 10:
                delta_list.append(60.0)
            elif atom_count < 20:
                delta_list.append(45.0)
            elif atom_count < 30:
                delta_list.append(30.0)
            elif atom_count < 40:
                delta_list.append(20.0)
            else:
                delta_list.append(15.0)
        return delta_list

    # Single-component fallback
    pool = molecules[0] if isinstance(molecules, (list, tuple)) and molecules else None
    if not pool:
        return DEFAULT
    mol = pool[0]
    try:
        atom_count = len(mol.mol)
    except Exception:
        atom_count = len(getattr(mol, "atoms", []))
    if atom_count <= 3:
        return 90.0
    elif atom_count < 10:
        return 60.0
    elif atom_count < 20:
        return 45.0
    elif atom_count < 30:
        return 30.0
    elif atom_count < 40:
        return 20.0
    else:
        return 15.0


def plot_id_vs_energy(code, energies, match_ids=None, match_energies=None, out_dir="qrs_plots", time_cost_s=None, coverage=None, n_conformers=None, energy_unit="kcal/mol"):
    """Save a plot of visited-structure ID vs energy for one QRS run."""
    if not energies:
        print(f"No energies collected for {code}; skipping plot.")
        return

    os.makedirs(out_dir, exist_ok=True)
    ids = list(range(1, len(energies) + 1))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(ids, energies, s=10, alpha=0.8, label="Visited")
    ax.plot(ids, energies, linewidth=0.8, alpha=0.6)
    if match_ids and match_energies:
        ax.scatter(
            match_ids,
            match_energies,
            s=70,
            marker="*",
            c="crimson",
            edgecolors="black",
            linewidths=0.6,
            zorder=3,
            label="Match",
        )
    ax.set_xlabel("Visited Structure ID")
    ax.set_ylabel(f"Energy ({energy_unit})")
    title_parts = [f"{code}: "]
    if n_conformers is not None:
        title_parts.append(f"conf: {n_conformers}")
    if time_cost_s is not None:
        title_parts.append(f"time: {time_cost_s:.2f} s")
    if coverage is not None:
        title_parts.append(f"coverage: {coverage}")
    if len(title_parts) > 1:
        ax.set_title(f"{title_parts[0]} ({', '.join(title_parts[1:])})")
    else:
        ax.set_title(title_parts[0])
    ax.grid(alpha=0.25, linestyle="--", linewidth=0.6)
    ax.legend(loc="best")

    ymin = min(energies)
    ymax = max(energies)
    if ymin == ymax:
        margin = max(abs(ymin) * 0.05, 1.0)
        ax.set_ylim(ymin - margin, ymax + margin)
    else:
        y_max = min(ymin + 30, ymax)
        ax.set_ylim(ymin - 0.25, y_max)
    fig.tight_layout()

    fig_path = os.path.join(out_dir, f"{code}.png")
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Saved plot: {fig_path}")


if __name__ == "__main__":
    db = database("pyxtal/database/test.db")
    out_dir = "Tests-0611"
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "qrs_results.csv")

    with open(csv_path, "w", newline="") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(
            [
                "code",
                "smiles",
                "sg",
                "vector_dim",
                "n_pregenerated_confs",
                "time_cost_s",
                "success_rate",
                "result",
                "coverage",
                "n_conformers",
            ]
        )

    for code in db.get_all_codes():
        #if code not in ['ACSALA']: continue
        #if code not in ['FUNZOE']: continue
        #if code not in ['XAFQON']: continue
        #if code not in ['ACEMID02']: continue
        row = db.get_row(code=code)
        ref_xtal = db.get_pyxtal(code=code)
        if ref_xtal.has_special_site():
            ref_xtal = ref_xtal.to_subgroup()

        vector_dim = get_vector_dimension(ref_xtal)
        ref_pmg = ref_xtal.to_pymatgen()
        ref_pmg.remove_species(["H"])  # Ignore H for matching robustness.

        print(f"\n=== {code} ===")
        print(ref_xtal)

        workdir = os.path.join(out_dir, row.csd_code)
        os.makedirs(workdir, exist_ok=True)
        sites = build_sites_from_reference(ref_xtal)

        # Pregenerate per-component conformer pools aligned with site.type.
        # QRS expects one pool per molecular component in row.mol_smi.split('.').
        smiles_parts = row.mol_smi.split(".")
        n_types = len(ref_xtal.numMols)
        if len(smiles_parts) != n_types:
            print(
                f"SMILES/component mismatch for {code}: "
                f"{len(smiles_parts)} smiles parts vs {n_types} crystal components; skipping."
            )
            continue

        type_wps = [[] for _ in range(n_types)]
        for site in ref_xtal.mol_sites:
            type_wps[site.type].append(site.wp)

        molecules = []
        n_pregen_total = 0
        for type_idx, smi in enumerate(smiles_parts):
            smi = smi.strip()
            if not smi:
                print(f"Empty SMILES for component {type_idx}; skipping {code}.")
                molecules = None
                break

            if smi in single_smiles:
                # Monoatomic ions (e.g., Cl-) have no conformer degrees of freedom.
                # Use one fixed molecule instead of running conformer generation.
                try:
                    m0 = pyxtal_molecule(smi + ".smi", fix=True)
                    _, valid = m0.get_orientations_in_wps(type_wps[type_idx])
                except Exception as exc:
                    print(
                        f"Failed to build single-component molecule for component {type_idx} ({smi}) in {code}: "
                        f"{exc}; skipping."
                    )
                    molecules = None
                    break

                if not valid:
                    print(
                        f"Single-component molecule {smi} has no valid orientation in component {type_idx}; "
                        f"skipping {code}."
                    )
                    molecules = None
                    break

                p_mols = [m0]
                print(f"Component {type_idx} ({smi}) single-species pool: 1")
            else:
                try:
                    p_mols = generate_molecules(
                        smi,
                        wps=type_wps[type_idx],
                        N_iter=20,
                        N_conf=200,
                        tol=0.5,
                    )
                except Exception as exc:
                    print(
                        f"Failed to pregenerate conformers for component {type_idx} ({smi}) in {code}: "
                        f"{exc}; skipping."
                    )
                    molecules = None
                    break

            if len(p_mols) == 0:
                print(f"No valid pregenerated conformers for component {type_idx} ({smi}); skipping.")
                molecules = None
                break

            p_mols = filter_similar_molecules(p_mols, rmsd_tol=0.5)
            print(f"Component {type_idx} ({smi}) unique conformers: {len(p_mols)}")
            if len(p_mols) == 0:
                print(f"All conformers filtered out for component {type_idx} ({smi}); skipping.")
                molecules = None
                break

            molecules.append(p_mols)
            n_pregen_total += len(p_mols)

        if molecules is None:
            continue
        print(f"Total pregenerated conformers across components: {n_pregen_total}")

        param_xml = os.path.join(workdir, "parameters.xml")
        if os.path.exists(param_xml):
            os.remove(param_xml)
        print(f"Initialized QRS workdir: {workdir}", ref_xtal.lattice)
        composition = [int(a) for a in ref_xtal.get_zprime()]
        # determine per-component delta angles and show them
        selected_deltas = select_delta_angle(molecules, composition)
        print(f"Selected delta_angle(s) for components: {selected_deltas}")

        qrs = QRS(
            smiles=row.mol_smi,
            workdir=workdir,
            sg=ref_xtal.group.hall_number,
            tag=row.csd_code.lower(),
            use_hall=True,
            lattice=ref_xtal.lattice,  # Fixed cell.
            composition = composition,
            molecules=molecules,
            sites=sites,
            N_gen=200,
            N_pop=96,
            N_cpu=48,
            cif="all.cif",
            skip_mlp=True,
            verbose=False,
            delta_length=1.2,
            delta_angle=selected_deltas,
        )

        t0 = perf_counter()
        success_rate = qrs.run(ref_pmg=ref_pmg)
        time_cost_s = perf_counter() - t0

        if success_rate is not None and success_rate > 0:
            print(f"Match found! Success rate: {success_rate}%")
        else:
            print("No match found within the given generations/population.")
        print(f"Time cost: {time_cost_s:.2f} s")
        visited_energies = get_visited_energies(qrs)
        match_ids, match_energies = get_match_points(qrs)
        coverage = None
        if hasattr(qrs, "sampler") and hasattr(qrs.sampler, "current") and qrs.sampler.total:
            pct = 100.0 * qrs.sampler.current / qrs.sampler.total
            coverage = f"{qrs.sampler.current}/{qrs.sampler.total} [{pct:.1f}%]"
        plot_id_vs_energy(
            code,
            visited_energies,
            match_ids=match_ids,
            match_energies=match_energies,
            out_dir=os.path.join(out_dir, "qrs_plots"),
            time_cost_s=time_cost_s,
            coverage=coverage,
            n_conformers=n_pregen_total,
            energy_unit="eV/atom" if not qrs.skip_mlp else "kcal/mol",
        )

        with open(csv_path, "a", newline="") as fcsv:
            writer = csv.writer(fcsv)
            writer.writerow(
                [
                    code,
                    row.mol_smi,
                    ref_xtal.group.number,
                    vector_dim,
                    n_pregen_total,
                    f"{time_cost_s:.2f}",
                    f"{success_rate:.4f}" if success_rate is not None else "",
                    "success" if success_rate is not None and success_rate > 0 else "failure",
                    coverage if coverage is not None else "",
                    n_pregen_total,
                ]
            )

    print(f"\nSaved summary CSV: {csv_path}")

