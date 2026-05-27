"""
Example: Run QRS with pregenerated molecular conformers.

This script iterates over entries in pyxtal/database/test.db, precomputes a pool
of pyxtal_molecule conformers via generate_molecules, and passes them to QRS.
With the lattice fixed and conformers pregenerated, QRS samples only Wyckoff-site
fractional coordinates and molecular orientations for each trial crystal.

Results are appended to Tests/qrs_results_pregen_mols.csv.
"""

import csv
import os
from time import perf_counter

from pyxtal.db import database
from pyxtal.molecule import generate_molecules
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


if __name__ == "__main__":
    db = database("pyxtal/database/test.db")
    os.makedirs("Tests", exist_ok=True)
    csv_path = "Tests/qrs_results_pregen_mols.csv"

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
            ]
        )

    for code in db.get_all_codes()[5:8]:
        row = db.get_row(code=code)
        ref_xtal = db.get_pyxtal(code=code)
        if ref_xtal.has_special_site():
            ref_xtal = ref_xtal.to_subgroup()

        vector_dim = get_vector_dimension(ref_xtal)
        ref_pmg = ref_xtal.to_pymatgen()
        ref_pmg.remove_species(["H"])  # Ignore H for matching robustness.

        print(f"\n=== {code} ===")
        print(ref_xtal)

        workdir = os.path.join("Tests", row.csd_code)
        os.makedirs(workdir, exist_ok=True)
        sites = build_sites_from_reference(ref_xtal)

        # Pregenerate molecular conformers that are compatible with the target WP.
        # We provide the resulting conformer pool to QRS so torsions are fixed per
        # sampled molecule choice and QRS only samples WP xyz + orientation DOFs.
        target_wps = [site.wp for site in ref_xtal.mol_sites]
        p_mols = generate_molecules(
            row.mol_smi,
            wps=target_wps,
            N_iter=8,
            N_conf=50,
            tol=0.5,
        )
        if len(p_mols) == 0:
            print("No valid pregenerated conformers; skipping this entry.")
            continue
        p_mols = filter_similar_molecules(p_mols, rmsd_tol=0.5)
        print(f"Unique conformers after filtering: {len(p_mols)}")
        if len(p_mols) == 0:
            print("All pregenerated conformers were filtered out; skipping this entry.")
            continue
        print(f"Pregenerated conformers: {len(p_mols)}")

        param_xml = os.path.join(workdir, "parameters.xml")
        if os.path.exists(param_xml):
            os.remove(param_xml)

        qrs = QRS(
            smiles=row.mol_smi,
            workdir=workdir,
            sg=ref_xtal.group.hall_number,
            tag=row.csd_code.lower(),
            use_hall=True,
            lattice=ref_xtal.lattice,  # Fixed cell.
            molecules=[p_mols],        # One molecular component with many conformers.
            sites=sites,
            N_gen=20,
            N_pop=50,
            N_cpu=2,
            cif="all.cif",
            skip_mlp=True,
            verbose=False,
            delta_length=1.5,
            delta_angle=45.0,
        )

        t0 = perf_counter()
        success_rate = qrs.run(ref_pmg=ref_pmg)
        time_cost_s = perf_counter() - t0

        if success_rate is not None and success_rate > 0:
            print(f"Match found! Success rate: {success_rate}%")
        else:
            print("No match found within the given generations/population.")
        print(f"Time cost: {time_cost_s:.2f} s")

        with open(csv_path, "a", newline="") as fcsv:
            writer = csv.writer(fcsv)
            writer.writerow(
                [
                    code,
                    row.mol_smi,
                    ref_xtal.group.number,
                    vector_dim,
                    len(p_mols),
                    f"{time_cost_s:.2f}",
                    f"{success_rate:.4f}" if success_rate is not None else "",
                ]
            )

    print(f"\nSaved summary CSV: {csv_path}")
