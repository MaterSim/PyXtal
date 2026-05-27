"""
Example: Run QRS crystal structure prediction with known cell parameters
and check whether the reference structure is recovered.

Iterates over all entries in pyxtal/database/test.db, fixes the lattice
for each, runs QRS, and appends results (code, SMILES, SG, success_rate)
to qrs_results.csv.
"""

import csv

from pyxtal.optimize import QRS
from pyxtal.db import database

db = database("pyxtal/database/test.db")

csv_path = "qrs_results.csv"
with open(csv_path, "w", newline="") as fcsv:
    writer = csv.writer(fcsv)
    writer.writerow(["code", "smiles", "sg", "success_rate"])

for code in db.get_all_codes()[3:]:
    # ── Molecule ──────────────────────────────────────────────────────────────────
    row = db.get_row(code=code)
    ref_xtal = db.get_pyxtal(code=code)
    if ref_xtal.has_special_site(): ref_xtal = ref_xtal.to_subgroup()
    ref_pmg  = ref_xtal.to_pymatgen()
    ref_pmg.remove_species(["H"])  # ignore H for matching since positions are less certain
    print(ref_xtal)

    sites = [[] for _ in range(len(ref_xtal.numMols))]
    for site in ref_xtal.mol_sites:
        sites[site.type].append(site.wp.get_label())

    # ── QRS setup ─────────────────────────────────────────────────────────────────
    qrs = QRS(
        smiles        = row.mol_smi,                # molecule as SMILES string
        workdir       = row.csd_code,               # working directory for this QRS run
        sg            = ref_xtal.group.hall_number, # space group number (P2_1/c = 81)
        tag           = row.csd_code.lower(),       # tag for output files
        use_hall      = True,                       # interpret sg as a Hall number
        lattice       = ref_xtal.lattice,           # fix the cell; only WP positions are sampled
        N_gen         = 10,                         # number of QRS generations
        N_pop         = 50,                         # structures per generation
        N_cpu         = 1,
        cif           = "all.cif",                  # save all relaxed structures
        skip_mlp      = True,                       # no machine-learning potential relaxation
        verbose       = False,
        sites         = sites,
        delta_length  = 1.5,                        # grid spacing for fractional coords (Å)
        delta_angle   = 60.0,                       # grid spacing for Euler/torsion angles (°)
    )

    # ── Run and check for match ───────────────────────────────────────────────────
    success_rate = qrs.run(ref_pmg=ref_pmg)

    if success_rate is not None and success_rate > 0:
        print(f"\nMatch found! Success rate: {success_rate}%")
    else:
        print("\nNo match found within the given generations/population.")

    with open(csv_path, "a", newline="") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow([code, row.mol_smi, ref_xtal.group.number,
                         f"{success_rate:.4f}" if success_rate is not None else ""])
