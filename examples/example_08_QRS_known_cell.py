"""
Example: Run QRS crystal structure prediction with known cell parameters
and check whether the reference structure is recovered.

Iterates over all entries in pyxtal/database/test.db, fixes the lattice
for each, runs QRS, and appends results (code, SMILES, SG, success_rate)
to qrs_results.csv.
"""

import csv
import os
from time import perf_counter

import matplotlib.pyplot as plt

from pyxtal.db import database
from pyxtal.optimize import QRS


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


def get_vector_dimension(xtal):
    """Return the total dimension of sampling vectors for this crystal."""
    return sum(len(site.get_bounds()) for site in xtal.mol_sites)


def plot_id_vs_energy(code, energies, match_ids=None, match_energies=None, out_dir="qrs_plots"):
    """Save a plot of visited-structure ID vs energy for one QRS run."""
    if not energies:
        print(f"No energies collected for {code}; skipping plot.")
        return

    os.makedirs(out_dir, exist_ok=True)
    ids = list(range(1, len(energies) + 1))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(ids, energies, s=20, alpha=0.8, label="Visited")
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
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title(f"{code}: visited structure ID vs energy")
    ax.grid(alpha=0.25, linestyle="--", linewidth=0.6)
    ax.legend(loc="best")
    fig.tight_layout()

    fig_path = os.path.join(out_dir, f"{code}_id_vs_energy.png")
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Saved plot: {fig_path}")

db = database("pyxtal/database/test.db")

csv_path = "qrs_results.csv"
with open(csv_path, "w", newline="") as fcsv:
    writer = csv.writer(fcsv)
    writer.writerow(["code", "smiles", "sg", "vector_dim", "time_cost_s", "success_rate"])

for code in db.get_all_codes()[2:]:
    # ── Molecule ──────────────────────────────────────────────────────────────────
    row = db.get_row(code=code)
    ref_xtal = db.get_pyxtal(code=code)
    if ref_xtal.has_special_site(): ref_xtal = ref_xtal.to_subgroup()
    vector_dim = get_vector_dimension(ref_xtal)
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
    t0 = perf_counter()
    success_rate = qrs.run(ref_pmg=ref_pmg)
    time_cost_s = perf_counter() - t0

    if success_rate is not None and success_rate > 0:
        print(f"\nMatch found! Success rate: {success_rate}%")
    else:
        print("\nNo match found within the given generations/population.")
    print(f"Time cost: {time_cost_s:.2f} s")

    match_ids, match_energies = get_match_points(qrs)
    plot_id_vs_energy(code, qrs.engs, match_ids, match_energies)

    with open(csv_path, "a", newline="") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow([code, row.mol_smi, ref_xtal.group.number, vector_dim,
                         f"{time_cost_s:.2f}",
                         f"{success_rate:.4f}" if success_rate is not None else ""])
