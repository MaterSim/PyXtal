"""
Example: Run QRS with pregenerated molecular conformers.

This script iterates over entries in pyxtal/database/test.db, precomputes a pool
of pyxtal_molecule conformers via generate_molecules, and passes them to QRS.
With the lattice fixed and conformers pregenerated, QRS samples only Wyckoff-site
fractional coordinates and molecular orientations for each trial crystal.

Results are appended to <out-dir>/qrs_results.csv (default: Tests-0611).
"""

import argparse
import csv
import os
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher

from pyxtal.constants import single_smiles
from pyxtal.db import database
from pyxtal.molecule import generate_molecules, get_inertia_tensor, pyxtal_molecule
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


def pregen_component_pool(smi, wps, args):
    """
    Build the conformer pool for one molecular component.

    Rigid molecules (no rotatable bonds) use a single fixed ``pyxtal_molecule``,
    matching the fast path that worked for cases like MERQUY in Tests-0719.
    """
    if smi in single_smiles:
        m0 = pyxtal_molecule(smi + ".smi", fix=True)
        _, valid = m0.get_orientations_in_wps(wps)
        if not valid:
            return None, "single-species"
        return [m0], "single-species"

    m0 = pyxtal_molecule(smi + ".smi", fix=True)
    if not m0.torsionlist:
        _, valid = m0.get_orientations_in_wps(wps)
        if not valid:
            return None, "rigid"
        return [m0], "rigid"

    pool = generate_molecules(
        smi,
        wps=wps,
        N_iter=20,
        N_conf=200,
        tol=0.5,
        torsion_extras=args.torsion_extras,
        verbose=True,
    )
    return pool, "flexible"


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

    Monoatomic species (1 atom) have no orientation DOFs; their entry is
    ``None`` so callers can omit a meaningless angle resolution.

    Rules (per multi-atom component):
      - max_atoms <= 3 -> 90.0
      - max_atoms < 10 -> 60.0
      - max_atoms < 20 -> 45.0
      - max_atoms < 30 -> 30.0
      - max_atoms < 40 -> 20.0
      - else           -> 15.0
    """
    DEFAULT = 15.0

    def _atom_count(mol):
        try:
            return len(mol.mol)
        except Exception:
            return len(getattr(mol, "atoms", []))

    def _delta_for_count(atom_count):
        if atom_count <= 1:
            return None
        if atom_count <= 3:
            return 90.0
        if atom_count < 10:
            return 60.0
        if atom_count < 20:
            return 45.0
        if atom_count < 30:
            return 30.0
        if atom_count < 40:
            return 20.0
        return 15.0

    if not molecules:
        if composition is None or len(composition) <= 1:
            return DEFAULT
        return [DEFAULT for _ in composition]

    if composition is not None and len(composition) > 1:
        delta_list = []
        for idx, _count in enumerate(composition):
            pool = molecules[idx] if idx < len(molecules) else None
            if not pool:
                delta_list.append(DEFAULT)
                continue
            delta_list.append(_delta_for_count(_atom_count(pool[0])))
        return delta_list

    pool = molecules[0] if isinstance(molecules, (list, tuple)) and molecules else None
    if not pool:
        return DEFAULT
    return _delta_for_count(_atom_count(pool[0]))


def _molecule_atom_count(mol):
    try:
        return len(mol.mol)
    except Exception:
        return len(getattr(mol, "atoms", []))


def _component_pool(molecules, idx):
    if not molecules or idx >= len(molecules):
        return None
    pool = molecules[idx]
    return pool if pool else None


def _clamp_delta(value, lo=10.0, hi=90.0):
    return max(lo, min(hi, float(value)))


def compute_molecular_shape_metrics(mol):
    """Estimate molecular shape anisotropy from the inertia tensor.

    Returns a dict with principal moments ``i1 <= i2 <= i3``, axis ratios,
    and a coarse shape label: ``monoatomic``, ``compact``, ``rod``, or ``planar``.
    """
    atom_count = _molecule_atom_count(mol)
    if atom_count <= 1:
        return {
            "shape": "monoatomic",
            "atom_count": atom_count,
            "ratios": (1.0, 1.0, 1.0),
            "anisotropy": 0.0,
        }

    coords = np.asarray(mol.mol.cart_coords, dtype=float)
    numbers = np.asarray([site.specie.Z for site in mol.mol], dtype=int)
    heavy = numbers > 1
    if heavy.sum() < 2:
        heavy = np.ones(len(numbers), dtype=bool)

    coords = coords[heavy]
    weights = numbers[heavy].astype(float)
    inertia = get_inertia_tensor(coords.copy(), weights=weights)
    moments = np.sort(np.linalg.eigvalsh(inertia))
    i1, i2, i3 = moments
    eps = max(float(i1), 1e-6)

    r21 = i2 / eps
    r31 = i3 / eps
    r32 = i3 / max(i2, 1e-6)
    anisotropy = 1.0 - (i1 / max(i3, 1e-6))

    if r31 < 2.0:
        shape = "compact"
    elif abs(r32 - 1.0) <= 0.25:
        shape = "rod"
    else:
        shape = "planar"

    return {
        "shape": shape,
        "atom_count": atom_count,
        "moments": (float(i1), float(i2), float(i3)),
        "ratios": (float(r21), float(r31), float(r32)),
        "anisotropy": float(anisotropy),
        "extents": tuple(float(x) for x in np.ptp(coords, axis=0)),
    }


def _molecular_extents(mol):
    """Cartesian spans (Å) along lab axes for one pyxtal_molecule."""
    coords = np.asarray(mol.mol.cart_coords, dtype=float)
    return np.ptp(coords, axis=0)


def _packing_strain_ratio(mol, lattice):
    """
    How tightly a molecule must wrap against the fixed cell.

    Anisotropic molecules (small min/max extent ratio) that exceed the
    shortest cell edge rely on orientation/PBC and need a looser soft-clash
    buffer. The aspect cutoff is 0.45 so rods such as LEVJON (aspect ~0.37)
    still use plate/rod strain, not only the sorted-edge aligned ratio.
    """
    extents = _molecular_extents(mol)
    cell = np.asarray(lattice.get_para()[:3], dtype=float)
    aspect = float(np.min(extents) / max(np.max(extents), 1e-6))
    plate_strain = float(np.max(extents) / max(np.min(cell), 1e-6))
    ext_sorted = np.sort(extents)[::-1]
    cell_sorted = np.sort(cell)[::-1]
    aligned_strain = float(np.max(ext_sorted / np.maximum(cell_sorted, 1e-6)))
    if aspect < 0.45:
        return max(plate_strain, aligned_strain)
    return aligned_strain


def select_soft_clash_buffer(molecules, lattice, composition=None):
    """
    Pick a soft-clash buffer from molecular packing vs the fixed cell.

    Unrelaxed grid packings for plate-like / rod-like molecules in short cells
    (extent exceeding the shortest axis) often fail Tol_matrix at buffer 0 even
    when CHARMM relaxes to a valid structure (e.g. MERQUY, LEVJON).
    """
    if not molecules or lattice is None:
        return 0.0, {"mode": "default", "strain": 0.0}

    if composition is None:
        composition = [1] * len(molecules)

    strains = []
    for idx in range(len(composition)):
        pool = _component_pool(molecules, idx)
        if not pool:
            continue
        mol = pool[0]
        strain = _packing_strain_ratio(mol, lattice)
        shape = compute_molecular_shape_metrics(mol)["shape"]
        numbers = np.asarray(mol.mol.atomic_numbers, dtype=int)
        n_heavy = int(np.sum(numbers > 1))
        strains.append((idx, strain, shape, n_heavy))

    if not strains:
        return 0.0, {"mode": "default", "strain": 0.0}

    idx, strain, shape, n_heavy = max(strains, key=lambda row: row[1])
    if strain < 1.15:
        buffer = 0.0
    elif strain < 1.45:
        buffer = -0.2
    elif strain < 1.65:
        buffer = -0.5
    elif strain < 1.85:
        buffer = -0.7
    else:
        buffer = -1.0

    # Bulky organics still soft-clash heavily on the raw grid even at modest
    # strain (many F/heavy contacts). Enforce a looser floor.
    if n_heavy >= 40:
        buffer = min(buffer, -0.5)
    elif n_heavy >= 30:
        buffer = min(buffer, -0.2)

    # Planar / rod Z'>=2 packs (e.g. CEQGEL) reject ~99.5% at buffer 0 even
    # when strain < 1.15, turning Gen-0 into a 5–10 min serial filter. A
    # −0.5 floor matches Tests-0810 (match in ~14s/gen) without going as
    # loose as −1.5.
    z_same = max(int(c) for c in composition) if composition else 1
    if shape in ("planar", "rod") and z_same >= 2 and n_heavy >= 16:
        buffer = min(buffer, -0.5)

    return buffer, {
        "mode": "auto",
        "component": idx,
        "shape": shape,
        "strain": strain,
        "n_heavy": n_heavy,
    }


def resolve_soft_clash_buffer_arg(value):
    """Return ``None`` for auto selection, otherwise a float buffer."""
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in ("auto", "none"):
        return None
    return float(value)


def select_delta_angle_by_shape(
    molecules,
    composition=None,
    *,
    linear_ratio=2.5,
    rod_spin_factor=2.0,
    tilt_factor=0.5,
    planar_normal_factor=2.0,
):
    """Pick direction-dependent delta angles from molecular shape.

    Uses ``select_delta_angle`` for a baseline, then for strictly anisotropic
    molecules returns per-component Euler resolutions ``[alpha, beta, gamma]``:

    - **rod**: coarse rotation around the long axis, finer tilts
    - **planar**: finer in-plane rotations, coarser out-of-plane tilt
    - **compact/monoatomic**: scalar baseline (``None`` for monoatomic)
    """
    if not molecules:
        return select_delta_angle(molecules, composition)

    if composition is None or len(composition) <= 1:
        composition = [1]

    delta_specs = []
    for idx in range(len(composition)):
        pool = _component_pool(molecules, idx)
        if not pool:
            delta_specs.append(15.0)
            continue

        mol = pool[0]
        metrics = compute_molecular_shape_metrics(mol)
        shape = metrics["shape"]
        if shape == "monoatomic":
            delta_specs.append(None)
            print(f"  component {idx}: shape=monoatomic, delta=None", flush=True)
            continue

        base = select_delta_angle([pool], None)
        if base is None:
            delta_specs.append(None)
            continue

        _, r31, _r32 = metrics["ratios"]
        if shape == "compact" or r31 < linear_ratio:
            delta_specs.append(base)
            print(
                f"  component {idx}: shape={shape}, base delta={base}",
                flush=True,
            )
            continue

        if shape == "rod":
            euler = [
                _clamp_delta(base * rod_spin_factor),
                _clamp_delta(base * tilt_factor),
                _clamp_delta(base * tilt_factor),
            ]
        else:  # planar
            euler = [
                _clamp_delta(base * tilt_factor),
                _clamp_delta(base * planar_normal_factor),
                _clamp_delta(base * tilt_factor),
            ]

        delta_specs.append(euler)
        print(
            f"  component {idx}: shape={shape}, ratios={metrics['ratios']}, "
            f"delta_angle=[alpha={euler[0]}, beta={euler[1]}, gamma={euler[2]}]",
            flush=True,
        )

    if len(composition) <= 1:
        return delta_specs[0]
    return delta_specs


def plot_id_vs_energy(
    code,
    energies,
    match_ids=None,
    match_energies=None,
    out_dir="qrs_plots",
    time_min=None,
    coverage=None,
    n_conformers=None,
    energy_unit="kcal/mol",
):
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
    if time_min is not None:
        title_parts.append(f"time: {time_min:.2f} min")
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
        # Default: zoom to ~30 kcal/mol above the best energy.
        # Expand if matches (or almost all points) would otherwise be clipped —
        # absolute FF energies can span >>30 kcal/mol (e.g. UJIRIO02).
        y_hi = min(ymin + 30.0, ymax)
        if match_energies:
            y_hi = max(y_hi, max(match_energies))
        n_vis = sum(1 for e in energies if e <= y_hi + 1e-9)
        if len(energies) >= 20 and n_vis < max(20, 0.05 * len(energies)):
            y_hi = float(np.percentile(energies, 95))
            if match_energies:
                y_hi = max(y_hi, max(match_energies))
            y_hi = min(y_hi, ymax)
        ax.set_ylim(ymin - 0.25, y_hi + 0.25)
    fig.tight_layout()

    fig_path = os.path.join(out_dir, f"{code}.png")
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Saved plot: {fig_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run QRS with pregenerated molecular conformers for crystals in test.db",
    )
    parser.add_argument(
        "--out-dir",
        default="Tests-0611",
        help="Output root directory for per-case QRS results (default: Tests-0611)",
    )
    parser.add_argument(
        "--nproc",
        default=48,
        type=int,
        help="Number of CPUs for QRS local optimization (default: 48)",
    )
    parser.add_argument(
        "--ngen",
        default=200,
        type=int,
        help="Number of generations for QRS (default: 200)",
    )
    parser.add_argument(
        "--npop",
        default=96,
        type=int,
        help="Number of populations for QRS (default: 96)",
    )
    parser.add_argument(
        "--delta-length",
        default=1.1,
        type=float,
        help="Length grid spacing in Angstrom for fractional Wyckoff coords (default: 1.1)",
    )
    parser.add_argument(
        "--delta-angle-mode",
        choices=["size", "shape"],
        default="size",
        help="How to choose delta_angle: by molecule size (default) or shape anisotropy",
    )
    parser.add_argument(
        "--shape-linear-ratio",
        default=2.5,
        type=float,
        help="I3/I1 threshold to treat a molecule as shape-anisotropic (default: 2.5)",
    )
    parser.add_argument(
        "--max-grid-product",
        default=0,
        type=float,
        help="Cap on QRS grid product prod(n_levels) (default: 0 = disabled; >0 enables)",
    )
    parser.add_argument(
        "--min-grid-per-gen",
        default=0,
        type=int,
        help=(
            "Minimum QRS grid points evaluated per generation even after "
            "N_pop is filled (default: 0 = stop once population is full; "
            "try 8192 for easy-filter multi-component cases like XAFQON)"
        ),
    )
    parser.add_argument(
        "--match-ltol",
        default=0.3,
        type=float,
        help="StructureMatcher lattice length tolerance (default: 0.3)",
    )
    parser.add_argument(
        "--match-stol",
        default=0.3,
        type=float,
        help="StructureMatcher site distance tolerance (default: 0.3)",
    )
    parser.add_argument(
        "--match-angle-tol",
        default=5.0,
        type=float,
        help="StructureMatcher angle tolerance in degrees (default: 5.0)",
    )
    parser.add_argument(
        "--max-rmsd",
        default=0.3,
        type=float,
        help="Maximum RMSD for structure matching (default: 0.3)",
    )
    parser.add_argument(
        "--min-matches",
        default=5,
        type=int,
        help=(
            "Quit early after this many matched structures "
            "(default: 5; set very large to disable)"
        ),
    )
    parser.add_argument(
        "--check-stable",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "After each local relaxation, deterministically perturb Wyckoff DOFs "
            "by +/- half the QRS grid spacing and re-relax to detect shallow minima "
            "(default: enabled; use --no-check-stable to disable)"
        ),
    )
    parser.add_argument(
        "--soft-clash",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Reject Tol_matrix soft clashes before FF relaxation "
            "(default: enabled; use --no-soft-clash for hard short-distance "
            "only, as in older QRS runs like Tests-0708)"
        ),
    )
    parser.add_argument(
        "--soft-clash-buffer",
        default="auto",
        help=(
            "Add this value in Angstrom to molecular Tol_matrix cutoffs "
            "(default: auto from molecular extent vs cell; negative loosens, "
            "positive tightens; use 'auto' or a float; ignored with "
            "--no-soft-clash)"
        ),
    )
    parser.add_argument(
        "--close-grid-cutoff",
        default="0",
        help=(
            "Skip conformer checks when Wyckoff-orbit site centers are too "
            "close. Use a float (Angstrom), 'auto' for "
            "alpha*0.5*(r_i+r_j) with pool min radii, or 0/off to disable "
            "(default: 0)"
        ),
    )
    parser.add_argument(
        "--close-grid-alpha",
        default=0.7,
        type=float,
        help=(
            "Scale factor for --close-grid-cutoff auto "
            "(default: 0.7; try 0.6-0.8)"
        ),
    )
    parser.add_argument(
        "--fix-translation",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Pin continuous origin-free axes on one general-position site "
            "using Group.get_free_axis() (default: enabled; no-op when the "
            "space group has no free axes)"
        ),
    )
    parser.add_argument(
        "--asu-clamp",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Restrict that site's fractional bounds to the conventional ASU "
            "box (default: enabled; removes discrete space-group copies)"
        ),
    )
    parser.add_argument(
        "--gauge-site",
        default=0,
        type=int,
        help="Molecular site index used for translational gauge / ASU (default: 0)",
    )
    parser.add_argument(
        "--order-identical-sites",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "For chemically identical Z'>=2 sites with matching xyz bounds, "
            "require strictly increasing translational grid index (j>i) to "
            "remove label-exchange double counting (default: enabled)"
        ),
    )
    parser.add_argument(
        "--max-mol-combos-per-grid",
        default=0,
        type=int,
        help=(
            "Optional cap on conformer combinations per grid point during "
            "filtering (default: 0 = full pool; set e.g. 32 to subsample)"
        ),
    )
    parser.add_argument(
        "--max-valid-per-grid",
        default=0,
        type=int,
        help=(
            "Optional cap on valid conformers kept per grid point "
            "(default: 0 = keep all that pass filtering)"
        ),
    )
    parser.add_argument(
        "--code",
        nargs="+",
        metavar="CSD_CODE",
        help="Run only these CSD codes from test.db (default: all codes)",
    )
    parser.add_argument(
        "--torsion-extras",
        choices=["auto", "on", "off"],
        default="auto",
        help=(
            "Pregen torsion extras (default auto): charged→raw grid; "
            "bulky+many rotors→anti-bias ±180°; bulky+few rotors→off; "
            "medium flexible→MMFF-refined; on=force grid; off=MMFF only"
        ),
    )
    args = parser.parse_args()

    matcher = StructureMatcher(
        ltol=args.match_ltol,
        stol=args.match_stol,
        angle_tol=args.match_angle_tol,
    )
    print(
        f"StructureMatcher: ltol={args.match_ltol}, stol={args.match_stol}, "
        f"angle_tol={args.match_angle_tol}",
    )

    db = database("pyxtal/database/test.db")
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "qrs_results.csv")

    if args.code:
        codes = list(args.code)
        unknown = [c for c in codes if c not in db.get_all_codes()]
        if unknown:
            parser.error(f"Unknown CSD code(s) in test.db: {', '.join(unknown)}")
    else:
        codes = db.get_all_codes()

    with open(csv_path, "w", newline="") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(
            [
                "code",
                "smiles",
                "sg",
                "vector_dim",
                "time_min",
                "success_rate",
                "result",
                "coverage",
                "n_conformers",
            ]
        )

    max_grid_product = args.max_grid_product if args.max_grid_product > 0 else None

    for code in codes:
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

            try:
                p_mols, pool_kind = pregen_component_pool(smi, type_wps[type_idx], args)
            except Exception as exc:
                print(
                    f"Failed to pregenerate conformers for component {type_idx} ({smi}) in {code}: "
                    f"{exc}; skipping."
                )
                molecules = None
                break

            if p_mols is None:
                print(
                    f"Component {type_idx} ({smi}) has no valid orientation in Wyckoff sites; "
                    f"skipping {code}."
                )
                molecules = None
                break

            if len(p_mols) == 0:
                print(f"No valid pregenerated conformers for component {type_idx} ({smi}); skipping.")
                molecules = None
                break

            p_mols = filter_similar_molecules(p_mols, rmsd_tol=0.5)
            print(
                f"Component {type_idx} ({smi}) {pool_kind} pool: {len(p_mols)} unique conformers"
            )
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

        if args.delta_angle_mode == "shape":
            print("Selecting direction-dependent delta_angle(s) from molecular shape:")
            selected_deltas = select_delta_angle_by_shape(
                molecules,
                composition,
                linear_ratio=args.shape_linear_ratio,
            )
        else:
            selected_deltas = select_delta_angle(molecules, composition)
        print(f"Selected delta_angle(s) for components: {selected_deltas}")

        if not args.soft_clash:
            soft_clash_buffer = 0.0
            print("Soft-clash filter disabled (--no-soft-clash); using hard short-distance check")
        else:
            manual_soft_buffer = resolve_soft_clash_buffer_arg(args.soft_clash_buffer)
            if manual_soft_buffer is None:
                soft_clash_buffer, buffer_info = select_soft_clash_buffer(
                    molecules,
                    ref_xtal.lattice,
                    composition,
                )
                heavy_txt = (
                    f", n_heavy={buffer_info['n_heavy']}"
                    if buffer_info.get("n_heavy") is not None
                    else ""
                )
                print(
                    "Selected soft_clash_buffer: "
                    f"{soft_clash_buffer:+.1f} Angstrom "
                    f"(component {buffer_info['component']}, "
                    f"shape={buffer_info['shape']}, "
                    f"packing_strain={buffer_info['strain']:.2f}"
                    f"{heavy_txt})"
                )
            else:
                soft_clash_buffer = manual_soft_buffer
                print(f"Using manual soft_clash_buffer: {soft_clash_buffer:+.1f} Angstrom")

        qrs = QRS(
            smiles=row.mol_smi,
            workdir=workdir,
            sg=ref_xtal.group.hall_number,
            tag=row.csd_code.lower(),
            use_hall=True,
            lattice=ref_xtal.lattice,
            composition=composition,
            molecules=molecules,
            sites=sites,
            N_gen=args.ngen,
            N_pop=args.npop,
            N_cpu=args.nproc,
            cif="all.cif",
            skip_mlp=True,
            verbose=False,
            delta_length=args.delta_length,
            delta_angle=selected_deltas,
            matcher=matcher,
            check_stable=args.check_stable,
            soft_clash_check=args.soft_clash,
            soft_clash_buffer=soft_clash_buffer,
            close_grid_cutoff=args.close_grid_cutoff,
            close_grid_alpha=args.close_grid_alpha,
            max_mol_combos_per_grid=(
                args.max_mol_combos_per_grid
                if args.max_mol_combos_per_grid > 0
                else None
            ),
            max_valid_per_grid=(
                args.max_valid_per_grid
                if args.max_valid_per_grid > 0
                else None
            ),
            max_grid_product=max_grid_product,
            min_grid_points=(
                args.min_grid_per_gen if args.min_grid_per_gen > 0 else None
            ),
            fix_translation=args.fix_translation,
            asu_clamp=args.asu_clamp,
            gauge_site_index=args.gauge_site,
            order_identical_sites=args.order_identical_sites,
            N_min_matches=args.min_matches,
        )
        t0 = perf_counter()
        success_rate = qrs.run(ref_pmg=ref_pmg, max_rmsd=args.max_rmsd)
        time_cost_min = (perf_counter() - t0) / 60

        if success_rate is not None and success_rate > 0:
            print(f"Match found! Success rate: {success_rate}%")
        else:
            print("No match found within the given generations/population.")
        print(f"Time cost: {time_cost_min:.2f} min")
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
            time_min=time_cost_min,
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
                    f"{time_cost_min:.1f}",
                    f"{success_rate:.4f}" if success_rate is not None else "",
                    "success" if success_rate is not None and success_rate > 0 else "failure",
                    coverage if coverage is not None else "",
                    n_pregen_total,
                ]
            )

    print(f"\nSaved summary CSV: {csv_path}")
