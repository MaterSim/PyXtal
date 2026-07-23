"""
Global Optimizer base on Quasi-Random Sampling
"""
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from threading import local
from time import time
from typing import TYPE_CHECKING
import multiprocessing as mp

import numpy as np
from scipy.stats import qmc
from pymatgen.analysis.structure_matcher import StructureMatcher

from pyxtal.optimize.base import GlobalOptimize
from pyxtal.representation import representation
from pyxtal.lattice import Lattice
from pyxtal.tolerance import Tol_matrix

if TYPE_CHECKING:
    from pyxtal.molecule import pyxtal_molecule

# Shared Tol_matrix for soft-clash checks (avoid rebuilding per call).
_MOL_TM = Tol_matrix(prototype="molecular")
# Conservative vdW reach beyond molecular radius for broad-phase pruning.
_VDW_REACH = 2.5


def identical_site_groups(composition):
    """
    Site-index groups that are chemically identical with Z'>=2.

    For ``composition=[2]`` returns ``[[0, 1]]``; for ``[1, 1]`` returns ``[]``
    (distinct species — no label-exchange symmetry).
    """
    groups = []
    site = 0
    for count in composition or []:
        n = int(count)
        if n >= 2:
            groups.append(list(range(site, site + n)))
        site += n
    return groups


def _leading_frac_xyz_indices(site_bounds, n_dof=None):
    """Indices of leading fractional xyz DOFs in one site's bounds."""
    limit = len(site_bounds) if n_dof is None else min(len(site_bounds), int(n_dof))
    idx = []
    for i in range(limit):
        lb, ub = site_bounds[i]
        if float(lb) >= -1e-8 and float(ub) <= 1.0 + 1e-8:
            idx.append(i)
            if len(idx) == 3:
                break
        else:
            break
    return idx


def _xyz_bounds_match(bounds_a, bounds_b, xyz_idx):
    """True when both sites share the same fractional xyz sampling box."""
    if len(xyz_idx) < 1:
        return False
    for i in xyz_idx:
        if i >= len(bounds_a) or i >= len(bounds_b):
            return False
        la, ua = bounds_a[i]
        lb, ub = bounds_b[i]
        if abs(float(la) - float(lb)) > 1e-12 or abs(float(ua) - float(ub)) > 1e-12:
            return False
    return True


def orderable_identical_site_pairs(composition, site_bounds_list):
    """
    Consecutive same-species site pairs with matching xyz bounds.

    ASU / gauge on only one site makes bounds differ; those pairs are skipped
    because label exchange is no longer a symmetry of the sampling domain.
    """
    pairs = []
    for group in identical_site_groups(composition):
        for a, b in zip(group, group[1:]):
            if a >= len(site_bounds_list) or b >= len(site_bounds_list):
                continue
            xyz_a = _leading_frac_xyz_indices(site_bounds_list[a])
            xyz_b = _leading_frac_xyz_indices(site_bounds_list[b])
            if xyz_a != xyz_b or not xyz_a:
                continue
            if _xyz_bounds_match(site_bounds_list[a], site_bounds_list[b], xyz_a):
                pairs.append((a, b, tuple(xyz_a)))
    return pairs


def _xyz_linear_bin_index(unit_vals, n_levels):
    """Linearize per-axis bin indices (row-major) for a shared xyz grid."""
    bins = [min(int(float(v) * int(n)), int(n) - 1) for v, n in zip(unit_vals, n_levels)]
    lin = 0
    for b, n in zip(bins, n_levels):
        lin = lin * int(n) + int(b)
    return lin


def identical_sites_xyz_out_of_order(
    sample_unit,
    seqs,
    site_bounds_list,
    composition,
    n_levels=None,
):
    """
    True if a same-species site pair violates translational ordering.

    For matching xyz grids, require strictly increasing linearized bin index
    (site ``j`` uses grid ``i+1`` or higher than site ``i``). Without discrete
    levels, fall back to lexicographic order on absolute fractional xyz.
    """
    pairs = orderable_identical_site_pairs(composition, site_bounds_list)
    if not pairs:
        return False

    sample_unit = np.asarray(sample_unit, dtype=float).ravel()
    offsets = []
    off = 0
    for seq in seqs:
        offsets.append(off)
        off += int(seq)

    for a, b, xyz_idx in pairs:
        oa, ob = offsets[a], offsets[b]
        ua = [sample_unit[oa + k] for k in xyz_idx]
        ub = [sample_unit[ob + k] for k in xyz_idx]

        use_bins = False
        if n_levels is not None:
            try:
                la = [int(n_levels[oa + k]) for k in xyz_idx]
                lb = [int(n_levels[ob + k]) for k in xyz_idx]
                use_bins = la == lb and all(n >= 1 for n in la)
            except Exception:
                use_bins = False

        if use_bins:
            if _xyz_linear_bin_index(ub, la) <= _xyz_linear_bin_index(ua, la):
                return True
            continue

        # Continuous fallback: compare absolute fractional coordinates.
        ba, bb = site_bounds_list[a], site_bounds_list[b]
        xa = tuple(
            float(ba[k][0]) + float(ua[i]) * (float(ba[k][1]) - float(ba[k][0]))
            for i, k in enumerate(xyz_idx)
        )
        xb = tuple(
            float(bb[k][0]) + float(ub[i]) * (float(bb[k][1]) - float(bb[k][0]))
            for i, k in enumerate(xyz_idx)
        )
        if xb <= xa:
            return True
    return False


def _site_type_indices(comp):
    """Map each molecular site index to its composition-pool index."""
    indices = []
    for type_idx, count in enumerate(comp):
        indices.extend([type_idx] * int(count))
    return indices


def _pool_max_radii(mol_pools):
    """Largest molecular radius per conformer pool (for conservative broad-phase)."""
    return [max(mol.radius for mol in pool) for pool in mol_pools]


def _pool_min_radii(mol_pools):
    """Smallest molecular radius per conformer pool (for guaranteed-clash pruning)."""
    return [min(mol.radius for mol in pool) for pool in mol_pools]


def _site_radii_from_pools(comp, pool_radii):
    """Expand per-component pool radii to a per-site list."""
    site_types = _site_type_indices(comp)
    return [float(pool_radii[t]) for t in site_types]


def resolve_close_grid_cutoff(value):
    """
    Normalize close-grid cutoff settings.

    Returns:
        ``None`` (disabled), ``"auto"`` (radius-based), or a positive float
        (flat Angstrom threshold).
    """
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in ("", "0", "0.0", "off", "none", "false", "no"):
            return None
        if text == "auto":
            return "auto"
        value = float(text)
    value = float(value)
    if value <= 0:
        return None
    return value


def _grid_must_clash(
    sites,
    close_grid_cutoff=None,
    site_radii=None,
    close_grid_alpha=0.7,
):
    """
    True when Wyckoff site centers are very close and the conformer loop
    can be skipped.

    ``close_grid_cutoff``:
      - ``None`` / disabled: never skip
      - positive float: skip if any orbit-aware center distance is below
        that flat Angstrom threshold
      - ``"auto"``: skip if any pair has
        ``d_orbit < alpha * 0.5 * (r_i + r_j)`` using ``site_radii``
        (pool *min* radii). The 0.5 factor matches molecular packing
        size used elsewhere (``radius * 0.5``).
    """
    if close_grid_cutoff is None:
        return False
    if len(sites) < 2:
        return False

    if close_grid_cutoff == "auto":
        if not site_radii or len(site_radii) != len(sites):
            return False
        alpha = float(close_grid_alpha)
        if alpha <= 0:
            return False
        for i in range(len(sites)):
            for j in range(i + 1, len(sites)):
                # Envelope radii overestimate hard size; pack with half-radii.
                threshold = alpha * 0.5 * (
                    float(site_radii[i]) + float(site_radii[j])
                )
                if sites[i].min_orbit_center_distance(sites[j]) < threshold:
                    return True
        return False

    threshold = float(close_grid_cutoff)
    if threshold <= 0:
        return False
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            if sites[i].min_orbit_center_distance(sites[j]) < threshold:
                return True
    return False


def _inter_pairs_for_pools(sites, comp, pool_max_radii):
    """
    Site pairs that may need an inter-molecular clash check.

    Pairs omitted here are separated even for the bulkiest conformer in each
    pool (using full Wyckoff-orbit center distances), so
    ``short_dist_with_wp2`` can be skipped for every conformer at this grid
    point (exact w.r.t. inter-site soft clashes).
    """
    site_types = _site_type_indices(comp)
    pairs = []
    for i in range(len(sites)):
        ri = pool_max_radii[site_types[i]]
        for j in range(i + 1, len(sites)):
            rj = pool_max_radii[site_types[j]]
            if sites[i].min_orbit_center_distance(sites[j]) <= ri + rj + _VDW_REACH:
                pairs.append((i, j))
    return pairs


def resolve_max_mol_combos_per_grid(mol_pools, max_mol_combos):
    """
    Resolve conformer cap for grid filtering.

    ``None`` or ``<= 0`` evaluates the full conformer product.  A positive
    value subsamples up to that many deterministic combinations per grid.
    """
    if mol_pools is None:
        return None
    total = int(np.prod([len(pool) for pool in mol_pools]))
    if max_mol_combos is None or int(max_mol_combos) <= 0:
        return None
    return min(int(max_mol_combos), total)


def _format_cell_for_print(cell, precision=2):
    """Convert a cell representation to a string with fixed numeric precision."""
    values = []
    for value in cell:
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, float):
            values.append(f"{value:.{precision}f}")
        else:
            values.append(str(value))
    return "[" + ", ".join(values) + "]"

def _qmc_scale(sample, lb, ub):
    """
    Scale unit-cube samples to [lb, ub], allowing zero-span (gauge-fixed) axes.

    SciPy ``qmc.scale`` requires strict ``l_bounds < u_bounds``. Translational
    gauge pinning sets ``lb == ub`` for free axes, so widen those by a tiny
    epsilon; with ``n_levels=1`` the snapped sample is 0.5 and maps back to
    the pinned midpoint.
    """
    lo = np.asarray(lb, dtype=float)
    hi = np.asarray(ub, dtype=float)
    if lo.shape != hi.shape:
        raise ValueError(f"Bound length mismatch: {lo.shape} vs {hi.shape}")
    # Widen any non-increasing interval just enough for SciPy.
    hi = np.where(hi <= lo, lo + 1e-12, hi)
    return qmc.scale(np.asarray(sample, dtype=float), lo, hi)


def generate_qrs_cell(sampler, cell_bounds, ref_volume, ltype):
    """
    A routine to generate quasi random samples for lattice and wp
    """
    # Sample cell parameters
    min_vol, max_vol = 0.75*ref_volume, 2.5*ref_volume
    lb = [b[0] for b in cell_bounds]
    ub = [b[1] for b in cell_bounds]
    count = 0
    while True:
        count += 1
        sample = _qmc_scale(sampler.random(), lb, ub)[0].tolist()
        lat = Lattice.from_1d_representation(sample, ltype)
        if min_vol < lat.volume < max_vol:
            #print(sample, ltype)
            return sample
        if count == 1000:
            raise ValueError("Cannot generate valid cell with 1000 attempts")

def _has_soft_clash(xtal, tm=None, buffer=0.0, inter_pairs=None):
    """
    True if molecular packing has soft vdW clashes (Tol_matrix).

    Hard check ``check_short_distances(r≈0.9)`` only catches nearly overlapping
    nuclei. Soft clashes (e.g. ~1.2–2.0 Å heavy contacts) still pass that cut,
    then explode in CHARMM (connectivity change → Wrong.xyz / calc.error).

    Args:
        inter_pairs: optional list of (i, j) site indices needing inter-site
            checks.  When empty, all inter-site checks are skipped (intra-site
            checks still run).  When ``None``, every pair is checked.
    """
    if tm is None:
        tm = _MOL_TM
    sites = xtal.mol_sites
    if inter_pairs:
        for i, j in inter_pairs:
            if not sites[i].short_dist_with_wp2(sites[j], tm=tm, buffer=buffer):
                return True
    for site in sites:
        if site.short_dist(buffer=buffer):
            return True
    if inter_pairs is None:
        n = len(sites)
        for i in range(n):
            for j in range(i + 1, n):
                if not sites[i].short_dist_with_wp2(sites[j], tm=tm, buffer=buffer):
                    return True
    return False


def _iter_mol_combos(mol_pools, max_combos=None, rng=None):
    """
    Yield molecule combinations from per-type pools.

    Full Cartesian product when ``max_combos`` is None / exceeds the product size;
    otherwise sample up to ``max_combos`` unique combinations (much cheaper with
    large conformer pools).
    """
    from itertools import product

    if not mol_pools:
        return
    sizes = [len(p) for p in mol_pools]
    total = int(np.prod(sizes))
    if max_combos is None or max_combos <= 0 or total <= max_combos:
        yield from product(*mol_pools)
        return

    if rng is None:
        rng = np.random.default_rng(0)
    seen = set()
    n_need = int(max_combos)
    # Cap attempts so we do not loop forever on tiny pools after collisions
    max_tries = max(n_need * 40, n_need + 10)
    tries = 0
    while len(seen) < n_need and tries < max_tries:
        tries += 1
        idxs = tuple(int(rng.integers(0, s)) for s in sizes)
        if idxs in seen:
            continue
        seen.add(idxs)
        yield tuple(mol_pools[i][j] for i, j in enumerate(idxs))


class _MolecularGridEvaluator:
    """Picklable per-worker state for one grid point + conformer product."""

    def __init__(
        self,
        cell,
        comp,
        smiles,
        mol_pools,
        probe_molecules,
        comb_per_grid,
        resolved_max_combos,
        pool_max_radii,
        pool_min_radii,
        soft_clash_check,
        soft_clash_buffer,
        close_grid_cutoff,
        close_grid_alpha,
        d_tol,
        max_valid_per_grid,
        reuse_grid_template,
    ):
        self.cell = cell
        self.comp = comp
        self.smiles = smiles
        self.mol_pools = mol_pools
        self.probe_molecules = probe_molecules
        self.comb_per_grid = comb_per_grid
        self.resolved_max_combos = resolved_max_combos
        self.pool_max_radii = pool_max_radii
        # Guaranteed-clash pruning must use the *smallest* pool radii; max
        # radii over-prune grids that could still pack with compact conformers.
        self.site_radii = _site_radii_from_pools(comp, pool_min_radii)
        self.soft_clash_check = soft_clash_check
        self.soft_clash_buffer = soft_clash_buffer
        self.close_grid_cutoff = close_grid_cutoff
        self.close_grid_alpha = float(close_grid_alpha)
        self.d_tol = d_tol
        self.max_valid_per_grid = max_valid_per_grid
        self.reuse_grid_template = reuse_grid_template
        self._template_local = local()

    def _prepare_probe_xtal(self, rep):
        """Build once per worker, then update only grid-dependent site DOFs."""
        if not self.reuse_grid_template:
            return rep.to_pyxtal(
                composition=self.comp, molecules=self.probe_molecules
            )

        state = getattr(self._template_local, "state", None)
        if state is None:
            probe_xtal = rep.to_pyxtal(
                composition=self.comp, molecules=self.probe_molecules
            )
            base_wps = [
                probe_xtal.group[int(values[0])]
                for values in rep.x[1:]
            ]
            state = (probe_xtal, base_wps)
            self._template_local.state = state

        probe_xtal, base_wps = state
        lattice_matrix = probe_xtal.lattice.matrix
        for site, base_wp, values in zip(
            probe_xtal.mol_sites, base_wps, rep.x[1:]
        ):
            position, wp, _ = base_wp.merge(
                np.asarray(values[1:4]),
                lattice_matrix,
                0.01,
                group=probe_xtal.group,
            )
            site.position = position
            site.wp = wp
            site.PBC = wp.PBC
            site._invalidate_geom_cache()
            if len(values) >= 7:
                site.update_orientation(values[4:7])
        return probe_xtal

    def evaluate(self, item):
        """
        Evaluate the complete conformer product for one grid point.

        Each worker owns its reusable probe crystal, so separate grid points can
        safely run concurrently while sharing the read-only conformer pools.
        """
        rep, site_chunks = item
        result = {
            "xtals": [],
            "attempted": 0,
            "skipped_special": 0,
            "skipped_short": 0,
            "skipped_soft": 0,
            "skipped_inter_broad": 0,
            "skipped_close_grid": 0,
            "special_signatures": [],
        }
        probe_xtal = self._prepare_probe_xtal(rep)
        if probe_xtal.has_special_site():
            result["attempted"] = self.comb_per_grid
            result["skipped_special"] = self.comb_per_grid
            for site_id, site in enumerate(probe_xtal.mol_sites):
                if (
                    site.wp.index > 0
                    and site_id < len(site_chunks)
                    and len(site_chunks[site_id]) >= 3
                ):
                    signature = tuple(np.round(site_chunks[site_id][:3], 8))
                    result["special_signatures"].append((site_id, signature))
            return result

        if _grid_must_clash(
            probe_xtal.mol_sites,
            close_grid_cutoff=self.close_grid_cutoff,
            site_radii=self.site_radii,
            close_grid_alpha=self.close_grid_alpha,
        ):
            result["attempted"] = self.comb_per_grid
            result["skipped_close_grid"] = self.comb_per_grid
            result["skipped_soft"] = self.comb_per_grid
            return result

        inter_pairs = _inter_pairs_for_pools(
            probe_xtal.mol_sites, self.comp, self.pool_max_radii
        )
        if not inter_pairs:
            result["skipped_inter_broad"] = self.comb_per_grid

        flat_chunks = np.concatenate([np.asarray(chunk) for chunk in site_chunks])
        seed_values = np.asarray(np.round(flat_chunks * 1e8), dtype=np.uint32)
        local_rng = np.random.default_rng(np.random.SeedSequence(seed_values))
        for chosen_molecules in _iter_mol_combos(
            self.mol_pools,
            max_combos=self.resolved_max_combos,
            rng=local_rng,
        ):
            result["attempted"] += 1
            site_idx = 0
            for type_idx, count in enumerate(self.comp):
                mol = chosen_molecules[type_idx]
                for _ in range(int(count)):
                    # Give each site its own molecule instance. Sharing one
                    # pool object across Z'>=2 sites is unsafe: xtal.copy()
                    # preserves aliases, and CHARMM site.update() mutates .mol.
                    probe_xtal.mol_sites[site_idx].update_molecule(
                        mol if int(count) == 1 else mol.copy(),
                        same_species=True,
                    )
                    site_idx += 1

            if self.soft_clash_check:
                if _has_soft_clash(
                    probe_xtal,
                    buffer=self.soft_clash_buffer,
                    inter_pairs=inter_pairs,
                ):
                    result["skipped_soft"] += 1
                    continue
            elif len(probe_xtal.check_short_distances(r=self.d_tol)) > 0:
                result["skipped_short"] += 1
                continue

            # Materialize an independent crystal because probe_xtal is mutated.
            xtal = probe_xtal.copy()
            result["xtals"].append((xtal, "QRandom"))
            if (
                self.max_valid_per_grid is not None
                and len(result["xtals"]) >= self.max_valid_per_grid
            ):
                break
        return result


_molecular_grid_worker_state = None


def _init_molecular_grid_worker(state):
    global _molecular_grid_worker_state
    _molecular_grid_worker_state = state


def _molecular_grid_worker(item):
    return _molecular_grid_worker_state.evaluate(item)


def generate_qrs_xtals(
    cell,
    wp_bounds,
    N_pop,
    smiles,
    comp,
    sampler_wp=None,
    d_tol=0.9,
    molecules=None,
    max_grid_attempts=None,
    min_grid_points=None,
    attempt_log_interval=None,
    log_grid_rejections=False,
    soft_clash_check=True,
    soft_clash_buffer=0.0,
    close_grid_cutoff=None,
    close_grid_alpha=0.7,
    max_mol_combos_per_grid=None,
    max_valid_per_grid=None,
    grid_workers=1,
    grid_executor="process",
    reuse_grid_template=True,
    order_identical_sites=True,
):
    """
    Get the qrs xtal samples

    Args:
        cell (list): [hall, a, b, c]
        wp_bounds (list): [[wp0], [wp1], ...]
        N_pop (int): number of valid candidates
        smiles (list): []
        comp (list): [1]
        sampler_wp: sampler
        d_tol (float): short distance tolerance value
        molecules: optional conformer pools aligned with smiles
        max_grid_attempts: optional cap on number of grid points sampled
        min_grid_points: keep sampling until at least this many grid points
            are evaluated, even after ``N_pop`` valid candidates are found.
            Extra valids are buffered on the sampler for later generations.
            ``None`` / ``<= 0`` restores the default early-stop-on-fill behavior.
        attempt_log_interval: optional periodic logging interval for attempt status
        log_grid_rejections: whether to log precheck-rejected grids
        soft_clash_check: reject Tol_matrix soft clashes before FF (default True)
        soft_clash_buffer: additive adjustment in Angstrom to Tol_matrix
            cutoffs. Negative values loosen the check; positive values tighten it.
        close_grid_cutoff: skip the conformer loop when site centers are too
            close. ``None`` disables; a positive float is a flat Angstrom
            threshold; ``"auto"`` uses ``alpha * 0.5 * (r_i + r_j)`` with
            pool *min* radii (orbit-aware; 0.5 matches packing size).
        close_grid_alpha: scale factor for ``close_grid_cutoff="auto"``
            (default 0.7).
        max_mol_combos_per_grid: max conformer combinations tried per grid point.
            ``None`` / ``<= 0`` evaluates the full conformer product.
        max_valid_per_grid: stop trying more conformers at a grid point once this
            many packings pass filtering (``None`` / ``<= 0`` keeps all valid).
        grid_workers: number of independent grid points evaluated concurrently.
            Use ``grid_executor='process'`` (default) for CPU-bound molecular
            clash checks; ``'thread'`` only helps when workers are I/O bound.
        grid_executor: ``'process'``, ``'thread'``, or ``'none'`` to force serial
            evaluation even when ``grid_workers > 1``.
        reuse_grid_template: reuse one mutable probe crystal per grid worker
            instead of rebuilding symmetry/lattice/site objects for every grid.
        order_identical_sites: for chemically identical Z'>=2 sites with the
            same xyz bounds, require strictly increasing translational grid
            index (or lex order) to remove label-exchange double counting.
    """
    #cell = [81, 11.38,  6.48, 11.24,  96.9]
    filter_t0 = time()
    xtals = []
    attempted = 0
    valid = 0
    skipped_special = 0
    skipped_short = 0
    skipped_soft = 0
    skipped_pattern = 0
    skipped_inter_broad = 0
    skipped_close_grid = 0
    skipped_order = 0
    max_combinations = None
    mol_pools = None
    if molecules is not None:
        mol_pools = []
        for mol_group in molecules:
            if isinstance(mol_group, (list, tuple)):
                mol_pools.append(list(mol_group))
            else:
                mol_pools.append([mol_group])
        max_combinations = int(np.prod([len(pool) for pool in mol_pools]))

    pool_min_radii = None
    if mol_pools is None:
        trimmed_wp_bounds = wp_bounds
    else:
        trimmed_wp_bounds = []
        site_idx = 0
        for type_idx, count in enumerate(comp):
            pool = mol_pools[type_idx]
            torsion_count = len(pool[0].torsionlist) if pool and pool[0].torsionlist is not None else 0
            for _ in range(int(count)):
                site_bounds = wp_bounds[site_idx]
                if torsion_count > 0:
                    site_bounds = site_bounds[:-torsion_count]
                trimmed_wp_bounds.append(site_bounds)
                site_idx += 1
    lb, ub = [], []
    seqs = []
    for wp_bound in trimmed_wp_bounds:
        lb += [b[0] for b in wp_bound]
        ub += [b[1] for b in wp_bound]
        seqs.append(len(wp_bound))

    order_pairs = (
        orderable_identical_site_pairs(comp, trimmed_wp_bounds)
        if order_identical_sites
        else []
    )
    n_levels = getattr(sampler_wp, "n_levels", None) if sampler_wp is not None else None

    if sampler_wp is None:
        sampler_wp = qmc.Sobol(d=len(lb), scramble=False)

    if hasattr(sampler_wp, "total") and hasattr(sampler_wp, "current"):
        remaining_grid_points = max(int(sampler_wp.total) - int(sampler_wp.current), 0)
        grid_budget = remaining_grid_points
    else:
        # Fallback for samplers without finite-state metadata.
        grid_budget = 2 ** max(int(np.log2(N_pop)) + 3, 10)

    if max_grid_attempts is not None:
        grid_budget = min(grid_budget, int(max_grid_attempts))

    min_grid = int(min_grid_points) if min_grid_points is not None else 0
    if min_grid < 0:
        min_grid = 0
    # Ensure the floor cannot exceed the remaining budget.
    min_grid = min(min_grid, grid_budget)

    # A parallel batch may produce more valid candidates than one generation
    # needs. Keep that ordered overflow on the persistent sampler so later
    # generations do not skip already evaluated grid/conformer candidates.
    # When min_grid is set we force fresh exploration every generation, so
    # stale pending would only block coverage — drop it.
    pending_xtals = list(getattr(sampler_wp, "_qrs_pending_xtals", []))
    sampler_wp._qrs_pending_xtals = []
    if min_grid > 0:
        pending_xtals = []
    elif pending_xtals:
        take = min(N_pop, len(pending_xtals))
        xtals.extend(pending_xtals[:take])
        valid += take
        sampler_wp._qrs_pending_xtals = pending_xtals[take:]

    if mol_pools is None:
        comb_per_grid = 1
        pool_max_radii = None
        pool_min_radii = None
        resolved_max_combos = None
    else:
        resolved_max_combos = resolve_max_mol_combos_per_grid(
            mol_pools, max_mol_combos_per_grid
        )
        if resolved_max_combos is None or resolved_max_combos <= 0:
            comb_per_grid = max_combinations
        else:
            comb_per_grid = min(max_combinations, int(resolved_max_combos))
        pool_max_radii = _pool_max_radii(mol_pools)
        pool_min_radii = _pool_min_radii(mol_pools)
    probe_molecules = [pool[0] for pool in mol_pools] if mol_pools is not None else None
    if max_valid_per_grid is None or max_valid_per_grid <= 0:
        max_valid_per_grid = None
    else:
        max_valid_per_grid = int(max_valid_per_grid)

    log_attempts = (
        attempt_log_interval is not None
        and int(attempt_log_interval) > 0
    )
    if log_attempts:
        attempt_log_interval = int(attempt_log_interval)

    close_grid_cutoff = resolve_close_grid_cutoff(close_grid_cutoff)
    grid_evaluator = None
    if mol_pools is not None:
        grid_evaluator = _MolecularGridEvaluator(
            cell=cell,
            comp=comp,
            smiles=smiles,
            mol_pools=mol_pools,
            probe_molecules=probe_molecules,
            comb_per_grid=comb_per_grid,
            resolved_max_combos=resolved_max_combos,
            pool_max_radii=pool_max_radii,
            pool_min_radii=pool_min_radii,
            soft_clash_check=soft_clash_check,
            soft_clash_buffer=soft_clash_buffer,
            close_grid_cutoff=close_grid_cutoff,
            close_grid_alpha=close_grid_alpha,
            d_tol=d_tol,
            max_valid_per_grid=max_valid_per_grid,
            reuse_grid_template=reuse_grid_template,
        )

    sampled_grid_points = 0
    bad_center_signatures = [set() for _ in seqs]

    def _finalize_stats():
        generate_qrs_xtals.last_attempted = attempted
        generate_qrs_xtals.last_valid = valid
        generate_qrs_xtals.last_skipped_special = skipped_special
        generate_qrs_xtals.last_skipped_short = skipped_short
        generate_qrs_xtals.last_skipped_soft = skipped_soft
        generate_qrs_xtals.last_skipped_pattern = skipped_pattern
        generate_qrs_xtals.last_skipped_inter_broad = skipped_inter_broad
        generate_qrs_xtals.last_skipped_close_grid = skipped_close_grid
        generate_qrs_xtals.last_skipped_order = skipped_order
        generate_qrs_xtals.last_combos_per_grid = comb_per_grid
        generate_qrs_xtals.last_sampled_grid_points = sampled_grid_points
        generate_qrs_xtals.last_possible = sampled_grid_points * comb_per_grid
        generate_qrs_xtals.last_buffered = len(sampler_wp._qrs_pending_xtals)
        generate_qrs_xtals.last_filter_seconds = time() - filter_t0

    def _need_more_grid():
        """Continue while population is short or the min-grid floor is unmet."""
        if sampled_grid_points >= grid_budget:
            return False
        if len(xtals) < N_pop:
            return True
        return sampled_grid_points < min_grid

    def _close_grid_starved():
        """
        True when close-grid pruning rejected every attempt and no valids
        were found after a substantial sample. Without this, an over-aggressive
        cutoff walks the entire remaining grid budget (~1e9) with no output.
        """
        if close_grid_cutoff is None or valid > 0 or attempted <= 0:
            return False
        floor = max(min_grid, 1024) if min_grid > 0 else 1024
        if sampled_grid_points < floor:
            return False
        return skipped_close_grid >= attempted

    def _trim_to_population():
        """Keep N_pop for this generation; buffer or subsample the rest."""
        nonlocal xtals
        if len(xtals) <= N_pop:
            return
        if min_grid > 0:
            # Evenly subsample across the explored block so later grids are
            # not discarded after paying for the min-grid floor.
            idx = np.linspace(0, len(xtals) - 1, N_pop, dtype=int)
            # Unique in case linspace collapses on tiny overflows.
            idx = list(dict.fromkeys(int(i) for i in idx))
            if len(idx) < N_pop:
                extras = [i for i in range(len(xtals)) if i not in idx]
                idx.extend(extras[: N_pop - len(idx)])
            xtals = [xtals[i] for i in idx[:N_pop]]
        else:
            sampler_wp._qrs_pending_xtals.extend(xtals[N_pop:])
            xtals = xtals[:N_pop]

    # Pending-only fill: skip further sampling unless a min-grid floor is set.
    if len(xtals) >= N_pop and min_grid <= 0:
        _finalize_stats()
        return xtals

    def _make_grid_rep(sample_wp):
        """Convert one unit-cube sample to a representation and site chunks."""
        scaled = _qmc_scale(np.asarray(sample_wp).reshape(1, -1), lb, ub)[0].tolist()
        x = [cell]
        prev = 0
        chunks = []
        for seq in seqs:
            chunk = scaled[prev:prev + seq]
            chunks.append(chunk)
            x.append([0] + chunk + [0])
            prev += seq
        return representation(x, smiles), chunks

    def _reject_unordered_identical_sites(sample_wp):
        """Skip swap-duplicate grids for chemically identical sites."""
        nonlocal attempted, skipped_order
        if not order_pairs:
            return False
        if identical_sites_xyz_out_of_order(
            sample_wp,
            seqs,
            trimmed_wp_bounds,
            comp,
            n_levels=n_levels,
        ):
            attempted += comb_per_grid if mol_pools is not None else 1
            skipped_order += comb_per_grid if mol_pools is not None else 1
            return True
        return False

    # Molecular conformer pools are the expensive path. Evaluate independent
    # grid points concurrently, but consume results in sampler order to retain
    # deterministic candidate ordering. Every submitted grid evaluates its full
    # conformer product, including the final grid of a generation.
    grid_workers = max(1, int(grid_workers))
    executor_kind = str(grid_executor).lower()
    if executor_kind not in ("process", "thread", "none"):
        executor_kind = "process"
    use_parallel = (
        mol_pools is not None
        and grid_workers > 1
        and executor_kind != "none"
    )

    def _consume_grid_results(results):
        nonlocal attempted, skipped_special, skipped_short, skipped_soft
        nonlocal skipped_inter_broad, skipped_close_grid, valid
        for result in results:
            attempted += result["attempted"]
            skipped_special += result["skipped_special"]
            skipped_short += result["skipped_short"]
            skipped_soft += result["skipped_soft"]
            skipped_inter_broad += result.get("skipped_inter_broad", 0)
            skipped_close_grid += result.get("skipped_close_grid", 0)
            valid += len(result["xtals"])
            xtals.extend(result["xtals"])
            for site_id, signature in result["special_signatures"]:
                bad_center_signatures[site_id].add(signature)

    if use_parallel:
        if executor_kind == "process":
            try:
                mp_ctx = mp.get_context("fork")
            except ValueError:
                mp_ctx = mp.get_context()
            pool_ctx = ProcessPoolExecutor(
                max_workers=grid_workers,
                initializer=_init_molecular_grid_worker,
                initargs=(grid_evaluator,),
                mp_context=mp_ctx,
            )
            map_fn = _molecular_grid_worker
        else:
            pool_ctx = ThreadPoolExecutor(max_workers=grid_workers)
            map_fn = grid_evaluator.evaluate

        close_grid_abort = False
        with pool_ctx as executor:
            while _need_more_grid():
                batch_size = min(grid_workers, grid_budget - sampled_grid_points)
                grids = []
                for _ in range(batch_size):
                    try:
                        sample_wp = sampler_wp.random()
                    except StopIteration:
                        break
                    sampled_grid_points += 1
                    if _reject_unordered_identical_sites(sample_wp[0]):
                        continue
                    rep, site_chunks = _make_grid_rep(sample_wp[0])

                    matched = False
                    for site_id, chunk in enumerate(site_chunks):
                        if len(chunk) >= 3:
                            signature = tuple(np.round(chunk[:3], 8))
                            if signature in bad_center_signatures[site_id]:
                                matched = True
                                break
                    if matched:
                        attempted += comb_per_grid
                        skipped_special += comb_per_grid
                        skipped_pattern += comb_per_grid
                    else:
                        grids.append((rep, site_chunks))

                if not grids:
                    if batch_size == 0:
                        break
                    continue

                _consume_grid_results(executor.map(map_fn, grids))

                if log_attempts:
                    print(
                        f"  Grid batch: sampled={sampled_grid_points}, "
                        f"attempted={attempted}, valid={valid}"
                    )
                if _close_grid_starved():
                    close_grid_abort = True
                    print(
                        f"QRS close-grid abort: sampled_grid={sampled_grid_points}, "
                        f"attempted={attempted}, valid=0, "
                        f"close_grid={skipped_close_grid} "
                        f"(cutoff={close_grid_cutoff!r}, alpha={close_grid_alpha})"
                    )
                    break

        _trim_to_population()
        _finalize_stats()
        if len(xtals) < N_pop and not close_grid_abort:
            print(
                f"QRS grid exhausted: sampled_grid={sampled_grid_points}, "
                f"attempted={attempted}, valid={valid}, "
                f"skipped_special={skipped_special}, skipped_short={skipped_short}, "
                f"soft={skipped_soft}, "
                f"pattern={skipped_pattern}, "
                f"inter_broad={skipped_inter_broad}, "
                f"close_grid={skipped_close_grid}, "
                f"order={skipped_order}"
            )
        return xtals

    while _need_more_grid():
        try:
            sample_wp = sampler_wp.random()
        except StopIteration:
            break
        sampled_grid_points += 1
        if _reject_unordered_identical_sites(sample_wp[0]):
            continue
        rep, site_chunks = _make_grid_rep(sample_wp[0])
        if mol_pools is None:
            attempted += 1
            xtal = rep.to_pyxtal(composition=comp)
            if xtal.has_special_site():
                skipped_special += 1
                continue
            # Soft vdW check dominates hard d_tol≈0.9 for organics; skip the
            # expensive to_pymatgen hard check when soft_clash is enabled.
            if soft_clash_check:
                if _has_soft_clash(xtal, buffer=soft_clash_buffer):
                    skipped_soft += 1
                    continue
            elif len(xtal.check_short_distances(r=d_tol)) > 0:
                skipped_short += 1
                continue

            valid += 1
            xtals.append((xtal, "QRandom"))
            if len(xtals) >= N_pop and sampled_grid_points >= min_grid:
                _trim_to_population()
                _finalize_stats()
                return xtals
        else:
            # Fast pattern pre-prune: if this site-center signature is known to
            # always collapse into a special WP, skip the entire conformer grid.
            matched_signature = None
            for site_id, chunk in enumerate(site_chunks):
                if len(chunk) < 3:
                    continue
                signature = tuple(np.round(chunk[:3], 8))
                if signature in bad_center_signatures[site_id]:
                    matched_signature = (site_id, signature)
                    break

            if matched_signature is not None:
                attempted += comb_per_grid
                skipped_special += comb_per_grid
                skipped_pattern += comb_per_grid
                if log_grid_rejections:
                    site_id, signature = matched_signature
                    print(
                        f"  Grid pattern-pruned: site={site_id}, center={signature}, "
                        f"skipped_combinations={comb_per_grid}"
                    )
                continue

            result = grid_evaluator.evaluate((rep, site_chunks))
            attempted += result["attempted"]
            skipped_special += result["skipped_special"]
            skipped_short += result["skipped_short"]
            skipped_soft += result["skipped_soft"]
            skipped_inter_broad += result.get("skipped_inter_broad", 0)
            skipped_close_grid += result.get("skipped_close_grid", 0)
            valid += len(result["xtals"])
            xtals.extend(result["xtals"])
            for site_id, signature in result["special_signatures"]:
                bad_center_signatures[site_id].add(signature)

            if log_attempts and attempted % attempt_log_interval < result["attempted"]:
                print(
                    f"  Grid complete: sampled={sampled_grid_points}, "
                    f"attempted={attempted}, valid={valid}"
                )
            if _close_grid_starved():
                print(
                    f"QRS close-grid abort: sampled_grid={sampled_grid_points}, "
                    f"attempted={attempted}, valid=0, "
                    f"close_grid={skipped_close_grid} "
                    f"(cutoff={close_grid_cutoff!r}, alpha={close_grid_alpha})"
                )
                _finalize_stats()
                return xtals
            if len(xtals) >= N_pop and sampled_grid_points >= min_grid:
                _trim_to_population()
                _finalize_stats()
                return xtals
        #else:
        #    print(rep, len(xtal.check_short_distances(r=0.6)))
    _trim_to_population()
    _finalize_stats()
    if len(xtals) < N_pop:
        print(
            f"QRS grid exhausted: sampled_grid={sampled_grid_points}, "
            f"attempted={attempted}, valid={valid}, "
            f"skipped_special={skipped_special}, skipped_short={skipped_short}, "
            f"skipped_soft={skipped_soft}, "
            f"skipped_pattern={skipped_pattern}, "
            f"skipped_inter_broad={skipped_inter_broad}, "
            f"skipped_close_grid={skipped_close_grid}, "
            f"order={skipped_order}"
        )
    return xtals


def expand_delta_angle_per_dof(delta_spec, composition, wp_bounds, default=30.0):
    """Expand per-component delta specs to a flat per-DOF angle resolution list.

    Each component entry may be:
      - ``None`` (monoatomic / no orientation DOFs; frac slots get ``default``)
      - a scalar (same resolution for all angle DOFs at that site)
      - a length-3 sequence ``[alpha, beta, gamma]`` for Euler angles
    """
    flat = []
    site_idx = 0
    comp_specs = delta_spec if isinstance(delta_spec, (list, tuple)) else [delta_spec]

    for comp_idx, count in enumerate(composition):
        if comp_idx < len(comp_specs):
            spec = comp_specs[comp_idx]
        elif comp_specs:
            spec = comp_specs[-1]
        else:
            spec = default

        for _ in range(int(count)):
            bounds = wp_bounds[site_idx]
            euler = None
            if spec is None:
                site_scalar = None
            elif isinstance(spec, (list, tuple)) and len(spec) == 3:
                euler = [float(d) for d in spec]
                site_scalar = None
            elif isinstance(spec, (list, tuple)) and len(spec) == 1:
                site_scalar = float(spec[0])
            else:
                site_scalar = float(spec)

            euler_idx = 0
            for lb, ub in bounds:
                span = ub - lb
                if abs(span - 1.0) < 1e-6:
                    flat.append(float(default))
                elif euler is not None:
                    da = euler[euler_idx] if euler_idx < len(euler) else euler[-1]
                    flat.append(float(da))
                    euler_idx += 1
                elif site_scalar is None:
                    # No orientation expected for this component; keep placeholder.
                    flat.append(float(default))
                else:
                    flat.append(float(site_scalar))
            site_idx += 1

    return flat


def _delta_angle_enabled(delta_angle, delta_length):
    """Return True when uneven-grid angle resolution should be used."""
    if delta_length > 0:
        return True
    if delta_angle is None:
        return False
    if isinstance(delta_angle, (int, float)):
        return delta_angle > 0
    if isinstance(delta_angle, (list, tuple)):
        def _positive(spec):
            if spec is None:
                return False
            if isinstance(spec, (list, tuple)):
                return any(_positive(x) for x in spec)
            return float(spec) > 0
        return any(_positive(d) for d in delta_angle)
    return False


def compute_wp_resolutions(wp_bounds, cell_lengths, delta_length=1.0, delta_angle=30.0):
    """
    Compute per-DOF grid resolution (number of levels) for an uneven QRS grid.

    Rules:
    - Fractional coordinate DOF (bounds within [0, 1], including ASU-shrunk
      intervals): n = max(1, int(edge * span / delta_length)).  Zero-span
      (gauge-fixed) axes get n = 1.
    - Angle DOF (any other bounds): n = max(1, round(span / delta_angle))

    Args:
        wp_bounds (list): per-site bounds from mol_site.get_bounds()
        cell_lengths (list): [a, b, c] in Angstroms
        delta_length (float): Angstrom resolution for coordinate dims
        delta_angle (float): degree resolution for angle dims

    Returns:
        list[int]: number of levels per DOF (flat, same ordering as wp_bounds)
    """
    n_levels = []
    total_dofs = sum(len(site_bounds) for site_bounds in wp_bounds)

    per_dof_delta = None
    if isinstance(delta_angle, (list, tuple, np.ndarray)):
        if len(delta_angle) == total_dofs:
            per_dof_delta = [None if d is None else float(d) for d in delta_angle]
        elif len(delta_angle) == len(wp_bounds):
            per_site_delta = [None if d is None else float(d) for d in delta_angle]
        else:
            try:
                first = delta_angle[0]
                fill = None if first is None else float(first)
                per_site_delta = [fill] * len(wp_bounds)
            except Exception:
                per_site_delta = [float(delta_angle)] * len(wp_bounds)
    else:
        per_site_delta = [float(delta_angle)] * len(wp_bounds)

    def _is_fractional_bound(lb, ub):
        """True for frac/ASU coords in [0, 1]; false for Euler/torsion angles."""
        return float(lb) >= -1e-8 and float(ub) <= 1.0 + 1e-8

    dof_idx = 0
    for site_idx, site_bounds in enumerate(wp_bounds):
        site_delta = per_site_delta[site_idx] if per_dof_delta is None else None
        coord_idx = 0
        for (lb, ub) in site_bounds:
            span = float(ub) - float(lb)
            if _is_fractional_bound(lb, ub):
                if span <= 1e-12:
                    n = 1  # gauge-fixed axis
                else:
                    edge = cell_lengths[coord_idx] if coord_idx < len(cell_lengths) else 1.0
                    n = max(1, int(edge * span / delta_length))
                coord_idx += 1
            else:  # angle DOF (Euler or torsion)
                if per_dof_delta is not None:
                    site_delta = per_dof_delta[dof_idx]
                if site_delta is None or float(site_delta) <= 0:
                    n = 1
                else:
                    n = max(1, int(round(span / float(site_delta))))
            dof_idx += 1
            n_levels.append(n)
    return n_levels


def apply_translational_gauge(
    wp_bounds,
    group,
    *,
    site_index=0,
    fix_free_axes=True,
    asu_clamp=True,
):
    """
    Restrict translational sampling on one general-position site.

    For molecular crystals on general Wyckoff sites (no special positions):

    - ``fix_free_axes``: pin continuous origin-free axes from
      ``group.get_free_axis()`` to the midpoint of that site's allowed
      interval (removes continuous gauge redundancy).
    - ``asu_clamp``: intersect that site's fractional bounds with the
      conventional ASU box from ``ASU_bounds`` (removes discrete space-group
      copies of the whole crystal).

    Args:
        wp_bounds: per-site bounds (mutated copy is returned)
        group: ``pyxtal.symmetry.Group`` for the crystal
        site_index: which molecular site to constrain (default: 0)
        fix_free_axes: enable continuous gauge pinning
        asu_clamp: enable ASU box restriction

    Returns:
        (new_wp_bounds, info_dict)
    """
    from pyxtal.constants import ASU_bounds

    bounds = [[list(b) for b in site] for site in wp_bounds]
    info = {
        "site_index": int(site_index),
        "sg_number": int(group.number),
        "free_axes": [],
        "fixed_axes": [],
        "asu_box": None,
        "asu_applied": False,
    }
    if not bounds or site_index < 0 or site_index >= len(bounds):
        return bounds, info

    site = bounds[site_index]
    # Leading DOFs with bounds in [0, 1] are fractional xyz (general position).
    frac_idx = [
        i for i, (lb, ub) in enumerate(site)
        if float(lb) >= -1e-8 and float(ub) <= 1.0 + 1e-8
    ]
    if len(frac_idx) < 1:
        return bounds, info

    # Map up to 3 fractional DOFs -> crystal axes 0,1,2 in encounter order.
    axis_for_dof = {frac_idx[k]: k for k in range(min(3, len(frac_idx)))}

    asu_lo = np.zeros(3)
    asu_hi = np.ones(3)
    if asu_clamp:
        box = np.asarray(ASU_bounds[group.number - 1], dtype=float)
        asu_lo = box[0::2]
        asu_hi = box[1::2]
        info["asu_box"] = [float(x) for x in box]
        for dof_i, axis in axis_for_dof.items():
            lb, ub = site[dof_i]
            site[dof_i] = [
                float(max(lb, asu_lo[axis])),
                float(min(ub, asu_hi[axis])),
            ]
            if site[dof_i][1] < site[dof_i][0]:
                # Degenerate / inconsistent: fall back to ASU segment.
                mid = 0.5 * (asu_lo[axis] + asu_hi[axis])
                site[dof_i] = [mid, mid]
        info["asu_applied"] = True

    free_axes = list(group.get_free_axis() or [])
    info["free_axes"] = [int(a) for a in free_axes]
    if fix_free_axes and free_axes:
        for dof_i, axis in axis_for_dof.items():
            if axis not in free_axes:
                continue
            lb, ub = site[dof_i]
            mid = 0.5 * (float(lb) + float(ub))
            site[dof_i] = [mid, mid]
            info["fixed_axes"].append(int(axis))

    bounds[site_index] = site
    return bounds, info


def budget_grid_levels(n_levels, max_product, min_level=1):
    """
    Reduce per-DOF grid levels so that prod(n_levels) <= max_product.

    Uses a uniform log-space scale first (preserves relative fineness), then
    decrements the currently largest remaining levels until the product fits.
    Levels never drop below ``min_level``.

    Args:
        n_levels: per-DOF bin counts
        max_product: maximum allowed product (None / <=0 disables budgeting)
        min_level: minimum bins per DOF (default 1)

    Returns:
        list[int]: budgeted levels (same length as input)
    """
    import math

    levels = [max(int(min_level), int(n)) for n in n_levels]
    if not levels or max_product is None or max_product <= 0:
        return levels

    def _log_prod(lv):
        return sum(math.log(x) for x in lv)

    log_budget = math.log(float(max_product))
    if _log_prod(levels) <= log_budget + 1e-12:
        return levels

    # Uniform scale in log space so relative resolution is roughly preserved.
    scale = math.exp((log_budget - _log_prod(levels)) / len(levels))
    levels = [max(int(min_level), int(round(n * scale))) for n in levels]

    # Trim residual overshoot from the largest dimensions first.
    while _log_prod(levels) > log_budget + 1e-12:
        candidates = [i for i, n in enumerate(levels) if n > min_level]
        if not candidates:
            break
        i_max = max(candidates, key=lambda j: levels[j])
        levels[i_max] -= 1

    return levels


def trim_wp_bounds_for_molecules(wp_bounds, composition, molecules):
    """Remove torsion DOFs from wp bounds when fixed conformer pools are provided."""
    if molecules is None:
        return wp_bounds

    mol_pools = []
    for mol_group in molecules:
        if isinstance(mol_group, (list, tuple)):
            mol_pools.append(list(mol_group))
        else:
            mol_pools.append([mol_group])

    trimmed_wp_bounds = []
    site_idx = 0
    for type_idx, count in enumerate(composition):
        pool = mol_pools[type_idx]
        torsion_count = len(pool[0].torsionlist) if pool and pool[0].torsionlist is not None else 0
        for _ in range(int(count)):
            site_bounds = wp_bounds[site_idx]
            if torsion_count > 0:
                site_bounds = site_bounds[:-torsion_count]
            trimmed_wp_bounds.append(site_bounds)
            site_idx += 1
    return trimmed_wp_bounds


class GridSampler:
    """
    Quantized quasi-random sampler for uneven discrete grids.

    Draws points from a Sobol sequence and snaps each coordinate to the centre
    of the nearest bin defined by n_levels (the number of grid levels per
    dimension).  This gives Sobol's multi-dimensional low-discrepancy coverage
    from the very first sample — every dimension varies immediately — while
    keeping all output values on a per-dimension discrete grid.

    Has the same .random() interface as qmc.Sobol, so it is a drop-in
    replacement.  self.current advances on every call; successive generation
    runs automatically draw the next block of the Sobol sequence.
    """

    def __init__(self, n_levels):
        self.n_levels = list(n_levels)
        self.d = len(n_levels)
        self.total = int(np.prod(n_levels)) if n_levels else 0
        self.current = 0
        self._sobol = qmc.Sobol(d=self.d, scramble=False)
        if self.d > 0 and self.total > 0:
            # Skip the all-0.5 Sobol point, which is a poor first seed for
            # symmetry-constrained molecular sites and can land on special positions.
            self._sobol.fast_forward(1)
            self.current = 1

    @property
    def exhausted(self):
        return self.current >= self.total

    def random(self):
        """Return the next quantized Sobol point as a (1, d) array in [0, 1]^d."""
        if self.exhausted:
            raise StopIteration("GridSampler exhausted")
        raw = self._sobol.random()[0]          # values in [0, 1)
        snapped = [(min(int(v * n), n - 1) + 0.5) / n
                   for v, n in zip(raw, self.n_levels)]
        self.current += 1
        return np.array([snapped])

    def __repr__(self):
        pct = 100.0 * self.current / self.total if self.total else 0.0
        return (f"GridSampler(d={self.d}, progress={self.current}/{self.total} [{pct:.1f}%])")


class QRS(GlobalOptimize):
    """
    Quasi-Random Sampling

    Args:
        smiles (str): smiles string
        workdir (str): path of working directory
        sg (int or list): space group number or list of spg numbers
        tag (string): job prefix
        ff_opt (bool): activate on the fly FF mode
        ff_style (str): automated force style (`gaff` or `openff`)
        ff_parameters (str or list): ff parameter xml file or list
        reference_file (str): path of reference xml data for FF training
        N_gen (int): number of generation (default: `10`)
        N_pop (int): number of populations (default: `10`)
        N_cpu (int): number of cpus for parallel calculation (default: `1`)
        cif (str): cif file name to store all structure information
        block: block mode
        num_block: list of blocks
        compositions: list of composition, (default is [1]*Num_mol)
        lattice (bool): whether or not supply the lattice
        torsions: list of torsion angle
        molecules (list): list of pyxtal_molecule objects
        sites (list): list of wp sites, e.g., [['4a']]
        use_hall (bool): whether or not use hall number (default: False)
        mlp (str): MACE | UMA | ANI
        skip_mlp (bool): whether or not use ani or not (default: True)
        output_mlp (bool): whether or not output the ANI relaxed structure (default: True)
        eng_cutoff (float): the cutoff energy for FF training
        E_max (float): maximum energy defined as an invalid structure
        verbose (bool): show more details
        use_mpi (bool): whether or not use mpi for parallel calculation
        delta_length (float): grid resolution in Å for fractional-coordinate DOF (0 = use Sobol)
        delta_angle (float): grid resolution in degrees for angle DOF (0 = use Sobol)
        max_grid_product (float or None): cap on prod(n_levels); levels are scaled down
            if the raw product exceeds this (default: None disables; set e.g. 1e9 to cap)
        N_min_matches (int): quit early after this many matches (default: 10)
    """

    def __init__(
        self,
        smiles: str,
        workdir: str,
        sg: int | list,
        tag: str = 'test',
        info: dict[any, any] | None = None,
        ff_opt: bool = False,
        ff_style: str = "openff",
        ff_parameters: str = "parameters.xml",
        reference_file: str = "references.xml",
        ref_criteria: dict[any, any] | None = None,
        N_gen: int = 10,
        N_pop: int = 10,
        N_cpu: int = 1,
        cif: str | None = None,
        block: list[any] | None = None,
        num_block: list[any] | None = None,
        composition: list[any] | None = None,
        lattice: Lattice | None = None,
        torsions: list[any] | None = None,
        molecules: list[pyxtal_molecule] | None = None,
        sites: list[any] | None = None,
        use_hall: bool = False,
        mlp: str = "MACE",
        skip_mlp: bool = True,
        output_mlp: bool = True,
        factor: float = 1.1,
        eng_cutoff: float = 5.0,
        E_max: float = 1e10,
        verbose: bool = False,
        max_time: float | None = None,
        matcher: StructureMatcher | None = None,
        early_quit: bool = False,
        check_stable: bool = True,
        use_mpi: bool = False,
        delta_length: float = 1.0,
        delta_angle: float = 60.0,
        max_grid_attempts: int | None = None,
        min_grid_points: int | None = None,
        max_grid_product: float | None = None,
        soft_clash_check: bool = True,
        soft_clash_buffer: float = 0.0,
        close_grid_cutoff: float | str | None = None,
        close_grid_alpha: float = 0.7,
        max_mol_combos_per_grid: int | None = None,
        max_valid_per_grid: int | None = None,
        fix_translation: bool | str = False,
        asu_clamp: bool = False,
        gauge_site_index: int = 0,
        order_identical_sites: bool = True,
        N_min_matches: int = 10,
    ):

        # POPULATION parameters:
        self.N_gen = N_gen # Number of lattice points
        self.N_pop = N_pop # Number of wp varieties
        self.verbose = verbose
        self.name = 'QRS'
        self.delta_length = delta_length
        self.delta_angle = delta_angle
        self.max_grid_attempts = max_grid_attempts
        self.min_grid_points = (
            int(min_grid_points)
            if min_grid_points is not None and int(min_grid_points) > 0
            else None
        )
        self.max_grid_product = max_grid_product
        self.soft_clash_check = bool(soft_clash_check)
        self.soft_clash_buffer = float(soft_clash_buffer)
        self.close_grid_cutoff = resolve_close_grid_cutoff(close_grid_cutoff)
        self.close_grid_alpha = float(close_grid_alpha)
        self.max_mol_combos_per_grid = max_mol_combos_per_grid
        self.max_valid_per_grid = max_valid_per_grid
        # Translational gauge / ASU: True/"auto" enables; False/"off" disables.
        if isinstance(fix_translation, str):
            self.fix_translation = fix_translation.strip().lower() not in (
                "0", "false", "off", "no", "none",
            )
        else:
            self.fix_translation = bool(fix_translation)
        self.asu_clamp = bool(asu_clamp)
        self.gauge_site_index = int(gauge_site_index)
        self.order_identical_sites = bool(order_identical_sites)
        self._gauge_info = None

        # initialize other base parameters
        GlobalOptimize.__init__(
            self,
            smiles,
            workdir,
            sg,
            tag,
            info,
            ff_opt,
            ff_style,
            ff_parameters,
            reference_file,
            ref_criteria,
            N_cpu,
            cif,
            block,
            num_block,
            composition,
            lattice,
            torsions,
            molecules,
            sites,
            use_hall,
            mlp,
            skip_mlp,
            output_mlp,
            factor,
            eng_cutoff,
            E_max,
            None, #random_state,
            max_time,
            matcher,
            early_quit,
            check_stable,
            use_mpi,
            N_min_matches=N_min_matches,
        )

        if self.rank == 0:
            strs = self.full_str()
            self.logging.info(strs)
            print(strs)

    def full_str(self):
        s = str(self)
        s += "\nMethod    : Deterministic Quasi-Random Sampling"
        s += f"\nGeneration: {self.N_gen:4d}"
        s += f"\nPopulation: {self.N_pop:4d}"
        if self.max_grid_attempts is not None:
            s += f"\nGrid cap   : {self.max_grid_attempts:4d}"
        if self.min_grid_points is not None:
            s += f"\nMin grid/gen: {self.min_grid_points:d}"
        if self.max_grid_product is not None and self.max_grid_product > 0:
            s += f"\nGrid budget: {self.max_grid_product:.3g}"
        if self.soft_clash_check:
            s += f"\nClash buffer: {self.soft_clash_buffer:+.3f} Angstrom"
        else:
            s += "\nClash filter: hard short-distance only (soft clash off)"
        if self.close_grid_cutoff == "auto":
            s += (
                f"\nClose-grid  : auto α={self.close_grid_alpha:.2f} "
                f"(orbit, min radii)"
            )
        elif self.close_grid_cutoff is not None:
            s += f"\nClose-grid  : {self.close_grid_cutoff:.1f} Angstrom (orbit)"
        if self.fix_translation or self.asu_clamp:
            s += (
                f"\nTrans gauge : fix_axes={self.fix_translation} "
                f"asu_clamp={self.asu_clamp} site={self.gauge_site_index}"
            )
        if self.order_identical_sites:
            s += "\nSite order : identical-xyz index j>i (same-species Z'>=2)"
        if self.molecules is not None:
            pools = [
                list(pool) if isinstance(pool, (list, tuple)) else [pool]
                for pool in self.molecules
            ]
            resolved = resolve_max_mol_combos_per_grid(
                pools, self.max_mol_combos_per_grid
            )
            total = int(np.prod([len(pool) for pool in pools]))
            if resolved is not None and resolved < total:
                s += f"\nCombos/grid: {resolved:d} of {total:d}"
            if self.max_valid_per_grid and self.max_valid_per_grid > 0:
                s += f"\nValid/grid cap: {int(self.max_valid_per_grid):d}"
        return s

    def _grid_workers(self):
        """Process-based grid parallelism; cap workers to limit fork/pickle overhead."""
        if self.molecules is None:
            return min(self.ncpu, 8)
        comb = int(np.prod([len(pool) for pool in self.molecules]))
        return min(self.ncpu, max(1, comb), 16)

    def _init_qrs_params(self):
        """
        Bootstrap hall_number, wp_bounds, ltype, and the Sobol sampler from a
        single randomly generated crystal using the fixed lattice.  Called once
        before the main generation loop when self.lattice is not None, so that
        generate_qrs_xtals can be used from gen 0 onwards.
        """
        from pyxtal import pyxtal
        from pyxtal.symmetry import Group
        from pyxtal.wyckoff_site import mol_site

        tmp = pyxtal(molecular=True)
        #print(self.smiles, self.sg, self.lattice); import sys; sys.exit()
        smiles = [smi + ".smi" for smi in self.smiles]
        sg = self.sg[0]
        wp = Group(sg, use_hall=True)[0] if self.use_hall else Group(sg)[0]
        mult = len(wp)
        numIons = [int(c * mult) for c in self.composition]

        for _ in range(200):
            tmp.from_random(3, sg, smiles, numIons,
                            lattice=self.lattice, use_hall=self.use_hall,
                            sites=self.sites, t_factor=0.2)  # skip structures with bad distances
            if tmp.valid:
                break
        self.hall_number = sg
        self.ltype = self.lattice.ltype
        self.wp_bounds = [site.get_bounds() for site in tmp.mol_sites]
        # Resolve international SG for ASU / free-axis tables.
        sg_group = Group(sg, use_hall=True) if self.use_hall else Group(sg)
        if self.fix_translation or self.asu_clamp:
            self.wp_bounds, self._gauge_info = apply_translational_gauge(
                self.wp_bounds,
                sg_group,
                site_index=self.gauge_site_index,
                fix_free_axes=self.fix_translation,
                asu_clamp=self.asu_clamp,
            )
            info = self._gauge_info
            parts = [f"SG {info['sg_number']}"]
            if info.get("asu_applied") and info.get("asu_box") is not None:
                box = info["asu_box"]
                parts.append(
                    f"ASU=[{box[0]:.3g},{box[1]:.3g}]x"
                    f"[{box[2]:.3g},{box[3]:.3g}]x"
                    f"[{box[4]:.3g},{box[5]:.3g}]"
                )
            if info.get("fixed_axes"):
                parts.append(f"fixed_axes={info['fixed_axes']}")
            elif info.get("free_axes") == []:
                parts.append("no_free_axes")
            self._log_info(
                f"Translational gauge on site {info['site_index']}: "
                + ", ".join(parts)
            )
        grid_wp_bounds = trim_wp_bounds_for_molecules(self.wp_bounds, self.composition, self.molecules)
        if _delta_angle_enabled(self.delta_angle, self.delta_length):
            # Uneven grid: per-dim resolution derived from cell lengths / angle range
            dl = self.delta_length if self.delta_length > 0 else 1.0
            da = self.delta_angle if self.delta_angle is not None else 30.0
            a, b, c = self.lattice.get_para()[:3]

            # If da is a per-component list (one value per composition entry),
            # expand it into a per-site or per-DOF list matching grid_wp_bounds.
            if isinstance(da, (list, tuple)) and len(da) == len(self.composition):
                has_euler_spec = any(
                    isinstance(entry, (list, tuple)) and len(entry) == 3
                    for entry in da
                )
                if has_euler_spec:
                    per_dof_da = expand_delta_angle_per_dof(
                        da, self.composition, grid_wp_bounds, default=30.0,
                    )
                else:
                    per_site_da = []
                    for comp_idx, cnt in enumerate(self.composition):
                        # None = monoatomic (no orientation); unused for angle bins.
                        val = da[comp_idx]
                        per_site_da.extend([val] * int(cnt))
                    per_dof_da = None
            else:
                per_site_da = da
                per_dof_da = None

            n_levels = compute_wp_resolutions(
                grid_wp_bounds,
                [a, b, c],
                dl,
                per_dof_da if per_dof_da is not None else per_site_da,
            )
            raw_product = int(np.prod(n_levels)) if n_levels else 0
            self._log_info(f"DOF grid levels: {n_levels} (product={raw_product})")
            if self.max_grid_product is not None and self.max_grid_product > 0:
                budgeted = budget_grid_levels(n_levels, self.max_grid_product)
                if budgeted != list(n_levels):
                    budgeted_product = int(np.prod(budgeted)) if budgeted else 0
                    self._log_info(
                        f"Grid levels to product <= {self.max_grid_product:.3g}: "
                        f"{budgeted} (product={budgeted_product})"
                    )
                    n_levels = budgeted
            self.sampler = GridSampler(n_levels)
        else:
            len_reps = sum(len(b) for b in grid_wp_bounds)
            self.sampler = qmc.Sobol(d=len_reps, scramble=False)
        #print(f"Bootstrap QRS parameters from random structure: hall={self.hall_number}, "
        #      f"ltype={self.ltype}, wp_bounds={self.wp_bounds}")

    def _log_info(self, msg):
        """Write to stdout and the workdir loginfo file (rank 0 only)."""
        if self.rank == 0:
            print(msg)
            self.logging.info(msg)

    def _run(self, pool=None):
        """
        The main code to run QRS prediction

        Returns:
            success_rate or None
        """
        self.ref_volumes = []
        success_rate = 0
        print(f"Rank {self.rank} starts QRS in {self.tag}")

        # When the lattice is fixed, bootstrap QRS parameters before gen 0 so
        # that the Sobol grid can be used for every generation (including gen 0).
        if self.rank == 0 and self.lattice is not None:
            self._init_qrs_params()

        for gen in range(self.N_gen):
            self.generation = gen
            cur_xtals = None

            if self.rank == 0:
                gen_hdr = f"\nGeneration {gen:d} starts"
                if self.lattice is not None and isinstance(self.sampler, GridSampler):
                    gen_hdr += f"  [{self.sampler}]"
                print(gen_hdr)
                self.logging.info(gen_hdr)
                t0 = time()

                if self.lattice is not None:
                    # Fixed-lattice mode: use QRS sampling for every generation.
                    # self.sampler is advanced progressively so each generation
                    # draws the next block of the Sobol sequence without repeats.
                    cell = [self.hall_number] + self.lattice.encode()
                    cur_xtals = generate_qrs_xtals(cell,
                                                   self.wp_bounds,
                                                   self.N_pop,
                                                   self.smiles,
                                                   self.composition,
                                                   self.sampler,
                                                   molecules=self.molecules,
                                                   max_grid_attempts=self.max_grid_attempts,
                                                   min_grid_points=self.min_grid_points,
                                                   attempt_log_interval=(50 if self.verbose else None),
                                                   log_grid_rejections=self.verbose,
                                                   soft_clash_check=self.soft_clash_check,
                                                   soft_clash_buffer=self.soft_clash_buffer,
                                                   close_grid_cutoff=self.close_grid_cutoff,
                                                   close_grid_alpha=self.close_grid_alpha,
                                                   max_mol_combos_per_grid=self.max_mol_combos_per_grid,
                                                   max_valid_per_grid=self.max_valid_per_grid,
                                                   grid_workers=self._grid_workers(),
                                                   order_identical_sites=self.order_identical_sites)
                    t_filter = time()
                    filter_sec = getattr(generate_qrs_xtals, "last_filter_seconds", t_filter - t0)
                    self._log_info(
                        f"Gen-{gen:d} filter time: {filter_sec:5.1f}s "
                        f"({len(cur_xtals)} candidates)"
                    )
                    #self._log_info(
                    #    f"Cell parameters in Gen-{gen:d}: "
                    #    f"{_format_cell_for_print(cell, precision=3)} {len(cur_xtals)}"
                    #)
                    if self.molecules is not None:
                        attempted = getattr(generate_qrs_xtals, "last_attempted", 0)
                        valid = getattr(generate_qrs_xtals, "last_valid", 0)
                        skipped_special = getattr(generate_qrs_xtals, "last_skipped_special", 0)
                        skipped_short = getattr(generate_qrs_xtals, "last_skipped_short", 0)
                        skipped_soft = getattr(generate_qrs_xtals, "last_skipped_soft", 0)
                        skipped_inter_broad = getattr(
                            generate_qrs_xtals, "last_skipped_inter_broad", 0
                        )
                        skipped_close_grid = getattr(
                            generate_qrs_xtals, "last_skipped_close_grid", 0
                        )
                        skipped_order = getattr(
                            generate_qrs_xtals, "last_skipped_order", 0
                        )
                        combos_per_grid = getattr(
                            generate_qrs_xtals, "last_combos_per_grid", None
                        )
                        sampled_grid = getattr(generate_qrs_xtals, "last_sampled_grid_points", 0)
                        possible = getattr(generate_qrs_xtals, "last_possible", None)
                        buffered = getattr(generate_qrs_xtals, "last_buffered", 0)
                        combo_str = (
                            f", combos/grid={combos_per_grid}"
                            if combos_per_grid is not None
                            else ""
                        )
                        broad_str = (
                            f", inter_broad={skipped_inter_broad}"
                            if skipped_inter_broad
                            else ""
                        )
                        close_str = (
                            f", close_grid={skipped_close_grid}"
                            if (skipped_close_grid or self.close_grid_cutoff is not None)
                            else ""
                        )
                        order_str = (
                            f", order={skipped_order}"
                            if (skipped_order or self.order_identical_sites)
                            else ""
                        )
                        #if possible is None:
                        self._log_info(
                            f"QRS in Gen-{gen:d}: "
                            f"attempt/valid/buffer/special/short/soft={attempted}/{valid}/{buffered}/{skipped_special}/{skipped_short}/{skipped_soft}"
                            f"{combo_str}{broad_str}{close_str}{order_str}"
                            f", sampled_grid={sampled_grid}"
                        )
                else:
                    # Variable-lattice mode: gen 0 uses random crystals to
                    # calibrate the cell bounds; gen > 0 uses QRS.
                    cur_xtals = [(None, "Random")] * self.N_pop
                    if gen > 0:
                        cell = generate_qrs_cell(self.sampler,
                                                 self.cell_bounds,
                                                 self.ref_volumes[-1],
                                                 self.ltype)
                        cell = [self.hall_number] + cell
                        cur_xtals = generate_qrs_xtals(cell,
                                                       self.wp_bounds,
                                                       self.N_pop,
                                                       self.smiles,
                                                       self.composition,
                                                       molecules=self.molecules,
                                                       max_grid_attempts=self.max_grid_attempts,
                                                       min_grid_points=self.min_grid_points,
                                                       attempt_log_interval=(50 if self.verbose else None),
                                                       log_grid_rejections=self.verbose,
                                                       soft_clash_check=self.soft_clash_check,
                                                       soft_clash_buffer=self.soft_clash_buffer,
                                                       close_grid_cutoff=self.close_grid_cutoff,
                                                       close_grid_alpha=self.close_grid_alpha,
                                                       max_mol_combos_per_grid=self.max_mol_combos_per_grid,
                                                       max_valid_per_grid=self.max_valid_per_grid,
                                                       grid_workers=self._grid_workers(),
                                                       order_identical_sites=self.order_identical_sites)
                        strs = f"Cell parameters in Gen-{gen:d}: "
                        print(strs, _format_cell_for_print(cell, precision=2), self.ref_volumes[-1], len(cur_xtals))
                        if self.molecules is not None:
                            attempted = getattr(generate_qrs_xtals, "last_attempted", 0)
                            valid = getattr(generate_qrs_xtals, "last_valid", 0)
                            skipped_special = getattr(generate_qrs_xtals, "last_skipped_special", 0)
                            skipped_short = getattr(generate_qrs_xtals, "last_skipped_short", 0)
                            skipped_soft = getattr(generate_qrs_xtals, "last_skipped_soft", 0)
                            skipped_inter_broad = getattr(
                                generate_qrs_xtals, "last_skipped_inter_broad", 0
                            )
                            skipped_close_grid = getattr(
                                generate_qrs_xtals, "last_skipped_close_grid", 0
                            )
                            skipped_order = getattr(
                                generate_qrs_xtals, "last_skipped_order", 0
                            )
                            combos_per_grid = getattr(
                                generate_qrs_xtals, "last_combos_per_grid", None
                            )
                            sampled_grid = getattr(generate_qrs_xtals, "last_sampled_grid_points", 0)
                            possible = getattr(generate_qrs_xtals, "last_possible", None)
                            buffered = getattr(generate_qrs_xtals, "last_buffered", 0)
                            combo_str = (
                                f", combos/grid={combos_per_grid}"
                                if combos_per_grid is not None
                                else ""
                            )
                            broad_str = (
                                f", skipped_inter_broad={skipped_inter_broad}"
                                if skipped_inter_broad
                                else ""
                            )
                            close_str = (
                                f", skipped_close_grid={skipped_close_grid}"
                                if skipped_close_grid
                                else ""
                            )
                            order_str = (
                                f", order={skipped_order}"
                                if (skipped_order or self.order_identical_sites)
                                else ""
                            )
                            if possible is None:
                                print(
                                    f"QRS trial stats in Gen-{gen:d}: sampled_grid={sampled_grid}, "
                                    f"attempted={attempted}, valid={valid}, buffered={buffered}, "
                                    f"skipped_special={skipped_special}, skipped_short={skipped_short}, "
                                    f"skipped_soft={skipped_soft}{combo_str}{broad_str}{close_str}{order_str}"
                                )
                            else:
                                print(
                                    f"QRS trial stats in Gen-{gen:d}: sampled_grid={sampled_grid}, "
                                    f"attempted={attempted}/{possible}, valid={valid}, buffered={buffered}, "
                                    f"skipped_special={skipped_special}, skipped_short={skipped_short}, "
                                    f"skipped_soft={skipped_soft}{combo_str}{broad_str}{close_str}{order_str}"
                                )


            # Broadcast
            if self.use_mpi: cur_xtals = self.comm.bcast(cur_xtals, root=0)

            # Fixed-lattice QRS: stop cleanly once the discrete grid is used up.
            # An empty population would crash gen_summary (N_pop argsort vs len 0).
            quit = False
            if self.rank == 0 and self.lattice is not None:
                n_valid = 0 if cur_xtals is None else len(cur_xtals)
                grid_done = (
                    isinstance(getattr(self, "sampler", None), GridSampler)
                    and self.sampler.exhausted
                )
                if n_valid == 0:
                    msg = (
                        f"Stopping QRS: no valid structures in Gen-{gen:d}"
                        + (" (grid exhausted)" if grid_done else "")
                    )
                    print(msg)
                    self.logging.info(msg)
                    if not hasattr(self, "N_struc"):
                        self.N_struc = 0
                    quit = True

            if self.use_mpi:
                quit = self.comm.bcast(quit, root=0)

            if quit:
                self.logging.info(f"Early Termination in Rank {self.rank}")
                return success_rate

            # Local optimization
            gen_results = self.local_optimization(cur_xtals, qrs=True, pool=pool)
            if self.rank == 0 and self.lattice is not None:
                t_relax = time()
                self._log_info(f"Gen-{gen:d} relax time: {t_relax - t_filter:5.1f}s")

            # Summary and Ranking
            quit = False
            if self.rank == 0:
                cur_xtals, matches, engs = self.gen_summary(t0,
                                        gen_results, cur_xtals)

                # update hist_best
                vols = []
                for id, (xtal, _) in enumerate(cur_xtals):
                    if xtal is not None:
                        vols.append(xtal.lattice.volume)

                # update best volume
                self.ref_volumes.append(np.array(vols).mean())
                if gen == 0 and self.lattice is None:
                    # Variable-lattice mode only: initialise cell and WP bounds
                    # from the best structure found in the random gen-0 batch.
                    best_xtal = cur_xtals[0][0]
                    self.cell_bounds = best_xtal.lattice.get_bounds(2.5, 25)
                    self.ltype = best_xtal.lattice.ltype
                    self.wp_bounds = [site.get_bounds() for site in best_xtal.mol_sites]
                    self.hall_number = best_xtal.group.hall_number
                    self.sampler = qmc.Sobol(d=len(self.cell_bounds), scramble=False)

                if self.ref_pmg is not None:
                    success_rate = self.success_count(cur_xtals, matches)

                    if self.early_termination(success_rate):
                        quit = True

                elif self.ref_pxrd is not None:
                    self.count_pxrd_match(cur_xtals, matches)

            if self.use_mpi:
                quit = self.comm.bcast(quit, root=0)
                self.comm.Barrier()
                self.logging.info(f"Gen {gen} Finish in Rank {self.rank}")

            if quit:
                self.logging.info(f"Early Termination in Rank {self.rank}")
                return success_rate

        return success_rate

if __name__ == "__main__":
    from pyxtal.representation import representation
    smiles = ["CC(=O)OC1=CC=CC=C1C(=O)O"]
    string = "82 11.43  6.49 11.19 83.31 1  0 0.77  0.57  0.53 48.55 24.31 145.9 -77.85 -4.40 170.9 0"
    rep3 = representation.from_string(string, smiles)
    xtal_ref = rep3.to_pyxtal()
    print(xtal_ref)

    # ── Reproducibility check when the lattice is fixed ───────────────────────
    # When a lattice is provided to QRS, the cell is frozen and only Wyckoff
    # positions are sampled via a Sobol sequence (scramble=False).
    # Because Sobol(scramble=False) is fully deterministic, two independent
    # invocations with the same cell / wp_bounds / N_pop must produce
    # exactly the same list of structures in the same order.
    print("\n--- Reproducibility check (fixed lattice) ---")
    cell = [xtal_ref.group.hall_number] + xtal_ref.lattice.encode()
    wp_bounds = [site.get_bounds() for site in xtal_ref.mol_sites]
    N_pop = 8
    comp = [1]

    xtals_a = generate_qrs_xtals(cell, wp_bounds, N_pop, smiles, comp,
                                  sampler_wp=qmc.Sobol(d=sum(len(b) for b in wp_bounds), scramble=False))
    xtals_b = generate_qrs_xtals(cell, wp_bounds, N_pop, smiles, comp,
                                  sampler_wp=qmc.Sobol(d=sum(len(b) for b in wp_bounds), scramble=False))

    assert len(xtals_a) == len(xtals_b), \
        f"Different number of valid structures: {len(xtals_a)} vs {len(xtals_b)}"

    all_match = True
    for i, ((xa, _), (xb, _)) in enumerate(zip(xtals_a, xtals_b)):
        ra = xa.get_1D_representation().to_string()
        rb = xb.get_1D_representation().to_string()
        if ra != rb:
            print(f"  Structure {i}: MISMATCH\n  run-1: {ra}\n  run-2: {rb}")
            all_match = False

    if all_match:
        print(f"  PASSED — both runs produced {len(xtals_a)} identical structures.")
