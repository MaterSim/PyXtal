"""
crystal plane class
"""

import itertools
from math import gcd

import numpy as np


def has_reduction(hkl):
    h, k, l = hkl
    gcf = gcd(h, gcd(k, l))
    if gcf > 1:
        # like [2, 2, 0]
        return True
    elif hkl[np.nonzero(np.array(hkl))[0][0]] < 0:
        # like [-2, 0, 0]
        return True
    return False


def reduced_hkl(hkl):
    h, k, l = hkl
    gcf = gcd(h, gcd(k, l))
    if gcf > 1:
        return [int(h / gcf), int(k / gcf), int(l / gcf)], gcf
    else:
        return [h, k, l], 1


def structure_factor(pos, hkl, total=True):
    coords = np.dot(pos, hkl)
    F = np.exp(-2 * np.pi * (1j) * coords)
    if total:
        return F.sum()
    else:
        return F


def get_dspacing(inv_matrix, hkl):
    return 1 / np.linalg.norm(inv_matrix.dot(np.array(hkl)))


class planes:
    """
    This is a database class to process crystal data

    Args:
        db_name: *.db format from ase database
        d_min: the minimum layer spacing
        cp_factor: skip non-close packed plane
    """

    def __init__(self, extent=6, d_min=1.5, cp_factor=0.5):
        self.d_min = d_min
        self.extent = extent
        self.cp_factor = cp_factor
        self.set_planes()
        # self.set_xtal()

    def set_xtal(self, xtal):
        self.xtal = xtal
        self.atoms = xtal.to_ase(center_only=True)
        self.cell_reciprocal = xtal.lattice.inv_matrix

    def set_planes(self):
        planes = list(itertools.product(
            range(-self.extent, self.extent + 1), repeat=3))
        planes = [hkl for hkl in planes if hkl !=
                  (0, 0, 0) and not has_reduction(hkl)]
        self.planes = planes

    def get_cp_factor(self, hkl):
        hkl, hkl_factor = reduced_hkl(hkl)
        hkl = np.array(hkl)
        dspacing = get_dspacing(self.cell_reciprocal, hkl)
        cp_factor = 0
        for n in range(1, 10):
            if dspacing / n > self.d_min:
                _factor = np.abs(self.get_structure_factor(
                    n * hkl)) / len(self.atoms)
                if _factor > cp_factor:
                    cp_factor = _factor
            else:
                break
        return cp_factor

    def search_close_packing_planes(self, N_max=10):
        """
        Search for the close-packed molecular plane for a given crystal

        Args:
            N_max: maximum number of multiples
        """
        cp_planes = []
        for hkl in self.planes:
            dspacing = get_dspacing(self.cell_reciprocal, hkl)
            for n in range(1, N_max):
                if dspacing / n > self.d_min:
                    hkl1 = n * np.array(hkl)
                    F = self.get_structure_factor(hkl1)
                    if np.abs(F) >= len(self.atoms) * self.cp_factor:
                        # Scan the plane with high density
                        plane = self.get_separation(hkl1)
                        if plane is not None:
                            cp_planes.append(plane)
                        break
                else:
                    # early termination for very short ranged hkl
                    break
        if len(cp_planes) > 0:
            cp_planes = sorted(cp_planes, key=lambda x: -x[-1][0])
        return cp_planes

    def get_structure_factor(self, hkl):
        return structure_factor(self.atoms.get_scaled_positions(), hkl)

    def get_separation(self, hkl):
        """
        Compute the separation for the given hkl plane

        Args:
            - hkl: three indices
        """
        hkl_reduced, hkl_factor = reduced_hkl(hkl)
        d_spacing = get_dspacing(self.cell_reciprocal, hkl_reduced)

        # Go over each molecule to compute the molecular slab
        slabs = []
        for mol_site in self.xtal.mol_sites:
            N_atoms = len(mol_site.numbers)
            center_frac = mol_site.position  # fractional
            xyz, species = mol_site._get_coords_and_species(unitcell=True)
            for id_mol, op in enumerate(mol_site.wp.ops):
                start, end = id_mol * N_atoms, (id_mol + 1) * N_atoms
                coords = xyz[start:end, :]  # frac
                center_frac = op.operate(mol_site.position)
                center_frac -= np.floor(center_frac)
                coords -= center_frac  # place the center to (0, 0, 0)
                center_hkl = np.dot(center_frac, hkl_reduced)
                center_hkl -= np.floor(center_hkl)
                coords_hkl = np.dot(coords, hkl_reduced)

                lower, upper = (
                    center_hkl + coords_hkl.min(),
                    center_hkl + coords_hkl.max(),
                )
                slabs.append([center_hkl, lower, upper])
        # if np.abs(hkl-np.array([0, 8, 0])).sum()==0:
        groups = self.group_slabs(slabs, 0.5 / hkl_factor)
        groups = self.group_slabs(groups, 0.5 / hkl_factor)
        # print(hkl); print(groups)
        separations = self.find_unique_separations(groups, d_spacing)
        return (hkl, d_spacing / abs(hkl_factor), separations)

    def group_slabs(self, slabs, tol):
        groups = []
        for slab in slabs:
            new = True
            center, lower, upper = slab
            for group in groups:
                # Center is within the group
                if (
                    group[1] <= center <= group[2]
                    or lower >= group[1]
                    and upper <= group[2]
                    or lower <= group[1]
                    and upper >= group[2]
                ):
                    new = False
                # to include
                elif lower - 1 <= group[1] and upper - 1 >= group[2] or lower - 1 >= group[1] and upper - 1 <= group[2]:
                    center, lower, upper = center - 1, lower - 1, upper - 1
                    new = False
                else:
                    dist = center - group[0]
                    shift = np.rint(dist)
                    if abs(dist - shift) < tol:
                        new = False
                        lower -= shift
                        upper -= shift

                if not new:
                    # Update group
                    if lower < group[1]:
                        group[1] = lower
                    if upper > group[2]:
                        group[2] = upper
                    center = (group[1] + group[2]) / 2
                    break
            if new:
                groups.append([center, lower, upper])
        return sorted(groups)

    def find_unique_separations(self, groups, d):
        groups.append(groups[0])
        groups[-1] = [g + 1 for g in groups[-1]]
        separations = []
        for i in range(len(groups) - 1):
            separations.append(d * (groups[i + 1][1] - groups[i][2]))
        separations = np.unique(-1 * np.array(separations).round(decimals=3))
        return -np.sort(separations)

    def gather(self, planes, tol=-0.1):
        output = []
        for _plane in planes:
            (hkl, d, separations) = _plane
            if separations[0] <= d:
                if separations[0] > tol:
                    for separation in separations:
                        if separation > tol:
                            p = plane(hkl, self.cell_reciprocal, separation)
                            print(p)
                            output.append(_plane)
                else:
                    break
        return output


class plane:
    """
    This simplest possible plane object
    """

    def __init__(self, hkl, cell_reciprocal, separation=None):
        self.indices = hkl
        self.cell_reciprocal = cell_reciprocal
        self.simple_indices, self.factor = reduced_hkl(hkl)
        self.d_spacing = get_dspacing(self.cell_reciprocal, self.indices)
        self.separation = separation

    def __str__(self):
        s = "({:2d} {:2d} {:2d})".format(*self.indices)
        s += f" Spacing: {self.d_spacing:6.3f}"
        s += f" Separation: {self.separation:6.3f}"
        return s


if __name__ == "__main__":

    import sys
    from pyxtal.db import database

    try:
        from ccdc import io
        from ccdc.particle import SlipPlanes

        csd = io.EntryReader("CSD")
    except:
        csd = None
        print("Cannot import ccdc to check")

    db = database("pyxtal/database/mech.db")

    if len(sys.argv) == 1:
        codes = ["ANLINB02", "ADIPAC", "CHEXDC"]
        for code in [c for c in db.codes if c not in codes]:
            # for code in codes:
            print("\n", code)
            p = planes()
            p.set_xtal(db.get_pyxtal(code))
            cp_planes = p.search_close_packing_planes()
            ps = p.gather(cp_planes, 0)

            # Crossvalidation with csd_python
            if csd is not None:
                splanes = SlipPlanes(csd.crystal(code))
                splanes = [
                    splane
                    for splane in splanes
                    if splane.repeat_distance > p.d_min and \
                    p.get_cp_factor(list(splane.orientation.hkl)) > p.cp_factor
                ]
                if len(splanes) > len(ps):
                    for splane in splanes:
                        print(
                            code,
                            splane.orientation,
                            round(splane.repeat_distance, 3),
                            round(splane.slab_separation, 3),
                        )
    else:
        data = [
            ("UCECAG01", [3, 0, -2]),
            ("PIDGOZ", [1, 0, 2]),
        ]
        p = planes(cp_factor=0.2)
        for code, hkl in data:
            print(hkl, tuple(hkl) in p.planes)
            p.set_xtal(db.get_pyxtal(code))
            print(code, hkl, p.get_separation(hkl))
