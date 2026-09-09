# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import numpy as np
import pymatgen.analysis.structure_matcher as sm
from ase.neighborlist import neighbor_list
from pymatgen.core import Structure

from pyxtal import pyxtal
from pyxtal.crystal import random_crystal
from pyxtal.lattice import Lattice
from pyxtal.symmetry import Hall, Wyckoff_position
from pyxtal.tolerance import Tol_matrix
from pyxtal.wyckoff_site import atom_site


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestDof(unittest.TestCase):
    def test_atomic(self):
        s = pyxtal()
        s.from_random(3, 225, ["C"], [8])
        ans = s.get_dof()
        assert s.lattice.dof == 1
        assert ans == 1


class TestAtomic3D(unittest.TestCase):

    def test_partial(self):
        cell = Lattice.from_para(7.8758, 7.9794, 5.6139, 90, 90, 90, ltype='orthorhombic')
        spg = 58
        elements = ['Al', 'Si', 'O']
        composition = [8, 4, 20]

        sites = [{"4e": [0.0000, 0.0000, 0.2418],
                  "4g": [0.1294, 0.6392, 0.0000],
                 },
                 {"4g": [0.2458, 0.2522, 0.0000]},
                 {"4g": [[0.4241, 0.3636, 0.0000], [0.5538, 0.2648, 0.0000]]},
                ]

        s = pyxtal()
        s.from_random(3, spg, elements, composition, lattice=cell, sites=sites)
        assert s.valid

    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(3, 225, ["C"], [4], 1.2, conventional=False)
        struc.to_file("tmp-3d.cif")
        os.remove("tmp-3d.cif")
        assert struc.valid


    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(3, 99, ["Ba", "Ti", "O"], [1, 1, 3], 1.2)
        assert struc.valid

        struc.from_random(3, 99, ["Ba", "Ti", "O"], [1, 1, 3], 1.2, use_asu=True)
        assert struc.valid

    def test_preassigned_sites(self):
        sites = [["1b"], ["1b"], ["2c", "1b"]]
        struc = pyxtal()
        struc.from_random(3, 99, ["Ba", "Ti", "O"], [1, 1, 3], 1.0, sites=sites)
        assert struc.valid

        struc = pyxtal()
        struc.from_random(3, 225, ["C"], [12], 1.0, sites=[["4a", "8c"]])
        assert struc.valid

        struc.from_random(3, 99, ["Ba", "Ti", "O"], [1, 1, 3], 1.0, sites=sites, use_asu=True)
        assert struc.valid

        struc.from_random(3, 225, ["C"], [12], 1.0, sites=[["4a", "8c"]], use_asu=True)
        assert struc.valid

    def test_read(self):
        # test reading xtal from cif
        for name in ["FAU", "NaSb3F10", "PVO", "lt_quartz"]:
            cif_file = cif_path + name + ".cif"
            pmg1 = Structure.from_file(cif_file, primitive=True)
            struc = pyxtal()
            struc.from_seed(seed=cif_file)
            pmg_struc = struc.to_pymatgen()
            assert sm.StructureMatcher().fit(pmg_struc, pmg1)

    def test_read_spglib(self):
        # test reading xtal from cif
        for name in ["FAU"]:
            cif_file = cif_path + name + ".cif"
            pmg1 = Structure.from_file(cif_file, primitive=True)
            struc = pyxtal()
            struc.from_seed(seed=cif_file, style="spglib")
            pmg_struc = struc.to_pymatgen()
            assert sm.StructureMatcher().fit(pmg_struc, pmg1)
        # more space groups
        for name in ["I41amd", "P4nmm", "Pmmn", "Pn3m", "Fd3", "Pn3"]:
            cif_file = cif_path + name + ".vasp"
            pmg1 = Structure.from_file(cif_file, primitive=True)
            struc = pyxtal()
            struc.from_seed(seed=cif_file, style="spglib")
            pmg_struc = struc.to_pymatgen()
            assert sm.StructureMatcher().fit(pmg_struc, pmg1)

    def test_read_by_HN(self):
        for name in ["aspirin"]:
            cif_file = cif_path + name + ".cif"
            pmg1 = Structure.from_file(cif_file, primitive=True)
            struc = pyxtal()
            for hn in Hall(14).hall_numbers:
                struc._from_pymatgen(pmg1, hn=hn)
                pmg_struc = struc.to_pymatgen()
                assert sm.StructureMatcher().fit(pmg_struc, pmg1)

    def test_from_tabular(self):
        xtal = pyxtal()
        rep = [116,10.5754,10.7203,4.47208,1.5705,2.6561,2.0943,1,0.4447,0.3762,0.7526]
        xtal.from_tabular_representation(rep, normalize=False)
        assert xtal.valid == False
        rep = [116,10.5754,10.7203,4.47208,1.5705,1.5705,1.5705,0,0.4447,0.3762,0.7526]
        xtal.from_tabular_representation(rep, normalize=False)
        assert xtal.valid == True

        rep0 = xtal.get_tabular_representation(discrete_cell=True,
                                               discrete=True,
                                               N_grids=100)
        assert(int(rep0[6]) == 50)

        xtal.from_tabular_representation(rep0,
                                         discrete_cell=True,
                                         discrete=True,
                                         N_grids=100)
        assert(int(xtal.lattice.get_para(degree=True)[-1]) == 90)

        reps = xtal.get_tabular_representations(N_wp=1,
                                                discrete_cell=True,
                                                discrete=True,
                                                N_grids=100)
        assert(len(reps)==8)

class TestDistanceTolerance(unittest.TestCase):
    """Every pair of atoms must clear the `Tol_matrix` entry of *that pair*.

    Regression test: `check_wp` used to be handed a single tolerance, the
    like-like one of the species being placed, and applied it to every pair.
    With `prototype="atomic"` the pair tolerance is `f * (r_A + r_B)`, so
    `f * 2 * r_new` is too small whenever the species being placed is the
    smaller of the two -- which lets the two overlap -- and too large whenever
    it is the larger, which rejects legal structures.
    """

    @staticmethod
    def worst_pair_ratio(struc, tm):
        """`min(d / tol(pair))` over pairs of distinct atoms, images included.

        Below 1.0 means at least one pair is closer than it was allowed to be.

        An atom against its own periodic image is excluded: those are governed
        by the cell, not by `check_wp`, and are not checked at all for an orbit
        of multiplicity 1 (`short_distances` has no pair to look at), so a
        lattice vector shorter than the like-like tolerance survives
        generation. That is a separate gap from the one this class covers.
        """
        atoms = struc.to_ase()
        numbers = atoms.numbers
        elements = sorted({int(n) for n in numbers})
        cutoff = max(tm.get_tol(a, b) for a in elements for b in elements)
        first, second, dist = neighbor_list("ijd", atoms, cutoff)
        distinct = first != second
        first, second, dist = first[distinct], second[distinct], dist[distinct]
        if len(dist) == 0:
            return np.inf
        tols = np.array([tm.get_tol(int(numbers[a]), int(numbers[b]))
                         for a, b in zip(first, second)])
        return float(np.min(dist / tols))

    def test_check_wp_rejects_a_contact_below_the_pair_tolerance(self):
        # Cs-O has to clear 2.04 A, O-O only 0.91 A. A Cs-O contact of 1.5 A
        # sits between the two: legal under O's like-like tolerance, illegal
        # under the pair's.
        tm = Tol_matrix(prototype="atomic", factor=1.3)
        assert tm.get_tol("O", "O") < 1.5 < tm.get_tol("Cs", "O")

        wp = Wyckoff_position.from_group_and_letter(1, "1a")
        cell = np.eye(3) * 10.0
        placed = atom_site(wp, [0.0, 0.0, 0.0], "Cs")
        candidate = atom_site(wp, [0.15, 0.0, 0.0], "O")

        class _Stub:
            tol_matrix = tm

        accepted = random_crystal.check_wp(
            _Stub(), [], [placed], cell, candidate, tm.get_tol("O", "O")
        )
        assert not accepted

    def test_generated_structures_honour_the_pair_tolerance(self):
        # Cs and O differ by 3.5x in covalent radius, and every site is a
        # general position, so a wrong tolerance shows up as a real overlap.
        # Both species orders are checked: the sites are placed one species at
        # a time, so a tolerance taken from the species being placed makes the
        # outcome depend on the order they are given in.
        tm = Tol_matrix(prototype="atomic", factor=1.3)
        for species, num_ions in ((["Cs", "O"], [1, 3]), (["O", "Cs"], [3, 1])):
            sites = [["1a"] * n for n in num_ions]
            for seed in range(10):
                struc = pyxtal()
                struc.from_random(3, 1, species, num_ions, sites=sites,
                                  tm=tm, random_state=seed)
                assert struc.valid
                ratio = self.worst_pair_ratio(struc, tm)
                assert ratio >= 1.0, (
                    f"{species}, seed {seed}: closest contact is {ratio:.3f} "
                    f"of the tolerance it was generated under"
                )


class TestAtomic2D(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(2, 20, ["C"], [4], 1.0, thickness=2.0)
        struc.to_file("tmp-2d.cif")
        os.remove("tmp-2d.cif")
        assert struc.valid

    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(2, 4, ["Mo", "S"], [2, 4], 1.0)
        assert struc.valid


class TestAtomic1D(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(1, 20, ["C"], [4], 1.0)
        struc.to_file("tmp-1d.cif")
        os.remove("tmp-1d.cif")
        assert struc.valid

    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(1, 4, ["Mo", "S"], [2, 4], 1.0)
        assert struc.valid


class TestCluster(unittest.TestCase):
    def test_multi_sites(self):
        struc = pyxtal()
        struc.from_random(0, 1, ["C"], [60], 1.0)
        assert struc.valid

        struc = pyxtal()
        struc.from_random(0, 3, ["C"], [60], 1.0)
        assert struc.valid

    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(0, "Ih", ["C"], [60], 1.0)
        assert struc.valid

    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(0, 4, ["Mo", "S"], [2, 4], 1.0)
        assert struc.valid

if __name__ == "__main__":
    unittest.main()
