# python -m unittest pyxtal/test_all.py
import unittest

import numpy as np
from pkg_resources import resource_filename
from pymatgen.core.structure import Molecule
import pymatgen.analysis.structure_matcher as sm
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.symmetry import Group, Wyckoff_position, get_wyckoffs
from pyxtal.wyckoff_site import WP_merge
from pyxtal.XRD import Similarity

cif_path = resource_filename("pyxtal", "database/cifs/aspirin.cif")
l0 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l1 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l2 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
l3 = Lattice.from_para(4.08, 7.13, 5.50, 90, 38, 90, ltype="monoclinic")
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_index(36, "4a")


class TestGroup(unittest.TestCase):
    def test_list_wyckoff_combinations(self):
        g = Group(64)
        a1, _ = g.list_wyckoff_combinations([4, 2])
        self.assertTrue(a1 is None)
        a2, _ = g.list_wyckoff_combinations([4, 8], quick=False) 
        self.assertTrue(len(a2) == 8)

class TestWP(unittest.TestCase):
    def test_wp(self):
        symbol = str(wp1.multiplicity) + wp1.letter
        self.assertTrue(symbol == "8b")
        symbol = str(wp2.multiplicity) + wp2.letter
        self.assertTrue(symbol == "4a")

    def test_merge(self):
        pt, wp, _ = WP_merge([0.05, 0.7, 0.24], l1.get_matrix(), wp1, 0.5)
        symbol = str(wp.multiplicity) + wp.letter
        self.assertTrue(symbol == "4a")
        pt, wp, _ = WP_merge([0.15, 0.7, 0.24], l1.get_matrix(), wp1, 0.5)
        symbol = str(wp.multiplicity) + wp.letter
        self.assertTrue(symbol == "8b")

    def test_get_wyckoff(self):
        for i in [1, 2, 229, 230]:
            get_wyckoffs(i)
            get_wyckoffs(i, organized=True)

    # to add test from string


class TestMolecular(unittest.TestCase):
    def test_single_specie(self):
        # print("test_h2o")
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=[["8b"]])
        struc.to_file()
        self.assertTrue(struc.valid)

        # test space group
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        # print(sga.get_space_group_symbol())
        self.assertTrue(sga.get_space_group_number() >= 36)
        # print(pmg_struc.frac_coords[:3])

        # test rotation
        ax = struc.mol_sites[0].orientation.axis
        struc.mol_sites[0].rotate(axis=[1, 0, 0], angle=90)
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        pmg_struc.to("cif", "1.cif")
        self.assertTrue(sga.get_space_group_symbol() == "Cmc2_1")
        # print(pmg_struc.frac_coords[:3])

    def test_sites(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [4])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol() == "Cmc2_1")

        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=[["4a", "4a"]])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol() == "Cmc2_1")

    def test_read(self):
        # test reading structure from external
        struc = pyxtal(molecular=True)
        struc.from_seed(seed=cif_path, molecule="aspirin")
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol() == "P2_1/c")
        # todo support reading ice

    def test_big_molecule(self):
        # print("test_big_molecule")
        for mol in ["ROY", "aspirin"]:
            struc = pyxtal(molecular=True)
            struc.from_random(3, 19, [mol], [4], 1.2)
            self.assertTrue(struc.valid)
            pair = struc.check_short_distances()
            if len(pair) > 0:
                print("short distances were detected")
                print(mol)
                print(pair)
            self.assertTrue(len(pair) == 0)

    def test_c60(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["C60"], [4], 1.0)
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        Li = Molecule(["Li"], [[0.0, 0.0, 0.0]])
        coords = [
            [0.000000, 0.000000, 0.000000],
            [1.200000, 1.200000, -1.200000],
            [1.200000, -1.200000, 1.200000],
            [-1.200000, 1.200000, 1.200000],
            [-1.200000, -1.200000, -1.200000],
        ]
        ps4 = Molecule(["P", "S", "S", "S", "S"], coords)

        for i in range(3):
            struc = pyxtal(molecular=True)
            struc.from_random(3, 10, [Li, ps4], [6, 2], 1.2, conventional=False)
            if struc.valid:
                self.assertTrue(len(struc.to_pymatgen()) == 16)

    def test_molecular_2d(self):
        # print("test_molecular_2d")
        struc = pyxtal(molecular=True)
        struc.from_random(2, 20, ["H2O"], [4], 1.0, conventional=False)
        cif = struc.to_file()
        self.assertTrue(struc.valid)

    def test_molecular_1d(self):
        struc = pyxtal(molecular=True)
        struc.from_random(1, 20, ["H2O"], [4], 1.0, conventional=False)
        cif = struc.to_file()
        self.assertTrue(struc.valid)
        # def test_space_groups(self):

    def test_preassigned_sites(self):
        sites = [["4a", "4a"]]
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=sites)
        self.assertTrue(struc.valid)

class TestAtomic3D(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(3, 225, ["C"], [4], 1.2, conventional=False)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(3, 99, ["Ba", "Ti", "O"], [1, 1, 3], 1.2)
        self.assertTrue(struc.valid)

    def test_preassigned_sites(self):
        sites = [["1b"], ["1b"], ["2c", "1b"]]
        struc = pyxtal()
        struc.from_random(3, 99, ["Ba", "Ti", "O"], [1, 1, 3], 1.0, sites=sites)
        self.assertTrue(struc.valid)

        struc = pyxtal()
        struc.from_random(3, 225, ["C"], [12], 1.0, sites=[["4a", "8c"]])
        self.assertTrue(struc.valid)

class TestAtomic2D(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(2, 20, ["C"], [4], 1.0, thickness=2.0)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(2, 4, ["Mo", "S"], [2, 4], 1.0)
        self.assertTrue(struc.valid)

class TestAtomic1D(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(1, 20, ["C"], [4], 1.0)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(1, 4, ["Mo", "S"], [2, 4], 1.0)
        self.assertTrue(struc.valid)


class TestCluster(unittest.TestCase):
    def test_multi_sites(self):
        struc = pyxtal()
        struc.from_random(0, 1, ["C"], [60], 1.0)
        self.assertTrue(struc.valid)

        struc = pyxtal()
        struc.from_random(0, 3, ["C"], [60], 1.0)
        self.assertTrue(struc.valid)

    def test_single_specie(self):
        struc = pyxtal()
        struc.from_random(0, "Ih", ["C"], [60], 1.0)
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = pyxtal()
        struc.from_random(0, 4, ["Mo", "S"], [2, 4], 1.0)
        self.assertTrue(struc.valid)


class TestLattice(unittest.TestCase):
    def test_para_matrix(self):
        self.assertTrue(np.allclose(l1.matrix, l2.matrix))

    def test_swap(self):
        l1.swap_axis(ids=[1, 0, 2])
        abc = l1.get_para()[:3]
        self.assertTrue(abc, np.array([9.13, 4.08, 5.50]))

    def test_optimize(self):
        l4, tran, _ = l3.optimize()
        self.assertTrue(abs(l4.beta-1.495907)<1e-4)

    def test_setpara(self):
        l0.set_para([5, 5, 5, 90, 90, 90])
        self.assertTrue(l0.a == 5)


class TestSymmetry(unittest.TestCase):
    def test_P21(self):
        strs = ["x, y, z", "-x, y+1/2, -z"]
        wyc, perm = Wyckoff_position.from_symops(strs)
        self.assertTrue(wyc.number == 4)

    def test_Pmn21(self):
        strs = ["x, y, z", "-x+1/2, -y, z+1/2", "-x, y, z", "x+1/2, -y, z+1/2"]
        wyc, perm = Wyckoff_position.from_symops(strs)
        self.assertTrue(wyc.number == 31)

    def test_P21a(self):
        strs = ["x, y, z", "-x, -y, -z", "-x+1/2, y+1/2, -z", "x+1/2, -y+1/2, z"]
        wyc, perm = Wyckoff_position.from_symops(strs)
        self.assertTrue(wyc.number == 14)

    def test_P21n(self):
        strs = [
            "x, y, z",
            "-x, -y, -z",
            "-x+1/2, y+1/2, -z+1/2",
            "x+1/2, -y+1/2, z+1/2",
        ]
        wyc, perm = Wyckoff_position.from_symops(strs)
        self.assertTrue(wyc.number == 14)

class TestSubgroup(unittest.TestCase):
    def test_cubic_cubic(self):
        sites = ['8a', '32e']
        numIons = int(sum([int(i[:-1]) for i in sites]))
        C1 = pyxtal()
        C1.from_random(3, 227, ['C'], [numIons], sites=[sites])
        pmg_s1 = C1.to_pymatgen()
        sga1 = SpacegroupAnalyzer(pmg_s1).get_space_group_symbol()

        C2s = C1.subgroup(eps=1e-4)
        for C2 in C2s:
            pmg_s2 = C2.to_pymatgen()
            sga2 = SpacegroupAnalyzer(pmg_s2).get_space_group_symbol()
            self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))

    def test_from_seed(self):
        from pymatgen import Lattice, Structure
        coords = [[0, 0, 0], [0.75,0.5,0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120,
                                  beta=90, gamma=60)
        struct = Structure(lattice, ["Si", "C"], coords)
        s1 = pyxtal()
        s1.from_seed(struct)
        s2 = s1.subgroup(once=True)
        pmg_s1 = s1.to_pymatgen()
        pmg_s2 = s2.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))
 
class TestPXRD(unittest.TestCase):
    def test_similarity(self):
        sites = ['8a']
        C1 = pyxtal()
        C1.from_random(3, 227, ['C'], [8], sites=[['8a']])
        xrd1 = C1.get_XRD()
        C2 = C1.subgroup(once=True, eps=1e-3)
        xrd2 = C1.get_XRD()
        p1 = xrd1.get_profile()
        p2 = xrd2.get_profile()
        s = Similarity(p1, p2, x_range=[15, 90])
        self.assertTrue( 0.9 <s.S <1.001)
     

        C2.apply_perturbation(1e-3, 1e-3)
        xrd3 = C2.get_XRD()
        p3 = xrd3.get_profile()
        s = Similarity(p1, p2, x_range=[15, 90])
        self.assertTrue( 0.95 <s.S <1.001)

class TestLoad(unittest.TestCase):
    def test_atomic(self):
        s1 = pyxtal()
        s1.from_random(3, 36, ['C', 'Si'], [4, 8])
        s2 = pyxtal()
        s2.load_dict(s1.save_dict())
        pmg_s1 = s1.to_pymatgen()
        pmg_s2 = s2.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))

    def test_molecular(self):
        s1 = pyxtal(molecular=True)
        s1.from_random(3, 36, ['H2O'], [4])
        s2 = pyxtal()
        s2.load_dict(s1.save_dict())
        pmg_s1 = s1.to_pymatgen()
        pmg_s2 = s2.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))
 

# class TestIO(unittest.TestCase):

if __name__ == "__main__":
    unittest.main()
