#python -m unittest pyxtal/test_all.py
import unittest
import numpy as np
from pyxtal.crystal import *
from pyxtal.molecular_crystal import *
from pyxtal.lattice import Lattice
from pyxtal.symmetry import Wyckoff_position, get_wyckoffs
from pyxtal.wyckoff_site import WP_merge
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

l1 = Lattice.from_matrix([[4.08,0,0],[0,9.13,0],[0,0,5.50]])
l2 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_index(36, '4a')

class TestWP(unittest.TestCase):
    def test_wp(self):
        symbol = str(wp1.multiplicity) + wp1.letter
        self.assertTrue(symbol=='8b')
        symbol = str(wp2.multiplicity) + wp2.letter
        self.assertTrue(symbol=='4a')

    def test_merge(self):
        pt, wp, _ = WP_merge([0.05, 0.7, 0.24], l1.get_matrix(), wp1, 0.5)
        symbol = str(wp.multiplicity) + wp.letter
        self.assertTrue(symbol=='4a')
        pt, wp, _ = WP_merge([0.15, 0.7, 0.24], l1.get_matrix(), wp1, 0.5)
        symbol = str(wp.multiplicity) + wp.letter
        self.assertTrue(symbol=='8b')

    def test_get_wyckoff(self):
        for i in [1, 2, 229, 230]:
            get_wyckoffs(i)
            get_wyckoffs(i, organized=True)

    # to add test from string

class TestMolecular(unittest.TestCase):

    def test_single_specie(self):
        #print("test_h2o")
        struc = molecular_crystal(36, ['H2O'], [2], 1.0)
        struc.to_file()
        self.assertTrue(struc.valid)

        #test space group
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol()=='Cmc2_1')
        #print(pmg_struc.frac_coords[:3])

        #test rotation
        ax = struc.mol_sites[0].orientation.axis
        struc.mol_sites[0].rotate(axis=[1,0,0], angle=90)
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol()=='Cmc2_1')
        #print(pmg_struc.frac_coords[:3])

        #test reading structure from external
        mol = struc.molecules[0].mol
        struc = molecular_crystal(36, ['H2O'], [2], 1.0, pmg_struc)
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol()=='Cmc2_1')
        #print(pmg_struc.frac_coords[:3])

    def test_big_molecule(self):
        #print("test_big_molecule")
        for mol in ['ROY', 'aspirin']:
            struc = molecular_crystal(19, ['ROY'], [4], 1.2)
            self.assertTrue(struc.valid)
            pair = struc.check_short_distances()
            if len(pair) > 0:
                print(mol, pair)
            self.assertTrue(len(pair)==0)
 
    def test_c60(self):
        struc = molecular_crystal(36, ['C60'], [2], 1.0)
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        Li = Molecule(['Li'], [[0.0,0.0,0.0]])
        coords = [[ 0.000000,  0.000000,  0.000000],
                  [ 1.200000,  1.200000, -1.200000],
                  [ 1.200000, -1.200000,  1.200000],
                  [-1.200000,  1.200000,  1.200000],
                  [-1.200000, -1.200000, -1.200000]]
        ps4 = Molecule(['P', 'S', 'S', 'S', 'S'], coords)
    
        for i in range(3):
            struc = molecular_crystal(10, [Li, ps4], [6, 2], 1.2)
            if struc.valid:
                self.assertTrue(len(struc.to_pymatgen()) == 16)

    def test_molecular_2d(self):
        #print("test_molecular_2d")
        struc = molecular_crystal_2D(20, ['H2O'], [4], 1.0)
        cif = struc.to_file()
        self.assertTrue(struc.valid)
        
        struc = molecular_crystal_1D(20, ['H2O'], [4], 1.0)
        cif = struc.to_file()
        self.assertTrue(struc.valid)
        #def test_space_groups(self):


class TestAtomic3D(unittest.TestCase):

    def test_single_specie(self):
        struc = random_crystal(225, ['C'], [4], 1.2)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = random_crystal(99, ['Ba','Ti','O'], [1,1,3], 1.2)
        self.assertTrue(struc.valid)

    def test_preassigned_sites(self):
        sites=[["1b"], ["1b"], ["2c", "1b"]]
        struc = random_crystal(99, ['Ba','Ti','O'], [1,1,3], 1.0, sites=sites)
        self.assertTrue(struc.valid)

        struc = random_crystal(225, ['C'], [3], 1.0, sites=[["4a", "8c"]])
        self.assertTrue(struc.valid)

    #def test_space_groups(self):


class TestAtomic2D(unittest.TestCase):

    def test_single_specie(self):
        struc = random_crystal_2D(20, ['C'], [4], 1.0, thickness=2.0)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = random_crystal_2D(4, ['Mo','S'], [2,4], 1.0)
        self.assertTrue(struc.valid)

    #def test_space_groups(self):
class TestAtomic1D(unittest.TestCase):

    def test_single_specie(self):
        struc = random_crystal_1D(20, ['C'], [4], 1.0)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = random_crystal_1D(4, ['Mo','S'], [2,4], 1.0)
        self.assertTrue(struc.valid)

class TestCluster(unittest.TestCase):

    def test_multi_sites(self):
        struc = random_cluster(1, ['C'], [60], 1.0)
        struc.to_file()
        self.assertTrue(struc.valid)

        struc = random_cluster(3, ['C'], [60], 1.0)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_single_specie(self):
        struc = random_cluster('Ih', ['C'], [60], 1.0)
        struc.to_file()
        self.assertTrue(struc.valid)

    def test_mutiple_species(self):
        struc = random_cluster(4, ['Mo','S'], [2,4], 1.0)
        self.assertTrue(struc.valid)

class TestLattice(unittest.TestCase):

    def test_para_matrix(self):
        self.assertTrue(np.allclose(l1.matrix, l2.matrix))

    def test_swap(self):
        l1.swap_axis(ids=[1,0,2])
        abc = l1.get_para()[:3]
        self.assertTrue(abc, np.array([9.13, 4.08, 5.50]))

#class TestOperation(unittest.TestCase):
#class TestIO(unittest.TestCase):

if __name__ == '__main__':
    unittest.main()
