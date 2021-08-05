# python -m unittest pyxtal/test_all.py
import unittest

from random import choice
import numpy as np
from pkg_resources import resource_filename
from pymatgen.core import Structure
from pymatgen.core import Lattice as pmg_Lattice
from pymatgen.core.structure import Molecule
import pymatgen.analysis.structure_matcher as sm
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.molecule import pyxtal_molecule
from pyxtal.symmetry import Group, Wyckoff_position, get_wyckoffs
from pyxtal.XRD import Similarity
from pyxtal.operations import get_inverse

cif_path = resource_filename("pyxtal", "database/cifs/")
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

class TestOptLat(unittest.TestCase):
    def test_atomic(self):
        c1 = pyxtal()
        c1.from_seed(cif_path+"LiCs.cif", backend='pyxtal')
        pmg1 = c1.to_pymatgen()

        c2 = c1.copy()
        c2.optimize_lattice(1)
        pmg2 = c2.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg2))

        c2.optimize_lattice(1)
        pmg2 = c2.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg2))
        
        c3 = pyxtal()
        c3.from_seed(cif_path+"LiCs.cif")
        pmg3 = c3.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg3))

    def test_molecular(self):
        c1 = pyxtal(molecular=True)
        c1.from_seed(seed=cif_path+"aspirin-c.cif", molecules=["aspirin"])
        pmg1 = c1.to_pymatgen()

        c2 = c1.copy()
        c2.optimize_lattice(1)
        pmg2 = c2.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg2))

    def test_molecular_trans(self):
        c1 = pyxtal(molecular=True)
        c1.from_seed(seed=cif_path+"aspirin.cif", molecules=["aspirin"])
        pmg1 = c1.to_pymatgen()
        c1.transform(trans=[[1,0,0],[0,1,0],[-1,0,1]])
        pmg2 = c1.to_pymatgen()
        c1.optimize_lattice()
        pmg3 = c1.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg2))
        self.assertTrue(sm.StructureMatcher().fit(pmg2, pmg3))


    def test_molecular_diag(self):
        sgs, diags = [5, 7, 8, 12, 13, 14], [True, False]
        for diag in diags:
            for i in range(40):
                sg = choice(sgs)
                c1 = pyxtal(molecular=True)
                c1.from_random(3, sg, ["aspirin"], diag=diag) 
                pmg1 = c1.to_pymatgen()
                c2 = c1.copy()
                c2.optimize_lattice(1)
                pmg2 = c2.to_pymatgen()
                self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg2))

    def test_molecular_nodiag(self):
        sgs, diag = [5, 7, 8, 12, 13, 14], False
        for i in range(40):
            sg = choice(sgs)
            c1 = pyxtal(molecular=True)
            c1.from_random(3, sg, ["aspirin"], diag=diag) 
            pmg1 = c1.to_pymatgen()
            c2 = c1.copy()
            c2.optimize_lattice(1)
            pmg2 = c2.to_pymatgen()
            self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg2))


class TestWP(unittest.TestCase):
    def test_wp(self):
        symbol = str(wp1.multiplicity) + wp1.letter
        self.assertTrue(symbol == "8b")
        symbol = str(wp2.multiplicity) + wp2.letter
        self.assertTrue(symbol == "4a")

    def test_merge(self):
        pt, wp, _ = wp1.merge([0.05, 0.7, 0.24], l1.get_matrix(), 0.5)
        symbol = str(wp.multiplicity) + wp.letter
        self.assertTrue(symbol == "4a")
        pt, wp, _ = wp1.merge([0.15, 0.7, 0.24], l1.get_matrix(), 0.5)
        symbol = str(wp.multiplicity) + wp.letter
        self.assertTrue(symbol == "8b")

    def test_get_wyckoff(self):
        for i in [1, 2, 229, 230]:
            get_wyckoffs(i)
            get_wyckoffs(i, organized=True)

class TestDof(unittest.TestCase):
    def test_atomic(self):
        s = pyxtal() 
        s.from_random(3, 225, ['C'], [8])
        ans = s.get_dof()
        self.assertTrue(s.lattice.dof == 1)
        self.assertTrue(ans == 1)

class TestMolecule(unittest.TestCase):
    def test_get_orientations_in_wp(self):
        m = pyxtal_molecule('Benzene')
        g = Group(61)
        self.assertTrue(len(m.get_orientations_in_wp(g[0])) == 1)
        self.assertTrue(len(m.get_orientations_in_wp(g[1])) == 1)
        self.assertTrue(len(m.get_orientations_in_wp(g[2])) == 1)

class TestMolecular(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], sites=[["8b"]])
        struc.to_file()
        self.assertTrue(struc.valid)

        # test space group
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        # print(sga.get_space_group_symbol())
        self.assertTrue(sga.get_space_group_number() >= 36)

        # test rotation
        ax = struc.mol_sites[0].orientation.axis
        struc.mol_sites[0].rotate(ax_vector=[1, 0, 0], angle=90)
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        #pmg_struc.to("cif", "1.cif")
        self.assertTrue(sga.get_space_group_symbol() == "Cmc2_1")
        # print(pmg_struc.frac_coords[:3])

    def test_sites(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 19, ["H2O"])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol() == "P2_12_12_1")

        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=[["4a", "4a"]])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol() == "Cmc2_1")

    def test_read(self):
        # test reading structure from external
        struc = pyxtal(molecular=True)
        struc.from_seed(seed=cif_path+"aspirin.cif", molecules=["aspirin"])
        self.assertTrue(struc.lattice.ltype == "monoclinic")
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        self.assertTrue(sga.get_space_group_symbol() == "P2_1/c")
        C = struc.subgroup_once(eps=0, H=4)
        pmg_s2 = C.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg_struc, pmg_s2))

    def test_big_molecule(self):
        # print("test_big_molecule")
        for mol in ["ROY", "aspirin"]:
            struc = pyxtal(molecular=True)
            struc.from_random(3, 19, [mol], factor=1.2)
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
        struc.from_random(2, 20, ["H2O"])
        self.assertTrue(struc.valid)

    def test_molecular_1d(self):
        struc = pyxtal(molecular=True)
        struc.from_random(1, 20, ["H2O"])
        self.assertTrue(struc.valid)

    def test_preassigned_sites(self):
        sites = [["4a", "4a"]]
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=sites)
        self.assertTrue(struc.valid)

    def test_special_sites(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 61, ["Benzene"], [4])
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

    def test_read(self):
        # test reading xtal from cif
        for name in ["FAU", "NaSb3F10", "PVO", "lt_quartz"]:
            cif_file = cif_path + name + ".cif"
            pmg1 = Structure.from_file(cif_file)
            struc = pyxtal()
            struc.from_seed(seed=cif_file)
            pmg_struc = struc.to_pymatgen()
            self.assertTrue(sm.StructureMatcher().fit(pmg_struc, pmg1))


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

        C3s = C1.subgroup(permutations={"C":"Si"}, H=216)

    def test_from_seed(self):
        coords = [[0, 0, 0], [0.75,0.5,0.75]]
        lattice = pmg_Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120,
                                  beta=90, gamma=60)
        struct = Structure(lattice, ["Si", "C"], coords)
        s1 = pyxtal()
        s1.from_seed(struct)
        s2 = s1.subgroup_once(eps=0)
        pmg_s1 = s1.to_pymatgen()
        pmg_s2 = s2.to_pymatgen()
        self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))
    
        #pmg_s1 = Structure.from_file(cif_path + "B28.vasp")
        #struc = pyxtal()
        #struc.from_seed(seed=cif_path + "B28.vasp")
        #pmg_s2 = struc.to_pymatgen()
        #self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))
        #permutation = {"B":"C"}
        #struc.subgroup_once(0.01, None, permutation, max_cell=2) 

    def test_molecules(self):
        for name in ["aspirin", "resorcinol", "coumarin", "HAHCOI", "xxvi",\
                     "WEXBOS", "MERQIM", "LAGNAL", "YICMOP", "LUFHAW", \
                     "JAPWIH", "AXOSOW01", "PAHYON01"]:
            cif = cif_path + name + ".cif"
            struc = pyxtal(molecular=True)
            struc.from_seed(seed=cif, molecules=[name])
            pmg_struc = struc.to_pymatgen()
            pmg_s1 = Structure.from_file(cif)
            self.assertTrue(sm.StructureMatcher().fit(pmg_struc, pmg_s1))

            Cs = struc.subgroup(eps=0, max_cell=1)
            for C in Cs:
                pmg_s2 = C.to_pymatgen()
                self.assertTrue(sm.StructureMatcher().fit(pmg_struc, pmg_s2))
 
    def test_hydrate(self):
        # glycine dihydrate
        cif = cif_path + "gdh.cif"
        struc = pyxtal(molecular=True)
        struc.from_seed(seed=cif, molecules=['Glycine-z', 'H2O'])
        pmg_struc = struc.to_pymatgen()
        pmg_s1 = Structure.from_file(cif)
        self.assertTrue(sm.StructureMatcher().fit(pmg_struc, pmg_s1))

    def test_special(self):
        cif = cif_path + '191.vasp'
        struc = pyxtal()
        struc.from_seed(seed=cif)
        for i in range(100):
            s = struc.subgroup_once(0.2, None, None, 't+k', 2)
 
class TestPXRD(unittest.TestCase):
    def test_similarity(self):
        sites = ['8a']
        C1 = pyxtal()
        C1.from_random(3, 227, ['C'], [8], sites=[['8a']])
        xrd1 = C1.get_XRD()
        C2 = C1.subgroup_once(eps=1e-3)
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

class Test_operations(unittest.TestCase):
    def test_inverse(self):
        coord0 = [0.35, 0.1, 0.4]
        coords = np.array([
			   [0.350,  0.100,  0.400],
			   [0.350,  0.100,  0.000],
			   [0.350,  0.100,  0.000],
			   [0.350,  0.000,  0.667],
			   [0.350,  0.000,  0.250],
			   [0.350,  0.350,  0.400],
			   [0.350,  0.350,  0.500],
			   [0.350,  0.350,  0.000],
			   [0.350,  0.350,  0.350],
			   [0.100,  0.100,  0.100],
			   [0.400,  0.400,  0.400],
			   [0.350,  0.000,  0.000],
			   [0.000,  0.100,  0.400],
			   [0.350,  0.000,  0.400],
			 ])
        xyzs = ['x,y,z',
                'x,y,0',
                'y,x,0',
                'x,0,2/3',
                '0,x,1/4',
                'x,x,z',
                'x,-x,1/2',
                '2x,x,0',
                '-2x,-0.5x,-x+1/4',
                '-2y,-0.5y,-y+1/4',
                '-2z,-0.5z,-z+1/4',
                '0,0,x',
                '-y/2+1/2,-z,0',
                '-z,-x/2+1/2,0',
                ]
        
        for i, xyz in enumerate(xyzs):
            op = SymmOp.from_xyz_string(xyz)
            inv_op = get_inverse(op)
            coord1 = op.operate(coord0)
            coord2 = inv_op.operate(coord1)
            self.assertTrue(np.allclose(coord2, coords[i], rtol=1e-2))
            #strs = "{:6.3f} {:6.3f} {:6.3f}".format(*coord0)
            #strs += "  {:12s}  ".format(op.as_xyz_string())
            #strs += "{:6.3f} {:6.3f} {:6.3f}".format(*coord1)
            #strs += "  {:12s}  ".format(inv_op.as_xyz_string())
            #strs += "{:6.3f} {:6.3f} {:6.3f}".format(*coord2)
            #print(strs)

    def test_swap_wp(self):
        g = Group(38)
        wp = g[4]
        wp1, trans = wp.swap_axis([1,0,2])
        
        g = Group(71)
        wp = g[5]
        wp1, trans = wp.swap_axis([0,2,1])
        wp1, trans = wp.swap_axis([1,2,0])
        wp1, trans = wp.swap_axis([2,1,0])

    def test_alternative(self):
        for name in ["BTO-Amm2", "lt_quartz", "GeF2", "lt_cristobalite", "PVO"]:
            s = pyxtal()
            s.from_seed(cif_path+name+'.cif')
            pmg_s1 = s.to_pymatgen()
            strucs = s.get_alternatives()
            for struc in strucs:
                pmg_s2 = struc.to_pymatgen()
                #if not sm.StructureMatcher().fit(pmg_s1, pmg_s2):
                #    print(struc)
                #    print(s)
                self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))

    def test_wyc_sets(self):
        for i in range(1,229):
            res = Group(i).get_alternatives()['No.']

if __name__ == "__main__":
    unittest.main()
