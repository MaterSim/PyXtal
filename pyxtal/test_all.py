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
from pyxtal.supergroup import supergroups, supergroup

#import warnings
#warnings.filterwarnings("ignore", category=DeprecationWarning)

cif_path = resource_filename("pyxtal", "database/cifs/")
l0 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l1 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l2 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
l3 = Lattice.from_para(4.08, 7.13, 5.50, 90, 38, 90, ltype="monoclinic")
l4 = Lattice.from_para( 71.3649,   9.1273,  10.0753,  90.0000,  20.7978,  90.0000, ltype="monoclinic")
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_index(36, "4a")

wp6 = Wyckoff_position.from_group_and_index(167, 0)
l6 = Lattice.from_para(9.647000, 9.647000, 7.281000, 90, 90, 120, ltype="rhombohedral")


class TestGroup(unittest.TestCase):
    def test_list_wyckoff_combinations(self):
        g = Group(64)
        a1, _ = g.list_wyckoff_combinations([4, 2])
        self.assertTrue(a1 is None)
        a2, _ = g.list_wyckoff_combinations([4, 8], quick=False) 
        self.assertTrue(len(a2) == 8)

    def test_spg_symmetry(self):
        N_polar, N_centro, N_chiral = 0, 0, 0
        for sg in range(1, 231):
            g = Group(sg, quick=True)
            pg, polar, centro, chiral =g.point_group, g.polar, g.inversion, g.chiral
            #pg, polar, centro, chiral = get_point_group(sg)
            if polar:
                N_polar += 1
            if centro:
                N_centro += 1
            if chiral:
                N_chiral += 1
        
        self.assertTrue(N_polar == 68)
        self.assertTrue(N_centro == 92)
        self.assertTrue(N_chiral == 65)

    def test_ferroelectric(self):
        pairs = [(4, 1), (187, 4), (222, 5)]
        for pair in pairs:
            (sg, N) = pair
            self.assertTrue(len(Group(sg, quick=True).get_ferroelectric_groups()) == N)

    def test_check_compatible(self):
        self.assertTrue(Group(225).check_compatible([64, 28, 24]) == (True, True))
        self.assertTrue(Group(227).check_compatible([8]) == (True, False))
        self.assertTrue(Group(227).check_compatible([4]) == (False, False))
        self.assertTrue(Group(19).check_compatible([6]) == (False, False))

    def test_search_supergroup_paths(self):
        paths = Group(59, quick=True).search_supergroup_paths(139, 2) 
        self.assertTrue(paths == [[71, 139], [129, 139], [137, 139]])

    def test_get_splitters(self):
        s = pyxtal()
        s.from_seed(cif_path+'3-G139.cif')
        g = Group(225)
        solutions = g.get_splitters_from_structure(s, 't')
        self.assertTrue(len(solutions)==3)

class TestSupergroup(unittest.TestCase):
    def test_make_pyxtal(self):
        data = {
                "NbO2": 141,
                "GeF2": 62,
                "lt_quartz": 180,
                "NiS-Cm": 160, #9b->2a+4b
               }
        for cif in data.keys():
            s = pyxtal()
            s.from_seed(cif_path+cif+'.cif')
            my = supergroup(s, G=data[cif])
            sols = my.search_supergroup(max_solutions=1)
            for sol in sols:
                struc_high = my.make_pyxtal_in_supergroup(sol)
                strucs = my.make_pyxtals_in_subgroup(sol, 3) 
                pmg1 = struc_high.to_pymatgen()
                pmg2 = strucs[-1].to_pymatgen()
                rms = sm.StructureMatcher().get_rms_dist(pmg1, pmg2)[0]
                print(cif, rms)
                self.assertTrue(rms < 1e-3)

    def test_by_group(self):
        data = {
                "BTO-Amm2": 221,
                "BTO": 221,
                "lt_cristobalite": 227,
                "NaSb3F10": 194,
               }
        for cif in data.keys():
            s = pyxtal()
            s.from_seed(cif_path+cif+'.cif')
            sup = supergroups(s, G=data[cif], show=False)
            self.assertFalse(sup.strucs is None)

    def test_by_path(self):
        data = {
                "MPWO": [59, 71, 139, 225],
               }
        for cif in data.keys():
            s = pyxtal()
            s.from_seed(cif_path+cif+'.cif')
            sup = supergroups(s, path=data[cif], show=False)
            self.assertFalse(sup.strucs is None)

    def test_long(self):
        paras = (
                 ['I4_132', 98], #1-3
                 ["P3_112", 5], #1-3
                 ["P6_422", 21], #1-3
                )
        for para in paras:
            name, H = para
            name = cif_path + name + ".cif"
            s = pyxtal()
            s.from_seed(name)
            pmg_s1 = s.to_pymatgen()
            G = s.group.number
            struc_h = s.subgroup_once(eps=0.05, H=H, mut_lat=False)
            sup = supergroups(struc_h, G=G, d_tol=0.2, show=False, max_per_G=500)

            if sup.strucs is None: print("Problem in ", name)
            self.assertFalse(sup.strucs is None)

            match = False
            for struc in sup.strucs:
                pmg_g = struc.to_pymatgen()
                if sm.StructureMatcher().fit(pmg_g, pmg_s1):
                    match = True
                    break
            if not match: print("Problem in ", name)
            self.assertTrue(match)

    def test_long2(self):
        paras = (
                 ["P4_332", 155], #1-4 splitting
                 ["Fd3", 70], #1-3
                 ["Pm3", 47], #1-3
                 ["Fd3m", 166],
                 ["R-3c", 15],
                 ["R32", 5],
                 ["R-3", 147], #1-3, k
                 ["P4_332", 96], #1-3
                )
        for para in paras:
            name, H = para
            name = cif_path + name + ".cif"
            s = pyxtal()
            s.from_seed(name)
            pmg_s1 = s.to_pymatgen()
            G = s.group.number
            struc_h = s.subgroup_once(eps=0, H=H, mut_lat=False)
            my = supergroup(struc_h, G=G)
            sols = my.search_supergroup(max_solutions=1)

            if len(sols) == 0: print("problem", name)
            self.assertTrue(len(sols) > 0)

            struc_high = my.make_pyxtal_in_supergroup(sols[0])
            struc_sub = my.make_pyxtals_in_subgroup(sols[0], 3)[-1]

            pmg_high = struc_high.to_pymatgen()
            pmg_sub = struc_sub.to_pymatgen()

            dist1 = sm.StructureMatcher().get_rms_dist(pmg_high, pmg_sub)[0]
            dist2 = sm.StructureMatcher().get_rms_dist(pmg_s1, pmg_high)[0]
            self.assertTrue(dist1 < 1e-3)
            self.assertTrue(dist2 < 1e-3)

    def test_multi(self):
        data = {
                "BTO": [123, 221],
                "lt_cristobalite": [98, 210, 227],
                "BTO-Amm2": [65, 123, 221],
                "NaSb3F10": [186, 194],
                "NaSb3F10": [176, 194],
                #"MPWO": [59, 71, 139, 225],
               }
        for cif in data.keys():
            s = pyxtal()
            s.from_seed(cif_path+cif+'.cif')
            sup = supergroups(s, path=data[cif], show=False, max_per_G=2500)
            strucs = sup.get_transformation()
            pmg_0, pmg_1 = s.to_pymatgen(), sup.strucs[-1].to_pymatgen()
            pmg_2, pmg_3 = strucs[0].to_pymatgen(), strucs[1].to_pymatgen()
            dist1 = sm.StructureMatcher().get_rms_dist(pmg_0, pmg_2)[0]
            dist2 = sm.StructureMatcher().get_rms_dist(pmg_1, pmg_3)[0]
            self.assertTrue(dist1 < 1e-3)
            self.assertTrue(dist2 < 1e-3)

    def test_similarity(self):
        paras = [
                ('0-G62', '2-G71'),
                ('0-G62', '3-G139'),
                #('0-G62', '4-G225'),
               ]
        for para in paras:
            (cif1, cif2) = para
            s1 = pyxtal()
            s2 = pyxtal()
            s1.from_seed('pyxtal/database/cifs/'+cif1+'.cif')
            s2.from_seed('pyxtal/database/cifs/'+cif2+'.cif')
            pmg_s2 = s2.to_pymatgen()

            strucs, _, _ = s2.get_transition(s1)

            if strucs is None:
                print("Problem between ", cif1, cif2) 
            else:
                struc_G_in_H = strucs[-1]
                s3 = pyxtal()
                s3.from_seed(struc_G_in_H.to_pymatgen(), tol=1e-3)
                pmg_s3 = s3.to_pymatgen()
                dist = sm.StructureMatcher().get_rms_dist(pmg_s2, pmg_s3)[0]
                self.assertTrue(dist < 1e-3)
                self.assertTrue(s3.group.number == s2.group.number)

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
        c1.mol_sites[0].get_mol_object(0)
        c1.mol_sites[0].get_mol_object(1)
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
            for i in range(20):
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
        for i in range(20):
            sg = choice(sgs)
            c1 = pyxtal(molecular=True)
            c1.from_random(3, sg, ["aspirin"], diag=diag) 
            pmg1 = c1.to_pymatgen()
            c2 = c1.copy()
            c2.optimize_lattice(1)
            pmg2 = c2.to_pymatgen()
            self.assertTrue(sm.StructureMatcher().fit(pmg1, pmg2))


class TestWP(unittest.TestCase):
    def test_wp_site_symm(self):
        data = [(143, 1, '3 . .'), 
                (230, 6, '. 3 2'), 
                (160, 1, '. . m'),
                (160, 2, '3 m .')]
        for d in data:
            (sg, i, symbol) = d
            wp = Group(sg)[i]
            wp.get_site_symmetry()
            self.assertTrue(wp.site_symm == symbol)

    def test_wp_label(self):
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

        wp = Group(167)[0]
        cell = np.diag([9, 9, 7])
        for pt in [[0.12, 0, 0.25], [0, 0.1316, 0.25]]:
            _, wpt, _ = wp.merge(pt, cell, 0.1)
            symbol = str(wpt.multiplicity) + wpt.letter
            self.assertTrue(symbol == "18e")

        for pt in [[0, 0, 3/4], [2/3, 1/3, 7/12]]:
            _, wpt, _ = wp.merge(pt, cell, 0.1)
            symbol = str(wpt.multiplicity) + wpt.letter
            self.assertTrue(symbol == "6a")

    def test_search_generator(self):
        wp = Group(167)[1]
        for pt in [[0, 0.13, 0.25], [0.13, 0, 0.25]]:
            _, dist = wp.search_generator(pt)
            self.assertTrue(dist<1e-3)

    def test_get_wyckoff(self):
        for i in [1, 2, 229, 230]:
            get_wyckoffs(i)
            get_wyckoffs(i, organized=True)

    def test_is_equivalent(self):
        g = Group(15)
        wp = g[0]
        a = [ 0.10052793,  0.12726851,  0.27405404]
        b = [-0.10052642, -0.12726848, -0.27405526]
        c = [ 0.60052642,  0.62726848,  0.27405526]
        d = [-0.60052642, -0.62726848, -0.27405526]
        e = [0, 2.54537267e-01, 0]
        self.assertTrue(wp.is_equivalent(a,b))
        self.assertTrue(wp.is_equivalent(b,c))
        self.assertTrue(wp.is_equivalent(d,a))
        self.assertFalse(wp.is_equivalent(a,e))

        wp = g[1]
        a = [ 0.00,  0.127,  0.254]
        b = [-0.01, -0.127, -0.250]
        self.assertTrue(wp.is_equivalent(a,b))

    def test_euclidean(self):

        def check_error(spg, pt, cell):
            p0 = np.dot(pt, cell.matrix)
            for sg in spg:
                wp = Group(sg)[0]
                for i in range(len(wp)):
                    op0 = wp[i]
                    p1 = op0.operate(pt)
                    
                    op1 = wp.get_euclidean_generator(cell.matrix, i)
                    if wp.euclidean:
                        p2 = np.dot(op1.operate(p0), cell.inv_matrix)
                    else:
                        p2 = np.dot(op1.apply_rotation_only(p0), cell.inv_matrix) 
                        p2 += op1.translation_vector

                    diff = p1-p2
                    diff -= np.round(diff)
                    if np.linalg.norm(diff) > 0.02:
                        #res = '{:2d} {:28s}'.format(i, op0.as_xyz_string())
                        #res += ' {:28s}'.format(op1.as_xyz_string())
                        #res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p1)
                        #res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p2)
                        #res += '{:6.3f} {:6.3f} {:6.3f}'.format(*diff)
                        #print(res)
                        return False
            return True


        pt = [0.1333, 0.1496, 0.969]
        
        cell = Lattice.from_para(9.395, 7.395, 8.350, 91, 101, 92, ltype='triclinic')
        self.assertTrue(check_error(range(1, 3), pt, cell))
        
        cell = Lattice.from_para(9.395, 7.395, 8.350, 90, 101, 90, ltype='monoclinic')
        self.assertTrue(check_error(range(3, 16), pt, cell))
        
        cell = Lattice.from_para(9.395, 7.395, 8.350, 90, 90, 90, ltype='orthorhombic')
        self.assertTrue(check_error(range(16, 74), pt, cell))
        
        cell = Lattice.from_para(9.395, 9.395, 8.350, 90, 90, 90, ltype='tetragonal')
        self.assertTrue(check_error(range(74, 143), pt, cell))
        
        cell = Lattice.from_para(9.395, 9.395, 8.350, 90, 90, 120, ltype='hexagonal')
        self.assertTrue(check_error(range(143, 195), pt, cell))
        
        cell = Lattice.from_para(9.395, 9.395, 9.395, 90, 90, 90, ltype='cubic')
        self.assertTrue(check_error(range(195, 231), pt, cell))


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
            struc.from_random(3, 19, [mol], factor=1.4)
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
        lat, tran, _ = l3.optimize()
        self.assertTrue(abs(lat.beta-1.495907)<1e-4)

    def test_multi_optimize(self):
        lat0 = l4
        for i in range(10):
            lat, tran, opt = lat0.optimize()
            if opt:
                lat0 = lat
            else:
                break
        self.assertTrue(abs(lat.beta-1.7201)<1e-4)

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

class TestNeighbour(unittest.TestCase):
    def test_packing(self):
        c = pyxtal(molecular=True)
        for data in [("aspirin", 14),
                     ("WEXBOS", 14),
                     ("MERQIM", 12),
                     ("LAGNAL", 16),
                     ("YICMOP", 14),
                     ("LUFHAW", 18),
                     ("coumarin", 14),
                     ("HAHCOI", 14),
                     ("JAPWIH", 14),
                     ("AXOSOW01", 14),
                     ("PAHYON01", 13),
                     ("xxvi", 15),
                     ("resorcinol", 14),
                    ]:
 
            (name, CN) = data
            c.from_seed(seed=cif_path+name+".cif", molecules=[name])
            ds, _, _, _, engs = c.get_neighboring_molecules(0, 1.5)
            #print(engs)
            #print(name, CN, len(ds))
            self.assertTrue(len(ds) == CN)
 
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
        self.assertTrue( 0.9 <s.value <1.001)
     

        C2.apply_perturbation(1e-3, 1e-3)
        xrd3 = C2.get_XRD()
        p3 = xrd3.get_profile()
        s = Similarity(p1, p2, x_range=[15, 90])
        self.assertTrue( 0.95 <s.value <1.001)

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
                self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))

    def test_wyc_sets(self):
        for i in range(1, 229):
            res = Group(i, quick=True).get_alternatives()['No.']

if __name__ == "__main__":
    unittest.main()
