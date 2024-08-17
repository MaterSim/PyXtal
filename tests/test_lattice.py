# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import numpy as np
import pymatgen.analysis.structure_matcher as sm
from pymatgen.core import Structure

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.symmetry import Group, Hall, Wyckoff_position


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestLattice(unittest.TestCase):
    def test_para_matrix(self):
        assert np.allclose(l01.matrix, l02.matrix)

    def test_swap(self):
        l01.swap_axis(ids=[1, 0, 2])
        abc = l01.get_para()[:3]
        assert abc, np.array([9.13, 4.08, 5.5])

    def test_optimize_once(self):
        l3 = Lattice.from_para(4.08, 7.13, 5.50, 90, 38, 90, ltype="monoclinic")
        lat, tran, _ = l3.optimize_once()
        assert abs(lat.beta - 1.495907) < 0.0001

    def test_optimize_multi(self):
        l4 = Lattice.from_para(71.364, 9.127, 10.075, 90.00, 20.80, 90.00, ltype="monoclinic")
        lat, _ = l4.optimize_multi(7)
        assert abs(lat.beta - 1.7201) < 0.01

    def test_setpara(self):
        l0 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
        l0.set_para([5, 5, 5, 90, 90, 90])
        assert l0.a == 5

    def test_search_transformation(self):
        l6 = Lattice.from_para(3.454, 3.401, 5.908, 90.00, 105.80, 90.00, ltype="monoclinic")
        l7 = Lattice.from_para(6.028, 3.419, 6.028, 90.00, 146.92, 90.00, ltype="monoclinic")
        l7, _ = l7.optimize_multi()
        trans, diff = l7.search_transformation(l6)
        l7 = l7.transform_multi(trans)
        assert np.abs(l7.matrix - l6.matrix).sum() < 0.25

    def test_is_valid_lattice(self):
        l8 = Lattice.from_para(3.454, 3.401, 5.908, 90.00, 105.80, 91.00, ltype="monoclinic")
        l9 = Lattice.from_para(3.454, 3.401, 5.908, 90.00, 105.80, 90.00, ltype="monoclinic")
        l10 = Lattice.from_para(3.454, 3.401, 5.908, 90.00, 90.00, 90.00, ltype="cubic")
        assert not l8.is_valid_lattice()
        assert l9.is_valid_lattice()
        assert not l10.is_valid_lattice()

    def test_from_1d_representation(self):
        lat = Lattice.from_1d_representation([5.09, 6.11], "trigonal")
        assert abs(lat.a - 5.09) < 0.001
        assert abs(lat.c - 6.11) < 0.001
        assert abs(lat.gamma - 2 / 3 * np.pi) < 0.001


class TestOptimizeLattice(unittest.TestCase):
    def test_atomic(self):
        c1 = pyxtal()
        c1.from_seed(cif_path + "LiCs.cif", backend="pyxtal")
        pmg1 = c1.to_pymatgen()

        c2 = c1.copy()
        c2.optimize_lattice(1)
        pmg2 = c2.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg1, pmg2)

        c2.optimize_lattice(1)
        pmg2 = c2.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg1, pmg2)

        c3 = pyxtal()
        c3.from_seed(cif_path + "LiCs.cif")
        pmg3 = c3.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg1, pmg3)

    def test_molecular_from_seed(self):
        c1 = pyxtal(molecular=True)
        c1.from_seed(seed=cif_path + "aspirin-c.cif", molecules=["aspirin"])
        c1.mol_sites[0].get_mol_object(0)
        c1.mol_sites[0].get_mol_object(1)
        pmg1 = c1.to_pymatgen()

        c2 = c1.copy()
        c2.optimize_lattice(1)
        pmg2 = c2.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg1, pmg2)

    def test_molecular_trans(self):
        c1 = pyxtal(molecular=True)
        c1.from_seed(seed=cif_path + "aspirin.cif", molecules=["aspirin"])
        pmg1 = c1.to_pymatgen()
        c1.transform(trans=[[1, 0, 0], [0, 1, 0], [-1, 0, 1]])
        pmg2 = c1.to_pymatgen()
        c1.optimize_lattice()
        pmg3 = c1.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg1, pmg2)
        assert sm.StructureMatcher().fit(pmg2, pmg3)

    def test_molecular_nonstd(self):
        sgs = [5, 7, 8, 9, 12, 13, 14, 15]
        c1 = pyxtal(molecular=True)
        for sg in sgs:
            hns = Hall(sg).hall_numbers
            for hn in hns[1:]:
                c1.from_random(3, hn, ["aspirin"], use_hall=True)
                pmg1 = c1.to_pymatgen()
                c2 = c1.copy()
                c2.optimize_lattice(1)
                pmg2 = c2.to_pymatgen()
                assert sm.StructureMatcher().fit(pmg1, pmg2)
                pmg3 = Structure.from_str(c1.to_file(), fmt="cif")
                assert sm.StructureMatcher().fit(pmg1, pmg3)

    def test_molecular_from_random(self):
        rng = np.random.default_rng(0)
        sgs = [5, 7, 8, 12, 13, 14]
        c1 = pyxtal(molecular=True)
        for _i in range(20):
            sg = rng.choice(sgs)
            c1.from_random(3, sg, ["aspirin"])
            pmg1 = c1.to_pymatgen()
            c2 = c1.copy()
            c2.optimize_lattice(1)
            pmg2 = c2.to_pymatgen()
            assert sm.StructureMatcher().fit(pmg1, pmg2)

    def test_transform(self):
        l = Lattice.from_para(48.005, 7.320, 35.864, 90.000, 174.948, 90.000, ltype="monoclinic")
        sites = [
            [0.0000, 0.0000, 0.1250],
            [0.2296, 0.7704, 0.5370],
            [0.2296, 0.2296, 0.5370],
            [0.5408, 0.0000, 0.6186],
            [0.0000, 0.0000, 0.3074],
        ]
        c = pyxtal()
        c.build(9, ["S"], [8], lattice=l, sites=[sites])
        pmg0 = c.to_pymatgen()
        c.optimize_lattice()
        pmg1 = c.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg0, pmg1)

    def test_build(self):
        l = np.array([48.005, 7.320, 35.864, 90.000, 174.948, 90.000])
        sites = [
            [0.0000, 0.0000, 0.1250],
            [0.2296, 0.7704, 0.5370],
            [0.2296, 0.2296, 0.5370],
            [0.5408, 0.0000, 0.6186],
            [0.0000, 0.0000, 0.3074],
        ]
        c = pyxtal()
        c.build(9, ["S"], [8], lattice=l, sites=[sites])
        pmg0 = c.to_pymatgen()
        l1 = c.lattice.matrix

        c.optimize_lattice()
        pmg1 = c.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg0, pmg1)

        c1 = pyxtal()
        c1.build(9, ["S"], [8], lattice=l1, sites=[sites])
        c1.optimize_lattice()
        pmg2 = c1.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg0, pmg2)

    def test_build_1D(self):
        lat = Lattice.from_para(6.8472, 6.8472, 3.3198, 90, 90, 90, "tetragonal")
        c1 = pyxtal()
        sites = [
            [
                ("4g", 0.9236547993389047, 0.0, 0.5),
                ("4f", 0.39177300078207977, 0.0, 0.0),
            ]
        ]
        c1.build(30, ["C"], [8], lat, sites=sites, dim=1)
        assert c1.valid

    def test_build_2D(self):
        lat = Lattice.from_para(6.8472, 6.8472, 3.3198, 90, 90, 90, "tetragonal")
        c1 = pyxtal()
        sites = [[("4c", 0.15, 0.50, 0.28), ("4c", 0.59, 0.85, 0.87)]]
        c1.build(30, ["C"], [8], lat, sites=sites, dim=2)
        assert c1.valid

    # def test_build_0D(self):
    #    lat = Lattice.from_para(2.85, 2.85, 2.85, 90, 90, 90, 'spherical')
    #    c1 = pyxtal()
    #    sites = [[('8b', 0.70, 0.70, 0.70)]]
    #    c1.build(30, ['C'], [8], lat, sites=sites, dim=0)
    #    self.assertTrue(c1.valid)

    def test_transforms(self):
        paras = [
            (5.0317, 19.2982, 5.8004, 90.0000, 122.2672, 90.0000, "monoclinic", 6),
            (5.0317, 19.2982, 5.8004, 90.0000, 57.7328, 90.0000, "monoclinic", 6),
            (9.0640, 8.3522, 5.2856, 90.0000, 103.9699, 90.0000, "monoclinic", 3),
            (9.4913, 8.5844, 5.3358, 90.0000, 110.0035, 90.0000, "monoclinic", 3),
            (3.7727, 6.7490, 7.6446, 64.4392, 77.9087, 75.7214, "triclinic", 1),
            (3.7297, 6.5421, 7.9915, 73.9361, 71.5867, 103.2590, "triclinic", 1),
            (5.5740, 7.0902, 11.4529, 96.5703, 90.0224, 108.1857, "triclinic", 2),
            (5.5487, 7.1197, 11.4349, 83.0223, 89.8479, 72.0003, "triclinic", 2),
            (5.7985, 30.6352, 7.6374, 90.0000, 112.2615, 90.0000, "monoclinic", 3),
            (5.8280, 30.5992, 7.6373, 90.0000, 112.6020, 90.0000, "monoclinic", 3),
        ]

        for i in range(int(len(paras) / 2)):
            (a1, b1, c1, alpha1, beta1, gamma1, ltype1, N) = paras[i * 2]
            (a2, b2, c2, alpha2, beta2, gamma2, ltype2, N) = paras[i * 2 + 1]
            l1 = Lattice.from_para(a1, b1, c1, alpha1, beta1, gamma1, ltype=ltype1)
            l2 = Lattice.from_para(a2, b2, c2, alpha2, beta2, gamma2, ltype=ltype2)
            trans, diffs = l2.search_transformations(l1)
            print(i, len(trans), N)
            assert len(trans) == N

            # print("\nTarget:", l1)
            # print("Start :", l2)
            # if len(trans) == N:
            # for tran, diff in zip(trans, diffs):
            #    l = l2.transform_multi(tran)
            #    strs = "Success:" + str(l) + " {:6.3f} {:6.3f} {:6.3f}".format(*diff)
            #    print(strs)

    def test_optlat_setting(self):
        paras = [
            (18.950, 10.914, 31.672, 90, 168.63, 90, "monoclinic"),
            (48.005, 7.320, 35.864, 90, 174.948, 90, "monoclinic"),
        ]
        sites = [0.3600, 0.7500, 0.6743]
        for para in paras:
            (a, b, c, alpha, beta, gamma, ltype) = para
            l = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
            for sg in range(3, 16):
                mult = Group(sg, quick=True)
                c = pyxtal()
                c.build(sg, ["S"], [mult], lattice=l, sites=[[sites]])
                pmg0 = c.to_pymatgen()
                c0 = c.copy()
                c0.optimize_lattice(standard=True)
                pmg1 = c0.to_pymatgen()
                c1 = c.copy()
                c1.optimize_lattice(standard=False)
                pmg2 = c1.to_pymatgen()
                d1 = sm.StructureMatcher().get_rms_dist(pmg0, pmg1)
                d2 = sm.StructureMatcher().get_rms_dist(pmg0, pmg2)
                assert sum(d1) + sum(d2) < 0.001

if __name__ == "__main__":
    unittest.main()
