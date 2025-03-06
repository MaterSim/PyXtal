# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import numpy as np
import pymatgen.analysis.structure_matcher as sm
from pymatgen.core import Lattice as pmg_Lattice
from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.operations import get_inverse
from pyxtal.symmetry import Group, Hall, Wyckoff_position


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestNeighbour(unittest.TestCase):
    def test_packing(self):
        c = pyxtal(molecular=True)
        for data in [
            ("aspirin", 14),
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
            c.from_seed(seed=cif_path + name + ".cif", molecules=[name])
            ds, _, _, _, engs = c.get_neighboring_molecules(0, 1.5)
            # print(engs)
            # print(name, CN, len(ds))
            assert len(ds) == CN


class TestSubgroup(unittest.TestCase):
    def test_cubic_cubic(self):
        sites = ["8a", "32e"]
        numIons = int(sum([int(i[:-1]) for i in sites]))
        C1 = pyxtal()
        C1.from_random(3, 227, ["C"], [numIons], sites=[sites])
        pmg_s1 = C1.to_pymatgen()
        # sga1 = SpacegroupAnalyzer(pmg_s1).get_space_group_symbol()

        C2s = C1.subgroup(eps=1e-5)
        for C2 in C2s:
            pmg_s2 = C2.to_pymatgen()
            # sga2 = SpacegroupAnalyzer(pmg_s2).get_space_group_symbol()
            # prevent some numerical error
            if not sm.StructureMatcher().fit(pmg_s1, pmg_s2):
                C12 = pyxtal()
                C12.from_seed(pmg_s2)
                pmg_12 = C12.to_pymatgen()
                assert sm.StructureMatcher().fit(pmg_s1, pmg_12)

            # self.assertTrue(sm.StructureMatcher().fit(pmg_s1, pmg_s2))

        C1.subgroup(perms={"C": "Si"}, H=216)

    def test_from_seed(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = pmg_Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)
        struct = Structure(lattice, ["Si", "C"], coords)
        s1 = pyxtal()
        s1.from_seed(struct)
        s2 = s1.subgroup_once(eps=0)
        pmg_s1 = s1.to_pymatgen()
        pmg_s2 = s2.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg_s1, pmg_s2)

    def test_molecules(self):
        for name in [
            "aspirin",
            "resorcinol",
            "coumarin",
            "HAHCOI",
            "xxvi",
            "WEXBOS",
            "MERQIM",
            "LAGNAL",
            "YICMOP",
            "LUFHAW",
            "JAPWIH",
            "AXOSOW01",
            "PAHYON01",
        ]:
            cif = cif_path + name + ".cif"
            struc = pyxtal(molecular=True)
            struc.from_seed(seed=cif, molecules=[name])
            pmg_struc = struc.to_pymatgen()
            pmg_s1 = Structure.from_file(cif, primitive=True)
            assert sm.StructureMatcher().fit(pmg_struc, pmg_s1)

            Cs = struc.subgroup(eps=0, max_cell=1)
            for C in Cs:
                pmg_s2 = C.to_pymatgen()
                assert sm.StructureMatcher().fit(pmg_struc, pmg_s2)

    def test_hydrate(self):
        # glycine dihydrate
        cif = cif_path + "gdh.cif"
        struc = pyxtal(molecular=True)
        struc.from_seed(seed=cif, molecules=["Glycine-z", "H2O"])
        pmg_struc = struc.to_pymatgen()
        pmg_s1 = Structure.from_file(cif, primitive=True)
        assert sm.StructureMatcher().fit(pmg_struc, pmg_s1)

    def test_special(self):
        cif = cif_path + "191.vasp"
        struc = pyxtal()
        struc.from_seed(seed=cif)
        for _i in range(100):
            struc.subgroup_once(0.2, None, None, "t+k", 2)


class TestLoad(unittest.TestCase):
    def test_atomic(self):
        s1 = pyxtal()
        s1.from_random(3, 36, ["C", "Si"], [4, 8])
        s2 = pyxtal()
        s2.load_dict(s1.save_dict())
        pmg_s1 = s1.to_pymatgen()
        pmg_s2 = s2.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg_s1, pmg_s2)

    def test_molecular(self):
        s1 = pyxtal(molecular=True)
        s1.from_random(3, 36, ["H2O"], [4])
        s2 = pyxtal()
        s2.load_dict(s1.save_dict())
        pmg_s1 = s1.to_pymatgen()
        pmg_s2 = s2.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg_s1, pmg_s2)


class TestPartial(unittest.TestCase):
    def test_Al2SiO5(self):
        cell = Lattice.from_para(7.8758, 7.9794, 5.6139, 90, 90, 90)
        spg = 58
        elements = ["Al", "Si", "O"]
        composition = [8, 4, 20]
        sites = [
            {
                "4e": [0.0000, 0.0000, 0.2418],
                "4g": [0.1294, 0.6392, 0.0000],
            },
            {"4g": [0.2458, 0.2522, 0.0000]},
            [],  # empty for oxygen
        ]

        s = pyxtal()
        s.from_random(3, spg, elements, composition, lattice=cell, sites=sites)
        assert s.valid

        sites2 = [
            {
                "4e": [0.0000, 0.0000, 0.2418],
                "4g": [0.1294, 0.6392, 0.0000],
            },
            {"4g": [0.2458, 0.2522, 0.0000]},
            {"4g": [0.4241, 0.3636, 0.0000]},  # partial info on O
        ]

        s = pyxtal()
        s.from_random(3, spg, elements, composition, lattice=cell, sites=sites2)
        assert s.valid


class Test_resort(unittest.TestCase):
    def test_molecule(self):
        rng = np.random.default_rng(0)
        # glycine dihydrate
        cif = cif_path + "gdh.cif"
        struc = pyxtal(molecular=True)
        struc.from_seed(seed=cif, molecules=["Glycine-z", "H2O"])
        N1 = len(struc.mol_sites)
        l = list(range(len(struc.mol_sites)))
        rng.shuffle(l)
        struc.mol_sites = [struc.mol_sites[i] for i in l]
        struc.resort()
        N2 = len(struc.mol_sites)
        assert N1 == N2

    def test_atom(self):
        cif = cif_path + "aspirin.cif"
        struc = pyxtal()
        rng = np.random.default_rng(0)
        struc.from_seed(seed=cif)
        N1 = len(struc.atom_sites)
        l = list(range(len(struc.atom_sites)))
        rng.shuffle(l)
        struc.atom_sites = [struc.atom_sites[i] for i in l]
        struc.resort()
        N2 = len(struc.atom_sites)
        assert N1 == N2

#class Test_rng(unittest.TestCase):
#    """
#    Test rng generators in two ways
#    1. random_state as a fixed integer
#    2. random_state as a generator
#    """
#    def test_rng_integer(self):
#        xtal = pyxtal(); xtal.from_random(3, 194, ['C'], [8], random_state=0)
#        xs = xtal.get_1d_rep_x()
#        assert np.sum((xs - np.array([4.679, 6.418, 0.943])**2)) < 1e-2
#
#        xtal = pyxtal(molecular=True)
#        xtal.from_random(3, 19, ['aspirin'], [4], random_state=0)
#        rep = xtal.get_1D_representation().x
#        d1 = np.array([115, 17.294, 15.077, 9.018])
#        d2 = np.array([0.677, 0.243, 0.612, 0.057, -1.194, 0.110])
#        assert np.sum((rep[0] - d1)**2) < 1e-2
#        assert np.sum((rep[1][1:-1] - d2)**2) < 1e-2
#
#
#    def test_rng_generator(self):
#        rng = np.random.default_rng(1)
#        xtal = pyxtal()
#        xtal.from_random(3, 194, ['C'], [8], random_state=rng)
#        xs = xtal.get_1d_rep_x()
#        assert np.sum((xs - np.array([7.442, 2.110])**2)) < 1e-2
#
#        xtal.from_random(3, 194, ['C'], [8], random_state=rng)
#        xs = xtal.get_1d_rep_x()
#        assert np.sum((xs - np.array([5.864, 5.809])**2)) < 1e-2
#
#        xtal = pyxtal(molecular=True)
#        xtal.from_random(3, 19, ['aspirin'], [4], random_state=rng)
#        rep = xtal.get_1D_representation().x
#        d1 = np.array([115, 14.207, 18.334, 9.028])
#        d2 = np.array([0.294, 0.627, 0.528, 157.052, -11.968, -171.851])
#        assert np.sum((rep[0] - d1)**2) < 1e-2
#        assert np.sum((rep[1][1:-1] - d2)**2) < 1e-2
#
#        xtal.from_random(3, 19, ['aspirin'], [4], random_state=rng)
#        rep = xtal.get_1D_representation().x
#        d1 = np.array([115, 12.763, 16.639, 11.073])
#        d2 = np.array([0.504, 0.127, 0.585, -21.523, -68.406, 152.839])
#        assert np.sum((rep[0] - d1)**2) < 1e-2
#        assert np.sum((rep[1][1:-1] - d2)**2) < 1e-2

class Test_operations(unittest.TestCase):
    def test_inverse(self):
        coord0 = [0.35, 0.1, 0.4]
        coords = np.array(
            [
                [0.350, 0.100, 0.400],
                [0.350, 0.100, 0.000],
                [0.350, 0.100, 0.000],
                [0.350, 0.000, 0.667],
                [0.350, 0.000, 0.250],
                [0.350, 0.350, 0.400],
                [0.350, 0.350, 0.500],
                [0.350, 0.350, 0.000],
                [0.350, 0.350, 0.350],
                [0.100, 0.100, 0.100],
                [0.400, 0.400, 0.400],
                [0.350, 0.000, 0.000],
                [0.000, 0.100, 0.400],
                [0.350, 0.000, 0.400],
            ]
        )
        xyzs = [
            "x,y,z",
            "x,y,0",
            "y,x,0",
            "x,0,2/3",
            "0,x,1/4",
            "x,x,z",
            "x,-x,1/2",
            "2x,x,0",
            "-2x,-0.5x,-x+1/4",
            "-2y,-0.5y,-y+1/4",
            "-2z,-0.5z,-z+1/4",
            "0,0,x",
            "-y/2+1/2,-z,0",
            "-z,-x/2+1/2,0",
        ]

        for i, xyz in enumerate(xyzs):
            op = SymmOp.from_xyz_str(xyz)
            inv_op = get_inverse(op)
            coord1 = op.operate(coord0)
            coord2 = inv_op.operate(coord1)
            assert np.allclose(coord2, coords[i], rtol=0.01)
            # strs = "{:6.3f} {:6.3f} {:6.3f}".format(*coord0)
            # strs += "  {:12s}  ".format(op.as_xyz_str())
            # strs += "{:6.3f} {:6.3f} {:6.3f}".format(*coord1)
            # strs += "  {:12s}  ".format(inv_op.as_xyz_str())
            # strs += "{:6.3f} {:6.3f} {:6.3f}".format(*coord2)
            # print(strs)

    def test_swap_wp(self):
        g = Group(38)
        wp = g[4]
        wp1, trans = wp.swap_axis([1, 0, 2])

        g = Group(71)
        wp = g[5]
        wp1, trans = wp.swap_axis([0, 2, 1])
        wp1, trans = wp.swap_axis([1, 2, 0])
        wp1, trans = wp.swap_axis([2, 1, 0])

    def test_alternative(self):
        for name in ["BTO-Amm2", "lt_quartz", "GeF2", "lt_cristobalite", "PVO"]:
            s = pyxtal()
            s.from_seed(cif_path + name + ".cif")
            pmg_s1 = s.to_pymatgen()
            strucs = s.get_alternatives()
            for struc in strucs:
                pmg_s2 = struc.to_pymatgen()
                assert sm.StructureMatcher().fit(pmg_s1, pmg_s2)

    def test_wyc_sets(self):
        for i in range(1, 229):
            Group(i, quick=True).get_alternatives()["No."]

    def test_trans(self):
        mats = [
            np.array([[1, 0, 1], [0, 1, 0], [0, 0, 1]]),
            np.array([[1, 0, -1], [0, 1, 0], [0, 0, 1]]),
        ]

        for n in range(3, 16):
            hns = Hall(n).hall_numbers
            for hn in hns:
                g = Group(hn, use_hall=True)
                for i in range(1, len(g)):
                    wp = g[i].copy()
                    for mat in mats:
                        wp.transform_from_matrix(mat)
                        c1, c2 = wp.update()
                        # print(hn, i, (c1 or c2))
                        assert c1 or c2

    # def test_image(self):
    #    from pyxtal.descriptor import spherical_image
    #    c1 = pyxtal(molecular=True)
    #    for name in ['benzene', 'aspirin', 'naphthalene']:
    #        c1.from_seed(seed=cif_path+name+".cif", molecules=[name])
    #        for model in ['contact', 'molecule']:
    #            sph = spherical_image(c1, model=model)
    #            sph.align()
    #            print(name, model)

    class TestSubstitution(unittest.TestCase):
        def test_substitute_1_2(self):
            data = [
                (227, ["8a"], [3.6], ["C"], 1),
                (
                    92,
                    ["4a", "8b"],
                    [5.0847, 7.0986, 0.2944, 0.0941, 0.2410, 0.8256],
                    ["Si", "O"],
                    1,
                ),
            ]
            for d in data:
                (spg, wps, rep, elements, N) = d
                s = pyxtal()
                xtals = s.from_spg_wps_rep(spg, wps, rep, elements)
                assert len(xtals) == N

        def test_criteria(self):
            criteria = {"CN": {"B": 4, "N": 4}, "cutoff": 1.9, "exclude_ii": True}
            xtals = pyxtal().substitute_1_2({"C": ["B", "N"]}, ratio=[1, 1], criteria=criteria)
            assert xtals[0].check_validity(criteria)


if __name__ == "__main__":
    unittest.main()
