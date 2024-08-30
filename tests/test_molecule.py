# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import pymatgen.analysis.structure_matcher as sm
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.molecule import pyxtal_molecule
from pyxtal.symmetry import Group, Wyckoff_position


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestMolecule(unittest.TestCase):
    def test_get_orientations_in_wp(self):
        m = pyxtal_molecule("Benzene")
        g = Group(61)
        assert len(m.get_orientations_in_wp(g[0])) == 1
        assert len(m.get_orientations_in_wp(g[1])) == 1
        assert len(m.get_orientations_in_wp(g[2])) == 1


class TestMolecular(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], sites=[["8b"]])
        struc.to_file("tmp-molecular.cif")
        os.remove("tmp-molecular.cif")
        assert struc.valid

        # test space group
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        # print(sga.get_space_group_symbol())
        assert sga.get_space_group_number() >= 36

        # test rotation
        struc.mol_sites[0].orientation.axis
        struc.mol_sites[0].rotate(ax_vector=[1, 0, 0], angle=90)
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        # pmg_struc.to("cif", "1.cif")
        assert sga.get_space_group_symbol() == "Cmc2_1"
        # print(pmg_struc.frac_coords[:3])

    def test_sites(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 19, ["H2O"])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        assert sga.get_space_group_symbol() == "P2_12_12_1"

        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=[["4a", "4a"]])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        assert sga.get_space_group_symbol() == "Cmc2_1"

    def test_sites_xyz(self):
        struc = pyxtal(molecular=True)
        sites = [{"4e": [0.77, 0.57, 0.53]}]
        lat = Lattice.from_para(11.43, 6.49, 11.19, 90, 83.31, 90, ltype="monoclinic")
        struc.from_random(3, 14, ["aspirin"], [4], lattice=lat, sites=sites)
        assert struc.valid

    def test_read(self):
        # test reading structure from external
        struc = pyxtal(molecular=True)
        struc.from_seed(seed=cif_path + "aspirin.cif", molecules=["aspirin"])
        assert struc.lattice.ltype == "monoclinic"
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        assert sga.get_space_group_symbol() == "P2_1/c"
        C = struc.subgroup_once(eps=0, H=4)
        pmg_s2 = C.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg_struc, pmg_s2)

    def test_big_molecule(self):
        # print("test_big_molecule")
        for mol in ["ROY", "aspirin"]:
            struc = pyxtal(molecular=True)
            struc.from_random(3, 19, [mol], factor=1.4)
            assert struc.valid
            pair = struc.check_short_distances()
            if len(pair) > 0:
                print("short distances were detected")
                print(mol)
                print(pair)
            assert len(pair) == 0

    def test_from_random_site(self):
        spg, wps, elements, numIons = 224, [["24j"]], ["C"], [24]
        for _i in range(10):
            c = pyxtal()
            c.from_random(3, spg, elements, numIons, sites=wps)
            pair = c.check_short_distances(0.5)
            assert len(pair) == 0

    def test_c60(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["C60"], [4], 1.0)
        assert struc.valid

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

        for _i in range(3):
            struc = pyxtal(molecular=True)
            struc.from_random(3, 10, [Li, ps4], [6, 2], 1.2, conventional=False)
            if struc.valid:
                assert len(struc.to_pymatgen()) == 16

    def test_distance(self):
        Ag_xyz = """1
        AgC2N2H
        Ag         4.30800        8.26300       -0.2200
        """

        C2N2H7_xyz = """12
        AgC2N2H
        H          5.95800        5.80600       -0.9530
        N          5.24100        6.16800       -1.1210
        N          2.23200        6.99000       -0.6820
        C          4.02900        5.47000       -1.5870
        C          2.78500        5.61100       -0.6610
        H          3.69300        5.63500       -2.5830
        H          4.17400        4.42800       -1.6990
        H          3.12700        5.50500        0.4260
        H          2.05000        5.01200       -0.9300
        H          1.96000        7.20500       -1.3860
        H          1.59400        6.99200       -0.0710
        """
        with open("Ag.xyz", "w") as f:
            f.write(Ag_xyz)
        with open("C2N2H7.xyz", "w") as f:
            f.write(C2N2H7_xyz)

        for _i in range(10):
            c = pyxtal(molecular=True)
            c.from_random(3, 9, ["Ag.xyz", "C2N2H7.xyz"], [12, 12])
            short_bonds = c.check_short_distances(r=1.1)
            assert len(short_bonds) == 0
        os.remove("Ag.xyz")
        os.remove("C2N2H7.xyz")

    def test_molecular_2d(self):
        # print("test_molecular_2d")
        struc = pyxtal(molecular=True)
        struc.from_random(2, 20, ["H2O"])
        assert struc.valid

    def test_molecular_1d(self):
        struc = pyxtal(molecular=True)
        struc.from_random(1, 20, ["H2O"])
        assert struc.valid

    def test_preassigned_sites(self):
        sites = [["4a", "4a"]]
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=sites)
        assert struc.valid

    def test_special_sites(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 61, ["Benzene"], [4])
        assert struc.valid

if __name__ == "__main__":
    unittest.main()
