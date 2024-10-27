# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import pymatgen.analysis.structure_matcher as sm
from pymatgen.core import Structure

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.symmetry import Hall, Wyckoff_position


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

    def test_preassigned_sites(self):
        sites = [["1b"], ["1b"], ["2c", "1b"]]
        struc = pyxtal()
        struc.from_random(3, 99, ["Ba", "Ti", "O"], [1, 1, 3], 1.0, sites=sites)
        assert struc.valid

        struc = pyxtal()
        struc.from_random(3, 225, ["C"], [12], 1.0, sites=[["4a", "8c"]])
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
