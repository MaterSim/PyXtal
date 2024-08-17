# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.symmetry import Wyckoff_position
from pyxtal.XRD import Similarity


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestPXRD(unittest.TestCase):
    def test_similarity(self):
        C1 = pyxtal()
        C1.from_random(3, 227, ["C"], [8], sites=[["8a"]])
        xrd1 = C1.get_XRD()
        C2 = C1.subgroup_once(eps=1e-3)
        xrd2 = C1.get_XRD()
        p1 = xrd1.get_profile()
        p2 = xrd2.get_profile()
        s = Similarity(p1, p2, x_range=[15, 90])
        assert 0.9 < s.value < 1.001

        C2.apply_perturbation(1e-3, 1e-3)
        xrd3 = C2.get_XRD()
        xrd3.get_profile()
        s = Similarity(p1, p2, x_range=[15, 90])
        assert 0.95 < s.value < 1.001

if __name__ == "__main__":
    unittest.main()
