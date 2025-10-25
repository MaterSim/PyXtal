# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

from pyxtal.lattice import Lattice
from pyxtal.symmetry import Group, Wyckoff_position

def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)

cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestSymmetry(unittest.TestCase):
    def test_from_symops_wo_grou(self):
        data = [
        (["x, y, z", "-x, y+1/2, -z"], 4, 6),
        (["x, y, z", "-x+1/2, -y, z+1/2", "-x, y, z", "x+1/2, -y, z+1/2"], 31, 155),
        (["x, y, z", "-x, -y, -z", "-x+1/2, y+1/2, -z", "x+1/2, -y+1/2, z"], 14, 83),
        (["x, y, z", "-x, -y, -z", "-x+1/2, y+1/2, -z+1/2", "x+1/2, -y+1/2, z+1/2"], 14, 82),
        ]
        for d in data:
            (strs, spg, hall) = d
            wp = Wyckoff_position.from_symops_wo_group(strs)
            assert wp.number == spg
            assert wp.hall_number == hall

    def test_from_symops(self):
        data = [
        (["x, y, z", "-x, y+1/2, -z"], 4, 6),
        (["x, y, z", "-x+1/2, -y, z+1/2", "-x, y, z", "x+1/2, -y, z+1/2"], 31, 155),
        ]
        for d in data:
            (strs, spg, hall) = d
            G = Group(spg)
            wp = Wyckoff_position.from_symops(strs, G)
            assert wp.number == spg
            assert wp.hall_number == hall


    def test_reflections(self):
        # 141: hkl: h+k+l=2n
        # hk0: h,k=2n

        data = [(230, 1, 1, 1, False),
                (229, 1, 1, 1, False),
                (228, 2, 1, 1, False),
                (228, 2, 2, 0, True),
                (227, 1, 1, 1, True),
                (227, 2, 2, 0, True),
                (226, 2, 2, 2, True),
                (226, 2, 2, 0, True),
                (225, 1, 1, 1, True),
                (225, 2, 2, 1, False),
                (142, 1, 1, 2, True),
                (141, 1, 1, 2, True),
                (141, 2, 2, 1, False),
                (62, 1, 1, 4, True),
                (62, 2, 2, 0, True),
                (36, 1, 1, 8, True),
                (36, 2, 2, 0, True),
                (18, 3, 0, 0, False),
                (17, 0, 0, 1, False),
                (14, 3, 0, 0, True),
                (14, 0, 0, 3, False),
                (12, 1, 2, 0, False),
                (2, 2, 2, 0, True),
                ]
        for d in data:
            (spg, h, k, l, expected) = d
            G = Group(spg)
            print(spg, h, k, l, G.is_valid_hkl(h, k, l), expected)
            assert G.is_valid_hkl(h, k, l) == expected

if __name__ == "__main__":
    unittest.main()
