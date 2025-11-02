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
                (227, 2, 2, 2, True),
                (228, 2, 2, 2, True),
                (226, 2, 2, 2, True),
                (226, 2, 2, 0, True),
                (225, 1, 1, 1, True),
                (225, 2, 2, 1, False),
                (142, 1, 1, 2, True),
                (141, 1, 1, 2, True),
                (141, 2, 2, 1, False),
                (127, 0, 0, 1, True),
                (128, 0, 0, 1, False),
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

    def test_site_symms(self):

        data = [
            (146, "3."),
            (147, "-1", "-1", "3..", "3..", "-3..", "-3.."),
            (148, "-1", "-1", "3.", "-3.", "-3."),
            (149, "..2", "..2", "3..", "3..", "3..", "3.2", "3.2", "3.2", "3.2", "3.2", "3.2"),
            (150, ".2.", ".2.", "3..", "3..", "32.", "32."),
            (151, "..2", "..2"),
            (152, ".2.", ".2."),
            (153, "..2", "..2"),
            (154, ".2.", ".2."),
            (155, ".2", ".2", "3.", "32", "32"),
            (156, ".m.", "3m.", "3m.", "3m."),
            (157, "..m", "3..", "3.m"),
            (158, "3..", "3..", "3.."),
            (159, "3..", "3.."),
            (160, ".m", "3m"),
            (161, "3."),
            (162, "..m", "..2", "..2", "3..", "..2/m", "..2/m", "3.m", "3.2", "3.2", "-3.m", "-3.m"),
            (163, "..2", "-1", "3..", "3..", "3.2", "3.2", "-3..", "3.2"),
            (164, ".m.", ".2.", ".2.", ".2/m.", ".2/m.", "3m.", "3m.", "-3m.", "-3m."),
            (165, ".2.", "-1", "3..", "3..", "-3..", "32."),
            (166, ".m", ".2", ".2", ".2/m", ".2/m", "3m", "-3m", "-3m"),
            (167, ".2", "-1", "3.", "-3.", "32"),
            (168, "2..", "3..", "6.."),
            (171, "2..", "2.."),
            (172, "2..", "2.."),
            (173, "3..", "3.."),
            (174, "m..", "m..", "3..", "3..", "3..", "-6..", "-6..", "-6..", "-6..", "-6.."),
            (175, "m..", "m..", "2..", "3..", "2/m..", "2/m..", "6..", "-6..", "-6..", "6/m..", "6/m.."),
            (176, "m..", "-1", "3..", "3..", "-6..", "-6..", "-3..", "-6.."),
            (177, "..2", "..2", ".2.", ".2.", "2..", "3..", "222", "222", "6..", "3.2", "3.2", "622", "622"),
            (178, "..2", ".2."),
            (179, "..2", ".2."),
            (180, "..2", "..2", ".2.", ".2.", "2..", "2..", "222", "222", "222", "222"),
            (181, "..2", "..2", ".2.", ".2.", "2..", "2..", "222", "222", "222", "222"),
            (182, "..2", ".2.", "3..", "3..", "3.2", "3.2", "3.2", "32."),
            (183, ".m.", "..m", "2mm", "3m.", "6mm"),
            (184, "2..", "3..", "6.."),
            (185, "..m", "3..", "3.m"),
            (186, ".m.", "3m.", "3m."),
            (187, ".m.", "m..", "m..", "mm2", "mm2", "3m.", "3m.", "3m.", "-6m2", "-6m2", "-6m2", "-6m2", "-6m2", "-6m2"),
            (188, "m..", "..2", "3..", "3..", "3..", "-6..", "3.2", "-6..", "3.2"),
            (189, "m..", "m..", "..m", "3..", "m2m", "m2m", "3.m", "-6..", "-6..", "-62m", "-62m"),
            (190, "m..", ".2.", "3..", "3..", "-6..", "-6..", "-6..", "32."),
            (191, "m..", "m..", ".m.", "..m", "mm2", "mm2", "m2m", "m2m", "2mm", "3m.", "mmm", "mmm", "6mm",
             "-6m2", "-6m2", "6/mmm", "6/mmm"),
            (192, "m..", "..2", ".2.", "2..", "3..", "2/m..", "222", "6..", "-6..", "3.2", "6/m..", "622"),
            (193, "..m", "m..", "..2", "3..", "m2m", "..2/m", "3.m", "3.2", "-6..", "-3.m", "-62m"),
            (194, ".m.", "m..", ".2.", "mm2", ".2/m.", "3m.", "3m.", "-6m2", "-6m2", "-6m2", "-3m."),
            ]
        for d in data:
            g = Group(d[0])
            for i in range(1, len(d)):
                wp = g[i]
                wp.get_site_symmetry()
                assert wp.site_symm == d[i]


if __name__ == "__main__":
    unittest.main()
