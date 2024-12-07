# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import numpy as np

from pyxtal.lattice import Lattice
from pyxtal.symmetry import Group, Wyckoff_position, get_wyckoffs
from pyxtal.wyckoff_site import atom_site


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestWP(unittest.TestCase):
    def test_wp_check_translation(self):
        pass

    def test_wp_site_symm(self):
        data = [
            (143, 1, "3.."),
            (150, 1, ".2."),
            (152, 1, ".2."),
            (154, 1, ".2."),
            (155, 1, ".2"),
            (160, 1, ".m"),
            (160, 2, "3m"),
            (164, 4, ".2/m."),
            (165, 1, ".2."),
            (177, 3, ".2."),
            (178, 2, ".2."),
            (180, 3, ".2."),
            (181, 3, ".2."),
            (230, 6, ".32"),
        ]
        for d in data:
            (sg, i, symbol) = d
            wp = Group(sg)[i]
            wp.get_site_symmetry()
            if wp.site_symm != symbol: print("\n========", wp.site_symm, d, "==========\n")
            assert wp.site_symm == symbol

    def test_wp_dof(self):
        for sg in range(1, 231):
            g = Group(sg)
            for wp in g:
                axs = wp.get_frozen_axis()
                assert wp.get_dof() + len(axs) == 3

    def test_wp_label(self):
        symbol = wp1.get_label()
        assert symbol == "8b"
        symbol = wp2.get_label()
        assert symbol == "4a"

    def test_merge(self):
        pt, wp, _ = wp1.merge([0.05, 0.7, 0.24], l01.matrix, 0.5)
        symbol = wp.get_label()
        assert symbol == "4a"
        pt, wp, _ = wp1.merge([0.15, 0.7, 0.24], l01.matrix, 0.5)
        symbol = wp.get_label()
        assert symbol == "8b"

        wp = Group(167)[0]
        cell = np.diag([9, 9, 7])
        for pt in [[0.12, 0, 0.25], [0, 0.1316, 0.25]]:
            _, wpt, _ = wp.merge(pt, cell, 0.1)
            symbol = wpt.get_label()
            assert symbol == "18e"

        for pt in [[0, 0, 3 / 4], [2 / 3, 1 / 3, 7 / 12]]:
            _, wpt, _ = wp.merge(pt, cell, 0.1)
            symbol = wpt.get_label()
            assert symbol == "6a"

    def test_search_generator(self):
        wp = Group(167)[1]
        for pt in [[0, 0.13, 0.25], [0.13, 0, 0.25]]:
            wp0 = wp.search_generator(pt)
            assert wp0 is not None

    def test_get_wyckoff(self):
        for i in [1, 2, 229, 230]:
            get_wyckoffs(i)
            get_wyckoffs(i, organized=True)

    def test_setting(self):
        wp = Group(14)[0]
        wp.transform_from_matrix()
        assert not wp.is_standard_setting()

    def test_standarize(self):
        pass

    def test_is_equivalent(self):
        g = Group(15)
        wp = g[0]
        a = [0.10052793, 0.12726851, 0.27405404]
        b = [-0.10052642, -0.12726848, -0.27405526]
        c = [0.60052642, 0.62726848, 0.27405526]
        d = [-0.60052642, -0.62726848, -0.27405526]
        e = [0, 2.54537267e-01, 0]
        assert wp.are_equivalent_pts(a, b)
        assert wp.are_equivalent_pts(b, c)
        assert wp.are_equivalent_pts(d, a)
        assert not wp.are_equivalent_pts(a, e)

        wp = g[1]
        a = [0.00, 0.127, 0.254]
        b = [-0.01, -0.127, -0.250]
        assert wp.are_equivalent_pts(a, b)

    def test_project(self):
        pt = np.array([0.0629, 0.1258, 0.25])
        g = Group(178)
        wp = g[1]
        xyz0 = wp.project(pt, np.eye(3))
        diff = np.sum((pt - xyz0) ** 2)
        assert diff < 1e-08

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

                    diff = p1 - p2
                    diff -= np.rint(diff)
                    if np.linalg.norm(diff) > 0.02:
                        # res = '{:2d} {:28s}'.format(i, op0.as_xyz_str())
                        # res += ' {:28s}'.format(op1.as_xyz_str())
                        # res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p1)
                        # res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p2)
                        # res += '{:6.3f} {:6.3f} {:6.3f}'.format(*diff)
                        # print(res)
                        return False
            return True

        pt = [0.1333, 0.1496, 0.969]

        cell = Lattice.from_para(9.395, 7.395, 8.350, 91, 101, 92, ltype="triclinic")
        assert check_error(range(1, 3), pt, cell)

        cell = Lattice.from_para(9.395, 7.395, 8.350, 90, 101, 90, ltype="monoclinic")
        assert check_error(range(3, 16), pt, cell)

        cell = Lattice.from_para(9.395, 7.395, 8.350, 90, 90, 90, ltype="orthorhombic")
        assert check_error(range(16, 74), pt, cell)

        cell = Lattice.from_para(9.395, 9.395, 8.350, 90, 90, 90, ltype="tetragonal")
        assert check_error(range(74, 143), pt, cell)

        cell = Lattice.from_para(9.395, 9.395, 8.350, 90, 90, 120, ltype="hexagonal")
        assert check_error(range(143, 195), pt, cell)

        cell = Lattice.from_para(9.395, 9.395, 9.395, 90, 90, 90, ltype="cubic")
        assert check_error(range(195, 231), pt, cell)


class Test_wyckoff_site(unittest.TestCase):
    def test_atom_site(self):
        """
        Test the search function
        """
        wp = Group(227)[-5]
        arr = np.array([0.1, 0.1, 0.1])
        for xyz in [
            [0.6501, 0.15001, 0.5999],
            [0.6501, 0.14999, 0.5999],
            [0.6499, 0.14999, 0.5999],
            [0.6499, 0.14999, 0.60001],
        ]:
            site = atom_site(wp, xyz, search=True)
            assert np.allclose(site.position, arr, rtol=0.001)

if __name__ == "__main__":
    unittest.main()
