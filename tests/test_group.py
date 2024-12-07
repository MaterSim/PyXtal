# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.symmetry import Group, Wyckoff_position
from pyxtal.util import generate_wp_lib


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestGroup(unittest.TestCase):
    def test_generate_wp_lib(self):
        wps = generate_wp_lib([227, 228], composition=[1, 2])
        assert len(wps) == 18
        wps = generate_wp_lib([227, 228], composition=[1, 1, 3])
        assert len(wps) == 9

    def test_list_wyckoff_combinations(self):
        g = Group(64)
        a1, _, _ = g.list_wyckoff_combinations([4, 2])
        assert len(a1) == 0
        a2, _, _ = g.list_wyckoff_combinations([4, 8], quick=False)
        assert len(a2) == 8

    def test_print_group_and_dof(self):
        for d in [(1, 6), (15, 4), (60, 3), (143, 2), (208, 1)]:
            (sg, dof_ref) = d
            g = Group(sg)
            dof = g.get_lattice_dof()
            assert dof == dof_ref

    def test_get_spg_symmetry_object(self):
        spg_list = [14, 36, 62, 99, 143, 160, 225, 230] #182, 191, 
        ans = [32, 18, 36, 21, 16, 19, 62, 62] # 24, 48
        for spg, num in zip(spg_list, ans):
            g = Group(spg)
            ss = g.get_spg_symmetry_object()
            matrix = ss.to_matrix_representation_spg()
            assert num == sum(sum(matrix))

    def test_short_path(self):
        g = Group(217)
        path = g.short_path_to_general_wp(7)
        assert path[-1][2] == 145

    def test_spg_symmetry(self):
        N_polar, N_centro, N_chiral = 0, 0, 0
        for sg in range(1, 231):
            g = Group(sg, quick=True)
            _pg, polar, centro, chiral = g.point_group, g.polar, g.inversion, g.chiral
            if polar:
                N_polar += 1
            if centro:
                N_centro += 1
            if chiral:
                N_chiral += 1
        assert N_polar == 68
        assert N_centro == 92
        assert N_chiral == 65

    def test_ferroelectric(self):
        pairs = [(4, 1), (187, 4), (222, 5)]
        for pair in pairs:
            (sg, N) = pair
            assert len(Group(sg, quick=True).get_ferroelectric_groups()) == N

    def test_check_compatible(self):
        assert Group(225).check_compatible([64, 28, 24]) == (True, True)
        assert Group(227).check_compatible([8]) == (True, False)
        assert Group(227).check_compatible([4]) == (False, False)
        assert Group(19).check_compatible([6]) == (False, False)

    def test_search_supergroup_paths(self):
        paths = Group(59, quick=True).search_supergroup_paths(139, 2)
        assert paths == [[71, 139], [129, 139], [137, 139]]

    def test_get_splitters(self):
        s = pyxtal()
        s.from_seed(cif_path + "3-G139.cif")
        g = Group(225)
        solutions = g.get_splitters_from_structure(s, "t")
        assert len(solutions) == 3

    def test_add_k_transitions(self):
        paras = [
            (189, 26, 6),
            (11, 6, 2),
            (62, 33, 3),
            (63, 33, 5),
        ]

        for p in paras:
            (g, h, n_path) = p
            gr = Group(g, quick=True)
            path = gr.search_subgroup_paths(h)[0]
            solutions = gr.add_k_transitions(path)
            assert len(solutions) == n_path

    def test_to_subgroup(self):
        s = pyxtal(molecular=True)
        s.from_seed(cif_path + "benzene.cif", ["benzene"])
        c = s.to_subgroup()
        assert c.valid

    def test_get_wyckoff_position_from_xyz(self):
        g = Group(5)
        pos = [
            ([0.0, 0.0, 0.0], "2a"),
            ([0.5, 0.5, 0.5], None),
            ([0.1, 0.0, 0.1], "4c"),
            ([0.5, 0.0, 0.0], None),
            ([0.0, 0.1, 0.0], "2a"),
            ([0.0, 0.5, 0.5], "2b"),
            ([0.1, 0.2, 0.3], "4c"),
        ]
        for d in pos:
            (p0, wp0) = d
            wp = g.get_wyckoff_position_from_xyz(p0)
            if wp is None:
                assert wp0 is None
            else:
                assert wp.get_label() == wp0

    def test_short_path_to_general_wp(self):
        data = [('t', 0, 141, ['16h']), ('t', 6, 122, ['16e'])]
        G = Group(227)
        assert G.short_path_to_general_wp(4) == data

if __name__ == "__main__":
    unittest.main()
