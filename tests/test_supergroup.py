# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import pymatgen.analysis.structure_matcher as sm

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.supergroup import supergroup, supergroups
from pyxtal.symmetry import Wyckoff_position


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")


class TestSupergroup(unittest.TestCase):
    def test_supergroup(self):
        """
        call supergroup from pyxtal
        """
        data = [
            ("NbO2", 141, 5),
        ]
        for d in data:
            (cif, g, N) = d
            s = pyxtal()
            s.from_seed(cif_path + cif + ".cif")
            strucs, sols = s.supergroup(g, 0.5)
            assert len(strucs) > 0

    def test_make_pyxtal(self):
        data = {
            # "NbO2": 141,
            "GeF2": 62,
            "lt_quartz": 180,
            "NiS-Cm": 160,  # 9b->2a+4b
        }
        for cif in data:
            s = pyxtal()
            s.from_seed(cif_path + cif + ".cif")
            my = supergroup(s, G=data[cif])
            sols = my.search_supergroup(max_solutions=1)
            for sol in sols:
                struc_high = my.make_pyxtal_in_supergroup(sol)
                strucs = my.make_pyxtals_in_subgroup(sol, 3)
                pmg1 = struc_high.to_pymatgen()
                pmg2 = strucs[-1].to_pymatgen()
                rms = sm.StructureMatcher().get_rms_dist(pmg1, pmg2)[0]
                print(cif, rms)
                assert rms < 0.001

    def test_by_group(self):
        data = {
            "BTO-Amm2": 221,
            "BTO": 221,
            #"lt_cristobalite": 227,
            "NaSb3F10": 194,
        }
        for cif in data:
            s = pyxtal()
            s.from_seed(cif_path + cif + ".cif")
            sup = supergroups(s, G=data[cif], show=False)
            assert sup.strucs is not None

    def test_by_path(self):
        data = {
            "MPWO": [59, 71, 139, 225],
        }
        for cif in data:
            s = pyxtal()
            s.from_seed(cif_path + cif + ".cif")
            sup = supergroups(s, path=data[cif], show=False)
            assert sup.strucs is not None

    def test_long(self):
        paras = (
            ["I4_132", 98],  # 1-3
            ["P3_112", 5],  # 1-3
            ["P6_422", 21],  # 1-3
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

            if sup.strucs is None:
                print("Problem in ", name)
            assert sup.strucs is not None

            match = False
            for struc in sup.strucs:
                pmg_g = struc.to_pymatgen()
                if sm.StructureMatcher().fit(pmg_g, pmg_s1):
                    match = True
                    break
            if not match:
                print("Problem in ", name)
            assert match

    def test_long2(self):
        paras = (
            ["P4_332", 155],  # 1-4 splitting
            #["Fd3", 70],  # 1-3
            ["Pm3", 47],  # 1-3
            ["Fd3m", 166],
            ["R-3c", 15],
            #["R32", 5],
            ["R-3", 147],  # 1-3, k
            ["P4_332", 96],  # 1-3
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

            if len(sols) == 0:
                print("problem", name)
            assert len(sols) > 0

            struc_high = my.make_pyxtal_in_supergroup(sols[0])
            struc_sub = my.make_pyxtals_in_subgroup(sols[0], 3)[-1]

            pmg_high = struc_high.to_pymatgen()
            pmg_sub = struc_sub.to_pymatgen()

            dist1 = sm.StructureMatcher().get_rms_dist(pmg_high, pmg_sub)[0]
            dist2 = sm.StructureMatcher().get_rms_dist(pmg_s1, pmg_high)[0]
            assert dist1 < 0.001
            assert dist2 < 0.001

    def test_multi(self):
        data = {
            "BTO": [123, 221],
            "lt_cristobalite": [98, 210, 227],
            #"BTO-Amm2": [65, 123, 221],
            # "NaSb3F10": [186, 194],
            # "NaSb3F10": [176, 194],
            # "MPWO": [59, 71, 139, 225],
        }
        for cif in data:
            s = pyxtal()
            s.from_seed(cif_path + cif + ".cif")
            sup = supergroups(s, path=data[cif], show=False, max_per_G=2500)
            strucs = sup.get_transformation()
            pmg_0, pmg_1 = s.to_pymatgen(), sup.strucs[-1].to_pymatgen()
            pmg_2, pmg_3 = strucs[0].to_pymatgen(), strucs[1].to_pymatgen()
            dist1 = sm.StructureMatcher().get_rms_dist(pmg_0, pmg_2)[0]
            dist2 = sm.StructureMatcher().get_rms_dist(pmg_1, pmg_3)[0]
            print(cif, dist1, dist2)
            if dist2 > 1e-3:
                print(pmg_1)
                print(pmg_3)
            assert dist1 < 0.001
            assert dist2 < 0.001

    def test_similarity(self):
        paras = [
            ("0-G62", "2-G71"),
            ("0-G62", "3-G139"),
            # ('0-G62', '4-G225'),
        ]
        for para in paras:
            (cif1, cif2) = para
            s1 = pyxtal()
            s2 = pyxtal()
            s1.from_seed(cif_path + cif1 + ".cif")
            s2.from_seed(cif_path + cif2 + ".cif")
            pmg_s2 = s2.to_pymatgen()

            strucs, _, _, _, _ = s2.get_transition(s1)

            if strucs is None:
                print("Problem between ", cif1, cif2)
            else:
                struc_G_in_H = strucs[-1]
                s3 = pyxtal()
                s3.from_seed(struc_G_in_H.to_pymatgen(), tol=1e-3)
                pmg_s3 = s3.to_pymatgen()
                dist = sm.StructureMatcher().get_rms_dist(pmg_s2, pmg_s3)[0]
                assert dist < 0.001
                assert s3.group.number == s2.group.number

    def test_get_disp_sets(self):
        s1 = pyxtal()
        s1.from_seed(cif_path + "dist_6_0.cif")
        s2 = pyxtal()
        s2.from_seed(cif_path + "dist_6_1.cif")
        _, _, _, d = s1.get_disps_sets(s2, 1.0)
        assert d < 0.15

        #s1 = pyxtal()
        #s1.from_seed(cif_path + "sim-0.vasp")
        #s2 = pyxtal()
        #s2.from_seed(cif_path + "sim-1.vasp")
        #_, _, _, d = s1.get_disps_sets(s2, 1.0, 0.3)
        #assert d < 0.02

if __name__ == "__main__":
    unittest.main()
