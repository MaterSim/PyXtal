# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

from pyxtal import pyxtal
from pyxtal.XRD import Similarity, check_pxrd_match

def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)

cif_path = resource_filename("pyxtal", "database/cifs/")

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

    def test_match(self):
        dia = pyxtal()
        dia.from_prototype('diamond')
        xrd = dia.get_XRD()
        pxrd1 = xrd.get_profile(res=0.25, user_kwargs={"FWHM": 0.25})
        assert check_pxrd_match(dia, pxrd1)
        sub = dia.subgroup_once(H=141, eps=0.025)
        sub = sub.subgroup_once(H=74, eps=0.01)
        sub = sub.subgroup_once(H=12, eps=0.01)
        assert check_pxrd_match(sub, pxrd1, s_tol=0.7)


if __name__ == "__main__":
    unittest.main()
