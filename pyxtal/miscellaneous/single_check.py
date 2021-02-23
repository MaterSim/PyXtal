"""
A script to systematicly check
1, read_cif
2, alternative setting
3, subgroup
4, supergroup
"""

from glob import glob
import numpy as np
import pymatgen.analysis.structure_matcher as sm
from pymatgen import Structure
from pyxtal import pyxtal
from pyxtal.supergroup import supergroups

name, H, gtype = "R-3", 148, 'k' #needs 1->3
name, H, gtype = "Pm3", 47, 't' #needs 1->3
name, H, gtype = "Fd3m", 166, 't'
name, H, gtype = "R32", 5, 't'
#name, H, gtype = "R-3c", 15, 't'
#name, H, gtype = "Pmma", 25, 't'
name = "pyxtal/miscellaneous/cifs/" + name + ".cif"
print(name)

s = pyxtal()
s.from_seed(name)
#print(s)
pmg_s1 = s.to_pymatgen()
G = s.group.number
for i in range(100):
    struc_h = s.subgroup_once(eps=0.15, H=H, group_type=gtype, mut_lat=False, max_cell=3)
    print(struc_h)
    print(s)
    sup = supergroups(struc_h, G=G, d_tol=0.2, show=True, max_per_G=500)
    if sup.strucs is not None:
        match = False
        for struc in sup.strucs:
            pmg_g = struc.to_pymatgen()
            if sm.StructureMatcher().fit(pmg_g, pmg_s1):
                match = True
                break
        if not match:
            print("==Wrong: cannot recover the original structure")
    else:
        print(sup)
        print("==Error: cannot generate structure")
        break
