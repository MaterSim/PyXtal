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

name, H, gtype = "Fd3m", 166, 't'
name, H, gtype = "R-3", 148, 'k' #needs 1->3
#Name, H, gtype = "Acam", 12, 't'
#Name, H, gtype = "Pcmb", 61, 'k'
#Name, H, gtype = "P2_122_1", 61, 'k'
#Name, H, gtype = "R32", 5, 't'
#Name, H, gtype = "I2mb", 8, 't'
#Name, H, gtype = "Acam", 57, 'k'
#Name, H, gtype = "P6_3mc", 36, 't'
#Name, H, gtype = "Ia-3d", 220, 't'
#Name, H, gtype = "R-3c", 15, 't'
#name, H, gtype = "I2_12_12_1", 17, 'k'
#name, H, gtype = "Fd-3m", 166, 't'
#name, H, gtype = "Pmma", 25, 't'
#name, H, gtype = "Pman", 60, 'k'
#name, H, gtype = "R32", 152, 'k'
#name, H, gtype = "P-42m", 16, 't'
name = "pyxtal/miscellaneous/cifs/" + name + ".cif"
print(name)

s = pyxtal()
s.from_seed(name)
pmg_s1 = s.to_pymatgen()
G = s.group.number
for i in range(100):
    struc_h = s.subgroup_once(eps=0.15, H=H, group_type=gtype, mut_lat=False, max_cell=3)
    print(struc_h)
    print(s)
    sup = supergroups(struc_h, G=G, d_tol=0.2, show=True, max_per_G=500)
    if sup.strucs is not None:
        pmg_g = sup.strucs[-1].to_pymatgen()
        if not sm.StructureMatcher().fit(pmg_g, pmg_s1):
            print("==================Error")
            break
    else:
        print(sup)
        print("==================wrong")
        break
