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

for i, name in enumerate(glob("pyxtal/miscellaneous/cifs/*.cif")[:160]):
    print(i, name)

    # 1, read from cif
    s = pyxtal()
    s.from_seed(name)
    pmg_s1 = s.to_pymatgen()
    pmg0 = Structure.from_file(name)
    G = s.group.number
    
    if not sm.StructureMatcher().fit(pmg_s1, pmg0):
        print("Error in reading cif")

    # 2, alternative setting
    strucs = s.get_alternatives()
    for i, struc in enumerate(strucs):
        pmg_s2 = struc.to_pymatgen()
        if not sm.StructureMatcher().fit(pmg_s1, pmg_s2):
            print("Error in alternative setting")
            print(s)
            print(struc)
            break
        
    # 3, subgroup
    for gtype in ['t', 'k']:
        valid = True
        try:
            struc_h = s.subgroup_once(eps=0, group_type=gtype, max_cell=3)
            H = struc_h.group.number
            pmg_h = struc_h.to_pymatgen()
            if not sm.StructureMatcher().fit(pmg_s1, pmg_h):
                print("Error in subgroup", gtype)
        except RuntimeError:
            print("no splitter skip", name)
            valid = False
        # 4, supergroup
        if valid:
            #print(G, H)
            if H != G and H in s.group.get_max_subgroup_numbers():
                struc_h = s.subgroup_once(eps=0.05, H=H, group_type=gtype)
                try: 
                    error = False
                    sup = supergroups(struc_h, G=G, d_tol=0.3)
                    if sup.strucs is not None:
                        pmg_g = sup.strucs[-1].to_pymatgen()
                        if not sm.StructureMatcher().fit(pmg_g, pmg_s1):
                            error = True
                    else:
                        error = True

                    if error:
                        print("Error in supergroup", G, '<-', H)
                except RuntimeError:
                    print("no splitter skip", name)
