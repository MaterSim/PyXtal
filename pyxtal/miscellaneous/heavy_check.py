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
from pymatgen.core import Structure
from pyxtal import pyxtal
from pyxtal.supergroup import supergroups

for i, name in enumerate(glob("pyxtal/miscellaneous/cifs/*.cif")):

    # 1, read from cif
    s = pyxtal()
    s.from_seed(name)
    pmg_s1 = s.to_pymatgen()
    pmg0 = Structure.from_file(name)
    G = s.group.number
    
    print(i, name, len(s.atom_sites))
    if len(s.atom_sites) <= 6:
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
                if H>2 and H != G and H in s.group.get_max_subgroup_numbers():
                    struc_h = s.subgroup_once(eps=0.05, H=H, group_type=gtype, mut_lat=False)
                    try: 
                        sup = supergroups(struc_h, G=G, d_tol=0.3, max_per_G=500)
                        if sup.strucs is not None:
                            match = False
                            for struc in sup.strucs:
                                pmg_g = struc.to_pymatgen()
                                if sm.StructureMatcher().fit(pmg_g, pmg_s1):
                                    match = True
                                    break
                            if not match:
                                print("Cannot recover the original structure", G, '<-', H)
                        else:
                            print("Error in supergroup", G, '<-', H)
                    except RuntimeError:
                        print("no splitter skip", name)
