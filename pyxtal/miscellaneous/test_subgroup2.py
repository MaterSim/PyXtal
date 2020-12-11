from pyxtal import pyxtal
from pyxtal.symmetry import Group
import pymatgen.analysis.structure_matcher as sm
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np

for G in range(143, 195):
    g = Group(G)
    letter = str(g[0].multiplicity) + g[0].letter
    C1 = pyxtal()
    C1.from_random(3, G, ['C'], [g[0].multiplicity], sites=[[letter]])
    #print(C1)
    pmg_s1 = C1.to_pymatgen()
    sga1 = SpacegroupAnalyzer(pmg_s1).get_space_group_symbol()

    # each subgroup
    #C2s = C1.subgroup(eps=0, group_type='t')
    try: 
        C2s = C1.subgroup(eps=0, group_type='k', max_cell=4)
        for C2 in C2s:
            #print(C2)
            pmg_s2 = C2.to_pymatgen()
            try:
                sga2 = SpacegroupAnalyzer(pmg_s2, symprec=1e-4).get_space_group_symbol()
            except:
                #print("unable to find the space group")
                sga2 = None
            print(G, C2.group.number, g.symbol, C2.group.symbol, sga1, sga2)
            if not sm.StructureMatcher().fit(pmg_s1, pmg_s2):
                print('WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW')
                print(C1)
                print(C2)
    except RuntimeError:
        pass
            #print(pmg_s1)
            #print(pmg_s2)
            #print(tran[i])
            #import sys; sys.exit()
