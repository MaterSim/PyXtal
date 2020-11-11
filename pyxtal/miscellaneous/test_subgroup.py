from pyxtal.crystal import random_crystal
import pymatgen.analysis.structure_matcher as sm
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


#G, fac = 197, 2 
#sites = ['6b', '8c']
#for H in [146, 23]:
#    for i in range(10):
#        numIons = int(sum([int(i[:-1]) for i in sites])/fac)
#        C1 = random_crystal(G, ['C'], [numIons], sites=[sites])
#        C2 = C1.subgroup(H=H)
#        #print(C1)
#        #print(C2)
#        pmg_s1 = C1.to_pymatgen()
#        pmg_s2 = C2.to_pymatgen()
#        sga1 = SpacegroupAnalyzer(pmg_s1).get_space_group_symbol()
#        sga2 = SpacegroupAnalyzer(pmg_s2).get_space_group_symbol()
#        print(i, sga1, sga2, sm.StructureMatcher().fit(pmg_s1, pmg_s2))
#        if not sm.StructureMatcher().fit(pmg_s1, pmg_s2):
#            import sys
#            sys.exit()
G, fac = 227, 4
sites = ['8a', '32e']
#for H in [166, 141, 203, 210, 216]:
for H in [203, 210, 216]:
    for i in range(10):
        numIons = int(sum([int(i[:-1]) for i in sites])/fac)
        C1 = random_crystal(G, ['C'], [numIons], sites=[sites])
        C2 = C1.subgroup(H=H)
        #print(C1)
        #print(C2)
        pmg_s1 = C1.to_pymatgen()
        pmg_s2 = C2.to_pymatgen()
        sga1 = SpacegroupAnalyzer(pmg_s1).get_space_group_symbol()
        sga2 = SpacegroupAnalyzer(pmg_s2).get_space_group_symbol()
        print(i, sga1, sga2, sm.StructureMatcher().fit(pmg_s1, pmg_s2))
        if not sm.StructureMatcher().fit(pmg_s1, pmg_s2):
            import sys
            sys.exit()


