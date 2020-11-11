from pyxtal.crystal import random_crystal
from pyxtal.symmetry import Group
import pymatgen.analysis.structure_matcher as sm
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
from pyxtal.lattice import cellsize

for G in range(1, 231):
    g = Group(G)
    subs = g.get_max_t_subgroup()
    indices = subs['index']
    hs = subs['subgroup']
    relations = subs['relations']
    tran = subs['transformation']
    letter = str(g[0].multiplicity) + g[0].letter
    print(G)
    C1 = random_crystal(G, ['C'], [int(g[0].multiplicity/cellsize(g))], sites=[[letter]])
    pmg_s1 = C1.to_pymatgen()
    sga1 = SpacegroupAnalyzer(pmg_s1).get_space_group_symbol()
    # each subgroup
    for i in range(len(relations)):
        H = hs[i]
        if H > 1:
            C2 = C1.subgroup(H=H, eps=0, idx=i)
            pmg_s2 = C2.to_pymatgen()
            try:
                sga2 = SpacegroupAnalyzer(pmg_s2, symprec=1e-4).get_space_group_symbol()
            except:
                #print("unable to find the space group")
                sga2 = None
            print(G, H, g.symbol, sga1, Group(H).symbol, sga2, i)
            if not sm.StructureMatcher().fit(pmg_s1, pmg_s2):
                print('WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW')
                print(C1)
                print(C2)
                #print(pmg_s1)
                #print(pmg_s2)
                #print(tran[i])
                #import sys; sys.exit()
    #    # each site
    #    wps.reverse()
    #    for j, wp in enumerate(wps):
    #        N_G = int(g[j].multiplicity * np.linalg.det(tran[i][:3,:3]))
    #        N_H = sum([int(w[:-1]) for w in wp])

    #        if abs(N_G - N_H)>0:
    #            strs = "problem in G:{:d} [{:d}] -> H:{:d} [{:d}]".format(G, N_G, H, N_H)
    #            print(strs, wp)
