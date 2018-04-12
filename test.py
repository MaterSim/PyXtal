from structure import random_crystal
from spglib import get_symmetry_dataset
from sys import exit
import profile
#from pymatgen.io.cif import CifWriter

#Test paramaters: [sg_min, sg_max], [species], [numIons]
params_list = [
[[120,121],[3],[8]]]
#[[2,230],[20,3],[4,8]],
#[[2,230],[20,3,8],[4,8,10]],
#[[2,15],[3],[8]],
#[[143, 167],[3],[8]],
#[[2,15],[3],[4]],
#[[143,167],[3],[4]]
#]

factor = 1.0
for params in params_list:
    sg_min, sg_max = params[0][0], params[0][1]
    species, numIons = params[1], params[2]
    for sg in range(sg_min, sg_max+1):
        profile.run("random_crystal(sg, species, numIons, factor)")
        #rand_crystal = random_crystal(sg, species, numIons, factor)

        #if rand_crystal.valid:
        #    rand_crystal.struct.to(fmt="poscar", filename = '1.vasp')
        #    ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)['number']
        #    print('Space group  requested: ', sg, 'generated', ans)
            #if ans < int(sg/1.2):
            #    print('something is wrong')
            #    exit(0)
            #    break
            #print(CifWriter(rand_crystal.struc, symprec=0.1).__str__())
            #print('Space group:', finder.get_space_group_symbol(), 'tolerance:', tol)
            #output wyckoff sites only

        #else: 
        #    print('something is wrong')
        #    #print(len(new_struct.frac_coords))
        #    #print(new_struct)
        #    #exit(0)
        #    break

