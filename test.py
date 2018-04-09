from structure import *

factor = 2.0
#Test paramaters: [sg_min, sg_max], [species], [numIons]
params_list = [
[[2,230],[3],[8]],
[[2,230],[20,3],[4,8]],
[[2,230],[20,3,8],[4,8,10]],
[[2,15],[3],[8]],
[[143, 167],[3],[8]],
[[2,15],[3],[4]],
[[143,167],[3],[4]]
]
for params in params_list:
    sg_min, sg_max = params[0][0], params[0][1]
    species, numIons = params[1], params[2]
    for sg in range(sg_min, sg_max+1):
        rand_crystal = random_crystal(sg, species, numIons, factor)
        new_struct, good_struc = rand_crystal.struct, rand_crystal.good_struct
        if good_struc:
            new_struct.to(fmt="poscar", filename = '1.vasp')
            cell = read_vasp('1.vasp')
            #ans = spglib.get_spacegroup(cell)
            ans = spglib.get_symmetry_dataset(cell, symprec=1e-1)['number']
            print('Space group  requested: ', sg, 'generated', ans)
            if ans < int(sg/1.2):
                print('something is wrong')
                break
            #print(CifWriter(new_struct, symprec=0.1).__str__())
            #print('Space group:', finder.get_space_group_symbol(), 'tolerance:', tol)
            #output wyckoff sites only

        else: 
            print('something is wrong')
            #print(len(new_struct.frac_coords))
            break
            #print(new_struct)
