from structure import random_crystal
from spglib import get_symmetry_dataset
from ase.io import read, write
from ase.calculators.vasp import Vasp
from random import randint
from vasprun import vasprun
import os

"""
This is a script to generate random crystals 
and then perform multiple steps of optimization with ASE-VASP 
"""
#Test paramaters: [sg_min, sg_max], [species], [numIons]
pstress = 50.0
maxvec = 25.0
minvec = 2.5
maxangle = 150
minangle = 30

setup = {'B': '_s'} #use soft pp to speed up calc 

def set_vasp(level=0, pstress=0.0000, setup=None):
    default0 = {'xc': 'pbe',
            'npar': 4,
            'kgamma': True,
            'lcharg': False,
            'lwave': False,
            'ibrion': 2,
            'nsw': 50,
            'potim': 0.02,
            'pstress': pstress,
            }
    if level==0:
        default1 = {'prec': 'low',
                'kspacing': 0.4,
                'isif': 2,
	        'ediff': 1e-2,
                }
    elif level==1:
        default1 = {'prec': 'normal',
                'kspacing': 0.3,
                'isif': 3,
	        'ediff': 1e-3,
                }
    elif level==2:
        default1 = {'prec': 'high',
                'kspacing': 0.2,
                'isif': 3,
	        'ediff': 1e-3,
                }
    elif level==3:
        default1 = {'prec': 'high',
                'encut': 600,
                'kspacing': 0.15,
                'isif': 3,
	        'ediff': 1e-4,
                }
    dict_vasp = dict(default0, **default1)
    return Vasp(**dict_vasp)

def good_lattice(struc):
    para = struc.get_cell_lengths_and_angles()
    if (max(para[:3])<maxvec) and (max(para[3:])<maxangle) and (min(para[3:])>minangle):
        return True
    else:
        return False

def optimize(struc, dir1):
    os.mkdir(dir1)
    os.chdir(dir1)
    
    # Step1: ISIF = 2
    struc.set_calculator(set_vasp(level=0))
    print(struc.get_potential_energy())

    # Step2: ISIF = 3
    struc = read('CONTCAR',format='vasp')
    struc.set_calculator(set_vasp(level=1))
    print(struc.get_potential_energy())

    # Step3: ISIF = 3 with high precision 
    struc = read('CONTCAR',format='vasp')
    if good_lattice(struc):
        struc.set_calculator(set_vasp(level=2, pstress=pstress))
        print(struc.get_potential_energy())
        struc = read('CONTCAR',format='vasp')

        if good_lattice(struc):
            struc.set_calculator(set_vasp(level=2, pstress=pstress))
            print(struc.get_potential_energy())
            struc = read('CONTCAR',format='vasp')
           
            if good_lattice(struc):
                struc.set_calculator(set_vasp(level=3, pstress=pstress))
                struc.get_potential_energy()
                result = vasprun().values
                print('%12.6f %6.2f' % (result['calculation']['energy_per_atom'], result['gap']))

species = ['B']
factor = 1.0
dir0 = os.getcwd()

for i in range(2000):
    os.chdir(dir0)
    numIons = [12]
    run = True
    while run:
        #numIons[0] = randint(8,16)
        sg = randint(2,230)
        #print(sg, species, numIons, factor)
        rand_crystal = random_crystal(sg, species, numIons, factor)
        if rand_crystal.valid:
            run = False
    try:
        struc = rand_crystal.struct
        label = str(struc.formula).replace(" ","")
        dir1 = str(i)+'-'+label
        struc.to(fmt="poscar", filename = '1.vasp')
        ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)['number']
        print('Space group  requested: ', sg, 'generated', ans)
        struc = read('1.vasp',format='vasp')
        optimize(struc, dir1)
    except:
        print('something is wrong')
