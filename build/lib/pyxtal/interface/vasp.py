from ase.io import read
from ase import Atoms
from ase.calculators.vasp import Vasp
from pyxtal.interface.util import symmetrize_cell, good_lattice, pymatgen2ase
import os, time

"""
A script to perform multistages vasp calculation
"""

def set_vasp(level=0, pstress=0.0000, setup=None):
    default0 = {'xc': 'pbe',
            'npar': 8,
            'kgamma': True,
            'lcharg': False,
            'lwave': False,
            'ibrion': 2,
            'pstress': pstress,
            'setups': setup,
            }
    if level==0:
        default1 = {'prec': 'low',
                'algo': 'normal',
                'kspacing': 0.4,
                'isif': 4,
                'ediff': 1e-2,
                'nsw': 10,
                'potim': 0.02,
                }
    elif level==1:
        default1 = {'prec': 'normal',
                'algo': 'normal',
                'kspacing': 0.3,
                'isif': 3,
                'ediff': 1e-3,
                'nsw': 25,
                'potim': 0.05,
                }
    elif level==2:
        default1 = {'prec': 'accurate',
                'kspacing': 0.2,
                'isif': 3,
                'ediff': 1e-3,
                'nsw': 50,
                'potim': 0.1,
                }
    elif level==3:
        default1 = {'prec': 'accurate',
                'encut': 600,
                'kspacing': 0.15,
                'isif': 3,
                'ediff': 1e-4,
                'nsw': 50,
                }
    elif level==4:
        default1 = {'prec': 'accurate',
                'encut': 600,
                'kspacing': 0.15,
                'isif': 3,
                'ediff': 1e-4,
                'nsw': 0,
                }

    dict_vasp = dict(default0, **default1)
    return Vasp(**dict_vasp)

def read_OUTCAR(path='OUTCAR'):
    """read time and ncores info from OUTCAR"""
    time = 0
    ncore = 0
    for line in open(path, 'r'):
        if line.rfind('running on  ') > -1:
            ncore = int(line.split()[2])
        elif line.rfind('Elapsed time ') > -1:
            time = float(line.split(':')[-1])

    return time, ncore

def single_optimize(struc, level, pstress, mode, setup):
    """single optmization"""
    struc = symmetrize_cell(struc, mode)
    struc.set_calculator(set_vasp(level, pstress, setup))
    energy = struc.get_potential_energy()
    print(energy)
    time, ncore = read_OUTCAR()
    struc = read('CONTCAR',format='vasp')
    return struc, energy, time

def optimize(struc, dir0, modes=['C','C','C','P','P'], pstress=0, setup=None):
    """multi optimization"""
    struc = pymatgen2ase(struc)
    os.mkdir(dir0)
    os.chdir(dir0)
    time0 = []
    struc0 = []
    energy0 = []
    for i, mode in enumerate(modes):
        struc, energy, time = single_optimize(struc, i, pstress, mode, setup)
        time0.append(time)
        struc0.append(struc)
        energy0.append(energy)
        if not good_lattice(struc):
            break
    return struc0, energy0, time0


