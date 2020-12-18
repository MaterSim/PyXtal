from pyxtal import pyxtal
from pyxtal.interface.util import good_lattice
from ase.calculators.vasp import Vasp
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

def single_optimize(struc, level, pstress, setup, dir0=None):
    """
    single optmization

    Args: 
        struc: pyxtal structure
        level: vasp calc level
        pstress: external pressure
        setup: vasp setup 
        dir0: calculation directory

    Returns:
        the structure, energy and time costs
    """
    cwd = os.getcwd()
    if dir0 is not None:
        if not os.path.exists(dir0):
            os.makedirs(dir0)
        os.chdir(dir0)

    ase_atoms = struc.to_ase()
    ase_atoms.set_calculator(set_vasp(level, pstress, setup))
    energy = ase_atoms.get_potential_energy()/len(ase_atoms)
    time, ncore = read_OUTCAR()
    struc = pyxtal()
    struc.from_seed('CONTCAR')
    struc.optimize_lattice()

    os.chdir(cwd)
    return struc, energy, time

def single_point(struc, setup=None, dir0=None):
    """
    single optmization

    Args: 
        struc: pyxtal structure
        level: vasp calc level
        pstress: external pressure
        setup: vasp setup 
        dir0: calculation directory

    Returns:
        the energy and forces
    """
    cwd = os.getcwd()
    if dir0 is not None:
        if not os.path.exists(dir0):
            os.makedirs(dir0)
        os.chdir(dir0)

    ase_atoms = struc.to_ase()
    ase_atoms.set_calculator(set_vasp(level=4, setup=setup))
    energy = ase_atoms.get_potential_energy()
    forces = ase_atoms.get_forces()
    os.chdir(cwd)
    return energy, forces


def optimize(struc, dir0, levels=[0,2,3], pstress=0, setup=None):
    """
    multi optimization

    Args:
        struc: pyxtal structure
        dir0: calculation directory
        levels: list of vasp calc levels
        pstress: external pressure
        setup: vasp setup 

    Returns:
        list of structures, energies and time costs
    """

    times = []
    strucs = []
    engs = []
    for level in levels:
        struc, e, t = single_optimize(struc, level, pstress, setup, dir0)
        times.append(t)
        strucs.append(struc)
        engs.append(e)
        # skip the structures with bad lattices
        if not good_lattice(struc):
            break
    return strucs, engs, times


