import signal
from time import time
import numpy as np
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter
from ase.optimize.fire import FIRE
import logging

import os
from mace.calculators import mace_mp
_cached_mace_mp = None


def get_calculator(calculator):
    global _cached_mace_mp

    if type(calculator) is str:
        if calculator == 'ANI':
            import torchani
            calc = torchani.models.ANI2x().ase()
        else:
            if _cached_mace_mp is None:
                _cached_mace_mp = mace_mp(model='small', dispersion=True, device='cpu')
            calc = _cached_mace_mp
    else:
        calc = calculator

    return calc

def ASE_relax(struc, calculator, opt_cell=False, step=500, fmax=0.1, logfile=None, max_time=10.0, label='ase'):
#def ASE_relax(struc, calculator, opt_cell=False, step=500, fmax=0.1, logfile='ase.log', max_time=10.0, label='ase'):
    """
    ASE optimizer
    Args:
        struc: ase atoms object
        calculator (str): 'ANI', 'MACE'
        step: optimization steps (int)
        max_time: float (minutes)
    """

    def handler(signum, frame):
        raise TimeoutError("Optimization timed out")

    step_init = min([30, int(step/2)])
    logger = logging.getLogger()
    max_time *= 60
    timeout = int(max_time)
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout)
    #logger.info(f"{label} start calculation")
    _fmax = 1e+5

    try:
    #if True:
        calc = get_calculator(calculator)
        struc.set_calculator(calc)
        struc.set_constraint(FixSymmetry(struc))
        if opt_cell:
            ecf = UnitCellFilter(struc)
            dyn = FIRE(ecf, a=0.1, logfile=logfile) if logfile is not None else FIRE(ecf, a=0.1)
        else:
            dyn = FIRE(struc, a=0.1, logfile=logfile) if logfile is not None else FIRE(struc, a=0.1)

        # Run relaxation
        dyn.run(fmax=fmax, steps=step_init)
        forces = dyn.optimizable.get_forces()
        _fmax = np.sqrt((forces ** 2).sum(axis=1).max())

        if _fmax < 1e+3 and step > step_init:
            dyn.run(fmax=fmax, steps=step-step_init)
            forces = dyn.optimizable.get_forces()
            _fmax = np.sqrt((forces ** 2).sum(axis=1).max())
            eng = struc.get_potential_energy() / len(struc)
            if _fmax > 100:
                logger.info(f"Warning {label} big stress {eng:.2f} / {_fmax:.2f}, skip")
                struc = None
            else:
                logger.info(f"{label} Success  {eng:.2f} / {_fmax:.2f}")
        else:
            logger.info(f"Warning {label} big stress {_fmax:.2f} for 20 steps, skip")
            struc = None
        signal.alarm(0)  # Cancel the alarm if finished within time

    except TimeoutError:
        logger.info(f"Warning {label} timed out after {timeout} seconds.")
        struc = None

    except TypeError:
        logger.info(f"Warning {label} spglib error in getting the lattice")
        struc = None
        signal.alarm(0)  # Cancel the alarm if finished within time

    tag = 'False' if struc is None else 'True'
    logger.info(f"Finishing {label} {tag}")
    #signal.alarm(0)  # Cancel the alarm
    return struc #, eng, _fmax

class ASE_optimizer:
    """
    This is a ASE optimizer to perform oragnic crystal structure optimization.
    We assume that the geometry has been well optimized by classical FF

    Args:
        struc: pyxtal object
        calculator (str): 'ANI', 'MACE'
        opt_lat (bool): to opt lattice or not
        log_file (str): output file
    """

    def __init__(self, struc, calculator='MACE', opt_lat=True, logfile=None):
        self.structure = struc
        self.calculator = get_calculator(calculator)
        self.opt_lat = opt_lat
        self.stress = None
        self.forces = None
        self.optimized = True
        self.positions = None
        self.cell = None
        self.cputime = 0
        self.logfile = logfile

    def run(self, steps=10):
        t0 = time()
        s = self.structure.to_ase(resort=False)
        s.set_constraint(FixSymmetry(s))
        s.set_calculator(self.calculator)#; print("Setup Fire")

        if not self.opt_lat:
            dyn = FIRE(s, a=0.1, logfile=self.logfile)#, force_consistent=False)
            dyn.run(fmax=0.1, steps=steps)
        else:
            #ecf = FrechetCellFilter(s)
            ecf = UnitCellFilter(s)
            dyn = FIRE(ecf, a=0.1, logfile=self.logfile)#, force_consistent=False)
            dyn.run(fmax=0.1, steps=steps)
            self.structure.lattice.set_matrix(s.get_cell())

        positions = s.get_scaled_positions()
        #try:
        if True:
            # s.write('../1.cif', format='cif')
            count = 0
            for _i, site in enumerate(self.structure.mol_sites):
                coords0, _ = site._get_coords_and_species(first=True)
                coords1 = positions[count : count + len(site.molecule.mol)]
                for j, coor in enumerate(coords1):
                    diff = coor - coords0[j]
                    diff -= np.round(diff)
                    abs_diff = np.dot(diff, s.get_cell())
                    # print(j, coor, coords0[j], diff, np.linalg.norm(abs_diff))
                    if abs(np.linalg.norm(abs_diff)) < 2.0:
                        coords1[j] = coords0[j] + diff
                    else:
                        print(coords1[j], coords1[j], np.linalg.norm(abs_diff))
                        import sys; sys.exit()

                site.update(coords1, self.structure.lattice)
                count += len(site.molecule.mol) * site.wp.multiplicity
            self.structure.optimize_lattice()
            self.structure.energy = s.get_potential_energy()
            self.cell = s.get_cell()
            # print(self.structure.lattice)
        #except:
        #    self.structure.energy = 10000
        #    self.optimized = False
        #    print("Structure is wrong after optimization")
        s.set_calculator()
        s.set_constraint()
        self.cputime = time() - t0


if __name__ == "__main__":
    import os, warnings
    from pyxtal.db import database
    warnings.filterwarnings("ignore")

    work_dir = "tmp"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    db = database("pyxtal/database/test.db")
    struc = db.get_pyxtal("ACSALA")

    calc = ASE_optimizer(struc)
    print(calc.structure.lattice)
    calc.run()
    print(calc.structure.energy)
    print(calc.structure.lattice)
