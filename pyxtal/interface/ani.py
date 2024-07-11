import os
from time import time
import numpy as np
import torchani

from ase.optimize.fire import FIRE
from ase.constraints import ExpCellFilter
from ase.spacegroup.symmetrize import FixSymmetry


def ANI_relax(struc, opt_cell=False, step=500, fmax=0.1, logfile=None, max_time=6.0):
    """
    ani optimizer
    Args:
    struc: ase atoms object
    step: optimization steps (int)
    max_time: float (minutes)
    """
    calc = torchani.models.ANI2x().ase()
    struc.set_calculator(calc)
    struc.set_constraint(FixSymmetry(struc)) 
    if opt_cell:
        ecf = ExpCellFilter(struc)
        if logfile is not None:
            dyn = FIRE(ecf, a=0.1, logfile=logfile)
        else:
            dyn = FIRE(ecf, a=0.1)
    else:
        if logfile is not None:
            dyn = FIRE(struc, a=0.1, logfile=logfile)
        else:
            dyn = FIRE(struc, a=0.1)
    
    # Run relaxation
    if step < 50:
        dyn.run(fmax=fmax, steps=step)
    else:
        t0 = time()
        dyn.run(fmax=fmax, steps=int(step/2))
        # If time is too long, only run half steps
        if (time()-t0)/60 < max_time/2 :
            dyn.run(fmax=fmax, steps=int(step/2))
    return struc


class ANI():
    """
    This is a calculator to perform oragnic crystal structure optimization in ANI
    We assume that the geometry has been well optimized by classical FF

    Args:
        struc: pyxtal object
        opt: 'conv', 'conp', 'single'
    """

    def __init__(self, struc, opt_lat=True, logfile=None):
        self.structure = struc
        self.opt_lat = opt_lat
        self.stress = None
        self.forces = None
        self.optimized = True
        self.positions = None
        self.cell = None
        self.group = self.structure.group
        self.cputime = 0
        self.calculator = torchani.models.ANI2x().ase()
        self.logfile = logfile

    def run(self):

        t0 = time()
        s = self.structure.to_ase(resort=False)
        s.set_constraint(FixSymmetry(s))
        s.set_calculator(self.calculator)
        dyn = FIRE(s, a=0.1, logfile=self.logfile)
        dyn.run(fmax=0.1, steps=10)
        #print(s)
        if self.opt_lat:
            ecf = ExpCellFilter(s)
            dyn = FIRE(ecf, a=0.1, logfile=self.logfile)
            dyn.run(fmax=0.1, steps=10)
            self.structure.lattice.set_matrix(s.get_cell())
        positions = s.get_scaled_positions()
        try:
            #s.write('../1.cif', format='cif')
            count = 0
            for i, site in enumerate(self.structure.mol_sites):
                coords0, _ = site._get_coords_and_species(first=True)
                coords1 = positions[count:count+len(site.molecule.mol)]
                for j, coor in enumerate(coords1):
                    diff = coor - coords0[j]
                    diff -= np.round(diff)
                    abs_diff = np.dot(diff, s.get_cell())
                    #print(j, coor, coords0[j], diff, np.linalg.norm(abs_diff))
                    if abs(np.linalg.norm(abs_diff)) < 2.0:
                        coords1[j] = coords0[j] + diff
                    else:
                        print(coords1[j], coords1[j], np.linalg.norm(abs_diff))
                        import sys; sys.exit()

                site.update(coords1, self.structure.lattice)
                count += len(site.molecule.mol)*site.wp.multiplicity
            self.structure.optimize_lattice()
            self.structure.energy = s.get_potential_energy()
            self.cell = s.get_cell()
            #print(self.structure.lattice)
        except:
            self.structure.energy = 10000
            self.optimized = False
            print("Structure is wrong after optimization")

        self.cputime = time()-t0

if __name__ == "__main__":

    from pyxtal.db import database
    import warnings
    warnings.filterwarnings("ignore")

    work_dir = "tmp"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    db = database('benchmarks/test.db')
    struc = db.get_pyxtal('ACSALA')

    calc = ANI(struc)
    print(calc.structure.lattice)
    calc.run()
    print(calc.structure.energy)
    print(calc.structure.lattice)