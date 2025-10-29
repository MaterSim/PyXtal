import torch

# --- Temporary fix for PyTorch 2.6 weights_only change ---
from torch.serialization import add_safe_globals
add_safe_globals([slice])   # allow 'slice' to be unpickled
# ----------------------------------------------------------

from ase.io import read, write
from fairchem.core import pretrained_mlip, FAIRChemCalculator
import signal
from time import time
import numpy as np
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter
from ase.optimize.fire import FIRE
import logging
from pyxtal.optimize import WFS, DFS, QRS
from pyxtal import pyxtal
from pyxtal.util import get_pmg_dist
import os
from mace.calculators import mace_mp
_cached_mace_mp = None
from mace.calculators import mace_off


def get_calculator(calculator):
    global _cached_mace_mp

    if type(calculator) is str:
        if calculator == 'ANI':
            import torchani
            calc = torchani.models.ANI2x().ase()

        elif calculator == 'MACE':
            if _cached_mace_mp is None:
                _cached_mace_mp = mace_mp(
                    model='small',
                    dispersion=True,
                    device='cpu'
                )
            calc = _cached_mace_mp

        elif calculator == 'MACEOFF':
            calc = mace_off(model='medium', device='cpu')
            
        elif calculator == 'FAIRChem':
            # Initialize FAIRChem
            predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
            calc = FAIRChemCalculator(predictor, task_name="omc")
            
        else:
            raise ValueError(f"Unknown calculator: {calculator}")
            
    else:
        calc = calculator

    return calc

class ASE_optimizer:
    """
    This is a ASE optimizer to perform oragnic crystal structure optimization.
    We assume that the geometry has been well optimized by classical FF.

    Args:
        struc: pyxtal object
        calculator (str): 'ANI', 'MACE'
        opt_lat (bool): to opt lattice or not
        log_file (str): output file
    """

    def __init__(self, struc, calculator="FAIRChem", opt_lat=True, logfile="ase_log_ebdc_local"):
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

    def run(self, fmax_target=0.01):
        t0 = time()
        s = self.structure.to_ase(resort=False)
        s.set_constraint(FixSymmetry(s))
        s.set_calculator(self.calculator)
    
        obj = UnitCellFilter(s) if self.opt_lat else s
        dyn = FIRE(obj, a=0.01, logfile=self.logfile)
    
        # <-- Key line: no 'steps' argument
        dyn.run(fmax=fmax_target)
    
        if self.opt_lat:
            self.structure.lattice.set_matrix(s.get_cell())
    
        positions = s.get_scaled_positions()
        count = 0
        for _i, site in enumerate(self.structure.mol_sites):
            coords0, _ = site._get_coords_and_species(first=True)
            coords1 = positions[count : count + len(site.molecule.mol)]
            for j, coor in enumerate(coords1):
                diff = coor - coords0[j]
                diff -= np.round(diff)
                abs_diff = np.dot(diff, s.get_cell())
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
    
        s.set_calculator()
        s.set_constraint()
        self.cputime = time() - t0
        self.optimized = bool(getattr(dyn, "converged", False))


#("/Users/mmukta/Downloads/HOF-EBDC_aka_ZJU-HOF-60.cif", "O=C(O)c2cc(C#Cc1cc(C(=O)O)cc(C(=O)O)c1)cc(C(=O)O)c2")
#("/Users/mmukta/Downloads/HOF-BDDC_aka_ZJU-HOF-62.cif", "O=C(O)c2cc(C#CC#Cc1cc(C(=O)O)cc(C(=O)O)c1)cc(C(=O)O)c2")
#("/Users/mmukta/Downloads/HOF-1a.cif", "Nc8nc(N)nc(c7ccc(C(c2ccc(c1nc(N)nc(N)n1)cc2)(c4ccc(c3nc(N)nc(N)n3)cc4)c6ccc(c5nc(N)nc(N)n5)cc6)cc7)n8")

if __name__ == "__main__":
    import os, warnings
    from pyxtal.db import database
    warnings.filterwarnings("ignore")

    work_dir = "tmp"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    #db = database("pyxtal/database/test.db")
    #struc = db.get_pyxtal("ACSALA")
    data = [
    ("/Users/mmukta/Desktop/Cocrystal/ebdc-local-maceOff1500.cif", "O=C(O)c2cc(C#Cc1cc(C(=O)O)cc(C(=O)O)c1)cc(C(=O)O)c2")
    ]

    for d in data:
        cif, smiles = d
        c = pyxtal(molecular=True)
        c.from_seed(cif, molecules=[smiles+'.smi'])
        pmg0 = c.to_pymatgen()
        if c.has_special_site():
            c1 = c.to_subgroup(); print(c1)
            pmg = c1.to_pymatgen()
            if get_pmg_dist(pmg0, pmg) > 0.1:
                print("The reference structure is not a valid subgroup.")
                m = c1.mol_sites[0]
                m.rotate(ax_id=2, angle=180)
                pmg = c1.to_pymatgen()
                print("Distance after flip", get_pmg_dist(pmg0, pmg))
            c = c1
            pmg = c.to_pymatgen()
        else:
            pmg = pmg0
    calc = ASE_optimizer(c)
    print(calc.structure.lattice)
    #calc.run(steps=1500)
    calc.run(fmax_target=0.1)
    print(calc.structure.energy)
    print(calc.structure.lattice)
    '''
    calc.structure.to_file("maceOff/ebdc-local.cif")
    from pymatgen.core import Structure
    # Load CIF
    structure = Structure.from_file("maceOff/ebdc-local.cif")
    # Get density in g/cm³
    print("Density:", structure.density, "g/cm³")
    total_molecules = 0
    print("Molecular sites and multiplicities:")
    for i, site in enumerate(c.mol_sites):
        print(f"Site {i+1}: Wyckoff multiplicity = {site.wp.multiplicity}")
        total_molecules += site.wp.multiplicity

    print(f"\nTotal molecules per unit cell: {total_molecules}")
    '''