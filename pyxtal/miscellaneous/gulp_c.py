from pyxtal.interface.gulp import GULP
from pyxtal.crystal import random_crystal
from ase import Atoms
import os
from spglib import get_symmetry_dataset

file = "C-POSCARs"
if os.path.exists(file):
    os.remove(file)
for i in range(10):
    #struc = random_crystal(194, ["C"], [4], 1.0)
    struc = random_crystal(227, ["C"], [2], 1.0)
    if struc.valid:
        calc = GULP(struc, ff="tersoff.lib")
        calc.run()
        s = Atoms(calc.sites, scaled_positions=calc.positions, cell=calc.cell)
        info = get_symmetry_dataset(s, symprec=1e-1)
        s.write("1.vasp", format='vasp', vasp5=True, direct=True)
        os.system("cat 1.vasp >> " + file)
        print("{:4d} {:8.3f} {:s}".format(i, calc.energy/8, info['international']))
        #print(calc.stress)
        #print(calc.cell)
