from pyxtal.interface.gulp import GULP
from pyxtal.crystal import random_crystal
from ase import Atoms
import os

file = "C-POSCARs"
if os.path.exists(file):
    os.remove(file)
for i in range(10):
    struc = random_crystal(19, ["C"], [16], 1.0)
    if struc.valid:
        calc = GULP(struc, ff="tersoff.lib")
        calc.run()
        s = Atoms(struc.sites, scaled_positions=calc.positions, cell=calc.cell)
        s.write("1.vasp", format='vasp', vasp5=True, direct=True)
        os.system("cat 1.vasp >> " + file)
        print(i, calc.energy)
        #print(calc.stress)
        #print(calc.cell)
