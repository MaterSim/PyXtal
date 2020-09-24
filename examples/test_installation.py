import pyxtal
print("PyXtal: ", pyxtal.__version__)

import ase
print("ase: ", ase.__version__)

from pyxtal.crystal import random_crystal_2D
print("Using PyXtal to generate structure")
struc = random_crystal_2D(75, ["C"], [4], thickness=0)
print("Convert PyXtal structure to ASE")
ase_struc = struc.to_ase()

from pyxtal.interface.lammpslib import run_lammpslib
from pyxtal.interface.gulp import single_optimize as gulp_opt
from lammps import lammps

# Set up lammps
import os
calc_folder = 'tmp' 
for folder in [calc_folder]:
    if not os.path.exists(folder):
        os.makedirs(folder)


lammps_name=''
comm=None
log_file= calc_folder + '/lammps.log'
cmd_args = ['-echo', 'log', '-log', log_file,
            '-screen', 'none', '-nocite']
lmp = lammps(lammps_name, cmd_args, comm)
parameters = ["mass * 1.",
              "pair_style tersoff",
              "pair_coeff * * SiCGe.tersoff C",
             ]


print("launch the LAMMPS calculator")
s, eng = run_lammpslib(ase_struc, lmp, parameters, method='opt', path=calc_folder)
print(eng)

print("launch the GULP calculator")
s, eng, time, error = gulp_opt(s, ff='tersoff.lib', symmetrize=False)
print(eng)
