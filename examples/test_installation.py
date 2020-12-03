import pyxtal
print("PyXtal: ", pyxtal.__version__)

import ase
print("ase: ", ase.__version__)

from pyxtal import pyxtal
print("Using PyXtal to generate structure")
struc = pyxtal()
struc.from_random(3, 75, ["C"], [8], 0.9)
print("Convert PyXtal structure to ASE")
ase_struc = struc.to_ase()

calc_folder = 'tmp'

print("launch the GULP calculator")
from pyxtal.interface.gulp import single_optimize as gulp_opt
s, eng, time, error = gulp_opt(ase_struc, ff='tersoff.lib', path=calc_folder, clean=False)
print(eng)

print("launch the LAMMPS calculator")
# Set up lammps
from pyxtal.interface.lammpslib import opt_lammpslib
from lammps import lammps

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
s = opt_lammpslib(ase_struc, lmp, parameters, path=calc_folder)

# todo: figure out the results
