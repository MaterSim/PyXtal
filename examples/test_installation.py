import pyxtal
print("PyXtal: ", pyxtal.__version__)

import ase
print("ase: ", ase.__version__)

from pyxtal import pyxtal
print("Using PyXtal to generate structure")
struc = pyxtal()
struc.from_random(3, 227, ["C"], [8], sites=[['8a']])
print(struc)
print("Convert PyXtal structure to ASE")
ase_struc = struc.to_ase()

calc_folder = 'tmp'
import os
if not os.path.exists(calc_folder):
    os.mkdir(calc_folder)

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

s = opt_lammpslib(ase_struc, lmp, parameters, path=calc_folder, opt_cell=True)
if abs(s.get_potential_energy())<1e-8:
    cell=s.get_cell()
    pos = s.get_scaled_positions()
    s.set_cell(cell*0.8)
    s.set_scaled_positions(pos)
    s = opt_lammpslib(ase_struc, lmp, parameters, path=calc_folder, opt_cell=True)

print("launch the GULP calculator")
from pyxtal.interface.gulp import single_optimize as opt_gulp
s, eng, time, error = opt_gulp(ase_struc, ff='tersoff.lib', path=calc_folder, clean=False)
print(eng*sum(s.numIons))
#s.to_ase().write('1.vasp', format='vasp', vasp5=True, direct=True)
