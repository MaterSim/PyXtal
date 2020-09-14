from pyxtal.molecular_crystal import molecular_crystal
from spglib import get_symmetry_dataset
from random import choice
#from pyxtal.interface.lammpslib import LAMMPSlib
from lammpslib import run_lammpslib, optimize_lammpslib
from pyxtal.interface.gulp import single_optimize as gulp_opt
from ase.db import connect
from lammps import lammps
import numpy as np
import os
import logging

log_file = "07results.log"
if os.path.exists(log_file):
    os.remove(log_file)

logging.basicConfig(format="%(asctime)s :: %(message)s", filename=log_file, level=logging.INFO)

calc_folder = '07-tmp' #store tmp files for lammps and lasp
out_folder = '07-out'  #store tmp files for lammps and lasp
for folder in [calc_folder , out_folder]:
    if not os.path.exists(folder):
        os.makedirs(folder)

lammps_name=''
comm=None
log_file= calc_folder + '/lammps.log'
cmd_args = ['-echo', 'log', '-log', log_file,
            '-screen', 'none', '-nocite']
lmp = lammps(lammps_name, cmd_args, comm)

parameters = [
    "atom_style full",
    "pair_style      lj/cut/tip4p/long 1 2 1 1 0.278072379 17.007",
    "bond_style      class2 ",
    "angle_style     harmonic",
    "kspace_style pppm/tip4p 0.0001",
    "read_data 07-tmp/data.lammps",
    "pair_coeff  2  2 0 0",
    "pair_coeff  1  2 0 0",
    "pair_coeff  1  1  0.000295147 5.96946",
]

for i in range(100):
    while True:
        sg, numIons = choice(range(3,231)), choice(range(4,6))
        struc = molecular_crystal(19, ["H2O"], [8]) 
        if struc.valid:
            break
    ase_struc = struc.to_ase(resort=False)

    #ase_struc = ase_struc.repeat((2,2,2))
    print("init", struc.group.symbol, struc.numMols, ase_struc.get_cell_lengths_and_angles())
    s, eng = run_lammpslib(ase_struc, lmp, parameters, molecule=True, method='opt', path=calc_folder)
    s = optimize_lammpslib(s, lmp, parameters, molecule=True, fmax=0.01, path=calc_folder)
    s, eng = run_lammpslib(s, lmp, parameters, molecule=True, method='opt', path=calc_folder)

    Eng = s.get_potential_energy() * 96 / len(s) * 3
    Vol = s.get_volume() / len(s) * 3
    try:
        spg = get_symmetry_dataset(s, symprec=1e-1)['international']
    except:
        spg = 1
    
    logging.info("{:4d} {:6.3f} eV/atom  {:6.3f} A^3  {:10s} --> {:10s}".format(i, Eng, Vol, struc.group.symbol, spg))
    print("{:4d} {:6.3f} eV/atom  {:6.3f} A^3  {:10s} --> {:10s}".format(i, Eng, Vol, struc.group.symbol, spg))
    permutation = np.argsort(s.numbers)
    s = s[permutation]
    s.write(out_folder + '/' + str(i)+'.vasp', format='vasp', vasp5=True, direct=True)
