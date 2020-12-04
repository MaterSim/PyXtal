from pyxtal import pyxtal
from pyxtal.interface.lammpslib import run_lammpslib as lmp_run 
from pyxtal.interface.lammpslib import opt_lammpslib as lmp_opt
from spglib import get_symmetry_dataset
from random import choice
from ase.db import connect
from lammps import lammps
import numpy as np
import os
import logging

log_file= '07-results.log'
logging.basicConfig(format="%(asctime)s| %(message)s", filename=log_file, level=logging.INFO)

calc_folder = '07-tmp' #store tmp files for lammps and lasp
for folder in [calc_folder]:
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
    "read_data " + calc_folder + "/data.lammps",
    "kspace_style pppm/tip4p 0.0001",
    "pair_coeff  2  2 0 0",
    "pair_coeff  1  2 0 0",
    "pair_coeff  1  1  0.000295147 5.96946",
]

filename = '07.db'
logfile = calc_folder + '/log'
with connect(filename) as db:
    for i in range(100):
        while True:
            sg, numIons = choice(range(3,231)), choice(range(4,12))
            struc = pyxtal(molecular=True)
            struc.from_random(3, sg, ["H2O"], [numIons], force_pass=True) 
            if struc.valid:
                break
        s = struc.to_ase(resort=False)
        #s, _ = lmp_run(s, lmp, parameters, molecule=True, method='opt', path=calc_folder)
        s = lmp_opt(s, lmp, parameters, logfile=logfile, molecule=True, fmax=0.01, path=calc_folder, opt_cell=False, a=0, steps=100)
        s = lmp_opt(s, lmp, parameters, logfile=logfile, molecule=True, fmax=0.01, path=calc_folder, opt_cell=True, a=0, steps=100)
        s = lmp_opt(s, lmp, parameters, logfile=logfile, molecule=True, fmax=0.01, path=calc_folder, opt_cell=True)
        s, _ = lmp_run(s, lmp, parameters, molecule=True, method='opt', path=calc_folder)

        Eng = s.get_potential_energy() * 96 / len(s) * 3
        Vol = s.get_volume() / len(s) * 3
        try:
            spg = get_symmetry_dataset(s, symprec=1e-1)['international']
        except:
            spg = 'P1'
        strs = "{:4d} {:6.3f} eV/atom {:4d} {:6.3f} A^3 {:10s}-->{:10s}".format(\
                i, Eng, numIons, Vol, struc.group.symbol, spg)
        logging.info(strs)
        print(strs)
        permutation = np.argsort(s.numbers)
        s = s[permutation]
