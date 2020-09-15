from pyxtal.crystal import random_crystal_2D
from pyxtal.interface.lammpslib import run_lammpslib as lmp_run 
from pyxtal.interface.lammpslib import optimize_lammpslib as lmp_opt
from pyxtal.interface.gulp import single_optimize as gulp_opt
from spglib import get_symmetry_dataset
import pymatgen.analysis.structure_matcher as sm
from pymatgen.io.ase import AseAtomsAdaptor
from random import choice
from ase.db import connect
from lammps import lammps
import numpy as np
import os
import logging

def new_struc(db, s, eng):
    s_pmg = AseAtomsAdaptor().get_structure(s)
    for row in db.select():
        ref = db.get_atoms(id=row.id)
        if abs(eng-row.ff_energy) < 1e-2:
            ref_pmg = AseAtomsAdaptor().get_structure(ref)
            if sm.StructureMatcher().fit(s_pmg, ref_pmg):
                return False
    return True

log_file = "06-results.log"
if os.path.exists(log_file):
    os.remove(log_file)

logging.basicConfig(format="%(asctime)s| %(message)s", filename=log_file, level=logging.INFO)


calc_folder = '06-tmp/' #store tmp files for lammps and lasp
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

strain = np.zeros([3,3]) 
strain[:2,:2] += 1
# Here we do plane groups
pgs = [3, 11, 12, 13, 23, 24, 25, 26, 49, 55, 56, 65, 69, 70, 73, 75]

filename = '06.db'
with connect(filename) as db:
    for i in range(100):
        while True:
            sg, numIons = choice(pgs), choice(range(4,24))
            struc = random_crystal_2D(sg, ["C"], [numIons], thickness=0)
            if struc.valid:
                ase_struc = struc.to_ase()
                break
        s, eng = lmp_run(ase_struc, lmp, parameters, method='opt', path=calc_folder)
        # relax the structure with multiple steps
        for m in range(2):
            s = lmp_opt(s, lmp, parameters, strain=strain, fmax=0.01, path=calc_folder)
            s, eng = lmp_run(s, lmp, parameters, method='opt', path=calc_folder)

        s, eng, _, error = gulp_opt(s, ff='tersoff.lib', path=calc_folder, symmetrize=False)
        if error:
            spg = get_symmetry_dataset(s, symprec=1e-1)['international']
            pos = s.get_positions()[:,2]
            thickness = max(pos) - min(pos) 
            new = new_struc(db, s, eng)
            strs = "{:4d} {:6.3f} eV/atom *{:2d} {:6.3f} {:10s}->{:10s} {:}".format(\
                    i, eng/len(s), len(s), thickness, struc.group.symbol, spg, new)
            logging.info(strs)
            print(strs)
            if new_struc(db, s, eng/len(s)):
                area = np.linalg.det(s.get_cell()[:2,:2])/len(s)
                kvp = {"spg": spg, 
                       "ff_energy": eng/len(s), 
                       "thickness": thickness, 
                       "area": area}
                db.write(s, key_value_pairs=kvp)

