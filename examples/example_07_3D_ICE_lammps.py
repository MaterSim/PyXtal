from pyxtal.crystal import Lattice
from pyxtal.molecular_crystal import molecular_crystal
from pyxtal.interface.lammpslib import LAMMPSlib
from pyxtal.interface.mushybox import mushybox
from lammps import lammps
from spglib import get_symmetry_dataset
from ase.optimize.fire import FIRE
from ase import Atoms
from ase.optimize import BFGS
from ase.build import sort
import logging
from random import randint
import pandas as pd
import os
import numpy as np

# set up for lammps
lammps_name = ""
comm = None
log_file = "lammps.log"
cmd_args = ["-echo", "log", "-log", log_file, "-screen", "none", "-nocite"]
lmp = lammps(lammps_name, cmd_args, comm)

logging.basicConfig(
    format="%(asctime)s :: %(message)s", filename="results.log", level=logging.INFO
)


parameters = [
    "atom_style full",
    "pair_style      lj/cut/tip4p/long 1 2 1 1 0.278072379 17.007",
    "bond_style      class2 ",
    "angle_style     harmonic",
    "kspace_style pppm/tip4p 0.0001",
    "read_data tmp/data.lammps",
    "pair_coeff  2  2 0 0",
    "pair_coeff  1  2 0 0",
    "pair_coeff  1  1  0.000295147 5.96946",
    "neighbor 2.0 bin",
    "min_style cg",
    "minimize 1e-6 1e-6 10000 10000",
]

# setup folders
calc_folder = "tmp"  # store tmp files for lammps and lasp
out_folder = "out"  # store tmp files for lammps and lasp
for folder in [calc_folder, out_folder]:
    if not os.path.exists(folder):
        os.makedirs(folder)


# para = Lattice.from_para(4.45, 7.70, 7.28, 90, 90, 90)

# Here we generate many random structures and optimize them
data = {  # "ID":[],
    "Sym": [],
    "Eng": [],
    "abc": [],
    "Volume": [],
}

for i in range(100):
    while True:
        sg = randint(16, 191)
        crystal = molecular_crystal(sg, ["H2O"], [4], 1.0)  # , lattice=para)
        if crystal.valid:
            struc = Atoms(
                crystal.spg_struct[2],
                cell=crystal.spg_struct[0],
                scaled_positions=crystal.spg_struct[1],
            )
            break
    # struc = struc.repeat((2,2,2))
    lammps = LAMMPSlib(lmp=lmp, lmpcmds=parameters, mol=True)
    struc.set_calculator(lammps)
    box = mushybox(struc)
    dyn = FIRE(box)
    dyn.run(fmax=0.01, steps=200)
    dyn = BFGS(box)
    dyn.run(fmax=0.01, steps=200)
    Eng = struc.get_potential_energy() * 96 / len(struc) * 3
    Vol = struc.get_volume() / len(struc) * 3
    stress = np.max(struc.get_stress())
    struc = sort(struc)
    try:
        spg = get_symmetry_dataset(struc, symprec=1e-1)["number"]
    except:
        spg = 1
    struc.write(out_folder + "/" + str(i) + ".vasp", format="vasp", vasp5=True)
    logging.info(
        "{:4d} Spg: {:4d} Eng: {:8.4f} Vol: {:8.4f} Stress: {:5.2f}".format(
            i, spg, Eng, Vol, stress
        )
    )
    abc = struc.get_cell_lengths_and_angles()
    abc = [int(i * 1000) / 1000 for i in abc]
    # data['ID'].append(i)
    data["Sym"].append(spg)
    data["Eng"].append(Eng)
    data["abc"].append(abc)
    data["Volume"].append(Vol)

df = pd.DataFrame(data)
df = df.sort_values(["Eng", "Volume"], ascending=[True, True])
print(df)
