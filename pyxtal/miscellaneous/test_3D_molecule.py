import os
from random import randint
from time import time

from pymatgen.io.cif import CifWriter
from pyxtal.molecular_crystal import molecular_crystal

mols = ["CH4", "H2O", "NH3", "urea", "benzene", "roy", "aspirin", "pentacene", "C60"]
filename = "out.cif"
if os.path.isfile(filename):
    os.remove(filename)

for mol in mols:
    for i in range(10):
        run = True
        while run:
            sg = randint(4, 191)
            start = time()
            rand_crystal = molecular_crystal(sg, [mol], [4], 1.0)
            if rand_crystal.valid:
                run = False
                print(
                    "Mol:{0:10s}  SG:{1:3d} Time:{2:4.2f} seconds, N_attempts:{3:4d} Vol: {4:6.2f}".format(
                        mol,
                        sg,
                        time() - start,
                        rand_crystal.numattempts,
                        rand_crystal.volume,
                    )
                )
                content = str(CifWriter(rand_crystal.struct, symprec=0.1))
                with open(filename, "a+") as f:
                    f.writelines(content)
