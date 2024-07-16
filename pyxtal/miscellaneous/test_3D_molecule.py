import os
from time import time

import numpy as np
from pymatgen.io.cif import CifWriter

from pyxtal.molecular_crystal import molecular_crystal

rng = np.random.default_rng(0)
mols = ["CH4", "H2O", "NH3", "urea", "benzene", "roy", "aspirin", "pentacene", "C60"]
filename = "out.cif"
if os.path.isfile(filename):
    os.remove(filename)

for mol in mols:
    for _i in range(10):
        run = True
        while run:
            sg = rng.integers(4, 191)
            start = time()
            rand_crystal = molecular_crystal(sg, [mol], [4], 1.0)
            if rand_crystal.valid:
                run = False
                print(
                    f"Mol:{mol:10s}  SG:{sg:3d} Time:{time() - start:4.2f} seconds, N_attempts:{rand_crystal.numattempts:4d} Vol: {rand_crystal.volume:6.2f}"
                )
                content = str(CifWriter(rand_crystal.struct, symprec=0.1))
                with open(filename, "a+") as f:
                    f.writelines(content)
