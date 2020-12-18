from pyxtal import pyxtal
from pyxtal.interface.vasp import optimize
from ase.db import connect
import os
from random import randint
from time import time
import warnings

warnings.filterwarnings("ignore")

"""
This is a script to 
1, generate random structures
2, perform multiple steps of optmization with ASE-VASP

Requirement:
You must have ASE installed and can call vasp from ASE
"""
N = 10
elements = {"C": [4, 6]}
levels = [0, 2, 3]
dir1 = "Calc"
filename = "C-VASP.db"

for i in range(N):
    t0 = time()
    while True:
        sg = randint(2, 230)
        species = []
        numIons = []
        for ele in elements.keys():
            species.append(ele)
            if len(elements[ele]) == 2:
                num = randint(elements[ele][0], elements[ele][1])
                numIons.append(num)
            else:
                numIons.append(elements[ele])

        crystal = pyxtal()
        crystal.from_random(3, sg, species, numIons, force_pass=True)

        if crystal.valid:
            break

    try:
        # relax the structures with levels 0, 2, 3
        strucs, energies, times = optimize(crystal, dir1, levels)

        if len(strucs) == len(levels): # successful calculation
            s = strucs[-1].to_ase()
            with connect(filename) as db:
                kvp = {"spg": strucs[-1].group.symbol, 
                       "dft_energy": energies[-1], 
                       "energy": energies[-1], 
                      }
                db.write(s, key_value_pairs=kvp)
            t = (time() - t0) / 60.0
            strs = "{:3d}".format(i)
            strs += " {:12s} -> {:12s}".format(crystal.group.symbol, strucs[-1].group.symbol)
            strs += " {:6.3f} eV/atom".format(energies[-1])
            strs += " {:6.2f} min".format(t)
            print(strs)

            #from pyxtal.interface.vasp import single_point
            #print(single_point(strucs[-1], dir0=dir1))
    except:
        print("vasp calculation is wrong")
