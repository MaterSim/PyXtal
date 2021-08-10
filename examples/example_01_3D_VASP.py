from pyxtal import pyxtal
from pyxtal.interface.vasp import optimize
from ase.db import connect
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
elements = {"C": [2, 4]}
levels = [0, 2] #, 3]
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

    #try:
    # relax the structures with levels 0, 2, 3
    struc, energy, cputime, error = optimize(crystal, path=dir1, levels=levels)

    if not error: # successful calculation
        s = struc.to_ase()
        with connect(filename) as db:
            kvp = {"spg": struc.group.symbol, 
                   "dft_energy": energy/len(s), 
                  }
            db.write(s, key_value_pairs=kvp)
        cputime /= 60.0
        strs = "{:3d}".format(i)
        strs += " {:12s} -> {:12s}".format(crystal.group.symbol, struc.group.symbol)
        strs += " {:6.3f} eV/atom".format(energy/len(s))
        strs += " {:6.2f} min".format(cputime)
        print(strs)

        #from pyxtal.interface.vasp import single_point
        #print(single_point(struc, path=dir1))
    #except:
    #    print("vasp calculation is wrong")
