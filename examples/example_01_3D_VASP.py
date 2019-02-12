from pyxtal.crystal import random_crystal
from pyxtal.interface.vasp import optimize
import os
from random import randint
from time import time

"""
This is a script to 
1, generate random structures
2, perform multiple steps of optmization with ASE-VASP
3, optionally, you can also save the configuration to json file for other purposes

Requirement:
You must have ASE installed and can call vasp from ASE
"""

N = 100
elements = {'C': [8, 16]} 
factor = 1.1

dir0 = os.getcwd()
t0 = time()

for i in range(N):
    os.chdir(dir0)
    run = True
    while run:
        sg = randint(2, 230)
        species = []
        numIons = []
        for ele in elements.keys():
            species.append(ele)
            if len(elements[ele])==2:
                numIons.append(randint(elements[ele][0], elements[ele][1]))
            else:
                numIons.append(elements[ele])

        crystal = random_crystal(sg, species, numIons, factor)

        if crystal.valid:
            struc = crystal.struct
            run = False

    print('SG requested: {0:3d} Vol: {1:6.2f} Time: {2:6.1f} mins'.\
            format(sg, struc.volume, (time()- t0)/60.0))
    try:
        dir1 = str(i) + '-' + str(struc.formula).replace(" ","")
        print(dir1)
        [struc, energy, time] = optimize(struc, dir1)
        if len(struc) == 5:
            dump_json()

    except:
        print('something is wrong')

