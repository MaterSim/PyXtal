from pyxtal.crystal import random_crystal
from pyxtal.interface.vasp import optimize

# from vasp import optimize
import json
from monty.serialization import MontyEncoder, MontyDecoder, loadfn
import os
from random import randint
from time import time
import shutil
import pandas as pd
from tabulate import tabulate
import warnings

warnings.filterwarnings("ignore")

pyxtal_verbosity = 0

"""
This is a script to 
1, generate random structures
2, perform multiple steps of optmization with ASE-VASP
3, optionally, you can also save the configuration to json file for other purposes

Requirement:
You must have ASE installed and can call vasp from ASE
"""


def dump_json(strucs, energies, times, json_file1, json_file2):
    struc_info = []
    for struc, energy, time in zip(strucs, energies, times):
        dic = {
            "formula": struc.get_chemical_formula(),
            "coords": struc.get_scaled_positions(),
            "lattice": struc.cell,
            "elements": struc.get_chemical_symbols(),
            "energy": energy,
            "time": time,
        }
        struc_info.append(dic)
    with open(json_file1, "a+") as f:
        json.dump(struc_info, f, cls=MontyEncoder, indent=2)
    with open(json_file2, "a+") as f:
        json.dump(struc_info[-1], f, cls=MontyEncoder, indent=2)


def write_summary(file1):
    content = loadfn(file1, cls=MontyDecoder)
    col_name = {
        "Formula": formula,
        "Energy": energy,
        "Time": time,
    }
    for dct in content:
        formula.append(dct["formula"])
        energy.append(dct["energy"] / len(dct["elements"]))
        time.append(dct["time"])
    df = pd.DataFrame(col_name)
    print(tabulate(df, headers="keys", tablefmt="psql"))


N = 2
elements = {"C": [4, 6]}
factor = 1.1
json1 = "train.json"
json2 = "csp.json"
modes = ["C", "P"]

dir0 = os.getcwd()
t0 = time()

for i in range(N):
    run = True
    while run:
        sg = randint(2, 230)
        species = []
        numIons = []
        for ele in elements.keys():
            species.append(ele)
            if len(elements[ele]) == 2:
                numIons.append(randint(elements[ele][0], elements[ele][1]))
            else:
                numIons.append(elements[ele])

        crystal = random_crystal(sg, species, numIons, factor)

        if crystal.valid:
            struc = crystal.struct
            run = False

    print(
        "SG requested: {0:3d} Vol: {1:6.2f} Time: {2:6.1f} mins".format(
            sg, struc.volume, (time() - t0) / 60.0
        )
    )
    dir1 = str(i) + "-" + str(struc.formula).replace(" ", "")
    [strucs, energies, times] = optimize(struc, dir1, modes=modes)

    os.chdir(dir0)
    if len(strucs) == len(modes):
        dump_json(strucs, energies, times, json1, json2)
    shutil.rmtree(dir1)

# write_summary(json2)
