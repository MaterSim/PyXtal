from pymatgen.core.operations import SymmOp
import pandas as pd
import numpy as np
from structure import *
from ast import literal_eval as eval

"""mylist = [None]
print("Adding space group:")
for sg in range(1, 231):
    print(sg)
    wyckoffs = get_wyckoff_positions(sg)
    mylist.append([])
    for x in wyckoffs:
        #x is a sorted list of positions
        for y in x:
            #y is a wyckoff position
            mylist[-1].append([])
            for z in y:
                #z is a SymmOp
                ops = site_symm(z, wyckoffs[0][0])
                oplist = []
                for op in ops:
                    oplist.append(op.as_xyz_string())
                mylist[-1][-1].append(oplist)
print(len(mylist))
df = pd.DataFrame(data=mylist)
df.to_csv("wyckoff_symmetry.csv")"""

wyckoff_df = pd.read_csv("wyckoff_list.csv")


def get_wyckoffs(sg):
    wyckoff_strings = eval(wyckoff_df["0"][sg])
    wyckoffs = []
    for x in wyckoff_strings:
        wyckoffs.append([])
        for y in x:
            wyckoffs[-1].append(SymmOp.from_xyz_string(y))
    return wyckoffs


print(type(get_wyckoffs(2)[0][0]))
print(get_wyckoffs(2))
