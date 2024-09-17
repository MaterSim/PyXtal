from juliacall import Main as jl
# Import the Julia package manager
jl.seval("using CrystalNets")
print("Success CrystalNets")

data = """
12,4.6997,7.1199,3.6592,1.5708,2.1366,1.5708,0,0.438,0.8075,0.6503,2,0.0,0.9344,0.5,-1,-1.0,-1.0,-1.0,-1,-1.0,-1.0,-1.0,-1,-1.0,-1.0,-1.0,-1,-1.0,-1.0,-1.0,-1,-1.0,-1.0,-1.0,-1,-1.0,-1.0,-1.0
43,10.1657,9.213,8.7536,1.5708,1.5708,1.5708,0,0.8343,0.3453,0.1644,0,0.3345,0.7476,0.8868,0,0.7693,0.9651,0.6999,0,0.8732,0.0387,0.0689,0,0.7571,0.9465,0.4255,0,0.9706,0.7866,0.5912,-1,-1.0,-1.0,-1.0,-1,-1.0,-1.0,-1.0
"""

import numpy as np
from pyxtal import pyxtal
from pyxtal.lego.builder import mof_builder

lines = data.strip().splitlines()
dat = [list(map(float, line.split(','))) for line in lines]

xtal = pyxtal(); xtal.from_prototype('graphite')
cif_file = xtal.to_pymatgen()

builder = mof_builder(['C'], [1], verbose=True)
builder.set_descriptor_calculator(mykwargs={'rcut': 2.0})
builder.set_reference_enviroments(cif_file)
builder.set_criteria(CN={'C': [3]})
print(builder)


xtals = []
for i in range(len(dat)):
    rep = dat[i]
    print(i, rep)
    if (len(rep)-7)%4 == 1: rep = rep[1:]
    xtal = pyxtal()
    xtal.from_tabular_representation(rep, normalize=False)
    if xtal.valid and len(xtal.atom_sites) > 0:
        xtals.append(xtal)

xtals_opt = builder.optimize_xtals(xtals,
                                   minimizers=[('Nelder-Mead', 100),
                                               ('L-BFGS-B', 400)],
                                   early_quit=0.15)
xtals = []
builder.db.update_row_topology(overwrite=True, prefix='tmp')
builder.db.clean_structures_spg_topology(dim=3)
