"""
global optimizer
ACSALA
81 11.38  6.48 11.24  96.9 1 0 0.23 0.43 0.03  -44.6   25.0   34.4  -76.6   -5.2  171.5 0 -70594.48
"""
# ACSALA has 13 variables (4 cell + 3 xyz + 3 ori + 3 torsions)
import os
from scipy.stats import qmc
from pyxtal.representation import representation
from pyxtal.optimize.common import optimizer
from pyxtal.db import database
from pyxtal.lattice import Lattice
import pymatgen.analysis.structure_matcher as sm

w_dir = "tmp"
if not os.path.exists(w_dir):
    os.makedirs(w_dir)

db = database("pyxtal/database/test.db")
code = "ACSALA"
row = db.get_row(code)
xtal0 = db.get_pyxtal(code)
if xtal0.has_special_site():
    xtal0 = xtal0.to_subgroup()
smiles = row.mol_smi; print(smiles)
pmg0 = xtal0.to_pymatgen()
pmg0.remove_species("H")

# Prepare prm
c_info = row.data["charmm_info"]
with open(w_dir + "/pyxtal.prm", "w") as prm:
    prm.write(c_info["prm"])
with open(w_dir + "/pyxtal.rtf", "w") as rtf:
    rtf.write(c_info["rtf"])


# Prepare the sampling.....
bounds = xtal0.get_bounds([2.5, 25])
lb = [b[0] for b in bounds]
ub = [b[1] for b in bounds]
sampler = qmc.Sobol(d=xtal0.get_dof(), scramble=False)
sample = sampler.random_base2(m=10)
sample = qmc.scale(sample, lb, ub)

l_dof = xtal0.lattice.dof
est_volume = sum(n*mol.volume for n, mol in zip(xtal0.numMols, xtal0.molecules))
min_vol, max_vol = 0.25*est_volume, 5.0*est_volume

for i, des in enumerate(sample):
    cell = [xtal0.group.hall_number] + des[:l_dof].tolist()
    lat = Lattice.from_1d_representation(cell[1:], xtal0.lattice.ltype)

    if min_vol < lat.volume < max_vol:
        #print('reset volume0', lat.encode(), lat.volume)
        coef = (est_volume*0.5/lat.volume)**(1/3)
        lat = lat.scale(coef)
        cell = [xtal0.group.hall_number] + lat.encode()
        #print('reset volume1', lat.encode(), coef, lat.volume)

        x = [cell, [0] + des[l_dof:].tolist() + [False]]

        rep0 = representation(x, [smiles])
        xtal = rep0.to_pyxtal()
        if len(xtal.check_short_distances(r=0.75)) == 0:
            res = optimizer(xtal, c_info, w_dir, skip_ani=True)
            if res is not None:
                xtal, eng = res["xtal"], res["energy"]/sum(xtal.numMols)
                if eng < 1e+3:
                    rep1 = xtal.get_1D_representation()
                    strs = f'Opt {i:4d} ' + rep1.to_string(eng) + f' Vol {lat.volume:7.1f} => {xtal.lattice.volume:7.1f}'
                    pmg1 = xtal.to_pymatgen()
                    pmg1.remove_species("H")
                    if sm.StructureMatcher().fit(pmg0, pmg1):
                        strs += '++++++'
                    print(strs)

    if i % 10000 == 0:
        print("Progress", i, lat.volume, xtal0.lattice.volume)
