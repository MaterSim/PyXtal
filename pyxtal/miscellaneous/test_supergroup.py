from pyxtal import pyxtal
from pkg_resources import resource_filename
from pyxtal.supergroup import supergroups
from time import time
data = {
        "GeF2": 62,
        "NiS-Cm": 160,
        "lt_quartz": 180,
        "BTO-Amm2": 221,
        "BTO": 221,
        "lt_cristobalite": 227,
        "NaSb3F10": 194,
        "NbO2": 141,
        #"MPWO": 225,
       }

cif_path = resource_filename("pyxtal", "database/cifs/")

for cif in data.keys():
    t0 = time()
    print("===============", cif, "===============")
    for tol in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]:
        s = pyxtal()
        s.from_seed(cif_path+cif+'.cif', tol=tol)
        print("tol", tol, s.group.number)

    s = pyxtal()
    s.from_seed(cif_path+cif+'.cif')
    if isinstance(data[cif], list):
        sup = supergroups(s, path=data[cif], show=True, max_per_G=2500)
    else:
        sup = supergroups(s, G=data[cif], show=False, max_per_G=2500)
    print(sup)
    print("{:6.3f} seconds".format(time()-t0))
