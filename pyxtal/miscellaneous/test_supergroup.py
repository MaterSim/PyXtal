import importlib.resources
from time import time

from pyxtal import pyxtal
from pyxtal.supergroup import supergroups

data = {
    "GeF2": 62,
    "NiS-Cm": 160,
    "lt_quartz": 180,
    "BTO-Amm2": 221,
    "BTO": 221,
    "lt_cristobalite": 227,
    "NaSb3F10": 194,
    "NbO2": 141,
    # "MPWO": 225,
}

with importlib.resources.as_file(importlib.resources.files("pyxtal") / "database" / "cifs") as path:
    cif_path = path

for cif in data:
    t0 = time()
    print("===============", cif, "===============")
    for tol in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]:
        s = pyxtal()
        s.from_seed(str(cif_path / f"{cif}.cif"), tol=tol)
        print("tol", tol, s.group.number)

    s = pyxtal()
    s.from_seed(str(cif_path / f"{cif}.cif"))
    if isinstance(data[cif], list):
        sup = supergroups(s, path=data[cif], show=True, max_per_G=2500)
    else:
        sup = supergroups(s, G=data[cif], show=False, max_per_G=2500)
    print(sup)
    print(f"{time() - t0:6.3f} seconds")
