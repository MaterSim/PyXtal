# %%

from pymatgen.core import Molecule

from pyxtal import pyxtal

Li = Molecule(["Li"], [[0.0, 0.0, 0.0]])
coords = [
    [0.000000, 0.000000, 0.000000],
    [1.200000, 1.200000, -1.200000],
    [1.200000, -1.200000, 1.200000],
    [-1.200000, 1.200000, 1.200000],
    [-1.200000, -1.200000, -1.200000],
]
ps4 = Molecule(["P", "S", "S", "S", "S"], coords)

for _i in range(3):
    struc = pyxtal(molecular=True)
    struc.from_random(3, 10, [Li, ps4], [6, 2], 1.2, conventional=False)
    if struc.valid:
        assert len(struc.to_pymatgen()) == 16
