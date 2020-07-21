from structure import *
from pymatgen.core.operations import SymmOp

print("---Site Symmetry for a point in a space group---")
sg = input("Space group Number (1-230): ")
point = input("Point (can have x,y,z variables): ")
gen_pos = get_wyckoff_positions(int(sg))[0][0]
point = SymmOp.from_xyz_string(point)
mylist = site_symm(point, gen_pos)
print("Found symmetry:")
for x in mylist:
    print(x.as_xyz_string())
