from pyxtal.crystal import Lattice
from pyxtal.molecular_crystal import molecular_crystal

lat0 = Lattice.from_para(21.339, 5.714, 8.900, 90, 90.529, 90, ltype='monoclinic')
for i in range(2):
    lat = lat0.swap_axis(random=True)
    lat = lat.swap_angle()
    struc = molecular_crystal(4, ["OFIXUX"], [2], lattice = lat)
    print(i, struc.valid)
