from pymatgen.core.operations import SymmOp
import pandas as pd
from pyxtal.crystal import *


site_symmetry = [None]
#site_symm is stored by space group number starting with 1 (site_symm[1] is P1))
print("Generating site symmetry to store in ")
print("Calculating space group:")

#Change 231 to 81 for 2D, 76 for 1D
for num in range(1, 231):
    print(num)
    site_symmetry.append([])

    #Replace get_wyckoffs with get_layer or get_rod
    wyckoffs = get_wyckoffs(num)

    gen_pos = wyckoffs[0]
    for wp in wyckoffs:
        site_symmetry[-1].append([])
        for point in wp:
            site_symmetry[-1][-1].append([])
            new = site_symm(point, gen_pos)
            for y in new:
                site_symmetry[-1][-1][-1].append(y.as_xyz_string())

#Rename this file
path = "wyckoff_symmetry_new.csv"

print("Saving file to "+path+" ...")
array = np.array(site_symmetry)
print(type(array[1][0][0][0]))
df = pd.DataFrame(data=array)
df.to_csv(path)
print("Done")
