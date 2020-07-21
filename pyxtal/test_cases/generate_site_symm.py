from pymatgen.core.operations import SymmOp
import pandas as pd
from pyxtal.crystal import *

"""
Generate the site symmetry csv files
"""

# EDIT
path = "point_symmetry_new.csv"
PBC = [0, 0, 0]
# Change 231 to 81 for 2D, 76 for 1D
maxn = 33

site_symmetry = [None]

# site_symm is stored by space group number starting with 1 (site_symm[1] is P1))
print("Generating site symmetry to store in ")
print("Calculating group:")


for num in range(1, maxn):
    print(num)
    site_symmetry.append([])

    # EDIT
    wyckoffs = get_point(num)

    gen_pos = wyckoffs[0]
    for wp in wyckoffs:
        site_symmetry[-1].append([])
        for point in wp:
            site_symmetry[-1][-1].append([])
            new = site_symm(point, gen_pos, PBC=PBC)
            for y in new:
                site_symmetry[-1][-1][-1].append(y.as_xyz_string())

print("Saving file to " + path + " ...")
array = np.array(site_symmetry)
if type(array[1][0][0][0]) != str:
    print("Error: data incorrectly stored.")
df = pd.DataFrame(data=array)
df.to_csv(path)
print("Done")
