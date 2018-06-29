from bug import *
from pymatgen.core.operations import SymmOp
import pandas as pd

site_symmetry = [None]
#site_symm is stored by space group number starting with 1 (site_symm[1] is P1))
print("Calculating space group:")
for sg in range(1, 231):
    print(sg)
    site_symmetry.append([])
    wyckoffs = get_wyckoffs(sg)
    gen_pos = wyckoffs[0]
    #Get site symmetry for every point in each wp, and store
    for wp in wyckoffs:
        site_symmetry[sg].append([])
        for point in wp:
            site_symmetry[sg][-1].append([])
            new = site_symm(point, gen_pos)
            for y in new:
                site_symmetry[sg][-1][-1].append(y.as_xyz_string())

print("Saving file...")
array = np.array(site_symmetry)
print(type(array[1][0][0][0]))
df = pd.DataFrame(data=array)
df.to_csv("wyckoff_symmetry_new.csv")
