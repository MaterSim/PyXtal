from pyxtal.symmetry import Group
import numpy as np
from pyxtal.lattice import Lattice

cell = Lattice.from_para(9.395000, 9.395000, 8.350000, 90, 90, 120, ltype='hexagonal')
pt = [0.1333, 0.1496, 0.969]
p0 = np.dot(pt, cell.matrix)


for sg in range(143, 195):
#for sg in [149]:
    wp = Group(sg)[0]
    print(sg, wp.symbol, wp.euclidean)
    for i in range(len(wp)):
        op0 = wp[i]
        rot, _ = wp.get_euclidean_rotation(i)
        #print(op0.affine_matrix[:3, :3].T)
        #print(rot)
        p1 = op0.operate(pt)
        p2 = np.dot(np.dot(p0, rot), cell.inv_matrix) + op0.translation_vector
        diff = p1-p2
        diff -= np.round(diff)
        if np.linalg.norm(diff) > 0.02:
            res = '{:2d} {:28s}'.format(i, op0.as_xyz_string())
            res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p1)
            res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p2)
            res += '{:6.3f} {:6.3f} {:6.3f}'.format(*diff)
            print(res)
