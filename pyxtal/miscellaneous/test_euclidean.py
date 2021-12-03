from pyxtal.symmetry import Group
import numpy as np
from pyxtal.lattice import Lattice


def test(spg, pt, cell):

    p0 = np.dot(pt, cell.matrix)
    
    for sg in spg:
        wp = Group(sg)[0]
        print(sg, wp.symbol, wp.euclidean)
        for i in range(len(wp)):
            op0 = wp[i]
            p1 = op0.operate(pt)
            
            op1 = wp.get_euclidean_rotation(cell.matrix, i)
            if wp.euclidean:
                p2 = np.dot(op1.operate(p0), cell.inv_matrix)
            else:
                p2 = np.dot(op1.apply_rotation_only(p0), cell.inv_matrix) 
                p2 += op1.translation_vector

            diff = p1-p2
            diff -= np.round(diff)
            if np.linalg.norm(diff) > 0.02:
                res = '{:2d} {:28s}'.format(i, op0.as_xyz_string())
                res += ' {:28s}'.format(op1.as_xyz_string())
                res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p1)
                res += '{:6.3f} {:6.3f} {:6.3f} -> '.format(*p2)
                res += '{:6.3f} {:6.3f} {:6.3f}'.format(*diff)
                print(res)


pt = [0.1333, 0.1496, 0.969]

cell = Lattice.from_para(9.395000, 7.395000, 8.350000, 91, 101, 92, ltype='triclinic')
test(range(1, 3), pt, cell)

cell = Lattice.from_para(9.395000, 7.395000, 8.350000, 90, 101, 90, ltype='monoclinic')
test(range(3, 16), pt, cell)

cell = Lattice.from_para(9.395000, 7.395000, 8.350000, 90, 90, 90, ltype='orthorhombic')
test(range(3, 16), pt, cell)

cell = Lattice.from_para(9.395000, 7.395000, 8.350000, 90, 90, 90, ltype='orthorhombic')
test(range(16, 74), pt, cell)

cell = Lattice.from_para(9.395000, 9.395000, 8.350000, 90, 90, 90, ltype='tetragonal')
test(range(74, 143), pt, cell)

cell = Lattice.from_para(9.395000, 9.395000, 8.350000, 90, 90, 120, ltype='hexagonal')
test(range(143, 195), pt, cell)

cell = Lattice.from_para(9.395000, 9.395000, 9.3950000, 90, 90, 90, ltype='cubic')
test(range(195, 231), pt, cell)
