from structure import get_wyckoff_symmetry
import numpy as np
from numpy import allclose
from numpy import isclose
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import generate_full_symmops
import pandas as pd
from copy import deepcopy
from matrix import *
from math import sqrt

"""site_symmetry = [None]
#site_symm is stored by space group number starting with 1 (site_symm[1] is P1))
print("Calculating space group:")
for sg in range(1, 231):
    print(sg)
    site_symmetry.append([])
    symmetry = get_wyckoff_symmetry(sg)
    #Get site symmetry for every point in each wp, and store
    for i, wp in enumerate(symmetry):
        site_symmetry[sg].append([])
        for j, point in enumerate(wp):
            site_symmetry[sg][-1].append([])
            symm_point_unedited = symmetry[i][j]
            found_params = []
            symm_point_partial = []
            symm_point = []
            orthogonal = True
            #Check if operations are orthogonal; 3-fold and 6-fold operations are not
            for op in symm_point_unedited:
                m1 = np.dot(op.rotation_matrix, np.transpose(op.rotation_matrix))
                m2 = np.dot(np.transpose(op.rotation_matrix), op.rotation_matrix)
                if ( not allclose(m1, np.identity(3)) ) or ( not allclose(m2, np.identity(3)) ):
                    rotation_order = 0
                    d = np.linalg.det(op.rotation_matrix)
                    e = np.linalg.eig(op.rotation_matrix)[1]
                    e = np.transpose(e)
                    for a in e:
                        if allclose(a, [1,0,0]) or allclose(a, [0,1,0]) or allclose(a, [0,0,1]):
                            op_axis = a
                        else:
                            print("Error: could not find axis")
                    if isclose(d, 1):
                        op_type = "rotation"
                    elif isclose(d, -1):
                        op_type = "improper rotation"
                    if op_type == "rotation":
                        newop = deepcopy(op)
                        for n in range(1, 7):
                            newop = (newop*op)
                            if allclose(newop.rotation_matrix, np.identity(3), rtol=1e-2):
                                rotation_order = n + 1
                                break
                    elif op_type == "improper rotation":
                        #We only want the order of the rotational part of op,
                        #So we multiply op by -1 for rotoinversions
                        op_1 = SymmOp.from_rotation_and_translation(op.rotation_matrix*-1,[0,0,0])
                        new_op = deepcopy(op_1)
                        for n in range(1, 7):
                            newop = (newop*op_1)
                            if allclose(newop.rotation_matrix, np.identity(3)):
                                rotation_order = n + 1
                                break
                    if rotation_order == 0:
                        print("Error: could not convert to orthogonal operation:")
                        print(op)
                        symm_point_partial.append(op)
                    else:
                        params = [rotation_order, op_type, op_axis]
                        found_params.append(params)
                else:
                    symm_point_partial.append(op)
            #Add 3-fold and 6-fold symmetry back in using absolute coordinates
            for params in found_params:
                order = params[0]
                op_type = params[1]
                if op_type == "rotation" and (order == 3 or order == 6):
                    m = aa2matrix(op_axis, (2*pi)/order)
                elif op_type == "improper rotation" and (order == 3 or order == 6):
                    m = aa2matrix(op_axis, (2*pi)/order)
                    m[2] *= -1
                symm_point_partial.append(SymmOp.from_rotation_and_translation(m, [0,0,0]))
            #print("Generating...")
            symm_point = generate_full_symmops(symm_point_partial, 1e-2)
            #print("Done")
            for op in symm_point:
                site_symmetry[sg][-1][-1].append(op.as_xyz_string())"""

P = SymmOp.from_rotation_and_translation(
    [[1, -0.5, 0], [0, sqrt(3) / 2, 0], [0, 0, 1]], [0, 0, 0]
)

site_symmetry = [None]
# site_symm is stored by space group number starting with 1 (site_symm[1] is P1))
print("Calculating space group:")
for sg in range(1, 231):
    print(sg)
    site_symmetry.append([])
    symm_sg = get_wyckoff_symmetry(sg)
    # Get site symmetry for every point in each wp, and store
    for i, symm_wp in enumerate(symm_sg):
        site_symmetry[sg].append([])
        for j, symm_point in enumerate(symm_wp):
            site_symmetry[sg][-1].append([])
            for op in symm_point:
                new_op = P * op * P.inverse
                new_op = new_op.from_rotation_and_translation(
                    new_op.rotation_matrix, [0, 0, 0]
                )
                site_symmetry[sg][-1][-1].append(new_op.as_xyz_string())


print("Saving file...")
array = np.array(site_symmetry)
df = pd.DataFrame(data=array)
df.to_csv("wyckoff_symmetry_molecular.csv")
print("File saved")
