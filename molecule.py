from pymatgen.core.structure import Molecule
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.symmetry.analyzer import generate_full_symmops
import numpy as np
from numpy import isclose
from numpy.linalg import eigh
from numpy.linalg import det
from copy import deepcopy
from math import fabs
from random import random
from random import choice as choose
from matrix import *
from structure import get_wyckoff_symmetry

#from ase.build import molecule

identity = np.array([[1,0,0],[0,1,0],[0,0,1]])
inversion = np.array([[-1,0,0],[0,-1,0],[0,0,-1]])

def get_inertia_tensor(mol):
    '''
    Calculate the symmetric inertia tensor for a Molecule object
    '''
    mo = mol.get_centered_molecule()
    # Initialize elements of the inertia tensor
    I11 = I22 = I33 = I12 = I13 = I23 = 0.0
    for i in range(len(mo)):
        x, y, z = mo.cart_coords[i]
        m = mo[i].specie.number
        I11 += m * (y ** 2 + z ** 2)
        I22 += m * (x ** 2 + z ** 2)
        I33 += m * (x ** 2 + y ** 2)
        I12 += -m * x * y
        I13 += -m * x * z
        I23 += -m * y * z
    return np.array([[I11, I12, I13],
                  [I12, I22, I23],
                  [I13, I23, I33]])

def reoriented_molecule(mol, nested=False):
    '''
    Return a molecule reoriented so that its principle axes
    are aligned with the identity matrix, and the matrix P
    used to rotate the molecule into this orientation
    '''
    def reorient(mol):
        new_mol = mol.get_centered_molecule()
        A = get_inertia_tensor(new_mol)
        #Store the eigenvectors of the inertia tensor
        P = np.transpose(eigh(A)[1])
        if det(P) < 0:
            P[0] *= -1
        #reorient the molecule
        P = SymmOp.from_rotation_and_translation(P,[0,0,0])
        new_mol.apply_operation(P)
        #Our molecule should never be inverted during reorientation.
        if det(P.rotation_matrix) < 0:
            print("Error: inverted reorientation applied.")
        return new_mol, P
    #If needed, recursively apply reorientation (due to numerical errors)
    iterations = 1
    max_iterations = 100
    new_mol, P = reorient(mol)
    while iterations < max_iterations:
        is_okay = True
        for i in range(3):
            for j in range(3):
                x = eigh(get_inertia_tensor(new_mol))[1][i][j]
                okay = True
                if i == j:
                    #Check that diagonal elements are 0 or 1
                    if (not isclose(x, 0)) and (not isclose(x, 1)):
                        okay = False
                else:
                    #Check that off-diagonal elements are 0
                    if (not isclose(x, 0)):
                        okay = False
        if okay is False:
            #If matrix is not diagonal with 1's and/or 0's, reorient
            new_mol, Q = reorient(new_mol)
            P = Q*P
            iterations += 1
        elif okay is True:
            break
    if iterations == max_iterations:
        print("Error: Could not reorient molecule after "+str(max_iterations)+" attempts")
        print(new_mol)
        print(get_inertia_tensor(new_mol))
        return False
    return new_mol, P

def get_symmetry(mol, already_oriented=False):
    '''
    Return a list of SymmOps for a molecule's point symmetry
    already_oriented: whether or not the principle axes of mol are already reoriented 
    '''
    if already_oriented == True:
        pga = PointGroupAnalyzer(mol)
    elif already_oriented == False:
        #Reorient the molecule
        oriented_mol, P = reoriented_molecule(mol)
        pga = PointGroupAnalyzer(oriented_mol)
    pg = pga.get_pointgroup()
    symm_m = []
    for op in pg:
        symm_m.append(op)
    #Handle linear molecules
    if '*' in pga.sch_symbol:
        #Add 12-fold  and reflections in place of ininitesimal rotation
        for axis in [[1,0,0],[0,1,0],[0,0,1]]:
            op = SymmOp.from_rotation_and_translation(aa2matrix(axis, pi/6), [0,0,0])
            if pga.is_valid_op(op):
                symm_m.append(op)
                #Any molecule with infinitesimal symmetry is linear;
                #Thus, it possess mirror symmetry for any axis perpendicular
                #To the rotational axis. pymatgen does not add this symmetry
                #for all linear molecules - for example, hydrogen
                if axis == [1,0,0]:
                    symm_m.append(SymmOp.from_xyz_string('x,-y,z'))
                    symm_m.append(SymmOp.from_xyz_string('x,y,-z'))
                elif axis == [0,1,0]:
                    symm_m.append(SymmOp.from_xyz_string('-x,y,z'))
                    symm_m.append(SymmOp.from_xyz_string('x,y,-z'))
                elif axis == [0,0,1]:
                    symm_m.append(SymmOp.from_xyz_string('-x,y,z'))
                    symm_m.append(SymmOp.from_xyz_string('x,-y,z'))
                #Generate a full list of SymmOps for the molecule's pointgroup
                symm_m = generate_full_symmops(symm_m, 1e-3)
                break
    #Handle non-linear molecules
    else:
        for op in pg:
            symm_m.append(op)
    #Reorient the SymmOps into mol's original frame
    if not already_oriented:
        new = []
        for op in symm_m:
            new.append(P*op*P.inverse)
        return new
    elif already_oriented:
        return symm_m

def orientation_in_wyckoff_position(mol, sg, index, randomize=True):
    '''
    Tests if a molecule meets the symmetry requirements of a Wyckoff position.
    If it does, return the rotation matrix needed. Otherwise, returns False.

    args:
        mol: pymatgen.core.structure.Molecule object. Orientation is arbitrary
        sg: the spacegroup to check
        index: the index of the Wyckoff position within the sg to check
        randomize: whether or not to apply a random rotation consistent with
            the symmetry requirements. 
    '''
    #Analyze the symmetry of the molecule and the Wyckoff position
    symm_m = get_symmetry(mol)
    symm_w = get_wyckoff_symmetry(sg)[index][0]
    #Store OperationAnalyzer objects for each SymmOp
    opa_w = []
    for op_w in symm_w:
        opa_w.append(OperationAnalyzer(op_w))
    opa_m = []
    for op_m in symm_m:
        opa_m.append(OperationAnalyzer(op_m))
    #Check for constraints from the Wyckoff symmetry...
    #If we find ANY two constraints (SymmOps with unique axes), the molecule's
    #point group MUST contain SymmOps which can be aligned to these particular
    #constraints. However, there may be multiple compatible orientations of the
    #molecule consistent with these constraints
    constraint1 = None
    constraint2 = None
    for i, op_w in enumerate(symm_w):
        if opa_w[i].axis is not None:
            constraint1 = opa_w[i]
            for j, op_w in enumerate(symm_w):
                if opa_w[j].axis is not None:
                    dot = np.dot(opa_w[i].axis, opa_w[j].axis)
                    if (not isclose(dot, 1)) and (not isclose(dot, -1)):
                        constraint2 = opa_w[j]
                        break
            break
    #Indirectly store the angle between the constraint axes
    if (constraint1 is not None
        and constraint2 is not None):
        dot_w = np.dot(constraint1.axis, constraint2.axis)
    #Store corresponding symmetry elements of the molecule
    constraints_m = []
    if constraint1 is not None:
        for i, op_m in enumerate(symm_m):
            if opa_m[i].is_conjugate(constraint1):
                constraints_m.append([opa_m[i], []])
                if constraint2 is not None:
                    for j, op_m in enumerate(symm_m):
                        if opa_m[j].is_conjugate(constraint2):
                            dot_m = np.dot(opa_m[i].axis, opa_m[j].axis)
                            #Ensure that the angles are equal
                            if isclose(dot_m, dot_w):
                                constraints_m[-1][1].append(opa_m[j])
    #Generate orientations consistent with the possible constraints
    orientations = []
    #Loop over molecular constraint sets
    for c1 in constraints_m:
        v1 = c1[0].axis
        v2 = constraint1.axis
        T = rotate_vector(v1, v2)
        #Loop over second molecular constraints
        if c1[1] == []:
            if randomize is True:
                angle = rand()*2*pi
                R = aa2matrix(v1, angle)
                T2 = np.dot(T, R)
                orientations.append(T2)
            else:
                orientations.append(T)
        for opa in c1[1]:
            v1 = np.dot(T, opa.axis)
            v2 = constraint2.axis
            R = rotate_vector(v1, v2)
            T2 = np.dot(T, R)
            orientations.append(T2)
    if constraints_m == []:
        if randomize is True:
            R = aa2matrix(1,1,random=True)
            orientations.append(R)
        else:
            orientation.append(np.identity())
    #Check each of the found orientations for consistency with the Wyckoff pos.
    #If consistent, put into an array of valid orientations
    allowed = []
    for o in orientations:
        o = SymmOp.from_rotation_and_translation(o,[0,0,0])
        mo = deepcopy(mol)
        mo.apply_operation(o)
        pga = PointGroupAnalyzer(mo)
        valid = True
        for op in symm_w:
            op_m = SymmOp.from_rotation_and_translation(op.rotation_matrix,[0,0,0])
            if not pga.is_valid_op(op_m):
                valid = False
        if valid:
            allowed.append(o)
    #Return the array of allowed orientations. If there are none, return False
    if allowed == []:
        return False
    else:
        return allowed

#Test Functionality
if __name__ == "__main__":
#---------------------------------------------------
    #Test cases: water, methane, and c60 via pymatgen
    h20 = Molecule.from_file('xyz/water.xyz')
    pga_h20 = PointGroupAnalyzer(h20)
    pg_h20 = pga_h20.get_pointgroup()
    c60 = Molecule.from_file('xyz/c60.xyz')
    pga_c60 = PointGroupAnalyzer(c60)
    pg_c60 = pga_c60.get_pointgroup()
    h2 = Molecule.from_file('xyz/hydrogen.xyz')
    pga_h2 = PointGroupAnalyzer(h2)
    pg_h2 = pga_h2.get_pointgroup()
    ch4 = Molecule.from_file('xyz/methane.xyz')
    pga_ch4 = PointGroupAnalyzer(ch4)
    pg_ch4 = pga_ch4.get_pointgroup()
    rand_mol = Molecule.from_file('xyz/random.xyz')
    pga_rand_mol = PointGroupAnalyzer(rand_mol)
    pg_rand_mol = pga_rand_mol.get_pointgroup()

    #To use:
    #orientation_in_wyckoff_position(mol, sg, WP's index in sg)
    #returns a list of orientations consistent with the WP.
    #We can choose any of these orientations at random using np.random.choice
    #To use an orientation, do mol.apply_operation(orientation)
    #More testing needs to be done to confirm the results.
    print(orientation_in_wyckoff_position(ch4, 221, 2, randomize=False))
