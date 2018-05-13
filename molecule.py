from pymatgen.core.structure import Molecule
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
#from ase.build import molecule
from matrix import *
import numpy as np
from numpy import isclose
from numpy.linalg import eigh
from numpy.linalg import det
from copy import deepcopy
from math import fabs
from random import random

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

def orientation_in_wyckoff_position(mol, wp, randomize=True):
    '''
    Tests if a molecule meets the symmetry requirements of a Wyckoff position.
    If it does, return the rotation matrix needed. Otherwise, returns False.

    args:
        mol: pymatgen.core.structure.Molecule object
        wp: a list of SymmOp objects corresponding to a Wyckoff position.
            Can be obtained from get_wyckoffs(sg, organized=False)[i]
        randomize: whether or not to apply a random rotation consistent with
            the symmetry requirements.
    '''
    

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

    '''for mol in [h20, c60, h2, ch4, rand_mol]:
        print((mol.formula))
        mol, P = reoriented_molecule(mol)
        print("Inertia tensor:")
        print(get_inertia_tensor(mol))
        print("Inertia tensor Eigenvalues:")
        print(eigh(get_inertia_tensor(mol))[0])
        print("Inertia tensor Eigenvectors:")
        print(eigh(get_inertia_tensor(mol))[1])'''

    mol = deepcopy(ch4)

    print('--------Initial Molecule--------')
    print(mol)
    print('Inertia tensor:')
    print(get_inertia_tensor(mol))
    print('Eigenvectors of inertia tensor:')
    print(eigh(get_inertia_tensor(mol))[1])

    v = [random(), random(), random()]
    a = random()*math.pi*2
    rot = aa2matrix(v, a)
    mol.apply_operation(SymmOp.from_rotation_and_translation(rot, [0,0,0]))

    print('--------After random rotation--------')
    print(mol)
    print('Inertia tensor:')
    print(get_inertia_tensor(mol))
    print('Eigenvectors of inertia tensor:')
    print(eigh(get_inertia_tensor(mol))[1])

    mol, P = reoriented_molecule(mol)

    print('--------After reorientation--------')
    print(mol)
    print('Inertia tensor:')
    print(get_inertia_tensor(mol))
    print('Eigenvectors of inertia tensor:')
    print(eigh(get_inertia_tensor(mol))[1])

    pga = PointGroupAnalyzer(mol)
    pg = pga_ch4.get_pointgroup()

    '''print('-----------Symmetry Operations:-----------')
    for op in pg:
        print(op)'''
