from pymatgen.core.structure import Molecule
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
#from ase.build import molecule
from matrix import *
import numpy as np
from numpy.linalg import eigh
from numpy.linalg import det
from copy import deepcopy
from math import fabs

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
    for x in [I11,I12,I13,I22,I23,I33]:
        x = np.round(x, decimals=4)
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
        P = eigh(A)[1]
        #reorient the molecule
        P = SymmOp.from_rotation_and_translation(P,[0,0,0])
        new_mol.apply_operation(P.inverse)
        return new_mol, P
    #If needed, recursively apply reorientation (due to numerical errors)
    iterations = 1
    max_iterations = 10
    new_mol, P = reorient(mol)
    while (not np.allclose(eigh(get_inertia_tensor(new_mol))[1], identity)
            and iterations < 10):
        new_mol, Q = reorient(new_mol)
        P = Q*P
        iterations += 1
    if iterations == max_iterations:
        print("Error: Could not reorient molecule after 10 attempts")
        return False
    return new_mol, P

#Test Functionality
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

mol = deepcopy(h20)

print('--------Initial Molecule--------')
print(mol)
print('Eigenvectors of inertia tensor:')
print(eigh(get_inertia_tensor(mol))[1])

mol, P = reoriented_molecule(mol)
print('--------After reorientation--------')
print(mol)
print('Eigenvectors of inertia tensor:')
print(eigh(get_inertia_tensor(mol))[1])

mol.apply_operation(SymmOp.from_rotation_and_translation(aa2matrix(),[0,0,0]))
print('--------After Random Rotation--------')
print(mol)
print('Eigenvectors of inertia tensor:')
print(eigh(get_inertia_tensor(mol))[1])

mol, P = reoriented_molecule(mol)
print('--------After reorientation--------')
print(mol)
print('Eigenvectors of inertia tensor:')
print(eigh(get_inertia_tensor(mol))[1])
