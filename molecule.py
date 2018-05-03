from pymatgen.core.structure import Molecule
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from ase.build import molecule
import numpy as np

def get_inertia_tensor(mol):
    #Calculate the inertia tensor for an ASE molecule
    com = mol.get_center_of_mass()
    positions = mol.get_positions()
    positions -= com  # translate center of mass to origin
    masses = mol.get_masses()

    # Initialize elements of the inertial tensor
    I11 = I22 = I33 = I12 = I13 = I23 = 0.0
    for i in range(len(mol)):
        x, y, z = positions[i]
        m = masses[i]

        I11 += m * (y ** 2 + z ** 2)
        I22 += m * (x ** 2 + z ** 2)
        I33 += m * (x ** 2 + y ** 2)
        I12 += -m * x * y
        I13 += -m * x * z
        I23 += -m * y * z

    return np.array([[I11, I12, I13],
                  [I12, I22, I23],
                  [I13, I23, I33]])
#Test Functionality
#---------------------------------------------------
from structure import *
#Test cases: water, methane, and c60 via pymatgen
'''h20 = Molecule.from_file('xyz/water.xyz')
ch4 = Molecule.from_file('xyz/methane.xyz')
c60 = Molecule.from_file('xyz/c60.xyz')
pga_h20 = PointGroupAnalyzer(h20)
pga_ch4 = PointGroupAnalyzer(ch4)
pga_c60 = PointGroupAnalyzer(c60)
pg_h20 = pga_h20.get_pointgroup()
pg_ch4 = pga_ch4.get_pointgroup()
pg_c60 = pga_c60.get_pointgroup()'''

#Test cases: water and methane via ASE
h2 = molecule('H2')
water = molecule('H2O')
methane = molecule('CH4')
c60 = molecule('C60')
newpos = [
[0., 0., 0.], #C
[0.5288, 0.161,  0.9359], #H
[ 0.2051,  0.824,  -0.6786], #H
[ 0.3345, -0.9314, -0.4496], #H
[-1.0685, -0.0537,  0.1921]] #H
methane.set_positions(newpos, apply_constraint=False)
print(get_inertia_tensor(methane))
