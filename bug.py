import numpy as np
from scipy.spatial.distance import cdist
from spglib import get_symmetry_dataset
from pymatgen.symmetry.groups import sg_symbol_from_int_number
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from pandas import read_csv
from numpy.random import random

wyckoff_df = read_csv("database/wyckoff_list.csv")
Euclidean_lattice = np.array([[1,0,0],[0,1,0],[0,0,1]])
wyckoff_symmetry_df = read_csv("database/wyckoff_symmetry.csv")
#Euclidean distance
def distance(xyz, lattice): 
    xyz = xyz - np.round(xyz)
    matrix = create_matrix()
    matrix += xyz
    matrix = np.dot(matrix, lattice)
    return np.min(cdist(matrix,[[0,0,0]]))       

def create_matrix():
    matrix = []
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                matrix.append([i,j,k])
    return np.array(matrix, dtype=float)

def get_wyckoff_symmetry(sg):
    '''
    Returns a list of Wyckoff position site symmetry for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a point within the WP
    3rd index: a site symmetry SymmOp of the point
    '''
    symmetry_strings = eval(wyckoff_symmetry_df["0"][sg])
    symmetry = []
    #Loop over Wyckoff positions
    for x in symmetry_strings:
        symmetry.append([])
        #Loop over points in WP
        for y in x:
            symmetry[-1].append([])
            #Loop over 
            for z in y:
                symmetry[-1][-1].append(SymmOp.from_xyz_string(z))
    return symmetry


def get_wyckoffs(sg, organized=False):
    '''
    Returns a list of Wyckoff positions for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a SymmOp object in the WP
    '''
    wyckoff_strings = eval(wyckoff_df["0"][sg])
    wyckoffs = []
    for x in wyckoff_strings:
        wyckoffs.append([])
        for y in x:
            wyckoffs[-1].append(SymmOp.from_xyz_string(y))
    if organized:
        wyckoffs_organized = [[]] #2D Array of WP's organized by multiplicity
        old = len(wyckoffs[0])
        for wp in wyckoffs:
            mult = len(wp)
            if mult != old:
                wyckoffs_organized.append([])
                old = mult
            wyckoffs_organized[-1].append(wp)
        return wyckoffs_organized
    else:
        return wyckoffs


def site_symm(point, gen_pos, tol=1e-3, lattice=Euclidean_lattice):
    '''
    Given gen_pos (a list of SymmOps), return the list of symmetry operations
    leaving a point (coordinate or SymmOp) invariant.
    '''
    #Convert point into a SymmOp
    if type(point) != SymmOp:
        point = SymmOp.from_rotation_and_translation([[0,0,0],[0,0,0],[0,0,0]], point)
    symmetry = []
    for op in gen_pos:
        is_symmetry = True
        #Calculate the effect of applying op to point
        difference = SymmOp(point.affine_matrix - (op*point).affine_matrix)
        #Check that the rotation matrix is unaltered by op
        if not np.allclose(difference.rotation_matrix, np.zeros((3,3)), rtol = 1e-3, atol = 1e-3):
            is_symmetry = False
        #Check that the displacement is less than tol
        displacement = difference.translation_vector
        if distance(displacement, lattice) > tol:
            is_symmetry = False
        if is_symmetry:
            '''The actual site symmetry's translation vector may vary from op by
            a factor of +1 or -1 (especially when op contains +-1/2).
            We record this to distinguish between special Wyckoff positions.
            As an example, consider the point (-x+1/2,-x,x+1/2) in position 16c
            of space group Ia-3(206). The site symmetry includes the operations
            (-z+1,x-1/2,-y+1/2) and (y+1/2,-z+1/2,-x+1). These operations are
            not listed in the general position, but correspond to the operations
            (-z,x+1/2,-y+1/2) and (y+1/2,-z+1/2,-x), respectively, just shifted
            by (+1,-1,0) and (0,0,+1), respectively.
            '''
            el = SymmOp.from_rotation_and_translation(op.rotation_matrix, op.translation_vector + np.round(displacement))
            symmetry.append(el)
    return symmetry

def check_wyckoff_position(points, sg, wyckoffs=None):
    '''
    Given a list of points, return index of Wyckoff position in space group.
    If no match found, returns False.

    Args:
        points: a list of 3d coordinates or SymmOps to check
        sg: the international space group number to check
        wyckoffs: a list of wyckoff positions obtained from get_wyckoffs.
    '''
    #TODO: Create function for assigning WP to a single point
    #QZ: I am not sure if this is really needed
    points = np.array(points)
    points = np.around((points*1e+10))/1e+10

    if wyckoffs == None:
        wyckoffs = get_wyckoffs(sg)
        gen_pos = wyckoffs[0]
    else:
        gen_pos = wyckoffs[0][0]
    w_symm_all = get_wyckoff_symmetry(sg)
    p_symm = []
    for x in points:
        p_symm.append(site_symm(x, gen_pos))

    '''#------------------debug------------------------#
    print('------ site symmetry from our code------------')
    for p in p_symm:
        print([p0.as_xyz_string() for p0 in p])
    print('------ site symmetry from database------------')
    for i, wp in enumerate(wyckoffs):
        w_symm = w_symm_all[i]
        if len(p_symm) == len(w_symm):
            print('---------------------------------')
            for w in w_symm:
                print([p0.as_xyz_string() for p0 in w])
    #------------------debug------------------------#'''

    for i, wp in enumerate(wyckoffs):
        w_symm = w_symm_all[i]
        if len(p_symm) == len(w_symm):
            temp = w_symm
            for p in p_symm:
                for w in temp:
                    if p == w:
                        temp.remove(w)
            if temp == []:
                return i
    return False


#It is stange that check_wyckoff returns false for the following set.
#should be 8j position for spg 97
'''
coor = np.array(
[[0.85540127, 0.35540127, 0.25],
 [0.14459873, 0.64459873, 0.25],
 [0.64459873, 0.85540127, 0.25],
 [0.35540127, 0.14459873, 0.25],
 [0.35540127, 0.85540127, 0.75],
 [0.64459873, 0.14459873, 0.75],
 [0.14459873, 0.35540127, 0.75],
 [0.85540127, 0.64459873, 0.75]]
)
print(check_wyckoff_position(coor, 97))

coor = np.array(
[[ 0.23631801, -0.23631801, -0.06786002],
 [-0.23631801,  0.23631801, -0.06786002],
 [ 0.23631801,  0.23631801, -0.06786002],
 [-0.23631801, -0.23631801, -0.06786002]]
)
print(check_wyckoff_position(coor, 99))
'''

#Test check_wyckoff_position() by plugging in random coordinates
#Passes using new wyckoff_list.csv and wyckoff_symmetry.csv
print("Calculating spacegroup:")
allpassed = True
for sg in range(1, 231):
    print(sg,)
    wyckoffs = get_wyckoffs(sg)
    for i, wp in enumerate(wyckoffs):
        xyz = [random(), random(), random()]
        coor = []
        for p in wp:
            coor.append(p.operate(xyz))
        passed = check_wyckoff_position(coor, sg)
        if passed == False and passed != 0:
            print("Failure for spacegroup "+str(sg)+" position # "+str(i))
            allpassed = False
if allpassed: print("All spacegroups passed.")


#This set is to check the numerical tolerance of check_wyckoff
#coor = np.array(
#    [[-2.77555756e-17, -2.77555756e-17,  1.29634884e+00],
#     [-5.00000000e-01, -5.00000000e-01,  7.96348839e-01]])
#print(check_wyckoff_position(coor, 79))
#coor = np.around(coor*1000)/1000
#print(check_wyckoff_position(coor, 79))

