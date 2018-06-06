'''
Program for generation of random crystal structures.
by Scott Fredericks and Qiang Zhu, Spring 2018
Args:
sg: space group number between 1 and 230,
specie: type of atoms
N : number of atoms in the primitive cell,

output:
a structure class

possibly output cif fileS
cif file with conventional setting
'''
import sys
from spglib import get_symmetry_dataset
from pymatgen.symmetry.groups import sg_symbol_from_int_number
from pymatgen.symmetry.analyzer import generate_full_symmops
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter

from optparse import OptionParser
from scipy.spatial.distance import cdist
import numpy as np
from random import uniform as rand
from random import choice as choose
from random import randint
from math import sqrt, pi, sin, cos, acos, fabs
from database.element import Element
from copy import deepcopy
from pandas import read_csv

import database.hall as hall
from matrix import OperationAnalyzer

#some optional libs
#from vasp import read_vasp
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#from os.path import isfile

#Define variables
#------------------------------
tol_m = 1.0 #seperation tolerance in Angstroms
max1 = 30 #Attempts for generating lattices
max2 = 30 #Attempts for a given lattice
max3 = 30 #Attempts for a given Wyckoff position
minvec = 2.0 #minimum vector length
ang_min = 30
ang_max = 150
Euclidean_lattice = np.array([[1,0,0],[0,1,0],[0,0,1]])
wyckoff_df = read_csv("database/wyckoff_list.csv")
wyckoff_symmetry_df = read_csv("database/wyckoff_symmetry.csv")
wyckoff_generators_df = read_csv("database/wyckoff_generators.csv")
letters = "abcdefghijklmnopqrstuvwxyzA"
#Define functions
#------------------------------

def gaussian(min, max, sigma=3.0):
    '''
    Choose a random number from a Gaussian probability distribution centered
    between min and max. sigma is the number of standard deviations that min
    and max are away from the center. Thus, sigma is also the largest possible
    number of standard deviations corresponding to the returned value. sigma=2
    corresponds to a 95.45% probability of choosing a number between min and max
    '''
    center = (max+min)*0.5
    delta = fabs(max-min)*0.5
    ratio = delta/sigma
    while True:
        x = np.random.normal(scale=ratio, loc=center)
        if x > min and x < max:
            return x
            
def random_matrix(width=1.0, unitary=False):
    '''
    Generate a random matrix with Gaussian elements. If unitary is True,
    normalize to determinant 1
    '''
    mat = np.zeros([3,3])
    determinant = 0
    while determinant == 0:
        for x in range(3):
            for y in range(3):
                mat[x][y] = np.random.normal(scale=width)
        determinant = np.linalg.det(mat)
    if unitary:
        new = mat / np.cbrt(np.linalg.det(mat))
        return new
    else: return mat

def random_shear_matrix(width=1.0, unitary=False):
    '''
    Generate a random symmetric shear matrix with Gaussian elements. If unitary
    is True, normalize to determinant 1
    '''
    mat = np.zeros([3,3])
    determinant = 0
    while determinant == 0:
        a, b, c = np.random.normal(scale=width), np.random.normal(scale=width), np.random.normal(scale=width)
        mat = np.array([[1,a,b],[a,1,c],[b,c,1]])
        determinant = np.linalg.det(mat)
    if unitary:
        new = mat / np.cbrt(np.linalg.det(mat))
        return new
    else: return mat

def random_vector(minvec=[0.,0.,0.], maxvec=[1.,1.,1.], width=0.35, unit=False):
    '''
    Generate a random vector for lattice constant generation. The ratios between
    x, y, and z of the returned vector correspond to the ratios between a, b,
    and c. Results in a Gaussian distribution of the natural log of the ratios.
    '''
    vec = np.array([np.exp(np.random.normal(scale=width)), np.exp(np.random.normal(scale=width)), np.exp(np.random.normal(scale=width))])
    if unit:
        return vec/np.linalg.norm(vec)
    else:
        return vec

def ss_string_from_ops(ops, sg, complete=False):
    '''
    Print the Hermann-Mauguin symbol for a site symmetry group, using a list of
    SymmOps as input. Note that the symbol does not necessarily refer to the
    x,y,z axes. For information on reading these symbols, see:
    http://en.wikipedia.org/wiki/Hermann-Mauguin_notation#Point_groups
    args:
        ops: a list of SymmOp objects representing the site symmetry
        sg: International number of the spacegroup. Used to determine which
            axes to show. For example, a 3-fold rotation in a cubic system is
            written as ".3.", whereas a 3-fold rotation in a trigonal system
            is written as "3.."
        complete: whether or not all symmetry operations in the group
            are present. If False, we generate the rest
    '''
    #Return the symbol for a single axis
    #Will be called later in the function
    def get_symbol(opas, order, has_reflection):
        #ops: a list of Symmetry operations about the axis
        #order: highest order of any symmetry operation about the axis
        #has_reflection: whether or not the axis has mirror symmetry
        if has_reflection is True:
            #rotations have priority
            for opa in opas:
                if opa.order == order and opa.type == "rotation":
                    return str(opa.rotation_order)+"/m"
            for opa in opas:
                if (opa.order == order and opa.type == "rotoinversion"
                    and opa.order != 2):
                    return "-"+str(opa.rotation_order)
            return "m"
        elif has_reflection is False:
            #rotoinversion has priority
            for opa in opas:
                if opa.order == order and opa.type == "rotoinversion":
                    return "-"+str(opa.rotation_order)
            for opa in opas:
                if opa.order == order and opa.type == "rotation":
                    return str(opa.rotation_order)
            return "."
    #Given a list of single-axis symbols, return the one with highest symmetry
    #Will be called later in the function
    def get_highest_symbol(symbols):
        symbol_list = ['.','2','m','-2','2/m','3','4','-4','4/m','-3','6','-6','6/m']
        max_index = 0
        for symbol in symbols:
            i = symbol_list.index(symbol)
            if i > max_index:
                max_index = i
        return symbol_list[max_index]
    if complete is False:
        ops = generate_full_symmops(ops, 1e-3)
    #Get OperationAnalyzer object for all ops
    opas = []
    for op in ops:
        opas.append(OperationAnalyzer(op))
    #Store the symmetry of each axis
    params = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
    has_inversion = False
    #Store possible symmetry axes for crystallographic point groups
    axes = [[1,0,0],[0,1,0],[0,0,1],
            [1,1,0],[0,1,1],[1,0,1],[1,-1,0],[0,1,-1],[1,0,-1],
            [1,1,1],[-1,1,1],[1,-1,1],[-1,-1,1]]
    for i, axis in enumerate(axes):
        axes[i] = axis/np.linalg.norm(axis)
    for opa in opas:
        if opa.type != "identity" and opa.type != "inversion":
            for i, axis in enumerate(axes):
                if np.isclose(abs(np.dot(opa.axis, axis)), 1):
                    params[i].append(opa)
        elif opa.type == "inversion":
            has_inversion = True
    #Determine how many high-symmetry axes are present
    n_axes = 0
    #Store the order of each axis
    orders = []
    #Store whether or not each axis has reflection symmetry
    reflections = []
    for axis in params:
        order = 1
        high_symm = False
        has_reflection = False
        for opa in axis:
            if opa.order >= 3:
                high_symm = True
            if opa.order > order:
                order = opa.order
            if opa.order == 2 and opa.type == "rotoinversion":
                has_reflection = True
        orders.append(order)
        if high_symm == True:
            n_axes += 1
        reflections.append(has_reflection)
    #Triclinic, monoclinic, orthorhombic
    #Positions in symbol refer to x,y,z axes respectively
    if sg >= 1 and sg <= 74:
        symbol = (get_symbol(params[0], orders[0], reflections[0])+
                get_symbol(params[1], orders[1], reflections[1])+
                get_symbol(params[2], orders[2], reflections[2]))
        if symbol != "...":
            return symbol
        elif symbol == "...":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    #Trigonal, Hexagonal, Tetragonal
    elif sg >= 75 and sg <= 194:
        #1st symbol: z axis
        s1 = get_symbol(params[2], orders[2], reflections[2])
        #2nd symbol: x or y axes (whichever have higher symmetry)
        s2x = get_symbol(params[0], orders[0], reflections[0])
        s2y = get_symbol(params[1], orders[1], reflections[1])
        s2 = get_highest_symbol([s2x, s2y])
        #3rd symbol: face-diagonal axes (whichever have highest symmetry)
        s3a = get_symbol(params[3], orders[3], reflections[3])
        s3b = get_symbol(params[4], orders[4], reflections[4])
        s3c = get_symbol(params[5], orders[5], reflections[5])
        s3d = get_symbol(params[6], orders[6], reflections[6])
        s3e = get_symbol(params[7], orders[7], reflections[7])
        s3f = get_symbol(params[8], orders[8], reflections[8])
        s3 = get_highest_symbol([s3a, s3b, s3c, s3d, s3e, s3f])
        symbol = s1 + s2 + s3
        if symbol != "...":
            return symbol
        elif symbol == "...":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    #Cubic
    elif sg >= 195 and sg <= 230:
        pass
        #1st symbol: x, y, and/or z axes (whichever have highest symmetry)
        s1x = get_symbol(params[0], orders[0], reflections[0])
        s1y = get_symbol(params[1], orders[1], reflections[1])
        s1z = get_symbol(params[2], orders[2], reflections[2])
        s1 = get_highest_symbol([s1x, s1y, s1z])
        #2nd symbol: body-diagonal axes (whichever has highest symmetry)
        s2a = get_symbol(params[9], orders[9], reflections[9])
        s2b = get_symbol(params[10], orders[10], reflections[10])
        s2c = get_symbol(params[11], orders[11], reflections[11])
        s2d = get_symbol(params[12], orders[12], reflections[12])
        s2 = get_highest_symbol([s2a, s2b, s2c, s2d])
        #3rd symbol: face-diagonal axes (whichever have highest symmetry)
        s3a = get_symbol(params[3], orders[3], reflections[3])
        s3b = get_symbol(params[4], orders[4], reflections[4])
        s3c = get_symbol(params[5], orders[5], reflections[5])
        s3d = get_symbol(params[6], orders[6], reflections[6])
        s3e = get_symbol(params[7], orders[7], reflections[7])
        s3f = get_symbol(params[8], orders[8], reflections[8])
        s3 = get_highest_symbol([s3a, s3b, s3c, s3d, s3e, s3f])
        symbol = s1 + s2 + s3
        if symbol != "...":
            return symbol
        elif symbol == "...":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    else:
        print("Error: invalid spacegroup number")
        return

def are_equal(op1, op2, allow_pbc=True, rtol=1e-3, atol=1e-3):
    #Check two SymmOps for equivalence
    #pbc=True means integer translations will be ignored
    m1 = op1.rotation_matrix
    m2 = op2.rotation_matrix
    #Check that rotations are equivalent
    if not np.allclose(m1, m2, rtol=rtol, atol=atol):
        return False
    v1 = op1.translation_vector
    v2 = op2.translation_vector
    if allow_pbc is False:
        #Check if translation vectors are equal
        if np.allclose(v1, v2, rtol=rtol, atol=atol):
            return True
        else: return False
    elif allow_pbc is True:
        #Check if translation vectors are equal up to integer difference
        difference = v1 - v2
        if np.allclose(difference, np.round(difference), rtol=rtol, atol=atol):
            return True
        else: return False

def create_matrix():
    matrix = []
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                matrix.append([i,j,k])
    return np.array(matrix, dtype=float)

#Euclidean distance
def distance(xyz, lattice): 
    xyz = xyz - np.round(xyz)
    matrix = create_matrix()
    matrix += xyz
    matrix = np.dot(matrix, lattice)
    return np.min(cdist(matrix,[[0,0,0]]))       

def check_distance(coord1, coord2, specie1, specie2, lattice):
    """
    check the distances between two set of atoms
    Args:
    coord1: multiple list of atoms e.g. [[0,0,0],[1,1,1]]
    specie1: the corresponding type of coord1, e.g. ['Na','Cl']
    coord2: a list of new atoms: [0.5, 0.5 0.5]
    specie2: the type of coord2: 'Cl'
    lattice: cell matrix
    """
    #add PBC
    coord2s = []
    matrix = create_matrix()
    for coord in coord2:
        for m in matrix:
            coord2s.append(coord+m)
    coord2 = np.array(coord2s)

    coord2 = np.dot(coord2, lattice)
    if len(coord1)>0:
        for coord, element in zip(coord1, specie1):
            coord = np.dot(coord, lattice)
            d_min = np.min(cdist(coord, coord2))
            tol = 0.5*(Element(element).covalent_radius + Element(specie2).covalent_radius)
            #print(d_min, tol)
            if d_min < tol:
                return False
        return True
    else:
        return True

def get_center(xyzs, lattice):
    """
    to find the geometry center of the clusters under PBC
    """
    matrix0 = create_matrix()
    xyzs -= np.round(xyzs)
    for atom1 in range(1,len(xyzs)):
        dist_min = 10.0
        for atom2 in range(0, atom1):
            #shift atom1 to position close to atom2
            #print(np.round(xyzs[atom1] - xyzs[atom2]))
            #xyzs[atom2] += np.round(xyzs[atom1] - xyzs[atom2])
            #print(xyzs[atom1] - xyzs[atom2])
            matrix = matrix0 + (xyzs[atom1] - xyzs[atom2])
            #print(matrix)
            matrix = np.dot(matrix, lattice)
            dists = cdist(matrix, [[0,0,0]])
            if np.min(dists) < dist_min:
                dist_min = np.min(dists)
                #print(dist_min, matrix[np.argmin(dists)]/4.72)
                matrix_min = matrix0[np.argmin(dists)]
        #print(atom1, xyzs[atom1], matrix_min, dist_min)
        xyzs[atom1] += matrix_min
    #print(xyzs)
    return xyzs.mean(0)

def para2matrix(cell_para, radians=True, format='lower'):
    """ 1x6 (a, b, c, alpha, beta, gamma) -> 3x3 representation -> """
    a = cell_para[0]
    b = cell_para[1]
    c = cell_para[2]
    alpha = cell_para[3]
    beta = cell_para[4]
    gamma = cell_para[5]
    if radians is not True:
        rad = pi/180.
        alpha *= rad
        beta *= rad
        gamma *= rad
    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)
    sin_alpha = np.sin(alpha)
    matrix = np.zeros([3,3])
    if format == 'lower':
        #Generate a lower-diagonal matrix
        c1 = c*cos_beta
        c2 = (c*(cos_alpha - (cos_beta * cos_gamma))) / sin_gamma
        matrix[0][0] = a
        matrix[1][0] = b * cos_gamma
        matrix[1][1] = b * sin_gamma
        matrix[2][0] = c1
        matrix[2][1] = c2
        matrix[2][2] = sqrt(c**2 - c1**2 - c2**2)
    elif format == 'symmetric':
        #TODO: allow generation of symmetric matrices
        pass
    elif format == 'upper':
        #Generate an upper-diagonal matrix
        a3 = a*cos_beta
        a2 = (a*(cos_gamma - (cos_beta * cos_alpha))) / sin_alpha
        matrix[2][2] = c
        matrix[1][2] = b * cos_alpha
        matrix[1][1] = b * sin_alpha
        matrix[0][2] = a3
        matrix[0][1] = a2
        matrix[0][0] = sqrt(a**2 - a3**2 - a2**2)
        pass
    return matrix

def matrix2para(matrix, radians=True):
    """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
    cell_para = np.zeros(6)
    #a
    cell_para[0] = np.linalg.norm(matrix[0])
    #b
    cell_para[1] = np.linalg.norm(matrix[1])
    #c
    cell_para[2] = np.linalg.norm(matrix[2])

    #alpha
    cell_para[3] = angle(matrix[1], matrix[2])
    #beta
    cell_para[4] = angle(matrix[0], matrix[2])
    #gamma
    cell_para[5] = angle(matrix[0], matrix[1])
    
    if not radians:
        #convert radians to degrees
        deg = 180./pi
        cell_para[3] *= deg
        cell_para[4] *= deg
        cell_para[5] *= deg
    return cell_para

def cellsize(sg):
    """
    Returns the number of duplications in the conventional lattice
    """
    symbol = sg_symbol_from_int_number(sg)
    letter = symbol[0]
    if letter == 'P':
    	return 1
    if letter in ['A', 'C', 'I']:
    	return 2
    elif letter in ['R']:
    	return 3
    elif letter in ['F']:
    	return 4
    else: return "Error: Could not determine lattice type"

def find_short_dist(coor, lattice, tol):
    """
    here we find the atomic pairs with shortest distance
    and then build the connectivity map
    """
    pairs=[]
    graph=[]
    for i in range(len(coor)):
        graph.append([])

    for i1 in range(len(coor)-1):
        for i2 in range(i1+1,len(coor)):
            dist = distance(coor[i1]-coor[i2], lattice)
            if dist <= tol:
                #dists.append(dist)
                pairs.append([i1,i2,dist])
    pairs = np.array(pairs)
    if len(pairs) > 0:
        #print('--------', dists <= (min(dists) + 0.1))
        d_min = min(pairs[:,-1]) + 1e-3
        sequence = [pairs[:,-1] <= d_min]
        #print(sequence)
        pairs = pairs[sequence]
        #print(pairs)
        #print(len(coor))
        for pair in pairs:
            pair0=int(pair[0])
            pair1=int(pair[1])
            #print(pair0, pair1, len(graph))
            graph[pair0].append(pair1)
            graph[pair1].append(pair0)

    return pairs, graph

def connected_components(graph):
    '''Given an undirected graph (a 2d array of indices), return a set of
    connected components, each connected component being an (arbitrarily
    ordered) array of indices which are connected either directly or indirectly.
    '''
    def add_neighbors(el, seen=[]):
        '''
        Find all elements which are connected to el. Return an array which
        includes these elements and el itself.
        '''
        #seen stores already-visited indices
        if seen == []: seen = [el]
        #iterate through the neighbors (x) of el
        for x in graph[el]:
            if x not in seen:
                seen.append(x)
                #Recursively find neighbors of x
                add_neighbors(x, seen)
        return seen

    #Create a list of indices to iterate through
    unseen = list(range(len(graph)))
    sets = []
    i = 0
    while (unseen != []):
        #x is the index we are finding the connected component of
        x = unseen.pop()
        sets.append([])
        #Add neighbors of x to the current connected component
        for y in add_neighbors(x):
            sets[i].append(y)
            #Remove indices which have already been found
            if y in unseen: unseen.remove(y)
        i += 1
    return sets

def merge_coordinate(coor, lattice, wyckoff, sg, tol):
    while True:
        pairs, graph = find_short_dist(coor, lattice, tol)
        if len(pairs)>0:
            if len(coor) > len(wyckoff[-1][0]):
                merged = []
                groups = connected_components(graph)
                for group in groups:
                    #print(coor[group])
                    #print(coor[group].mean(0))
                    merged.append(get_center(coor[group], lattice))
                merged = np.array(merged)
                #if check_wyckoff_position(merged, sg, wyckoff) is not False:
                if check_wyckoff_position(merged, sg, exact_translation=False) is False:
                    '''print('something is wrong')
                    print(coor)
                    print(merged)
                    print(sg)'''
                    #exit(0)
                    return coor, False
                else:
                    coor = merged

            else:#no way to merge
                #print('no way to Merge, FFFFFFFFFFFFFFFFFFFFFFF----------------')
                return coor, False
        else: 
            return coor, True

def estimate_volume(numIons, species, factor=2.0):
    volume = 0
    for numIon, specie in zip(numIons, species):
        volume += numIon*4/3*pi*Element(specie).covalent_radius**3
    return factor*volume

def generate_lattice(sg, volume, minvec=tol_m, minangle=pi/6, max_ratio=10.0, maxattempts = 100):
    """
    generate the lattice according to the space group symmetry and number of atoms
    if the space group has centering, we will transform to conventional cell setting
    If the generated lattice does not meet the minimum angle and vector requirements,
    we try to generate a new one, up to maxattempts times

    args:
        sg: International number of the space group
        volume: volume of the lattice
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
    """
    maxangle = pi-minangle
    for n in range(maxattempts):
        #Triclinic
        if sg <= 2:
            #Derive lattice constants from a random matrix
            mat = random_shear_matrix(width=0.2)
            a, b, c, alpha, beta, gamma = matrix2para(mat)
            x = sqrt(1-cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*(cos(alpha)*cos(beta)*cos(gamma)))
            vec = random_vector()
            abc = volume/x
            xyz = vec[0]*vec[1]*vec[2]
            a = vec[0]*np.cbrt(abc)/np.cbrt(xyz)
            b = vec[1]*np.cbrt(abc)/np.cbrt(xyz)
            c = vec[2]*np.cbrt(abc)/np.cbrt(xyz)
        #Monoclinic
        elif sg <= 15:
            alpha, gamma  = pi/2, pi/2
            beta = gaussian(minangle, maxangle)
            x = sin(beta)
            vec = random_vector()
            xyz = vec[0]*vec[1]*vec[2]
            abc = volume/x
            a = vec[0]*np.cbrt(abc)/np.cbrt(xyz)
            b = vec[1]*np.cbrt(abc)/np.cbrt(xyz)
            c = vec[2]*np.cbrt(abc)/np.cbrt(xyz)
        #Orthorhombic
        elif sg <= 74:
            alpha, beta, gamma = pi/2, pi/2, pi/2
            x = 1
            vec = random_vector()
            xyz = vec[0]*vec[1]*vec[2]
            abc = volume/x
            a = vec[0]*np.cbrt(abc)/np.cbrt(xyz)
            b = vec[1]*np.cbrt(abc)/np.cbrt(xyz)
            c = vec[2]*np.cbrt(abc)/np.cbrt(xyz)
        #Tetragonal
        elif sg <= 142:
            alpha, beta, gamma = pi/2, pi/2, pi/2
            x = 1
            vec = random_vector()
            c = vec[2]/(vec[0]*vec[1])*np.cbrt(volume/x)
            a = b = sqrt((volume/x)/c)
        #Trigonal/Rhombohedral/Hexagonal
        elif sg <= 194:
            alpha, beta, gamma = pi/2, pi/2, pi/3*2
            x = sqrt(3.)/2.
            vec = random_vector()
            c = vec[2]/(vec[0]*vec[1])*np.cbrt(volume/x)
            a = b = sqrt((volume/x)/c)
        #Cubic
        else:
            alpha, beta, gamma = pi/2, pi/2, pi/2
            s = (volume) ** (1./3.)
            a, b, c = s, s, s
        #Check that lattice meets requirements
        maxvec = (a*b*c)/(minvec**2)
        if(a>minvec and b>minvec and c>minvec
        and a<maxvec and b<maxvec and c<maxvec
        and alpha>minangle and beta>minangle and gamma>minangle
        and alpha<maxangle and beta<maxangle and gamma<maxangle
        and a/b<max_ratio and a/c<max_ratio and b/c<max_ratio
        and b/a<max_ratio and c/a<max_ratio and c/b<max_ratio):
            return np.array([a, b, c, alpha, beta, gamma])
    #If maxattempts tries have been made without success
    print("Error: Could not generate lattice after "+str(n+1)+" attempts")
    return

def filter_site(v): #needs to explain
    #Adjusts coordinates to be greater than 0 and less than 1
    w = v
    for i in range(len(w)):
    	while w[i]<0: w[i] += 1
    	while w[i]>=1: w[i] -= 1
    return w

def choose_wyckoff(wyckoffs, number):
    """
    choose the wyckoff sites based on the current number of atoms
    rules 
    1, the newly added sites is equal/less than the required number.
    2, prefer the sites with large multiplicity
    """
    if rand(0,1)>0.5: #choose from high to low
        for wyckoff in wyckoffs:
            if len(wyckoff[0]) <= number:
                return choose(wyckoff)
        return False
    else:
        good_wyckoff = []
        for wyckoff in wyckoffs:
            if len(wyckoff[0]) <= number:
                for w in wyckoff:
                    good_wyckoff.append(w)
        if len(good_wyckoff) > 0:
            return choose(good_wyckoff)
        else:
            return False


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

def get_wyckoff_symmetry(sg, molecular=False):
    '''
    Returns a list of Wyckoff position site symmetry for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a point within the WP
    3rd index: a site symmetry SymmOp of the point
    molecular: whether or not to return the Euclidean point symmetry operations
        If True, cuts off translational part of operation, and converts non-orthogonal
        (3-fold and 6-fold rotation) operations to pure rotations
    '''
    P = SymmOp.from_rotation_and_translation([[1,-.5,0],[0,sqrt(3)/2,0],[0,0,1]], [0,0,0])
    symmetry_strings = eval(wyckoff_symmetry_df["0"][sg])
    symmetry = []
    convert = False
    if molecular is True:
        if sg >= 143 and sg <= 194:
            convert = True
    #Loop over Wyckoff positions
    for x in symmetry_strings:
        symmetry.append([])
        #Loop over points in WP
        for y in x:
            symmetry[-1].append([])
            #Loop over ops
            for z in y:
                op = SymmOp.from_xyz_string(z)
                if convert is True:
                    #Convert non-orthogonal trigonal/hexagonal operations
                    op = P*op*P.inverse
                if molecular is False:
                    symmetry[-1][-1].append(op)
                elif molecular is True:
                    op = SymmOp.from_rotation_and_translation(op.rotation_matrix,[0,0,0])
                    symmetry[-1][-1].append(op)
    return symmetry

def get_wyckoff_generators(sg):
    '''
    Returns a list of Wyckoff generators for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a generator for the WP
    '''
    generators_strings = eval(wyckoff_generators_df["0"][sg])
    generators = []
    #Loop over Wyckoff positions
    for wp in generators_strings:
        generators.append([])
        #Loop over points in WP
        for op in wp:
            generators[-1].append(SymmOp.from_xyz_string(op))
    return generators

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
        difference = SymmOp((op*point).affine_matrix - point.affine_matrix)
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
            el = SymmOp.from_rotation_and_translation(op.rotation_matrix, op.translation_vector - np.round(displacement))
            symmetry.append(el)
    return symmetry

def check_wyckoff_position(points, sg, wyckoffs=None, exact_translation=False):
    '''
    Given a list of points, return index of Wyckoff position in space group.
    If no match found, returns False.

    Args:
        points: a list of 3d coordinates or SymmOps to check
        sg: the international space group number to check
        wyckoffs: a list of Wyckoff positions obtained from get_wyckoffs.
        exact_translation: whether we require two SymmOps to have exactly equal
            translational components. If false, translations related by +-1
            are considered equal
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
    new_points = []
    #
    if exact_translation == False:
        for p in points:
            new_points.append(p - np.floor(p))
        points = new_points
    w_symm_all = get_wyckoff_symmetry(sg)
    p_symm = []
    #If exact_translation is false, store WP's which might be a match
    possible = []
    for x in points:
        p_symm.append(site_symm(x, gen_pos))
    for i, wp in enumerate(wyckoffs):
        w_symm = w_symm_all[i]
        if len(p_symm) == len(w_symm):
            temp = deepcopy(w_symm)
            for p in p_symm:
                for w in temp:
                    if exact_translation:
                        if p == w:
                            temp.remove(w)
                    elif not exact_translation:
                        temp2 = deepcopy(w)
                        for op_p in p:
                            for op_w in w:
                                #Check that SymmOp's are equal up to some integer translation
                                if are_equal(op_w, op_p, allow_pbc=True):
                                    temp2.remove(op_w)
                        if temp2 == []:
                            temp.remove(w)
            if temp == []:
                #If we find a match with exact translations
                if exact_translation:
                    return i
                elif not exact_translation:
                    possible.append(i)
        #If no matching WP's are found
    if len(possible) == 0:
        return False
    #If exactly one matching WP is found
    elif len(possible) == 1:
        return possible
    #If multiple WP's are found
    else:
        #TODO: add a way to differentiate between possible WP's
        #print("Warning: multiple Wyckoff positions found")
        return possible

class random_crystal():
    def __init__(self, sg, species, numIons, factor):
        #numIons *= cellsize(sg) #Should not be called twice
        #volume = estimate_volume(numIons, species, factor)
        #wyckoffs = get_wyckoffs(sg, organized=True) #2D Array of Wyckoff positions organized by multiplicity
        
        #Necessary input
        numIons = np.array(numIons) #must convert it to np.array
        self.factor = factor
        self.numIons0 = numIons
        self.sg = sg
        self.species = species
        self.Msgs()
        self.numIons = numIons * cellsize(self.sg)
        self.volume = estimate_volume(self.numIons, self.species, self.factor)
        self.wyckoffs = get_wyckoffs(self.sg, organized=True) #2D Array of Wyckoff positions organized by multiplicity
        self.generate_crystal()


    def Msgs(self):
        self.Msg1 = 'Error: the number is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'

    def check_compatible(self):
        """
        check if the number of atoms is compatible with the wyckoff positions
        needs to improve later
        """
        N_site = [len(x[0]) for x in self.wyckoffs]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for numIon in self.numIons:
            #Check that the number of ions is a multiple of the smallest Wyckoff position
            if numIon % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numIon
                    for x in self.wyckoffs:
                        for wp in x:
                            removed = False
                            while remaining >= len(wp) and wp not in removed_wyckoffs:
                                #Check if WP has at least one degree of freedom
                                op = wp[0]
                                remaining -= len(wp)
                                if np.allclose(op.rotation_matrix, np.zeros([3,3])):
                                    removed_wyckoffs.append(wp)
                                    removed = True
                                else:
                                    has_freedom = True
                    if remaining != 0:
                        return False
        if has_freedom:
            return True
        else:
            #print("Warning: Wyckoff Positions have no degrees of freedom.")
            return 0

    def generate_crystal(self, max1=max1, max2=max2, max3=max3):
        """the main code to generate random crystal """
        #Check the minimum number of degrees of freedom within the Wyckoff positions
        degrees = self.check_compatible()
        if degrees is False:
            print(self.Msg1)
            self.struct = None
            self.valid = False
            return
        else:
            if degrees is 0:
                max1 = 5
                max2 = 5
                max3 = 5
            #Calculate a minimum vector length for generating a lattice
            minvector = max(max(2.0*Element(specie).covalent_radius for specie in self.species), tol_m)
            for cycle1 in range(max1):
                #1, Generate a lattice
                cell_para = generate_lattice(self.sg, self.volume, minvec=minvector)
                cell_matrix = para2matrix(cell_para)
                if abs(self.volume - np.linalg.det(cell_matrix)) > 1.0: 
                    print('Error, volume is not equal to the estimated value: ', self.volume, ' -> ', np.linalg.det(cell_matrix))
                    print('cell_para:  ', cell_para)
                    sys.exit(0)

                coordinates_total = [] #to store the added coordinates
                sites_total = []      #to store the corresponding specie
                good_structure = False

                for cycle2 in range(max2):
                    coordinates_tmp = deepcopy(coordinates_total)
                    sites_tmp = deepcopy(sites_total)
                    
            	    #Add specie by specie
                    for numIon, specie in zip(self.numIons, self.species):
                        numIon_added = 0
                        tol = max(0.5*Element(specie).covalent_radius, tol_m)

                        #Now we start to add the specie to the wyckoff position
                        for cycle3 in range(max3):
                            #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                            ops = choose_wyckoff(self.wyckoffs, numIon-numIon_added) 
                            if ops is not False:
            	    	    #Generate a list of coords from ops
                                point = np.random.random(3)
                                #print('generating new points:', point)
                                coords = np.array([op.operate(point) for op in ops])
                                #merge_coordinate if the atoms are close
                                coords_toadd, good_merge = merge_coordinate(coords, cell_matrix, self.wyckoffs, self.sg, tol)
                                if good_merge:
                                    coords_toadd -= np.floor(coords_toadd) #scale the coordinates to [0,1], very important!
                                    #print('existing: ', coordinates_tmp)
                                    if check_distance(coordinates_tmp, coords_toadd, sites_tmp, specie, cell_matrix):
                                        coordinates_tmp.append(coords_toadd)
                                        sites_tmp.append(specie)
                                        numIon_added += len(coords_toadd)
                                    if numIon_added == numIon:
                                        coordinates_total = deepcopy(coordinates_tmp)
                                        sites_total = deepcopy(sites_tmp)
                                        break
                        if numIon_added != numIon:
                            break  #need to repeat from the 1st species

                    if numIon_added == numIon:
                        #print(self.Msg6)
                        good_structure = True
                        break
                    else: #reset the coordinates and sites
                        coordinates_total = []
                        sites_total = []

                if good_structure:
                    final_coor = []
                    final_site = []
                    final_number = []
                    final_lattice = cell_matrix
                    for coor, ele in zip(coordinates_total, sites_total):
                        for x in coor:
                            final_coor.append(x)
                            final_site.append(ele)
                            final_number.append(Element(ele).z)

                    self.lattice = final_lattice                    
                    self.coordinates = np.array(final_coor)
                    self.sites = final_site                    
                    self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                    self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                    self.valid = True
                    return
        if degrees == 0: print("Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2

if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-s", "--spacegroup", dest="sg", metavar='sg', default=206, type=int,
            help="desired space group number: 1-230, e.g., 206")
    parser.add_option("-e", "--element", dest="element", default='Li', 
            help="desired elements: e.g., Li", metavar="element")
    parser.add_option("-n", "--numIons", dest="numIons", default=16, 
            help="desired numbers of atoms: 16", metavar="numIons")
    parser.add_option("-v", "--volume", dest="factor", default=2.0, type=float, 
            help="volume factors: default 2.0", metavar="factor")

    (options, args) = parser.parse_args()    
    element = options.element
    number = options.numIons
    numIons = []

    if element.find(',') > 0:
        system = element.split(',')
        for x in number.split(','):
            numIons.append(int(x))
    else:
        system = [element]
        numIons = [int(number)]
    for i in range(100):
        numIons0 = np.array(numIons)
        sg = options.sg
        rand_crystal = random_crystal(options.sg, system, numIons0, options.factor)

        if rand_crystal.valid:
            #pymatgen style
            rand_crystal.struct.to(fmt="poscar", filename = '1.vasp')

            #spglib style structure called cell
            ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)['number']
            print('Space group  requested: ', sg, 'generated', ans)

            #print(CifWriter(new_struct, symprec=0.1).__str__())
            #print('Space group:', finder.get_space_group_symbol(), 'tolerance:', tol)
            #output wyckoff sites only

        else: 
            print('something is wrong')
            #print(len(new_struct.frac_coords))
            break
            #print(new_struct)
