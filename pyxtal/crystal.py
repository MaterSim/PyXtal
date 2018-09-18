"""
Module for generation of random atomic crystals with symmetry constraints. A
pymatgen- or spglib-type structure object is created, which can be saved to a
.cif file. Options (preceded by two dashes) are provided for command-line usage
of the module:  

    spacegroup (-s): the international spacegroup number (between 1 and 230)
        to be generated. In the case of a 2D crystal (using option '-d 2'),
        this will instead be the layer group number (between 1 and 80).
        Defaults to 36.  

    element (-e): the chemical symbol of the atom(s) to use. For multiple
        molecule types, separate entries with commas. Ex: "C", "H, O, N".
        Defaults to Li  

    numIons (-n): the number of atoms in the PRIMITIVE unit cell
        (For P-type spacegroups, this is the same as the number of molecules in
        the conventional unit cell. For A, B, C, and I-centered spacegroups,
        this is half the number of the conventional cell. For F-centered unit
        cells, this is one fourth the number of the conventional cell.).
        For multiple atom types, separate entries with commas.
        Ex: "8", "1, 4, 12". Defaults to 16  

    factor (-f): the relative volume factor used to generate the unit cell.
        Larger values result in larger cells, with atoms spaced further apart.
        If generation fails after max attempts, consider increasing this value.
        Defaults to 2.0  

    verbosity (-v): the amount of information which should be printed for each
        generated structure. For 0, only prints the requested and generated
        spacegroups. For 1, also prints the contents of the generated pymatgen
        structure. Defaults to 0  

    attempts (-a): the number of structures to generate. Note: if any of the
        attempts fail, the number of generated structures will be less than this
        value. Structures will be output to separate cif files. Defaults to 1  

    outdir (-o): the file directory where cif files will be output to.
        Defaults to "out"

    dimension (-d): 3 for 3D, or 2 for 2D, 1 for 1D. If 2D, generates a 2D
        crystal using a layer group number instead of a space group number. For
        1D, we use a Rod group number.  

    thickness (-t): The thickness, in Angstroms, to use when generating a
        2D crystal. Note that this will not necessarily be one of the lattice
        vectors, but will represent the perpendicular distance along the non-
        periodic direction. For 1D crystals, we use this value
        as the cross-sectional area of the crystal. Defaults to None  
"""

import sys
from time import time
#from pkg_resources import resource_string
from pkg_resources import resource_filename
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
from copy import deepcopy
from pandas import read_csv

from pyxtal.database.element import Element
import pyxtal.database.hall as hall
from pyxtal.database.layergroup import Layergroup
from pyxtal.operations import OperationAnalyzer
from pyxtal.operations import angle
from pyxtal.operations import random_vector
from pyxtal.operations import are_equal
from pyxtal.operations import random_shear_matrix


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
#Matrix for a Euclidean metric
Euclidean_lattice = np.array([[1,0,0],[0,1,0],[0,0,1]])

wyckoff_df = read_csv(resource_filename("pyxtal", "database/wyckoff_list.csv"))
wyckoff_symmetry_df = read_csv(resource_filename("pyxtal", "database/wyckoff_symmetry.csv"))
wyckoff_generators_df = read_csv(resource_filename("pyxtal", "database/wyckoff_generators.csv"))
layer_df = read_csv(resource_filename("pyxtal", "database/layer.csv"))
layer_symmetry_df = read_csv(resource_filename("pyxtal", "database/layer_symmetry.csv"))
layer_generators_df = read_csv(resource_filename("pyxtal", "database/layer_generators.csv"))
rod_df = read_csv(resource_filename("pyxtal", "database/rod.csv"))
rod_symmetry_df = read_csv(resource_filename("pyxtal", "database/rod_symmetry.csv"))
rod_generators_df = read_csv(resource_filename("pyxtal", "database/rod_generators.csv"))

#Define functions
#------------------------------

def filtered_coords(coords, PBC=[1, 2, 3]):
    """
    Given a list of 3d fractional coordinates or a single 3d point, transform
    all coordinates to less than 1 and greater than 0. If one axis is not
    periodic, does not transform the coordinates along that axis. For example,
    for the point [1.2,1.6, -.4] with periodicity along the x and z axes, but
    not the y axis (PBC=[1, 3]), the function would return [0.2, 1.6, 0.6].

    Args:
        coords: a list or array of real 3d vectors, or a single real 3d vector
        PBC: the axes, if any, which are periodic. 1, 2, and 3 correspond
            to x, y, and z repectively

    Returns:
        a new list of coordinates (or single point) with values scaled between
        0 and 1, except for values on the non-periodic axis
    """
    coords0 = np.array(coords)
    new_coords = coords0 - np.floor(coords0)
    if PBC is not None:
        if len(new_coords.shape) == 2:
            for x in range(0, 3):
                if x +1 not in PBC:
                    new_coords[:,x] = coords0[:,x]
        elif len(new_coords.shape) == 1:
            for x in range(0, 3):
                if x +1 not in PBC:
                    new_coords[x] = coords0[x]
        else:
            print("Warning: invalid array dimensions for filtered_coords. Shape: "+str(new_coords.shape))
            return coords
    return new_coords

def gaussian(min, max, sigma=3.0):
    """
    Choose a random number from a Gaussian probability distribution centered
    between min and max. sigma is the number of standard deviations that min
    and max are away from the center. Thus, sigma is also the largest possible
    number of standard deviations corresponding to the returned value. sigma=2
    corresponds to a 95.45% probability of choosing a number between min and
    max.

    Args:
        min: the minimum acceptable value
        max: the maximum acceptable value
        sigma: the number of standard deviations between the center and min or max

    Returns:
        a value chosen randomly between min and max
    """
    center = (max+min)*0.5
    delta = fabs(max-min)*0.5
    ratio = delta/sigma
    while True:
        x = np.random.normal(scale=ratio, loc=center)
        if x > min and x < max:
            return x
            
def letter_from_index(index, sg):
    """
    Given a Wyckoff position's index within a spacegroup, return its number
    and letter e.g. '4a'

    Args:
        index: a single integer describing the WP's index within the
            spacegroup (0 is the general position)
        sg: the international spacegroup number
   
    Returns:
        the Wyckoff letter corresponding to the Wyckoff position (for example,
        for position 4a, the function would return 'a')
    """
    letters = "abcdefghijklmnopqrstuvwxyzA"
    wyckoffs = get_wyckoffs(sg)
    length = len(wyckoffs)
    return letters[length - 1 - index]

def index_from_letter(letter, sg):
    """
    Given the Wyckoff letter, returns the index of a Wyckoff position within
    the spacegroup

    Args:
        letter: The wyckoff letter
        sg: the internationl spacegroup number

    Returns:
        a single index specifying the location of the Wyckoff position within
        the spacegroup (0 is the general position)
    """
    letters = "abcdefghijklmnopqrstuvwxyzA"
    wyckoffs = get_wyckoffs(sg)
    length = len(wyckoffs)
    return length - 1 - letters.index(letter)

def jk_from_i(i, olist):
    """
    Given an organized list (Wyckoff positions or orientations), determine the
    two indices which correspond to a single index for an unorganized list.
    Used mainly for organized Wyckoff position lists, but can be used for other
    lists organized in a similar way

    Args:
        i: a single index corresponding to the item's location in the
            unorganized list
        olist: the organized list

    Returns:
        [j, k]: two indices corresponding to the item's location in the
            organized list
    """
    num = -1
    found = False
    for j , a in enumerate(olist):
        for k , b in enumerate(a):
            num += 1
            if num == i:
                return [j, k]
    print("Error: Incorrect Wyckoff position list or index passed to jk_from_i")
    return None

def i_from_jk(j, k, olist):
    """
    Inverse operation of jk_from_i: gives one list index from 2

    Args:
        j, k: indices corresponding to the location of an element in the
            organized list
        olist: the organized list of Wyckoff positions or molecular orientations

    Returns:
        i: one index corresponding to the item's location in the
            unorganized list    
    """
    num = -1
    for x, a in enumerate(olist):
        for y, b in enumerate(a):
            num += 1
            if x == j and y == k:
                return num
    print("Error: Incorrect Wyckoff position list or index passed to jk_from_i")
    return None

def ss_string_from_ops(ops, sg, complete=True):
    """
    Print the Hermann-Mauguin symbol for a site symmetry group, using a list of
    SymmOps as input. Note that the symbol does not necessarily refer to the
    x,y,z axes. For information on reading these symbols, see:
    http://en.wikipedia.org/wiki/Hermann-Mauguin_notation#Point_groups

    Args:
        ops: a list of SymmOp objects representing the site symmetry
        sg: International number of the spacegroup. Used to determine which
            axes to show. For example, a 3-fold rotation in a cubic system is
            written as ".3.", whereas a 3-fold rotation in a trigonal system is
            written as "3.."
        complete: whether or not all symmetry operations in the group
            are present. If False, we generate the rest

    Returns:
        a string representing the site symmetry. Ex: "2mm"
    """
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
    #Return whether or not two axes are symmetrically equivalent
    #It is assumed that both axes possess the same symbol
    #Will be called within combine_axes
    def are_symmetrically_equivalent(index1, index2):
        axis1 = axes[index1]
        axis2 = axes[index2]
        condition1 = False
        condition2 = False
        #Check for an operation mapping one axis onto the other
        for op in ops:
            if condition1 is False or condition2 is False:
                new1 = op.operate(axis1)
                new2 = op.operate(axis2)
                if np.isclose(abs(np.dot(new1, axis2)), 1):
                    condition1 = True
                if np.isclose(abs(np.dot(new2, axis1)), 1):
                    condition2 = True
        if condition1 is True and condition2 is True:
            return True
        else:
            return False
    #Given a list of axis indices, return the combined symbol
    #Axes may or may not be symmetrically equivalent, but must be of the same
    #type (x/y/z, face-diagonal, body-diagonal)
    #Will be called for mid- and high-symmetry crystallographic point groups
    def combine_axes(indices):
        symbols = {}
        for index in deepcopy(indices):
            symbol = get_symbol(params[index],orders[index],reflections[index])
            if symbol == ".":
                indices.remove(index)
            else:
                symbols[index] = symbol
        if indices == []:
            return "."
        #Remove redundant axes
        for i in deepcopy(indices):
            for j in deepcopy(indices):
                if j > i:
                    if symbols[i] == symbols[j]:
                        if are_symmetrically_equivalent(i, j):
                            if j in indices:
                                indices.remove(j)
        #Combine symbols for non-equivalent axes
        new_symbols = []
        for i in indices:
            new_symbols.append(symbols[i])
        symbol = ""
        while new_symbols != []:
            highest = get_highest_symbol(new_symbols)
            symbol += highest
            new_symbols.remove(highest)
        if symbol == "":
            print("Error: could not combine site symmetry axes.")
            return
        else:
            return symbol
    #Generate needed ops
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
            [1,1,1],[-1,1,1],[1,-1,1],[1,1,-1]]
    for i, axis in enumerate(axes):
        axes[i] = axis/np.linalg.norm(axis)
    for opa in opas:
        if opa.type != "identity" and opa.type != "inversion":
            found = False
            for i, axis in enumerate(axes):
                if np.isclose(abs(np.dot(opa.axis, axis)), 1):
                    found = True
                    params[i].append(opa)
            #Store uncommon axes for trigonal and hexagonal lattices
            if found is False:
                axes.append(opa.axis)
                #Check that new axis is not symmetrically equivalent to others
                unique = True
                for i, axis in enumerate(axes):
                    if i != len(axes)-1:
                        if are_symmetrically_equivalent(i, len(axes)-1):
                            unique = False
                if unique is True:
                    params.append([opa])
                elif unique is False:
                    axes.pop()
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
        s2 = combine_axes([0,1])
        #3rd symbol: face-diagonal axes (whichever have highest symmetry)
        s3 = combine_axes(list(range(3, len(axes))))
        symbol = s1+" "+s2+" "+s3
        if symbol != ". . .":
            return symbol
        elif symbol == ". . .":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    #Cubic
    elif sg >= 195 and sg <= 230:
        pass
        #1st symbol: x, y, and/or z axes (whichever have highest symmetry)
        s1 = combine_axes([0,1,2])
        #2nd symbol: body-diagonal axes (whichever has highest symmetry)
        s2 = combine_axes([9,10,11,12])
        #3rd symbol: face-diagonal axes (whichever have highest symmetry)
        s3 = combine_axes([3,4,5,6,7,8])
        symbol = s1+" "+s2+" "+s3
        if symbol != ". . .":
            return symbol
        elif symbol == ". . .":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    else:
        print("Error: invalid spacegroup number")
        return

def create_matrix(PBC=[1,2,3]):
    """
    Used for calculating distances in lattices with periodic boundary
    conditions. When multiplied with a set of points, generates additional
    points in cells adjacent to and diagonal to the original cell

    Args:
        PBC: The axes with periodic boundary conditions.
            Ex: PBC=[2,3] cancels periodic boundary conditions along the x axis

    Returns:
        A numpy array of matrices which can be multiplied by a set of
        coordinates
    """
    matrix = []
    i_list = [-1, 0, 1]
    j_list = [-1, 0, 1]
    k_list = [-1, 0, 1]
    if 1 not in PBC:
        i_list = [0]
    if 2 not in PBC:
        j_list = [0]
    if 3 not in PBC:
        k_list = [0]
    for i in i_list:
        for j in j_list:
            for k in k_list:
                matrix.append([i,j,k])
    return np.array(matrix, dtype=float)

def distance(xyz, lattice, PBC=[1,2,3]):
    """
    Returns the Euclidean distance from the origin for a fractional
    displacement vector. Takes into account the lattice metric and periodic
    boundary conditions, including up to one non-periodic axis.
    
    Args:
        xyz: a fractional 3d displacement vector. Can be obtained by
            subtracting one fractional vector from another
        lattice: a 3x3 matrix describing a unit cell's lattice vectors
        PBC: the axes, if any, which are periodic. 1, 2, and 3 correspond
            to x, y, and z respectively.

    Returns:
        a scalar for the distance of the point from the origin
    """
    xyz = filtered_coords(xyz, PBC=PBC)
    matrix = create_matrix(PBC=PBC)
    matrix += xyz
    matrix = np.dot(matrix, lattice)
    return np.min(cdist(matrix,[[0,0,0]]))       

def check_distance(coord1, coord2, specie1, specie2, lattice, PBC=[1,2,3], d_factor=1.0):
    """
    Check the distances between two set of atoms. The first set is
    generally larger than the second. Distances between coordinates within the
    first set are not checked, and distances between coordinates within the
    second set are not checked. Only distances between points from different
    sets are checked.

    Args:
        coord1: multiple lists of fractional coordinates e.g. [[[.1,.6,.4]
            [.3,.8,.2]],[[.4,.4,.4],[.3,.3,.3]]]
        coord2: a list of new fractional coordinates e.g. [[.7,.8,.9],
            [.4,.5,.6]]
        specie1: a list of atomic symbols for coord1. Ex: ['C', 'O']
        specie2: the atomic symbol for coord2. Ex: 'Li'
        lattice: matrix describing the unit cell vectors
        PBC: the axes, if any, which are periodic. 1, 2, and 3 correspond
            to x, y, and z respectively.
        d_factor: the tolerance is multiplied by this amount. Larger values
            mean atoms must be farther apart

    Returns:
        a bool for whether or not the atoms are sufficiently far enough apart
    """
    #add PBC
    coord2s = []
    matrix = create_matrix(PBC=PBC)
    for coord in coord2:
        for m in matrix:
            coord2s.append(coord+m)
    coord2 = np.array(coord2s)

    coord2 = np.dot(coord2, lattice)
    if len(coord1)>0:
        for coord, element in zip(coord1, specie1):
            coord = np.dot(coord, lattice)
            d_min = np.min(cdist(coord, coord2))
            tol = d_factor*0.5*(Element(element).covalent_radius + Element(specie2).covalent_radius)
            if d_min < tol:
                return False
        return True
    else:
        return True

def get_center(xyzs, lattice, PBC=[1,2,3]):
    """
    Finds the geometric centers of the clusters under periodic boundary
    conditions.

    Args:
        xyzs: a list of fractional coordinates
        lattice: a matrix describing the unit cell
        PBC: the axes, if any, which are periodic. 1, 2, and 3 correspond
            to x, y, and z respectively.

    Returns:
        x,y,z coordinates for the center of the input coordinate list
    """
    matrix0 = create_matrix(PBC=PBC)
    xyzs -= np.round(xyzs)
    matrix_min = [0,0,0]
    for atom1 in range(1,len(xyzs)):
        dist_min = 10.0
        for atom2 in range(0, atom1):
            #shift atom1 to position close to atom2
            matrix = matrix0 + (xyzs[atom1] - xyzs[atom2])
            matrix = np.dot(matrix, lattice)
            dists = cdist(matrix, [[0,0,0]])
            if np.min(dists) < dist_min:
                dist_min = np.min(dists)
                matrix_min = matrix0[np.argmin(dists)]
        xyzs[atom1] += matrix_min
    center = xyzs.mean(0)
    '''for a in range(1, 4):
        if a not in PBC:
            if abs(center[a-1])<1e-4:
                center[a-1] = 0.5'''
    return center

def para2matrix(cell_para, radians=True, format='lower'):
    """
    Given a set of lattic parameters, generates a matrix representing the
    lattice vectors

    Args:
        cell_para: a 1x6 list of lattice parameters [a, b, c, alpha, beta,
            gamma]. a, b, and c are the length of the lattice vectos, and
            alpha, beta, and gamma are the angles between these vectors. Can
            be generated by matrix2para
        radians: if True, lattice parameters should be in radians. If False,
            lattice angles should be in degrees
        format: a string ('lower', 'symmetric', or 'upper') for the type of
            matrix to be output

    Returns:
        a 3x3 matrix representing the unit cell. By default (format='lower'),
        the a vector is aligined along the x-axis, and the b vector is in the
        y-z plane
    """
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

def Add_vacuum(lattice, coor, vacuum=10.0, dim = 2):
    """
    Adds space above and below a 2D crystal. This allows for treating the
    structure as a 3D crystal during energy optimization

    Args:
        lattice: the lattice matrix of the crystal
        coor: the relative coordinates of the crystal
        vacuum: the amount of space, in Angstroms, to add above and below
        dim: the axis (0,1,2)->(x,y,z) along which to add space

    Returns:
        lattice, coor: The transformed lattice and coordinates after the
            vacuum space is added
    """
    old = lattice[dim, dim]
    new = old + vacuum
    coor[:,dim] = coor[:,dim]*old/new
    coor[:,dim] = coor[:,dim] - np.mean(coor[:,dim]) + 0.5
    lattice[dim, dim] = new
    return lattice, coor

def Permutation(lattice, coor, PB):
    """
    Permutes a list of coordinates. Not currently implemented.
    """
    para = matrix2para(lattice)
    para1 = deepcopy(para)
    coor1 = deepcopy(coor)
    for axis in [0,1,2]:
        para1[axis] = para[PB[axis]-1]
        para1[axis+3] = para[PB[axis]+2]
        coor1[:,axis] = coor[:,PB[axis]-1]
    return para2matrix(para1), coor1

def matrix2para(matrix, radians=True):
    """
    Given a 3x3 matrix representing a unit cell, outputs a list of lattice
    parameters.

    Args:
        matrix: a 3x3 array or list, where the first, second, and third rows
            represent the a, b, and c vectors respectively
        radians: if True, outputs angles in radians. If False, outputs in
            degrees

    Returns:
        a 1x6 list of lattice parameters [a, b, c, alpha, beta, gamma]. a, b,
        and c are the length of the lattice vectos, and alpha, beta, and gamma
        are the angles between these vectors
    """
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
    Returns the number of duplicate atoms in the conventional lattice (in
    contrast to the primitive cell). Based on the type of cell centering (P,
    A, C, I, R, or F)

    Args:
        sg: the international spacegroup number

    Returns:
        a number between 1 and 4
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

def find_short_dist(coor, lattice, tol, PBC=[1,2,3]):
    """
    Given a list of fractional coordinates, finds pairs which are closer
    together than tol, and builds the connectivity map

    Args:
        coor: a list of fractional 3-dimensional coordinates
        lattice: a matrix representing the crystal unit cell
        tol: the distance tolerance for pairing coordinates
        PBC: the axes, if any, which are periodic. 1, 2, and 3 correspond
            to x, y, and z respectively.
    
    Returns:
        pairs, graph: (pairs) is a list whose entries have the form [index1,
        index2, distance], where index1 and index2 correspond to the indices
        of a pair of points within the supplied list (coor). distance is the
        distance between the two points. (graph) is a connectivity map in the
        form of a list. Its first index represents a point within coor, and
        the second indices represent which point(s) it is connected to.
    """
    pairs=[]
    graph=[]
    for i in range(len(coor)):
        graph.append([])

    for i1 in range(len(coor)-1):
        for i2 in range(i1+1,len(coor)):
            dist = distance(coor[i1]-coor[i2], lattice, PBC=PBC)
            if dist <= tol:
                #dists.append(dist)
                pairs.append([i1,i2,dist])
    pairs = np.array(pairs)
    if len(pairs) > 0:
        d_min = min(pairs[:,-1]) + 1e-3
        sequence = [pairs[:,-1] <= d_min]
        #Avoid Futurewarning
        #pairs1 = deepcopy(pairs)
        #pairs = pairs1[sequence]
        pairs = pairs[tuple(sequence)]
        for pair in pairs:
            pair0=int(pair[0])
            pair1=int(pair[1])
            graph[pair0].append(pair1)
            graph[pair1].append(pair0)

    return pairs, graph

def connected_components(graph):
    """
    Given an undirected graph (a 2d array of indices), return a set of
    connected components, each connected component being an (arbitrarily
    ordered) array of indices which are connected either directly or
    indirectly.

    Args:
        graph: a list reprenting the connections between points. The first index
            represents a point, and the 2nd indices represent the points to
            which the first point is connected. Can be generated by
            find_short_dist

    Returns:
        a list of connected components. The first index denotes a separate
        connected component. The second indices denote the points within the
        connected component which are connected to each other
    """
    def add_neighbors(el, seen=[]):
        """
        Find all elements which are connected to el. Return an array which
        includes these elements and el itself.
        """
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

def merge_coordinate(coor, lattice, wyckoffs, w_symm_all, tol, PBC=[1,2,3]):
    """
    Given a list of fractional coordinates, merges them within a given
    tolerance, and checks if the merged coordinates satisfy a Wyckoff
    position. Used for merging general Wyckoff positions into special Wyckoff
    positions within the random_crystal (and its derivative) classes.

    Args:
        coor: a list of fractional coordinates
        lattice: a 3x3 matrix representing the unit cell
        wyckoffs: an unorganized list of Wyckoff positions to check
        w_symm_all: A list of Wyckoff site symmetry obtained from
            get_wyckoff_symmetry
        tol: the cutoff distance for merging coordinates
        PBC: the axes, if any, which are periodic. 1, 2, and 3 correspond
            to x, y, and z respectively.

    Returns:
        coor, index: (coor) is the new list of fractional coordinates after
        merging, and index is a single index of the Wyckoff position within
        the spacegroup. If merging is unsuccesful, or no index is found,
        returns the original coordinates and False
    """
    while True:
        pairs, graph = find_short_dist(coor, lattice, tol, PBC=PBC)
        index = None
        if len(pairs)>0:
            if len(coor) > len(wyckoffs[-1]):
                merged = []
                groups = connected_components(graph)
                for group in groups:
                    merged.append(get_center(coor[group], lattice, PBC=PBC))
                merged = np.array(merged)
                index = check_wyckoff_position(merged, wyckoffs, w_symm_all, exact_translation=False, PBC=PBC)
                if index is False:
                    return coor, False
                else:
                    coor = merged

            else:#no way to merge
                return coor, False
        else:
            if index is None:
                index = check_wyckoff_position(coor, wyckoffs, w_symm_all, exact_translation=False, PBC=PBC)
            return coor, index

def estimate_volume(numIons, species, factor=1.0):
    """
    Estimates the volume of a unit cell based on the number and types of ions.
    Assumes each atom takes up a sphere with radius equal to its covalent bond
    radius.

    Args:
        numIons: a list of the number of ions for each specie
        species: a corresponding list for the specie of each type of ion. Each
            element in the list should be a string for the atomic symbol
        factor: an optional factor to multiply the result by. Larger values
            allow more space between atoms
    
    Returns:
        a float value for the estimated volume
    """
    volume = 0
    for numIon, specie in zip(numIons, species):
        r = rand(Element(specie).covalent_radius, Element(specie).vdw_radius)
        volume += numIon*4/3*pi*r**3
    return factor*volume

def generate_lattice(sg, volume, minvec=tol_m, minangle=pi/6, max_ratio=10.0, maxattempts = 100):
    """
    Generates a lattice (3x3 matrix) according to the space group symmetry and
    number of atoms. If the spacegroup has centering, we will transform to
    conventional cell setting. If the generated lattice does not meet the
    minimum angle and vector requirements, we try to generate a new one, up to
    maxattempts times.

    Args:
        sg: International number of the space group
        volume: volume of the conventional unit cell
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
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
        if minvec < maxvec:
            #Check minimum Euclidean distances
            smallvec = min(a*cos(max(beta, gamma)), b*cos(max(alpha, gamma)), c*cos(max(alpha, beta)))
            if(a>minvec and b>minvec and c>minvec
            and a<maxvec and b<maxvec and c<maxvec
            and smallvec < minvec
            and alpha>minangle and beta>minangle and gamma>minangle
            and alpha<maxangle and beta<maxangle and gamma<maxangle
            and a/b<max_ratio and a/c<max_ratio and b/c<max_ratio
            and b/a<max_ratio and c/a<max_ratio and c/b<max_ratio):
                return np.array([a, b, c, alpha, beta, gamma])
    #If maxattempts tries have been made without success
    print("Error: Could not generate lattice after "+str(n+1)+" attempts for volume ", volume)
    return

def generate_lattice_2D(num, volume, thickness=None, minvec=tol_m, minangle=pi/6, max_ratio=10.0, maxattempts = 100):
    """
    Generates a lattice (3x3 matrix) according to the spacegroup symmetry and
    number of atoms. If the layer group has centering, we will use the
    conventional cell setting. If the generated lattice does not meet the
    minimum angle and vector requirements, we try to generate a new one, up to
    maxattempts times.
    Note: The monoclinic layer groups have different unique axes. Groups 3-7
        have unique axis c, while 8-18 have unique axis a. We use non-periodic
        axis c for all layer groups.

    Args:
        num: International number of the space group
        volume: volume of the lattice
        thickness: 3rd-dimensional thickness of the unit cell. If set to None,
            a thickness is chosen automatically
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    #Store the non-periodic axis
    NPA = 3
    #Set the unique axis for monoclinic cells
    if num in range(3, 8): unique_axis = "c"
    elif num in range(8, 19): unique_axis = "a"
    maxangle = pi-minangle
    for n in range(maxattempts):
        abc = np.ones([3])
        if thickness is None:
            v = random_vector()
            thickness1 = np.cbrt(volume)*(v[0]/(v[0]*v[1]*v[2]))
        else:
            thickness1 = thickness
        abc[NPA-1] = thickness1
        alpha, beta, gamma  = pi/2, pi/2, pi/2
        #Triclinic
        if num <= 2:
            mat = random_shear_matrix(width=0.2)
            a, b, c, alpha, beta, gamma = matrix2para(mat)
            x = sqrt(1-cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*(cos(alpha)*cos(beta)*cos(gamma)))
            abc[NPA-1] = abc[NPA-1]/x #scale thickness by outer product of vectors
            ab = volume/(abc[NPA-1]*x)
            ratio = a/b
            if NPA == 3:
                abc[0] = sqrt(ab*ratio)
                abc[1] = sqrt(ab/ratio)
            elif NPA == 2:
                abc[0] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)
            elif NPA == 1:
                abc[1] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)

        #Monoclinic
        elif num <= 18:
            a, b, c = random_vector()
            if unique_axis == "a":
                alpha = gaussian(minangle, maxangle)
                x = sin(alpha)
            elif unique_axis == "b":
                beta = gaussian(minangle, maxangle)
                x = sin(beta)
            elif unique_axis == "c":
                gamma = gaussian(minangle, maxangle)
                x = sin(gamma)
            ab = volume/(abc[NPA-1]*x)
            ratio = a/b
            if NPA == 3:
                abc[0] = sqrt(ab*ratio)
                abc[1] = sqrt(ab/ratio)
            elif NPA == 2:
                abc[0] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)
            elif NPA == 1:
                abc[1] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)

        #Orthorhombic
        elif num <= 48:
            vec = random_vector()
            if NPA == 3:
                ratio = abs(vec[0]/vec[1]) #ratio a/b
                abc[1] = sqrt(volume/(thickness1*ratio))
                abc[0] = abc[1]* ratio
            elif NPA == 2:
                ratio = abs(vec[0]/vec[2]) #ratio a/b
                abc[2] = sqrt(volume/(thickness1*ratio))
                abc[0] = abc[2]* ratio
            elif NPA == 1:
                ratio = abs(vec[1]/vec[2]) #ratio a/b
                abc[2] = sqrt(volume/(thickness1*ratio))
                abc[1] = abc[2]* ratio

        #Tetragonal
        elif num <= 64:
            if NPA == 3:
                abc[0] = abc[1] = sqrt(volume/thickness1)
            elif NPA == 2:
                abc[0] = abc[1]
                abc[2] = volume/(abc[NPA-1]**2)
            elif NPA == 1:
                abc[1] = abc[0]
                abc[2] = volume/(abc[NPA-1]**2)

        #Trigonal/Rhombohedral/Hexagonal
        elif num <= 80:
            gamma = pi/3*2
            x = sqrt(3.)/2.
            if NPA == 3:
                abc[0] = abc[1] = sqrt((volume/x)/abc[NPA-1])
            elif NPA == 2:
                abc[0] = abc[1]
                abc[2] = (volume/x)(thickness1**2)
            elif NPA == 1:
                abc[1] = abc[0]
                abc[2] = (volume/x)/(thickness1**2)

        para = np.array([abc[0], abc[1], abc[2], alpha, beta, gamma])

        a, b, c = abc[0], abc[1], abc[2]
        maxvec = (a*b*c)/(minvec**2)
        if minvec < maxvec:
            smallvec = min(a*cos(max(beta, gamma)), b*cos(max(alpha, gamma)), c*cos(max(alpha, beta)))
            if(a>minvec and b>minvec and c>minvec
            and a<maxvec and b<maxvec and c<maxvec
            and smallvec < minvec
            and alpha>minangle and beta>minangle and gamma>minangle
            and alpha<maxangle and beta<maxangle and gamma<maxangle
            and a/b<max_ratio and a/c<max_ratio and b/c<max_ratio
            and b/a<max_ratio and c/a<max_ratio and c/b<max_ratio):
                return para

    #If maxattempts tries have been made without success
    print("Error: Could not generate lattice after "+str(n+1)+" attempts")
    return

def generate_lattice_1D(num, volume, area=None, minvec=tol_m, minangle=pi/6, max_ratio=10.0, maxattempts = 100):
    """
    Generates a lattice (3x3 matrix) according to the spacegroup symmetry and
    number of atoms. If the spacegroup has centering, we will transform to
    conventional cell setting. If the generated lattice does not meet the
    minimum angle and vector requirements, we try to generate a new one, up to
    maxattempts times.
    Note: The monoclinic Rod groups have different unique axes. Groups 3-7
        have unique axis a, while 8-12 have unique axis c. We use periodic
        axis c for all Rod groups.

    Args:
        num: number of the Rod group
        volume: volume of the lattice
        area: cross-sectional area of the unit cell in Angstroms squared. If
            set to None, a value is chosen automatically
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    #Store the periodic axis
    PA = 3
    #Set the unique axis for monoclinic cells
    if num in range(3, 8): unique_axis = "a"
    elif num in range(8, 13): unique_axis = "c"
    maxangle = pi-minangle
    for n in range(maxattempts):
        abc = np.ones([3])
        if area is None:
            v = random_vector()
            thickness1 = np.cbrt(volume)*(v[0]/(v[0]*v[1]*v[2]))
        else:
            thickness1 = volume/area
        abc[PA-1] = thickness1
        alpha, beta, gamma  = pi/2, pi/2, pi/2
        #Triclinic
        if num <= 2:
            mat = random_shear_matrix(width=0.2)
            a, b, c, alpha, beta, gamma = matrix2para(mat)
            x = sqrt(1-cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*(cos(alpha)*cos(beta)*cos(gamma)))
            abc[PA-1] = abc[PA-1]/x #scale thickness by outer product of vectors
            ab = volume/(abc[PA-1]*x)
            ratio = a/b
            if PA == 3:
                abc[0] = sqrt(ab*ratio)
                abc[1] = sqrt(ab/ratio)
            elif PA == 2:
                abc[0] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)
            elif PA == 1:
                abc[1] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)

        #Monoclinic
        elif num <= 12:
            a, b, c = random_vector()
            if unique_axis == "a":
                alhpa = gaussian(minangle, maxangle)
                x = sin(alpha)
            elif unique_axis == "b":
                beta = gaussian(minangle, maxangle)
                x = sin(beta)
            elif unique_axis == "c":
                gamma = gaussian(minangle, maxangle)
                x = sin(gamma)
            ab = volume/(abc[PA-1]*x)
            ratio = a/b
            if PA == 3:
                abc[0] = sqrt(ab*ratio)
                abc[1] = sqrt(ab/ratio)
            elif PA == 2:
                abc[0] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)
            elif PA == 1:
                abc[1] = sqrt(ab*ratio)
                abc[2] = sqrt(ab/ratio)

        #Orthorhombic
        elif num <= 22:
            vec = random_vector()
            if PA == 3:
                ratio = abs(vec[0]/vec[1]) #ratio a/b
                abc[1] = sqrt(volume/(thickness1*ratio))
                abc[0] = abc[1]* ratio
            elif PA == 2:
                ratio = abs(vec[0]/vec[2]) #ratio a/b
                abc[2] = sqrt(volume/(thickness1*ratio))
                abc[0] = abc[2]* ratio
            elif PA == 1:
                ratio = abs(vec[1]/vec[2]) #ratio a/b
                abc[2] = sqrt(volume/(thickness1*ratio))
                abc[1] = abc[2]* ratio

        #Tetragonal
        elif num <= 41:
            if PA == 3:
                abc[0] = abc[1] = sqrt(volume/thickness1)
            elif PA == 2:
                abc[0] = abc[1]
                abc[2] = volume/(abc[PA-1]**2)
            elif PA == 1:
                abc[1] = abc[0]
                abc[2] = volume/(abc[PA-1]**2)

        #Trigonal/Rhombohedral/Hexagonal
        elif num <= 75:
            gamma = pi/3*2
            x = sqrt(3.)/2.
            if PA == 3:
                abc[0] = abc[1] = sqrt((volume/x)/abc[PA-1])
            elif PA == 2:
                abc[0] = abc[1]
                abc[2] = (volume/x)(thickness1**2)
            elif PA == 1:
                abc[1] = abc[0]
                abc[2] = (volume/x)/(thickness1**2)

        para = np.array([abc[0], abc[1], abc[2], alpha, beta, gamma])

        a, b, c = abc[0], abc[1], abc[2]
        maxvec = (a*b*c)/(minvec**2)
        if minvec < maxvec:
            smallvec = min(a*cos(max(beta, gamma)), b*cos(max(alpha, gamma)), c*cos(max(alpha, beta)))
            if(a>minvec and b>minvec and c>minvec
            and a<maxvec and b<maxvec and c<maxvec
            and smallvec < minvec
            and alpha>minangle and beta>minangle and gamma>minangle
            and alpha<maxangle and beta<maxangle and gamma<maxangle
            and a/b<max_ratio and a/c<max_ratio and b/c<max_ratio
            and b/a<max_ratio and c/a<max_ratio and c/b<max_ratio):
                return para

    #If maxattempts tries have been made without success
    print("Error: Could not generate lattice after "+str(n+1)+" attempts")
    return

def choose_wyckoff(wyckoffs, number):
    """
    Choose a Wyckoff position to fill based on the current number of atoms
    needed to be placed within a unit cell
    Rules:
        1) The new position's multiplicity is equal/less than (number).
        2) We prefer positions with large multiplicity.

    Args:
        wyckoffs: an organized list of Wyckoff positions
        number: the number of atoms still needed in the unit cell

    Returns:
        a single index for the Wyckoff position. If no position is found,
        returns False
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

def get_wyckoffs(sg, organized=False, PBC=[1,2,3]):
    """
    Returns a list of Wyckoff positions for a given space group. Has option to
    organize the list based on multiplicity (this is used for
    random_crystal.wyckoffs) For an unorganized list:

    1st index: index of WP in sg (0 is the WP with largest multiplicity)

    2nd index: a SymmOp object in the WP

    For an organized list:

    1st index: specifies multiplicity (0 is the largest multiplicity)

    2nd index: corresponds to a Wyckoff position within the group of equal
        multiplicity.

    3nd index: corresponds to a SymmOp object within the Wyckoff position

    You may switch between organized and unorganized lists using the methods
    i_from_jk and jk_from_i. For example, if a Wyckoff position is the [i]
    entry in an unorganized list, it will be the [j][k] entry in an organized
    list.

    Args:
        sg: the international spacegroup number
        organized: whether or not to organize the list based on multiplicity
        PBC: a list of periodic axes (1,2,3)->(x,y,z)
    
    Returns: 
        a list of Wyckoff positions, each of which is a list of SymmOp's
    """
    if PBC != [1,2,3]:
        for a in range(1, 4):
            if a not in PBC:
                coor = [0,0,0]
                coor[a-1] = 0.5
        coor = np.array(coor)

    wyckoff_strings = eval(wyckoff_df["0"][sg])
    wyckoffs = []
    for x in wyckoff_strings:
        if PBC != [1,2,3]:
            op = SymmOp.from_xyz_string(x[0])
            coor1 = op.operate(coor)
            invalid = False
            for a in range(1, 4):
                if a not in PBC:
                    if abs(coor1[a-1]-0.5) < 1e-2:
                        pass
                    else:
                        #invalid wyckoffs for layer group
                        invalid = True
            if invalid == False:
                wyckoffs.append([])
                for y in x:
                    wyckoffs[-1].append(SymmOp.from_xyz_string(y))
        else:
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

def get_layer(num, organized=False):
    """
    Returns a list of Wyckoff positions for a given 2D layer group. Has
    option to organize the list based on multiplicity (this is used for
    random_crystal_2D.wyckoffs) For an unorganized list:

    1st index: index of WP in layer group (0 is the WP with largest multiplicity)

    2nd index: a SymmOp object in the WP

    For an organized list:

    1st index: specifies multiplicity (0 is the largest multiplicity)

    2nd index: corresponds to a Wyckoff position within the group of equal
        multiplicity.

    3nd index: corresponds to a SymmOp object within the Wyckoff position

    You may switch between organized and unorganized lists using the methods
    i_from_jk and jk_from_i. For example, if a Wyckoff position is the [i]
    entry in an unorganized list, it will be the [j][k] entry in an organized
    list.

    For layer groups with more than one possible origin, origin choice 2 is
    used.

    Args:
        num: the international layer group number
        organized: whether or not to organize the list based on multiplicity
    
    Returns: 
        a list of Wyckoff positions, each of which is a list of SymmOp's
    """
    wyckoff_strings = eval(layer_df["0"][num])
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

def get_rod(num, organized=False):
    """
    Returns a list of Wyckoff positions for a given 1D Rod group. Has option to
    organize the list based on multiplicity (this is used for
    random_crystal_1D.wyckoffs) For an unorganized list:

    1st index: index of WP in layer group (0 is the WP with largest multiplicity)

    2nd index: a SymmOp object in the WP

    For an organized list:

    1st index: specifies multiplicity (0 is the largest multiplicity)

    2nd index: corresponds to a Wyckoff position within the group of equal
        multiplicity.

    3nd index: corresponds to a SymmOp object within the Wyckoff position

    You may switch between organized and unorganized lists using the methods
    i_from_jk and jk_from_i. For example, if a Wyckoff position is the [i]
    entry in an unorganized list, it will be the [j][k] entry in an organized
    list.

    For Rod groups with more than one possible setting, setting choice 1
    is used.

    Args:
        num: the international Rod group number
        organized: whether or not to organize the list based on multiplicity
    
    Returns: 
        a list of Wyckoff positions, each of which is a list of SymmOp's
    """
    wyckoff_strings = eval(rod_df["0"][num])
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

def get_wyckoff_symmetry(sg, PBC=[1,2,3], molecular=False):
    """
    Returns a list of Wyckoff position site symmetry for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a point within the WP
    3rd index: a site symmetry SymmOp of the point

    Args:
        sg: the international spacegroup number
        PBC: a list of periodic axes (1,2,3)->(x,y,z)
        molecular: whether or not to return the Euclidean point symmetry
            operations. If True, cuts off translational part of operation, and
            converts non-orthogonal operations (3-fold and 6-fold rotations)
            to (orthogonal) pure rotations. Should be used when dealing with
            molecular crystals

    Returns:
        a 3d list of SymmOp objects representing the site symmetry of each
        point in each Wyckoff position
    """
    if PBC != [1,2,3]:
        coor = [0,0,0]
        for a in range(1,4):
            if a not in PBC:
                coor[a-1] = 0.5
        coor = np.array(coor)
    wyckoffs = get_wyckoffs(sg, PBC=PBC)

    P = SymmOp.from_rotation_and_translation([[1,-.5,0],[0,sqrt(3)/2,0],[0,0,1]], [0,0,0])
    symmetry_strings = eval(wyckoff_symmetry_df["0"][sg])
    symmetry = []
    convert = False
    if molecular is True:
        if sg >= 143 and sg <= 194:
            convert = True
    #Loop over Wyckoff positions
    for x, w in zip(symmetry_strings, wyckoffs):
        if PBC != [1,2,3]:
            op = w[0]
            coor1 = op.operate(coor)
            invalid = False
            for a in range(1,4):
                if a not in PBC:
                    if abs(coor1[a-1]-0.5) < 1e-2:
                        pass
                    else:
                        invalid = True
            if invalid == False:
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
        else:
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

def get_layer_symmetry(num, molecular=False):
    """
    Returns a list of Wyckoff position site symmetry for a given space group.
    1st index: index of WP in group (0 is the WP with largest multiplicity)
    2nd index: a point within the WP
    3rd index: a site symmetry SymmOp of the point

    Args:
        num: the layer group number
        molecular: whether or not to return the Euclidean point symmetry
            operations. If True, cuts off translational part of operation, and
            converts non-orthogonal operations (3-fold and 6-fold rotations)
            to (orthogonal) pure rotations. Should be used when dealing with
            molecular crystals

    Returns:
        a 3d list of SymmOp objects representing the site symmetry of each
        point in each Wyckoff position
    """

    P = SymmOp.from_rotation_and_translation([[1,-.5,0],[0,sqrt(3)/2,0],[0,0,1]], [0,0,0])
    symmetry_strings = eval(layer_symmetry_df["0"][num])
    symmetry = []
    convert = False
    if molecular is True:
        if num >= 65:
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

def get_rod_symmetry(num, molecular=False):
    """
    Returns a list of Wyckoff position site symmetry for a given Rod group.
    1st index: index of WP in group (0 is the WP with largest multiplicity)
    2nd index: a point within the WP
    3rd index: a site symmetry SymmOp of the point

    Args:
        num: the Rod group number
        molecular: whether or not to return the Euclidean point symmetry
            operations. If True, cuts off translational part of operation, and
            converts non-orthogonal operations (3-fold and 6-fold rotations)
            to (orthogonal) pure rotations. Should be used when dealing with
            molecular crystals

    Returns:
        a 3d list of SymmOp objects representing the site symmetry of each
        point in each Wyckoff position
    """

    P = SymmOp.from_rotation_and_translation([[1,-.5,0],[0,sqrt(3)/2,0],[0,0,1]], [0,0,0])
    symmetry_strings = eval(rod_symmetry_df["0"][num])
    symmetry = []
    convert = False
    if molecular is True:
        if num >= 42:
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

def get_wyckoff_generators(sg, PBC=[1,2,3], molecular=False):
    """
    Returns a list of Wyckoff generators for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a generator for the WP
    This function is useful for rotating molecules based on Wyckoff position,
    since special Wyckoff positions only encode positional information, but not
    information about the orientation. The generators for each Wyckoff position
    form a subset of the spacegroup's general Wyckoff position.
    
    Args:
        sg: the international spacegroup number
        PBC: a list of periodic axes (1,2,3)->(x,y,z)
        molecular: whether or not to return the Euclidean point symmetry
            operations. If True, cuts off translational part of operation, and
            converts non-orthogonal operations (3-fold and 6-fold rotations)
            to (orthogonal) pure rotations. Should be used when dealing with
            molecular crystals
    
    Returns:
        a 2d list of SymmOp objects which can be used to generate a Wyckoff position given a
        single fractional (x,y,z) coordinate
    """
    if PBC != [1,2,3]:
        coor = [0,0,0]
        for a in range(1,4):
            if a not in PBC:
                coor[a-1] = 0.5
        coor = np.array(coor)
    wyckoffs = get_wyckoffs(sg, PBC=PBC)

    P = SymmOp.from_rotation_and_translation([[1,-.5,0],[0,sqrt(3)/2,0],[0,0,1]], [0,0,0])
    generator_strings = eval(wyckoff_generators_df["0"][sg])
    generators = []
    convert = False
    if molecular is True:
        if sg >= 143 and sg <= 194:
            convert = True
    #Loop over Wyckoff positions
    for x, w in zip(generator_strings, wyckoffs):
        if PBC != [1,2,3]:
            op = w[0]
            coor1 = op.operate(coor)
            invalid = False
            for a in range(1,4):
                if a not in PBC:
                    if abs(coor1[a-1]-0.5) < 1e-2:
                        pass
                    else:
                        invalid = True
            if invalid == False:
                generators.append([])
                #Loop over ops
                for y in x:
                    op = SymmOp.from_xyz_string(y)
                    if convert is True:
                        #Convert non-orthogonal trigonal/hexagonal operations
                        op = P*op*P.inverse
                    if molecular is False:
                        generators[-1].append(op)
                    elif molecular is True:
                        op = SymmOp.from_rotation_and_translation(op.rotation_matrix,[0,0,0])
                        generators[-1].append(op)
        else:
            generators.append([])
            for y in x:
                op = SymmOp.from_xyz_string(y)
                if convert is True:
                    #Convert non-orthogonal trigonal/hexagonal operations
                    op = P*op*P.inverse
                if molecular is False:
                    generators[-1].append(op)
                elif molecular is True:
                    op = SymmOp.from_rotation_and_translation(op.rotation_matrix,[0,0,0])
                    generators[-1].append(op)
    return generators

def get_layer_generators(num, molecular=False):
    """
    Returns a list of Wyckoff generators for a given layer group.
    1st index: index of WP in group (0 is the WP with largest multiplicity)
    2nd index: a generator for the WP
    This function is useful for rotating molecules based on Wyckoff position,
    since special Wyckoff positions only encode positional information, but not
    information about the orientation. The generators for each Wyckoff position
    form a subset of the group's general Wyckoff position.
    
    Args:
        num: the layer group number
        molecular: whether or not to return the Euclidean point symmetry
            operations. If True, cuts off translational part of operation, and
            converts non-orthogonal operations (3-fold and 6-fold rotations)
            to (orthogonal) pure rotations. Should be used when dealing with
            molecular crystals
    
    Returns:
        a 2d list of SymmOp objects which can be used to generate a Wyckoff position given a
        single fractional (x,y,z) coordinate
    """

    P = SymmOp.from_rotation_and_translation([[1,-.5,0],[0,sqrt(3)/2,0],[0,0,1]], [0,0,0])
    generator_strings = eval(layer_generators_df["0"][num])
    generators = []
    convert = False
    if molecular is True:
        if num >= 65:
            convert = True
    #Loop over Wyckoff positions
    for x in generator_strings:
        generators.append([])
        #Loop over ops
        for y in x:
            op = SymmOp.from_xyz_string(y)
            if convert is True:
                #Convert non-orthogonal trigonal/hexagonal operations
                op = P*op*P.inverse
            if molecular is False:
                generators[-1].append(op)
            elif molecular is True:
                op = SymmOp.from_rotation_and_translation(op.rotation_matrix,[0,0,0])
                generators[-1].append(op)
    return generators

def get_rod_generators(num, molecular=False):
    """
    Returns a list of Wyckoff generators for a given Rod group.
    1st index: index of WP in group (0 is the WP with largest multiplicity)
    2nd index: a generator for the WP
    This function is useful for rotating molecules based on Wyckoff position,
    since special Wyckoff positions only encode positional information, but not
    information about the orientation. The generators for each Wyckoff position
    form a subset of the group's general Wyckoff position.
    
    Args:
        num: the Rod group number
        molecular: whether or not to return the Euclidean point symmetry
            operations. If True, cuts off translational part of operation, and
            converts non-orthogonal operations (3-fold and 6-fold rotations)
            to (orthogonal) pure rotations. Should be used when dealing with
            molecular crystals
    
    Returns:
        a 2d list of SymmOp objects which can be used to generate a Wyckoff position given a
        single fractional (x,y,z) coordinate
    """

    P = SymmOp.from_rotation_and_translation([[1,-.5,0],[0,sqrt(3)/2,0],[0,0,1]], [0,0,0])
    generator_strings = eval(rod_generators_df["0"][num])
    generators = []
    convert = False
    if molecular is True:
        if num >= 42:
            convert = True
    #Loop over Wyckoff positions
    for x in generator_strings:
        generators.append([])
        #Loop over ops
        for y in x:
            op = SymmOp.from_xyz_string(y)
            if convert is True:
                #Convert non-orthogonal trigonal/hexagonal operations
                op = P*op*P.inverse
            if molecular is False:
                generators[-1].append(op)
            elif molecular is True:
                op = SymmOp.from_rotation_and_translation(op.rotation_matrix,[0,0,0])
                generators[-1].append(op)
    return generators

def site_symm(point, gen_pos, tol=1e-3, lattice=Euclidean_lattice, PBC=[1,2,3]):
    """
    Given a point and a general Wyckoff position, return the list of symmetry
    operations leaving the point (coordinate or SymmOp) invariant. The returned
    SymmOps are a subset of the general position. The site symmetry can be used
    for determining the Wyckoff position for a set of points, or for
    determining the valid orientations of a molecule within a given Wyckoff
    position.

    Args:
        point: a 1x3 coordinate or SymmOp object to find the symmetry of. If a
            SymmOp is given, the returned symmetries must also preserve the
            point's orientaion
        gen_pos: the general position of the spacegroup. Can be obtained using
            get_wyckoffs(sg)[0], where sg is the desired spacegroup number
        tol:
            the numberical tolerance for determining equivalent positions and
            orientations.
        lattice:
            a 3x3 matrix representing the lattice vectors of the unit cell
        PBC: a list of periodic axes (1,2,3)->(x,y,z)

    Returns:
        a list of SymmOp objects which leave the given point invariant
    """
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
        if distance(displacement, lattice, PBC=PBC) > tol:
            is_symmetry = False
        if is_symmetry:
            """The actual site symmetry's translation vector may vary from op by
            a factor of +1 or -1 (especially when op contains +-1/2).
            We record this to distinguish between special Wyckoff positions.
            As an example, consider the point (-x+1/2,-x,x+1/2) in position 16c
            of space group Ia-3(206). The site symmetry includes the operations
            (-z+1,x-1/2,-y+1/2) and (y+1/2,-z+1/2,-x+1). These operations are
            not listed in the general position, but correspond to the operations
            (-z,x+1/2,-y+1/2) and (y+1/2,-z+1/2,-x), respectively, just shifted
            by (+1,-1,0) and (0,0,+1), respectively.
            """
            el = SymmOp.from_rotation_and_translation(op.rotation_matrix, op.translation_vector - np.round(displacement))
            symmetry.append(el)
    return symmetry

def find_generating_point(coords, generators, PBC=[1,2,3]):
    """
    Given a set of coordinates and Wyckoff generators, return the coord which
    can be used to generate the others. This is useful for molecular Wyckoff
    positions, for which the orientation, and not just the position, is
    needed for each point in the Wyckoff position. Thus, we need to know which
    coordinates to use for x, y, and z, so that rotations can be applied
    correctly using the Wyckoff geneators

    Args:
        coords: a list of fractional coordinates corresponding to a Wyckoff
            position
        generators: the list of Wyckoff generators for the Wyckoff position.
            Can be obtained from get_wyckoff_generators
        PBC: a list of periodic axes (1,2,3)->(x,y,z)
    
    Returns:
        a fractional coordinate [x, y, z] corresponding to the first listed
        point in the Wyckoff position
     """
    for coord in coords:
        tmp_c = deepcopy(coords)
        tmp_c = filtered_coords(tmp_c, PBC=PBC)
        generated = list(gen.operate(coord) for gen in generators)
        generated = filtered_coords(generated, PBC=PBC)
        index_list1 = list(range(len(tmp_c)))
        index_list2 = list(range(len(generated)))
        if len(generated) != len(tmp_c):
            print("Warning: coordinate and generator lists have unequal length.")
            print("In check_wyckoff_position.find_generating_point:")
            print("len(coords): "+str(len(coords))+", len(generators): "+str(len(generators)))
            return None
        for index1, c1 in enumerate(tmp_c):
            for index2, c2 in enumerate(generated):
                if np.allclose(c1, c2, atol=.001, rtol=.001):
                    if index1 in index_list1 and index2 in index_list2:
                        index_list1.remove(index1)
                        index_list2.remove(index2)
                        break
        if index_list2 == []:
            return coord
    #If no valid coordinate is found
    return None

def check_wyckoff_position(points, wyckoffs, w_symm_all, exact_translation=False, PBC=[1,2,3]):
    """
    Given a list of points, returns a single index of a matching Wyckoff
    position in the space group. Checks the site symmetry of each supplied
    point against the site symmetry for each point in the Wyckoff position.

    Args:
        points: a list of 3d coordinates or SymmOps to check
        wyckoffs: an unorganized list of Wyckoff positions obtained from
            get_wyckoffs, get_layer, or get_rod
        w_symm_all: a list of site symmetry operations obtained from
            get_wyckoff_symmetry, get_layer_symmetry, or get_rod_symmetry
        exact_translation: whether we require two SymmOps to have exactly equal
            translational components. If false, translations related by +-1
            are considered equal. If points have been directly generated from
            a Wyckoff position, we may set this to True. Otherwise, leave False
        PBC: a list of periodic axes (1,2,3)->(x,y,z)

    Returns:
        a single index for the Wyckoff position within the sg. If no matching
        WP is found, returns False
    """
    points = np.array(points)
    gen_pos = wyckoffs[0]

    new_points = []
    if exact_translation == False:
        new_points = filtered_coords(points, PBC=PBC)
    else:
        new_points = deepcopy(points)

    p_symm = []
    #If exact_translation is false, store WP's which might be a match
    possible = []
    for x in new_points:
        p_symm.append(site_symm(x, gen_pos, PBC=PBC))
    
    for i, wp in enumerate(wyckoffs):
        w_symm = w_symm_all[i]
        if len(p_symm) == len(w_symm) and len(wp) == len(points):
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
                                if are_equal(op_w, op_p, PBC=PBC):
                                    temp2.remove(op_w)
                        if temp2 == []:
                            temp.remove(w)
            if temp == []:
                #If we find a match with exact translations
                if exact_translation:
                    generators = wyckoffs[i]
                    p = find_generating_point(points, generators, PBC=PBC)
                    if p is not None:
                        return i
                    else:
                        return False
                elif not exact_translation:
                    possible.append(i)
        #If no matching WP's are found
    if len(possible) == 0:
        return False
    #If exactly one matching WP is found
    elif len(possible) == 1:
        i = possible[0]
        generators = wyckoffs[i]
        p = find_generating_point(points, generators, PBC=PBC)
        if p is not None:
            return i
        else:
            return False
    #If multiple WP's are found
    else:
        #Check that points are generated from generators
        for i in possible:
            generators = wyckoffs[i]
            p = find_generating_point(points, generators, PBC=PBC)
            if p is not None:
                return i
        return False

def verify_distances(coordinates, species, lattice, factor=1.0, PBC=[1,2,3]):
    """
    Checks the inter-atomic distance between all pairs of atoms in a crystal.

    Args:
        coordinates: a 1x3 list of fractional coordinates
        species: a list of atomic symbols for each coordinate
        lattice: a 3x3 matrix representing the lattice vectors of the unit cell
        factor: a tolerance factor for checking distances. A larger value means
            atoms must be farther apart
        PBC: a list of periodic axes (1,2,3)->(x,y,z)
    
    Returns:
        True if no atoms are too close together, False if any pair is too close
    """
    for i, c1 in enumerate(coordinates):
        specie1 = species[i]
        for j, c2 in enumerate(coordinates):
            if j > i:
                specie2 = species[j]
                diff = np.array(c2) - np.array(c1)
                d_min = distance(diff, lattice, PBC=PBC)
                tol = factor*0.5*(Element(specie1).covalent_radius + Element(specie2).covalent_radius)
                if d_min < tol:
                    return False
    return True

class random_crystal():
    """
    Class for storing and generating atomic crystals based on symmetry
    constraints. Given a spacegroup, list of atomic symbols, the stoichiometry,
    and a volume factor, generates a random crystal consistent with the
    spacegroup's symmetry. This crystal is stored as a pymatgen struct via
    self.struct
    
    Args:
        sg: the international spacegroup number
        species: a list of atomic symbols for each ion type
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell)
        factor: a volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between atoms
    """
    def __init__(self, sg, species, numIons, factor):
        
        #Necessary input
        numIons = np.array(numIons) #must convert it to np.array
        self.factor = factor
        """The supplied volume factor for the unit cell."""
        self.numIons0 = numIons
        self.sg = sg
        """The international spacegroup number of the crystal."""
        self.species = species
        """A list of atomic symbols for the types of atoms in the crystal."""
        self.Msgs()
        """A list of warning messages to use during generation."""
        self.numIons = numIons * cellsize(self.sg)
        """The number of each type of atom in the CONVENTIONAL cell"""
        self.volume = estimate_volume(self.numIons, self.species, self.factor)
        """The volume of the generated unit cell"""
        self.wyckoffs = get_wyckoffs(self.sg)
        """The Wyckoff positions for the crystal's spacegroup."""
        self.wyckoffs_organized = get_wyckoffs(self.sg, organized=True)
        """The Wyckoff positions for the crystal's spacegroup. Sorted by
        multiplicity."""
        self.w_symm = get_wyckoff_symmetry(self.sg)
        """The site symmetry of the Wyckoff positions"""
        self.generate_crystal()


    def Msgs(self):
        """
        Define a set of error and warning message if generation fails.

        Returns:
            nothing
        """
        self.Msg1 = 'Error: the number is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'

    def check_compatible(self):
        """
        Checks if the number of atoms is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        N_site = [len(x[0]) for x in self.wyckoffs_organized]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for numIon in self.numIons:
            #Check that the number of ions is a multiple of the smallest Wyckoff position
            if numIon % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs_organized[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numIon
                    for x in self.wyckoffs_organized:
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
            #Wyckoff Positions have no degrees of freedom
            return 0

    def generate_crystal(self, max1=max1, max2=max2, max3=max3):
        """
        The main code to generate a random atomic crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
        """
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
                if cell_para is None:
                    break
                else:
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
                                ops = choose_wyckoff(self.wyckoffs_organized, numIon-numIon_added) 
                                if ops is not False:
            	        	    #Generate a list of coords from ops
                                    point = np.random.random(3)
                                    coords = np.array([op.operate(point) for op in ops])
                                    #Merge coordinates if the atoms are close
                                    coords_toadd, good_merge = merge_coordinate(coords, cell_matrix, self.wyckoffs, self.w_symm, tol)
                                    if good_merge is not False:
                                        coords_toadd -= np.floor(coords_toadd) #scale the coordinates to [0,1], very important!
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
                        """A 3x3 matrix representing the lattice of the unit
                        cell."""                 
                        self.coordinates = np.array(final_coor)
                        """The fractional coordinates for each molecule in the
                        final structure"""
                        self.sites = final_site
                        """A list of atomic symbols corresponding to the type
                        of atom for each site in self.coordinates"""
                        self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                        """A pymatgen.core.structure.Structure object for the
                        final generated crystal."""
                        self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                        """A list of information describing the generated
                        crystal, which may be used by spglib for symmetry
                        analysis."""
                        self.valid = True
                        """Whether or not a valid crystal was generated."""
                        return
        if degrees == 0: print("Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2

class random_crystal_2D():
    """
    A 2d counterpart to random_crystal. Generates a random atomic crystal based
    on a 2d layer group instead of a 3d spacegroup. Note that each layer group
    is equal to a corresponding 3d spacegroup, but without periodicity in one
    direction. The generated pymatgen structure can be accessed via self.struct

    Args:
        number: the layer group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
        species: a list of atomic symbols for each ion type
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell)
        thickness: the thickness, in Angstroms, of the unit cell in the 3rd
            dimension (the direction which is not repeated periodically)
        factor: a volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between atoms
    """
    def __init__(self, number, species, numIons, thickness, factor):

        self.number = number
        self.lgp = Layergroup(number)
        """The number (between 1 and 80) for the crystal's layer group."""
        self.sg = self.lgp.sgnumber
        """The number (between 1 and 230) for the international spacegroup."""
        numIons = np.array(numIons) #must convert it to np.array
        self.factor = factor
        """"The volume factor used to generate the unit cell."""
        self.thickness = thickness
        """the thickness, in Angstroms, of the unit cell in the 3rd
        dimension."""
        self.numIons0 = numIons
        self.species = species
        """A list of atomic symbols for the types of atoms in the crystal."""
        self.PBC = [1,2]
        """The periodic axes of the crystal."""
        self.PB = self.lgp.permutation[3:6] 
        #TODO: add docstring
        self.P = self.lgp.permutation[:3] 
        #TODO: add docstring
        self.Msgs()
        """A list of warning messages to use during generation."""
        self.numIons = numIons * cellsize(self.sg)
        """The number of each type of atom in the CONVENTIONAL cell"""
        self.volume = estimate_volume(self.numIons, self.species, self.factor)
        """The volume of the generated unit cell"""
        self.wyckoffs = get_layer(self.number)
        """The Wyckoff positions for the crystal's spacegroup."""      
        self.wyckoffs_organized = get_layer(self.number, organized=True)
        """The Wyckoff positions for the crystal's spacegroup. Sorted by
        multiplicity."""
        self.w_symm = get_layer_symmetry(self.number)
        """A list of site symmetry operations for the Wyckoff positions, obtained
            from get_wyckoff_symmetry."""
        self.generate_crystal()


    def Msgs(self):
        """
        Define a set of error and warning message if generation fails.
        """
        self.Msg1 = 'Error: the number is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'

    def check_compatible(self):
        """
        Checks if the number of atoms is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        N_site = [len(x[0]) for x in self.wyckoffs_organized]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for numIon in self.numIons:
            #Check that the number of ions is a multiple of the smallest Wyckoff position
            if numIon % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs_organized[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numIon
                    for x in self.wyckoffs_organized:
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
            #Wyckoff Positions have no degrees of freedom
            return 0

    def generate_crystal(self, max1=max1, max2=max2, max3=max3):
        """
        The main code to generate a random atomic crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
        """
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
                cell_para = generate_lattice_2D(self.number, self.volume, thickness=self.thickness, minvec=minvector)
                if cell_para is None:
                    break
                cell_matrix = para2matrix(cell_para)
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
                            ops = choose_wyckoff(self.wyckoffs_organized, numIon-numIon_added)
                            if ops is not False:
                	    	    #Generate a list of coords from ops
                                point = np.random.random(3)
                                for a in range(1, 4):
                                    if a not in self.PBC:
                                        point[a-1] -= 0.5
                                coords = np.array([op.operate(point) for op in ops])
                                coords_toadd, good_merge = merge_coordinate(coords, cell_matrix, self.wyckoffs, self.w_symm, tol, PBC=self.PBC)
                                if good_merge is not False:
                                    coords_toadd = filtered_coords(coords_toadd, PBC=self.PBC) #scale the coordinates to [0,1], very important!
                                    if check_distance(coordinates_tmp, coords_toadd, sites_tmp, specie, cell_matrix, PBC=self.PBC):
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
                    final_coor = np.array(final_coor)
                    #final_lattice, final_coor = Permutation(final_lattice, final_coor, self.PB)
                    final_lattice, final_coor = Add_vacuum(final_lattice, final_coor)
                    self.lattice = final_lattice
                    """A 3x3 matrix representing the lattice of the unit
                    cell."""                        
                    self.coordinates = final_coor
                    """The fractional coordinates for each molecule in the
                    final structure"""
                    self.sites = final_site  
                    """A list of atomic symbols corresponding to the type
                    of atom for each site in self.coordinates"""                  
                    self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                    """A pymatgen.core.structure.Structure object for the
                    final generated crystal."""                    
                    self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                    """A list of information describing the generated
                    crystal, which may be used by spglib for symmetry
                    analysis."""                    
                    self.valid = True
                    """Whether or not a valid crystal was generated."""
                    return
        if degrees == 0: print("Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2

class random_crystal_1D():
    """
    A 1d counterpart to random_crystal. Generates a random atomic crystal based
    on a 1d Rod group instead of a 3d spacegroup. Note that each layer group
    is equal to a corresponding 3d spacegroup, but without periodicity in one
    direction. The generated pymatgen structure can be accessed via self.struct

    Args:
        number: the Rod group number between 1 and 75. NOT equal to the
            international space group number, which is between 1 and 230
        species: a list of atomic symbols for each ion type
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell)
        area: the effective cross-sectional area, in Angstroms squared, of the
            unit cell
        factor: a volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between atoms
    """
    def __init__(self, number, species, numIons, area, factor):

        self.number = number
        numIons = np.array(numIons) #must convert it to np.array
        self.factor = factor
        """"The volume factor used to generate the unit cell."""
        self.area = area
        """the effective cross-sectional area, in Angstroms squared, of the
        unit cell."""
        self.numIons0 = numIons
        self.species = species
        """A list of atomic symbols for the types of atoms in the crystal."""
        self.PBC = [3]
        """The periodic axis of the crystal."""
        self.Msgs()
        """A list of warning messages to use during generation."""
        self.numIons = numIons #cellsize always == 1
        """The number of each type of atom in the CONVENTIONAL cell"""
        self.volume = estimate_volume(self.numIons, self.species, self.factor)
        """The volume of the generated unit cell"""
        self.wyckoffs = get_rod(self.number)
        """The Wyckoff positions for the crystal's spacegroup."""      
        self.wyckoffs_organized = get_rod(self.number, organized=True)
        """The Wyckoff positions for the crystal's spacegroup. Sorted by
        multiplicity."""
        self.w_symm = get_rod_symmetry(self.number)
        """A list of site symmetry operations for the Wyckoff positions, obtained
            from get_rod_symmetry."""
        self.generate_crystal()


    def Msgs(self):
        """
        Define a set of error and warning message if generation fails.
        """
        self.Msg1 = 'Error: the number is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'

    def check_compatible(self):
        """
        Checks if the number of atoms is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        N_site = [len(x[0]) for x in self.wyckoffs_organized]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for numIon in self.numIons:
            #Check that the number of ions is a multiple of the smallest Wyckoff position
            if numIon % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs_organized[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numIon
                    for x in self.wyckoffs_organized:
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
            #Wyckoff Positions have no degrees of freedom
            return 0

    def generate_crystal(self, max1=max1, max2=max2, max3=max3):
        """
        The main code to generate a random atomic crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
        """
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
                cell_para = generate_lattice_1D(self.number, self.volume, area=self.area, minvec=minvector)
                if cell_para is None:
                    break
                cell_matrix = para2matrix(cell_para)
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
                            ops = choose_wyckoff(self.wyckoffs_organized, numIon-numIon_added)
                            if ops is not False:
                	    	    #Generate a list of coords from ops
                                point = np.random.random(3)
                                for a in range(1, 4):
                                    if a not in self.PBC:
                                        if self.number < 46:
                                            point[a-1] -= 0.5
                                        elif self.number >= 46:
                                            point[a-1] *= 1./sqrt(3.)
                                coords = np.array([op.operate(point) for op in ops])
                                coords_toadd, good_merge = merge_coordinate(coords, cell_matrix, self.wyckoffs, self.w_symm, tol, PBC=self.PBC)
                                if good_merge is not False:
                                    coords_toadd = filtered_coords(coords_toadd, PBC=self.PBC) #scale the coordinates to [0,1], very important!
                                    if check_distance(coordinates_tmp, coords_toadd, sites_tmp, specie, cell_matrix, PBC=self.PBC):
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
                    final_coor = np.array(final_coor)
                    #TODO: Implement Add_vacuum for 1D lattices
                    self.lattice = final_lattice
                    """A 3x3 matrix representing the lattice of the unit
                    cell."""                        
                    self.coordinates = final_coor
                    """The fractional coordinates for each molecule in the
                    final structure"""
                    self.sites = final_site  
                    """A list of atomic symbols corresponding to the type
                    of atom for each site in self.coordinates"""                  
                    self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                    """A pymatgen.core.structure.Structure object for the
                    final generated crystal."""                    
                    self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                    """A list of information describing the generated
                    crystal, which may be used by spglib for symmetry
                    analysis."""                    
                    self.valid = True
                    """Whether or not a valid crystal was generated."""
                    return
        if degrees == 0: print("Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2


if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    import os
    parser = OptionParser()
    parser.add_option("-s", "--spacegroup", dest="sg", metavar='sg', default=36, type=int,
            help="desired space group number (1-230) or layer group number (1-80), e.g., 36")
    parser.add_option("-e", "--element", dest="element", default='Li', 
            help="desired elements: e.g., Li", metavar="element")
    parser.add_option("-n", "--numIons", dest="numIons", default=16, 
            help="desired numbers of atoms: 16", metavar="numIons")
    parser.add_option("-f", "--factor", dest="factor", default=3.0, type=float, 
            help="volume factor: default 3.0", metavar="factor")
    parser.add_option("-v", "--verbosity", dest="verbosity", default=0, type=int, 
            help="verbosity: default 0; higher values print more information", metavar="verbosity")
    parser.add_option("-a", "--attempts", dest="attempts", default=1, type=int, 
            help="number of crystals to generate: default 1", metavar="attempts")
    parser.add_option("-o", "--outdir", dest="outdir", default="out", type=str, 
            help="Directory for storing output cif files: default 'out'", metavar="outdir")
    parser.add_option("-d", "--dimension", dest="dimension", metavar='dimension', default=3, type=int,
            help="desired dimension: (3, 2, or 1 for 3d, 2d, or 1D respectively)")
    parser.add_option("-t", "--thickness", dest="thickness", metavar='thickness', default=None, type=float,
            help="Thickness, in Angstroms, of a 2D crystal: default 2.0")

    (options, args) = parser.parse_args()
    sg = options.sg
    dimension = options.dimension
    if dimension == 3:
        if sg < 1 or sg > 230:
            print("Invalid space group number. Must be between 1 and 230.")
            sys.exit(0)
    elif dimension == 2:
        if sg < 1 or sg > 80:
            print("Invalid layer group number. Must be between 1 and 80.")
            sys.exit(0)
    elif dimension == 1:
        if sg < 1 or sg > 75:
            print("Invalid Rod group number. Must be between 1 and 75.")
            sys.exit(0)
    elif dimension == 0:
        print("0d clusters cannot currently be generated. Use dimension 1, 2, or 3.")
        sys.exit(0)
    else:
        print("Invalid dimension. Use dimension 1, 2, or 3.")
        sys.exit(0)

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

    factor = options.factor
    if factor < 0:
        print("Error: Volume factor must be greater than 0.")
        sys.exit(0)

    verbosity = options.verbosity
    attempts = options.attempts
    outdir = options.outdir
    dimension = options.dimension
    thickness = options.thickness

    try:
        os.mkdir(outdir)
    except: pass

    filecount = 1 #To check whether a file already exists
    for i in range(attempts):
        numIons0 = np.array(numIons)
        sg = options.sg
        start = time()
        if dimension == 3:
            rand_crystal = random_crystal(options.sg, system, numIons0, factor)
            sg1 = sg
        elif dimension == 2:
            rand_crystal = random_crystal_2D(options.sg, system, numIons0, thickness, factor)
            sg1 = rand_crystal.sg
        elif dimension == 1:
            rand_crystal = random_crystal_1D(options.sg, system, numIons0, thickness, factor)
            sg1 = "?"
        end = time()
        timespent = np.around((end - start), decimals=2)

        if rand_crystal.valid:
            #Output a cif file
            written = False
            try:
                comp = str(rand_crystal.struct.composition)
                comp = comp.replace(" ", "")
                cifpath = outdir + '/' + comp + "_" + str(filecount) + '.cif'
                while os.path.isfile(cifpath):
                    filecount += 1
                    cifpath = outdir + '/' + comp + "_" + str(filecount) + '.cif'
                CifWriter(rand_crystal.struct, symprec=0.1).write_file(filename = cifpath)
                written = True
            except: pass
            #POSCAR output
            #rand_crystal.struct.to(fmt="poscar", filename = '1.vasp')

            #spglib style structure called cell
            ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)['number']
            print('Space group  requested:', sg1, ' generated:', ans)
            if written is True:
                print("    Output to "+cifpath)
            else:
                print("    Could not write cif file.")

            #Print additional information about the structure
            if verbosity > 0:
                print("Time required for generation: " + str(timespent) + "s")
                print(rand_crystal.struct)


        #If generation fails
        else: 
            print('something is wrong')
            print('Time spent during generation attempt: ' + str(timespent) + "s")
