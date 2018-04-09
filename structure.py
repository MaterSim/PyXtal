'''Program for generation of random crystal structures.
by Scott Fredericks, Spring 2018
Given a space group number between 1 and 230,
and a number N of atoms in the primitive cell,
produces a crystal structure with random atomic coordinates.
Outputs a cif file with conventional setting'''
import spglib
from vasp import read_vasp
from pymatgen.symmetry.groups import sg_symbol_from_int_number
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from optparse import OptionParser
from scipy.spatial.distance import cdist
import numpy as np
from os.path import isfile
from random import uniform as rand
from random import choice as choose
from random import randint
from math import sqrt, pi, sin, cos, acos, fabs
import database.make_sitesym as make_sitesym
import database.hall as hall
from database.element import Element
from copy import deepcopy
#Define variables
#------------------------------
tol_m = 1.0 #seperation tolerance in Angstroms
max1 = 30 #Attempts for generating lattices
max2 = 30 #Attempts for a given lattice
max3 = 30 #Attempts for a given Wyckoff position
minvec = 2.0 #minimum vector length
ang_min = 30
ang_max = 150
#Define functions
#------------------------------
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
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                np.append(coord2, coord2+[i,j,k])
    
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
    for atom1 in range(1,len(xyzs)):
        dist_min = 10.0
        for atom2 in range(0, atom1):
            #shift atom1 to position close to atom2
            #print(np.round(xyzs[atom1] - xyzs[atom2]))
            xyzs[atom2] += np.round(xyzs[atom1] - xyzs[atom2])
            #print(xyzs[atom1] - xyzs[atom2])
            matrix = matrix0 + (xyzs[atom1] - xyzs[atom2])
            matrix = np.dot(matrix, lattice)
            dists = cdist(matrix, [[0,0,0]])
            if np.min(dists) < dist_min:
                dist_min = np.min(dists)
                #print(dist_min, matrix[np.argmin(dists)]/4.72)
                matrix_min = matrix0[np.argmin(dists)]
        #print(atom1, xyzs[atom1], matrix_min, dist_min)
        xyzs[atom1] += matrix_min
    return xyzs.mean(0)


#generate random coordinate
def rand_coords(xyz_string): 
    """
    args:
    xyz_string: 0, y, 1/4
    return: random numbers for the places where x, y, z is present
    """
    #xyz_string = ops[0].as_xyz_string()
    xyz = []
    x,y,z = np.random.random(3)
    for content in xyz_string.strip('()').split(','):
        if content.find('x')>=0:
            xyz.append(x)
        elif content.find('y')>=0:
            xyz.append(y)
        elif content.find('z')>=0:
            xyz.append(z)
        elif content.find('/')>=0:
            tmp = content.split('/')
            xyz.append(float(tmp[0])/float(tmp[1]))
        else:
            xyz.append(float(content))

    return xyz

def para2matrix(cell_para):
    """ 1x6 (a, b, c, alpha, beta, gamma) -> 3x3 representation -> """
    matrix = np.zeros([3,3])
    matrix[0][0] = cell_para[0]
    matrix[1][0] = cell_para[1]*cos(cell_para[5])
    matrix[1][1] = cell_para[1]*sin(cell_para[5])
    matrix[2][0] = cell_para[2]*cos(cell_para[4])
    matrix[2][1] = cell_para[2]*cos(cell_para[3])*sin(cell_para[4])
    matrix[2][2] = sqrt(cell_para[2]**2 - matrix[2][0]**2 - matrix[2][1]**2)
    
    return matrix

def matrix2para(matrix):
    """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
    cell_para = np.zeros(6)
    cell_para[0] = np.linalg.norm(matrix[0])
    cell_para[1] = np.linalg.norm(matrix[1])
    cell_para[2] = np.linalg.norm(matrix[2])

    cell_para[5] = angle(matrix[0], matrix[1])
    cell_para[4] = angle(matrix[0], matrix[2])
    cell_para[3] = angle(matrix[1], matrix[2])

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

def merge_coordinate(coor, lattice, wyckoff, tol):
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
                #print('Merging----------------')
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

def generate_lattice(sg, volume):
    """
    generate the lattice according to the space group symmetry and number of atoms
    if the space group has centering, we will transform to conventional cell setting
    """
    minvec2 = minvec*minvec
    alpha = np.radians(rand(ang_min, ang_max))
    beta  = np.radians(rand(ang_min, ang_max))
    gamma = np.radians(rand(ang_min, ang_max))
    #Triclinic
    if sg <= 2:
        x = sqrt(fabs(1. - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2.*cos(alpha)*cos(beta)*cos(gamma)))
        volume = volume/x
        a = rand(minvec, volume/minvec2)
        b = rand(minvec, volume/(minvec*a))
        c = volume/(a*b)
    #Monoclinic
    elif sg <= 15:
        alpha, gamma  = pi/2, pi/2
        x = sin(beta)
        volume = volume/x
        a = rand(minvec, volume/(minvec2))
        b = rand(minvec, volume/(minvec*a))
        c = volume/(a*b)
    #Orthorhombic
    elif sg <= 74:
        alpha, beta, gamma = pi/2, pi/2, pi/2
        a = rand(minvec, volume/minvec2)
        b = rand(minvec, volume/(minvec*a))
        c = volume/(a*b)
    #Tetragonal
    elif sg <= 142:
        alpha, beta, gamma = pi/2, pi/2, pi/2
        c = rand(minvec, volume/minvec2)
        a = sqrt(volume/c)
        b = a
    #Trigonal/Rhombohedral/Hexagonal
    elif sg <= 194:
        alpha, beta, gamma = pi/2, pi/2, pi/3*2
        x = sqrt(3.)/2.
        volume = volume/x
        c = rand(minvec, volume/minvec2)
        a = sqrt(volume/c)
        b = a
    #Cubic
    else:
        alpha, beta, gamma = pi/2, pi/2, pi/2
        s = (volume) ** (1./3.)
        a, b, c = s, s, s
    return np.array([a, b, c, alpha, beta, gamma])

def check_compatible(numIons, wyckoff):
    """
    check if the number of atoms is compatible with the wyckoff positions
    needs to improve later
    """
    N_site = [len(x[0]) for x in wyckoff]
    for numIon in numIons:
        if numIon % N_site[-1] > 0:
            return False
    return True

def filter_site(v): #needs to explain
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
    for wyckoff in wyckoffs:
        if len(wyckoff[0]) <= number:
            return choose(wyckoff)
    return False

def get_wyckoff_positions(sg):
    """find the wyckoff positions
    Args:
        sg: space group number (19)
    Return: a list containg the operation matrix sorted by multiplicity
    4a
        [Rot:
        [[ 1.  0.  0.]
        [ 0.  1.  0.]
        [ 0.  0.  1.]]
        tau
        [ 0.  0.  0.], 
        Rot:
        [[-1.  0.  0.]
         [ 0. -1.  0.]
         [ 0.  0.  1.]]
        tau
        [ 0.5  0.   0.5], 
        Rot:
        [[-1.  0.  0.]
         [ 0.  1.  0.]
         [ 0.  0. -1.]]
        tau
        [ 0.   0.5  0.5], 
        Rot:
        [[ 1.  0.  0.]
         [ 0. -1.  0.]
         [ 0.  0. -1.]]
        tau
        [ 0.5  0.5  0. ]]]
    """
    array = []
    hall_number = hall.hall_from_hm(sg)
    wyckoff_positions = make_sitesym.get_wyckoff_position_operators('database/Wyckoff.csv', hall_number)
    for x in wyckoff_positions:
    	temp = []
    	for y in x:
    	    temp.append(SymmOp.from_rotation_and_translation(list(y[0]), filter_site(y[1]/24)))
    	array.append(temp)
        
    i = 0
    wyckoffs_organized = [[]] #2D Array of Wyckoff positions organized by multiplicity
    old = len(array[0])
    for x in array:
        mult = len(x)
        if mult != old:
            wyckoffs_organized.append([])
            i += 1
            old = mult
        wyckoffs_organized[i].append(x)
    return wyckoffs_organized

class generate_crystal(Structure):
    def __init__(self, sg, species, numIons, factor):
        numIons *= cellsize(sg)
        volume = estimate_volume(numIons, species, factor)
        wyckoffs = get_wyckoff_positions(sg) #2D Array of Wyckoff positions organized by multiplicity
        Msg1 = 'Error: the number is incompatible with the wyckoff sites choice'
        Msg2 = 'Error: failed in the cycle of generating structures'
        Msg3 = 'Warning: failed in the cycle of adding species'
        Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        Msg5 = 'Finishing: added the specie'
        Msg6 = 'Finishing: added the whole structure'

        if check_compatible(numIons, wyckoffs) is False:
            print(Msg1)
        else:
            for cycle1 in range(max1):
                #1, Generate a lattice
                cell_para = generate_lattice(sg, volume)
                cell_matrix = para2matrix(cell_para)
                coordinates_total = [] #to store the added coordinates
                sites_total = []      #to store the corresponding specie
                good_structure = False

                for cycle2 in range(max2):
                    coordinates_tmp = deepcopy(coordinates_total)
                    sites_tmp = deepcopy(sites_total)
                    
            	#Add specie by specie
                    for numIon, specie in zip(numIons, species):
                        numIon_added = 0
                        tol = max(0.5*Element(specie).covalent_radius, tol_m)
                        #Now we start to add the specie to the wyckoff position
                        #print(wyckoffs[-1][0][0].as_xyz_string)
                        for cycle3 in range(max3):

                            #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                            ops = choose_wyckoff(wyckoffs, numIon-numIon_added) 
                            if ops is not False:
            	    	    #Generate a list of coords from ops
                                #point = rand_coords()  #ops[0].as_xyz_string()) 
                                point = np.random.random(3)
                                #print('generating new points:', point)
                                coords = np.array([op.operate(point) for op in ops])
                                #merge_coordinate if the atoms are close
                                coords_toadd, good_merge = merge_coordinate(coords, cell_matrix, wyckoffs, tol)
                                if good_merge:
                                    coords_toadd -= np.floor(coords_toadd) #scale the coordinates to [0,1], very important!
                                    #print('existing: ', coordinates_tmp)
                                    if check_distance(coordinates_tmp, coords_toadd, sites_tmp, specie, cell_matrix):
                                        coordinates_tmp.append(coords_toadd)
                                        sites_tmp.append(specie)
                                        numIon_added += len(coords_toadd)
                                    if numIon_added == numIon:
                                        #print(Msg5)
                                        coordinates_total = deepcopy(coordinates_tmp)
                                        sites_total = deepcopy(sites_tmp)
                                        break

                    if numIon_added == numIon:
                        print(Msg6)
                        good_structure = True
                        break
                    elif cycle2+1 == max2:
                        #print(coordinates_total)
                        print(Msg3)
                if good_structure:
                    final_coor = []
                    final_site = []
                    final_lattice = cell_matrix
                    for coor, ele in zip(coordinates_total, sites_total):
                        for x in coor:
                            final_coor.append(x)
                            final_site.append(ele)

                    if len(final_coor) > 48:
                        self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                        self.good_struct = False
                        return
                    self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                    self.good_struct = True
                    return

        self.struct = Msg2
        self.good_struct = False
        return

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
        sg = randint(2,230)
        #new_struct, good_struc = generate_crystal(options.sg, system, numIons0, options.factor)
        rand_crystal = generate_crystal(sg, system, numIons0, options.factor)
        new_struct, good_struc = rand_crystal.struct, rand_crystal.good_struct
        if good_struc:
            new_struct.to(fmt="poscar", filename = '1.vasp')
            cell = read_vasp('1.vasp')
            #ans = spglib.get_spacegroup(cell)
            ans = spglib.get_symmetry_dataset(cell, symprec=1e-1)['number']
            print('Space group  requested: ', sg, 'generated', ans)
            if ans < int(sg/1.2):
                print('something is wrong')
                break
            #print(CifWriter(new_struct, symprec=0.1).__str__())
            #print('Space group:', finder.get_space_group_symbol(), 'tolerance:', tol)
            #output wyckoff sites only

        else: 
            print('something is wrong')
            #print(len(new_struct.frac_coords))
            break
            #print(new_struct)
