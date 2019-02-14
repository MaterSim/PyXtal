from pyxtal.crystal import random_cluster
from random import randint
from time import time
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from pyxtal.molecule import PointGroupAnalyzer
from pymatgen import Molecule
from pyxtal.database.collection import Collection
import numpy as np
import warnings
warnings.filterwarnings("ignore")

pyxtal_verbosity = 0

"""
This is a script to 
1, generate random clusters
2, perform optimization
"""
def LJ(pos):
    """
    Calculate the total energy
    input:
    pos: N*3 array which represents the atomic positions
    output
    E: the total energy
    """
    N_atom = int(len(pos)/3)
    pos = np.reshape(pos, (N_atom, 3))
    distance = cdist(pos, pos, 'euclidean')
    ids = np.triu_indices(N_atom)
    distance = distance[ids]
    distance = distance[distance > 0]
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    return np.sum(4*(1/r12 - 1/r6))
            
def single_optimize(pos):
    N_atom = int(len(pos)/3)
    res = minimize(LJ, pos, method='CG', tol=1e-4)
    pos = np.reshape(res.x, (N_atom, 3))
    energy = res.fun
    return energy, pos

def parse_symmetry(pos):
    mol = Molecule(['C']*len(pos), pos)
    try:
        symbol = PointGroupAnalyzer(mol).sch_symbol
    except:
        symbol = 'N/A'
    return symbol

N_attempts = 100
numIons = 12
factor = 1.1
ref = Collection('clusters')[str(numIons)]
print('The reference energy for LJ {0:3d} is {1:12.3f}, pointgroup: {2:4s}'.format(numIons, ref['energy'], ref['pointgroup']))

N_success = 0
t0 = time()
for i in range(N_attempts):
    run = True
    while run:
        pg = randint(2, 32)
        cluster = random_cluster(pg, ['Mo'], [numIons], factor)
        if cluster.valid:
            run = False

    pg1 = parse_symmetry(cluster.coordinates)
    pos = cluster.coordinates.flatten()
    [energy, pos] = single_optimize(pos)
    pg2 = parse_symmetry(pos)
    if abs(energy-ref['energy']) <1e-3:
        N_success += 1
        print('PG requested: {0:4s} relaxed: {1:4s} Energy: {2:12.3f} Time: {3:6.1f} mins ++++++'.format(pg1, pg2, energy, (time()- t0)/60.0))
    else:
        print('PG requested: {0:4s} relaxed: {1:4s} Energy: {2:12.3f} Time: {3:6.1f} mins'.format(pg1, pg2, energy, (time()- t0)/60.0))

print('Hit the ground state {0:4d} times out of {1:4d} attempts'.format(N_success, N_attempts))
