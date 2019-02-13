from pyxtal.crystal import random_cluster
from random import randint
from time import time
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
import numpy as np
import requests
"""
This is a script to 
1, generate random clusters
2, perform optimization
"""
def get_pos_from_url(N=7, address='http://doye.chem.ox.ac.uk/jon/structures/LJ/points/'):
    url_address = address + str(N)
    data_str = requests.get(url_address).text
    return parse_url_text(data_str)

def parse_url_text(data_str):
    x_array = []
    text = data_str.split('\n')
    for line in text:
        [x_array.append(float(i)) for i in line.split()]
    return np.array(x_array)

def reference(N=7):
    pos = get_pos_from_url(N)
    energy = LJ(pos)
    pos = np.reshape(pos, (N, 3))
    return energy, pos

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


N_attempts = 100
numIons = 13 
factor = 1.1
ref, pos = reference(numIons)
print('The reference energy for LJ {0:3d} is {1:12.3f}'.format(numIons, ref))

N_success = 0
t0 = time()
for i in range(N_attempts):
    run = True
    while run:
        pg = randint(2, 32)
        cluster = random_cluster(pg, ['Mo'], [numIons], factor)
        if cluster.valid:
            run = False

    pos = cluster.coordinates.flatten()
    [energy, pos] = single_optimize(pos)
    print('PG requested: {0:3d} Energy: {1:12.3f} Time: {2:6.1f} mins'.\
            format(pg, energy, (time()- t0)/60.0))
    if abs(energy-ref) <1e-3:
        N_success += 1

print('Hit the ground state {0:4d} times out of {1:4d} attempts'.format(N_success, N_attempts))
