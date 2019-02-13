from pyxtal.crystal import random_cluster
from random import randint
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from pyxtal.molecule import PointGroupAnalyzer
from pymatgen import Molecule
from pyxtal.database.collection import Collection
from time import time
import numpy as np
import warnings
warnings.filterwarnings("ignore")

"""
This is a script to 
1, generate random clusters
2, perform optimization
"""
def LJ(pos, dim, mu=0.1):
    """
    Calculate the total energy
    Args:
    pos: 1D array with N*dim numbers representing the atomic positions
    dim: dimension of the hyper/normal space
    output
    E: the total energy with punishing function
    """
    N_atom = int(len(pos)/dim)
    pos = np.reshape(pos, (N_atom, dim))
    
    distance = cdist(pos, pos, 'euclidean')
    ids = np.triu_indices(N_atom)
    distance = distance[ids]
    distance = distance[distance > 0]
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    Eng = np.sum(4*(1/r12 - 1/r6))
    if dim > 3:
        pos0 = pos
        norm = 0
        for i in range(3,dim):
            diff = pos[:, i] - np.mean(pos[:, i])
            norm += np.linalg.norm(diff)
        Eng += mu*norm
    return Eng
            
def single_optimize(pos, dim=3, kt=0.5, mu=0.1):
    """
    perform optimization for a given cluster
    Args: 
    pos: N*dim0 array representing the atomic positions
    dim: dimension of the hyper/normal space
    kt: perturbation factors

    output:
    energy: optmized energy
    pos: optimized positions
    """
    N_atom = len(pos)
    diff = dim - np.shape(pos)[1]
    # if the input pos has less dimensions, we insert a random array for the extra dimension
    # if the input pos has more dimensions, we delete the array for the extra dimension
    if diff > 0:
        pos = np.hstack((pos, 0.5*(np.random.random([N_atom, diff])-0.5) ))
    elif diff < 0:
        pos = pos[:, :dim]

    pos = pos.flatten()
    res = minimize(LJ, pos, args=(dim, mu), method='CG', tol=1e-4)
    pos = np.reshape(res.x, (N_atom, dim))
    energy = res.fun
    return energy, pos


def parse_symmetry(pos):
    mol = Molecule(['C']*len(pos), pos)
    try:
        symbol = PointGroupAnalyzer(mol).sch_symbol
    except:
        symbol = 'N/A'
    return symbol


class LJ_prediction():
    """
    A class to perform global optimization on LJ clusters
    Args:

    Attributes:

    """

    def __init__(self, numIons):
        self.numIons = numIons
        ref = Collection('clusters')[str(numIons)]
        print('\nReference for LJ {0:3d} is {1:12.3f} eV, PG: {2:4s}'.\
                format(numIons, ref['energy'], ref['pointgroup']))
        self.reference = ref

    def generate_cluster(self):
        run = True
        while run:
            pg = randint(2, 32)
            cluster = random_cluster(pg, ['Mo'], [self.numIons], 1.0)
            if cluster.valid:
                run = False
        return cluster.coordinates
 
    def predict(self, dim=3, maxN=100, ncpu=2):

        print('\nPerforming random search at {0:d}D space\n'.format(dim))
        cycle = range(maxN)
        if ncpu > 1:
            from multiprocessing import Pool
            from functools import partial

            with Pool(ncpu) as p:
                func = partial(self.relaxation, dim)
                res = p.map(func, cycle)
                p.close()
                p.join()
        else:
            res=[]
            for i in cycle:
                res.append(self.relaxation(dim, i))
        
        N_success = 0
        for dct in res:
            if dct['ground']:
                N_success +=1
        print('\nHit the ground state {0:4d} times out of {1:4d} attempts\n'.\
                format(N_success, maxN))

    def relaxation(self, dim, ind):
        pos = self.generate_cluster()
        pg1 = parse_symmetry(pos)
        if dim == 3:
            [energy, pos] = single_optimize(pos, 3)
        else:
            do = True
            while do:
                [energy1, pos1] = single_optimize(pos, 3)
                [energy2, pos2] = single_optimize(pos1, dim)
                [energy3, pos3] = single_optimize(pos2, 3)
                #print(energy1, energy2, energy3)
                if abs(energy3-energy1) < 1e-3 or energy3 > energy1:
                    pos = pos1
                    energy = energy1
                    do = False
                    #print('stop')
                else:
                    pos = pos3
        if abs(energy-self.reference['energy']) <1e-3:
            ground = True
        else:
            ground = False

        pg2 = parse_symmetry(pos)
        res = {'pos': pos,
               'energy': energy,
               'pg_init': pg1,
               'pg_finial': pg2,
               'ground': ground,
               'id': ind,
               }
        if ground:
            print('ID: {0:4d} PG initial: {1:4s} relaxed: {2:4s} Energy: {3:12.3f} ++++++'.\
                   format(ind, pg1, pg2, energy))
        elif ind%10 == 0:
            print('ID: {0:4d} PG initial: {1:4s} relaxed: {2:4s} Energy: {3:12.3f}'.\
                   format(ind, pg1, pg2, energy))
 
        return res

if __name__ == "__main__":

    lj_run = LJ_prediction(13)
    t0 = time()
    lj_run.predict(dim=4, maxN=10, ncpu=2)
    print('time: {0:6.2f} seconds'.format(time()-t0))
    lj_run.predict(dim=4, maxN=10, ncpu=2)
    print('time: {0:6.2f} seconds'.format(time()-t0))
    #lj_run.predict(4, 10)
