from scipy.optimize import minimize
from scipy.spatial.distance import pdist, cdist
import numpy as np
from time import time

"""
This is a script to 
1, generate random clusters
2, perform optimization
"""
def LJ(pos, dim=3):
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
    
    distance = pdist(pos) 
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    Eng = np.sum(4*(1/r12 - 1/r6))

    return Eng

def LJ_force(pos, dim=3):
    """
    Calculate the total energy
    Args:
        pos: 1D array with N*dim numbers representing the atomic positions
        dim: dimension of the hyper/normal space

    output
        F: 1D force array
    """
    N_atom = int(len(pos)/dim)
    pos = np.reshape(pos,[N_atom, dim])
    force = np.zeros([N_atom, dim])
    for i, pos0 in enumerate(pos):
        pos1 = pos.copy()
        pos1 = np.delete(pos1, i, 0)
        distance = cdist([pos0], pos1)
        r = pos1 - pos0
        r2 = np.power(distance, 2)
        r6 = np.power(r2, 3)
        r12 = np.power(r6, 2)
        force[i] = np.dot((48/r12-24/r6)/r2, r)
    return force.flatten()

def single_optimize(pos):
    """
    perform optimization for a given cluster

    Args: 
        pos: N*dim0 array representing the atomic positions
        dim: dimension of the hyper/normal space

    output:
        energy: optmized energy
        pos: optimized positions
    """
    pos = pos.flatten()
    fun, jac = LJ, LJ_force
    res = minimize(fun, pos, jac=jac, method='CG', tol=1e-3)
    pos = np.reshape(res.x, (int(len(pos)/3), 3))
    energy = res.fun
    return energy, pos

if __name__ == "__main__":
    N = 38
    L = 10
    for i in range(20):
        t0 = time()
        pos0 = L*np.random.random_sample((N*3,))
        eng, pos = single_optimize(pos0)
        print(i, eng, time()-t0)
