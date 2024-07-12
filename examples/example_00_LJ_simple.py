"""This script:
1. Generates random clusters
2. Performs optimization using the Lennard-Jones potential
"""

from time import time

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import pdist


def lennard_jones_potential(pos, dim=3):
    """Calculate the total Lennard-Jones potential energy.

    Args:
        pos (np.ndarray): 1D array with N*dim numbers representing atomic positions
        dim (int): Dimension of the space (default: 3)

    Returns:
        float: Total energy with punishing function
    """
    N_atom = len(pos) // dim
    pos = pos.reshape(N_atom, dim)

    distance = pdist(pos)
    r6 = distance**6
    r12 = r6**2
    return np.sum(4 * (1 / r12 - 1 / r6))


def lennard_jones_force(pos, dim=3):
    """Calculate the Lennard-Jones forces.

    Args:
        pos (np.ndarray): 1D array with N*dim numbers representing atomic positions
        dim (int): Dimension of the space (default: 3)

    Returns:
        np.ndarray: 1D force array
    """
    N_atom = len(pos) // dim
    pos = pos.reshape(N_atom, dim)
    force = np.zeros_like(pos)

    for i, pos0 in enumerate(pos):
        pos1 = np.delete(pos, i, axis=0)
        r = pos1 - pos0
        distance = np.linalg.norm(r, axis=1)
        r2 = distance**2
        r6 = r2**3
        r12 = r6**2
        force[i] = np.sum((48 / r12 - 24 / r6)[:, np.newaxis] * r / r2[:, np.newaxis], axis=0)

    return force.flatten()


def optimize_cluster(pos):
    """Perform optimization for a given cluster.

    Args:
        pos (np.ndarray): N*dim array representing atomic positions

    Returns:
        tuple: (optimized energy, optimized positions)
    """
    pos = pos.flatten()
    res = minimize(lennard_jones_potential, pos, jac=lennard_jones_force, method="CG", tol=1e-3)
    optimized_pos = res.x.reshape(-1, 3)
    optimized_energy = res.fun
    return optimized_energy, optimized_pos


if __name__ == "__main__":
    N_ATOMS = 38
    BOX_SIZE = 10
    N_ITERATIONS = 20

    rng = np.random.default_rng()

    for i in range(N_ITERATIONS):
        t0 = time()
        initial_pos = BOX_SIZE * rng.random((N_ATOMS * 3,))
        energy, optimized_pos = optimize_cluster(initial_pos)
        print(f"Iteration {i}: Energy = {energy:.6f}, Time = {time() - t0:.6f} s")
