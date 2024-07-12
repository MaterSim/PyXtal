import numpy as np
from scipy.spatial.distance import cdist, pdist

"""
LJ energy and force functions
"""


def LJ(pos, dim, mu=0.1, shift=False):
    """
    Calculate the total energy
    Args:
    pos: 1D array with N*dim numbers representing the atomic positions
    dim: dimension of the hyper/normal space
    output
    E: the total energy with punishing function
    """
    N_atom = int(len(pos) / dim)
    pos = np.reshape(pos, (N_atom, dim))

    distance = pdist(pos)
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    Eng = np.sum(4 * (1 / r12 - 1 / r6))

    if dim > 3:
        norm = 0
        for i in range(3, dim):
            diff = pos[:, i] - np.mean(pos[:, i]) if shift else pos[:, i]
            norm += np.sum(np.power(diff, 2))
        Eng += 0.5 * mu * norm
    return Eng


def LJ_force(pos, dim, mu=0.1, shift=False):
    N_atom = int(len(pos) / dim)
    pos = np.reshape(pos, [N_atom, dim])
    force = np.zeros([N_atom, dim])
    for i, pos0 in enumerate(pos):
        pos1 = pos.copy()
        pos1 = np.delete(pos1, i, 0)
        distance = cdist([pos0], pos1)
        r = pos1 - pos0
        r2 = np.power(distance, 2)
        r6 = np.power(r2, 3)
        r12 = np.power(r6, 2)
        force[i] = np.dot((48 / r12 - 24 / r6) / r2, r)
        # force from the punish function mu*sum([x-mean(x)]^2)
        if dim > 3:
            for j in range(3, dim):
                if shift:
                    force[i, j] += mu * (pos[i, j] - np.mean(pos[:, j]))
                else:
                    force[i, j] += mu * pos[i, j]
    return force.flatten()
