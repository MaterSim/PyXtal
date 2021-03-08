from pyxtal.crystal import random_cluster
from copy import deepcopy
from optparse import OptionParser
from random import randint, choice
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, cdist
from pyxtal.molecule import PointGroupAnalyzer
from pymatgen.core import Molecule
from pyxtal.database.collection import Collection
from time import time
import numpy as np
import matplotlib.pyplot as plt
import warnings
import sys

plt.style.use("bmh")
warnings.filterwarnings("ignore")

"""
This is a script to 
1, generate random clusters
2, perform optimization
"""


def LJ(pos, dim, method=1, mu=0.1):
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
            if method == 1:
                diff = pos[:, i]
            else:
                diff = pos[:, i] - np.mean(pos[:, i])
            norm += np.sum(np.multiply(diff, diff))
        Eng += 0.5 * mu * norm
    return Eng


def LJ_force(pos, dim, method=1, mu=0.1):
    N_atom = int(len(pos) / dim)
    pos = np.reshape(pos, [N_atom, dim])
    force = np.zeros([N_atom, dim])
    for i, pos0 in enumerate(pos):
        pos1 = deepcopy(pos)
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
                if method == 1:
                    force[i, j] += mu * pos[i, j]  # - np.mean(pos[:, j]))
                else:
                    force[i, j] += mu * (pos[i, j] - np.mean(pos[:, j]))
    return force.flatten()


def single_optimize(pos, dim=3, method=1, mu=0.1):
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
        # pos = np.hstack((pos, 0.5*(np.random.random([N_atom, diff])-0.5) ))
        pos = np.hstack((pos, np.random.uniform(-1, 1, (N_atom, diff))))
    elif diff < 0:
        pos = pos[:, :dim]

    pos = pos.flatten()
    res = minimize(LJ, pos, args=(dim, method, mu), jac=LJ_force, method="CG", tol=1e-3)
    pos = np.reshape(res.x, (N_atom, dim))
    energy = res.fun
    return energy, pos


def hyper_optimize(pos, dim, method=1, mu=0.1):
    """
    hyperspatial optimization 
    """
    [energy1, pos1] = single_optimize(pos, dim)

    while True:
        if dim == len(pos1[0]):
            disp = np.random.uniform(-1.0, 1, (len(pos), dim))
            pos2 = pos1 + disp
            [energy3, pos3] = single_optimize(pos2, dim)
        else:
            [energy2, pos2] = single_optimize(pos1, dim, method, mu)
            [energy3, pos3] = single_optimize(pos2, 3)

        if energy3 - energy1 > 1e-4:
            pos = pos1
            energy = energy1
            break
        else:
            pos1 = pos3
            energy1 = energy3
    return energy1, pos1


def parse_symmetry(pos):
    mol = Molecule(["C"] * len(pos), pos)
    try:
        symbol = PointGroupAnalyzer(mol, tolerance=0.1).sch_symbol
    except:
        symbol = "N/A"
    return symbol


class LJ_prediction:
    """
    A class to perform global optimization on LJ clusters
    Args:

    Attributes:

    """

    def __init__(self, numIons):
        self.numIons = numIons
        ref = Collection("clusters")[str(numIons)]
        print(
            "\nReference for LJ {0:3d} is {1:12.3f} eV, PG: {2:4s}".format(
                numIons, ref["energy"], ref["pointgroup"]
            )
        )
        self.reference = ref
        self.time0 = time()

    def generate_cluster(self, pgs=range(2, 33)):
        run = True
        while run:
            pg = choice(pgs)
            cluster = random_cluster(pg, ["Mo"], [self.numIons], 1.0)
            if cluster.valid:
                run = False
        return cluster.cart_coords

    def predict(self, dim=3, maxN=100, ncpu=2, pgs=range(2, 33), method=1):

        print("\nPerforming random search at {0:d}D space\n".format(dim))
        cycle = range(maxN)
        if ncpu > 1:
            from multiprocessing import Pool
            from functools import partial

            with Pool(ncpu) as p:
                func = partial(self.relaxation, dim, pgs, method)
                res = p.map(func, cycle)
                p.close()
                p.join()
        else:
            res = []
            for i in cycle:
                res.append(self.relaxation(dim, pgs, method, i))

        N_success = 0
        for dct in res:
            if dct["ground"]:
                N_success += 1
        print(
            "\nHit the ground state {0:4d} times out of {1:4d} attempts\n".format(
                N_success, maxN
            )
        )
        return res

    def relaxation(self, dim, pgs, method, ind):
        pos = self.generate_cluster(pgs)
        pg1 = parse_symmetry(pos)
        print("initial geometry:")
        print(pos)
        print("initial energy: ", LJ(pos.flatten(), 3))
        if dim == 3:
            [energy, pos] = single_optimize(pos, 3)
            energy = [energy]
        else:
            [energy1, pos1] = single_optimize(pos, 3)
            [energy2, pos2] = hyper_optimize(pos1, 3)
            # [energy2, pos2] = hyper_optimize(pos1, dim, method=1, mu=0.1)
            # [energy3, pos3] = hyper_optimize(pos1, dim, method=1, mu=0.2)
            # [energy4, pos4] = hyper_optimize(pos1, dim, method=2, mu=0.1)
            # [energy5, pos5] = hyper_optimize(pos1, dim, method=2, mu=0.2)
            # [energy6, pos6] = hyper_optimize(pos1, 3)
            energy = [energy1, energy2]  # , energy3, energy4, energy5, energy6]
            pos = [pos1, pos2]  # , pos3, pos4, pos5, pos6]
        if min(energy) - self.reference["energy"] < 1e-3:
            ground = True
        elif min(energy) > -10:
            print("high energy" + str(energy) + " after relaxation")
            print(pos1)
            sys.exit()
        else:
            ground = False

        pg2 = parse_symmetry(pos)
        res = {
            "pos_init": pos,
            "energy": energy,
            "pg_init": pg1,
            "pg_finial": pg2,
            "ground": ground,
            "id": ind,
        }

        if ground:
            print(
                "ID: {0:4d} PG initial: {1:4s} relaxed: {2:4s} Energy: {3:12.3f} Time: {4:6.1f} ++++++".format(
                    ind, pg1, pg2, energy[-1], (time() - self.time0) / 60
                )
            )
        elif ind % 2 == 0:
            print(
                "ID: {0:4d} PG initial: {1:4s} relaxed: {2:4s} Energy: {3:12.3f} Time: {4:6.1f} ".format(
                    ind, pg1, pg2, energy[-1], (time() - self.time0) / 60
                )
            )
        return res


if __name__ == "__main__":
    # -------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option(
        "-d",
        "--dimension",
        dest="dim",
        metavar="dim",
        default=3,
        type=int,
        help="dimension, 3 or higher",
    )
    parser.add_option(
        "-n",
        "--numIons",
        dest="numIons",
        default=16,
        type=int,
        help="desired numbers of atoms: 16",
    )
    parser.add_option(
        "-m",
        "--max",
        dest="max",
        default=100,
        type=int,
        help="maximum number of attempts",
    )
    parser.add_option(
        "-p",
        "--proc",
        dest="proc",
        default=1,
        type=int,
        help="number of processors, default 1",
    )
    parser.add_option(
        "-f",
        "--func",
        dest="func",
        default=1,
        type=int,
        help="penalty function, default 1: mu*r^2",
    )

    (options, args) = parser.parse_args()

    N = options.numIons  # 38
    maxN = options.max  # 1000
    dim = options.dim  # 4
    ncpu = options.proc
    method = options.func

    lj_run = LJ_prediction(N)
    eng_min = lj_run.reference["energy"]
    t0 = time()
    results1 = lj_run.predict(dim=4, maxN=maxN, ncpu=ncpu, pgs=[1])
    print("time: {0:6.2f} seconds".format(time() - t0))
    engs = []
    for dct in results1:
        engs.append(dct["energy"])
    engs = np.array(engs)
    for i in range(1, len(engs[0])):
        print("method " + str(i), np.mean(engs[:, i] - engs[:, 0]))
    grounds = []
    for i in range(len(engs[0])):
        eng_tmp = engs[:, i]
        grounds.append(len(eng_tmp[eng_tmp < (eng_min + 1e-3)]))

