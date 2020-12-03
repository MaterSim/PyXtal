from pyxtal import pyxtal
from copy import deepcopy
from optparse import OptionParser
from random import randint, choice
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, cdist
from pyxtal.molecule import PointGroupAnalyzer
from pymatgen import Molecule
from pyxtal.database.collection import Collection
from time import time
import numpy as np
import matplotlib.pyplot as plt

plt.switch_backend("agg")
import logging
import warnings

warnings.filterwarnings("ignore")
logging.basicConfig(
    format="%(asctime)s :: %(message)s", filename="results.log", level=logging.INFO
)
plt.style.use("bmh")

cluster_factor = 0.7

"""
This is a script to 
1, generate random clusters
2, perform optimization 
3, compare the efficiency of different algos (CG, BFGS)
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
    N_atom = int(len(pos) / dim)
    pos = np.reshape(pos, (N_atom, dim))

    distance = pdist(pos)
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    Eng = np.sum(4 * (1 / r12 - 1 / r6))

    if dim > 3:
        norm = 0
        for i in range(3, dim):
            # diff = pos[:, i] - np.mean(pos[:, i])
            diff = pos[:, i]
            norm += np.sum(np.multiply(diff, diff))
        Eng += 0.5 * mu * norm
    return Eng


def LJ_force(pos, dim, mu=0.1):
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
                # force[i, j] += mu*(pos[i, j] - np.mean(pos[:, j]))
                force[i, j] += mu * pos[i, j]  # - np.mean(pos[:, j]))
    return force.flatten()


def single_optimize(pos, dim=3, kt=0.5, mu=0.1, method="CG"):
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
        pos = np.hstack((pos, 0.5 * (np.random.random([N_atom, diff]) - 0.5)))
    elif diff < 0:
        pos = pos[:, :dim]

    pos = pos.flatten()
    res = minimize(LJ, pos, args=(dim, mu), jac=LJ_force, method=method, tol=1e-3)
    pos = np.reshape(res.x, (N_atom, dim))
    energy = res.fun
    return energy, pos


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
            cluster = pyxtal()
            cluster.from_random(0, pg, ["Mo"], [self.numIons], factor=cluster_factor)
            if cluster.valid:
                run = False
        return cluster._get_coords_and_species(absolute=True)[0]

    def predict(self, dim=3, maxN=100, ncpu=2, pgs=range(2, 33), method="CG"):

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
        pos0 = self.generate_cluster(pgs)
        pg0 = parse_symmetry(pos0)
        pos1 = pos0.copy()
        pos2 = pos0.copy()
        [energy1, pos1] = single_optimize(pos1, 3, method="CG")
        [energy2, pos2] = single_optimize(pos2, 3, method="BFGS")
        ground1, ground2 = False, False
        if abs(energy1 - self.reference["energy"]) < 1e-3:
            ground1 = True
        if abs(energy2 - self.reference["energy"]) < 1e-3:
            ground2 = True

        pg1 = parse_symmetry(pos1)
        pg2 = parse_symmetry(pos2)
        res = {
            "pos": [pos0, pos1, pos2],
            "energy": [energy1, energy2],
            "pg_init": pg0,
            "pg_finial": [pg1, pg2],
            "ground": [ground1, ground2],
            "id": ind,
        }
        if ground1 or ground2:
            logging.info(
                "ID: {:4d} PG initial: {:4s} relaxed: {:4s} {:4s} Energy: {:9.3f} {:9.3f}++++++++".format(
                    ind, pg0, pg1, pg2, energy1, energy2
                )
            )
        # elif ind%100 == 0:
        print(
            "ID: {:4d} PG initial: {:4s} relaxed: {:4s} {:4s} Energy: {:9.3f} {:9.3f}".format(
                ind, pg0, pg1, pg2, energy1, energy2
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
        "--factor",
        dest="cluster_factor",
        default=0.7,
        type=float,
        help="distance checking factor, default 0.7",
    )
    parser.add_option(
        "-v",
        "--volume",
        dest="volume_factor",
        default=1.0,
        type=float,
        help="volume factor, default 1.0",
    )

    (options, args) = parser.parse_args()

    N = options.numIons  # 38
    maxN = options.max  # 1000
    dim = options.dim  # 4
    ncpu = options.proc
    cluster_factor = options.cluster_factor
    volume_factor = options.volume_factor

    lj_run = LJ_prediction(N)
    eng_min = lj_run.reference["energy"]
    t0 = time()
    results1 = lj_run.predict(dim=dim, maxN=maxN, ncpu=ncpu, pgs=[1], method="CG")
    logging.info("-------------calculation is complete------------")
    ground1, ground2 = 0, 0
    engs = []
    for dct in results1:
        if dct["ground"][0]:
            ground1 += 1
        if dct["ground"][1]:
            ground2 += 1
        engs.append(dct["energy"])
    engs = np.array(engs)
    eng1 = engs[:, 0]
    eng2 = engs[:, 1]
    bins = np.linspace(eng_min - 0.1, eng_min + 20, 100)
    label1 = "3d CG: " + str(ground1) + "/" + str(len(eng1))
    label2 = "3d BFGS: " + str(ground2) + "/" + str(len(eng2))
    print(label1)
    print(label2)
    plt.hist(eng1, bins, alpha=0.5, label=label1)
    plt.hist(eng2, bins, alpha=0.5, label=label2)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Counts")
    plt.legend()
    eng_min_str = "{0:.2f}".format(eng_min)
    plt.title("LJ" + str(N) + " Ground state: " + eng_min_str)
    plt.savefig("LJ" + str(N) + "-" + str(maxN) + "samples-hist.png")
    plt.close()
    plt.scatter(eng1, eng2)
    plt.xlabel("Energy from CG optimization")
    plt.ylabel("Energy from BFGS optimization")
    plt.xlim([eng_min - 0.1, -20])
    plt.ylim([eng_min - 0.1, -20])
    plt.savefig("LJ" + str(N) + "-" + str(maxN) + "samples-scatter.png")
    plt.close()

