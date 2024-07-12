"""This is a script to:
1, generate random clusters.
2, perform optimization.
3, compare the efficiency of different algos (CG, BFGS).
"""

import logging
import warnings
from optparse import OptionParser
import numpy as np

from pyxtal import pyxtal
from pyxtal.database.collection import Collection
from pyxtal.optimize.myscipy_optimize import (
    _minimize_bfgs,
    _minimize_cg,
    _minimize_tpgd,
)
from pyxtal.potentials.LJ_cluster import LJ, LJ_force

warnings.filterwarnings("ignore")
logging.basicConfig(format="%(asctime)s :: %(message)s", filename="results.log", level=logging.INFO)


def single_optimize(pos, dim=3, kt=0.5, mu=0.1, beta=1.001, shift=False, method="mycg"):
    """Perform optimization for a given cluster.

    Args:
        pos: N*dim0 array representing the atomic positions
        dim: dimension of the hyper/normal space
        kt: perturbation factors
        mu: the weight for the punishing function
        beta: the factor for the punishing function
        shift: whether to shift the cluster to the center
        method: the optimization method

    output:
        energy: optmized energy
        pos: optimized positions
        niter: number of iterations
    """
    N_atom = len(pos)
    pos = pos.flatten()
    res = _minimize_tpgd(LJ, pos, args=(dim, mu, shift), jac=LJ_force, beta=beta, gtol=1e-4, maxiter=50)
    niter = res.nit
    pos = res.x
    if method == "mycg":
        res = _minimize_cg(LJ, pos, args=(dim, mu, shift), jac=LJ_force, beta=beta, gtol=1e-4)
    elif method == "mybfgs":
        res = _minimize_bfgs(LJ, pos, args=(dim, mu, shift), jac=LJ_force, beta=beta, gtol=1e-4)
    elif method == "mytpgd":
        res = _minimize_tpgd(LJ, pos, args=(dim, mu, shift), jac=LJ_force, beta=beta, gtol=1e-4)
    niter += res.nit
    energy = res.fun
    pos = np.reshape(res.x, (N_atom, dim))
    return energy, pos, niter


class LJ_prediction:
    """A class to perform global optimization on LJ clusters."""

    def __init__(self, numIons):
        """Initialize the class with the number of ions."""
        self.numIons = numIons
        ref = Collection("clusters")[str(numIons)]
        print(f"\nReference for LJ {numIons:3d} is {ref["energy"]:12.3f} eV, PG: {ref["pointgroup"]:4s}")

        self.reference = ref

    def generate_cluster(self, pgs: tuple[int, int], cluster_factor=0.6, seed=None):
        """Generate a random cluster with a given point group and number of ions."""
        rng = np.random.default_rng(seed)
        run = True
        while run:
            pg = rng.integers(pgs)
            cluster = pyxtal()
            cluster.from_random(0, pg, ["H"], [self.numIons], cluster_factor)
            if cluster.valid:
                run = False
        return cluster._get_coords_and_species(absolute=True)[0]

    def predict(self, maxN: int = 100, ncpu: int = 2, pgs: tuple[int, int] = (2, 33)):
        """Predict the energy of LJ clusters."""
        cycle = range(maxN)
        if ncpu > 1:
            from functools import partial
            from multiprocessing import Pool

            with Pool(ncpu) as p:
                func = partial(self.relaxation, pgs)
                res = p.map(func, cycle)
                p.close()
                p.join()
        else:
            res = [self.relaxation(pgs) for _ in cycle]

        return res

    def relaxation(self, pgs: tuple[int, int], **kwargs):
        """Perform relaxation for a given cluster."""
        pos = self.generate_cluster(pgs)
        res = []
        for method in ["mycg", "mybfgs", "mytpgd"]:
            pos1 = pos.copy()
            energy1, pos1, it1 = single_optimize(pos1, method=method)
            print(f"Optmization {method:10s}  3D: {energy1:10.4f}  {it1:6d} ")
            res.append([energy1, it1])
        res = np.array(res)
        return res.flatten()


if __name__ == "__main__":
    # -------------------------------- Options -------------------------
    parser = OptionParser()
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

    (options, args) = parser.parse_args()

    N = options.numIons
    maxN = options.max
    ncpu = options.proc

    lj_run = LJ_prediction(N)
    res = lj_run.predict(maxN=maxN, ncpu=ncpu, pgs=[1])
    res = np.array(res)
    eng_min = lj_run.reference["energy"]
    eng_cg3, iter_cg3 = res[:, 0], res[:, 1]
    eng_bf3, iter_bf3 = res[:, 2], res[:, 3]
    eng_gd3, iter_gd3 = res[:, 4], res[:, 5]
    iter_min, iter_max = np.min(res[:, [1, 3, 5]]), np.max(res[:, [1, 3, 5]])
    ground_cg3 = len(eng_cg3[eng_cg3 - 1e-2 < eng_min])
    ground_bf3 = len(eng_bf3[eng_bf3 - 1e-2 < eng_min])
    ground_gd3 = len(eng_gd3[eng_gd3 - 1e-2 < eng_min])

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    plt.switch_backend("agg")
    plt.style.use("bmh")

    e_max = eng_min * 0.87
    bins_eng = np.linspace(eng_min - 0.1, e_max, 100)
    bins_iter = np.linspace(iter_min, iter_max, 100)
    fig = plt.figure(figsize=(10, 8))

    gs = GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    label1 = "CG: " + str(ground_cg3) + "/" + str(len(eng_cg3))
    label2 = "BFGS: " + str(ground_bf3) + "/" + str(len(eng_bf3))
    label3 = "TPGD: " + str(ground_gd3) + "/" + str(len(eng_gd3))
    ax1.hist(eng_cg3, bins_eng, alpha=0.5, label=label1)
    ax1.hist(eng_bf3, bins_eng, alpha=0.5, label=label2)
    ax1.hist(eng_gd3, bins_eng, alpha=0.5, label=label3)
    ax1.set_xlim([eng_min - 0.1, e_max])
    ax1.set_xlabel("Energy (eV)")
    ax1.set_ylabel("Counts")
    ax1.legend()

    label1 = "CG: " + str(ground_cg3) + "/" + str(len(eng_cg3))
    label2 = "BFGS: " + str(ground_bf3) + "/" + str(len(eng_bf3))
    label3 = "TPGD: " + str(ground_gd3) + "/" + str(len(eng_gd3))
    ax2.hist(iter_cg3, bins_iter, alpha=0.5, label=label1)
    ax2.hist(iter_bf3, bins_iter, alpha=0.5, label=label2)
    ax2.hist(iter_gd3, bins_iter, alpha=0.5, label=label3)
    ax2.set_xlabel("# of iterations")
    ax2.set_ylabel("Counts")
    ax2.legend()

    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    x = np.linspace(eng_min, e_max, 10)

    ax3.plot(x, x, "k-", lw=1)
    ax3.scatter(eng_cg3, eng_bf3, c=iter_cg3 - iter_bf3, label="CG .v.s BFGS", s=5)
    ax3.set_xlabel("CG (eV)")
    ax3.set_ylabel("BFGS (eV)")
    ax3.set_xlim([eng_min - 0.1, e_max])
    ax3.set_ylim([eng_min - 0.1, e_max])
    ax3.legend()

    ax4.plot(x, x, "k-", lw=1)
    ax4.scatter(eng_cg3, eng_gd3, c=iter_cg3 - iter_gd3, label="CG .v.s TPGD", s=5)
    ax4.set_xlabel("CG (eV)")
    ax4.set_ylabel("BFGS (eV)")
    ax4.set_xlim([eng_min - 0.1, e_max])
    ax4.set_ylim([eng_min - 0.1, e_max])
    ax4.legend()
    plt.tight_layout()
    plt.savefig("LJ" + str(N) + "-" + str(maxN) + ".png")
