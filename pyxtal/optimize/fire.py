import numpy as np
from time import time
from scipy.spatial.distance import pdist, cdist
from copy import deepcopy
import random

"""
Scripts to compute LJ energy and force
"""


def LJ(pos):
    """
    Calculate the total energy
    """
    distance = pdist(pos)
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    Eng = np.sum(4 * (1 / r12 - 1 / r6))
    return Eng


def LJ_force(pos):
    N_atom = len(pos)
    force = np.zeros([N_atom, 3])
    for i, pos0 in enumerate(pos):
        pos1 = deepcopy(pos)
        pos1 = np.delete(pos1, i, 0)
        distance = cdist([pos0], pos1)
        r = pos1 - pos0
        r2 = np.power(distance, 2)
        r6 = np.power(r2, 3)
        r12 = np.power(r6, 2)
        force[i] = np.dot((48 / r12 - 24 / r6) / r2, r)
    return force


class FIRE:
    def __init__(
        self,
        xstruct,
        E,
        F,
        dt=0.1,
        maxmove=0.2,
        dtmax=1.0,
        Nmin=10,
        finc=1.1,
        fdec=0.5,
        astart=0.1,
        fa=0.99,
        a=0.1,
    ):
        """Parameters:
        An optimization engine based on the fire algorithm
        1, impose the symmetry operation
        2, enable the box condition
        """

        self.dt = dt
        self.maxmove = maxmove
        self.dtmax = dtmax
        self.Nmin = Nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart
        self.fa = fa
        self.a = a
        self.nsteps = 0
        self.Nsteps = 0
        self.get_energy = E
        self.get_force = F
        self.trajectory = {
            "position": [],
            "energy": [],
            "force": [],
            "v": [],
            "time": [],
        }
        self.time0 = time()
        self.xstruct = xstruct
        self.pos = xstruct.cart_coords
        self.get_energy = E
        self.get_force = F
        self.initialize()

    def initialize(self):
        """
        initilaize the positions, energy, force and veolcity for the 1st step
        """
        self.energy = self.get_energy(self.pos)
        self.force = self.get_force(self.pos)
        self.v = np.zeros((len(self.pos), 3))
        self.trajectory["position"].append(self.pos)
        self.trajectory["energy"].append(self.energy)
        self.trajectory["force"].append(self.force)
        self.trajectory["v"].append(self.v)
        self.trajectory["time"].append(time() - self.time0)

    def update(self):
        self.energy = self.get_energy(self.pos)
        self.force = self.get_force(self.pos)
        self.trajectory["position"].append(self.pos)
        self.trajectory["v"].append(self.v)
        self.trajectory["energy"].append(self.energy)
        self.trajectory["force"].append(self.force)
        self.trajectory["time"].append(time() - self.time0)
        fmax = np.max(np.abs(self.force.flatten()))
        print(
            "Step {0:4d} Eng: {1:12.4f} Fmax: {2:12.4f} Time: {3:6.2f} seconds".format(
                self.nsteps, self.energy, fmax, time() - self.time0
            )
        )

    def step(self):
        f = self.force
        vf = np.vdot(f, self.v)
        if vf > 0.0:
            self.v = (1.0 - self.a) * self.v + self.a * f / np.sqrt(
                np.vdot(f, f)
            ) * np.sqrt(np.vdot(self.v, self.v))
            if self.Nsteps > self.Nmin:
                self.dt = min(self.dt * self.finc, self.dtmax)
                self.a *= self.fa
            self.Nsteps += 1
        else:
            self.v[:] *= 0.0
            self.a = self.astart
            self.dt *= self.fdec
            self.Nsteps = 0
        self.v += self.dt * f
        dr = self.dt * self.v  # needs to impose constraints
        normdr = np.sqrt(np.vdot(dr, dr))
        if normdr > self.maxmove:
            dr = self.maxmove * dr / normdr

        # print(self.force)
        # print(self.v)
        # print(dr)
        self.pos = self.pos - dr

        # Apply symmetrization after dr is subtracted
        # self.pos = self.symmetrized_coords(self.pos)

        self.update()

    def run(self, max_steps=1000):
        while self.nsteps < max_steps:
            if self.check_convergence():
                break
            self.step()
            self.nsteps += 1

    def check_convergence(self):
        converged = False
        if np.max(self.force.flatten()) < 1e-3:
            if self.nsteps > 0:
                if self.trajectory["energy"][-1] - self.trajectory["energy"][-2] < 1e-3:
                    converged = True
        return converged

    def symmetrized_coords(self, coords):
        gen_coords = []
        gen_ops = []
        start_index = 0
        new_coords = []
        for ws in self.xstruct.wyckoff_sites:
            gen_coord = coords[start_index]
            gen_coord = ws.wp[0].operate(gen_coord)
            wp_coords = apply_ops(gen_coord, ws.wp.generators_m)
            if new_coords == []:
                new_coords = wp_coords
            else:
                new_coords = np.vstack([new_coords, wp_coords])
            start_index += ws.multiplicity
        return new_coords


# from pyxtal.interface.LJ import LJ, LJ_force
from pyxtal.database.collection import Collection

np.random.seed(10)

# pos = np.array(Collection('clusters')['20']['position']['data'])
from pyxtal.crystal import *

pg = random.choice(range(1, 57))
c = random_cluster(pg, ["C"], [20], 1.0)
pos = c.cart_coords

print("Point group " + c.group.symbol)
print("Initial energy: " + str(LJ(pos)))
print("Initial coords:")
print(pos)

pos += 0.5 * np.random.uniform(-1, 1, (len(pos), 3))
np.savetxt("1.txt", pos)
dyn = FIRE(c, LJ, LJ_force)
dyn.run(1000)
