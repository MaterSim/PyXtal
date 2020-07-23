import numpy as np
from time import time
from scipy.spatial.distance import pdist, cdist
from copy import deepcopy
from pyxtal.molecule import PointGroupAnalyzer
import random

"""
Scripts to compute LJ energy and force
"""


def check_struct_group(crystal, group, dim=3, tol=1e-2):
    # Supress pymatgen/numpy complex casting warnings
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        """Given a pymatgen structure, group number, and dimension, return
        whether or not the structure matches the group number."""
        if type(crystal) == random_crystal:
            lattice = struct.lattice.matrix
            if dim != 0:
                old_coords = deepcopy(crystal.struct.frac_coords)
                old_species = deepcopy(crystal.struct.atomic_numbers)
            elif dim == 0:
                old_coords = deepcopy(crystal.cart_coords)
                old_species = deepcopy(crystal.species)
        else:
            lattice = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            old_coords = np.array(crystal)
            old_species = ["C"] * len(old_coords)

        from pyxtal.symmetry import distance
        from pyxtal.symmetry import filtered_coords
        from copy import deepcopy

        PBC = [1, 1, 1]

        # Obtain the generators for the group
        if dim == 3:
            from pyxtal.symmetry import get_wyckoffs

            generators = get_wyckoffs(group)[0]

        elif dim == 2:
            from pyxtal.symmetry import get_layer

            generators = get_layer(group)[0]
            PBC = [1, 1, 0]
        elif dim == 1:
            from pyxtal.symmetry import get_rod

            generators = get_rod(group)[0]
            PBC = [0, 0, 1]
        elif dim == 0:
            from pyxtal.symmetry import Group

            generators = Group(group, dim=0)[0]
            PBC = [0, 0, 0]

        # TODO: Add check for lattice symmetry

        # Apply SymmOps to generate new points
        # old_coords = filtered_coords(struct.frac_coords,PBC=PBC)

        new_coords = []
        new_species = []
        for i, point in enumerate(old_coords):
            for j, op in enumerate(generators):
                if j != 0:
                    new_coords.append(op.operate(point))
                    new_species.append(old_species[i])
        # new_coords = filtered_coords(new_coords,PBC=PBC)

        # Check that all points in new list are still in old
        failed = False
        i_list = list(range(len(new_coords)))
        for i, point1 in enumerate(new_coords):
            found = False
            for j, point2 in enumerate(old_coords):
                if new_species[i] == old_species[j]:
                    difference = filtered_coords(point2 - point1, PBC=PBC)
                    if distance(difference, lattice, PBC=PBC) <= tol:
                        found = True
                        break
            if found is False:
                failed = True
                break

        if failed is False:
            return True
        else:
            return False


def parse_symmetry(pos):
    mol = Molecule(["C"] * len(pos), pos)
    try:
        symbol = PointGroupAnalyzer(mol, tolerance=0.1).sch_symbol
    except:
        symbol = "N/A"
    return symbol


def LJ_1d(pos, dim=3):
    N_atom = int(len(pos) / dim)
    pos = np.reshape(pos, (N_atom, dim))
    distance = pdist(pos)
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    Eng = np.sum(4 * (1 / r12 - 1 / r6))
    return Eng


def LJ_force_1d(pos, dim=3):
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
    return force.flatten()


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
        crystal,
        E,
        F,
        symmetrize=False,
        e_tol=1e-5,
        f_tol=1e-2,
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
        self.symmetrize = symmetrize
        self.f_tol = f_tol
        self.e_tol = e_tol
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
        try:
            self.xstruct = crystal
            self.pos = crystal.cart_coords
        except:
            self.pos = crystal
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
        print("Initial Energy: {:12.4f}".format(self.energy))

    def update(self, freq=500):
        self.energy = self.get_energy(self.pos)
        self.force = self.get_force(self.pos)
        self.trajectory["position"].append(self.pos)
        self.trajectory["v"].append(self.v)
        self.trajectory["energy"].append(self.energy)
        self.trajectory["force"].append(self.force)
        self.trajectory["time"].append(time() - self.time0)
        fmax = np.max(np.abs(self.force.flatten()))

        if self.nsteps % freq == 0:
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

        # Symmetrize the force
        if self.symmetrize:
            f = self.symmetrized_coords(f)

        self.v += self.dt * f
        dr = self.dt * self.v  # needs to impose constraints
        normdr = np.sqrt(np.vdot(dr, dr))
        if normdr > self.maxmove:
            dr = self.maxmove * dr / normdr

        # print(self.force)
        # print(self.v)
        # print(dr)
        self.pos = self.pos - dr

        # Symmetrize positions
        # if self.symmetrize:
        #    self.pos = self.symmetrized_coords(self.pos)

        self.update()

    def run(self, max_steps=1000):
        self.max_steps = max_steps
        while self.nsteps < max_steps:
            if self.check_convergence():
                break
            self.step()
            self.nsteps += 1

        fmax = np.max(np.abs(self.force.flatten()))
        print(
            "Finish at Step {0:4d} Eng: {1:12.4f} Fmax: {2:12.4f} Time: {3:6.2f} seconds".format(
                self.nsteps, self.energy, fmax, time() - self.time0
            )
        )

    def check_convergence(self):
        """
        There should be two options to terminate the optimization before it reaches the maximum number of cycles
        1, by energy difference compared to the previous step: self.e_tol
        2, by the residual force: self.f_tol
        Will implement both options later.
        """
        converged = False
        if np.max(self.force.flatten()) < self.f_tol:
            if self.nsteps > 0:
                if self.trajectory["energy"][-1] - self.trajectory["energy"][-2] < 1e-3:
                    converged = True
        return converged

    def symmetrized_coords(self, coords):
        start_index = 0
        new_coords = []
        for ws in self.xstruct.wyckoff_sites:
            gen_coord = coords[start_index]
            wp_coords = apply_ops(gen_coord, ws.wp)
            if len(new_coords) == 0:
                new_coords = wp_coords
            else:
                new_coords = np.vstack([new_coords, wp_coords])
            start_index += ws.multiplicity
        return new_coords


# Set global variables
dim = 0

# from pyxtal.interface.LJ import LJ, LJ_force
from pyxtal.database.collection import Collection

np.random.seed(10)
pos = np.array(Collection("clusters")["20"]["position"]["data"])

# Generate a random cluster with random point group
from pyxtal.crystal import *
from copy import deepcopy

pg = random.choice(range(1, 57))
c = random_cluster(pg, ["C"], [20], 0.7)
pos = c.cart_coords

# maxr = 0
# for p in pos:
#    r = np.linalg.norm(p)
#    if r > maxr:
#        maxr = r
# print("maxr: "+str(maxr))

print("Point group: ", c.group.symbol, parse_symmetry(pos))
print("Initial energy: " + str(LJ(pos)))
# print(LJ(pos))
# pos += 0.5*np.random.uniform(-1, 1, (len(pos), 3))
# np.savetxt('1.txt', pos)

print("\n optmization without symmetry constraints")
c0 = deepcopy(c)
dyn1 = FIRE(c0, LJ, LJ_force, f_tol=1e-2, dt=0.2, maxmove=2.0)
dyn1.run(2000)
print("symmetry after optimization: ", parse_symmetry(dyn1.pos))
print("Initial point group preserved: ", check_struct_group(dyn1.pos, pg, dim=dim))

print("\n optmization under symmetry constraints")
c1 = deepcopy(c)
dyn2 = FIRE(c1, LJ, LJ_force, symmetrize=True, f_tol=1e-2, dt=0.2, maxmove=2.0)
dyn2.run(2000)
print("symmetry after optimization: ", parse_symmetry(dyn2.pos))
print("Initial point group preserved: ", check_struct_group(dyn2.pos, pg, dim=dim))

from scipy.optimize import minimize

print("\n optmization with scipy CG")

c1 = deepcopy(c)
pos = c1.cart_coords.flatten()
res = minimize(LJ_1d, pos, jac=LJ_force_1d, method="CG", tol=1e-3)
print(res.fun)
pos = res.x
print(
    "symmetry after optimization: ",
    parse_symmetry(np.reshape(pos, (int(len(pos) / 3), 3))),
)
print(
    "Initial point group preserved: ",
    check_struct_group(np.reshape(pos, (int(len(pos) / 3), 3)), pg, dim=dim),
)

c1 = deepcopy(c)
pos = c1.cart_coords.flatten()
print("\n optmization with scipy BFGS")
pos = pos.flatten()

res = minimize(LJ_1d, pos, jac=LJ_force_1d, method="BFGS", tol=1e-3)
print(res.fun)
pos = res.x
print(
    "symmetry after optimization: ",
    parse_symmetry(np.reshape(pos, (int(len(pos) / 3), 3))),
)
print(
    "Initial point group preserved: ",
    check_struct_group(np.reshape(pos, (int(len(pos) / 3), 3)), pg, dim=dim),
)

