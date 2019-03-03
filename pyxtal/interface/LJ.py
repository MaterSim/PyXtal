import numpy as np
from time import time
from pymatgen.core.lattice import Lattice

eV2GPa = 160.217
GPa2eV = 1.0/eV2GPa

"""
Scripts to compute LJ energy and force
"""

def get_neighbors(struc, i, rcut):
    """
    script to get the neighbors info for an atom in the struct
    """
    lat = struc.lattice_matrix
    center = struc.cart_coords[i]
    c1 = struc.frac_coords[i]
    fcoords = struc.frac_coords
    r, dists, inds, images = Lattice(lat).get_points_in_sphere(
            fcoords, center=center, r=rcut, zip_results=False)
    ids = np.where(dists>0.1)
    r = r[ids] 
    r = np.dot(r, lat) - center
    return dists[ids]**2, r

class LJ():

    def __init__(self, struc, press=1e-4, epsilon=0.01, sigma=3.4, rcut=8.0):
        self.struc = struc
        self.press = press
        self.epsilon = epsilon
        self.sigma = sigma
        self.rcut = rcut
        pos = struc.cart_coords
        lat = struc.lattice_matrix
        volume = np.linalg.det(lat)
        #print('PV term', press*volume)


        # initiate the energy, force, stress
        energy = 0
        force = np.zeros([len(pos), 3])
        stress = np.zeros([3, 3])
        sigma6 = sigma**6
        sigma12 = sigma6*sigma6

        # calculating the energy/force, needs to update sigma/epsilons
        for i, pos0 in enumerate(pos):
            [r2, r] = get_neighbors(struc, i, rcut)
            r6 = np.power(r2, 3)     
            r12 = np.power(r6, 2)
            energy += np.sum(4.0*epsilon*(sigma12/r12 - sigma6/r6))
            f = (24*epsilon*(2.0*sigma12/r12-sigma6/r6)/r2)[:, np.newaxis] * r
            force[i] = f.sum(axis=0)
            stress += np.dot(f.T, r)

        self.energy = 0.5*energy    
        self.enthalpy = self.energy + press*volume*GPa2eV
        self.force = force
        self.stress = -0.5*stress/volume*eV2GPa

class FIRE():
    def __init__(self, crystal, symmetrize=False, e_tol=1e-5, f_tol=1e-2, dt=0.1, maxmove=0.1, dtmax=1.0, Nmin=10, finc=1.1, fdec=0.5,
                 astart=0.1, fa=0.99, a=0.1):
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
        self.trajectory = {"position": [],
                           "energy": [],
                           "force": [],
                           "v": [],
                           "time": [],
                          }
        self.time0 = time()
        self.initialize()

    def initialize(self):
        """
        initilaize the positions, energy, force and veolcity for the 1st step
        """
        self.crystal = crystal
        self.energy = LJ(self.crystal).energy
        self.force = LJ(self.crystal).force
        self.stress = LJ(self.crystal).stress
        self.fmax = np.max(np.abs(np.vstack((self.stress, self.force)).flatten()))
        self.v = np.zeros((len(self.force)+3, 3))
        print('Initial Energy: {:12.4f}'.format(self.energy))
        print('Initial stress: {:12.4f}'.format(self.stress[0,0]))

    def update(self, freq=10):
        if (self.nsteps % freq == 0):
            print('Step {0:4d} Eng: {1:12.4f} Fmax: {2:12.4f} Vol: {3:6.2f}'.
                format(self.nsteps, self.energy, self.fmax, self.volume))

    def step(self):
        f = np.vstack((self.stress, self.force))
        pos = np.vstack((self.crystal.lattice_matrix, self.crystal.cart_coords))

        vf = np.vdot(f, self.v)
        if vf > 0.0:
            self.v = (1.0 - self.a) * self.v + self.a * f / np.sqrt(
                np.vdot(f, f)) * np.sqrt(np.vdot(self.v, self.v))
            if self.Nsteps > self.Nmin:
                self.dt = min(self.dt * self.finc, self.dtmax)
                self.a *= self.fa
            self.Nsteps += 1
        else:
            self.v[:] *= 0.0
            self.a = self.astart
            self.dt *= self.fdec
            self.Nsteps = 0

        #Symmetrize the force
        if self.symmetrize:
            f = self.symmetrized_coords(f)

        self.v += self.dt * f
        dr = self.dt * self.v  #needs to impose constraints
        normdr = np.sqrt(np.vdot(dr, dr))
        if normdr > self.maxmove:
            dr = self.maxmove * dr / normdr

        #print(self.force)
        #print(self.v)
        pos = pos - dr
        self.crystal.lattice_matrix = pos[:3, :]
        self.crystal.cart_coords = pos[3:, :]
        self.crystal.frac_coords = np.dot(pos[3:, :], np.linalg.inv(pos[:3, :]))
        self.energy = LJ(self.crystal).energy
        self.force = LJ(self.crystal).force
        self.stress = LJ(self.crystal).stress
        self.fmax = np.max(np.abs(np.vstack((self.stress, self.force)).flatten()))
        self.volume = np.linalg.det(pos[:3, :])
        #Symmetrize positions
        #if self.symmetrize:
        #    self.pos = self.symmetrized_coords(self.pos)
    
        self.update()

    def run(self, max_steps=1000):
        self.max_steps = max_steps
        while self.nsteps < max_steps:
            if self.check_convergence():
                break
            self.step()
            self.nsteps += 1

        print('Finish at Step {0:4d} Eng: {1:12.4f} Fmax: {2:12.4f} Time: {3:6.2f} seconds'.
                format(self.nsteps, self.energy, self.fmax, time()-self.time0))

    def check_convergence(self):
        """
        There should be two options to terminate the optimization before it reaches the maximum number of cycles
        1, by energy difference compared to the previous step: self.e_tol
        2, by the residual force: self.f_tol
        Will implement both options later.
        """
        converged = False
        if self.fmax < self.f_tol:
            if self.nsteps > 0:
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
            if len(new_coords) == 0:
                new_coords = wp_coords
            else:
                new_coords = np.vstack([new_coords, wp_coords])
            start_index += ws.multiplicity
        return new_coords


from pyxtal.crystal import random_crystal
crystal = random_crystal(19, ['C'], [4], 1.0)
crystal.lattice_matrix = 4.641*np.eye(3)
crystal.frac_coords = np.array([[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])
crystal.cart_coords = np.dot(crystal.frac_coords, crystal.lattice_matrix)
if crystal.valid:
    test = LJ(crystal, epsilon=0.01, sigma=3.40, rcut=8.0)
    print('energy:', test.energy)
    print('enthalpy:', test.enthalpy)
    print('stress:', np.diag(test.stress))

    dyn1 = FIRE(crystal, f_tol=1e-3, dt=0.2, maxmove=2.0)
    dyn1.run(500)
    print('energy: ', LJ(crystal, epsilon=0.01, sigma=3.40).energy)
    print('enthalpy: ', LJ(crystal, epsilon=0.01, sigma=3.40).enthalpy)
    print('stress', np.diag(LJ(crystal, epsilon=0.01, sigma=3.40).stress))
    print('lattice:', np.diag(crystal.lattice_matrix))
    #d2, d = get_neighbors(crystal, 0, 3.0)
    #print(np.sqrt(d2))
