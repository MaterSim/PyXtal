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
    """
    LJ model for 3D crystals, maybe extended to 0D, 1D, 2D later
    """
    def __init__(self, epsilon=1.0, sigma=1.0, rcut=8.0):
        """
        passing the parameter to LJ model
        - epsilon
        - sigma
        - rcut
        """

        self.epsilon = epsilon
        self.sigma = sigma
        self.rcut = rcut

    def calc(self, struc, press=1e-4):
        pos = struc.cart_coords
        lat = struc.lattice_matrix
        volume = np.linalg.det(lat)

        # initiate the energy, force, stress
        energy = 0
        force = np.zeros([len(pos), 3])
        stress = np.zeros([3, 3])
        sigma6 = self.sigma**6
        sigma12 = sigma6*sigma6

        # calculating the energy/force, needs to update sigma/epsilons
        for i, pos0 in enumerate(pos):
            [r2, r] = get_neighbors(struc, i, self.rcut)
            r6 = np.power(r2, 3)     
            r12 = np.power(r6, 2)
            energy += np.sum(4.0*self.epsilon*(sigma12/r12 - sigma6/r6))
            f = (24*self.epsilon*(2.0*sigma12/r12-sigma6/r6)/r2)[:, np.newaxis] * r
            force[i] = f.sum(axis=0)
            stress += np.dot(f.T, r)

        energy = 0.5*energy    
        enthalpy = energy + press*volume*GPa2eV
        force = force
        stress = -0.5*stress/volume*eV2GPa

        return energy, enthalpy, force, stress


class FIRE():
    """
    FIRE algorithm for structure optimization
    """
    def __init__(self, struc, model, symmetrize=False, e_tol=1e-5, f_tol=1e-2, 
                 dt=0.1, maxmove=0.1, dtmax=1.0, Nmin=10, finc=1.1, fdec=0.5,
                 astart=0.1, fa=0.99, a=0.1):
        """
        Parameters for FIRE:
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
        self.time0 = time()
        self.model = model
        self.struc = struc
        self.initialize()

    def initialize(self):
        """
        initilaize the positions, energy, force for the 1st step
        """
        self.energy, self.enthalpy, self.force, self.stress = self.model.calc(self.struc)
        self.fmax = np.max(np.abs(np.vstack((self.stress, self.force)).flatten()))
        self.v = np.zeros((len(self.force)+3, 3))
        #print('Initial Energy: {:12.4f}'.format(self.energy))

    def update(self, freq=50):
        self.energy, self.enthalpy, self.force, self.stress = self.model.calc(self.struc)
        self.fmax = np.max(np.abs(np.vstack((self.stress, self.force)).flatten()))
        if (self.nsteps % freq == 0):
            print('Step {0:4d} Enthalpy: {1:12.4f} Fmax: {2:12.4f} Vol: {3:6.2f}'.
                format(self.nsteps, self.enthalpy, self.fmax, self.volume))

    def step(self):
        f = np.vstack((self.stress, self.force))
        pos = np.vstack((self.struc.lattice_matrix, self.struc.cart_coords))

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
            # f[:3, :] is the gradient for force, need to symmetrize it as well
            # f[:3, :] = 
            f[3:, :] = self.symmetrized_coords(f[3:, :])
        #f[0,1:] = 0
        #f[1,0] = 0
        #f[1,2] = 0
        #f[2,:2] = 0

        self.v += self.dt * f
        dr = self.dt * self.v  #needs to impose constraints
        normdr = np.sqrt(np.vdot(dr, dr))
        if normdr > self.maxmove:
            dr = self.maxmove * dr / normdr

        #print(self.force)
        #print(self.v)
        pos = pos - dr
        self.struc.lattice_matrix = pos[:3, :]
        self.struc.cart_coords = pos[3:, :]
        self.struc.frac_coords = np.dot(pos[3:, :], np.linalg.inv(pos[:3, :]))
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

        print('Finish at Step {0:4d} Enthalpy: {1:12.4f} Fmax: {2:12.4f} Time: {3:6.2f} seconds'.
                format(self.nsteps, self.enthalpy, self.fmax, time()-self.time0))

    def check_convergence(self):
        """
        There should be two options to terminate the optimization 
        before it reaches the maximum number of cycles
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
from spglib import get_symmetry_dataset

for i in range(10):
    crystal = random_crystal(19, ['C'], [4], 1.0)
    if crystal.valid:
        test = LJ(epsilon=0.01, sigma=3.40, rcut=8.0)
        struc = (crystal.lattice_matrix, crystal.frac_coords, [6]*4)
        eng, enth, force, stress = test.calc(crystal)
        sg =  get_symmetry_dataset(struc)['number']
        print('\nBefore relaxation Space group: {:4d}  Energy: {:12.4}  Enthalpy: {:12.4}\n'.format(sg, eng, enth))
    
        dyn1 = FIRE(crystal, test, f_tol=1e-5, dt=0.2, maxmove=0.2)
        dyn1.run(500)
        eng, enth, force, stress = test.calc(crystal)
        struc = (dyn1.struc.lattice_matrix, dyn1.struc.frac_coords, [6]*4)
        sg =  get_symmetry_dataset(struc, symprec=0.1)['number']
        print('\nAfter relaxation  Space group: {:4d}  Energy: {:12.4}  Enthalpy: {:12.4}'.format(sg, eng, enth))
        # right now, it seems structures goes to either HCP of FCC after relaxation, which is expected for 3D LJ system
        # need to compare with other code to see if the energy is correct
