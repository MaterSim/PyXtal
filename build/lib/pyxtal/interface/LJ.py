import numpy as np
import sys
from time import time
from pymatgen.core.lattice import Lattice
from pyxtal.operations import *
import logging

logging.basicConfig(filename='test.log', level=logging.DEBUG)
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
            logging.debug('Step {0:4d} Enthalpy: {1:12.4f} Fmax: {2:12.4f} Vol: {3:6.2f}'.
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

        self.v += self.dt * f
        dr = self.dt * self.v  #needs to impose constraints
        normdr = np.sqrt(np.vdot(dr, dr))
        if normdr > self.maxmove:
            dr = self.maxmove * dr / normdr
        #print('frac0', np.dot(pos[3:,:], np.linalg.inv(pos[:3,:])))

        #Symmetrize the force
        if self.symmetrize:
            # f[:3, :] is the gradient for force, need to symmetrize it as well
            dr[:3, :] = self.symmetrized_stress(dr[:3, :])
            dr[3:, :] = np.dot(self.struc.frac_coords, dr[:3, :]) 
            f_frac = np.dot(dr[3:,:], np.linalg.inv(pos[:3,:] - dr[:3, :]))
            f_frac = self.symmetrized_force(f_frac)
            dr[3:,:] += np.dot(f_frac, pos[:3,:] - dr[:3, :])

        #pos = pos - dr
        #Symmetrize positions
        #if self.symmetrize:
        #    pos = self.symmetrized_coords(pos)

        #print(self.force)
        #print(self.v)
        self.struc.lattice_matrix = pos[:3, :] - dr[:3, :]
        self.struc.cart_coords = pos[3:, :] - dr[3:, :]
        self.struc.frac_coords = np.dot(self.struc.cart_coords, np.linalg.inv(self.struc.lattice_matrix))
        self.volume = np.linalg.det(self.struc.lattice_matrix)
   
        #sg = get_symmetry_dataset((pos[:3, :], self.struc.frac_coords, [6]*4), symprec=0.1)['number']
        #print(sg)
        #if self.symmetrize and sg <19:
        #    print(sg)
        #    print('dr\n', dr)
        #    print('pos\n', pos)
        #    print('frac\n', self.struc.frac_coords)
        #    print('force\n', f_frac)
        #    sys.exit()
        self.update()

    def run(self, max_steps=1000):
        self.max_steps = max_steps
        while self.nsteps < max_steps:
            if self.check_convergence():
                break
            self.step()
            self.nsteps += 1

        logging.info('Finish at Step {0:4d} Enthalpy: {1:12.4f} Fmax: {2:12.4f} Time: {3:6.2f} seconds'.
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
        start_index = 0
        new_coords = []
        for ws in self.struc.wyckoff_sites:
            #Get coordinates associated with WP
            original_coords = coords[start_index:start_index+ws.multiplicity]
            #Re-generate the points from the first generating point
            gen_coord = coords[start_index]
            wp_coords0 = apply_ops(gen_coord, ws.wp.generators_m)
            #Calculate the PBC translation applied to each point
            translations = original_coords - wp_coords0
            frac_translations = np.dot(translations, np.linalg.inv(self.struc.lattice_matrix))
            frac_translations = np.round(frac_translations)
            cart_translations = np.dot(frac_translations, self.struc.lattice_matrix)
            #Subtract translations from original coords
            translated_coords = original_coords - cart_translations
            #Apply inverse WP operations to get generating_points
            inverse_ops = ws.wp.inverse_generators_m
            generating_points = apply_ops_diagonal(translated_coords, inverse_ops)
            #Find the average of the generating points
            average_point = np.sum(generating_points, axis=0) / ws.multiplicity
            #Convert to fractional coordinates
            average_point = np.dot(average_point, np.linalg.inv(self.struc.lattice_matrix))
            #Apply the WP operations to the averaged generating point
            wp_coords = apply_ops(average_point, ws.wp)
            #Convert back to Cartesian coordintes
            wp_coords = np.dot(wp_coords, self.struc.lattice_matrix)
            #Re-apply the cartesian translations
            wp_coords = wp_coords + cart_translations
            if len(new_coords) == 0:
                new_coords = wp_coords
            else:
                new_coords = np.vstack([new_coords, wp_coords])
            start_index += ws.multiplicity
        return new_coords

    def symmetrized_force(self, coords):
        start_index = 0
        new_coords = []
        for ws in self.struc.wyckoff_sites:
            #Get coordinates associated with WP
            original_coords = coords[start_index:start_index+ws.multiplicity]
            #Re-generate the points from the first generating point
            gen_coord = coords[start_index]
            wp_coords0 = apply_ops(gen_coord, ws.wp.generators_m)
            #Apply inverse WP operations to get generating_points
            inverse_ops = ws.wp.inverse_generators_m
            generating_points = apply_ops_diagonal(wp_coords0, inverse_ops)
            #Find the average of the generating points
            average_point = np.sum(generating_points, axis=0) / ws.multiplicity
            #Convert to fractional coordinates
            average_point = np.dot(average_point, np.linalg.inv(self.struc.lattice_matrix))
            #For forces, we do not apply translational WP operations, so we truncate the ops
            matrices = np.array([op.affine_matrix[:3,:3] for op in ws.wp])
            #Apply the truncated WP operations to the averaged generating point
            wp_coords = np.dot(average_point, matrices)
            #Convert back to Cartesian coordintes
            wp_coords = np.dot(wp_coords, self.struc.lattice_matrix)
            if len(new_coords) == 0:
                new_coords = wp_coords
            else:
                new_coords = np.vstack([new_coords, wp_coords])
            start_index += ws.multiplicity
        #if np.sum(np.abs(coords-new_coords))>1e-1:
        #    print(coords)
        #    print(new_coords)
        #    print(coords-new_coords)
        #    sys.exit()
        return new_coords

    def symmetrized_stress(self, stress):
        #Make the matrix lower-diagonal
        for (i,j) in [(1,0),(2,1),(2,2)]:
            stress[i][j] = stress[i][j] + stress[j][i]
            stress[j][i] = 0
        #Normalize the matrix
        snm = self.struc.lattice.stress_normalization_matrix
        m2 = np.multiply(stress, snm)
        #Normalize the on-diagonal elements
        indices = self.struc.lattice.stress_indices
        if len(indices) == 2:
            total = 0
            for index in indices:
                total += stress[index]        
            value = np.sqrt(total)
            for index in inices:
                m2[index] = value
        elif len(indices) == 3:
            total = 0
            for index in indices:
                total += stress[index]        
            value = np.cbrt(total)
            for index in inices:
                m2[index] = value
        return m2

from pyxtal.crystal import random_crystal
from spglib import get_symmetry_dataset

for i in range(10):
    crystal = random_crystal(11, ['C'], [4], 1.0)
    if crystal.valid:
        crystal1 = deepcopy(crystal)
        test = LJ(epsilon=0.01, sigma=3.40, rcut=8.0)
        struc = (crystal1.lattice_matrix, crystal1.frac_coords, [6]*4)
        eng, enth, force, stress = test.calc(crystal1)
        sg =  get_symmetry_dataset(struc)['number']
        print('\nBefore relaxation Space group:            {:4d}   Energy: {:12.4}  Enthalpy: {:12.4}'.format(sg, eng, enth))
    
        dyn1 = FIRE(crystal1, test, f_tol=1e-5, dt=0.2, maxmove=0.2) #, symmetrize=True)
        dyn1.run(500)
        eng, enth, force, stress = test.calc(crystal1)
        struc = (dyn1.struc.lattice_matrix, dyn1.struc.frac_coords, [6]*4)
        sg =  get_symmetry_dataset(struc, symprec=0.1)['number']
        print('After relaxation without symm Space group: {:4d}  Energy: {:12.4}  Enthalpy: {:12.4}'.format(sg, eng, enth))

        dyn1 = FIRE(crystal, test, f_tol=1e-5, dt=0.2, maxmove=0.2, symmetrize=True)
        dyn1.run(500)
        eng, enth, force, stress = test.calc(crystal)
        struc = (dyn1.struc.lattice_matrix, dyn1.struc.frac_coords, [6]*4)
        sg =  get_symmetry_dataset(struc, symprec=0.02)['number']
        if sg is None:
            sg = 0
        print('After relaxation with symm Space group:    {:4d}  Energy: {:12.4}  Enthalpy: {:12.4}'.format(sg, eng, enth))
        #dyn1 = FIRE(crystal, test, f_tol=1e-5, dt=0.2, maxmove=0.2)
        #struc = (dyn1.struc.lattice_matrix, dyn1.struc.frac_coords, [6]*4)
        #sg =  get_symmetry_dataset(struc, symprec=0.02)['number']
        #print('After relaxation without symm Space group: {:4d}  Energy: {:12.4}  Enthalpy: {:12.4}'.format(sg, eng, enth))
        # right now, it seems structures goes to either HCP of FCC after relaxation, which is expected for 3D LJ system
        # need to compare with other code to see if the energy is correct
