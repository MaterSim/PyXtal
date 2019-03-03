import numpy as np
from time import time
from pymatgen.core.lattice import Lattice

unit=6.242   #unit from GPa*m^3 to eV
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

    def __init__(self, struc, press=1e-4, epsilon=1.0, sigma=1.0, rcut=3.0):
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

        # calculating the energy/force, needs to update sigma/epsilons
        for i, pos0 in enumerate(pos):
            [r2, r] = get_neighbors(struc, i, rcut)
            r6 = np.power(r2, 3)     
            r12 = np.power(r6, 2)
            energy += np.sum(4*epsilon*(sigma**12/r12 - sigma**6/r6))
            f = (24*epsilon*(2*sigma**12/r12-sigma**6/r6)/r2)[:, np.newaxis] * r
            force[i] = f.sum(axis=0)
            stress += np.dot(f.T, r)
            ids = np.argsort(r2)
            #print(r[ids])
            #print(len(r), np.linalg.norm(np.dot( (24*(2/r12-1/r6)/r2), r)))

        self.energy = 0.5*energy    
        self.enthalpy = self.energy + press*volume*unit
        self.force = force
        self.stress = -0.5*stress

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
        self.v = np.zeros((len(self.force)+3, 3))
        self.trajectory['energy'].append(self.energy)
        self.trajectory['force'].append(self.force)
        self.trajectory['v'].append(self.v)
        self.trajectory['time'].append(time()-self.time0)
        print('Initial Energy: {:12.4f}'.format(self.energy))

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
        #print(dr)
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
        if np.max(self.force.flatten()) < self.f_tol and np.max(self.stress.flatten()) < self.f_tol :
            if self.nsteps > 0:
                #if self.trajectory['energy'][-1] - self.trajectory['energy'][-2] <1e-3:
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
crystal = random_crystal(225, ['C'], [4], 1.0)
if crystal.valid:
    test = LJ(crystal)
    print('Energy', test.energy)
    print('Enthalpy', test.enthalpy)
    print('Force', test.force)
    print('Stress', test.stress)
    print(np.linalg.det(crystal.lattice_matrix))
    #print(crystal.struct)
    dyn1 = FIRE(crystal, f_tol=1e-3, dt=0.2, maxmove=2.0)
    dyn1.run(500)
    print(LJ(crystal).force)
    print(LJ(crystal).stress)
    print(LJ(crystal).energy)
    print(LJ(crystal).enthalpy)
    print(crystal.lattice_matrix)
