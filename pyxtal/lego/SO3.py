from __future__ import division
import numpy as np
from scipy.special import sph_harm, spherical_in
from ase import Atoms
from ase.neighborlist import NeighborList


class SO3:
    '''
    A class to generate the SO3 power spectrum components
    based off of the Gaussian atomic neighbor density function
    defined in "On Representing Atomic Environments".

    args:
        nmax: int, degree of radial expansion
        lmax: int, degree of spherical harmonic expansion
        rcut: float, cutoff radius for neighbor calculation
        alpha: float, gaussian width parameter
        weight_on: bool, if True, the neighbors with different type will be counted as negative
    '''

    def __init__(self, nmax=3, lmax=3, rcut=3.5, alpha=2.0,
                 weight_on=False):
        # populate attributes
        self.nmax = nmax
        self.lmax = lmax
        self.rcut = rcut
        self.alpha = alpha
        self._type = "SO3"
        self.cutoff_function = 'cosine'
        self.weight_on = weight_on
        self.ncoefs = self.nmax*(self.nmax+1)//2*(self.lmax+1)
        self.tril_indices = np.tril_indices(self.nmax, k=0)
        self.ls = np.arange(self.lmax+1)
        self.norm = np.sqrt(2*np.sqrt(2)*np.pi/np.sqrt(2*self.ls+1))
        self.keys = ['keys', '_nmax', '_lmax', '_rcut', '_alpha',
                     '_cutoff_function', 'weight_on',
                     'ncoefs', 'ls', 'norm', 'tril_indices', 'args']
        self.args = (self.nmax, self.lmax, self.rcut, self.alpha, self._cutoff_function)

    def __str__(self):
        s = "SO3 descriptor with Cutoff: {:6.3f}".format(self.rcut)
        s += " lmax: {:d}, nmax: {:d}, alpha: {:.3f}\n".format(self.lmax, self.nmax, self.alpha)
        return s

    def __repr__(self):
        return str(self)

    def load_from_dict(self, dict0):
        self.nmax = dict0["nmax"]
        self.lmax = dict0["lmax"]
        self.rcut = dict0["rcut"]
        self.alpha = dict0["alpha"]
        self.derivative = dict0["derivative"]

    def save_dict(self):
        """
        save the model as a dictionary in json
        """
        dict = {"nmax": self.nmax,
                "lmax": self.lmax,
                "rcut": self.rcut,
                "alpha": self.alpha,
                "derivative": self.derivative,
                "_type": "SO3",
               }
        return dict

    @property
    def nmax(self):
        return self._nmax

    @nmax.setter
    def nmax(self, nmax):
        if isinstance(nmax, int) is True:
            if nmax < 1:
                raise ValueError('nmax must be greater than or equal to 1')
            if nmax > 11:
                raise ValueError('nmax > 11 yields complex eigenvalues')
            self._nmax = nmax
        else:
            raise ValueError('nmax must be an integer')

    @property
    def lmax(self):
        return self._lmax

    @lmax.setter
    def lmax(self, lmax):
        if isinstance(lmax, int) is True:
            if lmax < 0:
                raise ValueError('lmax must be greater than or equal to zero')
            elif lmax > 32:
                raise NotImplementedError('support a maxmimum l=32 for spherical harmonics')
            self._lmax = lmax
        else:
            raise ValueError('lmax must be an integer')

    @property
    def rcut(self):
        return self._rcut

    @rcut.setter
    def rcut(self, rcut):
        if isinstance(rcut, float) or isinstance(rcut, int):
            if rcut <= 0:
                raise ValueError('rcut must be greater than zero')
            self._rcut = rcut
        else:
            raise ValueError('rcut must be a float')

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, alpha):
        if isinstance(alpha, float) or isinstance(alpha, int):
            if alpha <= 0:
                raise ValueError('alpha must be greater than zero')
            self._alpha = alpha
        else:
            raise ValueError('alpha must be a float')

    @property
    def cutoff_function(self):
        return self._cutoff_function

    @cutoff_function.setter
    def cutoff_function(self, cutoff_function):
        self._cutoff_function = Cosine

    def clear_memory(self):
        '''
        Clears all non essential attributes for the calculator
        '''
        attrs = list(vars(self).keys())
        for attr in attrs:
            if attr not in self.keys:
                delattr(self, attr)
        return

    def init_atoms(self, atoms, atom_ids=None):
        """
        initilize atoms related attributes
        """
        self._atoms = atoms
        self.natoms = len(atoms)
        self.build_neighbor_list(atom_ids)

    def compute_p(self, atoms, atom_ids=None, return_CN=False):
        """
        Compute the powerspectrum function

        Args:
            atoms: ase atoms object
            atom_ids: optional list of atomic indices

        Returns:
            p array (N, M)
        """

        if atom_ids is None: atom_ids = range(len(atoms))
        self.init_atoms(atoms, atom_ids)
        plist = np.zeros((len(atom_ids), self.ncoefs), dtype=np.float64)
        dist = np.zeros((len(atom_ids), 2))

        if len(self.neighborlist) > 0:
            cs = compute_cs(self.neighborlist, *self.args)
            cs *= self.atomic_weights[:, np.newaxis, np.newaxis, np.newaxis]
            cs = np.einsum('inlm,l->inlm', cs, self.norm)

            # Get r_ij and compute C*np.conj(C)
            for _i, i in enumerate(atom_ids):
                centers = self.neighbor_indices[:, 0] == i
                CN = len(self.neighborlist[centers])
                if CN > 0:
                    ctot = cs[centers].sum(axis=0)
                    P = np.einsum('ijk,ljk->ilj', ctot, np.conj(ctot)).real
                    plist[_i] = P[self.tril_indices].flatten()
                    dist[_i, 0] += np.linalg.norm(self.neighborlist[centers], axis=1).max()
                    dist[_i, 1] += CN
        if return_CN:
            return plist, dist
        else:
            return plist


    def compute_dpdr(self, atoms, atom_ids=None):
        """
        Compute the powerspectrum function

        Args:
            atoms: ase atoms object
            atom_ids: optional list of atomic indices

        Returns:
            dpdr array (N, N, M, 3) and p array (N, M)
        """

        if atom_ids is None: atom_ids = range(len(atoms))
        self.init_atoms(atoms, atom_ids)
        p_list = np.zeros((len(atom_ids), self.ncoefs), dtype=np.float64)
        dp_list = np.zeros((len(atom_ids), self.natoms, self.ncoefs, 3), dtype=np.float64)

        if len(self.neighborlist) > 0:
            # get expansion coefficients and derivatives
            cs, dcs = compute_dcs(self.neighborlist, *self.args)

            # weight cs and dcs
            cs *= self.atomic_weights[:, np.newaxis, np.newaxis, np.newaxis]
            dcs *= self.atomic_weights[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
            cs = np.einsum('inlm,l->inlm', cs, self.norm)
            dcs = np.einsum('inlmj,l->inlmj', dcs, self.norm)
            #print('cs, dcs', self.neighbor_indices, cs.shape, dcs.shape)

            # cs: (N_ij, n, l, m)     => P (N_i, N_des)
            # dcs: (N_ij, n, l, m, 3) => dP (N_i, N_j, N_des, 3)
            # (n, l, m) needs to be merged to 1 dimension

            for _i, i in enumerate(atom_ids):
                # find atoms for which i is the center
                centers = self.neighbor_indices[:, 0] == i

                if len(self.neighbor_indices[centers]) > 0:
                    # total up the c array for the center atom
                    ctot = cs[centers].sum(axis=0) #(n, l, m)

                    # power spectrum P = c*c_conj
                    # eq_3 (n, n', l) eliminate m
                    P = np.einsum('ijk, ljk->ilj', ctot, np.conj(ctot)).real
                    p_list[_i] += P[self.tril_indices].flatten()

                    # gradient of P for each neighbor, eq_26
                    # (N_ijs, n, n', l, 3)
                    # dc * c_conj + c * dc_conj
                    dP = np.einsum('wijkn,ljk->wiljn', dcs[centers], np.conj(ctot))
                    dP += np.conj(np.transpose(dP, axes=[0, 2, 1, 3, 4]))
                    dP = dP.real

                    #print("shape of P/dP", P.shape, dP.shape)#; import sys; sys.exit()

                    # QZ: to check
                    ijs = self.neighbor_indices[centers]
                    for _id, j in enumerate(ijs[:, 1]):
                        tmp = dP[_id][self.tril_indices].flatten().reshape(self.ncoefs, 3)
                        dp_list[_i, j, :, :] += tmp
                        dp_list[_i, i, :, :] -= tmp

        return dp_list, p_list

    def compute_dpdr_5d(self, atoms, atom_ids=None):
        """
        Compute the powerspectrum function with respect to supercell

        Args:
            atoms: ase atoms object

        Returns:
            dpdr array (N, N, M, 3, 27) and p array (N, M)
        """
        if atom_ids is None: atom_ids = range(len(atoms))
        self.init_atoms(atoms, atom_ids)
        p_list = np.zeros((len(atom_ids), self.ncoefs), dtype=np.float64)
        dp_list = np.zeros((len(atom_ids), self.natoms, self.ncoefs, 3, 27), dtype=np.float64)

        if len(self.neighborlist) > 0:
            # get expansion coefficients and derivatives
            cs, dcs = compute_dcs(self.neighborlist, *self.args)

            # weight cs and dcs
            cs *= self.atomic_weights[:, np.newaxis, np.newaxis, np.newaxis]
            dcs *= self.atomic_weights[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
            cs = np.einsum('inlm,l->inlm', cs, self.norm)
            dcs = np.einsum('inlmj,l->inlmj', dcs, self.norm)
            #print('cs, dcs', self.neighbor_indices, cs.shape, dcs.shape)

            # cs: (N_ij, n, l, m)     => P (N_i, N_des)
            # dcs: (N_ij, n, l, m, 3) => dP (N_i, N_j, N_des, 3)
            # (n, l, m) needs to be merged to 1 dimension
            neigh_ids = np.arange(len(self.neighbor_indices))
            for _i, i in enumerate(atom_ids):
                # find atoms for which i is the center
                pair_ids = neigh_ids[self.neighbor_indices[:, 0] == i]
                if len(pair_ids) > 0:
                    ctot = cs[pair_ids].sum(axis=0) #(n, l, m)
                    dctot = dcs[pair_ids].sum(axis=0)
                    # power spectrum P = c * c_conj
                    P = np.einsum('ijk, ljk->ilj', ctot, np.conj(ctot)).real
                    p_list[_i] += P[self.tril_indices].flatten()

                    # loop over each pair
                    for pair_id in pair_ids:
                        (_, j, cell_id) = self.neighbor_indices[pair_id]
                        # map from (x, y, z) to (0, 27)

                        # (N_ijs, n, n', l, 3)
                        # dp = dc * c_conj + c * dc_conj
                        dP = np.einsum('ijkn, ljk->iljn', dcs[pair_id], np.conj(ctot))
                        dP += np.einsum('ijkn, ljk->iljn', np.conj(dcs[pair_id]), ctot)

                        dP = dP.real[self.tril_indices].flatten().reshape(self.ncoefs, 3)
                        #print(cs[pair_id].shape, dcs[pair_id].shape, dP.shape)

                        dp_list[_i, j, :, :, cell_id] += dP
                        dp_list[_i, i, :, :, 13] -= dP

        return dp_list, p_list


    def calculate(self, atoms, atom_ids=None, derivative=False):
        '''
        API for Calculating the SO(3) power spectrum components of the
        smoothened atomic neighbor density function

        Args:
            atoms: an ASE atoms object corresponding to the desired
                   atomic arrangement
            derivative: bool, whether to calculate the gradient of not
        '''
        p_list = None
        dp_list = None
        if derivative:
            dp_list, p_list = self.compute_dpdr(atoms, atom_ids)
        else:
            p_list = self.compute_p(atoms, atom_ids)

        x = {'x': p_list,
             'dxdr': dp_list,
             'elements': list(atoms.symbols)}
        self.clear_memory()
        return x


    def build_neighbor_list(self, atom_ids=None):
        '''
        Builds a neighborlist for the calculation of bispectrum components for
        a given ASE atoms object given in the calculate method.
        '''
        atoms = self._atoms
        cell_matrix = atoms.get_cell()
        VECTORS = np.array([[x1, y1, z1] for x1 in range(-1, 2) for y1 in range(-1, 2) for z1 in range(-1, 2)])
        neighbors = []
        neighbor_indices = []
        atomic_weights = []

        #if True: #atom_ids is None:
        if atom_ids is None:
            atom_ids = range(len(atoms))

            cutoffs = [self.rcut/2]*len(atoms)
            nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
            nl.update(atoms)

            for i in atom_ids:
                # get center atom position vector
                center_atom = atoms.positions[i]
                # get indices and cell offsets for each neighbor
                indices, offsets = nl.get_neighbors(i)
                #print(indices); import sys; sys.exit()

                for j, offset in zip(indices, offsets):
                    (x, y, z) = offset
                    cell_id = (x+1) * 9 + (y+1) * 3 + z + 1
                    pos = atoms.positions[j] + offset@cell_matrix - center_atom
                    if np.sum(np.abs(pos)) < 1e-3:
                        # to skip self
                        if cell_id == 13:
                            continue
                        # to prevent division by zero
                        else:
                            pos += 1e-3
                    neighbors.append(pos)
                    if self.weight_on and atoms[j].number != atoms[i].number:
                        factor = -1
                    else:
                        factor = 1
                    atomic_weights.append(factor*atoms[j].number)
                    neighbor_indices.append([i, j, cell_id])
        else:
            # A short cut version if we only compute the neighbors for a few atoms
            ref_pos = np.repeat(atoms.positions[:, :, np.newaxis], 27, axis=2)
            cell_shifts = np.dot(VECTORS, cell_matrix)
            ref_pos += cell_shifts.T[np.newaxis, :, :]
            ref_ids = np.arange(len(atoms))
            for i in atom_ids:
                center_atom = atoms.positions[i]

                # Compute distances relative to the center_atom across all cells
                dists = np.linalg.norm(ref_pos - center_atom[:, np.newaxis], axis=1)

                # Create mask for distances that are within the cutoff and greater than 1e-2
                mask = (dists > 1e-2) & (dists < self.rcut)

                # Use the mask to find neighbors in all cells at once
                valid_atoms, valid_cells = np.where(mask)
                # Append the neighbors and weights
                for j, cell_id in zip(valid_atoms, valid_cells):
                    neighbors.append(ref_pos[j, :, cell_id] - center_atom)
                    factor = -1 if self.weight_on and atoms[j].number != atoms[i].number else 1
                    atomic_weights.append(factor * atoms[j].number)
                    neighbor_indices.append([i, j, cell_id])
                #for cell_id in range(27):
                #    dists = np.linalg.norm(ref_pos[:, :, cell_id] - center_atom, axis=1)
                #    mask = (dists > 1e-2) & (dists < self.rcut)
                #    for j in ref_ids[mask]:
                #        neighbors.append(ref_pos[j, :, cell_id] - center_atom)
                #        factor = -1 if self.weight_on and atoms[j].number != atoms[i].number else 1
                #        atomic_weights.append(factor*atoms[j].number)
                #        neighbor_indices.append([i, j, cell_id])

        self.neighborlist = np.array(neighbors, dtype=np.float64)
        self.atomic_weights = np.array(atomic_weights, dtype=np.int64)
        self.neighbor_indices = np.array(neighbor_indices, dtype=np.int64)

def Cosine(Rij, Rc, derivative=False):
    # Rij is the norm
    if not derivative:
        result = 0.5 * (np.cos(np.pi * Rij / Rc) + 1.)
    else:
        result = -0.5 * np.pi / Rc * np.sin(np.pi * Rij / Rc)
    return result

def W(nmax):
    arr = np.zeros((nmax,nmax), np.float64)
    for alpha in range(1, nmax+1, 1):
        temp1 = (2*alpha+5)*(2*alpha+6)*(2*alpha+7)
        for beta in range(1, alpha+1, 1):
            temp2 = (2*beta+5)*(2*beta+6)*(2*beta+7)
            arr[alpha-1, beta-1] = np.sqrt(temp1*temp2)/(5+alpha+beta)/(6+alpha+beta)/(7+alpha+beta)
            arr[beta-1, alpha-1] = arr[alpha-1, beta-1]

    sinv = np.linalg.inv(arr)
    eigvals, V = np.linalg.eig(sinv)
    sqrtD = np.diag(np.sqrt(eigvals))
    arr[:,:] = np.dot(np.dot(V, sqrtD), np.linalg.inv(V))
    return arr

def phi(r, alpha, rcut):
    '''
    See g below
    '''
    return (rcut-r)**(alpha+2)/np.sqrt(2*rcut**(2*alpha+7)/(2*alpha+5)/(2*alpha+6)/(2*alpha+7))

def g(r, n, nmax, rcut, w):

    Sum = 0.0
    for alpha in range(1, nmax+1):
        Sum += w[n-1, alpha-1]*phi(r, alpha, rcut)

    return Sum

def GaussChebyshevQuadrature(nmax, lmax):
    NQuad = (nmax+lmax+1) * 10 #
    quad_array = np.zeros(NQuad, dtype=np.float64)
    for i in range(1, NQuad+1):
        # roots of Chebyshev polynomial of degree N
        x = np.cos((2*i-1)*np.pi/2/NQuad)
        quad_array[i-1] = x
    return quad_array, np.pi/NQuad

def compute_cs(pos, nmax, lmax, rcut, alpha, cutoff):
    """
    Compute expansion coefficients for a system based on the input positions.

    This function calculates the expansion coefficients for a set of atomic positions using
    Gauss-Chebyshev quadrature, spherical Bessel functions, and spherical harmonics. It is
    typically used in models that require high-dimensional projections of atomic environments.

    Args:
        pos (numpy.ndarray): An array of atomic positions (N x 3) where N is the number of atoms.
        nmax (int): Maximum radial quantum number used in the expansion.
        lmax (int): Maximum angular momentum quantum number used in the expansion.
        rcut (float): Cutoff radius for interactions and the radial expansion.
        alpha (float): Gaussian decay factor applied to the radial functions.
        cutoff (callable): A function to compute cutoff values for the radial distances.

    Returns:
        numpy.ndarray: A 4D array of expansion coefficients with shape (N_neighbors, nmax, lmax+1, 2*lmax+1),
        where `N_neighbors` is the number of neighbor atoms, `nmax` is the number of radial terms, and
        `lmax` and `m` correspond to angular momentum quantum numbers.
    """

    # 1. Init Overlap matrix, distances and quadrature points
    w = W(nmax)
    Ris = np.linalg.norm(pos, axis=1)  # (N_neighbors)
    GCQuadrature, weight = GaussChebyshevQuadrature(nmax, lmax)  # (Nquad)
    weight *= rcut / 2
    Quadrature = rcut / 2 * (GCQuadrature + 1)

    # 2. Rdial G functions
    Gs = np.zeros((nmax, len(Quadrature)), dtype=np.float64)  # (nmax, Nquad)
    for n in range(1, nmax + 1):
        Gs[n - 1, :] = g(Quadrature, n, nmax, rcut, w)

    Quad_Squared = Quadrature ** 2
    Gs *= Quad_Squared * np.exp(-alpha * Quad_Squared) * np.sqrt(1 - GCQuadrature**2) * weight

    # 3. Bessel functions
    BesselArgs = 2 * alpha * np.outer(Ris, Quadrature)  # (N_neighbors x Nquad)
    Bessels = np.zeros((len(Ris), len(Quadrature), lmax + 1), dtype=np.float64)  # (N_neighbors x Nquad x lmax+1)
    for l in range(lmax + 1):
        Bessels[:, :, l] = spherical_in(l, BesselArgs)

    # 4. Integration over radial coordinates using Einstein summation
    integral_array = np.einsum('ij,kjl->kil', Gs, Bessels)  # (N_neighbors x nmax x lmax+1)

    # 5. Gaussian factors for each atom and apply the cutoff function
    exparray = 4 * np.pi * np.exp(-alpha * Ris**2)  # (N_neighbors)
    cutoff_array = cutoff(Ris, rcut)
    np.multiply(exparray, cutoff_array, out=exparray)

    # 6. Compute the spherical harmonics for each atom
    # From Cartesian coordinates to spherical coordinates
    thetas = pos[:, 2] / Ris[:]
    thetas = np.arccos(thetas)
    phis = np.arctan2(pos[:, 1], pos[:, 0])
    msize = 2 * lmax + 1
    ylms = np.zeros((len(Ris), lmax + 1, msize), dtype=np.complex128)

    #for l in range(lmax + 1):
    #    for m in range(-l, l + 1):
    #        midx = msize // 2 + m
    #        ylms[:, l, midx] = sph_harm(m, l, phis, thetas)
    l_vals = np.repeat(np.arange(lmax + 1), [2 * l + 1 for l in range(lmax + 1)])
    m_vals = np.concatenate([np.arange(-l, l + 1) for l in range(lmax + 1)])
    midx_vals = msize // 2 + m_vals
    ylms[:, l_vals, midx_vals] = sph_harm(m_vals, l_vals, phis[:, np.newaxis], thetas[:, np.newaxis])


    # 7. Multiply Ylm and Gaussian factors
    Y_mul_innerprod = np.einsum('ijk,ilj->iljk', ylms, integral_array)
    C = np.einsum('i,ijkl->ijkl', exparray, Y_mul_innerprod)

    return C

def compute_dcs(pos, nmax, lmax, rcut, alpha, cutoff):
    """
    Compute exapnsion coefficients

    Args:
        pos:
        nmax (int):
        lmax (int):
        rcut (float):
        alpha (float):
        cutoff (callable):

    Returns:
        c(N_ij, nmax, lmax+1, 2lmax+1)
        dc(N_ij, nmax, lmax+1, 2lmax+1, 3) for each x,y,z
    """
    # compute the overlap matrix
    w = W(nmax)

    # get the norm of the position vectors
    Ris = np.linalg.norm(pos, axis=1) # (Nneighbors)

    # get unit vectors
    upos = pos/Ris[:,np.newaxis]

    # initialize Gauss Chebyshev Quadrature
    GCQuadrature, weight = GaussChebyshevQuadrature(nmax,lmax) #(Nquad)
    weight *= rcut/2
    # transform from (-1,1) to (0, rcut)
    Quadrature = rcut/2*(GCQuadrature+1)

    # compute the arguments for the bessel functions
    BesselArgs = 2*alpha*np.outer(Ris,Quadrature)#(Nneighbors x Nquad)

    # initalize the arrays for the bessel function values
    # and the G function values
    Bessels = np.zeros((len(Ris), len(Quadrature), lmax+1), dtype=np.float64) #(Nneighbors x Nquad x lmax+1)
    Gs = np.zeros((nmax, len(Quadrature)), dtype=np.float64) # (nmax, nquad)
    dBessels = np.zeros((len(Ris), len(Quadrature), lmax+1), dtype=np.float64) #(Nneighbors x Nquad x lmax+1)

    # compute the g values
    for n in range(1,nmax+1,1):
        Gs[n-1,:] = g(Quadrature, n, nmax, rcut,w)*weight

    # compute the bessel values
    for l in range(lmax+1):
        Bessels[:,:,l] = spherical_in(l, BesselArgs)
        dBessels[:,:,l] = spherical_in(l, BesselArgs, derivative=True)

    #(Nneighbors x Nquad x lmax+1) unit vector here
    gradBessels = np.einsum('ijk,in->ijkn',dBessels,upos)
    gradBessels *= 2*alpha

    # multiply with r for the integral
    gradBessels = np.einsum('ijkn,j->ijkn',gradBessels,Quadrature)

    # mutliply the terms in the integral separate from the Bessels
    Quad_Squared = Quadrature**2
    Gs *= Quad_Squared * np.exp(-alpha*Quad_Squared) * np.sqrt(1-GCQuadrature**2)

    # perform the integration with the Bessels
    integral_array = np.einsum('ij,kjl->kil', Gs, Bessels) # (Nneighbors x nmax x lmax+1)

    grad_integral_array = np.einsum('ij,kjlm->kilm', Gs, gradBessels)# (Nneighbors x nmax x lmax+1, 3)

    # compute the gaussian for each atom
    exparray = 4*np.pi*np.exp(-alpha*Ris**2) # (Nneighbors)

    gradexparray = (-2*alpha*Ris*exparray)[:,np.newaxis]*upos

    cutoff_array = cutoff(Ris, rcut)

    grad_cutoff_array = np.einsum('i,in->in',cutoff(Ris, rcut, True), upos)

    # get the spherical coordinates of each atom
    thetas = np.arccos(pos[:,2]/Ris[:])
    phis = np.arctan2(pos[:,1], pos[:,0])

    # the size changes temporarily for the derivative
    # determine the size of the m axis
    Msize = 2*(lmax+1)+1
    msize = 2*lmax + 1
    # initialize an array for the spherical harmonics and gradients
    #(Nneighbors, l, m, *3*)
    ylms = np.zeros((len(Ris), lmax+1+1, Msize), dtype=np.complex128)
    gradylms = np.zeros((len(Ris), lmax+1, msize, 3), dtype=np.complex128)
    # compute the spherical harmonics
    for l in range(lmax+1+1):
        for m in range(-l,l+1,1):
            midx = Msize//2 + m
            ylms[:,l,midx] = sph_harm(m, l, phis, thetas)


    for l in range(1, lmax+1):
        for m in range(-l, l+1, 1):
            midx = msize//2 + m
            Midx = Msize//2 + m
            # get gradient with recpect to spherical covariant components
            xcov0 = -np.sqrt(((l+1)**2-m**2)/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx]/Ris

            if abs(m) <= l-1:
                xcov0 += np.sqrt((l**2-m**2)/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx]/Ris


            xcovpl1 = -np.sqrt((l+m+1)*(l+m+2)/2/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx+1]/Ris

            if abs(m+1) <= l-1:
                xcovpl1 -= np.sqrt((l-m-1)*(l-m)/2/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx+1]/Ris


            xcovm1 = -np.sqrt((l-m+1)*(l-m+2)/2/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx-1]/Ris

            if abs(m-1) <= l-1:
                xcovm1 -= np.sqrt((l+m-1)*(l+m)/2/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx-1]/Ris

            #transform the gradient to cartesian
            gradylms[:,l,midx,0] = 1/np.sqrt(2)*(xcovm1-xcovpl1)
            gradylms[:,l,midx,1] = 1j/np.sqrt(2)*(xcovm1+xcovpl1)
            gradylms[:,l,midx,2] = xcov0

    # index ylms to get rid of extra terms for derivative
    ylms = ylms[:, 0:lmax+1, 1:1+2*lmax+1]
    # multiply the spherical harmonics and the radial inner product
    Y_mul_innerprod = np.einsum('ijk,ilj->iljk', ylms, integral_array)
    # multiply the gradient of the spherical harmonics with the radial inner get_radial_inner_product
    dY_mul_innerprod = np.einsum('ijkn,ilj->iljkn', gradylms, integral_array)
    # multiply the spherical harmonics with the gradient of the radial inner get_radial_inner_product
    Y_mul_dinnerprod = np.einsum('ijk,iljn->iljkn', ylms, grad_integral_array)
    # multiply the gaussians into the expression with 4pi
    C = np.einsum('i,ijkl->ijkl', exparray, Y_mul_innerprod)
    # multiply the gradient of the gaussian with the other terms
    gradexp_mul_y_inner = np.einsum('in,ijkl->ijkln', gradexparray, Y_mul_innerprod)
    # add gradient of inner product and spherical harmonic terms
    gradHarmonics_mul_gaussian = np.einsum('ijkln,i->ijkln', dY_mul_innerprod+Y_mul_dinnerprod, exparray)
    dC = gradexp_mul_y_inner + gradHarmonics_mul_gaussian
    dC *= cutoff_array[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
    dC += np.einsum('in,ijkl->ijkln', grad_cutoff_array, C)
    C *= cutoff_array[:, np.newaxis, np.newaxis, np.newaxis]
    return C, dC

if  __name__ == "__main__":
    from optparse import OptionParser
    from ase.io import read
    import time
    # ---------------------- Options ------------------------
    parser = OptionParser()
    parser.add_option("-c", "--crystal", dest="structure",
                      help="crystal from file, cif or poscar, REQUIRED",
                      metavar="crystal")

    parser.add_option("-r", "--rcut", dest="rcut", default=3.0, type=float,
                      help="cutoff for neighbor calcs, default: 3.0"
                      )

    parser.add_option("-l", "--lmax", dest="lmax", default=2, type=int,
                      help="lmax, default: 2"
                      )

    parser.add_option("-n", "--nmax", dest="nmax", default=2, type=int,
                      help="nmax, default: 2"
                      )

    parser.add_option("-a", "--alpha", dest="alpha", default=2.0, type=float,
                      help="cutoff for neighbor calcs, default: 2.0"
                      )

    parser.add_option("-f", dest="der", default=True,
                      action='store_false',help='derivative flag')

    (options, args) = parser.parse_args()

    if options.structure is None:
        from ase.build import bulk
        test = bulk('Si', 'diamond', a=5.459, cubic=True)
    else:
        test = read(options.structure, format='vasp')

    cell = test.get_cell()
    cell[0,1] += 0.5
    test.set_cell(cell)

    lmax = options.lmax
    nmax = options.nmax
    rcut = options.rcut
    alpha = options.alpha
    der = options.der

    start1 = time.time()
    f = SO3(nmax=nmax, lmax=lmax, rcut=rcut, alpha=alpha)
    start2 = time.time()
    print(f)
    p = f.compute_p(test, atom_ids=[0]); print('from P', p)
    x = f.calculate(test, derivative=der)
    print('x', x['x'])
    dp, p = f.compute_dpdr(test); print('from dP', p)
    dp, p = f.compute_dpdr_5d(test); print('from dP5d', p)

    (x, spg, wps) = ([13.45493847, 0.98, 0.04, 0.86, 0.08, 0.94, 0.14, 0.22, 0.66, 0.22, 0.74, 0.28, 0.2, 0.18, 0.68], 226, ['192j', '192j', '192j', '96i', '192j'])

    from pyxtal import pyxtal
    xtal = pyxtal()
    xtal.from_spg_wps_rep(spg, wps, x, ['C']*len(wps))
    atoms = xtal.to_ase()
    p1 = f.compute_p(atoms)[0]; print("from ase", p1)
    p2 = f.compute_p(atoms, atom_ids=[0])[0]; print("from self", p2)
    print(np.sum((p1-p2)**2))
    assert(np.sum((p1-p2)**2)<1e-4)

