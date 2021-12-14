import numpy as np
from scipy.special import sph_harm

def _qlm(center, neighbors, weights, l=4):
    '''
    Calculates the complex vector associated with an atomic site and 
    one of its neighbors

    Args:
        center: center xyz coordinate
        neighbors: a list of neighboring xyz coordinates
        weights: a list of weights for each neighbor
        l:  free integer quantum number
     Returns:
         q: numpy array(complex128), the complex vector qlm normalized
            by the number of nearest neighbors
    '''
    # initiate variable as a complex number
    q = np.zeros(2*l+1, dtype=np.complex128)
    neighbors_count = len(neighbors)

    for i, m in enumerate(range(-l, l+1)):
        for j, neighbor in enumerate(neighbors):
            # find the position vector of the site/neighbor pair
            r_vec = neighbor - center
            r_mag = np.linalg.norm(r_vec)
            theta = np.arccos(r_vec[2] / r_mag)
            if abs((r_vec[2] / r_mag) - 1.0) < 10.**(-8.):
                theta = 0.0
            elif abs((r_vec[2] / r_mag) + 1.0) < 10.**(-8.):
                theta = np.pi

            # phi
            if r_vec[0] < 0.:
                phi = np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif 0. < r_vec[0] and r_vec[1] < 0.:
                phi = 2 * np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif 0. < r_vec[0] and 0. <= r_vec[1]:
                phi = np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] == 0. and 0. < r_vec[1]:
                phi = 0.5 * np.pi
            elif r_vec[0] == 0. and r_vec[1] < 0.:
                phi = 1.5 * np.pi
            else:
                phi = 0.
            '''
            calculate the spherical harmonic associated with
            the neighbor and add to q
            '''
            q[i] += weights[j]*sph_harm(m, l, phi, theta)
    # normalize by number of neighbors
    return q / neighbors_count


def get_q4_and_q6(center, neighbors, weights=None):
    '''
    Computes the Steinhardt orientation order parameters q4 and q6 

    Args:
        center: center xyz coordinate
        neighbors: a list of neighboring xyz coordinates
        weights: a list of weights for each neighbor
        l:  free integer quantum number

    Returns:
         q: numpy array(complex128), the complex vector qlm normalized
            by the number of nearest neighbors
    '''
    if weights is None:
        weights = np.ones(len(neighbors), dtype=np.float)
    elif weights == 'auto': #scaled by r^6
        ds = np.linalg.norm(neighbors-center, axis=1)
        d2 = ds**6
        inv_a2 = 1/d2
        weights = inv_a2/np.sum(inv_a2) *len(ds) #intensity
    l = 4
    qlms = _qlm(center, neighbors, weights, l)
    dot = float(np.sum(qlms*np.conjugate(qlms)))
    factor = (4 * np.pi) / (2*l + 1)
    q4 = np.sqrt((4 * np.pi)/(2*l+1) * dot)

    l = 6
    qlms = _qlm(center, neighbors, weights, l)
    dot = float(np.sum(qlms*np.conjugate(qlms)))
    q6 = np.sqrt((4 * np.pi)/(2*l+1) * dot)

    return q4, q6


if __name__ == '__main__':

    from ase.build import bulk
    from ase.neighborlist import NeighborList

    data = [
            ('Si', 'diamond'),
            ('Fe', 'bcc'),
            ('Cu', 'fcc'),
            ('Mg', 'hcp'),
           ]
    
    for d in data:
        (element, packing) = d
        crystal = bulk(element, packing)
        center = crystal.get_positions()[0]
        cell = crystal.get_cell()
        for cutoff in [2.4, 2.8, 3.0, 3.6, 4.2, 5.0]:
            rs = [cutoff/2]*len(crystal)
            nl = NeighborList(rs, self_interaction=False, bothways=True, skin=0.0)
            nl.update(crystal)
            indices, offsets = nl.get_neighbors(0)
            if len(indices) > 0:
                neighbors = np.zeros([len(indices),3])
                for i in range(len(indices)):
                    pos, offset = crystal.positions[indices[i]], offsets[i]
                    neighbors[i, :] = pos + np.dot(offset, cell)
                q4, q6 = get_q4_and_q6(center, neighbors, 'auto')
                print(element, packing, cutoff, len(neighbors), q4, q6)
            

#import numpy as np
#
#class descriptor():
#    """
#    A class to handle the crystal descriptor
#
#    Args:
#        xtal:
#        cutoff: 
#        factor:
#        size:
#    """
#
#    def __init__(self, struc, cutoff=5.0, factor=1.5, size=201):
#        self.struc = struc
#        self.cutoff = cutoff
#        self.factor = factor
#        self.size = size
#
#        #initialize x, y
#        self.x = np.linspace(0, cutoff, size)
#        self.y = np.zeros([len(self.struc.mol_sites), size])
#        self.res = cutoff/(size-1)
#        self.compute()
#        #self.smear()
#
#    def compute(self):
#        """
#        compute the intensities
#
#        Args:
#            struc: pyxtal object
#        """
#        for i, site in enumerate(self.struc.mol_sites):
#            coord0, specie0 = site._get_coords_and_species(absolute=True, first=True)
#            res = self.struc.get_neighboring_molecules(i, self.factor, self.cutoff)
#            (min_ds, neighs, comps) = res
#            for neigh in neighs:
#                ds = cdist(coord0, neigh)
#                ds = ds[ds<self.cutoff]
#                d2 = ds**6
#                inv_a2 = 1/d2
#                intensities = inv_a2/np.sum(inv_a2) #intensity
#                for d, intensity in zip(ds, intensities):
#                    j = round(d/self.res)
#                    self.y[i, j] += intensity
#
#    def plot(self, filename=None, figsize=(10,8), smears=[0.1, 0.5, 1.0]):
#
#        import matplotlib.pyplot as plt
#
#        plt.figure(figsize=figsize)
#
#        x_min, x_max = 0, self.cutoff + 0.5
#        plt.gca()
#
#        y_smears = []
#        for s in smears:
#            y_smears.append(self.smear(s))
#
#        for id in range(self.y.shape[0]):
#            label = 'Site_{:d} (Raw)'.format(id)
#            plt.plot(self.x, self.y[id,:], label=label)
#            for j, s in enumerate(smears):
#                label = 'Site_{:d} ({:4.1f})'.format(id, s)
#                plt.plot(self.x, y_smears[j][id,:], label=label)
#        plt.legend()
#        plt.xlim([x_min, x_max])
#        plt.xlabel('R($\AA$)')
#        plt.ylabel('Intensity')
#
#        if filename is None:
#            plt.show()
#        else:
#            plt.savefig(filename)
#            plt.close()
#            
#    def smear(self, sigma=0.8):
#        """
#        Apply Gaussian smearing to spectrum y value.
#        Args:
#            sigma: Std dev for Gaussian smear function
#        """
#        from scipy.ndimage.filters import gaussian_filter1d
#        y = np.zeros([len(self.struc.mol_sites), self.size]) 
#        for id in range(self.y.shape[0]):
#            y[id] = gaussian_filter1d(self.y[id], sigma)
#        return y
#
#
#if __name__ == "__main__":
#
#    from pyxtal import pyxtal
#    from pyxtal.XRD import Similarity
#    from pyxtal.representation import representation
#    from scipy.optimize import minimize
#    import warnings
#    warnings.filterwarnings("ignore")
#
#
#    #for code in ['ACSALA', 'ACSALA13']:
#    #    c = pyxtal(molecular=True)
#    #    c.from_CSD(code)
#    #    des = descriptor(c)
#    #    des.plot(code+'.png')
#    #    print(np.sum(des.y))
#
#    #c1 = pyxtal(molecular=True); c1.from_CSD('ACSALA'); des=descriptor(c1); f1=des.smear(1.0)
#    #c2 = pyxtal(molecular=True); c2.from_CSD('ACSALA13'); des=descriptor(c2); f2=des.smear(1.0)
#
#    #smiles = ['CC(=O)Oc1ccccc1C(O)=O']
#    #rep1 = c1.get_1D_representation(); rep3 = representation(rep1.x, smiles); c3=rep3.to_pyxtal()
#    #rep2 = c2.get_1D_representation(); rep4 = representation(rep2.x, smiles); c4=rep4.to_pyxtal()
#
#    #des=descriptor(c3); f3=des.smear(1.0)    
#    #des=descriptor(c4); f4=des.smear(1.0)    
#
#    #f1 = [des.x, f1[0]]
#    #f2 = [des.x, f2[0]]
#    #f3 = [des.x, f3[0]]
#    #f4 = [des.x, f4[0]]
#
#    #S12 = Similarity(f1, f2, l=0.1); print(S12.value)
#    #S13 = Similarity(f1, f3, l=0.1); print(S13.value)
#    #S23 = Similarity(f2, f3, l=0.1); print(S23.value)
#    #S24 = Similarity(f2, f4, l=0.1); print(S24.value)
#    #S34 = Similarity(f2, f4, l=0.1); print(S34.value)
#
#    #
#    #c1.to_file('ACSALA.cif')
#    #c2.to_file('ACSALA13.cif')
#    #c3.to_file('rep-A.cif')
#    #c3.to_file('rep-B.cif')
#    ##print(S.value)
#    ##print(c)
#    ##print(c.mol_sites[0].molecule.smile)
#    ##rep = c.get_1D_representation()
#    #print(rep1.to_string())
#    #print(rep3.to_string())
#    #print(rep2.to_string())
#    #print(rep4.to_string())
#
#    #print(c1)
#    #print(c3)
#    #print(c1.mol_sites[0].molecule.mol.cart_coords[:3,:])
#    #print(c3.mol_sites[0].molecule.mol.cart_coords[:3,:])
#    #d = np.zeros([len(f1[0]), 2])
#    #print(f1[0])
#    #d[:,0] = f1[0]
#    #d[:,1] = f1[1]
#    #np.savetxt('ref.txt', d)
#
#    def fun(x, f0, smiles):
#        vec = [[14, 0, 11.45,  6.60, 11.39, 95.5], [0.22, 0.41, 0.03, x[0], x[1], x[2], -85.5, -2.4, -1.7, 1]]
#        rep = representation(vec, smiles)
#        c = rep.to_pyxtal()
#        des = descriptor(c)
#        f = des.smear(1.0)
#        f1 = [des.x, f[0]]
#        S = Similarity(f0, f1, l=0.1)
#        #print(x)
#        print(S.value)
#        return -S.value
#
#    # 14 0 11.45  6.60 11.39  95.5 1 0.22 0.41 0.03  -42.7   26.4   32.0  -85.5   -2.4   -1.7 1
#    # 14 0 12.10  6.49 11.32 111.5 1 0.27 0.42 0.92 -135.9   26.1  -32.1   84.0    5.3   -1.5 1 
#    smiles = ['CC(=O)Oc1ccccc1C(O)=O']
#    x0 = [-42.7, 26.4, 32.0]
#    f0 = np.loadtxt('ref.txt')
#    f0 = [f0[:,0], f0[:, 1]]
#    res = minimize(fun, x0, method='Powell', args=(f0, smiles), tol=1e-8, options={'maxiter':1000, 'disp': True}) 
#    #res = minimize(fun, x0, method='CG', args=(f0, smiles), tol=1e-6, options={'maxiter':1000, 'disp': True}) #opt w/o gradient
#    print(res.x)
#    x = res.x
#    vec = [[14, 0, 11.45,  6.60, 11.39, 95.5], [0.22, 0.41, 0.03, x[0], x[1], x[2], -85.5, -2.4, -1.7, 1]]
#    rep = representation(vec, smiles)
#    c = rep.to_pyxtal()
#    print(c)
#    c.to_file('ans.cif')
#
#    #res = minimize(fun, x0, method='Powell', args=(f0), tol=1e-4, options={'maxiter':20}) #opt w/o gradient


