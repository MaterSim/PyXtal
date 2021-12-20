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
    from ase import Atoms
    from ase.neighborlist import NeighborList

    def calculate(element, packing, a, covera, cutoff):
        if packing in ['bct']:
            crystal = bulk(element, packing, a=a, c=a*covera)
        elif packing in ['tetragonal']:
            crystal = Atoms(element, cell=(a, a, a*covera), pbc=True,
                          scaled_positions=[(0, 0, 0)])
        elif packing in ['diamond']:
            crystal = Atoms(8*element, cell=(a, a, a*covera), pbc=True,
                          scaled_positions=[(0, 0, 0), (0.25, 0.25, 0.25),
                                            (0, 0.5, 0.5), (0.25, 0.75, 0.75),
                                            (0.5, 0, 0.5), (0.75, 0.25, 0.75),
                                            (0.5, 0.5, 0), (0.75, 0.75, 0.25)])
        else:
            crystal = bulk(element, packing, a=a, covera=covera)
        center = crystal.get_positions()[0]
        cell = crystal.get_cell()
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

            return q4, q6
        else:
            return None

    _as = np.linspace(4.0, 8.0, 30) #[0, 4.0, 5.0, 6.0, 7.0]
    Npt = 80
    cut = _as[-1] + 1
    data = []
    for v in np.linspace(4, 16, Npt):
        data.append(('Xe', 'hcp', np.sqrt(v/3)))
    hcps = []
    for i, d in enumerate(data):
        (element, packing, covera) = d
        for j, a in enumerate(_as):
            q4, q6 = calculate(element, packing, a, covera, cut)
            hcps.append([q4, q6])

    data = []
    for v in np.linspace(1, 8, Npt):
        data.append(('Xe', 'bct', np.sqrt(v/4)))
    bcts = [] 
    for i, d in enumerate(data):
        (element, packing, covera) = d
        for j, a in enumerate(_as):
            q4, q6 = calculate(element, packing, a, covera, cut)
            bcts.append([q4, q6])

    data = []
    for v in np.linspace(1, 8, Npt):
        data.append(('Xe', 'tetragonal', np.sqrt(v/4)))
    ts = []
    for i, d in enumerate(data):
        (element, packing, covera) = d
        for j, a in enumerate(_as):
            q4, q6 = calculate(element, packing, a, covera, cut)
            ts.append([q4, q6])

    data = []
    for v in np.linspace(1, 8, Npt):
        data.append(('Xe', 'diamond', np.sqrt(v/4)))
 
    ds = []
    for i, d in enumerate(data):
        (element, packing, covera) = d
        for j, a in enumerate(_as+3):
            q4, q6 = calculate(element, packing, a, covera, cut)
            ds.append([q4, q6])


    hcps = np.array(hcps)
    bcts = np.array(bcts)
    ds = np.array(ds)
    ts = np.array(ts)
    import matplotlib.pyplot as plt
    plt.figure(figsize=[12, 12])
    plt.scatter(hcps[:,0], hcps[:,1], label='hcp')
    plt.scatter(bcts[:,0], bcts[:,1], label='bct')
    plt.scatter(ts[:,0], ts[:,1], label='sc')
    plt.scatter(ds[:,0], ds[:,1],  label='diamond')
    #for j, a in enumerate(_as):
    #    plt.plot(hcps[j,:,0], hcps[j,:,1], '--', label='hcp'+str(a))
    #    plt.plot(bcts[j,:,0], bcts[j,:,1], '-d', label='bct'+str(a))
    #    plt.plot(ts[j,:,0], ts[j,:,1], '-*', label='tetra'+str(a))
    #    plt.plot(ds[j,:,0], ds[j,:,1], '-s', label='diamond'+str(a))
    plt.legend()
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xlabel('q4')
    plt.ylabel('q6')
    plt.savefig('0.png')

            

