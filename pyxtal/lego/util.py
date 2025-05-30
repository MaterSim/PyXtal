import numpy as np
from pyxtal.symmetry import Group
from pyxtal import pyxtal

def get_input_from_letters(spg, sites, dim=3):
    """
    Prepare the input from letters

    Args:
        spg (int): 1-230
        sites (list): [['4a'], ['4b']] or [[0], [1]]
        dim (int): 0, 1, 2, 3

    Returns:
        space group instance + wp_list + dof
    """
    g = Group(spg, dim=dim)
    wp_combo = []
    dof = g.get_lattice_dof()
    for site in sites:
        wp_specie = []
        for s in site:
            if type(s) in [int, np.int64]:
                wp = g[s]
            else:
                wp = g.get_wyckoff_position(s)
            wp_specie.append(wp)
            dof += wp.get_dof()
        wp_combo.append(wp_specie)
    return g, wp_combo, dof

def create_trajectory(dumpfile, trjfile, modes=['Init', 'Iter'], dim=3):
    """
    create the ase trajectory from the dump file.
    This is used to analyze the performance of structure optimization
    """

    from ase.io import write

    keyword_start = 'Space Group'  # \n'
    keyword_end = 'END'  # \n'

    with open(dumpfile, 'r') as f:
        lines = f.readlines()
        starts, ends = [], []
        for i in range(len(lines)):
            l = lines[i].lstrip()
            if l.startswith(keyword_start):
                starts.append(i)
            elif l.startswith(keyword_end):
                ends.append(i)
    # print(len(starts), len(ends))
    assert (len(starts) == len(ends))

    atoms = []
    for i, j in zip(starts, ends):
        # parse each run
        xtal_strs = lines[i:j]
        spg = int(xtal_strs[0].split(':')[1])
        l_type = Group(spg).lattice_type
        elements = []
        numIons = []
        wps = []
        for count, tmp in enumerate(xtal_strs[1:]):
            if tmp.startswith('Element'):
                tmp1 = tmp.split(':')
                tmp2 = tmp1[1].split()
                elements.append(tmp2[0])
                numIons.append(int(tmp2[1]))
                wps.append(tmp2[2:])
            else:
                line_struc = count + 1
                break

        # print(spg, elements, numIons, wps)
        g, wps, dof = get_input_from_letters(spg, wps, dim)
        for line_number in range(line_struc, len(xtal_strs)):
            tmp0 = xtal_strs[line_number].split(':')
            tag = tmp0[0]
            if tag in modes:
                tmp = tmp0[1].split()
                sim = float(tmp[0])
                x = [float(t) for t in tmp[1:]]

                # QZ: to fix
                xtal = pyxtal()
                xtal.from_1d_rep(x, wps, numIons, l_type, elements, dim)

                struc = xtal.to_ase()
                struc.info = {'time': line_number-line_struc, 'fun': sim}
                atoms.append(struc)
                # print(xtal.lattice)
        write(filename=trjfile, images=atoms, format='extxyz')


def calculate_dSdx(x, xtal, des_ref, f, eps=1e-4, symmetry=True, verbose=False):
    """
    Compute dSdx via the chain rule of dSdr*drdx.
    This routine will serve as the jacobian for scipy.minimize

    Args:
        x (float): input variable array to describe a crystal structure
        xtal (instance): pyxtal object
        des_ref (array): reference environment
        f (callable): function to compute the environment
        symmetry (bool): whether or not turn on symmetry
        verbose (bool): output more information
    """
    xtal.update_from_1d_rep(x)
    ids = [0]
    weights = []
    for site in xtal.atom_sites:
        ids.append(site.wp.multiplicity + ids[-1])
        weights.append(site.wp.multiplicity)
    ids = ids[:-1]
    weights = np.array(weights, dtype=float)
    atoms = xtal.to_ase()
    dPdr, P = f.compute_dpdr_5d(atoms, ids)
    dPdr = np.einsum('i, ijklm -> ijklm', weights, dPdr)

    # Compute dSdr [N, M] [N, N, M, 3, 27] => [N, 3, 27]
    # print(P.shape, des_ref.shape, dPdr.shape)
    dSdr = np.einsum("ik, ijklm -> jlm", 2*(P - des_ref), dPdr)

    # Get supercell positions
    ref_pos = np.repeat(atoms.positions[:, :, np.newaxis], 27, axis=2)
    cell_shifts = np.dot(VECTORS, atoms.cell)
    ref_pos += cell_shifts.T[np.newaxis, :, :]

    # Compute drdx via numerical func
    drdx = np.zeros([len(atoms), 3, 27, len(x)])

    xtal0 = xtal.copy()
    for i in range(len(x)):
        x0 = x.copy()
        x0[i] += eps
        # Maybe expensive Reduce calls, just need positions
        xtal0.update_from_1d_rep(x0)
        atoms = xtal0.to_ase()

        # Get supercell positions
        pos = np.repeat(atoms.positions[:, :, np.newaxis], 27, axis=2)
        cell_shifts = np.dot(VECTORS, atoms.cell)
        pos += cell_shifts.T[np.newaxis, :, :]
        drdx[:, :, :, i] += (pos - ref_pos)/eps

    # [N, 3, 27] [N, 3, 27, H] => H
    dSdx = np.einsum("ijk, ijkl -> l", dSdr, drdx)
    return dSdx

def calculate_S(x, xtal, des_ref, f, return_rdf=False):
    """
    Optimization function used for structure generation

    Args:
        x (float): input variable array to describe a crystal structure
        xtal (instance): pyxtal object
        des_ref: reference environment
        f (callable): function to compute the enivronment
        return_rdf (bool): whether to return the radial distribution function

    Returns:
        obj (float): objective function value
    """
    if xtal.lattice is None:
        des = np.zeros(des_ref.shape)
        weights = 1
    else:
        xtal.update_from_1d_rep(x)  # ; print(xtal)
        if xtal.lattice is not None:
            ids = [0]
            weights = []
            for site in xtal.atom_sites:
                ids.append(site.wp.multiplicity + ids[-1])
                weights.append(site.wp.multiplicity)
            ids = ids[:-1]
            weights = np.array(weights, dtype=float)
            atoms = xtal.to_ase(resort=False)
            #des = f.calculate(atoms)['x'][ids]
            des = f.compute_p(atoms, ids, return_rdf=return_rdf)
        else:
            des = np.zeros(des_ref.shape)
            weights = 1

    sim = np.sum((des-des_ref)**2, axis=1)  # /des.shape[1]
    #if ref_CN is not None:
    #    coefs = CNs[:, 1] - ref_CN
    #    coefs[coefs < 0] = 0
    #    penalty = coefs * (f.rcut + 0.01 - CNs[:, 0])
    #    if penalty.sum() > 0: print('debug', penalty.sum(), x[:3])
    #    sim += 5 * coefs * penalty
    obj = np.sum(sim * weights) / sum(xtal.numIons)
    #print(des-des_ref); print(sim); print(obj)#; import sys; sys.exit()
    return obj
