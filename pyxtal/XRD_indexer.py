"""
Module for PXRD indexing and lattice parameter estimation.
"""
import numpy as np
from pyxtal.symmetry import get_bravais_lattice, get_lattice_type, generate_possible_hkls

def get_cell_params(bravais, hkls, two_thetas, wave_length=1.54184):
    """
    Calculate cell parameters for a given set of hkls.

    Args:
        bravais (int): Bravais lattice type (0-13)
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values
        wave_length: X-ray wavelength, default is Cu K-alpha
    """
    hkls = np.array(hkls)
    two_thetas = np.array(two_thetas)

    # Convert to d-spacings
    thetas = np.radians(two_thetas / 2)
    d_spacings = wave_length / (2 * np.sin(thetas))

    cells = []
    ltype = get_lattice_type(bravais)
    if ltype == 6:  # cubic, only need a
        h_sq_sum = np.sum(hkls**2, axis=1)
        cells = d_spacings * np.sqrt(h_sq_sum)
        #for m in range(len(cells)): print(cells[m], hkls[m], d_spacings[m], two_thetas[m])
        mask = cells < 50.0
        cells = cells[mask]
        cells = np.reshape(cells, [len(cells), 1])
        hkls_out = hkls[mask]

    elif ltype == 5:  #  hexagonal, need a and c
        # need two hkls to determine a and c
        len_solutions = len(hkls) // 2
        ds = (2 * np.sin(thetas) / wave_length)**2
        A = np.zeros([len(hkls), 2])
        A[:, 0] = 4/3 * (hkls[:, 0] ** 2 + hkls[:, 0] * hkls[:, 1] + hkls[:, 1] ** 2)
        A[:, 1] = hkls[:, 2] ** 2
        B = np.reshape(ds, [len_solutions, 2])
        A = np.reshape(A, [len_solutions, 2, 2])#; print(A.shape, B.shape)
        xs = np.linalg.solve(A, B)#; print(xs); import sys; sys.exit()
        mask1 = np.all(xs[:, :] > 0, axis=1)
        hkls_out = np.reshape(hkls, (len_solutions, 6))
        hkls_out = hkls_out[mask1]
        xs = xs[mask1]
        cells = np.sqrt(1/xs)
        mask2 = cells[:, 0] < 50.0
        cells = cells[mask2]
        hkls_out = hkls_out[mask2]

    elif ltype == 4:  # tetragonal, need a and c
        # need two hkls to determine a and c
        len_solutions = len(hkls) // 2
        ds = (2 * np.sin(thetas) / wave_length)**2
        A = np.zeros([len(hkls), 2])
        A[:, 0] = hkls[:, 0] ** 2 + hkls[:, 1] ** 2
        A[:, 1] = hkls[:, 2] ** 2
        B = np.reshape(ds, [len_solutions, 2])
        A = np.reshape(A, [len_solutions, 2, 2])#; print(A.shape, B.shape)
        xs = np.linalg.solve(A, B)#; print(xs); import sys; sys.exit()
        mask1 = np.all(xs[:, :] > 0, axis=1)
        hkls_out = np.reshape(hkls, (len_solutions, 6))
        hkls_out = hkls_out[mask1]
        xs = xs[mask1]
        cells = np.sqrt(1/xs)
        mask2 = np.all(cells[:, :2] < 50.0, axis=1)
        cells = cells[mask2]
        hkls_out = hkls_out[mask2]

    elif ltype == 3:  # orthorhombic, need a, b, c
        # need three hkls to determine a, b, c
        len_solutions = len(hkls) // 3
        ds = (2 * np.sin(thetas) / wave_length)**2
        A = np.zeros([len(hkls), 3])
        A[:, 0] = hkls[:, 0] ** 2
        A[:, 1] = hkls[:, 1] ** 2
        A[:, 2] = hkls[:, 2] ** 2
        B = np.reshape(ds, [len_solutions, 3])
        A = np.reshape(A, [len_solutions, 3, 3])
        xs = np.linalg.solve(A, B)#; print(xs); import sys; sys.exit()
        mask1 = np.all(xs[:, :] > 0, axis=1)
        hkls_out = np.reshape(hkls, (len_solutions, 9))
        hkls_out = hkls_out[mask1]
        xs = xs[mask1]
        cells = np.sqrt(1/xs)
        mask2 = np.all(cells[:, :3] < 65.0, axis=1)
        cells = cells[mask2]
        hkls_out = hkls_out[mask2]

    elif ltype == 2:  # monoclinic, need a, b, c, beta
        # need four hkls to determine a, b, c, beta
        len_solutions = len(hkls) // 4
        thetas = np.radians(two_thetas/2)
        ds = (2 * np.sin(thetas) / wave_length)**2
        # hkls (4*N, 3) A=> N, 4, 4
        A = np.zeros([len(hkls), 4])
        A[:, 0] = hkls[:, 0] ** 2
        A[:, 1] = hkls[:, 1] ** 2
        A[:, 2] = hkls[:, 2] ** 2
        A[:, 3] = hkls[:, 0] * hkls[:, 2]
        B = np.reshape(ds, [len_solutions, 4])
        A = np.reshape(A, [len_solutions, 4, 4])#; print(A.shape, B.shape)
        xs = np.linalg.solve(A, B)#; print(xs); import sys; sys.exit()
        mask1 = np.all(xs[:, :3] > 0, axis=1)
        hkls_out = np.reshape(hkls, (len_solutions, 12))#;print(hkls.shape, mask1.shape, A.shape)
        hkls_out = hkls_out[mask1]
        xs = xs[mask1]

        cos_betas = -xs[:, 3] / (2 * np.sqrt(xs[:, 0] * xs[:, 2]))
        masks = np.abs(cos_betas) <= 1/np.sqrt(2)
        xs = xs[masks]
        hkls_out = hkls_out[masks]

        cos_betas = cos_betas[masks]
        sin_betas = np.sqrt(1 - cos_betas ** 2)
        cells = np.zeros([len(xs), 4])
        cells[:, 1] = np.sqrt(1/xs[:, 1])
        cells[:, 3] = np.degrees(np.arccos(cos_betas))
        cells[:, 0] = np.sqrt(1/xs[:, 0]) / sin_betas
        cells[:, 2] = np.sqrt(1/xs[:, 2]) / sin_betas

        # force angle to be less than 90
        mask = cells[:, 3] > 90.0
        cells[mask, 3] = 180.0 - cells[mask, 3]

        mask2 = np.all(cells[:, :3] < 50.0, axis=1)
        cells = cells[mask2]
        hkls_out = hkls_out[mask2]
        #print(cells)
    else:
        msg = "Only cubic, tetragonal, hexagonal, and orthorhombic systems are supported."
        raise NotImplementedError(msg)

    return cells, hkls_out

def get_d_hkl_from_cell(bravais, cells, h, k, l):
    """
    Estimate the maximum hkl indices to consider based on the cell parameters and maximum 2theta.

    Args:
        bravais (int): Bravais lattice type (1-15)
        cells: cell parameters
        h: h index
        k: k index
        l: l index
    """
    ltype = get_lattice_type(bravais)
    if ltype == 6:  # cubic
        d = cells[:, 0] / np.sqrt(h**2 + k**2 + l**2)
    elif ltype == 5:  # hexagonal
        a, c = cells[:, 0], cells[:, 1]
        d = 1 / np.sqrt((4/3) * (h**2 + h*k + k**2) / a**2 + l**2 / c**2)
    elif ltype == 4:  # tetragonal
        a, c = cells[:, 0], cells[:, 1]
        d = 1 / np.sqrt((h**2 + k**2) / a**2 + l**2 / c**2)
    elif ltype == 3:  # orthorhombic    
        a, b, c = cells[:, 0], cells[:, 1], cells[:, 2]
        d = 1 / np.sqrt(h**2 / a**2 + k**2 / b**2 + l**2 / c**2)
    elif ltype == 2:  # monoclinic
        a, b, c, beta = cells[:, 0], cells[:, 1], cells[:, 2], np.radians(cells[:, 3])
        sin_beta = np.sin(beta)
        d = 1 / np.sqrt((h**2 / (a**2 * sin_beta**2)) + (k**2 / b**2) + (l**2 / (c**2 * sin_beta**2)) -
                (2 * h * l * np.cos(beta) / (a * c * sin_beta**2)))
    else:
        raise NotImplementedError("triclinic systems are not supported.")
    return d

def calc_two_theta_from_cell(bravais, hkls, cells, wave_length=1.54184):
    """
    Calculate expected 2theta values from hkls and cell parameters.

    Args:
        bravais (int): Bravais lattice type (1-15)
        hkls: hkl indices (np.array)
        cells: cell parameters
        wave_length: X-ray wavelength, default is Cu K-alpha
    """
    h, k, l = hkls[:, 0], hkls[:, 1], hkls[:, 2]
    ltype = get_lattice_type(bravais)
    if ltype == 6:  # cubic
        a = cells[0]
        d = a / np.sqrt(h**2 + k**2 + l**2)#; print('ddddd', d)
    elif ltype == 5:  # hexagonal
        a, c = cells[0], cells[1]
        d = 1 / np.sqrt((4/3) * (h**2 + h*k + k**2) / a**2 + l**2 / c**2)
    elif ltype == 4:  # tetragonal
        a, c = cells[0], cells[1]
        d = 1 / np.sqrt((h**2 + k**2) / a**2 + l**2 / c**2)
    elif ltype == 3:  # orthorhombic
        a, b, c = cells[0], cells[1], cells[2]
        d = 1 / np.sqrt(h**2 / a**2 + k**2 / b**2 + l**2 / c**2)
    elif ltype == 2:  # monoclinic
        a, b, c, beta = cells[0], cells[1], cells[2], np.radians(cells[3])
        sin_beta = np.sin(beta)
        d = 1 / np.sqrt((h**2 / (a**2 * sin_beta**2)) + (k**2 / b**2) + (l**2 / (c**2 * sin_beta**2)) -
                (2 * h * l * np.cos(beta) / (a * c * sin_beta**2)))
    else:
        raise NotImplementedError("triclinic systems are not supported.")

    # Handle cases where sin_theta > 1
    sin_theta = wave_length / (2 * d)
    valid = sin_theta <= 1#; print(d[~valid])
    two_thetas = 2 * np.degrees(np.arcsin(sin_theta[valid]))
    two_thetas = np.round(two_thetas, decimals=3)
    two_thetas, ids = np.unique(two_thetas, return_index=True)
    return two_thetas, hkls[valid][ids]

def get_seeds(bravais, hkls, two_thetas):
    """
    Select all possible seed hkls from the provided hkls based on the space group.

    Args:
        bravais (int): Bravais lattice type (0-13)
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values

    Returns:
        seed_hkls: a list of hkls for seeding
        seed_thetas: corresponding two_thetas
    """
    seed_hkls = []
    seed_thetas = []
    # For cubic, we can select any hkl.
    ltype = get_lattice_type(bravais)
    if ltype == 6:  # cubic
        seed_hkls, seed_thetas = hkls, two_thetas
    # For hexagonal/tetragonal, select two hkls and avoid two (h,k,0) or (0,0,l) cases
    elif ltype in [4, 5]:
        # loop all possible pairs and exclude invalid ones
        for i in range(len(hkls)-1):
            h1, k1, l1 = hkls[i]
            for j in range(i+1, len(hkls)):
                h2, k2, l2 = hkls[j]
                if (l1 == 0 and l2 == 0):
                    continue
                elif (h1 == 0 and k1 == 0) and (h2 == 0 and k2 == 0):
                    continue
                seed_hkls.append((h1, k1, l1))
                seed_hkls.append((h2, k2, l2))
                seed_thetas.append(two_thetas[i])
                seed_thetas.append(two_thetas[j])

    # For orthororhombic, select three hkls and avoid (h,0,0), (0,k,0), (0,0,l)
    elif ltype == 3:
        for i in range(len(hkls)-1):
            h1, k1, l1 = hkls[i]
            for j in range(i+1, len(hkls)):
                h2, k2, l2 = hkls[j]
                for k in range(j+1, len(hkls)):
                    h3, k3, l3 = hkls[k]
                    if (h1 == 0 and h2 == 0 and h3 == 0):
                        continue
                    elif (k1 == 0 and k2 == 0 and k3 == 0):
                        continue
                    elif (l1 == 0 and l2 == 0 and l3 == 0):
                        continue
                    # Note may still have something like (111), (222), (333)
                    seed_hkls.append((h1, k1, l1))
                    seed_hkls.append((h2, k2, l2))
                    seed_hkls.append((h3, k3, l3))
                    seed_thetas.append(two_thetas[i])
                    seed_thetas.append(two_thetas[j])
                    seed_thetas.append(two_thetas[k])
    elif ltype == 2:  # monoclinic
        for i in range(len(hkls)-1):
            h1, k1, l1 = hkls[i]
            for j in range(i+1, len(hkls)):
                h2, k2, l2 = hkls[j]
                for k in range(j+1, len(hkls)):
                    h3, k3, l3 = hkls[k]
                    for m in range(k+1, len(hkls)):
                        h4, k4, l4 = hkls[m]
                        if (h1 == 0 and h2 == 0 and h3 == 0 and h4 == 0):
                            continue
                        elif (k1 == 0 and k2 == 0 and k3 == 0 and k4 == 0):
                            continue
                        elif (l1 == 0 and l2 == 0 and l3 == 0 and l4 == 0):
                            continue
                        seed_hkls.append((h1, k1, l1))
                        seed_hkls.append((h2, k2, l2))
                        seed_hkls.append((h3, k3, l3))
                        seed_hkls.append((h4, k4, l4))
                        seed_thetas.append(two_thetas[i])
                        seed_thetas.append(two_thetas[j])
                        seed_thetas.append(two_thetas[k])
                        seed_thetas.append(two_thetas[m])
    else:
        raise NotImplementedError("mono/tri-clinic systems are not supported.")

    return seed_hkls, seed_thetas


def get_unique_thetas(xrd, bravais):
    ltype = get_lattice_type(bravais)
    if ltype == 6:  # cubic
        unique_thetas = xrd.pxrd[:20,0]
    else:
        # deal with (001, 002, 003) for spg < 195
        int_counts = []
        for jj in range(1, 5):
            ratio = np.sin(np.radians(xrd.pxrd[jj, 0] / 2)) / np.sin(np.radians(xrd.pxrd[0, 0] / 2))
            if ratio > 1.1 and abs(ratio - round(ratio)) < 0.01:
                int_counts.append(jj)
        # in case two are related
        if 1 not in int_counts:
            for jj in range(2, 10):
                ratio = np.sin(np.radians(xrd.pxrd[jj, 0] / 2)) / np.sin(np.radians(xrd.pxrd[1, 0] / 2))
                if ratio > 1.1 and abs(ratio - round(ratio)) < 0.01:
                    int_counts.append(jj)

        # if more than three ratios are close to integers, only keep the first peak
        if len(int_counts) > 0:
            unique_thetas = []
            for jj in range(min(15, len(xrd.pxrd))):
                if jj == 0 or jj not in int_counts:
                    unique_thetas.append(xrd.pxrd[jj, 0])
            unique_thetas = np.array(unique_thetas)
            print('processed', unique_thetas)
        else:
            unique_thetas = xrd.pxrd[:min([15, len(xrd.pxrd)]), 0]

    # Remove two consecutive peaks that are very close
    filtered_thetas = [unique_thetas[0]]
    for jj in range(1, len(unique_thetas)):
        if abs(unique_thetas[jj] - unique_thetas[jj-1]) > 0.02:
            filtered_thetas.append(unique_thetas[jj])
    unique_thetas = np.array(filtered_thetas)
    print("Final thetas:", unique_thetas)
    return unique_thetas


def get_cell_from_multi_hkls(bravais, hkls, two_thetas, long_thetas=None, wave_length=1.54184,
                             tolerance=0.1, use_seed=True, trial_hkls=None):
    """
    Estimate the cell parameters from multiple (hkl, two_theta) inputs.
    The idea is to use the Bragg's law and the lattice spacing formula to estimate the lattice parameters.
    It is possible to have mislabelled hkls, so we need to run multiple trials and select the best one.

    Args:
        bravais (int): Bravais lattice type (0-13)
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values
        long_thetas: array of  all observed 2theta values
        wave_length: X-ray wavelength, default is Cu K-alpha
        tolerance: tolerance for matching 2theta values, default is 0.1 degrees
        use_seed: whether to use seed hkls for initial cell estimation
        trial_hkls: pre-generated trial hkls to speed up calculation

    Returns:
        cells: list of solutions
    """
    if long_thetas is None: long_thetas = two_thetas

    if use_seed:
        seed_hkls, seed_thetas = get_seeds(bravais, hkls, two_thetas)#; print(seed_hkls, seed_thetas)
        cells, hkls = get_cell_params(bravais, seed_hkls, seed_thetas, wave_length)#; print(cells)
    else:
        cells, hkls = get_cell_params(bravais, hkls, two_thetas, wave_length)#; print(cells)

    if trial_hkls is None:
        trial_hkls = generate_possible_hkls(bravais, 100, 100, 100)
    #cells = np.array(cells)
    if len(cells) == 0: return []
    # keep cells up to 4 decimal places
    if bravais <= 3:
        cells[:, -1] = np.round(cells[:, -1], decimals=2)
        cells[:, :3] = np.round(cells[:, :3], decimals=4)
    elif bravais > 12:
        cells = np.round(cells, decimals=5)

    _, unique_ids = np.unique(cells, axis=0, return_index=True)
    hkls = hkls[unique_ids]#; print(cells)  # remove duplicates
    cells = cells[unique_ids]

    # get the maximum h from assuming the cell[-1] is (h00)
    d_100s = get_d_hkl_from_cell(bravais, cells, 1, 0, 0)
    d_010s = get_d_hkl_from_cell(bravais, cells, 0, 1, 0)
    d_001s = get_d_hkl_from_cell(bravais, cells, 0, 0, 1)
    theta_100s = 2*np.degrees(np.arcsin(wave_length / (2 * d_100s)))
    theta_010s = 2*np.degrees(np.arcsin(wave_length / (2 * d_010s)))
    theta_001s = 2*np.degrees(np.arcsin(wave_length / (2 * d_001s)))
    h_maxs = np.array(long_thetas[-1] / theta_100s, dtype=int); h_maxs[h_maxs > 100] = 100
    k_maxs = np.array(long_thetas[-1] / theta_010s, dtype=int); k_maxs[k_maxs > 100] = 100
    l_maxs = np.array(long_thetas[-1] / theta_001s, dtype=int); l_maxs[l_maxs > 100] = 100

    solutions = []
    for i, cell in enumerate(cells):
        if cell.min() <= 2.0: continue  # skip too small cells
        h_max, k_max, l_max = h_maxs[i], k_maxs[i], l_maxs[i]
        mask = (trial_hkls[:,0] <= h_max) & (trial_hkls[:,1] <= k_max) & (trial_hkls[:,2] <= l_max)
        test_hkls = trial_hkls[mask]
        #print(len(test_hkls), h_max, k_max, l_max, trial_hkls.shape)
        exp_thetas, exp_hkls = calc_two_theta_from_cell(bravais, test_hkls, cell, wave_length)
        if len(exp_thetas) == 0: continue

        errors_matrix = np.abs(long_thetas[:, np.newaxis] - exp_thetas[np.newaxis, :])
        within_tolerance = errors_matrix < tolerance
        has_obs_match = np.any(within_tolerance, axis=1)
        ids_matched = np.where(has_obs_match)[0]

        if len(ids_matched) == len(long_thetas):
            # Get the obs. peaks
            matched_peaks = []
            for id in ids_matched:
                errors = errors_matrix[id]
                obs_theta = long_thetas[id]
                hkl_id = np.argmin(errors)
                error = errors[hkl_id]
                matched_peaks.append((exp_hkls[hkl_id], obs_theta, error))

            mis_obs_match = np.any(within_tolerance, axis=0)
            ids_mis_matched = np.where(~mis_obs_match)[0]

            mis_matched_peaks = []
            for id in ids_mis_matched:
                hkl = exp_hkls[id]
                theta = exp_thetas[id]
                if theta < long_thetas[-1] and abs(hkl).max() < 3:
                    mis_matched_peaks.append((hkl, theta))

            if len(mis_matched_peaks) <= 15:
                solutions.append({
                    'cell': cell,
                    'matched_peaks': matched_peaks,
                    'mis_matched_peaks': mis_matched_peaks,
                    'id': hkls[i],
                })
                #if len(mis_matched_peaks) > 0:
                #    print(cell, mis_matched_peaks)#; import sys; sys.exit()
    
    return solutions


if __name__ == "__main__":
    from pyxtal import pyxtal
    from itertools import combinations
    from time import time
    np.set_printoptions(precision=4, suppress=True)

    xtal = pyxtal()
    data = []
    for cif in [
        #'pyxtal/database/cifs/JVASP-97915.cif', # Fm-3m, 0.9s
        #'pyxtal/database/cifs/JVASP-86205.cif', # Im-3 204, 0.1s
        #'pyxtal/database/cifs/JVASP-28634.cif', # P3m1, 0.1s
        'pyxtal/database/cifs/JVASP-85365.cif', # P4/mmm, 0.6s
        'pyxtal/database/cifs/JVASP-62168.cif', # Pnma, 33s
        #'pyxtal/database/cifs/JVASP-98225.cif', # P21/c, 14s
        #'pyxtal/database/cifs/JVASP-50935.cif', # Pm, 10s
        #'pyxtal/database/cifs/JVASP-28565.cif', # Cm, 100s
        #'pyxtal/database/cifs/JVASP-36885.cif', # Cm, 100s
        #'pyxtal/database/cifs/JVASP-42300.cif', # C2, 178s
        #'pyxtal/database/cifs/JVASP-47532.cif', # P2/m,
        #'pyxtal/database/cifs/JVASP-25063.cif', # bad slab
        #'pyxtal/database/cifs/JVASP-119184.cif', # Imm2 44
        #'pyxtal/database/cifs/JVASP-141590.cif', # R-3m 166
        #'pyxtal/database/cifs/JVASP-45907.cif', # R-3m 166
        #'pyxtal/database/cifs/JVASP-119739.cif', # C2/m 12
        ]:
        t0 = time()
        xtal.from_seed(cif)
        xrd = xtal.get_XRD(thetas=[0, 120], SCALED_INTENSITY_TOL=0.5)
        cell_ref = np.sort(np.array(xtal.lattice.encode()))
        long_thetas = xrd.pxrd[:20, 0]
        spg = xtal.group.number
        bravais = get_bravais_lattice(spg)
        print("\n", cif, xtal.lattice, xtal.group.symbol, spg)
        print(xrd.by_hkl(N_max=20))

        # Get the a list of hkl guesses and sort them by d^2
        if bravais > 10:
            N_add, N_batch = 3, 5
        elif bravais > 2:
            N_add, N_batch = 5, 20
        else:
            N_add, N_batch = 8, 20

        if bravais >= 7:
            guesses = xtal.group.generate_hkl_guesses(2, 3, 10, max_square=101, total_square=110, verbose=True)
        elif spg >= 3:
            guesses = xtal.group.generate_hkl_guesses(2, 2, 5, max_square=29, total_square=40, verbose=True)
            #guesses = xtal.group.generate_hkl_guesses(2, 2, 2, max_square=29, total_square=40, verbose=True)
        else:
            if bravais == 2:  # monoclinic-S
                guesses = xtal.group.generate_hkl_guesses(3, 3, 5, max_square=29, total_square=40, verbose=True)
            else:
                guesses = xtal.group.generate_hkl_guesses(4, 3, 4, max_square=29, total_square=35, verbose=True)
                #guesses = xtal.group.generate_hkl_guesses(3, 3, 3, max_square=15, total_square=36, verbose=True)

        guesses = np.array(guesses)
        print("Total guesses:", len(guesses))
        sum_squares = np.sum(guesses**2, axis=(1,2))
        sorted_indices = np.argsort(sum_squares)
        guesses = guesses[sorted_indices]
        if len(guesses) > 500000: guesses = guesses[:500000]
        #guesses = np.array([[[0, 0, 2], [1, 0, 1], [0, 1, 3]]])
        #guesses = np.array([[[2, 0, 0], [1, 1, 0], [0, 0, 2], [2, 0, -2]]])
        #guesses = np.array([[[0, 0, -1], [1, 1, 0], [1, 1, -1], [0, 2, -5]]])

        # Check the quality of each (hkl, 2theta) solutions
        cell2 = np.sort(np.array(xtal.lattice.encode()))
        if spg <= 15 and cell2[3] > 90: cell2[3] = 180 - cell2[3]

        trial_hkls = generate_possible_hkls(bravais, 100, 100, 100)

        # Try each combination of n peaks from the first n+1 peaks
        unique_thetas = get_unique_thetas(xrd, bravais)
        n_peaks = len(guesses[0])
        N = min(n_peaks + N_add, len(unique_thetas))
        available_peaks = unique_thetas[:N]; print(available_peaks)

        thetas = []
        for peak_combo in combinations(range(n_peaks + N_add), n_peaks):
            thetas.extend(available_peaks[list(peak_combo)])
        N_thetas = len(thetas) // n_peaks
        thetas = np.array(thetas)
        thetas = np.tile(thetas, N_batch)
        found = False
        d2 = 0
        cell_all = []
        for i in range(len(guesses)//N_batch + 1):
            if i == len(guesses)//N_batch:
                N_batch = len(guesses) - N_batch * i
                if N_batch == 0:
                    break
                else:
                    thetas = thetas[:N_thetas * n_peaks * N_batch]
            hkls_t = np.tile(guesses[N_batch*i:N_batch*(i+1)], (1, N_thetas, 1))
            hkls_t = np.reshape(hkls_t, (-1, 3))#, order='F')
            solutions = get_cell_from_multi_hkls(bravais, hkls_t, thetas, long_thetas,
                                                 tolerance=0.25,
                                                 use_seed=False,
                                                 trial_hkls=trial_hkls)
            if i % 1000 == 0:
                print(f"Processed {N_batch*i}/{d2}, found {len(cell_all)} cells.")

            for sol in solutions:
                cell1 = np.sort(np.array(sol['cell']))
                N_mis = len(sol['mis_matched_peaks'])
                guess = sol['id']
                d2 = np.sum(guess**2)
                if len(cell_all) == 0:
                    cell_all = np.array([np.sort(np.array(sol['cell']))])
                    print(f"Guess: {guess}, {d2}/{N_mis}/0 -> {cell1}")
                else:
                    diffs = np.sum((cell_all - cell1)**2, axis=1)
                    if len(cell_all[diffs < 0.1]) == 0:
                        print(f"Guess: {guess}, {d2}/{N_mis}/{len(cell_all)} -> {cell1}")
                        cell_all = np.vstack((cell_all, cell1))
                        if len(cell_all) >= 100 or np.max(cell1) > 100: # stop if find very big cells
                            found = True
                            break

            if found:
                break
        t1 = time()
        data.append((cif, bravais, d2, i*N_batch, len(cell_all), t1-t0))

    for d in data:
        print(d)
"""
('pyxtal/database/cifs/JVASP-97915.cif', 225, 11, 1, 0.9944178674744656, 0.8724310398101807)
('pyxtal/database/cifs/JVASP-86205.cif', 204, 4, 1, 0.9999808799880138, 0.10700225830078125)
('pyxtal/database/cifs/JVASP-28634.cif', 156, 2, 1, 0.9999999999999456, 0.12210202217102051)
('pyxtal/database/cifs/JVASP-85365.cif', 123, 16, 11, 0.999999999999978, 0.5885009765625)
('pyxtal/database/cifs/JVASP-62168.cif', 62, 34, 123, 0.9999999999999891, 32.97966909408569)
('pyxtal/database/cifs/JVASP-98225.cif', 14, 14, 3, 0.9999291925203678, 35.536349296569824)
('pyxtal/database/cifs/JVASP-50935.cif', 6, 6, 25, 0.9998350565457528, 11.00560998916626)
('pyxtal/database/cifs/JVASP-28565.cif', 8, 35, 3824, 0.9995855322890919, 490.50669598579407)

('pyxtal/database/cifs/JVASP-97915.cif', 225, 11, 1, 0.9942656034459425, 0.9070873260498047)
('pyxtal/database/cifs/JVASP-86205.cif', 204, 4, 1, 0.9997322519800832, 0.09690594673156738)
('pyxtal/database/cifs/JVASP-28634.cif', 156, 2, 1, 0.9997539680387753, 0.07663106918334961)
('pyxtal/database/cifs/JVASP-85365.cif', 123, 16, 9, 0.9997570762997596, 0.3383021354675293)
('pyxtal/database/cifs/JVASP-62168.cif', 62, 34, 15, 0.9997220764747486, 13.73779296875)
('pyxtal/database/cifs/JVASP-98225.cif', 14, 14, 1, 0.9997723881589121, 21.532269716262817)
('pyxtal/database/cifs/JVASP-50935.cif', 6, 6, 25, 0.9997130006576738, 10.676042079925537)
('pyxtal/database/cifs/JVASP-28565.cif', 8, 35, 2885, 0.9995837967820846, 218.40485000610352)
('pyxtal/database/cifs/JVASP-36885.cif', 6, 5, 25, 0.9997727973911672, 10.984601259231567)
('pyxtal/database/cifs/JVASP-42300.cif', 5, 25, 1, 0.9993642422227232, 84.00993585586548)
"""
