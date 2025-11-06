"""
Module for PXRD indexing and lattice parameter estimation.
"""
import numpy as np

def generate_possible_hkls(h_max, k_max, l_max, level=2):
    """
    Generate reasonable hkl indices within a cutoff for different crystal systems.

    Args:
        h_max: maximum absolute value for h
        k_max: maximum absolute value for k
        l_max: maximum absolute value for l
        level: level of indexing (0 for triclinic; 1 for monoclinic; 2 for orthorhombic or higher)
    """
    if level == 3: # orthorhombic or higher
        base_signs = [(1, 1, 1)]
    elif level == 2:  # hexagonal (110) (1-10)
        base_signs = [(1, 1, 1), (1, -1, 1)]
    elif level == 1: # monoclinic, baxis unique, (101) (10-1)
        base_signs = [(1, 1, 1), (1, 1, -1)]
    else:
        base_signs = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
                      (1, -1, -1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]
    # Create meshgrid for all h, k, l combinations
    h_vals, k_vals, l_vals = np.meshgrid(
        np.arange(h_max + 1),
        np.arange(k_max + 1),
        np.arange(l_max + 1),
        indexing='ij'
    )

    # Flatten to get all combinations
    h_flat = h_vals.flatten()
    k_flat = k_vals.flatten()
    l_flat = l_vals.flatten()

    # Filter out (0,0,0)
    non_zero_mask = (h_flat**2 + k_flat**2 + l_flat**2) > 0
    h_flat = h_flat[non_zero_mask]
    k_flat = k_flat[non_zero_mask]
    l_flat = l_flat[non_zero_mask]

    # Apply all sign combinations vectorized
    all_hkls = []
    for signs in base_signs:
        sh = signs[0] * h_flat
        sk = signs[1] * k_flat
        sl = signs[2] * l_flat
        hkls_with_signs = np.column_stack([sh, sk, sl])
        all_hkls.append(hkls_with_signs)

    return np.vstack(all_hkls)

def get_cell_params(spg, hkls, two_thetas, wave_length=1.54184):
    """
    Calculate cell parameters for a given set of hkls.

    Args:
        spg (int): space group number
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
    if spg >= 195:  # cubic, only need a
        h_sq_sum = np.sum(hkls**2, axis=1)
        cells = d_spacings * np.sqrt(h_sq_sum)
        #for m in range(len(cells)): print(cells[m], hkls[m], d_spacings[m], two_thetas[m])
        mask = cells < 50.0
        cells = cells[mask]
        cells = np.reshape(cells, [len(cells), 1])
        hkls_out = hkls[mask]

    elif 143 <= spg <= 194:  #  hexagonal, need a and c
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

    elif 75 <= spg <= 142:  # tetragonal, need a and c
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

    elif 16 <= spg <= 74:  # orthorhombic, need a, b, c
        # need three hkls to determine a, b, c
        len_solutions = len(hkls) // 3
        ds = (2 * np.sin(thetas) / wave_length)**2
        A = np.zeros([len(hkls), 3])
        A[:, 0] = hkls[:, 0] ** 2
        A[:, 1] = hkls[:, 1] ** 2
        A[:, 2] = hkls[:, 2] ** 2
        B = np.reshape(ds, [len_solutions, 3])
        A = np.reshape(A, [len_solutions, 3, 3])#; print(A.shape, B.shape)
        xs = np.linalg.solve(A, B)#; print(xs); import sys; sys.exit()
        mask1 = np.all(xs[:, :] > 0, axis=1)
        hkls_out = np.reshape(hkls, (len_solutions, 9))
        hkls_out = hkls_out[mask1]
        xs = xs[mask1]
        cells = np.sqrt(1/xs)
        mask2 = np.all(cells[:, :3] < 50.0, axis=1)
        cells = cells[mask2]
        hkls_out = hkls_out[mask2]

    elif 3 <= spg <= 15:  # monoclinic, need a, b, c, beta
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

def get_d_hkl_from_cell(spg, cells, h, k, l):
    """
    Estimate the maximum hkl indices to consider based on the cell parameters and maximum 2theta.

    Args:
        spg (int): space group number
        cells: cell parameters
        h: h index
        k: k index
        l: l index
    """
    if spg >= 195:  # cubic
        d = cells[:, 0] / np.sqrt(h**2 + k**2 + l**2)
    elif spg >= 143:  # hexagonal
        a, c = cells[:, 0], cells[:, 1]
        d = 1 / np.sqrt((4/3) * (h**2 + h*k + k**2) / a**2 + l**2 / c**2)
    elif spg >= 75:  # tetragonal
        a, c = cells[:, 0], cells[:, 1]
        d = 1 / np.sqrt((h**2 + k**2) / a**2 + l**2 / c**2)
    elif spg >= 16:  # orthorhombic
        a, b, c = cells[:, 0], cells[:, 1], cells[:, 2]
        d = 1 / np.sqrt(h**2 / a**2 + k**2 / b**2 + l**2 / c**2)
    elif spg >= 3:  # monoclinic
        a, b, c, beta = cells[:, 0], cells[:, 1], cells[:, 2], np.radians(cells[:, 3])
        sin_beta = np.sin(beta)
        d = 1 / np.sqrt((h**2 / (a**2 * sin_beta**2)) + (k**2 / b**2) + (l**2 / (c**2 * sin_beta**2)) -
                (2 * h * l * np.cos(beta) / (a * c * sin_beta**2)))
    else:
        raise NotImplementedError("triclinic systems are not supported.")
    return d

def calc_two_theta_from_cell(spg, hkls, cells, wave_length=1.54184):
    """
    Calculate expected 2theta values from hkls and cell parameters.

    Args:
        spg (int): space group number
        hkls: hkl indices (np.array)
        cells: cell parameters
        wave_length: X-ray wavelength, default is Cu K-alpha
    """
    h, k, l = hkls[:, 0], hkls[:, 1], hkls[:, 2]
    if spg >= 195:  # cubic
        a = cells[0]
        d = a / np.sqrt(h**2 + k**2 + l**2)#; print('ddddd', d)
    elif spg >= 143:  # hexagonal
        a, c = cells[0], cells[1]
        d = 1 / np.sqrt((4/3) * (h**2 + h*k + k**2) / a**2 + l**2 / c**2)
    elif spg >= 75:  # tetragonal
        a, c = cells[0], cells[1]
        d = 1 / np.sqrt((h**2 + k**2) / a**2 + l**2 / c**2)
    elif spg >= 16:  # orthorhombic
        a, b, c = cells[0], cells[1], cells[2]
        d = 1 / np.sqrt(h**2 / a**2 + k**2 / b**2 + l**2 / c**2)
    elif spg >= 3:  # monoclinic
        a, b, c, beta = cells[0], cells[1], cells[2], np.radians(cells[3])
        sin_beta = np.sin(beta)
        d = 1 / np.sqrt((h**2 / (a**2 * sin_beta**2)) + (k**2 / b**2) + (l**2 / (c**2 * sin_beta**2)) -
                (2 * h * l * np.cos(beta) / (a * c * sin_beta**2)))
    else:
        raise NotImplementedError("triclinic systems are not supported.")
    sin_theta = wave_length / (2 * d)
    # Handle cases where sin_theta > 1
    valid = sin_theta <= 1#; print(d[~valid])
    two_thetas = 2 * np.degrees(np.arcsin(sin_theta[valid]))
    two_thetas = np.round(two_thetas, decimals=3)
    two_thetas, ids = np.unique(two_thetas, return_index=True)
    return two_thetas, hkls[valid][ids]

def get_seeds(spg, hkls, two_thetas):
    """
    Select all possible seed hkls from the provided hkls based on the space group.

    Args:
        spg (int): space group number
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values

    Returns:
        seed_hkls: a list of hkls for seeding
        seed_thetas: corresponding two_thetas
    """
    seed_hkls = []
    seed_thetas = []
    # For cubic, we can select any hkl.
    if spg >= 195:  # cubic
        seed_hkls, seed_thetas = hkls, two_thetas
    # For hexagonal/tetragonal, select two hkls and  avoid two (h,k,0) or (0,0,l) cases
    elif spg >= 75:  # hexagonal or tetragonal
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

    # For tetragonal, two hkls are needed to determine a and c.
    elif spg >= 16:  # tetragonal
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
    elif spg >= 3:  # monoclinic
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


def get_cell_from_multi_hkls(spg, hkls, two_thetas, long_thetas=None, wave_length=1.54184,
                             tolerance=0.1, use_seed=True, min_score=0.999):
    """
    Estimate the cell parameters from multiple (hkl, two_theta) inputs.
    The idea is to use the Bragg's law and the lattice spacing formula to estimate the lattice parameters.
    It is possible to have mislabelled hkls, so we need to run multiple trials and select the best one.

    Args:
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values
        long_thetas: array of  all observed 2theta values
        spg (int): space group number
        wave_length: X-ray wavelength, default is Cu K-alpha
        tolerance: tolerance for matching 2theta values, default is 0.1 degrees
        use_seed: whether to use seed hkls for initial cell estimation
        min_score: threshold score for consideration

    Returns:
        cells: list of solutions
    """
    if long_thetas is None: long_thetas = two_thetas
    #test_hkls_array = np.array(generate_possible_hkls(max_h=max_h))
    if spg > 194 or 15 < spg < 143:
        level = 3  # orthorhombic or higher
    elif 142 < spg < 195:
        level = 2  # hexagonal
    elif 2 < spg < 16:
        level = 1  # monoclinic
    else:
        level = 0  # triclinic

    if use_seed:
        seed_hkls, seed_thetas = get_seeds(spg, hkls, two_thetas)#; print(seed_hkls, seed_thetas)
        cells, hkls = get_cell_params(spg, seed_hkls, seed_thetas, wave_length)#; print(cells)
    else:
        cells, hkls = get_cell_params(spg, hkls, two_thetas, wave_length)#; print(cells)

    #cells = np.array(cells)
    if len(cells) == 0: return []
    # keep cells up to 4 decimal places
    if spg < 16:
        cells[:, -1] = np.round(cells[:, -1], decimals=2)
        cells[:, :3] = np.round(cells[:, :3], decimals=4)
    elif spg > 194:
        cells = np.round(cells, decimals=5)

    _, unique_ids = np.unique(cells, axis=0, return_index=True)
    hkls = hkls[unique_ids]#; print(cells)  # remove duplicates
    cells = cells[unique_ids]

    # get the maximum h from assuming the cell[-1] is (h00)
    d_100s = get_d_hkl_from_cell(spg, cells, 1, 0, 0)
    d_010s = get_d_hkl_from_cell(spg, cells, 0, 1, 0)
    d_001s = get_d_hkl_from_cell(spg, cells, 0, 0, 1)
    theta_100s = 2*np.degrees(np.arcsin(wave_length / (2 * d_100s)))
    theta_010s = 2*np.degrees(np.arcsin(wave_length / (2 * d_010s)))
    theta_001s = 2*np.degrees(np.arcsin(wave_length / (2 * d_001s)))
    h_maxs = np.array(long_thetas[-1] / theta_100s, dtype=int); h_maxs[h_maxs > 100] = 100
    k_maxs = np.array(long_thetas[-1] / theta_010s, dtype=int); k_maxs[k_maxs > 100] = 100
    l_maxs = np.array(long_thetas[-1] / theta_001s, dtype=int); l_maxs[l_maxs > 100] = 100

    solutions = []
    for i, cell in enumerate(cells):
        test_hkls = np.array(generate_possible_hkls(h_max=h_maxs[i],
                                                    k_max=k_maxs[i],
                                                    l_max=l_maxs[i],
                                                    level=level))
        exp_thetas, exp_hkls = calc_two_theta_from_cell(spg, test_hkls, cell, wave_length)
        if len(exp_thetas) == 0: continue

        errors_matrix = np.abs(long_thetas[:, np.newaxis] - exp_thetas[np.newaxis, :])
        within_tolerance = errors_matrix < tolerance
        has_match = np.any(within_tolerance, axis=1)
        best_errors = np.min(errors_matrix, axis=1)
        #print('best errors', exp_thetas); import sys; sys.exit()

        # Filter to only those with valid matches
        valid_matches = has_match & (best_errors < tolerance)
        matched_peaks = []  # (index, hkl, obs_theta, error)
        valid_peak_indices = np.where(valid_matches)[0]

        for peak_idx in valid_peak_indices:
            obs_theta = long_thetas[peak_idx]
            error = best_errors[peak_idx]
            matched_peaks.append((peak_idx, obs_theta, error))
            #print(error)

        # Score this solution
        n_matched = len(matched_peaks)
        coverage = n_matched / len(long_thetas)
        avg_error = np.mean([match[-1] for match in matched_peaks])
        consistency_score = 1.0 / (1.0 + avg_error)  
        score = coverage * consistency_score
        #unmatches = exp_thetas[~within_tolerance.all(axis=0)]
        #mask = (unmatches > long_thetas[0]) & (unmatches < long_thetas[-1])
        #unmatches = exp_hkls[mask]

        if score > min_score:
            solutions.append({
                'cell': cell,
                'n_matched': n_matched,
                'score': score,
                'id': hkls[i],
                #'unmatched_thetas': unmatches,
            })
            #print(cell, len(unmatches), unmatches)

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
        #'pyxtal/database/cifs/JVASP-85365.cif', # P4/mmm, 0.6s
        #'pyxtal/database/cifs/JVASP-62168.cif', # Pnma, 33s
        #'pyxtal/database/cifs/JVASP-98225.cif', # P21/c, 14s
        #'pyxtal/database/cifs/JVASP-50935.cif', # Pm, 10s
        #'pyxtal/database/cifs/JVASP-28565.cif', # Cm, 100s
        #'pyxtal/database/cifs/JVASP-36885.cif', # Cm, 100s
        #'pyxtal/database/cifs/JVASP-42300.cif', # C2, 178s
        'pyxtal/database/cifs/JVASP-47532.cif', # P2/m,
        ]:
        t0 = time()
        xtal.from_seed(cif)
        xrd = xtal.get_XRD(thetas=[0, 120], SCALED_INTENSITY_TOL=0.5)
        cell_ref = np.sort(np.array(xtal.lattice.encode()))
        long_thetas = xrd.pxrd[:15, 0]
        spg = xtal.group.number
        print("\n", cif, xtal.lattice, xtal.group.symbol, xtal.group.number)
        print(xrd.by_hkl(N_max=10))

        # Get the a list of hkl guesses and sort them by d^2
        if spg >= 195:
            min_score, N_add, N_batch = 0.96, 3, 5
        elif spg > 15:
            min_score, N_add, N_batch = 0.999, 5, 20
        else:
            min_score, N_add, N_batch = 0.999, 8, 20

        if spg >= 16:
            guesses = xtal.group.generate_hkl_guesses(2, 2, 5, max_square=29, total_square=40, verbose=True)
        else:
            if spg in [5, 8, 12, 15]:
                guesses = xtal.group.generate_hkl_guesses(3, 3, 5, max_square=29, total_square=40, verbose=True)
            else:
                guesses = xtal.group.generate_hkl_guesses(3, 3, 4, max_square=29, total_square=35, verbose=True)
                #guesses = xtal.group.generate_hkl_guesses(3, 3, 3, max_square=15, total_square=36, verbose=True)

        guesses = np.array(guesses)
        print("Total guesses:", len(guesses))
        sum_squares = np.sum(guesses**2, axis=(1,2))
        sorted_indices = np.argsort(sum_squares)
        guesses = guesses[sorted_indices]
        if len(guesses) > 500000: guesses = guesses[:500000]
        #guesses = np.array([[[2, 0, 0], [1, 1, 0], [0, 1, 1], [0, 0, 2]]])
        #guesses = np.array([[[2, 0, 0], [1, 1, 0], [0, 0, 2], [2, 0, -2]]])
        #guesses = np.array([[[0, 0, -1], [1, 1, 0], [1, 1, -1], [0, 2, -5]]])

        # Check the quality of each (hkl, 2theta) solutions
        cell2 = np.sort(np.array(xtal.lattice.encode()))
        if spg <= 15 and cell2[3] > 90: cell2[3] = 180 - cell2[3]
        cells_all = np.reshape(cell2, (1, len(cell2)))

        # Try each combination of n peaks from the first n+1 peaks
        n_peaks = len(guesses[0])
        N = min(n_peaks + N_add, len(xrd.pxrd))
        available_peaks = xrd.pxrd[:N, 0]

        thetas = []
        for peak_combo in combinations(range(n_peaks + N_add), n_peaks):
            thetas.extend(available_peaks[list(peak_combo)])
        N_thetas = len(thetas) // n_peaks
        thetas = np.array(thetas)
        thetas = np.tile(thetas, N_batch)
        found = False
        d2 = 0
        for i in range(len(guesses)//N_batch + 1):
            if i == len(guesses)//N_batch:
                N_batch = len(guesses) - N_batch * i
                if N_batch == 0:
                    break
                else:
                    thetas = thetas[:N_thetas * n_peaks * N_batch]
            hkls_t = np.tile(guesses[N_batch*i:N_batch*(i+1)], (1, N_thetas, 1))
            hkls_t = np.reshape(hkls_t, (-1, 3))#, order='F')
            solutions = get_cell_from_multi_hkls(spg, hkls_t, thetas, long_thetas, min_score=min_score, use_seed=False)
            if i % 1000 == 0:
                print(f"Processed {N_batch*(i)}/{d2}, found {len(cells_all)-1} cells.")

            for sol in solutions:
                cell1 = np.sort(np.array(sol['cell']))

                # Check if it is a new solution
                diffs = np.sum((cells_all - cell1)**2, axis=1)
                guess = sol['id']
                score = sol['score']
                d2 = np.sum(guess**2)
                if len(cells_all[diffs < 0.1]) == 0:
                    print(f"Guess: {guess}, {d2}/{len(cells_all)-1} -> {cell1}, {score:.6f}")
                    cells_all = np.vstack((cells_all, cell1))

                # Early stopping for getting high-quality solutions
                if diffs[0] < 0.1:
                    print(f"Guess: {guess}, {d2}/{len(cells_all)-1} -> {cell1}, {score:.6f}")
                    print("High score, exiting early.")
                    found = True
                    break
            if found:
                break
        t1 = time()
        data.append((cif, spg, d2, len(cells_all), score, t1-t0))

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
