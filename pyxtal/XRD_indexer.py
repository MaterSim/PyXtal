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

    cell_values = []
    if spg >= 195:  # cubic, only need a
        h_sq_sum = np.sum(hkls**2, axis=1)
        a_values = d_spacings * np.sqrt(h_sq_sum)
        cell_values = [[a] for a in a_values if 0 < a < 50]

    elif 143 <= spg <= 194:  #  hexagonal, need a and c
        # need two hkls to determine a and c
        len_solutions = len(hkls) // 2
        for i in range(len_solutions):
            (h1, k1, l1), (h2, k2, l2) = hkls[2*i], hkls[2*i + 1]
            theta1 = np.radians(two_thetas[2*i] / 2)
            theta2 = np.radians(two_thetas[2*i + 1] / 2)
            d1 = wave_length / (2 * np.sin(theta1))
            d2 = wave_length / (2 * np.sin(theta2))
            A = np.array([[4/3 * (h1**2 + h1*k1 + k1**2), l1**2],
                          [4/3 * (h2**2 + h2*k2 + k2**2), l2**2]])
            b = np.array([1/d1**2, 1/d2**2])
            try:
                x = np.linalg.solve(A, b)
                if x[0] <= 0 or x[1] <= 0:
                    continue
                a = np.sqrt(1/x[0])
                c = np.sqrt(1/x[1])
                if a > 50: continue
                cell_values.append((a, c))
            except np.linalg.LinAlgError:
                continue

    elif 75 <= spg <= 142:  # tetragonal, need a and c
        # need two hkls to determine a and c
        len_solutions = len(hkls) // 2
        for i in range(len_solutions):
            (h1, k1, l1), (h2, k2, l2) = hkls[2*i], hkls[2*i + 1]
            theta1 = np.radians(two_thetas[2*i] / 2)
            theta2 = np.radians(two_thetas[2*i + 1] / 2)
            d1 = wave_length / (2 * np.sin(theta1))
            d2 = wave_length / (2 * np.sin(theta2))
            A = np.array([[h1**2 + k1**2, l1**2],
                          [h2**2 + k2**2, l2**2]])
            b = np.array([1/d1**2, 1/d2**2])
            try:
                x = np.linalg.solve(A, b)
                if x[0] <= 0 or x[1] <= 0:
                    continue
                a = np.sqrt(1/x[0])
                c = np.sqrt(1/x[1])
                if a > 50: continue
                cell_values.append((a, c))
            except np.linalg.LinAlgError:
                continue
    elif 16 <= spg <= 74:  # orthorhombic, need a, b, c
        # need three hkls to determine a, b, c
        len_solutions = len(hkls) // 3
        for i in range(len_solutions):
            (h1, k1, l1), (h2, k2, l2), (h3, k3, l3) = hkls[3*i], hkls[3*i + 1], hkls[3*i + 2]
            theta1 = np.radians(two_thetas[3*i] / 2)
            theta2 = np.radians(two_thetas[3*i + 1] / 2)
            theta3 = np.radians(two_thetas[3*i + 2] / 2)
            d1 = wave_length / (2 * np.sin(theta1))
            d2 = wave_length / (2 * np.sin(theta2))
            d3 = wave_length / (2 * np.sin(theta3))
            A = np.array([[h1**2, k1**2, l1**2],
                          [h2**2, k2**2, l2**2],
                          [h3**2, k3**2, l3**2]])
            b = np.array([1/d1**2, 1/d2**2, 1/d3**2])
            try:
                x = np.linalg.solve(A, b)
                if x[0] <= 0 or x[1] <= 0 or x[2] <= 0:
                    continue
                a = np.sqrt(1/x[0])
                b = np.sqrt(1/x[1])
                c = np.sqrt(1/x[2])
                if a > 50 or b > 50 or c > 50: continue
                cell_values.append((a, b, c))
            except np.linalg.LinAlgError:
                continue

    elif 3 <= spg <= 15:  # monoclinic, need a, b, c, beta
        # need four hkls to determine a, b, c, beta
        N = 4
        len_solutions = len(hkls) // N
        for i in range(len_solutions):
            hkls_sub = hkls[N*i:N*i+N]
            thetas_sub = np.radians(two_thetas[N*i:N*i+N]/2)
            d_sub = (1/wave_length / (2 * np.sin(thetas_sub)))**2

            # Non-linear system; use numerical methods
            from scipy.optimize import minimize
            def objective(params):
                a, b, c, beta = params
                sin_beta2 = np.sin(beta)**2
                cos_beta = np.cos(beta)
                h, k, l = hkls_sub[:, 0], hkls_sub[:, 1], hkls_sub[:, 2]
                d_inv_sq = (h**2 / (a**2 * sin_beta2)) + (k**2 / b**2) + \
                                (l**2 / (c**2 * sin_beta2)) - \
                                (2 * h * l * cos_beta / (a * c * sin_beta2))
                return np.sum((d_sub - d_inv_sq)**2)

            initial_guess = [2.0, 2.0, 2.0, np.pi/2]
            bounds = [(1.8, 50), (1.8, 50), (1.8, 50), (np.pi/4, np.pi*3/4)]
            result = minimize(objective, initial_guess, bounds=bounds)
            if result.success and result.fun < 1e-5:
                a, b, c, beta = result.x
                cell_values.append((a, b, c, np.degrees(beta)))
                #print(result.x, result.fun, hkls_sub, thetas_sub)
    else:
        msg = "Only cubic, tetragonal, hexagonal, and orthorhombic systems are supported."
        raise NotImplementedError(msg)

    return cell_values

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
        d = a / np.sqrt(h**2 + k**2 + l**2)
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
    valid = sin_theta <= 1
    thetas = np.zeros_like(sin_theta)
    thetas[valid] = np.arcsin(sin_theta[valid])
    two_thetas = 2 * np.degrees(thetas)
    two_thetas[~valid] = np.nan  # Mark invalid values as NaN
    return two_thetas

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
                             tolerance=0.05, use_seed=True, min_score=0.999):
    """
    Estimate the cell parameters from multiple (hkl, two_theta) inputs.
    The idea is to use the Bragg's law and the lattice spacing formula to estimate the lattice parameters.
    It is possible to have mislabelled hkls, so we need to run multiple trials and select the best one.

    Args:
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values
        long_thetas: list of all observed 2theta values
        spg (int): space group number
        wave_length: X-ray wavelength, default is Cu K-alpha
        tolerance: tolerance for matching 2theta values, default is 0.05 degrees
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
        cells = get_cell_params(spg, seed_hkls, seed_thetas, wave_length)#; print(cells)
    else:
        cells = get_cell_params(spg, hkls, two_thetas, wave_length)#; print(cells)

    cells = np.array(cells)
    if len(cells) == 0: return []
    cells = np.unique(cells, axis=0)#; print(cells)  # remove duplicates

    # get the maximum h from assuming the cell[-1] is (h00)
    d_100s = get_d_hkl_from_cell(spg, cells, 1, 0, 0)
    d_010s = get_d_hkl_from_cell(spg, cells, 0, 1, 0)
    d_001s = get_d_hkl_from_cell(spg, cells, 0, 0, 1)
    theta_100s = np.degrees(np.arcsin(wave_length / (2 * d_100s)))
    theta_010s = np.degrees(np.arcsin(wave_length / (2 * d_010s)))
    theta_001s = np.degrees(np.arcsin(wave_length / (2 * d_001s)))#; print(len(cells))
    h_maxs = np.array(long_thetas[-1] / theta_100s, dtype=int); h_maxs[h_maxs > 100] = 100
    k_maxs = np.array(long_thetas[-1] / theta_010s, dtype=int); k_maxs[k_maxs > 100] = 100
    l_maxs = np.array(long_thetas[-1] / theta_001s, dtype=int); l_maxs[l_maxs > 100] = 100

    solutions = []
    for i, cell in enumerate(cells):
        test_hkls = np.array(generate_possible_hkls(h_max=h_maxs[i],
                                                    k_max=k_maxs[i],
                                                    l_max=l_maxs[i],
                                                    level=level))
        expected_thetas = calc_two_theta_from_cell(spg, test_hkls, cell, wave_length)
        # Filter out None values
        valid_mask = expected_thetas != None
        valid_thetas = expected_thetas[valid_mask]
        # Now try to index all other peaks using this 'a'

        if len(valid_thetas) > 0:
            valid_thetas = np.array(valid_thetas, dtype=float)
            matched_peaks = []  # (index, hkl, obs_theta, error)

            for peak_idx, obs_theta in enumerate(long_thetas):
                best_match = None
                best_error = float('inf')
                errors = np.abs(obs_theta - valid_thetas)
                within_tolerance = errors < tolerance

                if np.any(within_tolerance):
                    min_idx = np.argmin(errors[within_tolerance])
                    valid_indices = np.where(within_tolerance)[0]
                    best_idx = valid_indices[min_idx]
                    best_error = errors[best_idx]
                    #best_match = (peak_idx, tuple(valid_hkls[best_idx]), obs_theta, best_error)
                    best_match = (peak_idx, obs_theta, best_error)

                #print(cell, peak_idx, best_match)
                if best_match is not None: matched_peaks.append(best_match)

        # Score this solution
        n_matched = len(matched_peaks)
        coverage = n_matched / len(long_thetas)
        avg_error = np.mean([match[-1] for match in matched_peaks])
        consistency_score = 1.0 / (1.0 + avg_error)  # lower error = higher score
        score = coverage * consistency_score

        if score > min_score:
            solutions.append({
                'cell': cell,
                'n_matched': n_matched,
                'n_total': len(long_thetas),
                'score': score,
            })

    return solutions


if __name__ == "__main__":
    from pyxtal import pyxtal
    from itertools import combinations
    xtal = pyxtal()
    xtal.from_seed('pyxtal/database/cifs/JVASP-62168.cif')
    xrd = xtal.get_XRD(thetas=[0, 120], SCALED_INTENSITY_TOL=0.5)
    cell_ref = np.sort(np.array(xtal.lattice.encode()))
    long_thetas = xrd.pxrd[:15, 0]
    spg = xtal.group.number
    print("\nTesting prototype:", xtal.lattice, xtal.group.symbol, xtal.group.number)
    print(xrd.by_hkl(N_max=10))

    # Get the a list of hkl guesses and sort them by d^2
    guesses = xtal.group.generate_hkl_guesses(2, 3, 5, max_square=29, total_square=40, verbose=True)
    guesses = np.array(guesses)
    print("Total guesses:", len(guesses))
    sum_squares = np.sum(guesses**2, axis=(1,2))
    sorted_indices = np.argsort(sum_squares)
    guesses = guesses[sorted_indices]

    # Check the quality of each (hkl, 2theta) solutions
    for guess in guesses[:]:
        found = False
        n_peaks = len(guess)

        # Try each combination of n peaks from the first n+1 peaks
        N_add = 5
        available_peaks = xrd.pxrd[:n_peaks + N_add, 0]

        thetas = []
        for peak_combo in combinations(range(n_peaks + N_add), n_peaks):
            thetas.extend(available_peaks[list(peak_combo)])
        hkls_t = np.tile(guess, (int(len(thetas)/len(guess)), 1))

        solutions = get_cell_from_multi_hkls(spg, hkls_t, thetas, long_thetas, use_seed=False)
        for sol in solutions:
            d2 = np.sum(guess**2)
            print("Guess:", guess.flatten(), d2, "->", sol['cell'], thetas[:len(guess)], "Score:", sol['score'])
            cell1 = np.sort(np.array(sol['cell']))
            diff = np.sum((cell1 - cell_ref)**2)
            if diff < 0.1:
                print("High score, exiting early.")
                found = True
                break
        if found:
            break


    #guesses = [
    #           [(1, 1, 1), (2, 2, 0), (3, 1, 1)],
    #           [(0, 0, 2), (1, 0, 0)],#, (1, 0, 1)],
    #           [(0, 0, 2), (1, 0, 1), (1, 1, 0)],
    #           [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
    #           [(2, 0, 0), (1, 0, 1), (2, 1, 0)],
    #           [(1, 0, 0), (1, 1, 1), (1, 0, 2)],
    #           [(1, 0, 1), (1, 1, 1), (3, 1, 1)],
    #           [(2, 2, 0), (3, 1, 1), (2, 2, 2)],
    #           [(1, 1, 1), (2, 0, 0), (2, 2, 0)],
    #           [(1, 1), (2, 0, 0), (2, 2, 0)], # test bad case
    #           [(0, 0, 1), (2, 0, -1), (2, 0, 1), (4, 0, 0), (1, 1, 0)],
    #        ]
    #for prototype in ["diamond", "graphite", "a-cristobalite", "olivine", "beta-Ga2O3"]:
    #    xtal.from_prototype(prototype)
    #    xrd = xtal.get_XRD(thetas=[0, 120], SCALED_INTENSITY_TOL=0.5)
    #    spg = xtal.group.number
    #    print("\nTesting prototype:", prototype, xtal.lattice, xtal.group.symbol, xtal.group.number)
    #    print(xrd.by_hkl(N_max=10))
    #    if True:
    #        if xtal.group.number >= 195:
    #            guesses = xtal.group.generate_hkl_guesses(3, max_square=16)
    #        elif xtal.group.number >= 75:
    #            guesses = xtal.group.generate_hkl_guesses(2, 2, 6, max_square=36)
    #        elif xtal.group.number >= 16:
    #            guesses = xtal.group.generate_hkl_guesses(3, 3, 3, max_square=20)
    #        else:
    #            guesses = xtal.group.generate_hkl_guesses(4, 2, 2, max_square=17, verbose=True)

    #    guesses = np.array(guesses)
    #    print("Total guesses:", len(guesses))
    #    sum_squares = np.sum(guesses**2, axis=(1,2))
    #    sorted_indices = np.argsort(sum_squares)
    #    guesses = guesses[sorted_indices]

    #    for guess in guesses:
    #        hkls = [tuple(hkl) for hkl in guess if len(hkl) == 3]
    #        if len(hkls)>3 and (hkls[0] == (1, 1, -1) and hkls[1] == (0, 0, -1) and hkls[2] == (2, 0, -1)):
    #            print("Trying guess:", hkls)
    #            run = True
    #        else:
    #            run = False
    #        n_peaks = len(guess)
    #        # select n peaks from first n+1 peaks, skipping a peak
    #        exit = False
    #        if n_peaks < len(xrd.pxrd):
    #            available_peaks = xrd.pxrd[:n_peaks+1, 0]
    #            # Try each combination of n peaks from the first n+1 peaks
    #            for peak_combo in combinations(range(n_peaks+1), n_peaks):
    #                theta = available_peaks[list(peak_combo)]
    #                result = get_cell_from_multi_hkls(spg, hkls, theta, xrd.pxrd[:15, 0])
    #                if result is not None and result['score'] > 0.95:
    #                    print("Guess:", hkls, "->", result['cell'], "Score:", result['score'],
    #                          "Matched peaks:", result['n_matched'], "/", result['n_total'])
    #                    if result['score'] > 0.992:
    #                        print("High score, exiting early.")
    #                        exit = True
    #                        break
    #        if exit:
    #            break
