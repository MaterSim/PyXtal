"""
Module for PXRD indexing and lattice parameter estimation.
"""
import numpy as np

def generate_possible_hkls(max_h=5, level=2):
    """
    Generate reasonable hkl indices within a cutoff for different crystal systems.

    Args:
        max_h: maximum absolute value for h, k, l
        level: level of indexing (0 for triclinic; 1 for monoclinic; 2 for orthorhombic or higher)
    """
    if level == 2:
        base_signs = [(1, 1, 1)]
    elif level == 1: # monoclinic, beta is not 90
        base_signs = [(1, 1, 1), (1, -1, 1)]
    else:
        base_signs = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
                      (1, -1, -1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]
    possible_hkls = []
    for h in range(0, max_h + 1):
        for k in range(0, max_h + 1):
            for l in range(0, max_h + 1):
                if h*h + k*k + l*l > 0:  # exclude (0,0,0)
                    # Add all permutations and sign variations
                    base_hkls = set()
                    for signs in base_signs:
                        sh, sk, sl = signs[0]*h, signs[1]*k, signs[2]*l
                        base_hkls.add((sh, sk, sl))
                    possible_hkls.extend(list(base_hkls))
    return list(set(possible_hkls))  # remove duplicates

def get_cell_params(spg, hkls, two_thetas, wave_length=1.54184):
    """
    Calculate cell parameters for a given set of hkls.

    Args:
        spg (int): space group number
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values
        wave_length: X-ray wavelength, default is Cu K-alpha
    """
    cell_values = []
    if spg >= 195:  # cubic, only need a
        for hkl, two_theta in zip(hkls, two_thetas):
            theta = np.radians(two_theta / 2)
            d = wave_length / (2 * np.sin(theta))
            cell_values.append([d * np.sqrt(hkl[0]**2 + hkl[1]**2 + hkl[2]**2)])
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
                a = np.sqrt(1/x[0])
                b = np.sqrt(1/x[1])
                c = np.sqrt(1/x[2])
                cell_values.append((a, b, c))
            except np.linalg.LinAlgError:
                continue
    elif 3 <= spg <= 15:  # monoclinic, need a, b, c, beta
        # need four hkls to determine a, b, c, beta
        len_solutions = len(hkls) // 4
        for i in range(len_solutions):
            (h1, k1, l1), (h2, k2, l2), (h3, k3, l3), (h4, k4, l4) = hkls[4*i], hkls[4*i + 1], hkls[4*i + 2], hkls[4*i + 3]
            theta1 = np.radians(two_thetas[4*i] / 2)
            theta2 = np.radians(two_thetas[4*i + 1] / 2)
            theta3 = np.radians(two_thetas[4*i + 2] / 2)
            theta4 = np.radians(two_thetas[4*i + 3] / 2)
            d1 = wave_length / (2 * np.sin(theta1))
            d2 = wave_length / (2 * np.sin(theta2))
            d3 = wave_length / (2 * np.sin(theta3))
            d4 = wave_length / (2 * np.sin(theta4))
            # Non-linear system; use numerical methods
            from scipy.optimize import minimize

            def objective(params):
                a, b, c, beta_deg = params
                beta = np.radians(beta_deg)
                eqs = []
                for (h, k, l), d_obs in zip([(h1, k1, l1), (h2, k2, l2), (h3, k3, l3), (h4, k4, l4)],
                                            [d1, d2, d3, d4]):
                    sin_beta_sq = np.sin(beta)**2
                    d_calc_inv_sq = (h**2 / (a**2 * sin_beta_sq)) + (k**2 / b**2) + \
                                    (l**2 / (c**2 * sin_beta_sq)) - \
                                    (2 * h * l * np.cos(beta) / (a * c * sin_beta_sq))
                    eqs.append((1/d_obs**2 - d_calc_inv_sq)**2)
                return sum(eqs)
            initial_guess = [2.0, 2.0, 2.0, 90.0]
            result = minimize(objective, initial_guess)#; print(result.x, result.fun)
            if result.success:
                a, b, c, beta = result.x
                cell_values.append((a, b, c, beta))
    else:
        msg = "Only cubic, tetragonal, hexagonal, and orthorhombic systems are supported."
        raise NotImplementedError(msg)

    return cell_values

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
                             tolerance=0.05, min_matched_peaks=2):
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
        min_matched_peaks: minimum number of matched peaks to consider a valid solution

    Returns:
        cells: estimated lattice parameter
    """
    if long_thetas is None:
        long_thetas = two_thetas
    all_possible_hkls = generate_possible_hkls(max_h=6)
    test_hkls_array = np.array(all_possible_hkls)

    best_solution = None
    best_score = 0

    seed_hkls, seed_thetas = get_seeds(spg, hkls, two_thetas)#; print(seed_hkls, seed_thetas)
    cells = get_cell_params(spg, seed_hkls, seed_thetas, wave_length)#; print(cells)
    cells = np.array(cells)
    cells = np.unique(cells, axis=0)#; print(cells)  # remove duplicates

    for cell in cells:
        # Now try to index all other peaks using this 'a'
        matched_peaks = []  # (index, hkl, obs_theta, error)
        for peak_idx, obs_theta in enumerate(long_thetas):
            best_match = None
            best_error = float('inf')

            # Try all possible hkls for this peak - vectorized version
            expected_thetas = calc_two_theta_from_cell(spg, test_hkls_array, cell, wave_length)

            # Filter out None values
            valid_mask = expected_thetas != None
            valid_thetas = expected_thetas[valid_mask]
            valid_hkls = test_hkls_array[valid_mask]

            if len(valid_thetas) > 0:
                valid_thetas = np.array(valid_thetas, dtype=float)
                errors = np.abs(obs_theta - valid_thetas)
                within_tolerance = errors < tolerance

                if np.any(within_tolerance):
                    min_idx = np.argmin(errors[within_tolerance])
                    valid_indices = np.where(within_tolerance)[0]
                    best_idx = valid_indices[min_idx]
                    best_error = errors[best_idx]
                    best_match = (peak_idx, tuple(valid_hkls[best_idx]), obs_theta, best_error)

            #print(cell, peak_idx, best_match)
            if best_match is not None:
                matched_peaks.append(best_match)

            # Score this solution
            n_matched = len(matched_peaks)
            if n_matched >= min_matched_peaks:
                coverage = n_matched / len(long_thetas)
                avg_error = np.mean([match[3] for match in matched_peaks])
                consistency_score = 1.0 / (1.0 + avg_error)  # lower error = higher score
                total_score = coverage * consistency_score

                if total_score > best_score:
                    best_score = total_score
                    best_solution = {
                        'cell': cell,
                        'n_matched': n_matched,
                        'n_total': len(long_thetas),
                        'score': total_score,
                        #'avg_error': avg_error,
                    }

    return best_solution


if __name__ == "__main__":
    from pyxtal import pyxtal

    xtal = pyxtal()
    guesses = [
               [(1, 1, 1), (2, 2, 0), (3, 1, 1)],
               [(0, 0, 2), (1, 0, 0), (1, 0, 1)],
               [(0, 0, 2), (1, 0, 1), (1, 1, 0)],
               [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
               [(2, 0, 0), (1, 0, 1), (2, 1, 0)],
               [(1, 0, 0), (1, 1, 1), (1, 0, 2)],
               [(1, 0, 1), (1, 1, 1), (3, 1, 1)],
               [(2, 2, 0), (3, 1, 1), (2, 2, 2)],
               [(1, 1, 1), (2, 0, 0), (2, 2, 0)],
               [(1, 1), (2, 0, 0), (2, 2, 0)], # test bad case
               [(0, 0, 1), (2, 0, -1), (2, 0, 1), (4, 0, 0), (1, 1, 0)],
            ]
    for prototype in ["diamond", "graphite", "a-cristobalite", "olivine", "beta-Ga2O3"]:
        xtal.from_prototype(prototype)
        xrd = xtal.get_XRD(thetas=[0, 120], SCALED_INTENSITY_TOL=0.5)
        print("\nTesting prototype:", prototype, xtal.lattice)
        print(xrd.by_hkl(N_max=5))
        spg = xtal.group.number
        for guess in guesses:
            n_peaks = len(guess)
            if spg < 16 and n_peaks < 4: continue
            theta = xrd.pxrd[:n_peaks, 0]
            hkls = [tuple(hkl) for hkl in guess if len(hkl) == 3]
            result = get_cell_from_multi_hkls(spg, hkls, theta, xrd.pxrd[:10, 0])
            if result is not None and result['score'] > 0.9:
                print("Guess:", guess, "->", result['cell'],)
