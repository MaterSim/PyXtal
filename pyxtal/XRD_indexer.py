"""
Module for PXRD indexing and lattice parameter estimation.
"""
from pyxtal import pyxtal
from pyxtal.lattice import Lattice
import numpy as np
from itertools import combinations
from pyxtal.symmetry import rf, Group, get_bravais_lattice, get_lattice_type, generate_possible_hkls
from pyxtal.database.element import Element
from collections import defaultdict

def find_wp_assignments(comp, ids, nums):
    """
    Assigns Wyckoff position IDs to a composition. This function can handle
    cases where a composition number is a sum of multiple Wyckoff multiplicities.

    Args:
        comp (list): The target composition counts, e.g., [18, 6].
        ids (list): The list of Wyckoff position IDs, e.g., [1, 1, 6, 8, 9].
        nums (list): The corresponding multiplicities for each ID, e.g., [8, 8, 4, 2, 2].

    Returns:
        list: A list of all possible valid assignments. Each assignment is a list
              of lists, where each inner list contains the WP IDs for an element.
    """

    # Pair IDs with their multiplicities and indices for unique tracking
    wp_info = list(zip(ids, nums, range(len(ids))))

    # --- Helper function to find all subsets that sum to a target ---
    def find_subsets(target, available_wps, start_index=0, current_subset=[]):
        if target == 0:
            yield current_subset
            return
        if target < 0 or start_index == len(available_wps):
            return

        for i in range(start_index, len(available_wps)):
            wp_id, wp_num, original_index = available_wps[i]

            # To handle duplicate numbers, we only proceed if this is the first
            # occurrence or if the previous identical element was not chosen.
            #if i > start_index and available_wps[i][1] == available_wps[i-1][1]:
            #    continue

            # Recurse with the current element included
            yield from find_subsets(
                target - wp_num,
                available_wps,
                i + 1,
                current_subset + [available_wps[i]]
            )

    # --- Main backtracking solver ---
    solutions = []

    def solve(comp_targets, current_assignment, available_wps):
        if not comp_targets:
            # Sort the assignment by the original index before storing
            sorted_assignment = sorted(current_assignment, key=lambda x: x[0])
            # Extract just the WP lists in the correct order
            final_solution = [part for index, part in sorted_assignment]
            solutions.append(final_solution)
            return

        original_index, target_comp = comp_targets[0]
        remaining_comp = comp_targets[1:]

        # Find all possible ways to form the target composition number
        possible_subsets = list(find_subsets(target_comp, available_wps))

        for subset in possible_subsets:
            # For each valid subset, create a new assignment
            new_assignment_part = [wp[0] for wp in subset] # Get just the IDs

            # Determine the remaining available WPs for the next recursive call
            used_indices = {wp[2] for wp in subset}
            next_available_wps = [wp for wp in available_wps if wp[2] not in used_indices]

            # Recurse, adding the original index along with the assignment part
            solve(remaining_comp, current_assignment + [(original_index, new_assignment_part)], next_available_wps)

    # Pair composition values with their original indices
    indexed_comp = list(enumerate(comp))
    # Sort by composition value (descending) to prune search space faster
    sorted_indexed_comp = sorted(indexed_comp, key=lambda x: x[1], reverse=True)

    solve(sorted_indexed_comp, [], wp_info)

    # --- Remove duplicate solutions ---
    # The backtracking algorithm can produce solutions that are identical in content
    # but were arrived at through different paths (e.g., choosing identical
    # input WPs in a different order). We remove these duplicates here.
    unique_solutions = []
    seen = set()
    for sol in solutions:
        # Convert the solution to a canonical, hashable representation
        # by sorting each inner list and converting the whole structure to tuples.
        canonical_sol = tuple(tuple(sorted(part)) for part in sol)
        if canonical_sol not in seen:
            seen.add(canonical_sol)
            unique_solutions.append(sol)

    return unique_solutions


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
                             tolerance=0.1, use_seed=True, trial_hkls=None,
                             max_mismatch=20):
    """
    Estimate the cell parameters from multiple (hkl, two_theta) inputs.
    The idea is to use the Bragg's law to estimate the lattice parameters.
    We need to run multiple trials and select the best one.

    Args:
        bravais (int): Bravais lattice type (0-13)
        hkls: list of (h, k, l) tuples
        two_thetas: list of 2theta values
        long_thetas: array of  all observed 2theta values
        wave_length: X-ray wavelength, default is Cu K-alpha
        tolerance: tolerance for matching 2theta values, default is 0.1 degrees
        use_seed: whether to use seed hkls for initial cell estimation
        trial_hkls: pre-generated trial hkls to speed up calculation
        max_mismatch: maximum number of mismatched peaks allowed

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
        exp_thetas, exp_hkls = calc_two_theta_from_cell(bravais, test_hkls, cell, wave_length)
        if len(exp_thetas) == 0: continue

        errors_matrix_raw = long_thetas[:, np.newaxis] - exp_thetas[np.newaxis, :]
        errors_matrix = np.abs(errors_matrix_raw)
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
                error = errors_matrix_raw[id, hkl_id]
                matched_peaks.append((exp_hkls[hkl_id], obs_theta, error))

            mis_obs_match = np.any(within_tolerance, axis=0)
            ids_mis_matched = np.where(~mis_obs_match)[0]

            mis_matched_peaks = []
            for id in ids_mis_matched:
                hkl = exp_hkls[id]
                theta = exp_thetas[id]
                if theta < long_thetas[-1] and abs(hkl).max() < 3:
                    mis_matched_peaks.append((hkl, theta))

            if len(mis_matched_peaks) <= max_mismatch:
                solutions.append({
                    'cell': cell,
                    'matched_peaks': matched_peaks,
                    'mis_matched_peaks': mis_matched_peaks,
                    'id': hkls[i],
                })
                #if len(mis_matched_peaks) > 0:
                #    print(cell, mis_matched_peaks)#; import sys; sys.exit()

    return solutions

def get_cell_from_thetas(spg, long_thetas, N_add=5, max_mismatch=20, theta_tol=0.25,
                         cell_tol=0.2, hkl_max=(2, 3, 10), max_square=20,
                         N_batch=20, verbose=False):
    """
    Estimate the cell parameters from multiple 2theta inputs.

    Args:
        spg (int): space group number (1-230)
        long_thetas: list of 2theta values
        N_add: number of extra peaks to consider beyond the number of hkls
        max_mismatch: maximum number of mismatched peaks allowed
        theta_tol: tolerance for matching 2theta values, default is 0.25 degrees
        cell_tol: tolerance for considering two cells as the same, default is 0.2
        N_batch: batch size for processing guesses
        verbose: whether to print verbose output

    Returns:
        cells: list of solutions
    """
    bravais = get_bravais_lattice(spg)
    trial_hkls = generate_possible_hkls(bravais, 50, 50, 50)
    unique_thetas = long_thetas
    group = Group(spg)
    h, k, l = hkl_max
    guesses = group.generate_hkl_guesses(h, k, l, max_square=max_square,
                                         verbose=True)
    guesses = np.array(guesses)
    print("Total guesses:", len(guesses), unique_thetas)
    sum_squares = np.sum(guesses**2, axis=(1,2))
    sorted_indices = np.argsort(sum_squares)
    guesses = guesses[sorted_indices]

    n_peaks = len(guesses[0])
    N = min([n_peaks + N_add, len(unique_thetas)])
    available_peaks = long_thetas[:N]

    thetas = []
    for peak_combo in combinations(range(N), n_peaks):
        thetas.extend(available_peaks[list(peak_combo)])
    N_thetas = len(thetas) // n_peaks
    thetas = np.array(thetas)
    thetas = np.tile(thetas, N_batch)

    results = []
    cell_all = []

    for i in range(len(guesses)//N_batch + 1):
        if i == len(guesses)//N_batch:
            N_batch = len(guesses) - N_batch * i
            if N_batch == 0:
                break
            else:
                thetas = thetas[:N_thetas * n_peaks * N_batch]
        hkls_t = np.tile(guesses[N_batch*i:N_batch*(i+1)], (1, N_thetas, 1))
        hkls_t = np.reshape(hkls_t, (-1, 3))
        sols = get_cell_from_multi_hkls(bravais, hkls_t, thetas, long_thetas,
                                        use_seed=False,
                                        trial_hkls=trial_hkls,
                                        tolerance=theta_tol,
                                        max_mismatch=max_mismatch)
        for sol in sols:
            guess, match, mis_match = sol['id'], sol['matched_peaks'], sol['mis_matched_peaks']
            if len(match) == len(long_thetas):
                cell1 = sol['cell'] #np.sort(np.array(sol['cell']))
                d2 = np.sum(guess**2)

                if len(cell_all) == 0:
                    cell_all = np.array([cell1])
                    if verbose:
                        print(f"Guess: {guess}, {d2}, {cell1}, {len(mis_match)}/{len(long_thetas)}")
                    results.append(sol)
                else:
                    diffs = np.sum((cell_all - cell1)**2, axis=1)
                    if len(cell_all[diffs < cell_tol]) == 0:
                        if verbose:
                            print(f"Guess: {guess}, {d2}, {cell1}, {len(mis_match)}/{len(long_thetas)}")
                        results.append(sol)
                        cell_all = np.vstack([cell_all, cell1])

    return results

class XtalManager:
    def __init__(self, spg, species, numIons, cell, WPs):
        """
        Crystal Manager is used to handle crystal structure related operations.

        Args:
            spg (int): Space group number
            WPs (list): Wyckoff positions
            cell (list): Cell parameters
        """
        self.spg = Group(spg)
        self.WPs = WPs
        self.cell = cell
        self.species = species
        self.numIons = numIons
        dof = 0
        sites = [[] for _ in range(len(species))]
        for i, wp in enumerate(WPs):
            for _wp in wp:
                dof += self.spg[_wp].get_dof()
                sites[i].append(self.spg[_wp].get_label())
        self.dof = dof
        self.sites = sites
        print(f"Space group: {spg}, Wyckoff positions: {sites}, DOF: {dof}")

    def generate_structure(self):
        """
        Generate the crystal structure from the Wyckoff positions and cell parameters.
        """
        xtal = pyxtal()
        xtal.from_random(3, self.spg, self.species, self.numIons,
                         lattice = self.cell, sites=self.sites, force_pass=True)
        return xtal

class WPManager:
    def __init__(self, spg, cell, composition={'Si': 1, 'O': 2}, max_wp=8, max_Z=8, max_dof=10):
        """
        WP Manager is used to infer likely Wyckoff positions from the given space group,
        cell, composition, and density constraint.

        Args:
            spg (int): Space group number
            cell (list): Cell parameters
            composition (dict): Elemental composition
            max_wp (int): Maximum number of Wyckoff positions to consider
            max_Z (int): Maximum Z value to consider for volume estimation
            max_dof (int): Maximum degrees of freedom to consider
        """
        from pandas import read_csv
        df = read_csv(rf("pyxtal", "database/spg_num_wps_raw.csv"))
        self.spg = spg
        self.df = df[df['spg'] == self.spg]
        self.cell = cell
        self.composition = composition
        self.comp = [composition[key] for key in composition.keys()]
        self.group = Group(spg)
        self.orders = self.group.get_orders()
        self.lattice = Lattice.from_1d_representation(cell, self.group.lattice_type)
        volume = self.lattice.volume
        self.max_wp = max_wp
        self.max_dof = max_dof
        vol = [0, 0]
        for el in composition.keys():
            sp = Element(el)
            vol1, vol2 = 0.8*sp.covalent_radius**3, 4*sp.covalent_radius**3 #sp.vdw_radius**3
            vol[0] += composition[el] * vol1 * np.pi * 4 / 3
            vol[1] += composition[el] * vol2 * np.pi * 4 / 3

        self.Zs = (int(np.round(volume / vol[1])), int(np.round(volume / vol[0])))
        if self.Zs[1] > max_Z:
            self.Zs = (self.Zs[0], max_Z)
        print(f"Estimated Z range: {self.Zs} {volume} {vol}")

    def get_wyckoff_positions_bak(self):
        """
        Infer possible Wyckoff position combinations based on the composition and Z range.
        """
        sols = []
        for Z in range(self.Zs[0], self.Zs[1]+1):
            comp = [n * Z for n in self.comp]
            wps, _, ids = self.group.list_wyckoff_combinations(comp, numWp=(0, self.max_wp))
            indices = sorted(range(len(ids)), key=lambda i: sum(len(x) for x in ids[i]))
            ids = [ids[i] for i in indices]
            wps = [wps[i] for i in indices]
            if len(ids) > 0:
                #print(f"Z={Z}: Found {len(ids)} Wyckoff position combinations.")
                wp_lists = []
                for id in ids:
                    # Check if the combination is alternative to existing ones
                    duplicate = False
                    tmp = [len(self.group)-1-item for sublist in id for item in sublist]
                    for order in self.orders:
                        tmp_list = order[tmp].tolist()
                        #print("Checking order:", wps[i], tmp, tmp_list)
                        if tmp_list in wp_lists:
                            duplicate = True
                            #print(wps[i], "is duplicate")
                            break
                    if not duplicate:
                        dof = [self.group[wp].get_dof() for sublist in id for wp in sublist]
                        wp_lists.append(tmp)
                        sols.append((self.spg, comp, self.lattice, id, len(tmp), sum(dof)))
                        #print("Added:", self.spg, comp, id)
            print(f"Z={Z}: Kept {len(sols)} Wyckoff position combinations.")
            # sort sols by DOF and number of WPs
            sols = sorted(sols, key=lambda x: (x[5], x[4]))
        return sols

    def get_wyckoff_positions(self):
        """
        Infer possible Wyckoff position combinations based on the composition and Z range.
        """
        sols = []
        for Z in range(self.Zs[0], self.Zs[1]+1):
            df_z = self.df[self.df['n_atoms'] == Z * sum(self.comp)]
            #print(df_z)
            if len(df_z) == 0: continue
            #print(f"Z={Z}: Found {len(df_z)} Wyckoff position combinations.")
            comp = [n * Z for n in self.comp]
            wp_lists = []
            for _, row in df_z.iterrows():
                ids = [int(x) for x in row['wps'].split('-')]
                nums = [self.group[id].multiplicity for id in ids]
                if len(ids) > self.max_wp: continue
                solutions = find_wp_assignments(comp, ids, nums)
                #print(comp, ids, nums, solutions)#; import sys; sys.exit()
                for sol in solutions:
                    dof = [self.group[w].get_dof() for wp in sol for w in wp]
                    #print([self.group[w].get_label() for wp in sol for w in wp])
                    if sum(dof) > self.max_dof: continue

                    duplicate = False
                    tmp = [len(self.group)-1-item for sublist in sol for item in sublist]
                    for order in self.orders:
                        tmp_list = order[tmp].tolist()
                        #print("Checking order:", wps[i], tmp, tmp_list)
                        if tmp_list in wp_lists:
                            duplicate = True
                            #print(wps[i], "is duplicate")
                            break
                    if not duplicate:
                        wp_lists.append(tmp)
                        sols.append((self.spg, comp, self.lattice, sol, len(sol), sum(dof)))
                        #print("Added:", self.spg, comp, ids)
            print(f"Z={Z}: Kept {len(sols)} Wyckoff position combinations.")
            # sort sols by DOF and number of WPs
            #import sys; sys.exit()
        sols = sorted(sols, key=lambda x: (x[5], x[4]))
        return sols



class CellManager:
    def __init__(self, spg, params, missing):
        # Store raw parameters
        self.raw_params = params
        # Sort dimensions immediately for consistent comparison (e.g. [5, 44] == [44, 5])
        #self.dims = np.sort(np.array(params))
        self.dims = np.array(params)
        self.missing = missing
        # Proxy for size (Area for 2D, Volume for 3D) used for sorting
        self.group = Group(spg)
        lattice = Lattice.from_1d_representation(self.dims, self.group.lattice_type)
        self.size_proxy = lattice.volume

    def is_supercell_of(self, other, tol=0.05):
        """Instance method: Check if 'self' is an integer multiple (supercell) of 'other'."""
        ratios = self.dims / other.dims

        # Check if ratios are close to integers (1, 2, 3...)
        is_integer = np.all(np.abs(ratios - np.round(ratios)) < tol)
        # Check if it is actually larger (at least one dimension is > 1.01x)
        is_larger = np.any(np.round(ratios) > 1.01)

        return is_integer and is_larger

    def is_similar_to(self, other, tol=0.04):
        """Instance method: Check if 'self' is nearly identical to 'other' (duplicate)."""
        diff = np.abs(self.dims - other.dims) / other.dims
        return np.all(diff < tol)

    def __repr__(self):
        return f"Cell: {self.dims}, Missing: {self.missing}"

    @classmethod
    def consolidate(cls, raw_data, merge_tol=0.04, supercell_tol=0.05, max_solutions=5, verbose=False):
        """
        Class method: Takes raw list of [dims, missing], instantiates objects,
        sorts, merges duplicates, removes supercells, and returns the clean list.
        """
        # 1. Instantiate objects
        solutions = [cls(d[0], d[1], d[2]) for d in raw_data]

        # 2. Sort: Primary = Fewest Missing (Quality), Secondary = Smallest Size (Parsimony)
        solutions.sort(key=lambda x: (x.size_proxy, x.missing))

        kept_solutions = []
        indices_to_skip = set()

        print(f"{'Status':<10} | {'Dims (Sorted)':<28} | {'Missing':<8} | {'Note'}")
        print("-" * 75)

        for i in range(len(solutions)):
            if i in indices_to_skip:
                continue

            base = solutions[i]
            kept_solutions.append(base)
            print(f"{'KEEP':<10} | {str(base.dims):<28} | {base.missing:<8} | Best candidate")
            if len(kept_solutions) >= max_solutions:
                print(f"Reached maximum of {max_solutions} solutions, stopping further consolidation.")
                break

            # Sweep through the rest of the list to find duplicates or supercells
            for j in range(i + 1, len(solutions)):
                if j in indices_to_skip:
                    continue

                candidate = solutions[j]

                # CHECK 1: Merge (Duplicate)
                if candidate.is_similar_to(base, tol=merge_tol):
                    if candidate.missing < base.missing:
                        if verbose:
                            print(f"{'  SWAP':<10} | {str(candidate.dims):<28} | {candidate.missing:<8} | Better than base, replacing")
                        # Replace the last added solution (which was the old base)
                        kept_solutions[-1] = candidate
                        # Update base for future comparisons in this inner loop
                        base = candidate
                    else:
                        # Otherwise, just drop the candidate
                        if verbose:
                            print(f"{'  MERGE':<10} | {str(candidate.dims):<28} | {candidate.missing:<8} | Similar to base (dropped)")
                    continue

                # CHECK 2: Drop (Supercell Artifact)
                if candidate.is_supercell_of(base, tol=supercell_tol):
                    if candidate.missing > base.missing + 2:
                        indices_to_skip.add(j)
                        ratio = np.round(candidate.dims / base.dims, 1)
                        if verbose:
                            print(f"{'  DROP':<10} | {str(candidate.dims):<28} | {candidate.missing:<8} | Supercell {ratio}")
        #if verbose:
        print(f"Consolidation: {len(kept_solutions)} unique solutions retained from {len(raw_data)} entries.")
        return kept_solutions

if __name__ == "__main__":
    data = [([9, 3], [2, 2, 2, 2], [3, 3, 3, 3], 1),
            ([4, 4, 4, 2, 4], [10, 10, 11, 12, 13], [4, 4, 4, 4, 2], 4),
           ]
    for d in data:
        comp, ids, nums, n = d
        sols = find_wp_assignments(comp, ids, nums)
        print('input: ', d, '\n', sols)

    if False:
        from pyxtal import pyxtal
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
                                                     tolerance=0.2,
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
