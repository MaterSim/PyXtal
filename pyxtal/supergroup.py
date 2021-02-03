"""
Module to search for the supergroup symmetry
"""

from copy import deepcopy
from random import sample
import itertools

import numpy as np
from scipy.optimize import minimize

from pymatgen.core.operations import SymmOp
import pymatgen.analysis.structure_matcher as sm

import pyxtal.symmetry as sym
from pyxtal.lattice import Lattice
from pyxtal.wyckoff_site import atom_site
from pyxtal.operations import apply_ops, get_inverse
from pyxtal.wyckoff_split import wyckoff_split

def new_solution(A, refs):
    """
    check if A is already in the reference solutions
    """
    for B in refs:
        match = True
        for a, b in zip(A, B):
            a.sort()
            b.sort()
            if a != b:
                match = False
                break
        if match:
            return False
    return True

def new_structure(struc, refs):
    """
    check if struc is already in the reference solutions
    """
    g1 = struc.group.number
    pmg1 = struc.to_pymatgen()
    for ref in refs:
        g2 = ref.group.number
        if g1 == g2:
            pmg2 = ref.to_pymatgen()
            if sm.StructureMatcher().fit(pmg1, pmg2):
                return False
    return True

def get_best_match(positions, ref, cell):
    """
    find the best match with the reference from a set of positions

    Args:
        positions: N*3 array
        ref: 1*3 array

    Returns:
        position: matched position
        id: matched id
    """

    diffs = positions - ref
    diffs -= np.round(diffs)
    diffs = np.dot(diffs, cell)
    dists = np.linalg.norm(diffs, axis=1)
    id = np.argmin(dists)
    return positions[id], dists[id]

def check_freedom(G, solutions):
    """
    check if the solutions are valid
    a special WP such as (0,0,0) cannot be occupied twice
    """
    valid_solutions = []
    G = sym.Group(G)
    for solution in solutions:
        sites = []
        for s in solution:
            sites.extend(s)
        if G.is_valid_combination(sites):
            valid_solutions.append(solution)
    return valid_solutions

def check_lattice(G, trans, struc, tol=1.0, a_tol=10):
    """
    check if the lattice mismatch is big
    used to save some computational cost
    """
    matrix = np.dot(trans, struc.lattice.get_matrix())
    l1 = Lattice.from_matrix(matrix)
    l2 = Lattice.from_matrix(matrix, ltype=sym.Group(G).lattice_type)
    (a1,b1,c1,alpha1,beta1,gamma1)=l1.get_para(degree=True)
    (a2,b2,c2,alpha2,beta2,gamma2)=l2.get_para(degree=True)
    abc_diff = np.abs(np.array([a2-a1, b2-b1, c2-c1])).max()
    ang_diff = np.abs(np.array([alpha2-alpha1, beta2-beta1, gamma2-gamma1])).max()
    if abc_diff > tol or ang_diff > a_tol:
        return False
    else:
        return True

def check_compatibility(G, relation, sites, elements):
    """
    find the compatible splitter to let the atoms of subgroup H fit the group G.

    Args:
        G: the target space group with high symmetry
        relation: a dictionary to describe the relation between G and H
    """
    G = sym.Group(G)

    #results = {}
    wyc_list = [(str(x.multiplicity)+x.letter) for x in G]
    wyc_list.reverse()

    good_splittings_list=[]

    # A lot of integer math below.
    # The goal is to find all the integer combinations of supergroup
    # wyckoff positions with the same number of atoms

    # each element is solved one at a time
    for i in range(len(elements)):

        site = np.unique(sites[i])
        site_counts = [sites[i].count(x) for x in site]
        possible_wyc_indices = []

        # the sum of all positions should be fixed.
        total_units = 0
        for j, x in enumerate(site):
            total_units += int(x[:-1])*site_counts[j]


        # collect all possible supergroup transitions
        # make sure all sites are included in the split
        for j, split in enumerate(relation):
            # print(j, split)
            if np.all([x in site for x in split]):
                possible_wyc_indices.append(j)
        # for the case of 173 ['2b'] -> 176
        # print(possible_wyc_indices) [2, 3, 5]

        # a vector to represent the possible combinations of positions
        # when the site is [6c, 2b]
        # the split from [6c, 6c] to [12i] will be counted as [2,0].
        # a entire split from [6c, 6c, 6c, 2b] will need [3, 1]

        possible_wycs = [wyc_list[x] for x in possible_wyc_indices]
        blocks = [np.array([relation[j].count(s) for s in site]) for j in possible_wyc_indices]
        block_units = [sum([int(x[:-1])*block[j] for j,x in enumerate(site)]) for block in blocks]

        # print(possible_wycs)  # ['2c', '2d', '4f']
        # print(blocks) # [array([1]), array([1]), array([2])]
        # print(block_units) # [2, 2, 4]

        # the position_block_units stores the total sum multiplicty
        # from the available G's wyckoff positions.
        # below is a brute force search for the valid combinations

        combo_storage = [np.zeros(len(block_units))]
        good_list = []
        # print(combo_storage)
        # print(block_units)
        # print(blocks)
        # print(possible_wycs)
        # print(total_units)
        # print(site_counts)
        while len(combo_storage)!=0:
            holder = []
            for j, x in enumerate(combo_storage):
                for k in range(len(block_units)):
                    trial = np.array(deepcopy(x)) # trial solution
                    trial[k] += 1
                    if trial.tolist() in holder:
                        continue
                    sum_units = np.dot(trial, block_units)
                    if sum_units > total_units:
                        continue
                    elif sum_units < total_units:
                        holder.append(trial.tolist())
                    else:
                        tester = np.zeros(len(site_counts))
                        for l, z in enumerate(trial):
                            tester += z*blocks[l]
                        if np.all(tester == site_counts):
                            G_sites = []
                            for l, number in enumerate(trial):
                                if number==0:
                                    continue
                                elif number==1:
                                    G_sites.append(possible_wycs[l])
                                else:
                                    for i in range(int(number)):
                                        G_sites.append(possible_wycs[l])
                            if G_sites not in good_list:
                                good_list.append(G_sites)
            combo_storage=holder

        if len(good_list)==0:
            # print("cannot find the valid split, quit the search asap")
            return None
        else:
            good_splittings_list.append(good_list)
    # if len(good_splittings_list[0])==1:
    #     print(good_splittings_list[0])
    return good_splittings_list

class supergroups():
    """
    Class to search for the feasible transition to a given super group

    Args:
        struc: pyxtal structure
        G (int): the desired super group number
    """

    def __init__(self, struc, G=None, path=None, d_tol=1.0):
        self.struc0 = struc
        if path is None:
            paths = self.search_paths(struc.group.number, G)
        else:
            paths = [path]

        self.strucs = None
        for p in paths:
            strucs = self.struc_along_path(p, d_tol)
            if len(strucs) == len(p) + 1:
                self.strucs = strucs
                self.path = [self.struc0.group.number] + p
                break

    def __str__(self):
        s = "\n===Transition to super group: "
        if self.strucs is None:
            s += "Unsuccessful, check your input"
        else:
            s += "{:d}".format(self.path[0])
            for i, p in enumerate(self.path[1:]):
                s += " -> {:d}[{:5.3f}]".format(p, self.strucs[i+1].disp)
            s += '\n'
            for struc in self.strucs:
                s += str(struc)
        return s

    def __repr__(self):
        return str(self)

    def search_paths(self, H, G):
        # enumerate all connections between H and G
        # paths = []
        """
        >>> from pyxtal.symmetry import Group
        >>> g = Group(227)
        >>> g.get_max_subgroup_numbers()
        [141, 141, 141, 166, 166, 166, 166, 203, 210, 216]
        """
        #G1 = sym.Group(G).get_max_subgroup_numbers()

        raise NotImplementedError

    def struc_along_path(self, path, d_tol=1.0):
        """
        search for the super group structure along a given path
        """
        strucs = []
        G_strucs = [self.struc0]
        for G in path:
            H = G_strucs[0].group.number
            if sym.get_point_group(G) == sym.get_point_group(H):
                group_type = 'k'
            else:
                group_type = 't'
            for G_struc in G_strucs:
                my = supergroup(G_struc, [G], group_type)
                solutions = my.search_supergroup(d_tol, max_per_G=500)
                new_G_strucs = my.make_supergroup(solutions, show_detail=False)
                if len(new_G_strucs) > 0:
                    strucs.append(G_struc)
                    G_strucs = new_G_strucs
                    break
            if len(new_G_strucs) == 0:
                print('cannot proceed from Spg {:d} to {:d}'.format(H, G))
                break
        # add the final struc
        if len(strucs) == len(path):
            strucs.append(new_G_strucs[0])
        return strucs

class supergroup():
    """
    Class to find the structure with supergroup symmetry

    Args:
        struc: pyxtal structure
        G: list of possible supergroup numbers, default is None
        group_type: `t` or `k`
    """
    def __init__(self, struc, G=None, group_type='t'):

        # initilize the necesary parameters
        self.group_type = group_type
        self.error = False

        # extract the supergroup information
        wyc_supergroups = struc.group.get_min_supergroup(group_type)
        if G is not None:
            self.wyc_supergroups = {}
            ids = [id for id, group in enumerate(wyc_supergroups['supergroup']) if group in G]
            if len(ids) == 0:
                self.error = True
            else:
                for key in wyc_supergroups:
                    self.wyc_supergroups[key] = [wyc_supergroups[key][id] for id in ids]
        else:
            self.wyc_supergroups = wyc_supergroups

        # list of all alternative wycsets
        strucs = struc.get_alternatives()
        for struc in strucs:
            # group the elements, sites, positions
            elements = []
            sites = []
            for at_site in struc.atom_sites:
                e = at_site.specie
                site = str(at_site.wp.multiplicity) + at_site.wp.letter
                if e not in elements:
                    elements.append(e)
                    sites.append([site])
                else:
                    id = elements.index(e)
                    sites[id].append(site)

            # search for the compatible solutions
            solutions = []
            for idx in range(len(self.wyc_supergroups['supergroup'])):
                G = self.wyc_supergroups['supergroup'][idx]
                relation = self.wyc_supergroups['relations'][idx]
                id = self.wyc_supergroups['idx'][idx]
                trans = np.linalg.inv(self.wyc_supergroups['transformation'][idx][:,:3])

                if check_lattice(G, trans, struc):
                    #print(G, relation)
                    results = check_compatibility(G, relation, sites, elements)
                    if results is not None:
                        sols = list(itertools.product(*results))
                        trials = check_freedom(G, sols)
                        sol = {'group': G, 'id': id, 'splits': trials}
                        solutions.append(sol)

            if len(solutions) > 0:
                # exit if one solution is found
                break

        if len(solutions) == 0:
            self.solutions = []
            self.error = True
            print("No compatible solution exists")
        else:
            #print(struc)
            self.sites = sites
            self.elements = elements
            self.struc = struc
            self.solutions = solutions
            self.cell = struc.lattice.matrix

    def search_supergroup(self, d_tol=1.0, max_per_G=2500):
        """
        search for valid supergroup transition

        Args:
            d_tol (float): tolerance for atomic displacement
            max_per_G (int): maximum number of possible solution for each G
        Returns:
            valid_solutions: dictionary
        """
        #self.d_tol = d_tol
        valid_solutions = []
        if len(self.solutions) > 0:
            # extract the valid
            for sols in self.solutions:
                G, id, sols = sols['group'], sols['id'], sols['splits']
                if len(sols) > max_per_G:
                    #print(len(sols))
                    sols = sample(sols, max_per_G)
                for sol in sols:
                    mae, disp, mapping, sp = self.get_displacement(G, id, sol, d_tol*1.1)
                    #print(G, sol, mae, disp)
                    if mae < d_tol:
                        valid_solutions.append((sp, mapping, disp, mae))
        return valid_solutions

    def make_supergroup(self, solutions, once=False, show_detail=True):
        """
        create supergroup structures based on the list of solutions

        Args:
            solutions: list of tuples (splitter, mapping, disp)
            show_detail (bool): print out the detail
        Returns:
            list of pyxtal structures
        """
        G_strucs = []

        if len(solutions) > 0:
            if once:
                disps = np.array([sol[-1] for sol in solutions])
                ID = np.argmin(disps)
                solutions = [solutions[ID]]
            for solution in solutions:
                (sp, mapping, disp, mae) = solution
                #print(mapping, disp, mae)
                G = sp.G.number
                lat1 = np.dot(sp.inv_R[:3,:3].T, self.struc.lattice.matrix)
                lattice = Lattice.from_matrix(lat1, ltype=sp.G.lattice_type)

                details = self.symmetrize(sp, mapping, disp)
                coords_G1, coords_G2, coords_H1, elements = details
                G_struc = self.struc.copy()
                G_struc.group = sp.G

                G_sites = []
                for i, wp in enumerate(sp.wp1_lists):
                    pos = coords_G1[i]
                    pos -= np.floor(pos)
                    pos1 = sym.search_matched_position(sp.G, wp, pos)
                    if pos1 is not None:
                        site = atom_site(wp, pos1, sp.elements[i])
                        G_sites.append(site)
                    else:
                        print(wp)
                        raise RuntimeError("cannot assign the right wp")

                G_struc.atom_sites = G_sites
                G_struc.source = 'supergroup {:6.3f}'.format(mae)
                G_struc.lattice = lattice
                G_struc.numIons *= round(np.abs(np.linalg.det(sp.R[:3,:3])))
                G_struc._get_formula()
                G_struc.disp = mae

                if new_structure(G_struc, G_strucs):
                    G_strucs.append(G_struc)
                    if show_detail:
                        details = self.symmetrize(sp, mapping, disp)
                        _, coords_G2, coords_H1, elements = details
                        self.print_detail(G, coords_H1, coords_G2, elements, disp)
                        print(G_struc)

        return G_strucs

    def get_displacement(self, G, split_id, solution, d_tol):
        """
        For a given solution, search for the possbile supergroup structure

        Args:
            G (int): supergroup number
            split_id (int): integer
            solution (list): e.g., [['2d'], ['6h'], ['2c', '6h', '12i']]
            d_tol (float): tolerance

        Returns:
            mae: mean absolute atomic displcement
            disp: overall cell translation
        """
        sites_G = []
        elements = []
        muls = []
        for i, e in enumerate(self.elements):
            sites_G.extend(solution[i])
            elements.extend([e]*len(solution[i]))
            muls.extend([int(sol[:-1]) for sol in solution[i]])

        # resort the sites_G by multiplicity
        ids = np.argsort(np.array(muls))
        elements = [elements[id] for id in ids]
        sites_G = [sites_G[id] for id in ids]
        #print(G, self.struc.group.number, sites_G)
        splitter = wyckoff_split(G, split_id, sites_G, self.group_type, elements)
        mappings = self.find_mapping(splitter)
        dists = []
        disps = []
        if len(mappings) > 0:
            for mapping in mappings:
                dist, disp, mask = self.symmetrize_dist(splitter, mapping, None, None, d_tol)
                #print("mask---disp", disp, mask, dist)
                dists.append(dist)
                disps.append(disp)
            dists = np.array(dists)
            mae = np.min(dists)
            id = np.argmin(dists)
            disp = disps[id]
            if 0.2 < mae < d_tol:
                # optimize further
                def fun(disp, mapping, splitter, mask):
                    return self.symmetrize_dist(splitter, mapping, disp, mask)[0]
                res = minimize(fun, disps[id], args=(mappings[id], splitter, mask),
                        method='Nelder-Mead', options={'maxiter': 20})
                if res.fun < mae:
                    mae = res.fun
                    disp = res.x
            return mae, disp, mappings[id], splitter
        else:
            return 1000, None, None, None

    def find_mapping(self, splitter, max_num=50):
        """
        search for all mappings for a given splitter

        Args:
            splitter: wyc_splitter object
            max_num (int): maximum number of atomic mapping

        Returns:
            unique solutions
        """
        solution_template = [[None]*len(wp2) for wp2 in splitter.wp2_lists]
        atom_sites_H = self.struc.atom_sites
        assigned_ids = []
        # look for unique assignment from sites_H to sites_G
        for i, wp2 in enumerate(splitter.wp2_lists):
            # choose the sites belong to the same element
            ele = splitter.elements[i]
            e_ids = [id for id, site in enumerate(atom_sites_H) if site.specie==ele]

            if len(wp2) == 1:
                ids = [id for id in e_ids if atom_sites_H[id].wp.letter == wp2[0].letter]
                if len(ids) == 1:
                    solution_template[i] = ids
                    assigned_ids.append(ids[0])

        # consider all permutations for to assign the rest atoms from H to G
        # https://stackoverflow.com/questions/65484940

        remaining_ids = [id for id in range(len(atom_sites_H)) if id not in assigned_ids]
        all_permutations = list(itertools.permutations(remaining_ids))
        unique_solutions = []
        if len(all_permutations)>max_num:
            all_permutations = sample(all_permutations, max_num)
        for permutation in all_permutations:
            permutation = list(permutation)
            solution = deepcopy(solution_template)
            for i, sol in enumerate(solution):
                valid = True
                if None in sol:
                    for j, id in enumerate(permutation[:len(sol)]):
                        if atom_sites_H[id].wp.letter == splitter.wp2_lists[i][j].letter:
                            solution[i][j] = id
                        else:
                            valid = False
                            break
                    if not valid:
                        break
                    else:
                        del permutation[:len(sol)]
            if valid and new_solution(solution, unique_solutions):
                unique_solutions.append(solution)
        return unique_solutions

    def symmetrize_dist(self, splitter, solution, disp=None, mask=None, d_tol=1.2):
        """
        For a given solution, search for the possbile supergroup structure

        Args:
            splitter: splitter object to specify the relation between G and H
            solution: list of sites in H, e.g., ['4a', '8b']
            disp: an overall shift from H to G, None or 3 vector
            mask: if need to freeze the direction
            d_tol: the tolerance in angstrom

        Returns:
            distortion
            cell translation
        """
        max_disps = []
        atom_sites_H = self.struc.atom_sites
        #n_atoms = sum([site.wp.multiplicity for site in atom_sites_H])

        # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
        for i, wp1 in enumerate(splitter.wp1_lists):
            if len(splitter.wp2_lists[i]) == 1:
                # symmetry info
                ops_H = splitter.H_orbits[i][0]  # ops for H
                matrix = splitter.G2_orbits[i][0][0].affine_matrix.copy() # op for G2

                change = True
                for mid in range(3):
                    # (x, 2x, 0) -> (x, y, z)
                    col = matrix[:3,mid]
                    if len(col[col==0])<2:
                        change=False
                        break
                if change:
                    # (x+1/4, 0, 0) -> (x, y, z)
                    # (z, 0, -x) -> (x1, y1, z1)
                    # transition between x and x+1/4 should be zero
                    for mid in range(3):
                        if np.linalg.norm(matrix[mid,:3])>0:
                            matrix[mid,:] = 0
                            matrix[mid,mid] = 1
                op_G2 = SymmOp(matrix)
                #print(change, op_G2.as_xyz_string(), ops_H[0].as_xyz_string())
                #import sys; sys.exit()

                # refine coord1 to find the best match on coord2
                coord = atom_sites_H[solution[i][0]].position
                coord0s = apply_ops(coord, ops_H) # possible coords in H
                dists = []
                for coord0 in coord0s:
                    coord2 = coord0.copy()
                    if disp is not None:
                        coord2 += disp
                    coord1 = op_G2.operate(coord2)
                    dist = coord1 - coord2
                    dist -= np.round(dist)
                    dist = np.dot(dist, self.cell)
                    dists.append(np.linalg.norm(dist))
                min_ID = np.argmin(np.array(dists))
                dist = dists[min_ID]
                coord2 = coord0s[min_ID].copy()
                if disp is None:
                    coord1 = op_G2.operate(coord2)
                    disp = (coord1 - coord2).copy()
                    # check if two are identical
                    mask = []
                    for m in range(3):
                        row1 = ops_H[0].affine_matrix[m]
                        row2 = op_G2.affine_matrix[m]
                        rot1 = row1[:3]
                        rot2 = row2[:3]
                        tran1 = row1[3]
                        tran2 = row2[3]
                        diff1 = np.sum((rot1-rot2)**2)
                        diff2 = tran1 - tran2
                        diff2 -= np.round(diff2)
                        diff2 *= diff2
                        if diff1 < 1e-3 and diff2 < 1e-3:
                            mask.append(m)
                elif dist < d_tol:
                    coord1 = op_G2.operate(coord2+disp)
                    if round(np.trace(op_G2.rotation_matrix)) in [1, 2]:
                        def fun(x, pt, ref, op):
                            pt[0] = x[0]
                            y = op.operate(pt)
                            diff = y - ref
                            diff -= np.round(diff)
                            diff = np.dot(diff, self.cell)
                            return np.linalg.norm(diff)

                        # optimize the distance by changing coord1
                        res = minimize(fun, coord1[0], args=(coord1, coord2, op_G2),
                                method='Nelder-Mead', options={'maxiter': 20})
                        coord1[0] = res.x[0]
                        coord1 = op_G2.operate(coord1)
                else:
                    return 10000, None, None

                if mask is not None:
                    disp[mask] = 0
                diff = coord1-(coord2+disp)
                diff -= np.round(diff)

                if np.linalg.norm(np.dot(diff, self.cell)) < d_tol:
                    max_disps.append(np.linalg.norm(np.dot(diff, self.cell)))
                else:
                    return 10000, None, None

            else:
                # symmetry operations
                ops_H1 = splitter.H_orbits[i][0]
                op_G21 = splitter.G2_orbits[i][0][0]
                ops_G22 = splitter.G2_orbits[i][1]
                #print(op_G21.as_xyz_string()) # (z, y, -x)
                #print(op_G22.as_xyz_string()) #
                #op_G2 = SymmOp(matrix)

                # atoms in H
                coord1 = atom_sites_H[solution[i][0]].position.copy()
                coord2 = atom_sites_H[solution[i][1]].position.copy()

                # refine coord1 to find the best match on coord2
                if disp is not None:
                    coord11 = coord1 + disp
                    coord22 = coord2 + disp
                else:
                    coord11 = coord1
                    coord22 = coord2

                if splitter.group_type == 'k':
                    # if it is t-type splitting
                    # the key is to restore the translation symmetry:
                    # e.g. (0.5, 0.5, 0.5), (0.5, 0, 0), .etc
                    # then find the best_match between coord11 and coord22
                    # regarding this vector, for example
                    # 8n in Immm -> 4f+4f in Pmmn
                    # x,y,0 -> -1/4+y,-1/4,-1/4+x -> 1/2+x1, 3/4, -z1
                    # -1/2+x,-1/2+y,-1/2 -> -1/4+y,-3/4,-x -> (x2, 1/4, z2 )
                    for op_G22 in ops_G22:
                        diff = (op_G22.rotation_matrix - op_G21.rotation_matrix).flatten()
                        if np.sum(diff**2) < 1e-3:
                            trans = op_G22.translation_vector - op_G21.translation_vector
                            break
                    trans -= np.round(trans)
                    coords11 = apply_ops(coord11, ops_H1)
                    coords11 += trans
                    tmp, dist = get_best_match(coords11, coord22, self.cell)
                    if dist > np.sqrt(2)*d_tol:
                        return 10000, None, mask

                    d = coord22 - tmp
                    d -= np.round(d)
                    max_disps.append(np.linalg.norm(np.dot(d/2, self.cell)))
                else:
                    op_G22 = ops_G22[0]
                    if round(np.trace(ops_H1[0].rotation_matrix)) == 0:
                        # if there is no freedom e.g. (0.5, 0, 0)
                        # no displacement to adjust
                        # G1/H1 and G2/H2 have to be identical
                        coord2 = op_G21.operate(coord2) #
                        coord_G1 = op_G21.operate(coord11)
                        coord_G2 = op_G22.operate(coord22)
                        dist1 = coord_G1 - coord11
                        dist2 = coord_G2 - coord22
                        dist1 -= np.round(dist1)
                        dist2 -= np.round(dist2)
                        dist1 = np.linalg.norm(dist1)
                        dist2 = np.linalg.norm(dist2)
                        if max([dist1, dist2]) > np.sqrt(2)*d_tol:
                            return 10000, None, mask
                        else:
                            max_disps.append(max([dist1, dist2]))

                    elif round(np.trace(ops_H1[0].rotation_matrix)) == 1:
                        # 24e in Fm-3m -> 4e+8h in I4/mmm
                        # x,0,0 -> x,x,0 -> x1,x1,0
                        # 0,0,-x -> 0,0,-x -> 0,0,z2
                        # only one freedom exists
                        coords11 = apply_ops(coord11, ops_H1)
                        for m, coord11 in enumerate(coords11):
                            coords11[m] = op_G22.operate(coord11)
                        tmp, dist = get_best_match(coords11, coord22, self.cell)
                        if dist > np.sqrt(2)*d_tol:
                            return 10000, None, mask

                        # recover the original position
                        try:
                            inv_op1 = get_inverse(op_G21)
                            inv_op2 = get_inverse(op_G22)
                            #print(op.as_xyz_string(), inv_op.as_xyz_string())
                        except:
                            print("Error in getting the inverse")
                            print(op_G21.as_xyz_string())
                            print(op_G22.as_xyz_string())
                            #import sys; sys.exit()

                        coord1 = inv_op2.operate(tmp)
                        coord2 = inv_op1.operate(coord2)
                        if disp is not None:
                            coord1 -= disp
                        d = coord22 - tmp
                        d -= np.round(d)

                        #coord22 -= d/2 #final coord2 after disp
                        #recover the displaced position
                        #coord11 = inv_op2.operate(coord22)
                        max_disps.append(np.linalg.norm(np.dot(d/2, self.cell)))

                    else:
                        # 8j in C2m -> 4e+4e in P21/c
                        # (z, y, -x) -> (x1, y1, z1)
                        # (z, 1/2+y, 1/2-x) -> (x2, y2, z2)
                        coord2 = op_G21.operate(coord2) #
                        coords11 = apply_ops(coord11, ops_H1)
                        # transform coords1 by symmetry operation
                        for m, coord11 in enumerate(coords11):
                            coords11[m] = op_G22.operate(coord11)
                        tmp, dist = get_best_match(coords11, coord22, self.cell)
                        #if wp1.letter == 'e': print(tmp, coord22, dist)

                        if dist > np.sqrt(2)*d_tol:
                            return 10000, None, mask

                        # recover the original position
                        try:
                            inv_op1 = get_inverse(op_G21)
                            inv_op2 = get_inverse(op_G22)
                            #print(op.as_xyz_string(), inv_op.as_xyz_string())
                        except:
                            print("Error in getting the inverse")
                            print(op_G21.as_xyz_string())
                            print(op_G22.as_xyz_string())
                            #import sys; sys.exit()
                        coord1 = inv_op2.operate(tmp)
                        coord2 = inv_op1.operate(coord2)
                        if disp is not None:
                            coord1 -= disp
                        d = coord22 - tmp
                        d -= np.round(d)

                        #coord22 -= d/2 #final coord2 after disp
                        # recover the displaced position
                        #coord11 = inv_op2.operate(coord22)

                        max_disps.append(np.linalg.norm(np.dot(d/2, self.cell)))
        #print(max_disps, mask)
        return max(max_disps), disp, mask

    def symmetrize(self, splitter, solution, disp):
        """
        For a given solution, search for the possbile supergroup structure

        Args:
            splitter: splitter object to specify the relation between G and H
            disp: an overall shift from H to G, None or 3 vector
            d_tol: the tolerance in angstrom

        Returns:
            coords_G1: coordinates in G
            coords_G2: coordinates in G under the subgroup setting
            coords_H1: coordinates in H
            elements: list of elements
        """

        atom_sites_H = self.struc.atom_sites
        coords_G1 = [] # position in G
        coords_G2 = [] # position in G on the subgroup bais
        coords_H1 = [] # position in H
        elements = []
        inv_rot = splitter.R[:3,:3] # inverse coordinate transformation
        inv_tran = splitter.R[:3,3] # needs to check
        # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
        for i, wp1 in enumerate(splitter.wp1_lists):

            if len(splitter.wp2_lists[i]) == 1:
                # symmetry info
                ops_H = splitter.H_orbits[i][0]  # ops for H
                matrix = splitter.G2_orbits[i][0][0].affine_matrix.copy() # ops for G2

                change = True
                for mid in range(3):
                    # (x, 2x, 0) -> (x, y, z)
                    col = matrix[:3,mid]
                    if len(col[col==0])<2:
                        change=False
                        break

                if change:
                    # (x+1/4, 0, 0) -> (x, y, z)
                    # (z, 0, -x) -> (x1, y1, z1)
                    # (x, 2x, 1/4) -> (x, y, 1/4) # will fail if x appears twice
                    # transition between x and x+1/4 should be zero
                    for mid in range(3):
                        if np.linalg.norm(matrix[mid,:3])>0:
                            matrix[mid,:] = 0
                            matrix[mid,mid] = 1
                else: #if splitter.G.number == 92:
                    # temp fix change from x-1/4, x, z to x, x+1/4, z
                    matrix[:2,3] -= matrix[0,3]

                op_G2 = SymmOp(matrix)
                # refine coord1 to find the best match on coord2
                coord = atom_sites_H[solution[i][0]].position
                coord0s = apply_ops(coord, ops_H) # possible coords in H
                dists = []
                for coord0 in coord0s:
                    coord2 = coord0.copy() + disp
                    coord1 = op_G2.operate(coord2) # coord in G
                    dist = coord1 - coord2
                    dist -= np.round(dist)
                    dist = np.dot(dist, self.cell)
                    dists.append(np.linalg.norm(dist))

                min_ID = np.argmin(np.array(dists))
                coord2 = coord0s[min_ID].copy()
                coord1 = op_G2.operate(coord2+disp)

                #if splitter.G.number == 92:
                #    print(change, op_G2.as_xyz_string(), "++++++++++++++")
                #    print(coord1, coord2, np.dot(inv_rot, coord1).T+inv_tran.T, "--------")

                if round(np.trace(op_G2.rotation_matrix)) in [1, 2]:
                    def fun(x, pt, ref, op):
                        pt[0] = x[0]
                        y = op.operate(pt)
                        diff = y - ref
                        diff -= np.round(diff)
                        diff = np.dot(diff, self.cell)
                        return np.linalg.norm(diff)

                    # optimize the distance by changing coord1
                    res = minimize(fun, coord1[0], args=(coord1, coord2, op_G2),
                            method='Nelder-Mead', options={'maxiter': 20})
                    coord1[0] = res.x[0]
                    coord1 = op_G2.operate(coord1)

                coords_G1.append(np.dot(inv_rot, coord1).T+inv_tran.T)
                coords_G2.append(coord1)
                coords_H1.append(coord2)
                elements.append(splitter.elements[i])

            else:
                # symmetry operations
                ops_H1 = splitter.H_orbits[i][0]
                #ops_H2 = splitter.H_orbits[i][0]
                op_G1  = splitter.G1_orbits[i][0][0]
                op_G21 = splitter.G2_orbits[i][0][0]
                ops_G22 = splitter.G2_orbits[i][1]

                # atoms in H
                coord1 = atom_sites_H[solution[i][0]].position.copy()
                coord2 = atom_sites_H[solution[i][1]].position.copy()

                if splitter.group_type == 'k':
                    # if it is t-type splitting
                    # the key is to restore the translation symmetry:
                    # e.g. (0.5, 0.5, 0.5), (0.5, 0, 0), .etc
                    # then find the best_match between coord11 and coord22
                    # regarding this vector, for example
                    # 8n in Immm -> 4f+4f in Pmmn
                    # x,y,0 -> -1/4+y,-1/4,-1/4+x -> 1/2+x1, 3/4, -z1
                    # -1/2+x,-1/2+y,-1/2 -> -1/4+y,-3/4,-x -> (x2, 1/4, z2 )
                    coord11 = coord1 + disp
                    coord22 = coord2 + disp

                    for op_G22 in ops_G22:
                        diff = (op_G22.rotation_matrix - op_G21.rotation_matrix).flatten()
                        if np.sum(diff**2) < 1e-3:
                            trans = op_G22.translation_vector - op_G21.translation_vector
                            break
                    trans -= np.round(trans)
                    coords11 = apply_ops(coord11, ops_H1)
                    coords11 += trans
                    tmp, dist = get_best_match(coords11, coord22, self.cell)

                    d = coord22 - tmp
                    d -= np.round(d)
                    coord22 -= d/2
                    coord11 += d/2

                    coords_G1.append(np.dot(inv_rot, coord22).T+inv_tran.T)
                    coords_G2.append(coord11)
                    coords_G2.append(coord22)
                    coords_H1.append(coord1)
                    coords_H1.append(coord2)

                else:
                    op_G22 = ops_G22[0]
                    if round(np.trace(ops_H1[0].rotation_matrix)) == 0:
                        # if there is no freedom e.g. (0.5, 0, 0)
                        # only check between
                        coord11 = coord1 + disp
                        coord22 = coord2 + disp
                        coords_G1.append(op_G1.operate(coord11))
                        coords_G2.append(op_G21.operate(coord11))
                        coords_G2.append(op_G22.operate(coord22))
                        coords_H1.append(coord1)
                        coords_H1.append(coord2)

                    elif round(np.trace(ops_H1[0].rotation_matrix)) == 1:
                        coord11 = coord1 + disp
                        coord22 = coord2 + disp
                        coords11 = apply_ops(coord11, ops_H1)
                        for m, coord in enumerate(coords11):
                            coords11[m] = op_G22.operate(coord)
                        tmp, dist = get_best_match(coords11, coord22, self.cell)

                        d = coord22 - tmp
                        d -= np.round(d)

                        coord22 -= d/2
                        coord11 += d/2

                        coords_G1.append(np.dot(inv_rot, coord22).T+inv_tran.T)
                        coords_G2.append(coord11)
                        coords_G2.append(coord22)
                        coords_H1.append(coord1)
                        coords_H1.append(coord2)

                    else:
                        coord2 = op_G21.operate(coord2) #
                        coord11 = coord1 + disp
                        coord22 = coord2 + disp
                        coords11 = apply_ops(coord11, ops_H1)
                        # transform coords1 by symmetry operation
                        for m, coord in enumerate(coords11):
                            coords11[m] = op_G22.operate(coord)
                        tmp, dist = get_best_match(coords11, coord22, self.cell)
                        # recover the original position
                        #inv_op2 = get_inverse(op_G22)
                        #coord1 = inv_op2.operate(tmp)
                        #coord1 -= disp
                        d = coord22 - tmp
                        d -= np.round(d)
                        coord22 -= d/2 #final coord2 after disp
                        coord11 += d/2

                        coords_G1.append(np.dot(inv_rot, coord22).T+inv_tran.T)
                        coords_G2.append(coord11)
                        coords_G2.append(coord22)
                        coords_H1.append(coord1)
                        coords_H1.append(coord2)
                        #print(splitter)
                elements.extend([splitter.elements[i]]*2)
        return coords_G1, coords_G2, coords_H1, elements

    def print_detail(self, G, coords_H1, coords_G1, elements, disp):
        """
        print out the details of tranformation
        """
        print("Valid structure", G)
        disps = []
        for x, y, ele in zip(coords_H1, coords_G1, elements):
            dis = y-(x+disp)
            dis -= np.round(dis)
            dis_abs = np.linalg.norm(dis.dot(self.cell))
            output = "{:2s} {:8.4f}{:8.4f}{:8.4f}".format(ele, *x)
            output += " -> {:8.4f}{:8.4f}{:8.4f}".format(*y)
            output += " -> {:8.4f}{:8.4f}{:8.4f} {:8.4f}".format(*dis, dis_abs)
            disps.append(dis_abs)
            print(output)
        print("cell: {:8.4f}{:8.4f}{:8.4f}, disp (A): {:8.4f}".format(*disp, max(disps)))



if __name__ == "__main__":

    from pyxtal import pyxtal

    data = {
            "NbO2": [141],
            "GeF2": [62],
            "lt_quartz": [180],
            "BTO": [123, 221],
            "NaSb3F10": [186,194],
            "BTO-Amm2": [65, 123, 221],
            "lt_cristobalite": [98, 210, 227],
            "MPWO": [59, 71, 139, 225],
           }
    # "PVO: 14-12-
    # "PPO: 15-12-?---bad example?
    cif_path = "pyxtal/database/cifs/"

    for cif in data.keys():
        print("===============", cif, "===============")
        s = pyxtal()
        s.from_seed(cif_path+cif+'.cif')
        sup = supergroups(s, path=data[cif])
        print(sup)
