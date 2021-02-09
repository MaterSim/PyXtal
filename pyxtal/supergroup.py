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
from pyxtal.operations import apply_ops
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

def find_mapping(atom_sites, splitter, max_num=720):
    """
    search for all mappings for a given splitter

    Args:
        atom_sites: list of wyc object
        splitter: wyc_splitter object
        max_num (int): maximum number of atomic mapping

    Returns:
        unique solutions
    """
    solution_template = [[None]*len(wp2) for wp2 in splitter.wp2_lists]
    assigned_ids = []
    # look for unique assignment from sites_H to sites_G
    for i, wp2 in enumerate(splitter.wp2_lists):
        # choose the sites belong to the same element

        ele = splitter.elements[i]
        e_ids = [id for id, site in enumerate(atom_sites) if site.specie==ele]
        wp_letters = set([wp.letter for wp in wp2])
        ids = [id for id in e_ids if atom_sites[id].wp.letter in wp_letters]
        if len(wp2) == len(ids):
            solution_template[i] = ids
            assigned_ids.extend(ids)
        else:
            if len(wp_letters) == len(wp2): # eg, ['4a', '2b']
                for j, wp in enumerate(wp2):
                    ids = [id for id in ids if atom_sites[id].wp.letter == wp_letters]
                    if len(ids) == 1:
                        solution_template[i][j] = ids[0]
                        assigned_ids.extend(ids)

    # print(assigned_ids, solution_template)
    # consider all permutations for to assign the rest atoms from H to G
    # https://stackoverflow.com/questions/65484940
    #print(solution_template)
    # splitting them by elements
    eles = set([site.specie for site in atom_sites])
    lists = []
    eles0 = []
    for ele in eles:
        ids = [id for id, site in enumerate(atom_sites) if id not in assigned_ids and site.specie==ele]
        # split them by sites e.g. ['4a', '4a'] will be no permutation
        if len(ids)>0:
            lists.append(list(itertools.permutations(ids)))
            eles0.append(ele)

    all_permutations = list(itertools.product(*lists))
    unique_solutions = []

    if len(all_permutations) > max_num:
        print("Warning: ignore some mapping: ", str(len(all_permutations)-max_num))
        all_permutations = sample(all_permutations, max_num)

    for permutation in all_permutations:
        solution = deepcopy(solution_template)
        valid = True
        for perm, ele in zip(permutation, eles0):
            perm = list(perm)
            for i, sol in enumerate(solution):
                if None in sol and splitter.elements[i]==ele:
                    for j, wp in enumerate(splitter.wp2_lists[i]):
                        if wp.letter == atom_sites[perm[j]].wp.letter:
                            solution[i][j] = perm[j]
                        else:
                            valid = False
                            break
                    if not valid:
                        break
                    else:
                        del perm[:len(splitter.wp2_lists[i])]
        if valid and new_solution(solution, unique_solutions):
            unique_solutions.append(solution)
    #print("mapping", len(unique_solutions), len(all_permutations))
    #import sys; sys.exit()
    return unique_solutions


def search_match_with_trans(rot, tran, pos1, pos2, cell=None, inverse=False):
    """
    apply symmetry operation on pos1 when it involves cell change.
    e.g., when the transformation is (a+b, a-b, c),
    trial translation needs to be considered to minimize the
    difference between the transformed pos1 and reference pos2
    """
    pos1 -= np.round(pos1)
    if np.linalg.det(rot) < 1:
        shifts = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])
    else:
        shifts = np.array([[0,0,0]])

    dists = []
    for shift in shifts:
        if inverse:
            res = np.dot(rot, pos1 + shift + tran.T)
        else:
            res = np.dot(rot, pos1 + shift) + tran.T

        diff = res - pos2
        diff -= np.round(diff)
        dists.append(np.linalg.norm(diff))

    dists = np.array(dists)
    dist = np.min(dists)
    shift = shifts[np.argmin(dists)]
    if inverse:
        pos = np.dot(rot, pos1 + shift + tran.T)
    else:
        pos = np.dot(rot, pos1 + shift) + tran.T

    diff = pos - pos2
    diff -= np.round(diff)

    if cell is not None:
        diff = np.dot(diff, cell)

    dist = np.linalg.norm(diff)

    return pos, dist


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
    matrix = np.dot(trans.T, struc.lattice.get_matrix())
    l1 = Lattice.from_matrix(matrix)
    l2 = Lattice.from_matrix(matrix, ltype=sym.Group(G).lattice_type)
    (a1,b1,c1,alpha1,beta1,gamma1)=l1.get_para(degree=True)
    (a2,b2,c2,alpha2,beta2,gamma2)=l2.get_para(degree=True)
    abc_diff = np.abs(np.array([a2-a1, b2-b1, c2-c1])).max()
    ang_diff = np.abs(np.array([alpha2-alpha1, beta2-beta1, gamma2-gamma1])).max()
    #print(l1, l2)
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

def search_paths(H, G, max_layers=5):
    """
    Search function throws away paths that take a roundabout. if
    path1:a>>e>>f>>g
    path2:a>>b>>c>>e>>f>>g
    path 2 will not be counted as there is already a shorter path from a>>e

    Args:
        H: starting structure IT group number
        G: final supergroup IT Group number
        max_layers: the number of supergroup calculations needed.

    Return:
        list of possible paths ordered from smallest to biggest


    """

    layers={}
    layers[0]={'groups':[G],'subgroups':[]}
    final=[]
    traversed=[]

    # searches for every subgroup of the the groups from the previous layer.
    # Stores the possible groups of each layer and their subgroups 
    # in a dictinoary to avoid redundant calculations.
    # Starts from G and goes down to H
    for l in range(1,max_layers+1):
        previous_layer_groups=layers[l-1]['groups']
        groups=[]
        subgroups=[]
        for g in previous_layer_groups:
            subgroup_numbers=np.unique(sym.Group(g).get_max_subgroup_numbers())

            # If a subgroup list has been found with H, will trace 
            # a path through the dictionary to build the path
            if H in subgroup_numbers:
                paths=[[H,g]]
                for j in reversed(range(l-1)):
                    holder=[]
                    for path in paths:
                        tail_number=path[-1]
                        indices=[]
                        for idx, numbers in enumerate(layers[j]['subgroups']):
                            if tail_number in numbers:
                                indices.append(idx)
                        for idx in indices:
                            holder.append(path+[layers[j]['groups'][idx]])
                    paths=deepcopy(holder)
                final.extend(paths)
                subgroups.append([])

            #will continue to generate a layer of groups if the path to H has not been found.
            else:
                subgroups.append(subgroup_numbers)
                [groups.append(x) for x in subgroup_numbers if (x not in groups) and (x not in traversed)]

        traversed.extend(groups)
        layers[l]={'groups':deepcopy(groups),'subgroups':[]}
        layers[l-1]['subgroups']=deepcopy(subgroups)
    return final

def new_path(path, paths):
    """
    check if struc is already in the reference solutions
    """
    for ref in paths:
        if path[:len(ref)] == ref:
            return False
    return True

class supergroups():
    """
    Class to search for the feasible transition to a given super group

    Args:
        struc: pyxtal structure
        G (int): the desired super group number
        path: the path to connect G and H, e.g, [62, 59, 74]
        d_tol (float): tolerance for largest atomic displacement
        show (bool): whether or not show the detailed process
    """

    def __init__(self, struc, G=None, path=None, d_tol=1.0, max_per_G=100, show=False):
        self.struc0 = struc
        self.show = show
        self.d_tol = d_tol
        self.max_per_G = max_per_G
        if path is None:
            paths = search_paths(struc.group.number, G, max_layers=5)
        else:
            paths = [path]
        
        print("{:d} paths will be checked".format(len(paths)))
        self.strucs = None
        failed_paths = []
        for i, p in enumerate(paths):
            status = "path{:2d}: {:s}, ".format(i, str(p))
            if new_path(p, failed_paths):
                strucs, w_path, valid = self.struc_along_path(p)
                status += "stops at: {:s}".format(str(w_path))
                if valid:
                    self.strucs = strucs
                    if len(strucs) > len(p):
                        self.path = [self.struc0.group.number] + p
                    else:
                        self.path = p
                    break
                else:
                    failed_paths.append(w_path)
            else:
                status += "skipped..."
            print(status)


    def __str__(self):
        s = "\nTransition to super group: "
        if self.strucs is None:
            s += "Unsuccessful, check your input"
        else:
            s += "{:d}".format(self.path[0])
            for i, p in enumerate(self.path[1:]):
                s += " -> {:d}[{:5.3f}]".format(p, self.strucs[i+1].disp)
            s += '\n'
            #for struc in self.strucs:
            #    s += str(struc)
        return s

    def __repr__(self):
        return str(self)

    def struc_along_path(self, path):
        """
        search for the super group structure along a given path
        """
        strucs = []
        G_strucs = [self.struc0]
        working_path = []
        for G in path:
            working_path.append(G)
            H = G_strucs[0].group.number
            if G != H:
                if sym.get_point_group(G) == sym.get_point_group(H):
                    group_type = 'k'
                else:
                    group_type = 't'

                for G_struc in G_strucs:
                    my = supergroup(G_struc, [G], group_type)
                    solutions = my.search_supergroup(self.d_tol, self.max_per_G)
                    new_G_strucs = my.make_supergroup(solutions, show_detail=self.show)
                    if len(new_G_strucs) > 0:
                        strucs.append(G_struc)
                        G_strucs = new_G_strucs
                        break
                if len(new_G_strucs) == 0:
                    break

        # add the final struc
        if len(new_G_strucs) > 0:
            ds = [st.disp for st in new_G_strucs]
            minID = np.argmin(np.array(ds))
            strucs.append(new_G_strucs[minID])
            valid = True
        else:
            valid = False
        return strucs, working_path, valid

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
                    print("Warning: ignore some solutions: ", len(sols)-max_per_G)
                    sols = sample(sols, max_per_G)
                #sols = [[['2f'], ['1a'], ['4n']]]
                #if G==59: sols=[(['2a'], ['4f'], ['2b'], ['4e', '4e', '4f'])]
                for sol in sols:
                    #print(sol)
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
                G = sp.G.number
                lat1 = np.dot(np.linalg.inv(sp.R[:3,:3]).T, self.struc.lattice.matrix)
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
                        print(">>>>>>>>>>>>>>")
                        print(pos)
                        print(wp)
                        print(">>>>>>>>>>>>>>")
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
        mappings = find_mapping(self.struc.atom_sites, splitter)
        dists = []
        disps = []
        masks = []
        if len(mappings) > 0:
            for mapping in mappings:
                dist, disp, mask = self.symmetrize_dist(splitter, mapping, None, None, d_tol)
                dists.append(dist)
                disps.append(disp)
                masks.append(mask)
            dists = np.array(dists)
            mae = np.min(dists)
            id = np.argmin(dists)
            disp = disps[id]
            mask = masks[id]
            if 0.2 < mae < d_tol:
                # optimize disp further
                if mask is None or len(mask)<3:
                    def fun(disp, mapping, splitter, mask):
                        return self.symmetrize_dist(splitter, mapping, disp, mask)[0]
                    res = minimize(fun, disps[id], args=(mappings[id], splitter, mask),
                            method='Nelder-Mead', options={'maxiter': 10})
                    if res.fun < mae:
                        mae = res.fun
                        disp = res.x
            return mae, disp, mappings[id], splitter
        else:
            print("bug in findding the mappings", solution)
            print(splitter.G.number, '->', splitter.H.number)

            return 1000, None, None, None


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


    def symmetrize_dist(self, splitter, mapping, disp=None, mask=None, d_tol=1.2):
        """
        For a given solution, search for the possbile supergroup structure

        Args:
            splitter: splitter object to specify the relation between G and H
            mapping: list of sites in H, e.g., ['4a', '8b']
            disp: an overall shift from H to G, None or 3 vector
            mask: if need to freeze the direction
            d_tol: the tolerance in angstrom

        Returns:
            distortion
            cell translation
        """
        cell = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.struc.lattice.matrix)
        max_disps = []
        atom_sites_H = self.struc.atom_sites
        rot = splitter.R[:3,:3] # inverse coordinate transformation
        tran = splitter.R[:3,3] # needs to check
        inv_rot = np.linalg.inv(rot)
        cell = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.struc.lattice.matrix)
        ops_G1  = splitter.G[0]
        #n_atoms = sum([site.wp.multiplicity for site in atom_sites_H])

        if mask is not None and disp is not None:
            disp[mask] = 0

        # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
        for i, wp1 in enumerate(splitter.wp1_lists):
            if len(splitter.wp2_lists[i]) == 1:
                op_G1 = splitter.G1_orbits[i][0][0]
                ops_H = splitter.H_orbits[i][0]
                base = atom_sites_H[mapping[i][0]].position.copy() 

                #choose the best coord1_H
                coord1s_H = apply_ops(base, ops_H)
                ds = []
                for coord1_H in coord1s_H:
                    if disp is not None:
                        coord1_G2 = coord1_H + disp
                    else:
                        coord1_G2 = coord1_H 

                    if np.linalg.det(rot) < 1:
                        shifts = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])
                    else:
                        shifts = np.array([[0,0,0]])

                    diffs = []
                    for shift in shifts:
                        res = np.dot(rot, coord1_G2+shift) + tran.T
                        tmp = sym.search_cloest_wp(splitter.G, wp1, op_G1, res)
                        diff = res - tmp
                        diff -= np.round(diff)
                        diffs.append(np.linalg.norm(diff))
                    diffs = np.array(diffs)
                    ds.append([np.min(diffs), np.argmin(diffs)])

                ds = np.array(ds)
                minID = np.argmin(ds[:,0])
                coord1_H = coord1s_H[minID]

                if disp is not None:
                    coord1_G2 = coord1_H + disp
                else:
                    coord1_G2 = coord1_H 

                coord1_G1 = np.dot(rot, coord1_G2+shifts[int(ds[minID, 1])]) + tran.T
                tmp = sym.search_cloest_wp(splitter.G, wp1, op_G1, coord1_G1)

                # initial guess on disp
                if disp is None:
                    coord1_G2, dist1 = search_match_with_trans(inv_rot, -tran, tmp, coord1_H, None, True)
                    diff = coord1_G2 - coord1_H
                    diff -= np.round(diff)
                    disp = diff.copy()
                    mask = []
                    for m in range(3):
                        if abs(diff[m])<1e-4:
                            mask.append(m)
                    dist = 0
                else:
                    coord1_G2, dist = search_match_with_trans(inv_rot, -tran, tmp, coord1_H+disp, self.cell, True)

                if dist < d_tol:
                    #if dist >0: print(dist, coord1_G2, coord1_H+disp)
                    max_disps.append(dist)
                else:
                    #import sys; sys.exit()
                    #print(wp1.letter, tmp, coord1_G2, coord1_H+disp, dist)
                    return 10000, None, None

            else:
                # assume zero shift, needs to check
                if disp is None:
                    disp = np.zeros(3)
                    mask = [0, 1, 2]

                # H->G2->G1
                if atom_sites_H[mapping[i][0]].wp.letter == splitter.wp2_lists[i][0].letter:
                    coord1_H = atom_sites_H[mapping[i][0]].position.copy()
                    coord2_H = atom_sites_H[mapping[i][1]].position.copy()
                else:
                    coord2_H = atom_sites_H[mapping[i][0]].position.copy()
                    coord1_H = atom_sites_H[mapping[i][1]].position.copy()

                #print("\n\n\nH", coord1_H, coord2_H)

                coord1_G2 = coord1_H + disp
                coord2_G2 = coord2_H + disp

                if splitter.group_type == 'k':
                    # For t-type splitting, restore the translation symmetry:
                    # e.g. (0.5, 0.5, 0.5), (0.5, 0, 0), .etc
                    # then find the best_match between coord1 and coord2, 

                    ops_H1 = splitter.H_orbits[i][0]
                    op_G21 = splitter.G2_orbits[i][0][0]
                    ops_G22 = splitter.G2_orbits[i][1]
                    for op_G22 in ops_G22:
                        diff = (op_G22.rotation_matrix - op_G21.rotation_matrix).flatten()
                        if np.sum(diff**2) < 1e-3:
                            trans = op_G22.translation_vector - op_G21.translation_vector
                            break
                    trans -= np.round(trans)
                    coords11 = apply_ops(coord1_G2, ops_H1)
                    coords11 += trans
                    tmp, dist = get_best_match(coords11, coord2_G2, self.cell)
                    if dist > np.sqrt(2)*d_tol:
                        #print("kkkk", dist, coords11, coord2_G2)
                        return 10000, None, mask
                    else:
                        d = coord2_G2 - tmp
                        d -= np.round(d)
                        max_disps.append(np.linalg.norm(np.dot(d/2, self.cell)))

                else:
                    #print("disp", disp)
                    #print("G2", coord1_G2, coord2_G2)
                    op_G12 = splitter.G1_orbits[i][1][0]
                    coord1_G1 = np.dot(rot, coord1_G2) + tran.T
                    coord2_G1 = np.dot(rot, coord2_G2) + tran.T
                    coord1_G1 -= np.round(coord1_G1)
                    coord2_G1 -= np.round(coord2_G1)
                    #print("G1", coord1_G1, coord2_G1, op_G12.as_xyz_string())
                    #print(splitter.G.number, splitter.wp1_lists[i].index, coord2_G1, op_G12.as_xyz_string())
                    coord2_G1 = sym.search_cloest_wp(splitter.G, wp1, op_G12, coord2_G1)
                    #print("G1(symm1)", coord1_G1, coord2_G1)
                    #import sys; sys.exit()
                    #find the best match
                    coords11 = apply_ops(coord1_G1, ops_G1)
                    tmp, dist = get_best_match(coords11, coord2_G1, cell)
                    #print(coords11)
                    #print("tmp_raw:", tmp, coord2_G1, "dist", dist)
                    tmp = sym.search_cloest_wp(splitter.G, wp1, op_G12, tmp)

                    # G1->G2->H
                    d = coord2_G1 - tmp
                    d -= np.round(d)
                    #print("dist", np.linalg.norm(d), "d", d, tmp, coord2_G1)

                    coord2_G1 -= d/2
                    coord1_G1 += d/2
                    #print("G1 (symm2)", coord1_G1, coord2_G1)
                    
                    coord1_G2, dist1 = search_match_with_trans(inv_rot, -tran, coord1_G1, coord1_H+disp, self.cell, True)
                    coord2_G2, dist2 = search_match_with_trans(inv_rot, -tran, coord2_G1, coord2_H+disp, self.cell, True)

                    max_disps.append(max([dist1, dist2]))
                    #print("1:", coord1_G2, coord1_H, dist1)
                    #print("2:", coord2_G2, coord2_H, dist2)
        #if max(max_disps)<0.1: print(max_disps)
        return max(max_disps), disp, mask


    def symmetrize(self, splitter, mapping, disp):
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

        cell = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.struc.lattice.matrix)
        atom_sites_H = self.struc.atom_sites
        coords_G1 = [] # position in G
        coords_G2 = [] # position in G on the subgroup bais
        coords_H1 = [] # position in H
        elements = []
        rot = splitter.R[:3,:3] # inverse coordinate transformation
        tran = splitter.R[:3,3] # needs to check
        inv_rot = np.linalg.inv(rot)
        ops_G1  = splitter.G[0]
        # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
        for i, wp1 in enumerate(splitter.wp1_lists):

            if len(splitter.wp2_lists[i]) == 1:
                op_G1 = splitter.G1_orbits[i][0][0]
                ops_H = splitter.H_orbits[i][0]
                base = atom_sites_H[mapping[i][0]].position.copy() 

                #choose the best coord1_H
                coord1s_H = apply_ops(base, ops_H)
                ds = []
                for coord1_H in coord1s_H:
                    coord1_G2 = coord1_H + disp
                    if np.linalg.det(rot) < 1:
                        shifts = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])
                    else:
                        shifts = np.array([[0,0,0]])

                    diffs = []
                    for shift in shifts:
                        res = np.dot(rot, coord1_G2+shift) + tran.T
                        tmp = sym.search_cloest_wp(splitter.G, wp1, op_G1, res)
                        diff = res - tmp
                        diff -= np.round(diff)
                        diffs.append(np.linalg.norm(diff))
                    diffs = np.array(diffs)
                    ds.append([np.min(diffs), np.argmin(diffs)])

                ds = np.array(ds)
                minID = np.argmin(ds[:,0])
                coord1_H = coord1s_H[minID]

                coord1_G2 = coord1_H + disp
                coord1_G1 = np.dot(rot, coord1_G2) + tran.T
                tmp = sym.search_cloest_wp(splitter.G, wp1, op_G1, coord1_G1)

                coords_G1.append(tmp)
                coord1_G2, dist = search_match_with_trans(inv_rot, -tran, tmp, coord1_H+disp, self.cell, True)
 
                coords_G2.append(coord1_G2)
                coords_H1.append(coord1_H)

            else:

                if atom_sites_H[mapping[i][0]].wp.letter == splitter.wp2_lists[i][0].letter:
                    coord1_H = atom_sites_H[mapping[i][0]].position.copy() 
                    coord2_H = atom_sites_H[mapping[i][1]].position.copy()
                else:
                    coord2_H = atom_sites_H[mapping[i][0]].position.copy() 
                    coord1_H = atom_sites_H[mapping[i][1]].position.copy()

                coord1_G2 = coord1_H + disp
                coord2_G2 = coord2_H + disp  

                if splitter.group_type == 'k':
                    ops_H1 = splitter.H_orbits[i][0]
                    op_G21 = splitter.G2_orbits[i][0][0]
                    ops_G22 = splitter.G2_orbits[i][1]

                    coord11 = coord1_H + disp
                    coord22 = coord2_H + disp

                    for op_G22 in ops_G22:
                        diff = (op_G22.rotation_matrix - op_G21.rotation_matrix).flatten()
                        if np.sum(diff**2) < 1e-3:
                            trans = op_G22.translation_vector - op_G21.translation_vector
                            break
                    trans -= np.round(trans)
                    coords11 = apply_ops(coord1_G2, ops_H1)
                    coords11 += trans
                    tmp, dist = get_best_match(coords11, coord2_G2, self.cell)

                    d = coord2_G2 - tmp
                    d -= np.round(d)
                    coord2_G2 -= d/2
                    coord1_G2 += d/2
                    coords_G1.append(np.dot(rot, coord2_G2).T + tran.T)

                else:
                    # H->G2->G1
                    op_G12 = splitter.G1_orbits[i][1][0]

                    coord1_G1 = np.dot(rot, coord1_G2) + tran.T
                    coord2_G1 = np.dot(rot, coord2_G2) + tran.T
                    coord1_G1 -= np.round(coord1_G1)
                    coord2_G1 -= np.round(coord2_G1)

                    coord2_G1 = sym.search_cloest_wp(splitter.G, wp1, op_G12, coord2_G1)
                    
                    #find the best match
                    coords11 = apply_ops(coord1_G1, ops_G1)
                    tmp, dist = get_best_match(coords11, coord2_G1, cell)
                    tmp = sym.search_cloest_wp(splitter.G, wp1, op_G12, tmp)
 
                    # G1->G2->H
                    d = coord2_G1 - tmp
                    d -= np.round(d)
                    coord2_G1 -= d/2
                    coord1_G1 += d/2
                    coords_G1.append(coord2_G1)

                    coord1_G2 = np.dot(inv_rot, coord1_G1 - tran.T)
                    coord2_G2 = np.dot(inv_rot, coord2_G1 - tran.T)

                coords_G2.append(coord1_G2)
                coords_G2.append(coord2_G2)
                coords_H1.append(coord1_H)
                coords_H1.append(coord2_H)

            elements.extend([splitter.elements[i]]*len(splitter.wp2_lists[i]))

        return coords_G1, coords_G2, coords_H1, elements

if __name__ == "__main__":

    from pyxtal import pyxtal
    from time import time
    data = {
            #"PVO": [12, 166],
            #"PPO": [12],
            #"BTO": [123, 221],
            #"lt_cristobalite": [98, 210, 227],
            #"MPWO": [59, 71, 139, 225],
            #"BTO-Amm2": [65, 123, 221],
            #"NaSb3F10": [186, 194],
            "GeF2": 62,
            "NbO2": 141,
            "NiS-Cm": 160,
            "lt_quartz": 180,
            "BTO-Amm2": 221,
            "BTO": 221,
            "MPWO": 225,
            "lt_cristobalite": 227,
            "NaSb3F10": 194,
           }
    cif_path = "pyxtal/database/cifs/"

    for cif in data.keys():
        t0 = time()
        print("===============", cif, "===============")
        s = pyxtal()
        s.from_seed(cif_path+cif+'.cif')
        if isinstance(data[cif], list):
            #sup = supergroups(s, path=data[cif], show=False, max_per_G=2500)
            sup = supergroups(s, path=data[cif], show=True, max_per_G=2500)
        else:
            sup = supergroups(s, G=data[cif], show=False, max_per_G=2500)
            #sup = supergroups(s, G=data[cif], show=True, max_per_G=2500)
        print(sup)
        print("{:6.3f} seconds".format(time()-t0))
