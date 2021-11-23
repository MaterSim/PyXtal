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

def find_mapping_per_element(sites1, sites2,max_num=720):
    """
    search for all mappings for a given splitter

    Args:
        sites1 (list): e.g., l layer ['4a', '8b', '4c']
        sites2 (list): e.g., 2 layers [['4a'], ['8b', '4c']]
        max_num (int): maximum number of atomic mapping

    Returns:
        unique solutions: e.g. 3 layers: [[[0], [1,2]]]
    """

    unique_letters=list(set(sites1))
    site1_letter_indices=[]
    for letter in unique_letters:
        site1_letter_indices.append([i for i, x in enumerate(sites1) if x==letter])
    site2_letter_bins=[]
    for lbin in sites2:
        site2_letter_bins.append([unique_letters.index(x) for x in lbin])

    combo_list=[]
    for s in site2_letter_bins:
        ls=list(set(s))
        rs=[s.count(r) for r in ls]
        p=[]
        for i, l in enumerate(ls):
            combo=itertools.combinations(site1_letter_indices[l],rs[i])
            combo=[list(x) for x in combo]
            p.append(deepcopy(combo))
        pr=p[0]
        for i in range(1,len(p)):
            pr=itertools.product(pr,p[i])
            pr=[sum(list(x),[]) for x in pr]
        combo_list.append(pr)
    unique_solutions=[[x] for x in combo_list[0]]
    for i in range(1,len(combo_list)):
        unique_solutions=[x+[y] for x in unique_solutions for y in combo_list[i] if len(set(sum(x,[])).intersection(y))==0]
    return unique_solutions


    #depreciated  mapping function
#     print('sites1=',sites1)
#     print('sites2=',sites2)
#     unique_solutions = []
#     solution_template = [[None]*len(site2) for site2 in sites2]
#     assigned_ids = []

#     # first identify the unique assignment
#     for i, site2 in enumerate(sites2):
#         wp_letters = set(site2)
#         if len(wp_letters) == len(site2): #site2: ['a', 'b'] or ['a']
#             for j, s2 in enumerate(site2):
#                 ids = [id for id, s1 in enumerate(sites1) if s1==s2]
#                 if len(ids) == 1:
#                     solution_template[i][j] = ids[0]
#                     assigned_ids.append(ids[0])
#         elif len(wp_letters)==1: #site2: ['a','a'] or ['a','a','a']
#             ids = [id for id, s1 in enumerate(sites1) if s1==list(wp_letters)[0]]
#             if len(ids) == len(site2):
#                 solution_template[i] = ids
#                 assigned_ids.extend(ids)
#         elif len(wp_letters)==2: #site2: ['a','a','b']
#             count = 0
#             for j in range(2):
#                 ids = [id for id, s1 in enumerate(sites1) if s1==list(wp_letters)[j]]
#                 if len(ids) == site2.count(list(wp_letters)[j]):
#                     solution_template[i][count:count+len(ids)] = ids
#                     assigned_ids.extend(ids)
#                     count += len(ids)
#             #raise NotImplementedError("unsupported:", site2)
#     #print(assigned_ids)

#     ids = [id for id, site in enumerate(sites1) if id not in assigned_ids]

#     all_permutations = list(itertools.permutations(ids))
#     if len(all_permutations) > max_num:
#         print("Warning: ignore some mapping: ", str(len(all_permutations)-max_num))
#         print(solution_template)
#         all_permutations = sample(all_permutations, max_num)

#     print('allpermuataions')
#     for x in all_permutations:
#         print(x)
#     #print(solution_template)
#     for perm in all_permutations:
#         solution = deepcopy(solution_template)
#         perm = list(perm)
#         valid = True
#         count = 0
#         for i, sol in enumerate(solution):
#             if None in sol:
#                 for j, s2 in enumerate(sites2[i]):
#                     if sol[j] is None:
#                         if s2 == sites1[perm[count]]:
#                             solution[i][j] = deepcopy(perm[count])
#                             count += 1
#                         else:
#                             valid = False
#                             break
#                 if not valid:
#                     break
#         if valid and new_solution(solution, unique_solutions):
#             unique_solutions.append(solution)
#     print('unique_solutions')
#     for x in unique_solutions:
#         print(x)
#     return unique_solutions

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
    eles = set([site.specie for site in atom_sites])

    # loop over the mapping for each element
    # then propogate the possible mapping via itertools.product
    lists = []
    for ele in eles:
        # ids of atom sites
        site_ids = [id for id, site in enumerate(atom_sites) if site.specie==ele]
        # ids to be assigned
        wp2_ids = [id for id, e in enumerate(splitter.elements) if e==ele]

        letters1 = [atom_sites[id].wp.letter for id in site_ids]
        letters2 = []
        for id in wp2_ids:
            wp2 = splitter.wp2_lists[id]
            letters2.append([wp.letter for wp in wp2])
        #print(ele, letters1, letters2)
        res = find_mapping_per_element(letters1, letters2, max_num=720)
        lists.append(res)
    mappings = list(itertools.product(*lists))
    # resort the mapping
    ordered_mappings = []
    for mapping in mappings:
        ordered_mapping = [None]*len(splitter.wp2_lists)
        for i, ele in enumerate(eles):
            site_ids = [id for id, site in enumerate(atom_sites) if site.specie==ele]
            count = 0
            for j, wp2 in enumerate(splitter.wp2_lists):
                if splitter.elements[j] == ele:
                    ordered_mapping[j] = [site_ids[m] for m in mapping[i][count]]
                    count += 1
        #print("res", ordered_mapping)
        ordered_mappings.append(ordered_mapping)

    #if len(ordered_mappings)==0: import sys; sys.exit()
    return ordered_mappings

def search_G1(G, rot, tran, pos, wp1, op):
    if np.linalg.det(rot) < 1:
        shifts = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])
    else:
        shifts = np.array([[0,0,0]])

    diffs = []
    coords = []
    for shift in shifts:
        res = np.dot(rot, pos + shift) + tran.T
        tmp = sym.search_cloest_wp(G, wp1, op, res)
        diff = res - tmp
        diff -= np.round(diff)
        dist = np.linalg.norm(diff)
        diffs.append(dist)
        coords.append(tmp)
        if dist < 1e-1:
            break

    diffs = np.array(diffs)
    minID = np.argmin(diffs)
    tmp = coords[minID]
    tmp -= np.round(tmp)
    return tmp, np.min(diffs)


def search_G2(rot, tran, pos1, pos2, cell=None, ortho=True):
    """
    apply symmetry operation on pos1 when it involves cell change.
    e.g., when the transformation is (a+b, a-b, c),
    trial translation needs to be considered to minimize the
    difference between the transformed pos1 and reference pos2
    """
    pos1 -= np.round(pos1)
    if np.linalg.det(rot) < 1:
        shifts = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])
    elif not ortho:
        shifts = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])
    else:
        shifts = np.array([[0,0,0]])
    shifts = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])

    dists = []
    for shift in shifts:
        res = np.dot(rot, pos1 + shift + tran.T)
        diff = res - pos2
        diff -= np.round(diff)
        dist = np.linalg.norm(diff)
        dists.append(dist)
        if dist < 1e-1:
            break
    #print("TTTTTTTTT", dists, np.linalg.det(rot))
    dists = np.array(dists)
    dist = np.min(dists)
    shift = shifts[np.argmin(dists)]
    pos = np.dot(rot, pos1 + shift + tran.T)

    diff = pos - pos2
    diff -= np.round(diff)

    if cell is not None:
        diff = np.dot(diff, cell)

    dist = np.linalg.norm(diff)

    return pos, dist

def find_xyz(G2_op, coord, quadrant=[0,0,0]):
    """
    Finds the x,y,z free parameter values for positions in the G_2 basis.

    Args:
        G2_op: a symmetry operation in G2
        coord: the coordinate that matches G2_op
        quadrant: a 3 item list (ex:[1,1,-1]) that contains information
                  on the orientation of the molecule

    Returns:
        G2_holder: The corresponding x,y,z parameters written in the G2 basis
    """
    if np.all(quadrant==[0,0,0]):
        for i,n in enumerate(coord):
            if n>=0.:
                quadrant[i]=1
            else:
                quadrant[i]=-1

    #prepare the rotation matrix and translation vector seperately
    G2_holder=[1,1,1]
    G2_op=np.array(G2_op.as_dict()['matrix'])
    rot_G2=G2_op[:3,:3].T
    tau_G2=G2_op[:3,3]
    b=coord-tau_G2
    for k in range(3):
        b[k]=b[k]%quadrant[k]

    #eliminate any unused free parameters in G2
    #The goal is to reduce the symmetry operations to be a full rank matrix
    #any free parameter that is not used has its spot deleted from the rotation matrix and translation vector
    for i,x in reversed(list(enumerate(rot_G2))):
        if set(x)=={0.}:
            G2_holder[i]=0
            rot_G2=np.delete(rot_G2,i,0)
            quadrant=np.delete(quadrant,i)

    #eliminate any leftover empty rows to have fulll rank matrix
    rot_G2=rot_G2.T
    for i,x in reversed(list(enumerate(rot_G2))):
        if set(x)=={0.}:
            rot_G2=np.delete(rot_G2,i,0)
            b=np.delete(b,i)
    while len(rot_G2)!=0 and len(rot_G2)!=len(rot_G2[0]):
        rot_G2=np.delete(rot_G2,len(rot_G2)-1,0)
        b=np.delete(b,len(b)-1)


    #Must come back later and add Schwarz Inequality check to elininate any dependent vectors
    #solves a linear system to find the free parameters
    if set(G2_holder)=={0.}:
        return np.array(G2_holder)

    else:
        try:
            G2_basis_xyz=np.linalg.solve(rot_G2,b)
            for i in range(len(quadrant)):
                G2_basis_xyz[i]=G2_basis_xyz[i]%quadrant[i]
            # print("!ST G2 HOLDER")
            for i in range(G2_holder.count(1)):
                G2_holder[G2_holder.index(1)]=G2_basis_xyz[i]
            # print('second G## holder')
            return np.array(G2_holder)

        except:
            raise RuntimeError('unable to find free parameters using this operation')


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
    #G = sym.Group(G)
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
    l2 = Lattice.from_matrix(matrix, ltype=G.lattice_type)
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
    #G = sym.Group(G)

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
    layers[0]={'groups':[G], 'subgroups':[]}
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
                paths=[[g]]
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
            #print(status)
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
            #print(status)


    def __str__(self):
        s = "\nTransition to super group: "
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
            #print(G, H)
            #if G != H:
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
        self.wyc_supergroups = struc.group.get_min_supergroup(group_type, G)

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
                G = sym.Group(self.wyc_supergroups['supergroup'][idx])
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
            #print("No compatible solution exists")
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
                #sols=[(['8c'], ['4a', '4b'], ['4b', '8c', '8c'])]
                for sol in sols:
                    #print(sol)
                    mae, disp, mapping, sp = self.get_displacement(G, id, sol, d_tol*1.1)
                    #print(G.number, sol, mae, disp)
                    if mae < d_tol:
                        valid_solutions.append((sp, mapping, disp, mae))
        return valid_solutions

    def make_supergroup(self, solutions, once=False, show_detail=False):
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
                        print(self.struc.group.number)
                        print(pos)
                        print(wp)
                        print(">>>>>>>>>>>>>>")
                        raise RuntimeError("cannot assign the right wp")

                G_struc.atom_sites = G_sites
                G_struc.source = 'supergroup {:6.3f}'.format(mae)
                G_struc.lattice = lattice
                G_struc.numIons *= int(round(np.abs(np.linalg.det(sp.R[:3,:3]))))
                G_struc._get_formula()
                G_struc.disp = mae

                if new_structure(G_struc, G_strucs):
                    G_strucs.append(G_struc)
                    if show_detail:
                        G = sp.G.number
                        self.print_detail(G, coords_H1, coords_G2, elements, disp)
                        print(G_struc)

        return G_strucs

    def get_displacement(self, G, split_id, solution, d_tol):
        """
        For a given solution, search for the possbile supergroup structure

        Args:
            G: group object
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
        if splitter.H.number <= 15 or 143<= splitter.H.number <= 194:
            ortho = False
        else:
            ortho = True
        cell = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.struc.lattice.matrix)
        max_disps = []
        atom_sites_H = self.struc.atom_sites
        rot = splitter.R[:3,:3] # inverse coordinate transformation
        tran = splitter.R[:3,3] # needs to check
        inv_rot = np.linalg.inv(rot)
        cell = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.struc.lattice.matrix)
        ops_G1  = splitter.G[0]
        #print(mapping)
        # if there involves wp transformation between frozen points, disp has to be zero
        for wp2 in splitter.wp2_lists:
            frozen = True
            for wp in wp2:
                if np.linalg.matrix_rank(wp[0].rotation_matrix) > 0:
                    frozen = False

            if frozen:
                disp = np.zeros(3)
                mask = [0, 1, 2]
                break

        if mask is not None and disp is not None:
            disp[mask] = 0

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
                    coord1_G1, diff = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G1)
                    ds.append(diff)

                ds = np.array(ds)
                minID = np.argmin(ds)
                coord1_H = coord1s_H[minID]

                if disp is not None:
                    coord1_G2 = coord1_H + disp
                else:
                    coord1_G2 = coord1_H

                tmp, _ = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G1)

                # initial guess on disp
                if disp is None:
                    coord1_G2, dist1 = search_G2(inv_rot, -tran, tmp, coord1_H, None, ortho)
                    diff = coord1_G2 - coord1_H
                    diff -= np.round(diff)
                    disp = diff.copy()
                    mask = []
                    for m in range(3):
                        if abs(diff[m])<1e-4:
                            mask.append(m)
                    dist = 0
                else:
                    coord1_G2, dist = search_G2(inv_rot, -tran, tmp, coord1_H+disp, self.cell, ortho)

                #print("--------", wp1.letter, tmp, coord1_G2, coord1_H+disp, dist)
                if dist < d_tol:
                    #if dist >0: print(dist, coord1_G2, coord1_H+disp)
                    max_disps.append(dist)
                else:
                    #import sys; sys.exit()
                    return 10000, None, None

            elif len(splitter.wp2_lists[i]) == 2:
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
                        #print("disp:", disp)
                        #print("G21:", coords11)
                        #print("G22:", coord2_G2)
                        #print("start: ", coord1_H, coord2_H)
                        #res = tmp - coord2_G2
                        #res -= np.round(res)
                        #print("kkkk", wp1.letter, dist, tmp, coord2_G2, "diff", res)
                        return 10000, None, mask
                    else:
                        d = coord2_G2 - tmp
                        d -= np.round(d)
                        max_disps.append(np.linalg.norm(np.dot(d/2, self.cell)))

                else:
                    #print("disp", disp)
                    #print("G2", coord1_G2, coord2_G2)
                    op_G11 = splitter.G1_orbits[i][0][0]
                    op_G12 = splitter.G1_orbits[i][1][0]
                    coord1_G1, _ = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G11)
                    coord2_G1, _ = search_G1(splitter.G, rot, tran, coord2_G2, wp1, op_G12)

                    #print("G1", coord1_G1, coord2_G1, op_G12.as_xyz_string())
                    #print(splitter.G.number, splitter.wp1_lists[i].index, coord2_G1, op_G12.as_xyz_string())
                    #coord2_G1 = sym.search_cloest_wp(splitter.G, wp1, op_G12, coord2_G1)
                    #print("G1(symm1)", coord1_G1, coord2_G1)
                    #import sys; sys.exit()
                    #find the best match
                    coords11 = apply_ops(coord1_G1, ops_G1)
                    tmp, dist = get_best_match(coords11, coord2_G1, cell)
                    #tmp = sym.search_cloest_wp(splitter.G, wp1, op_G12, tmp)

                    # G1->G2->H
                    d = coord2_G1 - tmp
                    d -= np.round(d)
                    #print("dist", np.linalg.norm(d), "d", d, tmp, coord2_G1)

                    coord2_G1 -= d/2
                    coord1_G1 += d/2

                    coord1_G2, dist1 = search_G2(inv_rot, -tran, coord1_G1, coord1_H+disp, self.cell, ortho)
                    coord2_G2, dist2 = search_G2(inv_rot, -tran, coord2_G1, coord2_H+disp, self.cell, ortho)


                    if max([dist1, dist2]) > np.sqrt(2)*d_tol:
                        #import sys; sys.exit()
                        return 10000, None, mask
                    else:
                        max_disps.append(max([dist1, dist2]))

            #for mergings greater than 2
            else:
                n=len(splitter.wp2_lists[i])
                # assume zero shift, needs to check
                if disp is None:
                    disp = np.zeros(3)
                    mask = [0, 1, 2]


                #organizes the numbers from the mapping to be in same order as the positions in splitter.wp2_lists
                #so that the positions from atoms_sites_H are in the correct assigned position.
                letters=[atom_sites_H[mapping[i][x]].wp.letter for x in range(n)]
                ordered_letter_index=[]
                for pos in splitter.wp2_lists[i]:
                    indice=letters.index(pos.letter)
                    ordered_letter_index.append(indice)
                    letters[indice]=0
                ordered_mapping=[mapping[i][x] for x in ordered_letter_index]

                coord_H=[atom_sites_H[ordered_mapping[x]].position.copy() for x in range(n)]

                #Finds the correct quadrant that the coordinates lie in to easily generate all possible_wycs
                #translations when trying to match
                quadrant=np.array(splitter.G2_orbits[i][0][0].as_dict()['matrix'])[:3,3]
                for k in range(3):
                    if quadrant[k]>=0.:
                        quadrant[k]=1
                    else:
                        quadrant[k]=-1
                coord_G2=[x+disp for x in coord_H]
                for j,x in enumerate(coord_G2):
                    for k in range(3):
                        coord_G2[j][k]=coord_G2[j][k]%quadrant[k]

                #uses 1st coordinate and 1st wyckoff position as starting example.
                #Finds the matching G2 operation bcased on nearest G1 search
                dist_list=[]
                coord_list=[]
                index=[]
                corresponding_ops=[]
                G2_xyz=[]
                for op in splitter.G1_orbits[i][0]:
                    coord,dist=search_G1(splitter.G,rot,tran,coord_G2[0],wp1,op)
                    dist_list.append(dist)
                    coord_list.append(coord)
                dist_list=np.array(dist_list)
                index.append(np.argmin(dist_list))
                corresponding_ops.append(splitter.G2_orbits[i][0][index[0]])

                #Finds the free parameters xyz in the G2 basis for this coordinate
                G2_xyz.append(find_xyz(corresponding_ops[0],coord_G2[0],quadrant))

                #systematically generates possible G2 positions to match the remainding coordinates
                #with. Also finds the corresponding G2 free parameters xyz for each coordinate

                for j in range(1,n):
                    possible_coords=[x.operate(G2_xyz[0]) for x in splitter.G2_orbits[i][j]]
                    corresponding_coord,_=get_best_match(possible_coords,coord_G2[j],cell)
                    index.append([np.all(x==corresponding_coord) for x in possible_coords].index(True))
                    corresponding_ops.append(splitter.G2_orbits[i][j][index[j]])
                    G2_xyz.append(find_xyz(corresponding_ops[j],coord_G2[j],quadrant))

                #Finds the average free parameters between all the coordinates as the best set of free
                #parameters that all coordinates must match
                final_xyz=np.sum(G2_xyz,axis=0)/n


                dist_list=[]
                for j in range(n):
                    G1_coord=splitter.G1_orbits[i][j][index[j]].operate(final_xyz)
                    G2_coord,dist=search_G2(inv_rot, -tran, G1_coord,coord_H[j]+disp,self.cell)
                    dist_list.append(dist)



                if max(dist_list) > np.sqrt(3)*d_tol:
                    return 10000, None, mask
                else:
                    max_disps.append(max(dist_list))

                #depreciated merging for specifically k_type 3->1 merge
                # if splitter.group_type == 'k':
                #     #This functionality uses the fact that, in a k type transitinos, atoms must be displaced such that
                #     #there is just 1 set of rotation matrices with n different trnaslation vectors
                #     #where n is the index of the splitting
                #     ops_H1 = splitter.H_orbits[i][0]
                #     op_G21 = splitter.G2_orbits[i][0][0]
                #
                #     ops_G22 = splitter.G2_orbits[i][1]
                #
                #     #tries to find the two translation vectors (along with [0,0,0]) specific to this
                #     #splitting.
                #     for op_G22 in ops_G22:
                #         diff = (op_G22.rotation_matrix - op_G21.rotation_matrix).flatten()
                #         if np.sum(diff**2) < 1e-3:
                #             trans2 = op_G22.translation_vector - op_G21.translation_vector
                #             trans3 = 2*deepcopy(trans2)
                #             for k in range(3):
                #                 trans2[k]=trans2[k]%quadrant[k]
                #                 trans3[k]=trans3[k]%quadrant[k]
                #             translations=[trans2,trans3]
                #             translations.append(np.array([0,0,0]))
                #             translations=np.array(translations)
                #             break
                #
                #
                #     #Generates all rotation matrices by applying every operation of the wyckoff position of the
                #     #1st coordinate to the 1st coordinate.
                #     coords11 = apply_ops(coord1_G2,ops_H1)
                #     for j in range(len(coords11)):
                #         for k in range(3):
                #             coords11[j][k]=coords11[j][k]%quadrant[k]
                #
                #     #Generates all the new possible positions using the translations
                #     possible_coord_G2=[]
                #     for coordinate in coords11:
                #         for translation in translations:
                #             t=coordinate+translation
                #             for j in range(3):
                #                 t[j]=t[j]%quadrant[j]
                #             possible_coord_G2.append(t)
                #
                #
                #     tmp12,dist12=get_best_match(possible_coord_G2,coord2_G2,self.cell)
                #     tmp13,dist13=get_best_match(possible_coord_G2,coord3_G2,self.cell)
                #     if max([dist12,dist13]) > np.sqrt(3)*d_tol:
                #         return 10000, None, mask
                #     else:
                #         max_disps.append(max([dist12,dist13]))

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
                    coord1_G1, diff = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G1)
                    ds.append(diff)

                ds = np.array(ds)
                minID = np.argmin(ds)
                coord1_H = coord1s_H[minID]

                # symmetrize coord_G1
                tmp, _ = search_G1(splitter.G, rot, tran, coord1_H+disp, wp1, op_G1)

                coords_G1.append(tmp)
                coord1_G2, _ = search_G2(inv_rot, -tran, tmp, coord1_H+disp, self.cell)
                #print("G2", wp1.letter, coord1_G2, tmp)
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
                op_G11 = splitter.G1_orbits[i][0][0]
                op_G12 = splitter.G1_orbits[i][1][0]

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
                    #print("G1:", tmp)
                    d = coord2_G2 - tmp
                    d -= np.round(d)
                    coord2_G2 -= d/2
                    coord1_G2 += d/2
                    coord1_G1, _ = search_G1(splitter.G, rot, tran, tmp, wp1, op_G11)
                    #print(coord1_G2, coord2_G2, np.dot(rot, tmp).T + tran.T)
                    coords_G1.append(coord1_G1)

                else:
                    # H->G2->G1
                    coord1_G1, _ = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G11)
                    coord2_G1, _ = search_G1(splitter.G, rot, tran, coord2_G2, wp1, op_G12)

                    #find the best match
                    coords11 = apply_ops(coord1_G1, ops_G1)
                    tmp, dist = get_best_match(coords11, coord2_G1, cell)
                    tmp = sym.search_cloest_wp(splitter.G, wp1, op_G12, tmp)

                    # G1->G2->H
                    d = coord2_G1 - tmp
                    d -= np.round(d)
                    coord2_G1 -= d/2
                    coord1_G1 += d/2
                    coords_G1.append(tmp)

                    coord1_G2, dist1 = search_G2(inv_rot, -tran, coord1_G1, coord1_H+disp, self.cell)
                    coord2_G2, dist2 = search_G2(inv_rot, -tran, coord2_G1, coord2_H+disp, self.cell)

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
                #"GeF2": 62,
                #"NiS-Cm": 160,
                #"lt_quartz": 180,
                #"BTO-Amm2": 221,
                #"BTO": 221,
                #"lt_cristobalite": 227,
                #"NaSb3F10": 194,
                "MPWO": 225,
                #"NbO2": 141,
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
            for i, struc in enumerate(sup.strucs):
                struc.to_file(str(i)+'-G'+str(struc.group.number)+'.cif')
