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
from pyxtal.operations import apply_ops, get_best_match
from pyxtal.wyckoff_split import wyckoff_split

ALL_SHIFTS = np.array([[0,0,0],[0,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,0],[1,0,1],[1,1,1]])

def write_poscars_intermediate(H_struc, G_struc, mapping, splitter, N_images=3):
    """
    Write the intermediate POSCARs betwee H and G structure,
    while assuming H is the maximum subgroup of G
    The key is to represent G in the subgroup representation
    Knowing the mapping, one just need to obtain the displacement

    Args:
        H_struc: PyXtal low symmetry structure
        G_struc: PyXtal high symmetry structure
        mapping: atomic mapping between H and G
        splitter: splitter object
        wyc_set: wyc_set transformation
        N_images: number of intermediate structures between H and G

    Return:
        a list of POSCARs
    """
    
    # 1, transform G_struc into a different set
    # 2, represent G in subgroup rep
    # 3, compute the displacement
    # 4, write a seriers of POSCAR by varying disp
    raise NotImplementedError

def write_poscars(H_struc, G_struc, mappings, splitters, wyc_sets, N_images=3):
    """
    Write the intermediate POSCARs betwee H and G structure,
    The key is to continuously change G to subgroup represenations with zero disp
    Finally, call `write_poscars_intermediate`

    Args:
        H_struc: PyXtal low symmetry structure
        G_strucs: a list of PyXtal high symmetry structures
        mapping: a list of atomic mappings 
        splitter: a list of splitter object
        wyc_set: a list of wyc_set transformation
        N_images: number of intermediate structures between H and G

    Return:
        a list of POSCARs
    """
    raise NotImplementedError

"""
def (H_struc, G_struc):

    1, paths between G and H
    2, for path in paths
            for s in s.get_alterntives()
                s.to_subgroup(p0)

        get all subgroup representations for any given wyc_set
        if match:
            break

    s = {}
    for i, p in enumerate(path):
        s[i] 

    for i, p in enumerate(path):
        for s[i] = None

    for s0 in s.get_alts():
    
        

https://www.educative.io/edpresso/how-to-implement-depth-first-search-in-python
#some depth first search

graph = {
    'A' : wyc_ids,
    'B' : ['D', 'E'],
    'C' : ['F'],
    'D' : [],
    'E' : ['F'],
    'F' : []
}

visited = set()
def dfs(visited, graph, node):
    if node not in visited:
        print (node)
        visited.add(node)
        for neighbour in graph[node]:
            dfs(visited, graph, neighbour)


visited = set()
def dfs(visited, )

# Driver Code
dfs(visited, graph, 'A')
"""

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

def new_path(path, paths):
    """
    check if struc is already in the reference solutions
    """
    for ref in paths:
        if path[:len(ref)] == ref:
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
        unique_solutions=[x+[y] for x in unique_solutions \
                for y in combo_list[i] if len(set(sum(x,[])).intersection(y))==0]
    return unique_solutions

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

    # resort the mapping
    mappings = list(itertools.product(*lists))
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
    """
    apply symmetry operation on pos1 when it involves cell change.
    e.g., when the transformation is (a+b, a-b, c),
    trial translation needs to be considered to minimize the
    difference between the transformed pos1 and reference pos2

    Args:
        G:
        rot:
        tran:
        pos:
        wp1:
        op:

    Return:
        tmp:
        diff:
    """

    if np.linalg.det(rot) < 1:
        shifts = ALL_SHIFTS
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


def search_G2(rot, tran, pos1, pos2, cell=None):
    """
    apply symmetry operation on pos1 when it involves cell change.
    e.g., when the transformation is (a+b, a-b, c),
    trial translation needs to be considered to minimize the
    difference between the transformed pos1 and reference pos2

    Args:
        rot:
        tran:
        pos1:
        pos2:
        cell:

    Return:
        pos:
        dist:
    """

    pos1 -= np.round(pos1)
    shifts = ALL_SHIFTS

    dists = []
    for shift in shifts:
        res = np.dot(rot, pos1 + shift + tran.T)
        diff = res - pos2
        diff -= np.round(diff)
        dist = np.linalg.norm(diff)
        dists.append(dist)
        if dist < 1e-1:
            break
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
    if np.all(quadrant == [0,0,0]):
        for i,n in enumerate(coord):
            if n>=0.:
                quadrant[i] = 1
            else:
                quadrant[i] = -1

    # prepare the rotation matrix and translation vector seperately
    G2_holder = [1,1,1]
    G2_op = np.array(G2_op.as_dict()['matrix'])
    rot_G2 = G2_op[:3,:3].T
    tau_G2 = G2_op[:3,3]
    b = coord - tau_G2
    for k in range(3):
        b[k] = b[k]%quadrant[k]

    # eliminate any unused free parameters in G2
    # The goal is to reduce the symmetry operations to be a full rank matrix
    # any free parameter that is not used has its spot deleted from the rotation and translation 
    for i,x in reversed(list(enumerate(rot_G2))):
        if set(x)=={0.}:
            G2_holder[i] = 0
            rot_G2 = np.delete(rot_G2,i,0)
            quadrant = np.delete(quadrant,i)

    # eliminate any leftover empty rows to have fulll rank matrix
    rot_G2=rot_G2.T
    for i, x in reversed(list(enumerate(rot_G2))):
        if set(x) == {0.}:
            rot_G2=np.delete(rot_G2,i,0)
            b=np.delete(b,i)
    while len(rot_G2) != 0 and len(rot_G2) != len(rot_G2[0]):
        rot_G2 = np.delete(rot_G2,len(rot_G2)-1,0)
        b = np.delete(b,len(b)-1)


    # Must come back later and add Schwarz Inequality check to elininate any dependent vectors
    # solves a linear system to find the free parameters
    if set(G2_holder) == {0.}:
        return np.array(G2_holder)

    else:
        try:
            G2_basis_xyz = np.linalg.solve(rot_G2,b)
            for i in range(len(quadrant)):
                G2_basis_xyz[i] = G2_basis_xyz[i]%quadrant[i]
            # print("!ST G2 HOLDER")
            for i in range(G2_holder.count(1)):
                G2_holder[G2_holder.index(1)] = G2_basis_xyz[i]
            # print('second G## holder')
            return np.array(G2_holder)

        except:
            raise RuntimeError('unable to find free parameters using this operation')

class supergroup():
    """
    Class to find the structure with supergroup symmetry

    Args:
        struc: pyxtal structure
        G: target supergroup number
    """
    def __init__(self, struc, G):

        # initilize the necesary parameters
        self.solutions = []
        self.error = True
        self.G = sym.Group(G)
        if self.G.point_group == struc.group.point_group:
            group_type = 'k'
        else:
            group_type = 't'
        self.group_type = group_type

        # list of all alternative wycsets
        strucs = struc.get_alternatives()
        for struc in strucs:
            solutions = self.G.get_splitters_from_structure(struc, group_type)
            if len(solutions) > 0:
                #print(struc)
                self.struc = struc
                self.elements, self.sites = struc._get_elements_and_sites()
                self.solutions = solutions
                self.cell = struc.lattice.matrix
                self.error = False
                break

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
                (id, sols) = sols
                if len(sols) > max_per_G:
                    print("Warning: ignore some solutions: ", len(sols)-max_per_G)
                    sols = sample(sols, max_per_G)
                    #sols=[(['8c'], ['4a', '4b'], ['4b', '8c', '8c'])]
                for sol in sols:
                    #print(sol)
                    mae, disp, mapping, sp = self.get_displacement(id, sol, d_tol*1.1)
                    #print(G.number, sol, mae, disp)
                    if mae < d_tol:
                        valid_solutions.append((sp, mapping, disp, mae))
        return valid_solutions

    def make_supergroup(self, solutions, show_detail=False):
        """
        create supergroup structures based on the list of solutions

        Args:
            solutions: list of tuples (splitter, mapping, disp)
            show_detail (bool): print out the detail

        Returns:
            list of pyxtal structures + mapping relation
        """
        G_strucs = []

        if len(solutions) > 0:

            for solution in solutions:
                (sp, mapping, translation, max_disp) = solution
                lat1 = np.dot(np.linalg.inv(sp.R[:3,:3]).T, self.struc.lattice.matrix)
                lattice = Lattice.from_matrix(lat1, ltype=sp.G.lattice_type)

                details = self.symmetrize(sp, mapping, translation)
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
                        print(self.struc.group.number); print(pos); print(wp)
                        print(">>>>>>>>>>>>>>")
                        raise RuntimeError("cannot assign the right wp")

                G_struc.atom_sites = G_sites
                G_struc.source = 'supergroup {:6.3f}'.format(max_disp)
                G_struc.lattice = lattice
                G_struc.numIons *= int(round(np.abs(np.linalg.det(sp.R[:3,:3]))))
                G_struc._get_formula()
                G_struc.disp = max_disp

                if new_structure(G_struc, G_strucs):
                    # store mapping information
                    G_struc.tag = {'mapping': mapping, 'splitter': sp}
                    G_strucs.append(G_struc)
                    if show_detail:
                        self.print_detail(coords_H1, coords_G2, elements, translation)
                        print(G_struc)
                        import sys; sys.exit()

        return G_strucs

    def get_displacement(self, split_id, solution, d_tol):
        """
        For a given solution, search for the possbile supergroup structure with
        minimum displacement by adusting `translation`.

        Args:
            split_id (int): integer
            solution (list): e.g., [['2d'], ['6h'], ['2c', '6h', '12i']]
            d_tol (float): tolerance

        Returns:
            max_disp: maximum atomic displcement
            translation: overall cell translation
        """
        sites_G = []
        elements = []
        muls = []
        for i, e in enumerate(self.elements):
            sites_G.extend(solution[i])
            elements.extend([e]*len(solution[i]))
            muls.extend([int(sol[:-1]) for sol in solution[i]])

        # resort the sites_G by multiplicity, needed the mask calculation
        ids = np.argsort(np.array(muls))
        elements = [elements[id] for id in ids]
        sites_G = [sites_G[id] for id in ids]

        splitter = wyckoff_split(self.G, split_id, sites_G, self.group_type, elements)
        mappings = find_mapping(self.struc.atom_sites, splitter)

        dists = []
        translations = []
        masks = []
        if len(mappings) > 0:
            mask = self.get_initial_mask(splitter)
            for mapping in mappings:
                dist, translation, mask = self.symmetrize_dist(splitter, mapping, mask, None, d_tol)
                dists.append(dist)
                translations.append(translation)
                masks.append(mask)
            dists = np.array(dists)
            max_disp = np.min(dists)
            id = np.argmin(dists)
            translation = translations[id]
            mask = masks[id]
            if 0.2 < max_disp < d_tol:
                # optimize disp further
                if mask is None or len(mask)<3:
                    def fun(translation, mapping, splitter, mask):
                        return self.symmetrize_dist(splitter, mapping, mask, translation)[0]
                    res = minimize(fun, translations[id], args=(mappings[id], splitter, mask),
                            method='Nelder-Mead', options={'maxiter': 10})
                    if res.fun < max_disp:
                        max_disp = res.fun
                        translation = res.x
            return max_disp, translation, mappings[id], splitter
        else:
            print("bug in findding the mappings", solution)
            print(splitter.G.number, '->', splitter.H.number)
            return 1000, None, None, None

    def get_initial_mask(self, splitter):
        for wp2 in splitter.wp2_lists:
            # if split into 2 sites
            if len(wp2) >= 2: 
                return [0, 1, 2]
            else:
                if wp2[0].get_dof() == 0:
                    return [0, 1, 2]
        return None


    def symmetrize_dist(self, splitter, mapping, mask, translation=None, d_tol=1.2):
        """
        For a given solution, search for the possbile supergroup structure
        based on a given `translation` and `mask`.

        Args:
            splitter: splitter object to specify the relation between G and H
            mapping: list of sites in H, e.g., ['4a', '8b']
            mask: if there is a need to freeze the direction
            translation: an overall shift from H to G, None or 3 vector
            d_tol: the tolerance in angstrom

        Returns:
            distortion
            cell translation
        """

        max_disps = []
        atom_sites_H = self.struc.atom_sites
        rot = splitter.R[:3,:3] # inverse coordinate transformation
        tran = splitter.R[:3,3] # needs to check
        inv_rot = np.linalg.inv(rot)
        cell = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.struc.lattice.matrix)
        ops_G1 = splitter.G[0]

        if mask is not None and translation is not None:
            translation[mask] = 0

        for i, wp1 in enumerate(splitter.wp1_lists):
            if len(splitter.wp2_lists[i]) == 1:
                op_G1 = splitter.G1_orbits[i][0][0]
                ops_H = splitter.H_orbits[i][0]
                base = atom_sites_H[mapping[i][0]].position.copy()

                # choose the best coord1_H
                coord1s_H = apply_ops(base, ops_H)
                ds = []
                for coord1_H in coord1s_H:
                    if translation is not None:
                        coord1_G2 = coord1_H + translation
                    else:
                        coord1_G2 = coord1_H
                    coord1_G1, diff = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G1)
                    ds.append(diff)

                ds = np.array(ds)
                minID = np.argmin(ds)
                coord1_H = coord1s_H[minID]

                if translation is not None:
                    coord1_G2 = coord1_H + translation
                else:
                    coord1_G2 = coord1_H

                tmp, _ = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G1)

                # initial guess on disp
                if translation is None:
                    coord1_G2, dist1 = search_G2(inv_rot, -tran, tmp, coord1_H, None)
                    diff = coord1_G2 - coord1_H
                    diff -= np.round(diff)
                    translation = diff.copy()
                    mask = []
                    for m in range(3):
                        if abs(diff[m])<1e-4:
                            mask.append(m)
                    dist = 0
                else:
                    coord1_H += translation
                    coord1_G2, dist = search_G2(inv_rot, -tran, tmp, coord1_H, self.cell)

                #print("--------", wp1.letter, tmp, coord1_G2, coord1_H+disp, dist)
                if dist < d_tol:
                    #if dist >0: print(dist, coord1_G2, coord1_H+disp)
                    max_disps.append(dist)
                else:
                    #import sys; sys.exit()
                    return 10000, None, None

            elif len(splitter.wp2_lists[i]) == 2:
                # H->G2->G1
                if atom_sites_H[mapping[i][0]].wp.letter == splitter.wp2_lists[i][0].letter:
                    coord1_H = atom_sites_H[mapping[i][0]].position.copy()
                    coord2_H = atom_sites_H[mapping[i][1]].position.copy()
                else:
                    coord2_H = atom_sites_H[mapping[i][0]].position.copy()
                    coord1_H = atom_sites_H[mapping[i][1]].position.copy()

                #print("\n\n\nH", coord1_H, coord2_H)

                coord1_G2 = coord1_H + translation
                coord2_G2 = coord2_H + translation

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
                    #print(splitter.wp1_lists[i].index, coord2_G1, op_G12.as_xyz_string())
                    #coord2_G1 = sym.search_cloest_wp(splitter.G, wp1, op_G12, coord2_G1)
                    #print("G1(symm1)", coord1_G1, coord2_G1)
                    #import sys; sys.exit()

                    # find the best match
                    coords11 = apply_ops(coord1_G1, ops_G1)
                    tmp, dist = get_best_match(coords11, coord2_G1, cell)
                    #tmp = sym.search_cloest_wp(splitter.G, wp1, op_G12, tmp)

                    # G1->G2->H
                    d = coord2_G1 - tmp
                    d -= np.round(d)
                    #print("dist", np.linalg.norm(d), "d", d, tmp, coord2_G1)

                    coord2_G1 -= d/2
                    coord1_G1 += d/2

                    coord1_H += translation
                    coord2_H += translation
                    _, dist1 = search_G2(inv_rot, -tran, coord1_G1, coord1_H, self.cell)
                    _, dist2 = search_G2(inv_rot, -tran, coord2_G1, coord2_H, self.cell)

                    if max([dist1, dist2]) > np.sqrt(2)*d_tol:
                        #import sys; sys.exit()
                        return 10000, None, mask
                    else:
                        max_disps.append(max([dist1, dist2]))

            else:
                # for mergings greater than 2
                # organizes the numbers in the order as the positions in splitter.wp2_lists
                # so that the positions from atoms_sites_H are in the correct assigned position.
                n = len(splitter.wp2_lists[i])
                letters=[atom_sites_H[mapping[i][x]].wp.letter for x in range(n)]
                ordered_letter_index=[]
                for pos in splitter.wp2_lists[i]:
                    indice = letters.index(pos.letter)
                    ordered_letter_index.append(indice)
                    letters[indice] = 0
                ordered_mapping = [mapping[i][x] for x in ordered_letter_index]

                coord_H=[atom_sites_H[ordered_mapping[x]].position.copy() for x in range(n)]

                # Finds the correct quadrant to easily generate all possible_wycs
                # translations when trying to match
                quadrant=np.array(splitter.G2_orbits[i][0][0].as_dict()['matrix'])[:3,3]
                for k in range(3):
                    if quadrant[k] >= 0.:
                        quadrant[k] = 1
                    else:
                        quadrant[k] = -1
                coord_G2 = [x+translation for x in coord_H]
                for j, x in enumerate(coord_G2):
                    for k in range(3):
                        coord_G2[j][k] = coord_G2[j][k]%quadrant[k]

                # uses 1st coordinate and 1st wyckoff position as starting example.
                # Finds the matching G2 operation bcased on nearest G1 search
                dist_list = []
                coord_list = []
                index = []
                corresponding_ops = []
                G2_xyz = []
                for op in splitter.G1_orbits[i][0]:
                    coord,dist = search_G1(splitter.G,rot,tran,coord_G2[0],wp1,op)
                    dist_list.append(dist)
                    coord_list.append(coord)
                dist_list=np.array(dist_list)
                index.append(np.argmin(dist_list))
                corresponding_ops.append(splitter.G2_orbits[i][0][index[0]])

                # Finds the free parameters xyz in the G2 basis for this coordinate
                G2_xyz.append(find_xyz(corresponding_ops[0],coord_G2[0],quadrant))

                # Systematically generates possible G2 positions to match the remaining coordinates
                # Also finds the corresponding G2 free parameters xyz for each coordinate

                for j in range(1, n):
                    possible_coords = [x.operate(G2_xyz[0]) for x in splitter.G2_orbits[i][j]]
                    corresponding_coord, _ = get_best_match(possible_coords,coord_G2[j],cell)
                    index.append([np.all(x==corresponding_coord) for x in possible_coords].index(True))
                    corresponding_ops.append(splitter.G2_orbits[i][j][index[j]])
                    G2_xyz.append(find_xyz(corresponding_ops[j],coord_G2[j],quadrant))

                # Finds the average free parameters between all the coordinates as the best set 
                # of free parameters that all coordinates must match
                final_xyz = np.sum(G2_xyz,axis=0)/n

                dist_list = []
                for j in range(n):
                    G1_coord = splitter.G1_orbits[i][j][index[j]].operate(final_xyz)
                    tmp = coord_H[j] + translation
                    G2_coord,dist = search_G2(inv_rot, -tran, G1_coord, tmp, self.cell)
                    dist_list.append(dist)

                if max(dist_list) > np.sqrt(3)*d_tol:
                    return 10000, None, mask
                else:
                    max_disps.append(max(dist_list))

        return max(max_disps), translation, mask

    def symmetrize(self, splitter, mapping, translation):
        """
        Symmetrize the structure (G) to supergroup symmetry (H)

        Args:
            splitter: splitter object to specify the relation between G and H
            mapping: atomic mapping between H and G
            translation: an overall shift from H to G, None or 3 vector

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

                # choose the best coord1_H
                coord1s_H = apply_ops(base, ops_H)
                ds = []
                for coord1_H in coord1s_H:
                    coord1_G2 = coord1_H + translation
                    coord1_G1, diff = search_G1(splitter.G, rot, tran, coord1_G2, wp1, op_G1)
                    ds.append(diff)

                ds = np.array(ds)
                minID = np.argmin(ds)
                coord1_H = coord1s_H[minID]
    
                # symmetrize coord_G1
                tmp, _ = search_G1(splitter.G, rot, tran, coord1_H+translation, wp1, op_G1)

                coords_G1.append(tmp)
                coord1_G2, _ = search_G2(inv_rot, -tran, tmp, coord1_H+translation, self.cell)
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

                coord1_G2 = coord1_H + translation
                coord2_G2 = coord2_H + translation
                op_G11 = splitter.G1_orbits[i][0][0]
                op_G12 = splitter.G1_orbits[i][1][0]

                if splitter.group_type == 'k':
                    ops_H1 = splitter.H_orbits[i][0]
                    op_G21 = splitter.G2_orbits[i][0][0]
                    ops_G22 = splitter.G2_orbits[i][1]

                    coord11 = coord1_H + translation
                    coord22 = coord2_H + translation

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

                    # find the best match
                    coords11 = apply_ops(coord1_G1, ops_G1)
                    tmp, dist = get_best_match(coords11, coord2_G1, cell)
                    tmp = sym.search_cloest_wp(splitter.G, wp1, op_G12, tmp)

                    # G1->G2->H
                    d = coord2_G1 - tmp
                    d -= np.round(d)
                    coord2_G1 -= d/2
                    coord1_G1 += d/2
                    coords_G1.append(tmp)

                    coord1_H += translation
                    coord2_H += translation
                    coord1_G2, dist1 = search_G2(inv_rot, -tran, coord1_G1, coord1_H, self.cell)
                    coord2_G2, dist2 = search_G2(inv_rot, -tran, coord2_G1, coord2_H, self.cell)

                coords_G2.append(coord1_G2)
                coords_G2.append(coord2_G2)
                coords_H1.append(coord1_H)
                coords_H1.append(coord2_H)

            elements.extend([splitter.elements[i]]*len(splitter.wp2_lists[i]))

        coords_G1 = np.array(coords_G1)
        coords_G2 = np.array(coords_G2)
        coords_H1 = np.array(coords_H1)
        return coords_G1, coords_G2, coords_H1, elements

    def print_detail(self, coords_H1, coords_G1, elements, translation):
        """
        print out the details of tranformation
        """
        print("Valid structure", self.G.number)
        disps = []
        for x, y, ele in zip(coords_H1, coords_G1, elements):
            dis = y - x - translation
            dis -= np.round(dis)
            dis_abs = np.linalg.norm(dis.dot(self.cell))
            output = "{:2s} {:8.4f}{:8.4f}{:8.4f}".format(ele, *x)
            output += " -> {:8.4f}{:8.4f}{:8.4f}".format(*y)
            output += " -> {:8.4f}{:8.4f}{:8.4f} {:8.4f}".format(*dis, dis_abs)
            disps.append(dis_abs)
            print(output)
        print("cell: {:8.4f}{:8.4f}{:8.4f}, disp (A): {:8.4f}".format(*translation, max(disps)))


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

    def __init__(self, struc, G=None, path=None, d_tol=1.0, max_per_G=100, max_layer=5, show=False):

        self.struc0 = struc
        self.show = show
        self.d_tol = d_tol
        self.max_per_G = max_per_G
        self.max_layer = max_layer

        if path is None:
            paths = struc.group.search_supergroup_paths(G, max_layer=max_layer)
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

    def write_poscar(self):
        """
        dump the poscar with the same order of atom mapping
        """
        file1, file2 = 'POSCAR_1', 'POSCAR_2'
        #do something


    def struc_along_path(self, path):
        """
        search for the super group structure along a given path
        
        Args: 
            - path: [59, 71, 139]

        Returns:
            - strucs: list of structures along the path
            - working_path:
            - valid: True or False
        """
        strucs = []
        G_strucs = [self.struc0]
        working_path = []
        for G in path:
            working_path.append(G)
            for G_struc in G_strucs:
                my = supergroup(G_struc, G)
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


if __name__ == "__main__":

    from pyxtal import pyxtal
    from time import time
    data = {
            #"PVO": [12, 166],
            #"PPO": [12],
            #"BTO": [123, 221],
            #"lt_cristobalite": [98, 210, 227],
            "MPWO": [59, 71, 139, 225],
            #"BTO-Amm2": [65, 123, 221],
            #"NaSb3F10": [186, 194],
            #"GeF2": 62,
            #"NiS-Cm": 160,
            #"lt_quartz": 180,
            #"BTO-Amm2": 221,
            #"BTO": 221,
            #"lt_cristobalite": 227,
            #"NaSb3F10": 194,
            #"MPWO": 225,
            #"NbO2": 141,
           }
    cif_path = "pyxtal/database/cifs/"

    for cif in data.keys():
        t0 = time()
        print("===============", cif, "===============")
        s = pyxtal()
        s.from_seed(cif_path+cif+'.cif')
        if isinstance(data[cif], list):
            sup = supergroups(s, path=data[cif], show=True, max_per_G=2500)
        else:
            sup = supergroups(s, G=data[cif], show=False, max_per_G=2500)
        print(sup)
        print("{:6.3f} seconds".format(time()-t0))
        #for i, struc in enumerate(sup.strucs):
        #    struc.to_file(str(i)+'-G'+str(struc.group.number)+'.cif')
