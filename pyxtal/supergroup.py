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

def find_mapping_per_element(sites1, sites2, max_num=720):
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
    For a given `pos`, apply the transformation due to (rot+tran),
    search for the cloest wp which satisfies the symmetry of op
    within the trial translation needs to be considered to minimize the
    difference between the transformed pos1 and reference pos2

    Args:
        G: the target space group object
        rot: rotation, 3*3 matrix
        tran: translation, 1*3 vector
        pos: starting position
        wp1: wyckoff symmetry
        op: one symmetry operation

    Return:
        the cloest position and the distance
    """

    #if np.linalg.det(rot) < 1:
    shifts = ALL_SHIFTS

    diffs = []
    coords = []
    # loop over all nearby translations
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
    # choose the one returns minimum difference
    diffs = np.array(diffs)
    minID = np.argmin(diffs)
    tmp = coords[minID]
    tmp -= np.round(tmp)
    return tmp, np.min(diffs)


def search_G2(rot, tran, pos1, pos2, cell=None):
    """

    trial translation needs to be considered to minimize the
    difference between the transformed pos1 and reference pos2

    Args:
        rot:
        tran:
        pos1: position in G1
        pos2: reference position in G2
        cell: 3*3 matrix

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
        G2_holder: x,y,z parameters written in the G2 basis
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

    def search_supergroup(self, d_tol=1.0, max_per_G=2500, max_solutions=None):
        """
        search for valid supergroup transition

        Args:
            d_tol (float): tolerance for atomic displacement
            max_per_G (int): maximum number of possible solution for each G
        Returns:
            solutions: list of solutions with small displacements
        """
        solutions = []
        if len(self.solutions) > 0:
            # extract the valid
            for sols in self.solutions:
                (id, sols) = sols
                if len(sols) > max_per_G:
                    print("Warning: ignore some solutions: ", len(sols)-max_per_G)
                    sols = sample(sols, max_per_G)
                    #sols=[(['8c'], ['4a', '4b'], ['4b', '8c', '8c'])]
                for i, sol in enumerate(sols):
                    #print("+++++++++++++++++++++", i, sol)
                    max_disp, trans, mapping, sp = self.calc_disps(id, sol, d_tol*1.1)
                    #print(G.number, sol, mae, disp)
                    if max_disp < d_tol:
                        solutions.append((sp, mapping, trans, max_disp))
                        if max_solutions is not None and len(solutions) == max_solutions:
                            break
        return solutions

    def make_supergroup(self, solutions, show_detail=False):
        """
        create unique supergroup structures from a list of solutions

        Args:
            solutions: list of tuples (splitter, mapping, translation, disp)
            show_detail (bool): print out the detail

        Returns:
            list of pyxtal structures
        """
        G_strucs = []
        for solution in solutions:
            (sp, mapping, translation, max_disp) = solution

            details = self.symmetrize(sp, mapping, translation)
            coords_G1, coords_G2, coords_H1, elements = details
            G_struc = self._make_pyxtal(sp, coords_G1)
            if new_structure(G_struc, G_strucs):
                if show_detail:
                    self.print_detail(solution, coords_H1, coords_G2, elements)
                G_struc.source = 'supergroup {:6.3f}'.format(max_disp)
                G_struc.disp = max_disp
                G_strucs.append(G_struc)
        return G_strucs


    def calc_disps(self, split_id, solution, d_tol):
        """
        For a given solution, compute the minimum displacement by adusting `translation`.

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
                #Maybe collect only valuable
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
            #if len(wp2) >= 2: 
            #    return [0, 1, 2]
            #else:
            for wp in wp2:
                if wp.get_dof() == 0:
                    return [0, 1, 2]
        #print(splitter); import sys; sys.exit()
        return None

    def get_coord_H(self, splitter, id, atom_sites_H, mapping):
        # number of split sites for a given WP
        n = len(splitter.wp2_lists[id])
        if n > 1:
            letters = [atom_sites_H[mapping[id][x]].wp.letter for x in range(n)]
            letters_wp = [wp.letter for wp in splitter.wp2_lists[id]]
            seq = [] #list(map(lambda x: letters.index(x), letters_wp))
            for l in letters_wp:
                index = letters.index(l)
                seq.append(index)
                letters[index] = 0

            ordered_mapping = [mapping[id][x] for x in seq]
            #print(letters, mapping[id], '->', letters_wp, ordered_mapping)
        else:
            ordered_mapping = mapping[id]
        coord_H = [atom_sites_H[ordered_mapping[x]].position.copy() for x in range(n)]
        return np.array(coord_H), ordered_mapping

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
        if mask is not None: 
            if translation is not None:
                translation[mask] = 0

        for i in range(len(splitter.wp1_lists)):

            n = len(splitter.wp2_lists[i])
            coord_H, _ = self.get_coord_H(splitter, i, self.struc.atom_sites, mapping)

            if n == 1:
                res = self.symmetrize_site_single(splitter, i, coord_H[0], translation)
                (dist, _tran, _mask) = res
                if translation is None:
                    translation = _tran
                    mask = _mask
            elif n == 2:
                if splitter.group_type == 'k':
                    dist = self.symmetrize_site_double_k(splitter, i, coord_H, translation)
                else:
                    dist = self.symmetrize_site_double_t(splitter, i, coord_H, translation)
            else:
                dist = self.symmetrize_site_multi(splitter, i, coord_H, translation)
    
            #strs = self.print_wp(splitter, i); print(strs, dist)

            if i == 0 and translation is None:
                translation = np.zeros(3)

            if dist < d_tol:
                max_disps.append(dist)
            else:
                return 10000, None, mask

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
        coords_G1 = [] # position in G
        coords_G2 = [] # position in G on the subgroup bais
        coords_H = [] # position in H
        elements = []

        # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
        for i, wp1 in enumerate(splitter.wp1_lists):

            n = len(splitter.wp2_lists[i])
            coord_H, _ = self.get_coord_H(splitter, i, self.struc.atom_sites, mapping)

            if n == 1:
                res = self.symmetrize_site_single(splitter, i, coord_H[0], translation, 0)       
            elif n==2:
                if splitter.group_type == 'k':
                    res = self.symmetrize_site_double_k(splitter, i, coord_H, translation, 0)
                else:
                    res = self.symmetrize_site_double_t(splitter, i, coord_H, translation, 0)
            else:
                res = self.symmetrize_site_multi(splitter, i, coord_H, translation, 0)


            coord_G1, coord_G2, coord_H = res
            coords_G1.append(coord_G1)
            coords_G2.extend(coord_G2)
            coords_H.extend(coord_H)
            elements.extend([splitter.elements[i]]*n)
        
            #self.print_wp(splitter, i)
            #print(coord_G1)
            #print(coord_G2)
 
        coords_G1 = np.array(coords_G1)
        coords_G2 = np.array(coords_G2)
        coords_H = np.array(coords_H)
        
        return coords_G1, coords_G2, coords_H, elements

    def print_wp(self, splitter, id):
        wp1 = splitter.wp1_lists[id]
        l = str(wp1.multiplicity) + wp1.letter + '->'
        for wp in splitter.wp2_lists[id]:
            l += str(wp.multiplicity) + wp.letter + ','
        strs = '{:2s}{:s} ID-{:d} {:s}'.format(splitter.elements[id], splitter.group_type, id, l)
        return strs

    def symmetrize_site_single(self, splitter, id, base, translation, run_type=1):
        """
        symmetrize one WP to another with higher symmetry 

        Args:
            splitter: splitter object
            id: index of splitter
            base: atomic position of site in H
            translation: 1*3 translation vector
            run_type: return distance or coordinates
        """
        # Some necessary items
        mask = []
        op_G1 = splitter.G1_orbits[id][0][0]
        ops_H = splitter.H_orbits[id][0]
        wp1 = splitter.wp1_lists[id]
        rot = splitter.R[:3,:3]
        tran = splitter.R[:3,3]
        inv_rot = np.linalg.inv(rot)

        # choose the best coord1_H
        coords_H = apply_ops(base, ops_H)

        ds = []
        for coord_H in coords_H:
            if translation is not None:
                coord_G2 = coord_H + translation
            else:
                coord_G2 = coord_H
            coord_G1, diff = search_G1(splitter.G, rot, tran, coord_G2, wp1, op_G1)
            ds.append(diff)

        ds = np.array(ds)
        minID = np.argmin(ds)
        coord_H = coords_H[minID]

        if run_type == 1:
            if translation is not None:
                coord_G2 = coord_H + translation
            else:
                coord_G2 = coord_H

            tmp, _ = search_G1(splitter.G, rot, tran, coord_G2, wp1, op_G1)

            # initial guess on disp
            if translation is None:
                coord_G2, dist1 = search_G2(inv_rot, -tran, tmp, coord_H, None)
                diff = coord_G2 - coord_H
                diff -= np.round(diff)
                translation = diff.copy()
                for m in range(3):
                    if abs(diff[m])<1e-4:
                        mask.append(m)
                dist = 0
            coord_G2, dist = search_G2(inv_rot, -tran, tmp, coord_H+translation, self.cell)

            return dist, translation, mask
        else:
            # symmetrize coord_G1
            tmp, _ = search_G1(splitter.G, rot, tran, coord_H+translation, wp1, op_G1)
            coord_G2, dist = search_G2(inv_rot, -tran, tmp, coord_H+translation, self.cell)
            #print('XXXXXXXXX', coord_H+translation, tmp, coord_G2, dist)
            return tmp, [coord_G2], [coord_H]

    def symmetrize_site_double_k(self, splitter, id, coord_H, translation, run_type=1):
        """
        symmetrize two WPs (wp_h1, wp_h2) to another wp_G with higher symmetry 

        Args:
            splitter: splitter object
            id: index of splitter
            coord_H: 2*3 coordinates
            translation: 1*3 transaltion vector 
            run_type: return distance or coordinates
        """
        # For k-type splitting, restore the translation symmetry:
        # e.g. (0.5, 0.5, 0.5), (0.5, 0, 0), .etc
        # then find the best_match between coord1 and coord2,

        if translation is None: translation = np.zeros(3)
        rot = splitter.R[:3,:3]
        tran = splitter.R[:3,3]
        inv_rot = np.linalg.inv(rot)

        wp1 = splitter.wp1_lists[id]          # wp_G
        ops_H1 = splitter.H_orbits[id][0]     # operations of wp_h1
        op_G21 = splitter.G2_orbits[id][0][0] # operation 1 of wp_h1 in subgroup
        ops_G22 = splitter.G2_orbits[id][1]   # operations of wp_h2 in subgroup

        coord1_H, coord2_H = coord_H[0], coord_H[1]
        coord1_G2, coord2_G2  = coord1_H + translation, coord2_H + translation
        
        # since rotation does not change, search for the closest match on rotation
        # then we can get the translation vector
        for op_G22 in ops_G22:
            diff = (op_G22.rotation_matrix - op_G21.rotation_matrix).flatten()
            if np.sum(diff**2) < 1e-3:
                trans = op_G22.translation_vector - op_G21.translation_vector
                break
        trans -= np.round(trans)
        coords11 = apply_ops(coord1_G2, ops_H1)
        coords11 += trans
        tmp, dist = get_best_match(coords11, coord2_G2, self.cell)

        # needed displacement 
        d = coord2_G2 - tmp
        d -= np.round(d)

        if run_type == 1:
            return np.linalg.norm(np.dot(d/2, self.cell))
        else:
            op_G11 = splitter.G1_orbits[id][0][0]
            coord2_G2 -= d/2
            coord1_G2 += d/2
            coord1_G1, _ = search_G1(splitter.G, rot, tran, tmp, wp1, op_G11)
            return coord1_G1, [coord1_G2, coord2_G2], coord_H

    def symmetrize_site_double_t(self, splitter, id, coord_H, translation, run_type=1):
        """
        symmetrize two WPs (wp_h1, wp_h2) to another wp_G with higher symmetry 
        assuming a zero translation

        Args:
            splitter: splitter object
            id: the id in the splitter
            coord_H: coordinates to work on
            translation: 1*3 transaltion vector 
            run_type: return distance or coordinates
        """

        if translation is None: translation = np.zeros(3)
        rot = splitter.R[:3,:3]
        tran = splitter.R[:3,3]
        inv_rot = np.linalg.inv(rot)
        cell_G = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.cell)
        wp1 = splitter.wp1_lists[id]          # wp_G
        ops_G11 = splitter.G1_orbits[id][0]   # operations of wp_h1 in subgroup
        ops_G12 = splitter.G1_orbits[id][1]   # operations of wp_h2 in subgroup
        ops_G1 = splitter.G[0]               # general operations of G

        coord1_H, coord2_H = coord_H[0], coord_H[1] #coordinates in H
        coord1_G2, coord2_G2  = coord1_H+translation, coord2_H+translation # in G


        # forward search for the best generator for wp_h1 and wp_h2 in subgroup
        coord1_G1, _ = search_G1(splitter.G, rot, tran, coord1_G2, wp1, ops_G11[0])
        coord2_G1, _ = search_G1(splitter.G, rot, tran, coord2_G2, wp1, ops_G12[0])

        # apply the operations in G 
        # find the position that is closest to coord2_G1 
        coords11 = apply_ops(coord1_G1, ops_G1)
        tmp, dist = get_best_match(coords11, coord2_G1, cell_G)

        #self.print_wp(splitter, id)
        if run_type == 1:
            return dist/2
        else:
            # G1->G2->H
            d = coord2_G1 - tmp
            d -= np.round(d)
            coord2_G1 -= d/2

            coords22 = apply_ops(coord2_G1, ops_G1)
            coord1_G1, dist = get_best_match(coords22, coord1_G1, cell_G)
            #print("in G", l1, coord1_G1, l2, coord2_G1)

            # backward search (G->H)
            coord1_G2, dist = search_G2(inv_rot, -tran, coord1_G1, coord1_G2, self.cell)
            coord2_G2, dist = search_G2(inv_rot, -tran, coord2_G1, coord2_G2, self.cell)
            #print("in G1", l1, coord1_G2, l2, coord2_G2, dist)
            
            return coord2_G1, np.array([coord1_G2, coord2_G2]), coord_H

    def symmetrize_site_multi(self, splitter, id, coord_H, translation, run_type=1):
        """
        symmetrize multiple WPs to another with higher symmetry 

        Args:
            splitter: splitter object
            id: the id in the splitter
            coord_H: coordinates to work on
            translation: 1*3 transaltion vector 
            run_type: return distance or coordinates
        """

        if translation is None: translation = np.zeros(3)
        
        n = len(splitter.wp2_lists[id])
        rot = splitter.R[:3,:3]
        tran = splitter.R[:3,3]
        inv_rot = np.linalg.inv(rot)
        cell_G = np.dot(np.linalg.inv(splitter.R[:3,:3]).T, self.cell)
        wp1 = splitter.wp1_lists[id]

        # Finds the correct quadrant to easily generate all possible_wycs
        # add translations when trying to match
        quadrant = np.array(splitter.G2_orbits[id][0][0].as_dict()['matrix'])[:3,3]
        for k in range(3):
            if quadrant[k] >= 0.:
                quadrant[k] = 1
            else:
                quadrant[k] = -1

        coord_G2 = coord_H + translation
        coord_G2 %= quadrant

        # uses 1st coordinate and 1st wyckoff position as starting example.
        # Finds the matching G2 operation based on the nearest G1 search
        dist_list = []
        coord_list = []
        index = []
        G2_xyz = np.zeros([n,3])
        corresponding_ops = []

        for op in splitter.G1_orbits[id][0]:
            coord, dist = search_G1(splitter.G, rot, tran, coord_G2[0], wp1, op)
            dist_list.append(dist)
            coord_list.append(coord)

        dist_list = np.array(dist_list)
        index.append(np.argmin(dist_list))
        corresponding_ops.append(splitter.G2_orbits[id][0][index[0]])

        # Finds the free parameters xyz in the G2 basis for this coordinate
        G2_xyz[0] += find_xyz(corresponding_ops[0], coord_G2[0], quadrant)

        # Systematically generates possible G2 positions to match the remaining coordinates
        # Also finds the corresponding G2 free parameters xyz for each coordinate
        for j in range(1, n):
            possible_coords = [x.operate(G2_xyz[0]) for x in splitter.G2_orbits[id][j]]
            corresponding_coord, _ = get_best_match(possible_coords, coord_G2[j], cell_G)
            index.append([np.all(x==corresponding_coord) for x in possible_coords].index(True))
            corresponding_ops.append(splitter.G2_orbits[id][j][index[j]])
            G2_xyz[j] += find_xyz(corresponding_ops[j], coord_G2[j], quadrant)
        #print(G2_xyz)

        # Finds the average free parameters between all the coordinates as the best set 
        # of free parameters that all coordinates must match
        final_xyz = np.mean(G2_xyz, axis=0)

        if run_type == 1:
            dist_list = []
            for j in range(n):
                d = np.dot(G2_xyz[j] - final_xyz, cell_G)
                dist_list.append(np.linalg.norm(d))
            return max(dist_list)
        else:
            coords_G1 = np.zeros([n, 3]) #xyz in G1
            coords_G2 = np.zeros([n, 3]) #xyz in G2
            #dist_list = []
            for j in range(n):
                coords_G1[j] = splitter.G1_orbits[id][j][index[j]].operate(final_xyz)
                tmp = coord_H[j] + translation
                coords_G2[j], dist = search_G2(inv_rot, -tran, coords_G1[j], tmp, self.cell)
                #dist_list.append(dist)
            #print("dist", dist)
            #print("G1", coords_G1)#; import sys; sys.exit()
            #print("G2", coords_G2)#; import sys; sys.exit()
            return coords_G1[0], coords_G2, coord_H

    def print_detail(self, solution, coords_H, coords_G, elements):
        """
        print out the details of tranformation
        """
        print("\nTransition: ", self.struc.group.number, '->', self.G.number)
        (sp, mapping, translation, max_disp) = solution
        print('Maximum displacement: {:6.3f}'.format(max_disp))
        print('Mapping:', mapping)

        count = 0
        disps = []
        for i, wp2 in enumerate(sp.wp2_lists):
            wp1 = sp.wp1_lists[i]
            strs = ''
            for wp in wp2:
                x, y, ele = coords_H[count], coords_G[count], elements[count]
                label = str(wp.multiplicity) + wp.letter + '->' 
                label += str(wp1.multiplicity) + wp1.letter 
                dis = y - x - translation
                dis -= np.round(dis)
                dis_abs = np.linalg.norm(dis.dot(self.cell))
                output = "{:2s}[{:8s}] {:8.4f}{:8.4f}{:8.4f}".format(ele, label, *x)
                output += " -> {:8.4f}{:8.4f}{:8.4f}".format(*y)
                output += " -> {:8.4f}{:8.4f}{:8.4f} {:8.4f}".format(*dis, dis_abs)
                count += 1
                disps.append(dis_abs)
                print(output)
        print("cell: {:8.4f}{:8.4f}{:8.4f}, disp (A): {:8.4f}".format(*translation, max(disps)))

    def make_pyxtals_in_subgroup(self, solution, N_images=5):
        """
        make the pyxtal according to the given solution

        Args: 
            - solution: a tuple of (sp, mapping, translation, max_disp)
            - N_images: number of images

        Return:
            a list of pyxtal structures in low symmetry
        """
        (sp, mapping, translation, max_disp) = solution
        details = self.symmetrize(sp, mapping, translation)
        _, coords_G2, coords_H1, elements = details
        
        #print(coords_G2)
        #print(coords_H1)
        #self.print_detail(solution, coords_H1, coords_G2, elements)
        # Get the list of atomic displacements
        disps = []
        count = 0
        for wp2 in sp.wp2_lists:
            for wp in wp2:
                x, y, ele = coords_H1[count], coords_G2[count], elements[count]
                disp = y - x - translation
                disp -= np.round(disp)
                disps.append(disp)
                count += 1
        # Create the PyXtals
        strucs = []
        disps = np.array(disps)
        disps /= (N_images-1)
        max_disp = np.max(np.linalg.norm(disps.dot(self.cell), axis=1))
        for i in range(N_images):
            coords = coords_H1 + i*disps + translation 
            #print('sub', i, coords_H1, translation, coords)
            struc = self._make_pyxtal(sp, coords, elements, 1, False)
            struc.source = 'supergroup {:d} {:6.3f}'.format(i, max_disp*i)
            strucs.append(struc)
        return strucs

    def make_pyxtal_in_supergroup(self, solution):
        """
        make the pyxtal according to the given solution

        Args: 
            - solution: a tuple of (sp, mapping, translation, max_disp)

        Return:
            a pyxtal structure in high symmetry
        """
        (sp, mapping, translation, max_disp) = solution
        details = self.symmetrize(sp, mapping, translation)
        coords_G1, coords_G2, coords_H1, elements = details
        struc = self._make_pyxtal(sp, coords_G1)
        struc.source = 'supergroup {:6.3f}'.format(max_disp)
        struc.disp = max_disp
        return struc

    def _make_pyxtal(self, sp, coords, elements=None, run_type=0, check=True):
        """
        create the pyxtal with high/low symmetries

        Args:
            sp: splitter object
            coords: coordinates for each WP
            run_type: 0: high symmetry, otherwise low symmetry

        Return:
            pyxtal structure
        """
        # Create the pyxtal
        struc = self.struc.copy()
        cell_G = np.dot(np.linalg.inv(sp.R[:3,:3]).T, self.cell)
        lattice_G = Lattice.from_matrix(cell_G, ltype=sp.G.lattice_type)

        # Collect the atom_sites
        G_sites = []
        if run_type == 0:
            for i, wp in enumerate(sp.wp1_lists):
                pos = coords[i]
                pos -= np.floor(pos)
                if check:
                    pos1 = sym.search_matched_position(sp.G, wp, pos)
                    if pos1 is not None:
                        site = atom_site(wp, pos1, sp.elements[i])
                        G_sites.append(site)
                    else:
                        print("Group:", self.struc.group.number)
                        print("Position:", pos)
                        print(wp)
                        raise RuntimeError("cannot assign the right wp")
                else:
                    pos1 = pos
            # Update space group and lattice
            struc.group = sp.G
            struc.lattice = lattice_G
        else:
            count = 0
            for wp2 in sp.wp2_lists:
                for wp in wp2:
                    pos = coords[count]
                    pos -= np.floor(pos)
                    pos1 = sym.search_matched_position(sp.H, wp, pos)
                    if pos1 is not None:
                        site = atom_site(wp, pos1, elements[count])
                        G_sites.append(site)
                        count += 1
                    else:
                        print("Position:", pos)
                        print(wp)
                        print(sp)
                        raise RuntimeError("cannot assign the right wp")

            cell_U = np.dot(sp.R[:3, :3].T, lattice_G.matrix)
            struc.lattice = Lattice.from_matrix(cell_U, ltype=sp.H.lattice_type)
            #for i, coord in enumerate(coords):
            #    struc.atom_sites[i].specie = elements[i]
            #    struc.atom_sites[i].update(coord)
        struc.atom_sites = G_sites
        struc._get_formula()
        return struc

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

    def write_cifs(self):
        """
        Dump the cif files in sequence
        """
        for i, struc in enumerate(self.strucs):
            struc.to_file(str(i)+'-G'+str(struc.group.number)+'.cif')


if __name__ == "__main__":

    from pyxtal import pyxtal
    from time import time
    data = {
            #"PVO": [12, 166],
            #"PPO": [12],
            #"BTO": [123, 221],
            #"lt_cristobalite": [98, 210, 227],
            #"BTO-Amm2": [65, 123, 221],
            #"NaSb3F10": [186, 194],
            "NaSb3F10": [176, 194],
            #"NbO2": 141,
            #"GeF2": 62,
            #"lt_quartz": 180,
            #"NiS-Cm": 160,
            #"BTO-Amm2": 221,
            #"BTO": 221,
            #"lt_cristobalite": 227,
            #"NaSb3F10": 194,
            #"MPWO": 225,
            #"MPWO": [59, 71, 139, 225],
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
            sup = supergroups(s, G=data[cif], show=True, max_per_G=2500)
        print(sup)
        strs = "====================================================="
        strs += "==============={:12.3f} seconds".format(time()-t0)
        print(strs)
        #my = supergroup(s, data[cif])
        #sols = my.search_supergroup(max_solutions=12)
        ##my.make_supergroup(sols, show_detail=True)
        #for sol in sols:
        #    struc_high = my.make_pyxtal_in_supergroup(sol)
        #    strucs = my.make_pyxtals_in_subgroup(sol)
        #    pmg1 = struc_high.to_pymatgen()
        #    pmg2 = strucs[-1].to_pymatgen()
        #    rms = sm.StructureMatcher().get_rms_dist(pmg1, pmg2)
        #    print('========================================', rms)
        #    if not sm.StructureMatcher().fit(pmg1, pmg2):
        #        print(struc_high)
        #        print(strucs[-1])
