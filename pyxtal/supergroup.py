import pyxtal.symmetry as sym
from pyxtal.operations import apply_ops
from pyxtal.wyckoff_split import wyckoff_split
import numpy as np
from copy import deepcopy
import itertools
from scipy.optimize import minimize

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

class supergroup():
    """
    Class to find the structure with supergroup symmetry

    Args:
        H: int (1-230), number of subgroup
        elements: list of elements ["Na", "Sb", "F", "F", "F", "F"]
        sites: the correpsonding wyckoff symbols ["2b", "6c", "6c", "6c", "6c", "2b"]
        idxs: index of which transformation to pick
        group_type: `t` or `k`
    """
    def __init__(self, struc, group_type='t'):

        # initilize the necesary parameters
        self.struc = struc
        self.wyc_supergroups = struc.group.get_min_supergroup(group_type)
        self.cell = struc.lattice.matrix
        self.group_type = group_type

        # group the elements, sites, positions
        self.elements = []
        self.sites = [] 
        for i, atom_site in enumerate(struc.atom_sites):
            e = atom_site.specie
            site = str(atom_site.wp.multiplicity) + atom_site.wp.letter
            if e not in self.elements:
                self.elements.append(e)
                self.sites.append([site])
            else:
                id = self.elements.index(e)
                self.sites[id].append(site)
        
        # perform actual search
        for idx in range(len(self.wyc_supergroups['supergroup'])):
            G = self.wyc_supergroups['supergroup'][idx]
            relation = self.wyc_supergroups['relations'][idx]
            id = self.wyc_supergroups['idx'][idx]
            results = self.check_compatibility(G, relation)
            if results is not None:
                solutions = list(itertools.product(*results))
                valid_solutions = self.check_freedom(G, solutions)
                for solution in valid_solutions:
                    print(G, solution)
                    self.make_supergroup(G, solution, id)

    def make_supergroup(self, G, solution, split_id):
        """
        For a given solution, search for the possbile supergroup structure

        Args: 
            - solution: e.g., [['2d'], ['6h'], ['2c', '6h', '12i']]
            - split_id: integer
        Returns:
            structure + displacement
        """
        # move it to the symmetrize function
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

        splitter = wyckoff_split(G, split_id, sites_G, self.group_type, elements)
        mappings = self.find_mapping(splitter)
        for maping in mappings:
            #disp = np.array([0.0, 0.0, 0.222222])
            res = self.symmetrize(splitter, mapping, disp=None)
            if res is not None:
                print("Valid structure in ", G)
                total_num = 0
                total_disp = 0
                for l, id in enumerate(ids):
                    ele = self.elements[l]
                    sites = np.array(self.sites[l])
                    for x, y, site in zip(coords_H1[l], coords_G[l], sites[id]):
                        num = int(site[:-1])
                        dis = y-(x+disp)
                        dis -= np.round(dis)
                        dis_abs = np.linalg.norm(dis.dot(self.cell))
                        output = "{:2s} {:7.3f}{:7.3f}{:7.3f}".format(ele, *x)
                        output += " -> {:7.3f}{:7.3f}{:7.3f}".format(*y)
                        output += " -> {:7.3f}{:7.3f}{:7.3f} {:7.3f}".format(*dis, dis_abs)
                        print(output)
                        total_disp += dis_abs*num
                        total_num += num
                mae = total_disp/total_num
                print("cell: {:7.3f}{:7.3f}{:7.3f}, disp (A): {:7.3f}".format(*disp, mae))
                #import sys; sys.exit()

    def find_mapping(self, splitter):
        """
        search for all mappings for a given splitter
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
        
    def symmetrize(self, splitter, solution, disp=None, d_tol=0.75):
        """
        For a given solution, search for the possbile supergroup structure

        Args: 
            - splitter: splitter object to specify the relation between G and H
            - disp: an overall shift from H to G, None or 3 vector
            - d_tol: the tolerance
        Returns:
            atom_site_G1 with the minimum displacements
            atom_site_G2 
        """
        atom_sites_H = self.struc.atom_sites

        # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
        for i, wp1 in enumerate(splitter.wp1_lists):

            if len(splitter.wp2_lists[i]) == 1:
                # one to one splitting, e.g., 2c->2d
                # this usually involves increase of site symmetry
               
                # symmetry info
                ops_H = splitter.H_orbits[i][0]  # ops for H
                ops_G2 = splitter.G2_orbits[i][0] # ops for G2
                
                # refine coord1 to find the best match on coord2
                coord = atom_sites_H[solution[i][0]].position
                coord0s = apply_ops(coord, ops_H) # possible coords in H
                dists = []
                for coord0 in coord0s:
                    coord2 = coord0.copy()
                    if disp is not None:
                        coord2 += disp 
                    coord1 = apply_ops(coord2, ops_G2)[0] # coord in G
                    dist = coord1 - coord2
                    dist -= np.round(dist)
                    dist = np.dot(dist, self.cell)
                    dists.append(np.linalg.norm(dist))
                min_ID = np.argmin(np.array(dists))

                dist = dists[min_ID]
                coord2 = coord0s[min_ID].copy()
                
                if disp is None:
                    coord1 = apply_ops(coord2, ops_G2)[0]
                    disp = (coord1 - coord2).copy()
                elif dist < d_tol:
                    coord1 = ops_G2[0].operate(coord2+disp)
                    if round(np.trace(ops_G2[0].rotation_matrix)) in [1, 2]:
                        def fun(x, pt, ref, op):
                            pt[0] = x[0]
                            y = op.operate(pt)
                            diff = y - ref
                            diff -= np.round(diff)
                            return np.linalg.norm(diff)

                        # optimize the distance by changing coord1
                        res = minimize(fun, coord1[0], args=(coord1, coord2, ops_G2[0]), 
                                method='Nelder-Mead', options={'maxiter': 20})
                        coord1[0] = res.x[0]
                        coord1 = ops_G2[0].operate(coord1)
                else:
                    return None
                        
            else: 
                # site symmetry does not change
                # if merge is to let two clusters gain additional symmetry (m, -1)
                # symmetry operations
                ops_H1 = splitter.H_orbits[i][0]
                ops_G22 = splitter.G2_orbits[i][1]

                coord1 = atom_sites_H[solution[i][0]].position.copy()
                coord2 = atom_sites_H[solution[i][0]].position.copy()
                # refine coord1 to find the best match on coord2
                if disp is not None:
                    coord11 = coord1 + disp
                    coord22 = coord2 + disp
                else:
                    coord11 = coord1
                    coord22 = coord2
                coords11 = apply_ops(coord11, ops_H1)
            
                # transform coords1 by symmetry operation
                op = ops_G22[0]
                inv_op = op.inverse
                for m, coord11 in enumerate(coords11):
                    coords11[m] = op.operate(coord11)
                tmp, dist = get_best_match(coords11, coord22, self.cell)
                if dist > np.sqrt(2)*d_tol:
                    return None
                else:
                    # recover the original position
                    coord1 = inv_op.operate(tmp)
                    if disp is not None:
                        coord1 -= disp
                    d = coord22 - tmp
                    d -= np.round(d)

                    coord22 -= d/2 #final coord2 after disp
                    # recover the displaced position
                    coord11 = inv_op.operate(coord22) 

                    coords_H1.append(coord1)
                    coords_H1.append(coord2)
                    coords_G.append(coord11)
                    coords_G.append(coord22)
                    
        return disp, coords_G, coords_H1, selected_ids

    def check_freedom(self, G, solutions):
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
                #if 'g' not in sites:
                #    print(sites)
        return valid_solutions

    def check_compatibility(self, G, relation):
        """
        find the compatible splitter too let the atoms of subgroup H fit the group G.

        Args:
            G: the target space group with high symmetry
            relation: a dictionary to describe the relation between G and H
        """

        G = sym.Group(G)

        results = {}
        wyc_list = [(str(x.multiplicity)+x.letter) for x in G]
        wyc_list.reverse()

        good_splittings_list=[]

        # A lot of integer math below. 
        # The goal is to find all the integer combinations of supergroup 
        # wyckoff positions with the same number of atoms

        # each element is solved one at a time
        for i, ele in enumerate(self.elements):

            site = np.unique(self.sites[i])
            site_counts = [self.sites[i].count(x) for x in site]
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

            while len(combo_storage)!=0:
                holder = []
                for j, x in enumerate(combo_storage):
                    for k in range(len(block_units)):
                        trial = np.array(deepcopy(x)) # trial solution
                        trial[k] += 1
                        sum_units = np.dot(trial, block_units)
                        if sum_units > total_units:
                            continue
                        elif sum_units < total_units:
                            holder.append(trial)
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

        return good_splittings_list

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

    diffs = positions - ref)
    diffs -= np.round(diffs)
    diffs = np.dot(diffs, cell)
    dists = np.linalg.norm(diffs, axis=1)
    id = np.argmin(dists)
    return positions[id], dists[id]

if __name__ == "__main__":

    from pyxtal import pyxtal
    
    s = pyxtal()
    #s.from_seed("pyxtal/database/cifs/BTO.cif")
    s.from_seed("pyxtal/database/cifs/NaSb3F10.cif")
    print(s)
    supergroup(s)
 
    
