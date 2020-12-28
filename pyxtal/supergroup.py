import pyxtal.symmetry as sym
from pyxtal.operations import apply_ops, get_inverse
from pyxtal.wyckoff_split import wyckoff_split
import numpy as np
from copy import deepcopy
import itertools
from scipy.optimize import minimize

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
    def __init__(self, H, elements, sites, positions, cell, group_type='t'):

        # initilize the necesary parameters
        self.H = sym.Group(H)
        self.cell = cell
        self.group_type = group_type
        self.wyc_supergroups=self.H.get_min_supergroup(group_type)
        if len(elements) != len(sites):
            raise ValueError("The lengths of elements and sites should be equal")

        # group the elements, sites, positions
        self.elements = []
        self.sites = [] 
        self.positions = [] 
        for i, e in enumerate(elements):
            if e not in self.elements:
                self.elements.append(e)
                self.sites.append([sites[i]])
                self.positions.append([positions[i]])
            else:
                id = self.elements.index(e)
                self.sites[id].append(sites[i])
                self.positions[id].append(positions[i])
        
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
        # find the mapping for each element
        coords_G = []
        coords_H1 = []
        ids = []
        #disp = np.array([0.0, 0.0, 0.222222])
        disp = None

        valid = True
        # loop over each site for the given elements
        for i in range(len(self.elements)):
            coords_H, sites_G, sites_H = self.positions[i], solution[i], self.sites[i]
            sites_G.reverse()
            splitter = wyckoff_split(G, split_id, sites_G, self.group_type)
            res = self.symmetrize(sites_G, sites_H, coords_H, splitter, disp)

            #if sites_G == ['i', 'h', 'c'] or sites_G == ['c', 'h', 'i']:
            #    import sys; sys.exit()

            if res is not None:
                (disp, coord_G, coord_H1, id) = res
                coords_G.append(coord_G)
                coords_H1.append(coord_H1)
                ids.append(id)
            else:
                #print("Skip this solution")
                valid = False
                break

        if valid:
            #disp = np.array(disps).mean(axis=0)
            #print(disp)
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
            #break

    def symmetrize(self, sites_G, sites_H, coords_H, splitter, disp=None, d_tol=0.15):
        """
        For a given solution, search for the possbile supergroup structure

        Args: 
            - sites_G: e.g., ['2c', '6h', '12i']
            - sites_H: e.g., ['2b', '6c', '6c', '6c']
            - coords_H: e.g., N*3 array, where N is the length of sites_H
        Returns:
            coords_G with a minimum displacement
        """

        # search for the mapping between sites_G and sites_H
        selected_ids = []
        coords_G = []
        coords_H1 = []
        axis = 2 # special axis

        for i, wp1 in enumerate(splitter.wp1_lists):
            # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
            mult_G = wp1.multiplicity
            mult_Hs = [int(site[:-1]) for site in sites_H]

            valid = False                  
            remaining_ids = [id for id in range(len(coords_H)) if id not in selected_ids]

            if len(splitter.wp2_lists[i]) == 1:
                # one to one splitting, e.g., 2c->2d
                # this usually involves increase of site symmetry
                # therefore, some position has to fixed
                # e.g. from 173 (P63) to 176 (P63/m), 
                # 1, x, y, z (1) -> x, y, 1/4 (m..)
                # 2, x, y, z (1) -> 1/2, 0, 0 (-1)
                # for 1, the translation on x, y is always 0
                # z can be moved to either 1/4 or -1/4
                # for 2, needs to consider all displacement
               
                # symmetry info
                ops_H = splitter.H_orbits[i][0]  # ops for H
                ops_G2 = splitter.G2_orbits[i][0] # ops for G2
                xyz1 = ops_H[0].as_xyz_string().split(',')
                xyz2 = ops_G2[0].as_xyz_string().split(',')
                # mask is a vector to control the freedom of x, y, z       
                mask = [False if xyz1[m]==xyz2[m] else True for m in range(3)]

                # chose the site with same multiplicity
                possible_ids = [j for j in remaining_ids if mult_G == mult_Hs[j]]
                for j in possible_ids:
                    coord = coords_H[j]
                    
                    # refine coord1 to find the best match on coord2
                    # coord1 = apply_ops(coord, ops_G2)[0] # coord in G
                    coord0s = apply_ops(coord, ops_H) # possible coords in H
                    
                    dists = []
                    for coord0 in coord0s:
                        coord2 = coord0.copy()
                        if disp is not None:
                            coord2 += disp 
                        coord1 = apply_ops(coord2, ops_G2)[0] # coord in G
                        dist = coord1 - coord2
                        dist -= np.round(dist)
                        dists.append(np.linalg.norm(dist))
                    min_ID = np.argmin(np.array(dists))

                    dist = dists[min_ID]
                    coord2 = coord0s[min_ID].copy()
                    
                    if disp is None:
                        coord1 = apply_ops(coord2, ops_G2)[0]
                        disp = (coord1 - coord2).copy()
                        valid = True
                    elif dist < d_tol:
                        coord1 = ops_G2[0].operate(coord2+disp)
                        if round(np.trace(ops_G2[0].rotation_matrix)) in [1, 2]:
                            # optimize the distance by changing coord1
                            
                            def fun(x, pt, ref, op):
                                pt[0] = x[0]
                                y = op.operate(pt)
                                diff = y - ref
                                diff -= np.round(diff)
                                #print("dddddddddd", y, diff)
                                return np.linalg.norm(diff)

                            res = minimize(fun, coord1[0], args=(coord1, coord2, ops_G2[0]), 
                                    method='Nelder-Mead', options={'maxiter': 20})
                            coord1[0] = res.x[0]
                            coord1 = ops_G2[0].operate(coord1)
                        valid = True
                            
                    #if disp is None:
                    #    # create cell translation
                    #    coord2 = coord0s[0]
                    #    disp = (coord1 - coord2).copy()
                    #    valid = True
                    #else:
                    #    dists = coord1 - disp - coord0s
                    #    dists -= np.round(dists)
                    #    for k, dist in enumerate(dists):
                    #        if np.sqrt((dist[mask]**2).sum()) < d_tol:
                    #            valid = True
                    #            coord2 = coord0s[k]
                    #            for ax in range(3):
                    #                if not mask[ax]:
                    #                    coord2[ax] = coord1[ax]
                    #            #print("break out---------------", np.sqrt((dist[mask]**2).sum()), wp1.letter)
                    #            break
                    if valid:
                        selected_ids.append(j)
                        coords_G.append(coord1)
                        coords_H1.append(coord2)
                        #print("=============pass", coord2)
                        break
            else: 
                # site symmetry does not change
                # if merge is to let two clusters gain additional symmetry (m, -1)

                # symmetry operations
                ops_H1 = splitter.H_orbits[i][0]
                ops_H2 = splitter.H_orbits[i][1]
                ops_G21 = splitter.G2_orbits[i][0]
                ops_G22 = splitter.G2_orbits[i][1]

                # choose two coords which satisfy the mutiplicity constraints
                sub_ids = list(itertools.combinations(remaining_ids, 2))
                possible_ids = [id for id in sub_ids if mult_Hs[id[0]]==len(ops_H1) \
                             and mult_Hs[id[1]]==len(ops_H2)]
                # search for the solution with minimum disp
                dists = []
                for ids in possible_ids: 
                    coord1, coord2 = coords_H[ids[0]], coords_H[ids[1]]
                    # refine coord1 to find the best match on coord2
                    coords1 = apply_ops(coord1, ops_H2)
                    tmp, match_id = get_best_match(coords1, coord2, axis)
                    diffs = (tmp - coord2)/2
                    diffs -= np.round(diffs)
                    diffs[axis] = 0
                    dists.append(np.sqrt(np.sum(diffs**2)))
                min_id = np.argmin(np.array(dists))

                if dists[min_id] < d_tol:
                    ids = possible_ids[min_id]
                    coord1, coord2 = coords_H[ids[0]], coords_H[ids[1]]
                    coords1 = apply_ops(coord1, ops_H2)
                    tmp, match_id = get_best_match(coords1, coord2, axis)
                    coord1 = coords1[match_id]
                    avg = (coord1+coord2)/2
                    # compute the overall shift along the special axis
                    
                    #print(np.sum(diffs**2), disp[axis],  avg[axis], "CCCCCCCCCCCCCCCCCCCC")
                    if disp is None:
                        shift = 0.75 - avg[axis]
                        disp = np.zeros(3)
                        disp[axis] = shift
                        valid = True
                    else:
                        # [x, y, z] -> [x, y, 1/2-z]
                        for pt in [0.25, 0.75]:
                            if abs(avg[axis]+disp[axis]-pt) < d_tol:
                                valid = True
                                shift = pt - avg[axis]
                                break
                if valid:
                    d2 = np.zeros(3)
                    mydisp = np.zeros(3)

                    d2[axis] = (coord1[axis]-coord2[axis])/2
                    mydisp[axis] = shift

                    coords_G.append(avg + mydisp + d2) 
                    coords_G.append(avg + mydisp - d2) 

                    coords_H1.append(coord1)
                    coords_H1.append(coord2)
                    selected_ids.extend(ids)

            if not valid:
                #print("Cannot make it, skip")
                return None   
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
                                        G_sites.append(possible_wycs[l][-1])
                                    else:
                                        for i in range(int(number)):
                                            G_sites.append(possible_wycs[l][-1])
                                if G_sites not in good_list:
                                    good_list.append(G_sites)

                combo_storage=holder

            if len(good_list)==0:
                # print("cannot find the valid split, quit the search asap")
                return None
            else:
                good_splittings_list.append(good_list)

        return good_splittings_list

def get_best_match(positions, ref, axis):
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
    diffs[:, axis] = 0
    dists = np.linalg.norm(diffs, axis=1)
    id = np.argmin(dists)
    return positions[id], id

if __name__ == "__main__":

    from pyxtal import pyxtal
    
    s = pyxtal()
    #s.from_seed("pyxtal/database/cifs/BTO.cif")
    s.from_seed("pyxtal/database/cifs/NaSb3F10.cif")
    print(s)
    species = []
    sites = []
    coords = []
    for site in s.atom_sites:
        sites.append(str(site.multiplicity)+site.wp.letter)
        coords.append(site.position)
        species.append(site.specie)
    print(sites)
    H = s.group.number
    supergroup(H, species, sites, coords, s.lattice.matrix)
 
    
