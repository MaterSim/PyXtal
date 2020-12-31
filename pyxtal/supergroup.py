import pyxtal.symmetry as sym
from pyxtal.lattice import Lattice
from pyxtal.wyckoff_site import atom_site
from pyxtal.operations import apply_ops, get_inverse
from pyxtal.wyckoff_split import wyckoff_split
import pymatgen.analysis.structure_matcher as sm
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


class supergroup():
    """
    Class to find the structure with supergroup symmetry

    Args:
        struc: pyxtal structure
        G: list of possible supergroup numbers, default is None 
        group_type: `t` or `k`
    """
    def __init__(self, struc, G=None, group_type='t', solution=None):

        # initilize the necesary parameters
        self.struc = struc
        self.cell = struc.lattice.matrix
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

        # group the elements, sites, positions
        self.elements = []
        self.sites = [] 
        for i, at_site in enumerate(struc.atom_sites):
            e = at_site.specie
            site = str(at_site.wp.multiplicity) + at_site.wp.letter
            if e not in self.elements:
                self.elements.append(e)
                self.sites.append([site])
            else:
                id = self.elements.index(e)
                self.sites[id].append(site)
        
        # search for the compatible solutions
        if solution is None:
            self.solutions = []
            for idx in range(len(self.wyc_supergroups['supergroup'])):
                G = self.wyc_supergroups['supergroup'][idx]
                relation = self.wyc_supergroups['relations'][idx]
                id = self.wyc_supergroups['idx'][idx]
                results = self.check_compatibility(G, relation)
                if results is not None:
                    solutions = list(itertools.product(*results))
                    trials = self.check_freedom(G, solutions)
                    sol = {'group': G, 'id': id, 'splits': trials}
                    self.solutions.append(sol)
        else: # load the solution
            raise NotImplementedError

        if len(self.solutions) == 0:
            self.error = True
            print("No compatible solution exists")

    def search_supergroup(self, d_tol=1.0):
        """
        search for valid supergroup transition
        
        Args:
            d_tol (float): tolerance

        Returns:
            valid_solutions: dictionary
        """
        self.d_tol = d_tol
        # extract the valid 
        valid_solutions = []
        for sols in self.solutions:
            G, id, sols = sols['group'], sols['id'], sols['splits']
            for sol in sols:
                mae, disp, mapping, sp = self.get_displacement(G, id, sol, d_tol*1.1)
                if mae < d_tol:
                    valid_solutions.append((sp, mapping, disp, mae))
        return valid_solutions

    def make_supergroup(self, solutions, show_detail=True):
        """
        create supergroup structures based on the list of solutions

        Args: 
            solutions: list of tuples (splitter, mapping, disp)

        Returns:
            list of pyxtal structures
        """
        G_strucs = []
        for solution in solutions:
            (sp, mapping, disp, mae) = solution
            G = sp.G.number
            details = self.symmetrize(sp, mapping, disp)
            coords_G1, coords_G2, coords_H1, elements, mults = details

            G_struc = self.struc.copy()
            G_struc.group = sp.G

            G_sites = []
            for i, wp in enumerate(sp.wp1_lists):
                site = atom_site(wp, coords_G1[i], sp.elements[i])
                G_sites.append(site)

            G_struc.atom_sites = G_sites
            G_struc.source = 'supergroup {:6.3f}'.format(mae) 

            lat1 = np.dot(sp.inv_R[:3,:3].T, self.struc.lattice.matrix)
            lattice = Lattice.from_matrix(lat1, ltype=sp.G.lattice_type)
            #lattice.reset_matrix() #make it has a regular shape
            #G_struc.lattice = self.struc.lattice.supergroup(sp.G.lattice_type)
            G_struc.lattice = lattice
            G_struc.numIons *= round(np.abs(np.linalg.det(sp.R[:3,:3])))
            G_struc._get_formula()

            if new_structure(G_struc, G_strucs):
                G_strucs.append(G_struc)
                if show_detail:
                    details = self.symmetrize(sp, mapping, disp)
                    _, coords_G2, coords_H1, elements, mults = details
                    self.print_detail(G, coords_H1, coords_G2, elements, mults, disp)
                    print(G_struc)
                    #print(sp.R)

        return G_strucs
            
    def get_displacement(self, G, split_id, solution, d_tol):
        """
        For a given solution, search for the possbile supergroup structure

        Args: 
            G: supergroup number 
            split_id: integer
            solution: e.g., [['2d'], ['6h'], ['2c', '6h', '12i']]
            d_tol: 
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
        for mapping in mappings:
            disp = None #np.array([0.0, 0.0, 0.222222])
            dist, disp = self.symmetrize_dist(splitter, mapping, disp, d_tol)
            dists.append(dist)
            disps.append(disp)
        dists = np.array(dists)
        mae = np.min(dists)
        id = np.argmin(dists)
        if (mae > 0.2) and (mae < d_tol):
            # optimize further
            def fun(disp, mapping, splitter):
                return self.symmetrize_dist(splitter, mapping, disp)[0]
            res = minimize(fun, disps[id], args=(mappings[id], splitter),
                    method='Nelder-Mead', options={'maxiter': 20})
            if res.fun < mae:
                mae = res.fun
                disp = res.x
        return mae, disp, mappings[id], splitter

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
        
    def symmetrize_dist(self, splitter, solution, disp=None, d_tol=1.2):
        """
        For a given solution, search for the possbile supergroup structure

        Args: 
            splitter: splitter object to specify the relation between G and H
            solution: list of sites in H, e.g., ['4a', '8b']
            disp: an overall shift from H to G, None or 3 vector
            d_tol: the tolerance in angstrom
        Returns:
            distortion
            cell translation
        """
        total_disp = 0
        atom_sites_H = self.struc.atom_sites
        n_atoms = sum([site.wp.multiplicity for site in atom_sites_H])
        #print("checking solution-----------------", solution)
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
                            diff = np.dot(diff, self.cell)
                            return np.linalg.norm(diff)

                        # optimize the distance by changing coord1
                        res = minimize(fun, coord1[0], args=(coord1, coord2, ops_G2[0]), 
                                method='Nelder-Mead', options={'maxiter': 20})
                        coord1[0] = res.x[0]
                        coord1 = ops_G2[0].operate(coord1)
                else:
                    return 10000, None
                
                diff = coord1-(coord2+disp)
                diff -= np.round(diff)
                total_disp += np.linalg.norm(np.dot(diff, self.cell))*len(ops_H)
                        
            else: 
                # symmetry operations
                ops_H1 = splitter.H_orbits[i][0]
                ops_G22 = splitter.G2_orbits[i][1]

                coord1 = atom_sites_H[solution[i][0]].position.copy()
                coord2 = atom_sites_H[solution[i][1]].position.copy()
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
                for m, coord11 in enumerate(coords11):
                    coords11[m] = op.operate(coord11)
                tmp, dist = get_best_match(coords11, coord22, self.cell)
                if dist > np.sqrt(2)*d_tol:
                    return 10000, None

                # recover the original position
                #print(op)
                #print(op.as_xyz_string())
                inv_op = get_inverse(op)
                coord1 = inv_op.operate(tmp)
                if disp is not None:
                    coord1 -= disp
                d = coord22 - tmp
                d -= np.round(d)

                coord22 -= d/2 #final coord2 after disp
                # recover the displaced position
                coord11 = inv_op.operate(coord22) 

                total_disp += np.linalg.norm(np.dot(d/2, self.cell))*len(ops_H1)*2

        return total_disp/n_atoms, disp

    def symmetrize(self, splitter, solution, disp):
        """
        For a given solution, search for the possbile supergroup structure

        Args: 
            splitter: splitter object to specify the relation between G and H
            disp: an overall shift from H to G, None or 3 vector
            d_tol: the tolerance in angstrom

        Returns:
            coords_G1
            coords_G2
            coords_H1
            elements
            mults
        """

        atom_sites_H = self.struc.atom_sites
        coords_G1 = [] # position in G
        coords_G2 = [] # position in G on the subgroup bais
        coords_H1 = [] # position in H
        elements = []
        mults = []
        inv_R = splitter.inv_R # inverse coordinate transformation
        # wp1 stores the wyckoff position object of ['2c', '6h', '12i']
        for i, wp1 in enumerate(splitter.wp1_lists):

            if len(splitter.wp2_lists[i]) == 1:
                # symmetry info
                ops_H = splitter.H_orbits[i][0]  # ops for H
                ops_G2 = splitter.G2_orbits[i][0] # ops for G2
                
                # refine coord1 to find the best match on coord2
                coord = atom_sites_H[solution[i][0]].position
                coord0s = apply_ops(coord, ops_H) # possible coords in H
                dists = []
                for coord0 in coord0s:
                    coord2 = coord0.copy()
                    coord2 += disp 
                    coord1 = apply_ops(coord2, ops_G2)[0] # coord in G
                    dist = coord1 - coord2
                    dist -= np.round(dist)
                    dist = np.dot(dist, self.cell)
                    dists.append(np.linalg.norm(dist))
                min_ID = np.argmin(np.array(dists))

                coord2 = coord0s[min_ID].copy()
                coord1 = ops_G2[0].operate(coord2+disp)
                if round(np.trace(ops_G2[0].rotation_matrix)) in [1, 2]:
                    def fun(x, pt, ref, op):
                        pt[0] = x[0]
                        y = op.operate(pt)
                        diff = y - ref
                        diff -= np.round(diff)
                        diff = np.dot(diff, self.cell)
                        return np.linalg.norm(diff)

                    # optimize the distance by changing coord1
                    res = minimize(fun, coord1[0], args=(coord1, coord2, ops_G2[0]), 
                            method='Nelder-Mead', options={'maxiter': 20})
                    coord1[0] = res.x[0]
                    coord1 = ops_G2[0].operate(coord1)

                coords_G1.append(np.dot(inv_R[:3,:3], coord1).T+inv_R[:3,3].T)
                #coords_G1.append(coord1)
                coords_G2.append(coord1)
                coords_H1.append(coord2)
                elements.append(splitter.elements[i])
                mults.append(len(ops_G2))
                        
            else: 
                # symmetry operations
                ops_H1 = splitter.H_orbits[i][0]
                ops_G22 = splitter.G2_orbits[i][1]

                coord1 = atom_sites_H[solution[i][0]].position.copy()
                coord2 = atom_sites_H[solution[i][1]].position.copy()
                # refine coord1 to find the best match on coord2
                coord11 = coord1 + disp
                coord22 = coord2 + disp
                coords11 = apply_ops(coord11, ops_H1)

                # transform coords1 by symmetry operation
                op = ops_G22[0]
                for m, coord11 in enumerate(coords11):
                    coords11[m] = op.operate(coord11)
                tmp, dist = get_best_match(coords11, coord22, self.cell)

                # recover the original position
                inv_op = get_inverse(op)
                coord1 = inv_op.operate(tmp)
                coord1 -= disp
                d = coord22 - tmp
                d -= np.round(d)

                coord22 -= d/2 #final coord2 after disp
                # recover the displaced position
                coord11 = inv_op.operate(coord22) 

                coords_G1.append(np.dot(inv_R[:3,:3], coord11).T+inv_R[:3,3].T)
                coords_G2.append(coord11)
                coords_G2.append(coord22)
                coords_H1.append(coord1)
                coords_H1.append(coord2)
                elements.extend([splitter.elements[i]]*2)
                mults.append(len(ops_H1))
                mults.append(len(ops_G22))
        return coords_G1, coords_G2, coords_H1, elements, mults

    def print_detail(self, G, coords_H1, coords_G1, elements, mults, disp):
        """
        print out the details of tranformation
        """
        print("Valid structure", G)
        total_disp = 0
        for x, y, ele, n in zip(coords_H1, coords_G1, elements, mults):
            dis = y-(x+disp)
            dis -= np.round(dis)
            dis_abs = np.linalg.norm(dis.dot(self.cell))
            output = "{:2s} {:7.3f}{:7.3f}{:7.3f}".format(ele, *x)
            output += " -> {:7.3f}{:7.3f}{:7.3f}".format(*y)
            output += " -> {:7.3f}{:7.3f}{:7.3f} {:7.3f}".format(*dis, dis_abs)
            total_disp += dis_abs*n
            print(output)
        mae = total_disp/sum(mults)
        print("cell: {:7.3f}{:7.3f}{:7.3f}, disp (A): {:7.3f}".format(*disp, mae))


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
            #else:
            #    print(solution)
            #    import sys; sys.exit()
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

    diffs = positions - ref
    diffs -= np.round(diffs)
    diffs = np.dot(diffs, cell)
    dists = np.linalg.norm(diffs, axis=1)
    id = np.argmin(dists)
    return positions[id], dists[id]

if __name__ == "__main__":

    from pyxtal import pyxtal
    
    s = pyxtal()
    #s.from_seed("pyxtal/database/cifs/BTO-P4mm.cif")
    #s.from_seed("pyxtal/database/cifs/NaSb3F10.cif")
    #s.from_seed("pyxtal/database/cifs/lt_cristobalite.cif")
    #s.from_seed("pyxtal/database/cifs/lt_quartz.cif")
    #s.from_seed("pyxtal/database/cifs/GeF2.cif")
    #s.from_seed("pyxtal/database/cifs/B28.cif")
    #s.from_seed("pyxtal/database/cifs/PPO.cif")
    s.from_seed("pyxtal/database/cifs/BTO-Amm2.cif")
    print(s)
    #my = supergroup(s, G=[165, 167])
    my = supergroup(s)
    solutions = my.search_supergroup(d_tol=0.8)
    G_strucs = my.make_supergroup(solutions)
    G_strucs[-1].to_ase().write('1.vasp', format='vasp', vasp5=True, direct=True)

    my = supergroup(G_strucs[-1])
    solutions = my.search_supergroup(d_tol=0.60)
    G_strucs = my.make_supergroup(solutions)
    G_strucs[-1].to_ase().write('2.vasp', format='vasp', vasp5=True, direct=True)

    my = supergroup(G_strucs[-1])
    solutions = my.search_supergroup(d_tol=0.60)
    G_strucs = my.make_supergroup(solutions)
    G_strucs[-1].to_ase().write('3.vasp', format='vasp', vasp5=True, direct=True)

