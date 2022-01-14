"""
Module for storing and accessing symmetry group information, including a Group class and
a Wyckoff_Position class. These classes are used for generation of random structures.
"""
# Imports
# ------------------------------
# Standard Libraries
import numpy as np
from pkg_resources import resource_filename as rf
from copy import deepcopy
import random
import itertools
import re

# External Libraries
from pymatgen.symmetry.analyzer import generate_full_symmops
from pandas import read_csv
from monty.serialization import loadfn

# PyXtal imports
from pyxtal.msg import printx
from pyxtal.operations import (
    SymmOp,
    apply_ops,
    filtered_coords,
    filtered_coords_euclidean,
    distance,
    distance_matrix,
    create_matrix,
    OperationAnalyzer,
    check_images,
)
from pyxtal.database.hall import hall_from_hm
from pyxtal.constants import letters
# ------------------------------ Constants ---------------------------------------

wyckoff_df = read_csv(rf("pyxtal", "database/wyckoff_list.csv"))
wyckoff_symmetry_df = read_csv(rf("pyxtal", "database/wyckoff_symmetry.csv"))
wyckoff_generators_df = read_csv(rf("pyxtal", "database/wyckoff_generators.csv"))
layer_df = read_csv(rf("pyxtal", "database/layer.csv"))
layer_symmetry_df = read_csv(rf("pyxtal", "database/layer_symmetry.csv"))
layer_generators_df = read_csv(rf("pyxtal", "database/layer_generators.csv"))
rod_df = read_csv(rf("pyxtal", "database/rod.csv"))
rod_symmetry_df = read_csv(rf("pyxtal", "database/rod_symmetry.csv"))
rod_generators_df = read_csv(rf("pyxtal", "database/rod_generators.csv"))
point_df = read_csv(rf("pyxtal", "database/point.csv"))
point_symmetry_df = read_csv(rf("pyxtal", "database/point_symmetry.csv"))
point_generators_df = read_csv(rf("pyxtal", "database/point_generators.csv"))
symbols = loadfn(rf("pyxtal", "database/symbols.json"))
t_subgroup = loadfn(rf("pyxtal",'database/t_subgroup.json'))
k_subgroup = loadfn(rf("pyxtal",'database/k_subgroup.json'))
wyc_sets = loadfn(rf("pyxtal",'database/wyckoff_sets.json'))
hex_cell = np.array([[1, -0.5, 0], [0, np.sqrt(3) / 2, 0], [0, 0, 1]])
# --------------------------- Group class -----------------------------
class Group:
    """
    Class for storing a set of Wyckoff positions for a symmetry group. See the documentation
    for details about settings.

    Examples
    --------
    >>> from pyxtal.symmetry import Group
    >>> g = Group(64)
    >>> g
    -- Spacegroup --# 64 (Cmce)--
    16g	site symm: 1
    8f	site symm: m..
    8e	site symm: .2.
    8d	site symm: 2..
    8c	site symm: -1
    4b	site symm: 2/m..
    4a	site symm: 2/m..

    one can access data with its attributes such as `symbol`, `number` and `Wyckoff_positions`

    >>> g.symbol
    'Cmce'
    >>> g.number
    64
    >>> g.Wyckoff_positions[0]
    Wyckoff position 16g in space group 64 with site symmetry 1
    x, y, z
    -x, -y+1/2, z+1/2
    -x, y+1/2, -z+1/2
    x, -y, -z
    -x, -y, -z
    x, y+1/2, -z+1/2
    x, -y+1/2, z+1/2
    -x, y, z
    x+1/2, y+1/2, z
    -x+1/2, -y+1, z+1/2
    -x+1/2, y+1, -z+1/2
    x+1/2, -y+1/2, -z
    -x+1/2, -y+1/2, -z
    x+1/2, y+1, -z+1/2
    x+1/2, -y+1, z+1/2
    -x+1/2, y+1/2, z

    We also provide several utilities functions, e.g.,
    one can search the possible wyckoff_combinations by a formula

    >>> g.list_wyckoff_combinations([4, 2])
    (None, False)
    >>> g.list_wyckoff_combinations([4, 8])
    ([[['4a'], ['8c']],
    [['4a'], ['8d']],
    [['4a'], ['8e']],
    [['4a'], ['8f']],
    [['4b'], ['8c']],
    [['4b'], ['8d']],
    [['4b'], ['8e']],
    [['4b'], ['8f']]],
    [False, True, True, True, False, True, True, True])

    or search the subgroup information

    >>> g.get_max_t_subgroup()['subgroup']
    [12, 14, 15, 20, 36, 39, 41]

    or check if a given composition is compatible with the list of wyckoff position

    >>> g = Group(225)
    >>> g.check_compatible([64, 28, 24])
    (True, True)

    or check the possible transition paths to a given supergroup

    >>> g = Group(59)
    >>> g.search_supergroup_paths(139, 2)
    [[71, 139], [129, 139], [137, 139]]



    Args:
        group: the group symbol or international number
        dim: the periodic dimension of the group
        quick: whether or not ignore the wyckoff information
    """

    def __init__(self, group, dim=3, quick=False):
        
        self.dim = dim
        names = ['Point', 'Rod', 'Layer', 'Space']
        self.header = "-- " + names[dim] + 'group --'
        self.symbol, self.number = get_symbol_and_number(group, dim)
        self.PBC, self.lattice_type = get_pbc_and_lattice(self.number, dim)
        self.alias = None

        if dim == 3:
            self.point_group, self.polar, self.inversion, self.chiral = get_point_group(self.number)

        if not quick:
            if dim == 3:
                if self.number in [5, 7, 8, 9, 12, 13, 14, 15]:
                    self.alias = self.symbol.replace("c","n")
                    #self.alias = self.alias.replace("C","I")
                self.hall_number = hall_from_hm(self.number)

            # Wyckoff positions, site_symmetry, generator
            self.wyckoffs = get_wyckoffs(self.number, dim=dim) 
            self.w_symm = get_wyckoff_symmetry(self.number, dim=dim)
            self.wyckoff_generators = get_generators(self.number, dim)

            wpdicts = [
                {
                    "index": i,
                    "letter": letter_from_index(i, self.wyckoffs, dim=self.dim),
                    "ops": self.wyckoffs[i],
                    "multiplicity": len(self.wyckoffs[i]),
                    "symmetry": self.w_symm[i],
                    "generators": self.wyckoff_generators[i],
                    "PBC": self.PBC,
                    "dim": self.dim,
                    "number": self.number,
                    "symbol": self.symbol,
                }
                for i in range(len(self.wyckoffs))
            ]

            # A list of Wyckoff_position objects, sorted by descending multiplicity
            self.Wyckoff_positions = [Wyckoff_position.from_dict(wpdict) for wpdict in wpdicts]

            # A 2D list of Wyckoff_position objects, grouped and sorted by multiplicity
            self.wyckoffs_organized = organized_wyckoffs(self)


    def __str__(self):
        try:
            return self.string
        except:
            s = self.header
            s += "# " + str(self.number) + " (" + self.symbol + ")--"

            for wp in self.Wyckoff_positions:
                ops = ss_string_from_ops(wp.symmetry[0], self.number, dim=self.dim)
                s += ("\n" + str(wp.multiplicity) + wp.letter + "\tsite symm: " + ops)
            self.string = s

            return self.string

    def __repr__(self):
        return str(self)

    def __iter__(self):
        yield from self.Wyckoff_positions

    def __getitem__(self, index):
        return self.get_wyckoff_position(index)

    def __len__(self):
        return len(self.wyckoffs)

    def get_ferroelectric_groups(self):
        """
        return the list of possible ferroelectric point groups
        """
        return para2ferro(self.point_group)

    def get_site_dof(self, sites):
        """
        compute the degree of freedom for each site
        Args:
            sites: list, e.g. ['4a', '8b'] or ['a', 'b']

        Returns:
            True or False
        """
        dof = np.zeros(len(sites))
        for i, site in enumerate(sites):
            if len(site) > 1:
                site = site[-1]
            id = len(self) - letters.index(site) - 1
            string = self[id].ops[0].as_xyz_string()
            dof[i] = len(set(re.sub('[^a-z]+', '', string)))

        return dof

    def is_valid_combination(self, sites):
        """
        check if the solutions are valid
        e.g., if a special WP with zero freedom (0,0,0) cannot be occupied twice

        Args:
            sites: list, e.g. ['4a', '8b'] or ['a', 'b']

        Returns:
            True or False
        """
        # remove the multiplicity:
        for i, site in enumerate(sites):
            if len(site) > 1:
                sites[i] = site[-1]
        
        for wp in self:
            letter = wp.letter
            if sites.count(letter)>1:
                freedom = np.trace(wp.ops[0].rotation_matrix) > 0
                if not freedom:
                    return False
        return True

    def list_wyckoff_combinations(self, numIons, quick=False):
        """
        List all possible wyckoff combinations for the given formula
        Note this is really design for a light weight calculation
        if the solution space is big, set quick as True

        Args:
            numIons: [12, 8]
            quick: Boolean, quickly generate some solutions

        Returns:
            Combinations: list of possible sites
            has_freedom: list of boolean numbers
        """

        numIons = np.array(numIons)

        basis = [] # [8, 4, 4]
        letters = [] # ['c', 'b', 'a']
        freedoms = [] # [False, False, False]

        # obtain the basis
        for i, wp in enumerate(self):
            mul = wp.multiplicity
            letter = wp.letter
            freedom = np.trace(wp.ops[0].rotation_matrix) > 0
            if mul <= max(numIons):
                if quick:
                    if mul in basis and freedom:
                        pass
                    #elif mul in basis and basis.count(mul) >= 3:
                    #    pass
                    else:
                        basis.append(mul)
                        letters.append(letter)
                        freedoms.append(freedom)
                else:
                    basis.append(mul)
                    letters.append(letter)
                    freedoms.append(freedom)

        basis = np.array(basis)

        # quickly exit
        if np.min(numIons) < np.min(basis):
            #print(numIons, basis)
            return None, False
        # odd and even
        elif np.mod(numIons, 2).sum()>0 and np.mod(basis, 2).sum()==0:
            #print("odd-even", numIons, basis)
            return None, False
        
        # obtain the maximum numbers for each basis
        max_solutions = np.floor(numIons[:,None]/basis)
        # reset the maximum to 1 if there is no freedom
        for i in range(len(freedoms)):
            if not freedoms[i]:
                max_solutions[:,i] = 1
        
        # find the integer solutions
        max_arrays = max_solutions.flatten()
        lists = []
        for a in max_arrays:
            d = int(a) + 1
            lists.append(list(range(d)))
        if len(lists) > 20:
            raise RuntimeError("Solution space is big, rerun it by setting quick as True")

        solutions = np.array(list(itertools.product(*lists)))
        big_basis = np.tile(basis, len(numIons))
        
        for i in range(len(numIons)):
            N = solutions[:,i*len(basis):(i+1)*len(basis)].dot(basis)
            solutions = solutions[N == numIons[i]]
            if len(solutions) == 0:
                #print("No solution is available")
                return None, False

        # convert the results to list
        combinations = []
        has_freedom = []
        for solution in solutions:
            res = solution.reshape([len(numIons), len(basis)])
            _com = []
            _free = []
            for i, numIon in enumerate(numIons):
                tmp = []
                bad_resolution = False
                frozen = []
                for j, b in enumerate(basis):
                    if not freedoms[j] and (res[:, j]).sum() > 1:
                        bad_resolution = True
                        break
                    else:
                        if res[i, j] > 0:
                            symbols = [str(b) + letters[j]] * res[i,j]
                            tmp.extend(symbols)
                            frozen.extend([freedoms[j]])
                if not bad_resolution:
                    _com.append(tmp)
                    _free.extend(frozen)
                       
            if len(_com) == len(numIons):
                combinations.append(_com)
                if True in _free:
                    has_freedom.append(True)
                else:
                    has_freedom.append(False)
 
        if len(combinations) == 0:
            #print("No solution is available")
            return None, False
        else:
            return combinations, has_freedom

    def get_wyckoff_position(self, index):
        """
        Returns a single Wyckoff_position object
        
        Args:
            index: the index of the Wyckoff position within the group
                The largest position is always 0

        Returns: a Wyckoff_position object
        """
        if type(index) == str:
            # Extract letter from number-letter combinations ("4d"->"d")
            for c in index:
                if c.isalpha():
                    letter = c
                    break
            index = index_from_letter(letter, self.wyckoffs, dim=self.dim)
        return self.Wyckoff_positions[index]

    def get_wyckoff_symmetry(self, index, molecular=False):
        """
        Returns the site symmetry symbol for the Wyckoff position

        Args:
            index: the index of the Wyckoff position within the group
            molecular: whether to use the Euclidean operations or not (for hexagonal groups)

        Returns: a Hermann-Mauguin style string for the site symmetry
        """
        if type(index) == str:
            # Extract letter from number-letter combinations ("4d"->"d")
            for c in index:
                if c.isalpha():
                    letter = c
                    break
            index = index_from_letter(letter, self.wyckoffs, dim=self.dim)
        if molecular:
            ops = self.w_symm_m[index][0]
        else:
            ops = self.w_symm[index][0]
        return ss_string_from_ops(ops, self.number, dim=self.dim)

    def get_alternatives(self):
        """
        Get the alternative settings as a dictionary
        """
        if self.dim == 3:
            return wyc_sets[str(self.number)]
        else:
            raise NotImplementedError("Only supports the subgroups for space group")

    def get_max_k_subgroup(self):
        """
        Returns the maximal k-subgroups as a dictionary
        """
        if self.dim == 3:
            return k_subgroup[str(self.number)]
        else:
            raise NotImplementedError("Only supports the subgroups for space group")

    def get_max_t_subgroup(self):
        """
        Returns the maximal t-subgroups as a dictionary
        """
        if self.dim == 3:
            return t_subgroup[str(self.number)]
        else:
            raise NotImplementedError("Only supports the subgroups for space group")

    def get_max_subgroup(self, H):
        if self.point_group == Group(H, quick=True).point_group:
            g_type = 'k'
            dicts = self.get_max_k_subgroup()
        else:
            g_type = 't'
            dicts = self.get_max_t_subgroup()
        return dicts, g_type

    def get_wp_list(self, reverse=False):
        """
        Get the reversed list of wps
        """
        wp_list = [(str(x.multiplicity)+x.letter) for x in self.Wyckoff_positions]
        if reverse: wp_list.reverse()
        return wp_list

    def get_splitters_from_structure(self, struc, group_type='t'):
        """
        Get the valid symmetry relations for a structure to transit to its supergroup
        e.g., 
        
        Args:
            - struc: pyxtal structure
            - group_type: `t` or `k`
        Returns:
            list of valid transitions [(id, (['4a'], ['4b'], [['4a'], ['4c']])]
        """

        if group_type=='t':
            dicts = self.get_max_t_subgroup()
        else:
            dicts = self.get_max_k_subgroup()

        # search for the compatible solutions
        solutions = []
        for i, sub in enumerate(dicts['subgroup']):
            if sub == struc.group.number:
                # extract the relation
                relation = dicts['relations'][i]
                trans = np.linalg.inv(dicts['transformation'][i][:,:3])
                if struc.lattice.check_mismatch(trans, self.lattice_type):
                    results = self.get_splitters_from_relation(struc, relation)
                    if results is not None:
                        sols = list(itertools.product(*results))
                        trials = self.get_valid_solutions(sols)
                        solutions.append((i, trials))
        return solutions

    def get_splitters_from_relation(self, struc, relation):        
        """
        Get the valid symmetry relations for a structure to transit to its supergroup
        e.g., 
        
        Args:
            - struc: pyxtal structure
            - group_type: `t` or `k`
        Returns:
            list of valid transitions 
        """

        elements, sites = struc._get_elements_and_sites()
        wp_list = self.get_wp_list(reverse=True)

        good_splittings_list=[]

        # search for all valid compatible relations
        # each element is solved one at a time
        for site in sites:
            # ['4a', '4a', '2b'] -> ['4a', '2b']
            _site = np.unique(site)
            _site_counts = [site.count(x) for x in _site]
            wp_indices = []

            # the sum of all positions should be fixed.
            total_units = 0
            for j, x in enumerate(_site):
                total_units += int(x[:-1])*_site_counts[j]

            # collect all possible supergroup transitions
            for j, split in enumerate(relation):
                if np.all([x in site for x in split]):
                    wp_indices.append(j)

            wps = [wp_list[x] for x in wp_indices]
            blocks = [np.array([relation[j].count(s) for s in _site]) for j in wp_indices]
            block_units = [sum([int(x[:-1])*block[j] for j, x in enumerate(_site)]) for block in blocks]

            # below is a brute force search for the valid combinations
            combo_storage = [np.zeros(len(block_units))]
            good_list = []
            while len(combo_storage) > 0:
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
                            tester = np.zeros(len(_site_counts))
                            for l, z in enumerate(trial):
                                tester += z*blocks[l]
                            if np.all(tester == _site_counts):
                                G_sites = []
                                for l, number in enumerate(trial):
                                    if number == 0:
                                        continue
                                    elif number == 1:
                                        G_sites.append(wps[l])
                                    else:
                                        for i in range(int(number)):
                                            G_sites.append(wps[l])
                                if G_sites not in good_list:
                                    good_list.append(G_sites)
                combo_storage=holder

            if len(good_list) == 0:
                return None
            else:
                good_splittings_list.append(good_list)
        return good_splittings_list


    def get_min_supergroup(self, group_type='t', G=None):
        """
        Returns the minimal supergroups as a dictionary
        """
        if self.dim == 3:
            dicts = {'supergroup': [],
                     'transformation': [],
                     'relations': [],
                     'idx': [],
                    }
            if G is None:
                sgs = range(1,231)
            else:
                sgs = G

            for sg in sgs:
                subgroups = None
                if group_type == 't':
                    if sg>self.number:
                        subgroups = Group(sg, quick=True).get_max_t_subgroup()
                else:
                    g1 = Group(sg)
                    if g1.point_group == self.point_group:
                        subgroups = Group(sg, quick=True).get_max_k_subgroup()
                if subgroups is not None:
                    for i, sub in enumerate(subgroups['subgroup']): 
                        if sub == self.number:
                            trans = subgroups['transformation'][i]
                            relation = subgroups['relations'][i]
                            dicts['supergroup'].append(sg)
                            dicts['transformation'].append(trans)
                            dicts['relations'].append(relation)
                            dicts['idx'].append(i)
            return dicts
        else:
            raise NotImplementedError("Only supports the supergroups for space group")

    def get_max_subgroup_numbers(self):
        """
        Returns the minimal supergroups as a dictionary
        """
        groups = []
        if self.dim == 3:
            k = k_subgroup[str(self.number)]['subgroup']
            t = t_subgroup[str(self.number)]['subgroup']
            return k+t
        else:
            raise NotImplementedError("Only supports the supergroups for space group")


    def get_lists(self, numIon, used_indices):
        """
        Compute the lists of possible mult/maxn/freedom/i_wp

        Args:
            numIon: integer number of atoms
            used_indices: a list of integer numbers
        """
        l_mult0 = []
        l_maxn0 = []
        l_free0 = []
        indices0 = []
        for i_wp, wp in enumerate(self):
            indices0.append(i_wp)
            l_mult0.append(len(wp))
            l_maxn0.append(numIon // len(wp))
            #check the freedom
            if np.allclose(wp[0].rotation_matrix, np.zeros([3, 3])):
                l_free0.append(False)
            else:
                l_free0.append(True)
        return self.clean_lists(numIon, l_mult0, l_maxn0, l_free0, indices0, used_indices)

    def get_lists_mol(self, numIon, used_indices, orientations):
        """
        Compute the lists of possible mult/maxn/freedom/i_wp

        Args:
            numIon: integer number of atoms
            used_indices: a list of integer numbers
            orientations: list of orientations
        """
        l_mult0 = []
        l_maxn0 = []
        l_free0 = []
        indices0 = []
        for i_wp, wp in enumerate(self):
            # Check that at least one valid orientation exists
            j, k = jk_from_i(i_wp, self.wyckoffs_organized)
            if len(orientations) > j and len(orientations[j]) > k:
                indices0.append(i_wp)
                l_mult0.append(len(wp))
                l_maxn0.append(numIon // len(wp))
                #check the freedom
                if np.allclose(wp[0].rotation_matrix, np.zeros([3, 3])):
                    l_free0.append(False)
                else:
                    l_free0.append(True)
        return self.clean_lists(numIon, l_mult0, l_maxn0, l_free0, indices0, used_indices)

    @staticmethod
    def clean_lists(numIon, l_mult0, l_maxn0, l_free0, indices0, used_indices):
        # Remove redundant multiplicities:
        l_mult = []
        l_maxn = []
        l_free = []
        indices = []
        for mult, maxn, free, i_wp in zip(l_mult0, l_maxn0, l_free0, indices0):
            if free:
                if mult not in l_mult:
                    l_mult.append(mult)
                    l_maxn.append(maxn)
                    l_free.append(True)
                    indices.append(i_wp)
            elif not free and i_wp not in used_indices:
                l_mult.append(mult)
                indices.append(i_wp)
                if mult <= numIon:
                    l_maxn.append(1)
                elif mult > numIon:
                    l_maxn.append(0)
                l_free.append(False)
        return l_mult, l_maxn, l_free, indices


    def check_compatible(self, numIons, valid_orientations=None):
        """
        Checks if the number of atoms is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.

        Args: 
            numIons: list of integers
            valid_orientations: list of possible orientations (molecule only)

        Returns:
            Compatible: True/False
            has_freedom: True/False
        """
        has_freedom = False #whether or not at least one degree of freedom exists
        used_indices = [] #wp's already used that don't have any freedom

        # Loop over species
        # Sort the specie from low to high so that the solution can be found ealier
        for id in np.argsort(numIons):
            # Get lists of multiplicity, maxn and freedom
            numIon = numIons[id]
            if valid_orientations is None:
                l_mult, l_maxn, l_free, indices = self.get_lists(numIon, used_indices)
            else:
                vo = valid_orientations[id]
                l_mult, l_maxn, l_free, indices = self.get_lists_mol(numIon, used_indices, vo)
            #print(numIon, l_mult, indices, l_maxn, l_free)

            # Loop over possible combinations
            p = 0  # Create pointer variable to move through lists
            # Store the number of each WP, used across possible WP combinations
            n0 = [0] * len(l_mult)
            n = deepcopy(n0)
            for i, mult in enumerate(l_mult):
                if l_maxn[i] != 0:
                    p = i
                    n[i] = l_maxn[i]
                    break
            p2 = p
            if n == n0:
                return False, False

            #print(numIon, n, n0, p)
            while True:
                num = np.dot(n, l_mult)
                dobackwards = False
                # The combination works: move to next species
                if num == numIon:
                    # Check if at least one degree of freedom exists
                    for val, free, i_wp in zip(n, l_free, indices):
                        if val > 0:
                            if free:
                                has_freedom = True
                            else:
                                used_indices.append(i_wp)
                    break
                # All combinations failed: return False
                if n == n0 and p >= len(l_mult) - 1:
                    #print('All combinations failed', numIon, n, n0)
                    return False, False
                # Too few atoms
                if num < numIon:
                    # Forwards routine
                    # Move p to the right and max out
                    if p < len(l_mult) - 1:
                        p += 1
                        n[p] = min((numIon - num) // l_mult[p], l_maxn[p])
                    elif p == len(l_mult) - 1:
                        # p is already at last position: trigger backwards routine
                        dobackwards = True
                # Too many atoms
                if num > numIon or dobackwards:
                    # Backwards routine
                    # Set n[p] to 0, move p backwards to non-zero, and decrease by 1
                    n[p] = 0
                    while p > 0 and p > p2:
                        p -= 1
                        if n[p] != 0:
                            n[p] -= 1
                            if n[p] == 0 and p == p2:
                                p2 = p + 1
                            break
            #print('used_indices', used_indices)
        if has_freedom:
            # All species passed: return True
            return True, True
        else:
            # All species passed, but no degrees of freedom: return 0
            return True, False

    def search_supergroup_paths(self, H, max_layer=5):
        """
        Search paths to transit to super group H. if
        - path1 is a>>e
        - path2 is a>>b>>c>>e
        path 2 will not be counted since path 1 exists
    
        Args:
            H: final supergroup number
            max_layer: the number of supergroup calculations needed.
    
        Returns:
            list of possible paths ordered from G to H
        """
    
        layers = {}
        layers[0] = {'groups': [H], 
                     'subgroups': [],
                   }
        final = []
        traversed = []
    
        # Searches for every subgroup of the the groups from the previous layer.
        # stores the possible groups of each layer and their subgroups
        # in a dictinoary to avoid redundant calculations.
        for l in range(1, max_layer+1):
            previous_layer_groups=layers[l-1]['groups']
            groups = []
            subgroups = []
            for g in previous_layer_groups:
                subgroup_numbers=np.unique(Group(g, quick=True).get_max_subgroup_numbers())
    
                # If a subgroup list has been found with H
                # trace a path through the dictionary to build the path
                if self.number in subgroup_numbers:
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
    
                # Continue to generate next layer if the path to H has not been found.
                else:
                    subgroups.append(subgroup_numbers)
                    for x in subgroup_numbers:
                        if (x not in groups) and (x not in traversed):
                            groups.append(x)
    
            traversed.extend(groups)
            layers[l] = {'groups': deepcopy(groups), 
                         'subgroups':[]}
            layers[l-1]['subgroups'] = deepcopy(subgroups)
        return final

    def search_subgroup_paths(self, G, max_layer=5):
        """
        Search paths to transit to subgroup H. if
        - path1 is a>>e
        - path2 is a>>b>>c>>e
        path 2 will not be counted since path 1 exists
    
        Args:
            G: final subgroup number
            max_layer: the number of supergroup calculations needed.
    
        Returns:
            list of possible paths ordered from H to G
        """
        paths = Group(G, quick=True).search_supergroup_paths(self.number, max_layer=max_layer)
        for p in paths:
            p.reverse()
            p.append(G)
        return paths

    def get_valid_solutions(self, solutions):
        """
        check if the solutions are valid
        a special WP such as (0,0,0) cannot be occupied twice

        Args:
            solutions: list of solutions about the distibution of WP sites

        Returns:
            the filtered solutions that are vaild
        """
        valid_solutions = []
        for solution in solutions:
            sites = []
            for s in solution:
                sites.extend(s)
            if self.is_valid_combination(sites):
                valid_solutions.append(solution)
        return valid_solutions

    def cellsize(self):
        """
        Returns the number of duplicate atoms in the conventional lattice (in
        contrast to the primitive cell). Based on the type of cell centering (P,
        A, C, I, R, or F)
        """
        if self.dim in [0, 1]:
            # Rod and point groups
            return 1
        elif self.dim == 2:
            # Layer groups
            if self.number in [10, 13, 18, 22, 26, 35, 36, 47, 48]:
                return 2
            else:
                return 1
        else:
            # space groups
            if self.symbol[0] == 'P':
                return 1 # P
            elif self.symbol[0] == 'R':
                return 3 # R
            elif self.symbol[0] == 'F':
                return 4 # F
            else:
                return 2  # A, C, I

    @classmethod
    def list_groups(cls, dim=3):
        """
        Function for quick print of groups and symbols
    
        Args:
            group: the group symbol or international number
            dim: the periodic dimension of the group
        """
    
        import pandas as pd
    
        keys = {
            3: "space_group",
            2: "layer_group",
            1: "rod_group",
            0: "point_group",
        }
        data = symbols[keys[dim]]
        df = pd.DataFrame(index=range(1, len(data) + 1), data=data, columns=[keys[dim]])
        pd.set_option("display.max_rows", len(df))
        print(df)

# --------------------------- Wyckoff Position class  -----------------------------

class Wyckoff_position:
    """
    Class for a single Wyckoff position within a symmetry group

    Examples
    --------
    >>> from pyxtal.symmetry import Wyckoff_position as wp
    >>> wp.from_group_and_index(19, 0)
    Wyckoff position 4a in space group 19 with site symmetry 1
    x, y, z
    -x+1/2, -y, z+1/2
    -x, y+1/2, -z+1/2
    x+1/2, -y+1/2, -z

    """

    def from_dict(dictionary):
        """
        Constructs a Wyckoff_position object using a dictionary. Used mainly by the
        Wyckoff class for constructing a list of Wyckoff_position objects at once
        """
        wp = Wyckoff_position()
        for key in dictionary:
            setattr(wp, key, dictionary[key])
        #wp.get_site_symmetry()
        wp.set_euclidean()
        return wp

    def get_dof(self):
        """
        Simply return the degree of freedom
        """
        return np.linalg.matrix_rank(self.ops[0].rotation_matrix)

    def get_frozen_axis(self):
        if self.index == 0:
            return []
        elif self.get_dof() == 0:
            return [0, 1, 2]
        else:
            if self.number >=75:
                if self.ops[0].rotation_matrix[2,2] == 1:
                    return [0, 1]
                else:
                    return [0, 1, 2]
            else:
                if self.get_dof() == 1:
                    if self.ops[0].rotation_matrix[2,2] == 1:
                        return [0, 1]
                    elif self.ops[0].rotation_matrix[1,1] == 1:
                        return [0, 2]
                    elif self.ops[0].rotation_matrix[0,0] == 1:
                        return [1, 2]
                else:
                    if self.ops[0].rotation_matrix[2,2] != 1:
                        return [2]
                    elif self.ops[0].rotation_matrix[1,1] != 1:
                        return [1]
                    elif self.ops[0].rotation_matrix[0,0] != 1:
                        return [0]

    def __str__(self, supress=False):
        if self.dim not in list(range(4)):
            return "invalid crystal dimension. Must be a number between 0 and 3."
        if not hasattr(self, "site_symm"): self.get_site_symmetry()   
        s = "Wyckoff position " + str(self.multiplicity) + self.letter + " in "
        if self.dim == 3:
            s += "space "
        elif self.dim == 2:
            s += "layer "
        elif self.dim == 1:
            s += "Rod "
        elif self.dim == 0:
            s += "Point group " + self.symbol
        if self.dim != 0:
            s += "group " + str(self.number)
        s += " with site symmetry " + self.site_symm

        if not supress:
            for op in self.ops:
                s += "\n" + op.as_xyz_string()

        self.string = s
        return self.string

    def __repr__(self):
        return str(self)

    def diagonalize_symops(self):
        """
        Obtain the symmetry in n representation for P2/c, Cc, P21/c, Pc, C2/c
        """
        ops = Group(self.number)[self.index]
        #Pn, Cn, P2/n, P21/n, C2/n
        if self.number in [7, 9, 13, 14, 15]:
            trans = np.array([[1,0,0],[0,1,0],[1,0,1]])
            for j, op in enumerate(ops):
                vec = op.translation_vector.dot(trans)
                vec -= np.floor(vec) 
                op1 = op.from_rotation_and_translation(op.rotation_matrix, vec)
                self.ops[j] = op1    
        #I2, Im, Ic, I2/m, I2/c
        elif self.number in [5, 8, 9, 12, 15]:
            trans = np.array([[1,0,1],[0,1,0],[1,0,1]])
            for j, op in enumerate(ops):
                vec = op.translation_vector.dot(trans)
                vec -= np.floor(vec) 
                op1 = op.from_rotation_and_translation(op.rotation_matrix, vec)
                self.ops[j] = op1  
        #In, I2/n

    def get_site_symm_wo_translation(self):
        ops = []
        for op in self.symmetry[0]:
            op = SymmOp.from_rotation_and_translation(op.rotation_matrix, [0, 0, 0])
            ops.append(op)
        return ops

    def equivalent_set(self, index):
        """
        Transform the wp to another equivalent set.
        Needs to update both wp and positions

        Args:
            transformation: index
        """
        if self.index > 0:
            G = Group(self.number)
            if len(G[index]) != len(G[self.index]):
                msg = "Spg {:d}, Invalid switch in Wyckoff Pos\n".format(self.number)
                msg += str(self)
                msg += "\n"+str(G[index])
                raise ValueError(msg)
            else:
                return G[index]
        return self

    def is_pure_translation(self, id):
        """
        Check if the operation is equivalent to pure translation
        """
        op = self.generators[id]
        diff = op.rotation_matrix - np.eye(3)
        if np.sum(diff.flatten()**2) < 1e-4:
            return True
        else:
            ops = self.get_site_symm_wo_translation()
            return (op in ops)

    def swap_axis(self, swap_id):
        """
        swap the axis may result in a new wp
        """
        if self.index > 0:
            perm_id = None
            _ops = [self.ops[0]]
            trans = [np.zeros(3)]
            if self.symbol[0] == "F":
                trans.append(np.array([0,0.5,0.5]))
                trans.append(np.array([0.5,0,0.5]))
                trans.append(np.array([0.5,0.5,0]))
            elif self.symbol[0] == "I":
                trans.append(np.array([0.5,0.5,0.5]))
            elif self.symbol[0] == "A":
                trans.append(np.array([0,0.5,0.5]))
            elif self.symbol[0] == "B":
                trans.append(np.array([0.5,0,0.5]))
            elif self.symbol[0] == "C":
                trans.append(np.array([0.5,0.5,0]))

            op_perm = swap_xyz_ops(_ops, swap_id)[0]
            for id, ops in enumerate(Group(self.number)):
                if len(ops) == len(self.ops):
                    for i, tran in enumerate(trans):
                        if i > 0:
                            # apply tran
                            op = op_translation(op_perm, tran)
                        else:
                            op = op_perm
                        #print(id, op.as_xyz_string(),tran)
                        if are_equivalent_ops(op, ops[0]):
                            perm_id = id
                            return Group(self.number)[id], tran
            if perm_id is None:
                raise ValueError("cannot swap", swap_id, self)
        return self, np.zeros(3)

    def from_symops(ops, group=None, permutation=True):
        """
        search Wyckoff Position by symmetry operations
        Now only supports space group symmetry

        Args:
            ops: a list of symmetry operations
            group: the space group number or Group object

        Returns:
            `Wyckoff_position`

        """
        if isinstance(ops[0], str):
            str1 = ops
        else:
            str1 = [op.as_xyz_string() for op in ops]

        str1 = [st.replace("-1/2","+1/2") for st in str1]

        N_sym = len(str1)
        # sometimes, we allow the permutation
        if permutation:
            permutations = [[0,1,2],[1,0,2],[2,1,0],[0,2,1]]
        else:
            permutations = [[0,1,2]]

        if group is None:
            groups = range(1,231)
        else:
            groups = [group]

        for i in groups:
            if len(groups) == 1:
                G_ops = i
            else:
                G_ops = Group(i)
            for wyc in G_ops:
                if len(wyc) == N_sym:
                    wyc.PBC = [1, 1, 1]
                    wyc.dim = 3
                    # Try permutation first
                    str2 = [op.as_xyz_string() for op in wyc.ops]
                    for perm in permutations:
                        str_perm = swap_xyz_string(str1, perm)
                        # Compare the pure rotation and then 
                        if set(str_perm) == set(str2):
                            return wyc, perm

                    # Try monoclinic space groups (P21/n, Pn, C2/n) 
                    if i in [7, 9, 13, 14, 15]:
                        trans = np.array([[1,0,0],[0,1,0],[1,0,1]])
                        str3 = []
                        for j, op in enumerate(wyc.ops):
                            vec = op.translation_vector.dot(trans)
                            vec -= np.floor(vec) 
                            op3 = op.from_rotation_and_translation(op.rotation_matrix, vec)
                            str3.append(op3.as_xyz_string().replace("-1/2","+1/2"))
                            #wyc.ops[j] = op3
                        if set(str3) == set(str1):
                            return wyc, trans
                    elif i in [5, 8, 9, 12, 15]:
                        trans = np.array([[1,0,1],[0,1,0],[1,0,1]])
                        str3 = []
                        for j, op in enumerate(wyc.ops):
                            vec = op.translation_vector.dot(trans)
                            vec -= np.floor(vec) 
                            op3 = op.from_rotation_and_translation(op.rotation_matrix, vec)
                            str3.append(op3.as_xyz_string().replace("-1/2","+1/2"))
                            #wyc.ops[j] = op3
                        if set(str3) == set(str1):
                            return wyc, trans
                elif len(wyc) < N_sym:
                    continue #break

        return None, None


    def from_group_and_index(group, index, dim=3, PBC=None):
        """
        Creates a Wyckoff_position using the space group number and index
        
        Args:
            group: the international number of the symmetry group
            index: the index or letter of the Wyckoff position within the group.
                0 is always the general position, and larger indeces represent positions
                with lower multiplicity. Alternatively, index can be the Wyckoff letter
                ("4a" or "6f")
            dim: the periodic dimension of the crystal
            PBC: the periodic boundary conditions
        """
        wp = Wyckoff_position()
        wp.dim = dim
        wp.number = group
        wp.PBC = PBC
        use_letter = False

        if np.issubdtype(type(index), np.integer):
            wp.index = index
        elif type(index) == str:
            use_letter = True
            # Extract letter from number-letter combinations ("4d"->"d")
            for c in index:
                if c.isalpha():
                    index = c
                    break
        if dim > 0:
            if dim == 3:
                if wp.PBC is None:
                    wp.PBC = [1, 1, 1]
            elif dim == 2:
                if wp.PBC is None:
                    wp.PBC = [1, 1, 0]
            elif dim == 1:
                if wp.PBC is None:
                    wp.PBC = [0, 0, 1]

            ops_all = get_wyckoffs(wp.number, dim=dim)
            if use_letter:
                wp.index = index_from_letter(index, ops_all, dim=dim)
                wp.letter = index
            else:
                wp.letter = letter_from_index(wp.index, ops_all, dim=dim)
            wp.ops = ops_all[wp.index]

            wp.multiplicity = len(wp.ops)
            wp.symmetry = get_wyckoff_symmetry(wp.number, dim)[wp.index]
            wp.generators = get_generators(wp.number, dim)[wp.index]

        elif dim == 0:
            # Generate a Group and retrieve Wyckoff position from it
            g = Group(group, dim=0)
            wp = g[index]
        #wp.generators = get_generators(wp.number, dim)[wp.index]
        #wp.get_site_symmetry()
        wp.set_euclidean()
        wp.symbol, _ = get_symbol_and_number(wp.number, wp.dim)
        return wp

    def wyckoff_from_generating_op(gen_op, gen_pos):
        """
        Given a general position and generating operation (ex: "x,0,0"), returns a
        Wyckoff_position object.
        
        Args:
            gen_op: a SymmOp into which the generating coordinate will be plugged
            gen_pos: a list of SymmOps representing the general position

        Returns:
            a list of SymmOps
        """

        def is_valid_generator(op):
            m = op.affine_matrix
            # Check for x,y+ax,z+bx+cy format
            # Make sure y, z are not referred to before second and third slots
            if (np.array([m[0][1], m[0][2], m[1][2]]) != 0.0).any():
                if (np.array([m[2][0], m[2][1], m[1][0]]) != 0.0).any():
                    return False
            for i in range(0, 3):
                # Check for translation
                if m[i][i] != 0 and m[i][3] != 0:
                    return False
                # Check that x, y, z are present in desired spot and positive
                if np.linalg.norm(m[:, i]) > 0:
                    if m[i][i] == 0:
                        return False
            return True

        def filter_zeroes(op):
            m = op.affine_matrix
            m2 = m
            for i, x in enumerate(m):
                for j, y in enumerate(x):
                    if np.isclose(y, 0, atol=1e-3):
                        m2[i][j] = 0
            return SymmOp(m2)

        Identity = SymmOp.from_xyz_string("x,y,z")
        # new_ops = [filter_zeroes(op*gen_op) for op in gen_pos]
        new_ops = [
            filter_zeroes(SymmOp(np.dot(op.affine_matrix, gen_op.affine_matrix)))
            for op in gen_pos
        ]

        # Remove redundant ops
        op_list = []
        for op1 in new_ops:
            passed = True
            for op2 in op_list:
                if np.allclose(op1.affine_matrix, op2.affine_matrix):
                    passed = False
                    break
            if passed is True:
                op_list.append(op1)

        # Check for identity op
        if Identity in op_list:
            index = op_list.index(Identity)
            op2 = op_list[0]
            op_list[index] = op2
            op_list[0] = Identity
            return op_list

        # Check for generating op
        found = False
        for op in op_list:
            if is_valid_generator(op):
                # Place generating op at beginning of list
                index = op_list.index(op)
                op2 = op_list[0]
                op_list[index] = op2
                op_list[0] = op
                found = True
                break
        if found is False:
            printx("Error: Could not find generating op.", priority=0)
            printx("gen_op: ")
            printx(str(gen_op))
            printx("gen_pos: ")
            printx(str(gen_pos))
            return op_list

        # Make sure first op is normalized
        m = op_list[0].affine_matrix
        x_factor = 1
        y_factor = 1
        z_factor = 1
        if m[0][0] != 0:
            if m[0][0] != 1:
                x_factor = m[0][0]
        if m[1][1] != 0:
            if m[1][1] != 1:
                y_factor = m[1][1]
        if m[2][2] != 0:
            if m[2][2] != 1:
                z_factor = m[2][2]
        new_list = []
        for op in op_list:
            m = op.affine_matrix
            m2 = m
            m2[:, 0] = m2[:, 0] / x_factor
            m2[:, 1] = m2[:, 1] / y_factor
            m2[:, 2] = m2[:, 2] / z_factor
            new_list.append(SymmOp(m2))

        return new_list

    def symmetry_from_wyckoff(wp, gen_pos):
        symm = []
        for op in wp:
            symm.append(site_symm(op, gen_pos))
        return symm

    def __iter__(self):
        yield from self.ops

    def __getitem__(self, index):
        return self.ops[index]

    def __len__(self):
        return self.multiplicity

    def get_site_symmetry(self):
        if self.euclidean:
            ops = self.get_euclidean_symmetries()
        else:
            ops = self.symmetry[0]
        self.site_symm = ss_string_from_ops(ops, self.number, dim=self.dim)


    def get_euclidean_symmetries(self):
        """
        return the symmetry operation object at the Euclidean space 

        Returns:
            list of pymatgen SymmOp object
        """
        ops = []
        for op in self.symmetry[0]:
            hat = SymmOp.from_rotation_and_translation(hex_cell, [0, 0, 0])
            ops.append(hat * op * hat.inverse)
        return ops

    def print_ops(self, ops=None):
        if ops is None:
            ops = self.ops
        for op in ops:
            print(op.as_xyz_string())

    def gen_pos(self):
        """
        Returns the general Wyckoff position
        """
        return self.Wyckoff_positions[0]

    def is_equivalent(self, pt1, pt2, cell=np.eye(3), tol=0.05):
        """
        Check two pts are equivalent
        """
        pt1, dist1 = self.search_generator(pt1, cell)
        pt2, dist2 = self.search_generator(pt2, cell)
        if dist1 > tol or dist2 > tol:
            return False
        else:
            pt1 = np.array(pt1); pt1 -= np.floor(pt1)
            pt2 = np.array(pt2); pt2 -= np.floor(pt2)
            pts = self.apply_ops(pt1); pts -= np.floor(pts)
            diffs = pt2 - pts
            diffs -= np.round(diffs)
            diffs = np.dot(diffs, cell)
            dists = np.linalg.norm(diffs, axis=1)
            #print(dists)
            if len(dists[dists<tol]) > 0:
                return True
            else:
                return False
        
    def merge(self, pt, lattice, tol, orientations=None):
        """
        Given a list of fractional coordinates, merges them within a given
        tolerance, and checks if the merged coordinates satisfy a Wyckoff
        position. 
    
        Args:
            pt: the originl point (3-vector)
            lattice: a 3x3 matrix representing the unit cell
            tol: the cutoff distance for merging coordinates
            orientations: the valid orientations for a given molecule. 
    
        Returns:
            pt: 3-vector after merge
            wp: a `pyxtal.symmetry.Wyckoff_position` object, If no match, returns False. 
            valid_ori: the valid orientations after merge
    
        """
        wp = deepcopy(self)
        PBC = wp.PBC
        group = Group(wp.number, wp.dim)
        pt = self.project(pt, lattice, PBC)
        coor = apply_ops(pt, wp)
        if orientations is None:
            valid_ori = None
        else:
            j, k = jk_from_i(wp.index, orientations)
            valid_ori = orientations[j][k]
        
        # Main loop for merging multiple times
        while True:
            # Check distances of current WP. If too small, merge
            dm = distance_matrix([coor[0]], coor, lattice, PBC=PBC)
            passed_distance_check = True
            x = np.argwhere(dm < tol)
            for y in x:
                # Ignore distance from atom to itself
                if y[0] == 0 and y[1] == 0:
                    pass
                else:
                    passed_distance_check = False
                    break
    
            # for molecular crystal, one more check
            if not check_images([coor[0]], [6], lattice, PBC=PBC, tol=tol):
                passed_distance_check = False
    
            if not passed_distance_check:
                mult1 = wp.multiplicity
                # Find possible wp's to merge into
                possible = []
                for i, wp0 in enumerate(group):
                    mult2 = wp0.multiplicity
                    # Check that a valid orientation exists
                    if orientations is not None:
                        res = jk_from_i(i, orientations)
                        if res is None:
                            continue
                        else:
                            j, k = res 
                            if orientations[j][k] == []:
                                continue
                            else:
                                valid_ori = orientations[j][k]
                    # factor = mult2 / mult1
    
                    if (mult2 < mult1) and (mult1 % mult2 == 0):
                        possible.append(i)
                if possible == []:
                    return None, False, valid_ori
                # Calculate minimum separation for each WP
                distances = []
                pts = []
                for i in possible:
                    wp = group[i]
                    p, d = wp.search_generator(pt.copy(), lattice)
                    distances.append(d)
                    pts.append(p)

                # Choose wp with shortest translation for generating point
                tmpindex = np.argmin(distances)
                index = possible[tmpindex]
                wp = group[index]
                pt = pts[tmpindex]
                coor = wp.apply_ops(pt)
            # Distances were not too small; return True
            else:
                return pt, wp, valid_ori

    def get_euclidean_generator(self, cell, idx=0):
        """
        return the symmetry operation object at the Euclidean space 

        Args:
            cell: 3*3 cell matrix
            idx: the index of wp generator

        Returns:
            pymatgen SymmOp object
        """
        op = self.generators[idx]
        if self.euclidean:
            hat = SymmOp.from_rotation_and_translation(cell.T, [0, 0, 0])
            op = hat * op * hat.inverse

        return op


    def set_euclidean(self):
        convert = False
        if self.dim == 3:
            if 143 <= self.number < 195:
                convert = True
        elif self.dim == 2:
            if self.number >= 65:
                convert = True
        elif self.dim == 1:
            if self.number >= 42:
                convert = True
        self.euclidean = convert

    def apply_ops(self, pt):
        """
        apply symmetry operation
        """
        return apply_ops(pt, self.ops)

    def search_generator(self, pt, lattice=np.eye(3)):
        """
        For a given special wp, (e.g., [(x, 0, 1/4), (0, x, 1/4)]),
        return the first position according to the symmetry operation
        
        Args:
            pt: 1*3 vector
            lattice: 3*3 matrix

        Returns:
            pt: the best matched pt
            diff: numerical difference

        """
        if self.index == 0: #general sites
            return pt, 0

        else:
            d = []

            if self.get_dof == 0: #fixed site like [0, 0, 0]
                pts = self.apply_ops(pt)
                for p0 in pts:
                    d.append(distance(p0, lattice, PBC=self.PBC))

            else: # sites like (x, 0, 0)
                wp0 = Group(self.number)[0]
                pts = wp0.apply_ops(pt)
                op = self.ops[0]
                for i, p0 in enumerate(pts):
                    coord = op.operate(p0)-p0
                    d.append(distance(coord, lattice, PBC=self.PBC))
            #print(d)
            d = np.array(d)
            return pts[np.argmin(d)], np.min(d)


    def project(self, point, cell=np.eye(3), PBC=[1, 1, 1], id=0):
        """
        Given a 3-vector and a Wyckoff position operator, 
        returns the projection onto the axis, plane, or point.
    
        >>> wp.project_point([0,0.3,0.1], 
        array([0. , 0.3, 0.1])
    
        Args:
            point: a 3-vector (numeric list, tuple, or array)
            cell: 3x3 matrix describing the unit cell vectors
            PBC: A periodic boundary condition list, 
                where 1 means periodic, 0 means not periodic.
                Ex: [1,1,1] -> full 3d periodicity, 
                    [0,0,1] -> periodicity along the z axis
    
        Returns:
            a transformed 3-vector (numpy array)
        """
        op = self.ops[id]
        rot = op.rotation_matrix
        trans = op.translation_vector
        point = np.array(point, dtype=float)
    
        def project_single(point, rot, trans):
            # move the point in the opposite direction of the translation
            point -= trans
            new_vector = np.zeros(3)
            # Loop over basis vectors of the symmetry element
            for basis_vector in rot.T:
                # b = np.linalg.norm(basis_vector)
                b = np.sqrt(basis_vector.dot(basis_vector))  # a faster version?
                #if not np.isclose(b, 0):
                if b > 1e-3:
                    new_vector += basis_vector * (np.dot(point, basis_vector) / (b ** 2))
            new_vector += trans
            return new_vector

        if PBC == [0, 0, 0]:
            return project_single(point, rot, trans)
        else:
            pt = filtered_coords(point)
            m = create_matrix(PBC=PBC)
            new_vectors = []
            distances = []
            for v in m:
                new_vector = project_single(pt, rot, trans+v)
                new_vectors.append(new_vector)
                tmp = (new_vector-point).dot(cell)
                distances.append(np.linalg.norm(tmp))
            i = np.argmin(distances)
            return filtered_coords(new_vectors[i], PBC=PBC)

    def get_all_positions(self, pos):
        """
        return the list of position from any single coordinate from wp.
        The position does not have to be the 1st number in the wp list

        >>> from pyxtal.symmetry import Group
        >>> wp2 = Group(62)[-1]
        >>> wp2
        Wyckoff position 4a in space group 62 with site symmetry -1
        0, 0, 0
        1/2, 0, 1/2
        0, 1/2, 0
        1/2, 1/2, 1/2
        >>> wp2.get_all_positions([1/2, 1/2, 1/2])
        array([[0. , 0. , 0. ],
               [0.5, 0. , 0.5],
               [0. , 0.5, 0. ],
               [0.5, 0.5, 0.5]])
        """

        pos0, _, _ = self.merge(pos, np.eye(3), 0.01)
        res = self.apply_ops(pos0)
        res -= np.floor(res)
        return res

# --------------------------- Wyckoff Position selection  -----------------------------

def choose_wyckoff(group, number=None, site=None, dim=3):
    """
    Choose a Wyckoff position to fill based on the current number of atoms
    needed to be placed within a unit cell
    Rules:
        0) use the pre-assigned list if this is provided
        1) The new position's multiplicity is equal/less than (number).
        2) We prefer positions with large multiplicity.

    Args:
        group: a pyxtal.symmetry.Group object
        number: the number of atoms still needed in the unit cell
        site: the pre-assigned Wyckoff sites (e.g., 4a)

    Returns:
        Wyckoff position. If no position is found, returns False
    """

    if site is not None:
        return Wyckoff_position.from_group_and_index(group.number, site, dim)
    else:
        wyckoffs_organized = group.wyckoffs_organized

        if random.uniform(0, 1) > 0.5:  # choose from high to low
            for wyckoff in wyckoffs_organized:
                if len(wyckoff[0]) <= number:
                    return random.choice(wyckoff)
            return False
        else:
            good_wyckoff = []
            for wyckoff in wyckoffs_organized:
                if len(wyckoff[0]) <= number:
                    for w in wyckoff:
                        good_wyckoff.append(w)
            if len(good_wyckoff) > 0:
                return random.choice(good_wyckoff)
            else:
                return False

def choose_wyckoff_molecular(group, number, site, orientations, general_site=True, dim=3):
    """
    Choose a Wyckoff position to fill based on the current number of molecules
    needed to be placed within a unit cell

    Rules:

        1) The new position's multiplicity is equal/less than (number).
        2) We prefer positions with large multiplicity.
        3) The site must admit valid orientations for the desired molecule.

    Args:
        group: a pyxtal.symmetry.Group object
        number: the number of molecules still needed in the unit cell
        orientations: the valid orientations for a given molecule. Obtained
            from get_sg_orientations, which is called within molecular_crystal

    Returns:
        Wyckoff position. If no position is found, returns False
    """
    wyckoffs = group.wyckoffs_organized

    if site is not None:
        return Wyckoff_position.from_group_and_index(group.number, site, dim)

    elif general_site or np.random.random() > 0.5:  # choose from high to low
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                good_wyckoff = []
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append(w)
                if len(good_wyckoff) > 0:
                    return random.choice(good_wyckoff)
        return False
    else:
        good_wyckoff = []
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append(w)
        if len(good_wyckoff) > 0:
            return random.choice(good_wyckoff)
        else:
            return False

# -------------------- quick utilities for symmetry conversion ----------------
def swap_xyz_string(xyzs, permutation):
    """
    Permutate the xyz string operation

    Args:
        xyzs: e.g. ['x', 'y+1/2', '-z']
        permuation: list, e.g., [0, 2, 1]

    Returns:
        the new xyz string after transformation
    """
    if permutation == [0,1,2]:
        return xyzs
    else:
        new = []
        for xyz in xyzs:
            tmp = xyz.replace(" ","").split(',')
            tmp = [tmp[it] for it in permutation]
            if permutation == [1,0,2]: #a,b
                tmp[0] = tmp[0].replace('y','x')
                tmp[1] = tmp[1].replace('x','y')
            elif permutation == [2,1,0]: #a,c
                tmp[0] = tmp[0].replace('z','x')
                tmp[2] = tmp[2].replace('x','z')
            elif permutation == [0,2,1]: #b,c
                tmp[1] = tmp[1].replace('z','y')
                tmp[2] = tmp[2].replace('y','z')
            elif permutation == [1,2,0]: #b,c
                tmp[0] = tmp[0].replace('y','x')
                tmp[1] = tmp[1].replace('z','y')
                tmp[2] = tmp[2].replace('x','z')
            elif permutation == [2,0,1]: #b,c
                tmp[0] = tmp[0].replace('z','x')
                tmp[1] = tmp[1].replace('x','y')
                tmp[2] = tmp[2].replace('y','z')
            new.append(tmp[0] + ", " + tmp[1] + ", " + tmp[2])

        return new

def swap_xyz_ops(ops, permutation):
    """
    change the symmetry operation by swaping the axes

    Args: 
        ops: SymmOp object
        permutation: list, e.g. [0, 1, 2]

    Returns:
        the new xyz string after transformation
    """
    if permutation == [0,1,2]:
        return ops
    else:
        new = []
        for op in ops:
            m = op.affine_matrix.copy()
            m[:3,:] = m[permutation, :]
            for row in range(3):
                # [0, y+1/2, 1/2] -> (0, y, 1/2)
                if np.abs(m[row,:3]).sum()>0:
                    m[row, 3] = 0
            m[:3,:3] = m[:3, permutation]
            new.append(SymmOp(m))
        return new

def op_transform(ops, affine_matrix):
    """
    x, y, z -> x+1/2, y+1/2, z
    0, 1/2, z -> 1/2, 0, z

    Args: 
        ops: SymmOp object
        permutation: list, e.g. [0, 1, 2]

    Returns:
        the new SymmOp object
    """
    matrix2 = affine_matrix.dot(ops.affine_matrix)
    return SymmOp(matrix2)

def op_translation(op, tran):
    m = op.affine_matrix.copy()
    m[:3,3] += tran
    for row in range(3):
        # [0, y+1/2, 1/2] -> (0, y, 1/2)
        if np.abs(m[row,:3]).sum()>0:
            m[row, 3] = 0
    return SymmOp(m)

def are_equivalent_ops(op1, op2, tol=1e-2):
    """
    check if two ops are equivalent
    """
    diff = op1.affine_matrix - op2.affine_matrix
    diff[:,3] -= np.round(diff[:,3])
    diff = np.abs(diff.flatten())
    if np.sum(diff) < tol:
        return True
    else:
        return False
        

def letter_from_index(index, group, dim=3):
    """
    Given a Wyckoff position's index within a spacegroup, return its number
    and letter e.g. '4a'

    Args:
        index: a single integer describing the WP's index within the
            spacegroup (0 is the general position)
        group: an unorganized Wyckoff position array or Group object (preferred)
        dim: the periodicity dimension of the symmetry group. Used for consideration
            of "o" Wyckoff positions in point groups. Not used if group is a Group
   
    Returns:
        the Wyckoff letter corresponding to the Wyckoff position (for example,
        for position 4a, the function would return 'a')
    """
    letters1 = letters
    # See whether the group has an "o" Wyckoff position
    checko = False
    if type(group) == Group and group.dim == 0:
        checko = True
    elif dim == 0:
        checko = True
    if checko is True:
        if len(group[-1]) == 1 and group[-1][0] == SymmOp.from_xyz_string("0,0,0"):
            # o comes before a
            letters1 = "o" + letters
    length = len(group)
    return letters1[length - 1 - index]


def index_from_letter(letter, group, dim=3):
    """
    Given the Wyckoff letter, returns the index of a Wyckoff position within
    the spacegroup

    Args:
        letter: The wyckoff letter
        group: an unorganized Wyckoff position array or Group object (preferred)
        dim: the periodicity dimension of the symmetry group. Used for consideration
            of "o" Wyckoff positions in point groups. Not used if group is a Group

    Returns:
        a single index specifying the location of the Wyckoff position within
        the spacegroup (0 is the general position)
    """
    letters1 = letters
    # See whether the group has an "o" Wyckoff position
    checko = False
    if type(group) == Group and group.dim == 0:
        checko = True
    elif dim == 0:
        checko = True
    if checko is True:
        if len(group[-1]) == 1 and group[-1][0] == SymmOp.from_xyz_string("0,0,0"):
            # o comes before a
            letters1 = "o" + letters
    length = len(group)
    return length - 1 - letters1.index(letter)


def jk_from_i(i, olist):
    """
    Given an organized list (Wyckoff positions or orientations), determine the
    two indices which correspond to a single index for an unorganized list.
    Used mainly for organized Wyckoff position lists, but can be used for other
    lists organized in a similar way

    Args:
        i: a single index corresponding to the item's location in the
            unorganized list
        olist: the organized list

    Returns:
        [j, k]: two indices corresponding to the item's location in the
            organized list
    """
    num = -1
    found = False
    for j, a in enumerate(olist):
        for k, b in enumerate(a):
            num += 1
            if num == i:
                return [j, k]
    return None


def i_from_jk(j, k, olist):
    """
    Inverse operation of jk_from_i: gives one list index from 2

    Args:
        j, k: indices corresponding to the location of an element in the
            organized list
        olist: the organized list of Wyckoff positions or molecular orientations

    Returns:
        i: one index corresponding to the item's location in the
            unorganized list    
    """
    num = -1
    for x, a in enumerate(olist):
        for y, b in enumerate(a):
            num += 1
            if x == j and y == k:
                return num
    return None


def ss_string_from_ops(ops, number, dim=3, complete=True):
    """
    Print the Hermann-Mauguin symbol for a site symmetry group, using a list of
    SymmOps as input. Note that the symbol does not necessarily refer to the
    x,y,z axes. For information on reading these symbols, see:
    http://en.wikipedia.org/wiki/Hermann-Mauguin_notation#Point_groups

    Args:
        ops: a list of SymmOp objects representing the site symmetry
        number: International number of the symmetry group. Used to determine which
            axes to show. For example, a 3-fold rotation in a cubic system is
            written as ".3.", whereas a 3-fold rotation in a trigonal system is
            written as "3.."
        dim: the dimension of the crystal. Also used to determine notation type
        complete: whether or not all symmetry operations in the group
            are present. If False, we generate the rest

    Returns:
        a string representing the site symmetry (e.g., `2mm`)
    """
    # TODO: Automatically detect which symm_type to use based on ops
    # Determine which notation to use
    symm_type = "high"
    if dim == 3:
        if number >= 1 and number <= 74:
            # Triclinic, monoclinic, orthorhombic
            symm_type = "low"
        elif number >= 75 and number <= 194:
            # Trigonal, Hexagonal, Tetragonal
            symm_type = "medium"
        elif number >= 195 and number <= 230:
            # cubic
            symm_type = "high"
    if dim == 2:
        if number >= 1 and number <= 48:
            # Triclinic, monoclinic, orthorhombic
            symm_type = "low"
        elif number >= 49 and number <= 80:
            # Trigonal, Hexagonal, Tetragonal
            symm_type = "medium"
    if dim == 1:
        if number >= 1 and number <= 22:
            # Triclinic, monoclinic, orthorhombic
            symm_type = "low"
        elif number >= 23 and number <= 75:
            # Trigonal, Hexagonal, Tetragonal
            symm_type = "medium"

    # TODO: replace sg with number, add dim variable
    # Return the symbol for a single axis
    # Will be called later in the function
    def get_symbol(opas, order, has_reflection):
        # ops: a list of Symmetry operations about the axis
        # order: highest order of any symmetry operation about the axis
        # has_reflection: whether or not the axis has mirror symmetry
        if has_reflection is True:
            # rotations have priority
            for opa in opas:
                if opa.order == order and opa.type == "rotation":
                    return str(opa.rotation_order) + "/m"
            for opa in opas:
                if (
                    opa.order == order
                    and opa.type == "rotoinversion"
                    and opa.order != 2
                ):
                    return "-" + str(opa.rotation_order)
            return "m"
        elif has_reflection is False:
            # rotoinversion has priority
            for opa in opas:
                if opa.order == order and opa.type == "rotoinversion":
                    return "-" + str(opa.rotation_order)
            for opa in opas:
                if opa.order == order and opa.type == "rotation":
                    return str(opa.rotation_order)
            return "."

    # Given a list of single-axis symbols, return the one with highest symmetry
    # Will be called later in the function
    def get_highest_symbol(symbols):
        symbol_list = [
            ".",
            "2",
            "m",
            "-2",
            "2/m",
            "3",
            "4",
            "-4",
            "4/m",
            "-3",
            "6",
            "-6",
            "6/m",
        ]
        max_index = 0
        use_list = True
        for j, symbol in enumerate(symbols):
            if symbol in symbol_list:
                i = symbol_list.index(symbol)
            else:
                use_list = False
                num_str = "".join(c for c in symbol if c.isdigit())
                i1 = int(num_str)
                if "m" in symbol or "-" in symbol:
                    if i1 % 2 == 0:
                        i = i1
                    elif i1 % 2 == 1:
                        i = i1 * 2
                else:
                    i = i1
            if i > max_index:
                max_j = j
                max_index = i
        if use_list is True:
            return symbol_list[max_index]
        else:
            return symbols[max_j]

    # Return whether or not two axes are symmetrically equivalent
    # It is assumed that both axes possess the same symbol
    # Will be called within combine_axes
    def are_symmetrically_equivalent(index1, index2):
        axis1 = axes[index1]
        axis2 = axes[index2]
        condition1 = False
        condition2 = False
        # Check for an operation mapping one axis onto the other
        for op in ops:
            if condition1 is False or condition2 is False:
                new1 = op.operate(axis1)
                new2 = op.operate(axis2)
                if np.isclose(abs(np.dot(new1, axis2)), 1):
                    condition1 = True
                if np.isclose(abs(np.dot(new2, axis1)), 1):
                    condition2 = True
        if condition1 is True and condition2 is True:
            return True
        else:
            return False

    # Given a list of axis indices, return the combined symbol
    # Axes may or may not be symmetrically equivalent, but must be of the same
    # type (x/y/z, face-diagonal, body-diagonal)
    # Will be called for mid- and high-symmetry crystallographic point groups
    def combine_axes(indices):
        symbols = {}
        for index in deepcopy(indices):
            symbol = get_symbol(params[index], orders[index], reflections[index])
            if symbol == ".":
                indices.remove(index)
            else:
                symbols[index] = symbol
        if indices == []:
            return "."
        # Remove redundant axes
        for i in deepcopy(indices):
            for j in deepcopy(indices):
                if j > i:
                    if symbols[i] == symbols[j]:
                        if are_symmetrically_equivalent(i, j):
                            if j in indices:
                                indices.remove(j)
        # Combine symbols for non-equivalent axes
        new_symbols = []
        for i in indices:
            new_symbols.append(symbols[i])
        symbol = ""
        while new_symbols != []:
            highest = get_highest_symbol(new_symbols)
            symbol += highest
            new_symbols.remove(highest)
        if symbol == "":
            printx("Error: could not combine site symmetry axes.", priority=1)
            return
        else:
            return symbol

    # Generate needed ops
    if not complete:
        ops = generate_full_symmops(ops, 1e-3)
    # Get OperationAnalyzer object for all ops
    opas = []
    for op in ops:
        opas.append(OperationAnalyzer(op))
    # Store the symmetry of each axis
    params = [[], [], [], [], [], [], [], [], [], [], [], [], []]
    has_inversion = False
    # Store possible symmetry axes for crystallographic point groups
    axes = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [0, 1, 1],
        [1, 0, 1],
        [1, -1, 0],
        [0, 1, -1],
        [1, 0, -1],
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ]
    for i, axis in enumerate(axes):
        axes[i] = axis / np.linalg.norm(axis)
    for opa in opas:
        if opa.type != "identity" and opa.type != "inversion":
            found = False
            for i, axis in enumerate(axes):
                if np.isclose(abs(np.dot(opa.axis, axis)), 1):
                    found = True
                    params[i].append(opa)
            # Store uncommon axes for trigonal and hexagonal lattices
            if found is False:
                axes.append(opa.axis)
                # Check that new axis is not symmetrically equivalent to others
                unique = True
                for i, axis in enumerate(axes):
                    if i != len(axes) - 1:
                        if are_symmetrically_equivalent(i, len(axes) - 1):
                            unique = False
                if unique is True:
                    params.append([opa])
                elif unique is False:
                    axes.pop()
        elif opa.type == "inversion":
            has_inversion = True
    # Determine how many high-symmetry axes are present
    n_axes = 0
    # Store the order of each axis
    orders = []
    # Store whether or not each axis has reflection symmetry
    reflections = []
    for axis in params:
        order = 1
        high_symm = False
        has_reflection = False
        for opa in axis:
            if opa.order >= 3:
                high_symm = True
            if opa.order > order:
                order = opa.order
            if opa.order == 2 and opa.type == "rotoinversion":
                has_reflection = True
        orders.append(order)
        if high_symm == True:
            n_axes += 1
        reflections.append(has_reflection)
    # Triclinic, monoclinic, orthorhombic
    # Positions in symbol refer to x,y,z axes respectively
    if symm_type == "low":
        symbol = (
            get_symbol(params[0], orders[0], reflections[0])
            + get_symbol(params[1], orders[1], reflections[1])
            + get_symbol(params[2], orders[2], reflections[2])
        )
        if symbol != "...":
            return symbol
        elif symbol == "...":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    # Trigonal, Hexagonal, Tetragonal
    elif symm_type == "medium":
        # 1st symbol: z axis
        s1 = get_symbol(params[2], orders[2], reflections[2])
        # 2nd symbol: x or y axes (whichever have higher symmetry)
        s2 = combine_axes([0, 1])
        # 3rd symbol: face-diagonal axes (whichever have highest symmetry)
        s3 = combine_axes(list(range(3, len(axes))))
        symbol = s1 + " " + s2 + " " + s3
        if symbol != ". . .":
            return symbol
        elif symbol == ". . .":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    # Cubic
    elif symm_type == "high":
        pass
        # 1st symbol: x, y, and/or z axes (whichever have highest symmetry)
        s1 = combine_axes([0, 1, 2])
        # 2nd symbol: body-diagonal axes (whichever has highest symmetry)
        s2 = combine_axes([9, 10, 11, 12])
        # 3rd symbol: face-diagonal axes (whichever have highest symmetry)
        s3 = combine_axes([3, 4, 5, 6, 7, 8])
        symbol = s1 + " " + s2 + " " + s3
        if symbol != ". . .":
            return symbol
        elif symbol == ". . .":
            if has_inversion is True:
                return "-1"
            else:
                return "1"
    else:
        printx("Error: invalid spacegroup number", priority=1)
        return


def organized_wyckoffs(group):
    """
    Takes a Group object or unorganized list of Wyckoff positions and returns
    a 2D list of Wyckoff positions organized by multiplicity.

    Args:
        group: a pyxtal.symmetry.Group object
    
    Returns:
        a 2D list of Wyckoff_position objects if group is a Group object.
        a 3D list of SymmOp objects if group is a 2D list of SymmOps
    """
    if type(group) == Group:
        wyckoffs = group.Wyckoff_positions
    else:
        wyckoffs = group
    wyckoffs_organized = [[]]  # 2D Array of WP's organized by multiplicity
    old = len(wyckoffs[0])
    for wp in wyckoffs:
        mult = len(wp)
        if mult != old:
            wyckoffs_organized.append([])
            old = mult
        wyckoffs_organized[-1].append(wp)
    return wyckoffs_organized


def symmetry_element_from_axis(axis):
    """
    Given an axis, returns a SymmOp representing a symmetry element on the axis.
    For example, the symmetry element for the vector (0,0,2) would be (0,0,z).
    
    Args:
        axis: a 3-vector representing the symmetry element

    Returns:
        a SymmOp object of form (ax, bx, cx), (ay, by, cy), or (az, bz, cz)
    """
    if len(axis) != 3:
        return
    # Vector must be non-zero
    if axis.dot(axis) < 1e-6:
        return
    v = np.array(axis) / np.linalg.norm(axis)
    # Find largest component (x, y, or z)
    abs_vals = [abs(a) for a in v]
    f1 = max(abs_vals)
    index1 = list(abs_vals).index(f1)
    # Initialize an affine matrix
    m = np.eye(4)
    m[:3] = [0.0, 0.0, 0.0, 0.0]
    # Set values for affine matrix
    m[:3, index1] = v
    return SymmOp(m)


def get_wyckoffs(num, organized=False, PBC=[1, 1, 1], dim=3):
    """
    Returns a list of Wyckoff positions for a given space group. Has option to
    organize the list based on multiplicity (this is used for
    random_crystal.wyckoffs) 
    
    For an unorganized list:

        - 1st index: index of WP in sg (0 is the WP with largest multiplicity)
        - 2nd index: a SymmOp object in the WP

    For an organized list:

        - 1st index: specifies multiplicity (0 is the largest multiplicity)
        - 2nd index: corresponds to a WP within the group of equal multiplicity.
        - 3nd index: corresponds to a SymmOp object within the Wyckoff position

    You may switch between organized and unorganized lists using the methods
    i_from_jk and jk_from_i. For example, if a Wyckoff position is the [i]
    entry in an unorganized list, it will be the [j][k] entry in an organized
    list.

    Args:
        num: the international group number
        dim: dimension [0, 1, 2, 3]
        organized: whether or not to organize the list based on multiplicity
        PBC: A periodic boundary condition list, 1 means periodic, 0 means not periodic.
            Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity along the z axis
    
    Returns: 
        a list of Wyckoff positions, each of which is a list of SymmOp's
    """
    if dim == 3:
        if PBC != [1, 1, 1]:
            coor = [0, 0, 0]
            for i, a in enumerate(PBC):
                if not a:
                    coor[i] = 0.5
            coor = np.array(coor)

        wyckoff_strings = eval(wyckoff_df["0"][num])
        wyckoffs = []
        for x in wyckoff_strings:
            if PBC != [1, 1, 1]:
                op = SymmOp.from_xyz_string(x[0])
                coor1 = op.operate(coor)
                invalid = False
                for i, a in enumerate(PBC):
                    if not a:
                        if not abs(coor1[i] - 0.5) < 1e-2:
                            # invalid wyckoffs for layer group
                            invalid = True
                if invalid is False:
                    wyckoffs.append([])
                    for y in x:
                        wyckoffs[-1].append(SymmOp.from_xyz_string(y))
            else:
                wyckoffs.append([])
                for y in x:
                    wyckoffs[-1].append(SymmOp.from_xyz_string(y))
        if organized:
            wyckoffs_organized = [[]]  # 2D Array of WP's organized by multiplicity
            old = len(wyckoffs[0])
            for wp in wyckoffs:
                mult = len(wp)
                if mult != old:
                    wyckoffs_organized.append([])
                    old = mult
                wyckoffs_organized[-1].append(wp)
            return wyckoffs_organized
        else:
            return wyckoffs

    else:
        if dim == 2:
            wyckoff_strings = eval(layer_df["0"][num])
        elif dim == 1:
            wyckoff_strings = eval(rod_df["0"][num])
        elif dim == 0:
            wyckoff_strings = eval(point_df["0"][num])

        wyckoffs = []
        for x in wyckoff_strings:
            wyckoffs.append([])
            for y in x:
                if dim == 0:
                    wyckoffs[-1].append(SymmOp(y))
                else:
                    wyckoffs[-1].append(SymmOp.from_xyz_string(y))
        if organized:
            wyckoffs_organized = [[]]  # 2D Array of WP's organized by multiplicity
            old = len(wyckoffs[0])
            for wp in wyckoffs:
                mult = len(wp)
                if mult != old:
                    wyckoffs_organized.append([])
                    old = mult
                wyckoffs_organized[-1].append(wp)
            return wyckoffs_organized
        else:
            return wyckoffs



def get_wyckoff_symmetry(num, dim=3, PBC=[1, 1, 1]):
    """
    Returns a list of Wyckoff position site symmetry for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a point within the WP
    3rd index: a site symmetry SymmOp of the point

    Args:
        sg: the international spacegroup number
        dim: 0, 1, 2, 3
        PBC: A periodic boundary condition list, Ex: [1,1,1] -> full 3d periodic

    Returns:
        a 3d list of SymmOp objects representing the site symmetry of each
        point in each Wyckoff position
    """
    if dim == 3:
        if PBC != [1, 1, 1]:
            coor = [0, 0, 0]
            for i, a in enumerate(PBC):
                if not a:
                    coor[i] = 0.5
            coor = np.array(coor)
        wyckoffs = get_wyckoffs(num, PBC=PBC)

        symmetry_strings = eval(wyckoff_symmetry_df["0"][num])
        symmetry = []
        # Loop over Wyckoff positions
        for x, w in zip(symmetry_strings, wyckoffs):
            if PBC != [1, 1, 1]:
                op = w[0]
                coor1 = op.operate(coor)
                invalid = False
                for i, a in enumerate(PBC):
                    if not a:
                        if not abs(coor1[i] - 0.5) < 1e-2:
                            invalid = True
                if not invalid:
                    symmetry.append([])
                    # Loop over points in WP
                    for y in x:
                        symmetry[-1].append([])
                        # Loop over ops
                        for z in y:
                            op = SymmOp.from_xyz_string(z)
                            symmetry[-1][-1].append(op)
            else:
                symmetry.append([])
                # Loop over points in WP
                for y in x:
                    symmetry[-1].append([])
                    # Loop over ops
                    for z in y:
                        op = SymmOp.from_xyz_string(z)
                        symmetry[-1][-1].append(op)
        return symmetry

    else:
        if dim == 2:
            symmetry_strings = eval(layer_symmetry_df["0"][num])
        elif dim == 1:
            symmetry_strings = eval(rod_symmetry_df["0"][num])
        elif dim == 0:
            symmetry_strings = eval(point_symmetry_df["0"][num])

        symmetry = []
        # Loop over Wyckoff positions
        for x in symmetry_strings:
            symmetry.append([])
            # Loop over points in WP
            for y in x:
                symmetry[-1].append([])
                # Loop over ops
                for z in y:
                    if dim == 0:
                        op = SymmOp(z)
                    else:
                        op = SymmOp.from_xyz_string(z)
                    symmetry[-1][-1].append(op)
        return symmetry
        

def get_generators(num, dim=3, PBC=[1, 1, 1]):
    """
    Returns a list of Wyckoff generators for a given space group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a generator for the WP
    This function is useful for rotating molecules based on Wyckoff position,
    since special Wyckoff positions only encode positional information, but not
    information about the orientation. The generators for each Wyckoff position
    form a subset of the spacegroup's general Wyckoff position.
    
    Args:
        num: the international spacegroup number
        dim: dimension
        PBC: A periodic boundary condition list, where 1 is periodic, 0 is nonperiodic.
            Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity along the z axis
    
    Returns:
        a 2d list of symmop objects [[wp0], [wp1], ... ]
    """

    generators = []
    if dim == 3:
        generator_strings = eval(wyckoff_generators_df["0"][num])

        # Loop over Wyckoff positions
        wyckoffs = get_wyckoffs(num, PBC=PBC)
        for x, w in zip(generator_strings, wyckoffs):
            if PBC != [1, 1, 1]:

                coor = [0, 0, 0]
                for i, a in enumerate(PBC):
                    if not a:
                        coor[i] = 0.5
                coor = np.array(coor)

                op = w[0]
                coor1 = op.operate(coor)
                valid = True
                for i, a in enumerate(PBC):
                    if not a and not abs(coor1[i] - 0.5) < 1e-2:
                        valid = False
                        break
                if valid:
                    generators.append([])
                    # Loop over ops
                    for y in x:
                        op = SymmOp.from_xyz_string(y)
                        generators[-1].append(op)
            else:
                generators.append([])
                for y in x:
                    op = SymmOp.from_xyz_string(y)
                    generators[-1].append(op)
        return generators

    else:
        if dim == 2:
            generator_strings = eval(layer_generators_df["0"][num])
        elif dim == 1:
            generator_strings = eval(rod_generators_df["0"][num])
        elif dim == 0:
            generator_strings = eval(point_generators_df["0"][num])

        # Loop over Wyckoff positions
        for x in generator_strings:
            generators.append([])
            # Loop over ops
            for y in x:
                if dim > 0:
                    op = SymmOp.from_xyz_string(y)
                else:
                    op = SymmOp(y)
                generators[-1].append(op)
        return generators


def site_symm(point, gen_pos, tol=1e-3, lattice=np.eye(3), PBC=None):
    """
    Given a point and a general Wyckoff position, return the list of symmetry
    operations leaving the point (coordinate or SymmOp) invariant. The returned
    SymmOps are a subset of the general position. The site symmetry can be used
    for determining the Wyckoff position for a set of points, or for
    determining the valid orientations of a molecule within a given Wyckoff
    position.

    Args:
        point: a 1x3 coordinate or SymmOp object to find the symmetry of. If a
            SymmOp is given, the returned symmetries must also preserve the
            point's orientaion
        gen_pos: the general position of the spacegroup. Can be a Wyckoff_position
            object or list of SymmOp objects.
        tol:
            the numberical tolerance for determining equivalent positions and
            orientations.
        lattice:
            a 3x3 matrix representing the lattice vectors of the unit cell
        PBC: A periodic boundary condition list, 1 means periodic, 0 means not periodic.
            Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity along the z axis.
            Need not be defined here if gen_pos is a Wyckoff_position object.

    Returns:
        a list of SymmOp objects which leave the given point invariant
    """
    if PBC is None:
        if type(gen_pos) == Wyckoff_position:
            PBC = gen_pos.PBC
        else:
            PBC = [1, 1, 1]
    # Convert point into a SymmOp
    if type(point) != SymmOp:
        point = SymmOp.from_rotation_and_translation(
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]], np.array(point)
        )
    symmetry = []
    for op in gen_pos:
        is_symmetry = True
        # Calculate the effect of applying op to point
        difference = SymmOp((op * point).affine_matrix - point.affine_matrix)
        # Check that the rotation matrix is unaltered by op
        if not np.allclose(
            difference.rotation_matrix, np.zeros((3, 3)), rtol=1e-3, atol=1e-3
        ):
            is_symmetry = False
        # Check that the displacement is less than tol
        displacement = difference.translation_vector
        if distance(displacement, lattice, PBC=PBC) > tol:
            is_symmetry = False
        if is_symmetry:
            """The actual site symmetry's translation vector may vary from op by
            a factor of +1 or -1 (especially when op contains +-1/2).
            We record this to distinguish between special Wyckoff positions.
            As an example, consider the point (-x+1/2,-x,x+1/2) in position 16c
            of space group Ia-3(206). The site symmetry includes the operations
            (-z+1,x-1/2,-y+1/2) and (y+1/2,-z+1/2,-x+1). These operations are
            not listed in the general position, but correspond to the operations
            (-z,x+1/2,-y+1/2) and (y+1/2,-z+1/2,-x), respectively, just shifted
            by (+1,-1,0) and (0,0,+1), respectively.
            """
            el = SymmOp.from_rotation_and_translation(
                op.rotation_matrix, op.translation_vector - np.round(displacement)
            )
            symmetry.append(el)
    return symmetry

def check_wyckoff_position(points, group, tol=1e-3):
    """
    Given a list of points, returns a single index of a matching Wyckoff
    position in the space group. Checks the site symmetry of each supplied
    point against the site symmetry for each point in the Wyckoff position.
    Also returns a point which can be used to generate the rest using the
    Wyckoff position operators
    Args:
        points: a list of 3d coordinates or SymmOps to check
        group: a Group object
        tol: the max distance between equivalent points

    Returns:
        index, p: index is a single index for the Wyckoff position within
        the sg. If no matching WP is found, returns False. point is a
        coordinate taken from the list points. When plugged into the Wyckoff
        position, it will generate all the other points.
    """
    points = np.array(points)
    wyckoffs = group.wyckoffs
    w_symm_all = group.w_symm
    PBC = group.PBC
    # new method
    # Store the squared distance tolerance
    t = tol ** 2
    # Loop over Wyckoff positions
    for i, wp in enumerate(wyckoffs):
        # Check that length of points and wp are equal
        if len(wp) != len(points):
            continue
        failed = False

        # Search for a generating point
        for p in points:
            failed = False
            # Check that point works as x,y,z value for wp
            xyz = filtered_coords_euclidean(wp[0].operate(p) - p, PBC=PBC)

            if xyz.dot(xyz) > t:
                continue
            # Calculate distances between original and generated points
            pw = np.array([op.operate(p) for op in wp])
            dw = distance_matrix(points, pw, None, PBC=PBC, metric="sqeuclidean")

            # Check each row for a zero
            for row in dw:
                num = (row < t).sum()
                if num < 1:
                    failed = True
                    break

            if failed:
                continue
            # Check each column for a zero
            for column in dw.T:
                num = (column < t).sum()
                if num < 1:
                    failed = True
                    break

            # Calculate distance between original and generated points
            ps = np.array([op.operate(p) for op in w_symm_all[i][0]])
            ds = distance_matrix([p], ps, None, PBC=PBC, metric="sqeuclidean")
            # Check whether any generated points are too far away
            num = (ds > t).sum()
            if num > 0:
                failed = True

            if failed:
                continue
            return i, p
    return False, None

def get_symbol_and_number(input_group, dim=3):
    """
    Function for quick conversion between symbols and numbers

    Args:
        input_group: the group symbol or international number
        dim: the periodic dimension of the group
    """

    keys = {
        3: "space_group",
        2: "layer_group",
        1: "rod_group",
        0: "point_group",
    }

    found = False
    lists = symbols[keys[dim]]
    number = None
    symbol = None
    if dim not in [0, 1, 2, 3]:
        raise ValueError("Dimension ({:d}) should in [0, 1, 2, 3] ".format(dim))

    if type(input_group) == str:
        for i, _symbol in enumerate(lists):
            if _symbol == input_group:
                number = i + 1
                symbol = input_group
                return symbol, number
        msg = "({:s}) not found in {:s} ".format(input_group, keys[dim])
        raise ValueError(msg)
    else:
        valid, msg = check_symmetry_and_dim(input_group, dim)
        if not valid:
            raise ValueError(msg)
        else:
            number = input_group
            symbol = lists[number - 1]
            return symbol, number

def check_symmetry_and_dim(number, dim=3):
    """
    check if it is a valid number for the given symmetry
    
    Args:
        number: int
        dim: 0, 1, 2, 3
    """
    valid = True
    msg = 'This is a valid group number'

    numbers = [56, 75, 80, 230]
    if dim not in [0, 1, 2, 3]:
        msg = "invalid dimension {:d}".format()
        valid = False
    else:
        max_num = numbers[dim]
        if number not in range(1, max_num+1):
            valid = False
            msg = "invalid symmetry group {:d}".format(number)
            msg += " in dimension {:d}".format(dim)
    return valid, msg

def get_pbc_and_lattice(number, dim):
    if dim == 3:
        PBC = [1, 1, 1]
        if number <= 2:
            lattice_type = "triclinic"
        elif number <= 15:
            lattice_type = "monoclinic"
        elif number <= 74:
            lattice_type = "orthorhombic"
        elif number <= 142:
            lattice_type = "tetragonal"
        elif number <= 194:
            lattice_type = "hexagonal"
        elif number <= 230:
            lattice_type = "cubic"
    elif dim == 2:
        PBC = [1, 1, 0]
        if number <= 2:
            lattice_type = "triclinic"
        elif number <= 18:
            lattice_type = "monoclinic"
        elif number <= 48:
            lattice_type = "orthorhombic"
        elif number <= 64:
            lattice_type = "tetragonal"
        elif number <= 80:
            lattice_type = "hexagonal"
    elif dim == 1:
        PBC = [0, 0, 1]
        if number <= 2:
            lattice_type = "triclinic"
        elif number <= 12:
            lattice_type = "monoclinic"
        elif number <= 22:
            lattice_type = "orthorhombic"
        elif number <= 41:
            lattice_type = "tetragonal"
        elif number <= 75:
            lattice_type = "hexagonal"
    elif dim == 0:
        PBC = [0, 0, 0]
        # "C1", "Ci", "D2", "D2h", "T", "Th",
        # "O", "Td", "Oh", "I", "Ih",
        if number in [1, 2, 6, 8, 28, 29, 30, 31, 32, 55, 56]:
            lattice_type = "spherical"
        else:
            lattice_type = "ellipsoidal"
    return PBC, lattice_type

def search_matched_position(G, wp, pos):
    """
    search generator for a special Wyckoff position

    Args:
        G: space group number or Group object
        wp: Wyckoff object
        pos: initial xyz position

    Return:
        pos1: the position that matchs the standard setting
    """
    wp0 = G[0]
    match = False

    for op in wp0:
        pos1 = op.operate(pos)
        pos0 = wp[0].operate(pos1)
        diff = pos1 - pos0
        diff -= np.round(diff)
        diff = np.abs(diff)
        #print(wp.letter, pos1, pos0, diff)
        if diff.sum()<1e-2:
            pos1 -= np.floor(pos1)
            match = True
            break
    if match:
        return pos1
    else:
        return None

def search_matched_positions(G, wp, pos):
    """
    search generator for a special Wyckoff position

    Args:
        G: space group number or Group object
        wp: Wyckoff object
        pos: initial xyz position

    Return:
        pos1: the position that matchs the standard setting
    """
    wp0 = G[0]
    coords = []

    for op in wp0:
        pos1 = op.operate(pos)
        pos0 = wp[0].operate(pos1)
        diff = pos1 - pos0
        diff -= np.round(diff)
        diff = np.abs(diff)
        #print(wp.letter, pos1, pos0, diff)
        if diff.sum()<1e-2:
            pos1 -= np.floor(pos1)
            coords.append(pos1)
    return coords

def search_cloest_wp(G, wp, op, pos):
    """
    For a given position, search for the cloest wp which
    satisfies the desired symmetry relation
    e.g., for pos (0.1, 0.12, 0.2) and op (x, x, z)
    the closest match is (0.11, 0.11, 0.2)

    Args:
        G: space group number or Group object
        wp: Wyckoff object
        op: symmetry operation belonging to wp 
        pos: initial xyz position

    Return:
        pos1: the position that matchs symmetry operation
    """
    #G = Group(wp.number)
    if np.linalg.matrix_rank(op.rotation_matrix) == 0:
        # fixed point (e.g, 1/2, 1/2, 1/2)
        return op.translation_vector
    elif np.linalg.matrix_rank(op.rotation_matrix) == 3:
        # fully independent, e.g., (x,y,z), (-x,y,z)
        return pos
    else:
        # check if this is already matched
        coords = search_matched_positions(G, wp, pos)
        if len(coords)>0:
            diffs = []
            for coord in coords:
                tmp = op.operate(coord)
                diff1 = tmp - pos
                diff1 -= np.round(diff1)
                dist = np.linalg.norm(diff1) 
                if dist < 1e-3:
                    return tmp
                else:
                    diffs.append(dist)
            minID = np.argmin(diffs)
            return op.operate(coords[minID])

        # if not match, search for the closest solution
        else:
            wp0 = G[0]
            # extract all possible xyzs
            all_xyz = apply_ops(pos, wp0)[1:]
            dists = all_xyz - pos
            dists -= np.round(dists)
            ds = np.linalg.norm(dists, axis=1)
            ids = np.argsort(ds)
            for id in ids:
                d = all_xyz[id] - pos
                d -= np.round(d)
                res = pos + d/2
                if search_matched_position(G, wp, res) is not None:
                    #print(ds[id], pos, res)
                    return res
            return op.operate(pos)


def get_point_group(number):
    """
    Parse the point group symmetry info from space group 
    http://img.chem.ucl.ac.uk/sgp/misc/pointgrp.htm
    Among 32(230) point(space) groups, there are
    10(68) polar groups,
    11(92) centrosymmetric groups,
    11(65) enantiomorphic groups

    Args:
        number: space group number

    Return:
        point group symbol
        polar: 1, 2, m, mm2, 3, 3m, 4, 4mm, 6, 6mm
        centrosymmetry: -1, 2/m, mmm, 4/m, 4/mmm, -3, -3m, 6/m, 6/mmm, m-3, m-3m
        enantiomorphic: 1, 2, 222, 4, 422, 3, 32, 6, 622, 23, 432
    """

    if number == 1:
        return '1', True, False, True
    elif number == 2:
        return '-1', False, True, False
    elif 3 <= number <= 5:
        return '2', True, False, True
    elif 6 <= number <= 9:
        return 'm', True, False, False
    elif 10 <= number <= 15:
        return '2/m', False, True, False
    elif 16 <= number <= 24:
        return '222', False, False, True
    elif 25 <= number <= 46:
        return 'mm2', True, False, False
    elif 47 <= number <= 74:
        return 'mmm', False, True, False
    elif 75 <= number <= 80:
        return '4', True, False, True
    elif 81 <= number <= 82:
        return '-4', False, False, False
    elif 83 <= number <= 88:
        return '4/m', False, True, False
    elif 89 <= number <= 98:
        return '422', False, False, True
    elif 99 <= number <= 110:
        return '4mm', True, False, False
    elif 111 <= number <= 122:
        return '-42m', False, False, False
    elif 123 <= number <= 142:
        return '4/mmm', False, True, False
    elif 143 <= number <= 146:
        return '3', True, False, True
    elif 147 <= number <= 148:
        return '-3', False, True, False
    elif 149 <= number <= 155:
        return '32', False, False, True
    elif 156 <= number <= 161:
        return '3m', True, False, False
    elif 162 <= number <= 167:
        return '-3m', False, True, False
    elif 168 <= number <= 173:
        return '6', True, False, True
    elif number == 174:
        return '-6', False, False, False
    elif 175 <= number <= 176:
        return '6/m', False, True, False
    elif 177 <= number <= 182:
        return '622', False, False, True
    elif 183 <= number <= 186:
        return '6mm', True, False, False
    elif 187 <= number <= 190:
        return '-62m', False, False, False
    elif 191 <= number <= 194:
        return '6/mmm', False, True, False
    elif 195 <= number <= 199:
        return '23', False, False, True
    elif 200 <= number <= 206:
        return 'm-3', False, True, False
    elif 207 <= number <= 214:
        return '432', False, False, True
    elif 215 <= number <= 220:
        return '-43m', False, False, False
    elif 221 <= number <= 230:
        return 'm-3m', False, True, False

def get_close_packed_groups(pg):
    """
    List the close packed groups based on the molcular symmetry
    Compiled from AIK Book, Table 2 P34

    Args:
        pg: point group symbol

    Return
        list of space group numbers
    """

    if pg == '1':
        return [1, 2, 4, 14, 19, 29, 33, 51, 54, 61, 62]
    elif pg == '2':
        return [1, 15, 18, 60]
    elif pg == 'm':
        return [1, 26, 36, 63, 64]
    elif pg == 'I':
        return [1, 2, 14, 15, 61]
    elif pg == 'mm':
        return [42, 51, 59]
    elif pg == '2/m':
        return [12, 54, 64]
    elif pg == '222':
        return [21, 22, 23, 68]
    elif pg == 'mmm':
        return [65, 69, 71]

def para2ferro(pg):
    """
    88 potential paraelectric-to-ferroelectric phase transitions
    https://journals.aps.org/prb/abstract/10.1103/PhysRevB.2.754
    https://pubs.rsc.org/en/content/articlelanding/2016/cs/c5cs00308c

    Args:
        paraelectric point group 

    Returns:
        list of ferroelectric point groups
    """
    #Triclinic: 1
    if pg == '-1': #2
        return ['1'] 
    #Monoclinic: 5
    elif pg in ['2', 'm']: #2
        return '1'
    elif pg == '2/m': #3
        return ['1', 'm', '2']
    #Orthorhombic: #7
    elif pg == '222': #2
        return ['1', '2'] 
    elif pg == 'mm2': #2
        return ['1', 'm'] 
    elif pg == 'mmm': #3
        return ['1', 'm', 'mm2']
    #Tetragonal: 20
    elif pg == '4': #1
        return ['1']
    elif pg == '-4': #2
        return ['1', '2']
    elif pg == '4/m': #3
        return ['1', '2', '4']
    elif pg == '422': #3
        #return ['1', '2(s)', '4']
        return ['1', '2', '4']
    elif pg == '4mm': #2
        return ['1', 'm']
    elif pg == '-42m': #4
        #return ['1', '2(s)', 'm', 'mm2']
        return ['1', '2', 'm', 'mm2']
    elif pg == '4/mmm': #5
        #return ['1', 'm(s)', 'm(p)', 'mm2(s)', '4mm']
        return ['1', 'm', 'mm2', '4mm']
    #Trigonal: 12
    elif pg == '3': #1
        return ['1']
    elif pg == '-3': #2
        return ['1', '3']
    elif pg == '32': #3
        return ['1', '2', '3']
    elif pg == '3m': #2
        return ['1', 'm']
    elif pg == '-3m': #4
        return ['1', '2', 'm', '3m']
    #Hexagonal: 22
    elif pg == '6': #1
        return ['1']
    elif pg == '-6': #3
        return ['1', 'm', '3']
    elif pg == '6/m': #3
        return ['1', 'm', '6']
    elif pg == '622': #3
        #return ['1', '2(s)', '6']
        return ['1', '2', '6']
    elif pg == '6mm': #2
        return ['1', '2']
    elif pg in ['-62m', '-6m2']: #5
        #return ['1', 'm(s)', 'm(p)', 'mm2', '3m']
        return ['1', 'm', 'mm2', '3m']
    elif pg == '6/mmm': #5
        #return ['1', 'm(s)', 'm(p)', 'mm2(s)', '6mm']
        return ['1', 'm', 'mm2', '6mm']
    #Cubic: 21
    elif pg == '23': #3
        return ['1', '2', '3']
    elif pg == 'm-3': #4
        return ['1', 'm', 'mm2', '3']
    elif pg == '432': #4
        #return ['1', '2(s)', '4', '3']
        return ['1', '2', '4', '3']
    elif pg == '-43m': #4
        return ['1', 'm', 'mm2', '3m']
    elif pg == 'm-3m': #6
        #return ['1', 'm(s)', 'm(p)', 'mm2', '4mm', '3m']
        return ['1', 'm', 'mm2', '4mm', '3m']

def get_all_polar_space_groups():
    ps, nps = [], []
    for i in range(1,231):
        g = Group(i, quick=True)
        if g.polar:
            ps.append(i)
        else:
            nps.append(i)
    return ps, nps
