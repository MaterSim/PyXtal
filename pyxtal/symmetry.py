"""
Module for storing & accessing symmetry group information, including
    - Group class
    - Wyckoff_Position class.
    - Hall class
"""

from __future__ import annotations

import functools
import importlib.util
import itertools
import operator
import os
import re
from copy import deepcopy
from ast import literal_eval

import numpy as np
from monty.serialization import loadfn
from numpy.random import Generator
from pandas import read_csv

from pyxtal.constants import all_sym_directions, hex_cell, letters
from pyxtal.operations import (
    OperationAnalyzer,
    SymmOp,
    apply_ops,
    check_images,
    create_matrix,
    distance,
    distance_matrix,
    filtered_coords,
    filtered_coords_euclidean,
)

def rf(package_name, resource_path):
    package_path = importlib.util.find_spec(
        package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)

"""
Properties for Lazy Loading
"""
class SymmetryData:

    _k_subgroup = None
    _t_subgroup = None

    def __init__(self):
        self._wyckoff_sg = None
        self._wyckoff_lg = None
        self._wyckoff_rg = None
        self._wyckoff_pg = None
        self._symmetry_sg = None
        self._symmetry_lg = None
        self._symmetry_rg = None
        self._symmetry_pg = None
        self._generator_sg = None
        self._generator_lg = None
        self._generator_rg = None
        self._generator_pg = None
        self._t_subgroup = None
        self._k_subgroup = None
        self._hall_table = None

    @classmethod
    def get_t_subgroup(cls):
        if cls._t_subgroup is None:
            cls._t_subgroup = loadfn(rf("pyxtal", "database/t_subgroup.json"))
        return cls._t_subgroup

    @classmethod
    def get_k_subgroup(cls):
        if cls._k_subgroup is None:
            cls._k_subgroup = loadfn(rf("pyxtal", "database/k_subgroup.json"))
        return cls._k_subgroup

    def get_wyckoff_sg(self):
        if self._wyckoff_sg is None:
            self._wyckoff_sg = read_csv(rf("pyxtal", "database/wyckoff_list.csv"))
        return self._wyckoff_sg

    def get_wyckoff_lg(self):
        if self._wyckoff_lg is None:
            self._wyckoff_lg = read_csv(rf("pyxtal", "database/layer.csv"))
        return self._wyckoff_lg

    def get_wyckoff_rg(self):
        if self._wyckoff_rg is None:
            self._wyckoff_rg = read_csv(rf("pyxtal", "database/rod.csv"))
        return self._wyckoff_rg

    def get_wyckoff_pg(self):
        if self._wyckoff_pg is None:
            self._wyckoff_pg = read_csv(rf("pyxtal", "database/point.csv"))
        return self._wyckoff_pg

    def get_symmetry_sg(self):
        if self._symmetry_sg is None:
            self._symmetry_sg = read_csv(rf("pyxtal", "database/wyckoff_symmetry.csv"))
        return self._symmetry_sg

    def get_symmetry_lg(self):
        if self._symmetry_lg is None:
            self._symmetry_lg = read_csv(rf("pyxtal", "database/layer_symmetry.csv"))
        return self._symmetry_lg

    def get_symmetry_rg(self):
        if self._symmetry_rg is None:
            self._symmetry_rg = read_csv(rf("pyxtal", "database/rod_symmetry.csv"))
        return self._symmetry_rg

    def get_symmetry_pg(self):
        if self._symmetry_pg is None:
            self._symmetry_pg = read_csv(rf("pyxtal", "database/point_symmetry.csv"))
        return self._symmetry_pg

    def get_generator_sg(self):
        if self._generator_sg is None:
            self._generator_sg = read_csv(rf("pyxtal", "database/wyckoff_generators.csv"))
        return self._generator_sg

    def get_generator_lg(self):
        if self._generator_lg is None:
            self._generator_lg = read_csv(rf("pyxtal", "database/layer_generators.csv"))
        return self._generator_lg

    def get_generator_rg(self):
        if self._generator_rg is None:
            self._generator_rg = read_csv(rf("pyxtal", "database/rod_generators.csv"))
        return self._generator_rg

    def get_generator_pg(self):
        if self._generator_pg is None:
            self._generator_pg = read_csv(rf("pyxtal", "database/point_generators.csv"))
        return self._generator_pg

    def get_hall_table(self):
        if self._hall_table is None:
            self._hall_table = read_csv(rf("pyxtal", "database/HM_Full.csv"), sep=",")
        return self._hall_table

# ------------------------------ Constants ---------------------------------------
SYMDATA = SymmetryData()
HALL_TABLE = SYMDATA.get_hall_table()

spglib_hall_numbers = [
    1,
    2,
    3,
    6,
    9,
    18,
    21,
    30,
    39,
    57,
    60,
    63,
    72,
    81,
    90,
    108,
    109,
    112,
    115,
    116,
    119,
    122,
    123,
    124,
    125,
    128,
    134,
    137,
    143,
    149,
    155,
    161,
    164,
    170,
    173,
    176,
    182,
    185,
    191,
    197,
    203,
    209,
    212,
    215,
    218,
    221,
    227,
    228,
    230,
    233,
    239,
    245,
    251,
    257,
    263,
    266,
    269,
    275,
    278,
    284,
    290,
    292,
    298,
    304,
    310,
    313,
    316,
    322,
    334,
    335,
    337,
    338,
    341,
    343,
    349,
    350,
    351,
    352,
    353,
    354,
    355,
    356,
    357,
    358,
    359,
    361,
    363,
    364,
    366,
    367,
    368,
    369,
    370,
    371,
    372,
    373,
    374,
    375,
    376,
    377,
    378,
    379,
    380,
    381,
    382,
    383,
    384,
    385,
    386,
    387,
    388,
    389,
    390,
    391,
    392,
    393,
    394,
    395,
    396,
    397,
    398,
    399,
    400,
    401,
    402,
    404,
    406,
    407,
    408,
    410,
    412,
    413,
    414,
    416,
    418,
    419,
    420,
    422,
    424,
    425,
    426,
    428,
    430,
    431,
    432,
    433,
    435,
    436,
    438,
    439,
    440,
    441,
    442,
    443,
    444,
    446,
    447,
    448,
    449,
    450,
    452,
    454,
    455,
    456,
    457,
    458,
    460,
    462,
    463,
    464,
    465,
    466,
    467,
    468,
    469,
    470,
    471,
    472,
    473,
    474,
    475,
    476,
    477,
    478,
    479,
    480,
    481,
    482,
    483,
    484,
    485,
    486,
    487,
    488,
    489,
    490,
    491,
    492,
    493,
    494,
    495,
    497,
    498,
    500,
    501,
    502,
    503,
    504,
    505,
    506,
    507,
    508,
    509,
    510,
    511,
    512,
    513,
    514,
    515,
    516,
    517,
    518,
    520,
    521,
    523,
    524,
    525,
    527,
    529,
    530,
]

# The map between standard space group and hall numbers
pyxtal_hall_numbers = [
    1,
    2,
    3,
    6,
    9,
    18,
    21,
    30,
    39,
    57,
    60,
    63,
    72,
    81,
    90,
    108,
    109,
    112,
    115,
    116,
    119,
    122,
    123,
    124,
    125,
    128,
    134,
    137,
    143,
    149,
    155,
    161,
    164,
    170,
    173,
    176,
    182,
    185,
    191,
    197,
    203,
    209,
    212,
    215,
    218,
    221,
    227,
    229,
    230,
    234,
    239,
    245,
    251,
    257,
    263,
    266,
    269,
    275,
    279,
    284,
    290,
    292,
    298,
    304,
    310,
    313,
    316,
    323,
    334,
    336,
    337,
    338,
    341,
    343,
    349,
    350,
    351,
    352,
    353,
    354,
    355,
    356,
    357,
    358,
    360,
    362,
    363,
    365,
    366,
    367,
    368,
    369,
    370,
    371,
    372,
    373,
    374,
    375,
    376,
    377,
    378,
    379,
    380,
    381,
    382,
    383,
    384,
    385,
    386,
    387,
    388,
    389,
    390,
    391,
    392,
    393,
    394,
    395,
    396,
    397,
    398,
    399,
    400,
    401,
    403,
    405,
    406,
    407,
    409,
    411,
    412,
    413,
    415,
    417,
    418,
    419,
    421,
    423,
    424,
    425,
    427,
    429,
    430,
    431,
    432,
    433,
    435,
    436,
    438,
    439,
    440,
    441,
    442,
    443,
    444,
    446,
    447,
    448,
    449,
    450,
    452,
    454,
    455,
    456,
    457,
    458,
    460,
    462,
    463,
    464,
    465,
    466,
    467,
    468,
    469,
    470,
    471,
    472,
    473,
    474,
    475,
    476,
    477,
    478,
    479,
    480,
    481,
    482,
    483,
    484,
    485,
    486,
    487,
    488,
    489,
    490,
    491,
    492,
    493,
    494,
    496,
    497,
    499,
    500,
    501,
    502,
    503,
    504,
    505,
    506,
    507,
    508,
    509,
    510,
    511,
    512,
    513,
    514,
    515,
    516,
    517,
    519,
    520,
    522,
    523,
    524,
    526,
    528,
    529,
    530,
]


# --------------------------- Hall class -----------------------------
class Hall:
    """
    Class for conversion between Hall and standard spacegroups
    http://cci.lbl.gov/sginfo/itvb_2001_table_a1427_hall_symbols.html

    Args:
        spg_num: interger number between 1 and 230
        style: spglib or pyxtal
        permutation: allow permutation or not
    """

    def __init__(self, spgnum, style="pyxtal", permutation=False):
        self.spg = spgnum
        if style == "pyxtal":
            self.hall_default = pyxtal_hall_numbers[spgnum - 1]
        else:
            self.hall_default = spglib_hall_numbers[spgnum - 1]
        self.hall_numbers = []
        self.hall_symbols = []
        self.Ps = []  # convertion from standard
        self.P1s = []  # inverse convertion to standard
        for id in range(len(HALL_TABLE["Hall"])):
            if HALL_TABLE["Spg_num"][id] == spgnum:
                include = True if permutation else HALL_TABLE["Permutation"][id] == 0
                if include:
                    self.hall_numbers.append(HALL_TABLE["Hall"][id])
                    self.hall_symbols.append(HALL_TABLE["Symbol"][id])
                    self.Ps.append(abc2matrix(HALL_TABLE["P"][id]))
                    self.P1s.append(abc2matrix(HALL_TABLE["P^-1"][id]))
            elif HALL_TABLE["Spg_num"][id] > spgnum:
                break
        if len(self.hall_numbers) == 0:
            msg = "hall numbers cannot be found, check input " + spgnum
            raise RuntimeError(msg)


# --------------------------- Group class -----------------------------
class Group:
    """
    Class for storing a set of Wyckoff positions for a symmetry group.
    See the documentation for details about settings.

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

    one can access data such as `symbol`, `number` and `Wyckoff_positions`

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
    ([], [], [])
    >>> g.list_wyckoff_combinations([4, 8])
    ([[['4a'], ['8c']],
    [['4a'], ['8d']],
    [['4a'], ['8e']],
    [['4a'], ['8f']],
    [['4b'], ['8c']],
    [['4b'], ['8d']],
    [['4b'], ['8e']],
    [['4b'], ['8f']]],
    [False, True, True, True, False, True, True, True],
    [[[6], [4]], [[6], [3]], [[6], [2]], [[6], [1]],
    [[5], [4]], [[5], [3]], [[5], [2]], [[5], [1]]]
    )

    or search the subgroup information

    >>> g.get_max_t_subgroup()['subgroup']
    [12, 14, 15, 20, 36, 39, 41]

    or check if a given composition is compatible with Wyckoff positions

    >>> g = Group(225)
    >>> g.check_compatible([64, 28, 24])
    (True, True)

    or check the possible transition paths to a given supergroup

    >>> g = Group(59)
    >>> g.search_supergroup_paths(139, 2)
    [[71, 139], [129, 139], [137, 139]]


    Args:
        group: the group symbol or international number
        dim (defult: 3): the periodic dimension of the group
        use_hall (default: False): whether or not use the hall number
        style (default: `pyxtal`): the choice of hall number ('pyxtal'/'spglib')
        quick (defaut: False): whether or not ignore the wyckoff information
    """

    def __init__(self, group, dim=3, use_hall=False, style="pyxtal", quick=False):
        self.string = None
        self.dim = dim
        names = ["Point", "Rod", "Layer", "Space"]
        self.header = "-- " + names[dim] + "group --"

        # Retrieve symbol and number for the group (avoid redundancy)
        if not use_hall:
            self.symbol, self.number = get_symbol_and_number(group, dim)
        else:
            self.symbol = HALL_TABLE["Symbol"][group - 1]
            self.number = HALL_TABLE["Spg_num"][group - 1]

        self.PBC, self.lattice_type = get_pbc_and_lattice(self.number, dim)

        if dim == 3:
            results = get_point_group(self.number)
            self.point_group, self.pg_number, self.polar, self.inversion, self.chiral = results

        # Lazy load Wyckoff positions and hall data unless quick=True
        if not quick:
            self._initialize_hall_data(group, use_hall, style, dim)
            self._initialize_wyckoff_data(dim)

    def _initialize_hall_data(self, group, use_hall, style, dim):
        """Initialize hall number and transformation matrices."""
        if dim == 3:
            if not use_hall:
                if style == "pyxtal":
                    self.hall_number = pyxtal_hall_numbers[self.number - 1]
                else:
                    self.hall_number = spglib_hall_numbers[self.number - 1]
            else:
                self.hall_number = group
            self.P = abc2matrix(HALL_TABLE["P"][self.hall_number - 1])
            self.P1 = abc2matrix(HALL_TABLE["P^-1"][self.hall_number - 1])
        else:
            self.hall_number, self.P, self.P1 = None, None, None

    def _initialize_wyckoff_data(self, dim):
        """Initialize Wyckoff positions and organize them."""
        # Wyckoff positions, site_symmetry, generator
        self.wyckoffs = get_wyckoffs(self.number, dim=dim)
        self.w_symm = get_wyckoff_symmetry(self.number, dim=dim)

        # Create dicts with relevant Wyckoff position data lazily
        wpdicts_gen = [
            {
                "index": i,
                "letter": letter_from_index(i, self.wyckoffs, dim=self.dim),
                "ops": self.wyckoffs[i],
                "multiplicity": len(self.wyckoffs[i]),
                "symmetry": self.w_symm[i],
                "PBC": self.PBC,
                "dim": self.dim,
                "number": self.number,
                "symbol": self.symbol,
                "P": self.P,
                "P1": self.P1,
                "hall_number": self.hall_number,
                "directions": self.get_symmetry_directions(),
                "lattice_type": self.lattice_type,
            }
            for i in range(len(self.wyckoffs))
        ]

        # A list of Wyckoff_positions sorted by descending multiplicity
        #self.Wyckoff_positions = []
        #for wpdict in wpdicts:
        #    wp = Wyckoff_position.from_dict(wpdict)
        #    self.Wyckoff_positions.append(wp)

        # Use a generator to avoid keeping the full list of dicts in memory
        self.Wyckoff_positions = [
            Wyckoff_position.from_dict(wpdict) for wpdict in wpdicts_gen
        ]

        # Organize wyckoffs by multiplicity
        self.wyckoffs_organized = organized_wyckoffs(self)

    def __str__(self):
        if self.string is not None:
            return self.string
        else:
            s = self.header
            s += "# " + str(self.number) + " (" + self.symbol + ")--"

            for wp in self.Wyckoff_positions:
                s += "\n" + wp.get_label()
                if not hasattr(wp, "site_symm"):
                    wp.get_site_symmetry()
                s += "\tsite symm: " + wp.site_symm
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
        """Compute the degree of freedom for each site.

        Args:
            sites (list): List of site labels, e.g. ['4a', '8b'] or ['a', 'b']

        Returns:
            numpy.ndarray: Array of degrees of freedom for each site
        """

        def extract_dof(site):
            site_letter = site[-1] if len(site) > 1 else site
            id = len(self) - letters.index(site_letter) - 1
            xyz_str = self[id].ops[0].as_xyz_str()
            return len(set(re.sub(r"[^a-z]+", "", xyz_str)))

        return np.array([extract_dof(site) for site in sites])

    def get_lattice_dof(self):
        """
        compute the degree of freedom for the lattice
        """
        if self.lattice_type in ["triclinic"]:
            dof = 6
        elif self.lattice_type in ["monoclinic"]:
            dof = 4
        elif self.lattice_type in ["orthorhombic"]:
            dof = 3
        elif self.lattice_type in ["tetragonal", "hexagonal", "trigonal"]:
            dof = 2
        else:
            dof = 1

        return dof

    def is_valid_combination(self, sites):
        """
        check if the solutions are valid. A special WP with zero freedom (0,0,0)
        cannot be occupied twice.

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
            if sites.count(letter) > 1:
                freedom = np.trace(wp.ops[0].rotation_matrix) > 0
                if not freedom:
                    return False
        return True

    def list_wyckoff_combinations(self, numIons, quick=False, numWp=(None, None), Nmax=10000000):
        """
        List all possible wyckoff combinations for the given formula. Note this
        is really designed for a light weight calculation. If the solution space
        is big, set quick as True.

        Args:
            numIons (list): [12, 8]
            quick (Boolean): quickly generate some solutions
            numWp (tuple): (min_wp, max_wp)
            Nmax: maximumly allowed combinations

        Returns:
            Combinations: list of possible sites
            has_freedom: list of boolean numbers
            indices: list of wp indices
        """

        numIons = np.array(numIons)
        (min_wp, max_wp) = numWp
        # Must be greater than the number of smallest wp multiplicity
        if numIons.min() < self[-1].multiplicity or max_wp is not None and sum(numIons) > self[0].multiplicity * max_wp:
            return [], [], []

        basis = []  # [8, 4, 4]
        letters = []  # ['c', 'b', 'a']
        freedoms = []  # [False, False, False]
        ids = []  # [2, 3, 4]
        # obtain the basis
        for i, wp in enumerate(self):
            mul = wp.multiplicity
            letter = wp.letter
            freedom = np.trace(wp.ops[0].rotation_matrix) > 0
            if mul <= max(numIons):
                if quick:
                    if mul in basis and freedom:
                        pass
                    # elif mul in basis and basis.count(mul) >= 3:
                    #    pass
                    else:
                        basis.append(mul)
                        letters.append(letter)
                        freedoms.append(freedom)
                else:
                    basis.append(mul)
                    letters.append(letter)
                    freedoms.append(freedom)
                    ids.append(i)

        basis = np.array(basis)

        # quickly exit
        if np.min(numIons) < np.min(basis):
            # print(numIons, basis)
            return [], []
        # odd and even
        elif np.mod(numIons, 2).sum() > 0 and np.mod(basis, 2).sum() == 0:
            # print("odd-even", numIons, basis)
            # return None, False
            return [], [], []

        # print("basis", basis)
        # print("numIons", numIons)
        # obtain the maximum numbers for each basis
        # reset the maximum to 1 if there is no freedom
        # find the integer solutions
        # reset solutions according to max_wp
        max_solutions = np.floor(numIons[:, None] / basis)
        for i in range(len(freedoms)):
            if not freedoms[i]:
                max_solutions[:, i] = 1

        if max_wp is not None:
            N_max = max_wp - (len(numIons) - 1)
            max_solutions[max_solutions > N_max] = N_max

        list_solutions = []
        for i, numIon in enumerate(numIons):
            lists = []
            prod = 1
            for a in max_solutions[i]:
                if prod <= Nmax:  # 10000000:
                    d = int(a) + 1
                    lists.append(list(range(d)))
                    prod *= d
                else:
                    # If the size is too big, we terminate it asap
                    lists.append([0])
                    # Terminate the list
                    # break
            # print(len(lists), prod)

            sub_solutions = np.array(list(itertools.product(*lists)))
            N = sub_solutions.dot(basis)
            sub_solutions = sub_solutions[numIon == N]
            list_solutions.append(sub_solutions.tolist())
            # print(i)
            # print(sub_solutions)#; import sys; sys.exit()
            if len(sub_solutions) == 0:
                return [], [], []

        # Gather all solutions and remove very large number solutions
        solutions = np.array(list(itertools.product(*list_solutions)))
        dim1 = solutions.shape[0]
        dim2 = np.prod(solutions.shape[1:])
        solutions = solutions.reshape([dim1, dim2])
        if max_wp is not None:
            solutions_total = solutions.sum(axis=1)
            solutions = solutions[solutions_total <= max_wp]
        if min_wp is not None:
            solutions_total = solutions.sum(axis=1)
            solutions = solutions[solutions_total >= min_wp]

        # convert the results to list
        combinations = []
        has_freedom = []
        indices = []
        for solution in solutions:
            res = solution.reshape([len(numIons), len(basis)])
            _com = []
            _free = []
            _ids = []

            # QZ: check what's going on
            for i, _ in enumerate(numIons):
                tmp = []
                bad_resolution = False
                frozen = []
                ids_in = []
                for j, b in enumerate(basis):
                    if not freedoms[j] and (res[:, j]).sum() > 1:
                        bad_resolution = True
                        break

                    if res[i, j] > 0:
                        symbols = [str(b) + letters[j]] * res[i, j]
                        tmp.extend(symbols)
                        frozen.extend([freedoms[j]])
                        ids_in.extend([ids[j]] * res[i, j])

                if not bad_resolution:
                    _com.append(tmp)
                    _free.extend(frozen)
                    _ids.append(ids_in)

            if len(_com) == len(numIons):
                combinations.append(_com)
                indices.append(_ids)
                if True in _free:
                    has_freedom.append(True)
                else:
                    has_freedom.append(False)

        return combinations, has_freedom, indices

    def get_spg_symmetry_object(self):
        """
        Generate the symmetry table for the given space group
        It only supports space group now!
        """

        if self.dim == 3:
            l_type, bravis = self.lattice_type, self.symbol[0]
            wp = self.get_wyckoff_position(0)

            ops = wp.get_euclidean_ops() if 143 <= self.number <= 194 else wp.ops

            if bravis in ["A", "B", "C", "I"]:
                ops = ops[: int(len(ops) / 2)]
            elif bravis == "R":
                ops = ops[: int(len(ops) / 3)]
            elif bravis == "F":
                ops = ops[: int(len(ops) / 4)]

            return site_symmetry(ops, l_type, bravis, True)
            # ss.to_beautiful_matrix_representation()
        else:
            raise ValueError("Only supports space group symmetry")

    def get_wyckoff_position(self, index):
        """
        Returns a single Wyckoff_position object.

        Args:
            index: the index of the Wyckoff position within the group.

        Returns: a Wyckoff_position object
        """
        if isinstance(index, str):
            # Extract letter from number-letter combinations ("4d"->"d")
            for c in index:
                if c.isalpha():
                    letter = c
                    break
            index = index_from_letter(letter, self.wyckoffs, dim=self.dim)
        return self.Wyckoff_positions[index]

    def get_wyckoff_position_from_xyz(self, xyz, decimals=4):
        """
        Returns a single Wyckoff_position object.

        Args:
            xyz: a trial [x, y, z] coordinate

        Returns: a Wyckoff_position object
        """
        xyz = np.round(np.array(xyz, dtype=float), decimals=decimals)
        xyz -= np.floor(xyz)

        for wp in self.Wyckoff_positions:
            pos = wp.apply_ops(xyz)
            pos -= np.floor(pos)
            is_present = np.any(np.all(pos == xyz, axis=1))
            if is_present and len(pos) == len(np.unique(pos, axis=0)):
                return wp

        print("Cannot find the suitable wp for the given input")
        return None

    def get_alternatives(self):
        """
        Get the alternative settings as a dictionary
        """
        if self.dim == 3:
            wyc_sets = loadfn(rf("pyxtal", "database/wyckoff_sets.json"))
            return wyc_sets[str(self.number)]
        else:
            msg = "Only supports the subgroups for space group"
            raise NotImplementedError(msg)

    @classmethod
    def _get_max_k_subgroup(cls, number=None):
        """
        Returns the maximal k-subgroups as a dictionary
        """
        if number is None:
            number = cls.number
        k_subgroup = SYMDATA.get_k_subgroup()
        return k_subgroup[str(number)]

    @classmethod
    def _get_max_t_subgroup(cls, number=None):
        """
        Returns the maximal t-subgroups as a dictionary
        """
        if number is None:
            number = cls.number
        t_subgroup = SYMDATA.get_t_subgroup()
        return t_subgroup[str(number)]

    def get_max_k_subgroup(self):

        return self._get_max_k_subgroup(self.number)

    def get_max_t_subgroup(self):

        return self._get_max_t_subgroup(self.number)

    def get_max_subgroup(self, H):
        """
        Returns the dicts for both t and k subgroup according the to
        trail group H

        Args:
            H (int): 1-230
        """
        #if self.point_group == Group(H, quick=True).point_group:
        if self.point_group == get_point_group(H):
            g_type = "k"
            dicts = self.get_max_k_subgroup()
        else:
            g_type = "t"
            dicts = self.get_max_t_subgroup()
        return dicts, g_type

    def get_wp_list(self, reverse=False):
        """
        Get the reversed list of wps
        """
        # wp_list = [(str(x.multiplicity)+x.letter) for x in self.Wyckoff_positions]
        wp_list = [(x.get_label()) for x in self.Wyckoff_positions]
        if reverse:
            wp_list.reverse()
        return wp_list

    def get_splitters_from_structure(self, struc, group_type="t"):
        """
        Get the valid symmetry relations for a structure to its supergroup
        e.g.,

        Args:
            - struc: pyxtal structure
            - group_type: `t` or `k`
        Returns:
            list of valid transitions [(id, (['4a'], ['4b'], [['4a'], ['4c']])]
        """
        if group_type == "t" :
            dicts = self.get_max_t_subgroup()
        else:
            dicts = self.get_max_k_subgroup()

        # search for the compatible solutions
        solutions = []
        for i, sub in enumerate(dicts["subgroup"]):
            if sub == struc.group.number:
                # extract the relation
                relation = dicts["relations"][i]
                trans = np.linalg.inv(dicts["transformation"][i][:, :3])
                if struc.lattice.check_mismatch(trans, self.lattice_type):
                    results = self.get_splitters_from_relation(struc, relation)
                    if results is not None:
                        sols = list(itertools.product(*results))
                        trials = self.get_valid_solutions(sols)
                        solutions.append((i, trials))
        return solutions

    def get_splitters_from_relation(self, struc, relation):
        """
        Get the valid symmetry relations for a structure to its supergroup
        e.g.,

        Args:
            - struc: pyxtal structure
            - group_type: `t` or `k`
        Returns:
            list of valid transitions
        """

        elements, sites = struc._get_elements_and_sites()
        wp_list = self.get_wp_list(reverse=True)

        good_splittings_list = []

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
                total_units += int(x[:-1]) * _site_counts[j]

            # collect all possible supergroup transitions
            for j, split in enumerate(relation):
                if np.all([x in site for x in split]):
                    wp_indices.append(j)

            wps = [wp_list[x] for x in wp_indices]
            blocks = [np.array([relation[j].count(s) for s in _site])
                      for j in wp_indices]
            block_units = [sum([int(x[:-1]) * block[j]
                               for j, x in enumerate(_site)]) for block in blocks]

            # below is a brute force search for the valid combinations
            combo_storage = [np.zeros(len(block_units))]
            good_list = []
            while len(combo_storage) > 0:
                holder = []
                for x in combo_storage:
                    for k in range(len(block_units)):
                        # trial solution
                        trial = np.array(deepcopy(x), dtype=int)
                        trial[k] += 1
                        if trial.tolist() in holder:
                            continue
                        sum_units = np.dot(trial, block_units)
                        if sum_units > total_units:
                            continue

                        if sum_units < total_units:
                            holder.append(trial.tolist())
                        else:
                            tester = np.zeros(len(_site_counts))
                            for l, z in enumerate(trial):
                                tester += z * blocks[l]
                            if np.all(tester == _site_counts):
                                G_sites = [wp for wp, num in zip(
                                    wps, trial) for _ in range(max(num, 1)) if num > 0]
                                if G_sites not in good_list:
                                    good_list.append(G_sites)
                combo_storage = holder

            if len(good_list) == 0:
                return None
            else:
                good_splittings_list.append(good_list)
        return good_splittings_list

    def get_min_supergroup(self, group_type="t", G=None):
        """
        Returns the minimal supergroups as a dictionary
        """
        if self.dim == 3:
            dicts = {
                "supergroup": [],
                "transformation": [],
                "relations": [],
                "idx": [],
            }
            sgs = range(1, 231) if G is None else G

            for sg in sgs:
                subgroups = None
                if group_type == "t":
                    if sg > self.number:
                        subgroups = Group._get_max_t_subgroup(sg)
                else:
                    #g1 = Group(sg)
                    if self.point_group == get_point_group(sg):
                        subgroups = Group._get_max_k_subgroup(sg)
                if subgroups is not None:
                    for i, sub in enumerate(subgroups["subgroup"]):
                        if sub == self.number:
                            trans = subgroups["transformation"][i]
                            relation = subgroups["relations"][i]
                            dicts["supergroup"].append(sg)
                            dicts["transformation"].append(trans)
                            dicts["relations"].append(relation)
                            dicts["idx"].append(i)
            return dicts
        else:
            msg = "Only supports the supergroups for space group"
            raise NotImplementedError(msg)

    @classmethod
    def _get_max_subgroup_numbers(cls, number, max_cell=9):
        """
        Returns the minimal supergroups as a dictionary
        """
        groups = []
        sub_k = Group._get_max_k_subgroup(number)
        sub_t = Group._get_max_t_subgroup(number)
        k = sub_k["subgroup"]
        t = sub_t["subgroup"]
        for i, n in enumerate(t):
            if np.linalg.det(sub_t["transformation"][i][:3, :3]) <= max_cell:
                groups.append(n)
        for i, n in enumerate(k):
            if np.linalg.det(sub_k["transformation"][i][:3, :3]) <= max_cell:
                groups.append(n)
        return groups

    def get_max_subgroup_numbers(self, max_cell=9):
        """
        Returns the minimal supergroups as a dictionary
        """
        groups = []
        if self.dim == 3:
            sub_k = self.get_max_k_subgroup()
            sub_t = self.get_max_t_subgroup()
            k = sub_k["subgroup"]
            t = sub_t["subgroup"]
            for i, n in enumerate(t):
                if np.linalg.det(sub_t["transformation"][i][:3, :3]) <= max_cell:
                    groups.append(n)
            for i, n in enumerate(k):
                if np.linalg.det(sub_k["transformation"][i][:3, :3]) <= max_cell:
                    groups.append(n)
            return groups
        else:
            msg = "Only supports the subgroups for space group"
            raise NotImplementedError(msg)

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
            # check the freedom
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
                # check the freedom
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
        has_freedom = False  # whether or not one degree of freedom exists
        used_indices = []  # wp's already used that don't have any freedom

        # Loop over species
        # Sort the specie from low to high so that the solution can be found ealier
        for id in np.argsort(numIons):
            # Get lists of multiplicity, maxn and freedom
            numIon = numIons[id]
            if valid_orientations is None:
                l_mult, l_maxn, l_free, indices = self.get_lists(
                    numIon, used_indices)
            else:
                vo = valid_orientations[id]
                l_mult, l_maxn, l_free, indices = self.get_lists_mol(
                    numIon, used_indices, vo)
            # print(numIon, l_mult, indices, l_maxn, l_free)

            # Loop over possible combinations
            p = 0  # Create pointer variable to move through lists
            # Store the number of each WP, used across possible WP combinations
            n0 = [0] * len(l_mult)
            n = deepcopy(n0)
            for i, _mult in enumerate(l_mult):
                if l_maxn[i] != 0:
                    p = i
                    n[i] = l_maxn[i]
                    break
            p2 = p
            if n == n0:
                return False, False

            # print(numIon, n, n0, p)
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
                    # print('All combinations failed', numIon, n, n0)
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
            # print('used_indices', used_indices)
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
        layers[0] = {
            "groups": [H],
            "subgroups": [],
        }
        final = []
        traversed = []

        # Searches for every subgroup of the the groups from the previous layer.
        # stores the possible groups of each layer and their subgroups
        # in a dictinoary to avoid redundant calculations.
        for l in range(1, max_layer + 1):
            previous_layer_groups = layers[l - 1]["groups"]
            groups = []
            subgroups = []
            for g in previous_layer_groups:
                #subgroup_numbers = np.unique(Group(g, quick=True).get_max_subgroup_numbers())
                subgroup_numbers = np.unique(Group._get_max_subgroup_numbers(g))

                # If a subgroup list has been found with H
                # trace a path through the dictionary to build the path
                if self.number in subgroup_numbers:
                    paths = [[g]]
                    for j in reversed(range(l - 1)):
                        holder = []
                        for path in paths:
                            tail_number = path[-1]
                            indices = [
                                idx for idx, numbers in enumerate(layers[j]["subgroups"]) if tail_number in numbers
                            ]
                            holder.extend([path + [layers[j]["groups"][idx]] for idx in indices])  # noqa: RUF005
                        paths = deepcopy(holder)
                    final.extend(paths)
                    subgroups.append([])

                # Continue to generate next layer if the path to H has not been found.
                else:
                    subgroups.append(subgroup_numbers)
                    groups.extend(
                        [x for x in subgroup_numbers if x not in groups and x not in traversed])

            traversed.extend(groups)
            layers[l] = {"groups": deepcopy(groups), "subgroups": []}
            layers[l - 1]["subgroups"] = deepcopy(subgroups)
        return final

    #def path_to_subgroup(self, H):
    #    """
    #    For a given a path, extract the
    #        a list of (g_types, subgroup_id, spg_number, wp_list (optional))
    #    """
    #    path_list = []
    #    paths = self.search_subgroup_paths(H)
    #    if len(paths) > 0:
    #        path = paths[0]
    #        g0 = Group(path[0], quick=True)
    #        for p in path[1:]:
    #            #g1 = Group(p, quick=True)
    #            #pg = get_point_group(p)
    #            if g0.point_group == get_point_group(p):
    #                g_type = "k"
    #                spgs = g0.get_max_k_subgroup()["subgroup"]
    #            else:
    #                g_type = "t"
    #                spgs = g0.get_max_t_subgroup()["subgroup"]
    #            for spg in spgs:
    #                if spg == p:
    #                    break
    #            path_list.append((g_type, id, p))
    #            g0 = g1
    #    return path_list

    def path_to_subgroup(self, H):
        """
        For a given a path, extract the
            a list of (g_types, subgroup_id, spg_number, wp_list (optional))
        """
        path_list = []
        paths = self.search_subgroup_paths(H)
        if len(paths) > 0:
            path = paths[0]
            sg0 = path[0]
            pg0 = get_point_group(path[0])
            #pg0 = Group(path[0], quick=True)
            for p in path[1:]:
                pg1 = get_point_group(p)
                if pg0 == pg1:
                    g_type = "k"
                    spgs = Group._get_max_k_subgroup(sg0)["subgroup"]
                else:
                    g_type = "t"
                    spgs = Group._get_max_t_subgroup(sg0)["subgroup"]
                for spg in spgs:
                    if spg == p:
                        break
                path_list.append((g_type, id, p))
                pg0 = pg1
        return path_list

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
        tmp = Group(G, quick=True)
        paths = tmp.search_supergroup_paths(self.number, max_layer=max_layer)
        for p in paths:
            p.reverse()
            p.append(G)
        return paths

    def add_k_transitions(self, path, n=1):
        """
        Adds additional k transitions to a subgroup path. ONLY n = 1 is
        supported for now. It will return viable additions in front of each
        group in the path.

        Args:
            path: a single result of search_subgroup_paths function
            n: number of extra k transitions to add to the given path
        Returns:
            a list of maximal subgroup chains with extra k type transitions
        """
        if n != 1:
            print("only 1 extra k type supported at this time")
            return None
        k_subgroup = SYMDATA.get_k_subgroup()
        t_subgroup = SYMDATA.get_t_subgroup()

        solutions = []
        for i in range(len(path[:-1])):
            g = path[i]
            h = path[i + 1]
            options = set(k_subgroup[str(g)]["subgroup"] +
                          t_subgroup[str(g)]["subgroup"])
            # print(g, h, options)
            for _g in options:
                ls = k_subgroup[str(_g)]["subgroup"] + \
                    t_subgroup[str(_g)]["subgroup"]
                if h in ls:
                    sol = deepcopy(path)
                    sol.insert(i + 1, _g)
                    solutions.append(sol)
        # https://stackoverflow.com/questions/2213923/removing-duplicates-from-a-list-of-lists
        solutions.sort()
        return [k for k, _ in itertools.groupby(solutions)]

    def path_to_general_wp(self, index=1, max_steps=1):
        """
        Find the path to transform the special wp into general site

        Args:
            index: the index of starting wp
            max_steps: the number of steps to search

        Return:
            a list of (g_types, subgroup_id, spg_number, wp_list (optional))
        """
        # label = [str(self[index].multiplicity) + self[index].letter]
        k_subgroup = SYMDATA.get_k_subgroup()
        t_subgroup = SYMDATA.get_t_subgroup()

        label = [self[index].get_label()]
        potential = [[(None, None, self.number, label)]]
        solutions = []

        for _step in range(max_steps):
            _potential = []
            for p in potential:
                tail = p[-1]
                tdict = t_subgroup[str(tail[2])]
                len_t = len(tdict["subgroup"])
                kdict = k_subgroup[str(tail[2])]
                len_k = len(kdict["subgroup"])
                _indexs = [ord(x[-1]) - 97 for x in tail[3]]
                next_steps = [
                    [
                        (
                            "t",
                            i,
                            tdict["subgroup"][i],
                            functools.reduce(
                                operator.iadd, [tdict["relations"][i][_index] for _index in _indexs], []),
                        )
                    ]
                    for i in range(len_t)
                ] + [
                    [
                        (
                            "k",
                            i,
                            kdict["subgroup"][i],
                            functools.reduce(
                                operator.iadd, [kdict["relations"][i][_index] for _index in _indexs], []),
                        )
                    ]
                    for i in range(len_k)
                    if kdict["subgroup"][i] != tail[2]
                ]
                _potential.extend([deepcopy(p) + n for n in next_steps])
            potential = _potential

            solutions.extend(
                [
                    deepcopy(p)[1:]
                    for p in potential
                    if (len(set(p[-1][3])) == 1 and p[-1][3][0][-1] == Wyckoff_position.from_group_and_index(p[-1][2], 0).letter)
                ]
            )
            potential = [
                p for p in potential if not (len(set(p[-1][3])) == 1 and p[-1][3][0][-1] == Wyckoff_position.from_group_and_index(p[-1][2], 0).letter)
            ]

        return solutions

    def short_path_to_general_wp(self, index=1, t_only=False):
        """
        Find a short path to turn the spcical wp to general position

        Args:
            index: index of the wp
            t_only: only consider t_spliting
        """
        for i in range(1, 5):
            paths = self.path_to_general_wp(index, max_steps=i)
            if len(paths) > 0:
                last_gs = np.array([p[-1][2] for p in paths])
                if t_only:
                    last_gs[last_gs > len(self[0])] = 0
                max_id = np.argmax(last_gs)
                return paths[max_id]
        return None

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
        contrast to the primitive cell). Based on the type of cell centering
        (P, A, C, I, R, or F)
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
            if self.symbol[0] == "P":
                return 1  # P
            elif self.symbol[0] == "R":
                return 3  # R
            elif self.symbol[0] == "F":
                return 4  # F
            else:
                return 2  # A, C, I

    def get_free_axis(self):
        """
        Get the free axis that can perform continus translation
        """
        number = self.number

        if number == 1:
            return [0, 1, 2]
        elif number == 2:
            return []
        elif 3 <= number <= 5:
            return [1]  # '2'
        elif 6 <= number <= 9:
            return [0, 2]  # 'm'
        elif 10 <= number <= 24:
            return []  # '2/m', '222'
        elif 25 <= number <= 46:
            return [2]  # 'mm2'
        elif 47 <= number <= 74:
            return []  # 'mmm'
        elif 75 <= number <= 80:
            return [2]  # '4'
        elif 81 <= number <= 98:
            return []  # '-4', '4/m', '422'
        elif 99 <= number <= 110:
            return [2]  # '4mm'
        elif 111 <= number <= 142:
            return []  # '-42m', '4/mmm'
        elif 143 <= number <= 146:
            return [2]  # '3'
        elif 147 <= number <= 155:
            return []  # '-3', '32'
        elif 156 <= number <= 161:
            return [2]  # '3m'
        elif 162 <= number <= 167:
            return []  # '-3m'
        elif 168 <= number <= 173:
            return [2]  # '6'
        elif 174 <= number <= 182:
            return []  # '-6', '6/m', '622'
        elif 183 <= number <= 186:
            return [2]  # '6mm'
        elif 187 <= number <= 194:
            return []  # '-62m', '6/mmm'
        elif 195 <= number <= 230:
            return []  # '23', 'm-3', '432', '-43m', 'm-3m',
        return None

    @classmethod
    def list_groups(cls, dim=3):
        """
        Function for quick print of groups and symbols.

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

        group_symbols = loadfn(rf("pyxtal", "database/symbols.json"))
        data = group_symbols[keys[dim]]
        index = range(1, len(data) + 1)
        df = pd.DataFrame(index=index, data=data, columns=[keys[dim]])
        pd.set_option("display.max_rows", len(df))
        print(df)

    def get_index_by_letter(self, letter):
        """
        get the wp object by the letter
        """
        if len(letter) > 1:
            letter = letter[-1]
        # print(letter); print(letters.index(letter))
        return len(self) - letters.index(letter) - 1

    def get_wp_by_letter(self, letter):
        """
        get the wp object by the letter
        """
        return self[self.get_index_by_letter(letter)]

    def get_symmetry_directions(self):
        """
        Table 2.1.3.1 from International Tables for Crystallography (2016).
        Vol. A, Chapter 2.1, pp. 142-174.
        including Primary, secondary and Tertiary
        """
        # if self.dim == 3:
        return get_symmetry_directions(self.lattice_type, self.symbol[0])
        # else:
        #    print("Cannot handle dimension other than 3")
        #    pass


#
# ----------------------- Wyckoff Position class  ------------------------


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

    # =============================Initialization===========================
    @classmethod
    def from_dict(cls, dictionary):
        """
        Constructs a Wyckoff_position object using a dictionary.
        """
        wp = cls()
        for key in dictionary:
            setattr(wp, key, dictionary[key])
        # wp.get_site_symmetry()
        wp.set_euclidean()
        # For nonstandard setting only
        if wp.P1 is not None and not identity_ops(wp.P1):
            wp.set_generators()
            wp.set_ops()
        return wp

    @classmethod
    def from_group_and_letter(cls, group, letter, dim=3, style="pyxtal", hn=None):
        """
        Creates a Wyckoff_position using the space group number and index

        Args:
            group: the international number of the symmetry group
            letter: e.g. 4a
            dim: the periodic dimension of the crystal
            style: 'pyxtal' or spglib, differing in the choice of origin
            hn: hall_number
        """

        for c in letter:
            if c.isalpha():
                letter = c
                break
        ops_all = get_wyckoffs(group, dim=dim)
        index = index_from_letter(letter, ops_all, dim=dim)
        if hn is not None:
            wp = cls.from_group_and_index(
                hn, index, dim, use_hall=True, wyckoffs=ops_all)
        else:
            wp = cls.from_group_and_index(
                group, index, dim, style=style, wyckoffs=ops_all)
        return wp

    @classmethod
    def from_group_and_index(cls, group, index, dim=3, use_hall=False, style="pyxtal", wyckoffs=None):
        """
        Creates a Wyckoff_position using the space group number and index

        Args:
            group: the international number of the symmetry group
            index: the index of the Wyckoff position within the group.
            dim: the periodic dimension of the crystal
            use_hall (default: False): whether or not use the hall number
            style (default: `pyxtal`): 'pyxtal' or 'spglib' for hall number
        """
        number, hall_number, P, P1 = group, None, None, None
        if not use_hall:
            symbol, number = get_symbol_and_number(group, dim)
        else:
            symbol = HALL_TABLE["Symbol"][group - 1]
            number = HALL_TABLE["Spg_num"][group - 1]
        pbc, lattice_type = get_pbc_and_lattice(number, dim)

        if dim == 3:
            PBC = [1, 1, 1]
            if not use_hall:
                if style == "pyxtal":
                    hall_number = pyxtal_hall_numbers[number - 1]
                else:
                    hall_number = spglib_hall_numbers[number - 1]
                    P =  abc2matrix(HALL_TABLE["P"][hall_number - 1])
                    P1 = abc2matrix(HALL_TABLE["P^-1"][hall_number - 1])
            else:
                hall_number = group
                P =  abc2matrix(HALL_TABLE["P"][hall_number - 1])
                P1 = abc2matrix(HALL_TABLE["P^-1"][hall_number - 1])
            directions = get_symmetry_directions(lattice_type, symbol[0])

        elif dim == 2:
            PBC = [1, 1, 0]
            directions = None
        elif dim == 1:
            PBC = [0, 0, 1]
            directions = None

        if wyckoffs is None:
            wyckoffs = get_wyckoffs(number, dim=dim)

        wpdict = {
            "index": index,
            "letter": letter_from_index(index, wyckoffs, dim=dim),
            "ops": wyckoffs[index],
            "multiplicity": len(wyckoffs[index]),
            "symmetry": get_wyckoff_symmetry(number, dim=dim)[index],
            # "generators": get_generators(number, dim=dim)[index],
            "PBC": PBC,
            "dim": dim,
            "number": number,
            "P": P,
            "P1": P1,
            "hall_number": hall_number,
            "symbol": symbol,
            "directions": directions,
            "lattice_type": lattice_type,
        }

        return cls.from_dict(wpdict)

    @classmethod
    def from_symops_wo_group(cls, ops):
        """
        search Wyckoff Position by symmetry operations
        Now only supports space group symmetry
        Assuming general position only

        Args:
            ops: a list of symmetry operations

        Returns:
            `Wyckoff_position`
        """
        _, spg_num = get_symmetry_from_ops(ops)
        wp = cls.from_group_and_index(spg_num, 0)
        if isinstance(ops[0], str):
            ops = [SymmOp.from_xyz_str(op) for op in ops]
        wp.ops = ops
        match_spg, match_hm = wp.update()
        # print("match_spg", match_spg, "match_hall", match_hm)
        return wp

    @classmethod
    def from_symops(cls, ops, G):
        """
        search Wyckoff Position by symmetry operations


        Args:
            ops: a list of symmetry operations
            G: the Group object

        Returns:
            `Wyckoff_position`

        """
        if isinstance(ops[0], str):
            ops = [SymmOp.from_xyz_str(op) for op in ops]

        for wp in G:
            if wp.has_equivalent_ops(ops):
                return wp

        if isinstance(ops[0], str):
            print(ops)
        else:
            for op in ops:
                print(op.as_xyz_str())
        raise RuntimeError("Cannot find the right wp")

    def from_index_quick(self, wyckoffs, index, P=None, P1=None):
        """
        A short cut to create the WP object from a given index
        ignore the site symmetry and generators
        Mainly used for the update function

        Args:
            wyckoffs: wyckoff position
            index: index of wp
            P: transformation matrix (rot + trans)
        """

        if P is None:
            P = self.P
            P1 = self.P1
        wpdict = {
            "index": index,
            "letter": letter_from_index(index, wyckoffs, dim=self.dim),
            "ops": wyckoffs[index],
            "multiplicity": len(wyckoffs[index]),
            "PBC": self.PBC,
            "dim": self.dim,
            "number": self.number,
            "P": P,
            "P1": P1,
            "hall_number": self.hall_number,
        }
        return Wyckoff_position.from_dict(wpdict)

    # =============================Fundamentals===========================
    def __str__(self, supress=False):
        if self.dim not in [0, 1, 2, 3]:
            return "invalid crystal dimension. Must be between 0 and 3."
        if not hasattr(self, "site_symm"):
            self.get_site_symmetry()
        s = "Wyckoff position " + self.get_label() + " in "
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
                s += "\n" + op.as_xyz_str()

        self.string = s
        return self.string

    def __repr__(self):
        return str(self)

    def __iter__(self):
        yield from self.ops

    def __getitem__(self, index):
        return self.ops[index]

    def __len__(self):
        return self.multiplicity

    def copy(self):
        """
        Simply copy the structure
        """
        return deepcopy(self)

    def save_dict(self):
        return {
            "group": self.number,
            "index": self.index,
            "dim": self.dim,
            # "transformation": self.get_transformation(),
        }

    @classmethod
    def load_dict(cls, dicts):
        g = dicts["group"]
        index = dicts["index"]
        dim = dicts["dim"]
        # , trans=trans)
        return Wyckoff_position.from_group_and_index(g, index, dim)

    # =============================Updates===========================
    def set_ops(self):
        self.ops = self.get_ops_from_transformation()

    def get_ops_from_transformation(self):
        """
        Get symmetry operation from the generators
        """
        # Get the position for the 1st site
        ops1 = []
        if self.index > 0:
            if self.P1 is not None and not identity_ops(self.P1):
                for op in self.ops:
                    rot_P = self.P[0].T
                    rot_Q = self.P1[0].T
                    tran_P = self.P[1]
                    # R = Q * R * P, suitable when P = {a-c, b, c}
                    tran = rot_Q.dot(op.translation_vector) - tran_P
                    rot = rot_Q.dot(op.rotation_matrix).dot(rot_P)
                    op0 = SymmOp.from_rotation_and_translation(rot, tran)
                    ops1.append(op0)
                    # print(op0.as_xyz_str())
                ops1 = trim_ops(ops1)
            else:
                op0 = self.ops[0]
        else:
            ops1 = self.generators
        return ops1

    def update(self):
        """
        update the spacegroup information if needed
        """
        match_spg, match_hall = False, False
        match_spg = self.update_index()
        if not match_spg:
            match_hall = self.update_hall()

        if not match_spg and not match_hall:
            print("match_spg", match_spg, "match_hall", match_hall)
            print(self)
            print(self.get_hm_symbol())
            raise RuntimeError("Cannot find the right hall_number")
        return match_spg, match_hall

    def update_hall(self, hall_numbers=None):
        """
        update the Hall number when the symmetry operation changes

        Args:
            hall_numbers: a list of numbers for consideration
        """
        # print("test", self)
        if hall_numbers is None:
            hall_numbers = Hall(self.number).hall_numbers

        candidates = self.process_ops()
        success = False
        for hall_number in hall_numbers:
            P =  abc2matrix(HALL_TABLE["P"][hall_number - 1])
            P1 = abc2matrix(HALL_TABLE["P^-1"][hall_number - 1])
            wyckoffs = get_wyckoffs(self.number, dim=self.dim)

            # Fist check the original index
            wp2 = self.from_index_quick(wyckoffs, self.index, P, P1)
            for ops in candidates:
                if wp2.has_equivalent_ops(ops):
                    success = True
                    # print("same letter") #; import sys; sys.exit()
                    break

            # Check other sites
            if not success:
                for i in range(len(wyckoffs)):
                    if i != self.index and len(wyckoffs[i]) == self.multiplicity:
                        wp2 = self.from_index_quick(wyckoffs, i, P, P1)
                        for ops in candidates:
                            if wp2.has_equivalent_ops(ops):
                                success = True
                                self.index = i
                                self.letter = wp2.letter
                                # print("new letter")
                                break
                    if success:
                        break
            if success:
                self.hall_number = hall_number
                self.P = wp2.P
                self.P1 = wp2.P1
                self.ops = wp2.ops

                return True
        return False

    def update_index(self):
        """
        check if needs to update the index due to lattice transformation
        """
        wyckoffs = get_wyckoffs(self.number, dim=self.dim)
        wp2 = self.from_index_quick(wyckoffs, self.index)
        if self.has_equivalent_ops(wp2):
            return True
        else:
            for i in range(len(wyckoffs)):
                if i != self.index and len(wyckoffs[i]) == self.multiplicity:
                    wp2 = self.from_index_quick(wyckoffs, i)
                    if self.has_equivalent_ops(wp2):
                        self.index = i
                        self.letter = wp2.letter
                        # adjust to normal
                        self.ops = wp2.ops
                        return True
        return False

    def transform_from_matrices(self, trans):
        """
        Args:
            trans: a list of transformation matrices
        """
        for tran in trans:
            self.transform_from_matrix(tran, False)

    def transform_from_matrix(self, trans=None, reset=True, update=False):
        """
        Transform the symmetry operation according to cell transformation
        Mostly needed when optimizing the lattice
        """

        if trans is None:
            if self.number in [7, 9, 13, 14, 15]:
                trans = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 1]])
            elif self.number in [5, 8, 9, 12]:
                trans = np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]])

        if 2 < self.number < 16:
            #ops = Group(self.number)[self.index] if reset else self.ops
            if reset:
                ops = Wyckoff_position.from_group_and_index(self.number, self.index)
            else:
                ops = self.ops

            for j, op in enumerate(ops):
                vec = op.translation_vector.dot(trans)
                vec -= np.floor(vec)
                op1 = op.from_rotation_and_translation(op.rotation_matrix, vec)
                self.ops[j] = op1

            if update:
                self.update_hall()

    def process_ops(self):
        """
        handle some annoying cases
        e.g., in I2, ['1/2, y, 1/2', '0, y+1/2, 0'] can be transfered to
        ['0, y, 0', '1/2, y+1/2, 1/2']
        """
        opss = [self.ops]
        if self.number in [5, 12] and self.index > 0:
            # replace y with y+1/2
            op2 = SymmOp.from_xyz_str("x, y+1/2, z")
            ops = [op2 * op for op in self.ops]
            opss.append(ops)

        if self.number in [13] and self.index > 0:
            op2 = SymmOp.from_xyz_str("x, -y, z")
            ops = [op2 * op for op in self.ops]
            opss.append(ops)

            # for op in ops: print('AAAA', op.as_xyz_str())
        return opss

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
                msg = f"Spg {self.number:d}, Invalid switch in Wyckoff Pos\n"
                msg += str(self)
                msg += "\n" + str(G[index])
                raise ValueError(msg)

            return G[index]
        return self

    # =============================Get functions===========================
    def get_site_symm_wo_translation(self):
        return [SymmOp.from_rotation_and_translation(op.rotation_matrix, [0, 0, 0]) for op in self.symmetry[0]]

    def get_site_symmetry_object(self, idx=0):
        ops = self.get_site_symm_ops(idx)#; print(self.number, self.index, self.letter)
        return site_symmetry(ops, self.lattice_type, self.symbol[0])

    def get_site_symmetry(self, idx=0):
        ss = self.get_site_symmetry_object(idx)
        # ss_string_from_ops(ops, self.number, dim=self.dim)
        self.site_symm = ss.name

    def get_site_symm_ops(self, idx=0):
        return self.get_euclidean_symmetries(idx) if self.euclidean else self.symmetry[idx]

    def get_hm_number(self, tol=1e-5):
        if self.index == 0:
            return get_symmetry_from_ops(self.ops, tol)[0]
        else:
            print(self)
            raise ValueError("input must be general position")

    def get_hm_symbol(self):
        """
        get Hermann-Mauguin symbol
        """
        return HALL_TABLE["Symbol"][self.hall_number - 1]

    def get_dof(self):
        """
        Simply return the degree of freedom
        """
        return np.linalg.matrix_rank(self.ops[0].rotation_matrix)

    def get_label(self):
        """
        get the string like 4a
        """
        return str(self.multiplicity) + self.letter

    def get_frozen_axis(self):
        if self.index == 0:
            return []
        elif self.get_dof() == 0:
            return [0, 1, 2]
        else:
            if self.number >= 75:
                # if self.ops[0].rotation_matrix[2,2] == 1:
                #    return [0, 1]
                # else:
                #    return [0, 1, 2]
                return [ax for ax in range(3) if self.ops[0].rotation_matrix[ax, ax] == 0]
            else:
                if self.get_dof() == 1:
                    if self.ops[0].rotation_matrix[2, 2] == 1:
                        return [0, 1]
                    elif self.ops[0].rotation_matrix[1, 1] == 1:
                        return [0, 2]
                    elif self.ops[0].rotation_matrix[0, 0] == 1:
                        return [1, 2]
                    return None
                else:
                    if self.ops[0].rotation_matrix[2, 2] != 1:
                        return [2]
                    elif self.ops[0].rotation_matrix[1, 1] != 1:
                        return [1]
                    elif self.ops[0].rotation_matrix[0, 0] != 1:
                        return [0]
                    return None

    def get_euclidean_symmetries(self, idx=0):
        """
        return the symmetry operation object at the Euclidean space

        Returns:
            list of pymatgen SymmOp object
        """
        if idx >= len(self.symmetry):
            raise ValueError(
                f"Cannot pick {idx:d} in {len(self.symmetry):d} operations")
        ops = []
        for op in self.symmetry[idx]:
            hat = SymmOp.from_rotation_and_translation(hex_cell, [0, 0, 0])
            ops.append(hat * op * hat.inverse)
        return ops

    def get_euclidean_ops(self):
        """
        return the symmetry operation object at the Euclidean space

        Returns:
            list of pymatgen SymmOp object
        """
        ops = [None] * len(self.ops)
        for i, op in enumerate(self.ops):
            hat = SymmOp.from_rotation_and_translation(hex_cell, [0, 0, 0])
            op_tmp = hat * op * hat.inverse
            ops[i] = op_tmp.from_rotation_and_translation(
                op_tmp.rotation_matrix, op.translation_vector)
            # ops[i].translation_vector = op.translation_vector

        return ops

    def get_euclidean_generator(self, cell, idx=0):
        """
        return the symmetry operation object at the Euclidean space

        Args:
            cell: 3*3 cell matrix
            idx: the index of wp generator

        Returns:
            pymatgen SymmOp object
        """
        if not hasattr(self, "generators"):
            self.set_generators()

        op = self.generators[idx]
        if self.euclidean:
            hat = SymmOp.from_rotation_and_translation(cell.T, [0, 0, 0])
            op = hat * op * hat.inverse

        return op

    def get_free_xyzs(self, pos, perturb=False, eps=0.1, random_state: int | None | Generator = None):
        """
        return the free xyz paramters from the given xyz position

        Args:
            pos (array): a 3-array to describe x, y, z
            perturb (bool): whether or not apply perturbation
            eps (float): the magnitude of perturbations

        Returns:
            free xyz array
        """
        if isinstance(random_state, Generator):
            random_state = random_state.spawn(1)[0]
        else:
            random_state = np.random.default_rng(random_state)

        # print(self.apply_ops(pos)[0])
        res = self.apply_ops(pos)[0]
        res = np.delete(res, self.get_frozen_axis())
        if perturb:
            res += eps * random_state.random(len(res)) - 0.5
        res -= np.floor(res)
        return res

    def get_position_from_free_xyzs(self, xyz):
        """
        generate the full xyz position from the free xyzs
        """
        pos = np.zeros(3)
        frozen = self.get_frozen_axis()
        count = 0
        for axis in range(3):
            if axis not in frozen:
                pos[axis] = xyz[count]
                count += 1
        pos = self.apply_ops(pos)[0]
        pos -= np.floor(pos)
        return pos

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

        pos0 = self.search_generator(pos)
        if pos0 is not None:
            res = self.apply_ops(pos0)
            res -= np.floor(res)
            return res
        else:
            return None

    # =============================Evaluations===========================
    def is_standard_setting(self):
        """
        Check if the symmetry operation follows the standard setting
        """
        G_ops = get_wyckoffs(self.number, dim=self.dim)
        for i, ops in enumerate(G_ops):
            if self.has_equivalent_ops(ops):
                self.ops = ops
                self.index = i
                self.letter = letter_from_index(i, G_ops, dim=self.dim)
                return True

        return False

    def has_equivalent_ops(self, wp2, tol=1e-3):
        """
        check if two wps are equivalent

        Args:
            wp2: wp object or list of operations
        """
        ops0 = wp2 if isinstance(wp2, list) else wp2.ops

        if len(ops0) == len(self.ops):
            count = 0
            for _i, op0 in enumerate(ops0):
                for _j, op1 in enumerate(self.ops):
                    diff0 = op0.translation_vector - op1.translation_vector
                    diff0 -= np.rint(diff0)
                    diff1 = op0.rotation_matrix - op1.rotation_matrix
                    if max([np.abs(diff0).sum(), np.abs(diff1).sum()]) < tol:
                        count += 1
            return count == len(ops0)
        else:
            return False

    def is_pure_translation(self, id):
        """
        Check if the operation is equivalent to pure translation
        """
        op = self.generators[id]
        diff = op.rotation_matrix - np.eye(3)
        if np.sum(diff.flatten() ** 2) < 1e-4:
            return True
        else:
            ops = self.get_site_symm_wo_translation()
            return op in ops

    def swap_axis(self, swap_id):
        """
        swap the axis may result in a new wp
        """
        if self.index > 0:
            perm_id = None
            _ops = [self.ops[0]]
            trans = [np.zeros(3)]
            if self.symbol[0] == "F":
                trans.append(np.array([0, 0.5, 0.5]))
                trans.append(np.array([0.5, 0, 0.5]))
                trans.append(np.array([0.5, 0.5, 0]))
            elif self.symbol[0] == "I":
                trans.append(np.array([0.5, 0.5, 0.5]))
            elif self.symbol[0] == "A":
                trans.append(np.array([0, 0.5, 0.5]))
            elif self.symbol[0] == "B":
                trans.append(np.array([0.5, 0, 0.5]))
            elif self.symbol[0] == "C":
                trans.append(np.array([0.5, 0.5, 0]))

            op_perm = swap_xyz_ops(_ops, swap_id)[0]
            for id, ops in enumerate(Group(self.number)):
                if len(ops) == len(self.ops):
                    for i, tran in enumerate(trans):
                        op = op_translation(
                            op_perm, tran) if i > 0 else op_perm
                        # print(id, op.as_xyz_str(),tran)
                        if are_equivalent_ops(op, ops[0]):
                            perm_id = id
                            return Group(self.number)[id], tran
            if perm_id is None:
                raise ValueError("cannot swap", swap_id, self)
        return self, np.zeros(3)

    def print_ops(self, ops=None):
        if ops is None:
            ops = self.ops
        for op in ops:
            print(op.as_xyz_str())

    def gen_pos(self):
        """
        Returns the general Wyckoff position
        """
        return self.ops[0]

    def are_equivalent_pts(self, pt1, pt2, cell=None, tol=0.05):
        """
        Check if two pts are equivalent
        """
        if cell is None:
            cell = np.eye(3)

        pt1 = self.search_generator(pt1, tol=tol)
        pt2 = self.search_generator(pt2, tol=tol)
        if pt1 is None or pt2 is None:
            return False
        else:
            pt1 = np.array(pt1)
            pt1 -= np.floor(pt1)
            pt2 = np.array(pt2)
            pt2 -= np.floor(pt2)
            pts = self.apply_ops(pt1)
            pts -= np.floor(pts)
            diffs = pt2 - pts
            diffs -= np.rint(diffs)
            diffs = np.dot(diffs, cell)
            dists = np.linalg.norm(diffs, axis=1)
            # print(dists)
            return len(dists[dists < tol]) > 0

    def distance_check(self, pt, lattice, tol):
        """
        Given a list of fractional coordinates, merges them within a given
        tolerance, and checks if the merged coordinates satisfy a Wyckoff
        position.

        Args:
            pt: the originl point (3-vector)
            lattice: a 3x3 matrix representing the unit cell
            tol: the cutoff distance for merging coordinates

        Returns:
            True or False
        """
        return not len(self.short_distances(pt, lattice, tol)) > 0

    def short_distances(self, pt, lattice, tol):
        """
        Given a list of fractional coordinates, merges them within a given
        tolerance, and checks if the merged coordinates satisfy a Wyckoff
        position.

        Args:
            pt: the originl point (3-vector)
            lattice: a 3x3 matrix representing the unit cell
            tol: the cutoff distance for merging coordinates

        Returns:
            a list of short distances
        """
        pt = self.project(pt, lattice, self.PBC)
        coor = self.apply_ops(pt)
        # coor -= np.round(coor)
        coor -= np.floor(coor)
        dm = distance_matrix([coor[0]], coor, lattice, PBC=self.PBC)[0][1:]
        # if len(dm[dm<tol]==0): print('+++++', pt, dm.shape, tol, dm[dm<tol], len(dm[dm<tol]))
        return dm[dm < tol]

    def merge(self, pt, lattice, tol, orientations=None, group=None):
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
            wp: a Wyckoff_position object, If no match, returns False.
            valid_ori: the valid orientations after merge

        """
        wp = self.copy()
        PBC = wp.PBC
        if group is None:
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
                    # print(dm); print(y)
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

                        j, k = res
                        if orientations[j][k] == []:
                            continue

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
                    # wp = group[i]
                    p, d = group[i].search_generator_dist(
                        pt.copy(), lattice, group)
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

    def set_generators(self):
        """
        set up generators, useful for many things
        """
        self.generators = get_generators(self.number, dim=self.dim)[self.index]
        if self.P is not None and not identity_ops(self.P):
            # self.print_ops(self.generators)
            ops = transform_ops(self.generators, self.P, self.P1)
            self.generators = ops
            # self.print_ops(ops)

    def set_euclidean(self):
        """
        For the hexagonal groups, need to consider the euclidean conversion
        """
        convert = False
        if self.dim == 3:
            if 143 <= self.number < 195:
                convert = True
        elif self.dim == 2:
            if self.number >= 65:
                convert = True
        elif self.dim == 1 and self.number >= 42:
            convert = True
        self.euclidean = convert

    def search_generator_dist(self, pt, lattice=None, group=None):
        """
        For a given special wp, (e.g., [(x, 0, 1/4), (0, x, 1/4)]),
        return the first position and distance

        Args:
            pt: 1*3 vector
            lattice: 3*3 matrix

        Returns:
            pt: the best matched pt
            diff: numerical difference

        """
        if lattice is None:
            lattice = np.eye(3)

        if self.index == 0:  # general sites
            return pt, 0

        if self.get_dof == 0:  # fixed site like [0, 0, 0]
            pts = self.apply_ops(pt)
            distances = [distance(p0, lattice, PBC=self.PBC) for p0 in pts]
        else:  # sites like (x, 0, 0)
            ops = group[0].ops if group is not None else get_wyckoffs(
                self.number, dim=self.dim)[0]
            pts = []
            distances = []
            for op in ops:
                pt0 = op.operate(pt)
                pt1 = self.ops[0].operate(pt0)
                coord = pt1 - pt0
                distances.append(distance(coord, lattice, PBC=self.PBC))
                pts.append(pt0)

        min_index = np.argmin(distances)
        return pts[min_index], np.min(distances)

    def search_generator(self, pos, ops=None, tol=1e-2, symmetrize=False):
        """
        search generator for a special Wyckoff position

        Args:
            pos: initial xyz position
            ops: list of symops
            tol: tolerance
            symmetrize (bool): apply symmetrization

        Return:
            pos1: the position that matchs the standard setting
        """

        if ops is None:
            ops = get_wyckoffs(self.number, dim=self.dim)[0]

        match = False
        for op in ops:
            pos1 = op.operate(pos)  #
            pos0 = self.ops[0].operate(pos1)
            diff = pos1 - pos0
            diff -= np.rint(diff)
            diff = np.abs(diff)
            # print(self.letter, "{:24s}".format(op.as_xyz_str()), pos, pos0, pos1, diff)
            if diff.sum() < tol:
                pos1 -= np.floor(pos1)
                match = True
                if symmetrize:
                    pos1 = pos0
                break

        if match:
            return pos1
        else:
            # print(pos, wp0, wp)
            return None

    def search_all_generators(self, pos, ops=None, tol=1e-2):
        """
        search generator for a special Wyckoff position

        Args:
            pos: initial xyz position
            ops: list of symops
            tol: tolerance

        Return:
            pos1: the position that matchs the standard setting
        """
        if ops is None:
            ops = get_wyckoffs(self.number, dim=self.dim)[0]

        coords = []
        for op in ops:
            pos1 = op.operate(pos)
            pos0 = self.ops[0].operate(pos1)
            diff = pos1 - pos0
            diff -= np.rint(diff)
            diff = np.abs(diff)
            # print(wp.letter, pos1, pos0, diff)
            if diff.sum() < tol:
                pos1 -= np.floor(pos1)
                coords.append(pos1)
        return coords

    def apply_ops(self, pt):
        """
        apply symmetry operation
        """
        return apply_ops(pt, self.ops)

    def project(self, point, cell=None, PBC=None, id=0):
        """
        Given a 3-vector and a Wyckoff position operator,
        returns the projection onto the axis, plane, or point.

        >>> wp.project_point([0,0.3,0.1],
        array([0. , 0.3, 0.1])

        Args:
            point: a 3-vector (numeric list, tuple, or array)
            cell: 3x3 matrix describing the unit cell vectors
            PBC: A periodic boundary condition list, where 1 means periodic, 0
                means not periodic. Ex: [1,1,1] -> full 3d periodicity, [0,0,1]
                -> 1d periodicity along the z axis

        Returns:
            a transformed 3-vector (numpy array)
        """
        if cell is None:
            cell = np.eye(3)

        # Must be different for hexcell
        if PBC is None:
            PBC = [1, 1, 1]
        op = self.get_euclidean_generator(
            cell, id) if self.euclidean else self.ops[id]

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
                # a faster version?
                b = np.sqrt(basis_vector.dot(basis_vector))
                # if not np.isclose(b, 0):
                if b > 1e-3:
                    new_vector += basis_vector * \
                        (np.dot(point, basis_vector) / (b**2))
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
                new_vector = project_single(pt, rot, trans + v)
                new_vectors.append(new_vector)
                tmp = (new_vector - point).dot(cell)
                distances.append(np.linalg.norm(tmp))
            i = np.argmin(distances)
            return filtered_coords(new_vectors[i], PBC=PBC)

    def to_discrete_grid(self, xyz, N_grids=50):
        """
        A function to convert (x, y, z) to a discrete grid
        """
        binwidth = 1.0 / N_grids
        x = int(xyz[0] // binwidth)
        y = int(xyz[1] // binwidth)
        z = int(xyz[2] // binwidth)
        return [x, y, z]

    def from_discrete_grid(self, xyz, N_grids=50):
        """
        A function to convert from a discrete grid to (x, y, z)
        """
        binwidth = 1.0 / N_grids
        x = binwidth * xyz[0]
        y = binwidth * xyz[1]
        z = binwidth * xyz[2]
        return [x, y, z]

# ----------------- Wyckoff Position selection  --------------------------


def choose_wyckoff(G, number=None, site=None, dim=3, random_state: int | None | Generator = None):
    """
    Choose a Wyckoff position to fill based on the current number of atoms
    needed to be placed within a unit cell
    Rules:
        0) use the pre-assigned list if this is provided
        1) The new position's multiplicity is equal/less than (number).
        2) We prefer positions with large multiplicity.

    Args:
        G: a pyxtal.symmetry.Group object
        number: the number of atoms still needed in the unit cell
        site: the pre-assigned Wyckoff sites (e.g., 4a)

    Returns:
        Wyckoff position. If no position is found, returns False
    """
    if isinstance(random_state, Generator):
        random_state = random_state.spawn(1)[0]
    else:
        random_state = np.random.default_rng(random_state)

    if site is not None:
        number = G.number
        hn = G.hall_number if G.hall_number is not None else None
        return Wyckoff_position.from_group_and_letter(number, site, dim, hn=hn)
    else:
        wyckoffs_organized = G.wyckoffs_organized

        if random_state.random() > 0.5:
            for wyckoff in wyckoffs_organized:
                if len(wyckoff[0]) <= number:
                    # NOTE wyckoff is a ragged list of lists
                    return wyckoff[random_state.choice(len(wyckoff))]
            return False
        else:
            good_wyckoff = [w for wyckoff in wyckoffs_organized if len(
                wyckoff[0]) <= number for w in wyckoff]

            if len(good_wyckoff) > 0:
                # NOTE good_wyckoff is a ragged list of lists
                return good_wyckoff[random_state.choice(len(good_wyckoff))]
            else:
                return False


def choose_wyckoff_mol(
    G: Group,
    number: int,
    site: str | None,
    orientations: list[list[list]],
    gen_site: bool = True,
    dim: int = 3,
    random_state: int | None | Generator = None,
) -> Wyckoff_position | bool:
    """
    Choose a Wyckoff position to fill based on the current number of molecules
    needed to be placed within a unit cell.

    Rules:
    1) The new position's multiplicity is equal/less than (number).
    2) We prefer positions with large multiplicity.
    3) The site must admit valid orientations for the desired molecule.

    Args:
        G: A pyxtal.symmetry.Group object.
        number: The number of molecules still needed in the unit cell.
        site: The specific Wyckoff site to use (if any).
        orientations: The valid orientations for a given molecule.
        gen_site: If True, consider only general Wyckoff positions.
        dim: Dimension of the space group.
        random_state: Seed for random number generation.

    Returns:
        Wyckoff position if found, False otherwise.
    """
    if isinstance(random_state, Generator):
        random_state = random_state.spawn(1)[0]
    else:
        random_state = np.random.default_rng(random_state)

    if site is not None:
        return Wyckoff_position.from_group_and_letter(G.number, site, dim, hn=G.hall_number)

    wyckoffs = G.wyckoffs_organized

    if gen_site or np.random.random() > 0.5:  # choose from high to low
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                good_wyckoffs = []
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoffs.append(w)
                if len(good_wyckoffs) > 0:
                    return good_wyckoffs[random_state.choice(len(good_wyckoffs))]
        return False
    else:
        good_wyckoffs = []
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoffs.append(w)
        if len(good_wyckoffs) > 0:
            return good_wyckoffs[random_state.choice(len(good_wyckoffs))]
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
    if permutation == [0, 1, 2]:
        return xyzs
    else:
        new = []
        for xyz in xyzs:
            tmp = xyz.replace(" ", "").split(",")
            tmp = [tmp[it] for it in permutation]
            if permutation == [1, 0, 2]:  # a,b
                tmp[0] = tmp[0].replace("y", "x")
                tmp[1] = tmp[1].replace("x", "y")
            elif permutation == [2, 1, 0]:  # a,c
                tmp[0] = tmp[0].replace("z", "x")
                tmp[2] = tmp[2].replace("x", "z")
            elif permutation == [0, 2, 1]:  # b,c
                tmp[1] = tmp[1].replace("z", "y")
                tmp[2] = tmp[2].replace("y", "z")
            elif permutation == [1, 2, 0]:  # b,c
                tmp[0] = tmp[0].replace("y", "x")
                tmp[1] = tmp[1].replace("z", "y")
                tmp[2] = tmp[2].replace("x", "z")
            elif permutation == [2, 0, 1]:  # b,c
                tmp[0] = tmp[0].replace("z", "x")
                tmp[1] = tmp[1].replace("x", "y")
                tmp[2] = tmp[2].replace("y", "z")
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
    if permutation == [0, 1, 2]:
        return ops
    else:
        new = []
        for op in ops:
            m = op.affine_matrix.copy()
            m[:3, :] = m[permutation, :]
            for row in range(3):
                # [0, y+1/2, 1/2] -> (0, y, 1/2)
                if np.abs(m[row, :3]).sum() > 0:
                    m[row, 3] = 0
            m[:3, :3] = m[:3, permutation]
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
    m[:3, 3] += tran
    for row in range(3):
        # [0, y+1/2, 1/2] -> (0, y, 1/2)
        if np.abs(m[row, :3]).sum() > 0:
            m[row, 3] = 0
    return SymmOp(m)


def are_equivalent_ops(op1, op2, tol=1e-2):
    """
    check if two ops are equivalent, assuming the same ordering
    """
    diff = op1.affine_matrix - op2.affine_matrix
    diff[:, 3] -= np.rint(diff[:, 3])
    diff = np.abs(diff.flatten())
    return np.sum(diff) < tol


def letter_from_index(index, group, dim=3):
    """
    Given a Wyckoff position's index within a spacegroup, return its number
    and letter e.g. '4a'

    Args:
        index: WP's index (0 is the general position)
        group: an unorganized Wyckoff position array or Group object (preferred)
        dim: the periodicity dimension of the symmetry group.

    Returns:
        the Wyckoff letter corresponding to the Wyckoff position (for example,
        for position 4a, the function would return 'a')
    """
    letters1 = letters
    # See whether the group has an "o" Wyckoff position
    checko = False
    if type(group) == Group and group.dim == 0 or dim == 0:
        checko = True
    if checko is True and len(group[-1]) == 1 and group[-1][0] == SymmOp.from_xyz_str("0,0,0"):
        # o comes before a
        letters1 = "o" + letters
    length = len(group)
    return letters1[length - 1 - index]


def index_from_letter(letter, group, dim=3):
    """
    Given the Wyckoff letter, returns the index of a Wyckoff position.

    Args:
        letter: The wyckoff letter
        group: an unorganized Wyckoff position array or Group object (preferred)
        dim: the periodicity dimension of the symmetry group.

    Returns:
        a single index specifying the location of the Wyckoff position.
    """
    letters1 = letters
    # See whether the group has an "o" Wyckoff position
    checko = False
    if isinstance(group, Group) and group.dim == 0 or dim == 0:
        checko = True
    if checko is True and len(group[-1]) == 1 and group[-1][0] == SymmOp.from_xyz_str("0,0,0"):
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
    for j, a in enumerate(olist):
        for k, _b in enumerate(a):
            num += 1
            if num == i:
                return [j, k]
    return None


class site_symmetry:
    """
    Derive the site symmetry group from symmetry operations
    site-symmetry group is indicated by an oriented symbol,
    which is a variation of the Hermann-Mauguin point-group
    symbol that provides information about the orientation
    of the symmetry elements. The constituents of the oriented
    symbol are ordered according to the symmetry directions of
    the corresponding crystal lattice (primary, secondary and tertiary)

    Args:
        ops: a list of SymmOp objects representing the site symmetry
        lattice_type (str): e.g., 'cubic'
        Bravis (str): 'P', 'R', 'A', 'B', 'C', 'F', 'I'
        parse_trans (bool): do space group or site

    Returns:
        a string representing the site symmetry (e.g., `2mm`)
    """

    def __init__(self, ops, lattice_type, Bravis="P", parse_trans=False):
        # self.G = G
        hexagonal = lattice_type in ["hexagonal", "trigonal"]

        self.parse_trans = parse_trans
        self.opas = [OperationAnalyzer(
            op, parse_trans, hexagonal) for op in ops]
        self.lattice_type = lattice_type
        self.directions = get_symmetry_directions(lattice_type, Bravis)

        if not parse_trans:
            self.symbols = ["1", "-1", "2", "m",
                            "3", "4", "-4", "-3", "6", "-6"]
        else:
            self.symbols = [
                "1",
                "-1",
                "2",
                "2_1",
                "m",
                "a",
                "b",
                "c",
                "n",
                "d",
                "3",
                "3_1",
                "3_2",
                "4",
                "-4",
                "4_1",
                "4_2",
                "4_3",
                "-3",
                "6",
                "6_1",
                "6_2",
                "6_3",
                "6_4",
                "6_5",
                "-6",
            ]
        self.set_table()

        if not parse_trans:
            self.set_full_hm_symbols(self.table)
            self.set_short_symbols()

    def to_one_hot(self):
        matrix = self.to_matrix_representation()
        one_hot_matrix = np.zeros([len(matrix), 13], dtype=int)
        for i, _axis in enumerate(all_sym_directions):
            symbol, id = self.get_highest_symmetry(matrix[i])
            one_hot_matrix[i, id] = 1
        return one_hot_matrix

    def to_matrix_representation_spg(self):
        """
        To create a binary matrix to represent the symmetry elements on each axis
        Translation is alos counted here.
        """

        matrix = np.zeros(
            [len(all_sym_directions), len(self.symbols)], dtype=int)
        # every direction must has identity symmetry
        matrix[:, 0] = 1
        self.inversion = False

        for opa in self.opas:
            if opa.type == "inversion":
                self.inversion = True

            elif opa.type != "identity":
                _ax0 = opa.axis / np.linalg.norm(opa.axis)

                store = False
                for i, ax in enumerate(all_sym_directions):
                    # ; print(opa.axis, ax, np.dot(_ax0, ax0))
                    ax0 = ax / np.linalg.norm(ax)
                    if np.isclose(abs(np.dot(_ax0, ax0)), 1):
                        store = True
                        break

                if store:
                    # Pure rotation
                    # print('add symmetry', i, ax, opa.type, opa.order)
                    if opa.symbol in self.symbols:
                        matrix[i, self.symbols.index(opa.symbol)] = 1
                    else:
                        print("To debug", opa.symbol, opa)
                        import sys

                        sys.exit()
                else:
                    raise ValueError("Cannot parse the axis",
                                     opa.axis, all_sym_directions)

            if self.inversion:
                matrix[:, 1] = 1  # if inversion is present
        return matrix

    def to_matrix_representation(self):
        """
        To create a 15 * 10 binary matrix to represent the
        symmetry elements on each axis
        #[1, -1, 2, m, 3, 2/m, 4, -3, 6, 4/m, -6, 6/m]
        [1, -1, 2, m, 3, 4, -3, 6, -6]
        """
        symbols = ["1", "-1", "2", "m", "3", "4", "-4", "-3", "6", "-6"]
        matrix = np.zeros([len(all_sym_directions), len(symbols)], dtype=int)
        # every direction must has identity symmetry
        matrix[:, 0] = 1
        self.inversion = False

        for opa in self.opas:
            if opa.type == "inversion":
                self.inversion = True

            elif opa.type != "identity":
                _ax0 = opa.axis / np.linalg.norm(opa.axis)

                for i, ax in enumerate(all_sym_directions):
                    store = False
                    # ; print(opa.axis, ax, np.dot(opa.axis, ax0))
                    ax0 = ax / np.linalg.norm(ax)
                    if np.isclose(abs(np.dot(_ax0, ax0)), 1):
                        store = True
                        break

                if store:
                    # Pure rotation
                    # print('add symmetry', i, ax, opa.type, opa.order)
                    if opa.type == "rotation":
                        if opa.order == 2:
                            matrix[i, 2] = 1
                        elif opa.order == 3:
                            matrix[i, 4] = 1
                        elif opa.order == 4:
                            matrix[i, 5] = 1
                        elif opa.order == 6:
                            matrix[i, 8] = 1
                        else:
                            raise RuntimeError(
                                "Unexpected rotation order", opa.order)
                    elif opa.type == "rotoinversion":
                        if opa.rotation_order == 2:
                            matrix[i, 3] = 1  # -2 is m
                        elif opa.rotation_order == 3:
                            matrix[i, 7] = 1  # -3
                        elif opa.rotation_order == 4:
                            matrix[i, 6] = 1  # -4
                        elif opa.rotation_order == 6:
                            matrix[i, 9] = 1
                        else:
                            raise RuntimeError(
                                "Unexpected rotinversion order", opa.order)
                else:
                    raise ValueError("Cannot parse the axis",
                                     opa.axis, all_sym_directions)

            if self.inversion:
                matrix[:, 1] = 1  # if inversion is present
        return matrix

    def set_table(self, skip=False):
        """
        Get the complete table representation.

        Args:
            skip (bool): whether or not skip 1 or -1 symmetry

        Returns:
            sorted table with (list of symmetry elements, symbols, order)
        """
        # Complete list of symmetry for one given axis
        # symbols = ['1', '-1', '2', 'm', '3', '4', '-4', '-3', '6', '-6']
        if self.lattice_type == "triclinic":
            skip = False

        matrix = self.to_matrix_representation_spg(
        ) if self.parse_trans else self.to_matrix_representation()

        tables = []
        for i, axis in enumerate(all_sym_directions):
            direction_id = find_axis_order(axis, self.directions)
            if direction_id is not None:
                num_symmetries = matrix[i, 2:].sum(
                ) if skip else matrix[i].sum()
                if num_symmetries > 0:
                    strs = "{:4d} ({:2d} {:2d} {:2d}): ".format(
                        direction_id, *axis)
                    for sym in matrix[i]:
                        strs += f"{sym:4d} "
                    # strs += "{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}".format(*matrix[i])
                    if not self.parse_trans:
                        symbol, _ = self.get_highest_symmetry(matrix[i])
                        strs += f"{symbol:>6s}"
                        tables.append((strs, symbol, direction_id))
                    else:
                        tables.append((strs, direction_id))
            # else:
            #    raise ValueError('Wrong input axis', axis, 'lattice_type', self.lattice_type)
        self.table = sorted(tables, key=lambda x: x[-1])

    def set_full_hm_symbols(self, tables):
        """
        Set the full hm symbols for each axis

        Args:
            tables: sorted table with (list of symmetry elements, symbols, order)

        Returns:
            a list of symmetry elements on {primary, secondary, tertiery} directions
        """
        hm_symbols = [[] for _ in range(len(self.directions))]
        # for row in tables: print(row)
        for row in tables:
            (_, symbol, direction_id) = row
            if symbol not in ["1", "-1"]:
                hm_symbols[direction_id].append(symbol)
            # print(hm_symbols, direction_id)

        for i, hm_symbol in enumerate(hm_symbols):
            if len(hm_symbol) == 0:
                hm_symbols[i] = ["."]
            # elif hm_symbol == ['1']:
            #    hm_symbols[i] = ['.']
        self.hm_symbols = hm_symbols

    def unique_symmetry(self, symbols, symmetry):
        return all(symbol in [".", symmetry] for symbol in symbols)

    def ref_symmetry(self, symbols, reference):
        return any(symbol in reference for symbol in symbols)

    def set_short_symbols(self):
        """
        Set short symbols from the Full symbols
        """
        # if hasattr(self, 'hm_symbols'):
        #    self.set_full_hm_symbols()

        self.symbols = []
        # print(self.hm_symbols)
        for hm_symbol in self.hm_symbols:
            if len(hm_symbol) == 1:
                # print('single', hm_symbol)
                self.symbols.append(hm_symbol[0])
            else:
                symbol = ""
                for hm in hm_symbol:
                    symbol += hm
                self.symbols.append(symbol)
                # print('multi', hm_symbol)

        # Some simplifications
        if self.lattice_type == "orthorhombic":
            if self.symbols == ["2/m", "2/m", "2/m"]:
                self.symbols = ["m", "m", "m"]

        elif self.lattice_type == "tetragonal":
            for i, symbol in enumerate(self.symbols):
                if symbol == "2/m2/m":
                    self.symbols[i] = "mm"
                elif symbol == "2/m":
                    if not self.unique_symmetry(self.symbols, "2/m"):
                        self.symbols[i] = "m"
                elif symbol == "m2":
                    self.symbols[i] = "2m"
            if self.symbols == ["4", "22", "22"]:
                self.symbols = ["4", "2", "2"]
            elif self.symbols == ["4", "mm", "mm"]:
                self.symbols = ["4", "m", "m"]
            elif self.symbols == ["-4", "22", "mm"]:
                self.symbols = ["-4", "2", "m"]
            elif self.symbols == ["-4", "mm", "22"]:
                self.symbols = ["-4", "m", "2"]
            elif self.symbols == ["2/m", "2/m", "2/m"]:
                self.symbols = ["m", "m", "m"]
            elif self.symbols == ["4/m", "mm", "mm"]:
                self.symbols = ["4/m", "m", "m"]

        elif self.lattice_type in ["trigonal", "hexagonal"]:
            for i, symbol in enumerate(self.symbols):
                if symbol in ["2/m2/m", "2/m2/m2/m", "mm", "mmm"]:
                    if not self.unique_symmetry(self.symbols, symbol):
                        self.symbols[i] = "m"
                elif symbol in ["22", "222"] and not self.unique_symmetry(self.symbols, symbol):
                    self.symbols[i] = "2"
            if self.symbols == ["2/m", "2/m", "2/m"]:
                self.symbols = ["m", "m", "m"]

        elif self.lattice_type == "cubic":
            for i, symbol in enumerate(self.symbols):
                if symbol in ["2/m2/m2/m"]:
                    if not self.unique_symmetry(self.symbols, symbol):
                        self.symbols[i] = "m"
                    else:
                        self.symbols[i] = "mmm"
                elif symbol == "mmm":
                    if not self.unique_symmetry(self.symbols, symbol):
                        self.symbols[i] = "m"
                elif symbol == "222":
                    # print(symbol, self.unique_symmetry(self.symbols, symbol))
                    if not self.unique_symmetry(self.symbols, symbol):
                        self.symbols[i] = "2"
                elif symbol in ["2222", "222222"]:
                    self.symbols[i] = "2"
                elif symbol in ["333", "3333"]:
                    self.symbols[i] = "3"
                elif symbol in ["-3-3-3-3"]:
                    self.symbols[i] = "-3"
                elif symbol == "444":
                    self.symbols[i] = "4"
                elif symbol == "-4-4-4":
                    self.symbols[i] = "-4"
                elif symbol == "422":
                    self.symbols[i] = "42"
                elif symbol == "-422":
                    self.symbols[i] = "-42"
                elif symbol == "2mm":
                    self.symbols[i] = "mm2"
                elif symbol in ["mmmm", "mmmmmm"]:
                    self.symbols[i] = "m"
                elif symbol in ["2m"]:
                    self.symbols[i] = "m2"
                elif symbol in ["4mm"]:
                    self.symbols[i] = "4m"
                elif symbol in ["-4mm"]:
                    self.symbols[i] = "-4m"
                elif symbol == "4/m2/m2/m":
                    self.symbols[i] = "4/mm"
                elif symbol in ["4/m4/m4/m", "2/m2/m2/m2/m2/m2/m"]:
                    self.symbols[i] = "m"
                elif symbol == "2/m2/m":
                    if self.ref_symmetry(self.symbols, ["4/mm"]):
                        self.symbols[i] = "m"
                    else:
                        self.symbols[i] = "mm"
                elif symbol == "2/m" and not self.unique_symmetry(self.symbols, symbol):
                    self.symbols[i] = "m"

            for i, symbol in enumerate(self.symbols):
                if symbol == "mm" and self.ref_symmetry(self.symbols, ["-42", "4m"]):
                    self.symbols[i] = "m"
                if symbol == "22" and self.ref_symmetry(self.symbols, ["42", "-4m"]):
                    self.symbols[i] = "2"

        #    #if self.symbols in [['4', '-3', '2'], ['-4', '-3', 'm']]:
        #    #    self.symbols = ['m', '-3', 'm']
        #    #if '222' in self.symbols:
        #    #    if len(self.opas) > 4:
        #    #        for i in range(len(self.symbols)):
        #    #            if self.symbols[i] == '222':
        #    #                self.symbols[i] = '2'#; print('Find ===')

        self.get_name()

    def get_name(self):
        if self.symbols in [[".", ".", "."], [".", "."], ["."]]:
            if self.inversion:
                self.name = "-1"
            else:
                self.name = "1"
        else:
            self.name = ""
            for symbol in self.symbols:
                # self.name += ' '
                for s in symbol:
                    self.name += s

    def to_beautiful_matrix_representation(self, skip=True):
        """
        A shortcut to check the representation

        Args:
            skip (bool): whether or not skip 1 or -1 symmetry
        """
        strs = "Order    Axis       "
        if self.parse_trans:
            symbols = [
                "1",
                "-1",
                "2",
                "2_1",
                "m",
                "a",
                "b",
                "c",
                "n",
                "d",
                "3",
                "3_1",
                "3_2",
                "4",
                "-4",
                "4_1",
                "4_2",
                "4_3",
                "-3",
                "6",
                "6_1",
                "6_2",
                "6_3",
                "6_4",
                "6_5",
                "-6",
            ]
        else:
            symbols = ["1", "-1", "2", "m", "3", "4", "-4", "-3", "6", "-6"]

        for symbol in symbols:
            strs += f"{symbol:<4s} "
        print(strs)
        if not hasattr(self, "table"):
            self.set_table(skip)
        for row in self.table:
            print(row[0])

    def get_highest_symmetry(self, row):
        # ['1']
        # ['1', '-1']
        # ['1', '2']
        # ['1', 'm']
        # ['1', '3']
        # ['1', '2', 'm', '2/m']
        # ['1', '2', '4']
        # ['1', '2', '-4']
        # ['1', '-1', '3', '-3']
        # ['1', '2', '3', '6']
        # ['1', 'm', '3', '-6']
        # ['1', '-1', '2', 'm', '4', '-4', '4/m']
        # ['1', '-1', '2', 'm', '3', '-3', '6', '-6', '6/m']
        ref_arrays = [
            (np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=int), "1"),
            (np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype=int), "-1"),
            (np.array([1, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype=int), "2"),
            (np.array([1, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype=int), "m"),
            (np.array([1, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype=int), "3"),
            (np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0], dtype=int), "2/m"),
            (np.array([1, 0, 1, 0, 0, 1, 0, 0, 0, 0], dtype=int), "4"),
            (np.array([1, 0, 1, 0, 0, 0, 1, 0, 0, 0], dtype=int), "-4"),
            (np.array([1, 1, 0, 0, 1, 0, 0, 1, 0, 0], dtype=int), "-3"),
            (np.array([1, 0, 1, 0, 1, 0, 0, 0, 1, 0], dtype=int), "6"),
            (np.array([1, 0, 0, 1, 1, 0, 0, 0, 0, 1], dtype=int), "-6"),
            (np.array([1, 1, 1, 1, 0, 1, 1, 0, 0, 0], dtype=int), "4/m"),
            (np.array([1, 1, 1, 1, 1, 0, 0, 1, 1, 1], dtype=int), "6/m"),
        ]

        for i, ref_array in enumerate(ref_arrays):
            if np.array_equal(row, ref_array[0]):
                return ref_array[1], i

        if self.lattice_type not in ["hexagonal", "trigonal"]:
            symbols = ["1", "-1", "2", "m", "3", "4", "-4", "-3", "6", "-6"]
            strs = [symbols[i] for i, x in enumerate(row) if x == 1]
            print(row)
            #raise ValueError("Incompatible symmetry list", strs)
        return ref_arrays[0][1], 0


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
    wyckoffs = group.Wyckoff_positions if type(group) == Group else group
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
        return None
    # Vector must be non-zero
    if axis.dot(axis) < 1e-6:
        return None
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


def get_wyckoffs(num, organized=False, dim=3):
    """
    Returns a list of Wyckoff positions for a given group. Has option to
    organize the list based on multiplicity (this is used for
    random_crystal.wyckoffs)

    For an unorganized list:

        - 1st index: index of WP in sg (0 is the WP with largest multiplicity)
        - 2nd index: a SymmOp object in the WP

    For an organized list:

        - 1st index: specifies multiplicity (0 is the largest multiplicity)
        - 2nd index: a WP within the group of equal multiplicity.
        - 3nd index: a SymmOp object within the Wyckoff position

    You may switch between organized and unorganized lists using the methods
    i_from_jk and jk_from_i. For example, if a Wyckoff position is the [i]
    entry in an unorganized list, it will be the [j][k] entry in an organized
    list.

    Args:
        num: the international group number
        dim: dimension [0, 1, 2, 3]
        organized: whether or not to organize the list based on multiplicity

    Returns:
        a list of Wyckoff positions, each of which is a list of SymmOp's
    """
    if dim == 3:
        df = SYMDATA.get_wyckoff_sg()
    elif dim == 2:
        df = SYMDATA.get_wyckoff_lg()
    elif dim == 1:
        df = SYMDATA.get_wyckoff_rg()
    elif dim == 0:
        df = SYMDATA.get_wyckoff_pg()

    # Convert the string from df into a list of wyckoff strings
    wyckoff_strings = literal_eval(df["0"][num])  # Use literal_eval instead of eval

    wyckoffs = []
    for x in wyckoff_strings:
        wyckoffs.append([])
        for y in x:
            if dim == 0:
                wyckoffs[-1].append(SymmOp(y))
            else:
                wyckoffs[-1].append(SymmOp.from_xyz_str(y))
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
        # Return Wyckoff positions without organization
        return wyckoffs

def get_wyckoff_symmetry(num, dim=3):
    """
    Returns a list of site symmetry for a given group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a point within the WP
    3rd index: a site symmetry SymmOp of the point

    Args:
        sg: the international spacegroup number
        dim: 0, 1, 2, 3

    Returns:
        a 3d list of SymmOp objects representing the site symmetry of each
        point in each Wyckoff position
    """
    if dim == 3:
        symmetry_df = SYMDATA.get_symmetry_sg()
    elif dim == 2:
        symmetry_df = SYMDATA.get_symmetry_lg()
    elif dim == 1:
        symmetry_df = SYMDATA.get_symmetry_rg()
    elif dim == 0:
        symmetry_df = SYMDATA.get_symmetry_pg()

    symmetry_strings = eval(symmetry_df["0"][num])

    symmetry = []
    # Loop over Wyckoff positions
    for x in symmetry_strings:
        symmetry.append([])
        # Loop over points in WP
        for y in x:
            symmetry[-1].append([])
            # Loop over ops
            for z in y:
                op = SymmOp(z) if dim == 0 else SymmOp.from_xyz_str(z)
                symmetry[-1][-1].append(op)
    return symmetry


def get_generators(num, dim=3):
    """
    Returns a list of Wyckoff generators for a given group.
    1st index: index of WP in sg (0 is the WP with largest multiplicity)
    2nd index: a generator for the WP
    This function is useful for rotating molecules based on Wyckoff position,
    since special Wyckoff positions only encode positional information, but not
    information about the orientation. The generators for each Wyckoff position
    form a subset of the spacegroup's general Wyckoff position.

    Args:
        num: the international spacegroup number
        dim: dimension

    Returns:
        a 2d list of symmop objects [[wp0], [wp1], ... ]
    """

    if dim == 3:
        generators_df = SYMDATA.get_generator_sg()
    elif dim == 2:
        generators_df = SYMDATA.get_generator_lg()
    elif dim == 1:
        generators_df = SYMDATA.get_generator_rg()
    elif dim == 0:
        generators_df = SYMDATA.get_generator_pg()

    generator_strings = eval(generators_df["0"][num])

    generators = []
    # Loop over Wyckoff positions
    for x in generator_strings:
        generators.append([])
        # Loop over ops
        for y in x:
            op = SymmOp.from_xyz_str(y) if dim > 0 else SymmOp(y)
            generators[-1].append(op)
    return generators


def site_symm(point, gen_pos, tol=1e-3, lattice=None, PBC=None):
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
    if lattice is None:
        lattice = np.eye(3)

    if PBC is None:
        PBC = gen_pos.PBC if type(gen_pos) == Wyckoff_position else [1, 1, 1]
    # Convert point into a SymmOp
    if type(point) != SymmOp:
        point = SymmOp.from_rotation_and_translation(
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]], np.array(point))
    symmetry = []
    for op in gen_pos:
        is_symmetry = True
        # Calculate the effect of applying op to point
        difference = SymmOp((op * point).affine_matrix - point.affine_matrix)
        # Check that the rotation matrix is unaltered by op
        if not np.allclose(difference.rotation_matrix, np.zeros((3, 3)), rtol=1e-3, atol=1e-3):
            is_symmetry = False
        # Check that the displacement is less than tol
        displacement = difference.translation_vector
        if distance(displacement, lattice, PBC=PBC) > tol:
            is_symmetry = False
        if is_symmetry:
            """
            The actual site symmetry's translation vector may vary from op by
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
                op.rotation_matrix, op.translation_vector - np.rint(displacement))
            symmetry.append(el)
    return symmetry


def check_wyckoff_position(points, group, tol=1e-3):
    """
    Given a list of points, returns a single index of a matching Wyckoff
    position in the space group. Checks the site symmetry of each supplied
    point against the site symmetry for each point in the Wyckoff position.
    Also returns a point which can be used to generate the rest using the
    Wyckoff position operators.

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
    t = tol**2
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
            dw = distance_matrix(points, pw, None, PBC=PBC,
                                 metric="sqeuclidean")

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
    Function for quick conversion between symbols and numbers.

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

    group_symbols = loadfn(rf("pyxtal", "database/symbols.json"))
    lists = group_symbols[keys[dim]]
    number = None
    symbol = None
    if dim not in [0, 1, 2, 3]:
        raise ValueError(f"Dimension ({dim:d}) should in [0, 1, 2, 3] ")

    if isinstance(input_group, str):
        for i, _symbol in enumerate(lists):
            if _symbol == input_group:
                number = i + 1
                symbol = input_group
                return symbol, number
        msg = f"({input_group:s}) not found in {keys[dim]:s} "
        raise ValueError(msg)

    valid, msg = check_symmetry_and_dim(input_group, dim)
    if not valid:
        raise ValueError(msg)

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
    msg = "This is a valid group number"

    numbers = [56, 75, 80, 230]
    if dim not in [0, 1, 2, 3]:
        msg = f"invalid dimension {dim:d}"
        valid = False
    else:
        max_num = numbers[dim]
        if number not in range(1, max_num + 1):
            valid = False
            msg = f"invalid symmetry group {number:d}"
            msg += f" in dimension {dim:d}"
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
        elif number <= 167:
            lattice_type = "trigonal"
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
        lattice_type = "spherical" if number in [
            1, 2, 6, 8, 28, 29, 30, 31, 32, 55, 56] else "ellipsoidal"
    return PBC, lattice_type


def search_cloest_wp(G, wp, op, pos):
    """
    For a given position, search for the cloest wp which satisfies the desired
    symmetry relation, e.g., for pos (0.1, 0.12, 0.2) and op (x, x, z) the
    closest match is (0.11, 0.11, 0.2)

    Args:
        G: space group number
        wp: Wyckoff object
        op: symmetry operation belonging to wp
        pos: initial xyz position

    Return:
        pos1: the position that matchs symmetry operation
    """
    # G = Group(wp.number)
    if np.linalg.matrix_rank(op.rotation_matrix) == 0:
        # fixed point (e.g, 1/2, 1/2, 1/2)
        return op.translation_vector
    elif np.linalg.matrix_rank(op.rotation_matrix) == 3:
        # fully independent, e.g., (x,y,z), (-x,y,z)
        return pos
    else:
        # check if this is already matched
        wp0 = G[0]
        coords = wp.search_all_generators(pos, wp0)
        if len(coords) > 0:
            diffs = []
            for coord in coords:
                tmp = op.operate(coord)
                diff1 = tmp - pos
                diff1 -= np.rint(diff1)
                dist = np.linalg.norm(diff1)
                if dist < 1e-3:
                    return tmp
                else:
                    diffs.append(dist)
            minID = np.argmin(diffs)
            return op.operate(coords[minID])

        # if not match, search for the closest solution
        else:
            # extract all possible xyzs
            all_xyz = apply_ops(pos, wp0)[1:]
            dists = all_xyz - pos
            dists -= np.rint(dists)
            ds = np.linalg.norm(dists, axis=1)
            ids = np.argsort(ds)
            for id in ids:
                d = all_xyz[id] - pos
                d -= np.rint(d)
                res = pos + d / 2
                if wp.search_generator(res, wp0) is not None:
                    # print(ds[id], pos, res)
                    return res
            return op.operate(pos)


def get_point_group(number):
    """
    Parse the point group symmetry info from space group. According to
    http://img.chem.ucl.ac.uk/sgp/misc/pointgrp.htm, among 32(230) point(space)
    groups, there are
        - 10(68) polar groups,
        - 11(92) centrosymmetric groups,
        - 11(65) enantiomorphic groups

    Args:
        number: space group number

    Return:
        point group symbol
        polar: 1, 2, m, mm2, 3, 3m, 4, 4mm, 6, 6mm
        centrosymmetry: -1, 2/m, mmm, 4/m, 4/mmm, -3, -3m, 6/m, 6/mmm, m-3, m-3m
        enantiomorphic: 1, 2, 222, 4, 422, 3, 32, 6, 622, 23, 432
    """

    if number == 1:
        return "1", 1, True, False, True
    elif number == 2:
        return "-1", 2, False, True, False
    elif 3 <= number <= 5:
        return "2", 3, True, False, True
    elif 6 <= number <= 9:
        return "m", 4, True, False, False
    elif 10 <= number <= 15:
        return "2/m", 5, False, True, False
    elif 16 <= number <= 24:
        return "222", 6, False, False, True
    elif 25 <= number <= 46:
        return "mm2", 7, True, False, False
    elif 47 <= number <= 74:
        return "mmm", 8, False, True, False
    elif 75 <= number <= 80:
        return "4", 9, True, False, True
    elif 81 <= number <= 82:
        return "-4", 10, False, False, False
    elif 83 <= number <= 88:
        return "4/m", 11, False, True, False
    elif 89 <= number <= 98:
        return "422", 12, False, False, True
    elif 99 <= number <= 110:
        return "4mm", 13, True, False, False
    elif 111 <= number <= 122:
        return "-42m", 14, False, False, False
    elif 123 <= number <= 142:
        return "4/mmm", 15, False, True, False
    elif 143 <= number <= 146:
        return "3", 16, True, False, True
    elif 147 <= number <= 148:
        return "-3", 17, False, True, False
    elif 149 <= number <= 155:
        return "32", 18, False, False, True
    elif 156 <= number <= 161:
        return "3m", 19, True, False, False
    elif 162 <= number <= 167:
        return "-3m", 20, False, True, False
    elif 168 <= number <= 173:
        return "6", 21, True, False, True
    elif number == 174:
        return "-6", 22, False, False, False
    elif 175 <= number <= 176:
        return "6/m", 23, False, True, False
    elif 177 <= number <= 182:
        return "622", 24, False, False, True
    elif 183 <= number <= 186:
        return "6mm", 25, True, False, False
    elif 187 <= number <= 190:
        return "-62m", 26, False, False, False
    elif 191 <= number <= 194:
        return "6/mmm", 27, False, True, False
    elif 195 <= number <= 199:
        return "23", 28, False, False, True
    elif 200 <= number <= 206:
        return "m-3", 29, False, True, False
    elif 207 <= number <= 214:
        return "432", 30, False, False, True
    elif 215 <= number <= 220:
        return "-43m", 31, False, False, False
    elif 221 <= number <= 230:
        return "m-3m", 32, False, True, False
    return None


def get_close_packed_groups(pg):
    """
    List the close packed groups based on the molecular symmetry.
    Compiled from AIK Book, Table 2 P34.

    Args:
        pg (str): Point group symbol.

    Returns:
        list or None: List of space group numbers, or None if not found.
    """
    close_packed_groups = {
        "1": [1, 2, 4, 14, 19, 29, 33, 51, 54, 61, 62],
        "2": [1, 15, 18, 60],
        "m": [1, 26, 36, 63, 64],
        "I": [1, 2, 14, 15, 61],
        "mm": [42, 51, 59],
        "2/m": [12, 54, 64],
        "222": [21, 22, 23, 68],
        "mmm": [65, 69, 71],
    }

    return close_packed_groups.get(pg)


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
    # Triclinic: 1
    if pg == "-1":  # 2
        return ["1"]
    # Monoclinic: 5
    elif pg in ["2", "m"]:  # 2
        return "1"
    elif pg == "2/m":  # 3
        return ["1", "m", "2"]
    # Orthorhombic: #7
    elif pg == "222":  # 2
        return ["1", "2"]
    elif pg == "mm2":  # 2
        return ["1", "m"]
    elif pg == "mmm":  # 3
        return ["1", "m", "mm2"]
    # Tetragonal: 20
    elif pg == "4":  # 1
        return ["1"]
    elif pg == "-4":  # 2
        return ["1", "2"]
    elif pg == "4/m":  # 3
        return ["1", "2", "4"]
    elif pg == "422":  # 3
        # return ['1', '2(s)', '4']
        return ["1", "2", "4"]
    elif pg == "4mm":  # 2
        return ["1", "m"]
    elif pg == "-42m":  # 4
        # return ['1', '2(s)', 'm', 'mm2']
        return ["1", "2", "m", "mm2"]
    elif pg == "4/mmm":  # 5
        # return ['1', 'm(s)', 'm(p)', 'mm2(s)', '4mm']
        return ["1", "m", "mm2", "4mm"]
    # Trigonal: 12
    elif pg == "3":  # 1
        return ["1"]
    elif pg == "-3":  # 2
        return ["1", "3"]
    elif pg == "32":  # 3
        return ["1", "2", "3"]
    elif pg == "3m":  # 2
        return ["1", "m"]
    elif pg == "-3m":  # 4
        return ["1", "2", "m", "3m"]
    # Hexagonal: 22
    elif pg == "6":  # 1
        return ["1"]
    elif pg == "-6":  # 3
        return ["1", "m", "3"]
    elif pg == "6/m":  # 3
        return ["1", "m", "6"]
    elif pg == "622":  # 3
        # return ['1', '2(s)', '6']
        return ["1", "2", "6"]
    elif pg == "6mm":  # 2
        return ["1", "2"]
    elif pg in ["-62m", "-6m2"]:  # 5
        # return ['1', 'm(s)', 'm(p)', 'mm2', '3m']
        return ["1", "m", "mm2", "3m"]
    elif pg == "6/mmm":  # 5
        # return ['1', 'm(s)', 'm(p)', 'mm2(s)', '6mm']
        return ["1", "m", "mm2", "6mm"]
    # Cubic: 21
    elif pg == "23":  # 3
        return ["1", "2", "3"]
    elif pg == "m-3":  # 4
        return ["1", "m", "mm2", "3"]
    elif pg == "432":  # 4
        # return ['1', '2(s)', '4', '3']
        return ["1", "2", "4", "3"]
    elif pg == "-43m":  # 4
        return ["1", "m", "mm2", "3m"]
    elif pg == "m-3m":  # 6
        # return ['1', 'm(s)', 'm(p)', 'mm2', '4mm', '3m']
        return ["1", "m", "mm2", "4mm", "3m"]
    return None


def get_all_polar_space_groups():
    ps, nps = [], []
    for i in range(1, 231):
        g = Group(i, quick=True)
        if g.polar:
            ps.append(i)
        else:
            nps.append(i)
    return ps, nps


def abc2matrix(abc):
    """
    convert the abc string representation to matrix
    Args:
        abc: string like 'a, b, c' or 'a+c, b, c' or 'a+1/4, b+1/4, c'

    Returns:
        4*4 affine matrix
    """
    rot_matrix = np.zeros((3, 3))
    trans = np.zeros(3)
    toks = abc.strip().replace(" ", "").lower().split(",")
    re_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([a-c])")
    re_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![a-c])")
    for i, tok in enumerate(toks):
        # build the rotation matrix
        for m in re_rot.finditer(tok):
            factor = -1.0 if m.group(1) == "-" else 1.0
            if m.group(2) != "":
                if m.group(3) != "":
                    factor *= float(m.group(2)) / float(m.group(3))
                else:
                    factor *= float(m.group(2))
            j = ord(m.group(4)) - 97
            try:
                rot_matrix[i, j] = factor
            except:
                print(abc)
                import sys

                sys.exit()

        # build the translation vector
        for m in re_trans.finditer(tok):
            factor = -1 if m.group(1) == "-" else 1
            num = float(m.group(2)) / float(m.group(3)
                                            ) if m.group(3) != "" else float(m.group(2))
            trans[i] = num * factor
    return (rot_matrix, trans)


def get_symmetry_from_ops(ops, tol=1e-5):
    """
    get the hall number from symmetry operations.

    Args:
        ops: tuple of (rotation, translation) or list of strings
        tol: tolerance
    """
    from spglib import get_hall_number_from_symmetry

    if isinstance(ops[0], str):
        ops = [SymmOp.from_xyz_str(op) for op in ops]
    rot = [op.rotation_matrix for op in ops]
    tran = [op.translation_vector for op in ops]
    hall_number = get_hall_number_from_symmetry(rot, tran, tol)
    spg_number = HALL_TABLE["Spg_num"][hall_number - 1]
    return hall_number, spg_number


def identity_ops(op):
    """
    check if the operation is the identity.
    """
    (rot, trans) = op
    return bool(np.allclose(rot, np.eye(3)) and np.sum(np.abs(trans)) < 0.001)


def transform_ops(ops, P, P1):
    """
    Transformation according to the P and P1 operations.

    Args:
        ops: list of symmtry ops
        P: transformation
        P1: inverse transformation
    """
    # print("}++++++++++++++++++++++", len(ops))
    rot_P = P[0].T
    rot_Q = P1[0].T
    tran_P = P[1]
    for i, op1 in enumerate(ops):
        # R = Q * R * P, suitable when P = {a-c, b, c}
        tran = rot_Q.dot(op1.translation_vector) - tran_P
        rot = rot_Q.dot(op1.rotation_matrix).dot(rot_P)
        op2 = SymmOp.from_rotation_and_translation(rot, tran)
        ops[i] = op2
        # print("{:25s} ==> {:25s} ==> {:24s}".format(s1, s2, s3))

    # in case of (x+1/2, y, z) as the first
    if np.linalg.norm(tran_P) > 1e-3:
        base = ops[0].translation_vector
        for i, op in enumerate(ops):
            inv = np.linalg.inv(op.affine_matrix)
            trans = inv[:3, 3] + base
            trans -= np.rint(trans)
            rot_ops = ops[i].rotation_matrix
            ops[i] = SymmOp.from_rotation_and_translation(rot_ops, trans)

    return ops


def trim_ops(ops):
    """
    Convert the operation to the simplest form. e.g.,
        - 'x+1/8, y+1/8, z+1/8' -> 'x, y, z'
        - '1/8 y+1/8 -y+1/8' -> '1/8, y, -y+1/4'

    Args:
        ops: A list of symmetry operations.

    Returns:
        A list of simplified symmetry operations.
    """

    def in_base(op, base):
        for b in base:
            if np.allclose(op, b[:3]) or np.allclose(op, -b[:3]):
                return b
        return None

    def process_rotation(rot, tran, base):
        for i in range(3):
            tmp = rot[i, :]
            if np.linalg.norm(tmp) > 1e-3:
                b = in_base(tmp, base)
                if b is None:
                    _base = np.zeros(4)
                    _base[:3] = tmp
                    _base[3] = tran[i]
                    base.append(_base)
                    tran[i] = 0
                else:
                    coef = next(tmp[j] / b[j]
                                for j in range(3) if abs(b[j]) > 0)
                    tran[i] -= coef * b[3]
        return tran

    base = []
    simplified_ops = []

    for op in ops:
        rot = op.rotation_matrix
        tran = op.translation_vector
        tran = process_rotation(rot, tran, base)

        simplified_ops.append(SymmOp.from_rotation_and_translation(rot, tran))

    return simplified_ops


def find_axis_order(axis, directions):
    for i, axes in enumerate(directions):
        for ax in axes:
            if ax == axis:
                return i
    return None


def get_symmetry_directions(lattice_type, symbol="P", unique_axis="b"):
    """
    Get the symmetry directions
    """
    if lattice_type == "monoclinic":
        if unique_axis == "b":
            return [[(0, 1, 0)]]
        elif unique_axis == "c":
            return [[(0, 0, 1)]]
        else:
            return [[(1, 0, 0)]]
    elif lattice_type == "orthorhombic":
        return [[(1, 0, 0)], [(0, 1, 0)], [(0, 0, 1)]]
    elif lattice_type == "tetragonal":
        return [[(0, 0, 1)], [(1, 0, 0), (0, 1, 0)], [(1, -1, 0), (1, 1, 0)]]
    elif lattice_type == "hexagonal" or (lattice_type == "trigonal" and symbol == "P"):
        return [
            [(0, 0, 1)],
            [(1, 0, 0), (0, 1, 0), (1, 1, 0)],
            [(1, -1, 0), (1, 2, 0), (2, 1, 0)],
        ]
    elif lattice_type == "trigonal" and symbol == "R":
        return [[(0, 0, 1)], [(1, 0, 0), (0, 1, 0), (1, 1, 0)]]
    # elif lattice_type == 'rhombohedral':
    #    return [[(0, 0, 1)],
    #            [(1, 0, 0), (0, 1, 0), (-1, -1, 0)]]
    elif lattice_type == "cubic":
        return [
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
            [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)],
            [(1, -1, 0), (1, 1, 0), (0, 1, -1), (0, 1, 1), (-1, 0, 1), (1, 0, 1)],
        ]
    else:
        return [[(0, 1, 0)]]


if __name__ == "__main__":
    print("Test pyxtal.wp.site symmetry")
    spg_list = [14, 36, 62, 99, 143, 160, 182, 191, 225, 230]
    for i in spg_list:
        g = Group(i)
        for wp in g:
            wp.get_site_symmetry()
            print(f"{wp.number:4d} {wp.get_label():10s} {wp.site_symm:10s}")

    print("Test pyxtal.wp.site symmetry representation")
    for i in spg_list:
        g = Group(i)
        for wp in g:
            if wp.index > 0:
                for idx in range(1):  # wp.multiplicity):
                    ss = wp.get_site_symmetry_object(idx)
                    print(
                        f"\n{wp.number:4d} {wp.get_label():10s} {ss.name:10s}",
                        ss.hm_symbols,
                    )
                    # ss.to_beautiful_matrix_representation(skip=True)
                    # print(ss.to_matrix_representation())
                    # print(ss.to_one_hot())

    print("Test pyxtal.wp.site space group")
    for i in spg_list:
        g = Group(i)
        print("\n", g.number, g.symbol)
        ss = g.get_spg_symmetry_object()
        ss.to_beautiful_matrix_representation()
        # matrix = ss.to_matrix_representation_spg()
        # print(matrix)
        # print(sum(sum(matrix)))
