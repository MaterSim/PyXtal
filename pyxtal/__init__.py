"""
PyXtal: A Python library for crystal structure generation and manipulation.

This module initializes the pyxtal package, providing essential imports and
utility functions.
"""

# Standard Libraries
import itertools
import json
from copy import deepcopy

# Third-Party Libraries
import numpy as np
from ase import Atoms
from pymatgen.core.structure import Molecule, Structure

# PyXtal Modules
from pyxtal.block_crystal import block_crystal
from pyxtal.crystal import random_crystal
from pyxtal.io import read_cif, structure_from_ext, write_cif
from pyxtal.lattice import Lattice
from pyxtal.molecule import pyxtal_molecule
from pyxtal.operations import SymmOp, apply_ops, get_inverse
from pyxtal.representation import representation, representation_atom
from pyxtal.symmetry import Group, Wyckoff_position
from pyxtal.tolerance import Tol_matrix
from pyxtal.version import __version__
from pyxtal.viz import display_atomic, display_cluster, display_molecular
from pyxtal.wyckoff_site import atom_site, mol_site
from pyxtal.wyckoff_split import wyckoff_split

# Uncomment the following line to set the package name
# name = "pyxtal"


def print_logo():
    """
    print out the logo and version information
    """

    print(
        r"""
             ______       _    _          _
            (_____ \     \ \  / /        | |
             _____) )   _ \ \/ / |_  ____| |
            |  ____/ | | | )  (|  _)/ _  | |
            | |    | |_| |/ /\ \ |_( (_| | |___
            |_|     \__  /_/  \_\___)__|_|_____)
                   (____/      """
    )
    print("\n")
    print("-------------------(version", __version__, ")--------------------\n")
    print("A Python package for random crystal generation")
    print("The source code is available at https://github.com/MaterSim/PyXtal")
    print("by Zhu's group at the University of North Carolina at Charlotte\n\n")


class pyxtal:
    """
    Class for handling atomic crystals based on symmetry constraints

    Examples
    --------

    To create a new structure instance

    >>> from pyxtal import pyxtal
    >>> struc = pyxtal()

    one can either use pyxtal to generate a random symmetric structure

    >>> struc.from_random(3, 227, ['C'], [8])

    or load a structure from a file or ASE atoms or Pymatgen structure object

    >>> struc.from_seed('diamond.cif')     # as a string
    >>> struc.from_seed(diamond_ase)       # as a ase atoms object
    >>> struc.from_seed(diamond_pymatgen)  # as a pymatgen structure object

    as long as the struc is created, you can check their symmetry as follows

    >>> struc.get_site_labels()
    {'C': ['8a']}
    >>> struc
    ------Crystal from random------
    Dimension: 3
    Composition: C8
    Group: Fd-3m (227)
    cubic lattice:   4.3529   4.3529   4.3529  90.0000  90.0000  90.0000
    Wyckoff sites:
         C @ [0.1250 0.1250 0.1250], WP:  8a, Site symmetry: -4 3 m

    The structure object can be easily manipulated via `apply_perturbtation`
    or `subgroup` function

    >>> struc2 = struc.subgroup_once(H=141)
    >>> struc2
    ------Crystal from Wyckoff Split------
    Dimension: 3
    Composition: C8
    Group: I41/amd (141)
    tetragonal lattice:   3.3535   3.3535   4.6461  90.0000  90.0000  90.0000
    Wyckoff sites:
             C @ [0.0000 0.2500 0.3750], WP:  4b, Site symmetry: -4 m 2

    Alternatively, one can easily compute its XRD via the `pyxtal.XRD` class

    >>> xrd = struc.get_XRD()
    >>> xrd
      2theta     d_hkl     hkl       Intensity  Multi
      32.706     2.738   [ 1  1  1]   100.00        8
      54.745     1.677   [ 2  2  0]    40.95       12
      65.249     1.430   [ 3  1  1]    20.65       24
      81.116     1.186   [ 4  0  0]     5.15        6
      90.236     1.088   [ 3  3  1]     8.24       24
     105.566     0.968   [ 4  2  2]    14.44       24
     115.271     0.913   [ 5  1  1]    10.03       24
     133.720     0.838   [ 4  4  0]     9.80       12
     148.177     0.802   [ 5  3  1]    28.27       48

    One can also try to get the transition path between two pyxtals
    that are symmetry related via the `get_transition` function

    >>> s1 = pyxtal()
    >>> s2 = pyxtal()
    >>> s1.from_seed("pyxtal/database/cifs/0-G62.cif") #structure with low symmetry
    >>> s2.from_seed("pyxtal/database/cifs/2-G71.cif") #structure with high symmetry
    >>> strucs, _, _, _, _ = s2.get_transition(s1) # get the transition from high to low
    >>> strucs
    [
    ------Crystal from Transition 0  0.000------
    Dimension: 3
    Composition: O24Mg4W4Pb8
    Group: Pnma (62)
    orthorhombic lattice:  11.6075   8.0526   5.8010  90.0000  90.0000  90.0000
    Wyckoff sites:
        Mg @ [ 0.8750  0.2500  0.7500], WP [4c] Site [.m.]
        Pb @ [ 0.6406  0.0053  0.7856], WP [8d] Site [1]
        W @ [ 0.6119  0.2500  0.2483], WP [4c] Site [.m.]
        O @ [ 0.6292  0.0083  0.2235], WP [8d] Site [1]
        O @ [ 0.4966  0.2500  0.0093], WP [4c] Site [.m.]
        O @ [ 0.5055  0.2500  0.4897], WP [4c] Site [.m.]
        O @ [ 0.7308  0.2500  0.9717], WP [4c] Site [.m.]
        O @ [ 0.7467  0.2500  0.4570], WP [4c] Site [.m.],
    ------Crystal from Transition 1  0.323------
    Dimension: 3
    Composition: O24Mg4W4Pb8
    Group: Pnma (62)
    orthorhombic lattice:  11.6020   8.0526   5.8038  90.0000  90.0000  90.0000
    Wyckoff sites:
        Mg @ [ 0.8750  0.2500  0.7500], WP [4c] Site [.m.]
        Pb @ [ 0.6250 -0.0053  0.7500], WP [8d] Site [1]
        W @ [ 0.6250  0.2500  0.2500], WP [4c] Site [.m.]
        O @ [ 0.6250  0.0083  0.2500], WP [8d] Site [1]
        O @ [ 0.5158  0.2500 -0.0068], WP [4c] Site [.m.]
        O @ [ 0.5158  0.2500  0.5068], WP [4c] Site [.m.]
        O @ [ 0.7342  0.2500  0.9932], WP [4c] Site [.m.]
        O @ [ 0.7342  0.2500  0.5068], WP [4c] Site [.m.]]

    Finally, the structure can be saved to different formats

    >>> struc.to_file('my.cif')
    >>> struc.to_file('my_poscar', fmt='poscar')

    or to Pymatgen/ASE structure object

    >>> pmg_struc = struc.to_pymatgen()
    >>> ase_struc = struc.to_ase()

    or to json file

    >>> struc.to_json('1.json')

    """

    def __init__(self, molecular=False, random_state=None):
        self.valid = False
        self.molecular = molecular
        self.standard_setting = True
        self.numIons = None
        self.numMols = None
        self.source = None
        self.formula = None
        self.species = None
        self.group = None
        self.lattice = None
        self.dim = 3
        self.factor = 1.0
        self.PBC = [1, 1, 1]
        self.random_state = np.random.default_rng(random_state)

        if molecular:
            self.molecules = []
            self.mol_sites = []
        else:
            self.atom_sites = []

    def __str__(self):
        if self.valid:
            s = f"\n------Crystal from {self.source:s}------"
            s += f"\nDimension: {self.dim}"
            s += f"\nComposition: {self.formula}"
            # if not self.standard_setting:
            if self.dim < 3:
                symbol = self.group.symbol
            else:
                if self.molecular:
                    if len(self.mol_sites) > 0:
                        symbol = self.mol_sites[0].wp.get_hm_symbol()
                else:
                    if len(self.atom_sites) > 0:
                        symbol = self.atom_sites[0].wp.get_hm_symbol()
            s += f"\nGroup: {symbol} ({self.group.number})"
            s += f"\n{self.lattice}"
            s += "\nWyckoff sites:"
            if self.molecular:
                for wyc in self.mol_sites:
                    s += f"\n\t{wyc}"
            else:
                for wyc in self.atom_sites:
                    s += f"\n\t{wyc}"
        else:
            s = "\nStructure not available."
        return s

    def __repr__(self):
        return str(self)

    def get_dof(self):
        """
        Get the number of dof for the given structures:
        """
        sites = self.mol_sites if self.molecular else self.atom_sites
        dof = 0
        for site in sites:
            dof += site.dof

        return self.lattice.dof + dof

    def get_bounds(self, vec=(2.0, 50.0), ang=(30, 150)):
        """
        Get the number of dof for the given structures:
        """
        bounds = self.lattice.get_bounds(vec[0], vec[1], ang[0], ang[1])
        sites = self.mol_sites if self.molecular else self.atom_sites
        for site in sites:
            bounds.extend(site.get_bounds())
        return bounds

    def get_site_labels(self):
        """
        Get the site_labels as a dictionary
        """
        if self.molecular:
            sites = self.mol_sites
            names = [site.molecule.name for site in sites]
        else:
            sites = self.atom_sites
            names = [site.specie for site in sites]

        dicts = {}
        for name, site in zip(names, sites):
            label = site.wp.get_label()
            if name not in dicts:
                dicts[name] = [label]
            else:
                dicts[name].append(label)
        return dicts

    def from_random(
        self,
        dim=3,
        group=None,
        species=None,
        numIons=None,
        factor=1.1,
        thickness=None,
        area=None,
        lattice=None,
        sites=None,
        conventional=True,
        t_factor=1.0,
        max_count=10,
        torsions=None,
        force_pass=False,
        block=None,
        num_block=None,
        seed=None,
        random_state=None,
        tm=None,
        use_hall=False,
    ):
        """
        The main function to generate random crystals.
        It calls `block_crysstal` or `random_crystal` internally.


        Args:
            dim (int): dimenion (0, 1, 2, 3)
            group (int): the group number (1-56, 1-75, 1-80, 1-230)
            species (list): a list of atomic symbols for each ion type, e.g., `["Ti", "O"]`
            numIons (list): a list of the number of each type of atom within the
                primitive cell (NOT the conventional cell), e.g., `[4, 2]`
            factor (optional): volume factor used to generate the crystal
            thickness:
            area (float):
            lattice (optional): `Lattice <pyxtal.lattice.Lattice.html>`_ to define the cell
            sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
            conventional (bool): whether or not use the conventional setting
            t_factor ():
            max_count ():
            tm (optional): `Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_ object
                to define the distances

        """
        prototype = "molecular" if self.molecular else "atomic"

        if tm is None:
            tm = Tol_matrix(prototype=prototype, factor=t_factor)

        count = 0
        quit = False

        while True:
            count += 1
            if self.molecular:
                struc = block_crystal(
                    dim,
                    group,
                    species,
                    numIons,
                    factor,
                    thickness=thickness,
                    area=area,
                    block=block,
                    num_block=num_block,
                    lattice=lattice,
                    torsions=torsions,
                    sites=sites,
                    conventional=conventional,
                    tm=tm,
                    seed=seed,
                    random_state=random_state,
                    use_hall=use_hall,
                )
            else:
                struc = random_crystal(
                    dim=dim,
                    group=group,
                    species=species,
                    numIons=numIons,
                    factor=factor,
                    thickness=thickness,
                    area=area,
                    lattice=lattice,
                    sites=sites,
                    conventional=conventional,
                    tm=tm,
                    use_hall=use_hall,
                    random_state=random_state,
                )
            if force_pass or struc.valid:
                quit = True
                break

            if count >= max_count:
                raise RuntimeError(
                    "long time to generate structure, check inputs")

        if quit:
            self.valid = struc.valid
            self.dim = dim
            try:
                self.lattice = struc.lattice
                if self.molecular:
                    self.numMols = struc.numMols
                    self.molecules = struc.molecules
                    self.mol_sites = struc.mol_sites
                    wp = self.mol_sites[0].wp
                    self.standard_setting = wp.is_standard_setting()
                else:
                    self.numIons = struc.numIons
                    self.species = struc.species
                    self.atom_sites = struc.atom_sites
                self.group = struc.group
                self.PBC = struc.PBC
                self.source = "random"
                self.factor = struc.factor
                self._get_formula()
            except:
                pass

            return struc

    def from_seed(
        self,
        seed,
        molecules=None,
        tol=1e-4,
        a_tol=5.0,
        ignore_HH=True,
        add_H=False,
        backend="pymatgen",
        style="pyxtal",
        hn=None,
        standard=False,
    ):
        """
        Load the seed structure from Pymatgen/ASE/POSCAR/CIFs
        Internally they will be handled by Pymatgen

        Args:
            seed: cif/poscar file or a Pymatgen Structure object
            molecules: a list of reference molecule (xyz file or Pyxtal molecule)
            tol: scale factor for covalent bond distance
            ignore_HH: whether or not ignore short H-H distance in molecules
            add_H: whether or not add H atoms
            backend: structure parser, default is pymatgen
            style: pyxtal for spglib
            standard: whether or not optimize lattice
        """

        if self.molecular:
            pmols = []
            for mol in molecules:
                if type(mol) == pyxtal_molecule:
                    pmols.append(mol)
                else:
                    pmols.append(pyxtal_molecule(mol, fix=True))
            # QZ: the default will not work for molecular H2, which is rare!
            struc = structure_from_ext(
                seed, pmols, ignore_HH=ignore_HH, add_H=add_H, hn=hn)
            self.mol_sites = struc.make_mol_sites()
            self.group = Group(struc.wyc.number)
            self.lattice = struc.lattice
            self.molecules = pmols
            self.numMols = struc.numMols
            self.standard_setting = True
            self.valid = True  # Need to add a check function
            if not standard:
                self.optimize_lattice()
        else:
            if isinstance(seed, dict):
                self.from_dict()
            elif isinstance(seed, Atoms):  # ASE atoms
                # from pymatgen.io.ase import AseAtomsAdaptor
                # pmg_struc = AseAtomsAdaptor.get_structure(seed)
                from pyxtal.util import ase2pymatgen

                pmg_struc = ase2pymatgen(seed)
                self._from_pymatgen(pmg_struc, tol, a_tol, style=style)
            elif isinstance(seed, Structure):  # Pymatgen
                self._from_pymatgen(seed, tol, style=style)
            elif isinstance(seed, str):
                if backend == "pymatgen":
                    pmg_struc = Structure.from_file(seed, primitive=True)
                    self._from_pymatgen(pmg_struc, tol, a_tol, style=style)
                else:
                    # Need to check
                    self.lattice, self.atom_sites = read_cif(seed)
                    wp = self.atom_sites[0].wp
                    self.group = Group(wp.hall_number, use_hall=True)
                    self.standard_setting = wp.is_standard_setting()
                    self.valid = True
        self.factor = 1.0
        self.source = "Seed"
        self.dim = 3
        self.PBC = [1, 1, 1]
        self._get_formula()

    def _from_pymatgen(self, struc, tol=1e-3, a_tol=5.0, style="pyxtal", hn=None):
        """
        Load structure from Pymatgen
        should not be used directly

        Args:
            struc: input pymatgen structure
            tol: symmetry tolerance
            a_tol: angle tolerance
            style: 'pyxtal' or spglib, differing in the choice of origin
            hn: hall_number
        """
        from pyxtal.util import get_symmetrized_pmg
        # import pymatgen.analysis.structure_matcher as sm

        self.valid = True
        try:
            sym_struc, number = get_symmetrized_pmg(
                struc, tol, a_tol, style, hn)
            # print(sym_struc)
            # import sys; sys.exit()
        except TypeError:
            print("Failed to load the Pymatgen structure")
        #    print(struc)
        #    self.valid = False

        if self.valid:
            d = sym_struc.composition.as_dict()
            species = list(d.keys())
            numIons = [int(d[ele]) for ele in species]
            self.numIons = numIons
            self.species = species
            if hn is None:
                self.group = Group(number, style=style)
            else:
                self.group = Group(hn, use_hall=True)
            # print(self.group[0]); import sys; sys.exit()
            matrix, ltype = sym_struc.lattice.matrix, self.group.lattice_type
            self.lattice = Lattice.from_matrix(matrix, ltype=ltype)
            atom_sites = []
            for i, site in enumerate(sym_struc.equivalent_sites):
                pos = site[0].frac_coords
                letter = sym_struc.wyckoff_symbols[i]
                # print(letter)
                wp = Wyckoff_position.from_group_and_letter(
                    number, letter, style=style, hn=hn)
                specie = site[0].specie.number
                # if wp.index>0: print(wp)
                pos1 = wp.search_generator(pos, self.group[0], tol=tol)
                if pos1 is not None:
                    atom_sites.append(atom_site(wp, pos1, specie))
                else:
                    pos1, wp, _ = self.group[0].merge(pos, matrix, 1e-3)
                    if pos1 is None:
                        print("Problem in ", site)
                        print(wp.symbol, wp.number, wp.letter, wp)
                        print("sym_struc_sites", sym_struc)
                        raise RuntimeError(
                            "Cannot extract the right mapping from spglib")
                        # break
                    else:
                        atom_sites.append(atom_site(wp, pos1, specie))
                        #print(pos, pos1, self.group[0])

            if len(atom_sites) != len(sym_struc.equivalent_sites):
                raise RuntimeError("Fail to extract the atom site")
            self.atom_sites = atom_sites
            # import pymatgen.analysis.structure_matcher as sm
            # self.dim = 3
            # self.PBC = [1, 1, 1]
            # pmg1 = self.to_pymatgen()
            # if not sm.StructureMatcher().fit(struc, pmg1):
            #    raise RuntimeError("The structure is inconsistent after conversion")

    def are_valid_numIons(self):
        """
        Check if the numIons are correct for debugging
        """

        N = [s.wp.multiplicity for s in self.atom_sites]
        # for s in self.atom_sites: print(s.wp.multiplicity, s.wp.get_label())
        if sum(N) != sum(self.numIons):
            print("Inconsistent numIons", sum(N), self.numIons)
            return False
        else:
            return True

    def check_H_coordination(self, r=1.12):
        """
        A function to check short if H is connected to more than one atom
        Mainly used for debug, powered by pymatgen

        Args:
            r: the given cutoff distances

        Returns:
            True or False
        """
        if self.dim > 0:
            pmg_struc = self.to_pymatgen()
            res = pmg_struc.get_all_neighbors(r)
            for i, neighs in enumerate(res):
                if pmg_struc.sites[i].specie.number == 1 and len(neighs) > 1:
                    return True
        else:
            raise NotImplementedError("Does not support cluster for now")
        return False

    def check_short_distances(self, r=0.7, exclude_H=True):
        """
        A function to check short distance pairs
        Mainly used for debug, powered by pymatgen

        Args:
            r: the given cutoff distances
            exclude_H: whether or not exclude the H atoms

        Returns:
            list of pairs within the cutoff
        """
        if self.dim > 0:
            pairs = []
            pmg_struc = self.to_pymatgen()
            if exclude_H:
                pmg_struc.remove_species("H")
            res = pmg_struc.get_all_neighbors(r)
            pairs = [
                [pmg_struc.sites[i].specie, n.specie, n.nn_distance] for i, neighs in enumerate(res) for n in neighs
            ]
        else:
            raise NotImplementedError("Does not support cluster for now")
        return pairs

    def check_short_distances_by_dict(self, dicts):
        """
        A function to check short distance pairs
        Mainly used for debug, powered by pymatgen

        Args:
            dicts: e.g., {"H-H": 1.0, "O-O": 2.0}

        Returns:
            N_pairs: number of atomic pairs within the cutoff
        """
        if self.dim > 0:
            N_pairs = 0
            r_cut = max([dicts[key] for key in dicts])
            pmg_struc = self.to_pymatgen()
            res = pmg_struc.get_all_neighbors(r_cut)
            for i, neighs in enumerate(res):
                ele1 = pmg_struc.sites[i].specie.value
                for n in neighs:
                    ele2 = n.specie.value
                    key1 = ele1 + "-" + ele2
                    key2 = ele2 + "-" + ele1
                    if (
                        key1 in dicts
                        and n.nn_distance < dicts[key1]
                        or (key1 != key2 and key2 in dicts and n.nn_distance < dicts[key2])
                    ):
                        N_pairs += 1
        else:
            raise NotImplementedError("Does not support cluster for now")
        return N_pairs

    def to_file(
        self,
        filename=None,
        fmt=None,
        permission="w",
        sym_num=None,
        header="from_pyxtal",
    ):
        """
        Creates a file with the given ame and type to store the structure.
        By default, creates cif files for crystals and xyz files for clusters.
        For other formats, Pymatgen is used

        Args:
            filename (string): the file path
            fmt (string): the file type (`cif`, `xyz`, etc.)
            permission (string): `w` or `a+`
            sym_num (int): number of sym_ops, None means writing all symops
            header (string): header

        Returns:
            Nothing. Creates a file at the specified path
        """
        if self.valid:
            if fmt is None:
                fmt = "xyz" if self.dim == 0 else "cif"

            if fmt == "cif":
                if self.dim == 3:
                    return write_cif(self, filename, header, permission, sym_num=sym_num)
                else:
                    pmg_struc = self.to_pymatgen()
                    if self.molecular:
                        pmg_struc.sort()
                return pmg_struc.to(fmt=fmt, filename=filename)
            else:
                pmg_struc = self.to_pymatgen()
                if self.molecular:
                    pmg_struc.sort()
                return pmg_struc.to(fmt=fmt, filename=filename)
        else:
            raise RuntimeError(
                "Cannot create file: structure did not generate")

    def supergroup(self, G=None, d_tol=1.0):
        """
        Generate a structure with higher symmetry

        Args:
            G: super space group number (list of integers)
            d_tol: maximum tolerance

        Returns:
            a list of pyxtal structures with minimum super group symmetries
        """

        from pyxtal.supergroup import supergroup

        my_super = supergroup(self, G=G)
        solutions = my_super.search_supergroup(d_tol=d_tol)
        return my_super.make_supergroup(solutions)

    def supergroups(self, G=None, d_tol=1.0):
        """
        Generate the structures with higher symmetry

        Args:
            G: super space group number (list of integers)
            d_tol: maximum tolerance

        Returns:
            a list of pyxtal structures with minimum super group symmetries
        """

        from pyxtal.supergroup import supergroups

        sup = supergroups(self, G=G, d_tol=d_tol)
        return sup.strucs

    def subgroup(
        self,
        perms=None,
        H=None,
        eps=0.05,
        idx=None,
        group_type="t",
        max_cell=4,
        min_cell=0,
        N_groups=None,
    ):
        """
        Generate a structure with lower symmetry

        Args:
            perms: e.g., {"Si": "C"}
            H: space group number (int)
            eps: pertubation term (float)
            idx: list
            group_type (string): `t`, `k` or `t+k`
            max_cell (float): maximum cell reconstruction
            min_cell (float): maximum cell reconstruction
            max_subgroups (int): maximum number of trial subgroups
        Returns:
            a list of pyxtal structures with lower symmetries
        """

        ids, sites, t_types, k_types = self._get_subgroup_ids(
            H, group_type, idx, max_cell, min_cell)
        # randomly choose a subgroup from the available list
        if N_groups is not None and len(ids) >= N_groups:
            ids = self.random_state.choice(ids, N_groups)
            # print('max_sub_group', len(idx), max_subgroups)

        valid_splitters = []
        bad_splitters = []
        for idx in ids:
            gtype = (t_types + k_types)[idx]
            if gtype == "k":
                idx -= len(t_types)
            splitter = wyckoff_split(
                G=self.group, wp1=sites, idx=idx, group_type=gtype)

            if not splitter.error:
                if perms is None:
                    if splitter.valid_split:
                        special = False
                        if self.molecular:
                            for i in range(len(self.mol_sites)):
                                for ops in splitter.H_orbits[i]:
                                    if len(ops) < len(splitter.H[0]):
                                        special = True
                                        break
                        if not special:
                            valid_splitters.append(splitter)
                        else:
                            bad_splitters.append(splitter)
                    else:
                        bad_splitters.append(splitter)
                else:
                    # apply permuation
                    if len(splitter.H_orbits) == 1:
                        if len(splitter.H_orbits[0]) > 1:
                            valid_splitters.append(splitter)
                        else:
                            bad_splitters.append(splitter)
                    else:
                        valid_splitters.append(splitter)

        if len(valid_splitters) == 0:
            # print("try do one more step")
            new_strucs = []
            for splitter in bad_splitters:
                trail_struc = self._subgroup_by_splitter(splitter, eps=eps)
                if trail_struc is not None:
                    new_strucs.extend(trail_struc.subgroup(
                        perms, group_type=group_type))
            return new_strucs
        else:
            # print(len(valid_splitters), "valid_splitters are present")
            new_strucs = []
            for splitter in valid_splitters:
                # print(splitter)
                if perms is None:
                    new_struc = self._subgroup_by_splitter(splitter, eps=eps)
                else:
                    new_struc = self._apply_substitution(splitter, perms)
                # if not new_struc.are_valid_numIons(): print(new_struc); import sys; sys.exit()
                new_strucs.append(new_struc)
            return new_strucs

    def subgroup_by_path(self, gtypes, ids, eps=0, mut_lat=False):
        """
        Generate a structure with lower symmetry (for atomic crystals only)

        Args:
            g_types: ['t', 'k']
            idx: list of ids for the splitter
            eps: degree of displacement
            mut_lat: mutate the lattice of not

        Returns:
            a pyxtal structure with lower symmetries
            list of splitters
        """
        struc = self.copy()
        splitters = []
        G = self.group
        for g_type, id in zip(gtypes, ids):
            _sites = struc.mol_sites if self.molecular else struc.atom_sites
            sites = [site.wp.index for site in _sites]
            # print(G.number, id, g_type, sites)
            splitter = wyckoff_split(G, wp1=sites, idx=id, group_type=g_type)
            struc = struc._subgroup_by_splitter(
                splitter, eps=eps, mut_lat=mut_lat)
            if struc is None:
                return None
            splitters.append(splitter)
            G = splitter.H
        return struc, splitters

    def subgroup_once(
        self,
        eps=0.1,
        H=None,
        perms=None,
        group_type="t",
        max_cell=4,
        min_cell=0,
        mut_lat=True,
        ignore_special=False,
        verbose=False,
    ):
        """
        Generate a structure with lower symmetry (for atomic crystals only)

        Args:
            perms: e.g., {"Si": "C"}
            H: space group number (int)
            idx: list
            group_type: `t` or `k`
            max_cell: maximum cell reconstruction (float)

        Returns:
            a pyxtal structure with lower symmetries
        """
        idx, sites, t_types, k_types = self._get_subgroup_ids(
            H, group_type, None, max_cell, min_cell)
        # print("Test subgroup_once", self.group.number, idx, sites)

        if idx is None or len(idx) == 0:
            if verbose:
                msg = "Cannot find valid splitter (likely due to large cellsize)"
                print(msg)
            return None

        # Try 100 times to see if a valid split can be found
        count = 0
        while count < 100:
            id = self.random_state.choice(idx)
            gtype = (t_types + k_types)[id]
            if gtype == "k":
                id -= len(t_types)
            # print(self.group.number, sites, id, gtype, idx)
            splitter = wyckoff_split(
                G=self.group.number, wp1=sites, idx=id, group_type=gtype)
            if not splitter.error:
                if perms is not None:
                    if len(splitter.H_orbits) == 1:
                        if len(splitter.H_orbits[0]) > 1:
                            return self._apply_substitution(splitter, perms)
                        else:
                            # print("try to find the next subgroup")
                            trail_struc = self._subgroup_by_splitter(
                                splitter, eps=eps, mut_lat=mut_lat)
                            if trail_struc is not None:
                                multiple = sum(
                                    trail_struc.numIons) / sum(self.numIons)
                                max_cell = max([1, max_cell / multiple])
                                ans = trail_struc.subgroup_once(
                                    eps, H, perms, group_type, max_cell)
                                if ans.group.number > 1:
                                    return ans
                    else:
                        return self._apply_substitution(splitter, perms)
                else:
                    # print('permuation')
                    if splitter.valid_split:
                        special = False
                        if self.molecular:
                            for i in range(len(self.mol_sites)):
                                for ops in splitter.H_orbits[i]:
                                    if len(ops) < len(splitter.H[0]):
                                        special = True
                                        break
                        if ignore_special:
                            return self._subgroup_by_splitter(splitter, eps=eps, mut_lat=mut_lat)
                        else:
                            if not special:
                                return self._subgroup_by_splitter(splitter, eps=eps, mut_lat=mut_lat)
                    else:
                        # print("try to find the next subgroup")
                        trail_struc = self._subgroup_by_splitter(
                            splitter, eps=eps, mut_lat=mut_lat)
                        if trail_struc is not None:
                            multiple = sum(trail_struc.numIons) / \
                                sum(self.numIons)
                            max_cell = max([1, max_cell / multiple])
                            ans = trail_struc.subgroup_once(
                                eps, H, None, group_type, max_cell)
                            if ans.group.number > 1:
                                return ans
            count += 1

        print("Cannot find the splitter")
        return None

    def _apply_substitution(self, splitter, perms):
        """
        Apply the substitution
        """
        try:
            new_struc = self._subgroup_by_splitter(splitter)
        except:
            print(self)
            print(splitter)
            print(len(splitter.H_orbits), len(
                splitter.G2_orbits), len(self.atom_sites))
            self._subgroup_by_splitter(splitter)

        # Create a list of tuples (site_id, original_specie) for sites that can be substituted
        site_info = [
            (site_id, site.specie) for site_id, site in enumerate(new_struc.atom_sites) if site.specie in perms
        ]

        # TODO range isn't inclusive of end number is this intended?
        N = self.random_state.choice(
            range(1, len(site_info))) if len(site_info) > 1 else 1
        sub_indices = self.random_state.choice(
            len(site_info), N, replace=False)

        for index in sub_indices:
            site_id, original_specie = site_info[index]
            new_struc.atom_sites[site_id].specie = perms[original_specie]

        new_struc._get_formula()
        return new_struc

    def _get_subgroup_ids(self, H, group_type, idx, max_cell, min_cell):
        """
        Generate the subgroup dictionary
        """

        # transform from p21/n to p21/n, need to fix later, wp.get_transformation_to_std()
        if not self.standard_setting:
            self.optimize_lattice(standard=True)

        if H is not None:
            group_type = "k" if Group(
                H, quick=True).point_group == self.group.point_group else "t"
        t_types = []
        k_types = []
        if group_type == "t":
            dicts = self.group.get_max_t_subgroup()
            t_types = ["t"] * len(dicts["subgroup"])
        elif group_type == "k":
            dicts = self.group.get_max_k_subgroup()
            k_types = ["k"] * len(dicts["subgroup"])
        else:
            dicts = self.group.get_max_t_subgroup()
            dict2 = self.group.get_max_k_subgroup()
            t_types = ["t"] * len(dicts["subgroup"])
            k_types = ["k"] * len(dict2["subgroup"])
            for key in dicts:
                dicts[key].extend(dict2[key])

        Hs = dicts["subgroup"]
        trans = dicts["transformation"]

        if idx is None:
            idx = []
            if not self.molecular or self.group.number > 142:
                for i, tran in enumerate(trans):
                    if min_cell <= np.linalg.det(tran[:3, :3]) <= max_cell:
                        idx.append(i)
            else:
                # for molecular crystals, assume the cell does not change
                for i, tran in enumerate(trans):
                    tran = np.abs(tran[:3, :3])
                    good = True
                    # only accepts trans like [a, b, c] [b, c, a]
                    if self.group.number != 5 and abs(abs(np.linalg.det(tran)) - 1) > 1e-3:
                        # print(self.group.number, np.linalg.det(tran))
                        good = False
                    if good:
                        # print(np.linalg.det(tran), tran)
                        idx.append(i)
        else:
            for id in idx:
                if id >= len(Hs):
                    raise ValueError(
                        "The idx exceeds the number of possible splits")

        if H is not None:
            idx = [id for id in idx if Hs[id] == H]

        # if len(idx) == 0:
        #    raise RuntimeError("Cannot find the splitter")

        struc_sites = self.mol_sites if self.molecular else self.atom_sites

        sites = [site.wp.get_label() for site in struc_sites]

        return idx, sites, t_types, k_types

    def _subgroup_by_splitter(self, splitter, eps=0.05, mut_lat=False):
        """
        Transform the crystal to subgroup symmetry from a splitter object

        Args:
            splitter: wyckoff splitter object
            eps (float): maximum atomic displacement in Angstrom
            mut_lat (bool): whether or not mutate the lattice
        """
        # print(splitter)
        lat1 = np.dot(splitter.R[:3, :3].T, self.lattice.matrix)

        multiples = np.linalg.det(splitter.R[:3, :3])
        new_struc = self.copy()
        new_struc.group = splitter.H
        try:
            lattice = Lattice.from_matrix(
                lat1, ltype=new_struc.group.lattice_type)
        except:
            self.optimize_lattice()
            lat1 = np.dot(splitter.R[:3, :3].T, self.lattice.matrix)
            try:
                lattice = Lattice.from_matrix(
                    lat1, ltype=new_struc.group.lattice_type)
            except:
                # print('problem with splitter, save it to bug.cif')
                # print(splitter)
                # print(self)
                # self.to_file('bug.cif')
                # import sys; sys.exit()
                return None
            # print(np.linalg.det(lat1))
            # print(self.lattice)
            # print(self.lattice.matrix)
            # print(splitter.R[:3,:3].T)
            # print(lat1); import sys; sys.exit()
        if mut_lat:
            try:
                lattice = lattice.mutate(degree=eps, frozen=False)  # True)
                # print('debug subgroup', mut_lat, eps, lattice)
            except:
                print("skip lattice mutation mostly for triclinic system")

        h = splitter.H
        split_sites = []
        if self.molecular:
            # below only works when the cell does not change
            # Fix when a, b, c swaps
            for i, site in enumerate(self.mol_sites):
                pos = site.position
                mol = site.molecule
                ori = site.orientation
                coord0 = np.dot(mol.mol.cart_coords, ori.matrix.T)

                wp1 = site.wp
                ori.reset_matrix(np.eye(3))
                id = 0
                for g1s, ops1, ops2 in zip(splitter.G1_orbits[i], splitter.G2_orbits[i], splitter.H_orbits[i]):
                    if site.wp.multiplicity == len(self.group[0]):
                        # general wyc
                        rot = g1s[0].affine_matrix[:3, :3].T
                    else:
                        # for special wyc, needs to get better treatment
                        op = wp1.get_euclidean_generator(
                            self.lattice.matrix, id)
                        rot = op.affine_matrix[:3, :3].T

                    # xyz in new lattice
                    # coord1 = np.dot(coord0, rot)
                    # coord1 = np.dot(coord1, splitter.inv_R[:3,:3].T)
                    # coord1 = np.array([np.dot(splitter.R[:3,:3].T, coord) for coord in coord1])
                    frac = np.dot(np.dot(coord0, self.lattice.inv_matrix), rot)
                    frac = np.dot(frac, splitter.inv_R[:3, :3].T)
                    coord1 = np.dot(frac, lattice.matrix)

                    _mol = mol.copy()
                    center = _mol.get_center(coord1)
                    _mol.reset_positions(coord1 - center)
                    # print('============'); print(_mol.mol.to(fmt='xyz'))
                    pos0 = apply_ops(pos, ops1)[0]  # ; print(pos0, pos)
                    pos0 -= np.floor(pos0)
                    dis = (np.random.sample(3) - 0.5).dot(self.lattice.matrix)
                    dis /= np.linalg.norm(dis)
                    pos0 += eps * dis * (np.random.random() - 0.5)
                    wp = Wyckoff_position.from_symops(ops2, h)
                    _site = mol_site(_mol, pos0, ori, wp, lattice)
                    _site.type = site.type
                    split_sites.append(_site)
                    id += wp.multiplicity
            new_struc.mol_sites = split_sites
            new_struc.numMols = [int(round(multiples * numMol))
                                 for numMol in self.numMols]

        else:
            for i, site in enumerate(self.atom_sites):
                pos = site.position
                for ops1, ops2 in zip(splitter.G2_orbits[i], splitter.H_orbits[i]):
                    pos0 = apply_ops(pos, ops1)[0]
                    pos0 -= np.floor(pos0)
                    dis = (np.random.sample(3) - 0.5).dot(self.lattice.matrix)
                    dis /= np.linalg.norm(dis)
                    pos0 += np.dot(eps * dis * (np.random.random() -
                                   0.5), self.lattice.inv_matrix)
                    wp = Wyckoff_position.from_symops(ops2, h)
                    split_sites.append(atom_site(wp, pos0, site.specie))

            new_struc.atom_sites = split_sites
            new_struc.numIons = [int(round(multiples * numIon))
                                 for numIon in self.numIons]
        new_struc.lattice = lattice
        new_struc.source = "subgroup"

        return new_struc

    def apply_perturbation(self, d_lat=0.05, d_coor=0.05, d_rot=1):
        """
        perturb the structure without breaking the symmetry

        Args:
            d_coor: magnitude of perturbation on atomic coordinates (in A)
            d_lat: magnitude of perturbation on lattice (in percentage)
        """

        self.lattice = self.lattice.mutate(degree=d_lat)

        if self.molecular:
            for site in self.mol_sites:
                site.perturbate(lattice=self.lattice.matrix,
                                trans=d_coor, rot=d_rot)
        else:
            for site in self.atom_sites:
                site.perturbate(lattice=self.lattice.matrix, magnitude=d_coor)

        self.source = "Perturbation"

    def copy(self):
        """
        Simply copy the structure
        """
        return deepcopy(self)

    def _get_coords_and_species(self, absolute=False, unitcell=True):
        """
        Extract the coordinates and species information

        Args:
            abosulte: if True, return the cartesian coords otherwise fractional

        Returns:
            total_coords (N*3 numpy array) and the list of species
        """
        species = []
        total_coords = None
        if self.molecular:
            for site in self.mol_sites:
                coords, site_species = site.get_coords_and_species(
                    absolute, unitcell=unitcell)
                species.extend(site_species)
                total_coords = coords if total_coords is None else np.append(
                    total_coords, coords, axis=0)
        else:
            for site in self.atom_sites:
                species.extend([site.specie] * site.multiplicity)
                total_coords = site.coords if total_coords is None else np.append(
                    total_coords, site.coords, axis=0)

            if absolute:
                total_coords = total_coords.dot(self.lattice.matrix)

        return total_coords, species

    def _get_formula(self):
        """
        A quick function to get the formula.
        """
        from pyxtal.database.element import Element

        formula = ""

        if self.molecular:
            numspecies = self.numMols
            species = [str(mol) for mol in self.molecules]
        else:
            specie_list = []
            for site in self.atom_sites:
                specie_list.extend([site.specie] * site.wp.multiplicity)
            if self.species is None:
                species = list(set(specie_list))
                self.species = species
            else:
                species = self.species

            numIons = np.zeros(len(species), dtype=int)
            for i, sp in enumerate(species):
                numIons[i] = specie_list.count(sp)
            self.numIons = numIons
            numspecies = self.numIons

        for i, s in zip(numspecies, species):
            specie = Element(s).short_name if isinstance(s, int) else s
            formula += f"{specie:s}{int(i):d}"

        self.formula = formula

    def get_zprime(self, integer=False):
        """
        Get zprime for molecular xtal
        """
        mult = len(self.group[0])
        comp = [c / mult for c in self.numMols]
        if integer:
            comp = [int(np.ceil(c)) for c in comp]
        return comp

    def get_1D_comp(self):
        """
        Get composition for 1d rep of molecular xtal
        """
        comp = [0] * len(self.molecules)
        for s in self.mol_sites:
            for i, m in enumerate(self.molecules):
                if s.molecule.name == m.name:
                    comp[i] += 1
        return comp

    def get_num_torsions(self):
        """
        Get number of torsions for molecular xtal
        """
        N_torsions = 0
        for s in self.mol_sites:
            N_torsions += len(s.molecule.torsionlist)
        return N_torsions

    def to_ase(self, resort=True, center_only=False, add_vaccum=True):
        """
        Export to ase Atoms object.

        Args:
            resort (bool): resort atoms
            center (bool): dump only molculer center for mol. xtal
            add_vaccum (bool): add vaccum for 0/1/2D systems
        """
        if self.valid:
            if self.dim > 0:
                lattice = self.lattice.copy()
                if self.molecular:
                    if center_only:
                        coords, species = [], []
                        for site in self.mol_sites:
                            _coords = site.wp.apply_ops(site.position)
                            _coords -= np.floor(_coords)
                            coords.extend(_coords.dot(self.lattice.matrix))
                            species.extend([site.type + 1] *
                                           site.wp.multiplicity)
                    else:
                        coords, species = self._get_coords_and_species(True)
                    if add_vaccum:
                        latt, coords = lattice.add_vacuum(
                            coords, frac=False, PBC=self.PBC)
                    else:
                        latt = lattice.matrix
                    atoms = Atoms(species, positions=coords,
                                  cell=latt, pbc=self.PBC)
                else:
                    coords, species = self._get_coords_and_species()
                    coords -= np.floor(coords)
                    if add_vaccum:
                        latt, coords = lattice.add_vacuum(coords, PBC=self.PBC)
                    else:
                        latt = lattice.matrix
                    atoms = Atoms(species, scaled_positions=coords,
                                  cell=latt, pbc=self.PBC)
                if resort:
                    permutation = np.argsort(atoms.numbers)
                    atoms = atoms[permutation]
                return atoms

            else:
                coords, species = self._get_coords_and_species(True)
                return Atoms(species, positions=coords)
        else:
            raise RuntimeError("No valid structure can be converted to ase.")

    def to_pymatgen(self, resort=True, shape="upper"):
        """
        Export to Pymatgen structure object.
        """

        if self.valid:
            if self.dim > 0:
                lattice = self.lattice.copy()
                # lattice.reset_matrix(shape)
                coords, species = self._get_coords_and_species()
                if resort:
                    permutation = sorted(
                        range(len(species)), key=species.__getitem__)
                    # permutation = np.argsort(species)
                    species = [species[id] for id in permutation]
                    coords = coords[permutation]
                # Add space above and below a 2D or 1D crystals
                latt, coords = lattice.add_vacuum(coords, PBC=self.PBC)
                return Structure(latt, species, coords)
            else:
                # Clusters are handled as large molecules
                coords, species = self._get_coords_and_species(True)
                return Molecule(species, coords)
        else:
            raise RuntimeError(
                "No valid structure can be converted to pymatgen.")

    def to_pyxtal_center(self):
        """
        Export to PyXtal object for molecular centers only.
        """

        if self.valid and self.molecular:
            new_struc = pyxtal()
            new_struc.lattice = self.lattice.copy()
            newsites = []
            for site in self.mol_sites:
                for i, mol in enumerate(self.molecules):
                    if mol.name == site.molecule.name:
                        break
                newsites.append(atom_site(site.wp, site.position, i + 1))
            new_struc.atom_sites = newsites
            new_struc.group = self.group
            new_struc.standard_setting = site.wp.is_standard_setting()
            new_struc.numIons = self.numMols
            new_struc.species = []
            for i in range(len(self.molecules)):
                new_struc.species.append(i + 1)
            new_struc.valid = True
            new_struc.factor = 1.0
            new_struc.source = "Mol. Center"
            new_struc.dim = self.dim
            new_struc.PBC = self.PBC
            new_struc._get_formula()
            return new_struc
        else:
            raise RuntimeError(
                "No valid structure can be converted to pymatgen.")

    def get_XRD(self, **kwargs):
        """
        Compute the PXRD object.

        ** kwargs include
            - wavelength (1.54184)
            - thetas [0, 180]
            - preferred_orientation: False
            - march_parameter: None
        """
        from pyxtal.XRD import XRD

        return XRD(self.to_ase(), **kwargs)

    def optimize_lattice(self, iterations=5, force=False, standard=False):
        """
        Optimize the lattice if the cell has a bad inclination angles
        We first optimize the angle to some good-looking range.
        If standard, force the structure to have the standard setting
        This only applies to monoclinic systems

        Args:
            iterations: maximum number of iterations
            force: whether or not do the early termination
            standard: True or False
        """
        for _i in range(iterations):
            lattice, trans, opt = self.lattice.optimize_once()
            # print(i, opt, self.lattice, "->", lattice)
            if force or opt:
                self.transform(trans)
            else:
                break

        # QZ: to check spglib
        if standard and 3 <= self.group.number <= 15 and not self.standard_setting:
            trans1 = self.lattice.get_permutation_matrices()
            trans2 = self.lattice.get_transformation_matrices()
            good_trans = None
            beta_diff = 90
            wp = self.mol_sites[0].wp if self.molecular else self.atom_sites[0].wp

            for tran1 in trans1:
                for tran2 in trans2:
                    _trans = [tran1, tran2]
                    wp0 = wp.copy()
                    lat0 = self.lattice.transform_multi(_trans)
                    wp0.transform_from_matrices(_trans)
                    beta_diff0 = abs(lat0.beta * 180 / np.pi - 90)
                    # print(wp0, wp0.is_standard_setting())
                    if wp0.is_standard_setting() and beta_diff0 < beta_diff:
                        good_trans = _trans
                        beta_diff = beta_diff0

            if good_trans is not None:
                for tran in good_trans:
                    self.transform(tran)
            else:
                print(self.lattice)
                if self.molecular:
                    print(self.mol_sites[0].wp)
                else:
                    print(self.atom_sites[0].wp)
                msg = "Cannot find the standard setting"
                raise RuntimeError(msg)

    def update_wyckoffs(self):
        """
        rescale the coordinates in the wyckoff site
        """

        sites = self.mol_sites if self.molecular else self.atom_sites

        for site in sites:
            site.update()

    def get_std_representation(self, trans):
        """
        Perform cell transformation so that the symmetry operations
        follow standard space group notation
        """
        pass

    def get_1D_representation(self, standard=False):
        """
        Get the 1D representation class for molecular crystals
        """
        if self.molecular:
            rep = representation.from_pyxtal(self)
        else:
            rep = representation_atom.from_pyxtal(self, standard=standard)
        return rep

    def transform(self, trans, lattice=None):
        """
        Perform cell transformation and symmetry operation

        Args:
            trans: 3*3 matrix
            lattice: pyxtal lattice object
            update: whether or not update each wp
        """
        # print(trans)
        if lattice is None:
            # print("perform cell transformation")
            lattice = self.lattice.transform(trans)

        lattice0 = lattice.copy()
        lattice0.reset_matrix()

        sites = self.mol_sites if self.molecular else self.atom_sites

        for j, site in enumerate(sites):
            pos_abs = np.dot(site.position, self.lattice.matrix)
            pos_frac = pos_abs.dot(lattice.inv_matrix)
            pos_frac -= np.floor(pos_frac)
            wp = site.wp.copy()
            wp.transform_from_matrix(trans, False, update=False)

            if self.molecular:
                # Obtain the transformed xyz
                xyz_abs, _ = site._get_coords_and_species(
                    absolute=True, first=True)
                xyz_frac = xyz_abs.dot(lattice.inv_matrix)
                xyz_abs = np.dot(xyz_frac, lattice0.matrix)

                mol = site.molecule
                ori = site.orientation
                ori.reset_matrix(np.eye(3))
                center = mol.get_center(xyz_abs)
                mol.reset_positions(xyz_abs - center)
                sites[j] = mol_site(mol, pos_frac, ori, wp,
                                    lattice0, site.type)
            else:
                sites[j] = atom_site(wp, pos_frac, site.specie)

        # update the hall number
        for i, site in enumerate(sites):
            if i == 0:
                match_spg, match_hall = site.wp.update()
                if not match_spg:
                    hall_numbers = [site.wp.hall_number]
                    # site.wp.update_hall(hall_numbers)
            else:
                if not match_spg:
                    site.wp.update_hall(hall_numbers)
                else:
                    site.wp.update_index()
        # reset the matrix
        wp0 = sites[0].wp
        self.lattice = lattice0
        self.standard_setting = wp0.is_standard_setting()
        self.group = Group(wp0.hall_number, use_hall=True)

    def to_json(self, filename="pyxtal.json"):
        """
        Save the model as a dictionary
        """
        from monty.json import MontyEncoder

        dict0 = self.save_dict()
        with open(filename, "w") as outfile:
            json.dump(dict0, outfile, cls=MontyEncoder)

    def from_json(self, filename):
        """
        Load the model from a json file
        """
        from monty.serialization import loadfn

        data = loadfn(filename)
        self.load_dict(data)

    def save_dict(self):
        """
        Save the model as a dictionary
        """
        if self.molecular:
            sites = [site.save_dict() for site in self.mol_sites]
        else:
            sites = [site.save_dict() for site in self.atom_sites]

        return {
            "lattice": self.lattice.matrix,
            "sites": sites,
            "group": self.group.number,
            "molecular": self.molecular,
            "numIons": self.numIons,
            "numMols": self.numMols,
            "factor": self.factor,
            "PBC": self.PBC,
            "formula": self.formula,
            "source": self.source,
            "dim": self.dim,
            "valid": self.valid,
        }

    def load_dict(self, dict0):
        """
        Load the structure from a dictionary
        """
        self.group = Group(dict0["group"], dict0["dim"])
        self.lattice = Lattice.from_matrix(
            dict0["lattice"], ltype=self.group.lattice_type)
        self.molecular = dict0["molecular"]
        self.factor = dict0["factor"]
        self.source = dict0["source"]
        self.dim = dict0["dim"]
        self.PBC = dict0["PBC"]
        self.numIons = dict0["numIons"]
        self.numMols = dict0["numMols"]
        self.valid = dict0["valid"]
        self.formula = dict0["formula"]

        if dict0["molecular"]:
            self.molecules = [None] * len(self.numMols)
            sites = [mol_site.load_dict(site) for site in dict0["sites"]]
            # TODO: this for loop makes repeated calls for duplicated molecules
            for msite in sites:
                if self.molecules[msite.type] is None:
                    self.molecules[msite.type] = msite.molecule
            self.mol_sites = sites
        else:
            sites = [atom_site.load_dict(site) for site in dict0["sites"]]
            self.atom_sites = sites

    def build(self, group, species, numIons, lattice, sites, tol=1e-2, dim=3, use_hall=False):
        """
        Build a atomic crystal based on the necessary input

        Args:
            group (int): 225
            species (list): ['Na', 'Cl']
            numIons (list): [4, 4]
            lattice: lattice object
            sites (list): [[{"4a": [0.0, 0.0, 0.0]}], [{"4b": [0.5, 0.5, 0.5]}]]
            dim (int):
            use_hll (bool):
        """

        from pyxtal.symmetry import choose_wyckoff

        if self.molecular:
            raise RuntimeError("Cannot support the molecular crystal")

        if type(group) == Group:
            self.group = group
        else:
            self.group = Group(group, dim=dim, use_hall=use_hall)

        # Lattica needs some special handling here
        if type(lattice) != Lattice:
            if type(lattice) == np.ndarray:
                ltype = self.group.lattice_type
                if len(lattice) == 3:
                    lattice = Lattice.from_matrix(lattice, ltype=ltype)
                elif len(lattice) == 6:  # cell para
                    [a, b, c, alpha, beta, gamma] = lattice
                    lattice = Lattice.from_para(
                        a, b, c, alpha, beta, gamma, ltype=ltype)
                else:
                    msg = "Cannot convert the input array to pyxtal.lattice.Lattice"
                    raise ValueError(msg, lattice)
            else:
                msg = "Cannot convert the input array to pyxtal.lattice.Lattice"
                raise ValueError(msg, lattice)
        # else:
        #    print(lattice, type(lattice), type(lattice)==Lattice, type(lattice)!=Lattice)
        #    msg = 'The input lattice needs to be a pyxtal.lattice.Lattice class'
        #    raise ValueError(msg, lattice)
        self.lattice = lattice
        self.dim = dim
        self.factor = 1.0
        self.PBC = self.group.PBC
        self.numIons = numIons
        self.species = species
        np.zeros(len(numIons), dtype=int)
        _sites = []

        if len(sites) != len(species):
            print(len(sites), len(species))
            raise RuntimeError("Inconsistency between sites and species")

        wp0 = self.group[0]
        for sp, wps in zip(species, sites):
            for wp in wps:
                if isinstance(wp, dict):  # dict
                    for pair in wp.items():
                        (key, pos) = pair
                        _wp = choose_wyckoff(self.group, site=key)
                        if _wp is not False:
                            pt = [0.0, 0.0, 0.0] if _wp.get_dof(
                            ) == 0 else _wp.get_all_positions(pos)[0]
                            _sites.append(atom_site(_wp, pt, sp))
                        else:
                            raise RuntimeError("Cannot interpret site", key)
                elif len(wp) == 4:  # tuple:
                    (key, x, y, z) = wp
                    _wp = choose_wyckoff(self.group, site=key, dim=dim)
                    #print('debug build', key, x, y, z, _wp.get_label())
                    if _wp is not False:
                        if _wp.get_dof() == 0:  # fixed pos
                            pt = [0.0, 0.0, 0.0]
                        else:
                            ans = _wp.get_all_positions([x, y, z])
                            pt = ans[0] if ans is not None else None
                            #print('debug build', ans, x, y, z)
                            # print('debug', ans)
                        if pt is not None:
                            _sites.append(atom_site(_wp, pt, sp))
                    else:
                        raise RuntimeError("Cannot interpret site", key)
                else:  # List of atomic coordinates
                    pt, _wp, _ = wp0.merge(wp, lattice.matrix, tol=tol)
                    _sites.append(atom_site(_wp, pt, sp))

        self.atom_sites = _sites
        self.standard_setting = True
        self.valid = True
        self.source = "Build"
        self._get_formula()

    def get_alternatives(self, include_self=True, same_letters=False, ref_lat=None, d_tol=2.0, f_tol=0.15):
        """
        Get alternative structure representations

        Args:
            include_self (bool): return the original structure
        Return:
            list of structures
        """
        if include_self:
            self.wyc_set_id = 0
            new_strucs = [self]
        else:
            new_strucs = []

        # the list of wyckoff indices in the original structure
        # e.g. [0, 2, 2, 4] -> [a, c, c, e]
        # ids = [len(self.group)-1-site.wp.index for site in self.atom_sites]

        wyc_sets = self.group.get_alternatives()
        No = len(wyc_sets["No."])
        letters = wyc_sets["Transformed WP"][0]
        if No > 1:
            # skip the first setting since it is identity
            for no in range(1, No):
                add = (wyc_sets["Transformed WP"][no] ==
                       letters) if same_letters else True
                if add:
                    new_struc = self._get_alternative(
                        wyc_sets, no, ref_lat, d_tol, f_tol)
                    if new_struc is not None:
                        new_strucs.append(new_struc)
        # print("Numbers===============", len(new_strucs)); import sys; sys.exit()
        return new_strucs

    def to_standard_setting(self):
        """
        A short cut to symmetrize the structure in the stardard setting
        """
        if self.molecular:
            pmg = self.to_pymatgen()
            self.from_seed(pmg, molecules=self.molecules, standard=True)

    def resort_species(self, species):
        """
        resort the atomic species

        Args:
            species: list of elements, e.g. ['Si', 'O']
        """
        sp1 = deepcopy(species)
        sp1.sort()
        sp2 = deepcopy(self.species)
        sp2.sort()
        if sp1 == sp2:
            self.species = species
            self.resort()
        else:
            ids = []
            for specie in species:
                for j, site in enumerate(self.atom_sites):
                    if site.specie == specie and j not in ids:
                        ids.append(j)
            self.atom_sites = [self.atom_sites[j] for j in ids]
            self.species = species
            self._get_formula()

    def resort(self):
        """
        A short cut to resort the sites by self.molecules or self.species
        """
        ids = []
        if self.molecular:
            for mol in self.molecules:
                for j, site in enumerate(self.mol_sites):
                    if site.molecule.name == mol.name and j not in ids:
                        ids.append(j)
            self.mol_sites = [self.mol_sites[j] for j in ids]
        else:
            for specie in self.species:
                for j, site in enumerate(self.atom_sites):
                    if site.specie == specie and j not in ids:
                        ids.append(j)
            self.atom_sites = [self.atom_sites[j] for j in ids]
            # print(self.atom_sites)

    def _get_alternative(self, wyc_sets, index, ref_lat=None, d_tol=2.0, f_tol=0.15):
        """
        Get alternative structure representations

        Args:
            wyc_sets: dictionary of `Coset Representative` and `Transformed WP`
            index: the index of target wyc_set
            ref_lat: a refernece lattice

        Returns:
            a new pyxtal structure after transformation
        """
        new_struc = self.copy()
        # xyz_string like 'x+1/4,y+1/4,z+1/4'
        xyz_string = wyc_sets["Coset Representative"][index]
        op = get_inverse(SymmOp.from_xyz_str(xyz_string))

        # transform lattice
        R = op.affine_matrix[:3, :3]  # rotation
        cell = self.lattice.matrix
        new_lat = Lattice.from_matrix(
            np.dot(R, cell), ltype=self.lattice.ltype)
        # matrix = new_lat.matrix
        if ref_lat is not None:
            d_tol1, f_tol1, a_tol1, switch = new_lat.get_diff(ref_lat)
            if (d_tol1 > d_tol and f_tol1 > f_tol) or (a_tol1 > 15.0) or switch:
                # print('bad setting', new_lat); print(ref_lat)
                return None

        # Lattice.from_matrix(matrix, ltype=self.group.lattice_type)
        new_struc.lattice = new_lat

        for i, site in enumerate(new_struc.atom_sites):
            id = len(self.group) - site.wp.index - 1
            letter = wyc_sets["Transformed WP"][index].split()[id]
            wp = Wyckoff_position.from_group_and_letter(
                self.group.number, letter)
            pos = op.operate(site.position)
            pos1 = wp.search_generator(pos, self.group[0])
            if pos1 is not None:
                new_struc.atom_sites[i] = atom_site(wp, pos1, site.specie)
            else:
                return None
                # print(pos)
                # print(wp)
                # raise RuntimeError("Cannot find the right pos")

        new_struc.source = f"Alt. Wyckoff Set [{index:d}]: {xyz_string:s}"
        new_struc.wyc_set_id = index

        return new_struc

    def _get_alternative_back(self, index):
        """
        Get alternative structure representations

        Args:
            index: the index of target wyc_set

        Returns:
            a new pyxtal structure after transformation
        """
        new_struc = self.copy()
        wyc_sets = self.group.get_alternatives()

        # xyz_string like 'x+1/4,y+1/4,z+1/4'
        xyz_string = wyc_sets["Coset Representative"][index]
        op = SymmOp.from_xyz_str(xyz_string)
        # op = get_inverse(SymmOp.from_xyz_str(xyz_string))
        letters = wyc_sets["Transformed WP"][0].split()
        letters1 = wyc_sets["Transformed WP"][index].split()

        for i, site in enumerate(new_struc.atom_sites):
            # id = len(self.group) - site.wp.index - 1
            letter1 = site.wp.letter
            letter = letters[letters1.index(letter1)]
            # print("transition", letter1, '->', letter)
            wp = Wyckoff_position.from_group_and_index(
                self.group.number, letter)
            pos = op.operate(site.position)
            pos1 = wp.search_generator(pos, self.group[0])
            if pos1 is not None:
                new_struc.atom_sites[i] = atom_site(wp, pos1, site.specie)
            else:
                print(pos)
                print(wp)
                raise RuntimeError("Cannot find the right pos")

        new_struc.source = f"Alt. Wyckoff Set [{index:d}]: {xyz_string:s}"
        new_struc.wyc_set_id = index

        # transform lattice
        R = op.affine_matrix[:3, :3]  # rotation
        matrix = np.dot(R, self.lattice.matrix)
        new_struc.lattice = Lattice.from_matrix(
            matrix, ltype=self.group.lattice_type)

        return new_struc

    def check_distance(self):
        """
        Check intermolecular distance for molecular crystal
        """
        if self.molecular:
            return all(ms.short_dist() for ms in self.mol_sites)
        return None

    def get_density(self):
        """
        Compute density
        """
        return self.to_pymatgen().density

    def has_special_site(self, species=None):
        """
        Check if the crystal has a special site
        """
        special = False
        if species is None:
            species = self.species

        sites = self.mol_sites if self.molecular else self.atom_sites

        for msite in sites:
            if self.molecular:
                if msite.wp.index > 0:
                    special = True
                    break
            else:
                if msite.specie in species and msite.wp.index > 0:
                    special = True
                    break
        return special

    def to_subgroup(self, path=None, t_only=True, iterate=False, species=None):
        """
        Transform a crystal with speical sites to subgroup
        represenatation with general sites

        Args:
            Path: list of path to get the general sites
            iterate (bool): whether or not do it iteratively
            species : to add
        """
        if not self.standard_setting:
            self.optimize_lattice(standard=True)
        if species is None:
            species = self.species

        # Compute the path is needed
        if path is None:
            if self.molecular:
                sites = self.mol_sites
                max_index = max([site.wp.index for site in sites])
            else:
                sites = self.atom_sites
                max_index = max(
                    [site.wp.index for site in sites if site.specie in species])
            # print([site.wp.index for site in sites])
            if self.molecular:
                path = self.group.short_path_to_general_wp(max_index, t_only)
            else:
                path = self.group.short_path_to_general_wp(max_index, t_only)
            # print(max_index, path)

        if path is not None:
            gtypes, ids = [], []
            for p in path:
                gtypes.append(p[0])
                ids.append(p[1])
            sub, _ = self.subgroup_by_path(gtypes, ids, eps=0)
            sub.optimize_lattice()
            sub.source = "subgroup"
        else:
            sub = self.copy()

        if iterate and sub.has_special_site():
            sub = sub.to_subgroup()

        return sub

    def show(self, **kwargs):
        """
        display the crystal structure
        """
        if self.molecular:
            return display_molecular(self, **kwargs)
        else:
            return display_atomic(self, **kwargs)

    def get_free_axis(self):
        """
        Check if the a, b, c axis have free parameters
        """
        free_axis = self.group.get_free_axis()

        for site in self.atom_sites:
            axis = site.wp.get_frozen_axis()
            for ax in axis:
                if ax in free_axis:
                    free_axis.remove(ax)
            if len(free_axis) == 0:
                break
        return free_axis

    def find_matched_lattice(self, ref_struc, d_tol=2.0, f_tol=0.15):
        """
        Compute the displacement w.r.t. the reference structure

        Args:
            ref_struc: reference pyxtal structure (assuming the same atomic ordering)
            d_tol: tolerence of mismatch in the absolute scale
            f_tol: tolerence of mismatch in the fractional scale

        Returns:
            ref_struc with matched lattice
        """
        ref_struc.optimize_lattice()
        l1 = self.lattice
        l2 = ref_struc.lattice
        # print(l1, l2)
        if self.group.number <= 15:
            # QZ: here we enumerate all possible transformations, maybe redundant
            trans_good, _ = l2.search_transformations(l1, d_tol, f_tol)
            # print(l1, l2, len(trans_good)); import sys; sys.exit()
            good_strucs = []

            for trans in trans_good:
                ref_struc0 = ref_struc.copy()
                for tran in trans:
                    ref_struc0.transform(tran)
                # print(ref_struc0, len(trans))
                if ref_struc0.standard_setting:
                    good_strucs.append(ref_struc0)
                # ============================== To remove
                # wp = ref_struc0.atom_sites[0].wp
                # pt = ref_struc0.atom_sites[0].position
                # if wp.is_standard_setting():
                #    good_strucs.append(ref_struc0)
                # else:
                #    valid, vector = wp.check_translation(pt)
                #    if valid:
                #        ref_struc0.translate(vector, reset_wp=True)
                #        ref_struc0.standard_setting = True
                #        good_strucs.append(ref_struc0)

            return good_strucs
        else:
            d_tol1, f_tol1, a_tol1, switch = l1.get_diff(l2)
            if d_tol1 > d_tol and f_tol1 > f_tol:
                return []
            else:
                return [ref_struc]

    def check_mapping(self, ref_struc):
        """
        Compute the displacement w.r.t. the reference structure

        Args:
            ref_struc: reference pyxtal structure (assuming the same atom order)

        Returns:
            True or False
        """
        orders = list(range(len(self.atom_sites)))
        atom_sites = self.atom_sites
        for site1 in atom_sites:
            match = False
            # search for the best match
            for i in orders:
                site2 = ref_struc.atom_sites[i]
                if site1.specie == site2.specie and site1.wp.index == site2.wp.index:
                    match = True
                    break
            if match:
                orders.remove(i)
            else:
                return False
        return True

    def get_disps_single(self, ref_struc, trans, d_tol=1.2):
        """
        Compute the displacement w.r.t. the reference structure

        Args:
            ref_struc: reference pyxtal structure (assuming the same atom order)
            trans: translation vector
            d_tol: tolerence of mismatch

        Returns:
            Atomic displacements in np.array
            translation:
        """
        cell1 = self.lattice.matrix
        disps = []
        orders = list(range(len(self.atom_sites)))

        atom_sites = self.atom_sites
        for site1 in atom_sites:
            match = False

            # search for the best match
            ds = []
            ids = []
            _disps = []
            for i in orders:
                site2 = ref_struc.atom_sites[i]
                if site1.specie == site2.specie and site1.wp.index == site2.wp.index:
                    disp, dist = site1.get_disp(site2.position, cell1, trans)
                    # strs = "{:2s} ".format(site1.specie)
                    # strs += "{:6.3f} {:6.3f} {:6.3f}".format(*site1.position)
                    # strs += " => {:6.3f} {:6.3f} {:6.3f} ".format(*site2.position)
                    # strs += "[{:6.3f} {:6.3f} {:6.3f}] {:6.3f}".format(*disp, dist)
                    # if dist < d_tol*1.2: strs += ' True'
                    # print(strs)

                    if dist < 0.3:
                        match = True
                        break

                    if dist < 1.2 * d_tol:
                        ds.append(dist)
                        ids.append(i)
                        _disps.append(disp)
                        # print("========", ds, ids, site1.specie, site2.specie)
            if match:
                disps.append(disp)
                orders.remove(i)
            else:
                if len(ds) > 0:
                    ds = np.array(ds)
                    id = np.argmin(ds)
                    disps.append(_disps[id])
                    orders.remove(ids[id])
                    match = True

            # print(match, site1, site2, trans, dist)
            if not match:
                return None, 5.0, False

        disps = np.array(disps)
        d = np.max(np.linalg.norm(disps.dot(cell1), axis=1))

        return disps, d, True

    def get_disps_optim(self, ref_struc, trans, d_tol):
        """
        Args:
            ref_struc: reference pyxtal structure (assuming the same atom order)
            trans: translation vector
            d_tol: tolerence of mismatch

        Returns:
            Atomic displacements in np.array
            translation:
        """

        from scipy.optimize import minimize

        def fun(tran, ref_struc, d_tol, axis):
            for i in range(3):
                if i not in axis:
                    tran[i] = 0
            disp, d, _ = self.get_disps_single(ref_struc, tran, d_tol)
            return d

        res = minimize(
            fun,
            trans,
            args=(ref_struc, d_tol, self.axis),
            method="Nelder-Mead",
            options={"maxiter": 10},
        )
        # print("Best_dist: {:6.3f} [{:6.3f} {:6.3f} {:6.3f}]".format(res.fun, *res.x))
        disp, d, _ = self.get_disps_single(ref_struc, res.x, d_tol)
        return disp, d, res.x

    def get_init_translations(self, ref_struc, tol=0.75):
        """
        Compute the displacement w.r.t. the reference structure

        Args:
            ref_struc: reference pyxtal structure (assuming the same atom order)

        Returns:
            list of possible translations
        """
        # print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        # print(ref_struc)

        axis = self.get_free_axis()

        translations = []
        # choose the one and avoid hydrogen is possible

        for specie in self.species:
            for site1 in self.atom_sites:
                if site1.specie == specie:
                    break
            for i in range(len(self.atom_sites)):
                site2 = ref_struc.atom_sites[i]
                if site1.specie == site2.specie and site1.wp.index == site2.wp.index:
                    trans0 = site1.get_translations(site2.position, axis)
                    translations.extend(trans0)

        # remove close translations
        good_translations = []
        for trans in translations:
            match = False
            for trans_ref in good_translations:
                diff = trans - trans_ref
                diff -= np.rint(diff)
                diff = np.dot(diff, self.lattice.matrix)
                if np.linalg.norm(diff) < tol:
                    match = True
                    break
            if not match:
                good_translations.append(trans)
        self.axis = axis
        return good_translations

    def is_duplicate(self, ref_strucs):
        """
        check if the structure is exactly the same
        """

        lat0 = self.lattice
        pos0 = self.atom_sites[0].position
        for ref_struc in ref_strucs:
            d_tol1, f_tol1, a_tol1, switch = ref_struc.lattice.get_diff(lat0)
            # print(d_tol1, f_tol1, a_tol1, switch); import sys; sys.exit()
            if (d_tol1 + a_tol1) < 1e-3 and not switch:
                tran = np.zeros(
                    [3]) if self.group.number > 15 else ref_struc.atom_sites[0].position - pos0
                disp, d, valid = self.get_disps_single(
                    ref_struc, -tran, d_tol=0.1)
                # print(self); print(ref_struc), print("=====", d); import sys; sys.exit()
                if d < 1e-3:
                    return True
        return False

    def get_disps_sets(self, ref_struc, d_tol, d_tol2=0.3, ld_tol=2.0, fd_tol=0.15, keep_lattice=False):
        """
        Compute the displacement w.r.t. a reference structure (considering all wycsets)

        Args:
            ref_struc: reference pyxtal structure (assuming the same atomic ordering)
            d_tol: maximally allowed atomic displacement
            d_tol2: displacement that allows early termination
            kepp_lattice: whether or not change the WP sets

        Returns:
            Atomic displacements in np.array
        """
        all_disps = []
        all_trans = []
        all_ds = []
        good_ref_strucs = []
        bad_ref_strucs = []
        # print(ld_tol, fd_tol)
        if keep_lattice:
            ref_strucs_matched = [ref_struc]
        else:
            ref_strucs_matched = self.find_matched_lattice(
                ref_struc, d_tol=ld_tol, f_tol=fd_tol)

        for _i, ref_struc_matched in enumerate(ref_strucs_matched):
            ref_strucs_alt = ref_struc_matched.get_alternatives(
                ref_lat=self.lattice, d_tol=ld_tol, f_tol=fd_tol)
            for _j, ref_struc_alt in enumerate(ref_strucs_alt):
                # initial setup
                d_min = 10.0
                disp_min = None
                trans = []

                # print('========================', ref_struc_alt)
                # print('=======', i, j, len(ref_strucs_matched), len(ref_strucs_alt))
                # must have the same wp letters and  different strucs
                if self.check_mapping(ref_struc_alt):
                    # if not ref_struc_alt.is_duplicate(good_ref_strucs+bad_ref_strucs):
                    trans = self.get_init_translations(ref_struc_alt)

                for _k, tran in enumerate(trans):
                    disp, d, valid = self.get_disps_single(
                        ref_struc_alt, tran, d_tol)
                    if valid and d > 0.3 and len(self.axis) > 0:
                        disp, d, tran = self.get_disps_optim(
                            ref_struc_alt, tran, d_tol)
                    # update
                    if d < d_min:
                        d_min = d
                        disp_min = disp
                        trans_min = tran

                    # strs = "\nlattice {:d} wyc {:d} trans {:d}".format(i, j, k)
                    # strs += "[{:6.3f} {:6.3f} {:6.3f}] {:6.3f}".format(*tran, d)
                    # print(strs)

                if d_min < d_tol2:  # Return it early
                    return disp_min, trans_min, ref_struc_alt, d_min
                elif d_min < d_tol:  # add to database
                    all_ds.append(d_min)
                    all_disps.append(disp_min)
                    all_trans.append(trans_min)
                    good_ref_strucs.append(ref_struc_alt)
                elif d_min < 5.1:  # add bad
                    bad_ref_strucs.append(ref_struc_alt)

        # choose the best
        # print("Good_candiates", len(good_ref_strucs), "Bad canidates", len(bad_ref_strucs))
        if len(all_ds) > 0:
            all_ds = np.array(all_ds)
            id = np.argmin(all_ds)
            if all_ds[id] < d_tol:
                return all_disps[id], all_trans[id], good_ref_strucs[id], all_ds[id]
            else:
                return None, None, None, None
        else:
            return None, None, None, None

        return None, None, None, None

    def _get_elements_and_sites(self):
        """
        Sometimes, the atoms are not arranged in order
        group the elements, sites

        Returns:
            elements: ['Si', 'O']
            sites: [['4b'], ['4a','4a']]
        """
        elements = []
        sites = []
        for at_site in self.atom_sites:
            e = at_site.specie
            site = at_site.wp.get_label()
            if e not in elements:
                elements.append(e)
                sites.append([site])
            else:
                id = elements.index(e)
                sites[id].append(site)
        return elements, sites

    def sort_sites_by_mult(self):
        mults = np.array([site.wp.multiplicity for site in self.atom_sites])
        seq = np.argsort(mults)
        self.atom_sites = [self.atom_sites[i] for i in seq]

    def sort_sites_by_numIons(self, seq=None):
        if seq is None:
            seq = np.argsort(self.numIons)

        self.atom_sites = [
            site for i in seq for site in self.atom_sites if self.species[i] == site.specie]

    def get_transition(self, ref_struc, d_tol=1.0, d_tol2=0.3, N_images=2, max_path=30, both=False):
        """
        Get the splitted wyckoff information along a given path:

        Args:
            ref_struc: structure with subgroup symmetry
            d_tol: maximally allowed atomic displacement
            d_tol2: displacement that allows early termination
            N_images: number of intermediate images
            max_path: maximum number of paths
            both: whether or not do interpolation along both sides

        Returns:
            - strucs:
            - displacements:
            - cell translation:
            - the list of space groups along the path
            - the list of splitters along the path
        """
        # ref_struc.sort_sites_by_numIons()
        # self.sort_sites_by_numIons()
        paths = self.group.search_subgroup_paths(ref_struc.group.number)
        if len(paths) == 0:
            print("No valid paths between the structure pairs")
            return None, None, None, None, None
        else:
            Skipped = len(paths) - max_path
            if Skipped > 0:
                paths = paths[:max_path]  # sample(paths, max_path)

            good_ds = []
            good_strucs = []
            good_disps = []
            good_paths = []
            good_trans = []
            good_splitters = []

            for p in paths:
                r = self.get_transition_by_path(
                    ref_struc, p, d_tol, d_tol2, N_images, both)
                (strucs, disp, tran, count, sps) = r
                if count == 0:
                    # prepare more paths to increase diversity
                    add_paths = self.group.add_k_transitions(p)
                    for p0 in add_paths:
                        r = self.get_transition_by_path(
                            ref_struc, p0, d_tol, d_tol2, N_images, both)
                        (strucs, disp, tran, count, sps) = r
                        if strucs is not None:
                            if strucs[-1].disp < d_tol2:  # stop
                                return strucs, disp, tran, p0, sps
                            else:
                                good_ds.append(strucs[-1].disp)
                                good_disps.append(disp)
                                good_paths.append(p0)
                                good_strucs.append(strucs)
                                good_trans.append(tran)
                                good_splitters.append(sps)
                else:
                    if strucs is not None:
                        if strucs[-1].disp < d_tol2:
                            return strucs, disp, tran, p, sps
                        else:
                            good_ds.append(strucs[-1].disp)
                            good_disps.append(disp)
                            good_paths.append(p)
                            good_strucs.append(strucs)
                            good_trans.append(tran)
                            good_splitters.append(sps)
                # Early stop
                if len(good_ds) > 5:
                    break
            if len(good_ds) > 0:
                # print("Number of candidate path:", len(good_ds))
                good_ds = np.array(good_ds)
                id = np.argmin(good_ds)
                return (
                    good_strucs[id],
                    good_disps[id],
                    good_trans[id],
                    good_paths[id],
                    good_splitters[id],
                )

            if Skipped > 0:
                print("Warning: ignore some solutions: ", Skipped)

            return None, None, None, p, None

    def get_transition_by_path(self, ref_struc, path, d_tol, d_tol2=0.5, N_images=2, both=False):
        """
        Get the splitted wyckoff information along a given path:

        Args:
            ref_struc: structure with subgroup symmetry
            path: a list of transition path
            d_tol: maximally allowed atomic displacement
            d_tol2: displacement that allows early termination
            N_images: number of intermediate images
            both: interpolation on both sides

        Returns:
            - strucs:
            - displacements:
            - cell translation:
        """
        import string

        # print("Searching the transition path.....", path)
        # Here we only check symbols
        count_match = 0
        elements0, sites_G = self._get_elements_and_sites()
        elements1, sites_H = ref_struc._get_elements_and_sites()
        # resort sites_H based on elements0
        seq = [elements1.index(x) for x in elements0]
        sites_H = [sites_H[i] for i in seq]
        numIons_H = [sum(int(l[:-1]) for l in site) for site in sites_H]

        # enumerate all possible solutions
        ids = []
        g_types = []
        G = self.group
        for p in path[1:]:
            dicts, g_type = G.get_max_subgroup(p)
            _ids = []
            for i, sub in enumerate(dicts["subgroup"]):
                tran = dicts["transformation"][i]
                if sub == p and np.linalg.det(tran[:3, :3]) <= 4:
                    _ids.append(i)
            ids.append(_ids)
            g_types.append(g_type)
            G = Group(p, quick=True)

        # print(ids)
        sols = list(itertools.product(*ids))

        # loop over all
        disps = []
        refs = []
        trans = []
        ds = []
        splitters = []

        for sol in sols:
            _sites = deepcopy(sites_G)
            G = self.group
            sol = list(sol)
            # mult = 1
            for p, s in zip(path[1:], sol):
                _sites0 = []
                dicts, _ = G.get_max_subgroup(p)
                relation = dicts["relations"][s]
                # tran = dicts['transformation'][s]; mult *= np.linalg.det(tran[:3,:3])
                # add site for each element
                for site in _sites:
                    _site = []
                    for label in site:
                        try:
                            index = string.ascii_lowercase.index(label[-1])
                        except ValueError:  # '8A'
                            index = 26
                        _site.extend(relation[index])
                    _sites0.append(_site)
                _sites = _sites0
                G = Group(p, quick=True)
            # print('============', path, mult)

            # match in sites and numbers
            match = True
            for i, site in enumerate(_sites):
                # sites
                if len(site) != len(sites_H[i]):
                    # print("bad sites", elements0[i], site, sites_H[i])
                    match = False
                    break
                # composition
                number = sum(int(l[:-1]) for l in site)
                if number != numIons_H[i]:
                    # print("bad number", site, number, numIons_H[i])
                    match = False
                    break
            # if int(mult) == 2: print(path, _sites0, match)
            # make subgroup
            if match:
                count_match += 1
                s, sps = self.subgroup_by_path(g_types, ids=sol, eps=0)
                if s is not None:
                    disp, tran, s, max_disp = ref_struc.get_disps_sets(
                        s, d_tol, d_tol2)
                    # import sys; sys.exit()
                    if disp is not None:
                        # early termination
                        if max_disp < d_tol2:
                            cell = s.lattice.matrix
                            strucs = ref_struc.make_transitions(
                                disp, cell, tran, N_images, both)
                            return strucs, disp, tran, count_match, sps
                        else:
                            disps.append(disp)
                            refs.append(s)
                            trans.append(tran)
                            ds.append(max_disp)
                            splitters.append(sps)

        if len(ds) > 0:
            ds = np.array(ds)
            id = np.argmin(ds)
            cell = refs[id].lattice.matrix
            tran = trans[id]
            disp = disps[id]
            sps = splitters[id]
            strucs = ref_struc.make_transitions(
                disp, cell, tran, N_images, both)
            return strucs, disp, tran, count_match, splitters

        else:
            return None, None, None, count_match, None

    def translate(self, trans, reset_wp=False):
        """
        move the atomic sites along a translation
        Note that this may change the structure

        Args:
            trans: 1*3 vector
        """
        for site in self.atom_sites:
            site.update(site.position + trans, reset_wp=reset_wp)

    def make_transitions(self, disps, lattice=None, translation=None, N_images=3, both=False):
        """
        make the pyxtals by following the atomic displacements

        Args:
            disps: N*3 atomic displacements in fractional coordinates
            lattice: 3*3 cell matrix (self.lattice.matrix if None)
            translation: overall translation
            N_images: number of images
            both: whether or not interpolar on both sides
        """
        N_images = max([N_images, 2])
        cell = self.lattice.matrix
        strucs = []

        disps = np.array(disps)
        disps /= N_images - 1
        max_disp = np.max(np.linalg.norm(disps.dot(cell), axis=1))
        l_disps = np.zeros([3, 3]) if lattice is None else (
            lattice - cell) / (N_images - 1)

        if translation is None:
            translation = np.zeros([3])

        i_list = range(2 * N_images - 1) if both else range(N_images)

        for i in i_list:
            struc = self.copy()
            for j, site in enumerate(struc.atom_sites):
                coord = site.position + i * disps[j] + translation
                struc.atom_sites[j].update(coord)
            struc.source = f"Transition {i:d} {max_disp * i:6.3f}"
            struc.disp = max_disp * i
            if i >= N_images:
                struc.lattice.set_matrix(
                    cell - (i - 2 * (N_images - 1)) * l_disps)
            else:
                struc.lattice.set_matrix(cell + i * l_disps)
            strucs.append(struc)
        return strucs

    def get_intermolecular_energy(self, factor=2.0, max_d=10.0):
        """
        For molecular crystals, get the intermolecular interactions from
        Gavezzotti, A., Filippini, G., J. Phys. Chem., 1994, 98 (18), 4831-4837

        Returns:
            Total energy
        """
        eng = 0
        for i in range(len(self.mol_sites)):
            res = self.get_neighboring_molecules(
                i, factor=factor, max_d=max_d, ignore_E=False)
            eng += np.array(res[-1]).sum()
        return eng

    def get_neighboring_molecules(self, site_id=0, factor=1.5, max_d=5.0, ignore_E=True):
        """
        For molecular crystals, get the neighboring molecules for a given WP

        Args:
            site_id: the index of reference site
            factor: factor of vdw tolerance
            max_d:
            ignore_E:

        Returns:
            min_ds: list of shortest distances
            neighs: list of neighboring molecular xyzs
            comps: list of molecular types
            Ps: list of [0, 1] to distinguish self and other molecules
            engs: list of energies from atom-atom potential
        """
        min_ds = []
        neighs = []
        comps = []
        Ps = []
        engs = []
        site0 = self.mol_sites[site_id]
        site0.get_ijk_lists()
        for id0, site1 in enumerate(self.mol_sites):
            if id0 == site_id:
                min_d0, neigh0, P, eng = site0.get_neighbors_auto(
                    factor, max_d, ignore_E)
            else:
                min_d0, neigh0, eng = site0.get_neighbors_wp2(
                    site1, factor, max_d, ignore_E)
                P = [1] * len(neigh0)
            comp = [site1.type] * len(min_d0)
            neighs.extend(neigh0)
            min_ds.extend(min_d0)
            Ps.extend(P)
            comps.extend(comp)
            engs.extend(eng)

        ids = np.argsort(min_ds) if engs[0] is None else np.argsort(engs)

        neighs = [neighs[i] for i in ids]
        comps = [comps[i] for i in ids]
        min_ds = [min_ds[i] for i in ids]
        Ps = [Ps[i] for i in ids]
        engs = [engs[i] for i in ids]
        return min_ds, neighs, comps, Ps, engs

    def get_spherical_images(self, **kwargs):
        """
        get the spherical image representation

        Args:
            model: either 'molecule' or 'contacts'

        Returns:
            the sph class
        """
        from pyxtal.descriptor import spherical_image

        return spherical_image(self, **kwargs)

    def get_neighboring_dists(self, site_id=0, factor=1.5, max_d=5.0):
        """
        For molecular crystals, get the neighboring molecules for a given WP

        Args:
            site_id: the index of reference site
            factor: factor of vdw tolerance
            max_d: maximum distances

        Returns:
            pairs: list of short contact pairs
            engs: list of energies from atom-atom potential
        """
        pairs = []
        engs = []
        dists = []

        site0 = self.mol_sites[site_id]
        site0.get_ijk_lists()
        for id0, site1 in enumerate(self.mol_sites):
            if id0 == site_id:
                _eng, _pair, _dist = site0.get_neighbors_auto(
                    factor, max_d, False, detail=True)
            else:
                _eng, _pair, _dist = site0.get_neighbors_wp2(
                    site1, factor, max_d, False, detail=True)
            pairs.extend([p[1] for p in _pair])
            engs.extend(_eng)
            dists.extend(_dist)
        return engs, pairs, dists

    def show_mol_cluster(
        self,
        id,
        factor=1.5,
        max_d=4.0,
        plot=True,
        ignore_E=False,
        cmap="YlGn",
        **kwargs,
    ):
        """
        display the local packing environment for a selected molecule
        """
        np.set_printoptions(precision=3)
        min_ds, neighs, comps, Ps, engs = self.get_neighboring_molecules(
            id, factor, max_d, ignore_E)
        print("Number of neighboring molecules", len(engs))
        print(np.array(engs))

        if plot:
            import matplotlib.pyplot as plt
            # from scipy.ndimage.filters import gaussian_filter1d

            plt.figure()
            x_min, x_max, size = 0, 2.5, 100
            res = (x_max - x_min) / size
            plt.gca()
            x = np.linspace(x_min, x_max, size)
            y = np.zeros([size, 2])
            for d, p in zip(min_ds, Ps):
                j = round(d / res)
                if p == 0:
                    y[j, 0] += 1
                else:
                    y[j, 1] += 1
            # y[:, 1] = gaussian_filter1d(y[:, 0], 0.2)
            plt.plot(x, y[:, 0], c="g", label=f"Self({int(sum(y[:, 0])):d})")
            if np.sum(y[:, 1]) > 1e-1:
                plt.plot(x, y[:, 1], c="c",
                         label=f"Other({int(sum(y[:, 1])):d})")

            cut = min_ds[-1]
            CN = len(min_ds)

            plt.axvline(x=cut, c="r", ls=":", label=f"Cutoff_{CN:d}")
            plt.legend()
            plt.xlim([x_min, x_max])
            plt.xlabel("R")
            plt.ylabel("Intensity")
            plt.show()

        site0 = self.mol_sites[id]
        coord0, specie0 = site0._get_coords_and_species(
            absolute=True, first=True)
        molecules = [Molecule(specie0, coord0)]

        for neigh, typ in zip(neighs, comps):
            specie0 = self.molecules[typ].mol.atomic_numbers
            molecules.append(Molecule(specie0, neigh))

        return display_cluster(molecules, self.lattice.matrix, engs, cmap, **kwargs)

    def substitute_1_2(
        self,
        dicts,
        ratio=None,
        group_type="t+k",
        max_cell=4,
        min_cell=0,
        max_wp=None,
        N_groups=None,
        N_max=10,
        criteria=None,
    ):
        """
        Derive the BC compounds from A via subgroup relation
        For example, from C to BN or from SiO2 to AlPO3.
        The key is to split A's wyckoff site to B and C by w.r.t. composition relation.
        If the current parent crystal doesn't allow splitting, look for the subgroup

        Args:
            dicts: e.g., {"F": "Cl"}
            ratio: e.g., [1, 1]
            group_type (string): `t`, `k` or `t+k`
            max_cell (float): maximum cell reconstruction
            min_cell (float): maximum cell reconstruction
            max_wp (int): maximum number of wp
            N_groups (int): maximum number of trial subgroups
            N_max (int): maximum number of structures
            criteria (dict): dictionary criteria

        Returns:
            A list of pyxtal structures
        """
        from pyxtal.util import new_struc_wo_energy

        if ratio is None:
            ratio = [1, 1]
        xtals = self._substitute_1_2(dicts, ratio)
        xtals = [
            xtal for xtal in xtals if criteria is None or xtal.check_validity(criteria)]

        if len(xtals) == 0:
            # print("\nCannot split, look for subgroup representation")
            subs = self.subgroup(
                group_type=group_type,
                max_cell=max_cell,
                min_cell=min_cell,
                N_groups=N_groups,
            )
            # print("Found {:d} subgroup representatios".format(len(subs)))
            for sub in subs:
                # print(sub)
                _xtals = sub._substitute_1_2(dicts, ratio)
                for _xtal in _xtals:
                    if max_wp is not None and len(_xtal.atom_sites) > max_wp:
                        continue
                    if new_struc_wo_energy(_xtal, xtals):
                        add = True

                        if criteria is not None and not _xtal.check_validity(criteria):
                            add = False

                        if add:
                            xtals.append(_xtal)
                            print("Add substitution", _xtal.get_xtal_string())
                            if len(xtals) == N_max:
                                break
                    # print('Add {:d} substitutions in subgroup {:d}'.format(len(_xtals), sub.group.number))
        else:
            print(
                f"Good representation ({len(xtals):d})", self.get_xtal_string())
        print(f"Found {len(xtals):d} substitutions in total")
        return xtals

    # , group_type='t', max_cell=4, min_cell=0):
    def _substitute_1_2(self, dicts, ratio=None):
        """
        Derive the BC compounds from A via subgroup relation
        For example, from C to BN or from SiO2 to AlPO3.
        The key is to split A's wyckoff site to B and C by w.r.t. composition relation.

        Args:
            dicts: e.g., {"F": "Cl"}
            ratio: e.g., [1, 1]
            group_type (string): `t`, `k` or `t+k`
            max_cell (float): maximum cell reconstruction
            min_cell (float): maximum cell reconstruction

        Returns:
            A list of pyxtal structures
        """
        from pyxtal.util import split_list_by_ratio

        if ratio is None:
            ratio = [1, 1]
        if len(dicts) > 1:
            raise ValueError("Only one item allowed in the input dicts", dicts)

        A, [B, C] = next(iter(dicts.items()))

        numbers = [
            site.wp.multiplicity for site in self.atom_sites if site.specie == A]
        solutions = split_list_by_ratio(numbers, ratio)

        # Output all possible substitutions
        xtals = []
        for group1, _group2 in solutions:
            xtal0 = self.copy()
            count = 0
            for site in xtal0.atom_sites:
                if site.specie == A:
                    if count in group1:
                        site.specie = B
                    else:
                        site.specie = C
                    count += 1
            xtal0.species = None
            xtal0._get_formula()  # update species
            xtals.append(xtal0)
        return xtals

        # idx, sites, t_types, k_types = self._get_subgroup_ids(H, group_type, None, max_cell, min_cell)
        # elements = [site.species for site in self.atom_sites]

        # valid_splitters = []
        # for id in idx:
        #    gtype = (t_types + k_types)[id]
        #    if gtype == 'k': id -= len(t_types)
        #    splitter = wyckoff_split(G=self.group, wp1=sites, idx=id, gtype, elements)
        #    # QZ: check if there is a redundancy in error and valid_split
        #    if not splitter.error and splitter.valid_split:
        #        for key in dicts.keys():
        #            # check whether the splitter can make valid decomposition
        #            # e.g., if we want to make a BN from C
        #            # The splitter from 8 to 4+4 is valid
        #            # The splitter from 8 to 6+2 is invalid
        #            numbers = []
        #            for i, wp2 in enumerate(splitter.wp2_lists):
        #                if sp.elements[i] == key:
        #                    for _wp2 in wp2:
        #                        numbers.append(_wp2.multiplicity)
        #            if not get_integer_solutions(numbers, ratio, quick=True):
        #                bad_splitter = True
        #                break
        #        if bad_splitter:
        #            continue
        #        else:
        #            valid_splitters.append(splitter)

        # new_strucs = []
        # for splitter in valid_splitters:
        #    #print(splitter)
        #    new_struc = self._apply_substitution(splitter, perms)
        #    new_strucs.append(new_struc)
        # return new_strucs

    def substitute(self, dicts):
        """
        A quick function to apply substitution

        Args:
            dicts: e.g., {"F": "Cl"}
        """
        pmg = self.to_pymatgen()
        pmg.replace_species(dicts)
        if self.molecular:
            for ele in dicts:
                smi = [m.smile.replace(ele, dicts[ele]) +
                       ".smi" for m in self.molecules]
            self.from_seed(pmg, smi)
        else:
            self.from_seed(pmg)

    def remove_species(self, species):
        """
        To remove the atom_sites by specie name
        Args:
            dicts: e.g., ["O"]
        """
        numIons = []
        new_species = []
        sites = []
        count = 0
        for site in self.atom_sites:
            sp = site.specie
            if sp not in species:
                num = site.wp.multiplicity
                sites.append(site)
                if sp in new_species:
                    id = new_species.index(sp)
                    numIons[id] += num
                else:
                    new_species.append(sp)
                    numIons.append(num)
            count += num
        struc = self.copy()
        struc.atom_sites = sites
        struc.numIons = numIons
        struc.species = new_species
        struc.resort()
        return struc

    def substitute_linear(self, dicts):
        """
        This is mainly designed for substitution between single atom and the linear block
        e.g., we want to make Zn(CN)2 from SiO2
        Args:
            dicts: e.g., {"Si": ["Zn"], "O": ["C","N"]}
        """
        from ase.neighborlist import neighbor_list

        if not hasattr(self, "cutoff"):
            self.set_cutoff()
            cutoff = self.cutoff

        atoms = self.to_ase(resort=False)
        (id1, id2, shifts) = neighbor_list("ijS", atoms, cutoff)

        self.lattice = self.lattice.scale(1.5)
        matrix = self.lattice.matrix
        numIons = []
        species = []
        sites = []
        count = 0
        for site in self.atom_sites:
            sp = site.specie
            if sp in dicts:
                if len(dicts[sp]) == 1:
                    e = dicts[sp][0]
                    sites.append(site.substitute_with_single(e))
                    if e in species:
                        id = species.index(e)
                        numIons[id] += site.wp.multiplicity
                    else:
                        species.append(e)
                        numIons.append(site.wp.multiplicity)
                else:
                    eles = dicts[sp]
                    ids = id2[id1 == count]  # index of ids
                    neighbors = atoms.positions[ids] + \
                        shifts[id1 == count].dot(atoms.cell)
                    direction = neighbors[1] - neighbors[0]
                    direction /= np.linalg.norm(direction)
                    s1, s2 = site.substitute_with_linear(
                        eles, direction, matrix)
                    for e in eles:
                        if e in species:
                            id = species.index(e)
                            numIons[id] += site.wp.multiplicity
                        else:
                            species.append(e)
                            numIons.append(site.wp.multiplicity)
                    sites.append(s1)
                    sites.append(s2)
            else:
                sites.append(site)
            count += site.wp.multiplicity

        struc = self.copy()
        struc.atom_sites = sites
        struc.numIons = numIons
        struc.species = species
        struc.resort()
        return struc

    def remove_water(self):
        """
        Remove water from hydrates
        """
        molecules = []
        numMols = []
        sites = []
        for i, m in enumerate(self.molecules):
            symbols = deepcopy(m.symbols)
            symbols.sort()
            if len(symbols) != 3 and symbols != ["H", "H", "O"]:
                molecules.append(m)
                numMols.append(self.numMols[i])

        for site in self.mol_sites:
            symbols1 = site.molecule.symbols
            for i, m in enumerate(molecules):
                symbols2 = m.symbols
                if len(symbols1) == len(symbols2) and symbols1 == symbols2:
                    site.type = i
                    sites.append(site)
                    # print("add sites", i, m)
                    break

        self.molecules = molecules
        self.numMols = numMols
        self.mol_sites = sites

    def set_cutoff(self, exclude_ii=False, value=None):
        """
        get the cutoff dictionary

        Args:
            exclude_ii (bool): whether or not exclude A-A distance
            value (float): cutoff distance
        """
        cutoff = {}
        tm = Tol_matrix(prototype="molecular")
        for i in range(len(self.species)):
            s1 = self.species[i]
            for j in range(i, len(self.species)):
                s2 = self.species[j]
                select = True
                if exclude_ii and s1 == s2:
                    select = False
                if select:
                    tuple_elements = (s1, s2)
                    if value is None:
                        cutoff[tuple_elements] = tm.get_tol(s1, s2)
                    else:
                        cutoff[tuple_elements] = value

        self.cutoff = cutoff

    def set_site_coordination(self, cutoff=None, verbose=False, exclude_ii=False):
        """
        Compute the coordination number from each atomic site
        """
        from ase.neighborlist import neighbor_list

        # if not hasattr(self, 'cutoff'):
        self.set_cutoff(exclude_ii, cutoff)
        my_cutoff = self.cutoff

        if verbose:
            print("\n The cutoff values for CN calculation are")
            print(my_cutoff)

        atoms = self.to_ase(resort=False)
        NL = neighbor_list("i", atoms, my_cutoff)
        coords = np.bincount(NL)

        count = 0
        for site in self.atom_sites:
            if count >= len(coords):
                site.coordination = 0
            else:
                site.coordination = coords[count]
            count += site.multiplicity

    def get_dimensionality(self, cutoff=None):
        """
        A quick wrapper to compute dimensionality from pymatgen
        https://pymatgen.org/pymatgen.analysis.dimensionality.html
        The dimensionality of the structure can be 1/2/3
        """
        from pymatgen.analysis.dimensionality import get_dimensionality_gorai

        if cutoff is None:
            if not hasattr(self, "cutoff"):
                self.set_cutoff()
            cutoff = self.cutoff

        return get_dimensionality_gorai(self.to_pymatgen(), bonds=cutoff)

    def from_CSD(self, csd_code):
        """
        Download the crystal from CCDC
        if csd_code is given, return the single pyxtal object
        if csd_family is given, do group analysis and ignore high pressure form

        Args:
            csd_code: e.g., ``ACSALA01``

        """
        from pymatgen.core.periodic_table import Element
        from pymatgen.io.cif import CifParser

        from pyxtal.msg import CSDError, ReadSeedError
        from pyxtal.util import get_struc_from__parser, process_csd_cif

        try:
            from ccdc import io
        except:
            msg = "No CSD-python api is available"
            raise CSDError(msg)
        try:
            from rdkit import Chem
        except:
            msg = "No rdkit is available"
            raise CSDError(msg)

        try:
            entry = io.EntryReader("CSD").entry(csd_code)
        except:
            msg = "Unknown CSD entry: " + csd_code
            raise CSDError(msg)

        if entry.has_3d_structure:
            smi = entry.molecule.smiles
            if smi is None:
                raise CSDError("No smile from CSD")

            if len(smi) > 350:
                raise CSDError(f"long smile {smi:s}")
            if Chem.MolFromSmiles(smi) is None:
                raise CSDError(f"problematic smiles: {smi:s}")

            cif = entry.to_string(format="cif")
            smiles = [s + ".smi" for s in smi.split(".")]

            # remove duplicates
            smiles = list(set(smiles))
            smi1 = ""
            for i, s in enumerate(smiles):
                smi1 += s[:-4]
                if i + 1 < len(smiles):
                    smi1 += "."

            self.tag = {
                "smiles": smi1,
                "csd_code": csd_code,
                "ccdc_number": entry.ccdc_number,
                # 'publication': entry.publication,
            }

            cif = process_csd_cif(cif)  # , remove_H=True)
            # print(cif)
            try:
                parser = CifParser.from_str(cif, occupancy_tolerance=2.0)
                pmg = get_struc_from__parser(parser)
                # pmg = Structure.from_str(cif, fmt='cif')
            except:
                print(cif)
                msg = "Problem in parsing CSD cif"
                raise CSDError(msg)

            organic = True
            for ele in pmg.composition.elements:
                if ele.symbol == "D":
                    pmg.replace_species({ele: Element("H")})
                elif ele.value not in "C Si H O N S F Cl Br I P".split(" "):
                    organic = False
                    break

            if not organic:
                msg = "Cannot handle the organometallic entry from CSD: "
                msg += entry.formula
                raise CSDError(msg)
            try:
                self.from_seed(pmg, smiles)
            except ReadSeedError:
                try:
                    # print("Add_H=============================================")
                    self.from_seed(pmg, smiles, add_H=True)
                except:
                    msg = f"unknown problems in Reading CSD {csd_code:s} {smi:s}"
                    raise CSDError(msg)
            except:
                msg = f"unknown problems in Reading CSD {csd_code:s} {smi:s}"
                raise CSDError(msg)
            self.source = "CSD: " + csd_code
        else:
            msg = csd_code + " does not have 3D structure"
            raise CSDError(msg)

        # check if the dumped cif is correct
        cif0 = self.to_file()
        try:
            Structure.from_str(cif0, fmt="cif")
        except:
            print(cif0)
            raise CSDError("Cif Error")
            # import sys; sys.exit()
        # ====================

        # check if the structure is consistent with the origin
        # pmg0 = self.to_pymatgen() #shape='lower')
        # if remove_H:
        #    print("REMOVE H")
        #    pmg.remove_species('H')
        #    pmg0.remove_species('H')
        # import pymatgen.analysis.structure_matcher as sm
        # if not sm.StructureMatcher().fit(pmg0, pmg):
        #    pmg = Structure.from_str(cif, fmt='cif')
        #    #pmg0 = Structure.from_str(self.to_file(), fmt='cif')
        #    pmg0 = Structure.from_str(pmg0.to(fmt='cif'), fmt='cif')
        #    print(cif)
        #    print(self)
        #    print(pmg0)
        #    pmg0.remove_species('H'); pmg0.remove_species('O')
        #    pmg.remove_species('H'); pmg.remove_species('O')
        #    #print(pmg0.to(fmt='cif'))
        #    print(sm.StructureMatcher().fit(pmg0, pmg))
        #    print(pmg) #reference
        #    print(pmg0) #pyxtal
        #    #print(cif)
        #    #print(self.to_file())
        #    print("Wrong", csd_code); import sys; sys.exit()


    def get_structure_factor(self, hkl, coeffs=None):
        """
        Compute the structure factor for a given hkl plane

        Args:
            - hkl: three indices hkl plane
            - coeffs (N lengths): Atomic scatering factors
        """
        coords, species = self._get_coords_and_species()
        g_dot_rs = np.dot(coords, np.array(hkl))
        F = np.exp(-2 * np.pi * (1j) * g_dot_rs)
        if coeffs is not None:
            F *= coeffs
        return F.sum()

    def to_molecular_xtal(self, molecules, oris=None, reflects=None):
        """
        Convert the atomic xtal to molecular xtal
        the input molecules must have the same length of the self.atom_sites
        """
        if not self.molecular:
            xtal = pyxtal(molecular=True)
            # updates mol_sites
            sites = []
            for i, site in enumerate(self.atom_sites):
                ori = [0, 0, 0] if oris is None else oris[i]
                reflect = False if reflects is None else reflects[i]
                sites.append(site.to_mol_site(
                    self.lattice, molecules[i], ori=ori, reflect=reflect))
            xtal.lattice = self.lattice
            xtal.group = self.group
            xtal.mol_sites = sites
            xtal.numMols = [sum([wp.multiplicity for wp in self.atom_sites])]
            xtal._get_formula()
            xtal.valid = True
            xtal.source = "Center"
            return xtal
        else:
            raise RuntimeError("The input must be molecular xtal")

    def to_atomic_xtal(self):
        """
        Convert the molecular
        """
        if self.molecular:
            xtal = pyxtal()
            sites = []
            for i, site in enumerate(self.mol_sites):
                sites.append(site.to_atom_site(i + 1))
            xtal.lattice = self.lattice
            xtal.group = self.group
            xtal.atom_sites = sites
            xtal._get_formula()
            xtal.valid = True
            xtal.source = "Mol. Center"
            return xtal
        else:
            raise RuntimeError("The input must be molecular xtal")

    def get_forcefield(self, ff_style="openff", code="lammps", chargemethod="am1bcc", parameters=None):
        """
        An interface to create forcefield for molecular simulation with
        - Charmm
        - LAMMPS
        Note that the current approach generates charge with its own conformer,
        need to provide the option to use the conformer from the given molecule

        Args:
            - ff_style: 'gaff' or 'openff'
            - code: 'lammps' or 'charmm'
            - charge_method: 'am1bcc', 'am1-mulliken', 'mmff94', 'gasteiger'
            - parameters: 1D-array of user-specified FF parameters

        Returns:
            An ase_atoms_objects with force field information
        """
        from pyocse.charmm import CHARMMStructure
        from pyocse.forcefield import forcefield
        from pyocse.interfaces.parmed import ParmEdStructure
        from pyocse.lmp import LAMMPSStructure

        if self.molecular:
            smiles = [mol.smile for mol in self.molecules]
            assert smiles[0] is not None

            # Lammps needs the xyz for all molecules
            if code == "lammps":
                # reset lammps cell
                from pyxtal.util import reset_lammps_cell

                atoms = self.to_ase(resort=False)
                atoms = reset_lammps_cell(atoms)
                converter = LAMMPSStructure
                n_mols = self.numMols

            # Charmm only needs the xyz for the asymmetric unit
            else:
                converter = CHARMMStructure
                coords, species = [], []
                checked = []
                for site in self.mol_sites:
                    if site.type not in checked:
                        _coords, _species = site._get_coords_and_species(
                            absolute=True, first=True)
                        coords.extend(_coords)
                        species.extend(_species)
                        checked.append(site.type)
                n_mols = [1] * len(self.molecules)
                atoms = Atoms(symbols=species, positions=coords)

            # print(len(atoms), n_mols)
            # Initialize the forcefield instance from pyocse
            ff = forcefield(smiles, style=ff_style, chargemethod=chargemethod)

            # update the parameters from the user
            if parameters is not None:
                ff.update_parameters(parameters)

            # Create the Parmed to handle FF parameters
            if sum(n_mols) == 1:
                pd_struc = ff.molecules[0].copy(cls=ParmEdStructure)
            else:
                from functools import reduce
                from operator import add

                mols = []
                for i, m in enumerate(n_mols):
                    mols += [ff.molecules[i] * m]
                pd_struc = reduce(add, mols)
            pd_struc.update(atoms)

            ase_with_ff = converter.from_structure(pd_struc)
            if code == "lammps":
                print("\n Write lmp.in and lmp.dat ......")
                ase_with_ff.write_lammps()
            else:
                print("\n Write charmm.rtf and charmm.prm ......")
                ase_with_ff.write_charmmfiles()
                # ase_with_ff.write_charmmfiles_by_parmd()

            return ase_with_ff
        else:
            raise NotImplementedError("\nNo support for atomic crystals")

    def update_from_1d_rep(self, x):
        """
        update the xtal from the 1d representation
        Group: I 41/a m d:2 (141)
        2.5019,   2.5019,   8.7534,  90.0000,  90.0000,  90.0000, tetragonal
        Wyckoff sites:
        C @ [ 0.0000  0.2500  0.4614], WP [8e] Site [2mm.]

        the above rep is [2.5019, 8.7514, 0.4614]

        Args:
            x (float): input variable array to describe a xtal
        """
        N = self.lattice.dof
        cell, pos = x[:N], x[N:]

        try:
            # update cell
            self.lattice.update_from_1d_representation(cell)
        except:
            self.lattice = None
            return
        #l_type = self.lattice.ltype
        #self.lattice = Lattice.from_1d_representation(cell, l_type)
        # print("lattice dof", N, cell, l_type, self.lattice)

        # update position
        start = 0
        for site in self.atom_sites:
            end = start + site.wp.get_dof()
            xyz = site.wp.get_position_from_free_xyzs(pos[start:end])
            site.update(xyz)
            start = end
        # xtal.to_ase().write('init.vasp', vasp5=True, format='vasp', direct=True)

    def get_1d_rep_x(self):
        """
        Get a simplified x representation for a crystal
        Group: I 41/a m d:2 (141)
        2.5019,   2.5019,   8.7534,  90.0000,  90.0000,  90.0000, tetragonal
        Wyckoff sites:
        C @ [ 0.0000  0.2500  0.4614], WP [8e] Site [2mm.]

        The above rep is [2.5019, 8.7514, 0.4614]
        """

        rep = self.get_1D_representation(standard=True)
        cell, xyzs = rep.x[0][1:], rep.x[1:]
        x = cell
        for xyz in xyzs:
            if len(xyz) > 2:
                x = np.hstack((x, xyz[2:]))
        return x

    def from_spg_wps_rep(self, spg, wps, x, elements=None):
        """
        An advanced way to build pyxtal from wyckoff position
        and compacted 1d x representation

        Args:
            spg (int): space group number 1-230
            wps (list): wyckoff representation (e.g., ['6a', '6a'])
            x (array): 1d representation
            elements (list): list of element string
        """
        if elements is None:
            elements = ["C"] * len(wps)
        group = Group(spg)
        sites = []
        for i, _wp in enumerate(wps):
            if type(_wp) is int:
                sites.append((elements[i], group[_wp]))
            else:
                letter = _wp[-1]
                for wp in group:
                    if wp.letter == letter:
                        sites.append((elements[i], wp))
                        break
        self.from_1d_rep(x, sites)

    def from_1d_rep(self, x, sites, dim=3):
        """
        An advanced way to build pyxtal from the 1d representation

        Args:
            x: 1d representation
            sites: list of (element, wp)
        """
        l_type = sites[0][1].lattice_type
        N = sum(Lattice.get_dofs(l_type))
        cell, pos = x[:N], x[N:]
        lat = Lattice.from_1d_representation(cell, l_type)

        species, numIons = [], []
        total_sites = []

        start = 0
        for site in sites:
            (specie, _wp) = site
            end = start + _wp.get_dof()
            xyz = pos[start:end]
            full_xyz = _wp.get_position_from_free_xyzs(xyz)
            [x, y, z] = full_xyz
            start = end
            label = _wp.get_label()
            if specie not in species:
                species.append(specie)
                numIons.append(_wp.multiplicity)
                total_sites.append([(label, x, y, z)])
            else:
                idx = species.index(specie)
                numIons[idx] += _wp.multiplicity
                total_sites[idx].append((label, x, y, z))
        # print(species, numIons,  total_sites)

        spg = sites[0][1].number
        try:
            self.build(spg, species, numIons, lat, total_sites, dim=dim)
        except:
            print("input x", x)
            print("input build", spg, numIons, lat, total_sites)
            raise ValueError("Problem in build")

    def check_validity(self, criteria, verbose=False):
        """
        Check if the xtal is valid for a given list of criteria.
        Currently support the following keywords
        - MIN_Density (float)
        - CN: coordination number (int)
        - Dimension: (int)

        Args:
            criteria (dict): keywords and the threshhold values
            verbose (bool): whether or not print out details

        Return:
            bool value regarding the validity
        """

        if len(criteria.keys()) == 0:
            return True

        if "exclude_ii" not in criteria:
            criteria["exclude_ii"] = False

        if "MIN_Density" in criteria:
            den1 = self.get_density()
            den2 = criteria["MIN_Density"]
            if den1 < den2:
                if verbose:
                    msg = f"===Invalid in MIN_density {den1:.2f}=>{den2:.2f}"
                    print(msg)
                return False
        if "CN" in criteria:
            options = [True, False] if criteria["exclude_ii"] else [False]

            coords = []
            for option in options:
                try:
                    self.set_site_coordination(
                        criteria["cutoff"], exclude_ii=option)
                    # print('exclude_ii', option, criteria['cutoff'], self)
                except:
                    if verbose:
                        print("=====Invalid in CN calculation")
                    return False
                coord = []
                for s in self.atom_sites:
                    ele = s.specie
                    cn1 = s.coordination
                    cn2 = criteria["CN"][ele]
                    # print(ele, cn1, option)
                    if cn1 not in cn2:
                        if verbose:
                            strs = f"=====Invalid CN {ele:s} [{cn1:d}=>{cn2[0]:d}]"
                            strs += ", exclude ii: " + str(option)
                            print(strs)
                        return False
                    coord.append(cn1)
                coords.append(coord)
            if len(coords) == 2 and coords[0] != coords[1]:
                return False

        if "Dimension" in criteria:
            try:
                # TODO: unclear what is the criteria_cutoff
                dim1 = self.get_dimensionality(criteria['cutoff'])
            except:
                dim1 = 3
            dim2 = criteria["Dimension"]
            if dim1 != dim2:
                if verbose:
                    print(f"=====Invalid in Dimension {dim1:d}=>{dim2:d}")
                return False

        return True

    def get_xtal_string(self, dicts=None, header=None):
        """
        get the xtal string

        Args:
            xtal: pyxtal instance
            sim (float): similarity metric
            status (bool): whether or not the structure is valid
            energy (float): the energy value
        """
        spg_num, spg_symbol = self.group.number, self.group.symbol
        density = self.get_density()
        dof = self.get_dof()
        N_atoms = sum(self.numIons)
        strs = header if header is not None else ""
        strs += f"*{N_atoms:3d} {dof:3d} {spg_num:4d} {spg_symbol:10s} {density:8.2f} "

        if dicts is not None:
            # Similarity, energy, status
            for key in dicts:
                value = dicts[key]
                if isinstance(value, (float, np.float64)):
                    strs += f" {value:8.3f} "
                elif isinstance(value, str):
                    strs += f" {value:24s} "
                elif isinstance(value, bool):
                    strs += f" {value!s:5s} "

        for s in self.atom_sites:
            strs += f"{s.wp.get_label():s} "
        return strs

    def get_tabular_representations(
        self,
        N_max=200,
        N_wp=8,
        normalize=False,
        perturb=False,
        eps=0.05,
        discrete=False,
        discrete_cell=False,
        N_grids=50,
    ):
        """
        Enumerate the equivalent representations for a give crystal.
        For example, 8a in space group 227 has 8 equivalent positions

        Args:
            N_max (int): the max number of samples from the enumeration
            N_wp (int): the maximum number of allowed WP sites
            normalize (bool): normalize all number to (0, 1)
            perturb (bool): whether or not perturb the variables
            eps (float): the degree of perturbation
            discrete (bool): to discretize xyz
            N_grids (int): number of grids used for discretization

        Returns:
            a list of equivalent 1D tabular representations
        """
        reps = []

        # To prevent the explosion of big multiplicity number
        if len(self.atom_sites) <= 1:
            min_wp = 192
        elif len(self.atom_sites) <= 2:
            min_wp = 100
        elif len(self.atom_sites) <= 4:
            min_wp = 20
        else:
            min_wp = int(np.power(100000.0, 1 / len(self.atom_sites)))

        sites_mul = [range(min([min_wp, site.wp.multiplicity]))
                     for site in self.atom_sites]
        ids = list(itertools.product(*sites_mul))
        if len(ids) > N_max:
            ids = self.random_state.choice(ids, N_max)

        #print(f"N_reps {len(ids)}/{N_max}: ", self.get_xtal_string())
        for sites_id in ids:
            rep = self.get_tabular_representation(
                    sites_id, normalize, N_wp, perturb,
                    eps=eps,
                    standardize=False,
                    discrete=discrete,
                    discrete_cell=discrete_cell,
                    N_grids=N_grids)
            reps.append(rep)
        return reps

    def get_tabular_representation(
        self,
        ids=None,
        normalize=False,
        N_wp=6,
        perturb=False,
        max_abc=50.0,
        max_angle=np.pi,
        eps=0.05,
        standardize=True,
        discrete=False,
        discrete_cell=False,
        N_grids=50,
    ):
        """
        Convert the xtal to a 1D reprsentation organized as
        (space_group_number, a, b, c, alpha, beta, gamma, wp1, wp2, ....)
        each wp is represented by 4 numbers (wp_id, x, y, z)

        Args:
            ids (list): index of generators for each wp
            normalize (bool): whether or not normalize the data to (0, 1)
            N_wp (int): number of maximum wp sites
            perturb (bool): whether or not perturb the variables
            max_abc (float): maximum abc length for normalization
            max_angle (float): maximum angle in radian for normalization
            eps (float): the degree of perturbation
            discrete (bool): to discretize xyz
            discrete_cell (bool): to discretize cell_para
            N_grids (int): number of grids used for discretization

        Returns:
            a 1D vector such as (141, 3, 3, 3, pi/2, pi/2, pi/2, wp_id, x, y, z)
        """

        if ids is None:
            ids = [0] * len(self.atom_sites)
        else:
            assert len(ids) == len(self.atom_sites)

        if len(ids) > N_wp:
            raise ValueError(
                "Cannot handle too many atom sites", len(ids), N_wp)

        N_cell, N_site = 7, N_wp * 4

        # If the number of wp is less then 6, assign -1
        rep = -np.ones(N_cell + N_site)  # 7 + 6 + 18 = 31
        rep[0] = self.group.number
        lattice = self.lattice.copy()
        if perturb:
            lattice = lattice.mutate(degree=eps)
        rep[1:7] = lattice.get_para()

        if normalize:
            rep[0] /= 230
            rep[1:4] /= max_abc
            rep[4:7] /= max_angle

        if discrete_cell:
            rep[1:4] = np.round(rep[1:4] / max_abc * N_grids)
            rep[4:7] = np.round(rep[4:7] / max_angle * N_grids)

        count = N_cell

        # Convert WP site
        for id, site in zip(ids, self.atom_sites):
            # Get the index
            rep[count] = site.wp.index
            if normalize:
                rep[count] /= len(self.group)

            # Convert xyz
            xyz = site.coords[id]  # position
            xyz -= np.floor(xyz)

            if discrete:
                xyz = site.wp.to_discrete_grid(xyz, N_grids)
            else:
                if standardize:
                    free_xyzs = site.wp.get_free_xyzs(
                        xyz, perturb=perturb, eps=eps)
                    xyz = site.wp.get_position_from_free_xyzs(free_xyzs)
            rep[count + 1: count + 4] = xyz
            count += 4
        return rep

    def from_tabular_representation(
        self,
        rep,
        max_abc=50.0,
        max_angle=180,
        normalize=False,
        tol=0.1,
        discrete=False,
        discrete_cell=False,
        N_grids=50,
        verbose=False,
    ):
        """
        Reconstruc xtal from 1d tabular_representation
        Currently assuming the elemental composition like carbon

        Args:
            rep: 1D array
            max_abc (float): maximum abc length for normalization
            max_angle (float): maximum angle in radian for normalization
            normalize (bool): whether normalize or not?
            tol (float): a tolerance value to assign the Wyckoff site
            discrete (bool): to discretize xyz
            N_grids (int): number of grids used for discretization
            verbose (bool): output detailed error
        """
        number = int(np.round(rep[0] * 230)) if normalize else int(rep[0])
        if 1 <= number <= 230:
            group = Group(number)
            [a, b, c, alpha, beta, gamma] = rep[1:7]
            if normalize:
                a, b, c = a * max_abc, b * max_abc, c * max_abc
                alpha, beta, gamma = (
                    alpha * max_angle,
                    beta * max_angle,
                    gamma * max_angle,
                )
            else:
                if discrete_cell:
                    a = a / N_grids *  max_abc
                    b = b / N_grids *  max_abc
                    c = c / N_grids *  max_abc
                    alpha = alpha / N_grids * max_angle
                    beta = beta / N_grids * max_angle
                    gamma = gamma / N_grids * max_angle
                    #print(a, b, c, alpha, beta, gamma)
                else:
                    alpha, beta, gamma = (
                        np.degrees(alpha),
                        np.degrees(beta),
                        np.degrees(gamma),
                    )
            try:
                lattice = Lattice.from_para(
                    a,
                    b,
                    c,
                    alpha,
                    beta,
                    gamma,
                    ltype=group.lattice_type,
                    force_symmetry=True,
                )
            except ValueError:
                print('Input lattice from rep is incorrect',
                      number, a, b, c, alpha, beta, gamma)
                self.valid = False
                return None

            # Convert WP site
            sites_info = np.reshape(rep[7:], (int((len(rep) - 7) / 4), 4))
            sites = []
            numIons = 0
            for site_info in sites_info:
                [id, x, y, z] = site_info
                # exclude data containing negative values
                if verbose:
                    print("site_info", site_info)

                if min(site_info) > -0.01:
                    wp_id = int(np.round(len(group) * id)
                                ) if normalize else int(id)
                    if wp_id >= len(group):
                        wp_id = -1
                    wp = group[wp_id]

                    # Conversion from discrete to continuous
                    if discrete:
                        #print('discrete', x, y, z)
                        [x, y, z] = wp.from_discrete_grid([x, y, z], N_grids)
                        #print('conversion', x, y, z)

                    # ; print(wp.get_label(), xyz)
                    xyz = wp.search_generator([x, y, z], tol=tol, symmetrize=True)
                    if xyz is not None:
                        #xyz, wp, _ = wp.merge(xyz, np.eye(3), 1e-3)
                        label = wp.get_label()
                        #print('after merge', x, y, z, label, xyz[0], xyz[1], xyz[2])
                        sites.append((label, xyz[0], xyz[1], xyz[2]))
                        numIons += wp.multiplicity
                        if verbose:
                            print('add wp', label, xyz)
                    else:
                        if verbose:
                            print("Cannot get wp in", x, y, z, wp.get_label())

            if len(sites) > 0:
                try:
                    # print(sites)
                    self.build(group, ["C"], [numIons], lattice, [sites])
                except:
                    print("Invalid Build", number, lattice, numIons, sites)
                    self.valid = False
            else:
                if verbose:
                    print("Empty sites in tabular_representation", rep)
                    print("parsed sites info", sites_info)
                self.valid = False
        else:
            if verbose:
                print("The input space group is invalid", rep[0])

    def get_Pearson_Symbol(self):
        """
        Based on https://en.wikipedia.org/wiki/Pearson_symbol
        """
        if self.group.number <= 2:
            l = "a"
        elif self.group.number <= 15:
            l = "m"
        elif self.group.number <= 74:
            l = "o"
        elif self.group.number <= 142:
            l = "t"
        elif self.group.number <= 194:
            l = "h"
        else:
            l = "c"

        b = "S" if self.group.symbol[0] in [
            "A", "B", "C"] else self.group.symbol[0]

        return l + b + str(sum(self.numIons))

    def resymmetrize(self, tol=1e-3):
        """
        Recheck the symmetry with a given tolerance value
        Useful to identify a higher space group symmetry

        Args:
            tol (float): the tolerance value for symmetry finder
        """

        xtal0 = pyxtal()
        xtal0.from_seed(self.to_pymatgen(), tol=tol)

        return xtal0

    def from_prototype(self, prototype):
        """
        Load a template structure based on a prototype name for quick testing.

        Args:
            prototype (str): Name of the crystal structure prototype to load.
                Supported ('diamond', 'graphite': 'a-cristobalite', .etc)

        Example:
            self.from_prototype('diamond')
            self.from_prototype('graphite')

        Notes:
            This method acts as a shortcut to quickly generate known crystals
            by loading the space group, Wyckoff positions, lattice parameters,
            and atomic species.
        """

        if prototype == 'graphite':
            self.from_spg_wps_rep(194, ['2c', '2b'], [2.46, 6.70])

        elif prototype == 'diamond':
            self.from_spg_wps_rep(227, ['8a'], [3.6])

        elif prototype == 'a-cristobalite':
            self.from_spg_wps_rep(92,
                                  ['4a', '8b'],
                                  [5.085, 7.099, 0.294, 0.094, 0.241, 0.826],
                                  ['Si', 'O'])
        elif prototype == 'b-cristobalite':
            self.from_spg_wps_rep(227,
                                  ['8a', '16c'],
                                  [7.46086],
                                  ['Si', 'O'])
        elif prototype == 'a-quartz':
            self.from_spg_wps_rep(152,
                                  ['3a', '6c'],
                                  [4.91, 5.43, 0.4689, 0.2692, 0.4134, 0.7849],
                                  ['Si', 'O'])
        elif prototype == 'b-quartz':
            self.from_spg_wps_rep(181,
                                  ['3c', '6j'],
                                  [5.06, 5.54, 0.2086],
                                  ['Si', 'O'])
        elif prototype in ['rocksalt', 'B1']:
            self.from_spg_wps_rep(225,
                                  ['4a', '4b'],
                                  [5.59],
                                  ['Na', 'Cl'])
        elif prototype in ['B2']:
            self.from_spg_wps_rep(221,
                                  ['1a', '1b'],
                                  [5.59],
                                  ['Cs', 'Cl'])
        else:
            raise ValueError("Cannot support the input prototype", prototype)

    def optimize_orientation_by_energy(self, max_iter=20, verbose=False):

        for i, site in enumerate(self.mol_sites):
            ms = [self.mol_sites[j] for j in range(len(self.mol_sites)) if j != i]
            site.optimize_orientation_by_energy(wps=ms, verbose=verbose)

    def get_orientation_energy(self):
        eng = 0
        for i, site in enumerate(self.mol_sites):
            ms = [self.mol_sites[j] for j in range(i+1, len(self.mol_sites))]
            eng += site.wp.multiplicity * site.get_energy(wps=ms)
        return eng

    def get_separations(self, hkls=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
        """
        Compute the separation for the given hkl plane
        """
        separations = []
        if self.molecular:
            from pyxtal.plane import planes
            p = planes()
            p.set_xtal(self)
            for hkl in hkls:
                (_, _, sep) = p.get_separation(hkl)
                separations.append(sep.max())

        return np.array(separations)

    def cut_lattice(self, max_separation=3.0, verbose=False):
        """
        An utility to reduce the empty spacing
        """
        seps = self.get_separations()
        ax = seps.argmax()
        if seps[ax] >= max_separation:
            cut = seps[ax] - max_separation
            # update coordinates
            for mol_site in self.mol_sites:
                mol_site.cut_lattice(ax, cut)
            self.lattice.update_para(ax, -cut)
            if verbose:
                print(f"Found large separation {ax} {seps[ax]:.2f}")
                print("Update lattice", self.lattice)

    def optimize_lattice_and_rotation(self, iterations=3, verbose=False):
        """
        Iteratively cut lattice and optimize the rotation
        """
        if self.molecular:
            for i in range(iterations):
                self.cut_lattice(2.0, verbose=verbose)
                if self.molecules[0].active_sites is not None:
                    self.optimize_orientation_by_energy(20, verbose=verbose)
