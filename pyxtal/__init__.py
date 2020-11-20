# Standard Libraries
from copy import deepcopy
from random import choice
import numpy as np

# PyXtal imports #avoid *
from pyxtal.version import __version__
from pyxtal.crystal import (
    random_cluster,
    random_crystal,
    random_crystal_1D,
    random_crystal_2D,
)
from pyxtal.symmetry import Group, Wyckoff_position
from pyxtal.wyckoff_site import atom_site
from pyxtal.wyckoff_split import wyckoff_split
from pyxtal.lattice import Lattice
from pyxtal.operations import apply_ops
from pyxtal.tolerance import Tol_matrix

name = "pyxtal"

def print_logo():
    """
    print out the logo and version information
    """

    print(
        """
             ______       _    _          _   
            (_____ \     \ \  / /        | |   
             _____) )   _ \ \/ / |_  ____| |  
            |  ____/ | | | )  (|  _)/ _  | | 
            | |    | |_| |/ /\ \ |_( (_| | |___
            |_|     \__  /_/  \_\___)__|_|_____)
                   (____/      """
    )
    print("\n")
    print("----------------------(version", __version__, ")----------------------\n")
    print("A Python package for random crystal generation")
    print("The source code is available at https://github.com/qzhu2017/pyxtal")
    print("Developed by Zhu's group at University of Nevada Las Vegas\n\n")

class pyxtal:
    """
    Class for handling atomic crystals based on symmetry constraints. 
    """

    def __init__(self):
        self.valid = False

    def __str__(self):
        if self.valid:
            s = "------Crystal from {:s}------".format(self.source)
            s += "\nComposition: {}".format(self.formula)
            s += "\nDimension: {}".format(self.dim)
            s += "\nGroup: {} ({})".format(self.group.symbol, self.group.number)
            s += "\n{}".format(self.lattice)
            s += "\nWyckoff sites:"
            for wyc in self.atom_sites:
                s += "\n\t{}".format(wyc)
        else:
            s += "\nStructure not available."
        return s

    def __repr__(self):
        return str(self)

    def from_random(
        self,
        dim = 3,
        group=None,
        species=None,
        numIons=None,
        factor=1.1,
        thickness = None,
        area = None,
        lattice=None,
        sites = None,
        tm=Tol_matrix(prototype="atomic"),
    ):
        count = 0
        while True:
            count += 1
            if dim == 3:
                struc = random_crystal(group, species, numIons, factor, lattice, sites, tm)
            elif dim == 2:
                struc = random_crystal_2D(group, species, numIons, factor, thickness, lattice, sites, tm)
            elif dim == 1:
                struc = random_crystal_1D(group, species, numIons, factor, area, lattice, sites, tm)
            else:
                struc = random_cluster(group, species, numIons, factor, lattice, sites, tm)

            if struc.valid:
                self.valid = True
                self.dim = dim
                self.lattice = struc.lattice
                self.numIons = struc.numIons
                self.group = struc.group
                self.atom_sites = struc.atom_sites
                self.PBC = struc.PBC
                self.source = 'random'
                self.formula = struc.formula
                self.factor = struc.factor
                self.number = struc.number
                break
            if count >= 10:
                raise RuntimeError("It takes long time to generate the structure, check inputs")

    def from_seed(self, seed):
        """
        Load the seed structure from Pymatgen/ASE/POSCAR/CIFs
        Internally they will be handled by Pymatgen
        """
        from ase import Atoms
        from pymatgen import Structure

        if isinstance(seed, dict):
            self.from_dict()
        elif isinstance(seed, Atoms): #ASE atoms
            from pymatgen.io.ase import AseAtomsAdaptor
            pmg_struc = AseAtomsAdaptor.get_structure(seed)
            self.from_pymatgen(pmg_struc)
        elif isinstance(seed, Structure): #Pymatgen
            self.from_pymatgen(seed)
        elif isinstance(seed, str):
            pmg_struc = Structure.from_file(seed)
            self.from_pymatgen(pmg_struc)

        formula = ""
        for i, s in zip(self.numIons, self.species):
            formula += "{:s}{:d}".format(s, int(i))
        self.formula = formula
        self.factor = 1.0
        self.number = self.group.number
        self.source = 'Seed'
        self.dim = 3
        self.PBC = [1, 1, 1]

    def from_pymatgen(self, structure):
        """
        Load the seed structure from Pymatgen/ASE/POSCAR/CIFs
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
        self.valid = True
        try:
            # needs to do it twice in order to get the conventional cell
            s = sga(structure)
            structure = s.get_refined_structure()
            s = sga(structure)
            sym_struc = s.get_symmetrized_structure()
            number = s.get_space_group_number()
        except:
            print("Failed to load the Pymatgen structure")
            self.valid = False

        if self.valid:
            d = sym_struc.composition.as_dict()
            species = [key for key in d.keys()]
            numIons = []
            for ele in species:
                numIons.append(int(d[ele]))
            self.numIons = numIons
            self.species = species
            self.group = Group(number)
            atom_sites = []
            for i, site in enumerate(sym_struc.equivalent_sites):
                pos = site[0].frac_coords
                wp = Wyckoff_position.from_group_and_index(number, sym_struc.wyckoff_symbols[i])
                specie = site[0].specie.number
                atom_sites.append(atom_site(wp, pos, specie))
            self.atom_sites = atom_sites
            self.lattice = Lattice.from_matrix(sym_struc.lattice.matrix, ltype=self.group.lattice_type)

    def check_short_distances(self, r=0.7, exclude_H = True):
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
                pmg_struc.remove_species('H')
            res = pmg_struc.get_all_neighbors(r)
            for i, neighs in enumerate(res):
                for n in neighs:
                    pairs.append([pmg_struc.sites[i].specie, n.specie, n.nn_distance])
        else:
            raise NotImplementedError("Does not support cluster for now")
        return pairs

    def to_file(self, filename=None, fmt=None, permission='w'):
        """
        Creates a file with the given filename and file type to store the structure.
        By default, creates cif files for crystals and xyz files for clusters.
        For other formats, Pymatgen is used

        Args:
            filename: the file path
            fmt: the file type (`cif`, `xyz`, etc.)
            permission: `w` or `a+`

        Returns:
            Nothing. Creates a file at the specified path
        """
        if self.valid:
            if fmt is None:
                if self.dim == 0:
                    fmt = 'xyz'
                else:
                    fmt = 'cif'

            if fmt == "cif":
                if self.dim == 3:
                    from pyxtal.io import write_cif
                    return write_cif(self, filename, "from_pyxtal", permission)
                else:
                    pmg_struc = self.to_pymatgen()
                return pmg_struc.to(fmt=fmt, filename=filename)
            else:
                pmg_struc = self.to_pymatgen()
                return pmg_struc.to(fmt=fmt, filename=filename)
        else:
            raise RuntimeError("Cannot create file: structure did not generate")


    def subgroup(self, H=None, eps=0.05, idx=None, once=False, group_type='t'):
        """
        generate a structure with lower symmetry

        Args:
            H: space group number (int)
            eps: pertubation term (float)
            idx: list
            once: generate only one structure, otherwise output all

        Returns:
            a list of pyxtal structures with lower symmetries
        """

        #randomly choose a subgroup from the available list
        if group_type == 't':
            Hs = self.group.get_max_t_subgroup()['subgroup']
        else:
            Hs = self.group.get_max_k_subgroup()['subgroup']

        if idx is None:
            idx = range(len(Hs))
        else:
            for id in idx:
                if id >= len(Hs):
                    raise ValueError("The idx exceeds the number of possible splits")
        
        if H is not None: 
            idx = [id for id in idx if Hs[idx] == H]

        if len(idx) == 0:
            raise ValueError("The space group H is incompatible with idx")

        sites = [str(site.wp.multiplicity)+site.wp.letter for site in self.atom_sites]
        valid_splitters = []
        bad_splitters = []
        for id in idx:
            splitter = wyckoff_split(G=self.group.number, wp1=sites, idx=id, group_type=group_type)
            if splitter.valid_split:
                valid_splitters.append(splitter)
            else:
                bad_splitters.append(splitter)

        if len(valid_splitters) == 0:
            # do one more step
            new_strucs = []
            for splitter in bad_splitters:
                trail_struc = self.subgroup_by_splitter(splitter)
                new_strucs.append(trail_struc.subgroup(once=True, group_type=group_type))
            return new_strucs
        else:
            if once:
                return self.subgroup_by_splitter(choice(valid_splitters), eps=eps)
            else:
                new_strucs = []
                for splitter in valid_splitters:
                    new_strucs.append(self.subgroup_by_splitter(splitter, eps=eps))
            return new_strucs

    def subgroup_by_splitter(self, splitter, eps=0.05):
        lat1 = np.dot(splitter.R[:3,:3].T, self.lattice.matrix)
        multiples = np.linalg.det(splitter.R[:3,:3])
        split_sites = []
        for i, site in enumerate(self.atom_sites):
            pos = site.position
            for ops1, ops2 in zip(splitter.G2_orbits[i], splitter.H_orbits[i]):
                pos0 = apply_ops(pos, ops1)[0]
                pos0 -= np.floor(pos0)
                pos0 += eps*(np.random.sample(3) - 0.5)
                wp, _ = Wyckoff_position.from_symops(ops2, group=splitter.H.number, permutation=False)
                split_sites.append(atom_site(wp, pos0, site.specie))
        new_struc = deepcopy(self)
        new_struc.group = splitter.H
        lattice = Lattice.from_matrix(lat1, ltype=new_struc.group.lattice_type)
        new_struc.lattice=lattice.mutate(degree=eps, frozen=True)
        new_struc.atom_sites = split_sites
        new_struc.numIons = [int(multiples*numIon) for numIon in self.numIons]
        new_struc.source = 'Wyckoff Split'
        
        return new_struc

    def show(self, **kwargs):
        """
        display the crystal structure
        """
        from pyxtal.viz import display_atomic
        return display_atomic(self, **kwargs)

    def copy(self):
        """
        simply copy the structure
        """
        return deepcopy(self)

    def _get_coords_and_species(self, absolute=False):
        """
        extract the coordinates and species information 

        Args:
            abosulte: if True, return the cartesian coords otherwise fractional

        Returns:
            total_coords (N*3 numpy array) and the list of species
        """
        species = []
        total_coords = None
        for site in self.atom_sites:
            species.extend([site.specie]*site.multiplicity)
            if total_coords is None:
                total_coords = site.coords
            else:
                total_coords = np.append(total_coords, site.coords, axis=0)

        if absolute:
            return total_coords.dot(self.lattice.matrix), species
        else:
            return total_coords, species

    def to_ase(self):
        """
        export to ase Atoms object
        """
        from ase import Atoms
        if self.valid:
            if self.dim > 0:
                lattice = self.lattice.copy()
                coords, species = self._get_coords_and_species()
                # Add space above and below a 2D or 1D crystals
                latt, coords = lattice.add_vacuum(coords, PBC=self.PBC)
                return Atoms(species, scaled_positions=coords, cell=latt, pbc=self.PBC)
            else:
                coords, species = self._get_coords_and_species(True)
                return Atoms(species, positions=coords)
        else:
            raise RuntimeError("No valid structure can be converted to ase.")

    def to_pymatgen(self):
        """
        export to Pymatgen structure object
        """
        from pymatgen.core.structure import Structure, Molecule

        if self.valid:
            if self.dim > 0:
                lattice = self.lattice.copy()
                coords, species = self._get_coords_and_species()
                # Add space above and below a 2D or 1D crystals
                latt, coords = lattice.add_vacuum(coords, PBC=self.PBC)
                return Structure(latt, species, coords)
            else:
                # Clusters are handled as large molecules
                coords, species = self._get_coords_and_species(True)
                return Molecule(species, coords)
        else:
            raise RuntimeError("No valid structure can be converted to pymatgen.")


