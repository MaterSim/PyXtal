# Standard Libraries
import os
import random
import numpy as np
from copy import deepcopy
from random import choice

# PyXtal imports #avoid *
from pyxtal.symmetry import Group, choose_wyckoff, Wyckoff_position
from pyxtal.wyckoff_site import atom_site, check_atom_sites, WP_merge
from pyxtal.wyckoff_split import wyckoff_split
from pyxtal.msg import printx
from pyxtal.tolerance import Tol_matrix
from pyxtal.lattice import Lattice, cellsize
from pyxtal.database.element import Element
from pyxtal.operations import apply_ops

# Define functions
# ------------------------------
class random_crystal:
    """
    Class for storing and generating atomic crystals based on symmetry
    constraints. Given a spacegroup, list of atomic symbols, the stoichiometry,
    and a volume factor, generates a random crystal consistent with the
    spacegroup's symmetry. 

    Args:
        group: the spacegroup number (1-230), or a 
            `pyxtal.symmetry.Group <pyxtal.symmetry.Group.html>`_ object
        species: a list of atomic symbols for each ion type, e.g., `["Ti", "O"]`
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell), e.g., `[4, 2]`
        factor (optional): volume factor used to generate the crystal
        sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
        lattice (optional): the `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_ 
            object to define the unit cell
        tm (optional): the `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_ 
            object to define the distances
        seed (optional): the cif/POSCAR file from user
    """

    def __init__(
        self,
        group=None,
        species=None,
        numIons=None,
        factor=1.1,
        lattice=None,
        sites = None,
        tm=Tol_matrix(prototype="atomic"),
        seed = None,
    ):

        self.dim = 3 #periodic dimensions of the crystal
        self.PBC = [1, 1, 1] #The periodic boundary axes of the crystal
        if seed is None:
            if type(group) != Group:
                group = Group(group, self.dim)
            self.sg = group.number #The international spacegroup number 
            self.seed = None
            self.init_common(species, numIons, factor, group, lattice, sites, tm)
        else:
            self.seed = seed
            self.from_seed()


    def from_seed(self):
        """
        Load the seed structure from Pymatgen/ASE/POSCAR/CIFs
        Internally they will be handled by Pymatgen
        """
        self.valid = True
        from ase import Atoms
        from pymatgen import Structure
        if isinstance(self.seed, Atoms): #ASE atoms
            from pymatgen.io.ase import AseAtomsAdaptor
            pmg_struc = AseAtomsAdaptor.get_structure(self.seed)
            self.from_pymatgen(pmg_struc)
        elif isinstance(self.seed, Structure): #Pymatgen
            self.from_pymatgen(self.seed)
        elif isinstance(self.seed, str):
            pmg_struc = Structure.from_file(self.seed)
            self.from_pymatgen(pmg_struc)

        formula = ""
        for i, s in zip(self.numIons, self.species):
            formula += "{:s}{:d}".format(s, int(i))
        self.formula = formula
        self.factor = 1.0
        self.number = self.group.number
        self.source = 'Seed'

    def from_pymatgen(self, structure):
        """
        Load the seed structure from Pymatgen/ASE/POSCAR/CIFs
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
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


    def init_common(self, species, numIons, factor, group, lattice, sites, tm):
        """
        Common init functionality for 0D-3D cases of random_crystal.
        """
        self.source = 'Random'
        self.valid = False
        # Check that numIons are integers greater than 0
        for num in numIons:
            if int(num) != num or num < 1:
                printx("Error: composition must be positive integers.", priority=1)
                return False
        if type(group) == Group:
            self.group = group
        else:
            self.group = Group(group, dim=self.dim)
        self.number = self.group.number
        """
        The international group number of the crystal:
        1-230 for 3D space groups
        1-80 for 2D layer groups
        1-75 for 1D Rod groups
        1-32 for crystallographic point groups
        None otherwise
        """

        # The number of attempts to generate the crystal
        # number of atoms
        # volume factor for the unit cell.
        # The number of atom in the PRIMITIVE cell
        # The number of each type of atom in the CONVENTIONAL cell.
        # A list of atomic symbols for the types of atoms
        # A list of warning messages

        self.numattempts = 0
        numIons = np.array(numIons)
        self.factor = factor
        self.numIons0 = numIons
        self.numIons = self.numIons0 * cellsize(self.group)
        formula = ""
        for i, s in zip(self.numIons, species):
            formula += "{:s}{:d}".format(s, int(i))
        self.formula = formula

        self.species = species
        self.Msgs()

        # Use the provided lattice
        if lattice is not None:
            self.lattice = lattice
            self.volume = lattice.volume
            # Make sure the custom lattice PBC axes are correct.
            if lattice.PBC != self.PBC:
                self.lattice.PBC = self.PBC
                printx("\n  Warning: converting custom lattice PBC to " + str(self.PBC))

        # Generate a Lattice instance based on a given volume estimation
        elif lattice is None:

            # Determine the unique axis
            if self.dim == 2:
                if self.number in range(3, 8):
                    unique_axis = "c"
                else:
                    unique_axis = "a"
            elif self.dim == 1:
                if self.number in range(3, 8):
                    unique_axis = "a"
                else:
                    unique_axis = "c"
            else:
                unique_axis = "c"

            self.volume = self.estimate_volume()

            if self.dim == 3 or self.dim == 0:
                self.lattice = Lattice(
                    self.group.lattice_type,
                    self.volume,
                    PBC=self.PBC,
                    unique_axis=unique_axis,
                )
            elif self.dim == 2:
                self.lattice = Lattice(
                    self.group.lattice_type,
                    self.volume,
                    PBC=self.PBC,
                    unique_axis=unique_axis,
                    # NOTE self.thickness is part of 2D class
                    thickness=self.thickness,
                )
            elif self.dim == 1:
                self.lattice = Lattice(
                    self.group.lattice_type,
                    self.volume,
                    PBC=self.PBC,
                    unique_axis=unique_axis,
                    # NOTE self.area is part of 1D class
                    area=self.area,
                )
        # Set the tolerance matrix for checking inter-atomic distances
        if type(tm) == Tol_matrix:
            self.tol_matrix = tm
        else:
            try:
                self.tol_matrix = Tol_matrix(prototype=tm)
            # TODO Remove bare except
            except:
                printx(
                    (
                        "Error: tm must either be a Tol_matrix object or "
                        "a prototype string for initializing one."
                    ),
                    priority=1,
                )
                self.valid = False
                return
        
        self.sites = {}
        for i, specie in enumerate(self.species):
            if sites is not None and sites[i] is not None:
                self.check_consistency(sites[i], self.numIons[i])
                self.sites[specie] = sites[i]
            else:
                self.sites[specie] = None
        # QZ: needs to check if it is compatible

        self.generate_crystal()

    def check_compatible(self, group, numIons):
        """
        Checks if the number of atoms is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
    
        NOTE Comprhys: Is degrees of freedom used symnomously with multiplicity?
        perhaps standardising to multiplicity would be clearer?
        """
        # Store whether or not at least one degree of freedom exists
        has_freedom = False
        # Store the wp's already used that don't have any freedom
        used_indices = []
        # Loop over species
        for numIon in numIons:
            # Get lists of multiplicity, maxn and freedom
            l_mult0 = []
            l_maxn0 = []
            l_free0 = []
            indices0 = []
            for i_wp, wp in enumerate(group):
                indices0.append(i_wp)
                l_mult0.append(len(wp))
                l_maxn0.append(numIon // len(wp))
                if np.allclose(wp[0].rotation_matrix, np.zeros([3, 3])):
                    l_free0.append(False)
                else:
                    l_free0.append(True)
            # Remove redundant multiplicities:
            l_mult = []
            l_maxn = []
            l_free = []
            indices = []
            for mult, maxn, free, i_wp in zip(l_mult0, l_maxn0, l_free0, indices0):
                if free is True:
                    if mult not in l_mult:
                        l_mult.append(mult)
                        l_maxn.append(maxn)
                        l_free.append(True)
                        indices.append(i_wp)
                elif free is False and i_wp not in used_indices:
                    l_mult.append(mult)
                    indices.append(i_wp)
                    if mult <= numIon:
                        l_maxn.append(1)
                    elif mult > numIon:
                        l_maxn.append(0)
                    l_free.append(False)
    
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
                return False
            while True:
                num = np.dot(n, l_mult)
                dobackwards = False
                # The combination works: move to next species
                if num == numIon:
                    # Check if at least one degree of freedom exists
                    for val, free, i_wp in zip(n, l_free, indices):
                        if val > 0:
                            if free is True:
                                has_freedom = True
                            elif free is False:
                                used_indices.append(i_wp)
                    break
                # All combinations failed: return False
                if n == n0 and p >= len(l_mult) - 1:
                    return False
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
                if num > numIon or dobackwards is True:
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
        # All species passed: return True
        if has_freedom is True:
            return True
        # All species passed, but no degrees of freedom: return 0
        elif has_freedom is False:
            return 0
    


    def check_consistency(self, site, numIon):
        num = 0
        for s in site:
            num += int(s[:-1])
        if numIon == num:
            return True
        else:
            msg = "\nThe requested number of atoms is inconsistent: " + str(site)
            msg += "\nfrom numIons: {:d}".format(numIon)
            msg += "\nfrom Wyckoff list: {:d}".format(num)
            raise ValueError(msg)

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
        if dim > 0:
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

    def Msgs(self):
        """
        Define a set of error and warning message if generation fails.

        Returns:
            nothing
        """
        self.Msg1 = (
            "Error: the stoichiometry is incompatible with the wyckoff sites choice"
        )
        self.Msg2 = "Error: failed in the cycle of generating structures"
        self.Msg3 = "Warning: failed in the cycle of adding species"
        self.Msg4 = "Warning: failed in the cycle of choosing wyckoff sites"
        self.Msg5 = "Finishing: added the specie"
        self.Msg6 = "Finishing: added the whole structure"
        self.Msg7 = "Error: invalid paramaters for initialization"

    def estimate_volume(self):
        """
        Estimates the volume of a unit cell based on the number and types of ions.
        Assumes each atom takes up a sphere with radius equal to its covalent bond
        radius.

        Returns:
            a float value for the estimated volume
        """
        volume = 0
        for numIon, specie in zip(self.numIons, self.species):
            r = random.uniform(
                Element(specie).covalent_radius, Element(specie).vdw_radius
            )
            volume += numIon * 4 / 3 * np.pi * r ** 3
        return self.factor * volume

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
            printx("Cannot create file: structure did not generate.", priority=1)

    def __str__(self):
        s = "------Crystal from {:s}------".format(self.source)
        s += "\nComposition: {}".format(self.formula)
        s += "\nDimension: {}".format(self.dim)
        s += "\nGroup: {} ({})".format(self.group.symbol, self.group.number)
        s += "\nVolume factor: {}".format(self.factor)
        s += "\n{}".format(self.lattice)
        if self.valid:
            s += "\nWyckoff sites:"
            for wyc in self.atom_sites:
                s += "\n\t{}".format(wyc)
        else:
            s += "\nStructure not generated."
        return s

    def __repr__(self):
        return str(self)

    def show(self, **kwargs):
        """
        display the crystal structure
        """
        from pyxtal.viz import display_atomic
        return display_atomic(self, **kwargs)

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



    def generate_crystal(self):
        """
        The main code to generate a random atomic crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

       """
        # Check the minimum number of degrees of freedom within the Wyckoff positions
        self.numattempts = 1
        degrees = self.check_compatible(self.group, self.numIons)
        if degrees is False:
            printx(self.Msg1, priority=1)
            self.valid = False
            return

        if degrees == 0:
            printx("Wyckoff positions have no degrees of freedom.", priority=2)
            # NOTE why do these need to be changed from defaults?
            self.lattice_attempts = 5
            self.coord_attempts = 5
            #self.wyckoff_attempts = 5
        else:
            self.lattice_attempts=40
            self.coord_attempts=10
            #self.wyckoff_attempts=10

        # Calculate a minimum vector length for generating a lattice
        # NOTE Comprhys: minvector never used?
        # minvector = max(self.tol_matrix.get_tol(s, s) for s in self.species)
        for cycle1 in range(self.lattice_attempts):
            self.cycle1 = cycle1

            # 1, Generate a lattice
            if self.lattice.allow_volume_reset:
                self.volume = self.estimate_volume()
                self.lattice.volume = self.volume
            self.lattice.reset_matrix()

            try:
                cell_matrix = self.lattice.get_matrix()
                if cell_matrix is None:
                    continue
            # TODO remove bare except
            except:
                continue

            # Check that the correct volume was generated
            if self.lattice.random:
                if self.dim != 0 and abs(self.volume - self.lattice.volume) > 1.0:
                    printx(
                        (
                            "Error, volume is not equal to the estimated value: "
                            "{} -> {} cell_para: {}"
                        ).format(self.volume, self.lattice.volume, self.lattice.get_para),
                        priority=0,
                    )
                    self.valid = False
                    return

            # to try to generate atomic coordinates
            for cycle2 in range(self.coord_attempts):
                self.cycle2 = cycle2
                output = self._generate_coords(cell_matrix)

                if output:
                    self.atom_sites = output
                    break

            if self.valid:
                return

        return

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
            printx("No valid structure can be converted to ase.", priority=1)

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
            printx("No valid structure can be converted to pymatgen.", priority=1)


    def _generate_coords(self, cell_matrix):
        """
        generate coordinates for random crystal
        """
        wyckoff_sites_list = []

        # generate coordinates for each ion type in turn
        for numIon, specie in zip(self.numIons, self.species):
            output = self._generate_ion_wyckoffs(
                numIon, specie, cell_matrix, wyckoff_sites_list
            )
            if output is not None:
                wyckoff_sites_list.extend(output)
            else:
                # correct multiplicity not achieved exit and start over
                return None

        # If numIon_added correct for all specie return structure
        self.valid = True
        return wyckoff_sites_list

    def _generate_ion_wyckoffs(self, numIon, specie, cell_matrix, wyks):
        """
        generates a set of wyckoff positions to accomodate a given number
        of ions

        Args:
            numIon: Number of ions to accomodate
            specie: Type of species being placed on wyckoff site
            cell_matrix: Matrix of lattice vectors
            wyks: current wyckoff sites

        Returns:
            Sucess:
                wyckoff_sites_tmp: list of wyckoff sites for valid sites
            Failue:
                None

        """
        numIon_added = 0
        tol = self.tol_matrix.get_tol(specie, specie)
        tol_matrix = self.tol_matrix
        wyckoff_sites_tmp = []

        # Now we start to add the specie to the wyckoff position
        sites_list = deepcopy(self.sites[specie]) # the list of Wyckoff site
        if sites_list is not None: 
            wyckoff_attempts = max(len(sites_list)*2, 10)
        else:
            # the minimum numattempts is to put all atoms to the general WPs
            min_wyckoffs = int(numIon/len(self.group.wyckoffs_organized[0][0]))
            wyckoff_attempts = max(2*min_wyckoffs, 10)

        cycle = 0
        while cycle < wyckoff_attempts:
            # Choose a random WP for given multiplicity: 2a, 2b
            if sites_list is not None:
                site = sites_list[0]
            else: # Selecting the merging 
                site = None
            
            wp = choose_wyckoff(self.group, numIon - numIon_added, site, self.dim)
            if wp is not False:
                # Generate a list of coords from ops
                mult = wp.multiplicity # remember the original multiplicity
                pt = self.lattice.generate_point()
                # Merge coordinates if the atoms are close
                pt, wp, _ = WP_merge(pt, cell_matrix, wp, tol)
                # For pure planar structure
                if self.dim == 2 and self.thickness is not None and self.thickness < 0.1:
                    pt[-1] = 0.5

               
                # If site the pre-assigned, do not accept merge
                if wp is not False:
                    if site is not None and mult != wp.multiplicity:
                        cycle += 1
                        continue 
                    # Use a Wyckoff_site object for the current site
                    new_site = atom_site(wp, pt, specie)

                    # Check current WP against existing WP's
                    passed_wp_check = True
                    for ws in wyckoff_sites_tmp + wyks:
                        if not check_atom_sites(new_site, ws, cell_matrix, tol_matrix):
                            passed_wp_check = False

                    if passed_wp_check:
                        if sites_list is not None:
                            sites_list.pop(0)
                        wyckoff_sites_tmp.append(new_site)
                        numIon_added += new_site.multiplicity 

                        # Check if enough atoms have been added
                        if numIon_added == numIon:
                            return wyckoff_sites_tmp

            cycle += 1
            self.numattempts += 1

        return None


class random_crystal_2D(random_crystal):
    """
    A 2d counterpart to random_crystal. Generates a random atomic crystal based
    on a 2d layer group instead of a 3d spacegroup. Note that each layer group
    is equal to a corresponding 3d spacegroup, but without periodicity in one
    direction. The generated pymatgen structure can be accessed via self.struct

    Args:
        group: the layer group number (1-80), or a 
            `pyxtal.symmetry.Group <pyxtal.symmetry.Group.html>`_ object
        species: a list of atomic symbols for each ion type, e.g., `["Ti", "O"]`
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell), e.g., `[4, 2]`
        factor (optional): volume factor used to generate the crystal
        thickness: the thickness, in Angstroms, of the unit cell in the 3rd
            dimension (the direction which is not repeated periodically)
        sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
        lattice (optional): the `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_ 
            object to define the unit cell
        tm (optional): the `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_ 
            object to define the distances
 
    """

    def __init__(
        self,
        group,
        species,
        numIons,
        factor=1.1,
        thickness=None,
        lattice=None,
        sites = None,
        tm=Tol_matrix(prototype="atomic"),
    ):
        self.dim = 2
        self.PBC = [1, 1, 0]

        if type(group) != Group:
            group = Group(group, self.dim)
        number = group.number  # The layer group number of the crystal
        self.thickness = thickness  # in Angstroms, in the 3rd dimenion of unit cell
        self.init_common(species, numIons, factor, number, lattice, sites, tm)


class random_crystal_1D(random_crystal):
    """
    A 1d counterpart to random_crystal. Generates a random atomic crystal based
    on a 1d Rod group instead of a 3d spacegroup. The generated pymatgen
    structure can be accessed via self.struct

    Args:
        group: the Rod group number (1-75), or a 
            `pyxtal.symmetry.Group <pyxtal.symmetry.Group.html>`_ object
        species: a list of atomic symbols for each ion type, e.g., `["Ti", "O"]`
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell), e.g., `[4, 2]`
        factor (optional): volume factor used to generate the crystal
        area: the effective cross-sectional area (A^2), of the unit cell
        sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
        lattice (optional): the `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_ 
            object to define the unit cell
        tm (optional): the `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_ 
            object to define the distances
 
    """

    def __init__(
        self,
        group,
        species,
        numIons,
        factor=1.1,
        area=None,
        lattice=None,
        sites = None,
        tm=Tol_matrix(prototype="atomic"),
    ):
        self.dim = 1
        self.PBC = [0, 0, 1]
        self.sg = None
        self.area = area  # the effective cross-sectional area, in A^2, of the unit cell.
        self.init_common(species, numIons, factor, group, lattice, sites, tm)


class random_cluster(random_crystal):
    """
    A 0d counterpart to random_crystal. Generates a random atomic cluster based
    on a 0d Point group instead of a 3d spacegroup. The generated pymatgen
    structure can be accessed via self.struct

    Args:
        group: the Schoenflies symbol for the point group (ex: `Oh, C5v, D3`)
            OR the number between 1-32 for a crystallographic point group,
            OR the `group <pyxtal.symmetry.Group.html>`_ object, see `wikipedia
            <https://en.wikipedia.org/wiki/Schoenflies_notation#Point_groups>`_
            for more information
        species: a list of atomic symbols for each ion type, e.g., `["Ti", "O"]`
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell), e.g., `[4, 2]`
        factor (optional): volume factor used to generate the crystal
        sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
        lattice (optional): the `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_ 
            object to define the unit cell
        tm (optional): the `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_ 
            object to define the distances
    """

    def __init__(
        self,
        group,
        species,
        numIons,
        factor=1.1,
        lattice=None,
        sites = None,
        tm=Tol_matrix(prototype="atomic", factor=0.7),
    ):
        # NOTE tol_m unused?
        tol_m = 0.1
        self.dim = 0
        self.PBC = [0, 0, 0]
        self.sg = None
        self.init_common(species, numIons, factor, group, lattice, sites, tm)
