# Standard Libraries
import os
import random
import numpy as np
from copy import deepcopy

# External Libraries
from pymatgen.core.structure import Structure, Molecule

# PyXtal imports #avoid *
from pyxtal.symmetry import Group, choose_wyckoff, check_wyckoff_position
from pyxtal.wyckoff_site import atom_site, check_atom_sites
from pyxtal.operations import (
    apply_ops,
    project_point,
    distance_matrix,
    distance,
)
from pyxtal.msg import printx
from pyxtal.tolerance import Tol_matrix
from pyxtal.lattice import Lattice, cellsize, add_vacuum
from pyxtal.database.element import Element

# Define functions
# ------------------------------


def merge_coordinate(coor, lattice, group, tol):
    """
    Given a list of fractional coordinates, merges them within a given
    tolerance, and checks if the merged coordinates satisfy a Wyckoff
    position. Used for merging general Wyckoff positions into special Wyckoff
    positions within the random_crystal (and its derivative) classes.

    Args:
        coor: a list of fractional coordinates
        lattice: a 3x3 matrix representing the unit cell
        group: a pyxtal.symmetry.Group object
        tol: the cutoff distance for merging coordinates

    Returns:
        coor, index, point: (coor) is the new list of fractional coordinates after
        merging. index is a single index for the Wyckoff position within
        the sg. If no matching WP is found, returns False. point is a 3-vector;
        when plugged into the Wyckoff position, it will generate all the other
        points.
    """
    coor = np.array(coor)
    # Get index of current Wyckoff position. If not one, return False
    index, point = check_wyckoff_position(coor, group)
    if index is False:
        return coor, False, None
    if point is None:
        printx("Error: Could not find generating point.", priority=1)
        printx("coordinates:")
        printx(str(coor))
        printx("Lattice: ")
        printx(str(lattice))
        printx("group: ")
        group.print_all()
        return coor, False, None
    PBC = group.PBC
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

        if passed_distance_check is False:
            mult1 = group[index].multiplicity
            # Find possible wp's to merge into
            possible = []
            for i, wp in enumerate(group):
                mult2 = wp.multiplicity
                # factor = mult2 / mult1
                if (mult2 < mult1) and (mult1 % mult2 == 0):
                    possible.append(i)
            if possible == []:
                return coor, False, None
            # Calculate minimum separation for each WP
            distances = []
            for i in possible:
                wp = group[i]
                projected_point = project_point(point, wp[0], lattice=lattice, PBC=PBC)
                # NOTE Comprhys: new_coor note used?
                new_coor = apply_ops(projected_point, wp)
                d = distance(point - projected_point, lattice, PBC=PBC)
                distances.append(np.min(d))
            # Choose wp with shortest translation for generating point
            tmpindex = np.argmin(distances)
            index = possible[tmpindex]
            newwp = group[index]
            projected_point = project_point(point, newwp[0], lattice=lattice, PBC=PBC)
            coor = apply_ops(projected_point, newwp)
            point = coor[0]
            index = newwp.index
        # Distances were not too small; return True
        else:
            return coor, index, point


def check_compatible(group, numIons):
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


class random_crystal:
    """
    Class for storing and generating atomic crystals based on symmetry
    constraints. Given a spacegroup, list of atomic symbols, the stoichiometry,
    and a volume factor, generates a random crystal consistent with the
    spacegroup's symmetry. This crystal is stored as a pymatgen struct via
    self.struct

    Args:
        group: the international spacegroup number, or a Group object
        species: a list of atomic symbols for each ion type, e.g. ["Ti", "O"]
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell), e.g.[4, 2]
        tm: the Tol_matrix object used to generate the crystal
        factor: a volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between atoms
        sites: pre-assigned wyckoff sites (e.g., [["4a"], ["2b"]])
        lattice: an optional Lattice object to use for the unit cell
    """

    def __init__(
        self,
        group,
        species,
        numIons,
        factor,
        lattice=None,
        sites = None,
        tm=Tol_matrix(prototype="atomic"),
    ):

        self.dim = 3
        """The number of periodic dimensions of the crystal"""
        if type(group) != Group:
            group = Group(group, self.dim)
        self.sg = group.number
        """The international spacegroup number of the crystal."""
        self.PBC = [1, 1, 1]
        """The periodic boundary axes of the crystal"""
        self.init_common(species, numIons, factor, group, lattice, sites, tm)

    def init_common(self, species, numIons, factor, group, lattice, sites, tm):
        """
        Common init functionality for 0D-3D cases of random_crystal.
        """
        self.valid = False
        # Check that numIons are integers greater than 0
        for num in numIons:
            if int(num) != num or num < 1:
                printx(
                    "Error: stoichiometry must consist of positive integers.", priority=1,
                )
                return False
        if type(group) == Group:
            """
            A pyxtal.symmetry.Group object storing information
            about the space/layer/Rod/point group,
            and its Wyckoff positions.
            """
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

        # The number of attempts to generate the crystal, max1*max2*max3.
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
                self.struct = None
                return
        
        self.sites = {}
        for i, specie in enumerate(self.species):
            if sites is not None and sites[i] is not None:
                self.check_consistency(sites[i], self.numIons[i])
                self.sites[specie] = sites[i]
            else:
                self.sites[specie] = None
        # QZ: needs to check if it is compatible

        self.lattice_attempts=40
        self.coord_attempts=10
        self.wyckoff_attempts=10

        self.generate_crystal()


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

    def to_file(self, fmt="cif", filename=None):
        """
        Creates a file with the given filename and file type to store the structure.
        By default, creates cif files for crystals and xyz files for clusters.
        By default, the filename is based on the stoichiometry.

        Args:
            fmt: the file type ('cif', 'xyz', etc.)
            filename: the file path

        Returns:
            Nothing. Creates a file at the specified path
        """
        if filename is None:
            given = False
        else:
            given = True
        if self.valid:
            if self.dim == 0:
                if filename is None:
                    filename = str(self.molecule.formula).replace(" ", "") + "." + fmt
            if self.dim != 0:
                if filename is None:
                    filename = str(self.struct.formula).replace(" ", "") + "." + fmt

            # Check if filename already exists
            # If it does, add a new number to end of filename

            if os.path.exists(filename):
                if given is False:
                    filename = filename[: (-len(fmt) - 1)]
                i = 1
                while True:
                    outdir = filename + "_" + str(i)
                    if given is False:
                        outdir += "." + fmt
                    if not os.path.exists(outdir):
                        break
                    i += 1
                    if i > 10000:
                        return "Could not create file: too many files already created."
            else:
                outdir = filename
            if self.dim == 0 and fmt == "xyz":
                self.molecule.to(fmt=fmt, filename=outdir)
            else:
                self.struct.to(fmt=fmt, filename=outdir)
            return "Output file to " + outdir
        elif self.valid is False:
            printx("Cannot create file: structure did not generate.", priority=1)

    def __str__(self):
        s = "------Random Crystal------"
        s += "\nComposition: {}".format(self.struct.formula)
        s += "\nDimension: {}".format(self.dim)
        s += "\nGroup: {} ({})".format(self.group.symbol, self.group.number)
        s += "\nVolume factor: {}".format(self.factor)
        s += "\n{}".format(self.lattice)
        if self.valid:
            s += "\nWyckoff sites:"
            for wyc in self.wyckoff_sites:
                s += "\n\t{}".format(wyc)
        else:
            s += "\nStructure not generated."
        return s

    def __repr__(self):
        return str(self)

    def print_all(self):
        """
        Prints useful information about the generated crystal.
        """
        print(str(self))

    def generate_crystal(self):
        """
        The main code to generate a random atomic crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

       """
        # Check the minimum number of degrees of freedom within the Wyckoff positions
        self.numattempts = 1
        degrees = check_compatible(self.group, self.numIons)
        if degrees is False:
            printx(self.Msg1, priority=1)
            self.struct = None
            self.valid = False
            return

        if degrees == 0:
            printx("Wyckoff positions have no degrees of freedom.", priority=2)
            # NOTE why do these need to be changed from defaults?
            self.lattice_attempts = 5
            self.coord_attempts = 5
            self.wyckoff_attempts = 5

        # Calculate a minimum vector length for generating a lattice
        # NOTE Comprhys: minvector never used?
        # minvector = max(self.tol_matrix.get_tol(s, s) for s in self.species)
        for cycle1 in range(self.lattice_attempts):
            self.cycle1 = cycle1

            # 1, Generate a lattice
            if self.lattice.allow_volume_reset is True:
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
            if self.lattice.random is True:
                if self.dim != 0 and abs(self.volume - self.lattice.volume) > 1.0:
                    printx(
                        (
                            "Error, volume is not equal to the estimated value: "
                            "{} -> {} cell_para: {}"
                        ).format(self.volume, self.lattice.volume, self.lattice.get_para),
                        priority=0,
                    )
                    self.valid = False
                    self.struct = None
                    return

            # to try to generate atomic coordinates
            for cycle2 in range(self.coord_attempts):
                self.cycle2 = cycle2
                output = self._generate_coords(cell_matrix)

                if output:
                    final_coords, final_species, wyckoff_sites = output
                    break
            if self.valid:
                final_number = [Element(ele).z for ele in final_species]
                final_coords = np.array(final_coords)
                cart_coords = np.dot(final_coords, cell_matrix)

                if self.dim != 0:
                    # Add space above and below a 2D or 1D crystals
                    cell_matrix, final_coords = add_vacuum(
                        cell_matrix, final_coords, PBC=self.PBC
                    )
                    self.struct = Structure(cell_matrix, final_species, final_coords)
                    self.spg_struct = (cell_matrix, final_coords, final_number)
                else:
                    self.species = final_species
                    # Clusters are handled as large molecules
                    self.molecule = Molecule(final_species, cart_coords)
                    # Calculate binding box
                    diffs = np.max(cart_coords, axis=0) - np.min(cart_coords, axis=0) + 10
                    # a, b, c = diffs[0], diffs[1], diffs[2]
                    self.struct = self.molecule.get_boxed_structure(*diffs)

                self.sites = final_species
                self.frac_coords = final_coords
                self.cart_coords = np.dot(final_coords, cell_matrix)
                self.lattice_matrix = cell_matrix
                self.wyckoff_sites = wyckoff_sites
                self.valid = True
                return

        self.valid = False
        self.struct = None
        return

    def _generate_coords(self, cell_matrix):
        """
        """
        coordinates_list = []  # to store the added coordinates
        species_list = []  # to store the corresponding specie
        wyckoff_sites_list = []

        # generate coordinates for each ion type in turn
        for numIon, specie in zip(self.numIons, self.species):
            output = self._generate_ion_wyckoffs(
                numIon, specie, cell_matrix, wyckoff_sites_list
            )

            if output:
                coords, species, wyks = output
                coordinates_list.extend(coords)
                species_list.extend(species)
                wyckoff_sites_list.extend(wyks)
            else:
                # correct multiplicity not achieved exit and start over
                return None

        # If numIon_added correct for all specie return structure
        self.valid = True
        return coordinates_list, species_list, wyckoff_sites_list

    def _generate_ion_wyckoffs(self, numIon, specie, cell_matrix, wyks):
        """
        generates a set of wyckoff positions to accomodate a given number
        of ions

        Args:
            numIon: Number of ions to accomodate
            specie: Type of species being placed on wyckoff site
            cell_matrix: Matrix of lattice vectors

        Returns:
            Sucess:
                coordinates_tmp: list of sets of co-ordinates for valid sites
                species_tmp: list of specie for valid sites
                wyckoff_sites_tmp: list of wyckoff sites for valid sites
            Failue:
                None

        """
        numIon_added = 0
        tol = self.tol_matrix.get_tol(specie, specie)

        coordinates_tmp = []
        species_tmp = []
        wyckoff_sites_tmp = []

        # Now we start to add the specie to the wyckoff position
        sites_list = deepcopy(self.sites[specie]) # the list of Wyckoff site
        if sites_list is not None: 
            self.wyckoff_attempts = len(sites_list)*2

        cycle3 = 0
        while cycle3 < self.wyckoff_attempts:

            self.cycle3 = cycle3
            # Choose a random WP for given multiplicity: 2a, 2b
            if sites_list is not None:
                site = sites_list[0]
            else: # Selecting the merging 
                site = None
            
            ops = choose_wyckoff(self.group, numIon - numIon_added, site, self.dim)
            if ops is not False:
                # Generate a list of coords from ops
                pt = self.lattice.generate_point()
                proj_pt = project_point(pt, ops[0], cell_matrix, self.PBC)
                coords = apply_ops(proj_pt, ops)


                # Merge coordinates if the atoms are close
                coords_toadd, wp_index, pt = merge_coordinate(
                    coords, cell_matrix, self.group, tol
                )
                if site is not None and len(coords_toadd) < len(coords):
                    continue # break the cycle if the merge happens
                #print(coords)
                #print(coords_toadd, wp_index, pt)
                #import sys
                #sys.exit()

                if wp_index is not False:
                    # Use a Wyckoff_site object for the current site
                    current_site = atom_site(self.group[wp_index], pt, specie)

                    # Check current WP against existing WP's
                    passed_wp_check = True
                    for ws in wyckoff_sites_tmp + wyks:
                        if not check_atom_sites(
                            current_site, ws, cell_matrix, self.tol_matrix,
                        ):
                            passed_wp_check = False

                    if passed_wp_check is True:
                        if sites_list is not None:
                            sites_list.pop(0)

                        # The current Wyckoff site passed; store it
                        if coordinates_tmp == []:
                            coordinates_tmp = current_site.coords
                        else:
                            coordinates_tmp = np.vstack(
                                [coordinates_tmp, current_site.coords]
                            )
                        species_tmp.extend([specie] * len(coords_toadd))
                        wyckoff_sites_tmp.append(current_site)
                        numIon_added += len(coords_toadd)

                        # Check if enough atoms have been added
                        if numIon_added == numIon:
                            return coordinates_tmp, species_tmp, wyckoff_sites_tmp

            cycle3 += 1
            self.numattempts += 1

        return None


class random_crystal_2D(random_crystal):
    """
    A 2d counterpart to random_crystal. Generates a random atomic crystal based
    on a 2d layer group instead of a 3d spacegroup. Note that each layer group
    is equal to a corresponding 3d spacegroup, but without periodicity in one
    direction. The generated pymatgen structure can be accessed via self.struct

    Args:
        group: the layer group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
            OR, a pyxtal.symmetry.Group object
        species: a list of atomic symbols for each ion type
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell)
        thickness: the thickness, in Angstroms, of the unit cell in the 3rd
            dimension (the direction which is not repeated periodically)
        factor: a volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between atoms
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        species,
        numIons,
        factor,
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
        group: the Rod group number between 1 and 75. NOT equal to the
            international space group number, which is between 1 and 230
            OR, a pyxtal.symmetry.Group object
        species: a list of atomic symbols for each ion type
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell)
        area: the effective cross-sectional area, in Angstroms squared, of the
            unit cell
        factor: a volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between atoms
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        species,
        numIons,
        factor,
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
        group: the Schoenflies symbol for the point group (ex: "Oh", "C5v", "D3")
            OR the number between 1-32 for a crystallographic point group,
            OR, a pyxtal.symmetry.Group object
            See:
            https://en.wikipedia.org/wiki/Schoenflies_notation#Point_groups
            for more information
        species: a list of atomic symbols for each ion type
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell)
        factor: a volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between atoms
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        species,
        numIons,
        factor,
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
