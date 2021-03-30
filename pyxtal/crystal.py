"""
Module for generating atomic crystals
"""

# Standard Libraries
import random
from copy import deepcopy
import numpy as np

# PyXtal imports #avoid *
from pyxtal.symmetry import Group, choose_wyckoff
from pyxtal.wyckoff_site import atom_site, WP_merge
from pyxtal.msg import printx
from pyxtal.tolerance import Tol_matrix
from pyxtal.lattice import Lattice, cellsize
from pyxtal.database.element import Element

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
        lattice (optional): `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_
            object to define the unit cell
        tm (optional): `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_
            object to define the distances
    """

    def __init__(
        self,
        group=None,
        species=None,
        numIons=None,
        factor=1.1,
        lattice=None,
        sites = None,
        conventional = True,
        tm=Tol_matrix(prototype="atomic"),
    ):

        self.dim = 3 #periodic dimensions of the crystal
        self.PBC = [1, 1, 1] #The periodic boundary axes of the crystal
        self.lattice_attempts = 0
        self.coord_attempts = 0

        if type(group) != Group:
            group = Group(group, self.dim)
        self.init_common(species, numIons, factor, group, lattice, sites, conventional, tm)

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
            s = "\nStructure not available."
        return s

    def __repr__(self):
        return str(self)


    def init_common(self, species, numIons, factor, group, lattice, sites, conventional, tm):
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
        if not conventional:
            mul = cellsize(self.group)
        else:
            mul = 1
        self.numIons = numIons * mul

        formula = ""
        for i, s in zip(self.numIons, species):
            formula += "{:s}{:d}".format(s, int(i))
        self.formula = formula

        self.species = species

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
        if has_freedom:
            # All species passed: return True
            return True
        else:
            # All species passed, but no degrees of freedom: return 0
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
            msg = "Warning: the stoichiometry is incompatible with wyckoff choice"
            printx(msg, priority=1)
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
                        if not new_site.check_with_ws2(ws, cell_matrix, tol_matrix):
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
    direction.

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
        lattice (optional): `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_
            object to define the unit cell
        tm (optional): `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_
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
        conventional = True,
        tm=Tol_matrix(prototype="atomic"),
    ):
        self.dim = 2
        self.PBC = [1, 1, 0]

        if type(group) != Group:
            group = Group(group, self.dim)
        number = group.number  # The layer group number of the crystal
        self.thickness = thickness  # in Angstroms, in the 3rd dimenion of unit cell
        self.init_common(species, numIons, factor, number, lattice, sites, conventional, tm)


class random_crystal_1D(random_crystal):
    """
    A 1d counterpart to random_crystal. Generates a random atomic crystal based
    on a 1d Rod group instead of a 3d spacegroup.

    Args:
        group: the Rod group number (1-75), or a
            `pyxtal.symmetry.Group <pyxtal.symmetry.Group.html>`_ object
        species: a list of atomic symbols for each ion type, e.g., `["Ti", "O"]`
        numIons: a list of the number of each type of atom within the
            primitive cell (NOT the conventional cell), e.g., `[4, 2]`
        factor (optional): volume factor used to generate the crystal
        area: the effective cross-sectional area (A^2), of the unit cell
        sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
        lattice (optional): `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_
            object to define the unit cell
        tm (optional): `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_
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
        conventional = True,
        tm=Tol_matrix(prototype="atomic"),
    ):
        self.dim = 1
        self.PBC = [0, 0, 1]
        self.area = area  # the effective cross-sectional area, in A^2, of the unit cell.
        self.init_common(species, numIons, factor, group, lattice, sites, conventional, tm)


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
        lattice (optional): `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_
            object to define the unit cell
        tm (optional): `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.Tol_matrix.html>`_
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
        self.dim = 0
        self.PBC = [0, 0, 0]
        self.init_common(species, numIons, factor, group, lattice, sites, False, tm)
