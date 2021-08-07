"""
Module for generating atomic crystals
"""

# Standard Libraries
import random
from copy import deepcopy
import numpy as np

# PyXtal imports #avoid *
from pyxtal.symmetry import Group, choose_wyckoff
from pyxtal.wyckoff_site import atom_site
from pyxtal.msg import Comp_CompatibilityError
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
        dim: dimenion (0, 1, 2, 3)
        group: the group number (1-56, 1-75, 1-80, 1-230)
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
        dim = 3,
        group = 227,
        species = ['C'],
        numIons = 8,
        factor = 1.1,
        thickness = None,
        area = None,
        lattice = None,
        sites = None,
        conventional = True,
        tm=Tol_matrix(prototype="atomic"),
    ):

        # Initialize
        self.source = 'Random'
        self.valid = False
        self.factor = factor

        # Dimesion
        self.dim = dim 
        self.area = area  # Cross-section area for 1D
        self.thickness = thickness # Thickness of 2D slab 

        #The periodic boundary condition
        if dim == 3:
            self.PBC = [1, 1, 1] 
        elif dim == 2:
            self.PBC = [1, 1, 0]
        elif dim == 1:
            self.PBC = [0, 0, 1]
        elif dim == 0:
            self.PBC = [0, 0, 0]

        # Symmetry group
        if type(group) == Group:
            self.group = group
        else:
            self.group = Group(group, dim=self.dim)
        self.number = self.group.number
        self.symbol = self.group.symbol

        # Composition
        numIons = np.array(numIons)
        if not conventional:
            mul = cellsize(self.group)
        else:
            mul = 1
        self.numIons = numIons * mul
        self.species = species

        # Tolerance matrix 
        if type(tm) == Tol_matrix:
            self.tol_matrix = tm
        else:
            self.tol_matrix = Tol_matrix(prototype=tm)

        # Wyckoff sites
        self.set_sites(sites) 

        # Lattice and coordinates
        compat, self.degrees = self._check_compatible()
        if not compat:
            self.valid = False
            msg = "Compoisition " + str(self.numIons) 
            msg += " not compatible with symmetry "
            msg += str(self.group.number) 
            raise Comp_CompatibilityError(msg)
        else:
            self.set_volume()
            self.set_lattice(lattice)
            self.set_crystal()

    def __str__(self):
        if self.valid:
            s = "------Crystal from {:s}------".format(self.source)
            s += "\nDimension: {}".format(self.dim)
            s += "\nGroup: {} ({})".format(self.symbol, self.number)
            s += "\n{}".format(self.lattice)
            s += "\nWyckoff sites:"
            for wyc in self.atom_sites:
                s += "\n\t{}".format(wyc)
        else:
            s = "\nStructure not available."
        return s

    def __repr__(self):
        return str(self)

    def set_sites(self, sites):
        """
        initialize Wyckoff sites

        Args:
            sites: list
        """
        # Symmetry sites
        self.sites = {}
        for i, specie in enumerate(self.species):
            if sites is not None and sites[i] is not None:
                self._check_consistency(sites[i], self.numIons[i])
                self.sites[specie] = sites[i]
            else:
                self.sites[specie] = None
 
    def set_volume(self):
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
        self.volume = self.factor * volume

    def set_lattice(self, lattice):
        """
        Generate the initial lattice
        """
        if lattice is not None:
            # Use the provided lattice
            self.lattice = lattice
            self.volume = lattice.volume
            # Make sure the custom lattice PBC axes are correct.
            if lattice.PBC != self.PBC:
                self.lattice.PBC = self.PBC
        else:
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

            # Generate a Lattice instance
            self.lattice = Lattice(
                self.group.lattice_type,
                self.volume,
                PBC=self.PBC,
                unique_axis=unique_axis,
                thickness=self.thickness,
                area=self.area,
                )

    def set_crystal(self):
        """
        The main code to generate a random atomic crystal. 
        If successful, `self.valid` is True 
       """
        self.numattempts = 0
        if not self.degrees:
            self.lattice_attempts = 5
            self.coord_attempts = 5
        else:
            self.lattice_attempts = 40
            self.coord_attempts = 10

        if not self.lattice.allow_volume_reset:
            self.lattice_attempts = 1

        for cycle1 in range(self.lattice_attempts):
            self.cycle1 = cycle1
            for cycle2 in range(self.coord_attempts):
                self.cycle2 = cycle2
                output = self._set_coords()

                if output:
                    self.atom_sites = output
                    break
            if self.valid:
                return
            else:
                self.lattice.reset_matrix()

        return

    def _set_coords(self):
        """
        generate coordinates for random crystal
        """
        wyks = []
        cell = self.lattice.get_matrix() 
        # generate coordinates for each ion type in turn
        for numIon, specie in zip(self.numIons, self.species):
            output = self._set_ion_wyckoffs(numIon, specie, cell, wyks)
            if output is not None:
                wyks.extend(output)
            else:
                # correct multiplicity not achieved exit and start over
                return None

        # If numIon_added correct for all specie return structure
        self.valid = True
        return wyks

    def _set_ion_wyckoffs(self, numIon, specie, cell, wyks):
        """
        generates a set of wyckoff positions to accomodate a given number
        of ions

        Args:
            numIon: Number of ions to accomodate
            specie: Type of species being placed on wyckoff site
            cellx: Matrix of lattice vectors
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
                pt, wp, _ = wp.merge(pt, cell, tol)
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
                        if not new_site.check_with_ws2(ws, cell, tol_matrix):
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

    def _check_compatible(self):
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
        for numIon in self.numIons:
            # Get lists of multiplicity, maxn and freedom
            l_mult0 = []
            l_maxn0 = []
            l_free0 = []
            indices0 = []
            for i_wp, wp in enumerate(self.group):
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
                return False, False
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
            return True, True
        else:
            # All species passed, but no degrees of freedom: return 0
            return True, False

    def _check_consistency(self, site, numIon):
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

