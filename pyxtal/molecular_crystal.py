"""
Module for generating molecular crystals
"""

# Standard Libraries
from copy import deepcopy
import numpy as np

from pyxtal.lattice import Lattice
from pyxtal.molecule import pyxtal_molecule
from pyxtal.msg import Comp_CompatibilityError, Symm_CompatibilityError, VolumeError
from pyxtal.symmetry import Group
from pyxtal.symmetry import choose_wyckoff_mol as wyc_mol
from pyxtal.tolerance import Tol_matrix
from pyxtal.wyckoff_site import mol_site


# Define functions
# ------------------------------
class molecular_crystal:
    """
    Class for storing and generating molecular crystals based on symmetry
    constraints. Based on the crystal.random_crystal class for atomic crystals.
    Given a spacegroup, list of molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints.

    Args:
        dim: dimenion (1, 2, 3)
        group: the group number (1-75, 1-80, 1-230)
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or the name of molecules from the built_in
            `database <pyxtal.database.collection.html>`_
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        factor: A volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between molecules
        lattice (optional): the `Lattice <pyxtal.lattice.Lattice.html>`_
            object to define the unit cell
        conventional (optional): count the atomic numbers in a conventional cell
        tm (optional): the `Tol_matrix <pyxtal.tolerance.tolerance.html>`_
            object to define the distances
        sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
        seed (optional): seeds
        use_hall: False
    """

    def __init__(
        self,
        dim,
        group,
        molecules,
        numMols,
        factor=1.1,
        thickness=None,
        area=None,
        lattice=None,
        torsions=None,
        tm=Tol_matrix(prototype="molecular"),
        sites=None,
        conventional=True,
        seed=None,
        random_state=None,
        use_hall=False,
    ):
        # Initialize
        self.source = "Random"
        self.valid = False
        self.factor = factor
        self.seed = seed
        self.random_state = np.random.default_rng(random_state)

        # Dimesion
        self.dim = dim
        self.area = area  # Cross-section area for 1D
        self.thickness = thickness  # Thickness of 2D slab

        # The periodic boundary condition
        if dim == 3:
            self.PBC = [1, 1, 1]
        elif dim == 2:
            self.PBC = [1, 1, 0]
        elif dim == 1:
            self.PBC = [0, 0, 1]

        # Symmetry Group
        if type(group) == Group:
            self.group = group
        else:
            self.group = Group(group, dim=self.dim, use_hall=use_hall)
        self.number = self.group.number
        self.hall_number = self.group.hall_number

        # Composition
        if numMols is None:
            numMols = [len(self.group[0])] * len(molecules)
            no_check_compability = True
        else:
            numMols = np.array(numMols)  # must convert it to np.array
            no_check_compability = False
        mul = self.group.cellsize() if not conventional else 1
        self.numMols = numMols * mul

        # Tolerance matrix
        if type(tm) == Tol_matrix:
            self.tol_matrix = tm
        else:
            self.tol_matrix = Tol_matrix(prototype=tm)

        # Wyckofff sites
        self.set_molecules(molecules, torsions)
        self.set_sites(sites)
        valid_orientation = self.set_orientations()

        if not valid_orientation:
            msg = "Molecular symmetry is compatible with WP site\n"
            msg += "".join(f"{mol}: {mol.pga.sch_symbol}, " for mol in self.molecules)
            raise Symm_CompatibilityError(msg)

        # Check the minimum dof within the Wyckoff positions
        if no_check_compability:
            compat, self.degrees = True, True
        else:
            compat, self.degrees = self.group.check_compatible(self.numMols, self.valid_orientations)

        if not compat:
            msg = f"Inincompatible compoisition {self.numMols} with symmetry {self.group.number}"
            raise Comp_CompatibilityError(msg)

        self.set_volume()
        self.set_lattice(lattice)
        self.set_crystal()

    def __str__(self):
        s = "------Random Molecular Crystal------"
        s += f"\nDimension: {self.dim}"
        # s += f"\nGroup: {self.symbol}"
        s += f"\nVolume factor: {self.factor}"
        s += f"\n{self.lattice}"
        if self.valid:
            s += "\nWyckoff sites:"
            s += "".join(f"\n\t{wyc}" for wyc in self.mol_sites)
        else:
            s += "\nStructure not generated."
        return s

    def __repr__(self):
        return str(self)

    def set_sites(self, sites):
        """
        Initialize and store symmetry sites for each molecule.

        This function processes a list of symmetry sites, validates them against the expected
        number of molecules, and stores them in the `self.sites` dictionary. Each entry in the
        `sites` list can be either a dictionary or another type (list, etc.), and if no valid
        site is provided for a molecule, `None` is assigned for that molecule's site.

        Args:
            sites (list): A list of sites corresponding to `self.molecules`. They can be:
                      - A dictionary of site information (keys represent Wyckoff letters or
                        other identifiers, and values are the corresponding information).
                      - A list or other type representing site information.
                      - None, if no symmetry site information is available for that molecule.
        """
        # Initialize the self.sites dictionary to store site information
        self.sites = {}

        # Iterate over the molecules and their corresponding sites
        for i, _mol in enumerate(self.molecules):
            if sites is not None and sites[i] is not None and len(sites[i]) > 0:
                self._check_consistency(sites[i], self.numMols[i])
                if isinstance(sites[i], dict):
                    self.sites[i] = []
                    for item in sites[i].items():
                        self.sites[i].append({item[0]: item[1]})
                else:
                    self.sites[i] = sites[i]
            else:
                self.sites[i] = None

    def set_molecules(self, molecules, torsions):
        """
        Initialize and store molecular information.

        This function processes a list of molecules and initializes each one
        as a `pyxtal_molecule` object. If torsions are provided, they are
        applied to the respective molecules during initialization.
        It stored information in the `self.molecules` attribute.

        Args:
            molecules (list): A list of molecules, where each entry can either be:
                          - A SMILES string or file path representing a molecule.
                          - An already-initialized `pyxtal_molecule` object.
            torsions (list): A list of torsion angles. The length of `torsions`
                         must be equal to the length of `molecules`, or None
                         can be provided if no torsions are needed.
        """

        # If no torsions are provided, initialize with None for each molecule
        if torsions is None:
            torsions = [None] * len(molecules)

        # Initialize the molecules list to store processed pyxtal_molecule objects
        self.molecules = []

        # Iterate over the molecules
        for i, mol in enumerate(molecules):
            # already a pyxtal_molecule object
            if isinstance(mol, pyxtal_molecule):
                p_mol = mol
            else:
                p_mol = pyxtal_molecule(mol,
                                        seed=self.seed,
                                        torsions=torsions[i],
                                        tm=self.tol_matrix,
                                        random_state=self.random_state)
            self.molecules.append(p_mol)

    def set_orientations(self):
        """
        Calculates the valid orientations for each Molecule and Wyckoff
        position. Returns a list with 4 indices:
            - index 1: the molecular prototype's index within self.molecules
            - index 2: the WP's 1st index (based on multiplicity)
            - index 3: the WP's 2nd index (within the group of same multiplicity)
            - index 4: the index of a valid orientation for the molecule/WP pair

        For example, self.valid_orientations[i][j][k] would be a list of valid
        orientations for self.molecules[i], in the Wyckoff position
        self.group.wyckoffs_organized[j][k]
        """
        valid_ori = False
        self.valid_orientations = []
        for i, pyxtal_mol in enumerate(self.molecules):
            self.valid_orientations.append([])
            for x in self.group.wyckoffs_organized:
                self.valid_orientations[-1].append([])
                for _j, wp in enumerate(x):
                    # Don't check the wp with high multiplicity
                    if len(wp) > self.numMols[i]:
                        allowed = []
                    else:
                        allowed = pyxtal_mol.get_orientations_in_wp(wp)
                        if len(allowed) > 0:
                            valid_ori = True
                    self.valid_orientations[-1][-1].append(allowed)
        return valid_ori

    def set_volume(self):
        """
        Given the molecular stoichiometry, estimate the volume for a unit cell.
        """
        volume = 0
        for numMol, mol in zip(self.numMols, self.molecules):
            volume += numMol * mol.volume
        self.volume = abs(self.factor * volume)

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
                raise ValueError("PBC is incompatible " + str(self.PBC))
        else:
            # Determine the unique axis
            if self.dim == 2:
                unique_axis = "c" if self.number in range(3, 8) else "a"
            elif self.dim == 1:
                unique_axis = "a" if self.number in range(3, 8) else "c"
            else:
                unique_axis = "c"

            # Generate a Lattice instance
            good_lattice = False
            for _cycle in range(10):
                try:
                    if self.group.number < 10:
                        coef = 1.0 * self.numMols[0] / self.group[0].multiplicity
                    elif 10 <= self.group.number <= 15:
                        coef = 2.0  # 2/m
                    elif 16 <= self.group.number <= 74:
                        coef = 1.5
                    else:
                        coef = 1.0

                    self.lattice = Lattice(
                        self.group.lattice_type,
                        self.volume,
                        PBC=self.PBC,
                        unique_axis=unique_axis,
                        thickness=self.thickness,
                        area=self.area,
                        min_special=coef * max([mol.get_max_length() for mol in self.molecules]),
                        random_state=self.random_state,
                    )
                    good_lattice = True
                    break
                except VolumeError:
                    self.volume *= 1.1
                    msg = "Warning: increase the volume by 1.1 times: "
                    msg += f"{self.volume:.2f}"
                    print(msg)

            if not good_lattice:
                msg = f"Volume estimation {self.volume:.2f} is very bad"
                msg += " with the given composition "
                msg += str(self.numMols)
                raise RuntimeError(msg)

    def set_crystal(self):
        """
        The main code to generate a random molecular crystal.
        If successful, `self.valid` is True
        """
        self.numattempts = 0
        if not self.degrees:
            self.lattice_attempts = 20
            self.coord_attempts = 3
            self.ori_attempts = 1
        else:
            self.lattice_attempts = 40
            self.coord_attempts = 30
            self.ori_attempts = 5

        if not self.lattice.allow_volume_reset:
            self.lattice_attempts = 1

        for cycle1 in range(self.lattice_attempts):
            self.cycle1 = cycle1
            for cycle2 in range(self.coord_attempts):
                self.cycle2 = cycle2
                output = self._set_coords()

                if output:
                    self.mol_sites = output
                    break
            if self.valid:
                return
            else:
                self.lattice.reset_matrix()

        print("Cannot generate crystal after max attempts.")

    def _set_coords(self):
        """
        generate coordinates for random crystal
        """

        mol_sites_total = []
        # Add molecules
        for i, numMol in enumerate(self.numMols):
            pyxtal_mol = self.molecules[i]
            valid_ori = self.valid_orientations[i]
            output = self._set_mol_wyckoffs(i, numMol, pyxtal_mol, valid_ori, mol_sites_total)
            if output is not None:
                mol_sites_total.extend(output)
            else:
                # correct multiplicity not achieved exit and start over
                return None

        self.valid = True
        return mol_sites_total

    def _set_mol_wyckoffs(self, id, numMol, pyxtal_mol, valid_ori, mol_wyks):
        """
        generates a set of wyckoff positions to accomodate a given number
        of molecules

        Args:
            id: molecular id
            numMol: Number of ions to accomodate
            pyxtal_mol: Type of species being placed on wyckoff site
            valid_ori: list of valid orientations
            mol_wyks: current wyckoff sites

        Returns:
            if sucess, wyckoff_sites_tmp: list of wyckoff sites for valid sites
            otherwise, None

        """
        numMol_added = 0
        mol_sites_tmp = []

        # Now we start to add the specie to the wyckoff position
        sites_list = deepcopy(self.sites[id])  # the list of Wyckoff site
        if sites_list is not None:
            self.wyckoff_attempts = max(len(sites_list) * 2, 10)
        else:
            # the minimum numattempts is to put all atoms to the general WPs
            min_wyckoffs = int(numMol / len(self.group.wyckoffs_organized[0][0]))
            self.wyckoff_attempts = max(2 * min_wyckoffs, 10)

        for _cycle in range(self.wyckoff_attempts):
            # Choose a random WP for given multiplicity: 2a, 2b, 2c
            if sites_list is not None and len(sites_list) > 0:
                site = sites_list[0]
            else:  # Selecting the merging
                site = None

            # NOTE: The molecular version return wyckoff indices, not ops
            diff = numMol - numMol_added

            if type(site) is dict:  # site with coordinates
                key = next(iter(site.keys()))
                wp = wyc_mol(self.group, diff, key, valid_ori, True, self.dim, self.random_state)
            else:
                wp = wyc_mol(self.group, diff, site, valid_ori, True, self.dim, self.random_state)

            if wp is not False:
                # Generate a list of coords from the wyckoff position
                mult = wp.multiplicity  # remember the original multiplicity

                pt = site[key] if type(site) is dict else self.lattice.generate_point()
                # merge coordinates if the atoms are close
                mtol = pyxtal_mol.radius * 0.5
                pt, wp, oris = wp.merge(pt, self.lattice.matrix, mtol, valid_ori, self.group)

                if wp is not False:
                    if site is not None and mult != wp.multiplicity:
                        continue
                    if self.dim == 2 and self.thickness is not None and self.thickness < 0.1:
                        pt[-1] = 0.5

                    ms0 = self._set_orientation(pyxtal_mol, pt, oris, wp)
                    if not ms0.short_dist():
                        # Check current WP against existing WP's
                        passed_wp_check = True
                        # print("Checking", ms0)
                        for ms1 in mol_wyks + mol_sites_tmp:
                            if not ms0.short_dist_with_wp2(ms1, tm=self.tol_matrix):
                                passed_wp_check = False
                            # else:
                            #    print("passing", ms1)

                        if passed_wp_check:
                            if sites_list is not None:
                                sites_list.pop(0)

                            ms0.type = id
                            mol_sites_tmp.append(ms0)
                            numMol_added += len(ms0.wp)

                            # We have enough molecules of the current type
                            if numMol_added == numMol:
                                return mol_sites_tmp
        return None


    def _set_orientation(self, pyxtal_mol, pt, oris, wp):
        """
        Generate valid orientations for a given molecule in a Wyckoff position.

        It tries to generate valid orientations for the molecule by:
        - Selecting a random orientation from a list of possible orientations.
        - Flipping the orientation to test different alignments.
        - Checking the smallest distance between atoms and ensuring it's valid.
        - Using the bisection method is refine the orientation.

        Args:
            pyxtal_mol: The pyxtal_molecule object representing the molecule.
            pt: Position of the molecule.
            oris: List of potential orientations.
            wp: Wyckoff position object representing the symmetry of the site.

        Returns:
            ms0: A valid `mol_site` object if an acceptable orientation is found.
             returns `None` if no valid orientation is found within the attempts.
        """

        # NOTE removing this copy causes tests to fail -> state not managed well
        ori = self.random_state.choice(oris).copy()
        ori.change_orientation(flip=True)

        # Create a mol_site object with the current orientation
        ms0 = mol_site(pyxtal_mol, pt, ori, wp, self.lattice)

        # Check if the current orientation results in valid distances
        if len(pyxtal_mol.mol) > 1 and ori.degrees > 0:
            if ms0.short_dist():
                ms0.optimize_orientation_by_dist(self.ori_attempts)
        return ms0

    def _check_consistency(self, site, numMol):
        """
        Check if N_mol is consistent with the symmetry constraints of the system.

        It verifies if the sum of molecules from the WP matches (`numMol`).
        Each Wyckoff site string in the `site` list includes a number that
        represents how many molecules are associated with that site.

        If a inconsistency is found, it raises a ValueError with a detailed message.

        Args:
            site (list of str): A list of strings for Wyckoff sites. Each string
                            contains a number (e.g., "3a", "4b") where the number
                            indicates how many molecules are at that site.
            numMol (int): The total number of molecules expected in the structure.

        Returns:
            bool: Returns `True` if the number of molecules matches `numMol`.
                Raises a ValueError if they do not match.
        """
        num = 0
        for s in site:
            num += int(s[:-1])
        if numMol == num:
            return True
        else:
            msg = "\nThe requested number of molecules is inconsistent: "
            msg += str(site)
            msg += f"\nfrom numMols: {numMol:d}"
            msg += f"\nfrom Wyckoff list: {num:d}"
            raise ValueError(msg)
