# Standard Libraries
import os
import random
import numpy as np
from copy import deepcopy

# External Libraries

# PyXtal imports
from pyxtal.constants import max1, max2, max3, max4
from pyxtal.msg import printx
from pyxtal.tolerance import Tol_matrix
from pyxtal.lattice import Lattice, cellsize
from pyxtal.io import write_cif, structure_from_ext
from pyxtal.database.element import Element
from pyxtal.wyckoff_site import mol_site, check_mol_sites, WP_merge
from pyxtal.molecule import pyxtal_molecule, orientation_in_wyckoff_position
from pyxtal.symmetry import Group, jk_from_i, choose_wyckoff_molecular



class molecular_crystal:
    """
    Class for storing and generating molecular crystals based on symmetry
    constraints. Based on the crystal.random_crystal class for atomic crystals.
    Given a spacegroup, list of molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. 

    Args:
        group: the spacegroup number (1-230), or a 
            `pyxtal.symmetry.Group <pyxtal.symmetry.Group.html>`_ object
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or the name of molecules from the built_in 
            `database <pyxtal.database.collection.html>`_
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        volume_factor: A volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between molecules
        allow_inversion: Whether or not to allow chiral molecules to be
            inverted. If True, the final crystal may contain mirror images of
            the original molecule. Unless the chemical properties of the mirror
            image are known, it is highly recommended to keep this value False
        orientations: Once a crystal with the same spacegroup and molecular
            stoichiometry has been generated, you may pass its
            valid_orientations attribute here to avoid repeating the
            calculation, but this is not required
        lattice (optional): the `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_ 
            object to define the unit cell
        tm (optional): the `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.tolerance.html>`_ 
            object to define the distances
        seed (optional): the cif file from user
        diag (optional): if use the nonstandart setting (P21/n, Pn, C2/n)?
    """

    def __init__(
        self,
        group,
        molecules,
        numMols,
        volume_factor=1.1,
        select_high=True,
        allow_inversion=True,
        orientations=None,
        lattice=None,
        tm=Tol_matrix(prototype="molecular"),
        seed = None,
        diag = False,
    ):

        self.dim = 3 # The number of periodic dimensions (1,2,3)
        self.PBC = [1, 1, 1]
        self.diag = diag

        if type(group) != Group:
            group = Group(group, self.dim)

        self.sg = group.number
        self.selec_high = select_high
        self.seed = seed

        self.init_common(
            molecules,
            numMols,
            volume_factor,
            select_high,
            allow_inversion,
            orientations,
            group,
            lattice,
            tm,
        )

    def init_common(
        self,
        molecules,
        numMols,
        volume_factor,
        select_high,
        allow_inversion,
        orientations,
        group,
        lattice,
        tm,
    ):
        # init functionality which is shared by 3D, 2D, and 1D crystals
        self.numattempts = 0 # number of attempts to generate the crystal.
        if type(group) == Group:
            self.group = group
            """A pyxtal.symmetry.Group object storing information about the space/layer
            /Rod/point group, and its Wyckoff positions."""
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
        self.Msgs()
        self.factor = volume_factor  # volume factor for the unit cell.
        numMols = np.array(numMols)  # must convert it to np.array
        self.numMols0 = numMols  # in the PRIMITIVE cell
        self.numMols = self.numMols0 * cellsize(self.group)  # in the CONVENTIONAL cell

        # boolean numbers
        self.allow_inversion = allow_inversion
        self.select_high = select_high

        # Set the tolerance matrix
        # The Tol_matrix object for checking inter-atomic distances within the structure.
        if type(tm) == Tol_matrix:
            self.tol_matrix = tm
        else:
            try:
                self.tol_matrix = Tol_matrix(prototype=tm)
            # TODO remove bare except
            except:
                msg = "Error: tm must either be a Tol_matrix object +\n"
                msg += "or a prototype string for initializing one."
                printx(msg, priority=1)
                self.valid = False
                return

        self.molecules = []  # A pyxtal_molecule objects,
        for mol in molecules:
            self.molecules.append(pyxtal_molecule(mol, self.tol_matrix))

        # if seeds, directly parse the structure from cif
        # At the moment, we only support one specie
        if self.seed is not None:
            seed = structure_from_ext(self.seed, self.molecules[0].mol)
            if seed.match():
                self.mol_sites = [seed.make_mol_site()]
                self.group = Group(seed.wyc.number)
                self.lattice = seed.lattice
                self.molecules = [pyxtal_molecule(seed.molecule)]
                self.diag = seed.diag
                self.valid = True # Need to add a check function
            else:
                raise ValueError("Cannot extract the structure from cif")
        """
        The valid orientations for each molecule and Wyckoff position.
        May be copied when generating a new molecular_crystal to save a
        small amount of time
        """

        if orientations is None:
            self.get_orientations()
        else:
            self.valid_orientations = orientations

        if self.seed is None:
            if lattice is not None:
                # Use the provided lattice
                self.lattice = lattice
                self.volume = lattice.volume
                # Make sure the custom lattice PBC axes are correct.
                if lattice.PBC != self.PBC:
                    self.lattice.PBC = self.PBC
                    printx("\n  Warning: converting custom lattice PBC to " + str(self.PBC))
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
                self.volume = self.estimate_volume()

                # Calculate the minimum, middle, and maximum box lengths for
                # the unit cell. make sure at least one non-overlapping orientation
                # exists for each molecule
                minls = []
                midls = []
                maxls = []
                for mol in self.molecules:
                    minls.append(mol.box.minl)
                    midls.append(mol.box.midl)
                    maxls.append(mol.box.maxl)

                # The Lattice object used to generate lattice matrices
                if self.dim == 3 or self.dim == 0:
                    self.lattice = Lattice(
                        self.group.lattice_type,
                        self.volume,
                        PBC=self.PBC,
                        unique_axis=unique_axis,
                        min_l=max(minls),
                        mid_l=max(midls),
                        max_l=max(maxls),
                    )
                elif self.dim == 2:
                    self.lattice = Lattice(
                        self.group.lattice_type,
                        self.volume,
                        PBC=self.PBC,
                        unique_axis=unique_axis,
                        min_l=max(minls),
                        mid_l=max(midls),
                        max_l=max(maxls),
                        thickness=self.thickness,
                    )
                elif self.dim == 1:
                    self.lattice = Lattice(
                        self.group.lattice_type,
                        self.volume,
                        PBC=self.PBC,
                        unique_axis=unique_axis,
                        min_l=max(minls),
                        mid_l=max(midls),
                        max_l=max(maxls),
                        area=self.area,
                    )


            self.generate_crystal()

    def estimate_volume(self):
        """
        Given the molecular stoichiometry, estimate the volume needed for a unit cell.

        Returns:
            the estimated volume (in cubic Angstroms) needed for the unit cell
        """
        volume = 0
        for numMol, mol in zip(self.numMols, self.molecules):
            volume += numMol * mol.volume
        return abs(self.factor * volume)

    def check_short_distances(self, r=1.0, exclude_H = True):
        """
        A function to check short distance pairs
        Mainly used for debug, powered by pymatgen

        Args:
            r: the given cutoff distances
            exclude_H: whether or not ignore the H atoms

        Returns:
            a list of pairs within the cutoff
        """
        pairs = []
        pmg_struc = self.to_pymatgen()
        if exclude_H:
            pmg_struc.remove_species('H')
        res = pmg_struc.get_all_neighbors(r)
        for i, neighs in enumerate(res):
            for n in neighs:
                pairs.append([pmg_struc.sites[i].specie, n.specie, n.nn_distance])
        return pairs


    def Msgs(self):
        self.Msg1 = (
            "Error: the stoichiometry is incompatible with the wyckoff sites choice"
        )
        self.Msg2 = "Error: failed in the cycle of generating structures"
        self.Msg3 = "Warning: failed in the cycle of adding species"
        self.Msg4 = "Warning: failed in the cycle of choosing wyckoff sites"
        self.Msg5 = "Finishing: added the specie"
        self.Msg6 = "Finishing: added the whole structure"
        self.Msg7 = "Error: invalid paramaters for initialization"

    def get_orientations(self):
        """
        Calculates the valid orientations for each Molecule and Wyckoff
        position. Returns a list with 4 indices:
            - index 1: the molecular prototype's index within self.molecules
            - index 2: the Wyckoff position's 1st index (based on multiplicity)
            - index 3: the WP's 2nd index (within the group of equal multiplicity)
            - index 4: the index of the valid orientation for the molecule/WP pair

        For example, self.valid_orientations[i][j][k] would be a list of valid
        orientations for self.molecules[i], in the Wyckoff position
        self.group.wyckoffs_organized[j][k]
        """
        self.valid_orientations = []
        for pyxtal_mol in self.molecules:
            self.valid_orientations.append([])
            wp_index = -1
            for i, x in enumerate(self.group.wyckoffs_organized):
                self.valid_orientations[-1].append([])
                for j, wp in enumerate(x):
                    wp_index += 1
                    allowed = orientation_in_wyckoff_position(
                        pyxtal_mol.mol,
                        wp,
                        already_oriented=True,
                        allow_inversion=self.allow_inversion,
                    )

                    if allowed is not False:
                        self.valid_orientations[-1][-1].append(allowed)
                    else:
                        self.valid_orientations[-1][-1].append([])

    def check_compatible(self, group, numMols, valid_orientations):
        """
        Checks if the number of molecules is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        # Store whether or not at least one degree of freedom exists
        has_freedom = False
        # Store the wp's already used that don't have any freedom
        used_indices = []
        # Loop over species
        for i_mol, numIon in enumerate(numMols):
            # Get lists of multiplicity, maxn and freedom
            l_mult0 = []
            l_maxn0 = []
            l_free0 = []
            indices0 = []
            for i_wp, wp in enumerate(group):
                # Check that at least one valid orientation exists
                j, k = jk_from_i(i_wp, group.wyckoffs_organized)
                if len(valid_orientations[i_mol][j][k]) > 0:
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
                    l_maxn.append(1)
                    l_free.append(False)
                    indices.append(i_wp)
            # Loop over possible combinations
            # Create pointer variable to move through lists
            p = 0
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
                print("n == n0", n, n0)
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
                                indices.append(i_wp)
                    break
                # All combinations failed: return False
                if n == n0 and p >= len(l_mult) - 1:
                    print("All combinations failed: return False")
                    return False
                # Too few atoms
                if num < numIon:
                    # Forwards routine
                    # Move p to the right and max out
                    if p < len(l_mult) - 1:
                        p += 1
                        n[p] = min((numIon - num) // l_mult[p], l_maxn[p])
                    else:
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

    def to_file(self, filename=None, fmt="cif", permission='w', **kwargs):
        """
        Creates a file with the given filename and file type to store the structure.
        By default, creates cif files for crystals and xyz files for clusters.
        By default, the filename is based on the stoichiometry.

        Args:
            filename: the file path
            fmt: the file type ('cif', 'xyz', etc.)
            permission: "w" or "a+"
            sym_num: "w" or "a+"

        Returns:
            Nothing. Creates a file at the specified path
        """
        if self.valid:
            if fmt == "cif":
                return write_cif(self, filename, "from_pyxtal", permission, **kwargs)
            else:
                pmg_struc = self.to_pymatgen()
                pmg_struc.sort()
                return pmg_struc.to(fmt=fmt, filename=filename)
        else:
            printx("Cannot create file: structure did not generate.", priority=1)

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
            total_coords: N*3 numpy array 
            species: N-length list, e.g. ["C", "C", ...]
        """

        species = []
        total_coords = None
        for site in self.mol_sites:
            coords, site_species = site.get_coords_and_species(absolute)
            species.extend(site_species)
            if total_coords is None:
                total_coords = coords
            else:
                total_coords = np.append(total_coords, coords, axis=0)

        return total_coords, species


    def to_ase(self, resort=True):
        """
        export to ase Atoms object
        """
        from ase import Atoms
        if self.valid:
            coords, species = self._get_coords_and_species(True)
            latt, coords = self.lattice.add_vacuum(coords, PBC=self.PBC)
            atoms = Atoms(species, positions=coords, cell=latt)
            if resort:
                permutation = np.argsort(atoms.numbers)
                atoms = atoms[permutation]
            return atoms
        else:
            printx("No valid structure can be converted to ase.", priority=1)

    def to_pymatgen(self):
        """
        export to Pymatgen structure object
        """
        from pymatgen.core.structure import Structure  

        if self.valid:
            coords, species = self._get_coords_and_species()
            # Add space above and below a 2D or 1D crystals
            latt, coords = self.lattice.add_vacuum(coords, PBC=self.PBC)
            return Structure(latt, species, coords)
        else:
            printx("No valid structure can be converted to pymatgen.", priority=1)


    def __str__(self):
        s = "------Random Molecular Crystal------"
        s += "\nDimension: " + str(self.dim)
        if self.diag and hasattr(self.group, 'alias'):
            symbol = self.group.alias
        else:
            symbol = self.group.symbol
        s += "\nGroup: " + symbol
        s += "\nVolume factor: " + str(self.factor)
        s += "\n" + str(self.lattice)
        if self.valid:
            s += "\nWyckoff sites:"
            for wyc in self.mol_sites:
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
        from pyxtal.viz import display_molecular
        return display_molecular(self, **kwargs)

    def generate_crystal(self, max1=max1, max2=max2, max3=max3, max4=max4):
        """
        The main code to generate a random molecular crystal. If successful,
        `self.valid` is True (False otherwise) 

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
            max4: the number of attempts for changing the molecular orientation
        """
        # Check the minimum number of degrees of freedom within the Wyckoff positions
        degrees = self.check_compatible(self.group, self.numMols, self.valid_orientations)
        if degrees is False:
            self.struct = None
            self.valid = False
            msg = "the space group is incompatible with the number of molecules"
            raise ValueError(msg)
            # return
        else:
            if degrees == 0:
                max1 = 20
                max2 = 3
                max3 = 1
                max4 = 1
            # Calculate a minimum vector length for generating a lattice
            all_lengths = []
            for mol in self.molecules:
                all_lengths.append(mol.box.minl)

            for cycle1 in range(max1):
                self.cycle1 = cycle1
                # 1, Generate a lattice
                if self.lattice.allow_volume_reset is True:
                    self.volume = self.estimate_volume()
                    self.lattice.volume = self.volume
                self.lattice.reset_matrix()
                #cell_matrix = self.lattice.get_matrix()

                for cycle2 in range(max2):
                    mol_sites_total = []
                    good_structure = False

                    self.cycle2 = cycle2
                    mol_sites_tmp = []

                    # Add molecules 
                    for i, numMol in enumerate(self.numMols):
                        pyxtal_mol = self.molecules[i]
                        valid_ori = self.valid_orientations[i]
                        numMol_added = 0

                        # Now we start to add the specie to the wyckoff position
                        for cycle3 in range(max3):
                            self.cycle3 = cycle3
                            self.cycle4 = 0
                            self.numattempts += 1
                            # Choose a random WP for given multiplicity: 2a, 2b, 2c
                            # NOTE: The molecular version return wyckoff indices, not ops
                            wp = choose_wyckoff_molecular(
                                self.group,
                                numMol - numMol_added,
                                valid_ori,
                                self.select_high,
                            )
                            if wp is not False:
                                # Generate a list of coords from the wyckoff position
                                mult = wp.multiplicity # remember the original multiplicity
                                pt = self.lattice.generate_point()

                                # merge coordinates if the atoms are close
                                mtol = mol.radius * 0.5
                                pt, wp, oris = WP_merge(pt, self.lattice.matrix, wp, mtol, valid_ori)

                                if wp is not False:
                                    # Use a Wyckoff_site object for the current site
                                    ori = random.choice(oris).copy()
                                    ori.change_orientation()
                                    ms0 = mol_site(pyxtal_mol, pt, ori, wp, self.lattice, self.diag)
                                    # Check distances within the WP
                                    if not ms0.check_distances():
                                        # Maximize the smallest distance for the general
                                        # positions if needed
                                        passed_ori = False
                                        if len(pyxtal_mol.mol) > 1 and ori.degrees > 0:
                                            # bisection method
                                            def fun_dist(angle, ori, mo, pt):
                                                ori0 = ori.copy()
                                                ori.change_orientation(angle)
                                                ms0 = mol_site(
                                                    mo,
                                                    pt,
                                                    ori,
                                                    wp,
                                                    self.lattice,
                                                    self.diag,
                                                )
                                                d = ms0.compute_distances()
                                                return d

                                            angle_lo = ori.angle
                                            angle_hi = angle_lo + np.pi
                                            fun_lo = fun_dist(angle_lo, ori, pyxtal_mol, pt)
                                            fun_hi = fun_dist(angle_hi, ori, pyxtal_mol, pt)
                                            fun = fun_hi
                                            for it in range(max4):
                                                if (fun > 0.8) & (ms0.check_distances()):
                                                    passed_ori = True
                                                    break
                                                angle = (angle_lo + angle_hi) / 2
                                                fun = fun_dist(angle, ori, pyxtal_mol, pt)
                                                #print('Bisection: ', it, fun)
                                                if fun_lo > fun_hi:
                                                    angle_hi, fun_hi = angle, fun
                                                else:
                                                    angle_lo, fun_lo = angle, fun
                                    else:
                                        passed_ori = True
                                    if not passed_ori:
                                        continue

                                    # Check distances with other WP's
                                    passed = True
                                    for ms1 in mol_sites_tmp:
                                        if not check_mol_sites(ms0, ms1, tm=self.tol_matrix):
                                            passed = False
                                            break

                                    if not passed:
                                        continue
                                    else:
                                        # Distance checks passed; store the new Wyckoff position
                                        mol_sites_tmp.append(ms0)
                                        numMol_added += len(ms0.wp)
                                        if numMol_added == numMol:
                                            # We have enough molecules of the current type
                                            mol_sites_total = deepcopy(mol_sites_tmp)
                                            break

                        if numMol_added != numMol:
                            break  # need to repeat from the 1st species

                    if numMol_added == numMol:
                        printx(self.Msg6, priority=3)
                        good_structure = True
                        break
                # placing molecules here
                if good_structure:
                    self.mol_sites = mol_sites_total
                    self.valid = True
                    return
                else:
                    printx("Failed final distance check.", priority=3)

        printx("Couldn't generate crystal after max attempts.", priority=1)
        if degrees == 0:
            printx("Note: Wyckoff positions have no degrees of freedom.", priority=2)
        self.valid = False

        return


class molecular_crystal_2D(molecular_crystal):
    """
    A 2d counterpart to molecular_crystal. Given a layer group, list of
    molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. This crystal is stored as a pymatgen struct via self.struct

    Args:
        group: the layer group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
            OR, a pyxtal.symmetry.Group object
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or the name of molecules from the built_in 
            `database <pyxtal.database.collection.html>`_
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        thickness: the thickness, in Angstroms, of the unit cell in the 3rd
            dimension (the direction which is not repeated periodically). A
            value of None causes a thickness to be chosen automatically. Note
            that this constraint applies only to the molecular centers; some
            atomic coordinates may lie outside of this range
        volume_factor: A volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between molecules
        allow_inversion: Whether or not to allow chiral molecules to be
            inverted. If True, the final crystal may contain mirror images of
            the original molecule. Unless the chemical properties of the mirror
            image are known, it is highly recommended to keep this value False
        orientations: Once a crystal with the same spacegroup and molecular
            stoichiometry has been generated, you may pass its
            valid_orientations attribute here to avoid repeating the
            calculation, but this is not required
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        molecules,
        numMols,
        volume_factor=1.1,
        select_high=True,
        allow_inversion=True,
        orientations=None,
        thickness=None,
        lattice=None,
        tm=Tol_matrix(prototype="molecular"),
        seed = None,
    ):

        self.dim = 2
        self.numattempts = 0
        self.seed = None
        if type(group) != Group:
            group = Group(group, self.dim)
        number = group.number  # The layer group number of the crystal."""
        self.diag = False
        self.thickness = thickness  # the thickness in Angstroms
        self.PBC = [1, 1, 0]
        self.init_common(
            molecules,
            numMols,
            volume_factor,
            select_high,
            allow_inversion,
            orientations,
            group,
            lattice,
            tm,
        )


class molecular_crystal_1D(molecular_crystal):
    """
    A 1d counterpart to molecular_crystal. Given a Rod group, list of
    molecule objects, molecular stoichiometry, volume factor, and area,
    generates a molecular crystal consistent with the given constraints.
    The crystal is stored as a pymatgen struct via self.struct

    Args:
        group: the Rod group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
            OR, a pyxtal.symmetry.Group object
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or the name of molecules from the built_in 
            `database <pyxtal.database.collection.html>`_
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        area: cross-sectional area of the unit cell in Angstroms squared. A
            value of None causes an area to be chosen automatically. Note that
            this constraint applies only to the molecular centers; some atomic
            coordinates may lie outside of this range
        volume_factor: A volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between molecules
        allow_inversion: Whether or not to allow chiral molecules to be
            inverted. If True, the final crystal may contain mirror images of
            the original molecule. Unless the chemical properties of the mirror
            image are known, it is highly recommended to keep this value False
        orientations: Once a crystal with the same spacegroup and molecular
            stoichiometry has been generated, you may pass its
            valid_orientations attribute here to avoid repeating the
            calculation, but this is not required
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        molecules,
        numMols,
        volume_factor=1.1,
        select_high=True,
        allow_inversion=False,
        orientations=None,
        area=None,
        lattice=None,
        tm=Tol_matrix(prototype="molecular"),
    ):
        self.dim = 1
        self.area = area  # the effective cross-sectional area in A^2
        self.diag = False
        self.PBC = [0, 0, 1]  # The periodic axes of the crystal (1,2,3)->(x,y,z)
        self.sg = None  # The international space group number, not rod groups
        self.seed = None
        self.init_common(
            molecules,
            numMols,
            volume_factor,
            select_high,
            allow_inversion,
            orientations,
            group,
            lattice,
            tm
        )
