"""
Module for generating molecular crystals
"""

# Standard Libraries
import random
from copy import deepcopy
import numpy as np

# PyXtal imports
from pyxtal.msg import printx
from pyxtal.tolerance import Tol_matrix
from pyxtal.lattice import Lattice, cellsize
from pyxtal.wyckoff_site import mol_site
from pyxtal.molecule import pyxtal_molecule
from pyxtal.symmetry import Group, jk_from_i
from pyxtal.symmetry import choose_wyckoff_molecular as wyc_mol
from pyxtal.msg import Comp_CompatibilityError, Symm_CompatibilityError

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
        lattice (optional): the `pyxtal.lattice.Lattice <pyxtal.lattice.Lattice.html>`_ 
            object to define the unit cell
        conventional (optional): count the number of atoms in the conventional cell
        tm (optional): the `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.tolerance.html>`_ 
            object to define the distances
        sites (optional): pre-assigned wyckoff sites (e.g., `[["4a"], ["2b"]]`)
        diag (optional): if use the nonstandart setting (P21/n, Pn, C2/n)?
        seed (optional): seeds
    """

    def __init__(
        self,
        dim,
        group,
        molecules,
        numMols,
        factor = 1.1,
        thickness = None,
        area = None,
        lattice = None,
        torsions = None,
        tm = Tol_matrix(prototype="molecular"),
        sites = None,
        conventional = True,
        diag = False,
        seed = None,
    ):

        # Initialize
        self.source = 'Random'
        self.valid = False
        self.factor = factor
        self.seed = seed

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

        # Symmetry Group
        if type(group) == Group:
            self.group = group
        else:
            self.group = Group(group, dim=self.dim)
        self.number = self.group.number
        self.diag = diag
        if self.diag and self.number not in [5, 7, 8, 9, 12, 13, 14, 15]:
            self.diag = False
        if self.diag and self.group.number in [7, 14, 15]:
            self.symbol = self.group.alias
        else:
            self.symbol = self.group.symbol

        # Composition
        if numMols is None:
            numMols = [len(self.group[0])] * len(molecules)
            no_check_compability = True
        else:
            numMols = np.array(numMols)  # must convert it to np.array
            no_check_compability = False
        if not conventional:
            mul = cellsize(self.group)
        else:
            mul = 1
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

        if valid_orientation:
            # Check the minimum dof within the Wyckoff positions
            if no_check_compability:
                compat, self.degrees = True, True
            else:
                compat, self.degrees = self._check_compatible()
            if not compat:
                msg = "Compoisition " + str(self.numMols) 
                msg += " not compatible with symmetry "
                msg += str(self.group.number) 
                raise Comp_CompatibilityError(msg)
            else:
                self.set_volume()
                self.set_lattice(lattice)
                self.set_crystal()
        else:
            msg = "Molecular symmetry is compatible with WP site\n"
            for mol in self.molecules:
                msg += str(mol) + ": "
                msg += mol.pga.sch_symbol
            raise Symm_CompatibilityError(msg)

    def __str__(self):
        s = "------Random Molecular Crystal------"
        s += "\nDimension: " + str(self.dim)
        s += "\nGroup: " + self.symbol
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

    def set_sites(self, sites):
        """
        initialize Wyckoff sites

        Args:
            sites: list
        """
        # Symmetry sites
        self.sites = {}
        for i, mol in enumerate(self.molecules):
            if sites is not None and sites[i] is not None:
                self._check_consistency(sites[i], self.numMols[i])
                self.sites[i] = sites[i]
            else:
                self.sites[i] = None
 
    def set_molecules(self, molecules, torsions):
        """
        Get molecular information

        Args:
            molecules: list of molecules
            torsions: list of torsions
        """
        if torsions is None:
            torsions = [None]*len(molecules)

        self.molecules = []  
        for i, mol in enumerate(molecules):
            # already a pyxtal_molecule object
            if isinstance(mol, pyxtal_molecule):
                p_mol = mol
            else:
                p_mol = pyxtal_molecule(mol, seed=self.seed, torsions=torsions[i], tm=self.tol_matrix)
            self.molecules.append(p_mol)
 
    def set_orientations(self):
        """
        Calculates the valid orientations for each Molecule and Wyckoff
        position. Returns a list with 4 indices:
            - index 1: the molecular prototype's index within self.molecules
            - index 2: the WP's 1st index (based on multiplicity)
            - index 3: the WP's 2nd index (within the group of equal multiplicity)
            - index 4: the index of the valid orientation for the molecule/WP pair

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
                for j, wp in enumerate(x):
                    #Don't check the wp with high multiplicity
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
        Given the molecular stoichiometry, estimate the volume needed for a unit cell.
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

        printx("Cannot generate crystal after max attempts.", priority=1)

    def _set_coords(self):
        """
        generate coordinates for random crystal
        """

        mol_sites_total = []
        # Add molecules 
        for i, numMol in enumerate(self.numMols):
            pyxtal_mol = self.molecules[i]
            valid_ori = self.valid_orientations[i]
            output = self._set_mol_wyckoffs(
                i, numMol, pyxtal_mol, valid_ori, mol_sites_total
            )
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
            numMol: Number of ions to accomodate
            pyxtal_mol: Type of species being placed on wyckoff site
            mol_wyks: current wyckoff sites

        Returns:
            if sucess, wyckoff_sites_tmp: list of wyckoff sites for valid sites
            otherwise, None

        """
        numMol_added = 0
        mol_sites_tmp = []

        # Now we start to add the specie to the wyckoff position
        sites_list = deepcopy(self.sites[id]) # the list of Wyckoff site
        if sites_list is not None:
            self.wyckoff_attempts = max(len(sites_list)*2, 10)
        else:
            # the minimum numattempts is to put all atoms to the general WPs
            min_wyckoffs = int(numMol/len(self.group.wyckoffs_organized[0][0]))
            self.wyckoff_attempts = max(2*min_wyckoffs, 10)

        for cycle in range(self.wyckoff_attempts):

            # Choose a random WP for given multiplicity: 2a, 2b, 2c
            if sites_list is not None:
                site = sites_list[0]
            else: # Selecting the merging 
                site = None

            # NOTE: The molecular version return wyckoff indices, not ops
            diff = numMol - numMol_added
            wp = wyc_mol(self.group, diff, site, valid_ori, True, self.dim)

            if wp is not False:
                # Generate a list of coords from the wyckoff position
                mult = wp.multiplicity # remember the original multiplicity
                pt = self.lattice.generate_point()

                # merge coordinates if the atoms are close
                mtol = pyxtal_mol.radius * 0.5
                pt, wp, oris = wp.merge(pt, self.lattice.matrix, mtol, valid_ori)

                if wp is not False:
                    if site is not None and mult != wp.multiplicity:
                        continue
                    if self.dim == 2 and self.thickness is not None and self.thickness < 0.1:
                        pt[-1] = 0.5 

                    ms0 = self._set_orientation(pyxtal_mol, pt, oris, wp)
                    if ms0 is not None:
                        # Check current WP against existing WP's  
                        passed_wp_check = True
                        for ms1 in mol_sites_tmp + mol_wyks:
                            if not ms0.check_with_ms2(ms1, tm=self.tol_matrix):
                                passed_wp_check = False
                        
                        if passed_wp_check:
                            if sites_list is not None:
                                sites_list.pop(0)

                            mol_sites_tmp.append(ms0)
                            numMol_added += len(ms0.wp)

                            # We have enough molecules of the current type
                            if numMol_added == numMol:
                                return mol_sites_tmp
        return None


    def _set_orientation(self, pyxtal_mol, pt, oris, wp): 
        """
        Generate good orientations
        """
        # Use a Wyckoff_site object for the current site
        self.numattempts += 1
        ori = random.choice(oris).copy()
        ori.change_orientation(flip=True)
        ms0 = mol_site(pyxtal_mol, pt, ori, wp, self.lattice, self.diag)
        # Check distances within the WP
        if ms0.check_distances():
            return ms0
        else:
            # Maximize the smallest distance for the general
            # positions if needed
            if len(pyxtal_mol.mol) > 1 and ori.degrees > 0:
                # bisection method
                def fun_dist(angle, ori, mo, pt):
                    # ori0 = ori.copy()
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
                for it in range(self.ori_attempts):
                    self.numattempts += 1
                    if (fun > 0.8) & (ms0.check_distances()):
                        return ms0
                    angle = (angle_lo + angle_hi) / 2
                    fun = fun_dist(angle, ori, pyxtal_mol, pt)
                    #print('Bisection: ', it, fun)
                    if fun_lo > fun_hi:
                        angle_hi, fun_hi = angle, fun
                    else:
                        angle_lo, fun_lo = angle, fun

        return None

    def _check_compatible(self):
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
        for i_mol, numIon in enumerate(self.numMols):
            # Get lists of multiplicity, maxn and freedom
            l_mult0 = []
            l_maxn0 = []
            l_free0 = []
            indices0 = []
            for i_wp, wp in enumerate(self.group):
                # Check that at least one valid orientation exists
                j, k = jk_from_i(i_wp, self.group.wyckoffs_organized)
                if len(self.valid_orientations[i_mol]) > j:
                    if len(self.valid_orientations[i_mol][j]) > k:
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
                if free:
                    if mult not in l_mult:
                        l_mult.append(mult)
                        l_maxn.append(maxn)
                        l_free.append(True)
                        indices.append(i_wp)
                #elif not free and i_wp not in used_indices:
                elif i_wp not in used_indices:
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
                #print("n == n0", n, n0)
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
                                indices.append(i_wp)
                    break
                # All combinations failed: return False
                if n == n0 and p >= len(l_mult) - 1:
                    #print("All combinations failed: return False")
                    return False, False
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
        if has_freedom:
            return True, True
        # All species passed, but no degrees of freedom: return 0
        else:
            return True, False

    def _check_consistency(self, site, numMol):
        """
        Check if the composition is consistent with symmetry
        """
        num = 0
        for s in site:
            num += int(s[:-1])
        if numMol == num:
            return True
        else:
            msg = "\nThe requested number of molecules is inconsistent: " + str(site)
            msg += "\nfrom numMols: {:d}".format(numMol)
            msg += "\nfrom Wyckoff list: {:d}".format(num)
            raise ValueError(msg)

