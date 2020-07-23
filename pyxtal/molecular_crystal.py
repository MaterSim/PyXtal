# Standard Libraries
import os
import random
import numpy as np
from copy import deepcopy

# External Libraries
from pymatgen.core.structure import Structure  # , Molecule

# from pymatgen.symmetry.analyzer import PointGroupAnalyzer

# PyXtal imports
# NOTE comprhys: surely importing constants to use as defaults
# is less clear than just specifying the defaults?
from pyxtal.constants import max1, max2, max3, max4
from pyxtal.msg import printx
from pyxtal.tolerance import Tol_matrix
from pyxtal.lattice import Lattice, cellsize, add_vacuum
from pyxtal.database.element import Element

# from pyxtal.database.collection import Collection
from pyxtal.wyckoff_site import mol_site, check_mol_sites
from pyxtal.molecule import pyxtal_molecule, orientation_in_wyckoff_position
from pyxtal.symmetry import (
    Group,  # choose_wyckoff,
    check_wyckoff_position,
    jk_from_i,
)
from pyxtal.operations import (
    apply_ops,
    check_images,
    project_point,
    distance,
    distance_matrix,
    filtered_coords,
)


# Define functions
# ------------------------------
def merge_coordinate_molecular(coor, lattice, group, tol, orientations):
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
        orientations: the valid orientations for a given molecule. Obtained
            from get_sg_orientations, which is called within molecular_crystal

    Returns:
        coor, index, point: (coor) is the new list of fractional coordinates after
        merging, and index is a single index of the Wyckoff position within
        the spacegroup. If merging is unsuccesful, or no index is found,
        returns the original coordinates and False. point is a 3-vector which can
        be plugged into the Wyckoff position to generate the rest of the points
    """
    # Get index of current Wyckoff position. If not one, return False
    index, point = check_wyckoff_position(coor, group)
    if index is False:
        return coor, False, None
    if point is None:
        return coor, False, None
    PBC = group.PBC
    # Main loop for merging multiple times
    while True:
        # Check distances of current WP. If too small, merge
        dm = distance_matrix([coor[0]], coor, lattice, PBC=PBC)
        passed_distance_check = True
        for i, x in enumerate(dm):
            for j, y in enumerate(x):
                if i != j and y < tol:
                    passed_distance_check = False
                    break
            if passed_distance_check is False:
                break
        # if check_images([coor[0]], ['C'], lattice, PBC=PBC, tol=tol) is False:
        if check_images([coor[0]], [6], lattice, PBC=PBC, tol=tol) is False:
            passed_distance_check = False
        if passed_distance_check is False:
            mult1 = group[index].multiplicity
            # Find possible wp's to merge into
            possible = []
            for i, wp in enumerate(group):
                mult2 = wp.multiplicity
                # Check that a valid orientation exists
                j, k = jk_from_i(i, orientations)
                if orientations[j][k] == []:
                    continue
                # Only allow smaller WP's that are an integer factor of the current one
                if (mult2 < mult1) and (mult1 % mult2 == 0):
                    possible.append(i)
            if possible == []:
                return coor, False, None
            # Calculate minimum separation for each WP
            distances = []
            for i in possible:
                wp = group[i]
                projected_point = project_point(point, wp[0], lattice=lattice, PBC=PBC)
                # NOTE new coor never used?
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


def choose_wyckoff_molecular(group, number, orientations, general_site_only=True):
    """
    Choose a Wyckoff position to fill based on the current number of molecules
    needed to be placed within a unit cell

    Rules:

        1) The new position's multiplicity is equal/less than (number).
        2) We prefer positions with large multiplicity.
        3) The site must admit valid orientations for the desired molecule.

    Args:
        group: a pyxtal.symmetry.Group object
        number: the number of molecules still needed in the unit cell
        orientations: the valid orientations for a given molecule. Obtained
            from get_sg_orientations, which is called within molecular_crystal

    Returns:
        a single index for the Wyckoff position. If no position is found,
        returns False
    """
    wyckoffs = group.wyckoffs_organized

    if general_site_only or np.random.random() > 0.5:  # choose from high to low
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                good_wyckoff = []
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append(w)
                if len(good_wyckoff) > 0:
                    return random.choice(good_wyckoff)
        return False
    else:
        good_wyckoff = []
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append(w)
        if len(good_wyckoff) > 0:
            return random.choice(good_wyckoff)
        else:
            return False


class molecular_crystal:
    """
    Class for storing and generating molecular crystals based on symmetry
    constraints. Based on the crystal.random_crystal class for atomic crystals.
    Given a spacegroup, list of molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. This crystal is stored as a pymatgen struct via self.struct

    Args:
        group: The international spacegroup number
            OR, a pyxtal.symmetry.Group object
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or give a string to convert (in which case fmt must be defined)
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
        check_atomic_distances: If True, checks the inter-atomic distances
            after each Wyckoff position is added. This requires slightly more
            time, but vastly improves accuracy. For approximately spherical
            molecules, or for large inter-molecular distances, this may be
            turned off
        fmt: Optional value for the input molecule string format. Used only
            when molecule values are strings
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        molecules,
        numMols,
        volume_factor,
        select_high=True,
        allow_inversion=False,
        orientations=None,
        check_atomic_distances=True,
        fmt="xyz",
        lattice=None,
        tm=Tol_matrix(prototype="molecular"),
    ):

        self.dim = 3
        """The number of periodic dimensions of the crystal"""
        # Necessary input
        self.PBC = [1, 1, 1]
        """The periodic axes of the crystal"""
        if type(group) != Group:
            group = Group(group, self.dim)
        self.sg = group.number
        self.selec_high = select_high
        self.init_common(
            molecules,
            numMols,
            volume_factor,
            select_high,
            allow_inversion,
            orientations,
            check_atomic_distances,
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
        check_atomic_distances,
        group,
        lattice,
        tm,
    ):
        """
        init functionality which is shared by 3D, 2D, and 1D crystals
        """
        self.numattempts = 0
        """The number of attempts needed to generate the crystal."""
        if type(group) == Group:
            self.group = group
            """A pyxtal.symmetry.Group object storing information about the space/layer
            /Rod/point group, and its Wyckoff positions."""
        else:
            self.group = Group(group, dim=self.dim)
        self.number = self.group.number
        """The international group number of the crystal:
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
        self.check_atomic_distances = check_atomic_distances
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
                self.struct = None
                return

        self.molecules = []  # A pyxtal_molecule objects,
        for mol in molecules:
            self.molecules.append(pyxtal_molecule(mol, self.tol_matrix))

        """
        The valid orientations for each molecule and Wyckoff position.
        May be copied when generating a new molecular_crystal to save a
        small amount of time
        """

        if orientations is None:
            self.get_orientations()
        else:
            self.valid_orientations = orientations

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
        index 1: the molecular prototype's index within self.molecules
        index 2: the Wyckoff position's 1st index (based on multiplicity)
        index 3: the WP's 2nd index (within the group of equal multiplicity)
        index 4: the index of the valid orientation for the molecule/WP pair

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

    def to_file(self, fmt=None, filename=None):
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
            if fmt is None:
                fmt = "cif"
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
            self.struct.to(fmt=fmt, filename=outdir)
            return "Output file to " + outdir
        elif self.valid:
            printx("Cannot create file: structure did not generate.", priority=1)

    def __str__(self):
        s = "------Random Molecular Crystal------"
        s += "\nDimension: " + str(self.dim)
        s += "\nGroup: " + self.group.symbol
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

    def print_all(self):
        print(str(self))

    def generate_crystal(self, max1=max1, max2=max2, max3=max3, max4=max4):
        """
        The main code to generate a random molecular crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
            max4: the number of attempts for changing the molecular orientation
        """
        atomic = self.check_atomic_distances
        # Check the minimum number of degrees of freedom within the Wyckoff positions
        degrees = self.check_compatible(self.group, self.numMols, self.valid_orientations)
        if degrees is False:
            self.struct = None
            self.valid = False
            raise ValueError(
                "the space group is incompatible with the number of molecules"
            )
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
            minvector = (
                max(all_lengths) * 1.2
            )  # Require slightly larger min vector for lattices

            for cycle1 in range(max1):
                self.cycle1 = cycle1
                # 1, Generate a lattice
                if self.lattice.allow_volume_reset is True:
                    self.volume = self.estimate_volume()
                    self.lattice.volume = self.volume
                self.lattice.reset_matrix()
                cell_matrix = self.lattice.get_matrix()

                for cycle2 in range(max2):
                    molecular_coordinates_total = []  # added molecular coordinates
                    molecular_sites_total = []  # corresponding molecular specie
                    coordinates_total = []  # added atomic coordinates
                    species_total = []  # corresponding atomic specie
                    wps_total = []  # corresponding Wyckoff position indices
                    points_total = []  # generating x,y,z points
                    mol_sites_total = []
                    good_structure = False

                    self.cycle2 = cycle2
                    molecular_coordinates_tmp = deepcopy(molecular_coordinates_total)
                    molecular_sites_tmp = deepcopy(molecular_sites_total)
                    coordinates_tmp = deepcopy(coordinates_total)
                    species_tmp = deepcopy(species_total)
                    wps_tmp = deepcopy(wps_total)
                    points_tmp = deepcopy(points_total)
                    mol_sites_tmp = []

                    # Add molecules specie by specie
                    for i, numMol in enumerate(self.numMols):
                        pyxtal_mol = self.molecules[i]
                        mol = pyxtal_mol.mol
                        symbols = pyxtal_mol.symbols
                        radius = pyxtal_mol.radius
                        tols_mat = pyxtal_mol.tols_matrix

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
                                point = self.lattice.generate_point()
                                projected_point = project_point(
                                    point, wp[0], lattice=cell_matrix, PBC=self.PBC
                                )

                                coords = apply_ops(projected_point, wp)
                                # merge coordinates if the atoms are close
                                if self.check_atomic_distances:
                                    mtol = radius * 0.5
                                else:
                                    mtol = radius * 2

                                (
                                    coords_toadd,
                                    good_merge,
                                    point,
                                ) = merge_coordinate_molecular(
                                    coords, cell_matrix, self.group, mtol, valid_ori
                                )
                                if good_merge is not False:
                                    wp_index = good_merge
                                    # scale the coordinates to [0,1], very important!
                                    coords_toadd = filtered_coords(
                                        coords_toadd, PBC=self.PBC
                                    )

                                    j, k = jk_from_i(
                                        wp_index, self.group.wyckoffs_organized
                                    )
                                    ori = random.choice(
                                        valid_ori[j][k]
                                    ).random_orientation()
                                    ms0 = mol_site(
                                        mol,
                                        point,
                                        symbols,
                                        ori,
                                        self.group[wp_index],
                                        cell_matrix,
                                        tols_mat,
                                        radius,
                                    )
                                    # Check distances within the WP

                                    if not ms0.check_distances(atomic):
                                        # Maximize the smallest distance for the general
                                        # positions if needed
                                        passed_ori = False
                                        if len(symbols) > 1 and ori.degrees > 0:
                                            # bisection method
                                            def fun_dist(angle, ori, mo, point):
                                                ori.change_orientation(angle)
                                                ms0 = mol_site(
                                                    mo,
                                                    point,
                                                    symbols,
                                                    ori,
                                                    self.group[wp_index],
                                                    cell_matrix,
                                                    tols_mat,
                                                    radius,
                                                )
                                                d = ms0.compute_distances()
                                                return d

                                            angle_lo = ori.angle
                                            angle_hi = angle_lo + np.pi
                                            fun_lo = fun_dist(angle_lo, ori, mol, point)
                                            fun_hi = fun_dist(angle_hi, ori, mol, point)
                                            fun = fun_hi
                                            for it in range(max4):
                                                if (fun > 0.8) & (
                                                    ms0.check_distances(atomic)
                                                ):
                                                    passed_ori = True
                                                    break
                                                angle = (angle_lo + angle_hi) / 2
                                                fun = fun_dist(angle, ori, mol, point)
                                                # print('Bisection: ', it, fun)
                                                if fun_lo > fun_hi:
                                                    angle_hi, fun_hi = angle, fun
                                                else:
                                                    angle_lo, fun_lo = angle, fun
                                    else:
                                        passed_ori = True
                                    if passed_ori is False:
                                        continue

                                    # Check distances with other WP's
                                    (
                                        coords_toadd,
                                        species_toadd,
                                    ) = ms0.get_coords_and_species()
                                    passed = True
                                    for ms1 in mol_sites_tmp:
                                        if not check_mol_sites(
                                            ms0, ms1, atomic, tm=self.tol_matrix
                                        ):
                                            passed = False
                                            break
                                    if not passed:
                                        continue
                                    else:
                                        # Distance checks passed; store the new Wyckoff position
                                        mol_sites_tmp.append(ms0)
                                        if len(coordinates_tmp) == 0:
                                            coordinates_tmp = coords_toadd
                                        else:
                                            coordinates_tmp = np.vstack(
                                                [coordinates_tmp, coords_toadd]
                                            )
                                        species_tmp += species_toadd
                                        numMol_added += len(coords_toadd) / len(symbols)
                                        if numMol_added == numMol:
                                            # We have enough molecules of the current type
                                            mol_sites_total = deepcopy(mol_sites_tmp)
                                            coordinates_total = deepcopy(coordinates_tmp)
                                            species_total = deepcopy(species_tmp)
                                            break

                        if numMol_added != numMol:
                            break  # need to repeat from the 1st species

                    if numMol_added == numMol:
                        printx(self.Msg6, priority=3)
                        good_structure = True
                        break
                    else:  # reset the coordinates and sites
                        molecular_coordinates_total = []
                        molecular_sites_total = []
                        wps_total = []
                # placing molecules here
                if good_structure:
                    final_lattice = cell_matrix
                    final_coor = []
                    final_site = []
                    final_number = []
                    self.mol_sites = []  # to regenerate the crystal

                    final_coor = deepcopy(coordinates_total)
                    final_site = deepcopy(species_total)
                    final_number = list(Element(ele).z for ele in species_total)
                    self.mol_sites = deepcopy(mol_sites_total)

                    final_coor = filtered_coords(final_coor, PBC=self.PBC)
                    lattice, coor = add_vacuum(final_lattice, final_coor, PBC=self.PBC)

                    self.lattice_matrix = lattice
                    self.frac_coords = np.array(coor)
                    self.cart_coords = np.dot(coor, lattice)
                    self.sites = final_site

                    # A pymatgen.core.structure.Structure object for
                    self.struct = Structure(lattice, self.sites, self.frac_coords)
                    # ASE/spglib structure
                    self.spg_struct = (lattice, self.frac_coords, final_number)
                    self.valid = True

                    return
                else:
                    printx("Failed final distance check.", priority=3)

        printx("Couldn't generate crystal after max attempts.", priority=1)
        if degrees == 0:
            printx("Note: Wyckoff positions have no degrees of freedom.", priority=2)
        self.struct = None
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
            or give a string to convert (in which case fmt must be defined)
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
        check_atomic_distances: If True, checks the inter-atomic distances
            after each Wyckoff position is added. This requires slightly more
            time, but vastly improves accuracy. For approximately spherical
            molecules, or for large inter-molecular distances, this may be
            turned off
        fmt: Optional value for the input molecule string format. Used only
            when molecule values are strings
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        molecules,
        numMols,
        volume_factor,
        select_high=True,
        allow_inversion=True,
        orientations=None,
        check_atomic_distances=True,
        fmt="xyz",
        thickness=None,
        lattice=None,
        tm=Tol_matrix(prototype="molecular"),
    ):

        self.dim = 2
        self.numattempts = 0
        if type(group) != Group:
            group = Group(group, self.dim)
        number = group.number  # The layer group number of the crystal."""
        self.thickness = thickness  # the thickness in Angstroms
        self.PBC = [1, 1, 0]
        self.init_common(
            molecules,
            numMols,
            volume_factor,
            select_high,
            allow_inversion,
            orientations,
            check_atomic_distances,
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
            or give a string to convert (in which case fmt must be defined)
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
        check_atomic_distances: If True, checks the inter-atomic distances
            after each Wyckoff position is added. This requires slightly more
            time, but vastly improves accuracy. For approximately spherical
            molecules, or for large inter-molecular distances, this may be
            turned off
        fmt: Optional value for the input molecule string format. Used only
            when molecule values are strings
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def __init__(
        self,
        group,
        molecules,
        numMols,
        volume_factor,
        select_high=True,
        allow_inversion=False,
        orientations=None,
        check_atomic_distances=True,
        fmt="xyz",
        area=None,
        lattice=None,
        tm=Tol_matrix(prototype="molecular"),
    ):
        self.dim = 1
        self.area = area  # the effective cross-sectional area in A^2
        self.PBC = [0, 0, 1]  # The periodic axes of the crystal (1,2,3)->(x,y,z)
        self.sg = None  # The international space group number, not rod groups
        self.init_common(
            molecules,
            numMols,
            volume_factor,
            select_high,
            allow_inversion,
            orientations,
            check_atomic_distances,
            group,
            lattice,
            tm,
        )
