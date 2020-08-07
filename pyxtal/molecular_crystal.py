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
from pyxtal.lattice import Lattice, cellsize, add_vacuum
from pyxtal.io import write_cif, structure_from_ext
from pyxtal.database.element import Element
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
        volume_factor=1.1,
        select_high=True,
        allow_inversion=False,
        orientations=None,
        check_atomic_distances=True,
        fmt="xyz",
        lattice=None,
        tm=Tol_matrix(prototype="molecular"),
        seed = None,
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
        self.seed = seed

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
        volume_factor=1.1,
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
        self.seed = None
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
            tm
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
        volume_factor=1.1,
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
        self.seed = None
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
            tm
        )
