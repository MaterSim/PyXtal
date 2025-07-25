"""
Module for generating and analyzing transformation operations. Several functions
for working with matrices are provided. The class OperationAnalyzer allows for
comparison between pymatgen.core.operations.SymmOp objects, and can be used to
identify conjugate operations. The orientation class can be used to identify
degrees of freedom for molecules in Wyckoff positions with certain symmetry
constraints.
"""

from copy import deepcopy

import numpy as np
from pymatgen.core import Element
from pymatgen.core.operations import SymmOp
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from pyxtal.constants import all_sym_directions, deg, hex_cell, rad
from pyxtal.tolerance import Tol_matrix


# ------------------------------
# Define functions
def check_distance(
    coord1,
    coord2,
    species1,
    species2,
    lattice,
    PBC=None,
    tm: Tol_matrix = None,
    d_factor=1.0,
):
    """
    Check the distances between two atom set. Only distances between points from
    different sets are checked.

    Args:
        coord1: a list of fractional coordinates e.g. [[.1,.6,.4]
            [.3,.8,.2]]
        coord2: a list of new fractional coordinates e.g. [[.7,.8,.9],
            [.4,.5,.6]]
        species1: a list of atomic species or numbers for coord1
        species2: a list of atomic species or numbers for coord2
        lattice: matrix describing the unit cell vectors
        PBC: A periodic boundary condition list,
            where 1 means periodic, 0 means not periodic.
            [1,1,1] -> full 3d periodicity,
            [0,0,1] -> periodicity along the z axis
        tm: a Tol_matrix object, or a string representing Tol_matrix
        d_factor: the tolerance is multiplied by this amount. Larger values
            mean atoms must be farther apart

    Returns:
        a bool for whether or not the atoms are sufficiently far enough apart
    """
    if tm is None:
        tm = Tol_matrix(prototype="atomic")

    # Check that there are points to compare
    if PBC is None:
        PBC = [1, 1, 1]
    if len(coord1) < 1 or len(coord2) < 1:
        return True

    # Create tolerance matrix from subset of tm
    tols = np.zeros((len(species1), len(species2)))
    for i1, specie1 in enumerate(species1):
        for i2, specie2 in enumerate(species2):
            tols[i1][i2] = tm.get_tol(specie1, specie2)

    # Calculate the distance between each i, j pair
    d = distance_matrix(coord1, coord2, lattice, PBC=PBC)

    return not (np.array(d) < np.array(tols)).any()


def verify_distances(coordinates, species, lattice, factor=1.0, PBC=None):
    """
    Checks the inter-atomic distance between all pairs of atoms in a crystal.

    Args:
        coordinates: a 1x3 list of fractional coordinates
        species: a list of atomic symbols for each coordinate
        lattice: a 3x3 matrix representing the lattice vectors of the unit cell
        factor: a tolerance factor for checking distances. A larger value means
            atoms must be farther apart
        PBC: A periodic boundary condition list, where 1 means periodic, 0 means
            not periodic. Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> 1d
            periodicity along the z axis

    Returns:
        True if no atoms are too close together, False if any pair is too close
    """
    if PBC is None:
        PBC = [1, 1, 1]
    for i, c1 in enumerate(coordinates):
        specie1 = species[i]
        for j, c2 in enumerate(coordinates):
            if j > i:
                specie2 = species[j]
                diff = np.array(c2) - np.array(c1)
                d_min = distance(diff, lattice, PBC=PBC)
                radius = Element(specie1).covalent_radius + Element(specie2).covalent_radius
                tol = factor * 0.5 * radius
                if d_min < tol:
                    return False
    return True


def check_images(
    coords,
    species,
    lattice,
    PBC=None,
    tm: Tol_matrix = None,
    tol=None,
    d_factor=1.0,
):
    """
    Given a set of (unfiltered) frac coordinates, checks if the periodic images
    are too close.

    Args:
        coords: a list of fractional coordinates
        species: the atomic species of each coordinate
        lattice: a 3x3 lattice matrix
        PBC: the periodic boundary conditions
        tm: a Tol_matrix object
        tol: a single override value for the distance tolerances
        d_factor: the tolerance is multiplied by this amount. Larger values
            mean atoms must be farther apart

    Returns:
        False if distances are too close. True if distances are not too close
    """
    if tm is None:
        tm = Tol_matrix(prototype="atomic")

    # If no PBC, there are no images to check
    if PBC is None:
        PBC = [1, 1, 1]

    if PBC == [0, 0, 0]:
        return True

    # Create image coords from given coords and PBC
    coords = np.array(coords)
    m = create_matrix(PBC=PBC, omit=True)

    new_coords = []
    for v in m:
        for v2 in coords + v:
            new_coords.append(v2)

    new_coords = np.array(new_coords)
    # Create a distance matrix
    dm = distance_matrix(coords, new_coords, lattice, PBC=[0, 0, 0])
    # Define tolerances
    if tol is None:
        tols = np.zeros((len(species), len(species)))
        for i, s1 in enumerate(species):
            for j, s2 in enumerate(species):
                if i <= j:
                    tols[i][j] = tm.get_tol(s1, s2)
                    tols[j][i] = tm.get_tol(s1, s2)
        tols2 = np.tile(tols, int(len(new_coords) / len(coords)))
        return not (dm < tols2).any()
    elif tol is not None:
        return not (dm < tol).any()
    return True


def distance(xyz, lattice, PBC=None):
    """
    Returns the Euclidean distance from the origin for a fractional
    displacement vector. Takes into account the lattice metric and periodic
    boundary conditions, including up to one non-periodic axis.

    Args:
        xyz: a fractional 3d displacement vector. Can be obtained by
            subtracting one fractional vector from another
        lattice: a 3x3 matrix describing a unit cell's lattice vectors
        PBC: A periodic boundary condition list, where 1 means periodic, 0 means
        not periodic. Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> 1d
        periodicity along the z axis

    Returns:
        a scalar for the distance of the point from the origin
    """
    if PBC is None:
        PBC = [1, 1, 1]
    xyz = filtered_coords(xyz, PBC=PBC)
    matrix = create_matrix(PBC=PBC)
    matrix += xyz
    matrix = np.dot(matrix, lattice)
    return np.min(np.linalg.norm(matrix, axis=1))


def distance_matrix(pts1, pts2, lattice, PBC=None, single=False, metric="euclidean"):
    """
    Returns the distances between two sets of fractional coordinates.
    Takes into account the lattice metric and periodic boundary conditions.

    Args:
        pts1: a list of fractional coordinates
        pts2: another list of fractional coordinates
        lattice: a 3x3 matrix describing a unit cell's lattice vectors
        PBC: A periodic boundary condition list, where 1 means periodic, 0 means
            not periodic. Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> 1d
            periodicity along the z axis
        single: return a scalor and matrix?
        metric: the metric to use with cdist, such as  'euclidean',
            'sqeuclidean', 'minkowski', and others

    Returns:
        a scalor or distance matrix
    """
    if PBC is None:
        PBC = [1, 1, 1]

    if PBC != [0, 0, 0]:
        l1 = filtered_coords(pts1, PBC=PBC)
        l2 = filtered_coords(pts2, PBC=PBC)
        l1 = np.dot(l1, lattice)
        l2 = np.dot(l2, lattice)
        matrix = create_matrix(PBC=PBC)
        matrix = np.dot(matrix, lattice)
        all_distances = np.zeros([len(matrix), len(l1), len(l2)])
        for i, v in enumerate(matrix):
            all_distances[i] += cdist(l1 + v, l2, metric)
        # m1 = np.array([(l1 + v) for v in matrix])
        # m1 = np.vstack([l1 + v for v in matrix])
        # all_distances = np.array([cdist(l, l2, metric) for l in m1])
        if single:
            return np.min(all_distances)
        else:
            # return np.apply_along_axis(np.min, 0, all_distances)
            return np.min(all_distances, axis=0)

    else:
        return distance_matrix_no_PBC(pts1, pts2, lattice, single, metric)


def distance_matrix_no_PBC(pts1, pts2, lattice, single=False, metric="euclidean"):
    """
    Returns the distances between two sets of fractional coordinates.
    Without periodic boundary conditions.

    Args:
        pts1: a list of fractional coordinates (N1*3)
        pts2: another list of fractional coordinates (N2*3)
        lattice: a 3x3 matrix describing a unit cell's lattice vectors
        single: return the minimum distance or the matrix
        metric: the metric to use with cdist. e.g. `euclidean`,
            `sqeuclidean`, `minkowski`, and others

    Returns:
        a scalor or distance matrix
    """

    l1 = np.dot(pts1, lattice)
    l2 = np.dot(pts2, lattice)
    d = cdist(l1, l2, metric)
    if single:
        return np.min(d)
    else:
        return d


def create_matrix(PBC=None, omit=False):
    """
    Used for calculating distances in lattices with periodic boundary
    conditions. When multiplied with a set of points, generates additional
    points in cells adjacent to and diagonal to the original cell

    Args:
        PBC: A periodic boundary condition list (1: periodic; 0: nonperiodic).

    Returns:
        A numpy array which can be multiplied by a set of coordinates
    """
    if PBC is None:
        PBC = [1, 1, 1]
    matrix = []
    i_list = [-1, 0, 1] if PBC[0] else [0]
    j_list = [-1, 0, 1] if PBC[1] else [0]
    k_list = [-1, 0, 1] if PBC[2] else [0]
    for i in i_list:
        for j in j_list:
            for k in k_list:
                if omit:
                    if [i, j, k] != [0, 0, 0]:
                        matrix.append([i, j, k])
                else:
                    matrix.append([i, j, k])
    return np.array(matrix, dtype=float)


def filtered_coords(coords, PBC=None):
    """
    Transform all coordinates to [0, 1] interval if PBC is allowed
    For example, [1.2, 1.6, -.4] becomes
    [0.2, 0.6, 0.6] when PBC=[1,1,1]
    [0.2, 1.6, 0.6] when PBC=[1,0,1]

    Args:
        coords: an array of real 3d vectors.
        PBC: A periodic boundary condition list (1: periodic; 0: nonperiodic).

    Returns:
        an array of filtered coords with the same shape as coords
    """

    if PBC is None:
        PBC = [1, 1, 1]
    if isinstance(coords, list):
        coords = np.array(coords)
    for i in range(3):
        if PBC[i] > 0:
            if len(coords.shape) > 1:
                coords[:, i] -= np.floor(coords[:, i])
            else:
                coords[i] -= np.floor(coords[i])
    return coords


def filtered_coords_euclidean(coords, PBC=None):
    """
    Given an array of fractional 3-vectors, filters coordinates to between 0 and
    1. Then, values which are greater than 0.5 are converted to 1 minus their
    value. This is used for converting displacement vectors with a Euclidean
    lattice.

    Args:
        coords: an array of real 3d vectors. The shape does not matter
        PBC: A periodic boundary condition list (1: periodic; 0: nonperiodic).

    Returns:
        an array of filtered coords with the same shape as coords
    """

    if PBC is None:
        PBC = [1, 1, 1]

    def filter_vector_euclidean(vector):
        for i, a in enumerate(PBC):
            if a:
                # QZ: check if this is equivalent to -= np.rint()
                vector[i] -= np.floor(vector[i])
                if vector[i] > 0.5:
                    vector[i] = 1 - vector[i]
        return vector

    return np.apply_along_axis(filter_vector_euclidean, -1, coords)


def get_inverse(op):
    """
    Given a SymmOp object, returns its inverse.

    Args:
        op: a Symmop object

    Returns:
        the inverse
    """
    matrix = op.affine_matrix.copy()
    # fill the matrix if it is ill conditioned
    # experimental
    if np.linalg.matrix_rank(matrix) < 4:
        for row in range(3):
            # fixed value
            if np.sum(matrix[row, :3] ** 2) < 1e-3:
                matrix[row, row] = 1
                matrix[row, 3] = 0

        if np.linalg.matrix_rank(matrix) == 3:
            # [-3x/2, -x/2, 1/4]
            # [0, x, 1/4]
            for rows in [[0, 1, 2], [1, 2, 0], [0, 2, 1]]:
                # m = (matrix[rows,:])[:,rows]
                # print(rows, m)
                if np.linalg.matrix_rank(matrix[rows[:2], :3]) != 2:
                    break
            id0, id1, id2 = rows[0], rows[1], rows[2]
            if matrix[id0, id1] == 0:
                matrix[id0, id1], matrix[id0, id0] = matrix[id0, id0], matrix[id0, id1]
                if np.linalg.matrix_rank(matrix) == 3:
                    matrix[id0, id1], matrix[id0, id2] = (
                        matrix[id0, id2],
                        matrix[id0, id1],
                    )
            else:
                matrix[id1, id0], matrix[id1, id1] = matrix[id1, id1], matrix[id1, id0]
                if np.linalg.matrix_rank(matrix) == 3:
                    matrix[id1, id0], matrix[id1, id2] = (
                        matrix[id1, id2],
                        matrix[id1, id0],
                    )

        elif np.linalg.matrix_rank(matrix) == 2:
            # -3x/2, -x/2, -x+1/4
            if np.sum(matrix[:, 0] ** 2) > 1e-3:
                matrix[1, 0], matrix[1, 1] = matrix[1, 1], matrix[1, 0]
                matrix[2, 0], matrix[2, 2] = matrix[2, 2], matrix[2, 0]
            elif np.sum(matrix[:, 1] ** 2) > 1e-3:
                matrix[0, 1], matrix[0, 0] = matrix[0, 0], matrix[0, 1]
                matrix[2, 1], matrix[2, 2] = matrix[2, 2], matrix[2, 1]
            else:
                matrix[0, 2], matrix[0, 0] = matrix[0, 0], matrix[0, 2]
                matrix[1, 2], matrix[1, 1] = matrix[1, 1], matrix[1, 2]
    return SymmOp(np.linalg.inv(matrix))


def get_inverse_ops(ops):
    """
    Given a inverse list of SymmOp objects

    Args:
        ops: a list of Symmop's

    Returns:
        a list of equal shape to ops, with the inverse operations
    """
    inverses = []
    for op in ops:
        if type(op) == SymmOp:
            inverses.append(op.inverse)
        else:
            inverses.append(get_inverse_ops(op))
    return inverses


def apply_ops(coord, ops):
    """
    Apply a list of SymmOps to a single 3-vector and return an array of
    the generated vectors. This is the inverse of SymmOp.operate_multi.

    Args:
        coord: a 3-vector (list or numpy array)
        ops: a list, tuple, or array of SymmOp objects

    Returns:
        an np array of floating-point 3-vectors
    """
    coord = np.array(coord)
    affine_point = np.concatenate([coord, np.ones(coord.shape[:-1] + (1,))], axis=-1)
    matrices = np.array([op.affine_matrix for op in ops])
    return np.inner(affine_point, matrices)[..., :-1]


def apply_ops_diagonal(coords, ops):
    """
    Given a list of coordinates and SymmOps, apply the ith op to the ith coord
    and return the list of transformed coordinates

    Args:
        coords: a list or array of 3-vectors

    Returns:
        a transformed numpy array of 3-vectors
    """
    coords = np.array(coords)
    affine_points = np.concatenate([coords, np.ones(coords.shape[:-1] + (1,))], axis=-1)
    matrices = np.array([op.affine_matrix for op in ops])
    return np.einsum("...ij,...j", matrices, affine_points)[:, :3]


def angle(v1, v2, radians=True):
    """
    Calculate the angle (in radians) between two vectors.

    Args:
        v1: a 1x3 vector
        v2: a 1x3 vector
        radians: whether to return angle in radians (default) or degrees

    Returns:
        the angle in radians between the two vectors
    """
    v1 = np.real(v1)
    v2 = np.real(v2)
    dot = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    # if np.isclose(dot, 1.0):
    if np.abs(dot - 1) < 1e-3:
        a = 0
    elif np.abs(dot + 1) < 1e-3:
        a = np.pi
    else:
        a = np.arccos(dot)

    if radians:
        return a
    else:
        return a * deg


def is_orthogonal(m, tol=0.001):
    """
    Check whether or not a 3x3 matrix is orthogonal. An orthogonal matrix has
    the property that it is equal to its transpose.

    Args:
        m: a 3x3 matrix (list or numpy array)
        tol: the numerical tolerance for checking if two matrices are equal

    Returns:
        True if the matrix is orthogonal, False if it is not
    """
    m1 = np.dot(m, np.transpose(m))
    m2 = np.dot(np.transpose(m), m)
    return not (not np.allclose(m1, np.identity(3), rtol=tol) or not np.allclose(m2, np.identity(3), rtol=tol))


def aa2matrix(axis, angle, radians=True, random=False):
    """
    Given an axis and an angle, return a 3x3 rotation matrix.
    Based on:
    https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle

    Args:
        axis: a vector about which to perform a rotation
        angle: the angle of rotation
        radians: whether the supplied angle is in radians (True)
            or in degrees (False)
        random: whether or not to choose a random rotation matrix. If True, the
            axis and angle are ignored, and a random orientation is generated

    Returns:
        a 3x3 numpy array representing a rotation matrix
    """
    # Convert to radians if necessary
    if not radians:
        angle *= rad
    # Allow for generation of random rotations
    if random:
        # a = np.random.random()
        axis = np.random.sample(3)
        angle = np.random.random() * np.pi * 2
    # Ensure axis is a unit vector
    axis = axis / np.linalg.norm(axis)
    # Define quantities which are reused
    x = np.real(axis[0])
    y = np.real(axis[1])
    z = np.real(axis[2])
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1 - c
    # Define the rotation matrix
    Q = np.zeros([3, 3])
    Q[0][0] = x * x * C + c
    Q[0][1] = x * y * C - z * s
    Q[0][2] = x * z * C + y * s
    Q[1][0] = y * x * C + z * s
    Q[1][1] = y * y * C + c
    Q[1][2] = y * z * C - x * s
    Q[2][0] = z * x * C - y * s
    Q[2][1] = z * y * C + x * s
    Q[2][2] = z * z * C + c
    return Q


def rotate_vector(v1, v2, rtol=1e-4):
    # TODO: Verify that multiplication order is correct
    # (matrix should come after vector in np.dot)
    """
    Rotates a vector v1 to v2 about an axis perpendicular to both. Returns the
    3x3 rotation matrix used to do so.

    Args:
        v1: a 1x3 vector (list or array) of floats
        v2: a 1x3 vector (list or array) of floats

    Returns:
        a 3x3 matrix corresponding to a rotation which
        can be applied to v1 to obtain v2
    """
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    dot = np.dot(v1, v2)
    # Handle collinear vectors
    if np.abs(dot - 1) < rtol:
        return np.identity(3)
    elif np.abs(dot + 1) < rtol:
        r = np.random.sample(3) #[np.random.random(), np.random.random(), np.random.random()]
        v3 = np.cross(v1, r)
        v3 /= np.linalg.norm(v3)
        # return aa2matrix(v3, np.pi)
        return Rotation.from_rotvec(np.pi * v3).as_matrix()
    theta = angle(v1, v2)
    v3 = np.cross(v1, v2)
    v3 /= np.linalg.norm(v3)
    # return aa2matrix(v3, theta)
    return Rotation.from_rotvec(theta * v3).as_matrix()


def are_equal(op1, op2, PBC=None, rtol=1e-3, atol=1e-3):
    """
    Check whether two SymmOp objects are equal up to some numerical tolerance.
    Allows for optional consideration of periodic boundary conditions. This
    option is useful for handling the periodicity of crystals.

    Args:
        op1: a SymmOp object
        op2: another SymmOp object
        PBC: A periodic boundary condition list.
        rtol: the relative numerical tolerance for equivalence
        atol: the absolute numerical tolerance for equivalence

    Returns:
        True if op1 and op2 are equivalent, False otherwise
    """
    # Check two SymmOps for equivalence
    # pbc=True means integer translations will be ignored
    if PBC is None:
        PBC = [1, 1, 1]
    m1 = op1.rotation_matrix
    m2 = op2.rotation_matrix
    # Check that rotations are equivalent
    if not np.allclose(m1, m2, rtol=rtol, atol=atol):
        return False
    v1 = op1.translation_vector
    v2 = op2.translation_vector

    difference = v2 - v1

    for i, a in enumerate(PBC):
        if a:
            difference[i] -= np.rint(difference[i])

    d = np.linalg.norm(difference)

    return abs(d) < rtol


class OperationAnalyzer(SymmOp):
    """
    Class for comparing operations. Stores rotation axis and angle, as well as
    the type and order of operation (identity, inversion, rotation, or
    rotoinversion). By default, takes a SymmOp as argument. This information can
    be accessed by calling print(object). The class methods is_conjugate and
    are_conjugate can be used to determine if two operations are conjugate to
    each other. That is, whether they represent the same angle of rotation and
    are either both inverted or both uninverted.

    Note: rotoinversions with odd-order rotational parts will have an over-all
        even order. For example, the order of (-3) is 6.

    Note: reflections are treated as rotoinversions of order 2.

    Args:
        SymmOp: a pymatgen.core.structure.SymmOp object to analyze
    """

    # ++++: include support for off-center operations
    # TODO: include support for shear and scaling operations
    # TODO: include support for matrix-column and axis-angle initialization
    def get_order(angle, rotoinversion=False, tol=1e-2):
        # Find the order of a rotation based on its angle
        found = False
        for n in range(1, 61):
            x = (n * angle) / (2.0 * np.pi)
            y = x - np.rint(x)
            if abs(y) <= tol:
                found = True
                break
        if found:
            # Double order of odd-rotation rotoinversions
            if rotoinversion:
                if n % 2 == 1:
                    return int(n * 2)
                else:
                    return int(n)
            else:
                return int(n)
        if not found:
            return "irrational"
        return None

    def __init__(self, op, parse_trans=False, hexagonal=False):
        if type(op) == deepcopy(SymmOp):
            # The numerical tolerance associated with op
            # The 4x4 affine matrix of the op
            # The 3x3 rotation (or rotoinversion) matrix
            # The determinant of self.m
            self.op = op
            self.tol = op.tol
            self.affine_matrix = op.affine_matrix
            self.m = op.rotation_matrix
            self.det = np.linalg.det(self.m)
            self.hexagonal = hexagonal

        elif (type(op) == np.ndarray) or (type(op) == np.matrix):
            if op.shape == (3, 3):
                self.op = SymmOp.from_rotation_and_translation(op, [0, 0, 0])
                self.m = self.op.rotation_matrix
                self.det = np.linalg.det(op)
        else:
            raise ValueError("Error: OperationAnalyzer requires a SymmOp or 3x3 array.")

        self.symbol = None
        self.parse_trans = parse_trans  # only for space group

        # If rotation matrix is not orthogonal
        if not is_orthogonal(self.m):
            if False: #self.hexagonal:
                #if self.op.as_xyz_string():
                #    print('debug')
                pass
            else:
                self.type = "general"
                self.axis, self.angle, self.order, self.rotation_order = (
                    None,
                    None,
                    None,
                    None,
                )
        # If rotation matrix is orthogonal
        else:
            # If determinant is positive, rotation
            if np.linalg.det(self.m) > 0:
                self.inverted = False
                rotvec = Rotation.from_matrix(self.m).as_rotvec()
                if np.sum(rotvec.dot(rotvec)) < 1e-6:
                    self.axis = None
                    self.angle = 0
                else:
                    self.angle = np.linalg.norm(rotvec)
                    self.axis = rotvec / self.angle
                    if self.hexagonal:
                        # print('convert hex', self.axis, np.dot(self.axis, hex_cell))
                        if np.linalg.norm(self.axis - np.array([1, 0, 0])) > 1e-4:
                            self.axis = np.dot(self.axis, hex_cell)
                    # parse symmetry direction
                    if self.parse_trans and not self.parse_axis():
                        self.axis *= -1
                        self.angle = 2 * np.pi - self.angle
                        # print('switch angle', self.angle)

                if np.isclose(self.angle, 0):
                    # Types: 'identity', 'inversion', 'rotation', or 'rotoinversion'.
                    # order: the number of times to get to origin
                    # rotation_order: 2, 3, 4, 6
                    self.type = "identity"
                    self.order = 1
                    self.rotation_order = 1
                    self.symbol = "1"
                else:
                    self.type = "rotation"
                    self.order = OperationAnalyzer.get_order(self.angle)
                    self.rotation_order = self.order
                    if self.parse_trans:
                        self.symbol = self.parse_screw_symmetry()
                        #print('debug', self.symbol)

            # If determinant is negative, rotoinversion (including reflection)
            elif np.linalg.det(self.m) < 0:
                self.inverted = True
                rotvec = Rotation.from_matrix(-1 * self.m).as_rotvec()
                if np.sum(rotvec.dot(rotvec)) < 1e-6:
                    self.axis = None
                    self.angle = 0
                else:
                    self.angle = np.linalg.norm(rotvec)
                    self.axis = rotvec / self.angle
                    if self.hexagonal:
                        if np.linalg.norm(self.axis - np.array([1, 0, 0])) > 1e-4:
                            #print('convert hex', self.axis, np.dot(self.axis, hex_cell))
                            self.axis = np.dot(self.axis, hex_cell)
                if np.isclose(self.angle, 0):
                    #self.axis = None
                    self.symbol = "-1"
                    self.type = "inversion"
                    self.order = 2
                    self.rotation_order = 1
                else:
                    # parse symmetry direction
                    # if self.parse_trans and not self.parse_axis():
                    #    self.axis *= -1
                    #    self.angle = 2*np.pi - self.angle
                    #    print('switch angle', self.angle)
                    self.axis *= -1
                    self.type = "rotoinversion"
                    self.order = OperationAnalyzer.get_order(self.angle, rotoinversion=True)
                    #print('debug', self.angle, self.axis, self.order)
                    self.rotation_order = OperationAnalyzer.get_order(self.angle, rotoinversion=False)
                    if self.parse_trans:
                        self.symbol = self.parse_glide_symmetry()
            elif np.linalg.det(self.m) == 0:
                self.type = "degenerate"
                self.axis, self.angle = None, None
        if self.symbol is None:
            if self.type == "rotoinversion":
                if self.rotation_order == 2:
                    self.symbol = 'm'
                elif self.rotation_order == 3:
                    self.symbol = '-3'
                elif self.rotation_order == 4:
                    self.symbol = '-4'
                elif self.rotation_order == 6:
                    self.symbol = '-6'
            elif self.type == 'rotation':
                self.symbol = str(self.order)

    def parse_screw_symmetry(self, tol=1e-2):
        """
        If the point group symmetry is rotation, parse the screw vector

        Returns:
            0, 2_1, 3_1, 3_2, 4_1, 4_2, 4_3, 6_1, 6_2, 6_3, 6_4, 6_5
        """
        # only count the translation on the given axis
        vec = self.translation_vector.copy()
        if np.isclose(abs(np.dot(self.axis, np.array([1, 0, 0]))), 1):
            vec[1] = 0
            vec[2] = 0
        elif np.isclose(abs(np.dot(self.axis, np.array([0, 1, 0]))), 1):
            vec[0] = 0
            vec[2] = 0
        elif np.isclose(abs(np.dot(self.axis, np.array([0, 0, 1]))), 1):
            vec[0] = 0
            vec[1] = 0
        else:
            vec = 0
        # No screw symmetry for other directios for tetragonal????

        if np.linalg.norm(vec) < tol:
            return str(self.order)  # 2, 3, 4, 6
        else:
            trans = np.sum(vec)
            if self.order == 2:
                return "2_1"
            elif self.order == 3:
                if abs(self.angle / trans - 2 * np.pi) < tol:
                    return "3_1"
                else:
                    return "3_2"
            elif self.order == 4:
                if abs(trans) < tol:  #
                    return "4_2"
                elif abs(self.angle / trans - 2 * np.pi) < tol:
                    return "4_1"
                elif abs(self.angle / trans - np.pi) < tol or abs(self.angle / trans - 3 * np.pi) < tol:
                    return "4_2"
                elif abs(self.angle / trans - 2 / 3 * np.pi) < tol or abs(self.angle / trans - 6 * np.pi) < tol:
                    return "4_3"
            elif self.order == 6:
                if abs(trans) < tol:  #
                    return "6_3"
                elif abs(self.angle / trans - 2 * np.pi) < tol:
                    return "6_1"
                elif abs(self.angle / trans - np.pi) < tol or abs(self.angle / trans - 5 / 2 * np.pi) < tol:
                    return "6_2"
                elif abs(self.angle / trans - 2 / 3 * np.pi) < tol or abs(self.angle / trans - 10 / 3 * np.pi) < tol:
                    return "6_3"
                elif abs(self.angle / trans - 1 / 2 * np.pi) < tol or abs(self.angle / trans - 5 * np.pi) < tol:
                    return "6_4"
                elif abs(self.angle / trans - 2 / 5 * np.pi) < tol or abs(self.angle / trans - 10 * np.pi) < tol:
                    return "6_5"

        #print("Cannot assign symbol", self.angle, trans)
        return "1" #None

    def parse_glide_symmetry(self, tol=1e-2):
        """
        If the point group symmetry is rotation, parse the screw vector
        Returns:
            m, a, b, c, n, d
        """
        if self.rotation_order > 2:
            return "-" + str(self.rotation_order)
        elif abs(self.angle - np.pi) > tol:
            #print("return m2", self.op)
            return "m"  # just indicate
        else:
            vec = self.translation_vector.copy()
            if np.isclose(abs(np.dot(self.axis, np.array([1, 0, 0]))), 1):
                vec[0] = 0
            elif np.isclose(abs(np.dot(self.axis, np.array([0, 1, 0]))), 1):
                vec[1] = 0
            elif np.isclose(abs(np.dot(self.axis, np.array([0, 0, 1]))), 1):
                vec[2] = 0

            if np.linalg.norm(vec) < tol:
                #print("return m1", self.op.as_xyz_str(), self.axis)
                return "m"
            else:
                if np.linalg.norm(vec - np.array([1 / 2, 0, 0])) < tol:
                    return "a"
                elif np.linalg.norm(vec - np.array([0, 1 / 2, 0])) < tol:
                    return "b"
                elif np.linalg.norm(vec - np.array([0, 0, 1 / 2])) < tol:
                    #print("return c1", self.op.as_xyz_str(), self.axis)
                    return "c"
                elif (
                    np.linalg.norm(vec - np.array([0, 1 / 2, 1 / 2])) < tol
                    or np.linalg.norm(vec - np.array([1 / 2, 0, 1 / 2])) < tol
                    or np.linalg.norm(vec - np.array([1 / 2, 1 / 2, 0])) < tol
                ):
                    return "n"
                elif np.linalg.norm(vec - np.array([1 / 2, 1 / 2, 1 / 2])) < tol:
                    if np.isclose(abs(np.dot(self.axis, np.array([0, -0.7071, 0.7071]))), 1) or np.isclose(
                        abs(np.dot(self.axis, np.array([-0.7071, 0, 0.7071]))), 1
                    ):
                        return "n"
                    else:
                        #print("return c2", self.op)
                        return "c"
                else:
                    return "d"

    def parse_axis(self):
        """
        parse if the axis follows the standard convention
        """
        ax = self.axis#; print('debug ax', ax)
        ax /= np.linalg.norm(ax)
        for direction in all_sym_directions:
            normed = direction / np.linalg.norm(direction)
            # print(direction, np.dot(direction, ax))
            if np.isclose(np.dot(normed, ax), 1):
                return True
            elif np.isclose(np.dot(normed, ax), -1):
                return False
        raise ValueError("Cannot find the symmetry direction", ax)

    def __str__(self):
        """
        A custom printing string for the object. The type, order, angle, and
        axis are printed. Converts values close to 0 to 0 for readability. Also
        only prints the real part of the axis.
        """
        # Avoid printing '-0.' instead of '0.'
        #print('type', self.type, self.op.as_xyz_str())
        if self.type not in ['inversion', 'identity'] and len(self.axis) == 3:
            for i, x in enumerate(self.axis):
                if np.isclose(x, 0):
                    self.axis[i] = 0.0
            self.axis = np.real(self.axis)
        #else:
        #    self.axis = None #np.array([1, 0, 0])
        return (
            f" Operation: {self.symbol:4s} {self.type}"
            + f"\n {self.op.as_xyz_str()} at: {self.axis}"
        )

    def is_conjugate(self, op2):
        """
        Returns whether or not another operation is conjugate (the same
        operation in a different reference frame). Rotations with the same order
        will not always return True. For example, a 5/12 and 1/12 rotation will
        not be considered conjugate.

        Args:
            op2: a SymmOp or OperationAnalyzer object to compare with

        Returns:
            True if op2 is conjugate to self.op, and False otherwise
        """
        if type(op2) != OperationAnalyzer:
            opa2 = OperationAnalyzer(op2)
            if opa2.type == self.type:
                if self.type == "rotation" or self.type == "rotoinversion":
                    ratio = self.angle / opa2.angle
                    if np.isclose(np.fabs(ratio), 1.0, atol=1e-2):
                        return True
                    return None
                elif self.type == "identity" or self.type == "inversion":
                    return True
                return None
            else:
                return False
        else:
            if op2.type == self.type:
                if self.type == "rotation" or self.type == "rotoinversion":
                    ratio = self.angle / op2.angle
                    if np.isclose(ratio, 1.0, atol=1e-2):
                        return True
                    return None
                elif self.type == "identity" or self.type == "inversion":
                    return True
                return None
            else:
                return False

    def are_conjugate(op1, op2):
        """
        Returns whether or not two operations are conjugate (the same
        operation in a different reference frame). Rotations with the same order
        will not always return True. For example, a 5/12 and 1/12 rotation will
        not be considered conjugate.

        Args:
            op1: a SymmOp or OperationAnalyzer object
            op2: a SymmOp or OperationAnalyzer object to compare with op1

        Returns:
            True if op2 is conjugate to op1, and False otherwise
        """
        if type(op1) != OperationAnalyzer:
            opa1 = OperationAnalyzer(op1)
        return opa1.is_conjugate(op2)


def find_ids(coords, ref, tol=1e-3):
    """
    find the refernce ids that can match
    """
    ids = []
    # print('ref', ref)
    for coord in coords:
        diffs = ref - coord
        diffs -= np.rint(diffs)
        norms = np.linalg.norm(diffs, axis=1)
        # print(norms, diffs)
        for i, norm in enumerate(norms):
            if norm < tol and i not in ids:
                ids.append(i)
                break
    return ids


def get_best_match(positions, ref, cell):
    """
    find the best match with the reference from a set of positions

    Args:
        positions: N*3 array
        ref: 1*3 array
        cell: cell matrix 3*3 array

    Returns:
        position: matched position
        id: matched id
    """
    diffs = positions - ref
    diffs -= np.rint(diffs)
    diffs = np.dot(diffs, cell)
    dists = np.linalg.norm(diffs, axis=1)
    id = np.argmin(dists)
    return positions[id], dists[id]


# Test Functionality
if __name__ == "__main__":
    # ----------------------------------------------------
    from pymatgen.symmetry.analyzer import generate_full_symmops

    rot = Rotation.from_rotvec(np.pi / 6 * np.array([1, 0, 0]))
    op = SymmOp.from_rotation_and_translation(rot.as_matrix(), np.array([1, 0, 0]))
    ops = [op]

    symm_m = generate_full_symmops(ops, 1e-3)
    for op in symm_m:
        opa = OperationAnalyzer(op)
        print(opa.order)
