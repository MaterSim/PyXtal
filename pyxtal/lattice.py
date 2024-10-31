# Standard Libraries
from __future__ import annotations

import numpy as np
from numpy.random import Generator

# PyXtal imports
from pyxtal.constants import deg, ltype_keywords, rad
from pyxtal.msg import VolumeError
from pyxtal.operations import angle, create_matrix


class Lattice:
    """
    Class for storing and generating crystal lattices. Allows for
    specification of constraint values. Lattice types include triclinic,
    monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic,
    spherical, and ellipsoidal. The last two are used for generating point
    group structures, and do not actually represent a parallelepiped lattice.

    Args:
        ltype: a string representing the type of lattice (from the above list)
        volume: the volume, in Angstroms cubed, of the lattice
        matrix: matrix in 3*3 form
        PBC: A periodic boundary condition list, where 1 is periodic,
            Ex: [1, 1, 1] -> 3d periodicity, [0, 0, 1] -> periodic at z axis
        kwargs: various values which may be defined. If none are defined,
            random ones will be generated. Values will be passed to
            generate_lattice. Options include:
            'area': The cross-sectional area (in Ang^2). Only for 1D crystals
            'thickness': The cell's thickness (in Angstroms) for 2D crystals
            'unique_axis': The unique axis for certain symmetry (and especially
                layer) groups. Because the symmetry operations are not also
                transformed, you should use the default values for random
                crystal generation
            'random': If False, keeps the stored values for the lattice geometry
                even upon applying reset_matrix. To alter the matrix, use
                set_matrix() or set_para
            'unique_axis': the axis ('a', 'b', or 'c') which is not symmetrically
                equivalent to the other two
            'min_l': the smallest allowed cell vector. The smallest vector must
                be larger than this.
            'mid_l': the second smallest allowed cell vector. The second
                smallest vector must be larger than this.
            'max_l': the third smallest allowed cell vector. The largest cell
                vector must be larger than this.
            'allow_volume_reset': a bool stating whether or not the volume
                should be reset during each crystal generation attempt
    """

    def __init__(
        self, ltype, volume=None, matrix=None, PBC=None, random_state: int | None | Generator = None, **kwargs
    ):
        # Set required parameters
        if PBC is None:
            PBC = [1, 1, 1]

        if ltype in ltype_keywords:
            self.ltype = ltype.lower()
        elif ltype is None:
            self.ltype = "triclinic"
        else:
            msg = "Invalid lattice type: " + ltype
            raise ValueError(msg)

        self.volume = float(volume)
        self.PBC = PBC
        self.dim = sum(PBC)
        self.kwargs = {}
        self.random = True

        if isinstance(random_state, Generator):
            self.random_state = random_state.spawn(1)[0]
        else:
            self.random_state = np.random.default_rng(random_state)

        # Set optional values
        self.allow_volume_reset = True
        for key, value in kwargs.items():
            if key in [
                "area",
                "thickness",
                "unique_axis",
                "random",
                "min_l",
                "mid_l",
                "max_l",
                "min_special",  # min special for mole
            ]:
                setattr(self, key, value)
                self.kwargs[key] = value
                if key == "allow_volume_reset" and value is False:
                    self.allow_volume_reset = False

        if not hasattr(self, "unique_axis"):
            self.unique_axis = "c"

        # Set stress normalization info
        if self.ltype == "triclinic":
            norm_matrix = np.ones([3, 3])

        elif self.ltype == "monoclinic":
            if self.PBC == [1, 1, 1]:
                norm_matrix = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 1]])
            else:
                if self.unique_axis == "a":
                    norm_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 1, 1]])
                elif self.unique_axis == "b":
                    norm_matrix = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 1]])
                elif self.unique_axis == "c":
                    norm_matrix = np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1]])

        elif self.ltype in [
            "orthorhombic",
            "tetragonal",
            "trigonal",
            "hexagonal",
            "cubic",
        ]:
            norm_matrix = np.eye(3)

        elif self.ltype in ["spherical", "ellipsoidal"]:
            norm_matrix = np.zeros([3, 3])

        self.stress_normalization_matrix = norm_matrix

        # Set info for on-diagonal stress symmetrization
        if self.ltype in ["tetragonal", "trigonal", "hexagonal"]:
            self.stress_indices = [(0, 0), (1, 1)]

        elif self.ltype in ["cubic"]:
            self.stress_indices = [(0, 0), (1, 1), (2, 2)]

        else:
            self.stress_indices = []

        # Set values for the matrix
        if matrix is None:
            self.reset_matrix()
        else:
            self.set_matrix(matrix)

        # Set tolerance
        if self.ltype in ["triclinic"]:
            self.a_tol = 15.0
        else:
            self.a_tol = 9.9
        self._get_dof()

    def _get_dof(self):
        """
        get the number of degree of freedom
        """
        if self.ltype in ["triclinic"]:
            self.dof = 6
        elif self.ltype in ["monoclinic"]:
            self.dof = 4
        elif self.ltype in ["orthorhombic"]:
            self.dof = 3
        elif self.ltype in ["tetragonal", "hexagonal", "trigonal"]:
            self.dof = 2
        else:
            self.dof = 1

    def get_bounds(self, min_vec=2.0, max_vec=50.0, min_ang=30, max_ang=150):
        """
        get the number of degree of freedom
        """
        if self.ltype in ["triclinic"]:
            self.bounds = [(min_vec, max_vec), (min_vec, max_vec), (min_vec, max_vec),
                           (min_ang, max_ang), (min_ang, max_ang), (min_ang, max_ang)]
        elif self.ltype in ["monoclinic"]:
            self.bounds = [(min_vec, max_vec), (min_vec, max_vec), (min_vec, max_vec),
                           (min_ang, max_ang)]
        elif self.ltype in ["orthorhombic"]:
            self.bounds = [(min_vec, max_vec), (min_vec, max_vec), (min_vec, max_vec)]
        elif self.ltype in ["tetragonal", "hexagonal", "trigonal"]:
            self.bounds = [(min_vec, max_vec), (min_vec, max_vec)]
        else:
            self.bounds = [(min_vec, max_vec)]
        return self.bounds

    @classmethod
    def get_dofs(self, ltype):
        """
        get the number of degree of freedom
        """
        # if ltype is None: ltype = self.ltype

        if ltype in ["triclinic"]:
            dofs = [3, 3]
        elif ltype in ["monoclinic"]:
            dofs = [3, 1]
        elif ltype in ["orthorhombic"]:
            dofs = [3, 0]
        elif ltype in ["tetragonal", "hexagonal", "trigonal"]:
            dofs = [2, 0]
        else:
            dofs = [1, 0]
        return dofs

    def copy(self):
        """
        simply copy the structure
        """
        from copy import deepcopy

        return deepcopy(self)

    def get_lengths(self):
        mat = create_matrix(self.PBC, True)
        mat = np.dot(mat, self.matrix)
        return mat, np.linalg.norm(mat, axis=1)

    def scale(self, factor=1.1):

        matrix = self.matrix
        return Lattice.from_matrix(matrix * factor, ltype=self.ltype)

    def get_permutation_matrices(self):
        """
        Return the possible permutation matrices that donot violate the symmetry
        """
        if self.ltype in ["monoclinic"]:  # permutation between a and c
            return np.array(
                [
                    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # self
                    [[0, 0, 1], [0, 1, 0], [1, 0, 0]],  # a-c
                ]
            )
        elif self.ltype in ["triclinic"]:
            return np.array(
                [
                    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # self
                    [[1, 0, 0], [0, 0, 1], [0, 1, 0]],  # b-c
                    [[0, 0, 1], [0, 1, 0], [1, 0, 0]],  # a-c
                    [[0, 1, 0], [1, 0, 0], [0, 0, 1]],  # a-b
                ]
            )
        else:
            return [np.eye(3)]

    def get_transformation_matrices(self):
        """
        Return possible transformation matrices that donot violate the symmetry
        """
        if self.ltype in ["monoclinic"]:
            return np.array(
                [
                    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, 0], [0, 1, 0], [1, 0, 1]],
                    [[1, 0, 0], [0, 1, 0], [-1, 0, 1]],
                    [[1, 0, 1], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, -1], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # change angle
                    # [[-1,0,0],[0,1,0],[0,0,1]], #change angle
                ]
            )

        elif self.ltype in ["triclinic"]:
            return np.array(
                [
                    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, 0], [0, 1, 0], [1, 0, 1]],
                    [[1, 0, 0], [0, 1, 0], [-1, 0, 1]],
                    [[1, 0, 1], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, -1], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, 0], [0, 1, 0], [0, 1, 1]],
                    [[1, 0, 0], [0, 1, 1], [0, 0, 1]],
                    [[1, 0, 0], [0, 1, 0], [0, -1, 1]],
                    [[1, 0, 0], [0, 1, -1], [0, 0, 1]],
                    [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
                    [[1, -1, 0], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, 0], [1, 1, 0], [0, 0, 1]],
                    [[1, 0, 0], [-1, 1, 0], [0, 0, 1]],
                    # [[-1,0,0],[0,-1,0],[0,0,1]],
                    # [[1,0,0],[0,-1,0],[0,0,-1]],
                    # [[-1,0,0],[0,1,0],[0,0,-1]],
                    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
                ]
            )
        else:
            return [np.eye(3)]

    def search_transformations(self, lat_ref, d_tol=1.0, f_tol=0.1):
        """
        search the closest match to the reference lattice object

        Args:
            lat_ref: reference lattice object
            d_tol: tolerance in angle
            f_tol:
            a_tol:

        Returns:
            a two steps of transformation matrix if the match is possible
        """
        # Find all possible permutation and transformation matrices
        trans1 = self.get_permutation_matrices()
        trans2 = self.get_transformation_matrices()
        tols = np.zeros([len(trans2) * len(trans1), 3])
        trans = []
        switchs = []

        count = 0
        for _i, tran1 in enumerate(trans1):
            lat0 = self.transform(tran1)
            for _j, tran2 in enumerate(trans2):
                tmp = np.dot(tran2, lat0.matrix)
                try:
                    # print("Start", np.linalg.det(tmp))
                    lat2 = Lattice.from_matrix(tmp, ltype=self.ltype)
                    # print("End", np.linalg.det(lat2.matrix))
                    d_tol1, f_tol1, a_tol1, switch = lat2.get_diff(lat_ref)
                    # print(d_tol1, f_tol1, a_tol1, switch)
                except:
                    d_tol1, f_tol1, a_tol1, switch = 10, 1.0, 90, None
                tols[count] = [d_tol1, f_tol1, a_tol1]
                trans.append([tran1, tran2])
                switchs.append(switch)
                count += 1

        trans_good = []
        tols_good = []
        for id in range(len(tols)):
            if (tols[id, 0] < d_tol or tols[id, 1] < f_tol) and tols[id, 2] < self.a_tol:
                if switchs[id]:
                    trans[id].extend([[[1, 0, 0], [0, -1, 0], [0, 0, -1]]])
                # print(tols[id], len(trans[id]))
                trans_good.append(trans[id])
                tols_good.append(tols[id])

        return trans_good, tols_good

    def search_transformation(self, lat_ref, d_tol=1.0, f_tol=0.1):
        """
        search the closest match to the reference lattice object

        Args:
            lat_ref: reference lattice object
            d_tol: tolerance in angle
            f_tol:
            a_tol:

        Returns:
            a two steps of transformation matrix if the match is possible
        """
        # Find all possible permutation and transformation matrices
        trans1 = self.get_permutation_matrices()
        trans2 = self.get_transformation_matrices()

        tols = np.zeros([len(trans2) * len(trans1) + 1, 3])
        trans = []
        switchs = []

        # Check it self
        d_tol1, f_tol1, a_tol1, switch = self.get_diff(lat_ref)
        tols[0] = [d_tol1, f_tol1, a_tol1]
        switchs.append(switch)
        trans.append([np.eye(3)])

        count = 0
        for _i, tran1 in enumerate(trans1):
            lat0 = self.transform(tran1)
            for _j, tran2 in enumerate(trans2):
                count += 1
                tmp = np.dot(tran2, lat0.matrix)
                try:
                    # print(i, j, self.ltype)
                    lat2 = Lattice.from_matrix(tmp, ltype=self.ltype)
                    d_tol1, f_tol1, a_tol1, switch = lat2.get_diff(lat_ref)
                    # print(d_tol1, f_tol1, a_tol1, switch)
                except:
                    d_tol1, f_tol1, a_tol1, switch = 10, 1.0, 90, None
                tols[count] = [d_tol1, f_tol1, a_tol1]
                trans.append([tran1, tran2])
                switchs.append(switch)

        # QZ: needs to figure out a better way to select the best
        rms = tols.sum(axis=1)
        ids = np.argsort(rms)
        id = ids[0]
        # print(tols, rms)
        # print(id, switchs[id])
        if abs(rms[ids[0]] - rms[ids[1]]) < 1e-3 and switchs[ids[0]] and not switchs[ids[1]]:
            id = ids[1]
            # print("change id 1", id)
        if id != 0 and abs(rms[0] - rms[id]) < 1.0:
            # print("change id 2", id, rms[0], rms[id])
            id = 0

        if (tols[id, 0] < d_tol or tols[id, 1] < f_tol) and tols[id, 2] < self.a_tol:
            if switchs[id]:
                trans[id].append([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
                return trans[id], tols[id]
            else:
                return trans[id], tols[id]
        else:
            # print("===================================Cannot match:", tols[id])
            return None, None

    def optimize_once(self, reset=False):
        """
        Optimize the lattice's inclination angles
        """
        opt = False
        trans = self.get_transformation_matrices()
        if len(trans) > 1:
            diffs = []
            for tran in trans:
                cell_new = np.dot(tran, self.matrix)
                try:
                    lat_new = Lattice.from_matrix(cell_new, ltype=self.ltype)
                    diffs.append(lat_new.get_worst_angle())
                except:
                    diffs.append(100)
            id = np.array(diffs).argmin()
            if id > 0 and diffs[id] < diffs[0] - 0.01:
                opt = True
                tran = trans[id]
                cell = np.dot(tran, self.matrix)
                lat = Lattice.from_matrix(cell, ltype=self.ltype, reset=reset)
                return lat, tran, opt
        return self, np.eye(3), opt

    def get_worst_angle(self):
        """
        return the worst inclination angle difference w.r.t 90 degree
        """
        return np.max(abs(np.array([self.alpha, self.beta, self.gamma]) - np.pi / 2))

    def optimize_multi(self, iterations=5):
        """
        Optimize the lattice if the cell has a bad inclination angles

        Args:
            iterations: maximum number of iterations
            force: whether or not do the early termination

        Returns:
            the optimized lattice
        """
        lattice = self
        trans_matrices = []
        for _i in range(iterations):
            lattice, trans, opt = lattice.optimize_once(reset=True)
            if opt:
                trans_matrices.append(trans)
            else:
                break
        return lattice, trans_matrices

    def standardize(self):
        """
        Force the angle to be smaller than 90 degree
        """
        change = False
        if self.ltype in ["monoclinic"]:
            if self.beta > np.pi / 2:
                self.beta = np.pi - self.beta
                change = True
        elif self.ltype in ["triclinic"]:
            if self.alpha > np.pi / 2:
                self.alpha = np.pi - self.alpha
                change = True
            if self.beta > np.pi / 2:
                self.beta = np.pi - self.beta
                change = True
            if self.gamma > np.pi / 2:
                self.gamma = np.pi - self.gamma
                change = True

        if change:
            para = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
            self.matrix = para2matrix(para)

    def transform(self, trans_mat=None, reset=False):
        """
        Optimize the lattice's inclination angles
        If reset is False, may return negative lattice
        """
        if trans_mat is None:
            trans_mat = np.eye(3)

        if not isinstance(trans_mat, np.ndarray):
            trans_mat = np.array(trans_mat)

        cell = np.dot(trans_mat, self.matrix)
        return Lattice.from_matrix(cell, ltype=self.ltype, reset=reset)

    def transform_multi(self, trans, reset=True):
        """
        Optimize the lattice's inclination angles
        """
        lat = self
        for tran in trans:
            lat = lat.transform(tran, reset)
        return lat

    def encode(self):
        a, b, c, alpha, beta, gamma = self.get_para(degree=True)
        if self.ltype in ["cubic"]:
            return [a]
        elif self.ltype in ["hexagonal", "trigonal", "tetragonal"]:
            return [a, c]
        elif self.ltype in ["orthorhombic"]:
            return [a, b, c]
        elif self.ltype in ["monoclinic"]:
            return [a, b, c, beta]
        else:
            return [a, b, c, alpha, beta, gamma]

    @classmethod
    def from_1d_representation(self, v, ltype):
        # print('test', v, type(v), len(v), v[0])
        if ltype == "triclinic":
            a, b, c, alpha, beta, gamma = v[0], v[1], v[2], v[3], v[4], v[5]
        elif ltype == "monoclinic":
            a, b, c, alpha, beta, gamma = v[0], v[1], v[2], 90, v[3], 90
        elif ltype == "orthorhombic":
            a, b, c, alpha, beta, gamma = v[0], v[1], v[2], 90, 90, 90
        elif ltype == "tetragonal":
            a, b, c, alpha, beta, gamma = v[0], v[0], v[1], 90, 90, 90
        elif ltype in ["trigonal", "hexagonal"]:
            a, b, c, alpha, beta, gamma = v[0], v[0], v[1], 90, 90, 120
        else:
            a, b, c, alpha, beta, gamma = v[0], v[0], v[0], 90, 90, 90
        try:
            return Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
        except:
            print(a, b, c, alpha, beta, gamma, ltype)

    def update_from_1d_representation(self, v):
        """
        Update the cell para and matrix from the 1d rep
        """
        if self.ltype == "triclinic":
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma = v[:6]
        elif self.ltype == "monoclinic":
            self.a, self.b, self.c, self.beta = v[:4]
        elif self.ltype == "orthorhombic":
            self.a, self.b, self.c = v[:3]
        elif self.ltype in ["tetragonal", "trigonal", "hexagonal"]:
            self.a, self.b, self.c = v[0], v[0], v[1]
        else:
            self.a, self.b, self.c = v[0], v[0], v[0]

        para = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        matrix = para2matrix(para)
        if matrix is not None:
            self.set_matrix(matrix)
        else:
            msg = f'error input {v} in update_from_1d_representation {self.ltype}'
            raise ValueError(msg)

    def mutate(self, degree=0.10, frozen=False):
        """
        Mutate the lattice object
        """
        rand = 1 + degree * (self.random_state.random(6) - 0.5)
        a0, b0, c0, alpha0, beta0, gamma0 = self.get_para()
        a = a0 * rand[0]
        b = b0 * rand[1]
        c = c0 * rand[2]
        alpha = np.degrees(alpha0 * rand[3])
        beta = np.degrees(beta0 * rand[4])
        gamma = np.degrees(gamma0 * rand[5])
        ltype = self.ltype

        if self.ltype in ["cubic"]:
            if frozen:
                lat = Lattice.from_para(a0, a0, a0, 90, 90, 90, ltype=ltype)
            else:
                lat = Lattice.from_para(a, a, a, 90, 90, 90, ltype=ltype)
        elif ltype in ["hexagonal", "trigonal"]:
            if frozen:
                lat = Lattice.from_para(a0, a0, c, 90, 90, 120, ltype=ltype)
            else:
                lat = Lattice.from_para(a, a, c, 90, 90, 120, ltype=ltype)
        elif ltype in ["tetragonal"]:
            if frozen:
                lat = Lattice.from_para(a0, a0, c, 90, 90, 90, ltype=ltype)
            else:
                lat = Lattice.from_para(a, a, c, 90, 90, 90, ltype=ltype)
        elif ltype in ["orthorhombic"]:
            lat = Lattice.from_para(a, b, c, 90, 90, 90, ltype=ltype)
        elif ltype in ["monoclinic"]:
            lat = Lattice.from_para(a, b, c, 90, beta, 90, ltype=ltype)
        elif ltype in ["triclinic"]:
            lat = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
        else:
            raise ValueError(f"ltype {ltype:s} is not supported")
        return lat

    def generate_para(self):
        if self.dim == 3:
            return generate_cellpara(self.ltype, self.volume, random_state=self.random_state, **self.kwargs)
        elif self.dim == 2:
            return generate_cellpara_2D(self.ltype, self.volume, random_state=self.random_state, **self.kwargs)
        elif self.dim == 1:
            return generate_cellpara_1D(self.ltype, self.volume, random_state=self.random_state, **self.kwargs)
        elif self.dim == 0:
            return generate_cellpara_0D(self.ltype, self.volume, random_state=self.random_state, **self.kwargs)
        return None

    def generate_matrix(self):
        """
        Generates a 3x3 matrix for a lattice based on the lattice type and volume
        """
        # Try multiple times in case of failure
        for _i in range(10):
            para = self.generate_para()
            if para is not None:
                return para2matrix(para)
        return None

    def get_matrix(self, shape="upper"):
        """
        Returns a 3x3 numpy array representing the lattice vectors.
        """
        return self.matrix

    def get_para(self, degree=False):
        """
        Returns a tuple of lattice parameters.
        """
        if degree:
            return (
                self.a,
                self.b,
                self.c,
                deg * self.alpha,
                deg * self.beta,
                deg * self.gamma,
            )
        else:
            return (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    def set_matrix(self, matrix=None):
        if matrix is not None:
            m = np.array(matrix)
            if np.shape(m) == (3, 3):
                self.matrix = m
                try:
                    self.inv_matrix = np.linalg.inv(m)
                except:
                    print(self.para)
                    print(matrix)
                    msg = "Error in getting the inv_matrix"
                    raise ValueError(msg)
            else:
                print(matrix)
                msg = "Error: matrix must be a 3x3 numpy array or list"
                raise ValueError(msg)
        else:
            self.reset_matrix()

        para = matrix2para(self.matrix)
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = para
        self.volume = np.linalg.det(self.matrix)

    def set_para(self, para=None, radians=False):
        if para is not None:
            if radians is False:
                para[3] *= rad
                para[4] *= rad
                para[5] *= rad
            self.set_matrix(para2matrix(para))
        else:
            self.set_matrix()

    def update_para(self, id, change):
        para = [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]
        para[id] += change
        self.set_matrix(para2matrix(para))

    def reset_matrix(self, shape="upper"):
        if self.random:
            success = False
            for _i in range(5):
                m = self.generate_matrix()
                if m is not None:
                    self.matrix = m
                    self.inv_matrix = np.linalg.inv(m)
                    [a, b, c, alpha, beta, gamma] = matrix2para(self.matrix)
                    self.a = a
                    self.b = b
                    self.c = c
                    self.alpha = alpha
                    self.beta = beta
                    self.gamma = gamma
                    success = True
                    break
            if not success:
                msg = "Cannot generate a good matrix"
                raise ValueError(msg)
        else:
            # a small utility to convert the cell shape
            para = matrix2para(self.matrix)
            self.matrix = para2matrix(para, format=shape)
            self.inv_matrix = np.linalg.inv(self.matrix)

    def set_volume(self, volume):
        if self.allow_volume_reset:
            self.volume = volume

    def swap_axis(self, random=False, ids=None):
        """
        For the lattice
        """
        # only applied to triclinic/monoclinic/orthorhombic
        if self.ltype in ["triclinic", "orthorhombic", "Orthorhombic"]:
            allowed_ids = [
                [0, 1, 2],
                [1, 0, 2],
                [0, 2, 1],
                [2, 1, 0],
                [1, 2, 0],
                [2, 0, 1],
            ]

        elif self.ltype in ["monoclinic"]:
            if abs(self.beta - 90 * rad) > 1e-3:
                allowed_ids = [[0, 1, 2], [2, 1, 0]]
            else:
                allowed_ids = [
                    [0, 1, 2],
                    [1, 0, 2],
                    [0, 2, 1],
                    [2, 1, 0],
                    [1, 2, 0],
                    [2, 0, 1],
                ]
        else:
            allowed_ids = [[0, 1, 2]]

        if random:
            ids = self.random_state.choice(allowed_ids)
        else:
            if ids not in allowed_ids:
                print(ids)
                raise ValueError("the above swap is not allowed in " + self.ltype)

        (a, b, c, alpha, beta, gamma) = self.get_para()
        alpha, beta, gamma = alpha * deg, beta * deg, gamma * deg
        if ids is None:
            return self
        elif ids == [1, 0, 2]:  # a->b
            return self.from_para(b, a, c, beta, alpha, gamma, self.ltype)
        elif ids == [2, 1, 0]:  # a->c
            return self.from_para(c, b, a, gamma, beta, alpha, self.ltype)
        elif ids == [0, 2, 1]:  # b-c
            return self.from_para(a, c, b, alpha, gamma, beta, self.ltype)
        elif ids == [2, 0, 1]:
            return self.from_para(c, a, b, gamma, alpha, beta, self.ltype)
        elif ids == [1, 2, 0]:
            return self.from_para(b, c, a, beta, gamma, alpha, self.ltype)
        else:
            return self

    def swap_angle(self, random=True, ids=None):
        # only applied to triclinic/monoclinic #/hexagonal
        """
        If the angle is not 90. There will be two equivalent versions
        e.g., 80 and 100.
        """
        if self.ltype in ["monoclinic"]:
            allowed_ids = ["beta", "No"]
        elif self.ltype in ["triclinic"]:
            allowed_ids = ["alpha", "beta", "gamma", "No"]
        else:
            allowed_ids = ["No"]

        if random:
            ids = self.random_state.choice(allowed_ids)
        else:
            if ids not in allowed_ids:
                print(ids)
                raise ValueError("the above swap is not allowed in " + self.ltype)

        (a, b, c, alpha, beta, gamma) = self.get_para()
        alpha, beta, gamma = alpha * deg, beta * deg, gamma * deg
        if ids is None:
            return self
        elif ids == "alpha":
            return self.from_para(a, b, c, 180 - alpha, beta, gamma, self.ltype)
        elif ids == "beta":
            return self.from_para(a, b, c, alpha, 180 - beta, gamma, self.ltype)
        elif ids == "gamma":
            return self.from_para(a, b, c, alpha, beta, 180 - gamma, self.ltype)
        else:
            return self

    def add_vacuum(self, coor, frac=True, vacuum=15, PBC=None):
        """
        Adds space above and below a 2D or 1D crystal.

        Args:
            coor: the relative coordinates of the crystal
            vacuum: the amount of space, in Angstroms, to add above and below
            PBC: A periodic boundary condition list,
                Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity
                along the z axis

        Returns:
            The transformed lattice and coordinates after the vacuum is added
        """
        if PBC is None:
            PBC = [0, 0, 0]
        matrix = self.matrix
        absolute_coords = np.dot(coor, matrix) if frac else coor

        for i, a in enumerate(PBC):
            if not a:
                ratio = 1 + vacuum / np.linalg.norm(matrix[i])
                matrix[i] *= ratio
                absolute_coords[:, i] += vacuum / 2
        coor = np.dot(absolute_coords, np.linalg.inv(matrix)) if frac else absolute_coords
        return matrix, coor

    def generate_point(self):
        if self.ltype in ["spherical", "ellipsoidal"]:
            # Choose a point within an octant of the unit sphere using Marsaglia's method
            # Generate 3D vector from normal distribution
            vec = self.random_state.normal(0, 1, 3)
            # Generate random radius
            r = self.random_state.random() ** (1 / 3)
            # Normalize and scale
            point = r * vec / np.linalg.norm(vec)
        else:
            point = self.random_state.random(3)

            for i, a in enumerate(self.PBC):
                if not a:
                    if self.ltype in ["hexagonal", "trigonal"]:
                        point[i] *= 1.0 / np.sqrt(3.0)
                    else:
                        point[i] -= 0.5
        return point

    @classmethod
    def from_para(
        self,
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        ltype="triclinic",
        radians=False,
        PBC=None,
        factor=1.0,
        force_symmetry=False,
        **kwargs,
    ):
        """
        Creates a Lattice object from 6 lattice parameters. Additional keyword
        arguments  are available. Unless specified by the keyword random=True,
        does not create a new matrix upon calling reset_matrix. This allows
        for generation of random crystals with a specific choice of unit cell.

        Args:
            a: The length (in Angstroms) of the unit cell vectors
            b: The length (in Angstroms) of the unit cell vectors
            c: The length (in Angstroms) of the unit cell vectors
            alpha: the angle (in degrees) between the b and c vectors
            beta: the angle (in degrees) between the a and c vectors
            gamma: the angle (in degrees) between the a and b vectors
            ltype: the lattice type ("cubic, tetragonal, etc."). Also available
                are "spherical", which confines generated points to lie within a
                sphere, and "ellipsoidal", which confines generated points to lie
                within an ellipse (oriented about the z axis)
            radians: whether or not to use radians (instead of degrees) for the
                lattice angles
            PBC: A periodic boundary condition list, where 1 means periodic,
                0 means not periodic.
                Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity along
                the z axis
            kwargs: various values which may be defined. If none are defined,
                random ones will be generated. Values will be passed to
                generate_lattice. Options include:
                area: The cross-sectional area (in Angstroms squared). Only used
                    to generate 1D crystals
                thickness: The unit cell's non-periodic thickness (in Angstroms).
                    Only used to generate 2D crystals
                unique_axis: The unique axis for certain symmetry (and especially
                    layer) groups. Because the symmetry operations are not also
                    transformed, you should use the default values for random
                    crystal generation
                random: If False, keeps the stored values for the lattice geometry
                    even upon applying reset_matrix. To alter the matrix,
                    use set_matrix() or set_para
                'unique_axis': the axis ('a', 'b', or 'c') which is unique
                'min_l': the smallest allowed cell vector.
                'mid_l': the second smallest allowed cell vector.
                'max_l': the third smallest allowed cell vector.

        Returns:
            a Lattice object with the specified parameters
        """
        if PBC is None:
            PBC = [1, 1, 1]
        try:
            cell_matrix = para2matrix((a, b, c, alpha, beta, gamma), radians=radians)
            cell_matrix *= factor
        except Exception as err:
            msg = "Error: invalid cell parameters for lattice."
            raise ValueError(msg) from err

        if force_symmetry:
            return Lattice.from_matrix(cell_matrix, ltype=ltype)
        else:
            volume = np.linalg.det(cell_matrix)
            # Initialize a Lattice instance
            l = Lattice(ltype, volume, PBC=PBC, **kwargs)
            l.a, l.b, l.c = factor * a, factor * b, factor * c
            l.alpha, l.beta, l.gamma = alpha * rad, beta * rad, gamma * rad
            l.matrix = cell_matrix
            l.inv_matrix = np.linalg.inv(cell_matrix)
            l.ltype = ltype
            l.volume = volume
            l.random = False
            l.allow_volume_reset = False
            return l

    @classmethod
    def from_matrix(
        self,
        matrix,
        reset=True,
        shape="upper",
        ltype="triclinic",
        PBC=None,
        **kwargs,
    ):
        """
        Creates a Lattice object from a 3x3 cell matrix. Additional keywords
        are available. Unless specified by the keyword random=True, does not
        create a new matrix upon calling reset_matrix. This allows for a random
        crystals with a specific choice of unit cell.

        Args:
            matrix: a 3x3 real matrix (numpy array or nested list) for the cell
            ltype: the lattice type ("cubic, tetragonal, etc."). Also can be
                - "spherical", confines points to lie within a sphere,
                - "ellipsoidal", points to lie within an ellipsoid (about z axis)
            PBC: A periodic boundary condition list, where 1 is periodic
                Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity at z.
            kwargs: various values which may be defined. Random ones if None
                Values will be passed to generate_lattice. Options include:
                `area: The cross-sectional area (in Ang^2) for 1D crystals
                `thickness`: The cell's thickness (in Ang) for 2D crystals
                `unique_axis`: The unique axis for layer groups.
                `random`: If False, keeps the stored values for the lattice
                geometry even applying reset_matrix. To alter the matrix,
                use `set_matrix()` or `set_para`
                'unique_axis': the axis ('a', 'b', or 'c') which is unique.
                'min_l': the smallest allowed cell vector.
                'mid_l': the second smallest allowed cell vector.
                'max_l': the third smallest allowed cell vector.

        Returns:
            a Lattice object with the specified parameters
        """
        if PBC is None:
            PBC = [1, 1, 1]

        m = np.array(matrix)
        if np.shape(m) != (3, 3):
            print(matrix)
            msg = "Error: matrix must be a 3x3 numpy array or list"
            raise ValueError(msg)

        [a, b, c, alpha, beta, gamma] = matrix2para(m)

        # symmetrize the lattice
        if reset:
            if ltype in ["cubic", "Cubic"]:
                a = b = c = (a + b + c) / 3
                alpha = beta = gamma = np.pi / 2
            elif ltype in ["hexagonal", "trigonal", "Hexagonal", "Trigonal"]:
                a = b = (a + b) / 2
                alpha = beta = np.pi / 2
                gamma = np.pi * 2 / 3
            elif ltype in ["tetragonal", "Tetragonal"]:
                a = b = (a + b) / 2
                alpha = beta = gamma = np.pi / 2
            elif ltype in ["orthorhombic", "Orthorhombic"]:
                alpha = beta = gamma = np.pi / 2
            elif ltype in ["monoclinic", "Monoclinic"]:
                alpha = gamma = np.pi / 2

            # reset matrix according to the symmetry
            m = para2matrix([a, b, c, alpha, beta, gamma], format=shape)

        # Initialize a Lattice instance
        volume = np.linalg.det(m)
        l = Lattice(ltype, volume, m, PBC=PBC, **kwargs)
        l.a, l.b, l.c = a, b, c
        l.alpha, l.beta, l.gamma = alpha, beta, gamma
        l.matrix = m
        l.inv_matrix = np.linalg.inv(m)
        l.ltype = ltype
        l.volume = volume
        l.random = False
        l.allow_volume_reset = False
        return l

    def is_valid_matrix(self):
        """
        check if the cell parameter is reasonable or not
        """

        try:
            paras = [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]
            para2matrix(paras)
            return True
        except:
            return False

    def is_valid_lattice(self, tol=1e-3):
        ltype = self.ltype.lower()

        def check_angles(angles):
            return all(abs(angle - np.pi / 2) <= tol for angle in angles)

        if ltype == "cubic":
            return (
                abs(self.a - self.b) <= tol
                and abs(self.a - self.c) <= tol
                and check_angles([self.alpha, self.beta, self.gamma])
            )

        elif ltype in ["hexagonal", "trigonal"]:
            return (
                abs(self.a - self.b) <= tol
                and check_angles([self.alpha, self.beta])
                and abs(self.gamma - 2 / 3 * np.pi) <= tol
            )

        elif ltype == "tetragonal":
            return abs(self.a - self.b) <= tol and check_angles([self.alpha, self.beta, self.gamma])

        elif ltype == "orthorhombic":
            return check_angles([self.alpha, self.beta, self.gamma])

        elif ltype == "monoclinic":
            return check_angles([self.alpha, self.gamma])

        return True

    def check_mismatch(self, trans, l_type, tol=1.0, a_tol=10):
        """
        check if the lattice mismatch is big after a transformation
        This is mostly used in supergroup function
        QZ: to fix ===============

        Args:
            trans: 3*3 matrix
            l_type: lattice_type like orthrhombic
            tol: tolerance in a, b, c
            a_tol: tolerance in alpha, beta, gamma

        Returns:
            True or False
        """
        matrix = np.dot(trans.T, self.matrix)
        l1 = Lattice.from_matrix(matrix)
        l2 = Lattice.from_matrix(matrix, ltype=l_type)
        (a1, b1, c1, alpha1, beta1, gamma1) = l1.get_para(degree=True)
        (a2, b2, c2, alpha2, beta2, gamma2) = l2.get_para(degree=True)
        abc_diff = np.abs(np.array([a2 - a1, b2 - b1, c2 - c1])).max()
        ang_diff = np.abs(np.array([alpha2 - alpha1, beta2 - beta1, gamma2 - gamma1])).max()
        return not (abc_diff > tol or ang_diff > a_tol)

    def get_diff(self, l_ref):
        """
        get the difference in length, angle, and check if switch is needed
        """
        (a1, b1, c1, alpha1, beta1, gamma1) = self.get_para(degree=True)
        (a2, b2, c2, alpha2, beta2, gamma2) = l_ref.get_para(degree=True)
        abc_diff = np.abs(np.array([a2 - a1, b2 - b1, c2 - c1])).max()
        abc_f_diff = np.abs(np.array([(a2 - a1) / a1, (b2 - b1) / b1, (c2 - c1) / c1])).max()
        ang_diff1 = abs(alpha1 - alpha2) + abs(beta1 - beta2) + abs(gamma1 - gamma2)
        ang_diff2 = abs(alpha1 - alpha2)
        ang_diff2 += abs(abs(beta1 - 90) - abs(beta2 - 90))
        ang_diff2 += abs(gamma1 - gamma2)
        # print(abc_diff, abc_f_diff, ang_diff1, ang_diff2, self.ltype)
        if ang_diff1 < ang_diff2 + 0.01:
            return abc_diff, abc_f_diff, ang_diff1, False
        else:
            if self.ltype == "monoclinic":
                return abc_diff, abc_f_diff, ang_diff2, True
            else:
                return abc_diff, abc_f_diff, ang_diff2, False

    def __str__(self):
        return (
            f"{self.a:8.4f}, {self.b:8.4f}, {self.c:8.4f}, {self.alpha * deg:8.4f}, "
            f"{self.beta * deg:8.4f}, {self.gamma * deg:8.4f}, {self.ltype!s:s}"
        )

    def __repr__(self):
        return str(self)

    def find_transition_to_orthoslab(self, c=(0, 0, 1), a=(1, 0, 0), m=5):
        """
        Create the slab model with an approximate orthogonal box shape
        """
        from pyxtal.plane import has_reduction

        tol = 1e-3
        direction = np.array(c)

        # find the simplest a-direction
        if np.dot(np.array(a), direction) < tol:
            a_hkl = np.array(a)
        else:
            a_hkls = []
            for h in range(-m, m + 1):
                for k in range(-m, m + 1):
                    for l in range(-m, m + 1):
                        hkl = np.array([h, k, l])
                        if (
                            ([h, k, l] != [0, 0, 0])
                            and (not has_reduction(hkl))
                            and (abs(np.dot(hkl, direction)) < tol)
                        ):
                            a_hkls.append(hkl)
            a_hkls = np.array(a_hkls)  # ; print(a_hkls)
            a_hkl = a_hkls[np.argmin(np.abs(a_hkls).sum(axis=1))]
        a_vector = np.dot(a_hkl, self.matrix)
        # print('a_hkl', a_hkl)

        # find the simplest b-direction
        b_hkl = None
        min_angle_ab = float("inf")
        for h in range(-m, m + 1):
            for k in range(-m, m + 1):
                for l in range(-m, m + 1):
                    hkl = np.array([h, k, l])
                    if ([h, k, l] != [0, 0, 0]) and (not has_reduction(hkl)) and (abs(np.dot(hkl, direction)) < tol):
                        vector = np.dot(hkl, self.matrix)
                        angle1 = angle(vector, a_vector, radians=False)
                        if abs(90 - angle1) < min_angle_ab:
                            min_angle_ab = abs(90 - angle1)
                            b_hkl = hkl
                            b_vector = vector

        # print('b_hkl', b_hkl, min_angle_ab)
        # change the sign
        if abs(angle(np.cross(a_hkl, b_hkl), direction)) > tol:
            b_hkl *= -1
            b_vector *= -1

        ## update the c_direction
        ab_plane = np.cross(a_vector, b_vector)  # ; print('ab_plane', ab_plane)
        c_hkl = None
        min_angle_c = float("inf")
        for h in range(-m, m + 1):
            for k in range(-m, m + 1):
                for l in range(-m, m + 1):
                    hkl = np.array([h, k, l])
                    if [h, k, l] != [0, 0, 0] and not has_reduction(hkl):
                        vector = np.dot(hkl, self.matrix)
                        angle1 = angle(vector, ab_plane, radians=False)
                        # print(hkl, angle)
                        if abs(angle1) < abs(min_angle_c):
                            min_angle_c = angle1
                            c_hkl = hkl

        # print(a_hkl, b_hkl, c_hkl)
        return np.vstack([a_hkl, b_hkl, c_hkl])

    def apply_transformation(self, trans):
        """
        Optimize the lattice's inclination angles
        """
        cell_new = np.dot(trans, self.matrix)
        return Lattice.from_matrix(cell_new)


def generate_cellpara(
    ltype,
    volume,
    minvec=1.2,
    minangle=np.pi / 6,
    max_ratio=10.0,
    min_special=None,
    maxattempts=100,
    random_state: None | int | Generator = None,
    **kwargs,
):
    """
    Generates the cell parameter (a, b, c, alpha, beta, gamma) according
    to the space group symmetry and number of atoms. If the spacegroup
    has centering, we will transform to conventional cell setting. If the
    generated lattice does not meet the minimum angle and vector
    requirements, we try to generate a new one, up to maxattempts times.

    Args:
        volume: volume of the conventional unit cell
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice
        kwargs: a dictionary of optional values. These include:
            'unique_axis': the axis ('a', 'b', or 'c') which is unique.
            'min_l': the smallest allowed cell vector.
            'mid_l': the second smallest allowed cell vector.
            'max_l': the third smallest allowed cell vector.

    Returns:
        a 6-length array representing the lattice of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    if isinstance(random_state, Generator):
        random_state = random_state.spawn(1)[0]
    else:
        random_state = np.random.default_rng(random_state)

    min_special = kwargs.get("min_special", min_special)  # ; print("min_special", min_special)
    maxangle = np.pi - minangle
    for _n in range(maxattempts):
        # Triclinic
        # if sg <= 2:
        if ltype == "triclinic":
            # Derive lattice constants from a random matrix
            mat = random_shear_matrix(width=0.2, random_state=random_state)
            a, b, c, alpha, beta, gamma = matrix2para(mat)
            x = np.sqrt(
                1
                - np.cos(alpha) ** 2
                - np.cos(beta) ** 2
                - np.cos(gamma) ** 2
                + 2 * (np.cos(alpha) * np.cos(beta) * np.cos(gamma))
            )
            vec = random_vector(random_state=random_state)
            abc = volume / x
            xyz = vec[0] * vec[1] * vec[2]
        # Monoclinic
        elif ltype in ["monoclinic"]:
            alpha, gamma = np.pi / 2, np.pi / 2
            beta = gaussian_random_variable(minangle, maxangle, random_state=random_state)
            x = np.sin(beta)
            vec = random_vector(random_state=random_state)
            xyz = vec[0] * vec[1] * vec[2]
            abc = volume / x
            xyz = vec[0] * vec[1] * vec[2]
        # Orthorhombic
        # elif sg <= 74:
        elif ltype in ["orthorhombic"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
            x = 1
            vec = random_vector(random_state=random_state)
            xyz = vec[0] * vec[1] * vec[2]
            abc = volume / x
        # Tetragonal
        # elif sg <= 142:
        elif ltype in ["tetragonal"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
            x = 1
            vec = random_vector(random_state=random_state)
            c = vec[2] / (vec[0] * vec[1]) * np.cbrt(volume / x)
            a = b = np.sqrt((volume / x) / c)
        # Trigonal/Rhombohedral/Hexagonal
        # elif sg <= 194:
        elif ltype in ["hexagonal", "trigonal"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 3 * 2
            x = np.sqrt(3.0) / 2.0
            vec = random_vector(random_state=random_state)
            c = vec[2] / (vec[0] * vec[1]) * np.cbrt(volume / x)
            a = b = np.sqrt((volume / x) / c)
        # Cubic
        # else:
        elif ltype in ["cubic"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
            s = (volume) ** (1.0 / 3.0)
            a, b, c = s, s, s

        # resort a/b/c if min_special is not None for mol. xtals
        if ltype in ["triclinic", "monoclinic", "orthorhombic"]:
            vec *= np.cbrt(abc) / np.cbrt(xyz)
            if min_special is not None:
                ax = random_state.choice([0, 1, 2])
                if vec[ax] < min_special:
                    coef = random_state.uniform(0.8, 1.2) * min_special / vec[ax]
                    for i in range(3):
                        if i == ax:
                            vec[i] *= coef
                        else:
                            vec[i] /= np.sqrt(coef)
            [a, b, c] = vec

        # Check that lattice meets requirements
        maxvec = (a * b * c) / (minvec**2)

        # Define limits on cell dimensions
        min_l = kwargs.get("min_l", minvec)
        mid_l = kwargs.get("mid_l", min_l)
        max_l = kwargs.get("max_l", mid_l)
        l_min = min(a, b, c)
        l_max = max(a, b, c)
        for x in (a, b, c):
            if x <= l_max and x >= l_min:
                l_mid = x
        if not (l_min >= min_l and l_mid >= mid_l and l_max >= max_l):
            continue

        if minvec < maxvec:
            # Check minimum Euclidean distances
            smallvec = min(
                a * np.cos(max(beta, gamma)),
                b * np.cos(max(alpha, gamma)),
                c * np.cos(max(alpha, beta)),
            )
            if (
                a > minvec
                and b > minvec
                and c > minvec
                and a < maxvec
                and b < maxvec
                and c < maxvec
                and smallvec < minvec
                and alpha > minangle
                and beta > minangle
                and gamma > minangle
                and alpha < maxangle
                and beta < maxangle
                and gamma < maxangle
                and a / b < max_ratio
                and a / c < max_ratio
                and b / c < max_ratio
                and b / a < max_ratio
                and c / a < max_ratio
                and c / b < max_ratio
            ):
                return np.array([a, b, c, alpha, beta, gamma])

    # If maxattempts tries have been made without success
    msg = f"lattice fails after {maxattempts:d} cycles"
    msg += f"for volume {volume:.2f}"
    raise VolumeError(msg)
    # return


def generate_cellpara_2D(
    ltype,
    volume,
    thickness=None,
    minvec=1.2,
    minangle=np.pi / 6,
    max_ratio=10.0,
    maxattempts=100,
    random_state: None | int | Generator = None,
    **kwargs,
):
    """
    Generates the cell parameter (a, b, c, alpha, beta, gamma) according
    to the layer group symmetry and number of atoms. If the layer group
    has centering, we will transform to conventional cell setting. If the
    generated lattice does not meet the minimum angle and vector
    requirements, we try to generate a new one, up to maxattempts times.

    Note: The monoclinic layer groups have different unique axes. Groups 3-7
        have unique axis c, while 8-18 have unique axis a. We use non-periodic
        axis c for all layer groups.

    Args:
        num: International number of the space group
        volume: volume of the lattice
        thickness: 3rd-dimensional thickness of the unit cell. If set to None,
            a thickness is chosen automatically
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice
        kwargs: a dictionary of optional values. These include:
            'unique_axis': the axis ('a', 'b', or 'c') which is unique.
            'min_l': the smallest allowed cell vector.
            'mid_l': the second smallest allowed cell vector.
            'max_l': the third smallest allowed cell vector.

    Returns:
        a 6-length representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    if isinstance(random_state, int):
        # NOTE if random_state is an integer make a Generator to ensure randomness
        # downstream that would be lost if integer seed used repeated
        random_state = np.random.default_rng(random_state)

    unique_axis = kwargs.get("unique_axis", "c")
    # Store the non-periodic axis
    NPA = 3
    # Set the unique axis for monoclinic cells
    # if num in range(3, 8): unique_axis = "c"
    # elif num in range(8, 19): unique_axis = "a"
    maxangle = np.pi - minangle
    for _n in range(maxattempts):
        abc = np.ones([3])
        if thickness is None:
            v = random_vector(random_state=random_state)
            thickness1 = np.cbrt(volume) * (v[0] / (v[0] * v[1] * v[2]))
        else:
            thickness1 = max([3.0, thickness])
        abc[NPA - 1] = thickness1
        alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
        # Triclinic
        # if num <= 2:
        if ltype == "triclinic":
            mat = random_shear_matrix(width=0.2, random_state=random_state)
            a, b, c, alpha, beta, gamma = matrix2para(mat)
            x = np.sqrt(
                1
                - np.cos(alpha) ** 2
                - np.cos(beta) ** 2
                - np.cos(gamma) ** 2
                + 2 * (np.cos(alpha) * np.cos(beta) * np.cos(gamma))
            )
            abc[NPA - 1] = abc[NPA - 1] / x  # scale thickness by outer product of vectors
            ab = volume / (abc[NPA - 1] * x)
            ratio = a / b
            if NPA == 3:
                abc[0] = np.sqrt(ab * ratio)
                abc[1] = np.sqrt(ab / ratio)
            elif NPA == 2:
                abc[0] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)
            elif NPA == 1:
                abc[1] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)

        # Monoclinic
        # elif num <= 18:
        elif ltype == "monoclinic":
            a, b, c = random_vector(random_state=random_state)
            if unique_axis == "a":
                alpha = gaussian_random_variable(minangle, maxangle, random_state=random_state)
                x = np.sin(alpha)
            elif unique_axis == "b":
                beta = gaussian_random_variable(minangle, maxangle, random_state=random_state)
                x = np.sin(beta)
            elif unique_axis == "c":
                gamma = gaussian_random_variable(minangle, maxangle, random_state=random_state)
                x = np.sin(gamma)
            ab = volume / (abc[NPA - 1] * x)
            ratio = a / b
            if NPA == 3:
                abc[0] = np.sqrt(ab * ratio)
                abc[1] = np.sqrt(ab / ratio)
            elif NPA == 2:
                abc[0] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)
            elif NPA == 1:
                abc[1] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)

        # Orthorhombic
        # elif num <= 48:
        elif ltype == "orthorhombic":
            vec = random_vector(random_state=random_state)
            if NPA == 3:
                ratio = abs(vec[0] / vec[1])  # ratio a/b
                abc[1] = np.sqrt(volume / (thickness1 * ratio))
                abc[0] = abc[1] * ratio
            elif NPA == 2:
                ratio = abs(vec[0] / vec[2])  # ratio a/b
                abc[2] = np.sqrt(volume / (thickness1 * ratio))
                abc[0] = abc[2] * ratio
            elif NPA == 1:
                ratio = abs(vec[1] / vec[2])  # ratio a/b
                abc[2] = np.sqrt(volume / (thickness1 * ratio))
                abc[1] = abc[2] * ratio

        # Tetragonal
        # elif num <= 64:
        elif ltype == "tetragonal":
            if NPA == 3:
                abc[0] = abc[1] = np.sqrt(volume / thickness1)
            elif NPA == 2:
                abc[0] = abc[1]
                abc[2] = volume / (abc[NPA - 1] ** 2)
            elif NPA == 1:
                abc[1] = abc[0]
                abc[2] = volume / (abc[NPA - 1] ** 2)

        # Trigonal/Hexagonal
        # elif num <= 80:
        elif ltype in ["hexagonal", "trigonal"]:
            gamma = np.pi / 3 * 2
            x = np.sqrt(3.0) / 2.0
            if NPA == 3:
                abc[0] = abc[1] = np.sqrt((volume / x) / abc[NPA - 1])
            elif NPA == 2:
                abc[0] = abc[1]
                abc[2] = (volume / x)(thickness1**2)
            elif NPA == 1:
                abc[1] = abc[0]
                abc[2] = (volume / x) / (thickness1**2)

        para = np.array([abc[0], abc[1], abc[2], alpha, beta, gamma])

        a, b, c = abc[0], abc[1], abc[2]
        maxvec = (a * b * c) / (minvec**2)

        # Define limits on cell dimensions
        min_l = kwargs.get("min_l", minvec)
        mid_l = kwargs.get("mid_l", min_l)
        max_l = kwargs.get("max_l", mid_l)
        l_min = min(a, b, c)
        l_max = max(a, b, c)
        for x in (a, b, c):
            if x <= l_max and x >= l_min:
                l_mid = x
        if not (l_min >= min_l and l_mid >= mid_l and l_max >= max_l):
            continue

        if minvec < maxvec:
            smallvec = min(
                a * np.cos(max(beta, gamma)),
                b * np.cos(max(alpha, gamma)),
                c * np.cos(max(alpha, beta)),
            )
            if (
                a > minvec
                and b > minvec
                and c > minvec
                and a < maxvec
                and b < maxvec
                and c < maxvec
                and smallvec < minvec
                and alpha > minangle
                and beta > minangle
                and gamma > minangle
                and alpha < maxangle
                and beta < maxangle
                and gamma < maxangle
                and a / b < max_ratio
                and a / c < max_ratio
                and b / c < max_ratio
                and b / a < max_ratio
                and c / a < max_ratio
                and c / b < max_ratio
            ):
                return para

    # If maxattempts tries have been made without success
    msg = f"Cannot get lattice after {maxattempts:d} cycles for volume {volume:.2f}"
    raise VolumeError(msg)


def generate_cellpara_1D(
    ltype,
    volume,
    area=None,
    minvec=1.2,
    minangle=np.pi / 6,
    max_ratio=10.0,
    maxattempts=100,
    random_state: None | int | Generator = None,
    **kwargs,
):
    """
    Generates a cell parameter (a, b, c, alpha, beta, gamma) according to
    the rod group symmetry and number of atoms. If the rod group has centering,
    we will transform to conventional cell setting. If the generated lattice
    does not meet the minimum angle and vector requirements, we try to
    generate a new one, up to maxattempts times.

    Note: The monoclinic Rod groups have different unique axes. Groups 3-7
        have unique axis a, while 8-12 have unique axis c. We use periodic
        axis c for all Rod groups.

    Args:
        num: number of the Rod group
        volume: volume of the lattice
        area: cross-sectional area of the unit cell in Angstroms squared. If
            set to None, a value is chosen automatically
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice
        kwargs: a dictionary of optional values. These include:
            'unique_axis': the axis ('a', 'b', or 'c') which is unique.
            'min_l': the smallest allowed cell vector.
            'mid_l': the second smallest allowed cell vector.
            'max_l': the third smallest allowed cell vector.

    Returns:
        a 6-length array representing the lattice of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    if isinstance(random_state, int):
        # NOTE if random_state is an integer make a Generator to ensure randomness
        # downstream that would be lost if integer seed used repeated
        random_state = np.random.default_rng(random_state)

    try:
        unique_axis = kwargs["unique_axis"]
    except:
        unique_axis = "a"
    # Store the periodic axis
    PA = 3
    # Set the unique axis for monoclinic cells
    # if num in range(3, 8): unique_axis = "a"
    # elif num in range(8, 13): unique_axis = "c"
    maxangle = np.pi - minangle
    for _n in range(maxattempts):
        abc = np.ones([3])
        if area is None:
            v = random_vector(random_state=random_state)
            thickness1 = np.cbrt(volume) * (v[0] / (v[0] * v[1] * v[2]))
        else:
            thickness1 = volume / area
        abc[PA - 1] = thickness1
        alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
        # Triclinic
        # if num <= 2:
        if ltype == "triclinic":
            mat = random_shear_matrix(width=0.2, random_state=random_state)
            a, b, c, alpha, beta, gamma = matrix2para(mat)
            x = np.sqrt(
                1
                - np.cos(alpha) ** 2
                - np.cos(beta) ** 2
                - np.cos(gamma) ** 2
                + 2 * (np.cos(alpha) * np.cos(beta) * np.cos(gamma))
            )
            abc[PA - 1] = abc[PA - 1] / x  # scale thickness by outer product of vectors
            ab = volume / (abc[PA - 1] * x)
            ratio = a / b
            if PA == 3:
                abc[0] = np.sqrt(ab * ratio)
                abc[1] = np.sqrt(ab / ratio)
            elif PA == 2:
                abc[0] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)
            elif PA == 1:
                abc[1] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)

        # Monoclinic
        # elif num <= 12:
        elif ltype == "monoclinic":
            a, b, c = random_vector(random_state=random_state)
            if unique_axis == "a":
                alpha = gaussian_random_variable(minangle, maxangle, random_state=random_state)
                x = np.sin(alpha)
            elif unique_axis == "b":
                beta = gaussian_random_variable(minangle, maxangle, random_state=random_state)
                x = np.sin(beta)
            elif unique_axis == "c":
                gamma = gaussian_random_variable(minangle, maxangle, random_state=random_state)
                x = np.sin(gamma)
            ab = volume / (abc[PA - 1] * x)
            ratio = a / b
            if PA == 3:
                abc[0] = np.sqrt(ab * ratio)
                abc[1] = np.sqrt(ab / ratio)
            elif PA == 2:
                abc[0] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)
            elif PA == 1:
                abc[1] = np.sqrt(ab * ratio)
                abc[2] = np.sqrt(ab / ratio)

        # Orthorhombic
        # lif num <= 22:
        elif ltype == "orthorhombic":
            vec = random_vector(random_state=random_state)
            if PA == 3:
                ratio = abs(vec[0] / vec[1])  # ratio a/b
                abc[1] = np.sqrt(volume / (thickness1 * ratio))
                abc[0] = abc[1] * ratio
            elif PA == 2:
                ratio = abs(vec[0] / vec[2])  # ratio a/b
                abc[2] = np.sqrt(volume / (thickness1 * ratio))
                abc[0] = abc[2] * ratio
            elif PA == 1:
                ratio = abs(vec[1] / vec[2])  # ratio a/b
                abc[2] = np.sqrt(volume / (thickness1 * ratio))
                abc[1] = abc[2] * ratio

        # Tetragonal
        # elif num <= 41:
        elif ltype == "tetragonal":
            if PA == 3:
                abc[0] = abc[1] = np.sqrt(volume / thickness1)
            elif PA == 2:
                abc[0] = abc[1]
                abc[2] = volume / (abc[PA - 1] ** 2)
            elif PA == 1:
                abc[1] = abc[0]
                abc[2] = volume / (abc[PA - 1] ** 2)

        # Trigonal/Rhombohedral/Hexagonal
        # elif num <= 75:
        elif ltype in ["hexagonal", "trigonal"]:
            gamma = np.pi / 3 * 2
            x = np.sqrt(3.0) / 2.0
            if PA == 3:
                abc[0] = abc[1] = np.sqrt((volume / x) / abc[PA - 1])
            elif PA == 2:
                abc[0] = abc[1]
                abc[2] = (volume / x)(thickness1**2)
            elif PA == 1:
                abc[1] = abc[0]
                abc[2] = (volume / x) / (thickness1**2)

        para = np.array([abc[0], abc[1], abc[2], alpha, beta, gamma])

        a, b, c = abc[0], abc[1], abc[2]
        maxvec = (a * b * c) / (minvec**2)

        # Define limits on cell dimensions
        min_l = kwargs.get("min_l", minvec)
        mid_l = kwargs.get("mid_l", min_l)
        max_l = kwargs.get("max_l", mid_l)
        l_min = min(a, b, c)
        l_max = max(a, b, c)
        for x in (a, b, c):
            if x <= l_max and x >= l_min:
                l_mid = x
        if not (l_min >= min_l and l_mid >= mid_l and l_max >= max_l):
            continue

        if minvec < maxvec:
            smallvec = min(
                a * np.cos(max(beta, gamma)),
                b * np.cos(max(alpha, gamma)),
                c * np.cos(max(alpha, beta)),
            )
            if (
                a > minvec
                and b > minvec
                and c > minvec
                and a < maxvec
                and b < maxvec
                and c < maxvec
                and smallvec < minvec
                and alpha > minangle
                and beta > minangle
                and gamma > minangle
                and alpha < maxangle
                and beta < maxangle
                and gamma < maxangle
                and a / b < max_ratio
                and a / c < max_ratio
                and b / c < max_ratio
                and b / a < max_ratio
                and c / a < max_ratio
                and c / b < max_ratio
            ):
                return para

    # If maxattempts tries have been made without success
    msg = f"Could not get lattice after {maxattempts:d} cycles for volume {volume:.2f}"
    raise VolumeError(msg)


def generate_cellpara_0D(
    ltype,
    volume,
    area=None,
    minvec=1.2,
    max_ratio=10.0,
    maxattempts=100,
    random_state: None | int | Generator = None,
    **kwargs,
):
    """
    Generates a cell parameter (a, b, c, alpha, beta, gamma) according to the
    point group symmetry and number of atoms. If the generated lattice does
    not meet the minimum angle and vector requirements, we try to generate
    a new one, up to maxattempts times.

    Args:
        num: number of the Rod group
        volume: volume of the lattice
        area: cross-sectional area of the unit cell in Angstroms squared. If
            set to None, a value is chosen automatically
        minvec: minimum allowed lattice vector length (among a, b, and c)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice
        kwargs: a dictionary of optional values. Only used for ellipsoidal
            lattices, which pass the value to generate_lattice. They include:
            'unique_axis': the axis ('a', 'b', or 'c') which is unique.
            'min_l': the smallest allowed cell vector.
            'mid_l': the second smallest allowed cell vector.
            'max_l': the third smallest allowed cell vector.

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    if isinstance(random_state, int):
        # NOTE if random_state is an integer make a Generator to ensure randomness
        # downstream that would be lost if integer seed used repeated
        random_state = np.random.default_rng(random_state)

    if ltype == "spherical":
        # Use a cubic lattice with altered volume
        a = b = c = np.cbrt((3 * volume) / (4 * np.pi))
        alpha = beta = gamma = 0.5 * np.pi
        return np.array([a, b, c, alpha, beta, gamma])
    if ltype == "ellipsoidal":
        # Use a matrix with only on-diagonal elements, with a = b
        alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
        x = (4.0 / 3.0) * np.pi
        for _numattempts in range(maxattempts):
            vec = random_vector(random_state=random_state)
            c = vec[2] / (vec[0] * vec[1]) * np.cbrt(volume / x)
            a = b = np.sqrt((volume / x) / c)
            if (a / c < 10.0) and (c / a < 10.0):
                return np.array([a, b, c, alpha, beta, gamma])

    # If maxattempts tries have been made without success
    msg = f"Cannot get lattice after {maxattempts:d} cycles for volume {volume:.2f}"
    raise VolumeError(msg)


def matrix2para(matrix, radians=True):
    """
    Given a 3x3 matrix representing a unit cell, outputs a list of lattice
    parameters.

    Args:
        matrix: a 3x3 array or list, where the first, second, and third rows
            represent the a, b, and c vectors respectively
        radians: if True, outputs angles in radians. If False, outputs in
            degrees

    Returns:
        a 1x6 list of lattice parameters [a, b, c, alpha, beta, gamma]. a, b,
        and c are the length of the lattice vectos, and alpha, beta, and gamma
        are the angles between these vectors (in radians by default)
    """
    cell_para = np.zeros(6)
    # a
    cell_para[0] = np.linalg.norm(matrix[0])
    # b
    cell_para[1] = np.linalg.norm(matrix[1])
    # c
    cell_para[2] = np.linalg.norm(matrix[2])
    # alpha
    cell_para[3] = angle(matrix[1], matrix[2])
    # beta
    cell_para[4] = angle(matrix[0], matrix[2])
    # gamma
    cell_para[5] = angle(matrix[0], matrix[1])

    if not radians:
        # convert radians to degrees
        deg = 180.0 / np.pi
        cell_para[3] *= deg
        cell_para[4] *= deg
        cell_para[5] *= deg
    return cell_para


# def para2matrix(cell_para, radians=True, format="lower"):
def para2matrix(cell_para, radians=True, format="upper"):
    """
    Given a set of lattic parameters, generates a matrix representing the
    lattice vectors

    Args:
        cell_para: a 1x6 list of lattice parameters [a, b, c, alpha, beta,
            gamma]. a, b, and c are the length of the lattice vectos, and
            alpha, beta, and gamma are the angles between these vectors. Can
            be generated by matrix2para
        radians: if True, lattice parameters should be in radians. If False,
            lattice angles should be in degrees
        format: a string ('lower', 'symmetric', or 'upper') for the type of
            matrix to be output

    Returns:
        a 3x3 matrix representing the unit cell. By default (format='lower'),
        the a vector is aligined along the x-axis, and the b vector is in the
        y-z plane
    """
    a = cell_para[0]
    b = cell_para[1]
    c = cell_para[2]
    alpha = cell_para[3]
    beta = cell_para[4]
    gamma = cell_para[5]
    if radians is not True:
        alpha *= rad
        beta *= rad
        gamma *= rad
    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)
    sin_alpha = np.sin(alpha)
    matrix = np.zeros([3, 3])

    if format == "lower":
        # Generate a lower-diagonal matrix
        c1 = c * cos_beta
        c2 = (c * (cos_alpha - (cos_beta * cos_gamma))) / sin_gamma
        matrix[0][0] = a
        matrix[1][0] = b * cos_gamma
        matrix[1][1] = b * sin_gamma
        matrix[2][0] = c1
        matrix[2][1] = c2
        matrix[2][2] = np.sqrt(c**2 - c1**2 - c2**2)
    elif format == "symmetric":
        # TODO: allow generation of symmetric matrices
        pass
    elif format == "upper":
        # Generate an upper-diagonal matrix
        a3 = a * cos_beta
        a2 = (a * (cos_gamma - (cos_beta * cos_alpha))) / sin_alpha
        matrix[2][2] = c
        matrix[1][2] = b * cos_alpha
        matrix[1][1] = b * sin_alpha
        matrix[0][2] = a3
        matrix[0][1] = a2
        tmp = a**2 - a3**2 - a2**2
        if tmp > 0:
            matrix[0][0] = np.sqrt(a**2 - a3**2 - a2**2)
        # elif abs(tmp) < 1e-5: #tmp is very close to 0
        #    matrix[0][0] = 0
        #    print(matrix)
        else:
            # print(tmp)
            return None
        # pass
    return matrix


def gaussian_random_variable(min, max, sigma=3.0, random_state: None | int | Generator = None):
    """
    Choose a random number from a Gaussian probability distribution centered
    between min and max. sigma is the number of standard deviations that min
    and max are away from the center. Thus, sigma is also the largest possible
    number of standard deviations corresponding to the returned value. sigma=2
    corresponds to a 95.45% probability of choosing a number between min and
    max.

    Args:
        min: the minimum acceptable value
        max: the maximum acceptable value
        sigma: the number of standard deviations between the center and min/max

    Returns:
        a value chosen randomly between min and max
    """
    if isinstance(random_state, Generator):
        random_state = random_state.spawn(1)[0]
    else:
        random_state = np.random.default_rng(random_state)

    center = (max + min) * 0.5
    delta = np.fabs(max - min) * 0.5
    ratio = delta / sigma
    while True:
        x = random_state.normal(scale=ratio, loc=center)
        if x > min and x < max:
            return x


def random_vector(minvec=None, maxvec=None, width=0.35, unit=False, random_state: None | int | Generator = None):
    """
    Generate a random vector for lattice constant generation. The ratios between
    x, y, and z of the returned vector correspond to the ratios between a, b,
    and c. Results in a Gaussian distribution of the natural log of the ratios.

    Args:
        minvec: the bottom-left-back minimum point which can be chosen
        maxvec: the top-right-front maximum point which can be chosen
        width: the width of the normal distribution to use when choosing values.
            Passed to np.random.normal
        unit: whether or not to normalize the vector to determinant 1

    Returns:
        a 1x3 numpy array of floats
    """
    if isinstance(random_state, Generator):
        random_state = random_state.spawn(1)[0]
    else:
        random_state = np.random.default_rng(random_state)

    # TODO these were not used in the original code, moved out of the defaults
    # and declared if None to avoid having lists as default arguments.
    if maxvec is None:
        maxvec = [1.0, 1.0, 1.0]
    if minvec is None:
        minvec = [0.0, 0.0, 0.0]

    vec = np.exp(random_state.normal(scale=width, size=3))
    if unit:
        return vec / np.linalg.norm(vec)
    else:
        return vec


def random_shear_matrix(width=1.0, unitary=False, random_state: None | int | Generator = None):
    """
    Generate a random symmetric shear matrix with Gaussian elements. If unitary
    is True, normalize to determinant 1

    Args:
        width: the width of the normal distribution to use when choosing values.
            Passed to np.random.normal
        unitary: whether or not to normalize the matrix to determinant 1

    Returns:
        a 3x3 numpy array of floats
    """
    if isinstance(random_state, Generator):
        random_state = random_state.spawn(1)[0]
    else:
        random_state = np.random.default_rng(random_state)

    mat = np.zeros([3, 3])
    determinant = 0
    while determinant == 0:
        a, b, c = random_state.normal(scale=width, size=3)
        mat = np.array([[1, a, b], [a, 1, c], [b, c, 1]])
        determinant = np.linalg.det(mat)
    if unitary:
        return mat / np.cbrt(np.linalg.det(mat))
    else:
        return mat
