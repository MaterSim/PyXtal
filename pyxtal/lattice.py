# Standard Libraries
import numpy as np
import random

# PyXtal imports
from pyxtal.msg import printx, VolumeError
from pyxtal.operations import angle, create_matrix
from pyxtal.constants import deg, rad

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
        PBC: A periodic boundary condition list, where 1 means periodic,
            0 means not periodic.
            Ex: [1, 1, 1] -> full 3d periodicity, [0, 0, 1] -> periodicity
            along the z axis
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

    def __init__(self, ltype, volume, PBC=[1, 1, 1], **kwargs):
        # Set required parameters
        if ltype in [
            "triclinic", "Triclinic",
            "monoclinic", "Monoclinic",
            "orthorhombic", "Orthorhombic",
            "tetragonal", "Tetragonal",
            "trigonal", "Trigonal",
            "hexagonal", "Hexagonal",
            "cubic", "Cubic",
            "spherical", "Spherical",
            "ellipsoidal", "Ellipsoidal",
        ]:
            self.ltype = ltype
        elif ltype == None:
            self.ltype = "triclinic"
        else:
            printx("Error: Invalid lattice type.", priority=1)
            return
        self.volume = float(volume)
        self.PBC = PBC
        self.dim = sum(PBC)
        self.kwargs = {}
        self.random = True
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
            ]:
                setattr(self, key, value)
                self.kwargs[key] = value
                if key == "allow_volume_reset":
                    if value == False:
                        self.allow_volume_reset = False
        try:
            self.unique_axis
        except:
            self.unique_axis = "c"
        # Set stress normalization info
        if self.ltype == "triclinic":
            self.stress_normalization_matrix = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        elif self.ltype == "monoclinic":
            if self.PBC == [1, 1, 1]:
                self.stress_normalization_matrix = np.array(
                    [[1, 0, 0], [0, 1, 0], [1, 0, 1]]
                )
            else:
                if self.unique_axis == "a":
                    self.stress_normalization_matrix = np.array(
                        [[1, 0, 0], [0, 1, 0], [0, 1, 1]]
                    )
                elif self.unique_axis == "b":
                    self.stress_normalization_matrix = np.array(
                        [[1, 0, 0], [0, 1, 0], [1, 0, 1]]
                    )
                elif self.unique_axis == "c":
                    self.stress_normalization_matrix = np.array(
                        [[1, 0, 0], [1, 1, 0], [0, 0, 1]]
                    )
        elif self.ltype in ["orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic",
                            "Orthorhombic", "Tetragonal", "Trigonal", "Hexagonal", "Cubic",]:
            self.stress_normalization_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        elif self.ltype in ["spherical", "ellipsoidal"]:
            self.stress_normalization_matrix = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        # Set info for on-diagonal stress symmetrization
        if self.ltype in ["tetragonal", "trigonal", "hexagonal", "rhombohedral",
                          "Tetragonal", "Trigonal", "Hexagonal", "Rhombohedral"]:
            self.stress_indices = [(0, 0), (1, 1)]
        elif self.ltype in ["cubic", "Cubic"]:
            self.stress_indices = [(0, 0), (1, 1), (2, 2)]
        else:
            self.stress_indices = []
        # Set values for the matrix
        self.reset_matrix()
        self._get_dof()

    def _get_dof(self):
        """
        get the number of degree of freedom
        """
        if self.ltype in ["triclinic", "Triclinic"]:
            self.dof = 6
        elif self.ltype in ["monoclinic", "Monoclinic"]:
            self.dof = 4
        elif self.ltype in ['orthorhombic', 'Orthorhombic']:
            self.dof = 3
        elif self.ltype in ['tetragonal', 'Tetragonal',
            'hexagonal', 'trigonal', 'Hexagonal', 'Trigonal']:
            self.dof = 2
        else:
            self.dof = 1

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

    def optimize(self):
        """
        Optimize the lattice's inclination angles
        """
        trans = np.zeros([1,3,3])
        trans[0] = np.eye(3)
        tmp = None
        opt = False
        if self.ltype in ["monoclinic", "Monoclinic"]:
            tmp = np.array([[[1,0,0],[0,1,0],[1,0,1]],
                           [[1,0,0],[0,1,0],[-1,0,1]],
                           [[1,0,1],[0,1,0],[0,0,1]],
                           [[1,0,-1],[0,1,0],[0,0,1]]])

        elif self.ltype in ["triclinic", "Triclinic"]:
            angles = np.array([self.alpha, self.beta, self.gamma])
            diff0s = np.abs(angles - np.pi/2)
            axis = np.argmax(diff0s)
            if axis == 1:
                tmp = np.array([[[1,0,0],[0,1,0],[1,0,1]],
                               [[1,0,0],[0,1,0],[-1,0,1]],
                               [[1,0,1],[0,1,0],[0,0,1]],
                               [[1,0,-1],[0,1,0],[0,0,1]]])
            elif axis == 0:
                tmp = np.array([[[1,0,0],[0,1,0],[0,1,1]],
                                [[1,0,0],[0,1,1],[0,0,1]],
                                [[1,0,0],[0,1,0],[0,-1,1]],
                                [[1,0,0],[0,1,-1],[0,0,1]]])
            else:
                tmp = np.array([[[1,1,0],[0,1,0],[0,0,1]],
                                [[1,-1,0],[0,1,0],[0,0,1]],
                                [[1,0,0],[1,1,0],[0,0,1]],
                                [[1,0,0],[-1,1,0],[0,0,1]]])
        if tmp is not None:
            trans = np.append(trans, tmp, axis=0)
            diffs = []
            for tran in trans:
                cell_new = np.dot(tran, self.matrix)
                try:
                    lat_new = Lattice.from_matrix(cell_new)
                    _, _, _, alpha, beta, gamma = lat_new.get_para()
                    diffs.append(np.max(abs(np.array([alpha, beta, gamma])-np.pi/2)))
                except:
                    diffs.append(100)
            id = np.array(diffs).argmin()
            tran = trans[id]
            cell = np.dot(tran, self.matrix)
            if id > 0:
                opt = True
            #print("cell:"); print(cell) #; import sys; sys.exit()
            return Lattice.from_matrix(cell, ltype=self.ltype, reset=False), tran, opt
        else:
            return self, np.eye(3), opt

    def transform(self, trans_mat=np.eye(3), reset=False):
        """
        Optimize the lattice's inclination angles
        """
        if type(trans_mat) == list:
            trans_mat = np.array(trans_mat)
        cell = np.dot(trans_mat, self.matrix)
        return Lattice.from_matrix(cell, ltype=self.ltype, reset=reset)

    def encode(self):
        a, b, c, alpha, beta, gamma = self.get_para(degree=True)
        if self.ltype in ['cubic', 'Cubic']:
            return [a]
        elif self.ltype in ['hexagonal', 'trigonal', 'Hexagonal', 'Trigonal', 'tetragonal', 'Tetragonal']:
            return [a, c]
        elif self.ltype in ['orthorhombic', 'Orthorhombic']:
            return [a, b, c]
        elif self.ltype in ['monoclinic', 'Monoclinic']:
            return [a, b, c, beta]
        else:
            return [a, b, c, alpha, beta, gamma]

    def mutate(self, degree=0.20, frozen=False):
        """
        mutate the lattice object
        """
        rand = 1 + degree*(np.random.sample(6)-0.5)
        a0, b0, c0, alpha0, beta0, gamma0 = self.get_para()
        a = a0*rand[0]
        b = b0*rand[1]
        c = c0*rand[2]
        alpha = np.degrees(alpha0*rand[3])
        beta = np.degrees(beta0*rand[4])
        gamma = np.degrees(gamma0*rand[5])
        ltype = self.ltype

        if self.ltype in ['cubic', 'Cubic']:
            if frozen:
                lat = Lattice.from_para(a0, a0, a0, 90, 90, 90, ltype=ltype)
            else:
                lat = Lattice.from_para(a, a, a, 90, 90, 90, ltype=ltype)
        elif ltype in ['hexagonal', 'trigonal', 'Hexagonal', 'Trigonal']:
            if frozen:
                lat = Lattice.from_para(a0, a0, c, 90, 90, 120, ltype=ltype)
            else:
                lat = Lattice.from_para(a, a, c, 90, 90, 120, ltype=ltype)
        #elif ltype in ['rhombohedral', 'Rhombohedral']:
        #    if forzen:
        #        lat = Lattice.from_para(a0, a0, a0, alpha0, alpha0, alpha0, ltype=ltype)
        #    else:
        #        lat = Lattice.from_para(a, a, a, alpha, alpha, alpha, ltype=ltype)
        elif ltype in ['tetragonal', 'Tetragonal']:
            if frozen:
                lat = Lattice.from_para(a0, a0, c, 90, 90, 90, ltype=ltype)
            else:
                lat = Lattice.from_para(a, a, c, 90, 90, 90, ltype=ltype)
        elif ltype in ['orthorhombic', 'Orthorhombic']:
            lat = Lattice.from_para(a, b, c, 90, 90, 90, ltype=ltype)
        elif ltype in ['monoclinic', 'Monoclinic']:
            lat = Lattice.from_para(a, b, c, 90, beta, 90, ltype=ltype)
        elif ltype in ['triclinic', 'Triclinic']:
            lat = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
        else:
            raise ValueError("ltype {:s} is not supported".format(ltype))
        return lat

    def generate_para(self):
        if self.dim == 3:
            return generate_lattice(self.ltype, self.volume, **self.kwargs)
        elif self.dim == 2:
            return generate_lattice_2D(self.ltype, self.volume, **self.kwargs)
        elif self.dim == 1:
            return generate_lattice_1D(self.ltype, self.volume, **self.kwargs)
        elif self.dim == 0:
            return generate_lattice_0D(self.ltype, self.volume, **self.kwargs)

    def generate_matrix(self):
        """
        Generates a 3x3 matrix for the lattice based on the lattice type and volume
        """
        # Try multiple times in case of failure
        for i in range(10):
            para = self.generate_para()
            if para is not None:
                return para2matrix(para)

    def get_matrix(self, shape='upper'):
        """
        Returns a 3x3 numpy array representing the lattice vectors.
        """
        try:
            return self.matrix
        # TODO remove bare except
        except:
            printx("Error: Lattice matrix undefined.", priority=1)
            return

    def get_para(self, degree=False):
        """
        Returns a tuple of lattice parameters.
        """
        if degree:
            return (self.a, self.b, self.c, deg*self.alpha, deg*self.beta, deg*self.gamma)
        else:
            return (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    def set_matrix(self, matrix=None):
        if matrix is not None:
            m = np.array(matrix)
            if np.shape(m) == (3, 3):
                self.matrix = m
                self.inv_matrix = np.linalg.inv(m)
            else:
                printx("Error: matrix must be a 3x3 numpy array or list", priority=1)
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

    def reset_matrix(self, shape='upper'):
        if self.random:
            for i in range(3):
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
                    break
        else:
            # a small utility to convert the cell shape
            para = matrix2para(self.matrix)
            self.matrix = para2matrix(para, format=shape)

    def set_volume(self, volume):
        if self.allow_volume_reset is True:
            self.volume = volume

    def swap_axis(self, random=False, ids=None):
        """
        For the lattice
        """
        # only applied to triclinic/monoclinic/orthorhombic
        if self.ltype in ["triclinic", "Triclinic", "orthorhombic", "Orthorhombic"]:
            allowed_ids = [[0,1,2],[1,0,2],[0,2,1],[2,1,0],[1,2,0],[2,0,1]]
        elif self.ltype == "monoclinic":
            if abs(self.beta-90*rad) > 1e-3:
                allowed_ids = [[0,1,2],[2,1,0]]
            else:
                allowed_ids = [[0,1,2],[1,0,2],[0,2,1],
                               [2,1,0],[1,2,0],[2,0,1]]
        else:
            allowed_ids = [[0,1,2]]

        if random:
            from random import choice
            ids = choice(allowed_ids)
        else:
            if ids not in allowed_ids:
                print(ids)
                raise ValueError("the above swap is not allowed in "+self.ltype)

        (a,b,c,alpha,beta,gamma) = self.get_para()
        alpha, beta, gamma = alpha*deg, beta*deg, gamma*deg
        if ids is None:
            return self
        elif ids == [1,0,2]: #a->b
            return self.from_para(b, a, c, beta, alpha, gamma, self.ltype)
        elif ids == [2,1,0]: #a->c
            return self.from_para(c, b, a, gamma, beta, alpha, self.ltype)
        elif ids == [0,2,1]: #b-c
            return self.from_para(a, c, b, alpha, gamma, beta, self.ltype)
        elif ids == [2,0,1]:
            return self.from_para(c, a, b, gamma, alpha, beta, self.ltype)
        elif ids == [1,2,0]:
            return self.from_para(b, c, a, beta, gamma, alpha, self.ltype)
        else:
            return self

    def swap_angle(self, random=True, ids=None):
        # only applied to triclinic/monoclinic #/hexagonal
        """
        If the angle is not 90. There will be two equivalent versions
        e.g., 80 and 100.
        """
        if self.ltype in ["monoclinic", "Monoclinic"]:
            allowed_ids = ["beta", "No"]
        elif self.ltype in ["triclinic", "Triclinic"]:
            allowed_ids = ["alpha", "beta", "gamma", "No"]
        else:
            allowed_ids = ["No"]

        if random:
            from random import choice
            ids = choice(allowed_ids)
        else:
            if ids not in allowed_ids:
                print(ids)
                raise ValueError("the above swap is not allowed in "+self.ltype)

        (a,b,c,alpha,beta,gamma) = self.get_para()
        alpha, beta, gamma = alpha*deg, beta*deg, gamma*deg
        if ids is None:
            return self
        elif ids == "alpha":
            return self.from_para(a, b, c, 180-alpha, beta, gamma, self.ltype)
        elif ids == "beta":
            return self.from_para(a, b, c, alpha, 180-beta, gamma, self.ltype)
        elif ids == "gamma":
            return self.from_para(a, b, c, alpha, beta, 180-gamma, self.ltype)
        else:
            return self

    def add_vacuum(self, coor, frac=True, vacuum=15, PBC=[0, 0, 0]):
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
        matrix = self.matrix
        if frac:
            absolute_coords = np.dot(coor, matrix)
        else:
            absolute_coords = coor

        for i, a in enumerate(PBC):
            if not a:
                ratio = 1 + vacuum/np.linalg.norm(matrix[i])
                matrix[i] *= ratio
                absolute_coords[:, i] += vacuum/2
        if frac:
            coor = np.dot(absolute_coords, np.linalg.inv(matrix))
        else:
            coor = absolute_coords
        return matrix, coor


    def generate_point(self):

        # point = np.random.RandomState().rand(3)
        # QZ: it was here because of multiprocess issue
        # https://github.com/numpy/numpy/issues/9650
        # now just fix it

        point = np.random.rand(3)
        if self.ltype in ["spherical", "ellipsoidal"]:
            # Choose a point within an octant of the unit sphere
            while point.dot(point) > 1:  # squared
                point = np.random.random(3)
            # Randomly flip some coordinates
            for index in range(len(point)):
                # Scale the point by the max radius
                if random.uniform(0, 1) < 0.5:
                    point[index] *= -1
        else:
            for i, a in enumerate(self.PBC):
                if not a:
                    if self.ltype in ["hexagonal", "trigonal", "rhombohedral"]:
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
        PBC=[1, 1, 1],
        factor=1.0,
        **kwargs
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
                random ones will be generated. Values will be passed to generate_lattice.
                Options include:
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
                'unique_axis': the axis ('a', 'b', or 'c') which is not symmetrically
                    equivalent to the other two
                'min_l': the smallest allowed cell vector. The smallest vector must
                    be larger than this.
                'mid_l': the second smallest allowed cell vector. The second
                    smallest vector must be larger than this.
                'max_l': the third smallest allowed cell vector. The largest cell
                    vector must be larger than this.

        Returns:
            a Lattice object with the specified parameters
        """
        try:
            cell_matrix = factor*para2matrix((a, b, c, alpha, beta, gamma), radians=radians)
        except:
            printx("Error: invalid cell parameters for lattice.", priority=1)
            return
        volume = np.linalg.det(cell_matrix)
        # Initialize a Lattice instance
        l = Lattice(ltype, volume, PBC=PBC, **kwargs)
        l.a, l.b, l.c = factor*a, factor*b, factor*c
        l.alpha, l.beta, l.gamma = alpha * rad, beta * rad, gamma * rad
        l.matrix = cell_matrix
        l.inv_matrix = np.linalg.inv(cell_matrix)
        l.ltype = ltype
        l.volume = volume
        l.random = False
        l.allow_volume_reset = False
        return l

    @classmethod
    def from_matrix(self, matrix, reset=True, shape='upper', ltype="triclinic", PBC=[1, 1, 1], **kwargs):
        """
        Creates a Lattice object from a 3x3 cell matrix. Additional keyword arguments
        are available. Unless specified by the keyword random=True, does not create a
        new matrix upon calling reset_matrix. This allows for generation of random
        crystals with a specific choice of unit cell.

        Args:
            matrix: a 3x3 real matrix (numpy array or nested list) describing the cell vectors
            ltype: the lattice type ("cubic, tetragonal, etc."). Also available are "spherical",
                which confines generated points to lie within a sphere, and "ellipsoidal", which
                confines generated points to lie within an ellipsoid (oriented about the z axis)
            PBC: A periodic boundary condition list, where 1 means periodic, 0 means not periodic.
                Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity along the z axis
            kwargs: various values which may be defined. If none are defined, random ones
                will be generated. Values will be passed to generate_lattice. Options include:
                area: The cross-sectional area (in Angstroms squared). Only used to generate 1D
                    crystals
                thickness: The unit cell's non-periodic thickness (in Angstroms). Only used to
                    generate 2D crystals
                unique_axis: The unique axis for certain symmetry (and especially layer) groups.
                    Because the symmetry operations are not also transformed, you should use the
                    default values for random crystal generation
                random: If False, keeps the stored values for the lattice geometry even applying
                    reset_matrix. To alter the matrix, use set_matrix() or set_para
                'unique_axis': the axis ('a', 'b', or 'c') which is not symmetrically
                    equivalent to the other two
                'min_l': the smallest allowed cell vector. The smallest vector must be larger
                    than this.
                'mid_l': the second smallest allowed cell vector. The second smallest vector
                    must be larger than this.
                'max_l': the third smallest allowed cell vector. The largest cell vector must
                    be larger than this.

        Returns:
            a Lattice object with the specified parameters
        """
        m = np.array(matrix)
        if np.shape(m) != (3, 3):
            printx("Error: Lattice matrix must be 3x3", priority=1)
            return
        [a, b, c, alpha, beta, gamma] = matrix2para(m)

        # symmetrize the lattice
        if reset:
            if ltype in ['cubic', 'Cubic']:
                a = b = c = (a+b+c)/3
                alpha = beta = gamma = np.pi/2
            elif ltype in ['hexagonal', 'trigonal', 'Hexagonal', 'Trigonal']:
                a = b = (a+b)/2
                alpha = beta = np.pi/2
                gamma = np.pi*2/3
            elif ltype in ['tetragonal', 'Tetragonal']:
                a = b = (a+b)/2
                alpha = beta = gamma = np.pi/2
            elif ltype in ['orthorhombic', 'Orthorhombic']:
                alpha = beta = gamma = np.pi/2
            elif ltype in ['monoclinic', 'Monoclinic']:
                alpha = gamma = np.pi/2
            
            # reset matrix according to the symmetry
            m = para2matrix([a, b, c, alpha, beta, gamma], format=shape)
        
        # Initialize a Lattice instance
        volume = np.linalg.det(m)
        l = Lattice(ltype, volume, PBC=PBC, **kwargs)
        l.a, l.b, l.c = a, b, c
        l.alpha, l.beta, l.gamma = alpha, beta, gamma
        l.matrix = m
        l.inv_matrix = np.linalg.inv(m)
        l.ltype = ltype
        l.volume = volume
        l.random = False
        l.allow_volume_reset = False
        return l


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
        matrix = np.dot(trans.T, self.get_matrix())
        l1 = Lattice.from_matrix(matrix)
        l2 = Lattice.from_matrix(matrix, ltype=l_type)
        (a1, b1, c1, alpha1, beta1, gamma1) = l1.get_para(degree=True)
        (a2, b2, c2, alpha2, beta2, gamma2) = l2.get_para(degree=True)
        abc_diff = np.abs(np.array([a2-a1, b2-b1, c2-c1])).max()
        ang_diff = np.abs(np.array([alpha2-alpha1, beta2-beta1, gamma2-gamma1])).max()
        if abc_diff > tol or ang_diff > a_tol:
            return False
        else:
            return True

    def __str__(self):
        s = "{:s}: {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}".format(
            str(self.ltype),
            self.a,
            self.b,
            self.c,
            self.alpha * deg,
            self.beta * deg,
            self.gamma * deg,
        )
        return s

    def __repr__(self):
        return str(self)

def generate_lattice(
    ltype,
    volume,
    minvec=1.2,
    minangle=np.pi / 6,
    max_ratio=10.0,
    maxattempts=100,
    **kwargs
):
    """
    Generates a lattice (3x3 matrix) according to the space group symmetry and
    number of atoms. If the spacegroup has centering, we will transform to
    conventional cell setting. If the generated lattice does not meet the
    minimum angle and vector requirements, we try to generate a new one, up to
    maxattempts times.

    Args:
        volume: volume of the conventional unit cell
        minvec: minimum allowed lattice vector length (among a, b, and c)
        minangle: minimum allowed lattice angle (among alpha, beta, and gamma)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice
        kwargs: a dictionary of optional values. These include:
            'unique_axis': the axis ('a', 'b', or 'c') which is not symmetrically
                equivalent to the other two
            'min_l': the smallest allowed cell vector. The smallest vector must be larger
                than this.
            'mid_l': the second smallest allowed cell vector. The second smallest vector
                must be larger than this.
            'max_l': the third smallest allowed cell vector. The largest cell vector must
                be larger than this.

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    maxangle = np.pi - minangle
    for n in range(maxattempts):
        # Triclinic
        # if sg <= 2:
        if ltype == "triclinic":
            # Derive lattice constants from a random matrix
            mat = random_shear_matrix(width=0.2)
            a, b, c, alpha, beta, gamma = matrix2para(mat)
            x = np.sqrt(
                1
                - np.cos(alpha) ** 2
                - np.cos(beta) ** 2
                - np.cos(gamma) ** 2
                + 2 * (np.cos(alpha) * np.cos(beta) * np.cos(gamma))
            )
            vec = random_vector()
            abc = volume / x
            xyz = vec[0] * vec[1] * vec[2]
            a = vec[0] * np.cbrt(abc) / np.cbrt(xyz)
            b = vec[1] * np.cbrt(abc) / np.cbrt(xyz)
            c = vec[2] * np.cbrt(abc) / np.cbrt(xyz)
        # Monoclinic
        elif ltype in ["monoclinic", "Monoclinic"]:
            alpha, gamma = np.pi / 2, np.pi / 2
            beta = gaussian(minangle, maxangle)
            x = np.sin(beta)
            vec = random_vector()
            xyz = vec[0] * vec[1] * vec[2]
            abc = volume / x
            a = vec[0] * np.cbrt(abc) / np.cbrt(xyz)
            b = vec[1] * np.cbrt(abc) / np.cbrt(xyz)
            c = vec[2] * np.cbrt(abc) / np.cbrt(xyz)
        # Orthorhombic
        # elif sg <= 74:
        elif ltype in ["orthorhombic", "Orthorhombic"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
            x = 1
            vec = random_vector()
            xyz = vec[0] * vec[1] * vec[2]
            abc = volume / x
            a = vec[0] * np.cbrt(abc) / np.cbrt(xyz)
            b = vec[1] * np.cbrt(abc) / np.cbrt(xyz)
            c = vec[2] * np.cbrt(abc) / np.cbrt(xyz)
        # Tetragonal
        # elif sg <= 142:
        elif ltype in ["tetragonal", "Tetragonal"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
            x = 1
            vec = random_vector()
            c = vec[2] / (vec[0] * vec[1]) * np.cbrt(volume / x)
            a = b = np.sqrt((volume / x) / c)
        # Trigonal/Rhombohedral/Hexagonal
        # elif sg <= 194:
        elif ltype in ["hexagonal", "trigonal", "rhombohedral"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 3 * 2
            x = np.sqrt(3.0) / 2.0
            vec = random_vector()
            c = vec[2] / (vec[0] * vec[1]) * np.cbrt(volume / x)
            a = b = np.sqrt((volume / x) / c)
        # Cubic
        # else:
        elif ltype in ["cubic", "Cubic"]:
            alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
            s = (volume) ** (1.0 / 3.0)
            a, b, c = s, s, s
        # Check that lattice meets requirements
        maxvec = (a * b * c) / (minvec ** 2)

        # Define limits on cell dimensions
        if "min_l" not in kwargs:
            min_l = minvec
        else:
            min_l = kwargs["min_l"]
        if "mid_l" not in kwargs:
            mid_l = min_l
        else:
            mid_l = kwargs["mid_l"]
        if "max_l" not in kwargs:
            max_l = mid_l
        else:
            max_l = kwargs["max_l"]
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
    msg = "Could not get lattice after {:d} cycles for volume {:.2f}".format(maxattempts, volume)
    raise VolumeError(msg)
    #return


def generate_lattice_2D(
    ltype,
    volume,
    thickness=None,
    minvec=1.2,
    minangle=np.pi / 6,
    max_ratio=10.0,
    maxattempts=100,
    **kwargs
):
    """
    Generates a lattice (3x3 matrix) according to the spacegroup symmetry and
    number of atoms. If the layer group has centering, we will use the
    conventional cell setting. If the generated lattice does not meet the
    minimum angle and vector requirements, we try to generate a new one, up to
    maxattempts times.
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
            'unique_axis': the axis ('a', 'b', or 'c') which is not symmetrically
                equivalent to the other two
            'min_l': the smallest allowed cell vector. The smallest vector must be larger
                than this.
            'mid_l': the second smallest allowed cell vector. The second smallest vector
                must be larger than this.
            'max_l': the third smallest allowed cell vector. The largest cell vector must
                be larger than this.

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    if "unique_axis" not in kwargs:
        unique_axis = "c"
    else:
        unique_axis = kwargs["unique_axis"]
    # Store the non-periodic axis
    NPA = 3
    # Set the unique axis for monoclinic cells
    # if num in range(3, 8): unique_axis = "c"
    # elif num in range(8, 19): unique_axis = "a"
    maxangle = np.pi - minangle
    for n in range(maxattempts):
        abc = np.ones([3])
        if thickness is None:
            v = random_vector()
            thickness1 = np.cbrt(volume) * (v[0] / (v[0] * v[1] * v[2]))
        else:
            thickness1 = max([3.0, thickness])
        abc[NPA - 1] = thickness1
        alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
        # Triclinic
        # if num <= 2:
        if ltype == "triclinic":
            mat = random_shear_matrix(width=0.2)
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
            a, b, c = random_vector()
            if unique_axis == "a":
                alpha = gaussian(minangle, maxangle)
                x = np.sin(alpha)
            elif unique_axis == "b":
                beta = gaussian(minangle, maxangle)
                x = np.sin(beta)
            elif unique_axis == "c":
                gamma = gaussian(minangle, maxangle)
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
            vec = random_vector()
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
                abc[2] = (volume / x)(thickness1 ** 2)
            elif NPA == 1:
                abc[1] = abc[0]
                abc[2] = (volume / x) / (thickness1 ** 2)

        para = np.array([abc[0], abc[1], abc[2], alpha, beta, gamma])

        a, b, c = abc[0], abc[1], abc[2]
        maxvec = (a * b * c) / (minvec ** 2)

        # Define limits on cell dimensions
        if "min_l" not in kwargs:
            min_l = minvec
        else:
            min_l = kwargs["min_l"]
        if "mid_l" not in kwargs:
            mid_l = min_l
        else:
            mid_l = kwargs["mid_l"]
        if "max_l" not in kwargs:
            max_l = mid_l
        else:
            max_l = kwargs["max_l"]
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
    msg = "Could not get lattice after {:d} cycles for volume {:.2f}".format(maxattempts, volume)
    raise VolumeError(msg)

def generate_lattice_1D(
    ltype,
    volume,
    area=None,
    minvec=1.2,
    minangle=np.pi / 6,
    max_ratio=10.0,
    maxattempts=100,
    **kwargs
):
    """
    Generates a lattice (3x3 matrix) according to the spacegroup symmetry and
    number of atoms. If the spacegroup has centering, we will transform to
    conventional cell setting. If the generated lattice does not meet the
    minimum angle and vector requirements, we try to generate a new one, up to
    maxattempts times.
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
            'unique_axis': the axis ('a', 'b', or 'c') which is not symmetrically
                equivalent to the other two
            'min_l': the smallest allowed cell vector. The smallest vector must be larger
                than this.
            'mid_l': the second smallest allowed cell vector. The second smallest vector
                must be larger than this.
            'max_l': the third smallest allowed cell vector. The largest cell vector must
                be larger than this.

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
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
    for n in range(maxattempts):
        abc = np.ones([3])
        if area is None:
            v = random_vector()
            thickness1 = np.cbrt(volume) * (v[0] / (v[0] * v[1] * v[2]))
        else:
            thickness1 = volume / area
        abc[PA - 1] = thickness1
        alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
        # Triclinic
        # if num <= 2:
        if ltype == "triclinic":
            mat = random_shear_matrix(width=0.2)
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
            a, b, c = random_vector()
            if unique_axis == "a":
                alhpa = gaussian(minangle, maxangle)
                x = np.sin(alpha)
            elif unique_axis == "b":
                beta = gaussian(minangle, maxangle)
                x = np.sin(beta)
            elif unique_axis == "c":
                gamma = gaussian(minangle, maxangle)
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
            vec = random_vector()
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
                abc[2] = (volume / x)(thickness1 ** 2)
            elif PA == 1:
                abc[1] = abc[0]
                abc[2] = (volume / x) / (thickness1 ** 2)

        para = np.array([abc[0], abc[1], abc[2], alpha, beta, gamma])

        a, b, c = abc[0], abc[1], abc[2]
        maxvec = (a * b * c) / (minvec ** 2)

        # Define limits on cell dimensions
        if "min_l" not in kwargs:
            min_l = minvec
        else:
            min_l = kwargs["min_l"]
        if "mid_l" not in kwargs:
            mid_l = min_l
        else:
            mid_l = kwargs["mid_l"]
        if "max_l" not in kwargs:
            max_l = mid_l
        else:
            max_l = kwargs["max_l"]
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
    msg = "Could not get lattice after {:d} cycles for volume {:.2f}".format(maxattempts, volume)
    raise VolumeError(msg)


def generate_lattice_0D(
    ltype, volume, area=None, minvec=1.2, max_ratio=10.0, maxattempts=100, **kwargs
):
    """
    Generates a lattice (3x3 matrix) according to the spacegroup symmetry and
    number of atoms. If the spacegroup has centering, we will transform to
    conventional cell setting. If the generated lattice does not meet the
    minimum angle and vector requirements, we try to generate a new one, up to
    maxattempts times.
    Note: The monoclinic Rod groups have different unique axes. Groups 3-7
        have unique axis a, while 8-12 have unique axis c. We use periodic
        axis c for all Rod groups.

    Args:
        num: number of the Rod group
        volume: volume of the lattice
        area: cross-sectional area of the unit cell in Angstroms squared. If
            set to None, a value is chosen automatically
        minvec: minimum allowed lattice vector length (among a, b, and c)
        max_ratio: largest allowed ratio of two lattice vector lengths
        maxattempts: the maximum number of attempts for generating a lattice
        kwargs: a dictionary of optional values. Only used for ellipsoidal
            lattices, which pass the value to generate_lattice. Possible values include:
            'unique_axis': the axis ('a', 'b', or 'c') which is not symmetrically
                equivalent to the other two
            'min_l': the smallest allowed cell vector. The smallest vector must be larger
                than this.
            'mid_l': the second smallest allowed cell vector. The second smallest vector
                must be larger than this.
            'max_l': the third smallest allowed cell vector. The largest cell vector must
                be larger than this.

    Returns:
        a 3x3 matrix representing the lattice vectors of the unit cell. If
        generation fails, outputs a warning message and returns empty
    """
    if ltype == "spherical":
        # Use a cubic lattice with altered volume
        a = b = c = np.cbrt((3 * volume) / (4 * np.pi))
        alpha = beta = gamma = 0.5 * np.pi
        return np.array([a, b, c, alpha, beta, gamma])
    if ltype == "ellipsoidal":
        # Use a matrix with only on-diagonal elements, with a = b
        alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
        x = (4.0 / 3.0) * np.pi
        for numattempts in range(maxattempts):
            vec = random_vector()
            c = vec[2] / (vec[0] * vec[1]) * np.cbrt(volume / x)
            a = b = np.sqrt((volume / x) / c)
            if (a / c < 10.0) and (c / a < 10.0):
                return np.array([a, b, c, alpha, beta, gamma])

    # If maxattempts tries have been made without success
    msg = "Could not get lattice after {:d} cycles for volume {:.2f}".format(maxattempts, volume)
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


#def para2matrix(cell_para, radians=True, format="lower"):
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
        matrix[2][2] = np.sqrt(c ** 2 - c1 ** 2 - c2 ** 2)
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
        matrix[0][0] = np.sqrt(a ** 2 - a3 ** 2 - a2 ** 2)
        #pass
    return matrix

def gaussian(min, max, sigma=3.0):
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
        sigma: the number of standard deviations between the center and min or max

    Returns:
        a value chosen randomly between min and max
    """
    center = (max + min) * 0.5
    delta = np.fabs(max - min) * 0.5
    ratio = delta / sigma
    while True:
        x = np.random.normal(scale=ratio, loc=center)
        if x > min and x < max:
            return x


def random_vector(minvec=[0.0, 0.0, 0.0], maxvec=[1.0, 1.0, 1.0], width=0.35, unit=False):
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
    vec = np.array(
        [
            np.exp(np.random.normal(scale=width)),
            np.exp(np.random.normal(scale=width)),
            np.exp(np.random.normal(scale=width)),
        ]
    )
    if unit:
        return vec / np.linalg.norm(vec)
    else:
        return vec



def random_shear_matrix(width=1.0, unitary=False):
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
    mat = np.zeros([3, 3])
    determinant = 0
    while determinant == 0:
        a, b, c = (
            np.random.normal(scale=width),
            np.random.normal(scale=width),
            np.random.normal(scale=width),
        )
        mat = np.array([[1, a, b], [a, 1, c], [b, c, 1]])
        determinant = np.linalg.det(mat)
    if unitary:
        new = mat / np.cbrt(np.linalg.det(mat))
        return new
    else:
        return mat

