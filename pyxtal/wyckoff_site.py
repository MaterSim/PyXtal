"""
Module for handling Wyckoff sites for both atom and molecule
"""

# Standard Libraries
import numpy as np
from scipy.spatial.transform import Rotation as R

# External Libraries
from pymatgen import Molecule

# PyXtal imports
from pyxtal.tolerance import Tol_matrix
from pyxtal.operations import apply_ops, distance, distance_matrix, project_point, filtered_coords, create_matrix
from pyxtal.symmetry import ss_string_from_ops as site_symm
from pyxtal.symmetry import Group
from pyxtal.database.element import Element
from pyxtal.constants import rad, deg
from pyxtal.msg import printx

class mol_site:
    """
    Class for storing molecular Wyckoff positions and orientations within
    the molecular_crystal class. Each mol_site object represenents an
    entire Wyckoff position, not necessarily a single molecule. This is the
    molecular version of Wyckoff_site

    Args:
        mol: a Pymatgen Molecule object
        position: the fractional 3-vector representing the generating molecule's position
        orientation: an Orientation object for the generating molecule
        wyckoff_position: a Wyckoff_position object
        lattice: a Lattice object for the crystal
        ellipsoid: an optional binding Ellipsoid object for checking distances.
        tm: a Tol_matrix object for distance checking
    """

    def __init__(
        self,
        mol,
        position,
        orientation,
        wyckoff_position,
        lattice,
        tols_matrix,
        radius,
        rotate_ref=True,
        ellipsoid=None,
    ):
        self.mol = mol
        """A Pymatgen molecule object"""
        self.position = wyckoff_position[0].operate(position)
        """Relative coordinates of the molecule's center within the unit cell"""
        self.orientation = orientation
        """The orientation object of the Mol in the first point in the WP"""
        self.ellipsoid = ellipsoid
        """A SymmOp representing the minimal ellipsoid for the molecule"""
        self.wp = wyckoff_position
        self.lattice = lattice
        self.inv_lattice = np.linalg.inv(lattice)
        """The crystal lattice in which the molecule resides"""
        self.multiplicity = self.wp.multiplicity
        """The multiplicity of the molecule's Wyckoff position"""
        self.PBC = wyckoff_position.PBC
        """The periodic axes"""
        self.symbols = [site.specie.value for site in self.mol.sites]
        self.numbers = self.mol.atomic_numbers
        self.tols_matrix = tols_matrix
        self.radius = radius
        if rotate_ref:
            tmp = self.mol.cart_coords
            tmp -= np.mean(tmp, axis=0)
            ax1 = self.get_principle_axes(tmp)
            self.coord0 = tmp.dot(ax1)
        else:
            self.coord0 = self.mol.cart_coords

    def __str__(self):

        if not hasattr(self, "site_symm"):
            self.site_symm = site_symm(
                self.wp.symmetry_m[0], self.wp.number, dim=self.wp.dim
            )
            #self.rotvec = self.orientation.r.as_rotvec()
            self.angles = self.orientation.r.as_euler('zxy', degrees=True)
        s = "{:} @ [{:6.4f} {:6.4f} {:6.4f}]  ".format(self.mol.formula.replace(" ",""), *self.position)
        s += "WP: {:2d}{:s}, ".format(self.wp.multiplicity, self.wp.letter)
        s += "Site symmetry {:} ==> Euler: ".format(self.site_symm)
        s += "{:6.3f} {:6.3f} {:6.3f}".format(*self.angles)
        return s

    def show(self, id=None, **kwargs):
        from pyxtal.viz import display_molecular_site
        return display_molecular_site(self, id, **kwargs)

    # NOTE appears deprecated?
    def get_ellipsoid(self):
        """
        Returns the bounding ellipsoid for the molecule. Applies the orientation
        transformation first.

        Returns:
            a re-orientated SymmOp representing the molecule's bounding ellipsoid
        """
        if self.ellipsoid is None:
            self.ellipsoid = find_ellipsoid(self.mol)
        e = self.ellipsoid
        # Apply orientation
        m = np.dot(e.rotation_matrix, self.orientation.get_matrix(angle=0))
        return SymmOp.from_rotation_and_translation(m, e.translation_vector)

    # NOTE appears deprecated?
    def get_ellipsoids(self):
        """
        Returns the bounding ellipsoids for the molecules in the WP. Includes the correct
        molecular centers and orientations.

        Returns:
            an array of re-orientated SymmOp's representing the molecule's bounding ellipsoids
        """
        # Get molecular centers
        centers0 = apply_ops(self.position, self.wp.generators)
        centers1 = np.dot(centers0, self.lattice)
        # Rotate ellipsoids
        e1 = self.get_ellipsoid()
        es = np.dot(self.wp.generators_m, e1)
        # Add centers to ellipsoids
        center_ops = [
            SymmOp.from_rotation_and_translation(Euclidean_lattice, c) for c in centers1
        ]
        es_final = []
        for e, c in zip(es, center_ops):
            es_final.append(e * c)
        return np.array(es_final)

    def _get_coords_and_species(self, absolute=False, add_PBC=False, first=False):
        """
        Used to generate coords and species for get_coords_and_species

        Args:
            absolute: whether or not to return absolute (Euclidean)
                coordinates. If false, return relative coordinates instead
            add_PBC: whether or not to add coordinates in neighboring unit cells, used for
                distance checking
            first: whether or not to extract the information from only the first site

        Returns:
            atomic coords: a numpy array of fractional coordinates for the atoms in the site
            species: a list of atomic species for the atomic coords
        """
        coord0 = self.coord0.dot(self.orientation.matrix.T)  #
        wp_atomic_sites = []
        wp_atomic_coords = None
        for point_index, op2 in enumerate(self.wp.generators):
            # Obtain the center in absolute coords
            center_relative = op2.operate(self.position)
            center_absolute = np.dot(center_relative, self.lattice)

            # Rotate the molecule (Euclidean metric)
            op2_m = self.wp.generators_m[point_index]
            rot = op2_m.affine_matrix[0:3][:, 0:3].T
            tau = op2_m.affine_matrix[0:3][:, 3]
            tmp = np.dot(coord0, rot) + tau
            # Add absolute center to molecule
            tmp += center_absolute
            tmp = tmp.dot(self.inv_lattice)
            if wp_atomic_coords is None:
                wp_atomic_coords = tmp
            else:
                wp_atomic_coords = np.append(wp_atomic_coords, tmp, axis=0)
            wp_atomic_sites.extend(self.symbols)

            if first:
                break

        if add_PBC is True:
            # Filter PBC of wp_atomic_coords
            wp_atomic_coords = filtered_coords(wp_atomic_coords, PBC=self.PBC)
            # Add PBC copies of coords
            m = create_matrix(PBC=self.PBC)
            # Move [0,0,0] PBC vector to first position in array
            m2 = [[0, 0, 0]]
            v0 = np.array([0.0, 0.0, 0.0])
            for v in m:
                if not (v == v0).all():
                    m2.append(v)
            new_coords = np.vstack([wp_atomic_coords + v for v in m2])
            wp_atomic_coords = new_coords

        if absolute:
            wp_atomic_coords = wp_atomic_coords.dot(self.lattice)

        return wp_atomic_coords, wp_atomic_sites

    def get_coords_and_species(self, absolute=False, add_PBC=False):
        """
        Lazily generates and returns the atomic coordinate and species for the
        Wyckoff position. Plugs the molecule into the provided orientation
        (with angle=0), and calculates the new positions.

        Args:
            absolute: whether or not to return absolute (Euclidean)
                coordinates. If false, return relative coordinates instead
            add_PBC: whether or not to add coordinates in neighboring unit cells, used for
                distance checking

        Returns:
            coords: a np array of 3-vectors.
            species: a list of atomic symbols, e.g. ['H', 'H', 'O', 'H', 'H', 'O']
        """
        if not absolute:
            # try to avoid repeating the computation
            try:
                return self.relative_coords, self.symbols_all
            except:
                self.relative_coords, self.symbols_all = self._get_coords_and_species(
                    absolute=absolute, add_PBC=add_PBC
                )
            return self.relative_coords, self.symbols_all
        else:
            try:
                return self.absolute_coords, self.symbols_all
            except:
                self.absolute_coords, self.symbols_all = self._get_coords_and_species(
                    absolute=absolute, add_PBC=add_PBC
                )
            return self.absolute_coords, self.symbols_all

    def get_centers(self, absolute=False):
        """
        Returns the fractional coordinates for the center of mass for each molecule in
        the Wyckoff position

        Returns:
            A numpy array of fractional 3-vectors
        """
        centers = apply_ops(self.position, self.wp.generators)
        # centers1 = filtered_coords(centers0, self.PBC)
        if absolute is False:
            return centers
        else:
            return np.dot(centers, self.lattice)

    def get_principle_axes(self, coords, adjust=False):
        """
        compute the principle axis
        """
        coords -= np.mean(coords, axis=0)
        Inertia = np.zeros([3,3])
        Inertia[0,0] = np.sum(coords[:,1]**2 + coords[:,2]**2)
        Inertia[1,1] = np.sum(coords[:,0]**2 + coords[:,2]**2)
        Inertia[2,2] = np.sum(coords[:,0]**2 + coords[:,1]**2)
        Inertia[0,1] = Inertia[1,0] = -np.sum(coords[:,0]*coords[:,1])
        Inertia[0,2] = Inertia[2,0] = -np.sum(coords[:,0]*coords[:,2])
        Inertia[1,2] = Inertia[2,1] = -np.sum(coords[:,1]*coords[:,2])
        _, matrix = np.linalg.eigh(Inertia)
        
        # search for the best direction
        if adjust:
            diffs = coords.dot(matrix) - self.coord0
            diffs = np.sqrt(np.sum(diffs**2,axis=0))/len(coords)
            for axis in range(3):
                if diffs[axis] > 0.05: #needs to check
                    matrix[:,axis] *= -1

            diffs = coords.dot(matrix) - self.coord0
            tol = np.sqrt(np.sum(diffs**2))/len(coords)
            if tol > 0.1:
                print("warining: molecular geometry changed")
                print(diffs)

        return matrix

    def get_Euler_angle(self):
        """
        To compute the Euler_angle for the given molecule
        """
        coord0 = self.coord0.dot(self.orientation.matrix.T)  #
        coord0 -= np.mean(coord0, axis=0)
        matrix = self.get_principle_axes(coord0, True)
        return R.from_matrix(matrix).as_euler('zxy', degrees=True)

    def translate(self, disp=np.array([0.0,0.0,0.0]), absolute=False):
        """
        To translate the molecule 
        Here we assume the molecule is free to rotate in SO(3)
        Needs to add the symmetry constraints later
        """
        disp = np.array(disp)
        if absolute:
            disp = disp.dot(self.inv_lattice)
        self.position += disp


    def rotate(self, axis=0, angle=180):
        """
        To rotate the molecule 
        Here we assume the molecule is free to rotate in SO(3)
        Needs to add the symmetry constraints later
        """
        p = self.orientation.r
        if type(axis) == int:
            coord0 = self.coord0.dot(self.orientation.r.as_matrix().T) 
            coord0 -= np.mean(coord0, axis=0)
            ax = self.get_principle_axes(coord0).T[axis]
        elif len(axis) == 3:
            ax = axis/np.linalg.norm(axis)

        q = R.from_rotvec(ax*rad*angle)
        o = q*p
        self.orientation.r = o 
        self.orientation.matrix = o.as_matrix()

    def get_mol_object(self, id=0):
        """
        make the pymatgen molecule object

        Args:
            id: the index of molecules in the given site

        Returns:
            a molecule object
        """
        coord0 = self.coord0.dot(self.orientation.matrix.T)  #
        # Obtain the center in absolute coords
        if id <= len(self.wp.generators):
            op = self.wp.generators[id]
            center_relative = op.operate(self.position)
            center_relative -= np.floor(center_relative)
            #print(center_relative)
            center_absolute = np.dot(center_relative, self.lattice)
            # Rotate the molecule (Euclidean metric)
            op_m = self.wp.generators_m[id]
            rot = op_m.affine_matrix[0:3][:, 0:3].T
            tau = op_m.affine_matrix[0:3][:, 3]
            tmp = np.dot(coord0, rot) + tau
            # Add absolute center to molecule
            tmp += center_absolute
            return Molecule(self.symbols, tmp)           
        else:
            raise ValueError("id is greater than the number of molecules")

    def update(self, coords, lattice=None, absolute=False):
        """
        After the geometry relaxation, the returned atomic coordinates
        maybe rescaled to [0, 1] bound. In this case, we need to refind
        the molecular coordinates according to the original neighbor list. 
        If the list does not change, we return the new coordinates
        otherwise, terminate the calculation.
        """
        from pymatgen.core.structure import Structure  
        from pyxtal.io import search_molecule_from_crystal
        from pyxtal.molecule import compare_mol_connectivity
        from openbabel import pybel, openbabel

        if lattice is not None:
            self.lattice = lattice
            self.inv_lattice = np.linalg.norm(lattice)
        if absolute:
            coords = coords.dot(self.inv_lattice)

        pmg = Structure(self.symbols, self.lattice, coords)
        coords, numbers = search_molecule_from_crystal(pmg, True)
        mol = Molecule(numbers, coords)
        match, _ = compare_mol_connectivity(mol, self.mol, True)
        if match:
            position = np.mean(coords, axis=0).dot(self.inv_lattice)
            position -= np.floor(position)
            self.position = position
            # orientation
            m1 = pybel.readstring('xyz', self.mol.to('xyz'))
            m2 = pybel.readstring('xyz', mol.to('xyz'))
            aligner = openbabel.OBAlign(True, False)
            aligner.SetRefMol(m1.OBMol)
            aligner.SetTargetMol(m2.OBMol)
            if aligner.Align():
                print("RMSD: ", aligner.GetRMSD())
                rot=np.zeros([3,3])
                for i in range(3):
                    for j in range(3):
                        rot[i,j] = aligner.GetRotMatrix().Get(i,j)
                if abs(np.linalg.det(rot) - 1) < 1e-2:
                    self.orientation.matrix = rot
                    self.orientation.r = R.from_matrix(rot)
                else:
                    raise ValueError("rotation matrix is wrong")
        else:
            raise ValueError("molecular connectivity changes! Exit")
        #todo check if connectivty changed
        


    def create_matrix(self):
        """
        Used for calculating distances in lattices with periodic boundary
        conditions. When multiplied with a set of points, generates additional
        points in cells adjacent to and diagonal to the original cell
        Returns:
            A numpy array of matrices which can be multiplied by a set of
            coordinates
        """
        matrix = []
        [a, b, c] = np.linalg.norm(self.lattice, axis=1)
        if a > 20 and self.radius<10:
            i_list = [0]
        elif a < 6.5:
            i_list = [-1,0,2]
        else:
            i_list = [-1, 0, 1]
            
        if b > 20 and self.radius<10:
            j_list = [0]
        elif b < 6.5:
            j_list = [-1, 0, 2]
        else:
            j_list = [-1, 0, 1]
            
        if c > 20 and self.radius<10:
            k_list = [0]
        elif c < 6.5:
            k_list = [-1, 0, 2]
        else:
            k_list = [-1, 0, 1]
        
        if not self.PBC[0]:
            i_list = [0]
        if not self.PBC[1]:
            j_list = [0]
        if not self.PBC[2]:
            k_list = [0]
        for i in i_list:
            for j in j_list:
                for k in k_list:
                    matrix.append([i, j, k])
        return np.array(matrix, dtype=float)

    def compute_distances(self):
        """
        compute if the atoms in the Wyckoff position are too close to each other
        or not. Does not check distances between atoms in the same molecule. Uses
        crystal.check_distance as the base code.

        Returns:
            True if the atoms are not too close together, False otherwise
        """
        m_length = len(self.symbols)
        # TODO: Use tm instead of tols lists
        # Get coords of WP with PBC
        coords, _ = self._get_coords_and_species()

        # Get coords of the generating molecule
        coords_mol = coords[:m_length]
        # Remove generating molecule's coords from large array
        coords = coords[m_length:]
        min_ds = []

        if self.PBC != [0, 0, 0]:
            # Check periodic images
            m = self.create_matrix()
            # Remove original coordinates
            m2 = []
            v0 = np.array([0.0, 0.0, 0.0])
            for v in m:
                if not (v == v0).all():
                    m2.append(v)
            if len(m2) > 0:
                coords_PBC = np.vstack([coords_mol + v for v in m2])
                d = distance_matrix(coords_mol, coords_PBC, self.lattice, [0, 0, 0], True)
                min_ds.append(d)
        if self.multiplicity > 1:
            # Check inter-atomic distances
            d = distance_matrix(coords_mol, coords, self.lattice, self.PBC, True)
            min_ds.append(d)
        return min(min_ds)

    def check_distances(self, atomic=True):
        """
        Checks if the atoms in the Wyckoff position are too close to each other
        or not. Does not check distances between atoms in the same molecule. Uses
        crystal.check_distance as the base code.

        Args:
            atomic: if True, checks inter-atomic distances. If False, checks ellipsoid
                overlap between molecules instead

        Returns:
            True if the atoms are not too close together, False otherwise
        """
        if atomic:
            m_length = len(self.symbols)
            # TODO: Use tm instead of tols lists
            # Get coords of WP with PBC
            coords, _ = self._get_coords_and_species()

            # Get coords of the generating molecule
            coords_mol = coords[:m_length]
            # Remove generating molecule's coords from large array
            coords = coords[m_length:]

            if self.PBC != [0, 0, 0]:
                # Check periodic images
                m = self.create_matrix()
                # Remove original coordinates
                m2 = []
                v0 = np.array([0.0, 0.0, 0.0])
                for v in m:
                    if not (v == v0).all():
                        m2.append(v)
                if len(m2) > 0:
                    coords_PBC = np.vstack([coords_mol + v for v in m2])
                    d = distance_matrix(coords_PBC, coords_mol, self.lattice, PBC=[0, 0, 0])
                    # only check if small distance is detected
                    if np.min(d) < np.max(self.tols_matrix):
                        tols = np.min(d.reshape([len(m2), m_length, m_length]), axis=0)
                        if (tols < self.tols_matrix).any():
                            return False

            if self.multiplicity > 1:
                # Check inter-atomic distances
                d = distance_matrix(coords, coords_mol, self.lattice, PBC=self.PBC)
                if np.min(d) < np.max(self.tols_matrix):
                    tols = np.min(
                        d.reshape([self.multiplicity - 1, m_length, m_length]), axis=0
                    )
                    if (tols < self.tols_matrix).any():
                        return False

            return True

            """New method - only checks some atoms/molecules"""
            # #Store length of molecule
            # m_length = len(self.mol)
            # #Get coordinates of center molecule and Wyckoff position
            # coords, species = self._get_coords_and_species(absolute=True)
            # coords_mol = coords[:m_length]

            # if self.PBC == [0,0,0]:
            #    #Check non-periodic Wyckoff positions
            #    if self.multiplicity == 1:
            #        return True
            #    coords_other = coords[m_length:]
            #    tols = np.repeat(self.tols_matrix, self.multiplicity-1, axis=1)
            #    d = cdist(coords_mol, coords_other)
            #    if (d<tols).any():
            #        return False
            #    else:
            #        return True

            # #Create PBC vectors
            # m = create_matrix(PBC=self.PBC)
            # ml = np.dot(m, self.lattice)
            # #Store the index of the (0,0,0) vector within ml
            # mid_index = len(ml) // 2

            # if self.multiplicity == 1:
            #    #Only check periodic images
            #    #Remove original coordinates
            #    m2 = []
            #    v0 = np.array([0.,0.,0.])
            #    for v in ml:
            #        if not (v==v0).all():
            #            m2.append(v)
            #    coords_PBC = np.vstack([coords_mol + v for v in m2])
            #    d = distance_matrix(coords_mol, coords_PBC, None, PBC=[0,0,0])
            #    tols = np.repeat(self.tols_matrix, len(m2), axis=1)
            #    if (d<tols).any():
            #        return False
            #    else:
            #        return True

            # #Generate centers of all molecules
            # centers = self.get_centers(absolute=True)
            # vectors = np.repeat(centers, len(ml), axis=0) + np.tile(ml,(len(centers),1)) - np.dot(self.position, self.lattice)
            # #Calculate distances between centers
            # distances = np.linalg.norm(vectors, axis=-1)
            # #Find which molecules need to be checked
            # indices_mol = np.where(distances < self.radius()*2)[0]
            # #Get indices of Wyckoff positions and PBC vectors
            # indices_wp = []
            # indices_pbc = []
            # indices_vector = []
            # for index in indices_mol:
            #    i_wp, i_pbc = divmod(index, len(ml))
            #    #Omit original center molecule
            #    if not (i_wp == 0 and i_pbc == mid_index):
            #        indices_wp.append(i_wp)
            #        indices_pbc.append(i_pbc)
            #        indices_vector.append(index)

            # if indices_wp == []:
            #    return True

            # #Get atomic positions of molecules with small separation vectors
            # original_coords = np.vstack([coords[index_wp*m_length:index_wp*m_length+m_length] for index_wp in indices_wp])
            # pbc_toadd = np.repeat(ml[indices_pbc], m_length, axis=0)
            # atomic_coords = original_coords + pbc_toadd
            # #Get inter-atomic tolerances
            # #tols = np.tile(self.get_tols_matrix(), len(indices_wp))
            # tols = np.tile(self.tols_matrix, len(indices_wp))
            # if m_length <= max_fast_mol_size:
            #    #Check all atomic pairs
            #    d = cdist(coords_mol, atomic_coords)

            #    """
            #    print("~~~~~~~~~~~~~~~~~~~~~~~")
            #    print("ml:", ml.shape)
            #    print(ml)
            #    print("centers:", centers.shape)
            #    print(centers)
            #    print("vectors:", vectors.shape)
            #    print(vectors)
            #    print("radius*2: ", self.get_radius()*2)
            #    print("distances:", distances.shape)
            #    print(distances)
            #    print("indices_mol:", len(indices_mol))
            #    print(indices_mol)
            #    print("indices_wp:", len(indices_wp))
            #    print(indices_wp)
            #    print("indices_pbc:", len(indices_pbc))
            #    print(indices_pbc)
            #    print("indices_vector:", len(indices_vector))
            #    print(indices_vector)
            #    print("coords_mol:", coords_mol.shape)
            #    print(coords_mol)
            #    print("coords:", coords.shape)
            #    print(coords)
            #    print("original_coords:", original_coords.shape)
            #    print(original_coords)
            #    print("pbc_toadd:", pbc_toadd.shape)
            #    print(pbc_toadd[:12])
            #    print("atomic_coords: ", atomic_coords.shape)
            #    print(atomic_coords)
            #    print("d:", d.shape)
            #    print(d)
            #    print("tols_matrix:", self.get_tols_matrix().shape)
            #    print(self.get_tols_matrix())
            #    print("tols:", tols.shape)
            #    print(tols)
            #    """

            #    if (d<tols).any():
            #        return False
            #    else:
            #        return True

            # elif m_length > max_fast_mol_size:
            #    #Get corresponding separation vectors
            #    new_vectors = np.repeat(vectors[indices_vector], m_length, axis=0)
            #    #Get atomic coordinates relative to molecular centers
            #    relative_atomic_coords = atomic_coords - new_vectors
            #    #Dot atomic coordinates with inter-molecular separation vectors
            #    dots = np.einsum('...j,...j', new_vectors, relative_atomic_coords)
            #    #Find where relative vectors point towards the original molecule
            #    new_indices = np.where(dots<0)[0]
            #    #Get new coordinates and tolerances for distance matrix
            #    new_atomic_coords = atomic_coords[new_indices]
            #    d = cdist(coords_mol, new_atomic_coords)
            #    tols2 = tols[:,new_indices]
            #    if (d<tols2).any():
            #        return False
            #    else:
            #        return True

        else:
            # Check molecular ellipsoid overlap
            if self.multiplicity == 1:
                return True
            es0 = self.get_ellipsoids()[1:]
            PBC_vectors = np.dot(create_matrix(PBC=self.PBC), self.lattice)
            PBC_ops = [
                SymmOp.from_rotation_and_translation(Euclidean_lattice, v)
                for v in PBC_vectors
            ]
            es1 = []
            for op in PBC_ops:
                es1.append(np.dot(es0, op))
            es1 = np.squeeze(es1)
            truth_values = np.vectorize(check_intersection)(es1, self.get_ellipsoid())
            if np.sum(truth_values) < len(truth_values):
                return False
            else:
                return True


def check_intersection(ellipsoid1, ellipsoid2):
    """
    Given SymmOp's for 2 ellipsoids, checks whether or not they overlap

    Args:
        ellipsoid1: a SymmOp representing the first ellipsoid
        ellipsoid2: a SymmOp representing the second ellipsoid

    Returns:
        False if the ellipsoids overlap.
        True if they do not overlap.
    """
    # Transform so that one ellipsoid becomes a unit sphere at (0,0,0)
    Op = ellipsoid1.inverse * ellipsoid2
    # We define a new ellipsoid by moving the sphere around the old ellipsoid
    M = Op.rotation_matrix
    a = 1.0 / (1.0 / np.linalg.norm(M[0]) + 1)
    M[0] = M[0] / np.linalg.norm(M[0]) * a
    b = 1.0 / (1.0 / np.linalg.norm(M[1]) + 1)
    M[1] = M[1] / np.linalg.norm(M[1]) * b
    c = 1.0 / (1.0 / np.linalg.norm(M[2]) + 1)
    M[2] = M[2] / np.linalg.norm(M[2]) * c
    p = Op.translation_vector
    # Calculate the transformed distance from the sphere's center to the new ellipsoid
    dsq = np.dot(p, M[0]) ** 2 + np.dot(p, M[1]) ** 2 + np.dot(p, M[2]) ** 2
    if dsq < 2:
        return False
    else:
        return True


def check_mol_sites(
    ms1, ms2, atomic=False, factor=1.0, tm=Tol_matrix(prototype="molecular")
):
    """
    Checks whether or not the molecules of two mol sites overlap. Uses
    ellipsoid overlapping approximation to check. Takes PBC and lattice
    into consideration.

    Args:
        ms1: a mol_site object
        ms2: another mol_site object
        atomic: if True, checks inter-atomic distances. If False, checks
            overlap between molecular ellipsoids
        factor: the distance factor to pass to check_distances. (only for
            inter-atomic distance checking)
        tm: a Tol_matrix object (or prototype string) for distance checking

    Returns:
        False if the Wyckoff positions overlap. True otherwise
    """
    if atomic is False:
        es0 = ms1.get_ellipsoids()
        PBC_vectors = np.dot(create_matrix(PBC=ms1.PBC), ms1.lattice)
        PBC_ops = [
            SymmOp.from_rotation_and_translation(Euclidean_lattice, v)
            for v in PBC_vectors
        ]
        es1 = []
        for op in PBC_ops:
            es1.append(np.dot(es0, op))
        es1 = np.squeeze(es1)
        truth_values = np.vectorize(check_intersection)(es1, ms2.get_ellipsoid())
        if np.sum(truth_values) < len(truth_values):
            return False
        else:
            return True

    elif atomic is True:
        # Get coordinates for both mol_sites
        c1, _ = ms1.get_coords_and_species()
        c2, _ = ms2.get_coords_and_species()

        # Calculate which distance matrix is smaller/faster
        m_length1 = len(ms1.numbers)
        m_length2 = len(ms2.numbers)
        wp_length1 = len(c1)
        wp_length2 = len(c2)
        size1 = m_length1 * wp_length2
        size2 = m_length2 * wp_length1

        # Case 1
        if size1 <= size2:
            coords_mol = c1[:m_length1]
            # Calculate tol matrix for species pairs
            tols = np.zeros((m_length1, m_length2))
            for i1, number1 in enumerate(ms1.numbers):
                for i2, number2 in enumerate(ms2.numbers):
                    tols[i1][i2] = tm.get_tol(number1, number2)
            tols = np.repeat(tols, ms2.multiplicity, axis=1)
            d = distance_matrix(coords_mol, c2, ms1.lattice, PBC=ms1.PBC)

        # Case 2
        elif size1 > size2:
            coords_mol = c2[:m_length2]
            # Calculate tol matrix for species pairs
            tols = np.zeros((m_length2, m_length1))
            for i1, number1 in enumerate(ms2.numbers):
                for i2, number2 in enumerate(ms1.numbers):
                    tols[i1][i2] = tm.get_tol(number1, number2)
            tols = np.repeat(tols, ms1.multiplicity, axis=1)
            d = distance_matrix(coords_mol, c1, ms1.lattice, PBC=ms1.PBC)

        # Check if distances are smaller than tolerances
        if (d < tols).any():
            return False
        return True


class atom_site:
    """
    Class for storing atomic Wyckoff positions with a single coordinate.

    Args:
        wp: a Wyckoff_position object
        coordinate: a fractional 3-vector for the generating atom's coordinate
        specie: an Element, element name or symbol, or atomic number of the atom
    """

    def __init__(self, wp, coordinate, specie=1):
        self.position = np.array(coordinate)
        self.specie = Element(specie).short_name
        self.multiplicity = wp.multiplicity
        self.wp = wp
        self.PBC = wp.PBC
        self.update_coords(coordinate)

    def __str__(self):
        if not hasattr(self, "site_symm"):
            self.site_symm = site_symm(
                self.wp.symmetry_m[0], self.wp.number, dim=self.wp.dim
            )

        s = "{:>2s} @ [{:6.4f} {:6.4f} {:6.4f}], ".format(self.specie, *self.position)
        s += "WP: {:2d}{:s}, ".format(self.wp.multiplicity, self.wp.letter)
        s += "Site symmetry: {:s}".format(self.site_symm)
        return s

    def update_coords(self, pos):
        """
        Used to generate coords from self.position
        """
        self.position = pos
        self.coords = apply_ops(pos, self.wp) 

    def __repr__(self):
        return str(self)


def check_atom_sites(ws1, ws2, lattice, tm, same_group=True):
    """
    Given two Wyckoff sites, checks the inter-atomic distances between them.

    Args:
        ws1: a Wyckoff_site object
        ws2: a different Wyckoff_site object (will always return False if
            two identical WS's are provided)
        lattice: a 3x3 cell matrix
        same_group: whether or not the two WS's are in the same structure.
            Default value True reduces the calculation cost

    Returns:
        True if all distances are greater than the allowed tolerances.
        False if any distance is smaller than the allowed tolerance
    """
    # Ensure the PBC values are valid
    if ws1.PBC != ws2.PBC:
        printx("Error: PBC values do not match between Wyckoff sites")
        return
    # Get tolerance
    tol = tm.get_tol(ws1.specie, ws2.specie)
    # Symmetry shortcut method: check only some atoms
    if same_group is True:
        # We can either check one atom in WS1 against all WS2, or vice-versa
        # Check which option is faster
        if ws1.multiplicity > ws2.multiplicity:
            coords1 = [ws1.coords[0]]
            coords2 = ws2.coords
        else:
            coords1 = [ws2.coords[0]]
            coords2 = ws1.coords
        # Calculate distances
        dm = distance_matrix(coords1, coords2, lattice, PBC=ws1.PBC)
        # Check if any distances are less than the tolerance
        if (dm < tol).any():
            return False
        else:
            return True
    # No symmetry method: check all atomic pairs
    else:
        dm = distance_matrix(ws1.coords, ws2.coords, lattice, PBC=ws1.PBC)
        # Check if any distances are less than the tolerance
        if (dm < tol).any():
            return False
        else:
            return True

def WP_merge_old(coor, lattice, group, tol):
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

    Returns:
        coor: the new list of fractional coordinates after merging. 
        index: a single index for the Wyckoff position within the sg. 
        If no matching WP is found, returns False. 
        point: is a 3-vector when plugged into the Wyckoff position,
    """
    coor = np.array(coor)

    # Get index of current Wyckoff position. If not one, return False
    index, point = check_wyckoff_position(coor, group)
    if index is False:
        return coor, False, None
    if point is None:
        printx("Error: Could not find generating point.", priority=1)
        printx("coordinates:")
        printx(str(coor))
        printx("Lattice: ")
        printx(str(lattice))
        printx("group: ")
        group.print_all()
        return coor, False, None
    PBC = group.PBC
    # Main loop for merging multiple times
    while True:
        # Check distances of current WP. If too small, merge
        dm = distance_matrix([coor[0]], coor, lattice, PBC=PBC)
        passed_distance_check = True
        x = np.argwhere(dm < tol)
        for y in x:
            # Ignore distance from atom to itself
            if y[0] == 0 and y[1] == 0:
                pass
            else:
                passed_distance_check = False
                break

        if passed_distance_check is False:
            mult1 = group[index].multiplicity
            # Find possible wp's to merge into
            possible = []
            for i, wp in enumerate(group):
                mult2 = wp.multiplicity
                # factor = mult2 / mult1
                if (mult2 < mult1) and (mult1 % mult2 == 0):
                    possible.append(i)
            if possible == []:
                return coor, False, None
            # Calculate minimum separation for each WP
            distances = []
            for i in possible:
                wp = group[i]
                projected_point = project_point(point, wp[0], lattice=lattice, PBC=PBC)
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

def WP_merge(pt, lattice, wp, tol):
    """
    Given a list of fractional coordinates, merges them within a given
    tolerance, and checks if the merged coordinates satisfy a Wyckoff
    position. Used for merging general Wyckoff positions into special Wyckoff
    positions within the random_crystal (and its derivative) classes.

    Args:
        pt: the originl point (3-vector)
        lattice: a 3x3 matrix representing the unit cell
        wp: a pyxtal.symmetry.Wyckoff_position object after merge
        tol: the cutoff distance for merging coordinates

    Returns:
        pt: 3-vector after merge
        wp: a pyxtal.symmetry.Wyckoff_position object
        If no matching WP is found, returns False. 
    """
    index = wp.index
    PBC = wp.PBC
    group = Group(wp.number, wp.dim)

    pt = project_point(pt, wp[0], lattice, PBC)
    coor = apply_ops(pt, wp)
    # Main loop for merging multiple times

    while True:
        # Check distances of current WP. If too small, merge
        dm = distance_matrix([coor[0]], coor, lattice, PBC=PBC)
        passed_distance_check = True
        x = np.argwhere(dm < tol)
        for y in x:
            # Ignore distance from atom to itself
            if y[0] == 0 and y[1] == 0:
                pass
            else:
                passed_distance_check = False
                break
        
        if not passed_distance_check:
            mult1 = group[index].multiplicity
            # Find possible wp's to merge into
            possible = []
            for i, wp0 in enumerate(group):
                mult2 = wp0.multiplicity
                # factor = mult2 / mult1
                if (mult2 < mult1) and (mult1 % mult2 == 0):
                    possible.append(i)
            if possible == []:
                return None, False
        
            # Calculate minimum separation for each WP
            distances = []
            for i in possible:
                wp = group[i]
                projected_point = project_point(pt, wp[0], lattice=lattice, PBC=PBC)
                d = distance(pt - projected_point, lattice, PBC=PBC)
                distances.append(np.min(d))
            # Choose wp with shortest translation for generating point
            tmpindex = np.argmin(distances)
            index = possible[tmpindex]
            wp = group[index]
            pt = project_point(pt, wp[0], lattice=lattice, PBC=PBC)
            coor = apply_ops(pt, wp)
        # Distances were not too small; return True
        else:
            return pt, wp


