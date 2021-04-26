"""
Module for handling Wyckoff sites for both atom and molecule
"""

# Standard Libraries
import numpy as np
from scipy.spatial.transform import Rotation as R

# External Libraries
from pymatgen.core import Molecule
from pymatgen.core.structure import Structure  

# PyXtal imports
from pyxtal.tolerance import Tol_matrix
from pyxtal.operations import (
    apply_ops, 
    check_images,
    distance, 
    distance_matrix, 
    project_point, 
    filtered_coords, 
    create_matrix,
    SymmOp,
)
from pyxtal.symmetry import Group, jk_from_i, Wyckoff_position
from pyxtal.symmetry import ss_string_from_ops as site_symm
from pyxtal.database.element import Element
from pyxtal.constants import rad, deg
from pyxtal.lattice import Lattice
from pyxtal.msg import printx

class mol_site:
    """
    Class for storing molecular Wyckoff positions and orientations within
    the molecular_crystal class. Each mol_site object represenents an
    entire Wyckoff position, not necessarily a single molecule. This is the
    molecular version of Wyckoff_site

    Args:
        mol: a `pyxtal_molecule <pyxtal.molecule.pyxtal_molecule.html>`_ object
        position: the 3-vector representing the generating molecule's position
        orientation: an `Orientation <pyxtal.molecule.Oreintation.html>`_ object 
        wp: a `Wyckoff_position <pyxtal.symmetry.Wyckoff_position.html>`_ object
        lattice: a `Lattice <pyxtal.lattice.Lattice>`_ object 
        diag: whether or not use the `n` representation
    """

    def __init__(self, mol, position, orientation, wp, lattice=None, diag=False):
        # describe the molecule
        self.molecule = mol
        self.diag = diag
        self.wp = wp
        self.position = position # fractional coordinate of molecular center
        self.orientation = orientation #pyxtal.molecule.orientation object 
        if isinstance(lattice, Lattice):
            self.lattice = lattice
        else:
            self.lattice = Lattice.from_matrix(lattice)
        self.PBC = self.wp.PBC
        self.mol = mol.mol # A Pymatgen molecule object
        self.symbols = mol.symbols #[site.specie.value for site in self.mol.sites]
        self.numbers = self.mol.atomic_numbers
        self.tols_matrix = mol.tols_matrix
        self.radius = mol.radius

        if self.diag:
            self.wp.diagonalize_symops()
            self.position = project_point(self.position, wp[0])

    def __str__(self):

        if not hasattr(self, "site_symm"):
            self.site_symm = site_symm(
                self.wp.symmetry_m[0], self.wp.number, dim=self.wp.dim
            )
        self.angles = self.orientation.r.as_euler('zxy', degrees=True)
        formula = self.mol.formula.replace(" ","")
        s = "{:} @ [{:7.4f} {:7.4f} {:7.4f}]  ".format(formula, *self.position)
        s += "WP: {:2d}{:s}, ".format(self.wp.multiplicity, self.wp.letter)
        s += "Site symmetry {:} ==> Euler: ".format(self.site_symm)
        s += "{:6.3f} {:6.3f} {:6.3f}".format(*self.angles)

        return s

    def _get_dof(self):
        """
        get the number of dof for the given structures:
        """
        freedom = np.trace(self.wp.ops[0].rotation_matrix) > 0
        self.dof = len(freedom[freedom==True])


    def save_dict(self):
        dict0 = {"position": self.position,
                 "number": self.wp.number,
                 "dim": self.wp.dim,
                 "index": self.wp.index,
                 "diag": self.diag,
                 "molecule": self.molecule.save_dict(),
                 "orientation": self.orientation.save_dict(),
                 "lattice": self.lattice.matrix,
                 "lattice_type": self.lattice.ltype,
                }
        if self.molecule.torsionlist is not None:
            xyz, _ = self._get_coords_and_species(absolute=True, first=True)
            dict0['rotors'] = molecule.get_torsion_angles(xyz)

        return dict0

    @classmethod
    def load_dict(cls, dicts):
        """
        load the sites from a dictionary
        """
        from pyxtal.molecule import pyxtal_molecule, Orientation

        mol = pyxtal_molecule.load_dict(dicts["molecule"])
        g = dicts["number"]
        index = dicts["index"]
        dim = dicts["dim"]
        position = dicts["position"]
        orientation = Orientation.load_dict(dicts['orientation'])
        wp = Wyckoff_position.from_group_and_index(g, index, dim)
        diag = dicts["diag"]
        lattice = Lattice.from_matrix(dicts["lattice"], ltype=dicts["lattice_type"])
        return cls(mol, position, orientation, wp, lattice, diag)

    def encode(self):
        """
        transform dict to 1D vector
        [x, y, z, or1, or2, or3, rotor1, rotor2, .etc]
        """
        xyz, _ = self._get_coords_and_species(absolute=True, first=True)
        xyz -= self.molecule.get_center(xyz)
        rotor = self.molecule.get_torsion_angles(xyz)
        ori, _, reflect = self.molecule.get_orientation(xyz)
        #print(self.molecule.mol)
        return list(self.position) + list(ori) + rotor + [reflect]
        
    def to_1D_dicts(self):
        """
        save the wp in 1D representation
        """
        xyz, _ = self._get_coords_and_species(absolute=True, first=True)
        dict0 = {"smile": self.molecule.smile}
        dict0["rotor"] = self.molecule.get_torsion_angles(xyz)
        dict0["orientation"], dict0["rmsd"], dict0["reflect"] = self.molecule.get_orientation(xyz)
        
        angs = dict0["rotor"] 
        rdkit_mol = self.molecule.rdkit_mol(self.molecule.smile)
        conf0 = rdkit_mol.GetConformer(0)
        #print(self.molecule.set_torsion_angles(conf0, angs))
        #import sys; sys.exit()

        #print("save matrix"); print(self.orientation.r.as_matrix())
        #print("save angle"); print(self.orientation.r.as_euler('zxy', degrees=True))
        #print("angle"); print(dict0["orientation"])
        dict0["center"] = self.position - np.floor(self.position) #self.molecule.get_center(xyz)
        dict0["number"] = self.wp.number
        dict0["index"] = self.wp.index
        dict0["dim"] = self.wp.dim
        dict0["diag"] = self.diag
        dict0["lattice"] = self.lattice.get_matrix()
        dict0["lattice_type"] = self.lattice.ltype

        return dict0

    @classmethod
    def from_1D_dicts(cls, dicts):
        from pyxtal.molecule import pyxtal_molecule, Orientation

        mol = pyxtal_molecule(mol=dicts['smile']+'.smi')
        rdkit_mol = mol.rdkit_mol(mol.smile)
        conf = rdkit_mol.GetConformer(0)
        #print("try")
        #print(conf.GetPositions()[:3])
        #print(dicts["rotor"])
        if dicts['reflect']:
            mol.set_torsion_angles(conf, dicts["rotor"], False)
        #    print(mol.set_torsion_angles(conf, dicts["rotor"], True))
        #    #import sys; sys.exit()
        xyz = mol.set_torsion_angles(conf, dicts["rotor"], dicts['reflect'])
        mol.reset_positions(xyz)
        g = dicts["number"]
        index = dicts["index"]
        dim = dicts["dim"]
        matrix = R.from_euler('zxy', dicts["orientation"], degrees=True).as_matrix()
        orientation = Orientation(matrix)
        #if dicts['reflect']:
        #    print('load'); print(xyz[:3])
        #    print("aaaaaaaaaaaaaa"); print(xyz[:3].dot(orientation.matrix.T))
        #    print("matrix"); print(orientation.matrix)
        wp = Wyckoff_position.from_group_and_index(g, index, dim)
        diag = dicts["diag"]
        lattice = Lattice.from_matrix(dicts["lattice"], ltype=dicts["lattice_type"])
        position = dicts["center"] #np.dot(dicts["center"], lattice.inv_matrix)

        return cls(mol, position, orientation, wp, lattice, diag)

    def show(self, id=None, **kwargs):
        from pyxtal.viz import display_molecular_site
        return display_molecular_site(self, id, **kwargs)

    def _get_coords_and_species(self, absolute=False, PBC=False, first=False, unitcell=False):
        """
        Used to generate coords and species for get_coords_and_species

        Args:
            absolute: whether or not to return absolute (Euclidean)
                coordinates. If false, return relative coordinates instead
            PBC: whether or not to add coordinates in neighboring unit cells, 
                used for distance checking
            first: whether or not to extract the information from only the first site
            unitcell: whether or not to move the molecular center to the unit cell

        Returns:
            atomic coords: a numpy array of fractional coordinates for the atoms in the site
            species: a list of atomic species for the atomic coords
        """
        coord0 = self.mol.cart_coords.dot(self.orientation.matrix.T)  #
        wp_atomic_sites = []
        wp_atomic_coords = None

        for point_index, op2 in enumerate(self.wp.ops):
            # Obtain the center in absolute coords
            center_relative = op2.operate(self.position)
            if unitcell:
                center_relative -= np.floor(center_relative)
            center_absolute = np.dot(center_relative, self.lattice.matrix)

            # Rotate the molecule (Euclidean metric)
            op2_m = self.wp.generators_m[point_index]
            rot = op2_m.affine_matrix[:3, :3].T
            if self.diag and self.wp.index > 0:
                tau = op2.translation_vector
            else:
                tau = op2_m.affine_matrix[:3, 3]
            tmp = np.dot(coord0, rot) + tau
            # Add absolute center to molecule
            tmp += center_absolute
            tmp = tmp.dot(self.lattice.inv_matrix)
            if wp_atomic_coords is None:
                wp_atomic_coords = tmp
            else:
                wp_atomic_coords = np.append(wp_atomic_coords, tmp, axis=0)
            wp_atomic_sites.extend(self.symbols)

            if first:
                break

        if PBC:
            # Filter PBC of wp_atomic_coords
            wp_atomic_coords = filtered_coords(wp_atomic_coords, PBC=self.PBC)
            # Add PBC copies of coords
            m = create_matrix(PBC=self.PBC, omit=True)
            # Move [0,0,0] PBC vector to first position in array
            m2 = [[0, 0, 0]]
            for v in m:
                m2.append(v)
            new_coords = np.vstack([wp_atomic_coords + v for v in m2])
            wp_atomic_coords = new_coords

        if absolute:
            wp_atomic_coords = wp_atomic_coords.dot(self.lattice.matrix)

        return wp_atomic_coords, wp_atomic_sites

    def get_coords_and_species(self, absolute=False, PBC=False, unitcell=False):
        """
        Lazily generates and returns the atomic coordinate and species for the
        Wyckoff position. Plugs the molecule into the provided orientation
        (with angle=0), and calculates the new positions.

        Args:
            absolute: whether or not to return absolute (Euclidean)
                coordinates. If false, return relative coordinates instead
            PBC: whether or not to add coordinates in neighboring unit cells, used for
                distance checking
            unitcell: whether or not to move the molecular center to the unit cell

        Returns:
            coords: a np array of 3-vectors.
            species: a list of atomic symbols, e.g. ['H', 'H', 'O', 'H', 'H', 'O']
        """
        return self._get_coords_and_species(absolute, PBC, unitcell=unitcell)

    #def get_centers(self, absolute=False):
    #    """
    #    Get the coordinates for the center of mass 

    #    Returns:
    #        A numpy array of fractional 3-vectors
    #    """
    #    centers = apply_ops(self.position, self.wp.ops)
    #    # centers1 = filtered_coords(centers0, self.PBC)
    #    if not absolute:
    #        return centers
    #    else:
    #        return np.dot(centers, self.lattice.matrix)

    def perturbate(self, lattice, trans=0.1, rot=5):
        """
        Random perturbation of the molecular site
        
        Args:
            lattice: lattice vectors
            trans: magnitude of tranlation vectors (default: 0.1 A)
            rot: magnitude of rotation degree (default: 5.0)
        """
        dis = (np.random.sample(3) - 0.5).dot(lattice)
        dis /= np.linalg.norm(dis)
        dis *= trans
        self.translate(dis, True)
        if rot == 'random':
            self.orientation.change_orientation()
        else:
            self.orientation.change_orientation(angle=rot/180*np.pi)
    
    def translate(self, disp=np.zeros(3), absolute=False):
        """
        To translate the molecule 
        """
        disp = np.array(disp)
        if absolute:
            disp = disp.dot(self.lattice.inv_matrix)
        position = self.position + disp
        self.position = project_point(position, self.wp[0])

    def rotate(self, ax_id=0, ax_vector=None, angle=180):
        """
        To rotate the molecule
        Args:
            ax_id: the principle axis id
            ax_vector (float): 3-vector to define the axis
            angle (float): angle to rotate
        """
        p = self.orientation.r

        if ax_vector is not None:
            ax = ax_vector/np.linalg.norm(ax_vector)
        else:
            xyz = self.mol.cart_coords.dot(p.as_matrix().T) 
            ax = self.molecule.get_principle_axes(xyz).T[axis]

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
        coord0 = self.mol.cart_coords.dot(self.orientation.matrix.T)  #
        # Obtain the center in absolute coords
        if id <= len(self.wp.generators):
            op = self.wp.generators[id]
            center_relative = op.operate(self.position)
            center_relative -= np.floor(center_relative)
            #print(center_relative)
            center_absolute = np.dot(center_relative, self.lattice.matrix)
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

    def update(self, coords, lattice=None, absolute=False, update_mol=True):
        """
        After the geometry relaxation, the returned atomic coordinates
        maybe rescaled to [0, 1] bound. In this case, we need to refind
        the molecular coordinates according to the original neighbor list. 
        If the list does not change, we return the new coordinates
        otherwise, terminate the calculation.
        """
        from pyxtal.molecule import compare_mol_connectivity, Orientation
        try:
            from openbabel import pybel, openbabel
        except:
            import pybel, openbabel
        if lattice is not None:
            self.lattice = lattice
        if not absolute:
            coords = coords.dot(self.lattice.matrix)
        #mol = Molecule(self.symbols, coords-np.mean(coords, axis=0))
        center = self.molecule.get_center(coords)
        mol = Molecule(self.symbols, coords-center)

        #match, _ = compare_mol_connectivity(mol, self.mol, True)
        match, _ = compare_mol_connectivity(mol, self.mol)
        if match:
            #position = np.mean(coords, axis=0).dot(self.lattice.inv_matrix)
            position = center.dot(self.lattice.inv_matrix)
            #position -= np.floor(position)
            self.position = position - np.floor(position)
            if update_mol:
                self.orientation = Orientation(np.eye(3))
                self.mol = mol
            else:
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
            import pickle
            with open('wrong.pkl', "wb") as f:
                pickle.dump([mol, self.mol], f)
                mol.to(filename='Wrong.xyz', fmt='xyz')
                self.mol.to(filename='Ref.xyz', fmt='xyz')
            raise ValueError("molecular connectivity changes! Exit")
        #todo check if connectivty changed
   
    def _find_gen_wyckoff_in_subgroup(self, groups=None):
        """
        Symmetry transformation
        group -> subgroup
        At the moment we only consider 
        for multiplicity 2: P-1, P21, P2, Pm and Pc
        to add: for multiplicity 4: P21/c, P212121
        Permutation is allowed
        """
        from pyxtal.symmetry import Wyckoff_position

        pos = self.position
        wp0 = self.wp
        pos0 = apply_ops(pos, wp0)

        if len(wp0) == 2:
            if self.diag: # P21/n -> Pn
                #print("----------P21n----------")
                wp1 = Wyckoff_position.from_group_and_index(7, 0)
                wp1.diagonalize_symops()
                axes = [[0,1,2],[2,1,0]]
                for ax in axes:
                    pos1 = apply_ops(pos[ax], wp1)
                    diff = (pos1[:, ax] - pos0)[1]
                    diff -= np.floor(diff)
                    if len(diff[diff==0]) >= 2:
                        #return wp1, ax, pos[ax] - 0.5*diff
                        return Wyckoff_position.from_group_and_index(7, 0), ax, pos[ax]-0.5*diff
            else:
                if groups is None:
                    groups = [4, 3, 6, 7, 2]
                if 15 < wp0.number < 71:
                    axes = [[0,1,2],[0,2,1],[1,0,2],[2,1,0]]
                elif wp0.number < 15:
                    axes = [[0,1,2],[2,1,0]]
                for group in groups:
                    wp1 = Wyckoff_position.from_group_and_index(group, 0)
                    for ax in axes:
                        pos1 = apply_ops(pos[ax], wp1)
                        diff = (pos1[:, ax] - pos0)[1]
                        diff -= np.floor(diff)
                        if len(diff[diff==0]) >= 2:
                            return wp1, ax, pos[ax]-0.5*diff

    def make_gen_wyckoff_site(self):

        if self.wp.index > 0:
            wp, ax, pos = self._find_gen_wyckoff_in_subgroup()
            ori = self.orientation.rotate_by_matrix(np.eye(3)[ax])
            lat = self.lattice.swap_axis(ids=ax)
            lat.ltype = Group(wp.number).lattice_type
            return mol_site(self.molecule, pos, ori, wp, lat, self.diag)
        else:
            print("This is already a general position")
            return self

    def _create_matrix(self):
        """
        Used for calculating distances in lattices with periodic boundary
        conditions. When multiplied with a set of points, generates additional
        points in cells adjacent to and diagonal to the original cell
        Returns:
            A numpy array of matrices which can be multiplied by a set of
            coordinates
        """
        matrix = []
        [a, b, c] = np.linalg.norm(self.lattice.matrix, axis=1)
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
                    if [i, j, k] != [0, 0, 0]:
                        matrix.append([i, j, k])
        #In case a,b,c are all greater than 20
        if len(matrix) == 0:
            matrix = [[1,0,0]]
        return np.array(matrix, dtype=float)

    def compute_distances(self):
        """
        compute if the atoms in the Wyckoff position are too close to each other
        or not. Does not check distances between atoms in the same molecule. Uses
        crystal.check_distance as the base code.

        Returns:
            minimum distances
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

        # Check periodic images
        m = self._create_matrix()
        coords_PBC = np.vstack([coords_mol + v for v in m])
        d = distance_matrix(coords_mol, coords_PBC, self.lattice.matrix, [0, 0, 0], True)
        if d < 0.9:
            return d
        else:
            min_ds.append(d)

        if self.wp.multiplicity > 1:
            # Check inter-atomic distances
            d = distance_matrix(coords_mol, coords, self.lattice.matrix, self.PBC, True)
            min_ds.append(d)
        return min(min_ds)

    def check_distances(self):
        """
        Checks if the atoms in the Wyckoff position are too close to each other
        or not. Does not check distances between atoms in the same molecule. Uses
        crystal.check_distance as the base code.

        Returns:
            True if the atoms are not too close together, False otherwise
        """
        m_length = len(self.symbols)
        coords, _ = self._get_coords_and_species()

        # Get coords of the generating molecule
        coords_mol = coords[:m_length]
        # Remove generating molecule's coords from large array
        coords = coords[m_length:]

        # Check periodic images
        m = self._create_matrix()
        coords_PBC = np.vstack([coords_mol + v for v in m])
        d = distance_matrix(coords_PBC, coords_mol, self.lattice.matrix, PBC=[0, 0, 0])
        # only check if small distance is detected
        if np.min(d) < np.max(self.tols_matrix):
            tols = np.min(d.reshape([len(m), m_length, m_length]), axis=0)
            if (tols < self.tols_matrix).any():
                return False

        if self.wp.multiplicity > 1:
            # Check inter-atomic distances
            d = distance_matrix(coords, coords_mol, self.lattice.matrix, PBC=self.PBC)
            if np.min(d) < np.max(self.tols_matrix):
                tols = np.min(
                    d.reshape([self.wp.multiplicity - 1, m_length, m_length]), axis=0
                )
                if (tols < self.tols_matrix).any():
                    return False

        return True



    def check_with_ms2(self, ms2, factor=1.0, tm=Tol_matrix(prototype="molecular")):
        """
        Checks whether or not the molecules of two mol sites overlap. Uses
        ellipsoid overlapping approximation to check. Takes PBC and lattice
        into consideration.
    
        Args:
            ms1: a mol_site object
            ms2: another mol_site object
            factor: the distance factor to pass to check_distances. (only for
                inter-atomic distance checking)
            tm: a Tol_matrix object (or prototype string) for distance checking
    
        Returns:
            False if the Wyckoff positions overlap. True otherwise
        """
        # Get coordinates for both mol_sites
        c1, _ = self.get_coords_and_species()
        c2, _ = ms2.get_coords_and_species()
    
        # Calculate which distance matrix is smaller/faster
        m_length1 = len(self.numbers)
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
            for i1, number1 in enumerate(self.numbers):
                for i2, number2 in enumerate(ms2.numbers):
                    tols[i1][i2] = tm.get_tol(number1, number2)
            tols = np.repeat(tols, ms2.wp.multiplicity, axis=1)
            d = distance_matrix(coords_mol, c2, self.lattice.matrix, PBC=self.PBC)
    
        # Case 2
        elif size1 > size2:
            coords_mol = c2[:m_length2]
            # Calculate tol matrix for species pairs
            tols = np.zeros((m_length2, m_length1))
            for i1, number1 in enumerate(ms2.numbers):
                for i2, number2 in enumerate(self.numbers):
                    tols[i1][i2] = tm.get_tol(number1, number2)
            tols = np.repeat(tols, self.wp.multiplicity, axis=1)
            d = distance_matrix(coords_mol, c1, self.lattice.matrix, PBC=self.PBC)
    
        # Check if distances are smaller than tolerances
        if (d < tols).any():
            return False
        return True


class atom_site:
    """
    Class for storing atomic Wyckoff positions with a single coordinate.

    Args:
        wp: a `Wyckoff_position <pyxtal.symmetry.Wyckoff_position.html> object 
        coordinate: a fractional 3-vector for the generating atom's coordinate
        specie: an Element, element name or symbol, or atomic number of the atom
        search: to search for the optimum position for special wyckoff site
    """

    def __init__(self, wp=None, coordinate=None, specie=1, diag=False, search=False):
        self.position = np.array(coordinate)
        self.specie = Element(specie).short_name
        self.diag = diag
        self.wp = wp
        if self.diag:
            self.wp.diagonalize_symops()
            #self.position = project_point(self.position, wp[0])

        self._get_dof()
        self.PBC = self.wp.PBC
        self.multiplicity = self.wp.multiplicity
        if search:
            self.search_position()
        self.update()

    def __str__(self):
        #number, dim = self.wp.number, self.wp.dim
        #self.site_symm = site_symm(self.wp.symmetry_m[0], number, dim=dim)

        if not hasattr(self, "site_symm"):
            self.site_symm = site_symm(
                self.wp.symmetry_m[0], self.wp.number, dim=self.wp.dim
            )
        s = "{:>2s} @ [{:7.4f} {:7.4f} {:7.4f}], ".format(self.specie, *self.position)
        s += "WP: {:2d}{:s}, ".format(self.wp.multiplicity, self.wp.letter)
        s += "Site symmetry: {:s}".format(self.site_symm)
        return s

    def __repr__(self):
        return str(self)

    def save_dict(self):
        dict0 = {"position": self.position,
                 "specie": self.specie,
                 "number": self.wp.number,
                 "dim": self.wp.dim,
                 "index": self.wp.index,
                 "PBC": self.wp.PBC,
                }
        return dict0

    def _get_dof(self):
        """
        get the number of dof for the given structures:
        """
        freedom = np.trace(self.wp.ops[0].rotation_matrix) > 0
        self.dof = len(freedom[freedom==True])


    @classmethod
    def load_dict(cls, dicts):
        """
        load the sites from a dictionary
        """
        g = dicts["number"]
        index = dicts["index"]
        dim = dicts["dim"]
        PBC = dicts["PBC"]
        position = dicts["position"]
        specie = dicts["specie"]
        wp = Wyckoff_position.from_group_and_index(g, index, dim, PBC)
        return cls(wp, position, specie)

    def perturbate(self, lattice, magnitude=0.1):
        """
        Random perturbation of the site
        
        Args:
            lattice: lattice vectors
            magnitude: the magnitude of displacement (default: 0.1 A)
        """
        dis = (np.random.sample(3) - 0.5).dot(lattice)
        dis /= np.linalg.norm(dis)
        dis *= magnitude
        pos = self.position + dis.dot(np.linalg.inv(lattice))
        self.update(pos)
 
    def search_position(self):
        """
        Sometimes, the initial posiition is not the proper generator
        Needs to find the proper generator
        """
        if self.wp.index > 0:
            pos = self.position
            coords = apply_ops(pos, Group(self.wp.number, self.wp.dim)[0])
            for coord in coords:
                ans = apply_ops(coord, [self.wp.ops[0]])[0]
                diff = coord - ans
                diff -= np.floor(diff)
                if np.sum(diff**2)<1e-4:
                    self.position = coord - np.floor(coord)
                    break
    def swap_axis(self, swap_id, shift=np.zeros(3)):
        """
        sometimes space groups like Pmm2 allows one to swap the a,b axes
        to get an alternative representation
        """
        self.position += shift
        self.position = self.position[swap_id]
        self.position -= np.floor(self.position)
        self.wp, _ = self.wp.swap_axis(swap_id)
        self.site_symm = site_symm(
            self.wp.symmetry_m[0], self.wp.number, dim=self.wp.dim
        )
        self.update()

    def shift_by_swap(self, swap_id):
        """
        check if a shift is needed during swap
        May occur for special WP in the I/A/B/C/F cases
        e.g., in space group 71 (Immm), the permutation
        4j(1/2,0,0.2) -> (0.2,0,1/2) -> 4f(0.7,1/2,0)
        it requires a shift of (0.5,0.5,0.5)
        """
        wp, shift = self.wp.swap_axis(swap_id)
        return shift

    def equivalent_set(self, tran, indices):
        """
        Transform the wp to another equivalent set.
        Needs to update both wp and positions

        Args:
            tran: affine matrix
            indices: the list of transformed wps 
        """
        self.position = SymmOp(tran).operate(self.position)
        self.position -= np.floor(self.position)
        self.wp = self.wp.equivalent_set(indices[self.wp.index]) #update the wp index
        self.site_symm = site_symm(
            self.wp.symmetry_m[0], self.wp.number, dim=self.wp.dim
        )
        self.update()


    def update(self, pos=None):
        """
        Used to generate coords from self.position
        """
        if pos is None:
            pos = self.position
        self.coords = apply_ops(pos, self.wp) 
        self.position = self.coords[0]

    def check_with_ws2(self, ws2, lattice, tm, same_group=True):
        """
        Given two Wyckoff sites, checks the inter-atomic distances between them.

        Args:
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
        if self.PBC != ws2.PBC:
            printx("Error: PBC values do not match between Wyckoff sites")
            return
        # Get tolerance
        tol = tm.get_tol(self.specie, ws2.specie)
        # Symmetry shortcut method: check only some atoms
        if same_group is True:
            # We can either check one atom in WS1 against all WS2, or vice-versa
            # Check which option is faster
            if self.multiplicity > ws2.multiplicity:
                coords1 = [self.coords[0]]
                coords2 = ws2.coords
            else:
                coords1 = [ws2.coords[0]]
                coords2 = self.coords
            # Calculate distances
            dm = distance_matrix(coords1, coords2, lattice, PBC=self.PBC)
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

def WP_merge(pt, lattice, wp, tol, orientations=None):
    """
    Given a list of fractional coordinates, merges them within a given
    tolerance, and checks if the merged coordinates satisfy a Wyckoff
    position. Used for merging general Wyckoff positions into special Wyckoff
    positions within the random_crystal (and its derivative) classes.

    Args:
        pt: the originl point (3-vector)
        lattice: a 3x3 matrix representing the unit cell
        wp: a `Wyckoff_position <pyxtal.symmetry.Wyckoff_position.html> object after merge
        tol: the cutoff distance for merging coordinates
        orientations: the valid orientations for a given molecule. Obtained
            from get_sg_orientations, which is called within molecular_crystal

    Returns:
        pt: 3-vector after merge
        wp: a `pyxtal.symmetry.Wyckoff_position` object, If no matching WP, returns False. 
        valid_ori: the valid orientations after merge

    """
    index = wp.index
    PBC = wp.PBC
    group = Group(wp.number, wp.dim)
    pt = project_point(pt, wp[0], lattice, PBC)
    coor = apply_ops(pt, wp)
    if orientations is None:
        valid_ori = None
    else:
        j, k = jk_from_i(index, orientations)
        valid_ori = orientations[j][k]
    
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

        # for molecular crystal, one more check
        if not check_images([coor[0]], [6], lattice, PBC=PBC, tol=tol):
            passed_distance_check = False

        if not passed_distance_check:
            mult1 = group[index].multiplicity
            # Find possible wp's to merge into
            possible = []
            for i, wp0 in enumerate(group):
                mult2 = wp0.multiplicity
                # Check that a valid orientation exists
                if orientations is not None:
                    j, k = jk_from_i(i, orientations)
                    if orientations[j][k] == []:
                        continue
                    else:
                        valid_ori = orientations[j][k]
                # factor = mult2 / mult1

                if (mult2 < mult1) and (mult1 % mult2 == 0):
                    possible.append(i)
            if possible == []:
                return None, False, valid_ori
            # Calculate minimum separation for each WP
            distances = []
            for i in possible:
                wp = group[i]
                projected_point = project_point(pt.copy(), wp[0], lattice=lattice, PBC=PBC)
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
            return pt, wp, valid_ori


