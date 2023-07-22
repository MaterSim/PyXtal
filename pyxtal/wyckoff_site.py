"""
Module for handling Wyckoff sites for both atom and molecule
"""

# Standard Libraries
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cdist
from copy import deepcopy

# External Libraries
from pymatgen.core import Molecule

# PyXtal imports
from pyxtal.tolerance import Tol_matrix
from pyxtal.operations import (
    check_images,
    distance_matrix,
    filtered_coords,
    create_matrix,
    SymmOp,
)
from pyxtal.symmetry import Group, Wyckoff_position
from pyxtal.database.element import Element
from pyxtal.constants import rad, deg
from pyxtal.lattice import Lattice

class atom_site:
    """
    Class for storing atomic Wyckoff positions with a single coordinate.

    Args:
        wp: a `Wyckoff_position <pyxtal.symmetry.Wyckoff_position.html> object
        coordinate: a fractional 3-vector for the generating atom's coordinate
        specie: an Element, element name or symbol, or atomic number of the atom
        search: to search for the optimum position for special wyckoff site
    """

    def __init__(self, wp=None, coordinate=None, specie=1, search=False):
        self.position = np.array(coordinate)
        self.specie = Element(specie).short_name
        self.wp = wp
        self.coordination = None

        self._get_dof()
        self.PBC = self.wp.PBC
        self.multiplicity = self.wp.multiplicity
        if search:
            self.search_position()
        self.update()

    def __str__(self):
        if not hasattr(self.wp, "site_symm"): self.wp.get_site_symmetry()
        s = "{:>2s} @ [{:7.4f} {:7.4f} {:7.4f}], ".format(self.specie, *self.position)
        s += "WP [{:}] ".format(self.wp.get_label())
        if self.coordination is not None:
            s += " CN [{:2d}] ".format(self.coordination)
        s += "Site [{:}]".format(self.wp.site_symm.replace(" ",""))

        return s

    def __repr__(self):
        return str(self)

    def copy(self):
        """
        Simply copy the structure
        """
        return deepcopy(self)

    def save_dict(self):
        dict0 = {"position": self.position,
                 "specie": self.specie,
                 "wp": self.wp.save_dict(),
                }
        return dict0

    def _get_dof(self):
        """
        get the number of dof for the given structures:
        """
        self.dof = self.wp.get_dof()
        #freedom = np.trace(self.wp.ops[0].rotation_matrix) > 0
        #self.dof = len(freedom[freedom==True])
        #self.dof = len(freedom[freedom==True])


    @classmethod
    def load_dict(cls, dicts):
        """
        load the sites from a dictionary
        """
        position = dicts["position"]
        specie = dicts["specie"]
        wp = Wyckoff_position.load_dict(dicts['wp'])
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
            wp0 = Group(self.wp.number, self.wp.dim)[0]
            pos = self.position
            coords = wp0.apply_ops(pos)
            for coord in coords:
                ans = self.wp.ops[0].operate(coord)
                diff = coord - ans
                diff -= np.floor(diff)
                if np.sum(diff**2)<1e-4:
                    self.position = coord - np.floor(coord)
                    break

    def encode(self):
        """
        transform dict to 1D vector
        [specie, wp.index, free x, y, z]
        """
        xyz = self.wp.get_free_xyzs(self.position)
        #print(self.wp.ops[0].rotation_matrix, self.wp.get_frozen_axis(), self.wp.get_dof())
        #print([self.specie, self.wp.index] + list(xyz))
        return [self.specie, self.wp.index] + list(xyz)
        

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
            self.wp.symmetry[0], self.wp.number, dim=self.wp.dim
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

    def update(self, pos=None, reset_wp=False):
        """
        Used to generate coords from self.position
        """
        if pos is None:
            pos = self.position
        if reset_wp:
            self.wp.ops = Group(self.wp.number)[self.wp.index].ops
        self.coords = self.wp.apply_ops(pos)
        self.position = self.coords[0]

    def get_translations(self, pos, axis):
        """
        return the displacement towards the reference positions

        Args:
            pos: reference position (1*3 vector)
            lattice: 3*3 matrix
            translation:
            axis:
        """
        #diffs0 = pos - self.coords
        diffs0 = self.wp.apply_ops(pos) - self.position
        diffs = diffs0.copy()
        diffs -= np.round(diffs)
        diffs[:, axis] = 0
        translations = diffs0 - diffs
        return translations

    def get_disp(self, pos, lattice, translation):
        """
        return the displacement towards the reference positions

        Args:
            pos: reference position (1*3 vector)
            lattice: 3*3 matrix
            translation:
        """
        coords = self.wp.apply_ops(pos)
        diffs = coords - (self.position + translation)
        #coords = self.wp.apply_ops(self.position + translation)
        #diffs = pos - coords

        diffs -= np.round(diffs)
        dists = np.linalg.norm(diffs.dot(lattice), axis=1)
        id = np.argmin(dists)

        #print("++++", id, dists[id], id, diffs[id], translation) #; import sys; sys.exit()
        return diffs[id], dists[id]

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
            raise ValueError("PBC values do not match between Wyckoff sites")
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

    def substitute_with_single(self, ele):
        """
        chemical substitution with another element

        Args:
            ele (str): e.g. 'Zn'
        """
        self.specie = ele
        return self

    def substitute_with_linear(self, eles, direction, lattice):
        """
        chemical substitution with another linear building block, e.g. CN

        Args:
            eles (str): e.g., ['C', 'N']
            neighbors: two nearest neiboring atom xyz
        """
        bond_length = 1.4
        #direction = neighbors[1] - neighbors[0]
        #direction = np.dot(direction, lattice)
        #direction /= np.linalg.norm(direction)

        # Get the fraction coordinates
        shift = np.dot(bond_length/2 * direction, np.linalg.inv(lattice))
        site1 = self.copy()
        site2 = self.copy()
        site1.specie = eles[0]
        site2.specie = eles[1]
        site1.update(site1.position + shift)
        site2.update(site2.position - shift)
        return site1, site2
        

class mol_site:
    """
    Class for storing molecular Wyckoff positions and orientations within
    the molecular_crystal class. Each mol_site object represenents an
    entire Wyckoff position, not necessarily a single molecule.
    This is the molecular version of Wyckoff_site

    Args:
        mol: a `pyxtal_molecule <pyxtal.molecule.pyxtal_molecule.html>`_ object
        position: the 3-vector representing the generating molecule's position
        orientation: an `Orientation <pyxtal.molecule.Oreintation.html>`_ object
        wp: a `Wyckoff_position <pyxtal.symmetry.Wyckoff_position.html>`_ object
        lattice: a `Lattice <pyxtal.lattice.Lattice>`_ object
        stype: integer number to specify the type of molecule
    """

    def __init__(self, mol, position, orientation, wp, lattice=None, stype=0):
        # describe the molecule
        self.molecule = mol
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
        self.type = stype

    def __str__(self):
        if not hasattr(self.wp, "site_symm"): 
            self.wp.get_site_symmetry()

        self.angles = self.orientation.r.as_euler('zxy', degrees=True)
        formula = self.mol.formula.replace(" ","")
        s = "{:12s} @ [{:7.4f} {:7.4f} {:7.4f}]  ".format(formula, *self.position)
        s += "WP [{:s}] ".format(self.wp.get_label())
        s += "Site [{:}]".format(self.wp.site_symm.replace(" ",""))
        if len(self.molecule.mol) > 1:
            s += " Euler [{:6.1f} {:6.1f} {:6.1f}]".format(*self.angles)

        return s

    def __repr__(self):
        return str(self)

    def _get_dof(self):
        """
        get the number of dof for the given structures:
        """
        freedom = np.trace(self.wp.ops[0].rotation_matrix) > 0
        self.dof = len(freedom[freedom==True])

    def save_dict(self):
        dict0 = {"position": self.position,
                 "wp": self.wp.save_dict(),
                 "molecule": self.molecule.save_str(),
                 "orientation": self.orientation.save_dict(),
                 "lattice": self.lattice.matrix,
                 "lattice_type": self.lattice.ltype,
                 "stype": self.type,
                }
        if self.molecule.torsionlist is not None:
            xyz, _ = self._get_coords_and_species(absolute=True, first=True)
            dict0['rotors'] = self.molecule.get_torsion_angles(xyz)

        return dict0

    @classmethod
    def load_dict(cls, dicts):
        """
        load the sites from a dictionary
        """
        from pyxtal.molecule import pyxtal_molecule, Orientation

        mol = pyxtal_molecule.load_str(dicts["molecule"])
        position = dicts["position"]
        orientation = Orientation.load_dict(dicts['orientation'])
        wp = Wyckoff_position.load_dict(dicts['wp'])
        lattice = Lattice.from_matrix(dicts["lattice"], ltype=dicts["lattice_type"])
        stype = dicts["stype"]
        return cls(mol, position, orientation, wp, lattice, stype)

    def encode(self):
        """
        transform dict to 1D vector
        [x, y, z, or1, or2, or3, rotor1, rotor2, .etc]
        """
        if len(self.molecule.mol)>1:
            xyz, _ = self._get_coords_and_species(absolute=True, first=True)
            #if len(xyz)==3: print("encode: \n", self.molecule.mol.cart_coords)
            rotor = self.molecule.get_torsion_angles(xyz)
            ori, _, reflect = self.molecule.get_orientation(xyz)
            return list(self.position) + list(ori) + rotor + [reflect]
        else:
            return list(self.position) + [0]

    def to_1D_dicts(self):
        """
        save the wp in 1D representation
        """
        xyz, _ = self._get_coords_and_species(absolute=True, first=True)
        dict0 = {"smile": self.molecule.smile}
        dict0["rotor"] = self.molecule.get_torsion_angles(xyz)
        dict0["orientation"], dict0["rmsd"], dict0["reflect"] = self.molecule.get_orientation(xyz)
        angs = dict0["rotor"]
        #rdkit_mol = self.molecule.rdkit_mol(self.molecule.smile)
        #conf0 = rdkit_mol.GetConformer(0)
        #print(self.molecule.set_torsion_angles(conf0, angs))
        #print("save matrix"); print(self.orientation.r.as_matrix())
        #print("save angle"); print(self.orientation.r.as_euler('zxy', degrees=True))
        #print("angle"); print(dict0["orientation"])
        dict0["center"] = self.position - np.floor(self.position) #self.molecule.get_center(xyz)
        dict0["number"] = self.wp.number
        dict0["index"] = self.wp.index
        dict0["PBC"] = self.wp.PBC
        dict0["dim"] = self.wp.dim
        dict0["lattice"] = self.lattice.matrix
        dict0["lattice_type"] = self.lattice.ltype

        return dict0

    @classmethod
    def from_1D_dicts(cls, dicts):
        from pyxtal.molecule import pyxtal_molecule, Orientation

        mol = pyxtal_molecule(mol=dicts['smile']+'.smi', fix=True)
        if len(mol.mol) > 1:
            if len(dicts['smile']) > 1:
                conf = mol.rdkit_mol().GetConformer(0)
                if dicts['reflect']:
                    mol.set_torsion_angles(conf, dicts["rotor"], False)
                xyz = mol.set_torsion_angles(conf, dicts["rotor"], dicts['reflect'])
            else:
                # for H2O, use the standard one
                xyz = np.array([[-0.00111384,  0.36313718,  0.        ],
                                [-0.82498189, -0.18196256,  0.        ],
                                [ 0.82609573, -0.18117463,  0.        ]])

            mol.reset_positions(xyz)
            matrix = R.from_euler('zxy', dicts["orientation"], degrees=True).as_matrix()
            orientation = Orientation(matrix)
        else:
            orientation = Orientation(np.eye(3))

        g = dicts["hn"]
        index = dicts["index"]
        dim = dicts["dim"]
        wp = Wyckoff_position.from_group_and_index(g, index, dim, dicts["PBC"])
        lattice = Lattice.from_matrix(dicts["lattice"], ltype=dicts["lattice_type"])
        position = dicts["center"] #np.dot(dicts["center"], lattice.inv_matrix)
        position, wp, _ = wp.merge(position, lattice.matrix, 0.01)

        return cls(mol, position, orientation, wp, lattice)

    def show(self, id=None, **kwargs):
        """
        display WP on the notebook
        """
        from pyxtal.viz import display_molecular_site
        return display_molecular_site(self, id, **kwargs)

    def _get_coords_and_species(self, absolute=False, PBC=False, first=False, unitcell=False):
        """
        Used to generate coords and species for get_coords_and_species

        Args:
            absolute: return absolute or relative coordinates
            PBC: whether or not to add coordinates in neighboring unit cells,
            first: whether or not to extract the information from only the first site
            unitcell: whether or not to move the molecular center to the unit cell

        Returns:
            atomic coords: a numpy array of atomic coordinates in the site
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
            #op2_m = self.wp.generators_m[point_index]
            op2_m = self.wp.get_euclidean_generator(self.lattice.matrix, point_index)
            rot = op2_m.affine_matrix[:3, :3].T
            # NOTE=====the euclidean_generator has wrong translation vectors,
            # but we don't care. This needs to be fixed later

            #if self.diag and self.wp.index > 0:
            #    tau = op2.translation_vector
            #else:
            #    tau = op2_m.translation_vector
            tmp = np.dot(coord0, rot) #+ tau

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
            absolute: return absolute or relative coordinates
            PBC: whether or not to add coordinates in neighboring unit cells
            unitcell: whether or not to move the molecule center to the unit cell

        Returns:
            coords: a np array of 3-vectors.
            species: a list of atomic symbols, e.g. ['H', 'H', 'O', 'H', 'H', 'O']
        """
        return self._get_coords_and_species(absolute, PBC, unitcell=unitcell)

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
        self.position = self.wp.project(position)

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
            ax = self.molecule.get_principle_axes(xyz).T[ax_id]

        q = R.from_rotvec(ax*rad*angle)
        o = q*p
        self.orientation.r = o
        self.orientation.matrix = o.as_matrix()

    #def is_compatible_symmetry(self, tol=0.3):
    #    """
    #    Check if the molecular symmetry matches the site symmetry
    #    """
    #    mol = self.molecule.mol
    #    if len(mol) == 1 or self.wp.index==0:
    #        return True
    #    else:
    #        pga = PointGroupAnalyzer(mol, tol)
    #        for op in self.wp.symmetry_m[0]:
    #            if not pga.is_valid_op(op):
    #                return False
    #        return True

    def is_valid_orientation(self):
        pass

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
        if not hasattr(self.wp, "generators"): self.wp.set_generators()

        if id <= len(self.wp.generators):
            #op = self.wp.generators[id]
            op = self.wp.get_euclidean_generator(self.lattice.matrix, id)
            center_relative = op.operate(self.position)
            center_relative -= np.floor(center_relative)
            center_absolute = np.dot(center_relative, self.lattice.matrix)

            # Rotate the molecule (Euclidean metric)
            #op_m = self.wp.generators_m[id]
            #rot = op_m.affine_matrix[0:3][:, 0:3].T
            #tau = op_m.affine_matrix[0:3][:, 3]
            op0 = self.wp.get_euclidean_generator(self.lattice.matrix, id)
            rot = op0.rotation_matrix.T
            tmp = np.dot(coord0, rot)
            # Add absolute center to molecule
            tmp += center_absolute
            return Molecule(self.symbols, tmp)
        else:
            raise ValueError("id is greater than the number of molecules")

    def show_molecule_in_box(self, id=0):
        """
        display molecule with box

        Args:
            id (int): molecular id
        """
        from pyxtal.viz import display_molecule
        mol = self.get_mol_object(id)
        cell, vertices, center = self.molecule.get_box_coordinates(mol.cart_coords)
        return display_molecule(mol, center, cell)

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
        # todo check if connectivty changed

    def _create_matrix(self, center=False, ignore=False):
        """
        Used for calculating distances in lattices with periodic boundary
        conditions. When multiplied with a set of points, generates additional
        points in cells adjacent to and diagonal to the original cell

        Returns:
            A numpy array of matrices which can be multiplied by a set of
            coordinates
        """
        abc = [self.lattice.a, self.lattice.b, self.lattice.c]

        # QZ: This should be based on the occupation of current molecule
        if hasattr(self, 'ijk_lists'):
            ijk_lists = self.ijk_lists
        else:
            ijk_lists = []
            for id in range(3):
                if self.PBC[id]:
                    if not ignore and abc[id] > 25 and self.radius < 10:
                        ijk_lists.append([0])
                    elif abc[id] < 7.0:
                        ijk_lists.append([-3, -2, -1, 0, 1, 2, 3])
                    else:
                        ijk_lists.append([-1, 0, 1])
                else:
                    ijk_lists.append([0])

        if center:
            matrix = [[0,0,0]]
        else:
            matrix = []
        for i in ijk_lists[0]:
            for j in ijk_lists[1]:
                for k in ijk_lists[2]:
                    if [i, j, k] != [0, 0, 0]:
                        matrix.append([i, j, k])
        # In case a,b,c are all greater than 20
        if len(matrix) == 0:
            matrix = [[1,0,0]]
        return np.array(matrix, dtype=float)


    def get_distances(self, coord1, coord2, m2=None, center=True, ignore=False):
        """
        Compute the distance matrix between the center molecule (m1 length) and
        neighbors (m2 length) within the PBC consideration (pbc)

        Args:
            coord1: fractional coordinates of the center molecule
            coord2: fractional coordinates of the reference neighbors
            m2: the length of reference molecule
            center: whether or not consider the self image for coord2

        Returns:
            distance matrix: [m1*m2*pbc, m1, m2]
            coord2 under PBC: [pbc, m2, 3]
        """
        m1 = len(coord1)
        if m2 is None: m2 = m1
        N2 = int(len(coord2)/m2)

        # peridoic images
        m = self._create_matrix(center, ignore) #PBC matrix
        coord2 = np.vstack([coord2 + v for v in m])

        # absolute xyz
        coord1 = np.dot(coord1, self.lattice.matrix)
        coord2 = np.dot(coord2, self.lattice.matrix)

        d = cdist(coord1, coord2)
        d = d.T.reshape([len(m)*N2, m1, m2])
        return d, coord2.reshape([len(m)*N2, m2, 3])


    def get_dists_auto(self, ignore=False):
        """
        Compute the distances between the periodic images

        Returns:
            a distance matrix (M, N, N)
            list of molecular xyz (M, N, 3)
        """
        m_length = len(self.numbers)
        coord1, _ = self._get_coords_and_species(first=True, unitcell=True)

        return self.get_distances(coord1, coord1, center=False, ignore=ignore)

    def get_dists_WP(self, ignore=False, id=None):
        """
        Compute the distances within the WP sites

        Returns:
            a distance matrix (M, N, N)
            list of molecular xyz (M, N, 3)
        """
        m_length = len(self.symbols)
        coords, _ = self._get_coords_and_species(unitcell=True)
        coord1 = coords[:m_length] #1st molecular coords
        if id is None:
            coord2 = coords[m_length:] #rest molecular coords
        else:
            coord2 = coords[m_length*(id):m_length*(id+1)] #rest molecular coords

        return self.get_distances(coord1, coord2, ignore=ignore)

    def get_min_dist(self):
        """
        Compute the minimum interatomic distance within the WP.

        Returns:
            minimum distance
        """
        # Self image
        ds, _ = self.get_dists_auto()
        min_dist = np.min(ds)

        if min_dist < 0.9:
            # terminate earlier
            return min_dist
        else:
            # Other molecules
            if self.wp.multiplicity > 1:
                ds, _ = self.get_dists_WP()
                if min_dist > np.min(ds):
                    min_dist = np.min(ds)
            return min_dist

    def short_dist(self):
        """
        Check if the atoms are too close within the WP.

        Returns:
            True or False
        """
        m_length = len(self.numbers)
        tols_matrix = self.tols_matrix
        # Check periodic images
        d, _ = self.get_dists_auto()
        if np.min(d) < np.max(tols_matrix):
            tols = np.min(d, axis=0)
            if (tols < tols_matrix).any():
                return False

        if self.wp.multiplicity > 1:
            d, _ = self.get_dists_WP()
            if np.min(d) < np.max(tols_matrix):
                tols = np.min(d, axis=0) #N*N matrix
                if (tols < tols_matrix).any():
                    return False

        return True

    def short_dist_with_wp2(self, wp2, tm=Tol_matrix(prototype="molecular")):
        """
        Check whether or not the molecules of two wp sites overlap. Uses
        ellipsoid overlapping approximation to check.

        Args:
            wp2: the 2nd wp sites
            tm: a Tol_matrix object (or prototype string) for distance checking
        Returns:
            True or False
        """

        # Get coordinates for both mol_sites
        c1, _ = self.get_coords_and_species()
        c2, _ = wp2.get_coords_and_species()
        m_length1 = len(self.numbers)
        m_length2 = len(wp2.numbers)

        # choose wp with bigger molecule as the center
        if len(c2) <= len(c1):
            coord1 = c1[:m_length1]
            coord2 = c2 #rest molecular coords
            tols_matrix = self.molecule.get_tols_matrix(wp2.molecule, tm)
            m2 = m_length2
        else:
            coord1 = c2[:m_length2]
            coord2 = c1
            tols_matrix = wp2.molecule.get_tols_matrix(self.molecule, tm)
            m2 = m_length1

        # compute the distance matrix
        d, _ = self.get_distances(coord1-np.floor(coord1), coord2-np.floor(coord2), m2)
        #print("short dist", len(c1), len(c2), d.min())
        
        if np.min(d) < np.max(tols_matrix):
            tols = np.min(d, axis=0)
            if (tols < tols_matrix).any():
                return False
        return True

    def get_neighbors_auto(self, factor=1.1, max_d=4.0, ignore_E=True, detail=False, etol=-5e-2):
        """
        Find the neigboring molecules

        Args:
            factor: volume factor
            max_d: maximum intermolecular distance
            ignore_E: 
            detail: show detailed energies

        Returns
            min_ds: list of shortest distances
            neighs: list of neighboring molecular xyzs
        """
        coord1, _ = self._get_coords_and_species(first=True, unitcell=True)
        tm = Tol_matrix(prototype="vdW", factor=factor)
        m_length = len(self.numbers)
        tols_matrix = self.molecule.get_tols_matrix(tm=tm)
        coef_matrix = None
        if not ignore_E:
            coef_matrix = self.molecule.get_coefs_matrix()
            if coef_matrix is not None:
                A = coef_matrix[:,:,0]
                B = coef_matrix[:,:,1]
                C = coef_matrix[:,:,2]

        min_ds = []
        neighs = []
        Ps = []
        engs = []
        pairs = []
        dists = []

        # Check periodic images
        d, coord2 = self.get_dists_auto(ignore=True)
        for i in range(d.shape[0]):
            if np.min(d[i]) < max_d and (d[i] < tols_matrix).any():
                if coef_matrix is not None:
                    if detail:
                        eng = A*np.exp(-B*d[i])-C/(d[i]**6)
                        ids = np.where(eng < etol)
                        for id in zip(*ids):
                            tmp1, tmp2 = coord1[id[0]], coord2[i][id[1]]
                            pairs.append((tmp1+tmp2)/2)
                            engs.append(eng[id])
                            dists.append(d[i][id])
                        #eng = eng.sum()
                    else:
                        eng0 = A*np.exp(-B*d[i])-C/(d[i]**6)
                        engs.append(eng0.sum())
                else:
                    engs.append(None)

                tmp = d[i]/tols_matrix
                _d = tmp[tmp < 1.0]
                id = np.argmin(tmp.flatten())
                d_min = d[i].flatten()[id]
                min_ds.append(min(_d)*factor)
                neighs.append(coord2[i])
                Ps.append(0)

        if self.wp.multiplicity > 1:
            for idx in range(1, self.wp.multiplicity):
                if self.wp.is_pure_translation(idx):
                    P = 0
                else:
                    P = 1
                d, coord2 = self.get_dists_WP(ignore=True, id=idx)
                for i in range(d.shape[0]):
                    if np.min(d[i])<max_d and (d[i] < tols_matrix).any():
                        if coef_matrix is not None:
                            if detail:
                                eng = A*np.exp(-B*d[i])-C/(d[i]**6)
                                ids = np.where(eng < etol)
                                for id in zip(*ids):
                                    tmp1, tmp2 = coord1[id[0]], coord2[i][id[1]]
                                    pairs.append((tmp1+tmp2)/2)
                                    engs.append(eng[id])
                                    dists.append(d[i][id])
                            else:
                                eng0 = A*np.exp(-B*d[i])-C/(d[i]**6)
                                engs.append(eng0.sum())
                        else:
                            engs.append(None)
                        tmp = d[i]/tols_matrix
                        _d = tmp[tmp < 1]
                        id = np.argmin(tmp.flatten())
                        d_min = d[i].flatten()[id]
                        min_ds.append(min(_d)*factor)
                        neighs.append(coord2[i])
                        Ps.append(P)
        if detail:
            return engs, pairs, dists
        else:
            return min_ds, neighs, Ps, engs

    def get_neighbors_wp2(self, wp2, factor=1.1, max_d=4.0, ignore_E=True, detail=False, etol=-5e-2):
        """
        Find the neigboring molecules from a 2nd wp site

        Returns
            min_ds: list of shortest distances
            neighs: list of neighboring molecular xyzs
        """

        tm=Tol_matrix(prototype="vdW", factor=factor)
        m_length1 = len(self.numbers)
        m_length2 = len(wp2.numbers)

        # Get coordinates for both mol_sites
        c1, _ = self.get_coords_and_species()
        c2, _ = wp2.get_coords_and_species()

        coord1 = c1[:m_length1]
        coord2 = c2 #rest molecular coords
        tols_matrix = self.molecule.get_tols_matrix(wp2.molecule, tm)
        coef_matrix = None
        if not ignore_E:
            coef_matrix = self.molecule.get_coefs_matrix(wp2.molecule)
            if coef_matrix is not None:
                A = coef_matrix[:,:,0]
                B = coef_matrix[:,:,1]
                C = coef_matrix[:,:,2]

        # compute the distance matrix
        d, coord2 = self.get_distances(coord1, coord2, m_length2, ignore=True)
        min_ds = []
        neighs = []
        engs = []
        dists = []
        pairs = []

        for i in range(d.shape[0]):
            if np.min(d[i])<max_d and (d[i] < tols_matrix).any():
                if coef_matrix is not None:
                    if detail:
                        eng = A*np.exp(-B*d[i])-C/(d[i]**6)
                        ids = np.where(eng < etol)
                        for id in zip(*ids):
                            tmp1, tmp2 = coord1[id[0]], coord2[i][id[1]]
                            pairs.append((tmp1+tmp2)/2)
                            engs.append(eng[id])
                            dists.append(d[i][id])
                    else:
                        eng0 = A*np.exp(-B*d[i])-C/(d[i]**6)
                        engs.append(eng0.sum())
                else:
                    engs.append(None)
                tmp = d[i]/tols_matrix
                _d = tmp[tmp < 1]
                id = np.argmin(tmp.flatten())
                d_min = d[i].flatten()[id]
                min_ds.append(min(_d)*factor)
                neighs.append(coord2[i])

        if detail:
            return engs, pairs, dists
        else:
            return min_ds, neighs, engs


    def get_ijk_lists(self, value=None):
        """
        Get the occupatation in the unit cell for the generating molecule
        This can be used to estimate the supercell size for finding neighbors

        Returns
            PBC
        """

        # Get coordinates for both mol_sites
        if value is None:
            ijk_lists = []
            c1, _ = self.get_coords_and_species(absolute=False)
            for id in range(3):
                if self.PBC[id]:
                    m1 = np.min(c1[:,id])
                    m2 = np.max(c1[:,id])
                    max_id = int(np.ceil(2*m2-m1))
                    min_id = int(np.floor(-m2))
                    ijk_lists.append(list(range(min_id-1, max_id+1)))
                else:
                    ijk_lists.append([0])
            self.ijk_lists = ijk_lists
        else:
            self.ijk_lists = value
