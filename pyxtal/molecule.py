"""
Module for handling molecules. Uses the pymatgen.core.structure.Molecule
class as a base. Has a function for reorienting molecules
(reoriented_molecule), and for calculating valid orientations within a Wyckoff
position based on symmetry (orientation_in_wyckoff_position).
The orientation class can be used to identify
degrees of freedom for molecules in Wyckoff positions with certain symmetry
constraints.

"""
# Imports
import os
from copy import deepcopy
from operator import itemgetter
from random import choice
import numpy as np
from scipy.spatial.transform import Rotation
import networkx as nx
# ------------------------------
# External Libraries
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer, generate_full_symmops
from pymatgen.core.bonds import CovalentBond

# PyXtal imports
from pyxtal.msg import printx
from pyxtal.tolerance import Tol_matrix
from pyxtal.database.element import Element
from pyxtal.operations import SymmOp, OperationAnalyzer, rotate_vector, angle
from pyxtal.database.collection import Collection

# Define functions
# ------------------------------
molecule_collection = Collection("molecules")

def cleaner(list_to_clean):
    """
    Remove duplicate torsion definion from a list of atom ind. tuples.
    """
    for_remove = []
    for x in reversed(range(len(list_to_clean))):
        for y in reversed(range(x)):
            ix1, ix2 = itemgetter(1)(list_to_clean[x]), itemgetter(2)(list_to_clean[x])
            iy1, iy2 = itemgetter(1)(list_to_clean[y]), itemgetter(2)(list_to_clean[y])
            if (ix1 == iy1 and ix2 == iy2) or (ix1 == iy2 and ix2 == iy1):
                for_remove.append(y)
    clean_list = [v for i, v in enumerate(list_to_clean)
                  if i not in set(for_remove)]
    return clean_list

def find_id_from_smile(smile):
    """
    Find the positions of rotatable bonds in the molecule.
    """
    from rdkit import Chem

    smarts_torsion="[*]~[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]~[*]"
    mol = Chem.MolFromSmiles(smile)
    pattern_tor = Chem.MolFromSmarts(smarts_torsion)
    torsion = list(mol.GetSubstructMatches(pattern_tor))
    return cleaner(torsion)

def dihedral(p):
    """
    dihedral from https://stackoverflow.com/questions/20305272
    """
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

class pyxtal_molecule:
    """
    Extended molecule class based on pymatgen.core.structure.Molecule
    The added features include:
    0, parse the input
    1, estimate volume/tolerance/radii
    2, find and store symmetry 
    3, get the principle axis 
    4, re-align the molecule 

    The molecule is always centered at (0, 0, 0).

    If the smile format is used, the center is defined as in
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolTransforms.html

    Otherwise, the center is just the mean of atomic positions

    Args:
        mol: a string to reprent the molecule
        tm: tolerance matrix
    """

    def __init__(self, mol=None, symmetrize=True, tm=Tol_matrix(prototype="molecular")):
        mo = None
        self.smile = None
        self.torsionlist = None
        if type(mol) == str:
            # Parse molecules: either file or molecule name
            tmp = mol.split(".")
            self.name = tmp[0]
            if len(tmp) > 1:
                # Load the molecule from the given file
                if tmp[-1] in ["xyz", "gjf", "g03", "json"]:
                    if os.path.exists(mol):
                        mo = Molecule.from_file(mol)
                    else:
                        raise NameError("{:s} is not a valid path".format(mol))
                elif tmp[-1] == 'smi':
                    self.smile = tmp[0]
                    symbols, xyz, self.torsionlist = self.rdkit_mol_init(tmp[0])
                    mo = Molecule(symbols, xyz)
                    symmetrize = False
                else:
                    raise NameError("{:s} is not a supported format".format(tmp[-1]))
            else:
                # print('\nLoad the molecule {:s} from collections'.format(mol))
                mo = molecule_collection[mol]
        elif hasattr(mol, "sites"):  # pymatgen molecule
            self.name = str(mol.formula)
            mo = mol

        if mo is None:
            msg = "Could not create molecules from given input: {:s}".format(mol)
            raise NameError(msg)

        self.props = mo.site_properties

        if len(mo) > 1:
            if symmetrize:
                pga = PointGroupAnalyzer(mo)
                mo = pga.symmetrize_molecule()["sym_mol"]
            mo = self.add_site_props(mo)

        self.mol = mo
        self.tm = tm
        self.get_box()
        self.volume = self.box.volume
        self.get_radius()
        self.get_symbols()
        self.get_tols_matrix()
        xyz = self.mol.cart_coords
        self.reset_positions(xyz-self.get_center(xyz))


    def __str__(self):
        return '[' + self.name + ']'

    def save_dict(self):
        """
        save the object as a dictionary
        """
        return self.mol.as_dict()

    def copy(self):
        """
        simply copy the structure
        """
        return deepcopy(self)

    @classmethod
    def load_dict(cls, dicts):
        """
        load the molecule from a dictionary
        """
        mol = Molecule.from_dict(dicts)
        return cls(mol)

    def swap_axis(self, ax):
        """
        swap the molecular axis
        """
        coords = self.mol.cart_coords[:, ax]
        mo = Molecule(self.symbols, coords)
        mo = self.add_site_props(mo)

        return pyxtal_molecule(mo, self.tm)

    def add_site_props(self, mo):
        """
        add site properties
        """
        if len(self.props) > 0:
            for key in self.props.keys():
                mo.add_site_property(key, self.props[key])
        return mo

    def get_box(self):
        """
        Given a molecule, find a minimum orthorhombic box containing it.
        Size is calculated using min and max x, y, and z values,
        plus the padding defined by the vdw radius
        For best results, call oriented_molecule first.

        Args:
            mol: a pymatgen Molecule object. Should be oriented along its principle axes.

        Returns:
            a Box object
        """
        mol, P = reoriented_molecule(self.mol)
        minx, miny, minz, maxx, maxy, maxz = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for p in mol:
            x, y, z = p.coords
            r = Element(p.species_string).vdw_radius
            if x - r < minx:
                minx = x - r
            if y - r < miny:
                miny = y - r
            if z - r < minz:
                minz = z - r
            if x + r > maxx:
                maxx = x + r
            if y + r > maxy:
                maxy = y + r
            if z + r > maxz:
                maxz = z + r
        self.box = Box(minx, maxx, miny, maxy, minz, maxz)
        self.axes = P

    def get_radius(self):
        """
        get the radius of a molecule
        """
        r_max = 0
        for coord, number in zip(self.mol.cart_coords, self.mol.atomic_numbers):
            radius = (
                np.sqrt(np.sum(coord * coord)) + self.tm.get_tol(number, number) * 0.5
            )
            if radius > r_max:
                r_max = radius
        self.radius = r_max
        # reestimate the radius if it has stick shape
        rmax = max([self.box.width,self.box.height,self.box.length])
        rmin = min([self.box.width,self.box.height,self.box.length])
        if rmax/rmin > 3 and rmax >12:
            self.radius = rmin

    def has_stick_shape(self):
        """
        check if the molecule is stick-like
        """
        sizes = [self.box.width,self.box.height,self.box.length]
        sizes.sort()
        if sizes[2]>15: #and sizes[2]/sizes[0]>2 and sizes[2]/sizes[1]>2:
            return True
        else:
            return False

    def get_symbols(self):
        self.symbols = [specie.name for specie in self.mol.species]

    def get_tols_matrix(self):
        """
        Returns: a 2D matrix which is used internally for distance checking.
        """
        numbers = self.mol.atomic_numbers
        tols = np.zeros((len(numbers), len(numbers)))
        for i1, number1 in enumerate(numbers):
            for i2, number2 in enumerate(numbers):
                tols[i1][i2] = self.tm.get_tol(number1, number2)
                # allow hydrogen bond
                if [number1, number2] in [[1,7], [1,8], [1,9], [7,1], [8,1], [9,1]]:
                    tols[i1][i2] *= 0.9
        if len(self.mol)==1:
            tols *= 0.8 # if only one atom, reduce the tolerance
        self.tols_matrix = tols

    def show(self):
        """
        show the molecule
        """
        from pyxtal.viz import display_molecules
        return display_molecules([self.mol])


    def rdkit_mol_init(self, smile):
        """
        initialize the mol xyz and torsion list
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem

        smarts_torsion="[*]~[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]~[*]"
        mol = Chem.MolFromSmiles(smile)
        pattern_tor = Chem.MolFromSmarts(smarts_torsion)
        torsion = list(mol.GetSubstructMatches(pattern_tor))
        torsionlist = cleaner(torsion)

        mol = Chem.AddHs(mol)
        symbols = []
        for id in range(mol.GetNumAtoms()):
            symbols.append(mol.GetAtomWithIdx(id).GetSymbol())

        AllChem.EmbedMultipleConfs(mol, numConfs=2, randomSeed=0xf00d)
        conf = mol.GetConformer(0)
        #print("Init: ", conf.GetPositions())
        xyz = self.align(conf)
        #print("Init: ", xyz[:3])
        return symbols, xyz, torsionlist

    def rdkit_mol(self, smile):
        """
        initialize the mol xyz and torsion list
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=3, randomSeed=0xf00d)
        #print("To_dict", mol.GetConformer(0).GetPositions())
        return mol

    def align(self, conf, reflect=False):
        """
        Align the molecule and return the xyz
        The default CanonicalizeConformer function may also include inversion
        """
        from rdkit.Chem import rdMolTransforms as rdmt

        #rotation
        #print("align: "); print(conf.GetPositions())
        trans = rdmt.ComputeCanonicalTransform(conf)
        if np.linalg.det(trans[:3,:3]) < 0:
            trans[:3,:3] *= -1

        if reflect:
            trans[:3,:3] *= -1
        #print(trans)
        rdmt.TransformConformer(conf, trans)
        
        #print("rot", conf.GetPositions()[:3])
        #translation
        pt = rdmt.ComputeCentroid(conf)
        center = np.array([pt.x, pt.y, pt.z])
        xyz = conf.GetPositions() - center
        #print("return", xyz[:3])
        return xyz

    def get_center(self, xyz):
        """
        get the molecular center for a transformed xyz
        """
        if self.smile is None:
            return np.mean(xyz, axis=0)
        else:
            # from rdkit
            from rdkit.Geometry import Point3D
            from rdkit.Chem import rdMolTransforms as rdmt

            conf1 = self.rdkit_mol(self.smile).GetConformer(0)
            for i in range(conf1.GetNumAtoms()):
                x, y, z = xyz[i]
                conf1.SetAtomPosition(i, Point3D(x,y,z))
            pt = rdmt.ComputeCentroid(conf1)
            return np.array([pt.x, pt.y, pt.z])

    def get_principle_axes(self, xyz):
        """
        get the principle axis for a rotated xyz
        """
        if self.smile is None:
            Inertia = get_inertia_tensor(xyz)
            _, matrix = np.linalg.eigh(Inertia)
            return matrix

        else:
            from rdkit.Geometry import Point3D
            from rdkit.Chem import rdMolTransforms as rdmt

            conf1 = self.rdkit_mol(self.smile).GetConformer(0)
            for i in range(len(self.mol)):
                x,y,z = xyz[i]
                conf1.SetAtomPosition(i,Point3D(x,y,z))

            return rdmt.ComputePrincipalAxesAndMoments(conf1)

    def get_torsion_angles(self, xyz):
        """
        get the torsion angles
        """
        from rdkit.Geometry import Point3D
        from rdkit.Chem import rdMolTransforms as rdmt
        mol = self.rdkit_mol(self.smile)
        conf = mol.GetConformer(0)

        for i in range(len(self.mol)):
            x,y,z = xyz[i]
            conf.SetAtomPosition(i,Point3D(x,y,z))

        angs = []
        for torsion in self.torsionlist:
            (i, j, k, l) = torsion
            angs.append(rdmt.GetDihedralDeg(conf, i, j, k, l))
        return angs
    
    def set_torsion_angles(self, conf, angles, reflect=False):
        """
        reset the torsion angles and update molecular xyz
        """
        from rdkit.Chem import rdMolTransforms as rdmt

        for id, torsion in enumerate(self.torsionlist):
            (i, j, k, l) = torsion
            rdmt.SetDihedralDeg(conf, i, j, k, l, angles[id])
        
        xyz = self.align(conf, reflect)
        return xyz

    def get_orientation(self, xyz, rtol=0.15):
        """
        get orientation, needs to check the tolerance
        """
        from rdkit.Geometry import Point3D
        from rdkit.Chem import rdMolAlign, RemoveHs, rdmolfiles, rdMolTransforms
        mol = self.rdkit_mol(self.smile)
        
        conf0 = mol.GetConformer(0)
        conf1 = mol.GetConformer(1)
        conf2 = mol.GetConformer(2)
        angs = self.get_torsion_angles(xyz)
        xyz0 = self.set_torsion_angles(conf0, angs) #conf0 with aligned
        xyz1 = self.set_torsion_angles(conf0, angs, True) #conf0 with aligned

        for i in range(len(self.mol)):
            x0,y0,z0 = xyz0[i]
            x1,y1,z1 = xyz1[i]
            x,y,z = xyz[i]
            conf0.SetAtomPosition(i,Point3D(x0,y0,z0))
            conf1.SetAtomPosition(i,Point3D(x,y,z))
            conf2.SetAtomPosition(i,Point3D(x1,y1,z1))

        mol = RemoveHs(mol)
        rmsd1, trans1 = rdMolAlign.GetAlignmentTransform(mol, mol, 1, 0)
        rmsd2, trans2 = rdMolAlign.GetAlignmentTransform(mol, mol, 1, 2)
        tol = rtol*mol.GetNumAtoms()

        #print(rmsd1, rmsd2)
        if rmsd1 < tol:
            trans = trans1[:3,:3].T
            r = Rotation.from_matrix(trans)
            return r.as_euler('zxy', degrees=True), rmsd1, False
        elif rmsd2 < tol:
            trans = trans2[:3,:3].T
            r = Rotation.from_matrix(trans)
            return r.as_euler('zxy', degrees=True), rmsd2, True
        else:
            print(rmsd1, rmsd2)
            #rdmolfiles.MolToXYZFile(mol, '1.xyz', 0)
            #rdmolfiles.MolToXYZFile(mol, '2.xyz', 1)
            #rdmolfiles.MolToXYZFile(mol, '3.xyz', 2)
            print(self.get_torsion_angles(xyz))   
            print(self.get_torsion_angles(xyz0))   
            print(self.get_torsion_angles(xyz1))   
            raise ValueError("Problem in conformer")


    def reset_positions(self, coors):
        """
        reset the coordinates
        """
        from pymatgen.core.sites import Site
        if len(coors) != len(self.mol._sites):
            raise ValueError("number of atoms is inconsistent!")
        else:
            for i, coor in enumerate(coors):
                _site = self.mol._sites[i]
                new_site = Site(_site.species, coor, properties=_site.properties)
                self.mol._sites[i] = new_site

    def apply_inversion(self):
        """
        reset the coordinates
        """
        from pymatgen.core.sites import Site
        xyz = self.mol.cart_coords
        center = self.get_center(xyz)
        xyz -= center
        xyz *= -1
        self.reset_positions(xyz)

class Box:
    """
    Class for storing the binding box for a molecule. Box is oriented along the x, y, and
    z axes.

    Args:
        minx: the minimum x value
        maxx: the maximum x value
        miny: the minimum y value
        maxy: the maximum y value
        minz: the minimum z value
        maxz: the maximum z value
    """

    def __init__(self, minx, maxx, miny, maxy, minz, maxz):
        self.minx = float(minx)
        self.maxx = float(maxx)
        self.miny = float(miny)
        self.maxy = float(maxy)
        self.minz = float(minz)
        self.maxz = float(maxz)

        self.width = float(abs(maxx - minx))
        self.length = float(abs(maxy - miny))
        self.height = float(abs(maxz - minz))

        self.minl = min(self.width, self.length, self.height)
        self.maxl = max(self.width, self.length, self.height)
        for x in (self.width, self.length, self.height):
            if x <= self.maxl and x >= self.minl:
                self.midl = x

        self.volume = float(self.width * self.length * self.height)


class Orientation:
    """
    Stores orientations for molecules based on vector constraints.
    Can be stored to regenerate orientations consistent with a given constraint
    vector, without re-calling orientation_in_wyckoff_position. Allows for
    generating orientations which differ only in their rotation about a given
    axis.

    Args:
        matrix: a 3x3 rotation matrix describing the orientation (and/or
            inversion) to store
        degrees: the number of degrees of freedom...
            0 - The orientation refers to a single rotation matrix
            1 - The orientation can be rotated about a single axis
            2 - The orientation can be any pure rotation matrix

        axis:
            an optional axis about which the orientation can rotate. Only used
            if degrees is equal to 1
    """

    def __init__(self, matrix=None, degrees=2, axis=None):
        self.matrix = np.array(matrix)
        self.degrees = degrees  # The number of degrees of freedom.
        if degrees == 1:
            if axis is None:
                raise ValueError("axis is required for orientation")
            else:
                axis /= np.linalg.norm(axis)
        self.axis = axis

        self.r = Rotation.from_matrix(self.matrix)  # scipy transform.Rotation class
        self.angle = None

    def __str__(self):
        s = "-------PyXtal.molecule.Orientation class----\n"
        s += "degree of freedom: {:d}\n".format(self.degrees)
        s += "Rotation matrix:\n"
        s += "{:6.3f} {:6.3f} {:6.3f}\n".format(*self.matrix[:,0])
        s += "{:6.3f} {:6.3f} {:6.3f}\n".format(*self.matrix[:,1])
        s += "{:6.3f} {:6.3f} {:6.3f}\n".format(*self.matrix[:,2])
        if self.axis is not None:
            s += "Rotation axis\n"
            s += "{:6.2f} {:6.2f} {:6.3f}\n".format(*self.axis)
        return s

    def reset_matrix(self, matrix):
        self.matrix = matrix
        self.r = Rotation.from_matrix(matrix)

    def __repr__(self):
        return str(self)

    def copy(self):
        return deepcopy(self)

    def save_dict(self):
        dict0 = {"matrix": self.matrix,
                 "degrees": self.degrees,
                 "axis": self.axis
                }
        return dict0

    @classmethod
    def load_dict(cls, dicts):
        matrix = dicts['matrix']
        degrees = dicts['degrees']
        axis = dicts['axis']
        return cls(matrix, degrees, axis)

    def change_orientation(self, angle="random", flip=False):
        """
        Allows for specification of an angle (possibly random) to
        rotate about the constraint axis.

        Args:
            angle: an angle to rotate about the constraint axis.
            If "random", chooses a random rotation angle.
            If self.degrees==2, chooses a random rotation matrix.
            If self.degrees==1, only apply on angle
            If self.degrees==0, no change

        """
        if self.degrees >= 1:
            # choose the axis
            if self.axis is None:
                axis = np.random.rand(3) - 0.5
                self.axis = axis / np.linalg.norm(axis)

            # parse the angle
            if angle == "random":
                angle = np.random.rand() * np.pi * 2
            self.angle = angle

            # update the matrix
            r1 = Rotation.from_rotvec(self.angle * self.axis)

            if self.degrees == 2 and flip:
                if np.random.random()>0.5:
                    ax = choice(['x','y','z'])
                    angle0 = choice([90, 180, 270])
                    r2 = Rotation.from_euler(ax, angle0, degrees=True)
                    r1 = r2*r1
            self.r = r1 * self.r
            self.matrix = self.r.as_matrix()

    def rotate_by_matrix(self, matrix, ignore_constraint=True):
        """
        rotate

        Args:
            matrix: 3*3 rotation matrix

        """
        if not ignore_constraint:
            if self.degrees == 0:
                raise ValueError("cannot rotate")
            elif self.degrees == 1:
                axis = self.axis
                vec = Rotation.from_matrix(matrix).as_rotvec()
                if angle(vec, self.axis) > 1e-2 and angle(vec, -self.axis) > 1e-2:
                    raise ValueError("must rotate along the given axis")
        else:
            axis = None

        matrix = matrix.dot(self.matrix)
        return Orientation(matrix, self.degrees, axis)

    def get_matrix(self, angle="random"):
        """
        Generate a 3x3 rotation matrix consistent with the orientation's
        constraints. Allows for specification of an angle (possibly random) to
        rotate about the constraint axis.

        Args:
            angle: an angle to rotate about the constraint axis. If "random",
                chooses a random rotation angle. If self.degrees==2, chooses a
                random 3d rotation matrix to multiply by. If the original matrix
                is wanted, set angle=0, or call self.matrix

        Returns:
            a 3x3 rotation (and/or inversion) matrix (numpy array)
        """
        if self.degrees == 2:
            if angle == "random":
                axis = np.random.sample(3)
                axis = axis / np.linalg.norm(axis)
                angle = np.random.random() * np.pi * 2
            else:
                axis = self.axis
            return Rotation.from_rotvec(angle * axis).as_matrix()

        elif self.degrees == 1:
            if angle == "random":
                angle = np.random.random() * np.pi * 2
            else:
                angle = self.angle
            return Rotation.from_rotvec(angle * self.axis).as_matrix()

        elif self.degrees == 0:
            return self.matrix

    def get_op(self): #, angle=None):
        """
        Generate a SymmOp object consistent with the orientation's
        constraints. Allows for specification of an angle (possibly random) to
        rotate about the constraint axis.

        Args:
            angle: an angle to rotate about the constraint axis. If "random",
                chooses a random rotation angle. If self.degrees==2, chooses a
                random 3d rotation matrix to multiply by. If the original matrix
                is wanted, set angle=0, or call self.matrix

        Returns:
            pymatgen.core.structure. SymmOp object
        """
        #if angle is not None:
        #    self.change_orientation(angle)
        return SymmOp.from_rotation_and_translation(self.matrix, [0, 0, 0])

    def random_orientation(self):
        """
        Applies random rotation (if possible) and returns a new orientation with
        the new base matrix.

        Returns:
            a new orientation object with a different base rotation matrix
        """

        self.change_orientation()
        return self

    def get_Euler_angles(self):
        """
        get the Euler angles
        """
        return self.r.as_euler('zxy', degrees=True)


def get_inertia_tensor(coords):
    """
    Calculate the symmetric inertia tensor for a Molecule
    the principal axes of symmetry.

    Args:
        coords: [N, 3] array of coordinates

    Returns:
        a 3x3 numpy array representing the inertia tensor
    """
    coords -= np.mean(coords, axis=0)
    Inertia = np.zeros([3,3])
    Inertia[0,0] = np.sum(coords[:,1]**2 + coords[:,2]**2)
    Inertia[1,1] = np.sum(coords[:,0]**2 + coords[:,2]**2)
    Inertia[2,2] = np.sum(coords[:,0]**2 + coords[:,1]**2)
    Inertia[0,1] = Inertia[1,0] = -np.sum(coords[:,0]*coords[:,1])
    Inertia[0,2] = Inertia[2,0] = -np.sum(coords[:,0]*coords[:,2])
    Inertia[1,2] = Inertia[2,1] = -np.sum(coords[:,1]*coords[:,2])

    return Inertia


def reoriented_molecule(mol): #, nested=False):
    """
    Reorient a molecule so that its principal axes are aligned with the
    identity matrix.

    Args:
        mol: a Molecule object

    Returns:
        new_mol: a reoriented copy of the original molecule.
        P: the 3x3 rotation matrix used to obtain it.
    """
    coords = mol.cart_coords
    numbers = mol.atomic_numbers
    coords -= np.mean(coords, axis=0)
    A = get_inertia_tensor(coords)
    # Store the eigenvectors of the inertia tensor
    P = np.linalg.eigh(A)[1]
    if np.linalg.det(P) < 0:
        P[0] *= -1
    coords = np.dot(coords, P)
    return Molecule(numbers, coords), P


def get_symmetry(mol, already_oriented=False):
    """
    Return a molecule's point symmetry.
    Note: for linear molecules, infinitessimal rotations are treated as 6-fold
    rotations, which works for 3d and 2d point groups.

    Args:
        mol: a Molecule object
        already_oriented: whether or not the principle axes of mol are already
            reoriented. Can save time if True, but is not required.

    Returns:
        a list of SymmOp objects which leave the molecule unchanged when applied
    """
    # For single atoms, we cannot represent the point group using a list of operations
    if len(mol) == 1:
        return []
    pga = PointGroupAnalyzer(mol)
    # Handle linear molecules
    if "*" in pga.sch_symbol:
        if not already_oriented:
            # Reorient the molecule
            oriented_mol, P = reoriented_molecule(mol)
            pga = PointGroupAnalyzer(oriented_mol)
        pg = pga.get_pointgroup()
        symm_m = []
        for op in pg:
            symm_m.append(op)
        # Add 12-fold  and reflections in place of ininitesimal rotation
        for axis in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
            # op = SymmOp.from_rotation_and_translation(aa2matrix(axis, np.pi/6), [0,0,0])
            m1 = Rotation.from_rotvec(np.pi / 6 * axis).as_matrix()
            op = SymmOp.from_rotation_and_translation(m1, [0, 0, 0])
            if pga.is_valid_op(op):
                symm_m.append(op)
                # Any molecule with infinitesimal symmetry is linear;
                # Thus, it possess mirror symmetry for any axis perpendicular
                # To the rotational axis. pymatgen does not add this symmetry
                # for all linear molecules - for example, hydrogen
                if axis == [1, 0, 0]:
                    symm_m.append(SymmOp.from_xyz_string("x,-y,z"))
                    symm_m.append(SymmOp.from_xyz_string("x,y,-z"))
                    #r = SymmOp.from_xyz_string("-x,y,-z")
                elif axis == [0, 1, 0]:
                    symm_m.append(SymmOp.from_xyz_string("-x,y,z"))
                    symm_m.append(SymmOp.from_xyz_string("x,y,-z"))
                    #r = SymmOp.from_xyz_string("-x,-y,z")
                elif axis == [0, 0, 1]:
                    symm_m.append(SymmOp.from_xyz_string("-x,y,z"))
                    symm_m.append(SymmOp.from_xyz_string("x,-y,z"))
                    #r = SymmOp.from_xyz_string("x,-y,-z")
                # Generate a full list of SymmOps for the molecule's pointgroup
                symm_m = generate_full_symmops(symm_m, 1e-3)
                break
        # Reorient the SymmOps into mol's original frame
        if not already_oriented:
            new = []
            for op in symm_m:
                new.append(P.inverse * op * P)
            return new
        elif already_oriented:
            return symm_m
    # Handle nonlinear molecules
    else:
        pg = pga.get_pointgroup()
        symm_m = []
        for op in pg:
            symm_m.append(op)
        return symm_m


def orientation_in_wyckoff_position(
    mol,
    wyckoff_position,
    randomize=True,
    exact_orientation=False,
    already_oriented=False,
    allow_inversion=True,
    rtol = 1e-2,
):
    """
    Tests if a molecule meets the symmetry requirements of a Wyckoff position,
    and returns the valid orientations.

    Args:
        mol: a Molecule object. Orientation is arbitrary
        wyckoff_position: a pyxtal.symmetry.Wyckoff_position object
        randomize: whether or not to apply a random rotation consistent with
            the symmetry requirements
        exact_orientation: whether to only check compatibility for the provided
            orientation of the molecule. Used within general case for checking.
            If True, this function only returns True or False
        already_oriented: whether or not to reorient the principle axes
            when calling get_symmetry. Setting to True can remove redundancy,
            but is not necessary
        allow_inversion: whether or not to allow chiral molecules to be
            inverted. Should only be True if the chemical and biological
            properties of the mirror image are known to be suitable for the
            desired application

    Returns:
        a list of operations.Orientation objects which can be applied to the
        molecule while allowing it to satisfy the symmetry requirements of the
        Wyckoff position. If no orientations are found, returns False.
    """
    # For single atoms, there are no constraints
    if len(mol) == 1:
        return [Orientation([[1, 0, 0], [0, 1, 0], [0, 0, 1]], degrees=2)]

    wyckoffs = wyckoff_position.ops
    w_symm = wyckoff_position.symmetry_m

    # Obtain the Wyckoff symmetry
    symm_w = w_symm[0]
    pga = PointGroupAnalyzer(mol)

    # Check exact orientation
    if exact_orientation is True:
        mo = deepcopy(mol)
        valid = True
        for op in symm_w:
            if not pga.is_valid_op(op):
                valid = False
        if valid is True:
            return True
        elif valid is False:
            return False

    # Obtain molecular symmetry, exact_orientation==False
    symm_m = get_symmetry(mol, already_oriented=already_oriented)
    # Store OperationAnalyzer objects for each molecular SymmOp
    chiral = True
    opa_m = []
    for op_m in symm_m:
        opa = OperationAnalyzer(op_m)
        opa_m.append(opa)
        if opa.type == "rotoinversion":
            chiral = False
        elif opa.type == "inversion":
            chiral = False

    # If molecule is chiral and allow_inversion is False,
    # check if WP breaks symmetry
    if chiral is True:
        if allow_inversion is False:
            for op in wyckoffs:
                if np.linalg.det(op.rotation_matrix) < 0:
                    printx(
                        "Warning: cannot place chiral molecule in spagegroup", priority=2,
                    )
                    return False

    # Store OperationAnalyzer objects for each Wyckoff symmetry SymmOp
    opa_w = []
    for op_w in symm_w:
        opa_w.append(OperationAnalyzer(op_w))

    # Check for constraints from the Wyckoff symmetry...
    # If we find ANY two constraints (SymmOps with unique axes), the molecule's
    # point group MUST contain SymmOps which can be aligned to these particular
    # constraints. However, there may be multiple compatible orientations of the
    # molecule consistent with these constraints
    constraint1 = None
    constraint2 = None
    for i, op_w in enumerate(symm_w):
        if opa_w[i].axis is not None:
            constraint1 = opa_w[i]
            for j, op_w in enumerate(symm_w):
                if opa_w[j].axis is not None:
                    dot = np.dot(opa_w[i].axis, opa_w[j].axis)
                    if (not np.isclose(dot, 1, rtol=rtol)) and (
                        not np.isclose(dot, -1, rtol=rtol)
                    ):
                        constraint2 = opa_w[j]
                        break
            break
    # Indirectly store the angle between the constraint axes
    if constraint1 is not None and constraint2 is not None:
        dot_w = np.dot(constraint1.axis, constraint2.axis)
    # Generate 1st consistent molecular constraints
    constraints_m = []
    if constraint1 is not None:
        for i, opa1 in enumerate(opa_m):
            if opa1.is_conjugate(constraint1):
                constraints_m.append([opa1, []])
                # Generate 2nd constraint in opposite direction
                extra = deepcopy(opa1)
                extra.axis = [opa1.axis[0] * -1, opa1.axis[1] * -1, opa1.axis[2] * -1]
                constraints_m.append([extra, []])

    # Remove redundancy for the first constraints
    list_i = list(range(len(constraints_m)))
    list_j = list(range(len(constraints_m)))
    copy = deepcopy(constraints_m)
    for i, c1 in enumerate(copy):
        if i in list_i:
            for j, c2 in enumerate(copy):
                if i > j and j in list_j and j in list_i:
                    # Check if axes are colinear
                    if np.isclose(np.dot(c1[0].axis, c2[0].axis), 1, rtol=rtol):
                        list_i.remove(j)
                        list_j.remove(j)
                    # Check if axes are symmetrically equivalent
                    else:
                        cond1 = False
                        # cond2 = False
                        for opa in opa_m:
                            if opa.type == "rotation":
                                op = opa.op
                                if np.isclose(
                                    np.dot(op.operate(c1[0].axis), c2[0].axis),
                                    1,
                                    rtol=5*rtol,
                                ):
                                    cond1 = True
                                    break
                        if cond1 is True:  # or cond2 is True:
                            list_i.remove(j)
                            list_j.remove(j)
    c_m = deepcopy(constraints_m)
    constraints_m = []
    for i in list_i:
        constraints_m.append(c_m[i])

    # Generate 2nd consistent molecular constraints
    valid = list(range(len(constraints_m)))
    if constraint2 is not None:
        for i, c in enumerate(constraints_m):
            opa1 = c[0]
            for j, opa2 in enumerate(opa_m):
                if opa2.is_conjugate(constraint2):
                    dot_m = np.dot(opa1.axis, opa2.axis)
                    # Ensure that the angles are equal
                    if abs(dot_m - dot_w) < 0.02 or abs(dot_m + dot_w) < 0.02:
                        constraints_m[i][1].append(opa2)
                        # Generate 2nd constraint in opposite direction
                        extra = deepcopy(opa2)
                        extra.axis = [
                            opa2.axis[0] * -1,
                            opa2.axis[1] * -1,
                            opa2.axis[2] * -1,
                        ]
                        constraints_m[i][1].append(extra)
            # If no consistent constraints are found, remove first constraint
            if constraints_m[i][1] == []:
                valid.remove(i)
    copy = deepcopy(constraints_m)
    constraints_m = []
    for i in valid:
        constraints_m.append(copy[i])

    # Generate orientations consistent with the possible constraints
    orientations = []
    # Loop over molecular constraint sets
    for c1 in constraints_m:
        v1 = c1[0].axis
        v2 = constraint1.axis
        T = rotate_vector(v1, v2)
        # If there is only one constraint
        if c1[1] == []:
            o = Orientation(T, degrees=1, axis=constraint1.axis)
            orientations.append(o)
        else:
            # Loop over second molecular constraints
            for opa in c1[1]:
                phi = angle(constraint1.axis, constraint2.axis)
                phi2 = angle(constraint1.axis, np.dot(T, opa.axis))
                if np.isclose(phi, phi2, rtol=rtol):
                    r = np.sin(phi)
                    c = np.linalg.norm(np.dot(T, opa.axis) - constraint2.axis)
                    theta = np.arccos(1 - (c ** 2) / (2 * (r ** 2)))
                    # R = aa2matrix(constraint1.axis, theta)
                    R = Rotation.from_rotvec(theta * constraint1.axis).as_matrix()
                    T2 = np.dot(R, T)
                    a = angle(np.dot(T2, opa.axis), constraint2.axis)
                    if not np.isclose(a, 0, rtol=rtol):
                        T2 = np.dot(np.linalg.inv(R), T)
                    o = Orientation(T2, degrees=0)
                    orientations.append(o)

    # Ensure the identity orientation is checked if no constraints are found
    if constraints_m == []:
        o = Orientation(np.identity(3), degrees=2)
        orientations.append(o)

    # Remove redundancy from orientations
    list_i = list(range(len(orientations)))
    list_j = list(range(len(orientations)))
    for i, o1 in enumerate(orientations):
        if i in list_i:
            for j, o2 in enumerate(orientations):
                if i > j and j in list_j and j in list_i:
                    # m1 = o1.get_matrix(angle=0)
                    # m2 = o2.get_matrix(angle=0)
                    m1 = o1.matrix
                    m2 = o2.matrix
                    new_op = SymmOp.from_rotation_and_translation(
                        np.dot(m2, np.linalg.inv(m1)), [0, 0, 0]
                    )
                    P = SymmOp.from_rotation_and_translation(np.linalg.inv(m1), [0, 0, 0])
                    old_op = P * new_op * P.inverse
                    if pga.is_valid_op(old_op):
                        list_i.remove(j)
                        list_j.remove(j)
    #copies = deepcopy(orientations)
    orientations_new = []
    for i in list_i:
        orientations_new.append(orientations[i])

    #Check each of the found orientations for consistency with the Wyckoff pos.
    #If consistent, put into an array of valid orientations
    allowed = []
    for o in orientations_new:
        if randomize is True:
            op = o.get_op()
        elif randomize is False:
            op = o.get_op() #do not change
        mo = deepcopy(mol)
        mo.apply_operation(op)
        if orientation_in_wyckoff_position(
                mo, wyckoff_position, exact_orientation=True,
                randomize=False, allow_inversion=allow_inversion
        ):
            allowed.append(o)
    if allowed == []:
        return False
    else:
        return allowed

def make_graph(mol, tol=0.2):
    """
    make graph object for the input molecule
    """
    #print("making graphs")
    G = nx.Graph()
    names = {}
    for i, site in enumerate(mol._sites):
        names[i] = site.specie.value

    for i in range(len(mol)-1):
        site1 = mol.sites[i]
        for j in range(i+1, len(mol)):
            site2 = mol.sites[j]
            #remove short X-H distances
            if names[i] == "H" and names[j]=="H":
                factor = -0.5
            elif [names[i], names[j]] in [["S","S"], ["S","O"], ["O","S"], ["F","O"], ["O","F"]]:
                factor = 0.05
            elif "H" in [names[i], names[j]]:
                factor = 0.5
            elif [names[i], names[j]] in [["N", "S"], ["S", "N"]]:
                factor = 1.25
            else:
                factor = 1.0
            try:
                if CovalentBond.is_bonded(site1, site2, factor*tol):
                    G.add_edge(i,j)
                    #if 'S' in [names[i], names[j]]: print(names[i], names[j], mol.get_distance(i, j))
            except ValueError:
                pass
    nx.set_node_attributes(G, names, 'name')

    return G

def compare_mol_connectivity(mol1, mol2, ignore_name=False):
    """
    Compare two molecules by connectivity
    """

    G1 = make_graph(mol1)
    G2 = make_graph(mol2)
    if ignore_name:
        GM = nx.isomorphism.GraphMatcher(G1, G2)
    else:
        fun = lambda n1, n2: n1['name'] == n2['name']
        GM = nx.isomorphism.GraphMatcher(G1, G2, node_match=fun)

    return GM.is_isomorphic(), GM.mapping
