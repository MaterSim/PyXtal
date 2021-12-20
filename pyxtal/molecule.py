"""
Module for handling molecules. 
"""
# Imports
import os
from pkg_resources import resource_filename
from copy import deepcopy
from operator import itemgetter
from random import choice, random
import numpy as np
from scipy.spatial.transform import Rotation
import networkx as nx

# External Libraries
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer, generate_full_symmops
from monty.serialization import loadfn

# PyXtal imports
from pyxtal.symmetry import Group
from pyxtal.tolerance import Tol_matrix
from pyxtal.database.element import Element
from pyxtal.operations import SymmOp, OperationAnalyzer, rotate_vector, angle
from pyxtal.database.collection import Collection
from pyxtal.msg import ConformerError

# Define functions
bonds = loadfn(resource_filename("pyxtal", "database/bonds.json"))
molecule_collection = Collection("molecules")

def find_id_from_smile(smile):
    """
    Find the positions of rotatable bonds in the molecule.
    """
    def cleaner(list_to_clean):
        """
        Remove duplicate torsion from a list of atom index tuples.
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

    if smile in ["Cl-", "F-", "Br-", "I-", "Li+", "Na+"]:
        return []
    else:
        from rdkit import Chem
        smarts_torsion1="[*]~[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]~[*]"
        smarts_torsion2="[*]~[^2]=[^2]~[*]" # C=C bonds
        #smarts_torsion2="[*]~[^1]#[^1]~[*]" # C-C triples bonds, to be fixed

        mol = Chem.MolFromSmiles(smile)
        pattern_tor1 = Chem.MolFromSmarts(smarts_torsion1)
        torsion1 = cleaner(list(mol.GetSubstructMatches(pattern_tor1)))
        pattern_tor2 = Chem.MolFromSmarts(smarts_torsion2)
        torsion2 = cleaner(list(mol.GetSubstructMatches(pattern_tor2)))
        tmp = cleaner(torsion1+torsion2)
        torsions = []
        for t in tmp:
            (i, j, k, l) = t
            b = mol.GetBondBetweenAtoms(j,k)
            if not b.IsInRing():
                torsions.append(t)
        return torsions

def generate_molecules(smile, wps=None, N_iter=4, N_conf=10, tol=0.5):
    """
    generate pyxtal_molecules from smiles codes.
    
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    torsionlist = find_id_from_smile(smile)

    def get_conformers(smile, seed):
        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)
        ps = AllChem.ETKDGv3()
        ps.randomSeed=seed
        ps.runeRmsThresh=tol
        AllChem.EmbedMultipleConfs(mol, max([4, 4*len(torsionlist)]), ps)
        return mol

    m0 = pyxtal_molecule(smile+'.smi', fix=True)
    _, valid = m0.get_orientations_in_wps(wps)
    if valid:
        mols = [m0]
    else:
        mols = []

    for i in range(N_iter):
        mol = get_conformers(smile, seed=i)
        res = AllChem.MMFFOptimizeMoleculeConfs(mol)
        for id, conf in enumerate(mol.GetConformers()):
            m = m0.copy()
            xyz = m.align(conf)
            m.reset_positions(xyz)
            m.get_symmetry(True)
            m.energy = res[id][1]
            _, valid = m.get_orientations_in_wps(wps)
            if valid:
                add = True
            else:
                add = False
            if add:
                match = False
                for mol in mols:
                    rms, _, _ = mol.get_rmsd(xyz)
                    if rms < tol:
                        match = True
                        break
                if not match:
                    mols.append(m)
            if len(mols) == N_conf:
                break
    #for m in mols:
    #    print(m.energy, m.pga.sch_symbol, len(torsionlist))
    return mols

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

    def __init__(self, 
                 mol=None, 
                 symmetrize=True, 
                 fix=False, 
                 torsions=None,
                 seed = None,
                 tm=Tol_matrix(prototype="molecular")):

        mo = None
        self.smile = None
        self.torsionlist = None
        if seed is None: seed = 0xf00d
        self.seed = seed

        # Parse molecules: either file or molecule name
        if type(mol) == str:
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
                    res = self.rdkit_mol_init(tmp[0], fix, torsions)
                    (symbols, xyz, self.torsionlist) = res
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

        # Molecule and symmetry analysis
        self.props = mo.site_properties
        if len(mo) > 1:
            if symmetrize:
                pga = PointGroupAnalyzer(mo, 0.5)
                mo = pga.symmetrize_molecule()["sym_mol"]
        self.mol = mo
        self.get_symmetry()

        self.tm = tm
        self.box = self.get_box()
        self.volume = self.box.volume
        self.get_symbols()
        self.get_radius()
        self.tols_matrix = self.get_tols_matrix()
        xyz = self.mol.cart_coords
        self.reset_positions(xyz-self.get_center(xyz))

    def __str__(self):
        return '[' + self.name + ']'

    def save_dict(self):
        """
        save the object as a dictionary
        """
        return self.mol.as_dict()

    @classmethod
    def load_dict(cls, dicts):
        """
        load the molecule from a dictionary
        """
        mol = Molecule.from_dict(dicts)
        return cls(mol)

    def copy(self):
        """
        simply copy the structure
        """
        return deepcopy(self)

    def swap_axis(self, ax):
        """
        swap the molecular axis
        """
        coords = self.mol.cart_coords[:, ax]
        mo = Molecule(self.symbols, coords)

        return pyxtal_molecule(mo, self.tm)

    #def add_site_props(self, mo):
    #    """
    #    add site properties
    #    """
    #    if len(self.props) > 0:
    #        for key in self.props.keys():
    #            mo.add_site_property(key, self.props[key])
    #    return mo

    def get_box(self, padding=None):
        """
        Given a molecule, find a minimum orthorhombic box containing it.
        Size is calculated using min and max x, y, and z values,
        plus the padding defined by the vdw radius
        For best results, call oriented_molecule first.

        Args:
            padding: float (default is 3.4 according to the vdw radius of C)

        Returns:
            box: a Box object
        """
        mol, P = reoriented_molecule(self.mol)
        xyz = mol.cart_coords
        dims = [0, 0, 0]
        for i in range(3):
            dims[i] = np.max(xyz[:,i]) - np.min(xyz[:,i])
            if padding is not None:
                dims[i] += padding
                dims[i] = max([dims[i], 2.0]) #for planar molecules
            else:
                ids = np.argsort(xyz[:, i])
                r = Element(mol[ids[0]].species_string).vdw_radius 
                r += Element(mol[ids[-1]].species_string).vdw_radius 
                dims[i] = max([dims[i]+r, 3.4]) #for planar molecules
        return Box(dims)

    def get_box_coordinates(self, xyz, padding=0, resolution=1.0):
        """
        create the points cloud to describe the molecular box

        Args:
            center: molecular center position
            orientation: orientation matrix
            padding: padding of the box
            resolution: float in angstrom

        Return:
            cell: box axis
            vertices: [N, 3] np.array, Cartesian coordinates to describe the box.
            center: box center
        """
        cell = self.get_principle_axes(xyz).T
        center = self.get_center(xyz) #, geometry=True)
        box = self.get_box(padding)
        #print(box)
        w, h, l = box.width, box.height, box.length      
        cell[0,:] *= l
        cell[1,:] *= w
        cell[2,:] *= h
        x_ = np.linspace(-1/2, 1/2, int(l/resolution)+1)
        y_ = np.linspace(-1/2, 1/2, int(w/resolution)+1)
        z_ = np.linspace(-1/2, 1/2, int(h/resolution)+1)
        
        #XY
        #print(len(x_), len(y_), len(z_))
        x, y = np.meshgrid(x_, y_, indexing='ij')
        size = len(x.flatten())
        xy = np.zeros([size*2, 3])
        xy[:size,0] = x.flatten()
        xy[size:,0] = x.flatten()
        xy[:size,1] = y.flatten()
        xy[size:,1] = y.flatten()
        xy[:size,2] = -0.5
        xy[size:,2] = 0.5
        #print(xy.shape)
        #print(xy)

        y, z = np.meshgrid(y_, z_, indexing='ij')
        size = len(y.flatten())
        yz = np.zeros([size*2, 3])
        yz[:size,1] = y.flatten()
        yz[size:,1] = y.flatten()
        yz[:size,2] = z.flatten()
        yz[size:,2] = z.flatten()
        yz[:size,0] = -0.5
        yz[size:,0] = 0.5

        x, z = np.meshgrid(x_, z_, indexing='ij')
        size = len(z.flatten())
        xz = np.zeros([size*2, 3])
        xz[:size,0] = x.flatten()
        xz[size:,0] = x.flatten()
        xz[:size,2] = z.flatten()
        xz[size:,2] = z.flatten()
        xz[:size,1] = -0.5
        xz[size:,1] = 0.5

        vertices = np.zeros([len(xy)+len(yz)+len(xz), 3])
        vertices[:len(xy),:] = xy
        vertices[len(xy):len(xy)+len(yz),:] = yz
        vertices[len(xy)+len(yz):,:] = xz
        vertices = vertices.dot(cell)
        vertices += center

        return cell, vertices, center

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
        rmax = max([self.box.width, self.box.height, self.box.length])
        rmin = min([self.box.width, self.box.height, self.box.length])
        if rmax/rmin > 3 and rmax >12:
            self.radius = rmin

    def get_symbols(self):
        self.symbols = [specie.name for specie in self.mol.species]

    def get_tols_matrix(self, mol2=None, tm=None):
        """
        Compute the 2D tolerance matrix between the current and other molecules

        Args:
            mol2: the 2nd pyxtal_molecule object
            tm: tolerance class

        Returns: 
            a 2D matrix which is used internally for distance checking.
        """
        if tm is None:
            tm = self.tm

        numbers1 = self.mol.atomic_numbers
        if mol2 is None:
            numbers2 = self.mol.atomic_numbers
        else:
            numbers2 = mol2.mol.atomic_numbers

        tols = np.zeros((len(numbers1), len(numbers2)))
        for i1, number1 in enumerate(numbers1):
            for i2, number2 in enumerate(numbers2):
                tols[i1, i2] = tm.get_tol(number1, number2)
                # allow hydrogen bond
                if [number1, number2] in [[1,7], [1,8], [1,9], [7,1], [8,1], [9,1]]:
                    tols[i1, i2] *= 0.9

        if len(self.mol)==1:
            tols *= 0.8 # if only one atom, reduce the tolerance
        return tols

    def get_coefs_matrix(self, mol2=None):
        """
        Get the Atom-Atom potential parameters
        E = A*exp(-B*r) - C*R^(-6)
        according to Gavezotti, Acc. Chem. Res., 27, 1994

        Args:
            mol2: the 2nd pyxtal_molecule object
            tm: tolerance class

        Returns: 
            a 3D matrix for computing the intermolecular energy
        """
        numbers1 = self.mol.atomic_numbers
        if mol2 is None:
            numbers2 = self.mol.atomic_numbers
        else:
            numbers2 = mol2.mol.atomic_numbers

        coefs = np.zeros([len(numbers1), len(numbers2), 3])
        for i1, number1 in enumerate(numbers1):
            for i2, number2 in enumerate(numbers2):
                if [number1, number2] in [[1,1]]:              #H-H
                    coefs[i1, i2, :] = [5774, 4.01, 26.1]
                elif [number1, number2] in [[1,6], [6,1]]:     #H-C
                    coefs[i1, i2, :] = [28870, 4.10, 113.]
                elif [number1, number2] in [[1,7], [7,1]]:     #H-N
                    coefs[i1, i2, :] = [54560, 4.52, 120.]
                elif [number1, number2] in [[1,8], [8,1]]:     #H-O
                    coefs[i1, i2, :] = [70610, 4.82, 105.]
                elif [number1, number2] in [[1,16],[16,1]]:    #H-S
                    coefs[i1, i2, :] = [64190, 4.03, 279.]
                elif [number1, number2] in [[1,17],[17,1]]:    #H-Cl
                    coefs[i1, i2, :] = [70020, 4.09, 279.]
                elif [number1, number2] in [[6,6]]:            #C-C
                    coefs[i1, i2, :] = [54050, 3.47, 578.]
                elif [number1, number2] in [[6,7], [7,6]]:     #C-N
                    coefs[i1, i2, :] = [117470, 3.86, 667.]
                elif [number1, number2] in [[6,8], [8,6]]:     #C-O
                    coefs[i1, i2, :] = [93950, 3.74, 641.]
                elif [number1, number2] in [[6,16], [16,6]]:   #C-S
                    coefs[i1, i2, :] = [126460, 3.41, 1504.]
                elif [number1, number2] in [[6,8], [8,6]]:     #C-Cl
                    coefs[i1, i2, :] = [93370, 3.52, 923.]
                elif [number1, number2] in [[7,7]]:            #N-N
                    coefs[i1, i2, :] = [87300, 3.65, 691.]
                elif [number1, number2] in [[7,8], [8,7]]:     #N-O
                    coefs[i1, i2, :] = [64190, 3.86, 364.]     
                elif [number1, number2] in [[8,8]]:            #O-O
                    coefs[i1, i2, :] = [46680, 3.74, 319.]     
                elif [number1, number2] in [[8,16], [16,8]]:   #O-S
                    coefs[i1, i2, :] = [110160, 3.63, 906.]     
                elif [number1, number2] in [[8,17], [17,8]]:   #O-Cl
                    coefs[i1, i2, :] = [80855, 3.63, 665.]   
                elif [number1, number2] in [[16, 16]]:         #S-S
                    coefs[i1, i2, :] = [259960, 3.52, 2571]  
                elif [number1, number2] in [[16,17], [17,16]]: #S-S
                    coefs[i1, i2, :] = [1800000, 3.52, 2000]   #made  
                elif [number1, number2] in [[17,17]]:          #Cl-Cl
                    coefs[i1, i2, :] = [140050, 3.52, 1385]  
                #todo add H-bonds
                #HB-O(amides):   3607810, 7.78, 238
                #HB-O(acids):    6313670, 8.75, 205
                #HB-O(alcohols): 4509750, 7.78, 298
                #HB-N(-N-H-N):   7215600, 7.78, 476
                #HB-N(-NH2-N):   1803920, 7.37, 165
                else:
                    return None
        return coefs

    def show(self):
        """
        show the molecule
        """
        from pyxtal.viz import display_molecules
        return display_molecules([self.mol])

    def show_box(self, center=np.zeros(3), orientation=None):
        """
        show the molecule
        """
        from pyxtal.viz import display_molecule
        return display_molecules(self.mol)


    def rdkit_mol(self, N_confs=1):
        """
        initialize the mol xyz and torsion list
        """
        from rdkit import Chem

        mol = Chem.MolFromMolBlock(self.rdkit_mb, removeHs=False)
        if N_confs > 1:
            conf = mol.GetConformer(0)
            for i in range(N_confs-1):
                mol.AddConformer(conf, True)
        return mol


    def rdkit_mol_init(self, smile, fix, torsions):
        """
        initialize the mol xyz and torsion list

        Args:
            smile: smile string
            fix: whether or not fix the seed
            torsions: None or list
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdMolTransforms as rdmt

        if smile not in ["Cl-", "F-", "Br-", "I-", "Li+", "Na+"]:
            torsionlist = find_id_from_smile(smile)
            mol = Chem.MolFromSmiles(smile)
            mol = Chem.AddHs(mol)
            symbols = []
            for id in range(mol.GetNumAtoms()):
                symbols.append(mol.GetAtomWithIdx(id).GetSymbol())
            if len(smile) > 100: #a tmp fix for KEKULN10
                AllChem.EmbedMultipleConfs(mol, numConfs=1, randomSeed=3)
                cid = 0
            else:
                ps = AllChem.ETKDGv3()
                ps.randomSeed=self.seed
                AllChem.EmbedMultipleConfs(mol, 1, ps)
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMultipleConfs(mol, 3, ps)

            res = AllChem.MMFFOptimizeMoleculeConfs(mol)
            engs = [c[1] for c in res]
            cid = engs.index(min(engs))
            self.rdkit_mb = Chem.MolToMolBlock(mol)
            self.energy = engs[cid]
            ref_conf = mol.GetConformer(cid) #always the reference molecule

            if fix or torsions is not None or len(torsionlist)==0:
                conf = ref_conf
            else:
                AllChem.EmbedMultipleConfs(mol, 
                                           numConfs=max([1, 4*len(torsionlist)]), 
                                           maxAttempts=200, 
                                           useRandomCoords=True, 
                                           pruneRmsThresh=0.5)
                N_confs = mol.GetNumConformers()
                conf_id = choice(range(N_confs))
                conf = mol.GetConformer(conf_id)
                #xyz = conf.GetPositions()
                #res = AllChem.MMFFOptimizeMoleculeConfs(mol)
                #print("Eng", res[conf_id])

            # set tosion angles from random or pre-defined values
            if torsions is not None:
                xyz = self.set_torsion_angles(conf, torsions, torsionlist=torsionlist)
            else:
                xyz = conf.GetPositions()
                xyz -= self.get_center(xyz)
        else:
            #single atom cation or anions
            symbols = ["Cl"]
            xyz = np.zeros([1,3])
            torsionlist = []
        return symbols, xyz, torsionlist

    def perturb_torsion(self, xyz):
        """
        slightly perturb the torsion
        """
        angs = self.get_torsion_angles(xyz, self.torsionlist)
        angs *= (1+0.1*np.random.uniform(-1., 1., len(angs))) 
        xyz = self.set_torsion_angles(conf, angs, torsionlist=self.torsionlist)
        xyz -= self.get_center(xyz)
        return xyz

    def align(self, conf, reflect=False, torsionlist=None):
        """
        Align the molecule and return the xyz
        The default CanonicalizeConformer function may also include inversion
        """
        from rdkit.Chem import rdMolTransforms as rdmt

        # Rotation
        if len(self.smile) > 1: 
            trans = rdmt.ComputeCanonicalTransform(conf)
            if abs(abs(np.linalg.det(trans))-1.0)>1e-1: 
                print("Bug in trans", np.linalg.det(trans))
                import sys; sys.exit()
            elif np.linalg.det(trans[:3,:3]) < 0: 
                trans[:3,:3] *= -1

            # add reflection if needed    
            if reflect: trans[:3,:3] *= -1
            rdmt.TransformConformer(conf, trans)

        # Translation
        pt = rdmt.ComputeCentroid(conf)
        center = np.array([pt.x, pt.y, pt.z])
        xyz = conf.GetPositions() - center

        # adjust cases like H2O
        if len(self.smile) == 1:
            xyz -= np.mean(xyz, axis=0)
            A = get_inertia_tensor(xyz)
            P = np.linalg.eigh(A)[1]
            if np.linalg.det(P) < 0:
                P[0] *= -1
            xyz = np.dot(xyz, P)
        return xyz

    def get_center(self, xyz, geometry=False):
        """
        get the molecular center for a transformed xyz
        """
        if geometry or self.smile is None:
            return np.mean(xyz, axis=0)
        else:
            if self.smile in ["Cl-", "F-", "Br-", "I-", "Li+", "Na+"]:
                return xyz[0]
            else:
                if len(self.smile) == 1:
                    #return xyz[0]
                    return np.mean(xyz, axis=0)
                else:
                    # from rdkit
                    from rdkit.Geometry import Point3D
                    from rdkit.Chem import rdMolTransforms as rdmt
                    
                    conf = self.rdkit_mol().GetConformer(0)
                    for i in range(len(xyz)):
                        x, y, z = xyz[i]
                        conf.SetAtomPosition(i, Point3D(x,y,z))
                    pt = rdmt.ComputeCentroid(conf)

                    return np.array([pt.x, pt.y, pt.z])

    def get_principle_axes(self, xyz, rdmt=True):
        """
        get the principle axis for a rotated xyz
        """
        if self.smile is None or len(self.smile)==1 or not rdmt: 
            Inertia = get_inertia_tensor(xyz)
            _, matrix = np.linalg.eigh(Inertia)
            return matrix
        else:
            from rdkit.Geometry import Point3D
            from rdkit.Chem import rdMolTransforms as rdmt

            conf1 = self.rdkit_mol().GetConformer(0)
            for i in range(len(self.mol)):
                x,y,z = xyz[i]
                conf1.SetAtomPosition(i,Point3D(x,y,z))

            return rdmt.ComputePrincipalAxesAndMoments(conf1)[0]

    def get_torsion_angles(self, xyz, torsionlist=None):
        """
        get the torsion angles
        """
        from rdkit.Geometry import Point3D
        from rdkit.Chem import rdMolTransforms as rdmt

        if torsionlist is None: torsionlist=self.torsionlist

        angs = []
        if len(torsionlist) > 0:
            conf = self.rdkit_mol().GetConformer(0)
            for i in range(len(xyz)):
                x,y,z = xyz[i]
                conf.SetAtomPosition(i,Point3D(x,y,z))

            for torsion in torsionlist:
                (i, j, k, l) = torsion
                angs.append(rdmt.GetDihedralDeg(conf, i, j, k, l))
        return angs
    
    def set_torsion_angles(self, conf, angles, reflect=False, torsionlist=None):
        """
        reset the torsion angles and update molecular xyz
        """
        from rdkit.Chem import rdMolTransforms as rdmt

        if torsionlist is None: torsionlist=self.torsionlist
        for id, torsion in enumerate(torsionlist):
            (i, j, k, l) = torsion
            rdmt.SetDihedralDeg(conf, i, j, k, l, angles[id])
        
        xyz = self.align(conf, reflect)
        return xyz

    def relax(self, xyz, align=False):
        """
        Relax the input xyz coordinates with rdit MMFF

        Args:
            xyz: 3D coordinates
            align: whether or not align the xyz

        Returns:
            xyz: new xyz
            eng: energy value
        """

        from rdkit.Chem import AllChem
        from rdkit.Geometry import Point3D

        mol = self.rdkit_mol(1)
        conf0 = mol.GetConformer(0)
        # reset the xyz
        for i in range(len(self.mol)):
            x,y,z = xyz[i]
            conf0.SetAtomPosition(i,Point3D(x,y,z))
        res = AllChem.MMFFOptimizeMoleculeConfs(mol)
        if align:
            xyz = self.align(conf0)
        else:   
            xyz = mol.GetConformer(0).GetPositions()
        return xyz, res[0][1]


    def get_rmsd(self, xyz):
        """
        Compute the rmsd with another 3D xyz coordinates

        Args:
            xyz: 3D coordinates

        Returns:
            rmsd:
            transition matrix:
            match:
        """

        from rdkit import Chem
        from rdkit.Geometry import Point3D
        from rdkit.Chem import rdMolAlign, RemoveHs

        mol = self.rdkit_mol(3)
        # 3 conformers for comparison
        conf0 = mol.GetConformer(0)
        conf1 = mol.GetConformer(1)  #reference+reflection
        conf2 = mol.GetConformer(2)  #trial xyz
        angs = self.get_torsion_angles(xyz)
        xyz0 = self.set_torsion_angles(conf0, angs) #conf0 with aligned
        xyz1 = self.set_torsion_angles(conf0, angs, True) #conf0 with aligned+reflect

        # reset the xyz
        for i in range(len(self.mol)):
            x0,y0,z0 = xyz0[i]
            x1,y1,z1 = xyz1[i]
            x,y,z = xyz[i]
            conf0.SetAtomPosition(i,Point3D(x0,y0,z0))
            conf1.SetAtomPosition(i,Point3D(x1,y1,z1))
            conf2.SetAtomPosition(i,Point3D(x,y,z))

        mol = RemoveHs(mol)
        rmsd1, trans1 = rdMolAlign.GetAlignmentTransform(mol, mol, 2, 0)
        rmsd2, trans2 = rdMolAlign.GetAlignmentTransform(mol, mol, 2, 1)
        #rdmolfiles.MolToXYZFile(mol, '1.xyz', 0)
        #rdmolfiles.MolToXYZFile(mol, '2.xyz', 1)
        #rdmolfiles.MolToXYZFile(mol, '3.xyz', 2)
       
        if rmsd1 <= rmsd2:
            return rmsd1, trans1, True
        else:
            return rmsd2, trans2, False

    def get_orientation(self, xyz, rtol=0.15):
        """
        For the given xyz, compute the orientation

        Args:
            xyz: molecular xyz
        """
        
        xyz -= self.get_center(xyz)

        if len(self.smile) > 1: # not in ["O", "o"]:
            rmsd, trans, reflect = self.get_rmsd(xyz)
            tol = rtol*len(xyz)
            
            if rmsd < tol:
                trans = trans[:3,:3].T
                r = Rotation.from_matrix(trans)
                return r.as_euler('zxy', degrees=True), rmsd, reflect
            else:
                msg = "Problem in conformer\n"
                msg += "{:5.2f} {:5.2f}\n".format(rmsd1, rmsd2)
                if len(self.torsionlist) > 0:
                    msg += str(self.get_torsion_angles(xyz)) + '\n'  
                    msg += str(self.get_torsion_angles(xyz0)) + '\n'  
                    msg += str(self.get_torsion_angles(xyz1)) + '\n'  
                raise ConformerError(msg)
        else:
            # the orientation of CH4, NH3, H2O
            ref = np.array([[-0.00111384,  0.36313718,  0.        ],
                            [-0.82498189, -0.18196256,  0.        ], 
                            [ 0.82609573, -0.18117463,  0.        ]])
            Inertia = get_inertia_tensor(xyz) 
            _, matrix = np.linalg.eigh(Inertia)

            ref0 = np.dot(xyz, matrix)
            # identify the rotation matrix
            libs = np.array([
                            [[1,1,1]],
                            [[-1,1,1]],
                            [[1,-1,1]],
                            [[1,1,-1]],
                            [[-1,-1,1]],
                            [[1,-1,-1]],
                            [[-1,1,-1]],
                            [[-1,-1,-1]]
                            ])
            dists = np.zeros(8)
            for i, lib in enumerate(libs):
                matrix0 = matrix*np.repeat(lib, 3, axis=0)
                res = np.dot(ref, np.linalg.inv(matrix0))
                dists[i] = np.sum((res-xyz)**2)
                #print(i, res)
            id = np.argmin(dists)
            matrix = matrix*np.repeat(libs[id], 3, axis=0)

            r = Rotation.from_matrix(np.linalg.inv(matrix).T)
            ang = r.as_euler('zxy', degrees=True)
            return ang, 0, False

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

    def get_symmetry(self, xyz=None, symmetrize=False, rtol=0.30):
        """
        Set the molecule's point symmetry.
            - pga: pymatgen.symmetry.analyzer.PointGroupAnalyzer object
            - pg: pyxtal.symmetry.Group object
            - symops: a list of SymmOp objects

        Args:
            symmetrize: boolean, whether or not symmetrize the coordinates
        """
        if xyz is None:
            mol = deepcopy(self.mol)
        else:
            mol = Molecule(self.symbols, xyz)

        if self.smile is not None:
            mol.remove_species('H')
            mol._spin_multiplicity = None #don't check spin

        if symmetrize:
            pga = PointGroupAnalyzer(mol, rtol, eigen_tolerance=1e-3)
            mol = pga.symmetrize_molecule()["sym_mol"]
        pga = PointGroupAnalyzer(mol, rtol, eigen_tolerance=1e-3)
        self.mol_no_h = mol

        #print(mol.to(fmt='xyz'), pga.sch_symbol)
        # For single atoms, we cannot represent the point group using a list of operations
        if len(mol) == 1:
            symm_m = []
            symbol = 'C1'
        else:
            symbol = pga.sch_symbol
            pg = pga.get_pointgroup()
            symm_m = [op for op in pg]
    
            if "*" in symbol: # linear molecules
                symbol = symbol.replace('*','6')
                # Add 12-fold  and reflections in place of ininitesimal rotation
                for i, axis in enumerate(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])):
                    # op = SymmOp.from_rotation_and_translation(aa2matrix(axis, np.pi/6), [0,0,0])
                    m1 = Rotation.from_rotvec(np.pi / 6 * axis).as_matrix()
                    op = SymmOp.from_rotation_and_translation(m1, [0, 0, 0])
                    if pga.is_valid_op(op):
                        symm_m.append(op)
                        # Any molecule with infinitesimal symmetry is linear;
                        # Thus, it possess mirror symmetry for any axis perpendicular
                        # To the rotational axis. pymatgen does not add this symmetry
                        # for all linear molecules - for example, hydrogen
                        if i == 0:
                            symm_m.append(SymmOp.from_xyz_string("x,-y,z"))
                            symm_m.append(SymmOp.from_xyz_string("x,y,-z"))
                            #r = SymmOp.from_xyz_string("-x,y,-z")
                        elif i == 1:
                            symm_m.append(SymmOp.from_xyz_string("-x,y,z"))
                            symm_m.append(SymmOp.from_xyz_string("x,y,-z"))
                            #r = SymmOp.from_xyz_string("-x,-y,z")
                        elif i == 2:
                            symm_m.append(SymmOp.from_xyz_string("-x,y,z"))
                            symm_m.append(SymmOp.from_xyz_string("x,-y,z"))
                            #r = SymmOp.from_xyz_string("x,-y,-z")
                        # Generate a full list of SymmOps for the molecule's pointgroup
                        symm_m = generate_full_symmops(symm_m, 1e-3)
                        break
        self.symops = symm_m
        self.pga = pga
        if symbol == 'S2': symbol = 'Ci'
        self.pg = Group(symbol, dim=0)

    def get_orientations_in_wps(self, wps=None, rtol=1e-2):
        """
        Compute the valid orientations from a given Wyckoff site symmetry.
    
        Args:
            wp: a pyxtal.symmetry.Wyckoff_position object
        Returns:
            a list of operations.Orientation objects 
        """
        if wps is None:
            return None, True
        else:
            valid = False
            valid_oris = []
            for wp in wps:
                allowed = self.get_orientations_in_wp(wp, rtol)
                if len(allowed) > 0:
                    valid = True
                valid_oris.append(allowed)
            return valid_oris, valid

    def get_orientations_in_wp(self, wp, rtol=1e-2):
        """
        Compute the valid orientations from a given Wyckoff site symmetry.
    
        Args:
            wp: a pyxtal.symmetry.Wyckoff_position object
        Returns:
            a list of pyxtal.molecule.Orientation objects 
        """
        # For single atoms, there are no constraints
        if len(self.mol) == 1 or wp.index == 0:
            return [Orientation([[1, 0, 0], [0, 1, 0], [0, 0, 1]], degrees=2)]
        # C1 molecule cannot take specical position
        elif wp.index > 1 and self.pga.sch_symbol == 'C1':
            return []

        symm_w = wp.get_site_symm_wo_translation() #symmetry without translation 
        # molecule has fewer symops
        if len(self.pg[0]) < len(symm_w):
            return []
    
        symm_m = self.symops 
        wyckoffs = wp.ops
        opa_m = []
        for op_m in symm_m:
            opa = OperationAnalyzer(op_m)
            opa_m.append(opa)
        
        # Store OperationAnalyzer objects for each Wyckoff symmetry SymmOp
        opa_w = []
        for op_w in symm_w:
            opa_w.append(OperationAnalyzer(op_w))
    
        # Check for constraints from the Wyckoff symmetry
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
                            if cond1:  # or cond2 is True:
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
                        if self.pga.is_valid_op(old_op):
                            list_i.remove(j)
                            list_j.remove(j)
        orientations_new = []
        for i in list_i:
            orientations_new.append(orientations[i])
    
        #Check each of the found orientations for consistency with the Wyckoff pos.
        #If consistent, put into an array of valid orientations
        #print("======", orientations_new)
        allowed = []
        for o in orientations_new:
            op = o.get_op()
            mo = deepcopy(self.mol_no_h)
            mo.apply_operation(op)
            #print(mo)
            if is_compatible_symmetry(mo, wp):
                allowed.append(o)
        return allowed

    def get_energy(self, xyz1, xyz2):
        """
        Get packing energy between two neighboring molecules
        """
        dists = cdist(xyz1-xyz2)



class Box:
    """
    Class for storing the binding box for a molecule. 

    Args:
        dims: [length, width, height]
    """

    def __init__(self, dims):
        self.length = dims[0] #float(abs(maxy - miny))
        self.width = dims[1] #float(abs(maxx - minx))
        self.height = dims[2] #float(abs(maxz - minz))
        self.volume = self.width * self.length * self.height

    def __str__(self):
        strs = "l: {:6.2f}, w: {:6.2f}, d: {:6.2f}".format(self.length, self.width, self.height)
        return strs

    def operate(self, rot=np.eye(3), center=np.zeros(3)):
        """
        Perform operation on the box:

        Args:
            rot: 3*3 rotation matrix
            center: center position
        """
        raise NotImplementedError

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


def get_inertia_tensor(coords, weights=None):
    """
    Calculate the symmetric inertia tensor for a Molecule
    the principal axes of symmetry.

    Args:
        coords: [N, 3] array of coordinates

    Returns:
        a 3x3 numpy array representing the inertia tensor
    """
    if weights is None: weights = np.ones(len(coords))
    coords -= np.mean(coords, axis=0)
    Inertia = np.zeros([3,3])
    Inertia[0,0] = np.sum(weights*coords[:,1]**2 + weights*coords[:,2]**2)
    Inertia[1,1] = np.sum(weights*coords[:,0]**2 + weights*coords[:,2]**2)
    Inertia[2,2] = np.sum(weights*coords[:,0]**2 + weights*coords[:,1]**2)
    Inertia[0,1] = Inertia[1,0] = -np.sum(weights*coords[:,0]*coords[:,1])
    Inertia[0,2] = Inertia[2,0] = -np.sum(weights*coords[:,0]*coords[:,2])
    Inertia[1,2] = Inertia[2,1] = -np.sum(weights*coords[:,1]*coords[:,2])

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
    if np.linalg.det(P) < 0: P[:,0] *= -1
    coords = np.dot(coords, P)
    return Molecule(numbers, coords), P


def is_compatible_symmetry(mol, wp):
    """
    Tests if a molecule meets the symmetry requirements of a Wyckoff position

    Args:
        mol: a pymatgen Molecule object.
        wp: a pyxtal.symmetry.Wyckoff_position object
    """
    # For single atoms, there are no constraints
    if len(mol) == 1 or wp.index == 0:
        return True
    pga = PointGroupAnalyzer(mol)
    for op in wp.get_site_symm_wo_translation(): #symmetry without translation
        #print("XXXX", pga.is_valid_op(op), op.as_xyz_string())
        if not pga.is_valid_op(op):
            return False
    return True


def make_graph(mol, tol=0.2):
    """
    make graph object for the input molecule
    """
    #print("making graphs")
    G = nx.Graph()
    names = {}
    for i, site in enumerate(mol._sites):
        names[i] = site.specie.value
        if names[i] not in ["C", "H", "O", "N", "S", "P", "Si", "F", "Cl", "Br", "I"]:
            raise ValueError(name[i]+' is not supported')

    for i in range(len(mol)-1):
        site1 = mol.sites[i]
        for j in range(i+1, len(mol)):
            site2 = mol.sites[j]
            key = "{:s}-{:s}".format(names[i], names[j])
            if site1.distance(site2) < bonds[key]:
                G.add_edge(i,j)
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
