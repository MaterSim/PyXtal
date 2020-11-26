"""
This module handles reading and write crystal files.
"""
from pyxtal.constants import deg, logo
import numpy as np
from pyxtal.symmetry import Group

def write_cif(struc, filename=None, header="", permission='w', sym_num=None, style='mp'):
    """
    Export the structure in cif format
    The default setting for _atom_site follows the materials project cif 

    Args:
        struc: pyxtal structure object
        filename: path of the structure file 
        header: additional information
        permission: write('w') or append('a+') to the given file
        sym_num: the number of symmetry operations, None means writing all symops
        style: "icsd" or "mp" (used in pymatgen)

    """
    if sym_num is None:
        l_type = struc.group.lattice_type
        symbol = struc.group.symbol
        number = struc.group.number
        G1 = struc.group.Wyckoff_positions[0]
    else: #P1 symmetry
        l_type = 'triclinic'
        symbol = 'P1'
        number = 1
        G1 = Group(1).Wyckoff_positions[0]

    if hasattr(struc, 'mol_sites'):
        sites = struc.mol_sites
        molecule = True
    else:
        sites = struc.atom_sites
        molecule = False

    change_set = False
    if number in [7, 14, 15]:
        if hasattr(struc, 'diag') and struc.diag:
            symbol = struc.group.alias 
            G1.diagonalize_symops()
            change_set = True
    
    lines = logo
    lines += 'data_' + header + '\n'
    if hasattr(struc, "energy"):
        lines += '#Energy: {:} eV/cell\n'.format(struc.energy/sum(struc.numMols))

    lines += "\n_symmetry_space_group_name_H-M '{:s}'\n".format(symbol)
    lines += '_symmetry_Int_Tables_number      {:>15d}\n'.format(number)
    lines += '_symmetry_cell_setting           {:>15s}\n'.format(l_type)

    a, b, c, alpha, beta, gamma = struc.lattice.get_para(degree=True)
    lines += '_cell_length_a        {:12.6f}\n'.format(a)
    lines += '_cell_length_b        {:12.6f}\n'.format(b)
    lines += '_cell_length_c        {:12.6f}\n'.format(c)
    lines += '_cell_angle_alpha     {:12.6f}\n'.format(alpha)
    lines += '_cell_angle_beta      {:12.6f}\n'.format(beta)
    lines += '_cell_angle_gamma     {:12.6f}\n'.format(gamma)

    lines += '\nloop_\n'
    lines += ' _symmetry_equiv_pos_site_id\n'
    lines += ' _symmetry_equiv_pos_as_xyz\n'

    if not change_set:
        wps = G1
    else:
        wps = sites[0].wp.ops
    for i, op in enumerate(wps):
        lines += "{:d} '{:s}'\n".format(i+1, op.as_xyz_string())

    lines += '\nloop_\n'
    lines += ' _atom_site_label\n'
    lines += ' _atom_site_type_symbol\n'
    lines += ' _atom_site_symmetry_multiplicity\n'
    if style == 'icsd':
        lines += ' _atom_site_Wyckoff_symbol\n'
    lines += ' _atom_site_fract_x\n'
    lines += ' _atom_site_fract_y\n'
    lines += ' _atom_site_fract_z\n'
    lines += ' _atom_site_occupancy\n'

    for site in sites:
        mul = site.wp.multiplicity
        letter = site.wp.letter
        if molecule:
            if sym_num is None:
                coords, species = site._get_coords_and_species(first=True)
            else:
                coords = None
                species = []
                for id in range(sym_num):
                    mol = site.get_mol_object(id)
                    tmp = mol.cart_coords.dot(site.lattice.inv_matrix)
                    if coords is None:
                        coords = tmp
                    else:
                        coords = np.append(coords, tmp, axis=0)
                    species.extend([s.value for s in mol.species])
                #coords, species = site._get_coords_and_species(ids=sym_num)
        else:
            coords, species = [site.position], [site.specie]
        for specie, coord in zip(species, coords):
            if style == 'mp':
                lines += '{:6s} {:6s} {:3d} {:12.6f}{:12.6f}{:12.6f} 1\n'.format(\
                    specie, specie, mul, *coord)
            else:
                lines += '{:6s} {:6s} {:3d} {:s} {:12.6f}{:12.6f}{:12.6f} 1\n'.format(\
                    specie, specie, mul, letter, *coord)
    lines +='#END\n\n'
    

    if filename is None:
        return lines
    else:
        with open(filename, permission) as f:
            f.write(lines)
        return

from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.bonds import CovalentBond
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pyxtal.wyckoff_site import mol_site, WP_merge
from pyxtal.molecule import pyxtal_molecule, Orientation, compare_mol_connectivity
from pyxtal.symmetry import Wyckoff_position, Group
from pyxtal.lattice import Lattice


class structure_from_ext():
    
    def __init__(self, struc, ref_mol=None, tol=0.2, relax_h=False):

        """
        extract the mol_site information from the give cif file 
        and reference molecule
    
        Args: 
            struc: cif/poscar file or a Pymatgen Structure object
            ref_mol: xyz file or a reference Pymatgen molecule object
            tol: scale factor for covalent bond distance
            relax_h: whether or not relax the position for hydrogen atoms in structure
        
    """
        if isinstance(ref_mol, str):
            ref_mol = Molecule.from_file(ref_mol)
        elif isinstance(ref_mol, Molecule):
            ref_mol = ref_mol
        else:
            print(type(ref_mol))
            raise NameError("reference molecule cannot be defined")

    
        if isinstance(struc, str):
            pmg_struc = Structure.from_file(struc)
        elif isinstance(struc, Structure):
            pmg_struc = struc
        else:
            print(type(struc))
            raise NameError("input structure cannot be intepretted")

        self.props = ref_mol.site_properties
        self.ref_mol = ref_mol.get_centered_molecule()
        self.tol = tol
        self.diag = False
        self.relax_h = relax_h

        sga = SpacegroupAnalyzer(pmg_struc)
        ops = sga.get_space_group_operations()
        self.wyc, perm = Wyckoff_position.from_symops(ops, sga.get_space_group_number())
        if self.wyc is not None:
            self.group = Group(self.wyc.number)
            if isinstance(perm, list):
                if perm != [0,1,2]:
                    lattice = Lattice.from_matrix(pmg_struc.lattice.matrix, self.group.lattice_type)
                    latt = lattice.swap_axis(ids=perm, random=False).get_matrix()
                    coor = pmg_struc.frac_coords[:, perm]
                    pmg_struc = Structure(latt, pmg_struc.atomic_numbers, coor)
            else:
                self.diag = True
                self.perm = perm

            coords, numbers = search_molecule_in_crystal(pmg_struc, self.tol)
            #coords -= np.mean(coords, axis=0)
            if self.relax_h:
                self.molecule = self.addh(Molecule(numbers, coords))
            else:
                self.molecule = Molecule(numbers, coords)
            self.pmg_struc = pmg_struc
            self.lattice = Lattice.from_matrix(pmg_struc.lattice.matrix, self.group.lattice_type)
            self.numMols = [len(self.wyc)]
        else:
            raise ValueError("Cannot find the space group matching the symmetry operation")

    def addh(self, mol):
        #if len(mol) < len(self.ref_mol):
        from pymatgen.io.babel import BabelMolAdaptor
        ad = BabelMolAdaptor(mol)
        ad.add_hydrogen()        
        ad.localopt()
        mol = ad.pymatgen_mol
        return mol


    def add_site_props(self, mo):
        if len(self.props) > 0:
            for key in self.props.keys():
                mo.add_site_property(key, self.props[key])
        return mo


    def make_mol_site(self, ref=False):
        if ref:
            mol = self.ref_mol
            ori = self.ori
        else:
            mol = self.molecule
            ori = Orientation(np.eye(3))
        mol = self.add_site_props(mol)
        pmol = pyxtal_molecule(mol, symmetrize=False)
        site = mol_site(pmol, self.position, ori, self.wyc, self.lattice, self.diag)
        return site


    def align(self):
        """
        compute the orientation wrt the reference molecule
        """
        try:
            from openbabel import pybel, openbabel
        except:
            import pybel, openbabel

        m1 = pybel.readstring('xyz', self.ref_mol.to('xyz'))
        m2 = pybel.readstring('xyz', self.molecule.to('xyz'))
        aligner = openbabel.OBAlign(True, False)
        aligner.SetRefMol(m1.OBMol)
        aligner.SetTargetMol(m2.OBMol)
        aligner.Align()
        print("RMSD: ", aligner.GetRMSD())
        rot=np.zeros([3,3])
        for i in range(3):
            for j in range(3):
                rot[i,j] = aligner.GetRotMatrix().Get(i,j)
        coord2 = self.molecule.cart_coords
        coord2 -= np.mean(coord2, axis=0)
        coord3 = rot.dot(coord2.T).T + np.mean(self.ref_mol.cart_coords, axis=0)
        self.mol_aligned = Molecule(self.ref_mol.atomic_numbers, coord3)
        self.ori = Orientation(rot)
   
    def match(self):
        """
        Check the two molecular graphs are isomorphic
        """
        match, mapping = compare_mol_connectivity(self.ref_mol, self.molecule)
        if not match:
            print(self.ref_mol.to("xyz"))
            print(self.molecule.to("xyz"))
            import pickle
            with open('wrong.pkl', "wb") as f:
                pickle.dump([self.ref_mol, self.molecule], f)

            return False
        else:
            # resort the atomic number for molecule 1
            order = [mapping[i] for i in range(len(self.ref_mol))]
            numbers = np.array(self.molecule.atomic_numbers)
            numbers = numbers[order].tolist()
            coords = self.molecule.cart_coords[order]
            position = np.mean(coords, axis=0).dot(self.lattice.inv_matrix)
            position -= np.floor(position)
            # check if molecule is on the special wyckoff position
            if len(self.pmg_struc)/len(self.molecule) < len(self.wyc):
                if self.diag:
                    #Transform it to the conventional representation
                    position = np.dot(self.perm, position).T
                position, wp, _ = WP_merge(position, self.lattice.matrix, self.wyc, 2.0)
                #print("After Mergey:---------------")
                #print(position)
                #print(wp)
                self.wyc = wp
            self.position = position
            self.molecule = Molecule(numbers, coords-np.mean(coords, axis=0))
            #self.align()
            return True

    def show(self, overlay=True):
        from pyxtal.viz import display_molecules
        if overlay:
            return display_molecules([self.ref_mol, self.mol_aligned])
        else:
            return display_molecules([self.ref_mol, self.molecule])


def search_molecule_in_crystal(struc, tol=0.2, keep_order=False, absolute=True):
    """
    This is a function to perform a search to find the molecule
    in a Pymatgen crystal structure

    Args:
        struc: Pymatgen Structure
        keep_order: whether or not use the orignal sequence
        absolute: whether or not output absolute coordindates

    Returns:
        coords: fractional coordinates
        numbers: atomic numbers
    """
    def check_one_layer(struc, sites0, visited):
        new_members = []
        for site0 in sites0:
            sites_add, visited = check_one_site(struc, site0, visited)
            new_members.extend(sites_add)
        return new_members, visited
    
    def check_one_site(struc, site0, visited):
        neigh_sites = struc.get_neighbors(site0, 3.0)
        ids = [m.index for m in visited]
        sites_add = []
        ids_add = []

        for site1 in neigh_sites:
            if (site1.index not in ids+ids_add):
                if CovalentBond.is_bonded(site0, site1, tol):
                    sites_add.append(site1)
                    ids_add.append(site1.index)
        if len(sites_add) > 0:
            visited.extend(sites_add)

        return sites_add, visited

    first_site = struc.sites[0]
    first_site.index = 0 #assign the index
    visited = [first_site] 
    ids = [0]

    n_iter, max_iter = 0, len(struc)
    while n_iter < max_iter:
        if n_iter == 0:
            new_sites, visited = check_one_site(struc, first_site, visited)
        else:
            new_sites, visited = check_one_layer(struc, new_sites, visited)
        n_iter += 1
        if len(new_sites)==0:
            break
    
    coords = [s.coords for s in visited]
    coords = np.array(coords)
    numbers = [s.specie.number for s in visited]
    
    if keep_order:
        ids = [s.index for s in visited]
        seq = np.argsort(ids)
        coords = coords[seq]
        numbers = numbers[seq]

    if not absolute:
        coords = coords.dot(struc.lattice.inv_matrix)
    return coords, numbers

#seed = structure_from_cif("254385.cif", "1.xyz")
#if seed.match():
#    print(seed.pmg_struc)

