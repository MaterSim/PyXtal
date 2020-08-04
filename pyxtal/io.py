"""
CIF in PyXtal format
"""
from pyxtal.constants import deg, logo
import numpy as np

def write_cif(struc, filename=None, header="", permission='w'):
    """
    Export the structure in cif format

    Args:
        struc: pyxtal structure object
        filename: path of the structure file 
        header: additional information
        permission: write('w') or append('a+') to the given file
    
    """

    l_type = struc.group.lattice_type
    symbol = struc.group.symbol
    number = struc.group.number
    change_set = False
    G1 = struc.group.Wyckoff_positions[0]
    if hasattr(struc, 'mol_sites'):
        sites = struc.mol_sites
        molecule = True
    else:
        sites = struc.atom_sites
        molecule = False

    if l_type == 'monoclinic':
        if G1 != sites[0].wp.generators:
            symbol = symbol.replace('c','n')
            change_set = True
    
    lines = logo
    lines += 'data_' + header + '\n'
    if hasattr(struc, "energy"):
        lines += '#Energy: {:} eV/cell\n'.format(struc.energy/sum(struc.numMols))

    lines += "\n_symmetry_space_group_name_H-M '{:s}'\n".format(symbol)
    lines += '_symmetry_Int_Tables_number      {:>15d}\n'.format(number)
    lines += '_symmetry_cell_setting           {:>15s}\n'.format(l_type)
    lines += '_cell_length_a        {:12.6f}\n'.format(struc.lattice.a)
    lines += '_cell_length_b        {:12.6f}\n'.format(struc.lattice.b)
    lines += '_cell_length_c        {:12.6f}\n'.format(struc.lattice.c)
    lines += '_cell_angle_alpha     {:12.6f}\n'.format(deg*struc.lattice.alpha)
    lines += '_cell_angle_beta      {:12.6f}\n'.format(deg*struc.lattice.beta)
    lines += '_cell_angle_gamma     {:12.6f}\n'.format(deg*struc.lattice.gamma)

    lines += '\nloop_\n'
    lines += ' _symmetry_equiv_pos_site_id\n'
    lines += ' _symmetry_equiv_pos_as_xyz\n'
    if not change_set:
        wps = G1
    else:
        wps = sites[0].wp.generators
    for i, op in enumerate(wps):
        lines += "{:d} '{:s}'\n".format(i+1, op.as_xyz_string())
    lines += '\nloop_\n'
    lines += ' _atom_site_label\n'
    lines += ' _atom_site_fract_x\n'
    lines += ' _atom_site_fract_y\n'
    lines += ' _atom_site_fract_z\n'

    for site in sites:
        if molecule:
            coords, species = site._get_coords_and_species(first=True)
        else:
            coords, species = [site.position], [site.specie]
        for specie, coord in zip(species, coords):
            lines += '{:6s}  {:12.6f}{:12.6f}{:12.6f}\n'.format(specie, *coord)
    lines +='#END\n\n'
    

    if filename is None:
        return lines
    else:
        with open(filename, permission) as f:
            f.write(lines)
        return

from pymatgen.io.cif import CifParser
from molmod.molecules import Molecule as mmm
from pymatgen import Molecule
from pymatgen.core.bonds import CovalentBond
import py3Dmol
from pyxtal.wyckoff_site import mol_site
from pyxtal.molecule import pyxtal_molecule, Orientation, compare_mol_connectivity
#from pyxtal.lattice import Lattice
from pyxtal.symmetry import Wyckoff_position

class structure_from_cif():
    
    def __init__(self, cif_file, ref_mol=None, tol=0.2):

        """
        extract the mol_site information from the give cif file 
        and reference molecule
    
        Args: 
            cif_file: file to store the experimental
            ref_mol: xyz file or a reference Pymatgen molecule object
            tol: scale factor for covalent bond distance
        
    """
        if isinstance(ref_mol, str):
            self.ref_mol = Molecule.from_file(ref_mol)
        elif isinstance(ref_mol, Molecule):
            self.ref_mol = ref_mol
        else:
            print(type(ref_mol))
            raise NameError("reference molecule cannot be defined")
        self.tol = tol
        
        self.cif = CifParser(cif_file)
        self.pmg_struc = self.cif.get_structures()[0]
        self.wyc = Wyckoff_position.from_symops(self.cif.symmetry_operations)  
        coords, numbers = search_molecule_in_crystal(self.pmg_struc, self.tol)
        self.molecule = Molecule(numbers, coords)
    
    def make_mol_site(self):
        mol = pyxtal_molecule(self.molecule)
        radius = mol.radius
        tols_matrix = mol.tols_matrix
        site = mol_site(self.molecule, 
                        self.position, 
                        Orientation(np.eye(3)),
                        self.wyc, 
                        self.lattice, #.matrix,
                        tols_matrix, 
                        radius,
                        rotate_ref = False,
                        )
        return site

   
    def match(self):
        """
        Check the two molecular graphs are isomorphic
        """
        self.G_ref = make_graph(self.ref_mol)
        self.G_cif = make_graph(self.molecule)
        fun = lambda n1, n2: n1['name'] == n2['name']
        match, mapping = compare_mol_connectivity(self.ref_mol, self.molecule)
        if not match:
            return False
        else:
            # resort the atomic number for molecule 1
            order = [mapping[i] for i in range(len(self.ref_mol))]
            numbers = np.array(self.molecule.atomic_numbers)
            numbers = numbers[order].tolist()
            coords = self.molecule.cart_coords[order]
            # Ideally, we want to 
            self.lattice = self.pmg_struc.lattice.matrix
            #Todo: figure out the molecular rotation
            #Needs to fix the lattice
            position = np.mean(coords, axis=0).dot(np.linalg.inv(self.lattice))
            position -= np.floor(position)
            position[0] = 0 #parse symmetry to adjust the position
            self.position = position
            self.molecule = Molecule(numbers, coords-np.mean(coords, axis=0))
 
            return True

    def show_overlay(self):
        from molmod import angstrom
        mol1 = mmm(self.ref_mol.atomic_numbers, self.ref_mol.cart_coords/angstrom)
        mol2 = mmm(self.molecule.atomic_numbers, self.molecule.cart_coords/angstrom)
        trans, coord, rmsd = mol1.rmsd(mol2)
        trans_mol = Molecule(self.molecule.atomic_numbers, (coord-5)*angstrom)
        #print(rmsd)

        view = py3Dmol.view()
        view.setStyle({'stick':{'colorscheme':'blueCarbon'}})
        view.addModel(self.molecule.to(fmt='xyz'), 'xyz')
        view.setStyle({'stick':{'colorscheme':'greenCarbon'}})
        view.addModel(self.ref_mol.to(fmt='xyz'), 'xyz')
        view.addModel(trans_mol.to(fmt='xyz'), 'xyz')
        view.setStyle({'stick':{'colorscheme':'redCarbon'}})
        return view.zoomTo()


def search_molecule_in_crystal(struc, tol=0.2, keep_order=False, absolute=True):
    """
    This is a function to perform a search to find the molecule
    in a Pymatgen crystal structure

    Args:
    struc: Pymatgen Structure
    keep_order: whether or not use the orignal sequence
    absolute: whether or not output absolute coordindates

    Returns:
    coords: frac
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
    numbers = [s.specie.number for s in visited]
    
    if keep_order:
        ids = [s.index for s in visited]
        seq = np.argsort(ids)
        coords = coords[seq]
        numbers = numbers[seq]

    if not absolute:
        coords = coords.dot(struc.lattice.inv)
    return coords, numbers



seed = structure_from_cif("254385.cif", "1.xyz")
print(seed.molecule)




