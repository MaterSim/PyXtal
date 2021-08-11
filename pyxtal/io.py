"""
This module handles reading and write crystal files.
"""
import numpy as np
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.bonds import CovalentBond
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pyxtal.wyckoff_site import atom_site, mol_site
from pyxtal.molecule import pyxtal_molecule, Orientation, compare_mol_connectivity
from pyxtal.symmetry import Wyckoff_position, Group
from pyxtal.lattice import Lattice
from pyxtal.util import get_symmetrized_pmg
from pyxtal.constants import deg, logo
from pyxtal.msg import ReadSeedError
from pkg_resources import resource_filename
from monty.serialization import loadfn

bonds = loadfn(resource_filename("pyxtal", "database/bonds.json"))

def write_cif(struc, filename=None, header="", permission='w', sym_num=None, style='mp'):
    """
    Export the structure in cif format
    The default setting for _atom_site follows the materials project cif 

    Args:
        struc: pyxtal structure object
        filename: path of the structure file 
        header: additional information
        permission: write(`w`) or append(`a+`) to the given file
        sym_num: the number of symmetry operations, None means writing all symops
        style: `icsd` or `mp` (used in pymatgen)

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
        if struc.molecular:
            eng = struc.energy/sum(struc.numMols)
        else:
            eng = struc.energy/sum(struc.numIons)
        lines += '#Energy: {:} eV/cell\n'.format(eng)

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
    lines += '_cell_volume          {:12.6f}\n'.format(struc.lattice.volume)
    #if struc.molecular:
    #    lines += '_cell_formula_units_Z     {:d}\n'.format(sum(struc.numMols))
    #else:
    #    lines += '_cell_formula_units_Z     {:d}\n'.format(sum(struc.numIons))

    lines += '\nloop_\n'
    lines += ' _symmetry_equiv_pos_site_id\n'
    lines += ' _symmetry_equiv_pos_as_xyz\n'

    if not change_set:
    #if change_set:
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
            lines += '{:6s} {:6s} {:3d} '.format(specie, specie, mul)
            if style != 'mp':
                lines += '{:s} '.format(letter)
            lines += '{:12.6f}{:12.6f}{:12.6f} 1\n'.format(*coord)
    lines +='#END\n\n'

    if filename is None:
        return lines
    else:
        with open(filename, permission) as f:
            f.write(lines)
        return

def read_cif(filename):
    """
    read the cif, mainly for pyxtal cif output
    Be cautious in using it to read other cif files

    Args:
        filename: path of the structure file 

    Return:
        pyxtal structure
    """
    species = []
    coords = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith('_symmetry_Int_Tables_number'):
                sg = int(line.split()[-1])
            elif line.startswith('_cell_length_a'):
                a = float(lines[i].split()[-1])
                b = float(lines[i+1].split()[-1])
                c = float(lines[i+2].split()[-1])
                alpha = float(lines[i+3].split()[-1])
                beta = float(lines[i+4].split()[-1])
                gamma = float(lines[i+5].split()[-1])
            elif line.startswith('_symmetry_cell_setting'):
                lat_type = line.split()[-1]
            elif line.startswith('_symmetry_space_group_name_H-M '):
                symbol = line.split()[-1]
                if eval(symbol) in ["Pn", "P21/n", "C2/n"]:
                    diag = True
                else:
                    diag = False

            elif line.find('_atom_site') >= 0:
                s = i
                while True:
                    s += 1
                    if lines[s].find('_atom_site') >= 0:
                        pass
                    elif len(lines[s].split()) <= 3:
                        break
                    else:
                        tmp = lines[s].split()
                        pos = [float(tmp[-4]), float(tmp[-3]), float(tmp[-2])]
                        species.append(tmp[0])
                        coords.append(pos)
                break

    wp0 = Group(sg)[0]
    lattice = Lattice.from_para(a, b, c, alpha, beta, gamma, lat_type)
    sites = []
    for specie, coord in zip(species, coords):
        pt, wp, _ = wp0.merge(coord, lattice.matrix, tol=0.1)
        sites.append(atom_site(wp, pt, specie, diag))
    return lattice, sites


class structure_from_ext():
    
    def __init__(self, struc, ref_mols, tol=0.2, ignore_HH=False, add_H=False):

        """
        extract the mol_site information from the give cif file 
        and reference molecule
    
        Args: 
            struc: cif/poscar file or a Pymatgen Structure object
            ref_mols: a list of reference molecule (xyz file or Pyxtal molecule)
            tol: scale factor for covalent bond distance
            ignore_HH: whether or not ignore the short H-H distance when checking molecule
            add_H: whether or not add the H atoms
        """

        for ref_mol in ref_mols:
            if isinstance(ref_mol, str):
                ref_mol = pyxtal_molecule(ref_mol, fix=True)
            elif isinstance(ref_mol, pyxtal_molecule):
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

        # reset the hydrogen position
        if add_H: pmg_struc.remove_species('H')

        self.ref_mols = ref_mols
        self.tol = tol
        self.diag = False
        self.add_H = add_H

        sym_struc, number = get_symmetrized_pmg(pmg_struc)
        group = Group(number)
        self.group = group
        self.wyc = group[0]
        self.perm = [0,1,2]

        molecules = search_molecules_in_crystal(sym_struc, self.tol, ignore_HH=ignore_HH)

        self.pmg_struc = sym_struc
        self.lattice = Lattice.from_matrix(sym_struc.lattice.matrix, ltype=group.lattice_type)
        self.resort(molecules)

    def resort(self, molecules):
        from pyxtal.operations import apply_ops, find_ids

        # filter out the molecular generators
        lat = self.pmg_struc.lattice.matrix
        inv_lat = self.pmg_struc.lattice.inv_matrix
        new_lat = self.lattice.matrix
        positions = np.zeros([len(molecules),3])
        for i in range(len(molecules)):
            positions[i, :] = np.dot(molecules[i].cart_coords.mean(axis=0), inv_lat) 

        ids = []  #id for the generator
        mults = [] #the corresponding multiplicities
        visited_ids = []
        for id, pos in enumerate(positions):
            if id not in visited_ids:
                ids.append(id)
                #print("pos", pos); print(self.wyc)
                centers = apply_ops(pos, self.wyc)
                tmp_ids = find_ids(centers, positions)
                visited_ids.extend(tmp_ids)
                mults.append(len(tmp_ids))
                #print("check", id, tmp_ids, mults)

        # add position and molecule, print("ids", ids, mults)
        N_sites = len(ids)
        self.numMols = [0] * len(self.ref_mols)
        self.positions = []
        self.p_mols = []
        self.wps = []
        ids_done = []
        for j, mol2_ref in enumerate(self.ref_mols):
            mol2 = mol2_ref.copy()
            if self.add_H: mol2.mol.remove_species("H")

            for i, id in enumerate(ids):
                mol1 = molecules[id]
                if id not in ids_done and len(mol2.mol) == len(mol1):
                    p_mol = mol2_ref.copy() # create p_mol
                    match, mapping = compare_mol_connectivity(mol2.mol, mol1)
                    if match:
                        self.numMols[j] += mults[i] 
                        if len(mol1) > 1:
                            # rearrange the order
                            order = [mapping[at] for at in range(len(mol1))]
                            xyz = mol1.cart_coords[order] 
                            # add hydrogen positions here
                            if self.add_H: xyz = self.add_Hydrogens(mol2.smile, xyz)
                            #print(xyz)
                            frac = np.dot(xyz, inv_lat) 
                            xyz = np.dot(frac, new_lat) 
                            center = p_mol.get_center(xyz) 
                            p_mol.reset_positions(xyz-center)
                            position = np.dot(center, np.linalg.inv(new_lat))
                        else:
                            position = np.dot(mol1.cart_coords[0], np.linalg.inv(new_lat))
                        position -= np.floor(position)
                        #print(position); print(lat)
                        #print(p_mol.mol.cart_coords[:10] + np.dot(position, new_lat))

                        # check if molecule is on the special wyckoff position
                        if mults[i] < len(self.wyc):
                            #Transform it to the conventional representation
                            if self.diag: position = np.dot(self.perm, position).T
                            #print("molecule is on the special wyckoff position")
                            position, wp, _ = self.wyc.merge(position, new_lat, 0.1)
                            self.wps.append(wp)
                            self.numMols[j] = len(wp)
                            #print("After Merge:---"); print(position); print(wp)
                        else:
                            self.wps.append(self.wyc)

                        self.positions.append(position)
                        self.p_mols.append(p_mol)
                        ids_done.append(id)

        if len(ids_done) < len(ids):
            for id in ids:
                if id not in ids_done:
                    msg = "This molecule cannot be matched to the reference molecule\n"
                    msg += molecules[id].to('xyz')
                    raise ReadSeedError(msg)

    def add_Hydrogens(self, smile, xyz):
        """
        add hydrogen for pymtagen molecule
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Geometry import Point3D

        #print(xyz)
        m1 = Chem.MolFromSmiles(smile)
        m2 = Chem.AddHs(m1)
        if len(smile) > 100:
            AllChem.EmbedMolecule(m2,randomSeed=3)
        else:
            AllChem.EmbedMolecule(m2,randomSeed=0xf00d)
        m2 = Chem.RemoveHs(m2)
        conf = m2.GetConformer(0)

        for i in range(conf.GetNumAtoms()):
            x, y, z = xyz[i]
            conf.SetAtomPosition(i, Point3D(x,y,z))
        #Chem.MolToPDBFile(m2, 'test.pdb')
        #conf = m2.GetConformer(0); print(conf.GetPositions())

        m1 = Chem.AddHs(m1)
        if len(smile) > 100:
            AllChem.EmbedMolecule(m1, randomSeed=3)
        else:
            AllChem.EmbedMolecule(m1, randomSeed=0xf00d)
        AllChem.UFFOptimizeMolecule(m1)
        m3 = AllChem.ConstrainedEmbed(m1,m2)
        conf = m3.GetConformer(0) #; print(conf.GetPositions())
        
        return conf.GetPositions()

    def make_mol_sites(self):
        """
        generate the molecular wyckoff sites
        """
        ori = Orientation(np.eye(3))
        sites = []
        for mol, pos, wp in zip(self.p_mols, self.positions, self.wps):
            site = mol_site(mol, pos, ori, wp, self.lattice, self.diag)
            #print(pos)
            #print(self.lattice.matrix)
            #print([a.value for a in site.molecule.mol.species])
            #print(site.molecule.mol.cart_coords)
            #print(site._get_coords_and_species(absolute=True)[0][:10])
            sites.append(site)
        return sites

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

    def show(self, overlay=True):
        from pyxtal.viz import display_molecules
        if overlay:
            return display_molecules([self.ref_mol, self.mol_aligned])
        else:
            return display_molecules([self.ref_mol, self.molecule])


def search_molecules_in_crystal(struc, tol=0.2, once=False, ignore_HH=True):
    """
    Function to perform to find the molecule in a Pymatgen structure

    Args:
        struc: Pymatgen Structure
        tol: tolerance value to check the connectivity
        once: search only one molecule or all molecules
        ignore_HH: whether or not ignore the short H-H distance when checking molecule

    Returns:
        molecules: list of pymatgen molecules
        positions: list of center positions
    """
    def check_one_layer(struc, sites0, visited):
        new_members = []
        for site0 in sites0:
            sites_add, visited = check_one_site(struc, site0, visited)
            new_members.extend(sites_add)
        return new_members, visited
    
    def check_one_site(struc, site0, visited, rmax=2.8):
        neigh_sites = struc.get_neighbors(site0, rmax)
        ids = [m.index for m in visited]
        sites_add = []
        ids_add = []
        pbc = isinstance(struc, Structure)

        for site1 in neigh_sites:
            if (site1.index not in ids+ids_add):
                try:
                    if CovalentBond.is_bonded(site0, site1, tol):
                        if pbc:
                            (d, image) = site0.distance_and_image(site1)
                        else:
                            d = site0.distance(site1)
                        key = "{:s}-{:s}".format(site1.specie.value, site0.specie.value)

                        #sometime the H-H short distance is not avoidable
                        if key == 'H-H': 
                            if not ignore_HH:
                                if pbc: site1.frac_coords += image
                                sites_add.append(site1)
                                ids_add.append(site1.index)
                        else:
                            if d < bonds[key]:
                                if pbc: site1.frac_coords += image
                                sites_add.append(site1)
                                ids_add.append(site1.index)
                except ValueError:
                    #QZ: use our own bond distance lib
                    if pbc:
                        (d, image) = site0.distance_and_image(site1)
                    else:
                        d = site0.distance(site1)

                    key = "{:s}-{:s}".format(site1.specie.value, site0.specie.value)
                    if d < bonds[key]:
                        if pbc: site1.frac_coords += image
                        sites_add.append(site1)
                        ids_add.append(site1.index)
                    
        if len(sites_add) > 0:
            visited.extend(sites_add)

        return sites_add, visited

    molecules = []
    visited_ids = []

    for id, site in enumerate(struc.sites):
        if id not in visited_ids:
            first_site = site
            visited = [first_site] 
            first_site.index = id
            n_iter, max_iter = 0, len(struc)-len(visited_ids)
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
            molecules.append(Molecule(numbers, coords))
            visited_ids.extend([s.index for s in visited])
            #print(molecules[-1].to(fmt='xyz')); import sys; sys.exit()
        if once and len(molecules) == 1:
            break

    return molecules

if __name__ == "__main__":

    from pyxtal.database.collection import Collection

    pmg = Structure.from_file('pyxtal/database/cifs/resorcinol.cif')
    mols = search_molecules_in_crystal(pmg, tol=0.2, once=False)
    print(len(mols))

    pmg = Collection("molecules")['xxv']
    mols = search_molecules_in_crystal(pmg, tol=0.2, once=False)
    print(len(mols))
