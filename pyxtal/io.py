"""
CIF in PyXtal format
"""
from pyxtal.constants import deg, logo

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

