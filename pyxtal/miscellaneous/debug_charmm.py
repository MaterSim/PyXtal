def parse_input(input_file):
    """
    Parse the input file to extract crystal parameters and symmetry operations
    """
    cell_params = []
    symops = ['x, y, z']

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('Crystal Build cutoff'):
                # Extract cell parameters
                num_ops = int(line.split()[-1])
                for i in range(num_ops):
                    next_line = next(f).strip()
                    print(next_line)
                    symops.append(next_line[1:-1])
            elif line.startswith('set shape'):
                for i in range(6):
                    next_line = next(f).strip()
                    if next_line.startswith('set'):
                        cell_params.append(float(next_line.split()[-1]))
    return cell_params, symops

def parse_crd(crd_file):
    atoms = []
    with open(crd_file, 'r') as f:
        lines = f.readlines()
        # Skip header lines
        num_atoms = int(lines[3].strip())
        # Parse atom coordinates
        for i in range(4, 4+num_atoms):
            line = lines[i]
            #print(line[41:])
            x = float(line[23:31])
            y = float(line[31:41])
            z = float(line[41:51])
            atom_type = line[15:19].split()[0].strip()#; print(atom_type, line[:-1])
            atoms.append({
            'type': atom_type,
            'coords': [x,y,z]
            })
    return atoms

def parse_pdb(pdb_file):
    """Extract cell parameters and atom coordinates from PDB format file"""
    cell_params = None
    atoms = []

    with open(pdb_file, 'r') as f:
        for line in f:
            if "REMARK CELL" in line:
                cell_params = [float(x) for x in line.split(':')[1].split()]
            elif line.startswith('ATOM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_type = line[12:16].split()[0].strip()
                atoms.append({
                    'type': atom_type,
                    'coords': [x,y,z]
                })

    return cell_params, atoms

def write_cif(atoms, cell_params, symops, spgnum, outfile):
    """Write structure to CIF format"""
    with open(outfile, 'w') as f:
        # Write cell parameters
        f.write("data_structure\n")
        f.write(f"_symmetry_Int_Tables_number {spgnum}\n")
        f.write(f"_symmetry_cell_setting monoclinic\n")
        f.write(f"_cell_length_a    {cell_params[0]}\n")
        f.write(f"_cell_length_b    {cell_params[1]}\n")
        f.write(f"_cell_length_c    {cell_params[2]}\n")
        f.write(f"_cell_angle_alpha {cell_params[3]}\n")
        f.write(f"_cell_angle_beta  {cell_params[4]}\n")
        f.write(f"_cell_angle_gamma {cell_params[5]}\n")

        f.write("\nloop_\n")
        f.write("_symmetry_equiv_pos_site_id\n")
        f.write("_symmetry_equiv_pos_as_xyz\n")
        for i, symop in enumerate(symops):
            f.write(f"{i} '{symop}'\n")

        # Write atom loop
        f.write("\nloop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")

        for atom in atoms:
            f.write(f"{atom['type']} {atom['coords'][0]:.6f} {atom['coords'][1]:.6f} {atom['coords'][2]:.6f}\n")

# Get cell parameters and atoms from PDB text
prefix = 'hof-1a-g0-p5'
cell_params, symops = parse_input(prefix+'charmm.in')
atoms = parse_crd(prefix+'pyxtal.crd')
write_cif(atoms, cell_params, symops, 9, 'ini.cif')
cell_params, atoms = parse_pdb(prefix+'result.pdb')
write_cif(atoms, cell_params, symops, 9, 'opt.cif')

from pyxtal import pyxtal
smiles = "Nc8nc(N)nc(c7ccc(C(c2ccc(c1nc(N)nc(N)n1)cc2)(c4ccc(c3nc(N)nc(N)n3)cc4)c6ccc(c5nc(N)nc(N)n5)cc6)cc7)n8"
c = pyxtal(molecular=True)
for cif in ['opt.cif', 'ini.cif']:
    c.from_seed(cif, molecules=[smiles+'.smi'])
    print(cif, c)
