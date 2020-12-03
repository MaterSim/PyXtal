"""ASE LAMMPS Calculator Library Version"""
from __future__ import print_function
import numpy as np
import ase.units
from ase.calculators.calculator import Calculator
from lammps import lammps
import os
from pkg_resources import resource_filename
from ase.optimize.fire import FIRE
from ase.optimize import LBFGS
from ase.constraints import ExpCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
#
def opt_lammpslib(struc, lmp, parameters=None, mask=None, logfile='-',
                  path='tmp', calc_type=None, lmp_file=None, molecule=False,
                  method='FIRE', fmax=0.01, opt_cell=False, a=0.1, steps=500):

    """
    mask: [1, 1, 1, 1, 1, 1], a, b, c, alpha, beta, gamma
          [1, 1, 0, 0, 0, 1]
    """

    if lmp_file is not None:
        lammps = LAMMPSlib(lmp=lmp, lmp_file=lmp_file, log_file='lammps.log', \
                           molecule=molecule, path=path)
    elif calc_type is not None:
        lammps = LAMMPSlib(lmp=lmp, lmpcmds=parameters, calc_type=calc_type, \
                           log_file='lammps.log', molecule=molecule, path=path)
    else:
        lammps = LAMMPSlib(lmp=lmp, lmpcmds=parameters, log_file='lammps.log', \
                           molecule=molecule, path=path)

    #check_symmetry(si, 1.0e-6, verbose=True)
    struc.set_calculator(lammps)
    struc.set_constraint(FixSymmetry(struc)) 
    if opt_cell:
        ecf = ExpCellFilter(struc, mask)
        if method == 'FIRE':
            dyn = FIRE(ecf, logfile=logfile, a=a)
        else:
            dyn = LBFGS(ecf, logfile=logfile)
    else:
        if method == 'FIRE':
            dyn = FIRE(struc, logfile=logfile)
        else:
            dyn = LBFGS(struc, logfile=logfile)
    dyn.run(fmax=fmax, steps=steps)
    return struc

def run_lammpslib(struc, lmp, parameters=None, path='tmp', \
                  calc_type=None, lmp_file=None, molecule=False, \
                  method='opt', temp=300, steps=10000, comp_z=None, min_style='cg'):

    if lmp_file is not None:
        lammps = LAMMPSlib(lmp=lmp, lmp_file=lmp_file, log_file='lammps.log', \
                           molecule=molecule, path=path)

    elif calc_type is not None:
        lammps = LAMMPSlib(lmp=lmp, lmpcmds=parameters, calc_type=calc_type, \
                           log_file='lammps.log', molecule=molecule, path=path)
    else:
        parameter0 = parameters + \
                    [
                    "neighbor 1.0 bin",
                    "neigh_modify  every 1  delay 1  check yes",
                    ]

        if method == 'md':
            parameter0 += [
                           "fix 1 all nvt temp " + str(temp) + ' '  + str(temp) + " 0.05",
                           "timestep  0.001",
                           "thermo 1000",
                           "run " + str(steps),
                           "reset_timestep 0",
                           f"min_style {min_style}",
                           f"minimize 1e-15 1e-15 {steps} {steps}",
                           "thermo 0",
                          ]
        elif method == 'opt':
            if comp_z is not None:
                parameter0 += [f'change_box all z scale {comp_z}']
            parameter0 += [
                           "timestep  0.001",
                           f"min_style {min_style}",
                           f"minimize 1e-15 1e-15 {steps} {steps}",
                          ]
        
        else:
            parameter0 += ['run 0']
        
        lammps = LAMMPSlib(lmp=lmp, lmpcmds=parameter0, molecule=molecule, \
                           lmp_file=lmp_file, log_file='lammps.log', path=path)
    
    struc.set_calculator(lammps)
    
    Eng = struc.get_potential_energy()
    
    return struc, Eng


def is_upper_triangular(arr, atol=1e-8):
    """test for upper triangular matrix based on numpy"""
    # must be (n x n) matrix
    assert len(arr.shape) == 2
    assert arr.shape[0] == arr.shape[1]
    return np.allclose(np.tril(arr, k=-1), 0., atol=atol) and \
        np.all(np.diag(arr) >= 0.0)


def convert_cell(ase_cell):
    """
    Convert a parallelepiped (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the bases are column vectors
    """
    cell = ase_cell.T
    if not is_upper_triangular(cell):
        tri_mat = np.zeros((3, 3))
        A = cell[:, 0]
        B = cell[:, 1]
        C = cell[:, 2]
        tri_mat[0, 0] = np.linalg.norm(A)
        Ahat = A / np.linalg.norm(A)
        AxBhat = np.cross(A, B) / np.linalg.norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = np.linalg.norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = np.linalg.norm(np.dot(C, AxBhat))
        # create and save the transformation for coordinates
        volume = np.linalg.det(ase_cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans /= volume
        coord_transform = np.dot(tri_mat, trans)
        return tri_mat, coord_transform
    else:
        return cell, None


class LAMMPSlib(Calculator):

    def __init__(self, lmp, lmpcmds=None, path='tmp', molecule=False,
                lmp_file=None, calc_type='single', *args, **kwargs):
        Calculator.__init__(self, *args, **kwargs)
        self.lmp = lmp
        self.lmp_file = lmp_file
        self.molecule = molecule
        self.calc_type = calc_type
        self.folder = path
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)
        self.lammps_data = self.folder+'/data.lammps'
        self.lammps_in = self.folder + '/in.lammps'
        self.ffpath = resource_filename("pyxtal", "potentials")
        self.lmpcmds = lmpcmds
        self.paras = []
        if lmpcmds is not None:
            for para in lmpcmds:
                if not self.molecule and para.find('pair_coeff') >= 0:
                    tmp = para.split()
                    filename = tmp[3]
                    filename1 = self.ffpath + '/' + filename 
                    para = para.replace(filename, filename1)
                self.paras.append(para)


    def calculate(self, atoms):
        """
        prepare lammps .in file and data file
        write_lammps_data(filename, self.atoms, )
        """
        boundary = ''
        for i in range(3):
            if atoms.pbc[i]:
                boundary += 'p ' 
            else:
                boundary += 'f '
        if boundary in ['f f p ', 'p p f ']: #needs some work later
            boundary = 'p p p '
        self.boundary = boundary
        if self.molecule:
            self.write_lammps_data_water(atoms)
        else:
            self.write_lammps_data(atoms)
        self.write_lammps_in()
        self.lmp.file(self.lammps_in)
        # Extract the forces and energy
        self.lmp.command('variable pxx equal pxx')
        self.lmp.command('variable pyy equal pyy')
        self.lmp.command('variable pzz equal pzz')
        self.lmp.command('variable pxy equal pxy')
        self.lmp.command('variable pxz equal pxz')
        self.lmp.command('variable pyz equal pyz')
        self.lmp.command('variable fx atom fx')
        self.lmp.command('variable fy atom fy')
        self.lmp.command('variable fz atom fz')
        self.lmp.command('variable pe equal pe')
        if self.calc_type.find('GB') >= 0:
            self.lmp.command('variable Etot equal c_eatoms')
            self.gb_energy = self.lmp.extract_variable("Etot", None, 0)
            #print('gb_energy from lammps: ', self.gb_energy)
            #print('update lammps GB energy')

        pos = np.array(
                [x for x in self.lmp.gather_atoms("x", 1, 3)]).reshape(-1, 3)
        
        self.energy = self.lmp.extract_variable('pe', None, 0) 
        #print('update lammps energy')

        xlo = self.lmp.extract_global("boxxlo", 1)
        xhi = self.lmp.extract_global("boxxhi", 1)
        ylo = self.lmp.extract_global("boxylo", 1)
        yhi = self.lmp.extract_global("boxyhi", 1)
        zlo = self.lmp.extract_global("boxzlo", 1)
        zhi = self.lmp.extract_global("boxzhi", 1)
        xy = self.lmp.extract_global("xy", 1)
        yz = self.lmp.extract_global("yz", 1)
        xz = self.lmp.extract_global("xz", 1)
        unitcell = np.array([[xhi-xlo, xy,  xz],
                               [0,  yhi-ylo,  yz],
                               [0,   0,  zhi-zlo]]).T

        stress = np.empty(6)
        stress_vars = ['pxx', 'pyy', 'pzz', 'pyz', 'pxz', 'pxy']

        for i, var in enumerate(stress_vars):
            stress[i] = self.lmp.extract_variable(var, None, 0)
        #print('update lammps stress')

        stress_mat = np.zeros((3, 3))
        stress_mat[0, 0] = stress[0]
        stress_mat[1, 1] = stress[1]
        stress_mat[2, 2] = stress[2]
        stress_mat[1, 2] = stress[3]
        stress_mat[2, 1] = stress[3]
        stress_mat[0, 2] = stress[4]
        stress_mat[2, 0] = stress[4]
        stress_mat[0, 1] = stress[5]
        stress_mat[1, 0] = stress[5]
        stress[0] = stress_mat[0, 0]
        stress[1] = stress_mat[1, 1]
        stress[2] = stress_mat[2, 2]
        stress[3] = stress_mat[1, 2]
        stress[4] = stress_mat[0, 2]
        stress[5] = stress_mat[0, 1]

        self.stress = -stress * 1e5 * ase.units.Pascal
        f = (np.array(self.lmp.gather_atoms("f", 1, 3)).reshape(-1,3) *
                (ase.units.eV/ase.units.Angstrom))
        #print('update lammps force')
        self.forces = f.copy()
        atoms.positions = pos.copy()
        atoms.cell = unitcell.copy()
        if self.molecule:
            atoms.positions *= 0.529
            atoms.cell *= 0.529
        self.atoms = atoms.copy()
        #self.atoms.info['GB_energy'] = self.gb_energy
        #print('update lammps all')

    def write_lammps_in(self):
        if self.lmp_file is not None:
            #print('cp ' + self.lmp_file + ' ' + self.folder + '/in.lammps')
            os.system('cp ' + self.lmp_file + ' ' + self.folder + '/in.lammps')
        else:
            with open(self.lammps_in, 'w') as fh:
                fh.write('clear\n')
                fh.write('box  tilt large\n')
                if self.molecule:
                    fh.write('units electron\n')
                else:
                    fh.write('units metal\n')
                if self.calc_type.find('GB') >= 0:
                    fh.write('boundary p p s\n')
                else:
                    fh.write('boundary {:s}\n'.format(self.boundary))
                fh.write('atom_modify sort 0 0.0\n') 
                if not self.molecule:
                    fh.write('read_data {:s}\n'.format(self.lammps_data))
                fh.write('\n### interactions\n')
                for para in self.paras:
                    fh.write("{:s}\n".format(para))
                fh.write('neighbor 1.0 bin\n')
                fh.write('neigh_modify  every 1  delay 1  check yes\n')
                if self.calc_type.find('GB') >= 0:
                    tmp = LAMMPS_collections().dict['GB_group']
                    for i in tmp:
                        fh.write(i)
                if self.calc_type.find('GB_opt') >= 0:
                    tmp = LAMMPS_collections().dict['GB_opt']
                    for i in tmp:
                        fh.write(i)

                elif self.calc_type.find('GB_md') >= 0:
                    tmp = LAMMPS_collections().dict['GB_opt']
                    tmp += LAMMPS_collections(temp=1000).dict['GB_md']
                    tmp += LAMMPS_collections().dict['GB_opt']
                    for i in tmp:
                        fh.write(i)
                
                if self.calc_type.find('GB') >= 0:
                    fh.write('thermo_style custom pe pxx pyy pzz pyz pxz pxy c_eatoms\n')
                else:
                    fh.write('thermo_style custom pe pxx pyy pzz pyz pxz pxy\n')
                fh.write('thermo_modify flush yes\n')
                fh.write('thermo 1\n')
                fh.write('run 1\n')

                #fh.write('print "__end_of_ase_invoked_calculation__"\n') 

    def write_lammps_data(self, atoms):
        atom_types = np.array(atoms.get_tags()) #[1]*len(atoms)
        if sum(atom_types) == 0:
            atom_types = np.array([1]*len(atoms))
        with open(self.lammps_data, 'w') as fh:
            comment = 'lammpslib autogenerated data file'
            fh.write(comment.strip() + '\n\n')
            fh.write('{:d} atoms\n'.format(len(atoms)))
            fh.write('{:d} atom types\n'.format(len(np.unique(atom_types))))
            cell, coord_transform = convert_cell(atoms.get_cell())
            #cell = atoms.get_cell()
            fh.write('\n')
            fh.write('{0:16.8e} {1:16.8e} xlo xhi\n'.format(0.0, cell[0, 0]))
            fh.write('{0:16.8e} {1:16.8e} ylo yhi\n'.format(0.0, cell[1, 1]))
            fh.write('{0:16.8e} {1:16.8e} zlo zhi\n'.format(0.0, cell[2, 2]))
            fh.write('{0:16.8e} {1:16.8e} {2:16.8e} xy xz yz\n'
                             ''.format(cell[0, 1], cell[0, 2], cell[1, 2]))

            fh.write('\n\nAtoms \n\n')
            for i, (typ, pos) in enumerate(
                    zip(atom_types, atoms.get_positions())):
                if coord_transform is not None:
                    pos = np.dot(coord_transform, pos.transpose())
                fh.write('{:4d} {:4d} {:16.8e} {:16.8e} {:16.8e}\n'
                         .format(i + 1, typ, pos[0], pos[1], pos[2]))

    def write_lammps_data_water(self, atoms):
        """
        Lammps input only for water model
        """
        atom_types = [1]*len(atoms)
        N_atom = len(atoms)
        N_mol = int(len(atoms)/3)
        N_bond = N_mol * 2
        N_angle = N_mol
        n_types = np.unique(atoms.numbers)
        lmp_types = np.zeros(N_atom, dtype=int)
        lmp_types[atoms.numbers==1] = 2
        lmp_types[atoms.numbers==8] = 1

        mol_types = np.zeros(N_atom, dtype=int)
        for i in range(N_mol):
            mol_types[i*3:(i+1)*3] = i+1

        with open(self.lammps_data, 'w') as fh:
            comment = 'lammpslib autogenerated data file'
            fh.write(comment.strip() + '\n\n')
            fh.write('{0} atoms\n'.format(N_atom))
            fh.write('{0} bonds\n'.format(N_bond))
            fh.write('{0} angles\n'.format(N_angle))

            fh.write('\n2 atom types\n')
            fh.write('1 bond types\n')
            fh.write('1 angle types\n')

            #cell = atoms.get_cell()/0.529
            cell, coord_transform = convert_cell(atoms.get_cell())
            cell /= 0.529
            fh.write('\n')
            fh.write('{0:16.8e} {1:16.8e} xlo xhi\n'.format(0.0, cell[0, 0]))
            fh.write('{0:16.8e} {1:16.8e} ylo yhi\n'.format(0.0, cell[1, 1]))
            fh.write('{0:16.8e} {1:16.8e} zlo zhi\n'.format(0.0, cell[2, 2]))
            fh.write('{0:16.8e} {1:16.8e} {2:16.8e} xy xz yz\n'
                             ''.format(cell[0, 1], cell[0, 2], cell[1, 2]))

            fh.write('\n\nMasses \n\n')
            fh.write('  1 15.9994\n')
            fh.write('  2  1.0000\n')

            fh.write('\n\nBond Coeffs \n\n')
            fh.write(' 1    1.78    0.2708585 -0.327738785 0.231328959\n')

            fh.write('\n\nAngle Coeffs \n\n')
            fh.write('  1    0.0700  107.400000')
            fh.write('\n\nAtoms \n\n')
            for i, (typ, mtyp, pos) in enumerate(
                    zip(lmp_types, mol_types, atoms.get_positions()/0.529)):
                if coord_transform is not None:
                    pos = np.dot(coord_transform, pos.transpose())
                #print(i, mtyp, typ)
                if typ==2:
                    fh.write('{0:4d} {1:4d} {2:4d}   0.5564 {3:16.8f} {4:16.8f} {5:16.8f}\n'
                         .format(i + 1, mtyp, typ, pos[0], pos[1], pos[2]))
                else:
                    fh.write('{0:4d} {1:4d} {2:4d}  -1.1128 {3:16.8f} {4:16.8f} {5:16.8f}\n'
                         .format(i + 1, mtyp, typ, pos[0], pos[1], pos[2]))

            fh.write('\nBonds \n\n')
            for i in range(N_mol):
                fh.write('{:4d} {:4d} {:4d} {:4d}\n'.format(i*2+1,1,i*3+1,i*3+2))
                fh.write('{:4d} {:4d} {:4d} {:4d}\n'.format(i*2+2,1,i*3+1,i*3+3))
                   
            fh.write('\nAngles \n\n')
            for i in range(N_angle):
                fh.write('{:4d} {:4d} {:4d} {:4d} {:4d}\n'.format(i+1,1,i*3+2,i*3+1,i*3+3))
                   

    def update(self, atoms):
        if not hasattr(self, 'atoms') or self.atoms != atoms:
            self.calculate(atoms)

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        if self.calc_type.find('GB') >= 0:
            return self.gb_energy
        else:
            return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()



class LAMMPS_collections():
    """
    A class to provide the lammps input file
    """
    def __init__(self, temp=500):
        self.dict = {}

        self.dict['GB_group'] = [
        "# ---------GB--------\n",
        "group lowlayer type 1\n",
        "group gbregion type 2\n",
        "group uplayer type 3\n",
        "#----Define compute quantities-----\n",
        "compute eng all pe/atom \n",
        "compute eatoms gbregion reduce sum c_eng\n"

        ]

        self.dict['GB_opt'] = [
        "# ---------- Run GB Minimization -----\n",
        "fix f1 lowlayer setforce 0.0 0.0 0.0\n",
        "fix f2 uplayer aveforce 0.0 0.0 0.0\n",
        "min_style cg\n",
        "#dump   dump_all all custom 100000 dump.xyz element x y z\n",
        "minimize 1e-15 1e-15 10000 10000\n",
        "run 0\n",
        "thermo 0\n",
        "thermo_style custom pe pxx c_eatoms\n",
        "thermo_modify flush yes\n",
        "unfix f1\n",
        "unfix f2\n",
        ]
        
        self.dict['GB_md'] = [
        "# ---------- Run GB MD ---------\n",
        "thermo 100\n",
        "velocity gbregion create {:d} {:d}\n".format(temp, temp),
        "fix 1 gbregion nvt temp {:d} {:d} 0.05\n".format(temp, temp),
        "timestep 0.001\n",
        "run 5000\n",
        "unfix 1\n",
        ]


if __name__ == '__main__':
    from ase.io import read
    from ase.build import bulk
    lammps_name=''
    comm=None
    log_file='lammps.log'
    cmd_args = ['-echo', 'log', '-log', log_file,
                '-screen', 'none', '-nocite']
    lmp = lammps(lammps_name, cmd_args, comm)

    struc = bulk('Si', 'diamond', cubic=True)
    struc.set_tags([1]*len(struc))
    parameters = ["mass * 1.0",
                  "pair_style tersoff",
                  "pair_coeff * * Si.tersoff Si",
                  "neighbor 1.0 bin",
                  "neigh_modify  every 1  delay 1  check yes",
                  ]
    lammps = LAMMPSlib(lmp=lmp, lmpcmds=parameters)
    struc.set_calculator(lammps)
    print('positions: ')
    print(struc.get_positions())
    print('energy: ')
    print(struc.get_potential_energy())
    print('force: ')
    print(struc.get_forces())
    print('stress: ')
    print(struc.get_stress())

    #dict = LAMMPS_collections().dict
    #for key in dict.keys():
    #    for string in dict[key]:
    #        print(string)
