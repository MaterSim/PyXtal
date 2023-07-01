import os
from ase.io import read
from ase.optimize.fire import FIRE
from ase.constraints import ExpCellFilter
from ase.spacegroup.symmetrize import FixSymmetry
from pyxtal.util import Kgrid
from ase.calculators.calculator import (FileIOCalculator, kpts2ndarray,
                                        kpts2sizeandoffsets)
from ase.units import Hartree, Bohr
import numpy as np


def make_Hamiltonian(skf_dir, atom_types, disp, kpts, write_band=False, use_omp=False):
    """
    Generate the DFTB Hamiltonian for DFTB+
    """
    if disp == 'D3': #from dftb manual
        dispersion = '''DftD3{
                      Damping = BeckeJohnson{
                      a1 = 0.5719
                      a2 = 3.6017
                      }
                      s6 = 1.0
                      s8 = 0.5883
                      }
                      '''

    elif disp == 'D30': #zero dampling
        dispersion = '''DftD3{
                      Damping = ZeroDamping{
                      sr6 = 0.746
                      alpah6 = 4.191
                      }
                      s6 = 1.0
                      s8 = 3.209
                      }
                     '''

    elif disp == 'D4':
        dispersion = '''DftD4{
                      s6 = 1
                      s8 = 0.6635015
                      s9 = 1
                      a1 = 0.5523240
                      a2 = 4.3537076
                        }
                     '''

    elif disp == 'MBD': #1.0 from J. Phys. Chem. Lett. 2018, 9, 399−405
        dispersion = 'MBD{\n\tKGrid = ' + str(kpts)[1:-1] + '\n\tBeta = 1.0}\n'

    elif disp == 'TS': #1.05 from J. Phys. Chem. Lett. 2018, 9, 399−405
        dispersion = '''TS{
                      Damping = 20.0
                      RangeSeparation = 1.0
                      }
                     '''

    elif disp == 'LJ':
        dispersion = 'LennardJones{Parameters = UFFParameters{}}'
    else:
        dispersion = None


    kwargs = {'Hamiltonian_SCC': 'yes',
              'Hamiltonian_SCCTolerance': 1e-06,
              'Hamiltonian_MaxSCCIterations': 1000,
              #'Hamiltonian_Mixer': 'DIIS{}', #Default is Broyden
              #'Hamiltonian_Dispersion': dispersion,
              'slako_dir': skf_dir,
              'Analysis_': '',
              'Analysis_WriteBandOut': 'No',
              'Analysis_MullikenAnalysis': 'No',
              'Analysis_CalculateForces': 'Yes',
             }
    if write_band: 
        kwargs['Analysis_WriteBandOut'] = 'Yes'
    if use_omp:
        kwargs['Parallel_'] = ''
        kwargs['Parallel_UseOmpThreads'] = 'Yes'
    if dispersion is not None:
        kwargs['Hamiltonian_Dispersion'] = dispersion

    if skf_dir.find('3ob') > 0: 
        calc_type = '3ob'
    elif skf_dir.find('mio') > 0: 
        calc_type = 'mio'
    elif skf_dir.find('pbc') > 0:
        calc_type = 'pbc'
    elif skf_dir.find('matsci') > 0:
        calc_type = 'matsci'

    #https://dftb.org/parameters/download/3ob/3ob-3-1-cc
    if calc_type == '3ob':
        kwargs['Hamiltonian_ThirdOrderFull'] = 'Yes'
        kwargs['Hamiltonian_HCorrection'] = 'Damping {\n\tExponent = 4.00\n\t}'
        HD = {"Br": -0.0573,
              "C":  -0.1492,
              "N":  -0.1535,
              "Ca": -0.0340, 
              "Na": -0.0454,
              "Cl": -0.0697, 
              "Zn": -0.03,
              "O":  -0.1575,
              "F":  -0.1623,
              "P":  -0.14,
              "H":  -0.1857, 
              "S":  -0.11,
              "I":  -0.0433, 
              "K":  -0.0339,
             }
        strs = '{'
        for ele in atom_types:
            if ele == 'H':
                kwargs['Hamiltonian_MaxAngularMomentum_H']='s'
            elif ele in ['Mg', 'C', 'N', 'Ca', 'Na', 'O', 'F', 'K']:
                kwargs['Hamiltonian_MaxAngularMomentum_'+ele]='p'
            elif ele in ['Br', 'Cl', 'P', 'S', 'I', 'Zn']:
                kwargs['Hamiltonian_MaxAngularMomentum_'+ele]='d'
            else:
                raise RuntimeError("3-ob-1 doesnot support", ele)
            strs +='\n\t'+ele+' = '+str(HD[ele])
        strs += '\n\t}'
        kwargs['Hamiltonian_HubbardDerivs'] = strs
    elif calc_type == 'pbc':
        #https://dftb.org/parameters/download/pbc/pbc-0-3-cc
        for ele in atom_types:
            if ele == 'H':
                kwargs['Hamiltonian_MaxAngularMomentum_H']='s'
            elif ele in ['C', 'O', 'N', 'F']:
                kwargs['Hamiltonian_MaxAngularMomentum_'+ele]='p'
            elif ele in ['Si', 'Fe']:
                kwargs['Hamiltonian_MaxAngularMomentum_'+ele]='d'
            else:
                raise RuntimeError("pbc-0-3 doesnot support", ele)
    elif calc_type in ['matsci', 'mio']:
         #https://dftb.org/parameters/download/pbc/pbc-0-3-cc
        for ele in atom_types:
            if ele == 'H':
                kwargs['Hamiltonian_MaxAngularMomentum_H']='s'
            elif ele in ['B', 'O', 'C', 'N']:
                kwargs['Hamiltonian_MaxAngularMomentum_'+ele]='p'
            elif ele in ['Si']:
                kwargs['Hamiltonian_MaxAngularMomentum_'+ele]='d'
            else:
                raise RuntimeError(calc_type, "doesnot support", ele)
        
    #DFTB2

    #pbc-0-3
    #matsci
    #ob2
    #pbc
    #print(calc_type, kwargs)
    return kwargs


def DFTB_relax(struc, skf_dir, opt_cell=False, step=500, \
               fmax=0.1, kresol=0.10, folder='tmp', disp='D3', \
               mask=None, symmetrize=True, logfile=None, use_omp=False):
    """
    DFTB optimizer based on ASE

    Args:
        struc: ase atoms object
        mode: [`single`, `relax`, `vc_relax`] (str)
        step: optimization steps (int)
        mask: apply constraints on strain
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
    cwd = os.getcwd()
    os.chdir(folder)

    kpts = Kgrid(struc, kresol)
    atom_types = set(struc.get_chemical_symbols())
    kwargs = make_Hamiltonian(skf_dir, atom_types, disp, kpts, use_omp=use_omp)

    calc = Dftb(label='test',
                atoms=struc,
                kpts=kpts,
                **kwargs,
                )
    struc.set_calculator(calc)

    # impose symmetry
    if symmetrize: struc.set_constraint(FixSymmetry(struc)) 

    # impose cell constraints
    if opt_cell:
        ecf = ExpCellFilter(struc, mask=mask)
        dyn = FIRE(ecf, logfile=logfile)
    else:
        dyn = FIRE(struc, logfile=logfile)
    try:
        dyn.run(fmax=fmax, steps=step)
        os.remove('dftb_pin.hsd')
        os.remove('geo_end.gen')
        os.remove('charges.bin')
    except:
        print("Problem in DFTB calculation")
        struc = None
    os.chdir(cwd)
    return struc

def DFTB_SCF(struc, skf_dir, kresol=0.10, folder='tmp', disp=None, filename=None):
    """
    DFTB SCF to get band structure

    Args:
        struc: ase atoms object
        skf_dir: 
        kresol: 
        filename: band structure output
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
    cwd = os.getcwd()
    os.chdir(folder)

    kpts = Kgrid(struc, kresol)
    atom_types = set(struc.get_chemical_symbols())
    kwargs = make_Hamiltonian(skf_dir, atom_types, disp, kpts, write_band=True)

    calc = Dftb(label='test',
                atoms=struc,
                kpts=kpts,
                **kwargs,
                )

    struc.set_calculator(calc)
    struc.get_potential_energy()
    eigvals = calc.read_eigenvalues()[0]
    ne = calc.read_electrons()
    nband = int(ne/2)
    vbm = eigvals[:, nband-1].max()
    cbm = eigvals[:, nband].min()
    gap = cbm - vbm
    #if filename is not None:
    #    # plot band structure
    os.chdir(cwd)
    return gap


class DFTB():
    """
    Modified DFTB calculator

    Args:
        struc: ase atoms object
        disp: dispersion (D3, D30, TS, MBD))
        step: optimization steps (int)
        kresol: kpoint resolution (float)
        folder: (str)
        label: (str)
        use_omp: True/False
    """

    def __init__(self, struc, skf_dir, 
                 disp = 'D3', 
                 kresol = 0.10, 
                 folder = 'tmp', 
                 label = 'test',
                 prefix = 'geo_final',
                 use_omp = False,
                ):

        self.struc = struc
        self.skf_dir = skf_dir
        self.folder = folder
        self.kresol = kresol
        self.disp = disp
        self.label = label
        self.kpts = Kgrid(struc, kresol)
        self.prefix = prefix
        self.use_omp = use_omp
        if not os.path.exists(self.folder):
           os.makedirs(self.folder)   

    def get_calculator(self, mode, step=500, ftol=1e-3, FixAngles=False, eVperA=True, md_params={}):
        """
        get the ase style calculator

        Args:
            mode: ['single', 'relax', 'vc_relax'] (str)
            step: relaxation steps (int)
            ftol: force tolerance (float)
            FixAngles: Fix angles between lattice vectors
            eVperA: unit in eV/A

        Returns:
            ase calculator
        """
        if eVperA: ftol *= 0.194469064593167E-01
        atom_types = set(self.struc.get_chemical_symbols())
        kwargs = make_Hamiltonian(self.skf_dir, atom_types, self.disp, self.kpts, use_omp=self.use_omp)
       
        if mode in ['relax', 'vc-relax']:
            #kwargs['Driver_'] = 'ConjugateGradient'
            kwargs['Driver_'] = 'GeometryOptimization'
            kwargs['Driver_Optimizer'] = 'FIRE {}'
            #kwargs['Driver_MaxForceComponent'] = ftol
            kwargs['Driver_MaxSteps'] = step
            kwargs['Driver_OutputPrefix'] = self.prefix
            kwargs['Driver_Convergence_'] = ''
            kwargs['Driver_Convergence_GradElem'] = ftol
            
            if mode == 'vc-relax':
                kwargs['Driver_MovedAtoms'] = "1:-1"
                kwargs['Driver_LatticeOpt'] = "Yes"
                if FixAngles:
                    kwargs['Driver_FixAngles'] = "Yes"

        elif mode in ['nve', 'nvt', 'npt']:
            # 1fs = 41.3 au
            # 1000K = 0.0031668 au
            dicts = {
                     'temperature': 300,
                     'pressure': 1e+5, #1atm
                     'timestep': 1,
                     'Thermostat': 'NoseHoover',
                     'MDRestartFrequency': 1000,
                     'band': 'No',
                    }
            dicts.update(md_params)

            kwargs['Analysis_WriteBandOut'] = dicts['band']
            kwargs['Driver_'] = 'VelocityVerlet'
            kwargs['Driver_Steps'] = step
            kwargs['Driver_TimeStep [fs]'] = dicts['timestep']
            kwargs['Driver_MDRestartFrequency'] = dicts['MDRestartFrequency']
            kwargs['Driver_MovedAtoms'] = "1:-1"
            kwargs['Driver_OutputPrefix'] = self.prefix

            if mode in ['nvt', 'npt']:
                kwargs['Driver_Thermostat_'] = dicts['Thermostat']
                kwargs['Driver_Thermostat_Temperature [Kelvin]'] = dicts['temperature']
                kwargs['Driver_Thermostat_CouplingStrength [cm^-1]'] = 3200

                if mode == 'npt':
                    kwargs['Driver_Barostat_'] = ''
                    kwargs['Driver_Barostat_Pressure [Pa]'] = dicts['pressure']
                    kwargs['Driver_Barostat_Timescale [ps]'] = 0.1

    
        calc = Dftb(label=self.label,
                    #run_manyDftb_steps=True,
                    atoms=self.struc,
                    kpts=self.kpts,
                    **kwargs,
                    )
        return calc

    def run(self, mode, step=500, ftol=1e-3, FixAngles=False, md_params={}):
        """
        execute the actual calculation
        """
        from time import time
        t0 = time()
        cwd = os.getcwd()
        os.chdir(self.folder)

        calc = self.get_calculator(mode, step, ftol, FixAngles, md_params=md_params)
        self.struc.set_calculator(calc)
        # self.struc.write('geo_o.gen', format='dftb')
        # execute the simulation
        calc.calculate(self.struc)
        if mode in ['relax', 'vc-relax']:
            final = read(self.prefix+'.gen')
        else:
            final = self.struc

        # get the final energy
        energy = self.struc.get_potential_energy()

        with open(self.label+'.out') as f:
            l = f.readlines()
            self.version = l[2]
        os.chdir(cwd)
        self.time = time() - t0
        return final, energy

class Dftb(FileIOCalculator):
    """ This module defines a FileIOCalculator for DFTB+
    
    http://www.dftbplus.org/
    http://www.dftb.org/
    
    Initial development: markus.kaukonen@iki.fi
    Modified by QZ to avoid the I/O load
    """

    if 'DFTB_COMMAND' in os.environ:
        command = os.environ['DFTB_COMMAND'] + ' > PREFIX.out'
    else:
        command = 'dftb+ > PREFIX.out'

    implemented_properties = ['energy', 'forces', 'stress']
    discard_results_on_any_change = True

    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='dftb', atoms=None, kpts=None,
                 slako_dir=None,
                 **kwargs):
        """
        All keywords for the dftb_in.hsd input file (see the DFTB+ manual)
        can be set by ASE. Consider the following input file block:

        >>> Hamiltonian = DFTB {
        >>>     SCC = Yes
        >>>     SCCTolerance = 1e-8
        >>>     MaxAngularMomentum = {
        >>>         H = s
        >>>         O = p
        >>>     }
        >>> }

        This can be generated by the DFTB+ calculator by using the
        following settings:

        >>> calc = Dftb(Hamiltonian_='DFTB',  # line is included by default
        >>>             Hamiltonian_SCC='Yes',
        >>>             Hamiltonian_SCCTolerance=1e-8,
        >>>             Hamiltonian_MaxAngularMomentum_='',
        >>>             Hamiltonian_MaxAngularMomentum_H='s',
        >>>             Hamiltonian_MaxAngularMomentum_O='p')

        In addition to keywords specific to DFTB+, also the following keywords
        arguments can be used:

        restart: str
            Prefix for restart file.  May contain a directory.
            Default is None: don't restart.
        ignore_bad_restart_file: bool
            Ignore broken or missing restart file. By default, it is an
            error if the restart file is missing or broken.
        label: str (default 'dftb')
            Prefix used for the main output file (<label>.out).
        atoms: Atoms object (default None)
            Optional Atoms object to which the calculator will be
            attached. When restarting, atoms will get its positions and
            unit-cell updated from file.
        kpts: (default None)
            Brillouin zone sampling:

            * ``(1,1,1)`` or ``None``: Gamma-point only
            * ``(n1,n2,n3)``: Monkhorst-Pack grid
            * ``dict``: Interpreted as a path in the Brillouin zone if
              it contains the 'path_' keyword. Otherwise it is converted
              into a Monkhorst-Pack grid using
              ``ase.calculators.calculator.kpts2sizeandoffsets``
            * ``[(k11,k12,k13),(k21,k22,k23),...]``: Explicit (Nkpts x 3)
              array of k-points in units of the reciprocal lattice vectors
              (each with equal weight)

        """

        if slako_dir is None:
            slako_dir = os.environ.get('DFTB_PREFIX', './')
            if not slako_dir.endswith('/'):
                slako_dir += '/'

        self.slako_dir = slako_dir

        self.default_parameters = dict(
            Hamiltonian_='DFTB',
            Hamiltonian_SlaterKosterFiles_='Type2FileNames',
            Hamiltonian_SlaterKosterFiles_Prefix=self.slako_dir,
            Hamiltonian_SlaterKosterFiles_Separator='"-"',
            Hamiltonian_SlaterKosterFiles_Suffix='".skf"',
            Hamiltonian_MaxAngularMomentum_='',
            Options_='',
            Options_WriteResultsTag='Yes',
            Options_WriteDetailedOut='No',

            )

        self.lines = None
        self.atoms = None
        self.atoms_input = None
        self.outfilename = 'dftb.out'

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms,
                                  **kwargs)

        # kpoint stuff by ase
        self.kpts = kpts
        self.kpts_coord = None

        if self.kpts is not None:
            initkey = 'Hamiltonian_KPointsAndWeights'
            mp_mesh = None
            offsets = None

            if isinstance(self.kpts, dict):
                if 'path' in self.kpts:
                    # kpts is path in Brillouin zone
                    self.parameters[initkey + '_'] = 'Klines '
                    self.kpts_coord = kpts2ndarray(self.kpts, atoms=atoms)
                else:
                    # kpts is (implicit) definition of
                    # Monkhorst-Pack grid
                    self.parameters[initkey + '_'] = 'SupercellFolding '
                    mp_mesh, offsets = kpts2sizeandoffsets(atoms=atoms,
                                                           **self.kpts)
            elif np.array(self.kpts).ndim == 1:
                # kpts is Monkhorst-Pack grid
                self.parameters[initkey + '_'] = 'SupercellFolding '
                mp_mesh = self.kpts
                offsets = [0.] * 3
            elif np.array(self.kpts).ndim == 2:
                # kpts is (N x 3) list/array of k-point coordinates
                # each will be given equal weight
                self.parameters[initkey + '_'] = ''
                self.kpts_coord = np.array(self.kpts)
            else:
                raise ValueError('Illegal kpts definition:' + str(self.kpts))

            if mp_mesh is not None:
                eps = 1e-10
                for i in range(3):
                    key = initkey + '_empty%03d' % i
                    val = [mp_mesh[i] if j == i else 0 for j in range(3)]
                    self.parameters[key] = ' '.join(map(str, val))
                    offsets[i] *= mp_mesh[i]
                    assert abs(offsets[i]) < eps or abs(offsets[i] - 0.5) < eps
                    # DFTB+ uses a different offset convention, where
                    # the k-point mesh is already Gamma-centered prior
                    # to the addition of any offsets
                    if mp_mesh[i] % 2 == 0:
                        offsets[i] += 0.5
                key = initkey + '_empty%03d' % 3
                self.parameters[key] = ' '.join(map(str, offsets))

            elif self.kpts_coord is not None:
                for i, c in enumerate(self.kpts_coord):
                    key = initkey + '_empty%09d' % i
                    c_str = ' '.join(map(str, c))
                    if 'Klines' in self.parameters[initkey + '_']:
                        c_str = '1 ' + c_str
                    else:
                        c_str += ' 1.0'
                    self.parameters[key] = c_str

    def write_dftb_in(self, outfile):
        """ Write the innput file for the dftb+ calculation.
            Geometry is taken always from the file 'geo_end.gen'.
        """

        outfile.write('Geometry = GenFormat { \n')
        outfile.write('    <<< "geo_end.gen" \n')
        outfile.write('} \n')
        outfile.write(' \n')

        params = self.parameters.copy()

        s = 'Hamiltonian_MaxAngularMomentum_'
        for key in params:
            if key.startswith(s) and len(key) > len(s):
                break
        # --------MAIN KEYWORDS-------
        previous_key = 'dummy_'
        myspace = ' '
        for key, value in sorted(params.items()):
            current_depth = key.rstrip('_').count('_')
            previous_depth = previous_key.rstrip('_').count('_')
            for my_backsclash in reversed(
                    range(previous_depth - current_depth)):
                outfile.write(3 * (1 + my_backsclash) * myspace + '} \n')
            outfile.write(3 * current_depth * myspace)
            if key.endswith('_') and len(value) > 0:
                outfile.write(key.rstrip('_').rsplit('_')[-1] +
                              ' = ' + str(value) + '{ \n')
            elif (key.endswith('_') and (len(value) == 0)
                  and current_depth == 0):  # E.g. 'Options {'
                outfile.write(key.rstrip('_').rsplit('_')[-1] +
                              ' ' + str(value) + '{ \n')
            elif (key.endswith('_') and (len(value) == 0)
                  and current_depth > 0):  # E.g. 'Hamiltonian_Max... = {'
                outfile.write(key.rstrip('_').rsplit('_')[-1] +
                              ' = ' + str(value) + '{ \n')
            elif key.count('_empty') == 1:
                outfile.write(str(value) + ' \n')
            else:
                outfile.write(key.rsplit('_')[-1] + ' = ' + str(value) + ' \n')
            previous_key = key
        current_depth = key.rstrip('_').count('_')
        for my_backsclash in reversed(range(current_depth)):
            outfile.write(3 * my_backsclash * myspace + '} \n')
        outfile.write('ParserOptions { \n')
        outfile.write('   IgnoreUnprocessedNodes = Yes  \n')
        outfile.write('} \n')

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore unit cell for molecules:
        if not atoms.pbc.any() and 'cell' in system_changes:
            system_changes.remove('cell')
        return system_changes

    def write_input(self, atoms, properties=None, system_changes=None):
        from ase.io import write
        FileIOCalculator.write_input(
            self, atoms, properties, system_changes)
        with open(os.path.join(self.directory, 'dftb_in.hsd'), 'w') as fd:
            self.write_dftb_in(fd)
        write(os.path.join(self.directory, 'geo_end.gen'), atoms,
              parallel=False)
        # self.atoms is none until results are read out,
        # then it is set to the ones at writing input
        self.atoms_input = atoms
        self.atoms = None

    def read_results(self):
        """ all results are read from results.tag file
            It will be destroyed after it is read to avoid
            reading it once again after some runtime error """

        with open(os.path.join(self.directory, 'results.tag'), 'r') as fd:
            self.lines = fd.readlines()

        self.atoms = self.atoms_input
        self.results['energy'] = float(self.lines[1])*Hartree
        forces = self.read_forces()
        self.results['forces'] = forces

        # stress stuff begins
        sstring = 'stress'
        have_stress = False
        stress = list()
        for iline, line in enumerate(self.lines):
            if sstring in line:
                have_stress = True
                start = iline + 1
                end = start + 3
                for i in range(start, end):
                    cell = [float(x) for x in self.lines[i].split()]
                    stress.append(cell)
        if have_stress:
            stress = -np.array(stress) * Hartree / Bohr**3
            self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]
        # stress stuff ends

        # calculation was carried out with atoms written in write_input
        os.remove(os.path.join(self.directory, 'results.tag'))

    def read_forces(self):
        """Read Forces from dftb output file (results.tag)."""

        # Initialise the indices so their scope
        # reaches outside of the for loop
        index_force_begin = -1
        index_force_end = -1
        
        # Force line indexes
        fstring = 'forces   '
        for iline, line in enumerate(self.lines):
            if line.find(fstring) >= 0:
                index_force_begin = iline + 1
                line1 = line.replace(':', ',')
                index_force_end = iline + 1 + \
                    int(line1.split(',')[-1])
                break
        gradients = []
        for j in range(index_force_begin, index_force_end):
            word = self.lines[j].split()
            gradients.append([float(word[k]) for k in range(0, 3)])
        gradients = np.array(gradients)* Hartree / Bohr

        return gradients 
    
    def read_eigenvalues(self):
        """ Read Eigenvalues from dftb output file (results.tag).
            Unfortunately, the order seems to be scrambled. """
        # Eigenvalue line indexes
        index_eig_begin = None
        for iline, line in enumerate(self.lines):
            fstring = 'eigenvalues   '
            if line.find(fstring) >= 0:
                index_eig_begin = iline + 1
                line1 = line.replace(':', ',')
                ncol, nband, nkpt, nspin = map(int, line1.split(',')[-4:])
                break
        else:
            return None

        # Take into account that the last row may lack
        # columns if nkpt * nspin * nband % ncol != 0
        nrow = int(np.ceil(nkpt * nspin * nband * 1. / ncol))
        index_eig_end = index_eig_begin + nrow
        ncol_last = len(self.lines[index_eig_end - 1].split())
        if ncol - ncol_last > 0:
            self.lines[index_eig_end - 1] = self.lines[index_eig_end - 1].replace('\n', '')
            self.lines[index_eig_end - 1] += ' 0.0 ' * (ncol - ncol_last)
            self.lines[index_eig_end - 1] += '\n'
        eig = np.loadtxt(self.lines[index_eig_begin:index_eig_end]).flatten()
        eig *= Hartree
        N = nkpt * nband
        eigenvalues = [eig[i * N:(i + 1) * N].reshape((nkpt, nband))
                       for i in range(nspin)]

        return eigenvalues

    def read_fermi_levels(self):
        """ Read Fermi level(s) from dftb output file (results.tag). """
        # Fermi level line indexes
        for iline, line in enumerate(self.lines):
            fstring = 'fermi_level   '
            if line.find(fstring) >= 0:
                index_fermi = iline + 1
                break
        else:
            return None

        fermi_levels = []
        words = self.lines[index_fermi].split()
        assert len(words) in [1, 2], 'Expected either 1 or 2 Fermi levels'

        for word in words:
            e = float(word)
            # In non-spin-polarized calculations with DFTB+ v17.1,
            # two Fermi levels are given, with the second one being 0,
            # but we don't want to add that one to the list
            if abs(e) > 1e-8:
                fermi_levels.append(e)

        return np.array(fermi_levels) * Hartree

    def read_electrons(self):
        """read number o electrons"""
        for iline, line in enumerate(self.lines):
            fstring = 'number_of_electrons'
            if line.find(fstring) >= 0:
                index_ele = iline + 1
                break
        return float(self.lines[index_ele].split('\n')[0])

if __name__ == '__main__':
    from ase.build import bulk

    skf_dir = os.environ['DFTB_PREFIX'] + 'pbc-0-3/'
    #skf_dir = os.environ['DFTB_PREFIX'] + '3ob-3-1/'

    struc = bulk('Si', 'diamond', cubic=True)
    struc.set_cell(1.1*struc.cell)
    for mode in ['single', 'relax', 'vc-relax', 'npt']:
        my = DFTB(struc, skf_dir)
        struc, energy = my.run(mode)
        res = "{:8s} ".format(mode)
        res += "{:8.4f} ".format(struc.cell[0,0])
        res += "{:8.4f} ".format(struc.cell[1,1])
        res += "{:8.4f} ".format(struc.cell[2,2])
        res += "{:12.4f}".format(energy)
        print(res)
    #gap = DFTB_SCF(struc, skf_dir)
    #print(gap)
