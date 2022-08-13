import os
from ase.io import read
from ase.calculators.dftb import Dftb
from ase.optimize.fire import FIRE
from ase.constraints import ExpCellFilter
from ase.spacegroup.symmetrize import FixSymmetry
from pyxtal.util import Kgrid

def make_Hamiltonian(skf_dir, atom_types, disp, kpts):
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


    kwargs = {'Hamiltonian_SCC': 'yes',
              'Hamiltonian_SCCTolerance': 1e-06,
              'Hamiltonian_MaxSCCIterations': 1000,
              'Hamiltonian_Mixer': 'DIIS{}',
              'Hamiltonian_Dispersion': dispersion,
              'slako_dir': skf_dir,
             }

    if skf_dir.find('3ob') > 0: 
        calc_type = '3ob'
    elif skf_dir.find('mio') > 0 or skf_dir.find('pbc') > 0:
        calc_type = 'mio'

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
    #elif calc_type == 'mio':
    #DFTB2

    #pbc-0-3
    #matsci
    #ob2
    #pbc
    #print(calc_type, kwargs)
    return kwargs


def DFTB_relax(struc, skf_dir, opt_cell=False, step=500, fmax=0.1, kresol=0.10, folder='tmp', disp='D3', logfile=None):
    """
    DFTB optimizer based on ASE

    Args:
        struc: ase atoms object
        mode: [`single`, `relax`, `vc_relax`] (str)
        step: optimization steps (int)
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
    cwd = os.getcwd()
    os.chdir(folder)

    kpts = Kgrid(struc, kresol)
    atom_types = set(struc.get_chemical_symbols())
    kwargs = make_Hamiltonian(skf_dir, atom_types, disp, kpts)

    calc = Dftb(label='test',
                atoms=struc,
                kpts=kpts,
                **kwargs,
                )

    struc.set_calculator(calc)
    struc.set_constraint(FixSymmetry(struc)) 
    if opt_cell:
        ecf = ExpCellFilter(struc)
        if logfile is not None:
            dyn = FIRE(ecf, logfile=logfile)
        else:
            dyn = FIRE(ecf)
    else:
        if logfile is not None:
            dyn = FIRE(struc, logfile=logfile)
        else:
            dyn = FIRE(struc)
    try:
        dyn.run(fmax=fmax, steps=step)
        os.remove('detailed.out')
        os.remove('dftb_pin.hsd')
        os.remove('geo_end.gen')
        os.remove('charges.bin')
    except:
        print("Problem in DFTB calculation")
    os.chdir(cwd)
    return struc

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
    """

    def __init__(self, struc, skf_dir, 
                 disp = 'D3', 
                 kresol = 0.10, 
                 folder = 'tmp', 
                 label = 'test',
                ):

        self.struc = struc
        self.skf_dir = skf_dir
        self.folder = folder
        self.kresol = kresol
        self.disp = disp
        self.label = label
        self.kpts = Kgrid(struc, kresol)
        if not os.path.exists(self.folder):
           os.makedirs(self.folder)   

    def get_calculator(self, mode, step=500, ftol=1e-3, eVperA=True):
        """
        get the ase style calculator

        Args:
            mode: ['single', 'relax', 'vc_relax'] (str)
            step: relaxation steps (int)
            ftol: force tolerance (float)
            eVperA: unit in eV/A

        Returns:
            ase calculator
        """
        if eVperA: ftol *= 0.194469064593167E-01
        atom_types = set(self.struc.get_chemical_symbols())
        kwargs = make_Hamiltonian(self.skf_dir, atom_types, self.disp, self.kpts)
       
        if mode in ['relax', 'vc-relax']:
            #kwargs['Driver_'] = 'ConjugateGradient'
            kwargs['Driver_'] = 'FIRE'
            kwargs['Driver_MaxForceComponent'] = ftol
            kwargs['Driver_MaxSteps'] = step
            
            if mode == 'vc-relax':
                kwargs['Driver_MovedAtoms'] = "1:-1"
                kwargs['Driver_LatticeOpt'] = "Yes"
    
        calc = Dftb(label=self.label,
                    #run_manyDftb_steps=True,
                    atoms=self.struc,
                    kpts=self.kpts,
                    **kwargs,
                    )
        return calc

    def run(self, mode, step=500, ftol=1e-3):
        """
        execute the actual calculation
        """
        from time import time
        t0 = time()
        cwd = os.getcwd()
        os.chdir(self.folder)

        calc = self.get_calculator(mode, step, ftol)
        self.struc.set_calculator(calc)
        self.struc.write('geo_o.gen', format='dftb')
        # execute the simulation
        calc.calculate(self.struc)
        try:
            final = read('geo_end.gen')
        except:
            print("Problem in reading the final structure", time()-t0)
            final = self.struc

        # get the final energy
        energy = self.struc.get_potential_energy()

        with open('test.out') as f:
            l = f.readlines()
            self.version = l[2]
        os.chdir(cwd)
        self.time = time() - t0
        return final, energy


if __name__ == '__main__':
    from ase.build import bulk

    skf_dir = os.environ['DFTB_PREFIX'] + 'pbc-0-3/'
    #skf_dir = os.environ['DFTB_PREFIX'] + '3ob-3-1/'

    struc = bulk('Si', 'diamond', cubic=True)
    struc.set_cell(1.1*struc.cell)
    for mode in ['single', 'relax', 'vc-relax']:
        my = DFTB(struc, skf_dir)
        struc, energy = my.run(mode)
        res = "{:8s} ".format(mode)
        res += "{:8.4f} ".format(struc.cell[0,0])
        res += "{:8.4f} ".format(struc.cell[1,1])
        res += "{:8.4f} ".format(struc.cell[2,2])
        res += "{:12.4f}".format(energy)
        print(res)
