"""
Database class
"""
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from ase.db import connect
from pyxtal import pyxtal
import pymatgen.analysis.structure_matcher as sm
from pyxtal.util import ase2pymatgen
from ase.calculators.calculator import CalculationFailed

def dftb_opt_par(ids, xtals, skf_dir, steps, folder, symmetrize, criteria):
    """
    Run GULP optimization in parallel for a list of atomic xtals
    Args:
        xtal: pyxtal instance
        ff (str): e.g., `reaxff`, `tersoff`
        path (str): path of calculation folder
        criteria (dicts): to check if the structure
    """
    cwd = os.getcwd()
    os.chdir(folder)
    results = []
    for id, xtal in zip(ids, xtals):
        res = dftb_opt_single(id, xtal, skf_dir, steps, symmetrize, criteria)
        (xtal, eng, status) = res
        if status:
            results.append((id, xtal, eng))
    os.chdir(cwd)
    return results

def dftb_opt_single(id, xtal, skf_dir, steps, symmetrize, criteria, kresol=0.05):
    """
    Single DFTB optimization for a given atomic xtal

    Args:
        id (int): id of the give xtal
        xtal: pyxtal instance
        skf_dir (str): path of skf files
        steps (int): number of relaxation steps
        criteria (dicts): to check if the structure
    """
    from pyxtal.interface.dftb import DFTB_relax, DFTB
    cwd = os.getcwd()
    atoms = xtal.to_ase(resort=False)
    eng = None
    stress = None
    try:
        if symmetrize:
            s = DFTB_relax(atoms, skf_dir, True, int(steps/2),
                           kresol=kresol*1.2, folder='.',
                           scc_iter = 100,
                           logfile='ase.log')
            s = DFTB_relax(atoms, skf_dir, True, int(steps/2),
                           kresol=kresol, folder='.',
                           scc_error = 1e-5,
                           scc_iter = 100,
                           logfile='ase.log')
            stress = np.sum(s.get_stress()[:3])/0.006241509125883258/3
        else:
            my = DFTB(atoms, skf_dir, kresol=kresol*1.5, folder='.',
                      scc_error=0.1, scc_iter=100)
            s, eng = my.run(mode='vc-relax', step=int(steps/2))
            my = DFTB(s, skf_dir, kresol=kresol, folder='.',
                      scc_error=1e-4, scc_iter=100)
            s, eng = my.run(mode='vc-relax', step=int(steps/2))
            s = my.struc
    except CalculationFailed:
        # This is due to covergence error in geometry optimization
        # Here we simply read the last energy
        my.calc.read_results()
        eng = my.calc.results['energy']
        s = my.struc
    except:
        s = None
        print("Problem in DFTB Geometry optimization", id)
        xtal.to_file('bug.cif')#; import sys; sys.exit()
    os.chdir(cwd)

    if s is not None:
        c = pyxtal(); c.from_seed(s)
        if eng is None:
            eng = s.get_potential_energy()/len(s)
        else:
            eng /= len(s)

        if criteria is not None:
            status = xtal.check_validity(criteria)
        else:
            status = True

        header = "{:4d}".format(id)
        dicts = {'validity': status, 'energy': eng}
        if stress is not None: dicts['stress'] = stress
        print(xtal.get_xtal_string(header=header, dicts=dicts))

        return xtal, eng, status
    else:
        return None, None, False

def gulp_opt_par(ids, xtals, ff, path, criteria):
    """
    Run GULP optimization in parallel for a list of atomic xtals
    Args:
        xtal: pyxtal instance
        ff (str): e.g., `reaxff`, `tersoff`
        path (str): path of calculation folder
        criteria (dicts): to check if the structure
    """
    results = []
    for id, xtal in zip(ids, xtals):
        res = gulp_opt_single(id, xtal, ff, path, criteria)
        (xtal, eng, status) = res
        if status:
            results.append((id, xtal, eng))
    return results


def gulp_opt_single(id, xtal, ff, path, criteria):
    """
    Single GULP optimization for a given atomic xtal

    Args:
        xtal: pyxtal instance
        ff (str): e.g., `reaxff`, `tersoff`
        path (str): path of calculation folder
        criteria (dicts): to check if the structure
    """
    from pyxtal.interface.gulp import single_optimize as gulp_opt

    xtal, eng, _, error = gulp_opt(xtal,
                                   ff=ff,
                                   label=str(id),
                                   #clean=False,
                                   path=path,
                                   symmetry=True)
    status = False
    if not error:
        if criteria is not None:
            status = xtal.check_validity(criteria)
        else:
            status = True
    if status:
        header = "{:4d}".format(id)
        dicts = {'validity': status, 'energy': eng}
        print(xtal.get_xtal_string(header=header, dicts=dicts))
    return xtal, eng, status


def make_entry_from_pyxtal(xtal):
    """
    make entry from the pyxtal object, assuming that
    the smiles/ccdc_number info is given

    Args:
        xtal: pyxtal object

    Returns:
        entry dictionary
    """
    from rdkit.Chem.Descriptors import ExactMolWt
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from rdkit import Chem

    if xtal.valid:
        url0 = "https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid="
        m = Chem.MolFromSmiles(xtal.tag['smiles'])
        mol_wt = ExactMolWt(m)
        mol_formula = CalcMolFormula(m)
        kvp = {
                "csd_code": xtal.tag['csd_code'],
                "mol_smi": xtal.tag['smiles'],
                "ccdc_number": xtal.tag['ccdc_number'],
                "space_group": xtal.group.symbol,
                "spg_num": xtal.group.number,
                "Z": sum(xtal.numMols),
                "Zprime": xtal.get_zprime()[0],
                "url": url0 + str(xtal.tag['ccdc_number']),
                "mol_formula": mol_formula,
                "mol_weight": mol_wt,
                "mol_name": xtal.tag['csd_code'],
                "l_type": xtal.lattice.ltype,
              }
        entry = (xtal.to_ase(), kvp, None)
        return entry
    else:
        return None

def make_entry_from_CSD_web(code, number, smiles, name=None):
    """
    make enetry dictionary from csd web https://www.ccdc.cam.ac.uk/structures

    Args:
        code: CSD style letter entry
        number: ccdc number
        smiles: the corresponding molecular smiles
        name: name of the compound
    """

    url0 = "https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid="
    #xtal = pyxtal(molecular=True)
    #
    #return make_entry_from_pyxtal(xtal)
    raise NotImplementedError('To do in future')


def make_entry_from_CSD(code):
    """
    make entry dictionary from CSD codes

    Args:
        code: a list of CSD codes
    """
    from pyxtal.msg import CSDError

    xtal = pyxtal(molecular=True)
    try:
        xtal.from_CSD(code)
        return make_entry_from_pyxtal(xtal)

    except CSDError as e:
        print("CSDError", code, e.message)
        return None

def make_db_from_CSD(dbname, codes):
    """
    make database from CSD codes

    Args:
        dbname: db file name
        codes: a list of CSD codes
    """
    # open
    db = database(dbname)

    # add structure
    for i, code in enumerate(codes):
        entry = make_entry_from_CSD(code)
        if entry is not None:
            db.add(entry)
            print(i, code)
    return db

class database():
    """
    This is a database class to process crystal data

    Args:
        db_name: *.db format from ase database
    """

    def __init__(self, db_name):
        self.db_name = db_name
        #if not os.path.exists(db_name):
        #    raise ValueError(db_name, 'doesnot exist')

        self.db = connect(db_name)
        self.get_all_codes()
        self.keys = ['csd_code', 'space_group', 'spg_num',
                     'Z', 'Zprime', 'url', 'mol_name',
                     'mol_smi', 'mol_formula', 'mol_weight',
                     'l_type',
                    ]
        self.calculators = ['charmm',
                            'gulp',
                            'ani',
                            'dftb_D3',
                            'dftb_TS',
                           ]

    def vacuum(self):
        self.db.vacuum()

    def get_all_codes(self):
        codes = []
        for row in self.db.select():
            if row.csd_code not in codes:
                codes.append(row.csd_code)
            else:
                print('find duplicate! remove', row.id, row.csd_code)
                self.db.delete([row.id])
        self.codes = codes

    def add(self, entry):
        (atom, kvp, data) = entry
        if kvp['csd_code'] not in self.codes:
            kvp0 = self.process_kvp(kvp)
            self.db.write(atom, key_value_pairs=kvp0, data=data)
            self.codes.append(kvp['csd_code'])

    def add_from_code(self, code):
        entry = make_entry_from_CSD(code)
        if entry is not None:
            self.add(entry)
        else:
            print("{:s} is not a valid entry".format(code))

    def process_kvp(self, kvp):
        kvp0 = {}
        for key in self.keys:
            if key in kvp:
                kvp0[key] = kvp[key]
            else:
                print('Error, cannot find ', key, ' from the input')
                return
        return kvp0

    def check_status(self, show=False):
        """
        Check the current status of each entry
        """
        ids = []
        for row in self.db.select():
            if len(row.data.keys())==len(self.calculators):
                ids.append(row.id)
                if show:
                    row_info = self.get_row_info(id=row.id)
                    self.view(row_info)
            else:
                print(row.csd_code) #, row.data['charmm_info']['prm'])
        return ids

    def copy(self, db_name, csd_codes):
        """
        copy the entries to another db

        Args:
            db_name: db file name
            csd_codes: list of codes
        """
        if db_name == self.db_name:
            raise RuntimeError('Cannot use the same db file for copy')
        with connect(db_name) as db:
            for csd_code in csd_codes:
                row_info = self.get_row_info(code=csd_code)
                (atom, kvp, data) = row_info
                db.write(atom, key_value_pairs=kvp, data=data)

    def view(self, row_info):
        """
        print the summary of benchmark results

        Args:
            row: row object
        """
        from pyxtal import pyxtal
        from pyxtal.util import ase2pymatgen
        from pyxtal.representation import representation
        from pyxtal.msg import ReadSeedError

        (atom, kvp, data) = row_info

        #Reference
        xtal = self.get_pyxtal(kvp['csd_code'])
        rep = xtal.get_1D_representation()
        print('\n', kvp['csd_code'], kvp['mol_smi'], xtal.lattice.volume)
        print(rep.to_string() + ' reference')

        #calcs
        for key in data.keys():
            calc = key[:-5]
            time = data[key]['time']

            rep = data[key]['rep']
            if type(rep[0]) is not list: rep = [rep]
            rep = representation(rep, kvp['mol_smi']).to_string()

            (dv, msd1, msd2) = data[key]['diff']
            strs = "{:s} {:8s} {:6.2f} {:6.3f}".format(rep, calc, time/60, dv)
            if msd1 is not None:
                strs += "{:6.3f}{:6.3f}".format(msd1, msd2)
            print(strs)


    def get_row_info(self, id=None, code=None):
        match = False
        if id is not None:
            for row in self.db.select(id=id):
                match = True
                break
        elif code is not None:
            for row in self.db.select(csd_code=code):
                match = True
                break
        if match:

            kvp = {}
            for key in self.keys:
                kvp[key] = row.key_value_pairs[key]

            data0 = {}
            for calc in self.calculators:
                key = calc + '_info'
                if key in row.data.keys():
                    data0[key] = row.data[key]

            atom = self.db.get_atoms(id=row.id)
            return (atom, kvp, data0)
        else:
            msg = 'cannot find the entry from ' + id + code
            raise RuntimeError(msg)

    def get_row(self, code):
        for row in self.db.select(csd_code=code):
            return row
        msg = 'cannot find the entry from ' + code
        raise RuntimeError(msg)

    def get_pyxtal(self, code):
        from pyxtal import pyxtal
        from pyxtal.util import ase2pymatgen
        from pyxtal.msg import ReadSeedError

        row = self.get_row(code)
        atom = self.db.get_atoms(id=row.id)
        #Reference
        pmg = ase2pymatgen(atom)
        smi = row.mol_smi
        smiles = smi.split('.')
        molecules = [smile+'.smi' for smile in smiles]

        xtal = pyxtal(molecular=True)
        try:
            xtal.from_seed(pmg, molecules=molecules)
        except ReadSeedError:
            xtal.from_seed(pmg, molecules=molecules, add_H=True)

        return xtal

    def compute(self, row, work_dir, skf_dir):
        if len(row.data.keys())<len(self.calculators):
            #not label information, run antechamber
            atom = self.db.get_atoms(id=row.id)
            if 'gulp_info' not in row.data.keys():
                pmg, c_info, g_info = get_parameters(row, atom)
                row.data = {'charmm_info': c_info,
                            'gulp_info': g_info
                           }
            else:
                pmg = ase2pymatgen(atom)

            data = compute(row, pmg, work_dir, skf_dir)
            self.db.update(row.id, data=data)
            print('updated the data for', row.csd_code)

class database_topology():
    """
    This is a database class to process atomic crystal data

    Args:
        db_name: *.db format from ase database
    """

    def __init__(self, db_name, ltol=0.05, stol=0.05, atol=3):
        #if not os.path.exists(db_name):
        #    raise ValueError(db_name, 'doesnot exist')

        self.db_name = db_name
        self.db = connect(db_name)
        self.keys = ['space_group_number',
                     'similarity0',
                     'similarity',
                     'density',
                     'dof',
                     'topology',
                     'topology_detail',
                     'dimension',
                     'wps',
                     'ff_energy',
                     'ff_lib',
                     'ff_relaxed',
                     'dftb_energy',
                     'dftb_relaxed',
                    ]
        self.matcher = sm.StructureMatcher(ltol=ltol, stol=stol, angle_tol=atol)

    def vacuum(self):
        self.db.vacuum()

    def get_pyxtal(self, id, use_ff):
        """
        Get pyxtal based on row_id, if use_ff, get pyxtal from the ff_relaxed file
        """
        from pyxtal import pyxtal
        from pyxtal.util import ase2pymatgen
        from pymatgen.core import Structure

        row = self.db.get(id) #; print(id, row.id)
        if use_ff and hasattr(row, 'ff_relaxed'):
            pmg = Structure.from_str(row.ff_relaxed, fmt='cif')
        else:
            atom = self.db.get_atoms(id=id)
            pmg = ase2pymatgen(atom)
        xtal = pyxtal()
        try:
            xtal.from_seed(pmg)
            if xtal is not None and xtal.valid:
                for key in self.keys:
                    if hasattr(row, key):
                        setattr(xtal, key, getattr(row, key))
            return xtal
        except:
            print('Cannot load the structure')

    def get_all_xtals(self):
        """
        Get all pyxtal instances from the current db
        """
        xtals = []
        for row in self.db.select():
            xtal = self.get_pyxtal(id=row.id)
            if xtal is not None:
                xtals.append(xtal)
        return xtals

    def add_xtal(self, xtal, kvp):
        """
        Add new xtal to the given db
        """
        spg_num = xtal.group.number
        density = xtal.get_density()
        dof = xtal.get_dof()
        wps = [s.wp.get_label() for s in xtal.atom_sites]
        _kvp = {"space_group_number": spg_num,
                "wps": str(wps),
                "density": density,
                "dof": dof,
               }
        kvp.update(_kvp)
        atoms = xtal.to_ase(resort=False)
        self.db.write(atoms, key_value_pairs=kvp)

    def add_strucs_from_db(self, db_file, check=False, tol=0.1, freq=50):
        """
        Add new structures from the given db_file

        Args:
            db_file (str): database path
            tol (float): tolerance in angstrom for symmetry detection
            freq (int): print frequency
        """
        print("\nAdding new strucs from {:s}".format(db_file))

        count = 0
        with connect(db_file) as db:
            for row in db.select():
                atoms = row.toatoms()
                xtal = pyxtal()
                try:
                    xtal.from_seed(atoms, tol=tol)
                except:
                    xtal = None
                if xtal is not None and xtal.valid:
                    if check:
                        add = self.check_new_structure(xtal)
                    else:
                        add = True
                    if add:
                        kvp = {}
                        for key in self.keys:
                            if hasattr(row, key):
                                kvp[key] = getattr(row, key)
                            elif key == 'space_group_number':
                                kvp[key] = xtal.group.number
                            elif key == 'density':
                                kvp[key] = xtal.get_density()
                            elif key == 'dof':
                                kvp[key] = xtal.get_dof()
                            elif key == 'wps':
                                kvp[key] == str(s.wp.get_label() for s in xtal.atom_sites)
                        self.db.write(atoms, key_value_pairs=kvp)
                        count += 1

                if count % freq == 0:
                    print("Adding {:4d} strucs from {:s}".format(count, db_file))


    def check_new_structure(self, xtal, same_group=True):
        """
        Check if the input xtal already exists in the db

        Args:
            xtal: pyxtal object
            same_group (bool): keep the same group or not
        """

        s_pmg = xtal.to_pymatgen()
        for row in self.db.select():
            if same_group:
                if row.space_group_number != xtal.group.number:
                    continue
            ref = self.db.get_atoms(id=row.id)
            ref_pmg = ase2pymatgen(ref)
            if self.matcher.fit(s_pmg, ref_pmg, symmetric=True):
                return False
        return True

    def clean_structures_spg_topology(self):
        """
        Clean up the db by removing the duplicate structures
        Here we check the follow criteria
            - same number of atoms
            - same space group
            - same topology
            - same wps

        Args:
            dtol (float): tolerance of density
            etol (float): tolerance of energy
        """

        unique_rows = []
        to_delete = []

        for row in self.db.select():
            unique = True
            for prop in unique_rows:
                (natoms, spg, wps, topology) = prop
                if natoms==row.natoms and spg == row.space_group_number \
                    and wps == row.wps:
                    if hasattr(row, 'topology'):
                        if row.topology == 'aaa':
                            if row.topology_detail == topology:
                                unique = False
                                break
                        elif row.topology == topology:
                            unique = False
                            break
            if unique:
                if hasattr(row, 'topology'):
                    unique_rows.append((row.natoms,
                                        row.space_group_number,
                                        row.wps,
                                        row.topology if row.topology !='aaa' else row.topology_detail))
                else:
                    unique_rows.append((row.natoms,
                                        row.space_group_number,
                                        row.wps,
                                        None))
            else:
                to_delete.append(row.id)
        print("The following structures were deleted", to_delete)
        self.db.delete(to_delete)


    def clean_structures(self, ids=(None, None), dtol=2e-3, etol=1e-3, criteria=None):
        """
        Clean up the db by removing the duplicate structures
        Here we check the follow criteria
            - same number of atoms
            - same density
            - same energy

        Args:
            dtol (float): tolerance of density
            etol (float): tolerance of energy
            criteria (dict): including
        """

        unique_rows = []
        to_delete = []
        ids, xtals = self.select_xtals(ids)
        for id, xtal in zip(ids, xtals):
            row = self.db.get(id)
            xtal = self.get_pyxtal(id)
            unique = True
            if criteria is not None:
                if not xtal.check_validity(criteria, True):
                    unique = False
                    print('Found unsatisfied criteria', row.id, row.space_group_number, row.wps)

                if unique:
                    if 'MAX_energy' in criteria and hasattr(row, 'ff_energy') \
                        and row.ff_energy > criteria['MAX_energy']:
                        unique = False
                        print('Unsatisfied energy', row.id, row.ff_energy, row.space_group_number, row.wps)
                if unique:
                    if 'MAX_similarity' in criteria and hasattr(row, 'similarity') \
                        and row.similarity > criteria['MAX_similarity']:
                        unique = False
                        print('Unsatisfied similarity', row.id, row.similarity, row.space_group_number, row.wps)
                if unique:
                    if 'BAD_topology' in criteria and hasattr(row, 'topology') \
                        and row.topology[:3] in criteria['BAD_topology']:
                        unique = False
                        print('Unsatisfied topology', row.id, row.topology, row.space_group_number, row.wps)
                if unique:
                    if 'BAD_dimension' in criteria and hasattr(row, 'dimension') \
                        and row.dimension in criteria['BAD_dimension']:
                        unique = False
                        print('Unsatisfied dimension', row.id, row.topology, row.space_group_number, row.wps)


            if unique:
                for prop in unique_rows:
                    (natoms, spg, wps, den, ff_energy) = prop
                    if natoms==row.natoms and spg==row.space_group_number and wps==row.wps:
                        if hasattr(row, 'ff_energy') and ff_energy is not None:
                            if abs(row.ff_energy-ff_energy) < etol:
                                unique = False
                                break
                        elif abs(den-row.density) < dtol:
                            unique = False
                            break

            if unique:
                if hasattr(row, 'ff_energy'):
                    unique_rows.append((row.natoms,
                                        row.space_group_number,
                                        row.wps,
                                        row.density,
                                        row.ff_energy))
                else:
                    unique_rows.append((row.natoms,
                                        row.space_group_number,
                                        row.wps,
                                        row.density,
                                        None))
            else:
                to_delete.append(row.id)
        print("The following structures were deleted", to_delete)
        self.db.delete(to_delete)

    def clean_structures_pmg(self, ids=(None, None), min_id=None, dtol=5e-2, criteria=None):
        """
        Clean up the db by removing the duplicate structures
        Here we check the follow criteria
            - same density
            - pymatgen check

        criteria should look like the following,
        {'CN': {'C': 3},
         'cutoff': 1.8,
         'MAX_energy': -8.00,
         #'MAX_similarity': 0.2,
         'BAD_topology': ['hcb'],
         'BAD_dimension': [0, 2],
        }

        Args:
            dtol (float): tolerance of density
            criteria (dict): including
        """

        unique_rows = []
        to_delete = []

        ids, xtals = self.select_xtals(ids)
        if min_id is None: min_id = min(ids)

        for id, xtal in zip(ids, xtals):
            row = self.db.get(id)
            xtal = self.get_pyxtal(id)
            unique = True

            if id > min_id and criteria is not None:
                if not xtal.check_validity(criteria, True):
                    unique = False
                    print('Found unsatisfied criteria', row.id, row.space_group_number, row.wps)

                if unique:
                    if 'MAX_energy' in criteria and hasattr(row, 'ff_energy') \
                        and row.ff_energy > criteria['MAX_energy']:
                        unique = False
                        print('Unsatisfied energy', row.id, row.ff_energy, row.space_group_number, row.wps)
                if unique:
                    if 'MAX_similarity' in criteria and hasattr(row, 'similarity') \
                        and row.similarity > criteria['MAX_similarity']:
                        unique = False
                        print('Unsatisfied similarity', row.id, row.similarity, row.space_group_number, row.wps)
                if unique:
                    if 'BAD_topology' in criteria and hasattr(row, 'topology') \
                        and row.topology[:3] in criteria['BAD_topology']:
                        unique = False
                        print('Unsatisfied topology', row.id, row.topology, row.space_group_number, row.wps)
                if unique:
                    if 'BAD_dimension' in criteria and hasattr(row, 'dimension') \
                        and row.dimension in criteria['BAD_dimension']:
                        unique = False
                        print('Unsatisfied dimension', row.id, row.topology, row.space_group_number, row.wps)


            if unique and id > min_id:
                for prop in unique_rows:
                    (rowid, den) = prop
                    if abs(den-row.density) < dtol:
                        ref_pmg = xtal.to_pymatgen()
                        s_pmg = ase2pymatgen(self.db.get_atoms(id=rowid))
                        if self.matcher.fit(s_pmg, ref_pmg):#, symmetric=True):
                            print('Found duplicate', row.id, row.space_group_number, row.wps)
                            unique = False
                            break
            if unique:
                unique_rows.append((row.id, row.density))
            else:
                to_delete.append(row.id)
        print("The following structures were deleted", to_delete)
        self.db.delete(to_delete)


    def select_xtals(self, ids, overwrite=False, attribute=None, use_ff=False):
        """
        Extract xtals based on attribute name.
        Mostly called by update_row_ff_energy or update_row_dftb_energy.

        Args:
            ids:
            overwrite:
        """
        (min_id, max_id) = ids
        if min_id is None: min_id = 1
        if max_id is None: max_id = self.db.count() + 10000

        ids, xtals = [], []
        for row in self.db.select():
            if overwrite or attribute is None or not hasattr(row, attribute):
                if min_id <= row.id <= max_id:
                    xtal = self.get_pyxtal(row.id, use_ff)
                    ids.append(row.id)
                    xtals.append(xtal)
                    if len(xtals) % 100 == 0:
                        print('Loading xtals from db', len(xtals))
        return ids, xtals

    def update_row_ff_energy(self, ff='reaxff', ids=(None, None), ncpu=1,
                             calc_folder='gulp_calc',
                             criteria=None,
                             overwrite=False):
        """
        Update row ff_energy with GULP calculator

        Args:
            ff (str): GULP force field library (e.g., 'reaxff', 'tersoff')
            ids (tuple): row ids e.g., (0, 100)
            ncpu (int): number of parallel processes
            calc_folder (str): temporary folder for GULP calculations
            overwrite (bool): remove the existing attributes
        """

        os.makedirs(calc_folder, exist_ok=True)
        ids, xtals = self.select_xtals(ids, overwrite, 'ff_energy')

        gulp_results = []

        # Serial or Parallel computation
        if ncpu == 1:
            for id, xtal in zip(ids, xtals):
                res = gulp_opt_single(id, xtal, ff, calc_folder, criteria)
                (xtal, eng, status) = res
                if status: gulp_results.append((id, xtal, eng))
        else:
            if len(ids) < ncpu: ncpu = len(ids)
            N_cycle = int(np.ceil(len(ids)/ncpu))
            print("\n# Parallel GULP optimizations", ncpu, N_cycle, len(ids))
            args_list = []

            for i in range(ncpu):
                id1 = i * N_cycle
                id2 = min([id1 + N_cycle, len(ids)])
                args_list.append((ids[id1:id2], xtals[id1:id2], ff, calc_folder, criteria))

            with ProcessPoolExecutor(max_workers=ncpu) as executor:
                results = [executor.submit(gulp_opt_par, *p) for p in args_list]
                for result in results:
                    gulp_results.extend(result.result())

        print("Wrap up the final results and update db", len(gulp_results))
        for gulp_result in gulp_results:
            (id, xtal, eng) = gulp_result
            if xtal is not None:
                self.db.update(id,
                               ff_energy=eng,
                               ff_lib=ff,
                               ff_relaxed=xtal.to_file())

    def update_row_dftb_energy(self, skf_dir, steps=500,
                               ids=(None, None),
                               use_ff = True,
                               ncpu=1,
                               calc_folder='dftb_calc',
                               criteria=None,
                               symmetrize=False,
                               overwrite=False):
        """
        Update row ff_energy with GULP calculator

        Args:
            skf_dir (str): GULP force field library (e.g., 'reaxff', 'tersoff')
            steps (int): relaxation steps
            ids (tuple): row ids e.g., (0, 100)
            ncpu (int): number of parallel processes
            calc_folder (str): temporary folder for GULP calculations
            symmetrize (bool): impose symmetry in optimization
            overwrite (bool): remove the existing attributes
        """

        os.makedirs(calc_folder, exist_ok=True)
        ids, xtals = self.select_xtals(ids, overwrite, 'dftb_energy', use_ff)

        dftb_results = []
        os.chdir(calc_folder)

        # Serial or Parallel computation
        if ncpu == 1:
            for id, xtal in zip(ids, xtals):
                res = dftb_opt_single(id, xtal, skf_dir, steps, symmetrize, criteria)
                (xtal, eng, status) = res
                if status: dftb_results.append((id, xtal, eng))
        else:
            # reset ncpus if the cpu is greater than the actual number of jobs
            if len(ids) < ncpu: ncpu = len(ids)
            N_cycle = int(np.ceil(len(ids)/ncpu))
            print("\n# Parallel DFTB optimizations", ncpu, N_cycle, len(ids))
            args_list = []

            for i in range(ncpu):
                id1 = i * N_cycle
                id2 = min([id1 + N_cycle, len(ids)])
                folder = self.get_label(i)
                os.makedirs(folder, exist_ok=True)
                args_list.append((ids[id1:id2],
                                  xtals[id1:id2],
                                  skf_dir,
                                  steps,
                                  folder,
                                  symmetrize,
                                  criteria))

            with ProcessPoolExecutor(max_workers=ncpu) as executor:
                results = [executor.submit(dftb_opt_par, *p) for p in args_list]
                for result in results:
                    dftb_results.extend(result.result())

        # Wrap up the final results and update db
        for dftb_result in dftb_results:
            (id, xtal, eng) = dftb_result
            self.db.update(id,
                           dftb_energy=eng,
                           dftb_relaxed=xtal.to_file())


    def update_row_topology(self, StructureType='Auto', overwrite=True):
        """
        Update row topology base on the CrystalNets.jl

        Args:
            StructureType (str): 'Zeolite', 'MOF' or 'Auto'
            overwrite (bool): remove the existing attributes
        """
        try:
            import juliacall
        except:
            print("Cannot load JuliaCall, no support on topology parser")

        def parse_topology(topology_info):
            """
            Obtain the dimensionality and topology name
            """
            dim = 0
            name = ""
            detail = 'None'
            for i, x in enumerate(topology_info):
                (d, n) = x
                if d > dim: dim = d
                tmp = n.split(',')[0]
                if tmp.startswith("UNKNOWN"):
                    detail = tmp[7:] #tuple(int(num) for num in tmp[7:].split())
                    tmp = 'aaa'
                elif tmp.startswith("unstable"):
                    tmp = 'unstable'
                name += tmp
                if i + 1 < len(topology_info): name += '-'
            return dim, name, detail

        jl = juliacall.newmodule("MOF_Builder")
        jl.seval("using CrystalNets")
        jl.CrystalNets.toggle_warning(False) # to disable warnings
        jl.CrystalNets.toggle_export(False) # to disable exports
        if StructureType == 'Zeolite':
            option = jl.CrystalNets.Options(structure=jl.StructureType.Zeolite)
        elif StructureType == 'MOF':
            option = jl.CrystalNets.Options(structure=jl.StructureType.MOF)
        else:
            option = jl.CrystalNets.Options(structure=jl.StructureType.Auto)

        for row in self.db.select():
            if overwrite or not hasattr(row, 'topology'):
                atoms = self.db.get_atoms(row.id)
                atoms.write('tmp.cif', format='cif')

                # Call crystalnet.jl
                result = jl.determine_topology('tmp.cif', option)
                #print(result)
                if len(result) > 1:
                    results = [x for x in result]
                else:
                    results = [result[0]]
                try:
                    topo = []
                    for res in results:
                        # topology for SingleNodes
                        name = str(res[0])
                        if res[1] > 1: name += '(' + str(res[1]) + ')'
                        genome = res[0][jl.Clustering.Auto]
                        dim = jl.ndims(jl.CrystalNets.PeriodicGraph(genome))
                        topo.append((dim, name))

                    # The maximum dimensionality and topology name
                    dim, name, detail = parse_topology(topo)
                except:
                    dim, name, detail = 3, 'error', 'error'
                print("Updating Topology", row.space_group_number, row.wps, dim, name, detail[:10])
                # Unknown will be labeled as aaa
                self.db.update(row.id, topology=name, dimension=dim, topology_detail=detail)
            else:
                print("Existing Topology", row.topology)

    def update_db_description(self):
        """
        update db description based on robocrys
        Call robocrys: https://github.com/hackingmaterials/robocrystallographer
        """
        from robocrys import StructureCondenser, StructureDescriber

        condenser = StructureCondenser()
        describer = StructureDescriber()

        for row in self.db.select():
            if not hasattr(row, 'description'):
                atoms = self.db.get_atoms(row.id)
                pmg = ase2pymatgen(atoms)
                try:
                    condensed_structure = condenser.condense_structure(pmg)
                    description = describer.describe(condensed_structure)
                except:
                    description = 'N/A'

                self.db.update(row.id, description=description)
                print("\n======Updating\n", description)
            else:
                print("\n======Existing\n", row.description)

    def export_structures(self, fmt='vasp', folder='mof_out', criteria=None,
                          sort_by='similarity', overwrite=True):
        """
        export structures from database according to the given criterion

        Args:
            fmt (str): 'vasp' or 'cif'
            folder (str): 'path of output folders'
            criteria (dict): check the validity with dict
            sort_by (str): sort by which attribute
            overwrite (bool): remove the existing folder
        """

        import shutil
        if not os.path.exists(folder):
            os.makedirs(folder)
        else:
            if overwrite:
                shutil.rmtree(folder)
                os.makedirs(folder)

        keys = ['space_group_number',
                'density',
                'dof',
                'similarity',
                'ff_energy',
                'topology',
               ]
        properties, atoms = [], []
        for row in self.db.select():
            spg = row.space_group_number
            den = row.density
            dof = row.dof
            sim = row.similarity if hasattr(row, 'similarity') else None
            top = row.topology if hasattr(row, 'topology') else None
            eng = row.ff_energy if hasattr(row, 'ff_energy') else None
            properties.append([row.id, spg, den, dof, sim, eng, top])
            atoms.append(self.db.get_atoms(id=row.id))

        if sort_by in keys:
            col = keys.index(sort_by) + 1
        else:
            print("supported attributes", keys)
            raise ValueError("Cannot sort by", sort_by)

        print("====Exporting {:} structures".format(len(atoms)))
        if len(atoms)>0:
            properties = np.array(properties)
            mids = np.argsort(properties[:, col])

            for mid in mids:
                [id, spg, den, dof, sim, eng, top] = properties[mid]
                id = int(id)
                spg = int(spg)
                sim = float(sim)
                den = float(den)
                dof = int(dof)
                if eng is not None: eng = float(eng)
                try:
                    xtal = pyxtal()
                    xtal.from_seed(atoms[mid])
                    number, symbol = xtal.group.number, xtal.group.symbol.replace('/','')
                    # convert to the desired subgroup representation if needed
                    if number != spg:
                        paths = xtal.group.path_to_subgroup(spg)
                        xtal = xtal.to_subgroup(paths)
                        number, symbol = xtal.group.number, xtal.group.symbol.replace('/','')

                    label = os.path.join(folder, '{:d}-{:d}-{:s}-N{:d}'.format(id, number, symbol, sum(xtal.numIons)))

                    if criteria is not None:
                        status = xtal.check_validity(criteria, True)
                    else:
                        status = True
                except:
                    status = False

                if status:
                    #if top is not None: print(top)
                    try:
                        xtal.set_site_coordination()
                        for s in xtal.atom_sites:
                            _l, _sp, _cn = s.wp.get_label(), s.specie, s.coordination
                            label += '-{:s}-{:s}{:d}'.format(_l, _sp, _cn)
                        label += '-S{:.3f}'.format(sim)
                    except:
                        print('Problem in setting site coordination')
                    if len(label) > 40: label = label[:40]

                    if den is not None: label += '-D{:.2f}'.format(abs(den))
                    if eng is not None: label += '-E{:.3f}'.format(abs(eng))
                    if top is not None: label += '-T{:s}'.format(top)
                    if sim is not None: label += '-S{:.2f}'.format(sim)

                    print("====Exporting:", label)
                    if fmt == 'vasp':
                        atoms[mid].write(label+'.vasp', format='vasp', vasp5=True, direct=True)
                    elif fmt == 'cif':
                        xtal.to_file(label+'.cif')
                else:
                    print("====Skippng:  ", label)

    def get_label(self, i):
        if i < 10:
            folder = f"cpu00{i}"
        elif i < 100:
            folder = f"cpu0{i}"
        else:
            folder = f"cpu0{i}"
        return folder


if __name__ == "__main__":
    # open
    if False:
        db = database('test.db')
        print("Total number of entries", len(db.codes))

        # view structure
        c = db.get_pyxtal('HXMTAM')
        print(c)
    if True:

        db = database_topology('../MOF-Builder/C-sp2/sp2-sacada-0506.db')
        #xtal = db.get_pyxtal(1)
        #print(xtal)
        #db.add_xtal(xtal, kvp={'similarity': 0.1})

        #db.update_row_ff_energy(ids=(0, 2), overwrite=True)
        #db.update_row_ff_energy(ncpu=2, ids=(2, 20), overwrite=True)
        # brew install coreutils to get timeout in maca
        #os.environ['ASE_DFTB_COMMAND'] = 'timeout 1m /Users/qzhu8/opt/dftb+/bin/dftb+ > PREFIX.out'
        os.environ['ASE_DFTB_COMMAND'] = '/Users/qzhu8/opt/dftb+/bin/dftb+ > PREFIX.out'
        skf_dir = '/Users/qzhu8/GitHub/MOF-Builder/3ob-3-1/'
        #db.update_row_dftb_energy(skf_dir, ncpu=1, ids=(0, 2), overwrite=True)
        db.update_row_dftb_energy(skf_dir, ncpu=1, ids=(17, 17), overwrite=True)
