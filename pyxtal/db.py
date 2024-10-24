"""
Database class
"""

import os, time
from concurrent.futures import ProcessPoolExecutor
import logging

import numpy as np
import pymatgen.analysis.structure_matcher as sm
from ase.calculators.calculator import CalculationFailed
from ase.db import connect

from pyxtal import pyxtal
from pyxtal.util import ase2pymatgen


def setup_worker_logger(log_file):
    """
    Set up the logger for each worker process.
    """
    logging.getLogger().handlers.clear()
    logging.basicConfig(format="%(asctime)s| %(message)s",
                        filename=log_file,
                        level=logging.INFO)

def call_opt_single(p):
    """
    Optimize a single structure and log the result.

    Args:
        p (tuple): A tuple where the first element is an identifier (id), and the
                   remaining elements are the arguments to pass to `opt_single`.

    Returns:
        tuple: A tuple (id, xtal, eng) where:
               - id (int): The identifier of the structure.
               - xtal: The optimized structure.
               - eng (float): The energy of the opt_structure, or None if it failed.

    Behavior:
        This function calls `opt_single` to perform the optimization of the structure
        associated with the given id.
    """
    #logger = logging.getLogger()
    #logger.info(f"ID: {p[0]} *{sum(p[1].numIons)}")
    myid = p[0]
    xtal, eng, status = opt_single(*p)
    return myid, xtal, eng


def opt_single(id, xtal, calc, *args):
    """
    Optimize a structure using the specified calculator.

    Args:
        id (int): Identifier of the structure to be optimized.
        xtal: Crystal structure object to be optimized.
        calc (str): The calculator to use ('GULP', 'DFTB', 'VASP', 'MACE').
        *args: Additional arguments to pass to the calculator function.

    Returns:
        tuple: The result of the optimization, which typically includes:
               - xtal: The optimized structure.
               - energy (float): The energy of the optimized structure.
               - status (bool): Whether the optimization was successful.

    Raises:
        ValueError: If an unsupported calculator is specified.
    """

    if calc == 'GULP':
        return gulp_opt_single(id, xtal, *args)
    elif calc == 'DFTB':
        return dftb_opt_single(id, xtal, *args)
    elif calc == 'VASP':
        return vasp_opt_single(id, xtal, *args)
    elif calc == 'MACE':
        return mace_opt_single(id, xtal, *args)
    else:
        raise ValueError("Cannot support this calcultor", calc)


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
    from pyxtal.interface.dftb import DFTB, DFTB_relax

    cwd = os.getcwd()
    atoms = xtal.to_ase(resort=False)
    eng = None
    stress = None
    try:
        if symmetrize:
            s = DFTB_relax(
                atoms,
                skf_dir,
                True,
                int(steps / 2),
                kresol=kresol,
                folder=".",
                scc_error=1e-5,
                scc_iter=100,
                logfile="ase.log",
            )
            stress = np.sum(s.get_stress()[:3]) / 0.006241509125883258 / 3
        else:
            my = DFTB(
                atoms,
                skf_dir,
                kresol=kresol * 1.5,
                folder=".",
                scc_error=0.1,
                scc_iter=100,
            )
            s, eng = my.run(mode="vc-relax", step=int(steps / 2))
            my = DFTB(s, skf_dir, kresol=kresol, folder=".",
                      scc_error=1e-4, scc_iter=100)
            s, eng = my.run(mode="vc-relax", step=int(steps / 2))
            s = my.struc
    except CalculationFailed:
        # This is due to covergence error in geometry optimization
        # Here we simply read the last energy
        my.calc.read_results()
        eng = my.calc.results["energy"]
        s = my.struc
    except:
        s = None
        print("Problem in DFTB Geometry optimization", id)
        xtal.to_file("bug.cif")  # ; import sys; sys.exit()
    os.chdir(cwd)

    status = False
    if s is not None:
        c = pyxtal()
        c.from_seed(s)
        if eng is None:
            eng = s.get_potential_energy() / len(s)
        else:
            eng /= len(s)

        status = process_xtal(id, xtal, eng, criteria)
        print(xtal.get_xtal_string())

        return xtal, eng, status
    else:
        return None, None, False


def vasp_opt_single(id, xtal, path, cmd, criteria):
    """
    Single VASP optimization for a given atomic xtal

    Args:
        id (int): id of the give xtal
        xtal: pyxtal instance
        path: calculation folder
        cmd: vasp command
        criteria (dicts): to check if the structure
    """
    from pyxtal.interface.vasp import optimize as vasp_opt
    cwd = os.getcwd()
    path += '/g' + str(id)
    status = False

    xtal, eng, _, error = vasp_opt(xtal,
                                   path,
                                   cmd=cmd,
                                   walltime="59m")
    if not error:
        status = process_xtal(id, xtal, eng, criteria)
    else:
        os.chdir(cwd)
    return xtal, eng, status


def gulp_opt_single(id, xtal, ff_lib, path, criteria):
    """
    Perform a single GULP optimization for a given crystal structure.

    Args:
        id (int): Identifier for the current structure.
        xtal: PyXtal instance representing the crystal to be optimized.
        ff_lib (str): Force field library for GULP, e.g., 'reaxff', 'tersoff'.
        path (str): Path to the folder where the calculation is stored.
        criteria (dict): Dictionary to check the validity of the opt_structure.

    Returns:
        tuple:
            - xtal: Optimized PyXtal instance.
            - eng (float): Energy of the optimized structure.
            - status (bool): Whether the optimization process is successful.

    Behavior:
        This function performs a GULP optimization using the force field.
        After the optimization, it checks the validity of the structure and
        attempts to remove the calculation folder if it is empty.
    """
    from pyxtal.interface.gulp import single_optimize as gulp_opt

    # Create the path for this specific structure
    path += '/g' + str(id)

    # Perform the optimization with GULP
    xtal, eng, _, error = gulp_opt(
        xtal,
        ff=ff_lib,
        label=str(id),
        path=path,
        symmetry=True,
    )

    # Default status to False, will be updated if successful
    status = False
    if not error:
        status = process_xtal(id, xtal, eng, criteria)
        try:
            os.rmdir(path)
        except:
            print("Folder is not empty", path)
    return xtal, eng, status


def mace_opt_single(id, xtal, step, criteria):
    """
    Perform a single MACE optimization for a given atomic crystal structure.

    Args:
        id (int): Identifier for the current structure.
        xtal: PyXtal instance representing the crystal structure.
        step (int): Maximum number of relaxation steps. Default is 250.
        criteria (dict): Dictionary to check the validity of the optimized structure.

    Returns:
        tuple:
            - xtal: Optimized PyXtal instance (or None if optimization failed).
            - eng (float): Energy/atom of the opt_structure (or None if it failed).
            - status (bool): Whether the optimization was successful.
    """
    from pyxtal.interface.ase_opt import ASE_relax as mace_opt

    logger = logging.getLogger()
    atoms = xtal.to_ase(resort=False)
    s = mace_opt(atoms,
                 'MACE',
                 opt_cell=True,
                 step=step,
                 max_time=9.0 * max([1, (len(atoms)/200)]),
                 label=str(id))
    if s is None:
        logger.info(f"mace_opt_single Failure {id}")
        return None, None, False

    try:
        xtal = pyxtal()
        xtal.from_seed(s)
        eng = s.get_potential_energy() / len(s)
        status = process_xtal(id, xtal, eng, criteria)
        logger.info(f"mace_opt_single Success {id}")
        return xtal, eng, status
    except:
        logger.info(f"mace_opt_single Bug {id}")
        return None, None, False


def process_xtal(id, xtal, eng, criteria):
    status = xtal.check_validity(
        criteria) if criteria is not None else True
    if status:
        header = f"{id:4d}"
        dicts = {"validity": status, "energy": eng}
        print(xtal.get_xtal_string(header=header, dicts=dicts))
    return status


def make_entry_from_pyxtal(xtal):
    """
    Generate an entry dictionary from a PyXtal object, assuming
    the SMILES and CCDC number information is provided.

    Args:
        xtal: PyXtal object (must contain the SMILES (`xtal.tag["smiles"]`)
        and CCDC number (`xtal.tag["ccdc_number"]`) in the `xtal.tag`.

    Returns:
        tuple: (ase_atoms, entry_dict, None)
            - ase_atoms: ASE Atoms object converted from the PyXtal structure.
            - entry_dict (dict): A dictionary containing information
            - None: Placeholder for future use (currently returns None).

    Structure of `entry_dict`:
        - "csd_code" (str): CSD code (if available) for the crystal structure.
        - "mol_smi" (str): SMILES representation of the molecule.
        - "ccdc_number" (str): CCDC identifier number.
        - "space_group" (str): Space group symbol of the crystal.
        - "spg_num" (int): Space group number.
        - "Z" (int): Number of molecules in the unit cell.
        - "Zprime" (float): Z' value of the crystal.
        - "url" (str): URL link to the CCDC database entry for the crystal.
        - "mol_formula" (str): Molecular formula of the structure.
        - "mol_weight" (float): Molecular weight of the structure.
        - "mol_name" (str): Name of the molecule, typically the CSD code.
        - "l_type" (str): Lattice type of the structure.

    Returns None if the PyXtal structure is invalid (i.e., `xtal.valid` is False).

    Example:
        entry = make_entry_from_pyxtal(xtal_instance)
        ase_atoms, entry_dict, _ = entry

    Notes:
        - The CCDC link is generated using the structure's CCDC number.
    """

    from rdkit import Chem
    from rdkit.Chem.Descriptors import ExactMolWt
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    if xtal.valid:
        url0 = "https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid="
        # Create RDKit molecule from SMILES string
        m = Chem.MolFromSmiles(xtal.tag["smiles"])

        # Calculate molecular weight and molecular formula using RDKit
        mol_wt = ExactMolWt(m)
        mol_formula = CalcMolFormula(m)

        # Create a dictionary containing information
        kvp = {
            "csd_code": xtal.tag["csd_code"],
            "mol_smi": xtal.tag["smiles"],
            "ccdc_number": xtal.tag["ccdc_number"],
            "space_group": xtal.group.symbol,
            "spg_num": xtal.group.number,
            "Z": sum(xtal.numMols),
            "Zprime": xtal.get_zprime()[0],
            "url": url0 + str(xtal.tag["ccdc_number"]),
            "mol_formula": mol_formula,
            "mol_weight": mol_wt,
            "mol_name": xtal.tag["csd_code"],
            "l_type": xtal.lattice.ltype,
        }
        # Return the ASE Atoms the entry dictionary, and None as a placeholder
        return (xtal.to_ase(), kvp, None)
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

    # xtal = pyxtal(molecular=True)
    #
    # return make_entry_from_pyxtal(xtal)
    raise NotImplementedError("To do in future")


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


class database:
    """
    This is a database class to process crystal data.

    Args:
        db_name: `*.db` format from ase database
    """

    def __init__(self, db_name):
        self.db_name = db_name
        # if not os.path.exists(db_name):
        #    raise ValueError(db_name, 'doesnot exist')

        self.db = connect(db_name, serial=True)
        self.codes = self.get_all_codes()
        self.keys = [
            "csd_code",
            "space_group",
            "spg_num",
            "Z",
            "Zprime",
            "url",
            "mol_name",
            "mol_smi",
            "mol_formula",
            "mol_weight",
            "l_type",
        ]
        self.calculators = [
            "charmm",
            "gulp",
            "ani",
            "dftb_D3",
            "dftb_TS",
        ]

    def vacuum(self):
        self.db.vacuum()

    def get_all_codes(self, group=None):
        """
        Get all codes
        """
        codes = []
        for row in self.db.select():
            if row.csd_code not in codes:
                if group is None:
                    codes.append(row.csd_code)
                else:
                    if row.group == group:
                        codes.append(row.csd_code)
            else:
                print("find duplicate! remove", row.id, row.csd_code)
                self.db.delete([row.id])
        return codes
        # self.codes = codes

    def add(self, entry):
        (atom, kvp, data) = entry
        if kvp["csd_code"] not in self.codes:
            kvp0 = self.process_kvp(kvp)
            self.db.write(atom, key_value_pairs=kvp0, data=data)
            self.codes.append(kvp["csd_code"])

    def add_from_code(self, code):
        entry = make_entry_from_CSD(code)
        if entry is not None:
            self.add(entry)
        else:
            print(f"{code:s} is not a valid entry")

    def process_kvp(self, kvp):
        kvp0 = {}
        for key in self.keys:
            if key in kvp:
                kvp0[key] = kvp[key]
            else:
                print("Error, cannot find ", key, " from the input")
                return None
        return kvp0

    def check_status(self, show=False):
        """
        Check the current status of each entry
        """
        ids = []
        for row in self.db.select():
            if len(row.data.keys()) == len(self.calculators):
                ids.append(row.id)
                if show:
                    row_info = self.get_row_info(id=row.id)
                    self.view(row_info)
            else:
                print(row.csd_code)  # , row.data['charmm_info']['prm'])
        return ids

    def copy(self, db_name, csd_codes):
        """
        copy the entries to another db

        Args:
            db_name: db file name
            csd_codes: list of codes
        """
        if db_name == self.db_name:
            raise RuntimeError("Cannot use the same db file for copy")
        with connect(db_name, serial=True) as db:
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
        from pyxtal.representation import representation

        (atom, kvp, data) = row_info

        # Reference
        xtal = self.get_pyxtal(kvp["csd_code"])
        rep = xtal.get_1D_representation()
        print("\n", kvp["csd_code"], kvp["mol_smi"], xtal.lattice.volume)
        print(rep.to_string() + " reference")

        # calcs
        for key in data:
            calc = key[:-5]
            time = data[key]["time"]

            rep = data[key]["rep"]
            if type(rep[0]) is not list:
                rep = [rep]
            rep = representation(rep, kvp["mol_smi"]).to_string()

            (dv, msd1, msd2) = data[key]["diff"]
            strs = f"{rep:s} {calc:8s} {time / 60:6.2f} {dv:6.3f}"
            if msd1 is not None:
                strs += f"{msd1:6.3f}{msd2:6.3f}"
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
                key = calc + "_info"
                if key in row.data:
                    data0[key] = row.data[key]

            atom = self.db.get_atoms(id=row.id)
            return (atom, kvp, data0)
        else:
            msg = "cannot find the entry from " + id + code
            raise RuntimeError(msg)

    def get_row(self, code):
        for row in self.db.select(csd_code=code):
            return row
        msg = "cannot find the entry from " + code
        raise RuntimeError(msg)

    def get_pyxtal(self, code):
        from pyxtal import pyxtal
        from pyxtal.msg import ReadSeedError
        from pyxtal.util import ase2pymatgen

        row = self.get_row(code)
        atom = self.db.get_atoms(id=row.id)
        # Reference
        pmg = ase2pymatgen(atom)
        smi = row.mol_smi
        smiles = smi.split(".")
        molecules = [smile + ".smi" for smile in smiles]

        xtal = pyxtal(molecular=True)
        try:
            xtal.from_seed(pmg, molecules=molecules)
        except ReadSeedError:
            xtal.from_seed(pmg, molecules=molecules, add_H=True)

        return xtal

    def compute(self, row, work_dir, skf_dir):
        if len(row.data.keys()) < len(self.calculators):
            # not label information, run antechamber
            atom = self.db.get_atoms(id=row.id)
            if "gulp_info" not in row.data:
                # pmg, c_info, g_info = get_parameters(row, atom)
                # row.data = {"charmm_info": c_info, "gulp_info": g_info}
                pass
            else:
                pmg = ase2pymatgen(atom)

            #data = compute(row, pmg, work_dir, skf_dir)
            #self.db.update(row.id, data=data)
            #print("updated the data for", row.csd_code)


class database_topology:
    """
    This is a database class to process atomic crystal data

    Args:
        db_name (str): `*.db` format from ase database
        rank (int): default 0
        size (int): default 1
        ltol (float): lattice tolerance
        stol (float): site tolerance
        atol (float): angle tolerance
        log_file (str): log_file
    """

    def __init__(self, db_name, rank=0, size=1, ltol=0.05, stol=0.05, atol=3,
                 log_file='db.log'):
        self.rank = rank
        self.size = size
        self.db_name = db_name
        self.db = connect(db_name, serial=True)
        self.keys = [
            "space_group_number",
            "pearson_symbol",
            "similarity0",
            "similarity",
            "density",
            "dof",
            "topology",
            "topology_detail",
            "dimension",
            "wps",
            "ff_energy",
            "ff_lib",
            "ff_relaxed",
            "mace_energy",
            "mace_relaxed",
            "dftb_energy",
            "dftb_relaxed",
            "vasp_energy",
            "vasp_relaxed",
        ]
        self.matcher = sm.StructureMatcher(
            ltol=ltol, stol=stol, angle_tol=atol)

        # Define logfile
        self.log_file = log_file
        logging.getLogger().handlers.clear()
        logging.basicConfig(format="%(asctime)s| %(message)s",
                            filename=self.log_file,
                            level=logging.INFO)
        self.logging = logging

    def vacuum(self):
        self.db.vacuum()

    def print_memory_usage(self):
        import psutil
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss / 1024 ** 2
        self.logging.info(f"Rank {self.rank} memory: {mem:.1f} MB")
        print(f"Rank {self.rank} memory: {mem:.1f} MB")

    def get_pyxtal(self, id, use_relaxed=None):
        """
        Get pyxtal based on row_id, if use_relaxed, get pyxtal from ff_relaxed

        Args:
            id (int): row id
            use_relaxed (str): 'ff_relaxed', 'vasp_relaxed'
        """
        from pymatgen.core import Structure
        from pyxtal import pyxtal
        from pyxtal.util import ase2pymatgen

        row = self.db.get(id)  # ; print(id, row.id)
        if use_relaxed is not None:
            if hasattr(row, use_relaxed):
                xtal_str = getattr(row, use_relaxed)
                pmg = Structure.from_str(xtal_str, fmt="cif")
            else:
                print(f"No {use_relaxed} attributes for structure", id)
                atom = self.db.get_atoms(id=id)
                pmg = ase2pymatgen(atom)
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
            print("Cannot load the structure")

    def get_all_xtals(self, include_energy=False):
        """
        Get all pyxtal instances from the current db
        """
        xtals = []
        for row in self.db.select():
            xtal = self.get_pyxtal(id=row.id)
            if xtal is not None:
                if include_energy and hasattr(row, 'vasp_energy'):
                    xtal.energy = row.vasp_energy
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
        _kvp = {
            "space_group_number": spg_num,
            "pearson_symbol": xtal.get_Pearson_Symbol(),
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
        print(f"\nAdding new strucs from {db_file:s}")

        count = 0
        with connect(db_file, serial=True) as db:
            for row in db.select():
                atoms = row.toatoms()
                xtal = pyxtal()
                try:
                    xtal.from_seed(atoms, tol=tol)
                except:
                    xtal = None
                if xtal is not None and xtal.valid:
                    add = self.check_new_structure(xtal) if check else True
                    if add:
                        kvp = {}
                        for key in self.keys:
                            if hasattr(row, key):
                                kvp[key] = getattr(row, key)
                            elif key == "space_group_number":
                                kvp[key] = xtal.group.number
                            elif key == "density":
                                kvp[key] = xtal.get_density()
                            elif key == "dof":
                                kvp[key] = xtal.get_dof()
                            elif key == "wps":
                                kvp[key] == str(s.wp.get_label()
                                                for s in xtal.atom_sites)
                            elif key == "pearson_symbol":
                                kvp[key] = xtal.get_Pearson_Symbol()

                        self.db.write(atoms, key_value_pairs=kvp)
                        count += 1

                if count % freq == 0:
                    print(f"Adding {count:4d} strucs from {db_file:s}")

    def check_new_structure(self, xtal, same_group=True):
        """
        Check if the input xtal already exists in the db

        Args:
            xtal: pyxtal object
            same_group (bool): keep the same group or not
        """

        s_pmg = xtal.to_pymatgen()
        for row in self.db.select():
            if same_group and row.space_group_number != xtal.group.number:
                continue
            ref = self.db.get_atoms(id=row.id)
            ref_pmg = ase2pymatgen(ref)
            if self.matcher.fit(s_pmg, ref_pmg, symmetric=True):
                return False
        return True

    def clean_structures_spg_topology(self, dim=None):
        """
        Clean up the db by removing the duplicate structures
        Here we check the follow criteria
            - same number of atoms
            - same space group
            - same topology
            - same wps

        Args:
            dim (int): wanted dimension
        """

        unique_rows = []
        to_delete = []

        for row in self.db.select():
            unique = True
            # Ignore unwanted dimension
            if dim is not None and hasattr(row, "dimension") and row.dimension != dim:
                # print(row.dimension, dim)
                unique = False
            else:
                for prop in unique_rows:
                    (natoms, spg, wps, topology) = prop
                    if (natoms == row.natoms and spg == row.space_group_number and wps == row.wps) and hasattr(
                        row, "topology"
                    ):
                        if row.topology == "aaa":
                            if row.topology_detail == topology:
                                unique = False
                                break
                        elif row.topology == topology:
                            unique = False
                            break
            if unique:
                if hasattr(row, "topology"):
                    unique_rows.append(
                        (
                            row.natoms,
                            row.space_group_number,
                            row.wps,
                            row.topology if row.topology != "aaa" else row.topology_detail,
                        )
                    )
                else:
                    unique_rows.append(
                        (row.natoms, row.space_group_number, row.wps, None))
            else:
                to_delete.append(row.id)
        if len(to_delete) > 0:
            print(len(to_delete), "structures were deleted", to_delete)
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
                    print(
                        "Found unsatisfied criteria",
                        row.id,
                        row.space_group_number,
                        row.wps,
                    )

                if unique and (
                    "MAX_energy" in criteria and hasattr(
                        row, "ff_energy") and row.ff_energy > criteria["MAX_energy"]
                ):
                    unique = False
                    print(
                        "Unsatisfied energy",
                        row.id,
                        row.ff_energy,
                        row.space_group_number,
                        row.wps,
                    )
                if unique and (
                    "MAX_similarity" in criteria
                    and hasattr(row, "similarity")
                    and row.similarity > criteria["MAX_similarity"]
                ):
                    unique = False
                    print(
                        "Unsatisfied similarity",
                        row.id,
                        row.similarity,
                        row.space_group_number,
                        row.wps,
                    )
                if unique and (
                    "BAD_topology" in criteria
                    and hasattr(row, "topology")
                    and row.topology[:3] in criteria["BAD_topology"]
                ):
                    unique = False
                    print(
                        "Unsatisfied topology",
                        row.id,
                        row.topology,
                        row.space_group_number,
                        row.wps,
                    )
                if unique and (
                    "BAD_dimension" in criteria
                    and hasattr(row, "dimension")
                    and row.dimension in criteria["BAD_dimension"]
                ):
                    unique = False
                    print(
                        "Unsatisfied dimension",
                        row.id,
                        row.topology,
                        row.space_group_number,
                        row.wps,
                    )

            if unique:
                for prop in unique_rows:
                    (natoms, spg, wps, den, ff_energy) = prop
                    if natoms == row.natoms and spg == row.space_group_number and wps == row.wps:
                        if hasattr(row, "ff_energy") and ff_energy is not None:
                            if abs(row.ff_energy - ff_energy) < etol:
                                unique = False
                                break
                        elif abs(den - row.density) < dtol:
                            unique = False
                            break

            if unique:
                if hasattr(row, "ff_energy"):
                    unique_rows.append(
                        (
                            row.natoms,
                            row.space_group_number,
                            row.wps,
                            row.density,
                            row.ff_energy,
                        )
                    )
                else:
                    unique_rows.append(
                        (row.natoms, row.space_group_number, row.wps, row.density, None))
            else:
                to_delete.append(row.id)
        print(len(to_delete), "structures were deleted", to_delete)
        self.db.delete(to_delete)

    def clean_structures_pmg(self, ids=(None, None), min_id=None, dtol=5e-2, criteria=None):
        """
        Clean up the db by removing the duplicate structures.
        Here we check the follow criteria same density and pymatgen matcher.
        The criteria should look like the following,

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
        if min_id is None:
            min_id = min(ids)

        for id, xtal in zip(ids, xtals):
            row = self.db.get(id)
            xtal = self.get_pyxtal(id)
            unique = True

            if id > min_id and criteria is not None:
                if not xtal.check_validity(criteria, True):
                    unique = False
                    print(
                        "Found unsatisfied criteria",
                        row.id,
                        row.space_group_number,
                        row.wps,
                    )

                if unique and (
                    "MAX_energy" in criteria and hasattr(
                        row, "ff_energy") and row.ff_energy > criteria["MAX_energy"]
                ):
                    unique = False
                    print(
                        "Unsatisfied energy",
                        row.id,
                        row.ff_energy,
                        row.space_group_number,
                        row.wps,
                    )
                if unique and (
                    "MAX_similarity" in criteria
                    and hasattr(row, "similarity")
                    and row.similarity > criteria["MAX_similarity"]
                ):
                    unique = False
                    print(
                        "Unsatisfied similarity",
                        row.id,
                        row.similarity,
                        row.space_group_number,
                        row.wps,
                    )
                if unique and (
                    "BAD_topology" in criteria
                    and hasattr(row, "topology")
                    and row.topology[:3] in criteria["BAD_topology"]
                ):
                    unique = False
                    print(
                        "Unsatisfied topology",
                        row.id,
                        row.topology,
                        row.space_group_number,
                        row.wps,
                    )
                if unique and (
                    "BAD_dimension" in criteria
                    and hasattr(row, "dimension")
                    and row.dimension in criteria["BAD_dimension"]
                ):
                    unique = False
                    print(
                        "Unsatisfied dimension",
                        row.id,
                        row.topology,
                        row.space_group_number,
                        row.wps,
                    )

            if unique and id > min_id:
                for prop in unique_rows:
                    (rowid, den) = prop
                    if abs(den - row.density) < dtol:
                        ref_pmg = xtal.to_pymatgen()
                        s_pmg = ase2pymatgen(self.db.get_atoms(id=rowid))
                        # , symmetric=True):
                        if self.matcher.fit(s_pmg, ref_pmg):
                            print(
                                "Found duplicate",
                                row.id,
                                row.space_group_number,
                                row.wps,
                            )
                            unique = False
                            break
            if unique:
                unique_rows.append((row.id, row.density))
            else:
                to_delete.append(row.id)
        print(len(to_delete), "structures were deleted", to_delete)
        self.db.delete(to_delete)

    def get_max_id(self):
        """
        Get the maximum row id
        """
        max_id = None
        for row in self.db.select():
            if max_id is None or row.id > max_id:
                max_id = row.id + 1
        return max_id

    def select_xtals(self, ids, N_atoms=(None, None), overwrite=False, attribute=None, use_relaxed=None):
        """
        Extract xtals based on attribute name.
        Mostly called by update_row_energy

        Args:
            ids:
            N_atoms:
            overwrite:
            atttribute:
            use_relaxed (str): 'ff_relaxed' or 'vasp_relaxed'
        """
        (min_id, max_id) = ids
        if min_id is None: min_id = 1
        if max_id is None: max_id = self.get_max_id()

        (min_atoms, max_atoms) = N_atoms
        if min_atoms is None: min_atoms = 1
        if max_atoms is None: max_atoms = 5000

        ids, xtals = [], []
        for row in self.db.select():
            if overwrite or attribute is None or not hasattr(row, attribute):
                if min_id <= row.id <= max_id and min_atoms < natoms <= max_atoms:
                    xtal = self.get_pyxtal(row.id, use_relaxed)
                    ids.append(row.id)
                    xtals.append(xtal)
                    if len(xtals) % 100 == 0:
                        print("Loading xtals from db", len(xtals))
        return ids, xtals

    def select_xtal(self, ids, N_atoms=(None, None), overwrite=False, attribute=None, use_relaxed=None):
        """
        Lazy extraction of select xtals

        Args:
            ids:
            N_atoms:
            overwrite:
            atttribute:
            use_relaxed (str): 'ff_relaxed' or 'vasp_relaxed'
        """
        (min_id, max_id) = ids
        if min_id is None: min_id = 1
        if max_id is None: max_id = self.get_max_id()

        (min_atoms, max_atoms) = N_atoms
        if min_atoms is None: min_atoms = 1
        if max_atoms is None: max_atoms = 5000

        ids, xtals = [], []
        for row in self.db.select():
            if overwrite or attribute is None or not hasattr(row, attribute):
                id, natoms = row.id, row.natoms
                if min_id <= id <= max_id and \
                    min_atoms < natoms <= max_atoms \
                    and id % self.size== self.rank:

                    xtal = self.get_pyxtal(id, use_relaxed)
                    yield id, xtal

    def update_row_energy(
        self,
        calculator='GULP',
        ids=(None, None),
        N_atoms=(None, None),
        ncpu=1,
        criteria=None,
        symmetrize=False,
        overwrite=False,
        write_freq=100,
        ff_lib='reaxff',
        steps=250,
        use_relaxed=None,
        cmd=None,
        calc_folder=None,
        skf_dir=None,
    ):
        """
        Update the row energy in the database for a given calculator.

        Args:
            calculator (str): 'GULP', 'MACE', 'VASP', 'DFTB'
            ids (tuple): A tuple specifying row IDs to update (e.g., (0, 100)).
            ncpu (int): number of parallel processes
            criteria (dict, optional): Criteria when selecting structures.
            symmetrize (bool): symmetrize the structure before calculation
            overwrite (bool): overwrite the existing energy attributes.
            write_freq (int): frequency to update db for ncpu=1
            ff_lib (str): Force field to use for GULP ('reaxff' by default).
            steps (int): Number of optimization steps for DFTB (default is 250).
            use_relaxed (str, optional): Use relaxed structures (e.g. 'ff_relaxed')
            cmd (str, optional): Command for VASP calculations
            calc_folder (str, optional): calc_folder for GULP/VASP calculations
            skf_dir (str, optional): Directory for DFTB potential files

        Functionality:
            Using the selected calculator, it updates the energy rows of the
            database. If `ncpu > 1`, run in parallel; otherwise in serial.

        Calculator Options:
            - 'GULP': Uses a force field (e.g., 'reaxff').
            - 'MACE': Uses the MACE calculator.
            - 'DFTB': Uses DFTB+ with symmetrization options.
            - 'VASP': Uses VASP, with a specified command (`cmd`).
        """

        label = calculator.lower() + "_energy"
        if calc_folder is None:
            calc_folder = calculator.lower() + "_calc"

        if calculator != 'MACE':
            #self.logging.info("make new folders", calc_folder, os.getpwd())
            os.makedirs(calc_folder, exist_ok=True)

        # Generate structures for calculation
        generator = self.select_xtal(ids, N_atoms, overwrite, label, use_relaxed)

        # Set up arguments for the chosen calculator
        args_up = []
        if calculator == 'GULP':
            args = [calculator, ff_lib, calc_folder, criteria]
            args_up = [ff_lib]
        elif calculator == 'MACE':
            args = [calculator, steps, criteria]
        elif calculator == 'DFTB':
            args = [calculator, skf_dir, steps, symmetrize, criteria]
        elif calculator == 'VASP':
            args = [calculator, calc_folder, cmd, criteria]
        else:
            raise ValueError(f"Unsupported calculator: {calculator}")

        # Perform calculation serially or in parallel
        self.logging.info(f"Rank-{self.rank} row_energy {calculator} {self.db_name}")
        if ncpu == 1:
            self.update_row_energy_serial(generator, write_freq, args, args_up)
        else:
            self.update_row_energy_mproc(ncpu, generator, args, args_up)
        self.logging.info(f"Rank-{self.rank} complete update_row_energy")

    def update_row_energy_serial(self, generator, write_freq, args, args_up):
        """
        Perform a serial update of row energies

        Args:
            generator (generator): Yielding tuples of (id, xtal), where:
                - `id` (int): Unique identifier for the structure.
                - `xtal` (object): pyxtal instance.
            write_freq (int): Frequency to update the database.
            args (list): Additional arguments to the function `opt_single`.
            args_up (list): Additional arguments for function `_update_db`.

        Functionality:
            It iterates over structures provided by `generator`,
            optimizes them using `opt_single`, and collects results that have
            converged (`status == True`). Once the number of results
            reaches `write_freq`, it updates the database.
        """
        results = []
        for id, xtal in generator:
            self.logging.info(f"Processing {id} {xtal.lattice} {args[0]}")
            print(f"Processing {id} {xtal.lattice} {args[0]}")
            res = opt_single(id, xtal, *args)
            (xtal, eng, status) = res
            if status:
                results.append((id, xtal, eng))
            if len(results) >= write_freq:
                self._update_db(results, args[0], *args_up)
                results = []
                self.print_memory_usage()
        if len(results) > 0:
            self._update_db(results, args[0], *args_up)

    def update_row_energy_mproc(self, ncpu, generator, args, args_up):
        """
        Perform parallel row energy updates by optimizing atomic structures.

        Args:
            ncpu (int): Number of CPUs to use for parallel processing.
            generator (generator): yielding tuples of (id, xtal), where:
                - `id` (int): Unique identifier for the structure.
                - `xtal` (object): pyxtal instance.
            args (list): Additional arguments passed to `call_opt_single`.
                - Typically includes a calculator or potential parameters.
            args_up (list): Additional arguments for function `_update_db`.

        Functionality:
            This function distributes the structures across multiple CPUs
            using `multiprocessing.Pool`. It creates chunks (based on `ncpu`),
            and process them in parallel by calling `call_opt_single`.
            Successful results are periodically written to the database.
            The function also prints memory usage after each database update.

        Parallelization Process:
            - The `Pool` is initialized with `ncpu` processes.
            - Structures are divided into chunks with the `chunkify` function.
            - Each chunk is processed by `call_opt_single` via the pool.
            - Successful results are periodically written to the database.
            - The pool is closed and joined after processing is complete.
        """
        from multiprocessing import Pool

        self.logging.info(f"Parallel optimizations {ncpu}")
        pool = Pool(processes=ncpu,
                    initializer=setup_worker_logger,
                    initargs=(self.log_file,))

        def chunkify(generator, chunk_size):
            chunk = []
            for item in generator:
                chunk.append(item)
                if len(chunk) == chunk_size:
                    yield chunk
                    chunk = []
            if chunk:
                yield chunk

        for chunk in chunkify(generator, ncpu*10):
            myargs = []
            for _id, xtal in chunk:
                if xtal is not None:
                    myargs.append(tuple([_id, xtal] + args))

            results = []
            self.logging.info(f"Start minicycle: {myargs[0][0]}-{myargs[-1][0]}")
            for result in pool.imap_unordered(call_opt_single,
                                              myargs,
                                              chunksize=1):
                if result is not None:
                    (myid, xtal, eng) = result
                    if eng is not None:
                        results.append(result)
                        numIons = sum(xtal.numIons)
                        count = len(results)
                        self.logging.info(f"Add {myid:4d} {eng:.3f} *{numIons} {count}")

                # Only do frequent update for slow calculator VASP
                if len(results) >= ncpu and args[0] == 'VASP':
                    self._update_db(results, args[0], *args_up)
                    self.logging.info(f"Finish minibatch: {len(results)}")
                    self.print_memory_usage()
                    results = []

            self.logging.info(f"Done  minicycle: {myargs[0][0]}-{myargs[-1][0]}")

            # After the loop, handle the remaining results
            if results:
                self.logging.info(f"Start  Update db: {len(results)}")
                self._update_db(results, args[0], *args_up)
                self.logging.info(f"Finish Update db: {len(results)}")

        pool.close()
        pool.join()

    def _update_db(self, results, calc, *args):
        """
        Update db with the calculation_results
        https://wiki.fysik.dtu.dk/ase/ase/db/db.html#writing-and-updating-many-rows-efficiently

        Args:
            results: list of (id, xtal, eng) tuples
            calc (str): calculator
        """
        #self.logging.info(f"====================Update db: {len(results)}")
        if calc == 'GULP':
            ff_lib = args[0]

        with self.db:
            for result in results:
                (id, xtal, eng) = result
                if xtal is not None:
                    if calc == 'GULP':
                        self.db.update(id,
                               ff_energy=eng,
                               ff_lib=ff_lib,
                               ff_relaxed=xtal.to_file())
                    elif calc == 'MACE':
                        self.db.update(id,
                               mace_energy=eng,
                               mace_relaxed=xtal.to_file())
                    elif calc == 'VASP':
                        self.db.update(id,
                                vasp_energy=eng,
                                vasp_relaxed=xtal.to_file())
                    elif calc == 'DFTB':
                        self.db.update(id,
                                dftb_energy=eng,
                                dftb_relaxed=xtal.to_file())
                    #self.logging.info(f'update_db_{calc}, {id}')

    def update_row_topology(self, StructureType="Auto", overwrite=True, prefix=None, ref_dim=3):
        """
        Update row topology base on the CrystalNets.jl

        Args:
            StructureType (str): 'Zeolite', 'MOF' or 'Auto'
            overwrite (bool): remove the existing attributes
            prefix (str): prefix for tmp cif file
            ref_dim (int): desired dimension
        """
        try:
            import juliacall
        except:
            raise RuntimeError(
                "Cannot load JuliaCall, Plz enable it before running")

        def parse_topology(topology_info):
            """
            Obtain the dimensionality and topology name
            """
            dim = 0
            name = ""
            detail = "None"
            for i, x in enumerate(topology_info):
                (d, n) = x
                if d > dim:
                    dim = d
                tmp = n.split(",")[0]
                if tmp.startswith("UNKNOWN"):
                    # tuple(int(num) for num in tmp[7:].split())
                    detail = tmp[7:]
                    tmp = "aaa"
                elif tmp.startswith("unstable"):
                    tmp = "unstable"
                name += tmp
                if i + 1 < len(topology_info):
                    name += "-"
            return dim, name, detail

        jl = juliacall.newmodule("MOF_Builder")
        jl.seval("using CrystalNets")
        jl.CrystalNets.toggle_warning(False)  # to disable warnings
        jl.CrystalNets.toggle_export(False)  # to disable exports
        if StructureType == "Zeolite":
            option = jl.CrystalNets.Options(structure=jl.StructureType.Zeolite)
        elif StructureType == "MOF":
            option = jl.CrystalNets.Options(structure=jl.StructureType.MOF)
        else:
            option = jl.CrystalNets.Options(structure=jl.StructureType.Auto)

        cif_file = prefix + ".cif" if prefix is not None else "tmp.cif"

        for row in self.db.select():
            if overwrite or not hasattr(row, "topology"):
                atoms = self.db.get_atoms(row.id)
                atoms.write(cif_file, format="cif", parallel=False)

                # Call crystalnet.jl
                result = jl.determine_topology(cif_file, option)
                # print(result)
                results = list(result) if len(result) > 1 else [result[0]]
                try:
                    topo = []
                    for res in results:
                        # topology for SingleNodes
                        name = str(res[0])
                        if res[1] > 1:
                            name += "(" + str(res[1]) + ")"
                        genome = res[0][jl.Clustering.Auto]
                        dim = jl.ndims(jl.CrystalNets.PeriodicGraph(genome))
                        topo.append((dim, name))

                    # The maximum dimensionality and topology name
                    dim, name, detail = parse_topology(topo)
                except:
                    dim, name, detail = 3, "error", "error"
                if dim == ref_dim:
                    print(
                        "Updating Topology",
                        row.space_group_number,
                        row.wps,
                        dim,
                        name,
                        detail[:10],
                    )
                # Unknown will be labeled as aaa
                self.db.update(row.id, topology=name,
                               dimension=dim, topology_detail=detail)
            # else:
            #    print("Existing Topology", row.topology)

    def update_db_description(self):
        """
        update db description based on robocrys
        Call robocrys: https://github.com/hackingmaterials/robocrystallographer
        """
        from robocrys import StructureCondenser, StructureDescriber

        condenser = StructureCondenser()
        describer = StructureDescriber()

        for row in self.db.select():
            if not hasattr(row, "description"):
                atoms = self.db.get_atoms(row.id)
                pmg = ase2pymatgen(atoms)
                try:
                    condensed_structure = condenser.condense_structure(pmg)
                    description = describer.describe(condensed_structure)
                except:
                    description = "N/A"

                self.db.update(row.id, description=description)
                print("\n======Updating\n", description)
            else:
                print("\n======Existing\n", row.description)

    def export_structures(
        self,
        fmt="vasp",
        folder="mof_out",
        criteria=None,
        sort_by="similarity",
        overwrite=True,
        cutoff=None,
        use_relaxed=None,
    ):
        """
        export structures from database according to the given criterion

        Args:
            fmt (str): 'vasp' or 'cif'
            folder (str): 'path of output folders'
            criteria (dict): check the validity with dict
            sort_by (str): sort by which attribute
            overwrite (bool): remove the existing folder
            cutoff (int): the maximum number of structures for export
        """

        import shutil

        if cutoff is None:
            cutoff = self.db.count()
        if not os.path.exists(folder):
            os.makedirs(folder)
        else:
            if overwrite:
                shutil.rmtree(folder)
                os.makedirs(folder)

        keys = [
            "id",
            "pearson_symbol",
            "space_group_number",
            "density",
            "dof",
            "similarity",
            "ff_energy",
            "vasp_energy",
            "mace_energy",
            "topology",
        ]
        properties = []
        for row in self.db.select():
            spg = row.space_group_number
            den = row.density
            dof = row.dof
            ps = row.pearson_symbol
            sim = float(row.similarity) if hasattr(
                row, "similarity") and row.similarity is not None else None
            top = row.topology if hasattr(row, "topology") else None
            ff_eng = float(row.ff_energy) if hasattr(
                row, "ff_energy") else None
            vasp_eng = float(row.vasp_energy) if hasattr(
                row, "vasp_energy") else None
            mace_eng = float(row.mace_energy) if hasattr(
                row, "mace_energy") else None
            properties.append([row.id, ps, spg, den, dof,
                              sim, ff_eng, vasp_eng, mace_eng,
                              top])

        dicts = {}
        for i, key in enumerate(keys):
            if properties[0][i] is not None:
                dicts[key] = [prop[i] for prop in properties]

        if sort_by in keys:
            col = keys.index(sort_by)  # + 1
        else:
            print("supported attributes", keys)
            raise ValueError("Cannot sort by", sort_by)

        print(f"====Exporting {len(properties)} structures")
        properties = [prop for prop in properties if prop[col] is not None]
        sorted_properties = sorted(properties, key=lambda x: x[col])

        for entry in sorted_properties[:cutoff]:
            [id, ps, spg, den, dof, sim, ff_eng, vasp_eng, mace_eng, top] = entry
            id = int(id)
            spg = int(spg)
            sim = float(sim)
            den = float(den)
            dof = int(dof)
            if vasp_eng is not None:
                eng = float(vasp_eng)
            elif mace_eng is not None:
                eng = float(mace_eng)
            elif ff_eng is not None:
                eng = float(ff_eng)
            else:
                eng = None
            if True:
            #try:
                xtal = self.get_pyxtal(id, use_relaxed)
                number, symbol = xtal.group.number, xtal.group.symbol.replace(
                    "/", "")
                # convert to the desired subgroup representation if needed
                #if number != spg:
                #    paths = xtal.group.path_to_subgroup(spg)
                #    xtal = xtal.to_subgroup(paths)
                #    number, symbol = (
                #        xtal.group.number,
                #        xtal.group.symbol.replace("/", ""),
                #    )

                label = os.path.join(
                    folder,
                    f"{id:d}-{xtal.get_Pearson_Symbol():s}-{number:d}-{symbol:s}",
                )

                status = xtal.check_validity(
                    criteria, True) if criteria is not None else True
            #except:
            #    status = False
            #    label = "Error"

            if status:
                try:
                    # if True:
                    xtal.set_site_coordination()
                    for s in xtal.atom_sites:
                        _l, _sp, _cn = s.wp.get_label(), s.specie, s.coordination
                        label += f"-{_l:s}-{_sp:s}{_cn:d}"
                    label += f"-S{sim:.3f}"
                    if len(label) > 40:
                        label = label[:40]
                except:
                    print("Problem in setting site coordination")

                if den is not None:
                    label += f"-D{abs(den):.2f}"
                if eng is not None:
                    label += f"-E{abs(eng):.3f}"
                if top is not None:
                    label += f"-T{top:s}"
                # if sim is not None: label += '-S{:.2f}'.format(sim)

                print("====Exporting:", label)
                if fmt == "vasp":
                    xtal.to_file(label + ".vasp", fmt="poscar")
                elif fmt == "cif":
                    xtal.to_file(label + ".cif")
            else:
                print("====Skippng:  ", label)

        return dicts

    def get_label(self, i):
        if i < 10:
            folder = f"cpu00{i}"
        elif i < 100:
            folder = f"cpu0{i}"
        else:
            folder = f"cpu0{i}"
        return folder

    def get_db_unique(self, db_name=None, prec=3):
        """
        Get a db file with only unique structures
        with the following identical attributes:
        (topology, ff_energy)

        Args:
            db_name (str): filename for the new db
            prec (int): ff_energy precision for the round number
        """

        print(f"The {self.db_name:s} has {self.db.count():d} strucs")
        if db_name is None:
            db_name = self.db_name[:-3] + "_unique.db"
        if os.path.exists(db_name):
            os.remove(db_name)

        unique_props = {}  # Using a dictionary to store unique properties
        for row in self.db.select():
            if hasattr(row, "topology") and hasattr(row, "ff_energy"):
                top, top_detail = row.topology, row.topology_detail
                dof, ff_energy = row.dof, round(row.ff_energy, prec)
                prop_key = (top, top_detail, ff_energy)
                # A dictionary lookup
                if prop_key in unique_props:
                    _id, _dof = unique_props[prop_key]
                    if dof < _dof:
                        print("Updating", row.id, top, ff_energy)
                        unique_props[prop_key] = (row.id, dof)
                    else:
                        print("Duplicate", row.id, top, ff_energy)
                else:
                    print("Adding", row.id, top, ff_energy)
                    unique_props[prop_key] = (row.id, dof)

        ids = [unique_props[key][0] for key in unique_props.keys()]
        with connect(db_name, serial=True) as db:
            for id in ids:
                row = self.db.get(id)
                kvp = {}
                for key in self.keys:
                    if hasattr(row, key):
                        kvp[key] = getattr(row, key)
                db.write(row.toatoms(), key_value_pairs=kvp)
        print(f"Created {db_name:s} with {db.count():d} strucs")
        return db.count()

    def check_overlap(self, reference_db, etol=2e-3, verbose=True):
        """
        Check the overlap w.r.t the reference database

        Args:
            reference_db (str): path of reference database
            etol (float): energy tolerence to distinguish the identical structure
            verbose (bool): whether or not print out details
        """

        db_ref = database_topology(reference_db, log_file = self.log_file)
        print(f"\nCurrent   database {self.db_name}: {self.db.count()}")
        print(f"Reference database {db_ref.db_name}: {db_ref.db.count()}")

        ref_data = []
        for row in db_ref.db.select():
            if hasattr(row, "topology") and hasattr(row, "ff_energy"):
                ref_data.append(
                    (row.topology, row.topology_detail, row.ff_energy))

        overlaps = []
        for row in self.db.select():
            if hasattr(row, "topology") and hasattr(row, "ff_energy"):
                for ref in ref_data:
                    (top, top_detail, ff_energy) = ref
                    if (
                        row.topology == top
                        and row.topology_detail == top_detail
                        and abs(row.ff_energy - ff_energy) < etol
                    ):
                        # strs = 'Find {:4d} {:6s}'.format(row.id, row.pearson_symbol)
                        # strs += ' {:12s} {:10.3f}'.format(row.topology, row.ff_energy)
                        # print(strs)
                        overlaps.append(
                            (
                                row.id,
                                row.pearson_symbol,
                                row.dof,
                                row.topology,
                                row.ff_energy,
                            )
                        )
                        break
        strs = f"\nThe number of overlap is: {len(overlaps):d}"
        strs += f"/{self.db.count():d}/{db_ref.db.count():d}"
        print(strs)
        sorted_overlaps = sorted(overlaps, key=lambda x: x[-1])
        if verbose:
            for entry in sorted_overlaps:
                print("{:4d} {:6s} {:4d} {:20s} {:10.3f}".format(*entry))

        return overlaps

    def print_info(self, excluded_ids=None, cutoff=100):
        """
        Print out the summary of the database based on the calculated energy
        Mostly used to quickly view the most interesting low-energy structures.
        Todo: show vasp_energy if available

        Args:
            excluded_ids (list): list of unwanted row ids
            cutoff (int): the cutoff value for the print
        """
        if excluded_ids is None:
            excluded_ids = []
        print(f"\nCurrent database {self.db_name}: {self.db.count()}")
        output = []
        for row in self.db.select():
            if row.id not in excluded_ids and hasattr(row, "topology") and hasattr(row, "ff_energy"):
                output.append(
                    (
                        row.id,
                        row.pearson_symbol,
                        row.dof,
                        row.topology,
                        row.ff_energy,
                    )
                )

        sorted_output = sorted(output, key=lambda x: x[-1])
        for entry in sorted_output[:cutoff]:
            print("{:4d} {:6s} {:4d} {:20s} {:10.3f}".format(*entry))

        strs = f"Showed structures: {len(sorted_output)}/{self.db.count()}"
        print(strs)

    def plot_histogram(self, prop, ax=None, filename=None, xlim=None, nbins=20):
        """
        Plot the histogram of a specified row property.

        Args:
            prop (str): The name of the property to plot (e.g., 'ff_energy').
            ax (matplotlib.axes.Axes, optional): Pre-existing axis to plot on.
                                                 If None, a new ax will be created.
            filename (str, optional): Path to save the plot (e.g., 'plot.png').
                                      If None, the plot will not be saved.
            xlim (tuple, optional): Limits for the x-axis (e.g., (0, 10)).
                                    If None, the x-axis will scale automatically.
            nbins (int, optional): Number of bins for the histogram. Default is 20.

        Returns:
            matplotlib.axes.Axes: The axis object with the histogram plotted.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            f, ax = plt.subplots()

        # Get the properties from the database
        props = self.get_properties(prop)

        # Check if there are values to plot
        if not props:
            raise ValueError(f"No rows contain the property '{prop}'.")

        ax.hist(props, nbins, density=True, alpha=0.75)

        # Set x-axis limits if provided
        if xlim is not None:
            ax.set_xlim(xlim)

        ax.set_xlabel(prop)

        # Save the plot if a filename is provided
        if filename is not None:
            plt.savefig(filename)

        return ax

    def get_properties(self, prop):
        """
        Retrieve a list of specific property values from the database rows.

        Args:
            prop (str): The property name to retrieve (e.g., 'ff_energy')

        Returns:
            list: A list of property values for rows that have the specified property.
                  If a row does not contain the property, it is ignored.

        Raises:
            Warning: If no rows in the database contain the specified property.
        """

        props = []

        # Loop through all rows in the database and collect the property values
        for row in self.db.select():
            if hasattr(row, prop):
                props.append(getattr(row, prop))

        # Print summary of rows
        name, count = self.db_name, self.db.count()
        print(f"Database {name} has {prop}: {len(props)}/{count}")

        # Warn if no properties were found
        if count == 0:
            raise Warning(
                f"No rows in the database contain the property '{prop}'.")

        return props


if __name__ == "__main__":
    # open
    if False:
        db = database("test.db")
        print("Total number of entries", len(db.codes))

        # view structure
        c = db.get_pyxtal("HXMTAM")
        print(c)
    if False:
        db = database_topology("../MOF-Builder/C-sp2/sp2-sacada-0506.db")
        # xtal = db.get_pyxtal(1)
        # print(xtal)
        # db.add_xtal(xtal, kvp={'similarity': 0.1})

        # db.update_row_ff_energy(ids=(0, 2), overwrite=True)
        # db.update_row_ff_energy(ncpu=2, ids=(2, 20), overwrite=True)
        # brew install coreutils to get timeout in maca
        # os.environ['ASE_DFTB_COMMAND'] = 'timeout 1m /Users/qzhu8/opt/dftb+/bin/dftb+ > PREFIX.out'
        os.environ["ASE_DFTB_COMMAND"] = "/Users/qzhu8/opt/dftb+/bin/dftb+ > PREFIX.out"
        skf_dir = "/Users/qzhu8/GitHub/MOF-Builder/3ob-3-1/"
        # db.update_row_dftb_energy(skf_dir, ncpu=1, ids=(0, 2), overwrite=True)
        db.update_row_dftb_energy(
            skf_dir, ncpu=1, ids=(17, 17), overwrite=True)

    db = database_topology("total.db")
    db.get_db_unique()
    db1 = database_topology("sp2_sacada.db")
    db1.get_db_unique()
    db = database_topology("total_unique.db")
    db.check_overlap("sp2_sacada_unique.db")
    db1.export_structures(folder="mof_out_sacada")
