"""
Database class
"""
import os
from ase.db import connect
from pyxtal import pyxtal
import pymatgen.analysis.structure_matcher as sm

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
                     'similarity',
                     'ff_energy',
                     'density',
                     'dof',
                     'topology',
                     'dimension',
                     'similarity0',
                     'wps'
                    ]
        self.matcher = sm.StructureMatcher(ltol=ltol, stol=stol, angle_tol=atol)

    def vacuum(self):
        self.db.vacuum()

    def get_pyxtal(self, id):
        from pyxtal import pyxtal
        from pyxtal.util import ase2pymatgen

        atom = self.db.get_atoms(id=id)
        pmg = ase2pymatgen(atom)
        xtal = pyxtal()
        try:
            xtal.from_seed(pmg)
            if xtal is not None and xtal.valid:
                for key in self.keys:
                    if hasattr(row, key):
                        setattr(xtal, key, getattr(row, key))
                xtals.append(xtal)
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

        return xtals

    def add_xtal(self, xtal, kvp):
        """
        Add new xtal to the given db
        """
        spg_num = xtal.group.number
        atoms = xtal.to_ase(resort=False)
        density = xtal.get_density()
        dof = xtal.get_dof()
        wps = [s.wp.get_label() for s in xtal.atom_sites]
        _kvp = {"space_group_number": spg_num,
                "wps": str(wps),
                "density": density,
                "dof": dof,
               }
        kvp.update(_kvp)
        self.db.write(atoms, key_value_pairs=kvp)

    def check_new_structure(self, xtal, same_group=True):
        """
        Check if the input xtal already exists in the db

        Args:
            xtal: pyxtal object
            same_group (bool): keep the same group or not
        """
        from pyxtal.util import ase2pymatgen

        s_pmg = xtal.to_pymatgen()
        for row in self.db.select():
            ref = self.db.get_atoms(id=row.id)
            ref_pmg = None
            if same_group:
                if row.space_group_number != xtal.group.number:
                    continue
            ref_pmg = ase2pymatgen(ref)
            if self.matcher.fit(s_pmg, ref_pmg, symmetric=True):
                return False
        return True

    def clean_structures(self, dtol=1e-3, etol=1e-3):
        """
        Clean up the db by removing the duplicate structures
        Here we check the follow criteria
            - same number of atoms
            - same density
            - same energy

        Args:
            dtol (float): tolerance of density
            etol (float): tolerance of energy
        """

        unique_rows = []
        to_delete = []

        for row in self.db.select():
            unique = True
            for prop in unique_rows:
                (natoms, den, ff_energy) = prop
                if natoms==row.natoms and abs(den-row.density) < dtol:
                    if hasattr(row, 'ff_energy'):
                        if abs(row.ff_energy-ff_energy) < etol:
                            unique = False
                            break
                    else:
                        unique = False
                        break
            if unique:
                unique_rows.append((row.natoms, row.density, row.ff_energy))
            else:
                to_delete.append(row.id)
        print("The following structures were deleted", to_delete)
        self.db.delete(to_delete)

    def update_row_topology(self):
        """
        Update row topology base on the CrystalNets.jl
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
            for i, x in enumerate(topology_info):
                (d, n) = x
                if d > dim: dim = d
                tmp = n.split(',')[0]
                if tmp.startswith("UNKNOWN"):
                    tmp = 'aaa'
                elif tmp.startswith("unstable"):
                    tmp = 'unstable'
                name += tmp
                if i + 1 < len(topology_info): name += '-'
            return dim, name


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
                    dim, name = parse_topology(topo)
                except:
                    dim, name = 3, 'error'
                print("Updating", row.spg, row.wps, dim, name)
                # Unknown will be labeled as aaa
                self.db.update(row.id, topology=name, dimension=dim)
            else:
                print("Existing", row.topology)

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

    def export_structures(self, fmt='vasp', folder='mof_out', check=False,
                          sort_by='similairty'):
        """
        export structures from database according to the given criterion

        Args:
            fmt (str): 'vasp' or 'cif'
            folder (str): 'path of output folders'
            check (bool): whether or not check the validity
            sort_by (str): sort by which attribute
        """

        if not os.path.exists(folder): os.makedirs(folder)
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
            atoms.append(db.get_atoms(id=row.id))

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
                xtal = pyxtal()
                try:
                    xtal.from_seed(atoms[mid])
                    number, symbol = xtal.group.number, xtal.group.symbol.replace('/','')
                    # convert to the desired subgroup representation if needed
                    if number != spg:
                        paths = xtal.group.path_to_subgroup(spg)
                        xtal = xtal.to_subgroup(paths)
                        number, symbol = xtal.group.number, xtal.group.symbol.replace('/','')

                    label = os.path.join(folder, '{:d}-{:d}-{:s}'.format(id, number, symbol))

                    if check:
                        status = xtal.check_validity(self.criteria)
                    else:
                        status = True
                except:
                    status = False

                if status:
                    if top is not None: print(top)
                    try:
                        xtal.set_site_coordination()
                        for s in xtal.atom_sites:
                            _l, _sp, _cn = s.wp.get_label(), s.specie, s.coordination
                            label += '-{:s}-{:s}{:d}'.format(_l, _sp, _cn)
                        label += '-S{:.3f}'.format(sim)
                    except:
                        print('Problem in setting site coordination')
                    if len(label)>40: label = label[:40]

                    if den is not None: label += '-D{:.2f}'.format(abs(den))
                    if eng is not None: label += '-E{:.3f}'.format(abs(eng))
                    if top is not None: label += '-{:s}'.format(top)

                    print("====Exporting:", label)
                    if fmt == 'vasp':
                        atoms[mid].write(label+'.vasp', format='vasp', vasp5=True, direct=True)
                    elif fmt == 'cif':
                        xtal.to_file(label+'.cif')
                else:
                    print("====Skippng:  ", label)


if __name__ == "__main__":
    # open
    db = database('test.db')
    print("Total number of entries", len(db.codes))

    # view structure
    c = db.get_pyxtal('HXMTAM')
    print(c)
