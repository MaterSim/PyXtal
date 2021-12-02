"""
some utilities
"""

from spglib import get_symmetry_dataset
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from pymatgen.core.structure import Structure
from ase import Atoms
from pyxtal.symmetry import Group
import re

def listToString(s): 
    # initialize an empty string
    #str1 = " " 
    str1 = "" 
    # return string  
    return (str1.join(s))


def pymatgen2ase(struc):
    """
    A short cut to convert between pymatgen to ase
    """
    atoms = Atoms(symbols = struc.atomic_numbers, cell = struc.lattice.matrix, pbc=True)
    atoms.set_scaled_positions(struc.frac_coords)
    return atoms

def ase2pymatgen(struc):
    """
    A short cut to convert between pymatgen to ase
    """
    lattice = struc.cell
    coordinates = struc.get_scaled_positions()
    species = struc.get_chemical_symbols()
    return Structure(lattice, species, coordinates)

def symmetrize_cell(struc, mode='C'):
    """
    symmetrize structure from pymatgen, and return the struc in conventional/primitive setting

    Args:
        struc: ase type
        mode: output conventional or primitive cell
    """
    P_struc = ase2pymatgen(struc)
    finder = sga(P_struc,symprec=0.06)
    if mode == 'C':
        P_struc = finder.get_conventional_standard_structure()
    else:
        P_struc = finder.get_primitive_standard_structure()

    return pymatgen2ase(P_struc)

def good_lattice(struc, maxvec=25.0, minvec=1.2, maxang=150, minang=30):
    """
    check if the lattice has a good shape

    Args:
        struc: pyxtal structure
    """

    para = struc.lattice.get_para(degree=True)
    if (max(para[:3])<maxvec) and (min(para[:3])>minvec)\
        and (max(para[3:])<maxang) and (min(para[3:])>minang):
        return True
    else:
        return False

def symmetrize(pmg, tol=1e-3, a_tol=5.0):
    """
    symmetrize the structure from spglib

    Args:
        pmg: pymatgen structure
        tol: tolerance

    Returns:
        pymatgen structure with symmetrized lattice
    """
    numbers = [site.species.elements[0].Z for site in pmg.sites]
    atoms = (pmg.lattice.matrix, pmg.frac_coords, numbers)
    dataset = get_symmetry_dataset(atoms, tol, angle_tolerance=a_tol)
    hn = Group(dataset['number'], 3).hall_number
    if hn != dataset['hall_number']:
        dataset = get_symmetry_dataset(atoms, tol, angle_tolerance=a_tol, hall_number=hn)
    cell = dataset['std_lattice']
    pos = dataset['std_positions']
    numbers = dataset['std_types']

    return Structure(cell, numbers, pos)

def get_symmetrized_pmg(pmg, tol=1e-3, a_tol=5.0):
    """
    Symmetrized Pymatgen structure
    A slight modification to ensure that the structure adopts the
    standard setting used in interational crystallography table

    Args:
        pmg: input pymatgen structure
        tol: symmetry tolerance
    """

    pmg = symmetrize(pmg, tol, a_tol=a_tol)
    s = sga(pmg, symprec=tol, angle_tolerance=a_tol)
    hn = Group(s.get_space_group_number()).hall_number
    # make sure that the coordinates are in standard setting
    if hn != s._space_group_data["hall_number"]:
        s._space_group_data = get_symmetry_dataset(s._cell, tol, angle_tolerance=a_tol, hall_number=hn)
    return s.get_symmetrized_structure(), s.get_space_group_number()

def extract_ase_db(db_file, id):
    """
    a short cut to extract the structural information
    from the ase db file by row id
    """

    from ase.db import connect
    from pyxtal import pyxtal
    import os

    if not os.path.exists("output"):
        os.makedirs("output")
    print("Dumping the structures to the folder output")
    with connect(db_file) as db:
        for id in ids:
            s = db.get_atoms(id=id)
            filename = "output/"+str(id)+".vasp"
            s.write(filename, format='vasp', direct=True, vasp5=True)
            my = pyxtal()
            my.from_seed(s)
            print(my)

def parse_cif(filename, header=False, spg=False, eng=False, csd=False, sim=False):
    """
    read structures from a cif (our own format with #END)
    Args:
        filename: string
        header: bool, whether or not return header
        spg: bool, whether or not return the spg
    """
    strings = []
    headers = []
    spgs = []
    engs = []
    csds = []
    sims = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        start = None
        end = None
        for i in range(len(lines)):
            if lines[i].find("data_") == 0:
                if sim:
                    sims.append(float(lines[i].split(':')[-1]))
                end = i
                if start is not None:
                    tmp = []
                    for l in lines[start:end-1]:
                        if len(re.findall(r"[0-9][B-C]", l))>0 or \
                        len(re.findall(r"[A-Z][0-9]\' [0-9]", l))>0:
                            #print(l) #; import sys; sys.exit()
                            continue
                        tmp.append(l)
                    cif = listToString(tmp)
                    strings.append(cif)
                start = i
                headers.append(lines[i])
            elif lines[i].find("_symmetry_Int_Tables_number") == 0:
                spgs.append(int(lines[i].split()[-1]))
            elif lines[i].find("#Energy") == 0:
                engs.append(float(lines[i].split()[1]))
            elif lines[i].find("_database_code") == 0:
                tmp = lines[i].split()[-1]
                csds.append(tmp[:-1])

        #Last one
        tmp = []
        for l in lines[start:]:
            if len(re.findall(r"[0-9][B-D]", l))>0 or \
            len(re.findall(r"[A-Z][0-9]\' [0-9]", l))>0:
                #print(l); 
                continue
            tmp.append(l)
        cif = listToString(tmp)
        strings.append(cif)

    if header:
        return strings, headers
    elif spg:
        return strings, spgs
    elif eng:
        return strings, engs
    elif csd:
        return strings, csds
    elif sim:
        return strings, sims
    else:
        return strings

def process_csd_cif(cif, remove_H=False):
    """
    process cif from CSD, sometimes it contains multiple 
    e.g., C2
    """
    lines = cif.split('\n')
    tmp = []
    for l in lines:
        if len(re.findall(r"[0-9][A-Z]", l))>0 or len(re.findall(r"[A-Z]\?", l))>0 or \
            len(re.findall(r"[0-9]\?", l))>0 or \
            len(re.findall(r"[A-Z][0-9]\' [0-9]", l))>0:
            #len(re.findall(r"[0-9]_[0-9]", l))>0 or \
            #print(l) #; import sys; sys.exit()
            continue
        else:
            if remove_H and len(re.findall(r" H ", l))>0:
                continue
            else:
                tmp.append(l+'\n')
    return listToString(tmp)

def get_similar_cids_from_pubchem(base, MaxRecords):

    """
    Args:
        base: PubChem CID of Starting chemical
        MaxRecords: Number of Similar Compounds

    Returns:
        List of the CIDs of PubChem compounds similar to the base compound.
    """

    import pubchempy as pcp

    if type(base) == int: base = str(base)
    cids = pcp.get_compounds(base, 
                             searchtype="similarity", 
                             MaxRecords=MaxRecords)
    results = []
    for x in cids:
        csd_codes = search_ccdc_structures(x.cid)
        if len(csd_codes)>0:
            d = {"cid": x.cid,
                 "smiles": x.canonical_smiles,
                 "name": x.iupac_name,
                 "csd_codes": csd_codes}
            results.append(d)
            print(d)
    return results

def search_csd_code_by_pubchem(cid):
    """
    Args:
        cid: PubChem cid

    Returns:
        CIDs that have CCDC crystal structure data
    """

    import urllib
    import json
    from monty.json import MontyDecoder

    url0 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'
    cid=str(cid)
    url = url0 + cid + '/JSON'
    csd_codes = []

    try:
        response = urllib.request.urlopen(url)
    except urllib.error.HTTPError:
        print("Problem in http connection", url)
        return None
    except urllib.error.URLError:
        print("Problem in parsing", url)
        return None

    try:
        contents = response.read()
        if contents.decode().find('CCDC') > 0:
            data = json.loads(contents, cls=MontyDecoder)
            if 'Section' in data['Record']['Section'][0].keys():
                if len(data['Record']['Section'][0]['Section']) == 3:
                    infos = data['Record']['Section'][0]['Section'][2]['Section'][0]['Information']
                    for info in infos:
                        csd_codes.append(info['Value']['StringWithMarkup'][0]['String'])
    except:
        print('Failed to parse json', url, '\n')

    return csd_codes

def search_csd_entries_by_code(code):
    """
    Args:
        code: CSD code, e.g., ACSALA

    Returns:
        list of csd ids
    """

    from ccdc.search import TextNumericSearch
    from ccdc.crystal import PackingSimilarity as PS

    def new_cryst(cryst, crysts, n_max):
        spg1 = cryst.spacegroup_number_and_setting[0]
        for ref in crysts:
            spg2 = ref.spacegroup_number_and_setting[0]
            if spg1 == spg2:
                h = PS().compare(cryst, ref)
                #print(cryst.identifier, ref.identifier, h.nmatched_molecules)
                if h is not None and h.nmatched_molecules == n_max:
                    return False
        return True

    n_max = PS().settings.packing_shell_size
    query = TextNumericSearch()
    query.add_identifier(code)
    hits = query.search()
    unique_crysts = []
    for hit in hits:
        if hit.entry.has_3d_structure and hit.entry.pressure is None:
            if new_cryst(hit.crystal, unique_crysts, n_max):
                unique_crysts.append(hit.crystal)
                #print(hit.entry.identifier, hit.entry.deposition_date)

    return [c.identifier for c in unique_crysts]


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "-f",
        dest="file",
        help="path of database file"
    )
    parser.add_argument(
        "-i",
        dest="id",
        help="index of the row",
    )

    options = parser.parse_args()
    ids = options.id
    if ids.find(",")>0:
        ids = [int(id) for id in ids.split(",")]
    else:
        ids = [int(ids)]
    extract_ase_db(options.file, ids)
