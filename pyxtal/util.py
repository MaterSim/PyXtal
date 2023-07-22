"""
some utilities
"""

import numpy as np
from spglib import get_symmetry_dataset
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from pymatgen.core.structure import Structure
from ase import Atoms
from pyxtal.symmetry import Hall
import re

def find_dir(dirs):
    """
    a short function to find the correct dir from a list
    """
    skf_dir = None
    for d in dirs:
        if os.path.isdir(d):
            return d
    raise RuntimeError("Cannot find the dirtory for dftb parameters")

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
    symmetrize structure from pymatgen, and return the struc in conventional or
    primitive setting.

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
    check if the lattice has a good shape.

    Args:
        struc: pyxtal structure
    """

    para = struc.lattice.get_para(degree=True)
    if (max(para[:3])<maxvec) and (min(para[:3])>minvec)\
        and (max(para[3:])<maxang) and (min(para[3:])>minang):
        return True
    else:
        return False

def symmetrize(pmg, tol=1e-3, a_tol=5.0, style='pyxtal', hn=None):
    """
    symmetrize the structure from spglib.

    Args:
        pmg: pymatgen structure
        tol: tolerance
        a_tol: angle tolerance
        style: 'pyxtal' or spglib, differing in the choice of origin
        hn: hall_number

    Returns:
        pymatgen structure with symmetrized lattice
    """
    numbers = [site.species.elements[0].Z for site in pmg.sites]
    atoms = (pmg.lattice.matrix, pmg.frac_coords, numbers)
    dataset = get_symmetry_dataset(atoms, tol, angle_tolerance=a_tol)
    if hn is None:
        hn = Hall(dataset['number'], style=style).hall_default
    if hn != dataset['hall_number']:
        dataset = get_symmetry_dataset(atoms, tol,
                                       angle_tolerance=a_tol,
                                       hall_number=hn)
    cell = dataset['std_lattice']
    pos = dataset['std_positions']
    numbers = dataset['std_types']

    return Structure(cell, numbers, pos)

def get_symmetrized_pmg(pmg, tol=1e-3, a_tol=5.0, style='pyxtal', hn=None):
    """
    Get the symmetrized Pymatgen structure. A slight modification to ensure that
    the structure adopts the standard setting according to the Interational
    Crystallography Table.

    Args:
        pmg: input pymatgen structure
        tol: symmetry tolerance
        a_tol: angle tolerance
        style: 'pyxtal' or spglib, differing in the choice of origin
        hn: hall_number

    Returns:
        pymatgen structure with symmetrized lattice
    """

    pmg = symmetrize(pmg, tol, a_tol=a_tol, style=style, hn=hn)
    s = sga(pmg, symprec=tol, angle_tolerance=a_tol)
    # make sure that the coordinates are in standard setting
    if hn is None:
        hn = Hall(s._space_group_data['number'], style=style).hall_default
    if hn != s._space_group_data["hall_number"]:
        s._space_group_data = get_symmetry_dataset(s._cell, tol,
                                                   angle_tolerance=a_tol,
                                                   hall_number=hn)
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


def get_struc_from__parser(p):
    """
    A utility to get the pymatgen structure from the CifParser
    Sometimes the cif structure may have repeated atom entries

    Args:
        p: pymatgen CifParser object

    Return:
        a single pymatgen structure
    """
    from pymatgen.util.coord import find_in_coord_list_pbc
    from collections import OrderedDict
    from pymatgen.core.periodic_table import get_el_sp
    from pymatgen.io.cif import str2float
    import numpy as np

    def get_matching_coord(coord, ops, atol=1e-4):
        keys = list(coord_to_species.keys())
        coords = np.array(keys)
        for op in ops:
            c = op.operate(coord)
            inds = find_in_coord_list_pbc(coords, c, atol=atol)
            if len(inds):
                return keys[inds[0]]
        return False


    for i, d in enumerate(p._cif.data.values()):
        ops = p.get_symops(d)
        coord_to_species = OrderedDict()
        d0 = {"_atom_site_label": [], 
              "_atom_site_fract_x": [],
              "_atom_site_fract_y": [],
              "_atom_site_fract_z": [],
             }
        for i in range(len(d["_atom_site_label"])):
            try:
                symbol = p._parse_symbol(d["_atom_site_type_symbol"][i])
            except KeyError:
                symbol = p._parse_symbol(d["_atom_site_label"][i])
        
            el = get_el_sp(symbol)
            x = str2float(d["_atom_site_fract_x"][i])
            y = str2float(d["_atom_site_fract_y"][i])
            z = str2float(d["_atom_site_fract_z"][i])
        
            coord = (x, y, z)
            match = get_matching_coord(coord, ops)
            if not match:
                d0['_atom_site_label'].append(el)
                d0["_atom_site_fract_x"].append(str(x))
                d0["_atom_site_fract_y"].append(str(y))
                d0["_atom_site_fract_z"].append(str(z))
                coord_to_species[coord] = el
        d.data['_atom_site_label'] = d0['_atom_site_label']
        d.data['_atom_site_fract_x'] = d0['_atom_site_fract_x']
        d.data['_atom_site_fract_y'] = d0['_atom_site_fract_y']
        d.data['_atom_site_fract_z'] = d0['_atom_site_fract_z']
    
        s = p._get_structure(d, primitive=False, symmetrized=False)
        return s

def Kgrid(struc, Kresol=0.10, dimension=3):
    """
    Assign kpoints based on the lattice
    """
    a, b, c, alpha, beta, gamma = struc.get_cell_lengths_and_angles()
    vol = struc.get_volume()
    dist    = np.zeros(3);
    dist[2] = np.abs(vol/(a*b*np.sin(np.radians(gamma))))
    dist[1] = np.abs(vol/(a*c*np.sin(np.radians(beta))))
    dist[0] = np.abs(vol/(b*c*np.sin(np.radians(alpha))))
    Kpoints = np.ceil(1./(dist*Kresol))
    if dimension == 2:
        Kpoints[-1] = 1;
    #print(a, b, c, alpha, beta, gamma)
    #print(vol/(a*b*np.sin(gamma)), a*b, np.sin(gamma))
    #print(vol/(a*c*np.sin(beta)), a*c, np.sin(beta))
    #print(vol/(b*c*np.sin(alpha)), b*c, np.sin(alpha))
    #print(Kpoints)
    #import sys; sys.exit()
    return Kpoints.astype(int)


def sort_by_dimer(atoms, N_mols, id=10, tol=4.0):
    """
    sort the ase atoms' xyz according to dimer
    so far only tested on aspirin

    Args:
        atoms: atoms object from pyxtal
        N_mols: number of molecules
        id: the refrence atom id
        tol: tolerence distance to check if it is a dimer
    """

    N_atoms = int(len(atoms)/N_mols)
    pos = atoms.get_scaled_positions()
    refs = pos[id:len(pos):N_atoms, :]
    #print(refs)

    # compuate the indices and shift
    orders = []
    shifts = []
    while len(orders) < N_mols:
        lefts = [i for i in range(N_mols) if i not in orders]
        i = lefts[0]
        orders.append(i)
        shifts.append(np.zeros(3))
        ref_i = refs[i]
        good = False
        for j in lefts[1:]:
            ref_j = refs[j]
            dist = ref_j - ref_i
            shift = np.round(dist)
            dist -= shift
            dist = np.linalg.norm(dist.dot(atoms.cell[:]))
            if dist < tol:
                orders.append(j)
                shifts.append(shift)
                good = True
                break
        if not good:
            raise RuntimeError('Cannot find match on molecule', i)
        else:
            print('get', i, j, dist, shift)

    pos0 = atoms.get_positions()
    pos1 = np.zeros([len(pos), 3])
    for i, id in enumerate(orders):
        s1, e1 = id*N_atoms, (id+1)*N_atoms
        s2, e2 = i*N_atoms, (i+1)*N_atoms
        pos1[s2:e2, :] += pos0[s1:e1, :] - shifts[i].dot(atoms.cell[:])

    atoms.set_positions(pos1)
    return atoms

def generate_wp_lib(spg_list, composition, 
                    num_wp=(None, None), 
                    num_fu=(None, None), 
                    num_dof=(None, None), 
                    N_max=1000):
    """
    Generate wps according to the composition constraint (e.g., SiO2)

    Args;
        - spg_list: list of space group choices
        - composition: chemical compositions [1, 2]
        - num_wp: (min_wp, max_wp)
        - num_fu: (min_fu, max_fu)
        - num_dof: (min_dof, max_dof)

    Returns:
        a list of wps [spg, ([wp1, ...], ... [wp1, ...]), dof]
    """

    from pyxtal.symmetry import Group

    composition = np.array(composition, dtype=int)
    (min_wp, max_wp) = num_wp
    (min_fu, max_fu) = num_fu
    (min_dof, max_dof) = num_dof

    if max_wp is None: max_wp = len(composition)
    if min_wp is None: min_wp = len(composition)
    if min_dof is None: min_dof = 1
    if max_dof is None: max_dof = 1000
    #print(max_wp, min_wp)
    wps = []
    for sg in spg_list:
        g = Group(sg)
        lat_dof = g.get_lattice_dof()
        # determine the upper and lower limit 
        if min_fu is None: min_fu = max([int(len(g[-1])/min(composition)), 1])
        if max_fu is None: max_fu = max([int(len(g[0])/max(composition)), 1])
        count = 0
        for i in range(max_fu, min_fu-1, -1):
            letters, _, wp_ids = g.list_wyckoff_combinations(
                    composition*i, max_wp=max_wp, 
                    min_wp=min_wp, Nmax=100000)
            for j, wp in enumerate(wp_ids):
                wp_dofs = 0
                num = 0
                for wp0 in wp:
                    for id in wp0:
                        wp_dofs += g[id].get_dof()
                        num += g[id].multiplicity
                #print(sg, wp, letters[j])
                num_dof = lat_dof + wp_dofs
                if min_dof <= num_dof <= max_dof:
                    wps.append((num, sg, wp, lat_dof + wp_dofs))
                    count += 1
            if count >= N_max:
                break
    return wps 
 

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
