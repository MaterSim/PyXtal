from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase import Atoms

"""
scripts to perform structure conformation
"""
def pymatgen2ase(struc):
    atoms = Atoms(symbols = struc.atomic_numbers, cell = struc.lattice.matrix)
    atoms.set_scaled_positions(struc.frac_coords)
    return atoms

def ase2pymatgen(struc):
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
    finder = SpacegroupAnalyzer(P_struc,symprec=0.06,angle_tolerance=5)
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

def symmetrize(pmg, tol=1e-3):
    """
    symmetrize the structure from spglib

    Args:
        pmg: pymatgen structure
        tol: tolerance

    Returns:
        pymatgen structure with symmetrized lattice
    """
    from spglib import get_symmetry_dataset

    atoms = (pmg.lattice.matrix, pmg.frac_coords, pmg.atomic_numbers)
    dataset = get_symmetry_dataset(atoms, tol)
    cell = dataset['std_lattice']
    pos = dataset['std_positions']
    numbers = dataset['std_types']
    pmg = Structure(cell, numbers, pos)
    return pmg


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
    
