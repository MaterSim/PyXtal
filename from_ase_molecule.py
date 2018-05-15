from molecule import *
from ase.build import molecule
from pymatgen import Molecule

def get_ase_mol(molname):
    """convert ase molecule to pymatgen style"""
    ase_mol = molecule(molname)
    pos = ase_mol.get_positions()
    symbols = ase_mol.get_chemical_symbols()
    return(Molecule(symbols, pos))

if __name__ == "__main__":
#---------------------------------------------------
    for name in ['H2', 'H2O', 'HCl', 'PH3', 'CH4', 'C60']:
        mol = get_ase_mol(name)
        pga = PointGroupAnalyzer(mol)
        print(name, ' has point group symmetry: ', pga.get_pointgroup())

        #orientation_in_wyckoff_position(mol, sg, WP's index in sg)
        #returns a list of orientations consistent with the WP's symmetry.
        #We can choose any of these orientations at random using np.random.choice
        #To use an orientation, do mol.apply_operation(orientation)
        #Spacegroup 16, index 6 has .2. symmetry
        allowed =  orientation_in_wyckoff_position(mol, 16, 6, randomize=True)
        if allowed is not False:
            print("Found "+str(len(allowed))+" orientations for ", name)
            for i, op in enumerate(allowed):
                mo = deepcopy(mol)
                mo.apply_operation(op)
                filename = 'xyz/' + name + '-' + str(i)+'.xyz'
                mo.to(fmt='xyz',filename=filename)
