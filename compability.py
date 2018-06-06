from molecule import *
from structure import *
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
    for name in ['H2O', 'CS2','CH4']:
        mol = get_ase_mol(name)
        pga = PointGroupAnalyzer(mol)

        #Symmetrize the molecule using pymatgen
        mol = pga.symmetrize_molecule()['sym_mol']
        pga = PointGroupAnalyzer(mol)

        print(name, ' has point group symmetry: ', pga.get_pointgroup())

        #Check if orders of rotation are detected correctly
        pg = pga.get_pointgroup()
        for op in pg:
            opa = OperationAnalyzer(op)
            if opa.order == 'irrational':
                print(opa)
            elif opa.order > 10:
                print(opa)

        for sg in range(1,230):
            symmetry = get_wyckoff_symmetry(sg, molecular=True)
            for index in range(1, len(symmetry)):
                ops=symmetry[index][0]
                allowed =  orientation_in_wyckoff_position(mol, sg, index, randomize=True)
                if allowed is not False:
                    print(name + ": found "+str(len(allowed))+ " orientations in " +
                            ' site symm: ' + ss_string_from_ops(ops, sg) + 
                            ' space group: ' + str(sg))
                    #for i, op in enumerate(allowed):
                    #    mo = deepcopy(mol)
                    #    mo.apply_operation(op)
                    #    filename = 'xyz/' + name + '-' + str(i)+'.xyz'
                    #    mo.to(fmt='xyz',filename=filename)
