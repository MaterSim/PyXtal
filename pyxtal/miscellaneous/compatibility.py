from pyxtal.molecule import *
from pyxtal.structure import *
from ase.build import molecule
from pymatgen.core import Molecule

letters = "abcdefghijklmnopqrstuvwxyzA"


def get_ase_mol(molname):
    """convert ase molecule to pymatgen style"""
    ase_mol = molecule(molname)
    pos = ase_mol.get_positions()
    symbols = ase_mol.get_chemical_symbols()
    return Molecule(symbols, pos)


if __name__ == "__main__":
    # ---------------------------------------------------
    for name in ["C60"]:  # ['H2O', 'CS2' ,'CH4']:
        mol = get_ase_mol(name)
        pga = PointGroupAnalyzer(mol)

        # Symmetrize the molecule using pymatgen
        mol = pga.symmetrize_molecule()["sym_mol"]
        pga = PointGroupAnalyzer(mol)

        print(name, " has point group symmetry: ", pga.get_pointgroup())

        # Check if orders of rotation are detected correctly
        pg = pga.get_pointgroup()
        for op in pg:
            opa = OperationAnalyzer(op)
            if opa.order == "irrational":
                print(opa)
            elif opa.order > 10:
                print(opa)

        """symm = get_symmetry(mol)
        opas = []
        for op in symm:
            opa = OperationAnalyzer(op)
            if opa.type == "rotation" and opa.order == 3:
                for op2 in symm:
                    opa2 = OperationAnalyzer(op2)
                    if opa2.type == "rotation" and opa2.order == 2:
                        if abs(np.dot(opa.axis, opa2.axis)) < .02:
                            print("=======")
                            print(np.dot(opa.axis,opa2.axis))
                            print(angle(opa.axis,opa2.axis))
                break"""

        for sg in range(142, 231):
            symmetry = get_wyckoff_symmetry(sg, molecular=True)
            for index in range(1, len(symmetry)):
                letter = letters[len(symmetry) - 1 - index]
                ops = symmetry[index][0]
                allowed = orientation_in_wyckoff_position(
                    mol, sg, index, randomize=True
                )
                if allowed is False:
                    print(
                        name
                        + ": found "
                        + "0"
                        + " orientations in "
                        + letter
                        + " site symm: "
                        + ss_string_from_ops(ops, sg)
                        + " space group: "
                        + str(sg)
                    )
                    # for i, op in enumerate(allowed):
                    #    mo = deepcopy(mol)
                    #    mo.apply_operation(op)
                    #    print(mo)
                    #    filename = 'xyz/' + name + '-' + str(i)+'.xyz'
                    #    mo.to(fmt='xyz',filename=filename)
