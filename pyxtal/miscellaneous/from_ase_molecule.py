from pyxtal.molecule import *
from ase.build import molecule
from pymatgen.core import Molecule


def get_ase_mol(molname):
    """convert ase molecule to pymatgen style"""
    ase_mol = molecule(molname)
    pos = ase_mol.get_positions()
    symbols = ase_mol.get_chemical_symbols()
    return Molecule(symbols, pos)


if __name__ == "__main__":
    # ---------------------------------------------------
    for name in ["H2", "H2O", "HCl", "CS2", "C2Cl4", "PH3", "CH4", "C6H6", "C60"]:
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

        # orientation_in_wyckoff_position(mol, sg, WP's index in sg)
        # returns a list of orientations consistent with the WP's symmetry.
        # We can choose any of these orientations at random using np.random.choice
        # To use an orientation, do mol.apply_operation(orientation)
        # Spacegroup 16, index 6 has .2. symmetry

        # check 2 fold rotation
        allowed = orientation_in_wyckoff_position(mol, 16, 6, randomize=True)
        if allowed is not False:
            print(
                "Found " + str(len(allowed)) + " orientations for ",
                name,
                " with site symm 2",
            )
            for i, op in enumerate(allowed):
                mo = deepcopy(mol)
                mo.apply_operation(op)
                filename = "xyz/" + name + "-" + str(i) + ".xyz"
                mo.to(fmt="xyz", filename=filename)

        # check reflection
        allowed = orientation_in_wyckoff_position(mol, 25, 2, randomize=True)
        if allowed is not False:
            print(
                "Found " + str(len(allowed)) + " orientations for ",
                name,
                " with site symm m",
            )
            for i, op in enumerate(allowed):
                mo = deepcopy(mol)
                mo.apply_operation(op)
                filename = "xyz/" + name + "-" + str(i) + ".xyz"
                mo.to(fmt="xyz", filename=filename)

        # check 3 fold rotation
        allowed = orientation_in_wyckoff_position(mol, 147, 4, randomize=True)
        if allowed is not False:
            print(
                "Found " + str(len(allowed)) + " orientations for ",
                name,
                " with site symm 3",
            )
            for i, op in enumerate(allowed):
                mo = deepcopy(mol)
                mo.apply_operation(op)
                filename = "xyz/" + name + "-" + str(i) + ".xyz"
                mo.to(fmt="xyz", filename=filename)

        # check -1
        allowed = orientation_in_wyckoff_position(mol, 2, 2, randomize=True)
        if allowed is not False:
            print(
                "Found " + str(len(allowed)) + " orientations for ",
                name,
                " with site symm -1",
            )
            for i, op in enumerate(allowed):
                mo = deepcopy(mol)
                mo.apply_operation(op)
                filename = "xyz/" + name + "-" + str(i) + ".xyz"
                mo.to(fmt="xyz", filename=filename)

        # check 2/m
        allowed = orientation_in_wyckoff_position(mol, 64, 6, randomize=True)
        if allowed is not False:
            print(
                "Found " + str(len(allowed)) + " orientations for ",
                name,
                " with site symm 2/m",
            )
            for i, op in enumerate(allowed):
                mo = deepcopy(mol)
                mo.apply_operation(op)
                filename = "xyz/" + name + "-" + str(i) + ".xyz"
                mo.to(fmt="xyz", filename=filename)

        # check 6
        allowed = orientation_in_wyckoff_position(mol, 168, 3, randomize=True)
        if allowed is not False:
            print(
                "Found " + str(len(allowed)) + " orientations for ",
                name,
                " with site symm 6",
            )
            for i, op in enumerate(allowed):
                mo = deepcopy(mol)
                mo.apply_operation(op)
                filename = "xyz/" + name + "-" + str(i) + ".xyz"
                mo.to(fmt="xyz", filename=filename)

