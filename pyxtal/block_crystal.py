"""
Module for generating crystals based on building blocks
"""

# Standard Libraries
import numpy as np

# External Libraries
from pymatgen.core import Molecule

from pyxtal.io import search_molecules_in_crystal

# PyXtal imports
from pyxtal.molecular_crystal import molecular_crystal as mol_xtal
from pyxtal.molecule import Orientation, compare_mol_connectivity, pyxtal_molecule
from pyxtal.wyckoff_site import mol_site


def block_crystal(
    dim,
    group,
    molecules,
    numMols,
    factor,
    thickness,
    area,
    block,
    num_block,
    lattice,
    torsions,
    sites,
    conventional,
    tm,
    seed,
    random_state,
    use_hall,
):
    # If block is None, directly generate mol. xtal.
    # Otherwise, generate crystal from building block
    if block is None:
        return mol_xtal(
            dim,
            group,
            molecules,
            numMols,
            factor,
            thickness=thickness,
            area=area,
            lattice=lattice,
            torsions=torsions,
            sites=sites,
            conventional=conventional,
            tm=tm,
            seed=seed,
            random_state=random_state,
            use_hall=use_hall,
        )

    else:
        p_mol = pyxtal_molecule(block)
        block_mols = search_molecules_in_crystal(p_mol.mol, tol=0.2, once=False)
        xtal_mols = [pyxtal_molecule(m, fix=True) for m in molecules]

        orders = []
        for m1 in xtal_mols:
            for m2 in block_mols:
                if len(m1.mol) == len(m2):
                    # match, mapping = compare_mol_connectivity(m2, m1.mol)
                    match, mapping = compare_mol_connectivity(m1.mol, m2)
                    # print(match, len(m1.mol), len(m2))
                    if match:
                        orders.append([mapping[at] for at in range(len(m2))])
                        break
                        # print(mapping)
                        # print("smile"); print(m1.mol.to('xyz'))
                        # print("reference"); print(m2.to('xyz'))
                        # print(orders[-1])
                        # import sys; sys.exit()

        if len(orders) != len(molecules):
            raise ValueError("Block is inconsistent with the molecules")

        # rearrange the order of block molecules
        numbers = []
        coords = np.zeros([len(p_mol.mol), 3])
        count = 0
        for order, m in zip(orders, block_mols):
            numbers.extend([m.atomic_numbers[o] for o in order])
            coords[count : count + len(m)] += m.cart_coords[order]
            count += len(m)
        mol = Molecule(numbers, coords)
        if num_block is not None:
            num_block = [num_block]
        # print(mol.to('xyz')); import sys; sys.exit()

        for i in range(10):
            struc = mol_xtal(
                dim,
                group,
                [mol],
                num_block,
                factor,
                thickness=thickness,
                area=area,
                lattice=lattice,
                torsions=torsions,
                sites=sites,
                conventional=conventional,
                tm=tm,
                random_state=random_state,
                use_hall=use_hall,
            )

            if struc.valid:
                break

        struc.molecules = xtal_mols
        struc.numMols = struc.numMols * len(xtal_mols)

        # XYZ from the block_sites
        b_site = struc.mol_sites[0]
        xyz, _ = b_site._get_coords_and_species(absolute=True, first=True)

        # Create mol sites
        ori = Orientation(np.eye(3))
        mol_sites = []

        count = 0
        for i in range(len(molecules)):
            mol = xtal_mols[i]
            m_xyz = xyz[count : count + len(mol.mol)]
            center = mol.get_center(m_xyz)
            mol.reset_positions(m_xyz - center)
            # print(mol.mol.to('xyz'))
            position = np.dot(center, np.linalg.inv(struc.lattice.matrix))
            position -= np.floor(position)
            m_site = mol_site(mol, position, ori, b_site.wp, struc.lattice)
            m_site.type = i
            mol_sites.append(m_site)
            count += len(mol.mol)

        struc.mol_sites = mol_sites

        return struc


if __name__ == "__main__":
    import pymatgen.analysis.structure_matcher as sm

    from pyxtal import pyxtal
    from pyxtal.representation import representation

    smiles = [
        "C1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])C(=O)O.smi",
        "CC1=CC2=C(C=C1)N3CC4=C(C=CC(=C4)C)N(C2)C3.smi",
    ]

    for block in [None, "xxv"]:
        print("block", block)
        for i in range(30):
            s = pyxtal(molecular=True)
            s.from_random(3, 14, smiles, block=block)

            rep2 = representation.from_pyxtal(s)
            strs = rep2.to_string()
            print(i, strs)

            rep3 = representation.from_string(strs, smiles)
            s1 = rep3.to_pyxtal()

            pmg1 = s.to_pymatgen()
            pmg2 = s1.to_pymatgen()
            pmg1.remove_species("H")
            pmg2.remove_species("H")
            if not sm.StructureMatcher().fit(pmg1, pmg2):
                print("Mismatch")
                print("S1", s.check_short_distances())
                print("S2", s1.check_short_distances())
                s.to_file("1.cif")
                s1.to_file("2.cif")
                import sys

                sys.exit()
