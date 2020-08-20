"""
Script to test the functionality of random_crystal_2D. For each layer
group (between 1 and 80), generates a 2d atomic crystal with the given
parameters. Most of the command line options for molecular_crystal.py
are also present here, with the addition of -t for the thickness of the
cell. For each layer group, the number of atoms used to generate the
crystal is equal to the multiplicity of the general position.
"""

if __name__ == "__main__":
    # -------------------------------- Options -------------------------
    from time import time
    from os import mkdir
    from pyxtal.crystal import *
    from pyxtal.database.layergroup import Layergroup

    parser = OptionParser()
    parser.add_option(
        "-e",
        "--atoms",
        dest="atoms",
        default="C",
        help="desired atoms: e.g., C",
        metavar="atoms",
    )
    parser.add_option(
        "-n",
        "--numions",
        dest="numions",
        default=12,
        help="desired numbers of ions: 12",
        metavar="numions",
    )
    parser.add_option(
        "-t",
        "--thickness",
        dest="thickness",
        default=3.0,
        type=float,
        help="volume factor: default 1.0",
        metavar="thickness",
    )
    parser.add_option(
        "-f",
        "--factor",
        dest="factor",
        default=4.0,
        type=float,
        help="volume factor: default 1.0",
        metavar="factor",
    )
    parser.add_option(
        "-v",
        "--verbosity",
        dest="verbosity",
        default=0,
        type=int,
        help="verbosity: default 0; higher values print more information",
        metavar="verbosity",
    )
    parser.add_option(
        "-a",
        "--attempts",
        dest="attempts",
        default=1,
        type=int,
        help="number of crystals to generate: default 1",
        metavar="attempts",
    )

    (options, args) = parser.parse_args()
    atoms = options.atoms
    number = options.numions
    verbosity = options.verbosity
    attempts = options.attempts

    numions = []
    if atoms.find(",") > 0:
        strings = atoms.split(",")
        system = []
        for atom in strings:
            system.append(atom)
        for x in number.split(","):
            numions.append(int(x))
    else:
        system = [atoms]
        numions = [int(number)]

    # Layergroup numbers to test
    numrange = list(range(1, 81))

    for num in numrange:
        print("---------------Layergroup " + str(num) + "---------------")
        for i in range(attempts):
            start = time()
            sg = Layergroup(num).sgnumber
            multiplicity = len(
                get_wyckoffs(sg)[0]
            )  # multiplicity of the general position
            rand_crystal = random_crystal_2D(
                num, system, [multiplicity], options.thickness, options.factor
            )
            end = time()
            timespent = np.around((end - start), decimals=2)
            if rand_crystal.valid:

                # spglib style structure called cell
                ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                sg = Layergroup(num).sgnumber
                if ans is not None:
                    print(
                        "Space group requested: " + str(sg) + " generated",
                        ans["number"],
                        "vol: ",
                        rand_crystal.volume,
                    )
                else:
                    print(
                        "Space group requested: "
                        + str(sg)
                        + " Could not calculate generated"
                    )

                # Print additional information about the structure
                if verbosity > 0:
                    print("Time required for generation: " + str(timespent) + "s")
                    print(rand_crystal.struct)

            # If generation fails
            else:
                print("something is wrong")
                print("Time spent during generation attempt: " + str(timespent) + "s")
