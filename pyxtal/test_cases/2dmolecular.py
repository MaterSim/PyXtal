"""
Script to test the functionality of molecular_crystal_2D. For each layer group
(between 1 and 80), generates a 2d molecular crystal with the given parameters.
Most of the command line options for molecular_crystal.py are also present here,
with the addition of -t for the thickness of the cell. For each layer group, the
number of molecules used to generate the crystal is equal to the multiplicity
of the general position.
"""

if __name__ == "__main__":
    # -------------------------------- Options -------------------------
    from os import mkdir
    from pyxtal.molecular_crystal import *
    from pyxtal.database.layergroup import Layergroup

    parser = OptionParser()
    parser.add_option(
        "-e",
        "--molecule",
        dest="molecule",
        default="H2O",
        help="desired molecules: e.g., H2O",
        metavar="molecule",
    )
    parser.add_option(
        "-n",
        "--numMols",
        dest="numMols",
        default=12,
        help="desired numbers of molecules: 12",
        metavar="numMols",
    )
    parser.add_option(
        "-t",
        "--thickness",
        dest="thickness",
        default=3.0,
        type=float,
        help="volume factor: default 3.0",
        metavar="thickness",
    )
    parser.add_option(
        "-f",
        "--factor",
        dest="factor",
        default=1.0,
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
    parser.add_option(
        "-o",
        "--outdir",
        dest="outdir",
        default="out",
        type=str,
        help="Directory for storing output cif files: default 'out'",
        metavar="outdir",
    )
    parser.add_option(
        "-c",
        "--checkatoms",
        dest="checkatoms",
        default="True",
        type=str,
        help="Whether to check inter-atomic distances at each step: default True",
        metavar="outdir",
    )
    parser.add_option(
        "-i",
        "--allowinversion",
        dest="allowinversion",
        default="False",
        type=str,
        help="Whether to allow inversion of chiral molecules: default False",
        metavar="outdir",
    )

    (options, args) = parser.parse_args()
    molecule = options.molecule
    number = options.numMols
    verbosity = options.verbosity
    attempts = options.attempts
    outdir = options.outdir
    if options.checkatoms == "True" or options.checkatoms == "False":
        checkatoms = eval(options.checkatoms)
    else:
        print("Invalid value for -c (--checkatoms): must be 'True' or 'False'.")
        checkatoms = True
    if options.allowinversion == "True" or options.allowinversion == "False":
        allowinversion = eval(options.allowinversion)
    else:
        print("Invalid value for -i (--allowinversion): must be 'True' or 'False'.")
        allowinversion = False

    numMols = []
    if molecule.find(",") > 0:
        strings = molecule.split(",")
        system = []
        for mol in strings:
            system.append(mol_from_collection(mol))
        for x in number.split(","):
            numMols.append(int(x))
    else:
        system = [mol_from_collection(molecule)]
        numMols = [int(number)]
    orientations = None

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
            numMols0 = np.array(numMols)
            rand_crystal = molecular_crystal_2D(
                num,
                system,
                [multiplicity],
                options.thickness,
                options.factor,
                orientations=orientations,
                check_atomic_distances=checkatoms,
                allow_inversion=allowinversion,
            )
            end = time()
            timespent = np.around((end - start), decimals=2)
            if rand_crystal.valid:
                """written = False
                try:
                    mkdir(outdir)
                except: pass
                try:
                    comp = str(rand_crystal.struct.composition)
                    comp = comp.replace(" ", "")
                    cifpath = outdir + '/' + comp + "_" + str(i+1) + '.cif'
                    CifWriter(rand_crystal.struct, symprec=0.1).write_file(filename = cifpath)
                    written = True
                except: pass"""

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
                        + " Could not calculate generated.***********"
                    )
                """if written is True:
                    print("    Output to "+cifpath)
                else:
                    print("    Could not write cif file.")"""

                # Print additional information about the structure
                if verbosity > 0:
                    print("Time required for generation: " + str(timespent) + "s")
                    print("Molecular Wyckoff positions:")
                    for ms in rand_crystal.mol_generators:
                        print(
                            str(ms.mol.composition)
                            + ": "
                            + str(ms.multiplicity)
                            + str(ms.letter)
                            + " "
                            + str(ms.position)
                        )
                if verbosity > 1:
                    print(rand_crystal.struct)

            # If generation fails
            else:
                print("something is wrong")
                print("Time spent during generation attempt: " + str(timespent) + "s")

