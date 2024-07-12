#!/usr/bin/env  python

import os
from argparse import ArgumentParser

import numpy as np

from pyxtal import print_logo, pyxtal
from pyxtal.symmetry import get_symbol_and_number

if __name__ == "__main__":
    # -------------------------------- Options -------------------------
    parser = ArgumentParser()
    parser.add_argument(
        "-s",
        "--symmetry",
        dest="sg",
        metavar="sg",
        default=36,
        type=str,
        help="desired symmetry, number or string, e.g., 36, Pbca, Ih",
    )
    parser.add_argument(
        "-e",
        "--element",
        dest="element",
        default="Li",
        help="desired elements: e.g., Li",
        metavar="element",
    )
    parser.add_argument(
        "-n",
        "--numIons",
        dest="numIons",
        default=16,
        help="desired numbers of atoms: 16",
        metavar="numIons",
    )
    parser.add_argument(
        "-f",
        "--factor",
        dest="factor",
        default=1.0,
        type=float,
        help="volume factor: default 1.0",
        metavar="factor",
    )
    parser.add_argument(
        "-a",
        "--attempts",
        dest="attempts",
        default=1,
        type=int,
        help="number of crystals to generate: default 1",
        metavar="attempts",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        default="out",
        type=str,
        help="Directory for storing output cif files: default 'out'",
        metavar="outdir",
    )
    parser.add_argument(
        "-d",
        "--dimension",
        dest="dimension",
        metavar="dimension",
        default=3,
        type=int,
        help="desired dimension: (3, 2, 1, 0): default 3",
    )
    parser.add_argument(
        "-t",
        "--thickness",
        dest="thickness",
        metavar="thickness",
        default=None,
        type=float,
        help="Thickness of a 2D crystal, or area of a 1D crystal, None generates a value automatically: default None",
    )

    parser.add_argument(
        "-m",
        "--molecular",
        dest="molecular",
        action="store_true",
        default=False,
        help="molecular? default: False",
    )

    parser.add_argument(
        "-c",
        "--conventional",
        dest="conventional",
        action="store_true",
        default=False,
        help="conventional setting? default: False",
    )

    print_logo()
    options = parser.parse_args()
    sg = options.sg
    dimension = options.dimension
    if isinstance(sg, str) and sg.isnumeric():
        sg = int(sg)
    symbol, sg = get_symbol_and_number(sg, dimension)

    element = options.element
    number = options.numIons
    numIons = []
    if element.find(",") > 0:
        system = element.split(",")
        for x in number.split(","):
            numIons.append(int(x))
    else:
        system = [element]
        numIons = [int(number)]

    factor = options.factor
    if factor < 0:
        raise ValueError(f"Volume factor {factor:.2f} must be greater than 0.")

    attempts = options.attempts
    outdir = options.outdir
    dimension = options.dimension
    thickness = options.thickness
    molecular = options.molecular
    conventional = options.conventional

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for i in range(attempts):
        numIons0 = np.array(numIons)
        rand_crystal = pyxtal(molecular=molecular)
        if dimension == 3:
            rand_crystal.from_random(3, sg, system, numIons0, factor, conventional=conventional)
        elif dimension == 2:
            rand_crystal.from_random(2, sg, system, numIons0, factor, thickness, conventional=conventional)
        elif dimension == 1:
            rand_crystal.from_random(1, sg, system, numIons0, factor, thickness, conventional=conventional)
        if dimension == 0:
            rand_crystal.from_random(0, sg, system, numIons0, factor)
        # Output a cif or xyz file
        outpath = options.outdir + "/" + str(i) + ".cif" if dimension > 0 else options.outdir + "/" + str(i) + ".xyz"
        rand_crystal.to_file(filename=outpath)

        print(rand_crystal)
