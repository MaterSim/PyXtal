#!/usr/bin/env  python
# encoding: utf-8

from pyxtal import print_logo
from pyxtal.symmetry import Group
from argparse import ArgumentParser

if __name__ == "__main__":
    # -------------------------------- Options -------------------------
    parser = ArgumentParser()
    parser.add_argument(
        "-s",
        "--symmetry",
        dest="sg",
        type=str,
        help="desired symmetry, number or string, e.g., 36, Pbca, Ih. if None, show all list of available groups",
    )
    parser.add_argument(
        "-d",
        "--dimension",
        dest="dimension",
        default=3,
        type=int,
        help="desired dimension: (3, 2, 1, 0): default 3",
    )

    print_logo()
    options = parser.parse_args()

    dimension = options.dimension

    if options.sg is not None:
        sg = options.sg
        if sg.isnumeric():
            sg = int(sg)
        print(Group(sg, dimension))
    else:
        Group.list_groups(dimension)
