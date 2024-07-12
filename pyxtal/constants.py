"""
Module to store the constants
"""

import numpy as np

from pyxtal.version import __version__

# Constants
# ------------------------------
tol_m = 0.3  # seperation tolerance in Angstroms
rad = np.pi / 180.0  # converting degrees to radians
deg = 180.0 / np.pi  # converting radians to degrees
pyxtal_verbosity = 1  # constant for printx function
# Max number of atoms per molecule before using fast distance check
max_fast_mol_size = 30
letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
ltype_keywords = [
    "triclinic",
    "Triclinic",
    "monoclinic",
    "Monoclinic",
    "orthorhombic",
    "Orthorhombic",
    "tetragonal",
    "Tetragonal",
    "trigonal",
    "Trigonal",
    "hexagonal",
    "Hexagonal",
    "cubic",
    "Cubic",
    "spherical",
    "Spherical",
    "ellipsoidal",
    "Ellipsoidal",
]
single_smiles = [
    "Cl-",
    "F-",
    "Br-",
    "I-",
    "Li+",
    "Na+",
    "Cs+",
    "Rb+",
    "[Cl-]",
    "[F-]",
    "[Br-]",
    "[I-]",
    "[Li+]",
    "[Na+]",
    "[Cs+]",
    "Rb+",
]
hex_cell = np.array([[1, -0.5, 0], [0, np.sqrt(3) / 2, 0], [0, 0, 1]])
all_sym_directions = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
    (1, -1, 0),
    (1, 1, 0),
    (0, 1, -1),
    (0, 1, 1),
    (-1, 0, 1),
    (1, 0, 1),
    (1, -2, 0),
    (2, -1, 0),
]

logo = rf"""#############################################################
#             ______       _    _          _                #
#            (_____ \     \ \  / /        | |               #
#             _____) )   _ \ \/ / |_  ____| |               #
#            |  ____/ | | | )  (|  _)/ _  | |               #
#            | |    | |_| |/ /\ \ |_( (_| | |___            #
#            |_|     \__  /_/  \_\___)__|_|_____)           #
#                   (____/                                  #
#---------------------(version {__version__:>8s})--------------------#
#       A Python package for random crystal generation      #
#       url: https://github.com/qzhu2017/pyxtal             #
#       @Zhu's group at University of Nevada Las Vegas      #
#############################################################
"""
