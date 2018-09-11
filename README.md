# PyXtal
Python code for the generation and study of crystal structures based on symmetry constraints.
Codeveloped by Qiang Zhu and Scott Fredericks at UNLV department of Physics.
Distributed under the MIT License.

## Current Features:
* Random generation of atomic (3d and 2d) crystals for a given space group and stoichiometry
* Random generation of molecular crystals (3d, with 2d under development), with consideration for each molecule's compatibility with the Wyckoff site symmetry
* Easy access to spacegroup and Wyckoff position information, including site symmetry operations and symbols

## Dependencies:
* [SciPy 1.0.1](https://www.scipy.org/install.html)
* [NumPy 1.14.3](https://www.scipy.org/scipylib/download.html)
* [Pandas 0.20.3](https://pandas.pydata.org/getpandas.html)
* [Pymatgen 2017.9.3](http://pymatgen.org/#getting-pymatgen)
* [SpgLib for Python 1.9.9.44](https://atztogo.github.io/spglib/python-spglib.html#installation)
* [openbabel 2.4.1](http://openbabel.org/wiki/Main_Page)
