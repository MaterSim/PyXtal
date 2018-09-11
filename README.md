# pyXtal
Python code for the generation and study of crystal structures based on symmetry constraints.
Codeveloped by Qiang Zhu and Scott Fredericks at UNLV department of Physics.
Distributed under the MIT License.

## Current Features:
* Random generation of atomic (3d and 2d) crystals for a given space group and stoichiometry
* Random generation of molecular crystals (3d, with 2d under development), with consideration for each molecule's compatibility with the Wyckoff site symmetry
* Easy access to spacegroup and Wyckoff position information, including site symmetry operations and symbols

## Dependencies:
* [SciPy](https://www.scipy.org/install.html)
* [NumPy](https://www.scipy.org/scipylib/download.html)
* [Pandas](https://pandas.pydata.org/getpandas.html)
* [Pymatgen](http://pymatgen.org/#getting-pymatgen)
* [SpgLib for Python](https://atztogo.github.io/spglib/python-spglib.html#installation)
* [openbabel](http://openbabel.org/wiki/Main_Page)
