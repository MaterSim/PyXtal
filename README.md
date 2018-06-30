# Crystallography
Python code for the generation and study of crystal structures based on symmetry constraints.
Codeveloped by Scott Fredericks and Qiang Zhu at UNLV department of Physics.
Distributed under the MIT License

## Current Features:
..*Random generation of crystals for a given space group and stoichiometry
...Can outputs cif files for generated crystals (via Pymatgen)
..*Generation of molecular crystals, with consideration for each molecule's compatibility with the Wyckoff site symmetry
...Allows input from xyz files for molecules (via Pymatgen)
...Allows generation of simple molecules from chemical symbol (via ASE)
..*Allows easy access to Wyckoff position information, including site symmetry operations and symbols

##Dependencies:
[SciPy](https://www.scipy.org/install.html)
[NumPy](https://www.scipy.org/scipylib/download.html)
[Pandas](https://pandas.pydata.org/getpandas.html)
[Pymatgen](http://pymatgen.org/#getting-pymatgen)
[SpgLib for Python](https://atztogo.github.io/spglib/python-spglib.html#installation)
[ASE](https://wiki.fysik.dtu.dk/ase/install.html) (Only needed for reading chemical symbols)
