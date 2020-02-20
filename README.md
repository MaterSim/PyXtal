<img src="https://raw.githubusercontent.com/qzhu2017/PyXtal/master/images/512px_type1.png" alt="PyXtal" width="300"/>

[![Documentation Status](https://readthedocs.org/projects/pyxtal/badge/?version=latest)](https://pyxtal.readthedocs.io/en/latest/?badge=latest)
[![Download Status](https://img.shields.io/pypi/dm/pyxtal)](https://pypi.org/project/pyxtal/)
[![Download Status](https://img.shields.io/pypi/pyversions/pyxtal)](https://pypi.org/project/pyxtal/)
[![Download Status](https://img.shields.io/pypi/v/pyxtal)](https://pypi.org/project/pyxtal/)
## Introduction
PyXtal is an open source Python package for generating random crystal structures based on symmetry constraints. The package allows for generation of both atomic and molecular crystals, with both general and special Wyckoff positions. These structures can be output to cif files for optimization and study. The package also allows users to easily access symmetry group information. See the [documentation](https://pyxtal.readthedocs.io/en/latest/) for information about installation and usage.

Codeveloped by Qiang Zhu and Scott Fredericks at department of Physics and Astronomy, University of Nevada Las Vegas.
Distributed under the MIT License.

## Current Features:
* Random generation of atomic (3D, 2D, and 1D) crystals and point group clusters for a given symmetry group and stoichiometry
* Random generation of (rigid) molecular crystals (3D, 2D, and 1D), with special Wyckoff positions. Molecules in special Wyckoff positions are automatically oriented to preserve the space group symmetry
* Structure output to cif or POSCAR files via pymatgen
* Easy access to symmetry group information, including Wyckoff positions, site symmetry operations, and point group symbols

## Dependencies:
* [SciPy 1.0.1](https://www.scipy.org/install.html)
* [NumPy 1.14.3](https://www.scipy.org/scipylib/download.html)
* [Pandas 0.20.3](https://pandas.pydata.org/getpandas.html)
* [Pymatgen 2017.9.3](http://pymatgen.org/#getting-pymatgen)
* [SpgLib for Python 1.9.9.44](https://atztogo.github.io/spglib/python-spglib.html#installation)

### Optional:
* [openbabel 2.4.1 (Python bindings)](http://openbabel.org/wiki/Main_Page) (allows for additional molecule file formats)

## Installation

```
$ pip install pyxtal
```
or 
```
$ pip install --upgrade git+https://github.com/qzhu2017/PyXtal.git@master
```
