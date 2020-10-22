<img src="https://raw.githubusercontent.com/qzhu2017/PyXtal/master/images/512px_type1.png" alt="PyXtal" width="300"/>

[![Documentation Status](https://readthedocs.org/projects/pyxtal/badge/?version=latest)](https://pyxtal.readthedocs.io/en/latest/?badge=latest)
[![Test Status](https://github.com/qzhu2017/PyXtal/workflows/tests/badge.svg)](https://github.com/qzhu2017/PyXtal/actions)
[![Download Status](https://img.shields.io/pypi/pyversions/pyxtal)](https://pypi.org/project/pyxtal/)
[![Download Status](https://img.shields.io/pypi/v/pyxtal)](https://pypi.org/project/pyxtal/)
[![Downloads](https://pepy.tech/badge/pyxtal)](https://pepy.tech/project/pyxtal)
[![DOI](https://zenodo.org/badge/128165891.svg)](https://zenodo.org/badge/latestdoi/128165891)
<img align="right" width="450" src="https://raw.githubusercontent.com/qzhu2017/PyXtal/master/images/water.gif">

## Content

- [Content](#content)
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Current Features](#current-features)
- [Installation](#installation)
- [Citation](#citation)
- [How to contribute?](#how-to-contribute)
  - [If you just want to use the code](#if-you-just-want-to-use-the-code)
  - [If you want to join the code development](#if-you-want-to-join-the-code-development)

## Introduction

PyXtal is an open source Python package which was initiated by [Qiang Zhu](http://qzhu2017.github.io) and Scott Fredericks at department of Physics and Astronomy, University of Nevada Las Vegas. The goal of PyXtal project is to develop a fundamental library to allow one to design the material structure with a certain symmetry constraint. So far, the package allows for generation/manipulation of crystals, with both general and special Wyckoff positions. These structures can exported to various structural formats for further study. See the [documentation](https://pyxtal.readthedocs.io/en/latest/) for information about installation and usage.

To contribute to this project, please check [How to contribute?](#how-to-contribute).

## Quick Start

Check the folloowing links to quickly understand how pyxtal works

- [Atomic crystal](https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/01_atomic_crystals.ipynb)
- [Molecular crystal](https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/02_molecular_crystals.ipynb)

## Current Features

- Random generation of atomic/molecular crystals in 3D, 2D, 1D and 0D crystals for a given symmetry group and stoichiometry
- Molecules in special Wyckoff positions are automatically oriented to preserve the space group symmetry
- Interfaces with Pymatgen and ASE for structural manipulation and analysis
- Easy access to symmetry group information (e.g., Wyckoff positions, site symmetry operations)
- X-ray diffraction analysis and [its online application](https://vxrd.physics.unlv.edu)
- Geometry optimization via different exploratory algorithms.

## Installation

To install the code, one just needs to do

```sh
pip install pyxtal
```

or

```sh
pip install --upgrade git+https://github.com/qzhu2017/PyXtal.git@master
```

## Citation

Fredericks S, Sayre D, Zhu Q\*(2019)
[PyXtal: a Python Library for Crystal Structure Generation and Symmetry Analysis](https://arxiv.org/pdf/1911.11123.pdf)

```bib
@article{pyxtal,
    title={PyXtal: a Python Library for Crystal Structure Generation and Symmetry Analysis},
    author={Scott Fredericks and Dean Sayre and Qiang Zhu},
    year={2019},
    eprint={1911.11123},
    archivePrefix={arXiv},
    primaryClass={cond-mat.mtrl-sci}
}
```

## How to contribute?

This is an open-source project. Its growth depends on the community. To contribute to PyXtal, you don't necessarily have to write the code. Any contributions from the following list will be helpful.

### If you just want to use the code

- Star the PyXtal project via GitHub and recommend it to your colleagues/friends
- Open an issue to report the bug or address your wishlist
- Suggestions to improve our documentation

### If you want to join the code development

- Fork the repository
- Suggest and implement new functions
- Send us the pull request
