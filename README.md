<img src="https://raw.githubusercontent.com/qzhu2017/PyXtal/master/images/512px_type1.png" alt="PyXtal" width="300"/>

[![Documentation Status](https://readthedocs.org/projects/pyxtal/badge/?version=latest)](https://pyxtal.readthedocs.io/en/latest/?badge=latest)
[![Test Status](https://github.com/qzhu2017/PyXtal/workflows/tests/badge.svg)](https://github.com/qzhu2017/PyXtal/actions)
[![Download Status](https://img.shields.io/pypi/pyversions/pyxtal)](https://pypi.org/project/pyxtal/)
[![Download Status](https://img.shields.io/pypi/v/pyxtal)](https://pypi.org/project/pyxtal/)
[![Downloads](https://pepy.tech/badge/pyxtal)](https://pepy.tech/project/pyxtal)
[![DOI](https://zenodo.org/badge/128165891.svg)](https://zenodo.org/badge/latestdoi/128165891)
[![Gitter](https://badges.gitter.im/PyXtal/community.svg)](https://gitter.im/PyXtal/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
<img align="right" width="450" src="https://raw.githubusercontent.com/qzhu2017/PyXtal/master/images/water.gif">


## Table of content
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Current Features](#current-features)
- [Installation](#installation)
- [Citation](#citation)
- [How to contribute?](#how-to-contribute)


## Introduction

PyXtal is an open source Python package which was initiated by [Qiang Zhu](http://qzhu2017.github.io) and Scott Fredericks at department of Physics and Astronomy, University of Nevada Las Vegas. The goal of PyXtal is to develop a fundamental library to allow one to design the material structure with a certain symmetry constraint. These structures can exported to various structural formats for further study. See the [documentation](https://pyxtal.readthedocs.io/en/latest/) for more information.

To contribute to this project, please check [How to contribute?](#how-to-contribute).

## Quick Start

- [Atomic crystal](https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/01_atomic_crystals.ipynb)
- [Molecular crystal](https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/02_molecular_crystals.ipynb)
- [XRD](https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/03_pxrd.ipynb)
- [Molecular Packing](https://nbviewer.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/05-crystal-packing.ipynb)

## Current Features

- Generation of atomic structures for a given symmetry and stoichiometry (0-3D)
- Generation of molecular crystals (1-3D) with the support of special Wyckoff positions.
- Structure manipulation via subgroup/supergroup symmetry relation
- Geometry optimization from built-in and external optimization methods.
- Internal support of ``cif`` file and many other formats via ``pymatgen`` or ``ASE``.
- Easy access to symmetry information (e.g., Wyckoff, site symmetry and international symbols).
- X-ray diffraction analysis and [its online application](https://vxrd.physics.unlv.edu)

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

Fredericks S, Parrish K, Sayre D, Zhu Q\*(2020)
[PyXtal: a Python Library for Crystal Structure Generation and Symmetry Analysis](https://www.sciencedirect.com/science/article/pii/S0010465520304057).
\[[arXiv link](https://arxiv.org/pdf/1911.11123.pdf)\]

```bib
@article{pyxtal,
title = "PyXtal: A Python library for crystal structure generation and symmetry analysis",
journal = "Computer Physics Communications",
volume = "261",
pages = "107810",
year = "2021",
issn = "0010-4655",
doi = "https://doi.org/10.1016/j.cpc.2020.107810",
url = "http://www.sciencedirect.com/science/article/pii/S0010465520304057",
author = "Scott Fredericks and Kevin Parrish and Dean Sayre and Qiang Zhu",
}
```

## How to contribute?

This is an open-source project. Its growth depends on the community. To contribute to PyXtal, you don't necessarily have to write the code. Any contributions from the following list will be helpful.

- [![Star on GitHub](https://img.shields.io/github/stars/qzhu2017/pyxtal.svg?style=social)](https://github.com/qzhu2017/pyxtal/stargazers)
 the PyXtal project and recommend it to your colleagues/friends
- Open an [![GitHub issues](https://img.shields.io/github/issues/qzhu2017/pyxtal.svg)](https://GitHub.com/qzhu2017/pyxtal/issues/) to report the bug or address your wishlist or improve our documentation
- [![GitHub forks](https://img.shields.io/github/forks/qzhu2017/pyxtal?style=social)](https://github.com/qzhu2017/PyXtal/network/members) the repository and send us the pull request
