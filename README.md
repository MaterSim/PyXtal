<img src="https://raw.githubusercontent.com/MaterSim/PyXtal/master/images/512px_type1.png" alt="PyXtal" width="300"/>

[![Documentation Status](https://readthedocs.org/projects/pyxtal/badge/?version=latest)](https://pyxtal.readthedocs.io/en/latest/?badge=latest)
[![Test Status](https://github.com/MaterSim/PyXtal/workflows/tests/badge.svg)](https://github.com/MaterSim/PyXtal/actions)
[![Download Status](https://img.shields.io/pypi/pyversions/pyxtal)](https://pypi.org/project/pyxtal/)
[![Download Status](https://img.shields.io/pypi/v/pyxtal)](https://pypi.org/project/pyxtal/)
[![Downloads](https://pepy.tech/badge/pyxtal)](https://pepy.tech/project/pyxtal)
[![DOI](https://zenodo.org/badge/128165891.svg)](https://zenodo.org/badge/latestdoi/128165891)
[![Gitter](https://badges.gitter.im/PyXtal/community.svg)](https://gitter.im/PyXtal/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
<img align="right" width="450" src="https://raw.githubusercontent.com/MaterSim/PyXtal/master/images/water.gif">


## Table of content
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Current Features](#current-features)
- [Installation](#installation)
- [Citation](#citation)
- [How to contribute?](#how-to-contribute)


## Introduction

PyXtal is an open source Python package which was initiated by [Qiang Zhu](http://qzhu2017.github.io) and Scott Fredericks. The goal of PyXtal is to develop a fundamental library to allow one to design the material structure with a certain symmetry constraint. These structures can exported to various structural formats for further study. See the [documentation](https://pyxtal.readthedocs.io/en/latest/) for more information.

To contribute to this project, please check [How to contribute?](#how-to-contribute).

## Quick Start

- [Atomic crystal](https://nbviewer.jupyter.org/github/MaterSim/PyXtal/blob/master/examples/tutorials_notebook/01_atomic_crystals.ipynb)
- [Molecular crystal](https://nbviewer.jupyter.org/github/MaterSim/PyXtal/blob/master/examples/tutorials_notebook/02_molecular_crystals.ipynb)
- [XRD](https://nbviewer.jupyter.org/github/MaterSim/PyXtal/blob/master/examples/tutorials_notebook/03_pxrd.ipynb)
- [Molecular Packing](https://nbviewer.org/github/MaterSim/PyXtal/blob/master/examples/tutorials_notebook/05-crystal-packing.ipynb)

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
pip install --upgrade git+https://github.com/MaterSim/PyXtal.git@master
```

If you want to add the Julia package install (required by the use of `pyxtal.lego` module), please use

```sh
export INSTALL_JULIA=1 && pip install pyxtal 
```

To check if the installation is successful, run the following script,
```python

from juliacall import Main as jl
jl.seval("using CrystalNets")
print("Success")
```

## Citation

### General PyXtal

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

### Site Symmetry Representation
Levy D, Panigrahi S-S, Kaba S-O, Zhu Q, Galkin M, Miret S, Ravanbakhsh S. (2024) 
[SymmCD: Symmetry-Preserving Crystal Generation with Diffusion Models AI for Accelerated Materials Design](https://openreview.net/forum?id=V7x2KZQn2v)
NeurIPS 2024, 

```bib
@inproceedings{
levy2024symmcd,
title={Symm{CD}: Symmetry-Preserving Crystal Generation with Diffusion Models},
author={Daniel Levy and Siba Smarak Panigrahi and S{\'e}kou-Oumar Kaba and Qiang Zhu and Mikhail Galkin and Santiago Miret and Siamak Ravanbakhsh},
booktitle={AI for Accelerated Materials Design - NeurIPS 2024},
year={2024},
url={https://openreview.net/forum?id=V7x2KZQn2v}
}
```


### Symmetry relation
Zhu Q, Kang B, Parrish K (2022). 
[Symmetry Relation Database and Its Application to Ferroelectric Materials Discovery](https://link.springer.com/article/10.1557/s43579-022-00268-4)

```bib
@article{zhu2022symmetry,
  title="Symmetry relation database and its application to ferroelectric materials discovery",
  author="Zhu, Qiang and Kang, Byungkyun and Parrish, Kevin",
  journal="MRS Communications",
  volume="12",
  number="5",
  pages="686--691",
  year="2022",
  doi="https://doi.org/10.1557/s43579-022-00268-4",
```

### Organic crystal packing motif
Zhu Q, Tang W-L, Hattori S. (2022). 
[Quantification of Crystal Packing Similarity from Spherical Harmonic Transform](https://pubs.acs.org/doi/10.1021/acs.cgd.2c00933)

```bib
@article{zhu2022quantification,
  title="Quantification of Crystal Packing Similarity from Spherical Harmonic Transform",
  author="Zhu, Qiang and Tang, Weilun and Hattori, Shinnosuke",
  journal="Crystal Growth \& Design",
  volume="22",
  number="12",
  pages="7308--7316",
  year="2022",
  doi="https://doi.org/10.1021/acs.cgd.2c00933",
}
```

## How to contribute?

This is an open-source project. Its growth depends on the community. To contribute to PyXtal, you don't necessarily have to write the code. Any contributions from the following list will be helpful.

- [![Star on GitHub](https://img.shields.io/github/stars/qzhu2017/pyxtal.svg?style=social)](https://github.com/MaterSim/pyxtal/stargazers)
 the PyXtal project and recommend it to your colleagues/friends
- Open an [![GitHub issues](https://img.shields.io/github/issues/qzhu2017/pyxtal.svg)](https://GitHub.com/MaterSim/pyxtal/issues/) to report the bug or address your wishlist or improve our documentation
- [![GitHub forks](https://img.shields.io/github/forks/qzhu2017/pyxtal?style=social)](https://github.com/MaterSim/PyXtal/network/members) the repository and send us the pull request
