.. PyXtal documentation master file, created by
   sphinx-quickstart on Mon Aug 27 10:19:01 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ../images/512px_type1.png
   :height: 512 px
   :width: 903 px
   :scale: 35 %
   :alt: PyXtal
   :align: center

Introduction
============

    PyXtal (pronounced "pie-crystal") is an open source Python library for the ab-initio generation of random crystal structures. It is available for use under the MIT license. Given a stoichiometry and space group, the user can quickly generate possible geometries, which can be output to .cif or .vasp files. The structure information can then be used in combination with various optimization methods and software, in order to determine the lowest-energy structure for a given compound. Currently, the software allows random generation of 3D, 2D, and 1D crystals, as well as point group clusters. Both atomic and molecular crystals can be generated; PyXtal will automatically check molecules for their symmetry compatibility with special Wyckoff positions. The software also allows access to symmetry information, including Wyckoff positions and site symmetry for a given space group. A basic tutorial is provided below for common functions. Additionally, documentation and source code are provided for individual modules. For more information about the project's development, see the GitHub page: https://github.com/qzhu2017/PyXtal

Version info
============
The current version is 0.1dev. Expect frequent updates.

Additional Documentation
========================

.. toctree::
   Installation
   Usage
   Settings
   :maxdepth: 2