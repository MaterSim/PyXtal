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

PyXtal
============

PyXtal (pronounced ``pie-crystal``) is an open source Python library for the ab-initio generation of random crystal structures. It has the following features:

- Generation of atomic crystals and clusters for a given symmetry and stoichiometry (0D, 1D, 2D and 3D)
- Generation of molecular crystals (1-3D), with the automatic determination of special Wyckoff positions. 
- Internal support of ``cif`` file and many other formats (e.g., ``xyz``, ``POSCAR``) via ``pymatgen`` and ``ASE``.
- Easy access to symmetry group information (e.g., Wyckoff, site symmetry, and point group symbols).
- Geometry optimization from built-in and external optimization methods.

The current version is ``0.0.6`` at `GitHub <https://github.com/qzhu2017/PyXtal>`_. It is available for use under the MIT license.

Expect updates upon request by `Qiang Zhu\'s group <http://www.physics.unlv.edu/~qzhu/index.html>`_ at University of Nevada Las Vegas.

PyXtal is an open source project. You are welcome to contribute it directly via the `GitHub platform <https://github.com/qzhu2017/PyXtal>`_ or send your comments and suggestions to the `developer <http://www.physics.unlv.edu/~qzhu/>`_.


A basic tutorial is provided below for common functions. Additionally, documentation and source code are provided for individual modules. 

Tutorial
========================

.. toctree::
   Installation
   COMMAND_MODE
   Usage
   Others
   Background
   Algorithm
   Settings
   :maxdepth: 2

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Modules
========================

.. toctree::
   :maxdepth: 4

   pyxtal
