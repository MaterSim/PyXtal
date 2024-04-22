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

PyXtal (pronounced ``pie-crystal``) is an open source Python library for the
ab-initio generation of random crystal structures. It has the following features:

- Generation of atomic structures for a given symmetry and stoichiometry (0-3D)
- Generation of molecular crystals (1-3D) with the support of special Wyckoff positions
- Internal support of ``cif`` file and many other formats via ``pymatgen`` or ``ASE``
- Easy access to symmetry information (e.g., Wyckoff, site symmetry and international symbols)
- X-ray diffraction analysis and its `online application <https://vxrd.physics.unlv.edu>`_
- Structural manipulation via symmetry relation (both subgroup and supergroup)
- Geometry optimization from built-in and external optimization methods

The current version is ``0.6.5`` at `GitHub <https://github.com/qzhu2017/PyXtal>`_.
It is available for use under the MIT license. Expect updates upon request by
`Qiang Zhu\'s group <https://qzhu2017.github.io>`_ at the
University of North Carolina at Charlotte.

PyXtal is an open source project. You are welcome to contribute it directly via
the `GitHub platform <https://github.com/qzhu2017/PyXtal>`_ or send your
comments and suggestions to the `developer <http://www.physics.unlv.edu/~qzhu/>`_.


A basic tutorial is provided below for common functions. Additionally, 
documentation and source code are provided for individual modules.

For an experienced Python user who are familiar with Jupyter notebook, you are
also encouraged to check out the following examples

- `Atomic crystal <https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/01_atomic_crystals.ipynb>`_
- `Molecular crystal <https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/02_molecular_crystals.ipynb>`_
- `XRD <https://nbviewer.jupyter.org/github/qzhu2017/PyXtal/blob/master/examples/tutorials_notebook/03_pxrd.ipynb>`_

Tutorial
========================

.. toctree::
   Installation
   COMMAND_MODE
   Usage
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
