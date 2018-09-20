.. PyXtal documentation master file, created by
   sphinx-quickstart on Mon Aug 27 10:19:01 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to PyXtal’s documentation!
==================================

    PyXtal is a Python library for the ab-initio generation of random crystal structures. It is made available under the MIT license. Given a stoichiometry and space group, the user can quickly generate possible geometries, which can be output to .cif or .vasp files. The structure information can then be used in combination with various optimization methods and software, in order to determine the lowest-energy structure for a given compound. Currently, the software allows random generation of 3D, 2D, and 1D crystals. Both atomic and molecular crystals can be generated; PyXtal will automatically check molecules for their symmetry compatibility with special Wyckoff positions. The software also allows access to symmetry information, including Wyckoff positions and site symmetry for a given space group. A basic tutorial is provided below for common functions. Additionally, documentation and source code are provided for individual modules. For more information about the project's development, see the GitHub page: https://github.com/qzhu2017/PyXtal

Dependencies
============

Versions indicated are those used during development. Other versions may be compatible, but have not yet been tested.

  * `SciPy 1.0.1 <https://www.scipy.org/install.html>`_  
  * `NumPy 1.14.3 <https://www.scipy.org/scipylib/download.html>`_  
  * `Pandas 0.20.3 <https://pandas.pydata.org/getpandas.html>`_  
  * `Pymatgen 2017.9.3 <http://pymatgen.org/#getting-pymatgen>`_  
  * `SpgLib for Python 1.9.9.44 <https://atztogo.github.io/spglib/python-spglib.html#installation>`_  
  * `Openbabel 2.4.1 <http://openbabel.org/wiki/Category:Installation>`_  

Note that for openbabel, you must install the C++ pacakge before installing the Python bindings. For Debian based systems, your distribution may already have installable packages:

``sudo apt-get install openbabel``  

``pip install openbabel``

Note that the openbabel Python bindings require swig to install:

``sudo apt-get install swig``  

For other systems, you must compile the openbabel bindings yourself. There are tutorials for this on the `openbabel website
<https://openbabel.org/docs/dev/UseTheLibrary/PythonInstall.html>`_, as well as in the `pymatgen documentation
<http://pymatgen.org/installation.html#openbabel-mac-os-x-tested-on-v2-3-2>`_.

Installation
============

To install PyXtal, first install all dependencies, then make a copy of the source code:

``git clone https://github.com/qzhu2017/pyxtal``

Then, inside of the downloaded directory, run

``python setup.py install``

This will install the module. The code can be used within Python via

.. code-block:: Python

  import pyxtal

Usage
=====

The package’s main modules are `crystal.py <pyxtal.crystal.html>`_ and `molecular_crystal.py. <pyxtal.molecular_crystal.html>`_ These modules contain functions for generating 3d atomic and molecular crystals, respectively. These modules can be run as scripts via

``python -m pyxtal.crystal``

or

``python -m pyxtal.molecular_crystal``

As scripts, these modules will generate 3d crystals with a given space group, stoichiometry, and size. There are various command line options for these modules, which are listed in their documentation: `crystal.py <pyxtal.crystal.html>`_, `molecular_crystal.py. <pyxtal.molecular_crystal.html>`_. To explore the full options, please use the -h flag:

``python -m pyxtal.molecular_crystal -h``

Alternatively, the functions can be called directly within Python, as illustrated in the following examples.

Understanding the Parameters
----------------------------

There are 230 types of space groups (with International number between 1 and 230), each of which has a different set of symmetry operations and shape of unit cell. For a list of space groups and their symmetry operations, see the free online utility `WYCKPOS <http://www.cryst.ehu.es/cryst/get_wp.html>`_, which is maintained by the Bilbao Crystallographic Server. For more information about space groups and crystals in general, consult `the International Tables of Crystallography <https://it.iucr.org/>`_. For 2D crystals, we instead use layer groups, and for 1D crystals, we use Rod groups. A list of these groups and their operations can be accessed using the `WPOS utility
<http://www.cryst.ehu.es/subperiodic/get_sub_wp.html>`_ on Bilbao.

Note that the input number of atoms is for the primitive unit cell, while the output uses the conventional unit cell. This means that for some space groups, the number of atoms in the output structure will be multiplied by 2, 3, or 4 compared to the input. Specifically, space groups beginning with “A”, “C” or “I” have twice as many atoms in the conventional cell, space groups beginning with “R” have 3 times as many, and space groups beginning with “F” have 4 times as many. Keep this in mind when choosing a range of stoichiometries and space groups to generate.

The function works by attempting to insert atoms into different Wyckoff positions within the unit cell. If a set of inserted atoms is too close together, it will be merged into a lower-multiplicity position. Because each Wyckoff position can only hold a set number of atoms, not all stoichiometries will work with all space groups. In this case, the package will output a warning message, and will not attempt to generate the crystal. Even if the stoichiometry is technically compatible, some parameter sets may still be difficult or impossible to generate. If the atoms are still too close together after a set number of attempts, the generation will fail, and random_crystal.valid will be set to False.

If this happens consistently, or if generation takes an unreasonably long time, consider increasing the volume factor. This creates a larger unit cell with more space between atoms, making it easier to find a possible structure. However, increasing the volume factor too much will influence which Wyckoff positions are chosen, and may thus give less accurate results. A combination of physical insight and trial/error should be used to determine an appropriate value. Typically, a value between 1 and 10 should suffice.

3D Atomic Crystals
------------------

The main function for generating and storing crystals is `crystal.random_crystal <pyxtal.crystal.html#pyxtal.crystal.random_crystal>`_. Here is a simple example of a carbon crystal:

.. code-block:: Python

  from pyxtal.crystal import random_crystal
  my_crystal = random_crystal(225, ['C'], [3], 8.0)

This would create a crystal structure with space group 225, 3 carbon atoms in the primitive cell, and a volume factor of 8. For stoichiometries with more than one type of atom, replace [‘C’] with a list of atomic symbols, and replace [3] with a list of numbers. For example,

.. code-block:: Python

  my_crystal = random_crystal(99, ['Ba','Ti','O'], [1,1,3], 8.0)

would create a random BaTiO3 crystal.

If the generation is successful, the value my_crystal.valid will be set to True; otherwise, it will be False. The geometric properties of the crystal are stored in my_crystal.struct, which is a `pymatgen.core.structure.Structure <http://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.Structure>`_ object. You can print my_crystal.struct directly to view the lattice vectors and atomic positions, or you can output the information using

.. code-block:: Python

  my_crystal.struct.to(filename='myfile',fmt='cif')

This will create a file called myfile.cif, which stores the structure information. Other file types are supported by changing the value of fmt to ‘poscar’, ‘cssr’, or ‘json’.

3D Molecular Crystals
---------------------

Molecular 3d crystals are generated in the same way as atomic 3d crystals, but atomic species are replaced with (rigid) molecules. The molecules may either be strings for the chemical composition, or `pymatgen.core.structure.Molecule <http://pymatgen.org/pymatgen.core.structure.html?highlight=class%20molecule#pymatgen.core.structure.Molecule>`_ objects. The generating class is `molecular_crystal.molecular_crystal <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal>`_:

.. code-block:: Python

  from pyxtal.molecular_crystal import molecular_crystal
  my_crystal = molecular_crystal(36, ['H2O'], [4], 8.0)

This would give a crystal with spacegroup 36, 4 water molecules in the primitive cell, and a volume factor of 8. As with atomic crystals, you may use lists as input for the (molecular) stoichiometry.

Because molecules are less symmetric than individual atoms, they may or may not fit within a given Wyckoff position. Furthermore, the molecule may need to be oriented in a certain direction to be compatible with a site. The molecular_crystal class handles this automatically, and only inserts molecules in positions and orientations for which the molecules are sufficiently symmetric.

Like atomic crystals, the atomic positions may be accessed with the struct attribute. However, for accessing the positions and orientations of the molecules themselves, there is an attribute called mol_generators. This provides a list of `mol_site <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.mol_site>`_ objects, which in turn give the type, location, Wyckoff position, and orientation of each molecule in the asymmetric unit. This can be used to generate the crystal using molecules instead of indivual atoms. Note that the coordinates here are fractional, and refer to the molecule’s center of mass.

The orientations stored in the mol_site class are members of the `operations.orientation <pyxtal.operations.html#pyxtal.operations.orientation>`_ class. A molecule in a Wyckoff position may be allowed to rotate about a certain axis, allowed to rotate freely, or may be rigidly constrained. This information is stored in the orientation class. To obtain a SymmOp which can be applied to the molecule, and which is consistent with the geometric constraints, call `orientation.get_op <pyxtal.operations.html#pyxtal.operations.orientation.get_op>`_. For a 3x3 matrix instead, call `orientation.get_matrix <pyxtal.operations.html#pyxtal.operations.orientation.get_matrix>`_. In either case, this will give a random rotation consistent with the degrees of freedom. To obtain the exact rotation used when generating the crystal (and avoid the random rotation), pass the parameter angle=0.

There are a few other parameters which may be passed to the class. See the `module documentation <pyxtal.molecular_crystal.html>`_ for details. Of particular importance is the variable allow_inversion=False. By default, chiral molecules will not be flipped or inverted while generating the crystal. This is because a chiral molecule’s mirror image may have different chemical properties, especially in a biological setting. But if the mirror images are acceptable for your application, you may use allow_inversion=True, which will allow more spacegroups to be generated. Note that this is only relevant if at least one molecule is chiral.

2D Atomic Crystals
------------------

To generate a 2d crystal, use the class `crystal.random_crystal_2D <pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_. For example,

.. code-block:: Python

  from pyxtal.crystal import random_crystal_2D
  my_crystal = random_crystal_2D(20, ['C'], [4], 2.0, 2.5)

would generate a 2d crystal with layer group 20, 4 carbon atoms in the primitive cell, a thickness of 2.0 Angstroms, and a volume factor of 2.5. As with the 3d case, for crystals with multiple atom types, you may replace [‘C’] and [4] with lists of the atomic symbols and amounts, respectively. The crystal will be periodic in two directions instead of 3. Thus, it is not recommended to store the structure in a .cif or POSCAR file without software specifically designed for handling 2d crystals. The axis of non-periodicity can be accessed via my_crystal.PBC; a value of 1, 2, or 3 corresponds to the x, y, or z axis, respectively.

Note that the layer group number is different from the international space group number, and ranges between 1 and 80. For a list of the layer groups and their symmetry operations, see `the International Tables of Crystallography, Volume E, part 4 <https://it.iucr.org/Eb/ch4o1v0001/contents/>`_.

Note that changing the volume factor will not change the thickness of the output crystal. So, if you are testing over a range of volume factors, consider how the shape of the unit cell will be affected, and change the thickness accordingly. If you are unsure what value to use for the thickness, using None will automatically generate a value for you.

2D Molecular Crystals  
---------------------

2d Molecular crystals are generated using the class `molecular_crystal.molecular_crystal_2D <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_2D>`_:

.. code-block:: Python

  from pyxtal.molecular_crystal import molecular_crystal_2D
  my_crystal = molecular_crystal_2D(20, ['H2O'], [4], 2.0, 2.5)

Here, the parameters correspond to those for `random_crystal_2D <pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_, except the atoms are again replaced with molecules. The additional options available for `molecular_crystal <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal>`_ are also available for `molecular_crystal_2D <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_2D>`_.

Because molecules have a certain thickness of their own, care should be used when choosing a thickness value. Currently, the thickness parameter only determines where the molecular centers of mass can be, so the final crystal may have individual atoms outside of this range.

1D Crystals
-----------

PyXtal also supports generation of 1D crystals using Rod groups (between 1 and 75). The corresponding classes are `crystal.random_crystal_1D
<pyxtal.crystal.html#pyxtal.crystal.random_crystal_1D>`_ and `molecular_crystal_1D
<pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_1D>`_. The parameters for these functions are the same as those for `random_crystal_2D
<pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_ and `molecular_crystal_2D <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_2D>`_. However, in place of the thickness of the unit cell, you should use the cross-sectional area of the unit cell (in Angstroms squared). Again, using an area of None will generate a value automatically.

Wyckoff Positions
-----------------

The package has several utilities for working with Wyckoff positions and their symmetries. Most of these functions are found within the `crystal <pyxtal.crystal.html>`_ module.

The main function is `get_wyckoffs(sg) <pyxtal.crystal.html#pyxtal.crystal.get_wyckoffs>`_, where sg is the international space group of interest. This will return a list of Wyckoff positions for the space group, with the 0th list element being the general Wyckoff position, and the remaining elements being the special Wyckoff positions. Each of the Wyckoff positions is, in turn, a list of symmetry operations (`pymatgen.core.operations.SymmOp <http://pymatgen.org/pymatgen.core.operations.html#pymatgen.core.operations.SymmOp>`_). These symmetry operations can be applied to 3d vectors using op.operate(vector), or can be composed together via multiplication: op3 = op1 * op2. Each SymmOp consists of a rotation matrix (op.rotation_matrix) and a translation vector (op.translation), and is represented by a 4x4 affine matrix (op.affine_matrix).

A Wyckoff position is typically denoted with a number-letter combination, depending on its multiplicity. For example, for space group 16 (P222), the general Wyckoff position is called “4u”. This is because the position has a multiplicity of 4, and the letters a-t are used by special Wyckoff positions. If you know the letter of a Wyckoff position, you can obtain its index within a get_wyckoffs(sg) list using the function `index_from_letter(letter, sg) <pyxtal.crystal.html#pyxtal.crystal.index_from_letter>`_. Likewise, if you know the index, you can obtain the letter using `letter_from_index(index, sg) <pyxtal.crystal.html#pyxtal.crystal.letter_from_index>`_. An example is provided below.

For a given space group, each Wyckoff position is a subgroup of the general Wyckoff position. Furthermore, the Wyckoff position occupied by a given structure is determined by the symmetry of that structure with respect to the crystal as a whole. Naturally, each Wyckoff position requires some symmetry for an atom or molecule to occupy it. This symmetry can be accessed using `get_wyckoff_symmetry(sg) <pyxtal.crystal.html#pyxtal.crystal.get_wyckoff_symmetry>`_. This returns a nested list, where the first index specifies a Wyckoff position, the second index specifies a point within that Wyckoff position, and the third index specifies a list of symmetry operations corresponding to that point. This list of operations can then be used to check whether a given molecule is consistent with a given Wyckoff position.

In order to make this information more human-readable, the function `ss_string_from_ops(ops, sg) <pyxtal.crystal.html#pyxtal.crystal.ss_string_from_ops>`_ can be used to generate a Hermann-Mauguin symbol for a given point group. Here, ops is a list of symmetry options, and sg is the international spacegroup number. So, for example, if we wanted to know the site symmetry of the Wyckoff position 2b in space group 35, we could use:

.. code-block:: Python

  >>> from pyxtal.crystal import *
  >>> sg = 35
  >>> index = index_from_letter('b', sg)
  >>> ops = get_wyckoff_symmetry(sg)[index][0]
  >>> ss_string_from_ops(ops, sg)
  'mm2'

In this example, ‘mm2’ denotes mirror planes accross the x and y axes, and a 2-fold rotation about the z axis. Note that in Hermann-Mauguin notation, the symbols do not always follow this x,y,z format. For more information on reading these symbols, see https://en.wikipedia.org/wiki/Hermann%E2%80%93Mauguin_notation.

Space group settings
--------------------

For the output 3D structures, PyXtal uses the conventional standard cell (the same as `Bilbao
<http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-def-choice>`_). This means unique axis b for monoclinic cells, the obverse triple hexagonal cell for rhombohedral groups, and origin choice 2 (0,0,0) when two origin choices are available.

For 2D structures, we use unique axis c for monoclinic layer groups 3-7, and unique axis a for layer groups 8-18. When two origin choices are available, we use origin choice 1.

For 1D structures, we use unique axis a for monoclinic Rod groups 3-7, and unique axis c for Rod groups 8-12. When two settings are available for a group, we use the 1st setting.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
