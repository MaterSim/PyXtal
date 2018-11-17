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
==================================

    PyXtal (pronounced "pie-crystal") is an open source Python library for the ab-initio generation of random crystal structures. It is available for use under the MIT license. Given a stoichiometry and space group, the user can quickly generate possible geometries, which can be output to .cif or .vasp files. The structure information can then be used in combination with various optimization methods and software, in order to determine the lowest-energy structure for a given compound. Currently, the software allows random generation of 3D, 2D, and 1D crystals, as well as point group clusters. Both atomic and molecular crystals can be generated; PyXtal will automatically check molecules for their symmetry compatibility with special Wyckoff positions. The software also allows access to symmetry information, including Wyckoff positions and site symmetry for a given space group. A basic tutorial is provided below for common functions. Additionally, documentation and source code are provided for individual modules. For more information about the project's development, see the GitHub page: https://github.com/qzhu2017/PyXtal

Dependencies
============

Versions indicated are those used during development. Other versions may be compatible, but have not yet been tested.

  * `SciPy 1.0.1 <https://www.scipy.org/install.html>`_  
  * `NumPy 1.14.3 <https://www.scipy.org/scipylib/download.html>`_  
  * `Pandas 0.20.3 <https://pandas.pydata.org/getpandas.html>`_  
  * `Pymatgen 2017.9.3 <http://pymatgen.org/#getting-pymatgen>`_  
  * `SpgLib for Python 1.9.9.44 <https://atztogo.github.io/spglib/python-spglib.html#installation>`_  

Optional:

  * `Openbabel 2.4.1 <http://openbabel.org/wiki/Category:Installation>`_  

Openbabel is not necessary, and only adds additional file format support for importing molecules. You must install the C++ pacakge before installing the Python bindings. For Debian based systems, your distribution may already have installable packages:

``sudo apt-get install openbabel``  

``pip install openbabel``

Also note that the openbabel Python bindings require swig to install:

``sudo apt-get install swig``  

Older version of swig (before 2.0) will not work. For other systems, you must compile the openbabel Python bindings yourself. There are tutorials for this on the `openbabel website
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

The current version is 0.1dev. Expect frequent updates.

Usage
=====

The package’s main modules are `crystal.py <pyxtal.crystal.html>`_ and `molecular_crystal.py. <pyxtal.molecular_crystal.html>`_ These modules contain functions for generating 3d atomic and molecular crystals, respectively. These modules can be run as scripts via

``python -m pyxtal.crystal``

or

``python -m pyxtal.molecular_crystal``

As scripts, these modules will generate 3d crystals with a given space group, stoichiometry, and size. There are various command line options for these modules, which are listed in their documentation: `crystal.py <pyxtal.crystal.html>`_, `molecular_crystal.py <pyxtal.molecular_crystal.html>`_. To explore the full options, please use the -h flag:

``python -m pyxtal.molecular_crystal -h``

Alternatively, the functions can be called directly within Python, as illustrated in the following examples.

Understanding the Parameters
----------------------------

There are 230 types of space groups (with International number between 1 and 230), each of which has a different set of symmetry operations and shape of unit cell. For a list of space groups and their symmetry operations, see the free online utility `WYCKPOS <http://www.cryst.ehu.es/cryst/get_wp.html>`_, which is maintained by the Bilbao Crystallographic Server. For more information about space groups and crystals in general, consult `the International Tables of Crystallography <https://it.iucr.org/>`_. For 2D crystals, we instead use layer groups, and for 1D crystals, we use Rod groups. For clusters, we use point groups. A list of these groups and their operations can be accessed using the `WPOS utility
<http://www.cryst.ehu.es/subperiodic/get_sub_wp.html>`_ on Bilbao.

Note that the input number of atoms is for the primitive unit cell, while the output uses the conventional unit cell. This means that for some space groups, the number of atoms in the output structure will be multiplied by 2, 3, or 4 compared to the input. Specifically, space groups beginning with “A”, “C” or “I” have twice as many atoms in the conventional cell, space groups beginning with “R” have 3 times as many, and space groups beginning with “F” have 4 times as many. Keep this in mind when choosing a range of stoichiometries and space groups to generate.

The function works by attempting to insert atoms into different Wyckoff positions within the unit cell. If a set of inserted atoms is too close together, it will be merged into a lower-multiplicity position. Because each Wyckoff position can only hold a set number of atoms, not all stoichiometries will work with all space groups. In this case, the package will output a warning message, and will not attempt to generate the crystal. Even if the stoichiometry is technically compatible, some parameter sets may still be difficult or impossible to generate. If the atoms are still too close together after a set number of attempts, the generation will fail, and random_crystal.valid will be set to False.

If this happens consistently, or if generation takes an unreasonably long time, consider increasing the volume factor. This creates a larger unit cell with more space between atoms, making it easier to find a possible structure. However, increasing the volume factor too much will influence which Wyckoff positions are chosen, and may thus give less accurate results. A combination of physical insight and trial/error should be used to determine an appropriate value. Typically, a value between 1.0 and 3.0 should suffice.

Lattices
~~~~~~~~

It is possible to supply your own unit cell lattice for a random crystal, via the `Lattice <pyxtal.crystal.html#pyxtal.crystal.Lattice>`_ class. You can define a lattice using either a 3x3 matrix, or using the lattice parameters:

.. code-block:: Python

  from pyxtal.crystal import Lattice
  l1 = Lattice.from_matrix([[4.08,0,0],[0,9.13,0],[0,0,5.50]])
  l2 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)

Here, both l1 and l2 describe the same lattice. In this case, it is an orthorhombic lattice with side lengths 4.08 Angstroms, 9.13 Angstroms, and 5.50 Angstroms, which is the unit cell for common water ice. The lattice parameters are, in order: (a, b, c, alpha, beta, gamma). a, b, and c are the lengths of the lattice vectors; alpha, beta, and gamma are the angles (in degrees) between these vectors. You can use a custom Lattice to generate a random_crystal or molecular_crystal:

.. code-block:: Python
 
  from pyxtal.molecular_crystal import molecular_crystal
  my_crystal = molecular_crystal(36, ['H2O'], [2], 1.0, lattice=l1)

This would generate a random water ice crystal, with space group 36, 4 molecules in the conventional cell (2 in the primitive cell), and using the lattice which we specified above. If you do not specify a lattice, a random one will be generated which is consistent with the chosen space group.

Note: For monoclinic layer groups, be careful when choosing the unique axis (see the "Group settings" section below).

Tolerance Matrices
~~~~~~~~~~~~~~~~~~

When generating random crystals, PyXtal performs inter-atomic distances checks to make sure the atoms are not too close together. By default, the covalent radius is used as a basis. However, the user may also define their own criteria using the `Tol_matrix <pyxtal.crystal.html#pyxtal.crystal.Tol_matrix>`_ class. To do this, initialize a Tol_matrix object using one of the built-in methods (see the Tol_matrix class documentation linked above for details):

.. code-block:: Python

  from pyxtal.crystal import Tol_matrix
  tol_m_1 = Tol_matrix(prototype="molecular", factor=2.0)
  tol_m_2 = Tol_matrix.from_radii(some_custom_list_of_atomic_radii)
  tol_m_3 = Tol_matrix.from_matrix(some_custom_2D_tolerance_matrix)

From here, you can alter the tolerance between certain inter-atomic pairs. Additionally, you can save and reload custom Tol_matrix objects for later use:

.. code-block:: Python

  >>> tol_m_1.set_tol('C', 'N', 2.1)
  >>> tol_m_1.set_tol(1, 3, 4.6)
  >>> tol_m_1.to_file("custom_matrix_file")
  'Output file to custom_matrix_file.npy'
  >>> reloaded_tol_matrix = Tol_matrix.from_file("custom_matrix_file.npy")
  >>> reloaded_tol_matrix.print_all()
  --Tol_matrix class object--
    Prototype: molecular
    Atomic radius type: covalent
    Radius scaling factor: 2.4
    Custom tolerance values:
      C, N: 2.1
      H, Li: 4.6

The Tol_matrix can now be passed to a random_crystal object:

.. code-block:: Python

  custom_tolerance_crystal = random_crystal(12, ['C','N'], [2,4], 1.0, tm=tol_m_1)

Alternatively, you can specify one of the preset tolerance matrices by passing a string to random_crystal or molecular_crystal. Possible values include "atomic", "molecular", or "metallic":

.. code-block:: Python

  metallic_crystal = random_crystal(12, ['Cu', 'Pd'], [2, 4], 1.0, tm="metallic")

By default, atomic crystals will use the average of the covalent radii between two atoms. Molecular crystals will use 1.2 times the sum of the covalent radii between two atoms. Using "metallic" will use the average of the metallic radius for metals, and the covalent radius for other atom types.

3D Atomic Crystals
------------------

The main function for generating and storing crystals is `crystal.random_crystal <pyxtal.crystal.html#pyxtal.crystal.random_crystal>`_. Here is a simple example of a carbon crystal:

.. code-block:: Python

  from pyxtal.crystal import random_crystal
  my_crystal = random_crystal(225, ['C'], [3], 1.0)

This would create a crystal structure with space group 225, 3 carbon atoms in the primitive cell, and a volume factor of 1.0. For stoichiometries with more than one type of atom, replace [‘C’] with a list of atomic symbols, and replace [3] with a list of numbers. For example,

.. code-block:: Python

  my_crystal = random_crystal(99, ['Ba','Ti','O'], [1,1,3], 1.0)

would create a random BaTiO3 crystal.

If the generation is successful, the value ``my_crystal.valid`` will be set to True; otherwise, it will be False. The geometric properties of the crystal are stored in ``my_crystal.struct``, which is a `pymatgen.core.structure.Structure <http://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.Structure>`_ object. You can print ``my_crystal.struct`` directly to view the lattice vectors and atomic positions, or you can output the information using

.. code-block:: Python

  my_crystal.to_file(filename)

By default, this will create a cif file storing the structure information. Other file types are supported by specifying the value fmt as ‘poscar’, ‘cssr’, or ‘json’.

3D Molecular Crystals
---------------------

Molecular 3d crystals are generated in the same way as atomic 3d crystals, but atomic species are replaced with (rigid) molecules.

The generating class is `molecular_crystal.molecular_crystal <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal>`_:

.. code-block:: Python
 
  from pyxtal.molecular_crystal import molecular_crystal
  my_crystal = molecular_crystal(36, ['H2O'], [2], 1.0)

This would give a crystal with spacegroup 36, 4 molecules in the conventional cell (2 in the primitive cell), and a volume factor of 1.0. As with atomic crystals, you may use lists as input for the (molecular) stoichiometry.

There are a few other parameters which may be passed to the class. See the `module documentation <pyxtal.molecular_crystal.html>`_ for details. Of particular importance is the variable allow_inversion=False. By default, chiral molecules will not be flipped or inverted while generating the crystal. This is because a chiral molecule’s mirror image may have different chemical properties, especially in a biological setting. But if the mirror images are acceptable for your application, you may use allow_inversion=True, which will allow more spacegroups to be generated. Note that this is only relevant if at least one of the imput molecules is chiral.

Working with Molecules
----------------------

There are 4 options for defining molecules within the molecular_crystal class. You may use a list with any of the following input types:

1) a pre-defined string for the chemical composition (currently supported: "C60", "H2O", "CH4", "NH3", "benzene", "naphthalene", "anthracene", "tetracene", "pentacene", "coumarin", "resorcinol", "benzamide", "aspirin", "ddt", "lindane", "glycine", "glucose", and "ROY"). This will load a molecule from PyXtal's database.

2) a `pymatgen.core.structure.Molecule <http://pymatgen.org/pymatgen.core.structure.html?highlight=class%20molecule#pymatgen.core.structure.Molecule>`_ object.

3) the path to a molecule file (as a string). This will generate a pymatgen Molecule object using the `from_file <http://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IMolecule.from_file>`_ method. Supported formats include .xyz, .gjf, .g03, .g09, .com, .inp, .out, and pymatgen's JSON serialized molecules.

4) a string representing the molecule. This will generate a pymatgen Molecule object using the `from_str <http://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IMolecule.from_str>`_ method. For this option, you must specify the string format (fmt) within the call to molecular_crystal. fmt must be one of: “xyz”, “gjf”, “g03”, or “json”.

For options 3 and 4, installing OpenBabel will allow additional file formats, but is not required.

Because molecules are less symmetric than individual atoms, they may or may not fit within a given Wyckoff position. Furthermore, the molecule may need to be oriented in a certain direction to be compatible with a site. The molecular_crystal class handles this automatically, and only inserts molecules in positions and orientations for which the molecules are sufficiently symmetric. Currently, PyXtal only works with rigid molecules; this simplifies the calculation of symmetry compatibility.

Like atomic crystals, the atomic positions may be accessed with the struct attribute, and stored using to_file(filename). However, for accessing the positions and orientations of the molecules themselves, there is an attribute called mol_generators. This provides a list of `mol_site <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.mol_site>`_ objects, which in turn give the type, location, Wyckoff position, and orientation of each molecule in the asymmetric unit. This can be used to generate the crystal using molecules instead of indivual atoms. Note that the coordinates here are fractional, and refer to the molecule’s center of mass.

The orientations stored in the mol_site class are members of the `operations.orientation <pyxtal.operations.html#pyxtal.operations.orientation>`_ class. A molecule in a Wyckoff position may be allowed to rotate about a certain axis, allowed to rotate freely, or may be rigidly constrained. This information is stored in the orientation class. To obtain a SymmOp which can be applied to the molecule, and which is consistent with the geometric constraints, call `orientation.get_op <pyxtal.operations.html#pyxtal.operations.orientation.get_op>`_. For a 3x3 matrix instead, call `orientation.get_matrix <pyxtal.operations.html#pyxtal.operations.orientation.get_matrix>`_. In either case, this will give a random rotation consistent with the degrees of freedom. To obtain the exact rotation used when generating the crystal (and avoid the random rotation), pass the parameter angle=0.

2D Atomic Crystals
------------------

To generate a 2d crystal, use the class `crystal.random_crystal_2D <pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_. For example,

.. code-block:: Python

  from pyxtal.crystal import random_crystal_2D
  my_crystal = random_crystal_2D(20, ['C'], [4], 2.0, 1.0)

would generate a 2d crystal with layer group 20, 4 carbon atoms in the primitive cell, a thickness of 2.0 Angstroms, and a volume factor of 1.0. As with the 3d case, for crystals with multiple atom types, you may replace [‘C’] and [4] with lists of the atomic symbols and amounts, respectively. The crystal will be periodic in two directions instead of 3. Thus, it is not recommended to store the structure in a .cif or POSCAR file without software specifically designed for handling 2d crystals. The axis of non-periodicity can be accessed via my_crystal.PBC; a value of 1, 2, or 3 corresponds to the x, y, or z axis, respectively.

Note that the layer group number is different from the international space group number, and ranges between 1 and 80. For a list of the layer groups and their symmetry operations, see `the International Tables of Crystallography, Volume E, part 4 <https://it.iucr.org/Eb/ch4o1v0001/contents/>`_.

Note that changing the volume factor will not change the thickness of the output crystal. So, if you are testing over a range of volume factors, consider how the shape of the unit cell will be affected, and change the thickness accordingly. If you are unsure what value to use for the thickness, using None will automatically generate a value for you.

2D Molecular Crystals  
---------------------

2d Molecular crystals are generated using the class `molecular_crystal.molecular_crystal_2D <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_2D>`_:

.. code-block:: Python

  from pyxtal.molecular_crystal import molecular_crystal_2D
  my_crystal = molecular_crystal_2D(20, ['H2O'], [4], 2.0, 1.0)

Here, the parameters correspond to those for `random_crystal_2D <pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_, except the atoms are again replaced with molecules. The additional options available for `molecular_crystal <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal>`_ are also available for `molecular_crystal_2D <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_2D>`_.

Because molecules have a certain thickness of their own, care should be used when choosing a thickness value. Currently, the thickness parameter only determines where the molecular centers of mass can be, so the final crystal may have individual atoms outside of this range.

1D Crystals
-----------

PyXtal also supports generation of 1D crystals using Rod groups (between 1 and 75). The corresponding classes are `crystal.random_crystal_1D
<pyxtal.crystal.html#pyxtal.crystal.random_crystal_1D>`_ and `molecular_crystal_1D
<pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_1D>`_. The parameters for these functions are the same as those for `random_crystal_2D
<pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_ and `molecular_crystal_2D <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_2D>`_. However, in place of the thickness of the unit cell, you should use the cross-sectional area of the unit cell (in Angstroms squared). Again, using an area of None will generate a value automatically.

Point Group Clusters
--------------------

PyXtal also supports generation of atomic clusters with point group symmetry. The corresponding class is `crystal.random_cluster <pyxtal.crystal.html#pyxtal.crystal.random_cluster>`_. As an example, the following code will generate a carbon cluster with 60 atoms and full icosohedral symmetry:

.. code-block:: Python

  from pyxtal.crystal import random_cluster
  my_cluster = molecular_crystal_2D('Ih', ['C'], [60], 1.0)

The parameters are the same as those for `random_crystal
<pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_. The resulting structure can be accessed either via a pymatgen Molecule object (my_cluster.molecule) or via a boxed pymatgen Structure object (my_cluster.struct)

The point group may be specified either by a number (only for the crystallographic point groups), or by a `Schoenflies symbol <https://en.wikipedia.org/wiki/Schoenflies_notation#Point_groups>`_ (ex: 'Ih', 'C*', 'D6h'). The numbers 1-32 correspond to the following crystallograpic point groups: 

+------------+------------+-----------+-----------+
| 1: C1      | 2: Ci      | 3: C2     | 4: Cs     |
+------------+------------+-----------+-----------+
| 5: C2h     | 6: D2      | 7: C2v    | 8: D2h    |
+------------+------------+-----------+-----------+
| 9: C4      | 10: S4     | 11: C4h   | 12: D4    |
+------------+------------+-----------+-----------+
| 13: C4v    | 14: D2d    | 15: D4h   | 16: C3    |
+------------+------------+-----------+-----------+
| 17: C3i    | 18: D3     | 19: C3v   | 20: D3d   |
+------------+------------+-----------+-----------+
| 21: C6     | 22: C3h    | 23: C6h   | 24: D6    |
+------------+------------+-----------+-----------+
| 25: C6v    | 26: D3h    | 27: D6h   | 28: T     |
+------------+------------+-----------+-----------+
| 29: Th     | 30: O      | 31: Td    | 32: Oh    |
+------------+------------+-----------+-----------+

For a list of Wyckoff positions, see the `Bilbao 3D WYCKPOS utility <http://www.cryst.ehu.es/cryst/point_wp.html>`_. The following finite noncrystallographic point groups are also available:

I, Ih, Cn, Cnh, Cnv, Sn, Cni, Dn, Dnh, Dnd.

where n should be replaced by an integer. I and Ih, which are the icosohedral and full icosohedral groups, are particularly useful (Buckminsterfullerene, for example has point group symmetry Ih). Finally, the infinite rotational and dihedral point groups are also available:

C*, C*v, C*h, C*v, D*h

However, only C* and C*h are needed, as the atomic positions will all lie along the z axis. These groups can thus be used for generating linear structures. C*h will have mirror symmetry, while C* will not.

Symmetry Groups
-----------------

The package makes working with symmetry groups simple. Useful information can be accessed directly through the `Group <pyxtal.crystal.html#pyxtal.symmetry.Group>`_ class:

  >>> from pyxtal.symmetry import Group
  >>> g = Group(45)
  >>> g
  -- Space group # 45 --
    8c	site symm: 1
    4b	site symm: ..2
    4a	site symm: ..2

Layer, Rod, and point groups can be accessed by passing the parameter dim=2, dim=1, or dim=0 respectively:

  >>> Group(5, dim=2)
  -- Layer group # 5 --
    2a	site symm: 1
  >>> Group(5, dim=1)
  -- Rod group # 5 --
    2a	site symm: 1
  >>> Group(5, dim=0)
  -- Point group 5 --
    4d	site symm: 1
    2c	site symm: m . .
    2b	site symm: 2 . .
    1a	site symm: 2/m . .

A Group instance contains the Wyckoff positions, site symmetry, and generators for the group. These are stored in the attributes (wyckoffs, w_symm, wyckoff_generators), respectively. Additionally, the Group class stores the lattice type (lattice_type), international number (number), symbol (symbol), and the periodic boundary conditions (PBC). Each group is divided into Wyckoff positions, which are sets of points which possess some subset of the complete group symmetry. Each Wyckoff position in the group has its own `Wyckoff_position <pyxtal.symmetry.html#pyxtal.symmetry.Wyckoff_position>`_ class object, which can be accessed with either a numerical index or the Wyckoff letter:

  >>> g[0]
  Wyckoff position 8c in space group 45 with site symmetry 1
  x, y, z
  -x, -y, z
  x+1/2, -y+1/2, z
  -x+1/2, y+1/2, z
  x+1/2, y+1/2, z+1/2
  -x+1/2, -y+1/2, z+1/2
  x+1, -y+1, z+1/2
  -x+1, y+1, z+1/2
  >>> g['b']
  Wyckoff position 4b in space group 45 with site symmetry ..2
  0, 1/2, z
  1/2, 0, z
  1/2, 1, z+1/2
  1, 1/2, z+1/2

A Wyckoff position is typically denoted with a number-letter combination, depending on its multiplicity. For example, for space group 45 (Iba2), the general Wyckoff position is called “8c”. This is because the position has a multiplicity of 8, and the letters a and b are used by special Wyckoff positions. Note that the naming convention is slightly different for point groups; a point group may have the special Wyckoff position 1o, which corresponds to the point (0,0,0). This is in contrast to the default name 1a.

Each Wyckoff position is further separated into individual operations ('-x,-y,z', '1,1/2,z+1/2', etc.). These are stored as `pymatgen.core.operations.SymmOp <http://pymatgen.org/pymatgen.core.operations.html#pymatgen.core.operations.SymmOp>`_ objects. These symmetry operations can be applied to 3d vectors using op.operate(vector), or can be composed together via multiplication: op3 = op1 * op2. Each SymmOp consists of a rotation matrix (op.rotation_matrix) and a translation vector (op.translation), and is represented by a 4x4 affine matrix (op.affine_matrix).

For a given symmetry group, each Wyckoff position is a subgroup of the general Wyckoff position. As a result, each Wyckoff position requires some point group symmetry for a molecule to occupy it. This symmetry can be accessed using g.w_symm. This returns a nested list, where the first index specifies a Wyckoff position, the second index specifies a point within that Wyckoff position, and the third index specifies a list of symmetry operations corresponding to that point. This list of operations can then be used to check whether a given molecule is consistent with a given Wyckoff position.

As displayed in the example above, the Wyckoff position 4b has site symmetry '..2'. In this example, ‘.’ denotes no symmetry about the x and y axes, and '2' denotes a 2-fold rotation about the z axis. Note that in Hermann-Mauguin notation, the symbols do not always follow this x,y,z format. For more information on reading these symbols, see https://en.wikipedia.org/wiki/Hermann%E2%80%93Mauguin_notation.

Group Settings
--------------------

For the output 3D structures, PyXtal uses the conventional standard cell (the same as `Bilbao
<http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-def-choice>`_). This means unique axis b for monoclinic cells, the obverse triple hexagonal cell for rhombohedral groups, and origin choice 2 (0,0,0) when two origin choices are available.

For 2D structures, we use unique axis c for monoclinic layer groups 3-7, and unique axis a for layer groups 8-18. When two origin choices are available, we use origin choice 1. We always choose c as the non-periodic axis.

For 1D structures, we use unique axis a for monoclinic Rod groups 3-7, and unique axis c for Rod groups 8-12. When two settings are available for a group, we use the 1st setting. We always choose c as the periodic axis.

For point group structures, we use unique axis c for all groups except the polyhedral groups T, Th, O, Td, Oh, I, and Ih. For all of these groups, we place the 2-fold rotation about the z axis and a 3-fold rotation about the (x,x,x) axis. For I and Ih, we use a 5-fold rotation about the axis (1,tau,0), where tau is the golden ratio ( tau = 0.5(sqrt(5)+1) ~ 1.618 ).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
