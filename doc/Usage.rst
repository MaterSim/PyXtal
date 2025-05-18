PyXtal as a library
===================

While the PyXtal can be used in the command mode, it can become much more
powerful with Python scripting. Here we describe the basic functionality of
PyXtal as a Python Library. This tutorial aims to cover the following contents:

- Built-in PyXtal tools
- Crystal structure generation
- Crystal structure manipulation


Available Tools in PyXtal
-------------------------

PyXtal includes the following functions:

- `Group <pyxtal.symmetry.html#pyxtal.symmetry.Group>`_
- `Wyckoff_position <pyxtal.symmetry.html#pyxtal.symmetry.Wyckoff_position>`_
- `pyxtal_molecule <pyxtal.molecule.html#pyxtal.molecule.pyxtal_molecule>`_
- `Lattice <pyxtal.lattice.html#pyxtal.lattice.Lattice>`_
- `Tol_matrix <pyxtal.tolerance.html#pyxtal.tolerance.Tol_matrix>`_

pyxtal.symmetry.Group
~~~~~~~~~~~~~~~~~~~~~

The Group package makes working with symmetry groups simple. Useful information
can be accessed directly as follows:

.. code-block:: Python

    >>> from pyxtal.symmetry import Group
    >>> g = Group(45)
    >>> g
    -- Space group # 45 --
      8c	site symm: 1
      4b	site symm: ..2
      4a	site symm: ..2
    >>> g.chiral   # check if the space group is enantiomorphic
    False
    >>> g.inversion #check if it has inversion symmetry
    False
    >>> g.polar #check if it is polar
    True


It is important to note that one space group may have multiple settings (see the
`Settings <Settings.html>`_ page for details). To avoid the ambiguity, Hall
introduced the explicit-origin space group notation. Following the Hall notation,
there exist 530 Concise space groups. The full list is available
`online <http://cci.lbl.gov/sginfo/itvb_2001_table_a1427_hall_symbols.html>`_.
In PyXtal, we also the initialization of space group according to Hall number.
Below shows an example to create the Group object of ``Fd-3m (227)``
with the choice 1 of origin.

.. code-block:: Python

    >>> g = Group(525, use_hall=True)
    >>> g.symbol
    'F d -3 m:1'
    >>> g[-1]
    Wyckoff position 8a in space group 227 with site symmetry -4 33 mm
    1/4, 1/4, 1/4
    ...
    1/2, 0, 1/2

For a comparison, we also show ``Fd-3m (227)`` in the standard setting of
with the choice 2 of origin. These two notations only differ in
which symmetry point is placed at (0,0,0).

.. code-block:: Python

    >>> g = Group(526, use_hall=True)
    >>> g.symbol
    'F d -3 m:2'
    >>> g[-1]
    Wyckoff position 8a in space group 227 with site symmetry -4 33 mm
    1/8, 1/8, 1/8
    ...
    3/8, 7/8, 3/8


If one wants to follow the spglib style to initialize the Group object, the
following way should work,

.. code-block:: Python

    >>> g = Group(227, style='spglib')
    >>> g.hall_number
    525


Layer, rod, and point groups can be accessed by passing the parameter ```dim=2``,
``dim=1``, or ``dim=0``, respectively.

.. code-block:: Python

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

A Group instance contains the Wyckoff positions, site symmetry, and generators
for the group. In addition, the Group class stores the lattice type
(``lattice_type``), international number (``number``), symbol (``symbol``),
and the periodic boundary conditions (``PBC``). Each group is divided into
Wyckoff positions, which are sets of points which possess some subset of the
complete group symmetry.


pyxtal.symmetry.Wyckoff_position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Wyckoff position is typically denoted with a number-letter combination,
depending on its multiplicity. For example, for space group ``Iba2 (45)``,
the general Wyckoff position is called ``8c``. This means the position has
a multiplicity of 8. The letters ``a`` and ``b`` are used by special Wyckoff
positions. Note that the name convention is different for point groups; a point
group may have the special Wyckoff position ``1o``, which corresponds to the point
(0,0,0). This is in contrast to the default name ``1a``. Each Wyckoff position
is further separated into individual operations
``('-x,-y,z', '1,1/2,z+1/2', etc.)``.

When a ``Group`` is defined, its ``Wyckoff_position`` can be accessed
with either a numerical index or letter.

.. code-block:: Python

    >>> g[0]
    Wyckoff position 8c in space group 45 with site symmetry 1
    x, y, z
    -x, -y, z
    ...
    x+1, -y+1, z+1/2
    -x+1, y+1, z+1/2
    >>> g['b']
    Wyckoff position 4b in space group 45 with site symmetry ..2
    0, 1/2, z
    ...
    1, 1/2, z+1/2

As displayed in the example above, the Wyckoff position ``4b`` has site symmetry
``..2``. In this example, ``.`` denotes no symmetry about the x and y axes, and
``2`` denotes a 2-fold rotation about the z axis in Hermann-Mauguin notation.
In each WP, the symmetry operations are stored as
`SymmOp <http://pymatgen.org/pymatgen.core.operations.html>`_ objects. These
symmetry operations can be applied to 3d vectors using ``op.operate``, or can be
composed together via multiplication: ``op3 = op1 * op2``. Each ``SymmOp``
consists of a rotation matrix (``op.rotation_matrix``) and a translation vector
(``op.translation_vector``), and is represented by a 4 x 4 affine
matrix (``op.affine_matrix``).

Alternatively, the WP can be initialized by itself.

.. code-block:: Python

    >>> from pyxtal.symmetry import Wyckoff_position as wp
    >>> wp.from_group_and_index(19, 0)
    Wyckoff position 4a in space group 19 with site symmetry 1
    x, y, z
    -x+1/2, -y, z+1/2
    -x, y+1/2, -z+1/2
    x+1/2, -y+1/2, -z




pyxtal.molecule.pyxtal_molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are three options for defining molecules within Pyxtal. First, you need to
import the ``pyxtal_molecule`` class,

.. code-block:: Python

    from pyxtal.molecule import pyxtal_molecule


1) From a pre-defined string for the chemical composition

.. code-block:: Python

    mol = pyxtal_molecule('H2O')

The list of supported molecules are accessible via:

.. code-block:: Python

    >>> pyxtal_molecule.list_molecules()
    >>> ['C60', 'Glycine-z', 'xxvi', 'xxv', 'BIPHEN', 'ANULEN',
    'QUPHEN', 'DBPERY', 'TBZPER', 'TBZPYR', 'YICMOP', 'MERQIM',
    'H2O', 'CH4', 'NH3', 'benzene', 'naphthalene', 'anthracene',
    'tetracene', 'Pentacene', 'coumarin', 'resorcinol', 'benzamide',
    'aspirin', 'ddt', 'lindane', 'Glycine', 'Glucose', 'ROY', 'LEFCIK',
    'OFIXUX', 'HAHCOI', 'JAPWIH', 'WEXBOS', 'LAGNAL', 'LUFHAW',
    'PAHYON01', 'AXOSOW01']


2) From a `Molecule <http://pymatgen.org/pymatgen.core.structure.html>`_ object.

.. code-block:: Python

    from pymatgen.core import Molecule

    xyz="""3
    Water molecule
    O          0.00000        0.00000        0.11779
    H          0.00000        0.75545       -0.47116
    H          0.00000       -0.75545       -0.47116
    """

    m = Molecule.from_str(xyz, fmt='xyz')
    mol = pyxtal_molecule(m)


    # Alternatively, one can load a xyz molecule file.
    # It will be converted to pymatgen.molecule and then passed to pyxtal.
    mol = pyxtal_molecule('h2o.xyz')


3) a smile string representing the molecule. For example, ``C1=CC=CC=C1.smi``
means a benzene molecule. Note that the `.smi` suffix must be included to
indicate that this is a smile string. In this case, **RDKit must be installed
to use this function.**. One can install RDKit by simply typing

``$ conda install -c conda-forge rdkit==2021.09.2``.

Note that the current code is designed for version no later than ``2021.09.2``.

.. code-block:: Python

    mol = pyxtal_molecule('CC(=O)NC1=CC=CC=C1C(=O)N.smi')


After the molecule is defined, its point group will also be parsed:

.. code-block:: Python

    mol = pyxtal_molecule('H2O')
    print(mol.pg)

::

    -- Pointgroup --# 7 (C2v)--
    4d	site symm: 1
    2c	site symm: m . .
    2b	site symm: m . .
    1a	site symm: mm2 . .


pyxtal.lattice.Lattice
~~~~~~~~~~~~~~~~~~~~~~

It is possible to supply your own unit cell lattice for a random crystal,
via the `pyxtal.lattice.Lattice <pyxtal.lattice.html>`_ class. You can define a lattice using
either a 3 x 3 matrix, or 6 cell parameters:

.. code-block:: Python

    from pyxtal.lattice import Lattice
    l1 = Lattice.from_matrix([[4.08,0,0],[0,9.13,0],[0,0,5.50]])
    l2 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)

Here, both ``l1`` and ``l2`` describe the same lattice.
In this case, it is an orthorhombic cell with lengths 4.08, 9.13, and 5.50 Å,
which is the unit cell for common water ice. The lattice parameters are,
in order: (a, b, c, :math:`\alpha, \beta, \gamma`).
a, b, and c are the lengths of the lattice vectors;
:math:`\alpha, \beta, \gamma` are the angles in degrees between these vectors.


pyxtal.tolerance.Tol_matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When generating random crystals, PyXtal performs inter-atomic distances checks
to make sure the atoms are not too close together. By default, the covalent
radius is used as a basis. However, the user may also define their own criteria
using the `pyxtal.tolerance.Tol_matrix <pyxtal.tolerance.html>`_ class.
To do this, initialize a ``Tol_matrix`` object using one of the built-in methods.

.. code-block:: Python

    from pyxtal.tolerance import Tol_matrix
    tol_m_1 = Tol_matrix(prototype="molecular", factor=2.0)
    tol_m_2 = Tol_matrix.from_radii(some_custom_list_of_atomic_radii)
    tol_m_3 = Tol_matrix.from_matrix(some_custom_2D_tolerance_matrix)



Crystal structure generation
----------------------------
PyXtal allows one to generate the crystal from either the existing structure or
from the scratch. First, One can always load an existing crystal from a given
file path. More importantly, PyXtal can generate the trial structure according
to the customized factors such as space group, cell parameters, partial
occupation. It also supports on handling different systems from atomic to
molecular, and from 1D to 3D.

Loading the existing structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming there is a file in either cif or VASP POSCAR format, one can just
load the structure by the ``from_seed`` function.

.. code-block:: Python

    from pyxtal import pyxtal
    my_crystal = pyxtal()

    my_crystal.from_seed(seed=struc_file, style='pyxtal')
    my_crystal.from_seed(seed=struc_file, style='spglib')

Note that the ``style`` flag allows one to represent the structure in different
space group settings.

For the molecular crystal, the molecular information must be provided as a list
(see the molecule section for more details).

.. code-block:: Python

    from pyxtal import pyxtal
    my_crystal = pyxtal(molecular=True)
    my_crystal.from_seed(seed=struc_file, molecule=['aspirin'])

In addition to the existing files in either cif or VASP POSCAR, pyxtal also
provides the interface with Pymatgen and ASE, which support a variety of
structure formats. Below we show a few working examples.

.. code-block:: Python

    from pyxtal import pyxtal
    c = pyxtal()

    # load the structure from ase
    from ase.io import read
    ase_atoms = read('1.cif', format='cif')
    c.from_seed(ase_atoms)

    # load the structure from pymatgen
    from pymatgen.core import Structure
    pmg = read('1.cif', format='cif')
    c.from_seed(pmg)

There is also a shortcut to allow the user to quickly make a well known structure prototype
via the ``from_prototype`` function,

.. code-block:: Python

	>>> from pyxtal import pyxtal
	>>> c = pyxtal()
	>>> c.from_prototype('graphite')
	>>> c

	------Crystal from Build------
	Dimension: 3
	Composition: C4
	Group: P 63/m m c (194)
  	2.4600,   2.4600,   6.7000,  90.0000,  90.0000, 120.0000, hexagonal

	Wyckoff sites:
	C @ [ 0.3333  0.6667  0.2500], WP [2c] Site [-622m2]
	C @ [ 0.0000  0.0000  0.2500], WP [2b] Site [-622m2]


Note that this function so far only supports `graphite`, `diamond`, `a-cristobalite`, `b-cristobalite`, `a-quartz`, `b-quartz`, `rocksalt`, `B1`, `B2`.

A more general approach is to call the ``build`` function,

.. code-block:: Python

	from pyxtal import pyxtal
	from pyxtal.lattice import Lattice

	c = pyxtal()
	l = Lattice.from_para(5.62, 5.62, 5.62, 90, 90, 90, ltype='Cubic')

	# define the sites as dictionary 
	sites = [{"4a": None}, # Na 
         	 {"4b": None}, # Cl
        	]
	c.build(225, ['Na', 'Cl'], [4, 4], lattice=l, sites=sites)
	print(c)

or the ``from_spg_wps_rep`` function

.. code-block:: Python

	from pyxtal import pyxtal
	c = pyxtal()
	c.from_spg_wps_rep(194, ['2c', '2b'], [2.46, 6.70])        # graphite
	c.from_spg_wps_rep(227, ['8a'], [3.6])                     # diamond
	c.from_spg_wps_rep(225, ['4a', '4b'], [5.59], ['Na', 'Cl'] # B1-rocksalt


Random 3D Atomic Crystals
~~~~~~~~~~~~~~~~~~~~~~~~~

PyXtal allows the user to generate random crystal structures with given symmetry
constraints. There are several parameters which can be specified, but only three
are necessary:

- the symmetry group,
- the types of atoms,
- the number of each atom in the primitive cell

Here is a simple example of a 3D carbon crystal:

.. code-block:: Python

    from pyxtal import pyxtal
    my_crystal = pyxtal()
    my_crystal.from_random(3, 225, ['C'], [12])

This would create a crystal structure with 3D structure with space group 225,
12 carbon atoms in the conventional cell. For stoichiometry with more than one
type of atom, replace ``[C]`` with a list of atomic symbols, and replace ``[12]``
with a list of numbers. For example,

.. code-block:: Python

    >>> my_crystal = pyxtal()
    >>> my_crystal.from_random(3, 99, ['Ba','Ti','O'], [1,1,3], random_state=10)
    >>> my_crystal

    ------Crystal from random------
    Dimension: 3
    Composition: Ba1Ti1O3
    Group: P 4 m m (99)
      4.0080,   4.0080,   9.3650,  90.0000,  90.0000,  90.0000, tetragonal
    Wyckoff sites:
    Ba @ [ 0.5000  0.5000  0.7043], WP [1b] Site [4mm]
    Ti @ [ 0.0000  0.0000  0.2199], WP [1a] Site [4mm]
     O @ [ 0.5000  0.0000  0.3652], WP [2c] Site [2mm.]
     O @ [ 0.0000  0.0000  0.6833], WP [1a] Site [4mm]

would create a random ``BaTiO3`` crystal. If the generation is successful, the value of ``my_crystal.valid`` will be set to ``True``; otherwise, it will be ``False``. Keep in mind that the you can pass an integer value to ``random_state`` to ensure the reproducibility. This feature would be useful for many different purposes. 

If you want to generate many random structures with a fixed random number generator, we recommend the following scripts by setting up a random number generator via ``np.random.default_rng``:

.. code-block:: Python

    from pyxtal import pyxtal
    import numpy as np
    rng = np.random.default_rng(0)

    xtals = []
    for i in range(10):
        xtal = pyxtal()
        xtal.from_random(3, 99, ['Ba','Ti','O'], [1,1,3], random_state=rng)
        xtals.append(xtal)

Similar setup can be applied to other modules as well.

Random 3D molecular crystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3D Molecular crystals are generated in the same way as atomic crystals,
but atomic species are replaced with (rigid) molecules. The following script
would give a crystal with space group 36, 4 molecules in the conventional
unit cell.

.. code-block:: Python

    my_crystal = pyxtal(molecular=True)
    my_crystal.from_random(3, 36, ['H2O'], [4])

    ------Random Molecular Crystal------
    Dimension: 3
    Group: Cmc21
    Volume factor: 1.0
    orthorhombic lattice:   5.6448   6.3389   4.4262  90.0000  90.0000  90.0000
    Wyckoff sites:
    H2 O1 @ [ 0.000  0.596  0.986]  Wyckoff letter:  4a, Site symmetry m.. ==> Rotvec: -0.343  0.000  0.000


For molecular crystals, it is possible that a structure is better represented in
a non-standard setting. PyXtal supports the generation of crystals from a
non-standard setting (as defined by the Hall number). Below compares how to
generate the crystals of :math:`P2_1/c` and :math:`P2_1/n`, which are both in
space group 14.

.. code-block:: Python

    >>> from pyxtal import pyxtal
    >>> c1 = pyxtal(molecular=True)
    >>> c1.from_random(3, 81, ["aspirin"], use_hall=True)
    >>> c1
    ------Crystal from random------
    Dimension: 3
    Composition: [aspirin]4
    Group: P 1 21/c 1 (14)
    12.6259,  15.1971,  12.3168,  90.0000,  84.2525,  90.0000, monoclinic
    Wyckoff sites:
	H8C9O4       @ [ 0.6281  0.9928  0.7032]  WP [4e] Site [1] Euler [  57.4  -46.9   89.8]

    >>> c1.from_random(3, 82, ["aspirin"], use_hall=True)
    >>> c1
    ------Crystal from random------
    Dimension: 3
    Composition: [aspirin]4
    Group: P 1 21/n 1 (14)
    16.4395,  16.5499,   9.4357,  90.0000, 113.6587,  90.0000, monoclinic
    Wyckoff sites:
	H8C9O4       @ [ 0.0181  0.6252  0.5789]  WP [4e] Site [1] Euler [-179.0   46.1  -63.9]


While the module is called molecular crystal, one can also use it to general MOF-like structures (e.g., `Na(N3)`) by setting molecular blocks as follows 

.. code-block:: Python

	from pyxtal import pyxtal
	xtal = pyxtal(molecular=True)
	from pyxtal.molecule import pyxtal_molecule
	m1 = pyxtal_molecule('Na+.smi')
	m2 = pyxtal_molecule('[N-]=[N+]=[N-].smi')

	for i, data in enumerate([(14, [4, 4]),
                          	  (36, [4, 4]),
                            ]):
    	(spg, numIons) = data
    	xtal.from_random(3, 14, [m1, m2], [4, 4])
    	xtal.to_file(f"{i}-{spg}.cif")


Random sub-periodic crystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PyXtal can also generate sub-periodic crystals. For example,

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(2, 20, ['C'], [4], thickness=2.0)

would generate a 2d crystal with

- layer group ``P2_122 (20)``,
- 4 carbon atoms in the conventional cell,
- a thickness of 2.0 Å.


The crystal will be periodic in two directions instead of three. PyXtal adds
10 Å of vacuum on the z axis (which is non-periodic). Note that the
layer group number is different from the space group number, and ranges between
1 and 80. By default, PyXtal will automatically generate a value for the
thickness of the unit cell, based on the volume. By specifying thickness value,
you override this behavior. So, if you are testing over a range of volume
factors, consider how the shape of the unit cell will be affected, and change
the thickness accordingly. Alternatively, you
may supply a custom Lattice object.

You can generate 1D crystals using Rod groups (between 1 and 75) and atomic
clusters with point group symmetry.

.. code-block:: Python

  1d = pyxtal()
  1d.from_random(1, 20, ['C'], [4])

  0d= pyxtal()
  0d.from_random(0, 'Ih', ['C'], [60])


The point group may be specified either by a number (only for the crystallographic
point groups), or by a symbol (see the `Settings <Settings.html>`_ page).


2D and 1D molecular crystals are also supported.

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(2, 20, ['H2O'], [4])
    my_crystal.from_random(1, 20, ['H2O'], [4])


Crystal structure Post Analysis
-------------------------------
After the crystal is built, PyXtal allows one to manipulate the structure in
different ways. The following script illustrate some useful functions.

.. code-block:: Python

    # create a random crystal
    c = pyxtal()
    c.from_random(3, 227, ['C'], [8])

    ------Crystal from random------
    Dimension: 3
    Composition: C8
    Group: F d -3 m:2 (227)
    4.9107,   4.9107,   4.9107,  90.0000,  90.0000,  90.0000, cubic
    Wyckoff sites:
	C @ [ 0.1250  0.1250  0.1250], WP [8a] Site [-433mm]

    # get a subgroup representation
    c.subgroup_once(H=141)
    ------Crystal from subgroup------
    Dimension: 3
    Composition: C8
    Group: I 41/a m d:2 (141)
    3.4724,   3.4724,   4.9667,  90.0000,  90.0000,  90.0000, tetragonal
    Wyckoff sites:
	C @ [ 0.0000  0.7500  0.1250], WP [4a] Site [-4mm2]

    # compute the pxrd
    >>> c.get_XRD()
      2theta     d_hkl     hkl       Intensity  Multi
      31.556     2.835   [ 1  1  1]   100.00        8
      52.723     1.736   [ 2  2  0]    42.05       12
      62.755     1.481   [ 3  1  1]    21.09       24
      77.799     1.228   [ 4  0  0]     5.08        6
      86.361     1.127   [ 3  3  1]     7.87       24
     100.543     1.002   [ 4  2  2]    12.92       24
     109.320     0.945   [ 5  1  1]     8.55       24
     125.261     0.868   [ 4  4  0]     7.45       12
     136.483     0.830   [ 5  3  1]    18.32       48
     166.319     0.776   [ 6  2  0]    58.30       24

In addition, the structure can be exported to a variety of formats for further
analysis and process.

.. code-block:: Python

    from pyxtal import pyxtal
    c = pyxtal()
    c.from_random(3, 225, ['C'], [16])

    # export the structure to pymatgen or ase.Atoms object.
    pmg_struc = c.to_pymatgen()
    ase_struc = c.to_ase()

    # ase.Atoms object supports a lot of methods for structural manipulation
    ase_struc *= 2             # create 2*2*2 supercell
    ase_struc *= [1, 2, 2]     # create 1*2*2 supercell

    # Export the structure into different formats
    ase_struc.write('1.vasp', format='vasp', vasp5=True, direct=True)
    ase_struc.write('1.xyz', format='extxyz')


For the molecular crystals, the atomic order will automatically adjusted
when converting when the structure is converted to `ASE Atoms` object.
If you want to keep the original order, just set ``resort=False``
when calling the ``to_ase()`` function.

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(3, 36, ['H2O'], [4], 1.0)
    xtal = my_crystal.to_ase(resort=False)
    print(xtal)

    Atoms(symbols='OH2OH2OH2OH2', pbc=True, cell=[[6.503138824544265, 0.0, 0.0], [3.0183112928813903e-16, 4.929276416649856, 0.0], [3.025303230945897e-16, 3.025303230945897e-16, 4.940695118057273]])

    ordered_xtal = my_crystal.to_ase()
    print(ordered_xtal)
    Atoms(symbols='H8O4', pbc=True, cell=[[6.503138824544265, 0.0, 0.0], [3.0183112928813903e-16, 4.929276416649856, 0.0], [3.025303230945897e-16, 3.025303230945897e-16, 4.940695118057273]])


Subgroup/supergroup manipulation
--------------------------------
Symmetry relation has been playing an important role in crystallography. PyXtal provides several utilities to allow one to conveniently explore the subgroup or supergroup symmetries.

Generating subgroup and supergroup xtals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: Python

    from pyxtal import pyxtal

    # load a graphite crystal
    xtal=pyxtal()
    xtal.from_spg_wps_rep(194, ['2c', '2b'], [2.46, 6.70])

    print("Derive subgroup graphite structures")
    sub_t_xtals = xtal.subgroup(eps=0.01, group_type='t')
    print("t_subgroup xtals", len(sub_t_xtals))
    sub_k_xtals = xtal.subgroup(eps=0.01, group_type='k', max_cell=9)
    print("k_subgroup xtals", len(sub_k_xtals))

Executing this above scripts will lead to the following output:

.. code-block:: Python

    Derive subgroup graphite structures
    t_subgroup xtals 10
    k_subgroup xtals 22


This way, you can easily find derivative crystals in the suboptimal representations. Conversely, it is also possible to identify the likely supergroup xtal. The following snippet codes can be used to design illustrate pyxtal functionalities.

.. code-block:: Python

    from pyxtal import pyxtal

    # load a graphite crystal and make the subgroup representation
    xtal=pyxtal()
    xtal.from_spg_wps_rep(194, ['2c', '2b'], [2.46, 6.70])
    xtal_sub = xtal.subgroup_once(H=164)

    # recheck the symmetry by varying the tolerance values
    for tol in [2e-1, 1e-1, 1e-2]:
        print("Refind the symmetry with tol", tol, xtal_sub.resymmetrize(tol).group.number)

    # make the supergroup symmetry
    xtals, solutions = xtal_sub.supergroup(G=194)
    print(xtals[0])
    print(solutions)


    # Output 
    Refind the symmetry with tol 0.2 194
    Refind the symmetry with tol 0.1 164
    Refind the symmetry with tol 0.01 164

    Recovered xtal with the supergroup symmetry

    ------Crystal from supergroup  0.014------
    Dimension: 3
    Composition: C4
    Group: P 63/m m c (194)
      2.4622,   2.4622,   6.7643,  90.0000,  90.0000, 120.0000, hexagonal
    Wyckoff sites:
    C @ [ 0.0000  0.0000  0.2500], WP [2b] Site [-622m2]
    C @ [ 0.3333  0.6667  0.2500], WP [2c] Site [-622m2]

    The mapping relation between sub/super group symmetries
    [(Wycokff split from 194 to 164

    2b -> 2c
    0, 0, 1/4           -> 0, 0, 1/4                -> 0, 0, z
    0, 0, 3/4           -> 0, 0, 3/4                -> 0, 0, -z

    2c -> 2d
    1/3, -1/3, 1/4      -> 1/3, 2/3, 1/4            -> 1/3, 2/3, z
    2/3, 1/3, 3/4       -> 2/3, 1/3, 3/4            -> 2/3, 1/3, -z
    , [1, 0], array([0.        , 0.        , 0.00206797]), 0, 0.014137159092037652)]



Chemical Substitution
~~~~~~~~~~~~~~~~~~~~~

In many cases, the crystal structures of mutlicompnent systems are strongly related to the structure of simple systems. For instance, the 1:1 ratio boron nitrides, as an isoelectronic analogue to carbon, exihibit very similar structural behaviors as compared to elemental carbon allotropes. Similarly, many of the known AlPO\ :sub:`4` polymorphs are related to SiO\ :sub:`2`. 

Inspired by these known correlation, PyXtal offers the `substitue_1_2 <pyxtal.html#pyxtal.substitue_1_2>`_  function to derive the BC compounds from A via subgroup relation (e.g., from C to BN or from SiO\ :sub:`2` to AlPO\ :sub:`4`). The key idea is to split A's Wyckoff sites to B and C according to the BC composition constraints. Unlike the random substitution, the Wyckoff position splitting strictly follows the group-subgroup relation. As such, the resulting compound retains a high space group symmetry from the parental structure. Below, we illustrate this function via a few examples.

Below is a script to make a 1:1 ratio BN compound from the diamond carbon allotrope.

.. code-block:: Python

    from pyxtal import pyxtal

    # load a diamond crystal
    xtal=pyxtal()
    xtal.from_spg_wps_rep(227, ['8a'], [3.6], ['C'])

    print("Derive symmetry related BN crystal from diamond")
    xtals = xtal.substitute_1_2({'C': ['B', 'N']})


    Derive symmetry related BN crystal from diamond
    Add substitution *  8   1  216 F-43m          3.53 4a 4d
    Found 1 substitutions in total


If you want to generate more BN crystals, you can first generate the subgroup representation and then apply the ``substitute_1_2`` function.

.. code-block:: Python

    print("Derive even more BN crystals from diamond subgroup representations")
    for sub in xtal.subgroup(group_type='t+k'):
        xtals = sub.substitute_1_2({'C': ['B', 'N']})

    Derive even more BN crystals from diamond subgroup representations
    Add substitution *  4   2  119 I-4m2          3.53 2d 2b
    Found 1 substitutions in total
    Add substitution *  4   2  119 I-4m2          3.53 2d 2b
    Found 1 substitutions in total
    Add substitution *  4   2  119 I-4m2          3.53 2d 2b
    Found 1 substitutions in total
    Add substitution *  6   4  160 R3m            3.53 3a 3a
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Found 3 substitutions in total
    Add substitution *  6   4  160 R3m            3.53 3a 3a
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Found 3 substitutions in total
    Add substitution *  6   4  160 R3m            3.53 3a 3a
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Found 3 substitutions in total
    Add substitution *  6   4  160 R3m            3.53 3a 3a
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Add substitution * 12   4  166 R-3m           3.53 6c 6c
    Found 3 substitutions in total
    Add substitution *  8   1  196 F23            3.53 4a 4d
    Found 1 substitutions in total
    Add substitution *  8   1  196 F23            3.53 4b 4d
    Found 1 substitutions in total


Additionally, you may want to generate a more homogeously substituted BN crystal in which each B(N) is 4-coordinated to N(B) from the diamond crystal. In this case, you can set a criteria dictionary to filter the unsatsified structure as follows.


.. code-block:: Python

    print("\nDerive a BN crystal in which each B(N) is 4-coordinated to N(B)")
    criteria = {'CN': {'B': 4, 'N': 4}, 'cutoff': 1.9, 'exclude_ii': True}
    xtals = sub.substitute_1_2({'C': ['B', 'N']}, ratio=[1, 1], criteria=criteria)

    Derive a BN crystal in which each B(N) is 4-coordinated to N(B)
    Add substitution *  8   1  216 F-43m          3.53 4a 4d
    Found 1 substitutions in total

Finally, this function is not limited to elemental-binary substitution. For example, 
below shows a case to convert an α-cristobalite SiO\ :sub:`2` structure to AlPO\ :sub:`4` by 
selectively substituting Si atoms:


.. code-block:: Python

    # load the alpha-cristobalite SiO2
    xtal_sio2 = pyxtal()
    xtal_sio2.from_spg_wps_rep(92, ['4a', '8b'], [5.0847, 7.0986, 0.2944, 0.0941, 0.2410, 0.8256], ['Si', 'O'])
    print(xtal_sio2)

    xtals_alpo4 = xtal_sio2.substitute_1_2({'Si': ['Al', 'P']})
    print(xtals_alpo4)


    ------Crystal from Build------
    Dimension: 3
    Composition: Si4O8
    Group: P 41 21 2 (92)
    5.0847,   5.0847,   7.0986,  90.0000,  90.0000,  90.0000, tetragonal
    Wyckoff sites:
	    Si @ [ 0.2944  0.2944  0.0000], WP [4a] Site [..2]
	     O @ [ 0.0941  0.2410  0.8256], WP [8b] Site [1]

    Add substitution * 24  11   20 C2221          2.21 4b 4a 8c 8c
    Found 1 substitutions in total
    [
    ------Crystal from subgroup------
    Dimension: 3
    Composition: O16P4Al4
    Group: C 2 2 21 (20)
      7.1909,   7.1909,   7.0986,  90.0000,  90.0000,  90.0000, orthorhombic
    Wyckoff sites:
    	Al @ [ 0.0000  0.7089  0.2500], WP [4b] Site [.2.]
    	 P @ [ 0.2077  0.0000  0.0000], WP [4a] Site [2..]
    	 O @ [ 0.9261  0.1658  0.5751], WP [8c] Site [1]
    	 O @ [ 0.8325  0.4278  0.8240], WP [8c] Site [1]]

Advanced examples in random structure generation
-------------------------------------------------

In addition to the required parameters, the user can provide additional
constraints.

Constraints on lattice and sites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, it is convenient to generate the crystal with partial information.
Below shows how to create a Al\ :sub:`2`\ SiO\ :sub:`5`\  crystal with a pre-assigned unit
cell and sites on ``8Al + 4Si + 4O``, and random coordinates on the 16 remaining
O atoms.

.. code-block:: Python

    from pyxtal.lattice import Lattice
    cell = Lattice.from_para(7.8758, 7.9794, 5.6139, 90, 90, 90, ltype='orthorhombic')
    spg = 58
    elements = ['Al', 'Si', 'O']
    composition = [8, 4, 20]

    sites = [{"4e": [0.0000, 0.0000, 0.2418],
              "4g": [0.1294, 0.6392, 0.0000],
             },
             {"4g": [0.2458, 0.2522, 0.0000]},
             {"4g": [0.4241, 0.3636, 0.0000]}, #partial information on O sites
            ]

    s = pyxtal()
    s.from_random(3, spg, elements, composition, lattice=cell, sites=sites)
    print(s)

    ------Crystal from random------
    Dimension: 3
    Composition: O20Si4Al8
    Group: Pnnm (58)
      7.8758,   7.9794,   5.6139,  90.0000,  90.0000,  90.0000, orthorhombic
    Wyckoff sites:
    Al @ [ 0.0000  0.0000  0.2418], WP [4e] Site [..2]
    Al @ [ 0.1294  0.6392  0.0000], WP [4g] Site [..m]
    Si @ [ 0.2458  0.2522  0.0000], WP [4g] Site [..m]
    O @ [ 0.4241  0.3636  0.0000], WP [4g] Site [..m]
    O @ [ 0.5538  0.2648  0.0000], WP [4g] Site [..m]
    O @ [ 0.0000  0.5000  0.6057], WP [4f] Site [..2]
    O @ [ 0.8809  0.5970  0.0786], WP [8h] Site [1]


Similarly, PyXtal allows the user to pre-assign the partial information (e.g.,
lattice, Wyckoff sites) before generating the crystals. A list of scripts is
shown below.

.. code-block:: Python

    s = pyxtal()
    # Generatation with minimum input
    s.from_random(from_random(3, 14, ['aspirin'], [4])

    # Add Lattice constraints
    from pyxtal.lattice import Lattice
    lat = Lattice.from_para(11.43, 6.49, 11.19, 90, 83.31, 90, ltype='monoclinic')
    s.from_random(3, 14, ['aspirin'], [4], lattice=lat)

    # Add sites constraints
    sites = [{"4e": [0.77, 0.57, 0.53]}]
    s.from_random(3, 14, ['aspirin'], [4], lattice=lat, sites=sites)

    # Crystal with 2 water molecules occupying two special wyckoff sites
    # This requires that the molecule is compatible with the site symmetry, be cautious!
    s.from_random(3, 36, ["H2O"], [8], sites=[["4a", "4a"]])


Random molecular crystal without calling pyxtal_molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you just want to generate a random molecular crystal, PyXtal will automatically
interpret the strings. Therefore, it is not necessary to call the
``pyxtal_molecule`` class. See a short example below.

.. code-block:: Python

    from pyxtal import pyxtal
    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, ['CC(=O)NC1=CC=CC=C1C(=O)N.smi'], [4])
    print(c1)


Random molecular crystal with constraints on torsion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the smiles string, one can specify the desired torsions

.. code-block:: Python

    from pyxtal import pyxtal

    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, ['CC(=O)NC1=CC=CC=C1C(=O)N.smi'], [4], torsions=[[-60.2, 1.7, 126.5]])
    print(c1)
    print("Torsions", c1.mol_sites[0].encode()[-4:-1])

::

    ------Crystal from random------
    Dimension: 3
    Composition: [CC(=O)NC1=CC=CC=C1C(=O)N]4
    Group: P21/c (14)
    monoclinic lattice:  19.2246  13.2842  10.1448  90.0000 113.3669  90.0000
    Wyckoff sites:
	    H10C9N2O2 @ [ 0.2497  0.4534  0.9597]  WP:  4e, Site symmetry 1 ==> Euler: -66.31  25.98 -37.99
    Torsions [-60.19971274864328, 1.6999253045986045, 126.50111998425088]



Symmetry Compatibility in Molecular Crystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the molecules with high point group symmetry, it is possible that the
molecule can occupy the special Wyckoff site. Different from other codes,
PyXtal offers an internal function to check if the molecular symmetry is
compatible with the Wyckoff site symmetry. Below is a short example to illustrate
the function.

.. code-block:: Python

    from pyxtal.symmetry import Group
    from pyxtal.molecule import pyxtal_molecule

    mol = pyxtal_molecule('H2O')
    sgs = [14, 36, 63]

    for sg in sgs:
        spg = Group(sg)
        for wp in spg.Wyckoff_positions:
            if len(mol.get_orientations_in_wp(wp)) > 0:
                print(wp.__str__(True))

If you run the above script, it is expected to return all the possible Wyckoff
sites that can host the H₂O molecule.

::

    Wyckoff position 4e in space group 14 with site symmetry 1
    Wyckoff position 8b in space group 36 with site symmetry 1
    Wyckoff position 4a in space group 36 with site symmetry m..
    Wyckoff position 16h in space group 63 with site symmetry 1
    Wyckoff position 8g in space group 63 with site symmetry ..m
    Wyckoff position 8f in space group 63 with site symmetry m..
    Wyckoff position 8e in space group 63 with site symmetry 2..
    Wyckoff position 4c in space group 63 with site symmetry m2m


1D Representation
-----------------

For the molecular crystal, PyXtal also provides a
`representation <pyxtal.representation.html>`_ class to handle the conversion
between Pyxtal and its 1D representation. With this module, one can represent the crystal into a 1D array.

.. code-block:: Python

    from pyxtal import pyxtal

    c1 = pyxtal(molecular=True)
    c1.from_seed('pyxtal/database/cifs/aspirin.cif', ['CC(=O)OC1=CC=CC=C1C(=O)O.smi'])
    rep = c1.get_1D_representation()
    print(rep.to_string())

::

    81 11.23  6.54 11.23  95.9 1 0.23 0.59 0.03   44.1  -25.2   32.5   82.9    2.8 -178.3 1

In the 1D string, the data is organized as follows

- Hall number (1-530)
- cell parameters: *a*, *b*, *c*, *α*, *β*, *γ*
- molecular site: center coordinates + orientation + torsions + inversion

Alternatively, one can read the structure from the 1D representation and smile string

.. code-block:: Python

    from pyxtal.representation import representation
    rep1 = representation(rep.x, ['CC(=O)OC1=CC=CC=C1C(=O)O'])
    xtal = rep1.to_pyxtal()
    print(xtal)


::

    ------Crystal from 1D rep.------
    Dimension: 3
    Composition: [CC(=O)OC1=CC=CC=C1C(=O)O]4
    Group: P 1 21/c 1 (14)
    11.2330,   6.5440,  11.2310,  90.0000,  95.8900,  90.0000, monoclinic
    Wyckoff sites:
	H8C9O4       @ [ 0.2252  0.5852  0.0308]  WP [4e] Site [1] Euler [  44.1  -25.2   32.5]


Database
--------

For molecular crystals, PyXtal provides a
`db <pyxtal.db.html>`_ class to handle store the database with additional information related to the Cambridge Crystallographic Database. **This function requires the access of `CSD Python-api <https://downloads.ccdc.cam.ac.uk/documentation/API/index.html>`_.**

To create a new database file (e.g., `test.db`),

.. code-block:: Python
    
    from pyxtal.db import make_db_from_CSD
    db = make_db_from_CSD('test.db', ['ACSALA', 'BENZEN', 'COUMAR01'])
    print("Initial list of codes", db.codes)
    db.add_from_code('NAPHTA')
    print("Updated list of codes", db.codes)

::

    0 ACSALA
    1 BENZEN
    2 COUMAR01
    Initial list of codes ['ACSALA', 'BENZEN', 'COUMAR01']
    Updated list of codes ['ACSALA', 'BENZEN', 'COUMAR01', 'NAPHTA']


To view the database file, 


``$ ase db test.db``

::

    csd_code|space_group|mol_smi              
    ACSALA  |P21/c      |CC(=O)Oc1ccccc1C(O)=O
    BENZEN  |Pbca       |c1ccccc1             
    COUMAR01|Pca21      |O=C1Oc2ccccc2C=C1    
    NAPHTA  |P21/c      |c1ccc2ccccc2c1       
    Rows: 4

To update some information,

.. code-block:: Python

    from pyxtal.db import database
    db = database('test.db')
    db.add_from_code('XATJOT')
    print("Updated list of codes", db.codes)
    row = db.get_row('XATJOT')
    print("Original smiles", row.mol_smi)
    db.db.update(row.id, mol_smi='[nH+]1cccc2cccnc12.OC(=O)/C=C/C(=O)[O-]')
    row = db.get_row('XATJOT')
    print("Update smiles", row.mol_smi)

::

    Updated list of codes ['ACSALA', 'BENZEN', 'COUMAR01', 'NAPHTA', 'XATJOT']
    Original smiles [nH+]1cccc2cccnc12.OC(=O)/C=C/C(=O)[O-]
    Update smiles [nH+]1cccc2cccnc12.OC(=O)/C=C/C(=O)[O-]


To access the pyxtal structure

.. code-block:: Python

    from pyxtal.db import database
    db = database('test.db')
    xtal = db.get_pyxtal('XATJOT')
    print(xtal)

::

    ------Crystal from Seed------
    Dimension: 3
    Composition: [[nH+]1cccc2cccnc12]4[OC(=O)/C=C/C(=O)[O-]]4
    Group: P c a 21 (29)
    23.5010,   3.7141,  12.6535,  90.0000,  90.0000,  90.0000, orthorhombic
    Wyckoff sites:
	    H7C8N2       @ [ 0.2272  0.3356  0.8232]  WP [4a] Site [1] Euler [   0.0    0.0    0.0]
	    H3C4O4       @ [ 0.5328  0.0993  0.0601]  WP [4a] Site [1] Euler [   0.0    0.0    0.0]


Space Group and Site Symmetry
-----------------------------

PyXtal provides a `site_symmetry <pyxtal.symmetry.html#pyxtal.symmetry.site_symmetry>`_ 
class to handle the conversion of site symmetry symbols and operations.


.. code-block:: Python
    
    from pyxtal import pyxtal
    c = pyxtal()
    c.from_seed('pyxtal/database/cifs/NaSb3F10.cif')
    for site in c.atom_sites:
        print(site)
        ss = site.wp.get_site_symmetry_object()
        ss.to_beautiful_matrix_representation()

::

    Na @ [ 0.3333  0.6667  0.0330], WP [2b] Site [3..]
    Order    Axis    1  -1   2   m   3   4  -4  -3   6  -6   Group
    0 ( 0  0  1):    1   0   0   0   1   0   0   0   0   0     3
    
    Sb @ [ 0.1163  0.3406  0.4500], WP [6c] Site [1]
    Order    Axis    1  -1   2   m   3   4  -4  -3   6  -6   Group
 
    F @ [ 0.9650  0.4560  0.4190], WP [6c] Site [1]
    Order    Axis    1  -1   2   m   3   4  -4  -3   6  -6   Group
 
    F @ [ 0.7960  0.1890  0.7060], WP [6c] Site [1]
    Order    Axis    1  -1   2   m   3   4  -4  -3   6  -6   Group
 
    F @ [ 0.8890  0.1180  0.3600], WP [6c] Site [1]
    Order    Axis    1  -1   2   m   3   4  -4  -3   6  -6   Group
 
    F @ [ 0.3333  0.6667  0.4550], WP [2b] Site [3..]
    Order    Axis    1  -1   2   m   3   4  -4  -3   6  -6   Group
    0 ( 0  0  1):   1   0   0   0   1   0   0   0   0   0     3


One can also access the matrix representation via the `to_matrix_representation <pyxtal.symmetry.html#pyxtal.symmetry.site_symmetry>`_ method.

.. code-block:: Python

    matrix = ss.to_matrix_representation()

This will result in a `15×10` array to represent the presence of 10 fundamental symmetry elements in 15 possible high symmetry crystallograph axes.
::
    
    # An example of 3-fold rotation symmetry on the (0 0 1) axis
    array([[1, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 1, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0]])   


The one-hot matrix representation `(15, 13)` can also be obtained via 

.. code-block:: Python

    one_hot = ss.to_one_hot()

::

    [[1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [0 0 0 0 1 0 0 0 0 0 0 0 0] # 3
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [0 0 1 0 0 0 0 0 0 0 0 0 0] # 2
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [0 0 1 0 0 0 0 0 0 0 0 0 0] # 2
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [0 0 1 0 0 0 0 0 0 0 0 0 0] # 2
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
     [1 0 0 0 0 0 0 0 0 0 0 0 0] # 1
    ]


Finally, the both ``Group`` and ``Wyckoff_position`` classes support the access of the symmetry element analysis via

.. code-block:: Python

    from pyxtal.symmetry import Group
    g = Group(14)
    ss_spg = g.get_spg_symmetry_object()
    print('Space group symmetry table', g.number, g.symbol)
    ss_spg.to_beautiful_matrix_representation()

    for wp in g:
        print('\nWp symmetry table', wp.get_label())
        ss_wp = wp.get_site_symmetry_object()
        ss_wp.to_beautiful_matrix_representation()

::

    Space group symmetry table 14 P21/c
    Order    Axis   1   -1  2   2_1  m   a   b   c   n   d   3  3_1  3_2  4  -4   
    0 ( 0  1  0):   1   1   0   1    0   0   0   1   0   0   0  0    0    0  0 
    
    Wp symmetry table 4e
    Order    Axis   1   -1  2   m    3   4   -4   Group
    0 ( 0  1  0):   1   0   0   0    0   0   0    1
    
    Wp symmetry table 2d
    Order    Axis   1   -1  2   m    3   4   -4   Group
    0 ( 0  1  0):   1   1   0   0    0   0   0    -1
    
    Wp symmetry table 2c
    Order    Axis   1   -1  2   m    3   4   -4   Group
    0 ( 0  1  0):   1   1   0   0    0   0   0    -1
    
    Wp symmetry table 2b
    Order    Axis   1   -1  2   m    3   4   -4   Group
    0 ( 0  1  0):   1   1   0   0    0   0   0    -1
    
    Wp symmetry table 2a
    Order    Axis   1   -1  2   m    3   4   -4   Group
    0 ( 0  1  0):   1   1   0   0    0   0   0    -1
    

More details on the symmetry table can be found in the `documentation <Wyckoff_positions.html>`_.
