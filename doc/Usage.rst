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

When a ``Group`` is defined, its Wyckoff position can be accessed
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

Alternatively, the WP can be initialized by itself.

.. code-block:: Python

    >>> from pyxtal.symmetry import Wyckoff_position as wp
    >>> wp.from_group_and_index(19, 0)
    Wyckoff position 4a in space group 19 with site symmetry 1
    x, y, z
    -x+1/2, -y, z+1/2
    -x, y+1/2, -z+1/2
    x+1/2, -y+1/2, -z


In each WP, the symmetry operations are stored as
`SymmOp <http://pymatgen.org/pymatgen.core.operations.html#pymatgen.core.operations.SymmOp>`_ objects.
These symmetry operations can be applied to 3d vectors using ``op.operate``
(vector), or can be composed together via multiplication:

``op3 = op1 * op2``.

Each ``SymmOp`` consists of a rotation matrix (``op.rotation_matrix``)
and a translation vector (``op.translation``), and is represented by a 4x4 affine
matrix (``op.affine_matrix``).

For a given symmetry group, each Wyckoff position is a subgroup of the general
Wyckoff position.

.. only:: comment

    As a result, each Wyckoff position requires some point group
    symmetry for a molecule to occupy it. This symmetry can be accessed using
    ``g.w_symm``. This returns a nested list, where the first index specifies a
    Wyckoff position, the second index specifies a point within that Wyckoff
    position, and the third index specifies a list of symmetry operations
    corresponding to that point. This list of operations can then be used to check
    whether a given molecule is consistent with a given Wyckoff position.


As displayed in the example above, the Wyckoff position ``4b`` has site symmetry
 ``..2``. In this example, ``.`` denotes no symmetry about the x and y axes, and
  ``2`` denotes a 2-fold rotation about the z axis in Hermann-Mauguin notation.
For more information on reading these symbols,
see https://en.wikipedia.org/wiki/Hermann%E2%80%93Mauguin_notation.


pyxtal.molecule.pyxtal_molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are four options for defining molecules within Pyxtal. First, you need to
import the
`pyxtal_molecule <pyxtal.molecule.html#pyxtal.molecule.pyxtal_molecule>`_ class,

.. code-block:: Python

    from pyxtal.molecule import pyxtal_molecule


You may use a list with any of the following input types:

1) From a pre-defined string for the chemical composition
(currently supported: ``C60``, ``H2O``, ``CH4``, ``NH3``, ``benzene``,
``naphthalene``, ``anthracene``, ``tetracene``, ``pentacene``, ``coumarin``,
``resorcinol``, ``benzamide``, ``aspirin``, ``ddt``, ``lindane``, ``glycine``,
``glucose``, and ``ROY``). This will load a molecule from PyXtal's database.

.. code-block:: Python

    mol = pyxtal_molecule('H2O')


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


3) From the path to a molecule file (as a string). This will generate a pymatgen
Molecule object as well. Supported formats include ``.xyz``.

.. code-block:: Python

    mol = pyxtal_molecule('h2o.xyz')


4) a smile string representing the molecule. For example, ``C1=CC=CC=C1.smi``
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


pyxtal.lattice.Lattices
~~~~~~~~~~~~~~~~~~~~~~~

It is possible to supply your own unit cell lattice for a random crystal,
via the `Lattice <pyxtal.crystal.html#pyxtal.crystal.Lattice>`_ class.
You can define a lattice using either a 3x3 matrix, or 6 cell parameters:

.. code-block:: Python

    from pyxtal.lattice import Lattice
    l1 = Lattice.from_matrix([[4.08,0,0],[0,9.13,0],[0,0,5.50]])
    l2 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)

Here, both ``l1`` and ``l2`` describe the same lattice.
In this case, it is an orthorhombic cell with lengths 4.08, 9.13, and 5.50 :math:`\AA`,
which is the unit cell for common water ice. The lattice parameters are,
in order: (a, b, c, :math:`\alpha, \beta, \gamma`).
a, b, and c are the lengths of the lattice vectors;
:math:`\alpha, \beta, \gamma` are the angles (in degrees) between these vectors.


pyxtal.tolerance.Tol_matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When generating random crystals, PyXtal performs inter-atomic distances checks
to make sure the atoms are not too close together. By default, the covalent
radius is used as a basis. However, the user may also define their own criteria
using the `Tol_matrix <pyxtal.crystal.html#pyxtal.crystal.Tol_matrix>`_ class.
To do this, initialize a ``Tol_matrix`` object using one of the built-in methods.

.. code-block:: Python

    from pyxtal.tolerance import Tol_matrix
    tol_m_1 = Tol_matrix(prototype="molecular", factor=2.0)
    tol_m_2 = Tol_matrix.from_radii(some_custom_list_of_atomic_radii)
    tol_m_3 = Tol_matrix.from_matrix(some_custom_2D_tolerance_matrix)

From here, you can alter the tolerance between certain inter-atomic pairs.
Additionally, you can save and reload custom Tol_matrix objects for later use:

.. code-block:: Python

    >>> tol_m_1.set_tol('C', 'N', 2.1)
    >>> tol_m_1.set_tol(1, 3, 4.6)
    >>> tol_m_1.to_file("custom_matrix_file")
    'Output file to custom_matrix_file.npy'
    >>> reloaded_tol_matrix = Tol_matrix.from_file("custom_matrix_file.npy")
    >>> print(reloaded_tol_matrix)
    --Tol_matrix class object--
      Prototype: molecular
      Atomic radius type: covalent
      Radius scaling factor: 2.4
      Custom tolerance values:
        C, N: 2.1
        H, Li: 4.6



Crystal structure generation
----------------------------

First, one can always load an existing crystal from a given file path. Assuming
there is a file.

.. code-block:: Python

    from pyxtal import pyxtal
    my_crystal = pyxtal()
    my_crystal.from_seed(seed=cif_file, style='spglib')


The most important feature of PyXtal is to generate the trial structure
according to customized factors such as space group, cell parameters, partial
occupation. Below we describe the supports on handling different systems from
atomic to molecular, and from 1D to 3D.


Atomic Crystals
~~~~~~~~~~~~~~~

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

    my_crystal = pyxtal()
    my_crystal.from_random(3, 99, ['Ba','Ti','O'], [1,1,3])
    my_crystal
    ------Random Crystal------
    Composition: Ba1 Ti1 O3
    Dimension: 3
    Group: P4mm (99)
    Volume factor: 1.0
    tetragonal lattice:   5.1029   5.1029   4.3018  90.0000  90.0000  90.0000
    Wyckoff sites:
    	Ba @ [0.5000 0.5000 0.3612], Wyckoff letter:  1b, Site symmetry: 4 m m
    	Ti @ [0.5000 0.5000 0.8701], Wyckoff letter:  1b, Site symmetry: 4 m m
    	 O @ [0.5000 0.0000 0.0823], Wyckoff letter:  2c, Site symmetry: 2 mm .
    	 O @ [0.5000 0.5000 0.8177], Wyckoff letter:  1b, Site symmetry: 4 m m

would create a random BaTiO3 crystal. If the generation is successful, the value
of ``my_crystal.valid`` will be set to ``True``;
otherwise, it will be ``False``.


Molecular Crystals
~~~~~~~~~~~~~~~~~~

3D Molecular crystals are generated in the same way as atomic crystals,
but atomic species are replaced with (rigid) molecules.

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


This would give a crystal with space group 36, 4 molecules in the conventional
unit cell.

For molecular crystals, it is possible that a structure is better represented in
a non-standard setting. PyXtal supports the generation of crystals from a
non-standard setting (as defined by the Hall number). Below compares how to
generate the crystals of P21/c and P21/n, which are both in space group 14.

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


There are a few other parameters which may be passed to the class. See the
`documentation <pyxtal.molecular_crystal.html>`_ for details. Of particular
importance is the variable allow_inversion=False. By default, chiral molecules
will not be flipped or inverted while generating the crystal. This is because
a chiral molecule's mirror image may have different chemical properties,
especially in a biological setting. But if the mirror images are acceptable for
your application, you may use allow_inversion=True, which will allow more space
groups to be generated.

The user may also define which orientations are allowed for each molecule
in each Wyckoff position. This is done by setting the orientations parameter.
By default, PyXtal will determine the valid orientations automatically using the
`get_orientations <pyxtal.molecular_crystal.html#molecular_crystal.get_orientations>`
_ function, which in turn calls the
`orientation_in_wyckoff_position <pyxtal.molecule.html#orientation_in_wyckoff_position>`_
function. Setting custom orientations will typically not be necessary, but may
be used to save time during generation; see the source code for more information.


2D Atomic Crystals
~~~~~~~~~~~~~~~~~~

PyXtal can also generate sub-periodic crystals. For example,

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(2, 20, ['C'], [4], thickness=2.0)

would generate a 2d crystal with

- layer group ``P2_122 (20)``,
- 4 carbon atoms in the conventional cell,
- a thickness of 2.0 Angstroms.

The crystal will be periodic in two directions instead of three. PyXtal adds
``10 Angstrom`` of vacuum on each side of the 2D lattice. The axis of
non-periodicity can be accessed via my_crystal.PBC; each axis will either
be 1 or 0, representing either periodicity or non-periodicity. For example,
PBC = [1,1,0] means that the x and y axes are periodic, while the z axis is
non-periodic. Note that the layer group number is different from the space group
number, and ranges between 1 and 80 (see settings). By default, PyXtal will
automatically generate a value for the thickness of the unit cell, based on the
volume. By specifying thickness value, you override this behavior. So, if you
are testing over a range of volume factors, consider how the shape of the unit
cell will be affected, and change the thickness accordingly. Alternatively, you
may supply a custom Lattice object, as described below.

You can generate 1D crystals using Rod groups (between 1 and 75). You should use
the cross-sectional area of the unit cell (in A^2). By default, PyXtal will
automatically generate a value for the area if one is not specified.

.. code-block:: Python

  1d = pyxtal()
  1d.from_random(1, 20, ['C'], [4])

PyXtal also supports generation of atomic clusters with point group symmetry.
As an example, the following code will generate a carbon cluster
with 60 atoms and full icosahedral symmetry:

.. code-block:: Python

  my_cluster = pyxtal()
  my_cluster.from_random(0, 'Ih', ['C'], [60])


The point group may be specified either by a number (only for the crystallographic
point groups), or by a
`Schoenflies <https://en.wikipedia.org/wiki/Schoenflies_notation#Point_groups>`_
symbol (ex: ``Ih``, ``C*``, ``D6h``).



2D and 1D molecular crystals are also supported.

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(2, 20, ['H2O'], [4])
    my_crystal.from_random(1, 20, ['H2O'], [4])


Crystal structure manipulation
------------------------------
To add.


Optional Parameters
-------------------

In addition to the four required parameters, the user can provide additional
constraints.


You can use a custom Lattice to generate a random_crystal or molecular_crystal:

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(3, 36, ['H2O'], [4], lattice=l1)


The Tol_matrix can be passed to a random_crystal object:

.. code-block:: Python

    crystal = pyxtal()
    crystal.from_random(3, 12, ['C','N'], [2,4], tm=tol_m_1)

By default, atomic crystals will use the average of the covalent radii between
two atoms. Molecular crystals will use 1.2 times the sum of the covalent radii
between two atoms. Using ``metallic`` will use the average of the metallic
radius for metals, and the ``covalent`` radius for other atom types.

Note: For monoclinic layer groups, be careful when choosing
the unique axis.

Examples with customized constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, it is convenient to generate the crystal with partial information.
Below shows how to create a Al2SiO5 crystal with a pre-assigned unit cell
and sites on ``8Al+4Si+4O``, and random coordinates on the 16 remaining O atoms.

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




Supports for Different File Formats
-----------------------------------
Once the structures are generated, they can be exported to a variety of formats
for further analysis. PyXtal offers three different formats.

Suppose we generated a carbon structure as follows,

.. code-block:: Python

    from pyxtal import pyxtal
    c = pyxtal()
    c.from_random(3, 225, ['C'], [16])

The `pyxtal` structure object can be conveniently converted to
`Pymatgen` or `ASE Atoms` object.

.. code-block:: Python

    ase_struc = c.to_ase()
    pmg_struc = c.to_pymatgen()


`ASE Atoms` object supports a lot of methods for structural manipulation
and file formats (`cif`, `poscar`, `extxyz`, .etc).

.. code-block:: Python

    ase_struc * 2
    Atoms(symbols='C128', pbc=True, cell=[[13.312249674597792, 0.0, 0.0], [8.151401976723291e-16, 13.312249674597792, 0.0], [8.151401976723291e-16, 8.151401976723291e-16, 13.312249674597792]])

    ase_struc * [1, 2, 2]
    Atoms(symbols='C64', pbc=True, cell=[[6.656124837298896, 0.0, 0.0], [8.151401976723291e-16, 13.312249674597792, 0.0], [8.151401976723291e-16, 8.151401976723291e-16, 13.312249674597792]])

    ase_struc.write('1.vasp', format='vasp', vasp5=True, direct=True)
    ase_struc.write('1.xyz', format='extxyz')


For the molecular crystals, the atomic order will automatically adjusted
when converting when the structure is converted to `ASE Atoms` object.
If you want to keep the original order, just set ``resort=False``
when you call the ``to_ase()`` function.

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(3, 36, ['H2O'], [4], 1.0)
    xtal = my_crystal.to_ase(resort=False)
    print(xtal)

    Atoms(symbols='OH2OH2OH2OH2', pbc=True, cell=[[6.503138824544265, 0.0, 0.0], [3.0183112928813903e-16, 4.929276416649856, 0.0], [3.025303230945897e-16, 3.025303230945897e-16, 4.940695118057273]])

    ordered_xtal = my_crystal.to_ase()
    print(ordered_xtal)
    Atoms(symbols='H8O4', pbc=True, cell=[[6.503138824544265, 0.0, 0.0], [3.0183112928813903e-16, 4.929276416649856, 0.0], [3.025303230945897e-16, 3.025303230945897e-16, 4.940695118057273]])


Molecule in PyXtal
------------------

Random molecular crystal from a customized pyxtal_molecule object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the molecule is defined, one can simply generate the molecular crystal by
the following

.. code-block:: Python

    from pyxtal import pyxtal
    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, [mol], [4])


Random molecular crystal without calling pyxtal_molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you just want to generate a random molecular crystal, Pyxtal will automatically
interpret the strings. Therefore, it is not necessary to call the
``pyxtal_molecule`` class. See a short example below.

.. code-block:: Python

    from pyxtal import pyxtal
    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, ['CC(=O)NC1=CC=CC=C1C(=O)N.smi'], [4])
    print(c1)


Constraints on torsion
~~~~~~~~~~~~~~~~~~~~~~

Using the smile string, one can specify the desired torsions

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
--------------------------------------------

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
sites that can host the H2O molecule.

::

    Wyckoff position 4e in space group 14 with site symmetry 1
    Wyckoff position 8b in space group 36 with site symmetry 1
    Wyckoff position 4a in space group 36 with site symmetry m..
    Wyckoff position 16h in space group 63 with site symmetry 1
    Wyckoff position 8g in space group 63 with site symmetry ..m
    Wyckoff position 8f in space group 63 with site symmetry m..
    Wyckoff position 8e in space group 63 with site symmetry 2..
    Wyckoff position 4c in space group 63 with site symmetry m2m


1D Representation (Experimental)
--------------------------------

For the molecular crystal, PyXtal also provides a
`representation <pyxtal.representation.html#pyxtal.representation.representation.>`_ class to handle the conversion between Pyxtal and its 1D representation. With this module, one can represent the crystal into a 1D array.

.. code-block:: Python

    from pyxtal import pyxtal
    from pyxtal.representation import representation

    c1 = pyxtal(molecular=True)
    print("\n1D string")
    c1.from_seed('pyxtal/database/cifs/aspirin.cif', ['CC(=O)NC1=CC=CC=C1C(=O)N.smi'])

::

    ------Crystal from Seed------
    Dimension: 3
    Composition: [CC(=O)OC1=CC=CC=C1C(=O)O]4
    Group: P21/c (14)
    monoclinic lattice:  11.2330   6.5440  11.2310  90.0000  95.8900  90.0000
    Wyckoff sites:
	H8C9O4 @ [ 0.2252  0.5852  0.0308]  WP:  4e, Site symmetry 1 ==> Euler:   0.00   0.00   0.00

    1D string
     14 0 11.23  6.54 11.23 95.89  0.23  0.59  0.03  130.3   24.9 -147.4   82.9    2.8 -178.3 0

In an 1D string, the data is organized as follows

- space group number (1-230)
- HM sequence
- cell parameter: ``a, b, c, alpha, beta, gamma``
- molecular site: fractional coordinates [``x, y, z``] + orientation [``ang_x,
ang_y, ang_z``] + torsions [``t1, t2, ...``]

Alternatively, one can read the structure from the 1D representation and smile string

.. code-block:: Python

    rep1 = representation(rep.x, ['CC(=O)OC1=CC=CC=C1C(=O)O'])
    xtal = rep1.to_pyxtal()
    print(xtal)


::

    ------Crystal from 1D rep.------
    Dimension: 3
    Composition: [CC(=O)OC1=CC=CC=C1C(=O)O]4
    Group: P21/c (14)
    monoclinic lattice:  11.2330   6.5440  11.2310  90.0000  95.8900  90.0000
    Wyckoff sites:
	H8C9O4 @ [ 0.2252  0.5852  0.0308]  WP:  4e, Site symmetry 1 ==> Euler: 130.31  24.91 -147.41
