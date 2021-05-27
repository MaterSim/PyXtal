Other features
===================

Working with Molecules
----------------------

There are 4 options for defining molecules within the molecular_crystal class. You may use a list with any of the following input types:

1) a pre-defined string for the chemical composition (currently supported: ``C60``, ``H2O``, ``CH4``, ``NH3``, ``benzene``, ``naphthalene``, ``anthracene``, ``tetracene``, ``pentacene``, ``coumarin``, ``resorcinol``, ``benzamide``, ``aspirin``, ``ddt``, ``lindane``, ``glycine``, ``glucose``, and ``ROY``). This will load a molecule from PyXtal's database.

2) a `pymatgen.core.structure.Molecule <http://pymatgen.org/pymatgen.core.structure.html?highlight=class%20molecule#pymatgen.core.structure.Molecule>`_ object.

3) the path to a molecule file (as a string). This will generate a pymatgen Molecule object using the `from_file <http://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IMolecule.from_file>`_ method. Supported formats include ``.xyz``.

4) a smile string representing the molecule. For example, 'C1=CC=CC=C1.smi' means a benzene molecule. Note that the `.smi` suffix must be included to indicate that this is a smile string. 

**RDKit must be installed to use this function.** One can install RDKit by simply typing 

::

    $ conda install -c rdkit rdkit
    
Below is an example to generate a random molecular crystal from smile string.

.. code-block:: Python

    from pyxtal import pyxtal
    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, ['CC(=O)NC1=CC=CC=C1C(=O)N.smi'], [4])
    print(c1)
    

One advantage of using the smile string is that one can also specify the desired torsions

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
    angles [-60.19971274864328, 1.6999253045986045, 126.50111998425088]


Because molecules are less symmetric than individual atoms, they may or may not fit within a given Wyckoff position. Furthermore, the molecule may need to be oriented in a certain direction to be compatible with a site. The `molecular_crystal class <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal>`_ handles this automatically, and only inserts molecules in positions and orientations for which the molecules are sufficiently symmetric. Currently, PyXtal only works with rigid molecules; this simplifies the calculation of symmetry compatibility.

Like atomic crystals, the atomic positions may be accessed with the struct attribute, and stored using to_file(filename). However, for accessing the positions and orientations of the molecules themselves, there is an attribute called mol_generators. This provides a list of `mol_site <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.mol_site>`_ objects, which in turn give the type, location, Wyckoff position, and orientation of each molecule in the asymmetric unit. This can be used to generate the crystal using molecules instead of indivual atoms. Note that the coordinates here are fractional, and refer to the molecule’s center of mass.

The orientations stored in the mol_site class are members of the `operations.orientation <pyxtal.operations.html#pyxtal.operations.orientation>`_ class. A molecule in a Wyckoff position may be allowed to rotate about a certain axis, allowed to rotate freely, or may be rigidly constrained. This information is stored in the orientation class. To obtain a SymmOp which can be applied to the molecule, and which is consistent with the geometric constraints, call `orientation.get_op <pyxtal.operations.html#pyxtal.operations.orientation.get_op>`_. For a 3x3 matrix instead, call `orientation.get_matrix <pyxtal.operations.html#pyxtal.operations.orientation.get_matrix>`_. In either case, this will give a random rotation consistent with the degrees of freedom. To obtain the exact rotation used when generating the crystal (and avoid the random rotation), pass the parameter angle=0.

Symmetry Groups and Wyckoff Positions
-------------------------------------

The package makes working with symmetry groups simple. Useful information can be accessed directly through the 
`Group <pyxtal.symmetry.html#yxtal.symmetry.Group>`_ class:

.. code-block:: Python

    >>> from pyxtal.symmetry import Group
    >>> g = Group(45)
    >>> g
    -- Space group # 45 --
      8c	site symm: 1
      4b	site symm: ..2
      4a	site symm: ..2

Layer, Rod, and point groups can be accessed by passing the parameter dim=2, dim=1, or dim=0 respectively:

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

A Group instance contains the Wyckoff positions, site symmetry, and generators for the group. These are stored in the attributes (``wyckoffs``, ``w_symm``, ``wyckoff_generators``), respectively. Additionally, the Group class stores the lattice type (``lattice_type``), international number (``number``), symbol (``symbol``), and the periodic boundary conditions (``PBC``). Each group is divided into Wyckoff positions, which are sets of points which possess some subset of the complete group symmetry. Each Wyckoff position in the group has its own `Wyckoff_position <pyxtal.symmetry.html#pyxtal.symmetry.Wyckoff_position>`_ class object, which can be accessed with either a numerical index or the Wyckoff letter:

.. code-block:: Python

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

A Wyckoff position is typically denoted with a number-letter combination, depending on its multiplicity. For example, for space group ``Iba2 (45)``, the general Wyckoff position is called ``8c``. This is because the position has a multiplicity of 8, and the letters a and b are used by special Wyckoff positions. Note that the naming convention is slightly different for point groups; a point group may have the special Wyckoff position 1o, which corresponds to the point (0,0,0). This is in contrast to the default name ``1a``.

Each Wyckoff position is further separated into individual operations ``('-x,-y,z', '1,1/2,z+1/2', etc.)``. These are stored as `pymatgen.core.operations.SymmOp <http://pymatgen.org/pymatgen.core.operations.html#pymatgen.core.operations.SymmOp>`_ objects. These symmetry operations can be applied to 3d vectors using ``op.operate`` (vector), or can be composed together via multiplication: ``op3 = op1 * op2``. Each ``SymmOp`` consists of a rotation matrix (``op.rotation_matrix``) and a translation vector (``op.translation``), and is represented by a 4x4 affine matrix (``op.affine_matrix``).

For a given symmetry group, each Wyckoff position is a subgroup of the general Wyckoff position. As a result, each Wyckoff position requires some point group symmetry for a molecule to occupy it. This symmetry can be accessed using ``g.w_symm``. This returns a nested list, where the first index specifies a Wyckoff position, the second index specifies a point within that Wyckoff position, and the third index specifies a list of symmetry operations corresponding to that point. This list of operations can then be used to check whether a given molecule is consistent with a given Wyckoff position.

As displayed in the example above, the Wyckoff position ``4b`` has site symmetry ``..2``. In this example, ``.`` denotes no symmetry about the x and y axes, and ``2`` denotes a 2-fold rotation about the z axis. Note that in Hermann-Mauguin notation, the symbols do not always follow this x,y,z format. For more information on reading these symbols, see https://en.wikipedia.org/wiki/Hermann%E2%80%93Mauguin_notation.
