Other features
===================

Molecule in pyxtal
------------------

There are 4 options for defining molecules within pyxtal. First, you need to import the `pyxtal_molecule <https://pyxtal.readthedocs.io/en/latest/pyxtal.molecule.html#pyxtal.molecule.pyxtal_molecule>`_ class,

.. code-block:: Python

    from pyxtal.molecule import pyxtal_molecule


You may use a list with any of the following input types:

1) From a pre-defined string for the chemical composition (currently supported: ``C60``, ``H2O``, ``CH4``, ``NH3``, ``benzene``, ``naphthalene``, ``anthracene``, ``tetracene``, ``pentacene``, ``coumarin``, ``resorcinol``, ``benzamide``, ``aspirin``, ``ddt``, ``lindane``, ``glycine``, ``glucose``, and ``ROY``). This will load a molecule from PyXtal's database.

.. code-block:: Python

    mol = pyxtal_molecule('H2O')


2) From a `pymatgen.core.structure.Molecule <http://pymatgen.org/pymatgen.core.structure.html?highlight=class%20molecule#pymatgen.core.structure.Molecule>`_ object.

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



3) From the path to a molecule file (as a string). This will generate a pymatgen Molecule object using the `from_file <http://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IMolecule.from_file>`_ method. Supported formats include ``.xyz``.

.. code-block:: Python

    mol = pyxtal_molecule('h2o.xyz')


4) a smile string representing the molecule. For example, 'C1=CC=CC=C1.smi' means a benzene molecule. Note that the `.smi` suffix must be included to indicate that this is a smile string. In this case, **RDKit must be installed to use this function.**. One can install RDKit by simply typing ``$ conda install -c conda-forge rdkit==2021.03.5``. Note that the current code is designed for version no later than ``2021.03.5``.

.. code-block:: Python

    mol = pyxtal_molecule('CC(=O)NC1=CC=CC=C1C(=O)N.smi')
	
	
After the molecule is defined, its point group will also be parsed, one can access this information by:

.. code-block:: Python
  
    mol = pyxtal_molecule('H2O')
    print(mol.pg)

::
    
    -- Pointgroup --# 7 (C2v)--
    4d	site symm: 1
    2c	site symm: m . .
    2b	site symm: m . .
    1a	site symm: mm2 . .


Generate the molecular crystal from a customized pyxtal_molecule object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the molecule is defined, one can simply generate the molecular crystal by the following

.. code-block:: Python

    from pyxtal import pyxtal
    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, [mol], [4])
    

Generate the molecular crystal from the user
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you just want to generate a random molecular crystal, pyxtal will automatically interpret the strings. Therefore, it is not necessary to call the ``pyxtal_molecule`` class. See a short example below.

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
    
	
Symmetry Compatibility
~~~~~~~~~~~~~~~~~~~~~~

For the molecules with high point group symmetry, it is possible that the molecule can occupy the special Wyckoff site. Different from other codes, PyXtal offers an internal function to check if the molecular symmetry is compatible with the Wyckoff site symmetry. Below is a short example to illustrate the function.

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
			
If you run the above script, it is expected to return all the possible Wyckoff sites that can host the H2O molecule.
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

For the molecular crystal, PyXtal also provides a `representation <pyxtal.representation.html#pyxtal.representation.representation.>`_ class to handle the conversion between Pyxtal and its 1D representation. With this module, one can represent the crystal into a 1D array.
    
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
- HM sequence (for monoclinic system like space group 14, 0 is ``P21/c``, 1 is ``P21/n``)
- cell parameter: ``a, b, c, alpha, beta, gamma`` (For othorhombic system, only a, b, c is specified)
- molecular site: fractional coordinates [``x, y, z``] + orientation [``ang_x, ang_y, ang_z``] + torsions [``t1, t2, ...``]

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


Symmetry Groups and Wyckoff Positions
-------------------------------------

The package makes working with symmetry groups simple. Useful information can be accessed directly through the 
`Group <pyxtal.symmetry.html#pyxtal.symmetry.Group>`_ class:

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
