PyXtal as a library
===================

While the PyXtal can be used in the command mode, it can become much more powerful with Python scripting. Here we describe the basic functionality of PyXtal as a Python Library.

Random Atomic Crystals
----------------------

PyXtal allows the user to generate random crystal structures with given symmetry constraints. The main class for this is the `crystal.random_crystal <pyxtal.crystal.html#pyxtal.crystal.random_crystal>`_ class. There are several parameters which can be specified, but only four are necessary: 

- the symmetry group, 
- the types of atoms, 
- the number of each atom in the primitive cell
- the volume factor. 
  
Here is a simple example of a 3D carbon crystal:

.. code-block:: Python

    from pyxtal import pyxtal
    my_crystal = pyxtal()
    my_crystal.from_random(3, 225, ['C'], [12])

This would create a crystal structure with 3D structure with space group 225, 12 carbon atoms in the conventional cell. For stoichiometries with more than one type of atom, replace ``[C]`` with a list of atomic symbols, and replace ``[12]`` with a list of numbers. For example,

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

would create a random BaTiO3 crystal.

If the generation is successful, the value ``my_crystal.valid`` will be set to ``True``; otherwise, it will be ``False``. 

2D Crystals
~~~~~~~~~~~

PyXtal can also generate subperiodic crystals. To generate a 2d crystal, use the class `crystal.random_crystal_2D <pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_. For example,

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(2, 20, ['C'], [4], thickness=2.0)

would generate a 2d crystal with 

- layer group ``P2_122 (20)``, 
- 4 carbon atoms in the conventional cell, 
- a thickness of 2.0 Angstroms. 
  
As with the 3d case, for crystals with multiple atom types, you may replace ``[C]`` and ``[4]`` with lists of the atomic symbols and amounts, respectively. The crystal will be periodic in two directions instead of 3. PyXtal adds ``10 Angstroms of vacuum`` on each side of the 2D lattice, so that optimization may be performed without altering the structure file. However, care should be taken when using the cif file for applications designed for 3D crystals. The axis of non-periodicity can be accessed via my_crystal.PBC; each axis will either be 1 or 0, representing either periodicity or non-periodicity. For example, PBC = [1,1,0] means that the x and y axes are periodic, while the z axis is non-periodic.

Note that the layer group number is different from the international space group number, and ranges between 1 and 80. For a list of the layer groups and their symmetry operations, see `the International Tables of Crystallography, Volume E, part 4 <https://it.iucr.org/Eb/ch4o1v0001/contents/>`_ or use the `pyxtal_symmetry utility <COMMAND_MODE.html#pyxtal-symmetry-utility>`_.

By default, PyXtal will automatically generate a value for the thickness of the unit cell, based on the volume. By specifying a value for thickness, you override this behavior. So, if you are testing over a range of volume factors, consider how the shape of the unit cell will be affected, and change the thickness accordingly. Alternatively, you may supply a custom Lattice object, as described below.

1D Crystals
~~~~~~~~~~~

You can generate 1D crystals using Rod groups (between 1 and 75). The corresponding class is `crystal.random_crystal_1D
<pyxtal.crystal.html#pyxtal.crystal.random_crystal_1D>`_. The parameters for this function are the same as those for `random_crystal_2D
<pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_. However, in place of the thickness of the unit cell, you should use the cross-sectional area of the unit cell (in Angstroms squared). Again, by default, PyXtal will automatically generate a value for the area if one is not specified.

Point Group Clusters
~~~~~~~~~~~~~~~~~~~~

PyXtal also supports generation of atomic clusters with point group symmetry. The corresponding class is `crystal.random_cluster <pyxtal.crystal.html#pyxtal.crystal.random_cluster>`_. As an example, the following code will generate a carbon cluster with 60 atoms and full icosohedral symmetry:

.. code-block:: Python

  my_cluster = pyxtal()
  my_cluster.from_random(0, 'Ih', ['C'], [60])


The point group may be specified either by a number (only for the crystallographic point groups), or by a `Schoenflies symbol <https://en.wikipedia.org/wiki/Schoenflies_notation#Point_groups>`_ (ex: ``Ih``, ``C*``, ``D6h``).

One can conveniently access the list of crystallographic point groups via the `Group <pyxtal.symmetry.html#yxtal.symmetry.Group>` class.

.. code-block:: Python

    >>> from pyxtal.symmetry import Group
    >>> g=Group.list_groups(dim=0)
   point_group
    1           C1
    2           Ci
    3           C2
    4           Cs
    5          C2h
    6           D2
    ...
    45         D8h
    46         D4d
    47         D5d
    48         D6d
    49         D7d
    50         D8d
    51          S6
    52          S8
    53         S10
    54         S12
    55           I
    56          Ih
    57          C*
    58         C*h


For a list of Wyckoff positions, see the `Bilbao 3D WYCKPOS utility <http://www.cryst.ehu.es/cryst/point_wp.html>`_. The following finite noncrystallographic point groups are also available:

``I, Ih, Cn, Cnh, Cnv, Sn, Cni, Dn, Dnh, Dnd.``

where n should be replaced by an integer. I and Ih, which are the icosohedral and full icosohedral groups, are particularly useful (Buckminsterfullerene, for example has point group symmetry Ih). Finally, the infinite rotational and dihedral point groups are also available:

``C*, C*v, C*h, D*, D*h``

However, only ``C*`` and ``C*h`` are needed, as the atomic positions will all lie along the z axis. 
These groups can thus be used for generating linear structures. ``C*h`` will have mirror symmetry, while ``C*`` will not.

Random Molecular Crystals
-------------------------

Molecular 3d crystals are generated in the same way as atomic 3d crystals, but atomic species are replaced with (rigid) molecules.

The generating class is `molecular_crystal.molecular_crystal <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal>`_:

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
      
This would give a crystal with spacegroup 36, 4 molecules in the conventional cell. As with atomic crystals, you may use lists as input for the (molecular) stoichiometry.

As with the random_crystal class, the molecular_crystal class has a `print_all <pyxtal.crystal.html#pyxtal.crystal.random_crystal.print_all>`_ function which shows useful information about the structure. In addition to the Wyckoff position and location, you can view the orientation angles for each molecule.


There are a few other parameters which may be passed to the class. See the `module documentation <pyxtal.molecular_crystal.html>`_ for details. Of particular importance is the variable allow_inversion=False. By default, chiral molecules will not be flipped or inverted while generating the crystal. This is because a chiral molecule's mirror image may have different chemical properties, especially in a biological setting. But if the mirror images are acceptable for your application, you may use allow_inversion=True, which will allow more spacegroups to be generated. Note that this is only relevant if at least one of the imput molecules is chiral.

The user may also define which orientations are allowed for each molecule in each Wyckoff position. This is done by setting the orientations parameter. By default, PyXtal will determine the valid orientations automatically using the `get_orientations <pyxtal.molecular_crystal.html#molecular_crystal.get_orientations>`_ function, which in turn calls the `orientation_in_wyckoff_position <pyxtal.molecule.html#orientation_in_wyckoff_position>`_ function. Setting custom orientations will typically not be necessary, but may be used to save time during generation; see the source code for more information.

2D/1D Molecular Crystals  
~~~~~~~~~~~~~~~~~~~~~

2d Molecular crystals are generated using the class `molecular_crystal.molecular_crystal_2D <pyxtal.molecular_crystal.html#pyxtal.molecular_crystal.molecular_crystal_2D>`_:

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(2, 20, ['H2O'], [4])
    my_crystal.from_random(1, 20, ['H2O'], [4])

Optional Parameters
-------------------

In addition to the four required parameters 

- symmetry group, 
- types of atom/molecules,
- number of atoms/molecules, 
- volume factor, 
  
the user can provide additional constraints:


Lattices
~~~~~~~~

It is possible to supply your own unit cell lattice for a random crystal, via the `Lattice <pyxtal.crystal.html#pyxtal.crystal.Lattice>`_ class. You can define a lattice using either a 3x3 matrix, or using the lattice parameters:

.. code-block:: Python

    from pyxtal.lattice import Lattice
    l1 = Lattice.from_matrix([[4.08,0,0],[0,9.13,0],[0,0,5.50]])
    l2 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)

Here, both ``l1`` and ``l2`` describe the same lattice. In this case, it is an orthorhombic lattice with side lengths 4.08, 9.13, and 5.50 Angstrom, which is the unit cell for common water ice. The lattice parameters are, in order: (a, b, c, :math:`\alpha, \beta, \gamma`). a, b, and c are the lengths of the lattice vectors; :math:`\alpha, \beta, \gamma` are the angles (in degrees) between these vectors. You can use a custom Lattice to generate a random_crystal or molecular_crystal:

.. code-block:: Python
 
    my_crystal = pyxtal()
    my_crystal.from_random(3, 36, ['H2O'], [4], lattice=l1)

This would generate a random water ice crystal, with 

- space group 36, 
- 4 molecules in the conventional cell (2 in the primitive cell)
- the lattice which we specified above. 
  
If you do not specify a lattice, a random one will be generated which is consistent with the chosen space group.

Note: For monoclinic layer groups, be careful when choosing the unique axis (see the `Settings <Settings.html>`_ page for details).

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

    crystal = pyxtal()
    crystal.from_random(3, 12, ['C','N'], [2,4], tm=tol_m_1)

By default, atomic crystals will use the average of the covalent radii between two atoms. Molecular crystals will use 1.2 times the sum of the covalent radii between two atoms. Using ``metallic`` will use the average of the metallic radius for metals, and the covalent radius for other atom types.

Supports for Different File Formats
-----------------------------------
Once the structures are generated, they can be exported to a variety of formats for further analysis. PyXtal offers there different mechanisms to manipulate the structureformats.

Suppose we generated a carbon structure as follows,

.. code-block:: Python

    from pyxtal import pyxtal
    c = pyxtal()
    c.from_random(3, 225, ['C'], [16])
    
The `pyxtal` structure object can be conveniently converted to `Pymatgen` or `ASE Atoms` object.

.. code-block:: Python

    ase_struc = c.to_ase()
    pmg_struc = c.to_pymatgen()


`ASE Atoms` object supports a lot of methods for structural manipulation and file formats (`cif`, `poscar`, `extxyz`, .etc).

.. code-block:: Python

    ase_struc * 2
    Atoms(symbols='C128', pbc=True, cell=[[13.312249674597792, 0.0, 0.0], [8.151401976723291e-16, 13.312249674597792, 0.0], [8.151401976723291e-16, 8.151401976723291e-16, 13.312249674597792]])
    
    
    ase_struc * [1, 2, 2]
    Atoms(symbols='C64', pbc=True, cell=[[6.656124837298896, 0.0, 0.0], [8.151401976723291e-16, 13.312249674597792, 0.0], [8.151401976723291e-16, 8.151401976723291e-16, 13.312249674597792]])
    
    
    ase_struc.write('1.vasp', format='vasp', vasp5=True, direct=True)
    ase_struc.write('1.xyz', format='extxyz')
    
    
For the molecular crytals, the atomic order will automatically adjusted when converting when the structure is converted to `ASE Atoms` object. If you want to keep the original order, just set ``resort=False`` when you call the ``to_ase()`` function.

.. code-block:: Python

    my_crystal = pyxtal()
    my_crystal.from_random(3, 36, ['H2O'], [4], 1.0)
    xtal = my_crystal.to_ase(resort=False)
    print(xtal)
    
    Atoms(symbols='OH2OH2OH2OH2', pbc=True, cell=[[6.503138824544265, 0.0, 0.0], [3.0183112928813903e-16, 4.929276416649856, 0.0], [3.025303230945897e-16, 3.025303230945897e-16, 4.940695118057273]])
    
    ordered_xtal = my_crystal.to_ase()
    print(ordered_xtal)
    Atoms(symbols='H8O4', pbc=True, cell=[[6.503138824544265, 0.0, 0.0], [3.0183112928813903e-16, 4.929276416649856, 0.0], [3.025303230945897e-16, 3.025303230945897e-16, 4.940695118057273]])
    
 
    
