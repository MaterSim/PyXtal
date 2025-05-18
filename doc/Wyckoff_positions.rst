Symmetry Representation
=======================

Space Group Symmetry
--------------------
In crystallography, the short Hermann-Mauguin notation is used to represent the symmetry of a crystal structure. Using this notation, the symmetry of a crystal structure can be described by a space group symbol, which is a combination of letters and numbers that represent the symmetry operations that can be applied to the structure. A full Hermann-Mauguin notation can be broken down into four parts:

1. A single letter (``P, A, B, C, I, R``) indicating the type of lattice
2. A symbol denoting the symmetry for the primary direction.
3. A symbol denoting the symmetry for the secondary direction.
4. A symbol denoting the symmetry for the tertiary direction. 

.. list-table:: Symmetry definitions in the Hermann-Mauguin notation.
  :header-rows: 1
  :widths: 7 7 7 7 30
  :align: left 

  * - System
    - Primary  
    - Secondary
    - Tertiary
    - Possible Symmetry Elements
  * - Triclinic
    - —
    - —
    - —
    - 1, -1.
  * - Monoclinic
    - [010]
    - —
    - —
    - 2, 2₁, m, c, 2/m, 2/c, 2₁/c.
  * - Orthorhombic
    - [100]
    - [010]
    - [001]
    - 2, 2₁, m, a, b, c, n, d.
  * - Tetragonal
    - [001]
    - [100]/[010]
    - [110]/[101]/...
    - 4ₙ(/m) @ [001].
  * - Tri/Hexagonal
    - [001]
    - [100]/[010]
    - [110]/[101]/..
    - 3ₙ/6ₙ, -3, -6 @ [001]
  * - Cubic
    - [100]/...
    - [111]/[1̄11]/...
    - [110]/[101]/...
    - 4ₙ @ [100], 3ₙ @ [111], 2ₙ @ [110]


For example, the space group 227 has a symbol of ``Fd-3m`` or ``F 4₁/d -3 2/m`` in the full Hermann-Mauguin notation. This means it has 

1. a face centered unit cell
2. a ``4₁/d`` (including ``1, -1, 2₁, d, 4₁, -4``) symmetry in the [100] family direction.
3. a ``-3`` (including ``1, -1, 3``) symmetry in the [111] family directions.
4. a ``2/m`` (including ``1, -1, 2, m``) symmetry in the [110] family directions.

Below shows the symmetry element analysis of the space group 227 (Fd-3m) on the 15 high symmetry crystallographic axes.

.. list-table:: Symmetry analysis of space group 227 (Fd-3m)
    :widths: 10 20 30 20
    :header-rows: 1

    * - Index 
      - Direction
      - Symmetry elements
      - Symbol  
    * - 0
      - (1, 0, 0)
      - 1, -1, 2₁, d, 4₁, -4
      - 4₁/d
    * - 1
      - (0, 1, 0)
      - 1, -1, 2₁, d, 4₁, -4
      - 4₁/d
    * - 2
      - (0, 0, 1)
      - 1, -1, 2₁, d, 4₁, -4
      - 4₁/d
    * - 3
      - (1, 1, 1)
      - 1, -1, 3
      - -3
    * - 4
      - (1, -1, -1)
      - 1, -1, 3
      - -3
    * - 5
      - (-1, 1, -1)
      - 1, -1, 3
      - -3
    * - 6
      - (-1, -1, 1)
      - 1, -1, 3
      - -3
    * - 7
      - (1, -1, 0)
      - 1, -1, 2, m
      - 2/m
    * - 8
      - (1, 1, 0)
      - 1, -1, 2, d
      - 2/d
    * - 9
      - (0, 1, -1)
      - 1, -1, 2, m
      - 2/m
    * - 10
      - (0, 1, 1)
      - 1, -1, 2, d
      - 2/d
    * - 11
      - (-1, 0, 1)
      - 1, -1, 2, m
      - 2/m
    * - 12
      - (1, 0, 1)
      - 1, -1, 2, d
      - 2/d
    * - 13
      - (1, -2, 0)
      - 1, -1
      - -1
    * - 14
      - (2, -1, 0)
      - 1, -1
      - -1


One can easily verify this via PyXtal as follows.

.. code-block:: python

    from pyxtal.symmetry import Group
    spg = Group(227)
    ss = spg.get_spg_symmetry_object()
    ss.to_beautiful_matrix_representation(skip=False)

>>> Order Axis   1  -1 2  2_1 m  a  b  c  n  d  3  3_1 3_2 4  4_1 4_2 4_3 -4   
  0 ( 1  0  0):  1  1  0  1   0  0  0  0  0  1  0  0   0   0  1   0   0   1  # 4_1/d
  0 ( 0  1  0):  1  1  0  1   0  0  0  0  0  1  0  0   0   0  1   0   0   1  # 4_1/d
  0 ( 0  0  1):  1  1  0  1   0  0  0  0  0  1  0  0   0   0  1   0   0   1  # 4_1/d
  1 ( 1  1  1):  1  1  0  0   0  0  0  0  0  0  1  0   0   0  0   0   0   0  # -3
  1 ( 1 -1 -1):  1  1  0  0   0  0  0  0  0  0  1  0   0   0  0   0   0   0  # -3
  1 (-1  1 -1):  1  1  0  0   0  0  0  0  0  0  1  0   0   0  0   0   0   0  # -3
  1 (-1 -1  1):  1  1  0  0   0  0  0  0  0  0  1  0   0   0  0   0   0   0  # -3
  2 ( 1 -1  0):  1  1  1  0   1  0  0  0  0  0  0  0   0   0  0   0   0   0  # 2/m
  2 ( 1  1  0):  1  1  1  0   0  0  0  0  0  1  0  0   0   0  0   0   0   0  # 2/m 
  2 ( 0  1 -1):  1  1  1  0   1  0  0  0  0  0  0  0   0   0  0   0   0   0  # 2/m
  2 ( 0  1  1):  1  1  1  0   0  0  0  0  0  1  0  0   0   0  0   0   0   0  # 2/m
  2 (-1  0  1):  1  1  1  0   1  0  0  0  0  0  0  0   0   0  0   0   0   0  # 2/m
  2 ( 1  0  1):  1  1  1  0   0  0  0  0  0  1  0  0   0   0  0   0   0   0  # 2/m

In the above output, we use a 1D 18-length array to represent the symmetry elements (``1, -1, 2, 2₁, m, a, b, c, n, d, 3, 3₁, 3₂, 4, 4₁, 4₂, 4₃, -4``). 

Space Group Symmetry's Digital Representation
---------------------------------------------

For the application of deep learning, the space group symmetry can be represented as a **14 lattice indices** + **15x18 matrix**. 

For the lattice index, we divde it into 14 possible cases (0-13) as shown below:

.. list-table::
    :widths: 15 35
    :header-rows: 1

    * - ID
      - Lattice Type
    * - 0 
      - Triclinic Primitive
    * - 1
      - Monoclinic Primitive 
    * - 2
      - Monoclinic Base Centered
    * - 3
      - Orthorhombic Primitive
    * - 4
      - Orthorhombic Base Centered
    * - 5
      - Orthorhombic Body Centered
    * - 6
      - Orthorhombic Face Centered
    * - 7
      - Tetragonal Primitive
    * - 8
      - Tetragonal Body Centered
    * - 9
      - Hexagonal Primitive
    * - 10
      - Hexagonal Rhombehedral 
    * - 11
      - Cubic Primitive
    * - 12
      - Cubic Body Centered
    * - 13
      - Cubic Face Centered

For the matrix, the rows represent 15 directions, including

    .. list-table::
        :widths: 15 35
        :header-rows: 1

        * - Index
          - Direction
        * - 0
          - (1, 0, 0)
        * - 1
          - (0, 1, 0)
        * - 2
          - (0, 0, 1)
        * - 3
          - (1, 1, 1)
        * - 4
          - (1, -1, -1)
        * - 5
          - (-1, 1, -1)
        * - 6
          - (-1, -1, 1)
        * - 7
          - (1, -1, 0)
        * - 8
          - (1, 1, 0)
        * - 9
          - (0, 1, -1)
        * - 10
          - (0, 1, 1)
        * - 11
          - (-1, 0, 1)
        * - 12
          - (1, 0, 1)
        * - 13
          - (1, 2, 0)
        * - 14
          - (2, 1, 0)


Ant the columns represent 18 symmetry elements. The value of each element in the matrix indicates the presence (1) or absence (0) of a symmetry operation. There exist a total of 48 combinational symmetries, 

.. list-table:: 
    :widths: 10 20 10
    :header-rows: 1

    * - Index
      - Symmetry elements  
      - Symbol
    * - 0
      - 1
      - 1
    * - 1  
      - 1, -1
      - ̄1
    * - 2
      - 1, 2
      - 2
    * - 3
      - 1, 2₁
      - 2₁
    * - 4
      - 1, m
      - m
    * - 5
      - 1, a
      - a
    * - 6
      - 1, b
      - b
    * - 7
      - 1, c
      - c
    * - 8
      - 1, n
      - n
    * - 9
      - 1, d
      - d
    * - 10
      - 1, 3
      - 3
    * - 11
      - 1, 3₁
      - 3₁
    * - 12
      - 1, 3₂
      - 3₂
    * - 13
      - 1, -1, 2, m
      - 2/m
    * - 14
      - 1, -1, 2, a
      - 2/a
    * - 15
      - 1, -1, 2, b
      - 2/b
    * - 16
      - 1, -1, 2, c
      - 2/c
    * - 17
      - 1, -1, 2, n
      - 2/n
    * - 18
      - 1, -1, 2, d
      - 2/d
    * - 19
      - 1, -1, 2₁, m
      - 2₁/m
    * - 20
      - 1, -1, 2₁, a
      - 2₁/a
    * - 21
      - 1, -1, 2₁, b
      - 2₁/b
    * - 22
      - 1, -1, 2₁, c
      - 2₁/c
    * - 23
      - 1, -1, 2₁, n
      - 2₁/n
    * - 24
      - 1, -1, 2₁, d
      - 2₁/d
    * - 25
      - 1, 2, 4
      - 4
    * - 26
      - 1, 2₁, 4₁
      - 4₁
    * - 27
      - 1, 2, 4₂
      - 4₂
    * - 28
      - 1, 2₁, 4₃
      - 4₃
    * - 29
      - 1, 2, -4
      - -4
    * - 30
      - 1, -1, 3
      - -3
    * - 31
      - 1, 2, 3
      - 6
    * - 32
      - 1, 2₁, 3₁
      - 6₁
    * - 33
      - 1, 2₁, 3₂
      - 6₅
    * - 34
      - 1, 2, 3₂
      - 6₂
    * - 35
      - 1, 2, 3₁
      - 6₄
    * - 36
      - 1, 2₁, 3
      - 6₃
    * - 37
      - 1, m, 3
      - -6
    * - 38
      - 1, -1, 2, m, 4, -4
      - 4/m
    * - 39
      - 1, -1, 2, n, 4, -4
      - 4/n
    * - 40
      - 1, -1, 2₁, a, 4₁, -4
      - 4₁/a
    * - 41
      - 1, -1, 2₁, b, 4₁, -4
      - 4₁/b
    * - 42
      - 1, -1, 2₁, c, 4₁, -4
      - 4₁/c
    * - 43
      - 1, -1, 2₁, d, 4₁, -4
      - 4₁/d
    * - 44
      - 1, -1, 2, m, 4₂, -4
      - 4₂/m
    * - 45
      - 1, -1, 2, m, 4₂, -4
      - 4₂/n
    * - 46
      - 1, -1, 2, m, 3
      - 6/m
    * - 47
      - 1, -1, 2₁, m, 3
      - 6₃/m


In PyXtal, this representation can be easily obtained via the ``get_spg_representation()`` method. The first element of the output is the lattice index, and the second element is a 15x48 matrix representing the symmetry elements.

.. code-block:: Python
    from pyxtal.symmetry import Group
    spg = Group(227)
    id, matrix = spg.get_spg_representation() 
    print(id)
    print(matrix)

    >>> 13 # lattice id
    >>> # one-hot encoding of 15*48 matrix to represent the space group Fd-3m
    [[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]])
 


Wyckoff Site Symmetry
---------------------
For each space group, the Wyckoff positions are defined by the symmetry of the site. The Wyckoff positions are labeled with a letter and a number, where the letter indicates the type of site and the number indicates the multiplicity of that site. Below, we show the Wyckoff site symmetry for the space group 227 (``Fd-3m``) as an example. 

.. code-block:: python

    from pyxtal.symmetry import Group
    spg = Group(227)
    print(spg)
    wp = spg[-1]   # Get the last Wyckoff position 8a
    ss.to_beautiful_matrix_representation() # List symmetry elements

>>> -- Spacegroup --# 227 (Fd-3m)--
192i  site symm: 1
96h	  site symm: ..2
96g	  site symm: ..m
48f	  site symm: 2.mm
32e	  site symm: .3m
16d	  site symm: .-3m
16c	  site symm: .-3m
8b	  site symm: -43m
8a	  site symm: -43m
  
>>> Order Axis     1    -1   2    m    3    4    -4   -3   6    -6   
  0 ( 1  0  0):    1    0    1    0    0    0    1    0    0    0     -4
  0 ( 0  1  0):    1    0    1    0    0    0    1    0    0    0     -4
  0 ( 0  0  1):    1    0    1    0    0    0    1    0    0    0     -4
  1 ( 1  1  1):    1    0    0    0    1    0    0    0    0    0      3
  1 ( 1 -1 -1):    1    0    0    0    1    0    0    0    0    0      3
  1 (-1  1 -1):    1    0    0    0    1    0    0    0    0    0      3
  1 (-1 -1  1):    1    0    0    0    1    0    0    0    0    0      3
  2 ( 1 -1  0):    1    0    0    1    0    0    0    0    0    0      m
  2 ( 1  1  0):    1    0    0    1    0    0    0    0    0    0      m
  2 ( 0  1 -1):    1    0    0    1    0    0    0    0    0    0      m
  2 ( 0  1  1):    1    0    0    1    0    0    0    0    0    0      m
  2 (-1  0  1):    1    0    0    1    0    0    0    0    0    0      m
  2 ( 1  0  1):    1    0    0    1    0    0    0    0    0    0      m

In space group 227, the Wyckoff position ``8a`` indicates that there are 8 equivalent sites in the unit cell, with a site symmetry of ``-43m``. Unlike the Hermann-Mauguin notation, the site symmetry does not count the translation symmetry. Hence, it does not include the screw axis (e.g., ``2₁, 3₁, 4₁, 6₁``) or glide plane symmetry (``a, b, c, n, d``). There are 7 fundamental point group symmetries (``1, -1, 2, m, 3, 4, 6, -6``) and 5 additional compound group symmetries (``-3, 6, 2/m, 4/m, 6/m``). For ``8a`` in space group 227, its site symmetry ``-43m`` includes

1. 4-fold rotation axis (``-4``) @ [100] family directions,
2. 3-fold rotation axis (``3``) @ [111] family directions,
3. 2-fold rotation axis (``m``) @ [110] family directions.

Site Symmetry's Digital Representation
---------------------------------------------

For the application of deep learning, the Wyckoff site symmetry can be represented as a 15x7 matrix, where the rows and columns represent the symmetry elements. The value of each element in the matrix indicates the presence or absence of a symmetry operation. For example, a value of 1 indicates that the symmetry operation is present, while a value of 0 indicates that it is absent. Given that there exist a total of 13 site symmetries (``1, -1, 2, m, 3, 4, 6, -6, -3, 6, 2/m, 4/m, 6/m``), it can be further converted to an one-hot encoding format of (15x13) matrix via PyXtal as follows:

.. code-block:: python

    from pyxtal.symmetry import Group
    spg = Group(227)
    rep = wp.get_site_symmetry_object().to_one_hot()

>>> # one-hot encoding of the 8a site symmetry of space group 227
[[1 0 0 0 0 0 0 0 0 0 0 0 0]
 [1 0 0 0 0 0 0 0 0 0 0 0 0]
 [1 0 0 0 0 0 0 0 0 0 0 0 0]
 [0 1 0 0 0 0 0 0 0 0 0 0 0]
 [0 1 0 0 0 0 0 0 0 0 0 0 0]
 [0 1 0 0 0 0 0 0 0 0 0 0 0]
 [0 1 0 0 0 0 0 0 0 0 0 0 0]
 [0 0 1 0 0 0 0 0 0 0 0 0 0]
 [0 0 1 0 0 0 0 0 0 0 0 0 0]
 [0 0 1 0 0 0 0 0 0 0 0 0 0]
 [0 0 1 0 0 0 0 0 0 0 0 0 0]
 [0 0 1 0 0 0 0 0 0 0 0 0 0]
 [0 0 1 0 0 0 0 0 0 0 0 0 0]
 [1 0 0 0 0 0 0 0 0 0 0 0 0]
 [1 0 0 0 0 0 0 0 0 0 0 0 0]]


Complete list of Wyckoff Site Symmetry Table
--------------------------------------------

Using PyXtal, you can easily access the Wyckoff site symmetry for all 230 space groups as follows:

.. code-block:: python

    from pyxtal.symmetry import Group
    for g in range(1, 231):
        spg = Group(g)
        for wp in spg:
            wp.get_site_symmetry()
            print(spg.number, spg.symbol, wp.get_label(), wp.site_symm)

The following table lists the space group number, symbol, Wyckoff label, and site symmetry for all 230 space groups.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Space Group Number
     - Space Group Symbol
     - Wyckoff Label
     - Site Symmetry
   * - 1
     - P1
     - 1a
     - 1
   * - 2
     - P-1
     - 2i
     - 1
   * - 2
     - P-1
     - 1h
     - -1
   * - 2
     - P-1
     - 1g
     - -1
   * - 2
     - P-1
     - 1f
     - -1
   * - 2
     - P-1
     - 1e
     - -1
   * - 2
     - P-1
     - 1d
     - -1
   * - 2
     - P-1
     - 1c
     - -1
   * - 2
     - P-1
     - 1b
     - -1
   * - 2
     - P-1
     - 1a
     - -1
   * - 3
     - P2
     - 2e
     - 1
   * - 3
     - P2
     - 1d
     - 2
   * - 3
     - P2
     - 1c
     - 2
   * - 3
     - P2
     - 1b
     - 2
   * - 3
     - P2
     - 1a
     - 2
   * - 4
     - P21
     - 2a
     - 1
   * - 5
     - C2
     - 4c
     - 1
   * - 5
     - C2
     - 2b
     - 2
   * - 5
     - C2
     - 2a
     - 2
   * - 6
     - Pm
     - 2c
     - 1
   * - 6
     - Pm
     - 1b
     - m
   * - 6
     - Pm
     - 1a
     - m
   * - 7
     - Pc
     - 2a
     - 1
   * - 8
     - Cm
     - 4b
     - 1
   * - 8
     - Cm
     - 2a
     - m
   * - 9
     - Cc
     - 4a
     - 1
   * - 10
     - P2/m
     - 4o
     - 1
   * - 10
     - P2/m
     - 2n
     - m
   * - 10
     - P2/m
     - 2m
     - m
   * - 10
     - P2/m
     - 2l
     - 2
   * - 10
     - P2/m
     - 2k
     - 2
   * - 10
     - P2/m
     - 2j
     - 2
   * - 10
     - P2/m
     - 2i
     - 2
   * - 10
     - P2/m
     - 1h
     - 2/m
   * - 10
     - P2/m
     - 1g
     - 2/m
   * - 10
     - P2/m
     - 1f
     - 2/m
   * - 10
     - P2/m
     - 1e
     - 2/m
   * - 10
     - P2/m
     - 1d
     - 2/m
   * - 10
     - P2/m
     - 1c
     - 2/m
   * - 10
     - P2/m
     - 1b
     - 2/m
   * - 10
     - P2/m
     - 1a
     - 2/m
   * - 11
     - P21/m
     - 4f
     - 1
   * - 11
     - P21/m
     - 2e
     - m
   * - 11
     - P21/m
     - 2d
     - -1
   * - 11
     - P21/m
     - 2c
     - -1
   * - 11
     - P21/m
     - 2b
     - -1
   * - 11
     - P21/m
     - 2a
     - -1
   * - 12
     - C2/m
     - 8j
     - 1
   * - 12
     - C2/m
     - 4i
     - m
   * - 12
     - C2/m
     - 4h
     - 2
   * - 12
     - C2/m
     - 4g
     - 2
   * - 12
     - C2/m
     - 4f
     - -1
   * - 12
     - C2/m
     - 4e
     - -1
   * - 12
     - C2/m
     - 2d
     - 2/m
   * - 12
     - C2/m
     - 2c
     - 2/m
   * - 12
     - C2/m
     - 2b
     - 2/m
   * - 12
     - C2/m
     - 2a
     - 2/m
   * - 13
     - P2/c
     - 4g
     - 1
   * - 13
     - P2/c
     - 2f
     - 2
   * - 13
     - P2/c
     - 2e
     - 2
   * - 13
     - P2/c
     - 2d
     - -1
   * - 13
     - P2/c
     - 2c
     - -1
   * - 13
     - P2/c
     - 2b
     - -1
   * - 13
     - P2/c
     - 2a
     - -1
   * - 14
     - P21/c
     - 4e
     - 1
   * - 14
     - P21/c
     - 2d
     - -1
   * - 14
     - P21/c
     - 2c
     - -1
   * - 14
     - P21/c
     - 2b
     - -1
   * - 14
     - P21/c
     - 2a
     - -1
   * - 15
     - C2/c
     - 8f
     - 1
   * - 15
     - C2/c
     - 4e
     - 2
   * - 15
     - C2/c
     - 4d
     - -1
   * - 15
     - C2/c
     - 4c
     - -1
   * - 15
     - C2/c
     - 4b
     - -1
   * - 15
     - C2/c
     - 4a
     - -1
   * - 16
     - P222
     - 4u
     - 1
   * - 16
     - P222
     - 2t
     - ..2
   * - 16
     - P222
     - 2s
     - ..2
   * - 16
     - P222
     - 2r
     - ..2
   * - 16
     - P222
     - 2q
     - ..2
   * - 16
     - P222
     - 2p
     - .2.
   * - 16
     - P222
     - 2o
     - .2.
   * - 16
     - P222
     - 2n
     - .2.
   * - 16
     - P222
     - 2m
     - .2.
   * - 16
     - P222
     - 2l
     - 2..
   * - 16
     - P222
     - 2k
     - 2..
   * - 16
     - P222
     - 2j
     - 2..
   * - 16
     - P222
     - 2i
     - 2..
   * - 16
     - P222
     - 1h
     - 222
   * - 16
     - P222
     - 1g
     - 222
   * - 16
     - P222
     - 1f
     - 222
   * - 16
     - P222
     - 1e
     - 222
   * - 16
     - P222
     - 1d
     - 222
   * - 16
     - P222
     - 1c
     - 222
   * - 16
     - P222
     - 1b
     - 222
   * - 16
     - P222
     - 1a
     - 222
   * - 17
     - P2221
     - 4e
     - 1
   * - 17
     - P2221
     - 2d
     - .2.
   * - 17
     - P2221
     - 2c
     - .2.
   * - 17
     - P2221
     - 2b
     - 2..
   * - 17
     - P2221
     - 2a
     - 2..
   * - 18
     - P21212
     - 4c
     - 1
   * - 18
     - P21212
     - 2b
     - ..2
   * - 18
     - P21212
     - 2a
     - ..2
   * - 19
     - P212121
     - 4a
     - 1
   * - 20
     - C2221
     - 8c
     - 1
   * - 20
     - C2221
     - 4b
     - .2.
   * - 20
     - C2221
     - 4a
     - 2..
   * - 21
     - C222
     - 8l
     - 1
   * - 21
     - C222
     - 4k
     - ..2
   * - 21
     - C222
     - 4j
     - ..2
   * - 21
     - C222
     - 4i
     - ..2
   * - 21
     - C222
     - 4h
     - .2.
   * - 21
     - C222
     - 4g
     - .2.
   * - 21
     - C222
     - 4f
     - 2..
   * - 21
     - C222
     - 4e
     - 2..
   * - 21
     - C222
     - 2d
     - 222
   * - 21
     - C222
     - 2c
     - 222
   * - 21
     - C222
     - 2b
     - 222
   * - 21
     - C222
     - 2a
     - 222
   * - 22
     - F222
     - 16k
     - 1
   * - 22
     - F222
     - 8j
     - 2..
   * - 22
     - F222
     - 8i
     - .2.
   * - 22
     - F222
     - 8h
     - ..2
   * - 22
     - F222
     - 8g
     - ..2
   * - 22
     - F222
     - 8f
     - .2.
   * - 22
     - F222
     - 8e
     - 2..
   * - 22
     - F222
     - 4d
     - 222
   * - 22
     - F222
     - 4c
     - 222
   * - 22
     - F222
     - 4b
     - 222
   * - 22
     - F222
     - 4a
     - 222
   * - 23
     - I222
     - 8k
     - 1
   * - 23
     - I222
     - 4j
     - ..2
   * - 23
     - I222
     - 4i
     - ..2
   * - 23
     - I222
     - 4h
     - .2.
   * - 23
     - I222
     - 4g
     - .2.
   * - 23
     - I222
     - 4f
     - 2..
   * - 23
     - I222
     - 4e
     - 2..
   * - 23
     - I222
     - 2d
     - 222
   * - 23
     - I222
     - 2c
     - 222
   * - 23
     - I222
     - 2b
     - 222
   * - 23
     - I222
     - 2a
     - 222
   * - 24
     - I212121
     - 8d
     - 1
   * - 24
     - I212121
     - 4c
     - ..2
   * - 24
     - I212121
     - 4b
     - .2.
   * - 24
     - I212121
     - 4a
     - 2..
   * - 25
     - Pmm2
     - 4i
     - 1
   * - 25
     - Pmm2
     - 2h
     - m..
   * - 25
     - Pmm2
     - 2g
     - m..
   * - 25
     - Pmm2
     - 2f
     - .m.
   * - 25
     - Pmm2
     - 2e
     - .m.
   * - 25
     - Pmm2
     - 1d
     - mm2
   * - 25
     - Pmm2
     - 1c
     - mm2
   * - 25
     - Pmm2
     - 1b
     - mm2
   * - 25
     - Pmm2
     - 1a
     - mm2
   * - 26
     - Pmc21
     - 4c
     - 1
   * - 26
     - Pmc21
     - 2b
     - m..
   * - 26
     - Pmc21
     - 2a
     - m..
   * - 27
     - Pcc2
     - 4e
     - 1
   * - 27
     - Pcc2
     - 2d
     - ..2
   * - 27
     - Pcc2
     - 2c
     - ..2
   * - 27
     - Pcc2
     - 2b
     - ..2
   * - 27
     - Pcc2
     - 2a
     - ..2
   * - 28
     - Pma2
     - 4d
     - 1
   * - 28
     - Pma2
     - 2c
     - m..
   * - 28
     - Pma2
     - 2b
     - ..2
   * - 28
     - Pma2
     - 2a
     - ..2
   * - 29
     - Pca21
     - 4a
     - 1
   * - 30
     - Pnc2
     - 4c
     - 1
   * - 30
     - Pnc2
     - 2b
     - ..2
   * - 30
     - Pnc2
     - 2a
     - ..2
   * - 31
     - Pmn21
     - 4b
     - 1
   * - 31
     - Pmn21
     - 2a
     - m..
   * - 32
     - Pba2
     - 4c
     - 1
   * - 32
     - Pba2
     - 2b
     - ..2
   * - 32
     - Pba2
     - 2a
     - ..2
   * - 33
     - Pna21
     - 4a
     - 1
   * - 34
     - Pnn2
     - 4c
     - 1
   * - 34
     - Pnn2
     - 2b
     - ..2
   * - 34
     - Pnn2
     - 2a
     - ..2
   * - 35
     - Cmm2
     - 8f
     - 1
   * - 35
     - Cmm2
     - 4e
     - m..
   * - 35
     - Cmm2
     - 4d
     - .m.
   * - 35
     - Cmm2
     - 4c
     - ..2
   * - 35
     - Cmm2
     - 2b
     - mm2
   * - 35
     - Cmm2
     - 2a
     - mm2
   * - 36
     - Cmc21
     - 8b
     - 1
   * - 36
     - Cmc21
     - 4a
     - m..
   * - 37
     - Ccc2
     - 8d
     - 1
   * - 37
     - Ccc2
     - 4c
     - ..2
   * - 37
     - Ccc2
     - 4b
     - ..2
   * - 37
     - Ccc2
     - 4a
     - ..2
   * - 38
     - Amm2
     - 8f
     - 1
   * - 38
     - Amm2
     - 4e
     - m..
   * - 38
     - Amm2
     - 4d
     - m..
   * - 38
     - Amm2
     - 4c
     - .m.
   * - 38
     - Amm2
     - 2b
     - mm2
   * - 38
     - Amm2
     - 2a
     - mm2
   * - 39
     - Aem2
     - 8d
     - 1
   * - 39
     - Aem2
     - 4c
     - .m.
   * - 39
     - Aem2
     - 4b
     - ..2
   * - 39
     - Aem2
     - 4a
     - ..2
   * - 40
     - Ama2
     - 8c
     - 1
   * - 40
     - Ama2
     - 4b
     - m..
   * - 40
     - Ama2
     - 4a
     - ..2
   * - 41
     - Aea2
     - 8b
     - 1
   * - 41
     - Aea2
     - 4a
     - ..2
   * - 42
     - Fmm2
     - 16e
     - 1
   * - 42
     - Fmm2
     - 8d
     - .m.
   * - 42
     - Fmm2
     - 8c
     - m..
   * - 42
     - Fmm2
     - 8b
     - ..2
   * - 42
     - Fmm2
     - 4a
     - mm2
   * - 43
     - Fdd2
     - 16b
     - 1
   * - 43
     - Fdd2
     - 8a
     - ..2
   * - 44
     - Imm2
     - 8e
     - 1
   * - 44
     - Imm2
     - 4d
     - m..
   * - 44
     - Imm2
     - 4c
     - .m.
   * - 44
     - Imm2
     - 2b
     - mm2
   * - 44
     - Imm2
     - 2a
     - mm2
   * - 45
     - Iba2
     - 8c
     - 1
   * - 45
     - Iba2
     - 4b
     - ..2
   * - 45
     - Iba2
     - 4a
     - ..2
   * - 46
     - Ima2
     - 8c
     - 1
   * - 46
     - Ima2
     - 4b
     - m..
   * - 46
     - Ima2
     - 4a
     - ..2
   * - 47
     - Pmmm
     - 8A
     - 1
   * - 47
     - Pmmm
     - 4z
     - ..m
   * - 47
     - Pmmm
     - 4y
     - ..m
   * - 47
     - Pmmm
     - 4x
     - .m.
   * - 47
     - Pmmm
     - 4w
     - .m.
   * - 47
     - Pmmm
     - 4v
     - m..
   * - 47
     - Pmmm
     - 4u
     - m..
   * - 47
     - Pmmm
     - 2t
     - mm2
   * - 47
     - Pmmm
     - 2s
     - mm2
   * - 47
     - Pmmm
     - 2r
     - mm2
   * - 47
     - Pmmm
     - 2q
     - mm2
   * - 47
     - Pmmm
     - 2p
     - m2m
   * - 47
     - Pmmm
     - 2o
     - m2m
   * - 47
     - Pmmm
     - 2n
     - m2m
   * - 47
     - Pmmm
     - 2m
     - m2m
   * - 47
     - Pmmm
     - 2l
     - 2mm
   * - 47
     - Pmmm
     - 2k
     - 2mm
   * - 47
     - Pmmm
     - 2j
     - 2mm
   * - 47
     - Pmmm
     - 2i
     - 2mm
   * - 47
     - Pmmm
     - 1h
     - mmm
   * - 47
     - Pmmm
     - 1g
     - mmm
   * - 47
     - Pmmm
     - 1f
     - mmm
   * - 47
     - Pmmm
     - 1e
     - mmm
   * - 47
     - Pmmm
     - 1d
     - mmm
   * - 47
     - Pmmm
     - 1c
     - mmm
   * - 47
     - Pmmm
     - 1b
     - mmm
   * - 47
     - Pmmm
     - 1a
     - mmm
   * - 48
     - Pnnn
     - 8m
     - 1
   * - 48
     - Pnnn
     - 4l
     - ..2
   * - 48
     - Pnnn
     - 4k
     - ..2
   * - 48
     - Pnnn
     - 4j
     - .2.
   * - 48
     - Pnnn
     - 4i
     - .2.
   * - 48
     - Pnnn
     - 4h
     - 2..
   * - 48
     - Pnnn
     - 4g
     - 2..
   * - 48
     - Pnnn
     - 4f
     - -1
   * - 48
     - Pnnn
     - 4e
     - -1
   * - 48
     - Pnnn
     - 2d
     - 222
   * - 48
     - Pnnn
     - 2c
     - 222
   * - 48
     - Pnnn
     - 2b
     - 222
   * - 48
     - Pnnn
     - 2a
     - 222
   * - 49
     - Pccm
     - 8r
     - 1
   * - 49
     - Pccm
     - 4q
     - ..m
   * - 49
     - Pccm
     - 4p
     - ..2
   * - 49
     - Pccm
     - 4o
     - ..2
   * - 49
     - Pccm
     - 4n
     - ..2
   * - 49
     - Pccm
     - 4m
     - ..2
   * - 49
     - Pccm
     - 4l
     - .2.
   * - 49
     - Pccm
     - 4k
     - .2.
   * - 49
     - Pccm
     - 4j
     - 2..
   * - 49
     - Pccm
     - 4i
     - 2..
   * - 49
     - Pccm
     - 2h
     - 222
   * - 49
     - Pccm
     - 2g
     - 222
   * - 49
     - Pccm
     - 2f
     - 222
   * - 49
     - Pccm
     - 2e
     - 222
   * - 49
     - Pccm
     - 2d
     - ..2/m
   * - 49
     - Pccm
     - 2c
     - ..2/m
   * - 49
     - Pccm
     - 2b
     - ..2/m
   * - 49
     - Pccm
     - 2a
     - ..2/m
   * - 50
     - Pban
     - 8m
     - 1
   * - 50
     - Pban
     - 4l
     - ..2
   * - 50
     - Pban
     - 4k
     - ..2
   * - 50
     - Pban
     - 4j
     - .2.
   * - 50
     - Pban
     - 4i
     - .2.
   * - 50
     - Pban
     - 4h
     - 2..
   * - 50
     - Pban
     - 4g
     - 2..
   * - 50
     - Pban
     - 4f
     - -1
   * - 50
     - Pban
     - 4e
     - -1
   * - 50
     - Pban
     - 2d
     - 222
   * - 50
     - Pban
     - 2c
     - 222
   * - 50
     - Pban
     - 2b
     - 222
   * - 50
     - Pban
     - 2a
     - 222
   * - 51
     - Pmma
     - 8l
     - 1
   * - 51
     - Pmma
     - 4k
     - m..
   * - 51
     - Pmma
     - 4j
     - .m.
   * - 51
     - Pmma
     - 4i
     - .m.
   * - 51
     - Pmma
     - 4h
     - .2.
   * - 51
     - Pmma
     - 4g
     - .2.
   * - 51
     - Pmma
     - 2f
     - mm2
   * - 51
     - Pmma
     - 2e
     - mm2
   * - 51
     - Pmma
     - 2d
     - .2/m.
   * - 51
     - Pmma
     - 2c
     - .2/m.
   * - 51
     - Pmma
     - 2b
     - .2/m.
   * - 51
     - Pmma
     - 2a
     - .2/m.
   * - 52
     - Pnna
     - 8e
     - 1
   * - 52
     - Pnna
     - 4d
     - 2..
   * - 52
     - Pnna
     - 4c
     - ..2
   * - 52
     - Pnna
     - 4b
     - -1
   * - 52
     - Pnna
     - 4a
     - -1
   * - 53
     - Pmna
     - 8i
     - 1
   * - 53
     - Pmna
     - 4h
     - m..
   * - 53
     - Pmna
     - 4g
     - .2.
   * - 53
     - Pmna
     - 4f
     - 2..
   * - 53
     - Pmna
     - 4e
     - 2..
   * - 53
     - Pmna
     - 2d
     - 2/m..
   * - 53
     - Pmna
     - 2c
     - 2/m..
   * - 53
     - Pmna
     - 2b
     - 2/m..
   * - 53
     - Pmna
     - 2a
     - 2/m..
   * - 54
     - Pcca
     - 8f
     - 1
   * - 54
     - Pcca
     - 4e
     - ..2
   * - 54
     - Pcca
     - 4d
     - ..2
   * - 54
     - Pcca
     - 4c
     - .2.
   * - 54
     - Pcca
     - 4b
     - -1
   * - 54
     - Pcca
     - 4a
     - -1
   * - 55
     - Pbam
     - 8i
     - 1
   * - 55
     - Pbam
     - 4h
     - ..m
   * - 55
     - Pbam
     - 4g
     - ..m
   * - 55
     - Pbam
     - 4f
     - ..2
   * - 55
     - Pbam
     - 4e
     - ..2
   * - 55
     - Pbam
     - 2d
     - ..2/m
   * - 55
     - Pbam
     - 2c
     - ..2/m
   * - 55
     - Pbam
     - 2b
     - ..2/m
   * - 55
     - Pbam
     - 2a
     - ..2/m
   * - 56
     - Pccn
     - 8e
     - 1
   * - 56
     - Pccn
     - 4d
     - ..2
   * - 56
     - Pccn
     - 4c
     - ..2
   * - 56
     - Pccn
     - 4b
     - -1
   * - 56
     - Pccn
     - 4a
     - -1
   * - 57
     - Pbcm
     - 8e
     - 1
   * - 57
     - Pbcm
     - 4d
     - ..m
   * - 57
     - Pbcm
     - 4c
     - 2..
   * - 57
     - Pbcm
     - 4b
     - -1
   * - 57
     - Pbcm
     - 4a
     - -1
   * - 58
     - Pnnm
     - 8h
     - 1
   * - 58
     - Pnnm
     - 4g
     - ..m
   * - 58
     - Pnnm
     - 4f
     - ..2
   * - 58
     - Pnnm
     - 4e
     - ..2
   * - 58
     - Pnnm
     - 2d
     - ..2/m
   * - 58
     - Pnnm
     - 2c
     - ..2/m
   * - 58
     - Pnnm
     - 2b
     - ..2/m
   * - 58
     - Pnnm
     - 2a
     - ..2/m
   * - 59
     - Pmmn
     - 8g
     - 1
   * - 59
     - Pmmn
     - 4f
     - .m.
   * - 59
     - Pmmn
     - 4e
     - m..
   * - 59
     - Pmmn
     - 4d
     - -1
   * - 59
     - Pmmn
     - 4c
     - -1
   * - 59
     - Pmmn
     - 2b
     - mm2
   * - 59
     - Pmmn
     - 2a
     - mm2
   * - 60
     - Pbcn
     - 8d
     - 1
   * - 60
     - Pbcn
     - 4c
     - .2.
   * - 60
     - Pbcn
     - 4b
     - -1
   * - 60
     - Pbcn
     - 4a
     - -1
   * - 61
     - Pbca
     - 8c
     - 1
   * - 61
     - Pbca
     - 4b
     - -1
   * - 61
     - Pbca
     - 4a
     - -1
   * - 62
     - Pnma
     - 8d
     - 1
   * - 62
     - Pnma
     - 4c
     - .m.
   * - 62
     - Pnma
     - 4b
     - -1
   * - 62
     - Pnma
     - 4a
     - -1
   * - 63
     - Cmcm
     - 16h
     - 1
   * - 63
     - Cmcm
     - 8g
     - ..m
   * - 63
     - Cmcm
     - 8f
     - m..
   * - 63
     - Cmcm
     - 8e
     - 2..
   * - 63
     - Cmcm
     - 8d
     - -1
   * - 63
     - Cmcm
     - 4c
     - m2m
   * - 63
     - Cmcm
     - 4b
     - 2/m..
   * - 63
     - Cmcm
     - 4a
     - 2/m..
   * - 64
     - Cmce
     - 16g
     - 1
   * - 64
     - Cmce
     - 8f
     - m..
   * - 64
     - Cmce
     - 8e
     - .2.
   * - 64
     - Cmce
     - 8d
     - 2..
   * - 64
     - Cmce
     - 8c
     - -1
   * - 64
     - Cmce
     - 4b
     - 2/m..
   * - 64
     - Cmce
     - 4a
     - 2/m..
   * - 65
     - Cmmm
     - 16r
     - 1
   * - 65
     - Cmmm
     - 8q
     - ..m
   * - 65
     - Cmmm
     - 8p
     - ..m
   * - 65
     - Cmmm
     - 8o
     - .m.
   * - 65
     - Cmmm
     - 8n
     - m..
   * - 65
     - Cmmm
     - 8m
     - ..2
   * - 65
     - Cmmm
     - 4l
     - mm2
   * - 65
     - Cmmm
     - 4k
     - mm2
   * - 65
     - Cmmm
     - 4j
     - m2m
   * - 65
     - Cmmm
     - 4i
     - m2m
   * - 65
     - Cmmm
     - 4h
     - 2mm
   * - 65
     - Cmmm
     - 4g
     - 2mm
   * - 65
     - Cmmm
     - 4f
     - ..2/m
   * - 65
     - Cmmm
     - 4e
     - ..2/m
   * - 65
     - Cmmm
     - 2d
     - mmm
   * - 65
     - Cmmm
     - 2c
     - mmm
   * - 65
     - Cmmm
     - 2b
     - mmm
   * - 65
     - Cmmm
     - 2a
     - mmm
   * - 66
     - Cccm
     - 16m
     - 1
   * - 66
     - Cccm
     - 8l
     - ..m
   * - 66
     - Cccm
     - 8k
     - ..2
   * - 66
     - Cccm
     - 8j
     - ..2
   * - 66
     - Cccm
     - 8i
     - ..2
   * - 66
     - Cccm
     - 8h
     - .2.
   * - 66
     - Cccm
     - 8g
     - 2..
   * - 66
     - Cccm
     - 4f
     - ..2/m
   * - 66
     - Cccm
     - 4e
     - ..2/m
   * - 66
     - Cccm
     - 4d
     - ..2/m
   * - 66
     - Cccm
     - 4c
     - ..2/m
   * - 66
     - Cccm
     - 4b
     - 222
   * - 66
     - Cccm
     - 4a
     - 222
   * - 67
     - Cmme
     - 16o
     - 1
   * - 67
     - Cmme
     - 8n
     - .m.
   * - 67
     - Cmme
     - 8m
     - m..
   * - 67
     - Cmme
     - 8l
     - ..2
   * - 67
     - Cmme
     - 8k
     - .2.
   * - 67
     - Cmme
     - 8j
     - .2.
   * - 67
     - Cmme
     - 8i
     - 2..
   * - 67
     - Cmme
     - 8h
     - 2..
   * - 67
     - Cmme
     - 4g
     - mm2
   * - 67
     - Cmme
     - 4f
     - .2/m.
   * - 67
     - Cmme
     - 4e
     - .2/m.
   * - 67
     - Cmme
     - 4d
     - 2/m..
   * - 67
     - Cmme
     - 4c
     - 2/m..
   * - 67
     - Cmme
     - 4b
     - 222
   * - 67
     - Cmme
     - 4a
     - 222
   * - 68
     - Ccce
     - 16i
     - 1
   * - 68
     - Ccce
     - 8h
     - ..2
   * - 68
     - Ccce
     - 8g
     - ..2
   * - 68
     - Ccce
     - 8f
     - .2.
   * - 68
     - Ccce
     - 8e
     - 2..
   * - 68
     - Ccce
     - 8d
     - -1
   * - 68
     - Ccce
     - 8c
     - -1
   * - 68
     - Ccce
     - 4b
     - 222
   * - 68
     - Ccce
     - 4a
     - 222
   * - 69
     - Fmmm
     - 32p
     - 1
   * - 69
     - Fmmm
     - 16o
     - ..m
   * - 69
     - Fmmm
     - 16n
     - .m.
   * - 69
     - Fmmm
     - 16m
     - m..
   * - 69
     - Fmmm
     - 16l
     - 2..
   * - 69
     - Fmmm
     - 16k
     - .2.
   * - 69
     - Fmmm
     - 16j
     - ..2
   * - 69
     - Fmmm
     - 8i
     - mm2
   * - 69
     - Fmmm
     - 8h
     - m2m
   * - 69
     - Fmmm
     - 8g
     - 2mm
   * - 69
     - Fmmm
     - 8f
     - 222
   * - 69
     - Fmmm
     - 8e
     - ..2/m
   * - 69
     - Fmmm
     - 8d
     - .2/m.
   * - 69
     - Fmmm
     - 8c
     - 2/m..
   * - 69
     - Fmmm
     - 4b
     - mmm
   * - 69
     - Fmmm
     - 4a
     - mmm
   * - 70
     - Fddd
     - 32h
     - 1
   * - 70
     - Fddd
     - 16g
     - ..2
   * - 70
     - Fddd
     - 16f
     - .2.
   * - 70
     - Fddd
     - 16e
     - 2..
   * - 70
     - Fddd
     - 16d
     - -1
   * - 70
     - Fddd
     - 16c
     - -1
   * - 70
     - Fddd
     - 8b
     - 222
   * - 70
     - Fddd
     - 8a
     - 222
   * - 71
     - Immm
     - 16o
     - 1
   * - 71
     - Immm
     - 8n
     - ..m
   * - 71
     - Immm
     - 8m
     - .m.
   * - 71
     - Immm
     - 8l
     - m..
   * - 71
     - Immm
     - 8k
     - -1
   * - 71
     - Immm
     - 4j
     - mm2
   * - 71
     - Immm
     - 4i
     - mm2
   * - 71
     - Immm
     - 4h
     - m2m
   * - 71
     - Immm
     - 4g
     - m2m
   * - 71
     - Immm
     - 4f
     - 2mm
   * - 71
     - Immm
     - 4e
     - 2mm
   * - 71
     - Immm
     - 2d
     - mmm
   * - 71
     - Immm
     - 2c
     - mmm
   * - 71
     - Immm
     - 2b
     - mmm
   * - 71
     - Immm
     - 2a
     - mmm
   * - 72
     - Ibam
     - 16k
     - 1
   * - 72
     - Ibam
     - 8j
     - ..m
   * - 72
     - Ibam
     - 8i
     - ..2
   * - 72
     - Ibam
     - 8h
     - ..2
   * - 72
     - Ibam
     - 8g
     - .2.
   * - 72
     - Ibam
     - 8f
     - 2..
   * - 72
     - Ibam
     - 8e
     - -1
   * - 72
     - Ibam
     - 4d
     - ..2/m
   * - 72
     - Ibam
     - 4c
     - ..2/m
   * - 72
     - Ibam
     - 4b
     - 222
   * - 72
     - Ibam
     - 4a
     - 222
   * - 73
     - Ibca
     - 16f
     - 1
   * - 73
     - Ibca
     - 8e
     - ..2
   * - 73
     - Ibca
     - 8d
     - .2.
   * - 73
     - Ibca
     - 8c
     - 2..
   * - 73
     - Ibca
     - 8b
     - -1
   * - 73
     - Ibca
     - 8a
     - -1
   * - 74
     - Imma
     - 16j
     - 1
   * - 74
     - Imma
     - 8i
     - .m.
   * - 74
     - Imma
     - 8h
     - m..
   * - 74
     - Imma
     - 8g
     - .2.
   * - 74
     - Imma
     - 8f
     - 2..
   * - 74
     - Imma
     - 4e
     - mm2
   * - 74
     - Imma
     - 4d
     - .2/m.
   * - 74
     - Imma
     - 4c
     - .2/m.
   * - 74
     - Imma
     - 4b
     - 2/m..
   * - 74
     - Imma
     - 4a
     - 2/m..
   * - 75
     - P4
     - 4d
     - 1
   * - 75
     - P4
     - 2c
     - 2..
   * - 75
     - P4
     - 1b
     - 4..
   * - 75
     - P4
     - 1a
     - 4..
   * - 76
     - P41
     - 4a
     - 1
   * - 77
     - P42
     - 4d
     - 1
   * - 77
     - P42
     - 2c
     - 2..
   * - 77
     - P42
     - 2b
     - 2..
   * - 77
     - P42
     - 2a
     - 2..
   * - 78
     - P43
     - 4a
     - 1
   * - 79
     - I4
     - 8c
     - 1
   * - 79
     - I4
     - 4b
     - 2..
   * - 79
     - I4
     - 2a
     - 4..
   * - 80
     - I41
     - 8b
     - 1
   * - 80
     - I41
     - 4a
     - 2..
   * - 81
     - P-4
     - 4h
     - 1
   * - 81
     - P-4
     - 2g
     - 2..
   * - 81
     - P-4
     - 2f
     - 2..
   * - 81
     - P-4
     - 2e
     - 2..
   * - 81
     - P-4
     - 1d
     - -4..
   * - 81
     - P-4
     - 1c
     - -4..
   * - 81
     - P-4
     - 1b
     - -4..
   * - 81
     - P-4
     - 1a
     - -4..
   * - 82
     - I-4
     - 8g
     - 1
   * - 82
     - I-4
     - 4f
     - 2..
   * - 82
     - I-4
     - 4e
     - 2..
   * - 82
     - I-4
     - 2d
     - -4..
   * - 82
     - I-4
     - 2c
     - -4..
   * - 82
     - I-4
     - 2b
     - -4..
   * - 82
     - I-4
     - 2a
     - -4..
   * - 83
     - P4/m
     - 8l
     - 1
   * - 83
     - P4/m
     - 4k
     - m..
   * - 83
     - P4/m
     - 4j
     - m..
   * - 83
     - P4/m
     - 4i
     - 2..
   * - 83
     - P4/m
     - 2h
     - 4..
   * - 83
     - P4/m
     - 2g
     - 4..
   * - 83
     - P4/m
     - 2f
     - 2/m..
   * - 83
     - P4/m
     - 2e
     - 2/m..
   * - 83
     - P4/m
     - 1d
     - 4/m..
   * - 83
     - P4/m
     - 1c
     - 4/m..
   * - 83
     - P4/m
     - 1b
     - 4/m..
   * - 83
     - P4/m
     - 1a
     - 4/m..
   * - 84
     - P42/m
     - 8k
     - 1
   * - 84
     - P42/m
     - 4j
     - m..
   * - 84
     - P42/m
     - 4i
     - 2..
   * - 84
     - P42/m
     - 4h
     - 2..
   * - 84
     - P42/m
     - 4g
     - 2..
   * - 84
     - P42/m
     - 2f
     - -4..
   * - 84
     - P42/m
     - 2e
     - -4..
   * - 84
     - P42/m
     - 2d
     - 2/m..
   * - 84
     - P42/m
     - 2c
     - 2/m..
   * - 84
     - P42/m
     - 2b
     - 2/m..
   * - 84
     - P42/m
     - 2a
     - 2/m..
   * - 85
     - P4/n
     - 8g
     - 1
   * - 85
     - P4/n
     - 4f
     - 2..
   * - 85
     - P4/n
     - 4e
     - -1
   * - 85
     - P4/n
     - 4d
     - -1
   * - 85
     - P4/n
     - 2c
     - 4..
   * - 85
     - P4/n
     - 2b
     - -4..
   * - 85
     - P4/n
     - 2a
     - -4..
   * - 86
     - P42/n
     - 8g
     - 1
   * - 86
     - P42/n
     - 4f
     - 2..
   * - 86
     - P42/n
     - 4e
     - 2..
   * - 86
     - P42/n
     - 4d
     - -1
   * - 86
     - P42/n
     - 4c
     - -1
   * - 86
     - P42/n
     - 2b
     - -4..
   * - 86
     - P42/n
     - 2a
     - -4..
   * - 87
     - I4/m
     - 16i
     - 1
   * - 87
     - I4/m
     - 8h
     - m..
   * - 87
     - I4/m
     - 8g
     - 2..
   * - 87
     - I4/m
     - 8f
     - -1
   * - 87
     - I4/m
     - 4e
     - 4..
   * - 87
     - I4/m
     - 4d
     - -4..
   * - 87
     - I4/m
     - 4c
     - 2/m..
   * - 87
     - I4/m
     - 2b
     - 4/m..
   * - 87
     - I4/m
     - 2a
     - 4/m..
   * - 88
     - I41/a
     - 16f
     - 1
   * - 88
     - I41/a
     - 8e
     - 2..
   * - 88
     - I41/a
     - 8d
     - -1
   * - 88
     - I41/a
     - 8c
     - -1
   * - 88
     - I41/a
     - 4b
     - -4..
   * - 88
     - I41/a
     - 4a
     - -4..
   * - 89
     - P422
     - 8p
     - 1
   * - 89
     - P422
     - 4o
     - .2.
   * - 89
     - P422
     - 4n
     - .2.
   * - 89
     - P422
     - 4m
     - .2.
   * - 89
     - P422
     - 4l
     - .2.
   * - 89
     - P422
     - 4k
     - ..2
   * - 89
     - P422
     - 4j
     - ..2
   * - 89
     - P422
     - 4i
     - 2..
   * - 89
     - P422
     - 2h
     - 4..
   * - 89
     - P422
     - 2g
     - 4..
   * - 89
     - P422
     - 2f
     - 222.\
   * - 89
     - P422
     - 2e
     - 222.\
   * - 89
     - P422
     - 1d
     - 422
   * - 89
     - P422
     - 1c
     - 422
   * - 89
     - P422
     - 1b
     - 422
   * - 89
     - P422
     - 1a
     - 422
   * - 90
     - P4212
     - 8g
     - 1
   * - 90
     - P4212
     - 4f
     - ..2
   * - 90
     - P4212
     - 4e
     - ..2
   * - 90
     - P4212
     - 4d
     - 2..
   * - 90
     - P4212
     - 2c
     - 4..
   * - 90
     - P4212
     - 2b
     - 2.22
   * - 90
     - P4212
     - 2a
     - 2.22
   * - 91
     - P4122
     - 8d
     - 1
   * - 91
     - P4122
     - 4c
     - ..2
   * - 91
     - P4122
     - 4b
     - .2.
   * - 91
     - P4122
     - 4a
     - .2.
   * - 92
     - P41212
     - 8b
     - 1
   * - 92
     - P41212
     - 4a
     - ..2
   * - 93
     - P4222
     - 8p
     - 1
   * - 93
     - P4222
     - 4o
     - ..2
   * - 93
     - P4222
     - 4n
     - ..2
   * - 93
     - P4222
     - 4m
     - .2.
   * - 93
     - P4222
     - 4l
     - .2.
   * - 93
     - P4222
     - 4k
     - .2.
   * - 93
     - P4222
     - 4j
     - .2.
   * - 93
     - P4222
     - 4i
     - 2..
   * - 93
     - P4222
     - 4h
     - 2..
   * - 93
     - P4222
     - 4g
     - 2..
   * - 93
     - P4222
     - 2f
     - 2.22
   * - 93
     - P4222
     - 2e
     - 2.22
   * - 93
     - P4222
     - 2d
     - 222.\
   * - 93
     - P4222
     - 2c
     - 222.\
   * - 93
     - P4222
     - 2b
     - 222.\
   * - 93
     - P4222
     - 2a
     - 222.\
   * - 94
     - P42212
     - 8g
     - 1
   * - 94
     - P42212
     - 4f
     - ..2
   * - 94
     - P42212
     - 4e
     - ..2
   * - 94
     - P42212
     - 4d
     - 2..
   * - 94
     - P42212
     - 4c
     - 2..
   * - 94
     - P42212
     - 2b
     - 2.22
   * - 94
     - P42212
     - 2a
     - 2.22
   * - 95
     - P4322
     - 8d
     - 1
   * - 95
     - P4322
     - 4c
     - ..2
   * - 95
     - P4322
     - 4b
     - .2.
   * - 95
     - P4322
     - 4a
     - .2.
   * - 96
     - P43212
     - 8b
     - 1
   * - 96
     - P43212
     - 4a
     - ..2
   * - 97
     - I422
     - 16k
     - 1
   * - 97
     - I422
     - 8j
     - ..2
   * - 97
     - I422
     - 8i
     - .2.
   * - 97
     - I422
     - 8h
     - .2.
   * - 97
     - I422
     - 8g
     - ..2
   * - 97
     - I422
     - 8f
     - 2..
   * - 97
     - I422
     - 4e
     - 4..
   * - 97
     - I422
     - 4d
     - 2.22
   * - 97
     - I422
     - 4c
     - 222.\
   * - 97
     - I422
     - 2b
     - 422
   * - 97
     - I422
     - 2a
     - 422
   * - 98
     - I4122
     - 16g
     - 1
   * - 98
     - I4122
     - 8f
     - .2.
   * - 98
     - I4122
     - 8e
     - ..2
   * - 98
     - I4122
     - 8d
     - ..2
   * - 98
     - I4122
     - 8c
     - 2..
   * - 98
     - I4122
     - 4b
     - 2.22
   * - 98
     - I4122
     - 4a
     - 2.22
   * - 99
     - P4mm
     - 8g
     - 1
   * - 99
     - P4mm
     - 4f
     - .m.
   * - 99
     - P4mm
     - 4e
     - .m.
   * - 99
     - P4mm
     - 4d
     - ..m
   * - 99
     - P4mm
     - 2c
     - 2mm.
   * - 99
     - P4mm
     - 1b
     - 4mm
   * - 99
     - P4mm
     - 1a
     - 4mm
   * - 100
     - P4bm
     - 8d
     - 1
   * - 100
     - P4bm
     - 4c
     - ..m
   * - 100
     - P4bm
     - 2b
     - 2.mm
   * - 100
     - P4bm
     - 2a
     - 4..
   * - 101
     - P42cm
     - 8e
     - 1
   * - 101
     - P42cm
     - 4d
     - ..m
   * - 101
     - P42cm
     - 4c
     - 2..
   * - 101
     - P42cm
     - 2b
     - 2.mm
   * - 101
     - P42cm
     - 2a
     - 2.mm
   * - 102
     - P42nm
     - 8d
     - 1
   * - 102
     - P42nm
     - 4c
     - ..m
   * - 102
     - P42nm
     - 4b
     - 2..
   * - 102
     - P42nm
     - 2a
     - 2.mm
   * - 103
     - P4cc
     - 8d
     - 1
   * - 103
     - P4cc
     - 4c
     - 2..
   * - 103
     - P4cc
     - 2b
     - 4..
   * - 103
     - P4cc
     - 2a
     - 4..
   * - 104
     - P4nc
     - 8c
     - 1
   * - 104
     - P4nc
     - 4b
     - 2..
   * - 104
     - P4nc
     - 2a
     - 4..
   * - 105
     - P42mc
     - 8f
     - 1
   * - 105
     - P42mc
     - 4e
     - .m.
   * - 105
     - P42mc
     - 4d
     - .m.
   * - 105
     - P42mc
     - 2c
     - 2mm.
   * - 105
     - P42mc
     - 2b
     - 2mm.
   * - 105
     - P42mc
     - 2a
     - 2mm.
   * - 106
     - P42bc
     - 8c
     - 1
   * - 106
     - P42bc
     - 4b
     - 2..
   * - 106
     - P42bc
     - 4a
     - 2..
   * - 107
     - I4mm
     - 16e
     - 1
   * - 107
     - I4mm
     - 8d
     - .m.
   * - 107
     - I4mm
     - 8c
     - ..m
   * - 107
     - I4mm
     - 4b
     - 2mm.
   * - 107
     - I4mm
     - 2a
     - 4mm
   * - 108
     - I4cm
     - 16d
     - 1
   * - 108
     - I4cm
     - 8c
     - ..m
   * - 108
     - I4cm
     - 4b
     - 2.mm
   * - 108
     - I4cm
     - 4a
     - 4..
   * - 109
     - I41md
     - 16c
     - 1
   * - 109
     - I41md
     - 8b
     - .m.
   * - 109
     - I41md
     - 4a
     - 2mm.
   * - 110
     - I41cd
     - 16b
     - 1
   * - 110
     - I41cd
     - 8a
     - 2..
   * - 111
     - P-42m
     - 8o
     - 1
   * - 111
     - P-42m
     - 4n
     - ..m
   * - 111
     - P-42m
     - 4m
     - 2..
   * - 111
     - P-42m
     - 4l
     - .2.
   * - 111
     - P-42m
     - 4k
     - .2.
   * - 111
     - P-42m
     - 4j
     - .2.
   * - 111
     - P-42m
     - 4i
     - .2.
   * - 111
     - P-42m
     - 2h
     - 2.mm
   * - 111
     - P-42m
     - 2g
     - 2.mm
   * - 111
     - P-42m
     - 2f
     - 222.\
   * - 111
     - P-42m
     - 2e
     - 222.\
   * - 111
     - P-42m
     - 1d
     - -42m
   * - 111
     - P-42m
     - 1c
     - -42m
   * - 111
     - P-42m
     - 1b
     - -42m
   * - 111
     - P-42m
     - 1a
     - -42m
   * - 112
     - P-42c
     - 8n
     - 1
   * - 112
     - P-42c
     - 4m
     - 2..
   * - 112
     - P-42c
     - 4l
     - 2..
   * - 112
     - P-42c
     - 4k
     - 2..
   * - 112
     - P-42c
     - 4j
     - .2.
   * - 112
     - P-42c
     - 4i
     - .2.
   * - 112
     - P-42c
     - 4h
     - .2.
   * - 112
     - P-42c
     - 4g
     - .2.
   * - 112
     - P-42c
     - 2f
     - -4..
   * - 112
     - P-42c
     - 2e
     - -4..
   * - 112
     - P-42c
     - 2d
     - 222.\
   * - 112
     - P-42c
     - 2c
     - 222.\
   * - 112
     - P-42c
     - 2b
     - 222.\
   * - 112
     - P-42c
     - 2a
     - 222.\
   * - 113
     - P-421m
     - 8f
     - 1
   * - 113
     - P-421m
     - 4e
     - ..m
   * - 113
     - P-421m
     - 4d
     - 2..
   * - 113
     - P-421m
     - 2c
     - 2.mm
   * - 113
     - P-421m
     - 2b
     - -4..
   * - 113
     - P-421m
     - 2a
     - -4..
   * - 114
     - P-421c
     - 8e
     - 1
   * - 114
     - P-421c
     - 4d
     - 2..
   * - 114
     - P-421c
     - 4c
     - 2..
   * - 114
     - P-421c
     - 2b
     - -4..
   * - 114
     - P-421c
     - 2a
     - -4..
   * - 115
     - P-4m2
     - 8l
     - 1
   * - 115
     - P-4m2
     - 4k
     - .m.
   * - 115
     - P-4m2
     - 4j
     - .m.
   * - 115
     - P-4m2
     - 4i
     - ..2
   * - 115
     - P-4m2
     - 4h
     - ..2
   * - 115
     - P-4m2
     - 2g
     - 2mm.
   * - 115
     - P-4m2
     - 2f
     - 2mm.
   * - 115
     - P-4m2
     - 2e
     - 2mm.
   * - 115
     - P-4m2
     - 1d
     - -4m2
   * - 115
     - P-4m2
     - 1c
     - -4m2
   * - 115
     - P-4m2
     - 1b
     - -4m2
   * - 115
     - P-4m2
     - 1a
     - -4m2
   * - 116
     - P-4c2
     - 8j
     - 1
   * - 116
     - P-4c2
     - 4i
     - 2..
   * - 116
     - P-4c2
     - 4h
     - 2..
   * - 116
     - P-4c2
     - 4g
     - 2..
   * - 116
     - P-4c2
     - 4f
     - ..2
   * - 116
     - P-4c2
     - 4e
     - ..2
   * - 116
     - P-4c2
     - 2d
     - -4..
   * - 116
     - P-4c2
     - 2c
     - -4..
   * - 116
     - P-4c2
     - 2b
     - 2.22
   * - 116
     - P-4c2
     - 2a
     - 2.22
   * - 117
     - P-4b2
     - 8i
     - 1
   * - 117
     - P-4b2
     - 4h
     - ..2
   * - 117
     - P-4b2
     - 4g
     - ..2
   * - 117
     - P-4b2
     - 4f
     - 2..
   * - 117
     - P-4b2
     - 4e
     - 2..
   * - 117
     - P-4b2
     - 2d
     - 2.22
   * - 117
     - P-4b2
     - 2c
     - 2.22
   * - 117
     - P-4b2
     - 2b
     - -4..
   * - 117
     - P-4b2
     - 2a
     - -4..
   * - 118
     - P-4n2
     - 8i
     - 1
   * - 118
     - P-4n2
     - 4h
     - 2..
   * - 118
     - P-4n2
     - 4g
     - ..2
   * - 118
     - P-4n2
     - 4f
     - ..2
   * - 118
     - P-4n2
     - 4e
     - 2..
   * - 118
     - P-4n2
     - 2d
     - 2.22
   * - 118
     - P-4n2
     - 2c
     - 2.22
   * - 118
     - P-4n2
     - 2b
     - -4..
   * - 118
     - P-4n2
     - 2a
     - -4..
   * - 119
     - I-4m2
     - 16j
     - 1
   * - 119
     - I-4m2
     - 8i
     - .m.
   * - 119
     - I-4m2
     - 8h
     - ..2
   * - 119
     - I-4m2
     - 8g
     - ..2
   * - 119
     - I-4m2
     - 4f
     - 2mm.
   * - 119
     - I-4m2
     - 4e
     - 2mm.
   * - 119
     - I-4m2
     - 2d
     - -4m2
   * - 119
     - I-4m2
     - 2c
     - -4m2
   * - 119
     - I-4m2
     - 2b
     - -4m2
   * - 119
     - I-4m2
     - 2a
     - -4m2
   * - 120
     - I-4c2
     - 16i
     - 1
   * - 120
     - I-4c2
     - 8h
     - ..2
   * - 120
     - I-4c2
     - 8g
     - 2..
   * - 120
     - I-4c2
     - 8f
     - 2..
   * - 120
     - I-4c2
     - 8e
     - ..2
   * - 120
     - I-4c2
     - 4d
     - 2.22
   * - 120
     - I-4c2
     - 4c
     - -4..
   * - 120
     - I-4c2
     - 4b
     - -4..
   * - 120
     - I-4c2
     - 4a
     - 2.22
   * - 121
     - I-42m
     - 16j
     - 1
   * - 121
     - I-42m
     - 8i
     - ..m
   * - 121
     - I-42m
     - 8h
     - 2..
   * - 121
     - I-42m
     - 8g
     - .2.
   * - 121
     - I-42m
     - 8f
     - .2.
   * - 121
     - I-42m
     - 4e
     - 2.mm
   * - 121
     - I-42m
     - 4d
     - -4..
   * - 121
     - I-42m
     - 4c
     - 222.\
   * - 121
     - I-42m
     - 2b
     - -42m
   * - 121
     - I-42m
     - 2a
     - -42m
   * - 122
     - I-42d
     - 16e
     - 1
   * - 122
     - I-42d
     - 8d
     - .2.
   * - 122
     - I-42d
     - 8c
     - 2..
   * - 122
     - I-42d
     - 4b
     - -4..
   * - 122
     - I-42d
     - 4a
     - -4..
   * - 123
     - P4/mmm
     - 16u
     - 1
   * - 123
     - P4/mmm
     - 8t
     - .m.
   * - 123
     - P4/mmm
     - 8s
     - .m.
   * - 123
     - P4/mmm
     - 8r
     - ..m
   * - 123
     - P4/mmm
     - 8q
     - m..
   * - 123
     - P4/mmm
     - 8p
     - m..
   * - 123
     - P4/mmm
     - 4o
     - m2m.
   * - 123
     - P4/mmm
     - 4n
     - m2m.
   * - 123
     - P4/mmm
     - 4m
     - m2m.
   * - 123
     - P4/mmm
     - 4l
     - m2m.
   * - 123
     - P4/mmm
     - 4k
     - m.2m
   * - 123
     - P4/mmm
     - 4j
     - m.2m
   * - 123
     - P4/mmm
     - 4i
     - 2mm.
   * - 123
     - P4/mmm
     - 2h
     - 4mm
   * - 123
     - P4/mmm
     - 2g
     - 4mm
   * - 123
     - P4/mmm
     - 2f
     - mmm.\
   * - 123
     - P4/mmm
     - 2e
     - mmm.\
   * - 123
     - P4/mmm
     - 1d
     - 4/mmm
   * - 123
     - P4/mmm
     - 1c
     - 4/mmm
   * - 123
     - P4/mmm
     - 1b
     - 4/mmm
   * - 123
     - P4/mmm
     - 1a
     - 4/mmm
   * - 124
     - P4/mcc
     - 16n
     - 1
   * - 124
     - P4/mcc
     - 8m
     - m..
   * - 124
     - P4/mcc
     - 8l
     - .2.
   * - 124
     - P4/mcc
     - 8k
     - .2.
   * - 124
     - P4/mcc
     - 8j
     - ..2
   * - 124
     - P4/mcc
     - 8i
     - 2..
   * - 124
     - P4/mcc
     - 4h
     - 4..
   * - 124
     - P4/mcc
     - 4g
     - 4..
   * - 124
     - P4/mcc
     - 4f
     - 222.\
   * - 124
     - P4/mcc
     - 4e
     - 2/m..
   * - 124
     - P4/mcc
     - 2d
     - 4/m..
   * - 124
     - P4/mcc
     - 2c
     - 422
   * - 124
     - P4/mcc
     - 2b
     - 4/m..
   * - 124
     - P4/mcc
     - 2a
     - 422
   * - 125
     - P4/nbm
     - 16n
     - 1
   * - 125
     - P4/nbm
     - 8m
     - ..m
   * - 125
     - P4/nbm
     - 8l
     - .2.
   * - 125
     - P4/nbm
     - 8k
     - .2.
   * - 125
     - P4/nbm
     - 8j
     - ..2
   * - 125
     - P4/nbm
     - 8i
     - ..2
   * - 125
     - P4/nbm
     - 4h
     - 2.mm
   * - 125
     - P4/nbm
     - 4g
     - 4..
   * - 125
     - P4/nbm
     - 4f
     - ..2/m
   * - 125
     - P4/nbm
     - 4e
     - ..2/m
   * - 125
     - P4/nbm
     - 2d
     - -42m
   * - 125
     - P4/nbm
     - 2c
     - -42m
   * - 125
     - P4/nbm
     - 2b
     - 422
   * - 125
     - P4/nbm
     - 2a
     - 422
   * - 126
     - P4/nnc
     - 16k
     - 1
   * - 126
     - P4/nnc
     - 8j
     - .2.
   * - 126
     - P4/nnc
     - 8i
     - .2.
   * - 126
     - P4/nnc
     - 8h
     - ..2
   * - 126
     - P4/nnc
     - 8g
     - 2..
   * - 126
     - P4/nnc
     - 8f
     - -1
   * - 126
     - P4/nnc
     - 4e
     - 4..
   * - 126
     - P4/nnc
     - 4d
     - -4..
   * - 126
     - P4/nnc
     - 4c
     - 222.\
   * - 126
     - P4/nnc
     - 2b
     - 422
   * - 126
     - P4/nnc
     - 2a
     - 422
   * - 127
     - P4/mbm
     - 16l
     - 1
   * - 127
     - P4/mbm
     - 8k
     - ..m
   * - 127
     - P4/mbm
     - 8j
     - m..
   * - 127
     - P4/mbm
     - 8i
     - m..
   * - 127
     - P4/mbm
     - 4h
     - m.2m
   * - 127
     - P4/mbm
     - 4g
     - m.2m
   * - 127
     - P4/mbm
     - 4f
     - 2.mm
   * - 127
     - P4/mbm
     - 4e
     - 4..
   * - 127
     - P4/mbm
     - 2d
     - m.mm
   * - 127
     - P4/mbm
     - 2c
     - m.mm
   * - 127
     - P4/mbm
     - 2b
     - 4/m..
   * - 127
     - P4/mbm
     - 2a
     - 4/m..
   * - 128
     - P4/mnc
     - 16i
     - 1
   * - 128
     - P4/mnc
     - 8h
     - m..
   * - 128
     - P4/mnc
     - 8g
     - ..2
   * - 128
     - P4/mnc
     - 8f
     - 2..
   * - 128
     - P4/mnc
     - 4e
     - 4..
   * - 128
     - P4/mnc
     - 4d
     - 2.22
   * - 128
     - P4/mnc
     - 4c
     - 2/m..
   * - 128
     - P4/mnc
     - 2b
     - 4/m..
   * - 128
     - P4/mnc
     - 2a
     - 4/m..
   * - 129
     - P4/nmm
     - 16k
     - 1
   * - 129
     - P4/nmm
     - 8j
     - ..m
   * - 129
     - P4/nmm
     - 8i
     - .m.
   * - 129
     - P4/nmm
     - 8h
     - ..2
   * - 129
     - P4/nmm
     - 8g
     - ..2
   * - 129
     - P4/nmm
     - 4f
     - 2mm.
   * - 129
     - P4/nmm
     - 4e
     - ..2/m
   * - 129
     - P4/nmm
     - 4d
     - ..2/m
   * - 129
     - P4/nmm
     - 2c
     - 4mm
   * - 129
     - P4/nmm
     - 2b
     - -4m2
   * - 129
     - P4/nmm
     - 2a
     - -4m2
   * - 130
     - P4/ncc
     - 16g
     - 1
   * - 130
     - P4/ncc
     - 8f
     - ..2
   * - 130
     - P4/ncc
     - 8e
     - 2..
   * - 130
     - P4/ncc
     - 8d
     - -1
   * - 130
     - P4/ncc
     - 4c
     - 4..
   * - 130
     - P4/ncc
     - 4b
     - -4..
   * - 130
     - P4/ncc
     - 4a
     - 2.22
   * - 131
     - P42/mmc
     - 16r
     - 1
   * - 131
     - P42/mmc
     - 8q
     - m..
   * - 131
     - P42/mmc
     - 8p
     - .m.
   * - 131
     - P42/mmc
     - 8o
     - .m.
   * - 131
     - P42/mmc
     - 8n
     - ..2
   * - 131
     - P42/mmc
     - 4m
     - m2m.
   * - 131
     - P42/mmc
     - 4l
     - m2m.
   * - 131
     - P42/mmc
     - 4k
     - m2m.
   * - 131
     - P42/mmc
     - 4j
     - m2m.
   * - 131
     - P42/mmc
     - 4i
     - 2mm.
   * - 131
     - P42/mmc
     - 4h
     - 2mm.
   * - 131
     - P42/mmc
     - 4g
     - 2mm.
   * - 131
     - P42/mmc
     - 2f
     - -4m2
   * - 131
     - P42/mmc
     - 2e
     - -4m2
   * - 131
     - P42/mmc
     - 2d
     - mmm.\
   * - 131
     - P42/mmc
     - 2c
     - mmm.\
   * - 131
     - P42/mmc
     - 2b
     - mmm.\
   * - 131
     - P42/mmc
     - 2a
     - mmm.\
   * - 132
     - P42/mcm
     - 16p
     - 1
   * - 132
     - P42/mcm
     - 8o
     - ..m
   * - 132
     - P42/mcm
     - 8n
     - m..
   * - 132
     - P42/mcm
     - 8m
     - .2.
   * - 132
     - P42/mcm
     - 8l
     - .2.
   * - 132
     - P42/mcm
     - 8k
     - 2..
   * - 132
     - P42/mcm
     - 4j
     - m.2m
   * - 132
     - P42/mcm
     - 4i
     - m.2m
   * - 132
     - P42/mcm
     - 4h
     - 2.mm
   * - 132
     - P42/mcm
     - 4g
     - 2.mm
   * - 132
     - P42/mcm
     - 4f
     - 2/m..
   * - 132
     - P42/mcm
     - 4e
     - 222.\
   * - 132
     - P42/mcm
     - 2d
     - -42m
   * - 132
     - P42/mcm
     - 2c
     - m.mm
   * - 132
     - P42/mcm
     - 2b
     - -42m
   * - 132
     - P42/mcm
     - 2a
     - m.mm
   * - 133
     - P42/nbc
     - 16k
     - 1
   * - 133
     - P42/nbc
     - 8j
     - ..2
   * - 133
     - P42/nbc
     - 8i
     - .2.
   * - 133
     - P42/nbc
     - 8h
     - .2.
   * - 133
     - P42/nbc
     - 8g
     - 2..
   * - 133
     - P42/nbc
     - 8f
     - 2..
   * - 133
     - P42/nbc
     - 8e
     - -1
   * - 133
     - P42/nbc
     - 4d
     - -4..
   * - 133
     - P42/nbc
     - 4c
     - 2.22
   * - 133
     - P42/nbc
     - 4b
     - 222.\
   * - 133
     - P42/nbc
     - 4a
     - 222.\
   * - 134
     - P42/nnm
     - 16n
     - 1
   * - 134
     - P42/nnm
     - 8m
     - ..m
   * - 134
     - P42/nnm
     - 8l
     - ..2
   * - 134
     - P42/nnm
     - 8k
     - ..2
   * - 134
     - P42/nnm
     - 8j
     - .2.
   * - 134
     - P42/nnm
     - 8i
     - .2.
   * - 134
     - P42/nnm
     - 8h
     - 2..
   * - 134
     - P42/nnm
     - 4g
     - 2.mm
   * - 134
     - P42/nnm
     - 4f
     - ..2/m
   * - 134
     - P42/nnm
     - 4e
     - ..2/m
   * - 134
     - P42/nnm
     - 4d
     - 2.22
   * - 134
     - P42/nnm
     - 4c
     - 222.\
   * - 134
     - P42/nnm
     - 2b
     - -42m
   * - 134
     - P42/nnm
     - 2a
     - -42m
   * - 135
     - P42/mbc
     - 16i
     - 1
   * - 135
     - P42/mbc
     - 8h
     - m..
   * - 135
     - P42/mbc
     - 8g
     - ..2
   * - 135
     - P42/mbc
     - 8f
     - 2..
   * - 135
     - P42/mbc
     - 8e
     - 2..
   * - 135
     - P42/mbc
     - 4d
     - 2.22
   * - 135
     - P42/mbc
     - 4c
     - 2/m..
   * - 135
     - P42/mbc
     - 4b
     - -4..
   * - 135
     - P42/mbc
     - 4a
     - 2/m..
   * - 136
     - P42/mnm
     - 16k
     - 1
   * - 136
     - P42/mnm
     - 8j
     - ..m
   * - 136
     - P42/mnm
     - 8i
     - m..
   * - 136
     - P42/mnm
     - 8h
     - 2..
   * - 136
     - P42/mnm
     - 4g
     - m.2m
   * - 136
     - P42/mnm
     - 4f
     - m.2m
   * - 136
     - P42/mnm
     - 4e
     - 2.mm
   * - 136
     - P42/mnm
     - 4d
     - -4..
   * - 136
     - P42/mnm
     - 4c
     - 2/m..
   * - 136
     - P42/mnm
     - 2b
     - m.mm
   * - 136
     - P42/mnm
     - 2a
     - m.mm
   * - 137
     - P42/nmc
     - 16h
     - 1
   * - 137
     - P42/nmc
     - 8g
     - .m.
   * - 137
     - P42/nmc
     - 8f
     - ..2
   * - 137
     - P42/nmc
     - 8e
     - -1
   * - 137
     - P42/nmc
     - 4d
     - 2mm.
   * - 137
     - P42/nmc
     - 4c
     - 2mm.
   * - 137
     - P42/nmc
     - 2b
     - -4m2
   * - 137
     - P42/nmc
     - 2a
     - -4m2
   * - 138
     - P42/ncm
     - 16j
     - 1
   * - 138
     - P42/ncm
     - 8i
     - ..m
   * - 138
     - P42/ncm
     - 8h
     - ..2
   * - 138
     - P42/ncm
     - 8g
     - ..2
   * - 138
     - P42/ncm
     - 8f
     - 2..
   * - 138
     - P42/ncm
     - 4e
     - 2.mm
   * - 138
     - P42/ncm
     - 4d
     - ..2/m
   * - 138
     - P42/ncm
     - 4c
     - ..2/m
   * - 138
     - P42/ncm
     - 4b
     - -4..
   * - 138
     - P42/ncm
     - 4a
     - 2.22
   * - 139
     - I4/mmm
     - 32o
     - 1
   * - 139
     - I4/mmm
     - 16n
     - .m.
   * - 139
     - I4/mmm
     - 16m
     - ..m
   * - 139
     - I4/mmm
     - 16l
     - m..
   * - 139
     - I4/mmm
     - 16k
     - ..2
   * - 139
     - I4/mmm
     - 8j
     - m2m.
   * - 139
     - I4/mmm
     - 8i
     - m2m.
   * - 139
     - I4/mmm
     - 8h
     - m.2m
   * - 139
     - I4/mmm
     - 8g
     - 2mm.
   * - 139
     - I4/mmm
     - 8f
     - ..2/m
   * - 139
     - I4/mmm
     - 4e
     - 4mm
   * - 139
     - I4/mmm
     - 4d
     - -4m2
   * - 139
     - I4/mmm
     - 4c
     - mmm.\
   * - 139
     - I4/mmm
     - 2b
     - 4/mmm
   * - 139
     - I4/mmm
     - 2a
     - 4/mmm
   * - 140
     - I4/mcm
     - 32m
     - 1
   * - 140
     - I4/mcm
     - 16l
     - ..m
   * - 140
     - I4/mcm
     - 16k
     - m..
   * - 140
     - I4/mcm
     - 16j
     - .2.
   * - 140
     - I4/mcm
     - 16i
     - ..2
   * - 140
     - I4/mcm
     - 8h
     - m.2m
   * - 140
     - I4/mcm
     - 8g
     - 2.mm
   * - 140
     - I4/mcm
     - 8f
     - 4..
   * - 140
     - I4/mcm
     - 8e
     - ..2/m
   * - 140
     - I4/mcm
     - 4d
     - m.mm
   * - 140
     - I4/mcm
     - 4c
     - 4/m..
   * - 140
     - I4/mcm
     - 4b
     - -42m
   * - 140
     - I4/mcm
     - 4a
     - 422
   * - 141
     - I41/amd
     - 32i
     - 1
   * - 141
     - I41/amd
     - 16h
     - .m.
   * - 141
     - I41/amd
     - 16g
     - ..2
   * - 141
     - I41/amd
     - 16f
     - .2.
   * - 141
     - I41/amd
     - 8e
     - 2mm.
   * - 141
     - I41/amd
     - 8d
     - .2/m.
   * - 141
     - I41/amd
     - 8c
     - .2/m.
   * - 141
     - I41/amd
     - 4b
     - -4m2
   * - 141
     - I41/amd
     - 4a
     - -4m2
   * - 142
     - I41/acd
     - 32g
     - 1
   * - 142
     - I41/acd
     - 16f
     - ..2
   * - 142
     - I41/acd
     - 16e
     - .2.
   * - 142
     - I41/acd
     - 16d
     - 2..
   * - 142
     - I41/acd
     - 16c
     - -1
   * - 142
     - I41/acd
     - 8b
     - 2.22
   * - 142
     - I41/acd
     - 8a
     - -4..
   * - 143
     - P3
     - 3d
     - 1
   * - 143
     - P3
     - 1c
     - 3..
   * - 143
     - P3
     - 1b
     - 3..
   * - 143
     - P3
     - 1a
     - 3..
   * - 144
     - P31
     - 3a
     - 1
   * - 145
     - P32
     - 3a
     - 1
   * - 146
     - R3
     - 9b
     - 1
   * - 146
     - R3
     - 3a
     - 3.\
   * - 147
     - P-3
     - 6g
     - 1
   * - 147
     - P-3
     - 3f
     - -1
   * - 147
     - P-3
     - 3e
     - -1
   * - 147
     - P-3
     - 2d
     - 3..
   * - 147
     - P-3
     - 2c
     - 3..
   * - 147
     - P-3
     - 1b
     - -3..
   * - 147
     - P-3
     - 1a
     - -3..
   * - 148
     - R-3
     - 18f
     - 1
   * - 148
     - R-3
     - 9e
     - -1
   * - 148
     - R-3
     - 9d
     - -1
   * - 148
     - R-3
     - 6c
     - 3.\
   * - 148
     - R-3
     - 3b
     - -3.
   * - 148
     - R-3
     - 3a
     - -3.
   * - 149
     - P312
     - 6l
     - 1
   * - 149
     - P312
     - 3k
     - ..2
   * - 149
     - P312
     - 3j
     - ..2
   * - 149
     - P312
     - 2i
     - 3..
   * - 149
     - P312
     - 2h
     - 3..
   * - 149
     - P312
     - 2g
     - 3..
   * - 149
     - P312
     - 1f
     - 322
   * - 149
     - P312
     - 1e
     - 322
   * - 149
     - P312
     - 1d
     - 322
   * - 149
     - P312
     - 1c
     - 322
   * - 149
     - P312
     - 1b
     - 322
   * - 149
     - P312
     - 1a
     - 322
   * - 150
     - P321
     - 6g
     - 1
   * - 150
     - P321
     - 3f
     - .2.
   * - 150
     - P321
     - 3e
     - .2.
   * - 150
     - P321
     - 2d
     - 3..
   * - 150
     - P321
     - 2c
     - 3..
   * - 150
     - P321
     - 1b
     - 32.\
   * - 150
     - P321
     - 1a
     - 32.\
   * - 151
     - P3112
     - 6c
     - 1
   * - 151
     - P3112
     - 3b
     - ..2
   * - 151
     - P3112
     - 3a
     - ..2
   * - 152
     - P3121
     - 6c
     - 1
   * - 152
     - P3121
     - 3b
     - .2.
   * - 152
     - P3121
     - 3a
     - .2.
   * - 153
     - P3212
     - 6c
     - 1
   * - 153
     - P3212
     - 3b
     - ..2
   * - 153
     - P3212
     - 3a
     - ..2
   * - 154
     - P3221
     - 6c
     - 1
   * - 154
     - P3221
     - 3b
     - .2.
   * - 154
     - P3221
     - 3a
     - .2.
   * - 155
     - R32
     - 18f
     - 1
   * - 155
     - R32
     - 9e
     - .2
   * - 155
     - R32
     - 9d
     - .2
   * - 155
     - R32
     - 6c
     - 3.\
   * - 155
     - R32
     - 3b
     - 32
   * - 155
     - R32
     - 3a
     - 32
   * - 156
     - P3m1
     - 6e
     - 1
   * - 156
     - P3m1
     - 3d
     - .m.
   * - 156
     - P3m1
     - 1c
     - 3m.
   * - 156
     - P3m1
     - 1b
     - 3m.
   * - 156
     - P3m1
     - 1a
     - 3m.
   * - 157
     - P31m
     - 6d
     - 1
   * - 157
     - P31m
     - 3c
     - .m.
   * - 157
     - P31m
     - 2b
     - 3..
   * - 157
     - P31m
     - 1a
     - 3mm
   * - 158
     - P3c1
     - 6d
     - 1
   * - 158
     - P3c1
     - 2c
     - 3..
   * - 158
     - P3c1
     - 2b
     - 3..
   * - 158
     - P3c1
     - 2a
     - 3..
   * - 159
     - P31c
     - 6c
     - 1
   * - 159
     - P31c
     - 2b
     - 3..
   * - 159
     - P31c
     - 2a
     - 3..
   * - 160
     - R3m
     - 18c
     - 1
   * - 160
     - R3m
     - 9b
     - .m
   * - 160
     - R3m
     - 3a
     - 3m
   * - 161
     - R3c
     - 18b
     - 1
   * - 161
     - R3c
     - 6a
     - 3.\
   * - 162
     - P-31m
     - 12l
     - 1
   * - 162
     - P-31m
     - 6k
     - .m.
   * - 162
     - P-31m
     - 6j
     - ..2
   * - 162
     - P-31m
     - 6i
     - ..2
   * - 162
     - P-31m
     - 4h
     - 3..
   * - 162
     - P-31m
     - 3g
     - .2/m.
   * - 162
     - P-31m
     - 3f
     - .2/m.
   * - 162
     - P-31m
     - 2e
     - 3mm
   * - 162
     - P-31m
     - 2d
     - 322
   * - 162
     - P-31m
     - 2c
     - 322
   * - 162
     - P-31m
     - 1b
     - -3m2/m
   * - 162
     - P-31m
     - 1a
     - -3m2/m
   * - 163
     - P-31c
     - 12i
     - 1
   * - 163
     - P-31c
     - 6h
     - ..2
   * - 163
     - P-31c
     - 6g
     - -1
   * - 163
     - P-31c
     - 4f
     - 3..
   * - 163
     - P-31c
     - 4e
     - 3..
   * - 163
     - P-31c
     - 2d
     - 322
   * - 163
     - P-31c
     - 2c
     - 322
   * - 163
     - P-31c
     - 2b
     - -3..
   * - 163
     - P-31c
     - 2a
     - 322
   * - 164
     - P-3m1
     - 12j
     - 1
   * - 164
     - P-3m1
     - 6i
     - .m.
   * - 164
     - P-3m1
     - 6h
     - .2.
   * - 164
     - P-3m1
     - 6g
     - .2.
   * - 164
     - P-3m1
     - 3f
     - .2/m.
   * - 164
     - P-3m1
     - 3e
     - .2/m.
   * - 164
     - P-3m1
     - 2d
     - 3m.
   * - 164
     - P-3m1
     - 2c
     - 3m.
   * - 164
     - P-3m1
     - 1b
     - -3m.
   * - 164
     - P-3m1
     - 1a
     - -3m.
   * - 165
     - P-3c1
     - 12g
     - 1
   * - 165
     - P-3c1
     - 6f
     - .2.
   * - 165
     - P-3c1
     - 6e
     - -1
   * - 165
     - P-3c1
     - 4d
     - 3..
   * - 165
     - P-3c1
     - 4c
     - 3..
   * - 165
     - P-3c1
     - 2b
     - -3..
   * - 165
     - P-3c1
     - 2a
     - 32.\
   * - 166
     - R-3m
     - 36i
     - 1
   * - 166
     - R-3m
     - 18h
     - .m
   * - 166
     - R-3m
     - 18g
     - .2
   * - 166
     - R-3m
     - 18f
     - .2
   * - 166
     - R-3m
     - 9e
     - .2/m
   * - 166
     - R-3m
     - 9d
     - .2/m
   * - 166
     - R-3m
     - 6c
     - 3m
   * - 166
     - R-3m
     - 3b
     - -3m
   * - 166
     - R-3m
     - 3a
     - -3m
   * - 167
     - R-3c
     - 36f
     - 1
   * - 167
     - R-3c
     - 18e
     - .2
   * - 167
     - R-3c
     - 18d
     - -1
   * - 167
     - R-3c
     - 12c
     - 3.\
   * - 167
     - R-3c
     - 6b
     - -3.
   * - 167
     - R-3c
     - 6a
     - 32
   * - 168
     - P6
     - 6d
     - 1
   * - 168
     - P6
     - 3c
     - 2..
   * - 168
     - P6
     - 2b
     - 3..
   * - 168
     - P6
     - 1a
     - 6..
   * - 169
     - P61
     - 6a
     - 1
   * - 170
     - P65
     - 6a
     - 1
   * - 171
     - P62
     - 6c
     - 1
   * - 171
     - P62
     - 3b
     - 2..
   * - 171
     - P62
     - 3a
     - 2..
   * - 172
     - P64
     - 6c
     - 1
   * - 172
     - P64
     - 3b
     - 2..
   * - 172
     - P64
     - 3a
     - 2..
   * - 173
     - P63
     - 6c
     - 1
   * - 173
     - P63
     - 2b
     - 3..
   * - 173
     - P63
     - 2a
     - 3..
   * - 174
     - P-6
     - 6l
     - 1
   * - 174
     - P-6
     - 3k
     - m..
   * - 174
     - P-6
     - 3j
     - m..
   * - 174
     - P-6
     - 2i
     - 3..
   * - 174
     - P-6
     - 2h
     - 3..
   * - 174
     - P-6
     - 2g
     - 3..
   * - 174
     - P-6
     - 1f
     - -6..
   * - 174
     - P-6
     - 1e
     - -6..
   * - 174
     - P-6
     - 1d
     - -6..
   * - 174
     - P-6
     - 1c
     - -6..
   * - 174
     - P-6
     - 1b
     - -6..
   * - 174
     - P-6
     - 1a
     - -6..
   * - 175
     - P6/m
     - 12l
     - 1
   * - 175
     - P6/m
     - 6k
     - m..
   * - 175
     - P6/m
     - 6j
     - m..
   * - 175
     - P6/m
     - 6i
     - 2..
   * - 175
     - P6/m
     - 4h
     - 3..
   * - 175
     - P6/m
     - 3g
     - 2/m..
   * - 175
     - P6/m
     - 3f
     - 2/m..
   * - 175
     - P6/m
     - 2e
     - 6..
   * - 175
     - P6/m
     - 2d
     - -6..
   * - 175
     - P6/m
     - 2c
     - -6..
   * - 175
     - P6/m
     - 1b
     - 6/m..
   * - 175
     - P6/m
     - 1a
     - 6/m..
   * - 176
     - P63/m
     - 12i
     - 1
   * - 176
     - P63/m
     - 6h
     - m..
   * - 176
     - P63/m
     - 6g
     - -1
   * - 176
     - P63/m
     - 4f
     - 3..
   * - 176
     - P63/m
     - 4e
     - 3..
   * - 176
     - P63/m
     - 2d
     - -6..
   * - 176
     - P63/m
     - 2c
     - -6..
   * - 176
     - P63/m
     - 2b
     - -3..
   * - 176
     - P63/m
     - 2a
     - -6..
   * - 177
     - P622
     - 12n
     - 1
   * - 177
     - P622
     - 6m
     - ..2
   * - 177
     - P622
     - 6l
     - ..2
   * - 177
     - P622
     - 6k
     - .2.
   * - 177
     - P622
     - 6j
     - .2.
   * - 177
     - P622
     - 6i
     - 2..
   * - 177
     - P622
     - 4h
     - 3..
   * - 177
     - P622
     - 3g
     - 22.\
   * - 177
     - P622
     - 3f
     - 22.\
   * - 177
     - P622
     - 2e
     - 6..
   * - 177
     - P622
     - 2d
     - 322
   * - 177
     - P622
     - 2c
     - 322
   * - 177
     - P622
     - 1b
     - 622
   * - 177
     - P622
     - 1a
     - 622
   * - 178
     - P6122
     - 12c
     - 1
   * - 178
     - P6122
     - 6b
     - .2.
   * - 178
     - P6122
     - 6a
     - .2.
   * - 179
     - P6522
     - 12c
     - 1
   * - 179
     - P6522
     - 6b
     - .2.
   * - 179
     - P6522
     - 6a
     - .2.
   * - 180
     - P6222
     - 12k
     - 1
   * - 180
     - P6222
     - 6j
     - .2.
   * - 180
     - P6222
     - 6i
     - .2.
   * - 180
     - P6222
     - 6h
     - .2.
   * - 180
     - P6222
     - 6g
     - .2.
   * - 180
     - P6222
     - 6f
     - 2..
   * - 180
     - P6222
     - 6e
     - 2..
   * - 180
     - P6222
     - 3d
     - 22.\
   * - 180
     - P6222
     - 3c
     - 22.\
   * - 180
     - P6222
     - 3b
     - 22.\
   * - 180
     - P6222
     - 3a
     - 22.\
   * - 181
     - P6422
     - 12k
     - 1
   * - 181
     - P6422
     - 6j
     - .2.
   * - 181
     - P6422
     - 6i
     - .2.
   * - 181
     - P6422
     - 6h
     - .2.
   * - 181
     - P6422
     - 6g
     - .2.
   * - 181
     - P6422
     - 6f
     - 2..
   * - 181
     - P6422
     - 6e
     - 2..
   * - 181
     - P6422
     - 3d
     - 22.\
   * - 181
     - P6422
     - 3c
     - 22.\
   * - 181
     - P6422
     - 3b
     - 22.\
   * - 181
     - P6422
     - 3a
     - 22.\
   * - 182
     - P6322
     - 12i
     - 1
   * - 182
     - P6322
     - 6h
     - .2.
   * - 182
     - P6322
     - 6g
     - .2.
   * - 182
     - P6322
     - 4f
     - 3..
   * - 182
     - P6322
     - 4e
     - 3..
   * - 182
     - P6322
     - 2d
     - 322
   * - 182
     - P6322
     - 2c
     - 322
   * - 182
     - P6322
     - 2b
     - 322
   * - 182
     - P6322
     - 2a
     - 32.\
   * - 183
     - P6mm
     - 12f
     - 1
   * - 183
     - P6mm
     - 6e
     - .m.
   * - 183
     - P6mm
     - 6d
     - .m.
   * - 183
     - P6mm
     - 3c
     - 2m.
   * - 183
     - P6mm
     - 2b
     - 3m.
   * - 183
     - P6mm
     - 1a
     - 6mm
   * - 184
     - P6cc
     - 12d
     - 1
   * - 184
     - P6cc
     - 6c
     - 2..
   * - 184
     - P6cc
     - 4b
     - 3..
   * - 184
     - P6cc
     - 2a
     - 6..
   * - 185
     - P63cm
     - 12d
     - 1
   * - 185
     - P63cm
     - 6c
     - .m.
   * - 185
     - P63cm
     - 4b
     - 3..
   * - 185
     - P63cm
     - 2a
     - 3mm
   * - 186
     - P63mc
     - 12d
     - 1
   * - 186
     - P63mc
     - 6c
     - .m.
   * - 186
     - P63mc
     - 2b
     - 3m.
   * - 186
     - P63mc
     - 2a
     - 3m.
   * - 187
     - P-6m2
     - 12o
     - 1
   * - 187
     - P-6m2
     - 6n
     - .m.
   * - 187
     - P-6m2
     - 6m
     - m..
   * - 187
     - P-6m2
     - 6l
     - m..
   * - 187
     - P-6m2
     - 3k
     - mm2
   * - 187
     - P-6m2
     - 3j
     - mm2
   * - 187
     - P-6m2
     - 2i
     - 3m.
   * - 187
     - P-6m2
     - 2h
     - 3m.
   * - 187
     - P-6m2
     - 2g
     - 3m.
   * - 187
     - P-6m2
     - 1f
     - -6m2m2
   * - 187
     - P-6m2
     - 1e
     - -6m2m2
   * - 187
     - P-6m2
     - 1d
     - -6m2m2
   * - 187
     - P-6m2
     - 1c
     - -6m2m2
   * - 187
     - P-6m2
     - 1b
     - -6m2m2
   * - 187
     - P-6m2
     - 1a
     - -6m2m2
   * - 188
     - P-6c2
     - 12l
     - 1
   * - 188
     - P-6c2
     - 6k
     - m..
   * - 188
     - P-6c2
     - 6j
     - ..2
   * - 188
     - P-6c2
     - 4i
     - 3..
   * - 188
     - P-6c2
     - 4h
     - 3..
   * - 188
     - P-6c2
     - 4g
     - 3..
   * - 188
     - P-6c2
     - 2f
     - -6..
   * - 188
     - P-6c2
     - 2e
     - 322
   * - 188
     - P-6c2
     - 2d
     - -6..
   * - 188
     - P-6c2
     - 2c
     - 322
   * - 188
     - P-6c2
     - 2b
     - -6..
   * - 188
     - P-6c2
     - 2a
     - 322
   * - 189
     - P-62m
     - 12l
     - 1
   * - 189
     - P-62m
     - 6k
     - m..
   * - 189
     - P-62m
     - 6j
     - m..
   * - 189
     - P-62m
     - 6i
     - .m.
   * - 189
     - P-62m
     - 4h
     - 3..
   * - 189
     - P-62m
     - 3g
     - m2m.
   * - 189
     - P-62m
     - 3f
     - m2m.
   * - 189
     - P-62m
     - 2e
     - 3mm
   * - 189
     - P-62m
     - 2d
     - -6..
   * - 189
     - P-62m
     - 2c
     - -6..
   * - 189
     - P-62m
     - 1b
     - -6mm2m
   * - 189
     - P-62m
     - 1a
     - -6mm2m
   * - 190
     - P-62c
     - 12i
     - 1
   * - 190
     - P-62c
     - 6h
     - m..
   * - 190
     - P-62c
     - 6g
     - .2.
   * - 190
     - P-62c
     - 4f
     - 3..
   * - 190
     - P-62c
     - 4e
     - 3..
   * - 190
     - P-62c
     - 2d
     - -6..
   * - 190
     - P-62c
     - 2c
     - -6..
   * - 190
     - P-62c
     - 2b
     - -6..
   * - 190
     - P-62c
     - 2a
     - 32.\
   * - 191
     - P6/mmm
     - 24r
     - 1
   * - 191
     - P6/mmm
     - 12q
     - m..
   * - 191
     - P6/mmm
     - 12p
     - m..
   * - 191
     - P6/mmm
     - 12o
     - .m.
   * - 191
     - P6/mmm
     - 12n
     - .m.
   * - 191
     - P6/mmm
     - 6m
     - mm2.
   * - 191
     - P6/mmm
     - 6l
     - mm2.
   * - 191
     - P6/mmm
     - 6k
     - m2m.
   * - 191
     - P6/mmm
     - 6j
     - m2m.
   * - 191
     - P6/mmm
     - 6i
     - 2m.
   * - 191
     - P6/mmm
     - 4h
     - 3m.
   * - 191
     - P6/mmm
     - 3g
     - 2/mm.
   * - 191
     - P6/mmm
     - 3f
     - 2/mm.
   * - 191
     - P6/mmm
     - 2e
     - 6mm
   * - 191
     - P6/mmm
     - 2d
     - -6m2m2
   * - 191
     - P6/mmm
     - 2c
     - -6m2m2
   * - 191
     - P6/mmm
     - 1b
     - 6/mm2/m
   * - 191
     - P6/mmm
     - 1a
     - 6/mm2/m
   * - 192
     - P6/mcc
     - 24m
     - 1
   * - 192
     - P6/mcc
     - 12l
     - m..
   * - 192
     - P6/mcc
     - 12k
     - .2.
   * - 192
     - P6/mcc
     - 12j
     - .2.
   * - 192
     - P6/mcc
     - 12i
     - 2..
   * - 192
     - P6/mcc
     - 8h
     - 3..
   * - 192
     - P6/mcc
     - 6g
     - 2/m..
   * - 192
     - P6/mcc
     - 6f
     - 22.\
   * - 192
     - P6/mcc
     - 4e
     - 6..
   * - 192
     - P6/mcc
     - 4d
     - -6..
   * - 192
     - P6/mcc
     - 4c
     - 322
   * - 192
     - P6/mcc
     - 2b
     - 6/m..
   * - 192
     - P6/mcc
     - 2a
     - 622
   * - 193
     - P63/mcm
     - 24l
     - 1
   * - 193
     - P63/mcm
     - 12k
     - .m.
   * - 193
     - P63/mcm
     - 12j
     - m..
   * - 193
     - P63/mcm
     - 12i
     - .2.
   * - 193
     - P63/mcm
     - 8h
     - 3..
   * - 193
     - P63/mcm
     - 6g
     - m2m.
   * - 193
     - P63/mcm
     - 6f
     - .2/m.
   * - 193
     - P63/mcm
     - 4e
     - 3mm
   * - 193
     - P63/mcm
     - 4d
     - 322
   * - 193
     - P63/mcm
     - 4c
     - -6..
   * - 193
     - P63/mcm
     - 2b
     - -3m2/m
   * - 193
     - P63/mcm
     - 2a
     - -6mm2m
   * - 194
     - P63/mmc
     - 24l
     - 1
   * - 194
     - P63/mmc
     - 12k
     - .m.
   * - 194
     - P63/mmc
     - 12j
     - m..
   * - 194
     - P63/mmc
     - 12i
     - .2.
   * - 194
     - P63/mmc
     - 6h
     - mm2.
   * - 194
     - P63/mmc
     - 6g
     - .2/m.
   * - 194
     - P63/mmc
     - 4f
     - 3m.
   * - 194
     - P63/mmc
     - 4e
     - 3m.
   * - 194
     - P63/mmc
     - 2d
     - -6m2m2
   * - 194
     - P63/mmc
     - 2c
     - -6m2m2
   * - 194
     - P63/mmc
     - 2b
     - -6m2m2
   * - 194
     - P63/mmc
     - 2a
     - -3m.
   * - 195
     - P23
     - 12j
     - 1
   * - 195
     - P23
     - 6i
     - 2..
   * - 195
     - P23
     - 6h
     - 2..
   * - 195
     - P23
     - 6g
     - 2..
   * - 195
     - P23
     - 6f
     - 2..
   * - 195
     - P23
     - 4e
     - .3.
   * - 195
     - P23
     - 3d
     - 222..
   * - 195
     - P23
     - 3c
     - 222..
   * - 195
     - P23
     - 1b
     - 23.\ 
   * - 195
     - P23
     - 1a
     - 23.\
   * - 196
     - F23
     - 48h
     - 1
   * - 196
     - F23
     - 24g
     - 2..
   * - 196
     - F23
     - 24f
     - 2..
   * - 196
     - F23
     - 16e
     - .3.
   * - 196
     - F23
     - 4d
     - 23.\
   * - 196
     - F23
     - 4c
     - 23.\
   * - 196
     - F23
     - 4b
     - 23.\
   * - 196
     - F23
     - 4a
     - 23.\
   * - 197
     - I23
     - 24f
     - 1
   * - 197
     - I23
     - 12e
     - 2..
   * - 197
     - I23
     - 12d
     - 2..
   * - 197
     - I23
     - 8c
     - .3.
   * - 197
     - I23
     - 6b
     - 222..
   * - 197
     - I23
     - 2a
     - 23.\
   * - 198
     - P213
     - 12b
     - 1
   * - 198
     - P213
     - 4a
     - .3.
   * - 199
     - I213
     - 24c
     - 1
   * - 199
     - I213
     - 12b
     - 2..
   * - 199
     - I213
     - 8a
     - .3.
   * - 200
     - Pm-3
     - 24l
     - 1
   * - 200
     - Pm-3
     - 12k
     - m..
   * - 200
     - Pm-3
     - 12j
     - m..
   * - 200
     - Pm-3
     - 8i
     - .3.
   * - 200
     - Pm-3
     - 6h
     - mm2..
   * - 200
     - Pm-3
     - 6g
     - mm2..
   * - 200
     - Pm-3
     - 6f
     - mm2..
   * - 200
     - Pm-3
     - 6e
     - mm2..
   * - 200
     - Pm-3
     - 3d
     - mmm..
   * - 200
     - Pm-3
     - 3c
     - mmm..
   * - 200
     - Pm-3
     - 1b
     - m-3.
   * - 200
     - Pm-3
     - 1a
     - m-3.
   * - 201
     - Pn-3
     - 24h
     - 1
   * - 201
     - Pn-3
     - 12g
     - 2..
   * - 201
     - Pn-3
     - 12f
     - 2..
   * - 201
     - Pn-3
     - 8e
     - .3.
   * - 201
     - Pn-3
     - 6d
     - 222..
   * - 201
     - Pn-3
     - 4c
     - .-3.
   * - 201
     - Pn-3
     - 4b
     - .-3.
   * - 201
     - Pn-3
     - 2a
     - 23.\
   * - 202
     - Fm-3
     - 96i
     - 1
   * - 202
     - Fm-3
     - 48h
     - m..
   * - 202
     - Fm-3
     - 48g
     - 2..
   * - 202
     - Fm-3
     - 32f
     - .3.
   * - 202
     - Fm-3
     - 24e
     - mm2..
   * - 202
     - Fm-3
     - 24d
     - 2/m..
   * - 202
     - Fm-3
     - 8c
     - 23.\
   * - 202
     - Fm-3
     - 4b
     - m-3.
   * - 202
     - Fm-3
     - 4a
     - m-3.
   * - 203
     - Fd-3
     - 96g
     - 1
   * - 203
     - Fd-3
     - 48f
     - 2..
   * - 203
     - Fd-3
     - 32e
     - .3.
   * - 203
     - Fd-3
     - 16d
     - .-3.
   * - 203
     - Fd-3
     - 16c
     - .-3.
   * - 203
     - Fd-3
     - 8b
     - 23.\
   * - 203
     - Fd-3
     - 8a
     - 23.\
   * - 204
     - Im-3
     - 48h
     - 1
   * - 204
     - Im-3
     - 24g
     - m..
   * - 204
     - Im-3
     - 16f
     - .3.
   * - 204
     - Im-3
     - 12e
     - mm2..
   * - 204
     - Im-3
     - 12d
     - mm2..
   * - 204
     - Im-3
     - 8c
     - .-3.
   * - 204
     - Im-3
     - 6b
     - mmm..
   * - 204
     - Im-3
     - 2a
     - m-3.
   * - 205
     - Pa-3
     - 24d
     - 1
   * - 205
     - Pa-3
     - 8c
     - .3.
   * - 205
     - Pa-3
     - 4b
     - .-3.
   * - 205
     - Pa-3
     - 4a
     - .-3.
   * - 206
     - Ia-3
     - 48e
     - 1
   * - 206
     - Ia-3
     - 24d
     - 2..
   * - 206
     - Ia-3
     - 16c
     - .3.
   * - 206
     - Ia-3
     - 8b
     - .-3.
   * - 206
     - Ia-3
     - 8a
     - .-3.
   * - 207
     - P432
     - 24k
     - 1
   * - 207
     - P432
     - 12j
     - ..2
   * - 207
     - P432
     - 12i
     - ..2
   * - 207
     - P432
     - 12h
     - 2..
   * - 207
     - P432
     - 8g
     - .3.
   * - 207
     - P432
     - 6f
     - 4..
   * - 207
     - P432
     - 6e
     - 4..
   * - 207
     - P432
     - 3d
     - 42.2
   * - 207
     - P432
     - 3c
     - 42.2
   * - 207
     - P432
     - 1b
     - 432
   * - 207
     - P432
     - 1a
     - 432
   * - 208
     - P4232
     - 24m
     - 1
   * - 208
     - P4232
     - 12l
     - ..2
   * - 208
     - P4232
     - 12k
     - ..2
   * - 208
     - P4232
     - 12j
     - 2..
   * - 208
     - P4232
     - 12i
     - 2..
   * - 208
     - P4232
     - 12h
     - 2..
   * - 208
     - P4232
     - 8g
     - .3.
   * - 208
     - P4232
     - 6f
     - 2.22
   * - 208
     - P4232
     - 6e
     - 2.22
   * - 208
     - P4232
     - 6d
     - 222..
   * - 208
     - P4232
     - 4c
     - .32
   * - 208
     - P4232
     - 4b
     - .32
   * - 208
     - P4232
     - 2a
     - 23.\
   * - 209
     - F432
     - 96j
     - 1
   * - 209
     - F432
     - 48i
     - 2..
   * - 209
     - F432
     - 48h
     - ..2
   * - 209
     - F432
     - 48g
     - ..2
   * - 209
     - F432
     - 32f
     - .3.
   * - 209
     - F432
     - 24e
     - 4..
   * - 209
     - F432
     - 24d
     - 2.22
   * - 209
     - F432
     - 8c
     - 23.\
   * - 209
     - F432
     - 4b
     - 432
   * - 209
     - F432
     - 4a
     - 432
   * - 210
     - F4132
     - 96h
     - 1
   * - 210
     - F4132
     - 48g
     - ..2
   * - 210
     - F4132
     - 48f
     - 2..
   * - 210
     - F4132
     - 32e
     - .3.
   * - 210
     - F4132
     - 16d
     - .32
   * - 210
     - F4132
     - 16c
     - .32
   * - 210
     - F4132
     - 8b
     - 23.\
   * - 210
     - F4132
     - 8a
     - 23.\
   * - 211
     - I432
     - 48j
     - 1
   * - 211
     - I432
     - 24i
     - ..2
   * - 211
     - I432
     - 24h
     - ..2
   * - 211
     - I432
     - 24g
     - 2..
   * - 211
     - I432
     - 16f
     - .3.
   * - 211
     - I432
     - 12e
     - 4..
   * - 211
     - I432
     - 12d
     - 2.22
   * - 211
     - I432
     - 8c
     - .32
   * - 211
     - I432
     - 6b
     - 42.2
   * - 211
     - I432
     - 2a
     - 432
   * - 212
     - P4332
     - 24e
     - 1
   * - 212
     - P4332
     - 12d
     - ..2
   * - 212
     - P4332
     - 8c
     - .3.
   * - 212
     - P4332
     - 4b
     - .32
   * - 212
     - P4332
     - 4a
     - .32
   * - 213
     - P4132
     - 24e
     - 1
   * - 213
     - P4132
     - 12d
     - ..2
   * - 213
     - P4132
     - 8c
     - .3.
   * - 213
     - P4132
     - 4b
     - .32
   * - 213
     - P4132
     - 4a
     - .32
   * - 214
     - I4132
     - 48i
     - 1
   * - 214
     - I4132
     - 24h
     - ..2
   * - 214
     - I4132
     - 24g
     - ..2
   * - 214
     - I4132
     - 24f
     - 2..
   * - 214
     - I4132
     - 16e
     - .3.
   * - 214
     - I4132
     - 12d
     - 2.22
   * - 214
     - I4132
     - 12c
     - 2.22
   * - 214
     - I4132
     - 8b
     - .32
   * - 214
     - I4132
     - 8a
     - .32
   * - 215
     - P-43m
     - 24j
     - 1
   * - 215
     - P-43m
     - 12i
     - ..m
   * - 215
     - P-43m
     - 12h
     - 2..
   * - 215
     - P-43m
     - 6g
     - 2.mm
   * - 215
     - P-43m
     - 6f
     - 2.mm
   * - 215
     - P-43m
     - 4e
     - .3m
   * - 215
     - P-43m
     - 3d
     - -42.m
   * - 215
     - P-43m
     - 3c
     - -42.m
   * - 215
     - P-43m
     - 1b
     - -43m
   * - 215
     - P-43m
     - 1a
     - -43m
   * - 216
     - F-43m
     - 96i
     - 1
   * - 216
     - F-43m
     - 48h
     - ..m
   * - 216
     - F-43m
     - 24g
     - 2.mm
   * - 216
     - F-43m
     - 24f
     - 2.mm
   * - 216
     - F-43m
     - 16e
     - .3m
   * - 216
     - F-43m
     - 4d
     - -43m
   * - 216
     - F-43m
     - 4c
     - -43m
   * - 216
     - F-43m
     - 4b
     - -43m
   * - 216
     - F-43m
     - 4a
     - -43m
   * - 217
     - I-43m
     - 48h
     - 1
   * - 217
     - I-43m
     - 24g
     - ..m
   * - 217
     - I-43m
     - 24f
     - 2..
   * - 217
     - I-43m
     - 12e
     - 2.mm
   * - 217
     - I-43m
     - 12d
     - -4..
   * - 217
     - I-43m
     - 8c
     - .3m
   * - 217
     - I-43m
     - 6b
     - -42.m
   * - 217
     - I-43m
     - 2a
     - -43m
   * - 218
     - P-43n
     - 24i
     - 1
   * - 218
     - P-43n
     - 12h
     - 2..
   * - 218
     - P-43n
     - 12g
     - 2..
   * - 218
     - P-43n
     - 12f
     - 2..
   * - 218
     - P-43n
     - 8e
     - .3.
   * - 218
     - P-43n
     - 6d
     - -4..
   * - 218
     - P-43n
     - 6c
     - -4..
   * - 218
     - P-43n
     - 6b
     - 222..
   * - 218
     - P-43n
     - 2a
     - 23.\
   * - 219
     - F-43c
     - 96h
     - 1
   * - 219
     - F-43c
     - 48g
     - 2..
   * - 219
     - F-43c
     - 48f
     - 2..
   * - 219
     - F-43c
     - 32e
     - .3.
   * - 219
     - F-43c
     - 24d
     - -4..
   * - 219
     - F-43c
     - 24c
     - -4..
   * - 219
     - F-43c
     - 8b
     - 23.\
   * - 219
     - F-43c
     - 8a
     - 23.\
   * - 220
     - I-43d
     - 48e
     - 1
   * - 220
     - I-43d
     - 24d
     - 2..
   * - 220
     - I-43d
     - 16c
     - .3.
   * - 220
     - I-43d
     - 12b
     - -4..
   * - 220
     - I-43d
     - 12a
     - -4..
   * - 221
     - Pm-3m
     - 48n
     - 1
   * - 221
     - Pm-3m
     - 24m
     - ..m
   * - 221
     - Pm-3m
     - 24l
     - m..
   * - 221
     - Pm-3m
     - 24k
     - m..
   * - 221
     - Pm-3m
     - 12j
     - m.m2
   * - 221
     - Pm-3m
     - 12i
     - m.m2
   * - 221
     - Pm-3m
     - 12h
     - mm2..
   * - 221
     - Pm-3m
     - 8g
     - .3m
   * - 221
     - Pm-3m
     - 6f
     - 4m.m
   * - 221
     - Pm-3m
     - 6e
     - 4m.m
   * - 221
     - Pm-3m
     - 3d
     - 4/mm.m
   * - 221
     - Pm-3m
     - 3c
     - 4/mm.m
   * - 221
     - Pm-3m
     - 1b
     - m-3m
   * - 221
     - Pm-3m
     - 1a
     - m-3m
   * - 222
     - Pn-3n
     - 48i
     - 1
   * - 222
     - Pn-3n
     - 24h
     - ..2
   * - 222
     - Pn-3n
     - 24g
     - 2..
   * - 222
     - Pn-3n
     - 16f
     - .3.
   * - 222
     - Pn-3n
     - 12e
     - 4..
   * - 222
     - Pn-3n
     - 12d
     - -4..
   * - 222
     - Pn-3n
     - 8c
     - .-3.
   * - 222
     - Pn-3n
     - 6b
     - 42.2
   * - 222
     - Pn-3n
     - 2a
     - 432
   * - 223
     - Pm-3n
     - 48l
     - 1
   * - 223
     - Pm-3n
     - 24k
     - m..
   * - 223
     - Pm-3n
     - 24j
     - ..2
   * - 223
     - Pm-3n
     - 16i
     - .3.
   * - 223
     - Pm-3n
     - 12h
     - mm2..
   * - 223
     - Pm-3n
     - 12g
     - mm2..
   * - 223
     - Pm-3n
     - 12f
     - mm2..
   * - 223
     - Pm-3n
     - 8e
     - .32
   * - 223
     - Pm-3n
     - 6d
     - -4m.2
   * - 223
     - Pm-3n
     - 6c
     - -4m.2
   * - 223
     - Pm-3n
     - 6b
     - mmm..
   * - 223
     - Pm-3n
     - 2a
     - m-3.
   * - 224
     - Pn-3m
     - 48l
     - 1
   * - 224
     - Pn-3m
     - 24k
     - ..m
   * - 224
     - Pn-3m
     - 24j
     - ..2
   * - 224
     - Pn-3m
     - 24i
     - ..2
   * - 224
     - Pn-3m
     - 24h
     - 2..
   * - 224
     - Pn-3m
     - 12g
     - 2.mm
   * - 224
     - Pn-3m
     - 12f
     - 2.22
   * - 224
     - Pn-3m
     - 8e
     - .3m
   * - 224
     - Pn-3m
     - 6d
     - -42.m
   * - 224
     - Pn-3m
     - 4c
     - .-3m
   * - 224
     - Pn-3m
     - 4b
     - .-3m
   * - 224
     - Pn-3m
     - 2a
     - -43m
   * - 225
     - Fm-3m
     - 192l
     - 1
   * - 225
     - Fm-3m
     - 96k
     - ..m
   * - 225
     - Fm-3m
     - 96j
     - m..
   * - 225
     - Fm-3m
     - 48i
     - m.m2
   * - 225
     - Fm-3m
     - 48h
     - m.m2
   * - 225
     - Fm-3m
     - 48g
     - 2.mm
   * - 225
     - Fm-3m
     - 32f
     - .3m
   * - 225
     - Fm-3m
     - 24e
     - 4m.m
   * - 225
     - Fm-3m
     - 24d
     - m.mm
   * - 225
     - Fm-3m
     - 8c
     - -43m
   * - 225
     - Fm-3m
     - 4b
     - m-3m
   * - 225
     - Fm-3m
     - 4a
     - m-3m
   * - 226
     - Fm-3c
     - 192j
     - 1
   * - 226
     - Fm-3c
     - 96i
     - m..
   * - 226
     - Fm-3c
     - 96h
     - ..2
   * - 226
     - Fm-3c
     - 64g
     - .3.
   * - 226
     - Fm-3c
     - 48f
     - 4..
   * - 226
     - Fm-3c
     - 48e
     - mm2..
   * - 226
     - Fm-3c
     - 24d
     - 4/m..
   * - 226
     - Fm-3c
     - 24c
     - -4m.2
   * - 226
     - Fm-3c
     - 8b
     - m-3.
   * - 226
     - Fm-3c
     - 8a
     - 432
   * - 227
     - Fd-3m
     - 192i
     - 1
   * - 227
     - Fd-3m
     - 96h
     - ..2
   * - 227
     - Fd-3m
     - 96g
     - ..m
   * - 227
     - Fd-3m
     - 48f
     - 2.mm
   * - 227
     - Fd-3m
     - 32e
     - .3m
   * - 227
     - Fd-3m
     - 16d
     - .-3m
   * - 227
     - Fd-3m
     - 16c
     - .-3m
   * - 227
     - Fd-3m
     - 8b
     - -43m
   * - 227
     - Fd-3m
     - 8a
     - -43m
   * - 228
     - Fd-3c
     - 192h
     - 1
   * - 228
     - Fd-3c
     - 96g
     - ..2
   * - 228
     - Fd-3c
     - 96f
     - 2..
   * - 228
     - Fd-3c
     - 64e
     - .3.
   * - 228
     - Fd-3c
     - 48d
     - -4..
   * - 228
     - Fd-3c
     - 32c
     - .-3.
   * - 228
     - Fd-3c
     - 32b
     - .32
   * - 228
     - Fd-3c
     - 16a
     - 23.\
   * - 229
     - Im-3m
     - 96l
     - 1
   * - 229
     - Im-3m
     - 48k
     - ..m
   * - 229
     - Im-3m
     - 48j
     - m..
   * - 229
     - Im-3m
     - 48i
     - ..2
   * - 229
     - Im-3m
     - 24h
     - m.m2
   * - 229
     - Im-3m
     - 24g
     - mm2..
   * - 229
     - Im-3m
     - 16f
     - .3m
   * - 229
     - Im-3m
     - 12e
     - 4m.m
   * - 229
     - Im-3m
     - 12d
     - -4m.2
   * - 229
     - Im-3m
     - 8c
     - .-3m
   * - 229
     - Im-3m
     - 6b
     - 4/mm.m
   * - 229
     - Im-3m
     - 2a
     - m-3m
   * - 230
     - Ia-3d
     - 96h
     - 1
   * - 230
     - Ia-3d
     - 48g
     - ..2
   * - 230
     - Ia-3d
     - 48f
     - 2..
   * - 230
     - Ia-3d
     - 32e
     - .3.
   * - 230
     - Ia-3d
     - 24d
     - -4..
   * - 230
     - Ia-3d
     - 24c
     - 2.22
   * - 230
     - Ia-3d
     - 16b
     - .32
   * - 230
     - Ia-3d
     - 16a
     - .-3.


References
----------

For the use of this function, please cite the following paper.

::

    @inproceedings{levy2024symmcd,
    title={Symm{CD}: Symmetry-Preserving Crystal Generation with Diffusion Models},
    author={Daniel Levy and Siba Smarak Panigrahi and S{\'e}kou-Oumar Kaba and Qiang Zhu and Mikhail Galkin and Santiago Miret and Siamak Ravanbakhsh},
    booktitle={AI for Accelerated Materials Design - NeurIPS 2024},
    year={2024},
    url={https://openreview.net/forum?id=V7x2KZQn2v}
    }
