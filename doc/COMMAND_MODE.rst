Run PyXtal executables
==============================

Currently, we provide several utilities to the users so that they can run the code from command line with Python scripting. 
They include:

- ``Pyxtal_symmetry``: a tool to access the symmetry information
- ``Pyxtal_atom``: a tool to generate atomic crystals
- ``Pyxtal_molecule``: a tool to generate molecular crystals
- ``Pyxtal_test``: a tool to test all modules

After a successfull `installation <Installation.html>`_, all of them can be accessed by invoking the ``-h`` command:

::

    $ pyxtal -h

                 ______       _    _          _   
                (_____ \     \ \  / /        | |   
                 _____) )   _ \ \/ / |_  ____| |  
                |  ____/ | | | )  (|  _)/ _  | | 
                | |    | |_| |/ /\ \ |_( ( | | |___
                |_|     \__  /_/  \_\___)_||_|_(___
                       (____/      
    
    
    ----------------------(version 0.0.1 )----------------------
    
    A Python package for random crystal generation
    The source code is available at https://github.com/qzhu2017/pyxtal
    Developed by Zhu's group at University of Nevada Las Vegas
    
    
    Usage: pyxtal [options]
    
    Options:
      -h, --help            show this help message and exit
      -s sg, --symmetry=sg  desired symmetry, number or string, e.g., 36, Pbca, Ih
      -e element, --element=element
                            desired elements: e.g., Li
      -n numIons, --numIons=numIons
                            desired numbers of atoms: 16
      -f factor, --factor=factor
                            volume factor: default 1.0
      -v verbosity, --verbosity=verbosity
                            verbosity: default 0; higher values print more
                            information
      -a attempts, --attempts=attempts
                            number of crystals to generate: default 1
      -o outdir, --outdir=outdir
                            Directory for storing output cif files: default 'out'
      -d dimension, --dimension=dimension
                            desired dimension: (3, 2, 1, 0): default 3
      -t thickness, --thickness=thickness
                            Thickness, in Angstroms, of a 2D crystal, or area of a
                            1D crystal, None generates a value automatically:
                            default None


Among them, ``pyxtal_test`` is mainly used for the internal test. In the following, we explain the rest utilites in detail.
   
PyXtal_symmetry utility
------------------------
``PyXtal_symmetry`` is a utility to handle the generation of moelcular crystals.

- `-d`, the dimension, e.g., ``3``, ``2``, ``1``, ``0``. The defult is ``3``.
- `-s`: the target symmetry (*space*, *layer*, *rod*, *point* group information), either by *string* (e.g., ``Ih``, ``Pbca``) and *integer* (``61``).

::
    
    $ pyxtal_symmetry -s 36

    -- Space group # 36 (Cmc2_1)--
    8b site symm: 1
      x, y, z
      -x, -y, z+1/2
      x, -y, z+1/2
      -x, y, z
      x+1/2, y+1/2, z
      -x+1/2, -y+1/2, z+1/2
      x+1/2, -y+1/2, z+1/2
      -x+1/2, y+1/2, z
    4a site symm: m..
      0, y, z
      0, -y, z+1/2
      1/2, y+1/2, z
      1/2, -y+1/2, z+1/2

::

    $ pyxtal_symmetry -s 20 -d 2
    
    -- Layer group # 20 (p2_122)--
    4d site symm: 1
      x, y, z
      x+1/2, -y, -z
      -x+1/2, y, -z
      -x, -y, z
    2c site symm: .2.
      1/4, y, 0
      3/4, -y, 0
    2b site symm: ..2
      0, 1/2, z
      1/2, 1/2, -z
    2a site symm: ..2
      0, 0, z
      1/2, 0, -z
 
if the ``-s`` tag is not given, it will output the list of all possible symmetry groups for the given dimension.

::

    $ pyxtal_symmetry -d 3
        space_group
    1            P1
    2           P-1
    3            P2
    4          P2_1
    5            C2
    6            Pm
    7            Pc
    8            Cm
    9            Cc
    10         P2/m
    11       P2_1/m
    12         C2/m
    13         P2/c
    14       P2_1/c
    15         C2/c
    16         P222
    17       P222_1
    18     P2_12_12
    19   P2_12_12_1
    20       C222_1
    ...
    ...
    212       P4332
    213      P4_132
    214      I4_132
    215       P-43m
    216       F-43m
    217       I-43m
    218       P-43n
    219       F-43c
    220       I-43d
    221       Pm-3m
    222       Pn-3n
    223       Pm-3n
    224       Pn-3m
    225       Fm-3m
    226       Fm-3c
    227       Fd-3m
    228       Fd-3c
    229       Im-3m
    230       Ia-3d

PyXtal_atom utility
--------------
``PyXtal_atom`` is a utility to handle the generation of atomic crystals.
Typically, four arguments are requried to describe the target structure:

- `-d`, the dimension, e.g., ``3``, ``2``, ``1``, ``0``.
- `-s`: the target symmetry (*space*, *layer*, *rod*, *point* group information), either by *string* (e.g., ``Ih``, ``Pbca``) and *integer* (``61``).
- `-e`: the list of elements, e.g., ``Si``, ``Si, O``
- `-n`: the number of atoms in the target primitive unit cell, e.g., ``12``, ``4, 8``. The size should be consistent with the ``-e`` tag.

For **symmetry group setting**, please refer to the `Group Setting page <Settings.html>`_.
**To our knowledge, PyXtal is perhaps the only open source code which can handle the crystal symmetry generation from 0 to 3 dimensional systems.**
Below we will introduce its capability in detail.

A quick example of C60
~~~~~~~~~~~~~~~~~~~~~~

Below is a quick example to generate a random ``C60`` clusters with icosahedral (``Ih``) symmetry. 

::

    $ pyxtal_atom -e C -n 60 -d 0 -s Ih

                 ______       _    _          _   
                (_____ \     \ \  / /        | |   
                 _____) )   _ \ \/ / |_  ____| |  
                |  ____/ | | | )  (|  _)/ _  | | 
                | |    | |_| |/ /\ \ |_( ( | | |___
                |_|     \__  /_/  \_\___)_||_|_(___
                       (____/      
    
    
    ----------------------(version 0.0.1 )----------------------
    
    A Python package for random crystal generation
    The source code is available at https://github.com/qzhu2017/pyxtal
    Developed by Zhu's group at University of Nevada Las Vegas
    
    
    Symmetry requested: 56(Ih), generated: Ih
    Output to out/C60.xyz


As described in the screen output, the run will generate a file called ``out/C60.xyz`` which stores the structural information about C60.
One can thus visualize via different third-party packages. For instance, below is the output from `VESTA <https://jp-minerals.org/vesta/en/>`_.

.. image:: ../images/C60.png
   :height: 763 px
   :width: 995 px
   :scale: 25 %
   :align: center

Note that this is a random process. So each time the structure is likely to be different.


3D crystals
~~~~~~~~~~~~~~~~~~~~~~
By default, ``-d`` tag is 3, which means to generate 3D crystal. Below is a quick example to generate a diamond like crystals for carbon.

::

    $ pyxtal_atom -e C -n 2 -s 227
    
    Symmetry requested: 227(Fd-3m), generated: Fd-3m
    Output to out/C8.cif


.. image:: ../images/C8-diamond.png
   :height: 763 px
   :width: 763 px
   :scale: 30 %
   :align: center

It is important to note that we specified ``2`` for ``-n`` tag, which means 2 carbon atoms in the primitivel unit cell. Because the space group ``Fd-3m (227)`` is *Face centered*, the resulting conventional unit cell with have ``2*4=8`` atoms.

2D and 1D crystals
~~~~~~~~~~~~~~~~~~~~~~
2D and 1D crystals need one more argument to specify the confinement. For 2D crystal, the ``thickness`` needs to be provided through ``-t`` tag in Angstrom. Below is an example fo generating a 2D MoS2 crystal.

::

    $ pyxtal_atom -e Mo,S -n 1,2 -s 77 -d 2 -t 2.4

    Symmetry requested: 77(p6mm), generated: P6mm
    Output to out/Mo1S2.cif


.. image:: ../images/MoS2.png
   :height: 763 px
   :width: 1263 px
   :scale: 30 %
   :align: center


PyXtal_molecule utility
------------------------

``PyXtal_molecule`` is a utility to handle the generation of moelcular crystals.

Molecular crystals occupying general Wyckoff positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example to generate of random crystal for a famours drug molecule ROY.

::

    $ pyxtal_molecule -e ROY -n 4 -s P2_12_12_1
    
    Symmetry requested: 19 (P2_12_12_1), generated: P2_12_12_1, vol: 2895.37 A^3
    Output to out/S4O8N12C48H36.cif
    
.. image:: ../images/ROY.png
   :height: 763 px
   :width: 963 px
   :scale: 30 %
   :align: center
    
Molecular crystals occupying special Wyckoff positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An import feature of PyXtal is that the program can automatically generate molecular crystals occupying special Wyckoff positions. 
This is very useful for molecules with high internal symmetry. During crystallization, these molecule can occupy some special Wyckoff positions as long as the site symmetry is compatible with the molecular symmetry. For instance, the space group ``Cmc_21`` has 4 symmetry operations (``mm2``) in its primitive cell. However, we can still generate a structure with 2 moleulces for C60 by placing them to the special Wycoff position. This will be automatically processed by our `internal algorithm <Algorithm.html#finding-valid-molecular-orientations>`_.

::

    $ pyxtal_molecule -e C60 -n 2 -s 36

.. image:: ../images/C60-x.png
   :height: 703 px
   :width: 683 px
   :scale: 50 %
   :align: center
 
How to define the molecules?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For the specification of molecule, please ref to the section of `Working with Molecules <Others.html#working-with-molecules>`_.


  
