Note: This page currently under construction

Background and Theory
=====================

Crystals and Structures
-----------------------

When studying solids, it is often useful to describe a material's structure at the atomic level. From this description one can (in theory) determine the material's physical properties, including mechanical strength, electrical and thermal conductivity, melting point, etc. Due to the near-infinite number of possible materials and atomic geometries, it is necessary to have a consistent mathematical framework. This is the job of crystallographers.

For an atomic structure, we can describe the geometry by specifying the type and position of each atom. This works alright for molecules, and is in fact how computers typically encode molecules. But for an ideal crystal, which is infinitely large, it is impossible to describe where each individual atom lies. Fortunately, because crystals are symmetrical, we can specify one part of the crystal, and then use the symmetry operations to "generate" the rest of the crystal. This would create a perfectly symmetrical structure which is infinitely large in size. Such objects do not exist in nature, but they are nevertheless useful for understanding small parts of real, imperfect crystals. So, we call this symmetrical object an ideal crystal.

Most inorganic materials are formed by many small (nearly) ideal crystals called grains. These grains may have different shapes, sizes, and orientations, but each grain has the same crystal structure at the small scale. If we can determine this crystal structure, it becomes possible to predict the way that the grains form and interact with each other. From this, we can go on to predict properties at larger and larger scales, and determine how useful a material will behave in different physical situations. For this reason, determining a material's small-scale crystal structure is absolutely essential for material science and engineering.

At different pressures and temperatures, a material may go through a solid phase transition, and take on a different crystal structure. 

Periodicity, Lattices, and Unit Cells
-------------------------------------

Formally, an ideal crystal is an atomic structure that is periodic in 3 dimensions. This means that when we translate the structure by a certain amount (in any one of 3 directions unique to the crystal), the crystal will look the same. This can be pictured in a few simple steps: 1) Define a small parallelepiped-shaped box. 2) Put atoms into the box (You can put as few or as many atoms as you like). 3) Make a copy of the box and place it adjacent to the original box. 4) Make a copy of the copy, and place that adjacent to the previous one, but along a different axis. 5) Repeat step 4 until you have filled all of space.

(TODO: Add image)

We say that the resulting object has translational symmetry, or that it is periodic. We can be more specific by defining the vectors of translational symmetry. For a given crystal, there are 3 such linearly independent vectors. These 3 vectors, placed into a matrix, define what is called the unit cell. Alternatively, we can define the unit cell using the lengths of each side of the box (usually called a, b, c), along with the angles between them (usually called alpha, beta, gamma). These 6 values are called the cell parameters. The unit cell is any parallepiped-shaped part of the crystal which can be used to generate the rest of the crystal through translations alone. Any unit cell which has the smallest possible volume is called the primitive cell.

Note: a given crystal can have multiple ways to define the primitive cell, and there is not always a clearly preferred choice. Consider a 2-dimensional square lattice. You could just as well define the lattice using parallelograms which run along the diagonal lines:

(TODO: Add image)

To avoid this confusion, there is a set of standards (defined in the `International Tables of Crystallography https://it.iucr.org/>`_) which is typically used. A cell based on these standards is called the conventional cell. In many cases, the conventional cell is not actually a primitive cell. Instead, the conventional cell may have extra atoms which exist in specific locations within the cell. So, the cell type is determined both by the cell parameters, and by any additional atomic sites within the cell. Below is a list of each type of conventional cell:

(TODO: Add image)

Symmetry Operations
-------------------

Translations are just one kind of symmetry operation which can be applied to a crystal. A