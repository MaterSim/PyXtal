Note: This page is currently under development.

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

We say that the resulting object has translational symmetry, or that it is periodic. We can be more specific by defining the vectors of translational symmetry. For a given crystal, there are 3 such linearly independent vectors. These 3 vectors, placed into a matrix, define what is called the unit cell. Alternatively, we can define the unit cell using the lengths of each side of the box (usually called a, b, c), along with the angles between them (usually called alpha, beta, gamma). These 6 values are called the cell parameters. The unit cell is any parallepiped-shaped part of the crystal which can be used to generate the rest of the crystal through translations alone. Any unit cell which has the smallest possible volume is called a primitive cell.

Note: a given crystal can have multiple ways to define a primitive cell, and there is not always a clearly preferred choice. Consider a 2-dimensional square lattice. You could just as well define the lattice using parallelograms which run along the diagonal lines:

(TODO: Add image)

To avoid this confusion, there is a set of standards (defined in the `International Tables of Crystallography <https://it.iucr.org/>`_) which is typically used. A cell based on these standards is called the conventional cell. In many cases, the conventional cell is not actually a primitive cell. Instead, the conventional cell may have extra atoms which exist in specific locations within the cell. So, the cell type is determined both by the cell parameters, and by any additional atomic sites within the cell. Below is a list of each type of conventional cell:

(TODO: Add image)

(TODO: Add section discussing the lattice types)

Symmetry Operations
-------------------

Translations are just one kind of transformation operation. More generally, we can perform any 3-dimensional transformation which preserves the lengths and angles between atoms. This means we can also apply rotations, reflections, and inversions, as well as any combination of these. Note that successive operations do not generally commute. That is, the order of operations determines the final outcome.

A symmetry operation is any transformation which leaves the original structure unchanged. In other words, if the structure looks the same before and after a transormation, then that transformation is a symmetry operation of the object. This includes the identity operation (doing nothing to the object), which means that every object has at least a trivial symmetry.

We can artificially split a transformation into two parts: the rotational/inversional part (given by a 3x3 matrix), and the translational part (given by a 3D vector, specifically a 3x1 column matrix). Often, we denote this as a matrix-column pair (P,p) or (P|p), where the capital letter P represents the rotation matrix, and the lowercase letter p represents the translation vector.

We can define the 3x3 rotation matrix by using 3 orthogonal unit vectors as the columns. The resulting matrix is orthogonal, meaning the determinant is either +1 or -1. If only a rotation is applied, then the determinant is +1, and if an inversion is applied, the determinant is -1. If an object has no symmetry operations with determinant -1, it is said to be chiral. In this case, the object's mirror image is different from the original, and cannot be rotated to match its twin. This is especially important for molecules with biochemical applications, since the mirror molecule may have a different effect.

Now, we can define how one operation is applied to another. We consider two operations: (P,p) and (Q,q). If we first apply (P,p), followed by (Q,q), then we get a new operation, which we will call (R,r): (Q,q)(P,p) = (R,r). Note that we "apply" operations from the left. Then, the relationships are:

R = Q*P

r = Q*p + q

where * denotes standard matrix multiplication. From this definition, we see that the rotation is always applied first, followed by the translation. This rule applies for multiple operations as well; with 3 operations (R,r)(Q,q)(P,p), we first apply (P,p), then (Q,q), then (R,r).

Alternatively, the matrix-column pair can be "combined" into a single 4x4 matrix. We simply place the vector to the right of the rotation matrix, place 0's on the bottom row, and place a 1 in the lower right-hand corner:

(TODO: insert matrix image)

This 4x4 matrix is called an affine transformation matrix. With it, we can apply operations using a single matrix multiplication operation. Although this may seem like just a mathematical trick, the affine matrix notation highlights the group structure of the transformations, as it allows translations and rotations to be placed on equal footing. Furthermore, we can use the additional dimension to represent time: the '1' value can be thought of as a single step forward in time, and thus we can define both rotational and translational reference frames (and equivalently, torques and forces) with a single 4x4 matrix. Objects which are (periodically) symmetric in time are called time crystals. Such objects have only recently been synthesized in the lab, and there is likely more research to be done. However, for most applications in crystallography, time is not a factor, and we consider only spatial symmetries.

Groups
------

Symmetry operations have several nice properties, and this allows certain sets of them to be classified as a mathematical object called a group. There are several simple and intuitive examples of groups, which we will discuss below. Formally, a group G is a set of mathematical objects (called elements) with 4 properties:

1) There is a binary operation (often denoted by *) which maps any two elements in the set onto a third element which is also in the set: A*B = C. The operation must be defined for every possible pair on the set, and must map onto an element which is inside of the set.

2) There must be exactly one identity element I which maps every element of the set onto itself: A*I = I*A = A for every A in G.

3) Every element A must have an inverse A^-1, such that multiplication by the inverse gives the identity: A*A^-1 = A^-1*A = I.

4) The operation * must be associative. That is, (A*B)*C = A*(B*C).

One of the simplest examples of a group is the additive group of real integers (Z,+). Here, the set is that of the integers (-1, 0, 1, ...), and the operation is addition. Here, the inverse of a number is just its negative. For example, the inverse of -2 is 2. One can easily verify that the 4 properties listed above hold true for this group. Similarly, we can consider the additive group of real numbers (R,+), or the additive group of complex numbers (C,+).

However, if we replace addition with multiplication, then we no longer have a group, because the element 0 does not have a multiplicitive inverse: any number multiplied by 0 is 0, but any number divided by 0 is undefined. We can fix this by considering the multiplicative group of all numbers except for 0. Or, equivalently, we can consider the multiplicitave group exp(x), where x is any complex number. Then, the inverse is defined as exp(-x), and the identity element is exp(0) = 1.

Interestingly, the  real numbers are a subset of the complex numbers, and yet both the complex numbers and the real numbers form groups in their own right. In this case, we call the real numbers a subgroup of the complex numbers. Likewise, we call the complex numbers a supergroup of the real numbers.

These are all examples of infinite groups, since there are infinitely many points in the complex plane. However, there also exist finite groups. For example, consider the permutation group of 3 objects (we'll call them a, b, and c). Here, our group elements are:

1: (a,b,c)
2: (a,c,b)
3: (b,a,c)
4: (b,c,a)
5: (c,a,b)
6: (c,b,a)

As you can see, there are only 6 elements in this group. Element (1) is the identity, as it represents keeping a, b, and c in their original order. Element (2) represents swapping b and c, element (3) represents swapping a and b, and so on.

Sometimes, it is inconvenient to list every member of a group. Instead, it is often possible to list only a few elements, which can be used to determine the other elements. These elements are called generators. For example, consider elements (2) and (3) in the permutation group shown above. We can define the remaining elements (1, 4, 5, and 6) as follows (with operations acting from the left):

2 * 2 = 1 : (a,c,b) * (a,c,b) = (a,b,c)

2 * 3 = 4 : (a,c,b) * (b,a,c) = (b,c,a)

3 * 4 = 6 : (b,a,c) * (b,c,a) = (c,b,a)

6 * 2 = 5 : (c,b,a) * (a,c,b) = (c,a,b)

Thus, we say that (2) and (3) are generators of the group. Typically, there is not a single "best" choice of generators for a group. We could just as easily have chosen (2) and (6) or (4) and (3) as our generators.

Symmetry Groups
---------------

One can verify that the four properties of groups listed above also hold for our 4x4 transormation matrices. Thus the set of all 3D transformations (with 4x4 matrix multiplication as our operation) forms a group. Because of this, the tools of group theory become available. 