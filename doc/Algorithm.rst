Note: This page is currently under development.

How PyXtal Works
================

PyXtal is a free, open source Python package intended to aid with crystal structure prediction. It is our aim to make our tools and theory freely available, in order to expand the knowledge and utilization of crystallographic symmetry. To this end, we outline the basic algorithms used by PyXtal. For more general information about crystallography, see the `Background and Theory <Background.html>`_ page.

Overview
--------

In order to generate random symmetrical crystals, PyXtal uses space groups and their Wyckoff positions as templates, then inserts atoms (or molecules) one Wyckoff position at a time, until the desired stoichiometry is reached. This ensures that the correct symmetry will be obtained, without the need to adjust atomic positions. Along the way, it checks that the atoms are sufficiently far apart. To generate a single random crystal, PyXtal performs roughly seven steps:

1) Get the input parameters from the user. This includes the crystal's stoichiometry, space group, and volume factor. Optionally, the desired lattice and allowed inter-atomic distances can also be defined (which is useful for testing at higher pressures). Otherwise, these parameters will be chosen automatically.

2) Check the compatibility of the stoichiometry with the space group. Because atoms lie in Wyckoff positions, and these can only have specific numbers of atoms in them, not every number of atoms will be able to fit into a unit cell without breaking the symmetry. To check this, we consider the full number of each type of atom, then reduce this number by the size of a Wyckoff position, beginning with the largest (general) position. If the number goes to exactly zero for each atom type, we say the stoichiometry is compatible with the space group. Additionally, we check the degrees of freedom of each Wyckoff position, so as to avoid placing multiple atoms in the same location. For molecular crystals, we also check whether the molecules can be symmetrically oriented into each Wyckoff position.

3) Generate a random lattice consistent with the space group. The cell parameters are based on both the crystal class (which determines the latice angles) and the stoichiometry (which determines the volume). Where there is some leeway for the lattice parameters, a value will be randomly chosen, with preference for more symmetrical values. If the user has defined a lattice, that will be used instead. For randomly generated lattices, we also check that the atoms or molecules can fit into the unit cell without extending outside of it.

4) Begin placing atoms into Wyckoff positions, one atomic specie at a time. First we check the multiplicity of the general Wyckoff position. If at least this number of atoms still needs to be added, then we place the atoms into the general position. If fewer atoms are needed, we instead place the atoms into a special Wyckoff position, beginning with the largest, then decreasing in multiplicity as needed. We choose a random vector between [0,0,0] and (1,1,1), and use this as the generating point for the atoms in the Wyckoff position.

5) Check that the atoms in the single Wyckoff position are not too close together. If the minimum distance is smaller than some tolerance (based on the atomic species), then we merge the Wyckoff position into a smaller special position. To do this, we first group atoms together based on the shortest distances between them, then replace the clusters with single atoms at the clusters' geometric centers. We check the Wyckoff position of the resulting cluster, then continue to merge as needed until the atoms are sufficiently far apart.

6) Check the inter-atomic distance between the newly generated Wyckoff position and the previously generated positions. If the distances are too close, we choose a new generating point within the same original Wyckoff position. For molecular crystals, we also attempt to reorient the molecule, so as to reduce the distance between atoms.

7) Continue to place Wyckoff positions and species into the crystal, until the proper stoichiometry is reached. If any step between 3-6 fails, we repeat it, up to a pre-determined maximum number of attempts. If this still fails, we go to the previous step and retry, up to a different maximum number of attempts. If we succeed, we store the information within the random_crystal class and set random_crystal.valid to True. If we fail after the maximum number of attempts, we output an error message and set random_crystal.valid to False.

If the generation succeeds, the resulting crystal can be (energetically) optimized and compared with other possible structures. After many attempts, the lowest-energy structure will be the one most likely to exist at the specified conditions. This can be used, for example, to determine the phase diagram of a solid without the need to perform any experiments.

Next, we explain the steps listed above in greater detail.

User Input
----------

Checking Compatibility
----------------------

Lattice Generation
------------------

Generation of Wyckoff Positions
-------------------------------

Checking Inter-atomic Distances
-------------------------------

Merging and Checking Wyckoff Positions
--------------------------------------

Checking Wyckoff Site Symmetry
------------------------------

Finding Valid Molecular Orientations
------------------------------------