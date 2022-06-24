Group Settings
==============

For the output 3D structures, PyXtal uses the conventional standard cell (same
as `Bilbao <http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-def-choice>`_).
Below are the links for each set.

- `Space group <http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-table>`_
- `Layer group <http://www.cryst.ehu.es/cgi-bin/subperiodic/programs/nph-sub_gen?subtype=layer&from=table>`_
- `Rod group <http://www.cryst.ehu.es/cgi-bin/subperiodic/programs/nph-sub_gen?subtype=rod&from=table>`_

One can conveniently access the list of crystallographic point groups, 1D
rod, 2D layer groups and 3D space groups by changing the ``dim`` flag.

.. code-block:: Python

    >>> from pyxtal.symmetry import Group
    >>> g=Group.list_groups(dim=0)
    point_group
    1           C1
    2           Ci
    3           C2
    ...
    56          Ih
    57          C*
    58         C*h

Space Group
-----------

By default, pyxtal follows the standard according to the Volume A of
International Tables for Crystallography.
They are defined as: unique axis b setting,
cell choice 1 for monoclinic groups, hexagonal axes setting for rhombohedral
groups, and origin choice 2 (origin in -1) for the centrosymmetric groups listed
with respect to two origins. The relation between the standard space group
and hall numbers are shown as follows,

.. code-block:: Python

	pyxtal_hall_numbers = [
	1,   2,   3,   6,   9,   18,  21,  30,  39,  57,  60,  63,  72,  81,  90,
	108, 109, 112, 115, 116, 119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
	155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 203, 209, 212, 215, 218,
	221, 227, 229, 230, 234, 239, 245, 251, 257, 263, 266, 269, 275, 279, 284,
	290, 292, 298, 304, 310, 313, 316, 323, 334, 336, 337, 338, 341, 343, 349,
	350, 351, 352, 353, 354, 355, 356, 357, 358, 360, 362, 363, 365, 366, 367,
	368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382,
	383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
	398, 399, 400, 401, 403, 405, 406, 407, 409, 411, 412, 413, 415, 417, 418,
	419, 421, 423, 424, 425, 427, 429, 430, 431, 432, 433, 435, 436, 438, 439,
	440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 452, 454, 455, 456, 457,
	458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
	475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489,
	490, 491, 492, 493, 494, 496, 497, 499, 500, 501, 502, 503, 504, 505, 506,
	507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 519, 520, 522, 523,
	524, 526, 528, 529, 530]


However, in some programs like Spglib, when converting from space group to Hall
numbers, the first description of the space-group type in International Tables
for Crystallography) is chosen. In this case, the Hall number ``525`` (instead
of ``526``) will be chosen for the space group ``227``.

.. code-block:: Python

	spglib_hall_numbers = [
	1,   2,   3,   6,   9,   18,  21,  30,  39,  57,  60,  63,  72,  81,  90,
	108, 109, 112, 115, 116, 119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
	155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 203, 209, 212, 215, 218,
	221, 227, 228, 230, 233, 239, 245, 251, 257, 263, 266, 269, 275, 278, 284,
	290, 292, 298, 304, 310, 313, 316, 322, 334, 335, 337, 338, 341, 343, 349,
	350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 361, 363, 364, 366, 367,
	368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382,
	383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
	398, 399, 400, 401, 402, 404, 406, 407, 408, 410, 412, 413, 414, 416, 418,
	419, 420, 422, 424, 425, 426, 428, 430, 431, 432, 433, 435, 436, 438, 439,
	440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 452, 454, 455, 456, 457,
	458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
	475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489,
	490, 491, 492, 493, 494, 495, 497, 498, 500, 501, 502, 503, 504, 505, 506,
	507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 520, 521, 523,
	524, 525, 527, 529, 530]


Layer Group
-----------

For 2D structures, we use unique axis c for monoclinic layer groups 3-7, and
unique axis a for layer groups 8-18. When two origin choices are available,
we use origin choice 1. We always choose c as the non-periodic axis.

Rod Group
---------

For 1D structures, we use unique axis a for monoclinic Rod groups 3-7, and
unique axis c for Rod groups 8-12. When two settings are available for a group,
we use the 1st setting. We always choose c as the periodic axis.


Point Group
-----------

For point group structures, we use unique axis c for all groups except the
polyhedral groups ``T, Th, O, Td, Oh, I, and Ih``. For all of these groups,
we place the 2-fold rotation about the z axis and a 3-fold rotation about the
(x,x,x) axis. For ``I`` and ``Ih``, we use a 5-fold rotation about the axis (1,
:math:`\tau`, 0), where :math:`\tau` is the golden ratio 1.618.

All supported point groups, listed by number:

+------------+------------+-----------+-----------+
| 1: C1      | 2: Ci      | 3: C2     | 4: Cs     |
+------------+------------+-----------+-----------+
| 5: C2h     | 6: D2      | 7: C2v    | 8: D2h    |
+------------+------------+-----------+-----------+
| 9: C4      | 10: S4     | 11: C4h   | 12: D4    |
+------------+------------+-----------+-----------+
| 13: C4v    | 14: D2d    | 15: D4h   | 16: C3    |
+------------+------------+-----------+-----------+
| 17: C3i    | 18: D3     | 19: C3v   | 20: D3d   |
+------------+------------+-----------+-----------+
| 21: C6     | 22: C3h    | 23: C6h   | 24: D6    |
+------------+------------+-----------+-----------+
| 25: C6v    | 26: D3h    | 27: D6h   | 28: T     |
+------------+------------+-----------+-----------+
| 29: Th     | 30: O      | 31: Td    | 32: Oh    |
+------------+------------+-----------+-----------+
| 33: C5     | 34: C7     | 35: C8    | 36: D5    |
+------------+------------+-----------+-----------+
| 37: D7     | 38: D8     | 39: C5v   | 40: C7v   |
+------------+------------+-----------+-----------+
| 41: C8v    | 42: C5h    | 43: D5h   | 44: D7h   |
+------------+------------+-----------+-----------+
| 45: D8h    | 46: D4d    | 47: D5d   | 48: D6d   |
+------------+------------+-----------+-----------+
| 49: D7d    | 50: D8d    | 51: S6    | 52: S8    |
+------------+------------+-----------+-----------+
| 53: S10    | 54: S12    | 55: I     | 56: Ih    |
+------------+------------+-----------+-----------+
| 57: C*     | 58: C*h    |           |           |
+------------+------------+-----------+-----------+

In addition to the 32 crystallographic point group , we add the following finite
non-crystallographic point groups:

``Cn, Cnh, Cnv, Sn, Cni, Dn, Dnh, Dnd.``

where n should be replaced by an integer. I and Ih, which are the icosahedral
and full icosahedral groups, are particularly useful (Buckminsterfullerene,
for example has point group symmetry ``Ih``). Finally, the infinite rotational
and dihedral point groups ``C*`` and ``C*h`` can be used for generating linear
structures. ``C*h`` will have mirror symmetry, while ``C*`` will not.
