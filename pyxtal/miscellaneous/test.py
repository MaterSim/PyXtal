"""
Test the functionality of ss_strin_from_ops
prints Hermann-Mauguin point group symbol for the site symmetry of each Wyckoff
position in each space group. Can be compared with WYCKPOS on the Bilbao server:
http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-table?from=getwp
Note that our notation is slightly different: if only one axis is present, Bilbao
sometimes forgoes placing dots in the symbol. We keep these dots to differentiate
between axes.
"""

from structure import *

letters = "abcdefghijklmnopqrstuvwxyzA"
for sg in range(1, 231):
    print("=====" + str(sg) + "=====")
    symmetry = get_wyckoff_symmetry(sg, molecular=True)
    length = len(symmetry)
    for i, wp in enumerate(symmetry):
        letter = letters[length - i - 1]
        ops = symmetry[i][0]
        print(letter + ": " + ss_string_from_ops(ops, sg))

