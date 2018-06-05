from structure import *

for sg in range(1, 195):
    print("====="+str(sg)+"=====")
    symmetry = get_wyckoff_symmetry(sg, molecular=True)
    for i, wp in enumerate(symmetry):
        ops = symmetry[i][0]
        print(ss_string_from_ops(ops, sg))

