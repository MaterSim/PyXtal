from pyxtal.structure import *

failed = []
sgs = range(2, 231)
for sg in sgs:
    print("calculating sg: " + str(sg))
    wyckoff_symmetry = get_wyckoff_symmetry(sg)
    for j, wp in enumerate(wyckoff_symmetry):
        seen = []
        for i, p in enumerate(wp):
            if i is 0:
                for op in p:
                    rot = SymmOp.from_rotation_and_translation(
                        op.rotation_matrix, [0, 0, 0]
                    )
                    if rot not in seen:
                        seen.append(rot)
            elif i is not 0:
                for op in p:
                    rot = SymmOp.from_rotation_and_translation(
                        op.rotation_matrix, [0, 0, 0]
                    )
                    if rot not in seen:
                        params = [sg]
                        if params not in failed:
                            failed.append(params)
for x in failed:
    print(x)
