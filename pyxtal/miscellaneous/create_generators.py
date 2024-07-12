from pandas import DataFrame

from pyxtal.crystal import *

fpath = "point_generators_new.csv"
PBC = [0, 0, 0]


def rounded(op):
    v1 = op.translation_vector
    v2 = v1 - np.floor(v1)
    if [1, 1, 1] == PBC:
        return SymmOp.from_rotation_and_translation(op.rotation_matrix, v2)
    elif [1, 1, 0] == PBC:
        return SymmOp.from_rotation_and_translation(
            op.rotation_matrix, [v2[0], v2[1], v1[2]]
        )
    elif [0, 0, 1] == PBC:
        return SymmOp.from_rotation_and_translation(
            op.rotation_matrix, [v1[0], v1[1], v2[2]]
        )
    elif [0, 0, 0] == PBC:
        return op
    return None


print("-------------------Creating generators-------------------")
generators = [None]
# Loop over spacegroups
for sg in range(1, 33):
    print("Calculating spacegroup: " + str(sg))
    sg_gen = []
    wyckoffs = get_point(sg)
    gen_pos = wyckoffs[0]
    if gen_pos == []:
        print(wyckoffs)
    # Loop over Wyckoff positions
    for _i, wp in enumerate(wyckoffs):
        wp_gen = [gen_pos[0].as_xyz_string()]
        first = wp[0]
        seen = [first]
        for op in gen_pos:
            if rounded(op * first) not in seen:
                seen.append(rounded(op * first))
                wp_gen.append(op.as_xyz_string())
        sg_gen.append(wp_gen)
    generators.append(sg_gen)

failed = False

"""print("-------------------Checking generators-------------------")
for sg in range(1, 231):
    print("Checking spacegroup: "+str(sg))
    wyckoffs = get_wyckoffs(sg)
    generators = get_wyckoff_generators(sg)
    for i, wp in enumerate(wyckoffs):
        temp = []
        for x in wp:
            temp.append(rounded(x))
        new = []
        for op in generators[i]:
            new.append(rounded(op*wp[0]))
        for op in new:
            if op in temp:
                temp.remove(op)
            else: print("error")
        if temp != []:
            print("Failed: sg="+str(sg)+", wp="+str(i))
            failed = True"""

if not failed:
    print("-------------------Writing to database-------------------")
    array = np.array(generators)
    df = DataFrame(data=array)
    df.to_csv(fpath)
