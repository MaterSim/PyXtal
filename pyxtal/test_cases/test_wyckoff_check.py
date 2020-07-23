from structure import *

allpassed = True
for sg in range(1, 231):
    print("Calculating spacegroup " + str(sg))
    wyckoffs = get_wyckoffs(sg)
    for index, wp in enumerate(wyckoffs):
        v = np.random.random(3)
        for i in range(3):
            if np.random.random() < 0.5:
                v[i] *= -1
        # v = SymmOp.from_rotation_and_translation(np.zeros([3,3]), v)
        points = []
        for p in wp:
            points.append(p.operate(v))
        for i, p in enumerate(points):
            for j in range(3):
                a = np.random.random()
                if a < 1 / 3:
                    points[i][j] += 1
                elif a < 2 / 3:
                    points[i][j] -= 1
        if check_wyckoff_position(points, sg) is not False:
            pass
        else:
            allpassed = False
            print("sg: " + str(sg) + ", index: " + str(index))
            print("points:")
            for p in points:
                print(p)
if allpassed is True:
    print("All spacegroups passed.")
