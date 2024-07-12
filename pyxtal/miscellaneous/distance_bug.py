import numpy as np

from pyxtal.crystal import distance_matrix, random_crystal


def find_short(dm):
    for x in dm:
        for y in x:
            if y > 0 and y < 1.07:
                print("  ~~~ Found short distance: " + str(y) + " ~~~")


n = 1000
rng = np.random.default_rng()
print("Testing " + str(n) + " structures...")


for i in range(n):
    if i % 10 == 0:
        print(str(i) + "/" + str(n) + " structures tested.")
    pg = rng.integers(2, 57)
    c = random_crystal(0, pg, ["Mo"], [16], 1.0)
    if c.valid:
        dm = distance_matrix(c.frac_coords, c.frac_coords, c.lattice.matrix, PBC=c.PBC)
        find_short(dm)
