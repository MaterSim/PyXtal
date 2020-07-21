from pyxtal.crystal import *


def check_cluster_distances(cluster, tol):
    try:
        dm = distance_matrix(
            cluster.cart_coords, cluster.cart_coords, Euclidean_lattice, PBC=[0, 0, 0]
        )
    except:
        dm = distance_matrix(
            cluster.coordinates, cluster.coordinates, Euclidean_lattice, PBC=[0, 0, 0]
        )
    for i, x in enumerate(dm):
        for j, y in enumerate(x):
            if i != j:
                if y < tol:
                    print("Found small distance: " + str(y))
                    return False
    return True


for i in range(20):
    print("Testing structure # " + str(i))
    pg = choose([1])
    c = random_cluster(pg, ["Mo"], [38], 1.0)
    tol = c.tol_matrix.get_tol("Mo", "Mo")
    if check_cluster_distances(c, tol):
        continue
    else:
        exit()

print("Success. No small distances found")
