from math import fabs
import numpy as np
from random import uniform
from numpy import arccos, cosh, sinh, cbrt, array, dot, pi, log
from numpy.linalg import norm, det, inv
from numpy.random import normal, random
import matplotlib.pyplot as plt
from timeit import timeit

from structure import *

rad = 360.0 / (2 * pi)
deg = 2 * pi / 360.0

identity = array([[1.0, 0, 0], [0, 1, 0], [0, 0, 1]])


def mat_a(a):
    return array([[1, 0, 0], [0, cosh(a), sinh(a)], [0, sinh(a), cosh(a)]])


def mat_b(b):
    return array([[cosh(b), 0, sinh(b)], [0, 1, 0], [sinh(b), 0, cosh(b)]])


def mat_c(c):
    return array([[cosh(c), sinh(c), 0], [sinh(c), cosh(c), 0], [0, 0, 1]])


def strain_matrix(a, b, c):
    a = mat_a(a)
    b = mat_b(b)
    c = mat_c(c)
    raw = a + b + c - 2 * identity
    return raw / cbrt(det(raw))


def shear_matrix(a, b, c):
    return array([[1, a, b], [a, 1, c], [b, c, 1]])


def random_strain():
    a, b, c = 0, 0, 0
    while (
        a < 30 * deg
        or b < 30 * deg
        or c < 30 * deg
        or a > 150 * deg
        or b > 150 * deg
        or c > 150 * deg
    ):
        mat = strain_matrix(normal(scale=0.1), normal(scale=0.1), normal(scale=0.1))
        a, b, c = alpha(mat), beta(mat), gamma(mat)
    return mat


def random_matrix(width=1.0, unitary=False):
    mat = np.zeros([3, 3])
    for x in range(3):
        for y in range(3):
            mat[x][y] = normal(scale=width)
    if unitary:
        new = mat / cbrt(det(mat))
        return new
    else:
        return mat


N_points = 10000
n_bins = 50

x = []
y = []
z = []

a = 0.3

"""for n in range(N_points):
    mat = generate_lattice(1, 100)
    x.append(mat[3]*rad)
    y.append(mat[4]*rad)
    z.append(mat[5]*rad)
    if n % 1000 == 0: print("...")
"""
for n in range(N_points):
    mat = generate_lattice(1, 100)
    x.append(log(mat[0] / mat[1]))
    y.append(log(mat[1] / mat[2]))
    z.append(log(mat[2] / mat[0]))
    if n % 1000 == 0:
        print("...")

fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=True)
plt.title("Lattice generation angles")
# We can set the number of bins with the `bins` kwarg
"""axs[0].hist(x, bins=n_bins)
axs[0].set_xlabel("alpha")
axs[0].set_ylabel("frequency")
axs[1].hist(y, bins=n_bins)
axs[1].set_xlabel("beta")
axs[2].hist(z, bins=n_bins)
axs[2].set_xlabel("gamma")
"""
axs[0].hist(x, bins=n_bins)
axs[0].set_xlabel("a/b")
axs[0].set_ylabel("frequency")
axs[1].hist(y, bins=n_bins)
axs[1].set_xlabel("b/c")
axs[2].hist(z, bins=n_bins)
axs[2].set_xlabel("c/a")

plt.show()
