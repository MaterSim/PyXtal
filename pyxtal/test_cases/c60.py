"""
Generates n C60 structures and outputs to cif files
"""
from pyxtal.crystal import *
from os import mkdir
from os.path import exists

# input variables
g = Group("Ih", dim=0)
species = ["C"]
numatoms = [60]
factor = 1.0
n = 10

# Generate n clusters
i = 1
outstructs = []
print("Generating " + str(n) + " structures...")
print("1")
while i <= n:
    c = random_cluster(g, species, numatoms, factor)
    if c.valid:
        outstructs.append(c.struct)
        i += 1
        if i <= n:
            print(i)

# Output files
print("Outputting cif files to c60_out")
if not exists("c60_out"):
    mkdir("c60_out")
for i, s in enumerate(outstructs):
    s.to(filename="c60_out/c60_" + str(i + 1) + ".cif")
