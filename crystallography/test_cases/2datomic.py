from crystallography.database.layergroup import Layergroup
from crystallography.crystal import random_crystal_2D

print("Testing generation of 2d atomic crystals")
for num in range(1,81):
    print("Layer group "+str(num)+"...")
    c12 = random_crystal_2D(num, ['C'], [12], 1.0, 1.0)
    if c12.valid is True:
        print("    Success!")
    else:
        print("    Failed generation----------------")
