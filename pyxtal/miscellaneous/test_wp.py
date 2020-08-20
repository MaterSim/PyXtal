import numpy as np
from pymatgen.core.operations import SymmOp
from sys import exit

data = np.load("wyckoff_symmetry.npy")
identity = data[1][0][0]
# loop sg's
for x in range(1, 231):
    seen = []
    # loop Wyckoff positions
    for i, y in enumerate(data[x]):
        # loop points
        for j, z in enumerate(y):
            if z not in seen:
                seen.append(z)
            else:
                if z != identity and z not in y:
                    print("Error. Seen:")
                    for x in seen:
                        print(x)
                    print("Current:")
                    print(z)
                    exit()
    print(x)
