from pyxtal.symmetry import Hall, Group
import numpy as np

for n in range(1, 231):
    h = Hall(n)#, permutation=True)
    for hn in h.hall_numbers:
        g = Group(hn, use_hall=True)
        print(g.number, g.hall_number, g.symbol)
        for wp in g:
            print(wp.letter)
            #print("=======")
            wp.set_generators()
            for i, op in enumerate(wp.generators):
                rot = op.rotation_matrix
                tran = op.translation_vector
                if abs(abs(np.linalg.det(rot))-1)>1e-2:
                    print("Error", rot)
                    import sys; sys.exit()
