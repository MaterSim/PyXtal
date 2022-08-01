from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pyxtal.XRD import XRD
from pyxtal import pyxtal
from pyxtal.util import ase2pymatgen
import numpy as np
from time import time
s = pyxtal()
thetas = (5, 50)
#thetas = (30, 33)

cifs = ['naphthalene', 'NaSb3F10', 'PAHYON01', 'GeF2']
for cif in cifs:
    s.from_seed('pyxtal/database/cifs/'+cif+'.cif')
    for rep in [1, 2, 3]: #, 4]: #, 5, 6]:
        a_struc = s.to_ase() * rep
        p_struc = ase2pymatgen(a_struc)
    
        xs = []
        ys = []
        for i in range(2):
            t0 = time()
            if i == 0:
                method = 'pyxtal'
                xrd0 = XRD(a_struc, thetas=thetas, per_N=3e+4)
                ids = np.where(xrd0.pxrd[:,-1]*100>1e-3)[0]
                xs.append(xrd0.pxrd[ids, 0])
                ys.append(xrd0.pxrd[ids,-1]*100)
                #print(xrd0)
            else:
                method = 'pymatgen'
                c = XRDCalculator()
                xrd0 = c.get_pattern(p_struc, two_theta_range=thetas)
                xs.append(xrd0.x)
                ys.append(xrd0.y)
            print("{:6d} {:8s} {:6.1f}".format(len(a_struc), method, time()-t0))
        if len(xs[0]) == len(xs[1]):
            xdiff = np.abs(xs[1]-xs[0]).max()
            ydiff = np.abs(ys[1]-ys[0]).max()
            if max([xdiff, ydiff]) > 1e-2:
                print('problem', xdiff, ydiff)
            #else:
            #    print(xs[0])
        else:
            if rep < 2:
                print('problem in calculation hkl')
                print('pyxtal')
                print(xs[0])
                print('pymatgen')
                print(xs[1])
            import sys; sys.exit()
