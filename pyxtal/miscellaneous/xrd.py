import numpy as np
from pyxtal.XRD import XRD
from ase.build import bulk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from time import time

cells = [1, 2, 3, 4] #, 5]
fig = plt.figure(figsize=(9.0, 2*len(cells)))
gs = gridspec.GridSpec(len(cells), 1, hspace=0)
thetas = [30, 80]
a = bulk('GaN', 'rocksalt', a=4.27, cubic=True) * 2

for _i, i in enumerate(cells):
    t0 = time()
    a0 = a*i
    permutation = np.argsort(-1*a0.numbers)
    a0 = a0[permutation]
    print(i, a0)

    ax0 = fig.add_subplot(gs[_i, 0])
    xrd = XRD(a0, thetas=thetas)
    xrd.plot_pxrd(ax=ax0, fontsize=12, res=0.01, fwhm=0.2, profile='gaussian', 
            legend="N={:d} {:.1f}s".format(len(a0), time()-t0))
    ax0.set_xlim(thetas)
    ax0.set_ylim([0, 0.99])
    if _i < len(cells) - 1:
        ax0.xaxis.set_visible(False)
    print(xrd)

fig.savefig('1.pdf')
