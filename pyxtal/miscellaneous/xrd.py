import numpy as np
from pyxtal.XRD import XRD
import warnings
warnings.filterwarnings("ignore")

from ase.build import bulk
from ase.io import read, lammpsrun
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from time import time

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", "--ncpu", dest="ncpu",
                  help="number of cpus, default: 1",
                  type=int,
                  default=1,
                  metavar="ncpu")
(options, args) = parser.parse_args()

# load calculator
cells = [1, 2, 3, 4] #, 5]
fig = plt.figure(figsize=(9.0, 2*len(cells)))
gs = gridspec.GridSpec(len(cells), 1, hspace=0)
thetas = [20, 90]
a = bulk('GaN', 'rocksalt', a=4.27, cubic=True) * 2

for _i, i in enumerate(cells):
    t0 = time()
    a0 = a*i
    permutation = np.argsort(-1*a0.numbers)
    a0 = a0[permutation]

    ax0 = fig.add_subplot(gs[_i, 0])
    xrd = XRD(a0, thetas=thetas, per_N=3e+4, ncpu=options.ncpu)
    xrd.plot_pxrd(ax=ax0, fontsize=12, res=0.01, fwhm=0.2, profile='gaussian', 
            legend="N={:d} {:.1f}s".format(len(a0), time()-t0))
    ax0.set_xlim(thetas)
    ax0.set_ylim([0, 0.99])
    if _i < len(cells) - 1:
        ax0.xaxis.set_visible(False)
    print(i, a0, time()-t0)
    print(xrd)

fig.savefig('test.pdf')
