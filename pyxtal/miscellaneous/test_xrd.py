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
ids = [0, 23, 24, 25, 26, 27, 28, 29]
pth = "META-16000-220000-1.0-0.5ps/META_" #0.dump"

fig = plt.figure(figsize=(9.0, 2.0*(len(ids)+1)))
gs = gridspec.GridSpec(len(ids)+1, 1, hspace=0)
thetas = [20, 90]
for _i, id in enumerate(ids):
    t0 = time()
    a0 = read(pth+str(id)+".dump", format='lammps-dump-text')
    permutation = np.argsort(-1*a0.numbers)
    a0 = a0[permutation]
    N = int(len(a0)/2)
    a0.set_atomic_numbers([31]*N + [7]*N)
    ax0 = fig.add_subplot(gs[_i, 0])
    xrd = XRD(a0, thetas=thetas, ncpu=options.ncpu)
    if _i == 0:
        legend = 'B4'
    else:
        legend = 'META_'+str(id)
        
    print(a0, time()-t0)    
    xrd.plot_pxrd(ax=ax0, res=0.01, fwhm=0.25, profile='gaussian', legend=legend)
    ax0.set_xlim(thetas)
    ax0.set_ylim([0, 0.99])
    ax0.xaxis.set_visible(False)

ax0 = fig.add_subplot(gs[-1, 0])
a0 = bulk('GaN', 'rocksalt', a=4.05, cubic=True)
xrd = XRD(a0, thetas=thetas)
legend = 'B1'
xrd.plot_pxrd(ax=ax0, res=0.01, fwhm=0.25, profile='gaussian', legend=legend)
ax0.set_xlim(thetas)
ax0.set_ylim([0, 0.99])


fig.savefig('16k.pdf')
