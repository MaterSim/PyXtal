from pyxtal.reciprocal import RECP
from pyxtal import pyxtal
import matplotlib.pyplot as plt

# Use a modern style
plt.style.use('seaborn-v0_8-darkgrid')

# Customize plot appearance
plt.rcParams.update({
    #'font.family': 'serif',
    'font.size': 14,
    'axes.labelsize': 15,
    'axes.titlesize': 15,
    'legend.fontsize': 14,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 2,
    'figure.dpi': 300
})

recp = RECP(dmax=10.0, nmax=10, lmax=10, rbasis='Bessel')
prototypes = ['diamond',
              'h-diamond',
              #'graphite',
              #'$\alpha$-boron',
              'a-quartz',
              'b-quartz']

fig, axs = plt.subplots(len(prototypes), 1, figsize=(6, 1.6*len(prototypes)))
data = []
for row, prototype in enumerate(prototypes):
    xtal = pyxtal()
    xtal.from_prototype(prototype)
    p, rdf = recp.compute(xtal.to_ase(), norm=True)
    # Plot p values
    if prototype == 'a-quartz':
        label = r'$\alpha$-quartz'
    elif prototype == 'b-quartz':
        label = r'$\beta$-quartz'
    else:
        label = prototype
    label += f" ({xtal.group.number}, {xtal.group.symbol})"
    axs[row].plot(p, label=label)
    axs[row].legend(loc=1)
    axs[row].set_ylabel('$P_{nl}$')

    if row == len(prototypes) - 1:
        axs[row].set_xlabel('Power Spectrum Index')
    else:
        axs[row].set_xticklabels([])
    axs[row].set_ylim(0, 1.0)
plt.tight_layout()
plt.savefig('p-demo.png', dpi=300)
plt.close()
