from pyxtal.xtal import RECP
from pyxtal import pyxtal
import numpy as np

xtal1 = pyxtal(); xtal1.from_prototype('diamond') #; xtal.to_file('dia.cif'); print(xtal)
xtal2 = pyxtal(); xtal2.from_prototype('cBN');
xtal3 = pyxtal(); xtal3.from_prototype('h-diamond') #; xtal.to_file('dia.cif'); print(xtal)
data = [(xtal1, 'dia'), (xtal2, 'cBN'), (xtal3, 'h-dia')]
for i in range(3, 5):
    xtal = pyxtal()
    xtal.from_seed(f'sub{i}.cif')
    data.append((xtal, f'sub{i}'))

data0 = []
recp = RECP(dmax=10.0, nmax=20, lmax=8)
for i, (xtal, label) in enumerate(data):
    d, rdf = recp.compute(xtal.to_ase(), norm=True)
    if i > 0:
        print(f"{label:10s} : {np.linalg.norm(d - data0[0][0]):12.4f}")
    data0.append((d, rdf, label))
recp.plot(data0, filename='reciprocal_dia.png')
