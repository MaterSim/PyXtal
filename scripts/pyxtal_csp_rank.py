import pymatgen as mg
import numpy as np
from pyxtal import pyxtal
from pyxtal.util import parse_cif
from optparse import OptionParser
import pymatgen.analysis.structure_matcher as sm

import warnings
warnings.filterwarnings("ignore")

def new_struc(xtal, xtals, max_num=100):
    """
    check if this is a new structure

    Args:
        xtal: input structure
        xtals: list of reference structures

    Return:
        `None` or the id of matched structure
    """
    if len(xtals) > max_num:
        start = len(xtals) - max_num
    else:
        start = 0
    sg1 = xtal.group.number
    pmg_s1 = xtal.to_pymatgen()
    pmg_s1.remove_species("H")
    vol1 = pmg_s1.lattice.volume

    for xtal2 in xtals[start:]:
        sg2 = xtal2.group.number
        if sg1 == sg2:
            pmg_s2 = xtal2.to_pymatgen()
            vol2 = pmg_s2.lattice.volume
            if abs(vol1-vol2)/vol1<5e-2:
                pmg_s2.remove_species("H")
                if sm.StructureMatcher().fit(pmg_s1, pmg_s2):
                    return False
    return True


parser = OptionParser()
parser.add_option("-f", "--cif", dest="cif",
                  help="cif file name, optional")
parser.add_option("-r", "--rank", dest="rank", default='energy',
                  help="ranking criteria: default is energy")
parser.add_option("-s", "--n1", dest="n1", type=int, default=0,
                  help="starting id, optional")
parser.add_option("-e", "--n2", dest="n2", type=int, default=-1,
                  help="ennding id, optional")
parser.add_option("-c", "--cut", dest="cut", type=int,
                  help="cutoff number, optional")

(options, args) = parser.parse_args()
rank = options.rank
n1 = options.n1
n2 = options.n2
output1 = 'Ranked.cif'

"""
Read the smile from the following contents
-------Global Crystal Structure Prediction------
smile     : CC(=O)OC1=CC=CC=C1C(=O)O
"""

with open(options.cif, 'r') as f:
    lines = f.readlines()
    smiles = []
    for l in lines:
        if 'smile' in l:
            smile_str = l.split(':')[1].strip()
            smiles = [smile_str + '.smi']
            break
print(smiles)

cifs, engs = parse_cif(options.cif, eng=True)
print("Total Number of Structures:", len(cifs))
if options.cut is None:
    cut = len(cifs)
else:
    cut = options.cut

if options.rank == 'energy':
    engs = np.array(engs)
    ids = np.argsort(engs)
else:
    cifs, sims = parse_cif(options.cif, sim=True)
    sims = np.array(sims)
    ids = np.argsort(-1*sims)
    sims = [sims[id] for id in ids]

cifs = [cifs[id] for id in ids]
engs = [engs[id] for id in ids]

if n2 == -1:
    cifs = cifs[n1:]
    engs = engs[n1:]
else:
    n2 = min(n2, len(cifs))
    cifs = cifs[n1:n2]
    engs = engs[n1:n2]
with open(output1, 'w') as f: f.write(l)

xtals = []
count = 0
with open(output1, 'a+') as f:
    for id, cif in enumerate(cifs):
        pmg = mg.core.Structure.from_str(cif, fmt='cif')
        try:
            xtal = pyxtal(molecular=True)
            xtal.from_seed(pmg, molecules = smiles)
            if new_struc(xtal, xtals, 100):
                xtals.append(xtal)
                spg = xtal.group.number
                den = xtal.get_density()
                eng = engs[id]
                label = f"{count}-d{den:.3f}-spg{spg}-e{eng:.3f}"
                f.writelines(xtal.to_file(header=label))
                print(f"{ids[id]:6d} {count:4d} {label}")
                count += 1
                if count == cut:
                    break
        except:
            print("Problem in reading")
