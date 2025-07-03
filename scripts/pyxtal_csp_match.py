import pymatgen as mg
import numpy as np
from pyxtal.util import parse_cif
from optparse import OptionParser
import pymatgen.analysis.structure_matcher as sm
from pyxtal import pyxtal
import warnings
warnings.filterwarnings("ignore")


parser = OptionParser()
parser.add_option("-f", "--cif", dest="cif",
                  help="cif file name, optional")
parser.add_option("-r", "--ref", dest="ref",
                  help="reference structure")
parser.add_option("-b", "--n1", dest="n1", type=int, default=0,
                  help="starting id, optional")
parser.add_option("-e", "--n2", dest="n2", type=int, default=-1,
                  help="ennding id, optional")
parser.add_option("-c", "--cut", dest="cut", type=int,
                  help="cutoff number, optional")
parser.add_option("--early_stop", dest="early_stop", 
                  action="store_true", default=False,
                  help="stop when the first match is found")

(options, args) = parser.parse_args()
n1 = options.n1
n2 = options.n2

with open(options.cif, 'r') as f:
    lines = f.readlines()
    smiles = []
    for l in lines:
        if 'smile' in l:
            smile_str = l.split(':')[1].strip()
            smiles = [smile_str + '.smi']
            break
print(smiles)

if options.ref is None:
    raise ValueError("Reference structure is required.")
else:
    pmg_ref = mg.core.Structure.from_file(options.ref)
    pmg_ref.remove_species("H")
    print(f"Reference Structure loaded from {options.ref}")

cifs, engs = parse_cif(options.cif, eng=True)
print("Total Number of Structures:", len(cifs))
if options.cut is None:
    cut = len(cifs)
else:
    cut = options.cut

engs = np.array(engs)
ids = np.argsort(engs)

cifs = [cifs[id] for id in ids]
engs = [engs[id] for id in ids]

if n2 == -1:
    cifs = cifs[n1:]
    engs = engs[n1:]
else:
    n2 = min(n2, len(cifs))
    cifs = cifs[n1:n2]
    engs = engs[n1:n2]

count = 0
xtal = pyxtal(molecular=True)
for id, cif in enumerate(cifs):
    pmg = mg.core.Structure.from_str(cif, fmt='cif')
    xtal.from_seed(pmg, molecules = smiles)
    print(f"Struc {id + n1:4d}: {xtal.group.number:3d} Eng: {engs[id]:.3f} eV, Den: {pmg.density:.3f} g/cm^3")
    pmg.remove_species("H")
    if abs(pmg.density-pmg_ref.density)<=0.05 and sm.StructureMatcher().fit(pmg, pmg_ref):
        print(f"Struc {ids[id]} is matched with reference structure.")
        if options.early_stop:
            break
