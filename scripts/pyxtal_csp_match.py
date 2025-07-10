import pymatgen as mg
import numpy as np
from pyxtal.util import parse_cif
from optparse import OptionParser
import pymatgen.analysis.structure_matcher as sm
from pyxtal import pyxtal
import warnings
warnings.filterwarnings("ignore")


parser = OptionParser()
parser.add_option("-f", "--cif", dest="cif", default="WFS-gaff.cif",
                  help="cif file name, optional")
parser.add_option("-r", "--ref", dest="ref",
                  help="reference structure")
parser.add_option("--emin", dest="emin", type=float, default=0,
                  help="minimum energy, optional")
parser.add_option("--emax", dest="emax", type=float, default=100,
                  help="maximum energy, optional")
parser.add_option("-c", "--cut", dest="cut", type=int,
                  help="cutoff number, optional")
parser.add_option("--early_stop", dest="early_stop",
                  action="store_true", default=False,
                  help="stop when the first match is found")

(options, args) = parser.parse_args()
matcher = sm.StructureMatcher(ltol=0.3, stol=0.3, angle_tol=5.0)

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
    xtal = pyxtal(molecular=True)
    xtal.from_seed(pmg_ref, molecules = smiles)
    print(f"Reference Structure loaded from {options.ref} {pmg_ref.density:.3f}")
    pmg_ref.remove_species("H")
    print(xtal)

cifs, engs = parse_cif(options.cif, eng=True)
print("Total Number of Structures:", len(cifs))
if options.cut is None:
    cut = len(cifs)
else:
    cut = options.cut

engs = np.array(engs)
ids = np.argsort(engs)

cifs = [cifs[id] for id in ids]
engs = engs[ids]
print(f"Min energy in eV: {engs.min()+options.emin/96.485:.4f} {engs.min()+options.emax/96.485}")
engs -= engs.min()  # Normalize energies to the lowest one
engs *= 96.485

# Find the id of energy that is between [options.emin, options.emax]
n1 = np.searchsorted(engs, options.emin, side='left')
n2 = np.searchsorted(engs, options.emax, side='right')
engs = engs[n1:n2]
cifs = [cifs[id] for id in range(n1, n2)]
ids = ids[n1:n2]

output1 = 'Matched.cif'
count = 0
xtal = pyxtal(molecular=True)
with open(output1, 'w') as f:
    for id, cif in enumerate(cifs):
        pmg = mg.core.Structure.from_str(cif, fmt='cif')
        try:
            xtal.from_seed(pmg, molecules = smiles)
            strs = f"Struc {ids[id]:6d}: {xtal.group.number:3d} {engs[id]:.3f} kJ/mol, {pmg.density:.3f} g/cm^3"
            pmg.remove_species("H")
            if abs(pmg.density-pmg_ref.density) <= 0.15:
                strs += '****'
                if matcher.fit(pmg, pmg_ref):
                    count += 1
                    strs += '+++++++++++'
                    spg = xtal.group.number
                    den = xtal.get_density()
                    eng = engs[id]
                    label = f"{count}-d{den:.3f}-spg{spg}-e{eng:.3f}"
                    f.writelines(xtal.to_file(header=label))
                    if options.early_stop:
                        break
            print(strs)
        except:
            continue

print(f"Found {count} matches")
