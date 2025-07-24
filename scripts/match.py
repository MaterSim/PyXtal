
import pymatgen as mg
import numpy as np
from pyxtal.util import parse_cif
from optparse import OptionParser
import pymatgen.analysis.structure_matcher as sm
from pyxtal import pyxtal
from pyxtal.XRD import Similarity
from pyxtal.optimize.base import GlobalOptimize
import warnings
warnings.filterwarnings("ignore")

parser = OptionParser()
parser.add_option("-f", dest="cif", default="WFS-gaff.cif", help="input cif file")
parser.add_option("-r", dest="ref", help="reference")
parser.add_option("-o", dest="out", default="Matched.cif", help="output")
parser.add_option("--emin", dest="emin", type=float, default=0,
                  help="minimum energy, default 0")
parser.add_option("--emax", dest="emax", type=float, default=100,
                  help="maximum energy, default 100")
parser.add_option("--early_stop", dest="early", action="store_true", default=False,
                  help="stop when the first match is found")
parser.add_option("--XRD", dest="xrd",
                  action="store_true", default=False,
                  help="Compare XRD")
parser.add_option("--smin", dest="smin", type=float, default=0.80,
                  help="min similarity for XRD")

(options, args) = parser.parse_args()
matcher = sm.StructureMatcher(ltol=0.3, stol=0.4, angle_tol=5.0)

with open(options.cif, 'r') as f:
    lines = f.readlines()
    smiles = []
    for l in lines:
        if 'smile' in l:
            smile_str = l.split(':')[1].strip()
            smiles = [s + '.smi' for s in smile_str.split('.')]
            break
print(smiles)
with open(options.out, 'w') as f: f.write(f'smiles: {smile_str}\n')

if options.ref is None:
    raise ValueError("Reference structure is required.")
else:
    pmg_ref = mg.core.Structure.from_file(options.ref)
    xtal = pyxtal(molecular=True)
    xtal.from_seed(pmg_ref, molecules = smiles)
    print(f"Reference Structure loaded from {options.ref} {pmg_ref.density:.3f}")
    pmg_ref.remove_species("H")
    print(xtal)

    if options.xrd:
        thetas = [0, 35.0]
        xrd = xtal.get_XRD(thetas=thetas)
        p_ref = xrd.get_profile(res=0.15, user_kwargs={"FWHM": 0.25})

cifs, engs = parse_cif(options.cif, eng=True)
print("Total Number of Structures:", len(cifs))
engs = np.array(engs)
ids = np.argsort(engs)
cifs = [cifs[id] for id in ids]
engs = engs[ids]
print(f"Min energy in eV: {engs.min()+options.emin/96.485:.4f} {engs.min()+options.emax/96.485:.4f}")
engs_norm = engs - engs.min()  # Normalize energies to the lowest one
engs_norm *= 96.485

# Find the id of energy that is between [options.emin, options.emax]
n1 = np.searchsorted(engs_norm, options.emin, side='left')
n2 = np.searchsorted(engs_norm, options.emax, side='right')
engs = engs[n1:n2]
engs_norm = engs_norm[n1:n2]
cifs = [cifs[id] for id in range(n1, n2)]
ids = ids[n1:n2]
count = 0
xtal = pyxtal(molecular=True)
for id, cif in enumerate(cifs):
    pmg = mg.core.Structure.from_str(cif, fmt='cif')
    try:
        xtal.from_seed(pmg, molecules = smiles)
        #xtal.energy = engs_norm[id]
        spg = xtal.group.number
        den = xtal.get_density()
        raw_eng = engs[id]
        norm_eng = engs_norm[id]
        match = False
        strs = f"Struc {ids[id]:6d}: {spg:3d} {raw_eng:.3f} kJ/mol, {den:.3f} g/cm^3, {norm_eng:.3f}"
        if options.xrd:
            p1 = xtal.get_XRD(thetas=thetas).get_profile(res=0.15, user_kwargs={"FWHM": 0.25})
            sim = Similarity(p1, p_ref, x_range=thetas).value
            if sim > options.smin: match = True
            strs += f'{sim:12.3f} in PXRD similarity'
        else:
            pmg.remove_species("H")
            if abs(pmg.density-pmg_ref.density) <= 0.35:
                strs += '****'
                if matcher.fit(pmg, pmg_ref):
                    match = True
                    
                    pmg_gen = xtal.to_pymatgen()
                    pmg_gen.remove_species("H")

                    ref_copy = pmg_ref.copy()
                    ref_copy.remove_species("H")
                    print(f"Generated structure sites: {len(pmg_gen.sites)}, Reference sites: {len(ref_copy.sites)}")
                    print(f"Species (gen): {[str(sp) for sp in pmg_gen.species]}")
                    print(f"Species (ref): {[str(sp) for sp in ref_copy.species]}")
                    rmsd = matcher.get_rms_dist( ref_copy, pmg_gen)
                    if rmsd is not None:
                        rms_lat, rms_cart = rmsd
                        print(f"RMSD – Lattice: {rms_lat:.3f} Å, Cartesian: {rms_cart:.3f} Å")
                    else:
                        print("RMSD calculation failed: structures not comparable")
                    
    except:
        continue

    if match:
        count += 1
        strs += '+++++++++++'
        label = f"{count}-d{den:.3f}-spg{spg}-e{norm_eng:.3f}"
        if options.xrd: label += f"-s{sim:.3f}"
        with open(options.out, 'a+') as f: f.writelines(xtal.to_file(header=label))
        if options.early: break
    print(strs)
print(f"Found {count} matches")

