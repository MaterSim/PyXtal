"""
Extra test cases with rdkit
"""
from pyxtal import pyxtal
import pymatgen.analysis.structure_matcher as sm
import warnings
warnings.filterwarnings("ignore")

s = pyxtal()
s.from_json('pyxtal/database/bug.json')
pmg1=s.to_pymatgen()

s.optimize_lattice(standard=True)
pmg2=s.to_pymatgen()

rms = sm.StructureMatcher().get_rms_dist(pmg1, pmg2)
if rms is None or sum(rms)>1e-4:
    print("Problem in optimizing molecular crystals")
    import sys; sys.exit()

for code in ["QQQCIG04", "AFUVAZ", "HXMTAM", "TRIZIN", "UREAXX02", "TCYETY02", "JAPCIM", "MAMFIT"]:
    s = pyxtal(molecular=True)
    s.from_CSD(code)
    if s.has_special_site():
        sub = s.to_subgroup()
        pmg0 = s.to_pymatgen()
        pmg1 = sub.to_pymatgen()
        rms = sm.StructureMatcher().get_rms_dist(pmg0, pmg1)
        print(code, s.group.number, sub.group.number)
        if rms is None or sum(rms)>1e-4:
            print("Problem ========================")
            break
