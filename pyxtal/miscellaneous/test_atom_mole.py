from pyxtal.db import database

db = database('pyxtal/database/mech.db')

c = db.get_pyxtal('DAHMUX')
rep = c.get_1D_representation()
a = c.to_atomic_xtal()

dicts = c.mol_sites[0].to_1D_dicts()
b = a.to_molecular_xtal([c.mol_sites[0].molecule],
                        reflects = [dicts['reflect']],
                        oris = [dicts['orientation']])
d = rep.to_pyxtal()


import pymatgen.analysis.structure_matcher as sm

pmg1 = c.to_pymatgen()
pmg2 = b.to_pymatgen()
pmg3 = d.to_pymatgen()
print(sm.StructureMatcher().get_rms_dist(pmg1, pmg2))
print(sm.StructureMatcher().get_rms_dist(pmg1, pmg3))
print(sm.StructureMatcher().get_rms_dist(pmg2, pmg3))
