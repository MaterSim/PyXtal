from pyxtal import pyxtal
from ase.io import read
from ase.spacegroup.symmetrize import prep_symmetry
from spglib import get_symmetry_dataset

#ans1 = get_symmetry_dataset(s, symprec=1e-2)
#print(ans1)

s = pyxtal()
s.from_seed('bug.vasp', tol=1e-2)
print(s)
#s1=s.subgroup(eps=0.1, group_type='t+k', max_cell=4)
#for a in s1:
#    print(a)
#s1=s.subgroup(eps=0.1, group_type='k', max_cell=4)
#for a in s1:
#    print(a)
#permutation = {"C":"Si", "Si":"C"}
#for i in range(100):
#    struc = s.subgroup_once(0.01, None, permutation, max_cell=1)
#    print(struc.group.number, struc.formula)

for i in range(100):
    struc = s.subgroup_once(0.2, None, None, 't+k', max_cell=2)
    print(struc.group.number, struc.formula)
#for i in range(1000):
#    struc = s.subgroup_with_substitution(permutation, once=True, max_cell=4)
#    print(struc)
