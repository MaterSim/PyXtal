from pyxtal import pyxtal

cif_path = "pyxtal/database/cifs/"
mol = ["ACBNZA01", "CC(=O)NC1=CC=CC=C1C(=O)N"]

print("read structure from seeds")
for i in range(10):
    c1 = pyxtal(molecular=True)
    c1.from_seed(cif_path+mol[0]+'.cif', molecules=[mol[1]+'.smi'])
    print(i, "angles", c1.mol_sites[0].encode()[-4:-1])

print("generate structure with prespecified torsions")
torsions=[[-60.2, 1.7, 126.5]]
for i in range(10):
    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, [mol[1]+'.smi'], [4], torsions=torsions)
    print(i, "angles", c1.mol_sites[0].encode()[-4:-1])

print("generate structure with random torsions")
torsions=None 
for i in range(100):
    c1 = pyxtal(molecular=True)
    c1.from_random(3, 14, [mol[1]+'.smi'], [4], torsions=torsions)
    print(i, "angles", c1.mol_sites[0].encode()[-4:-1])
