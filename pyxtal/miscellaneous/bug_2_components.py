from pyxtal import pyxtal

Ag_xyz = """1
AgC2N2H
Ag         4.30800        8.26300       -0.2200
"""

C2N2H7_xyz = """12
AgC2N2H
H          5.95800        5.80600       -0.9530
N          5.24100        6.16800       -1.1210
N          2.23200        6.99000       -0.6820
C          4.02900        5.47000       -1.5870
C          2.78500        5.61100       -0.6610
H          3.69300        5.63500       -2.5830
H          4.17400        4.42800       -1.6990
H          3.12700        5.50500        0.4260
H          2.05000        5.01200       -0.9300
H          1.96000        7.20500       -1.3860
H          1.59400        6.99200       -0.0710
"""
with open('Ag.xyz', "w") as f:
    f.write(Ag_xyz)
with open('C2N2H7.xyz', "w") as f:
    f.write(C2N2H7_xyz)

#for xyz in ['Ag.xyz', 'C2N2H7.xyz']:
#    with open(xyz, 'r') as f:
#        print('filename', xyz, '\n')
#        lines = f.readlines()
#        for l in lines: print(l[:-2])
count = 0
for i in range(100):
    c = pyxtal(molecular=True)
    c.from_random(3, 9, ['Ag.xyz', 'C2N2H7.xyz'], [12, 12])
    short_bonds = c.check_short_distances(r=1.2)
    if len(short_bonds) > 0:
        print(c)
        print(i, len(short_bonds), short_bonds[0])
        c.to_file('bug-'+str(i)+'.cif')
        c.to_ase().write('bug-'+str(i)+'.vasp', format='vasp', vasp5=True, direct=True)
        count += 1
print("\nTotal failures: ", count)
