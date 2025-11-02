from pyxtal.symmetry import Group
import numpy as np

fcs = []
bcs = []
acs = []
ccs = []
#Table 2.2.13.2 from International Tables for Crystallography, Vol. A
screw_21a = []
screw_41a = []
screw_42a = []
screw_43a = []
screw_21b = []
screw_41b = []
screw_42b = []
screw_43b = []
screw_21c = []
screw_41c = []
screw_42c = []
screw_43c = []
screw_31c = []
screw_32c = []
screw_61c = []
screw_62c = []
screw_63c = []
screw_64c = []
screw_65c = []

b_glide_a = [] # 0kl, k=2n
c_glide_a = [] # 0kl, l=2n
n_glide_a = [] # 0kl, k+l=2n
d_glide_a = [] # 0kl, k-l=4n
c_glide_b = [] # h0l, l=2n
a_glide_b = [] # h0l, h=2n
n_glide_b = [] # h0l, h+l=2n
d_glide_b = [] # h0l, h+l=4n
a_glide_c = [] # hk0, h=2n
b_glide_c = [] # hk0, k=2n
n_glide_c = [] # hk0, h+k=2n
d_glide_c = [] # hk0, h+k=4n
cn_glide_110 = [] # hhl, l=2n
d_glide_110 = [] # hhll=2n, 2h+l=4n
an_glide_011 = [] # hkk, h=2n
d_glide_011 = [] # hkk, h=2n, 2k+h=4n
bn_glide_101 = [] # hkh, k=2n
d_glide_101 = [] # hkh, k=2n, 2h+k=4n


for i in range(1, 231):
    g = Group(i)
    if g.symbol[0] == 'F':
        fcs.append(i)
    elif g.symbol[0] == 'I':
        bcs.append(i)
    elif g.symbol[0] == 'C':
        ccs.append(i)
    elif g.symbol[0] == 'A':
        acs.append(i)
    ss = g.get_spg_symmetry_object()
    matrix = ss.to_matrix_representation() # 15 * 18
    # 18 columns: 1  -1 2  2_1 m  a  b  c  n  d  3  3_1 3_2 4  4_1 4_2 4_3 -4
    # 15 row: (100) ..... https://pyxtal.readthedocs.io/en/latest/Symmetry_representation.html

    # srew axis
    if matrix[0, 3] == 1: screw_21a.append(i)
    if matrix[0, -4] == 1: screw_41a.append(i)
    if matrix[0, -3] == 1: screw_42a.append(i)
    if matrix[0, -2] == 1: screw_43a.append(i)
    if matrix[1, 3] == 1: screw_21b.append(i)
    if matrix[1, -4] == 1: screw_41b.append(i)
    if matrix[1, -3] == 1: screw_42b.append(i)
    if matrix[1, -2] == 1: screw_43b.append(i)
    if matrix[2, 3] == 1: screw_21c.append(i)
    if matrix[2, -4] == 1: screw_41c.append(i)
    if matrix[2, -3] == 1: screw_42c.append(i)
    if matrix[2, -2] == 1: screw_43c.append(i)
    if matrix[2, 11] == 1: screw_31c.append(i)
    if matrix[2, 12] == 1: screw_32c.append(i)

    if matrix[2, 3] == 1 and matrix[2, 10] == 1: # 2_1, 3 => 6_3 axis
        screw_63c.append(i)
    if matrix[2, 2] == 1 and matrix[2, 12] == 1: # 2, 3_2 => 6_2 axis
        screw_62c.append(i)
    if matrix[2, 2] == 1 and matrix[2, 11] == 1: # 2, 3_1
        screw_64c.append(i)
    if matrix[2, 3] == 1 and matrix[2, 11] == 1: # 2_1, 3_1
        screw_61c.append(i)
    if matrix[2, 3] == 1 and matrix[2, 12] == 1: # 2_1, 3_2
        screw_65c.append(i)

    # glide planes
    if matrix[0, 6] == 1: b_glide_a.append(i)
    if matrix[0, 7] == 1: c_glide_a.append(i)
    if matrix[0, 8] == 1: n_glide_a.append(i)
    if matrix[0, 9] == 1: d_glide_a.append(i)
    if matrix[1, 5] == 1: a_glide_b.append(i)
    if matrix[1, 7] == 1: c_glide_b.append(i)
    if matrix[1, 8] == 1: n_glide_b.append(i)
    if matrix[1, 9] == 1: d_glide_b.append(i)
    if matrix[2, 5] == 1: a_glide_c.append(i)
    if matrix[2, 6] == 1: b_glide_c.append(i)
    if matrix[2, 8] == 1: n_glide_c.append(i)
    if matrix[2, 9] == 1: d_glide_c.append(i)
    # 110
    if 1 in [matrix[8, 7], matrix[8, 8]]: # c, n, d
        cn_glide_110.append(i)
    if 1 == matrix[8, 9]: # c, n, d
        d_glide_110.append(i)
    # 011
    if 1 in [matrix[10, 5], matrix[10, 8]]:
        an_glide_011.append(i)
    elif 1 == matrix[10, 9]: # a, n, d
        d_glide_011.append(i)
    # 101
    if 1 in [matrix[12, 6], matrix[12, 8]]:
        bn_glide_101.append(i)
    elif 1 == matrix[12, 9]:
        d_glide_101.append(i)

dicts = {"fcs (all odd/even)": fcs,
         "bcs (h+k+l=2n)": bcs,
         "acs (k+l=2n)": acs,
         "ccs (h+k=2n)": ccs,
         "screw_21a (h00), h=2n": screw_21a,
         "screw_41a (h00), h=4n": screw_41a,
         "screw_42a (h00), h=2n": screw_42a,
         "screw_43a (h00), h=4n": screw_43a,
         "screw_21b (0k0), k=2n": screw_21b,
         "screw_41b (0k0), k=4n": screw_41b,
         "screw_42b (0k0), k=2n": screw_42b,
         "screw_43b (0k0), k=4n": screw_43b,
         "screw_21c (00l), l=2n": screw_21c,
         "screw_41c (00l), l=4n": screw_41c,
         "screw_42c (00l), l=2n": screw_42c,
         "screw_43c (00l), l=4n": screw_43c,
         "screw_31c (00l), l=3n": screw_31c,
         "screw_32c (00l), l=3n": screw_32c,
         "screw_61c (00l), l=6n": screw_61c,
         "screw_62c (00l), l=3n": screw_62c,
         "screw_63c (00l), l=2n": screw_63c,
         "screw_64c (00l), l=3n": screw_64c,
         "screw_65c (00l), l=6n": screw_65c,
         "b_glide_a (0kl), k=2n": b_glide_a,
         "c_glide_a (0kl), l=2n": c_glide_a,
         "n_glide_a (0kl), k+l=2n": n_glide_a,
         "d_glide_a (0kl), k+l=4n": d_glide_a,
         "a_glide_b (h0l), h=2n": a_glide_b,
         "c_glide_b (h0l), l=2n": c_glide_b,
         "n_glide_b (h0l), h+l=2n": n_glide_b,
         "d_glide_b (h0l), h+l=4n": d_glide_b,
         "a_glide_c (0kl), h=2n": a_glide_c,
         "b_glide_c (0kl), l=2n": b_glide_c,
         "n_glide_c (0kl), h+l=2n": n_glide_c,
         "d_glide_c (0kl), l=2n": d_glide_c,
         "cn_glide_110 (hhl), l=2n": cn_glide_110,
         "an_glide_011 (hkk), h=2n": an_glide_011,
         "bn_glide_101 (hkh), k=2n": bn_glide_101,
         "d_glide_110 (hhl), l=2n 2h+l=4n": d_glide_110,
         "d_glide_011 (hkk), h=2n 2k+h=4n": d_glide_011,
         "d_glide_101 (hkh), k=2n 2h+k=4n": d_glide_101,
}
for key in dicts:
    print('\n', key)
    for item in dicts[key]:
        spg = Group(item)
        print(f'{item}: {spg.symbol}')

for key in dicts:
    print(key, dicts[key])
