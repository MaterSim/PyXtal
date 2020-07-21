"""Module for converting between hall and Hermann-Mauguin spacegroup numbers"""
hall_numbers = [
    0,
    1,
    2,
    3,
    6,
    9,
    18,
    21,
    30,
    39,
    57,
    60,
    63,
    72,
    81,
    90,
    108,
    109,
    112,
    115,
    116,
    119,
    122,
    123,
    124,
    125,
    128,
    134,
    137,
    143,
    149,
    155,
    161,
    164,
    170,
    173,
    176,
    182,
    185,
    191,
    197,
    203,
    209,
    212,
    215,
    218,
    221,
    227,
    229,
    230,
    234,
    239,
    245,
    251,
    257,
    263,
    266,
    269,
    275,
    279,
    284,
    290,
    292,
    298,
    304,
    310,
    313,
    316,
    323,
    334,
    336,
    337,
    338,
    341,
    343,
    349,
    350,
    351,
    352,
    353,
    354,
    355,
    356,
    357,
    358,
    360,
    362,
    363,
    365,
    366,
    367,
    368,
    369,
    370,
    371,
    372,
    373,
    374,
    375,
    376,
    377,
    378,
    379,
    380,
    381,
    382,
    383,
    384,
    385,
    386,
    387,
    388,
    389,
    390,
    391,
    392,
    393,
    394,
    395,
    396,
    397,
    398,
    399,
    400,
    401,
    403,
    405,
    406,
    407,
    409,
    411,
    412,
    413,
    415,
    417,
    418,
    419,
    421,
    423,
    424,
    425,
    427,
    429,
    430,
    431,
    432,
    433,
    435,
    436,
    438,
    439,
    440,
    441,
    442,
    443,
    444,
    446,
    447,
    448,
    449,
    450,
    452,
    454,
    455,
    456,
    457,
    458,
    460,
    462,
    463,
    464,
    465,
    466,
    467,
    468,
    469,
    470,
    471,
    472,
    473,
    474,
    475,
    476,
    477,
    478,
    479,
    480,
    481,
    482,
    483,
    484,
    485,
    486,
    487,
    488,
    489,
    490,
    491,
    492,
    493,
    494,
    496,
    497,
    499,
    500,
    501,
    502,
    503,
    504,
    505,
    506,
    507,
    508,
    509,
    510,
    511,
    512,
    513,
    514,
    515,
    516,
    517,
    519,
    520,
    522,
    523,
    524,
    526,
    528,
    529,
    530,
]
# array position (starting with 1) represent international #
# Default Descriptions (consistent with cryst.ehu.eh):
# unique axis b and cell choice 1 (b1 when there are a1, b1, c1)
# origin choice 2
# obverse setting triple hexagonal axes for R space groups
# assumes unique axis b for (a or b or c) Hall symbols

# TODO: incorporate non-standard setting transformations
def hall_from_hm(num):
    """
    Returns a space group's Hall number from the international number

    Args:
        num: the international space group number

    Returns:
        the Hall number
    """
    return hall_numbers[num]


"""
List of Hermann-Mauguin to Hall symbol conversions:
HM #    HM Symbol      Hall Symbol(different from wyckoff.csv)
 1        P 1            P 1
 2        P -1          -P 1
 3:b      P 1 2 1        P 2y
 3:c      P 1 1 2        P 2
 3:a      P 2 1 1        P 2x
 4:b      P 1 21 1       P 2yb
 4:c      P 1 1 21       P 2c
 4:a      P 21 1 1       P 2xa
 5:b1     C 1 2 1        C 2y
 5:b2     A 1 2 1        A 2y
 5:b3     I 1 2 1        I 2y
 5:c1     A 1 1 2        A 2
 5:c2     B 1 1 2        B 2
 5:c3     I 1 1 2        I 2
 5:a1     B 2 1 1        B 2x
 5:a2     C 2 1 1        C 2x
 5:a3     I 2 1 1        I 2x
 6:b      P 1 m 1        P -2y
 6:c      P 1 1 m        P -2
 6:a      P m 1 1        P -2x
 7:b1     P 1 c 1        P -2yc
 7:b2     P 1 n 1        P -2yac
 7:b3     P 1 a 1        P -2ya
 7:c1     P 1 1 a        P -2a
 7:c2     P 1 1 n        P -2ab
 7:c3     P 1 1 b        P -2b
 7:a1     P b 1 1        P -2xb
 7:a2     P n 1 1        P -2xbc
 7:a3     P c 1 1        P -2xc
 8:b1     C 1 m 1        C -2y
 8:b2     A 1 m 1        A -2y
 8:b3     I 1 m 1        I -2y
 8:c1     A 1 1 m        A -2
 8:c2     B 1 1 m        B -2
 8:c3     I 1 1 m        I -2
 8:a1     B m 1 1        B -2x
 8:a2     C m 1 1        C -2x
 8:a3     I m 1 1        I -2x
 9:b1     C 1 c 1        C -2yc
 9:b2     A 1 n 1        A -2yac
 9:b3     I 1 a 1        I -2ya
 9:-b1    A 1 a 1        A -2ya
 9:-b2    C 1 n 1        C -2ybc
 9:-b3    I 1 c 1        I -2yc
 9:c1     A 1 1 a        A -2a
 9:c2     B 1 1 n        B -2bc
 9:c3     I 1 1 b        I -2b
 9:-c1    B 1 1 b        B -2b
 9:-c2    A 1 1 n        A -2ac
 9:-c3    I 1 1 a        I -2a
 9:a1     B b 1 1        B -2xb
 9:a2     C n 1 1        C -2xbc
 9:a3     I c 1 1        I -2xc
 9:-a1    C c 1 1        C -2xc
 9:-a2    B n 1 1        B -2xbc
 9:-a3    I b 1 1        I -2xb
10:b      P 1 2/m 1     -P 2y
10:c      P 1 1 2/m     -P 2
10:a      P 2/m 1 1     -P 2x
11:b      P 1 21/m 1    -P 2yb
11:c      P 1 1 21/m    -P 2c
11:a      P 21/m 1 1    -P 2xa
12:b1     C 1 2/m 1     -C 2y
12:b2     A 1 2/m 1     -A 2y
12:b3     I 1 2/m 1     -I 2y
12:c1     A 1 1 2/m     -A 2
12:c2     B 1 1 2/m     -B 2
12:c3     I 1 1 2/m     -I 2
12:a1     B 2/m 1 1     -B 2x
12:a2     C 2/m 1 1     -C 2x
12:a3     I 2/m 1 1     -I 2x
13:b1     P 1 2/c 1     -P 2yc
13:b2     P 1 2/n 1     -P 2yac
13:b3     P 1 2/a 1     -P 2ya
13:c1     P 1 1 2/a     -P 2a
13:c2     P 1 1 2/n     -P 2ab
13:c3     P 1 1 2/b     -P 2b
13:a1     P 2/b 1 1     -P 2xb
13:a2     P 2/n 1 1     -P 2xbc
13:a3     P 2/c 1 1     -P 2xc
14:b1     P 1 21/c 1    -P 2ybc
14:b2     P 1 21/n 1    -P 2yn
14:b3     P 1 21/a 1    -P 2yab
14:c1     P 1 1 21/a    -P 2ac
14:c2     P 1 1 21/n    -P 2n
14:c3     P 1 1 21/b    -P 2bc
14:a1     P 21/b 1 1    -P 2xab
14:a2     P 21/n 1 1    -P 2xn
14:a3     P 21/c 1 1    -P 2xac
15:b1     C 1 2/c 1     -C 2yc
15:b2     A 1 2/n 1     -A 2yac
15:b3     I 1 2/a 1     -I 2ya
15:-b1    A 1 2/a 1     -A 2ya
15:-b2    C 1 2/n 1     -C 2ybc
15:-b3    I 1 2/c 1     -I 2yc
15:c1     A 1 1 2/a     -A 2a
15:c2     B 1 1 2/n     -B 2bc
15:c3     I 1 1 2/b     -I 2b
15:-c1    B 1 1 2/b     -B 2b
15:-c2    A 1 1 2/n     -A 2ac
15:-c3    I 1 1 2/a     -I 2a
15:a1     B 2/b 1 1     -B 2xb
15:a2     C 2/n 1 1     -C 2xbc
15:a3     I 2/c 1 1     -I 2xc
15:-a1    C 2/c 1 1     -C 2xc
15:-a2    B 2/n 1 1     -B 2xbc
15:-a3    I 2/b 1 1     -I 2xb
16        P 2 2 2        P 2 2
17        P 2 2 21       P 2c 2
17:cab    P 21 2 2       P 2a 2a
17:bca    P 2 21 2       P 2 2b
18        P 21 21 2      P 2 2ab
18:cab    P 2 21 21      P 2bc 2
18:bca    P 21 2 21      P 2ac 2ac
19        P 21 21 21     P 2ac 2ab
20        C 2 2 21       C 2c 2
20:cab    A 21 2 2       A 2a 2a
20:bca    B 2 21 2       B 2 2b
21        C 2 2 2        C 2 2
21:cab    A 2 2 2        A 2 2
21:bca    B 2 2 2        B 2 2
22        F 2 2 2        F 2 2
23        I 2 2 2        I 2 2
24        I 21 21 21     I 2b 2c
25        P m m 2        P 2 -2
25:cab    P 2 m m        P -2 2
25:bca    P m 2 m        P -2 -2
26        P m c 21       P 2c -2
26:ba-c   P c m 21       P 2c -2c
26:cab    P 21 m a       P -2a 2a
26:-cba   P 21 a m       P -2 2a
26:bca    P b 21 m       P -2 -2b
26:a-cb   P m 21 b       P -2b -2
27        P c c 2        P 2 -2c
27:cab    P 2 a a        P -2a 2
27:bca    P b 2 b        P -2b -2b
28        P m a 2        P 2 -2a
28:ba-c   P b m 2        P 2 -2b
28:cab    P 2 m b        P -2b 2
28:-cba   P 2 c m        P -2c 2
28:bca    P c 2 m        P -2c -2c
28:a-cb   P m 2 a        P -2a -2a
29        P c a 21       P 2c -2ac
29:ba-c   P b c 21       P 2c -2b
29:cab    P 21 a b       P -2b 2a
29:-cba   P 21 c a       P -2ac 2a
29:bca    P c 21 b       P -2bc -2c
29:a-cb   P b 21 a       P -2a -2ab
30        P n c 2        P 2 -2bc
30:ba-c   P c n 2        P 2 -2ac
30:cab    P 2 n a        P -2ac 2
30:-cba   P 2 a n        P -2ab 2
30:bca    P b 2 n        P -2ab -2ab
30:a-cb   P n 2 b        P -2bc -2bc
31        P m n 21       P 2ac -2
31:ba-c   P n m 21       P 2bc -2bc
31:cab    P 21 m n       P -2ab 2ab
31:-cba   P 21 n m       P -2 2ac
31:bca    P n 21 m       P -2 -2bc
31:a-cb   P m 21 n       P -2ab -2
32        P b a 2        P 2 -2ab
32:cab    P 2 c b        P -2bc 2
32:bca    P c 2 a        P -2ac -2ac
33        P n a 21       P 2c -2n
33:ba-c   P b n 21       P 2c -2ab
33:cab    P 21 n b       P -2bc 2a
33:-cba   P 21 c n       P -2n 2a
33:bca    P c 21 n       P -2n -2ac
33:a-cb   P n 21 a       P -2ac -2n
34        P n n 2        P 2 -2n
34:cab    P 2 n n        P -2n 2
34:bca    P n 2 n        P -2n -2n
35        C m m 2        C 2 -2
35:cab    A 2 m m        A -2 2
35:bca    B m 2 m        B -2 -2
36        C m c 21       C 2c -2
36:ba-c   C c m 21       C 2c -2c
36:cab    A 21 m a       A -2a 2a
36:-cba   A 21 a m       A -2 2a
36:bca    B b 21 m       B -2 -2b
36:a-cb   B m 21 b       B -2b -2
37        C c c 2        C 2 -2c
37:cab    A 2 a a        A -2a 2
37:bca    B b 2 b        B -2b -2b
38        A m m 2        A 2 -2
38:ba-c   B m m 2        B 2 -2
38:cab    B 2 m m        B -2 2
38:-cba   C 2 m m        C -2 2
38:bca    C m 2 m        C -2 -2
38:a-cb   A m 2 m        A -2 -2
39        A b m 2        A 2 -2c
39:ba-c   B m a 2        B 2 -2c
39:cab    B 2 c m        B -2c 2
39:-cba   C 2 m b        C -2b 2
39:bca    C m 2 a        C -2b -2b
39:a-cb   A c 2 m        A -2c -2c
40        A m a 2        A 2 -2a
40:ba-c   B b m 2        B 2 -2b
40:cab    B 2 m b        B -2b 2
40:-cba   C 2 c m        C -2c 2
40:bca    C c 2 m        C -2c -2c
40:a-cb   A m 2 a        A -2a -2a
41        A b a 2        A 2 -2ac
41:ba-c   B b a 2        B 2 -2bc
41:cab    B 2 c b        B -2bc 2
41:-cba   C 2 c b        C -2bc 2
41:bca    C c 2 a        C -2bc -2bc
41:a-cb   A c 2 a        A -2ac -2ac
42        F m m 2        F 2 -2
42:cab    F 2 m m        F -2 2
42:bca    F m 2 m        F -2 -2
43        F d d 2        F 2 -2d
43:cab    F 2 d d        F -2d 2
43:bca    F d 2 d        F -2d -2d
44        I m m 2        I 2 -2
44:cab    I 2 m m        I -2 2
44:bca    I m 2 m        I -2 -2
45        I b a 2        I 2 -2c
45:cab    I 2 c b        I -2a 2
45:bca    I c 2 a        I -2b -2b
46        I m a 2        I 2 -2a
46:ba-c   I b m 2        I 2 -2b
46:cab    I 2 m b        I -2b 2
46:-cba   I 2 c m        I -2c 2
46:bca    I c 2 m        I -2c -2c
46:a-cb   I m 2 a        I -2a -2a
47        P m m m       -P 2 2
48:1      P n n n:1      P 2 2 -1n
48:2      P n n n:2     -P 2ab 2bc
49        P c c m       -P 2 2c
49:cab    P m a a       -P 2a 2
49:bca    P b m b       -P 2b 2b
50:1      P b a n:1      P 2 2 -1ab
50:2      P b a n:2     -P 2ab 2b
50:1cab   P n c b:1      P 2 2 -1bc
50:2cab   P n c b:2     -P 2b 2bc
50:1bca   P c n a:1      P 2 2 -1ac
50:2bca   P c n a:2     -P 2a 2c
51        P m m a       -P 2a 2a
51:ba-c   P m m b       -P 2b 2
51:cab    P b m m       -P 2 2b
51:-cba   P c m m       -P 2c 2c
51:bca    P m c m       -P 2c 2
51:a-cb   P m a m       -P 2 2a
52        P n n a       -P 2a 2bc
52:ba-c   P n n b       -P 2b 2n
52:cab    P b n n       -P 2n 2b
52:-cba   P c n n       -P 2ab 2c
52:bca    P n c n       -P 2ab 2n
52:a-cb   P n a n       -P 2n 2bc
53        P m n a       -P 2ac 2
53:ba-c   P n m b       -P 2bc 2bc
53:cab    P b m n       -P 2ab 2ab
53:-cba   P c n m       -P 2 2ac
53:bca    P n c m       -P 2 2bc
53:a-cb   P m a n       -P 2ab 2
54        P c c a       -P 2a 2ac
54:ba-c   P c c b       -P 2b 2c
54:cab    P b a a       -P 2a 2b
54:-cba   P c a a       -P 2ac 2c
54:bca    P b c b       -P 2bc 2b
54:a-cb   P b a b       -P 2b 2ab
55        P b a m       -P 2 2ab
55:cab    P m c b       -P 2bc 2
55:bca    P c m a       -P 2ac 2ac
56        P c c n       -P 2ab 2ac
56:cab    P n a a       -P 2ac 2bc
56:bca    P b n b       -P 2bc 2ab
57        P b c m       -P 2c 2b
57:ba-c   P c a m       -P 2c 2ac
57:cab    P m c a       -P 2ac 2a
57:-cba   P m a b       -P 2b 2a
57:bca    P b m a       -P 2a 2ab
57:a-cb   P c m b       -P 2bc 2c
58        P n n m       -P 2 2n
58:cab    P m n n       -P 2n 2
58:bca    P n m n       -P 2n 2n
59:1      P m m n:1      P 2 2ab -1ab
59:2      P m m n:2     -P 2ab 2a
59:1cab   P n m m:1      P 2bc 2 -1bc
59:2cab   P n m m:2     -P 2c 2bc
59:1bca   P m n m:1      P 2ac 2ac -1ac
59:2bca   P m n m:2     -P 2c 2a
60        P b c n       -P 2n 2ab
60:ba-c   P c a n       -P 2n 2c
60:cab    P n c a       -P 2a 2n
60:-cba   P n a b       -P 2bc 2n
60:bca    P b n a       -P 2ac 2b
60:a-cb   P c n b       -P 2b 2ac
61        P b c a       -P 2ac 2ab
61:ba-c   P c a b       -P 2bc 2ac
62        P n m a       -P 2ac 2n
62:ba-c   P m n b       -P 2bc 2a
62:cab    P b n m       -P 2c 2ab
62:-cba   P c m n       -P 2n 2ac
62:bca    P m c n       -P 2n 2a
62:a-cb   P n a m       -P 2c 2n
63        C m c m       -C 2c 2
63:ba-c   C c m m       -C 2c 2c
63:cab    A m m a       -A 2a 2a
63:-cba   A m a m       -A 2 2a
63:bca    B b m m       -B 2 2b
63:a-cb   B m m b       -B 2b 2
64        C m c a       -C 2bc 2
64:ba-c   C c m b       -C 2bc 2bc
64:cab    A b m a       -A 2ac 2ac
64:-cba   A c a m       -A 2 2ac
64:bca    B b c m       -B 2 2bc
64:a-cb   B m a b       -B 2bc 2
65        C m m m       -C 2 2
65:cab    A m m m       -A 2 2
65:bca    B m m m       -B 2 2
66        C c c m       -C 2 2c
66:cab    A m a a       -A 2a 2
66:bca    B b m b       -B 2b 2b
67        C m m a       -C 2b 2
67:ba-c   C m m b       -C 2b 2b
67:cab    A b m m       -A 2c 2c
67:-cba   A c m m       -A 2 2c
67:bca    B m c m       -B 2 2c
67:a-cb   B m a m       -B 2c 2
68:1      C c c a:1      C 2 2 -1bc
68:2      C c c a:2     -C 2b 2bc
68:1ba-c  C c c b:1      C 2 2 -1bc
68:2ba-c  C c c b:2     -C 2b 2c
68:1cab   A b a a:1      A 2 2 -1ac
68:2cab   A b a a:2     -A 2a 2c
68:1-cba  A c a a:1      A 2 2 -1ac
68:2-cba  A c a a:2     -A 2ac 2c
68:1bca   B b c b:1      B 2 2 -1bc
68:2bca   B b c b:2     -B 2bc 2b
68:1a-cb  B b a b:1      B 2 2 -1bc
68:2a-cb  B b a b:2     -B 2b 2bc
69        F m m m       -F 2 2
70:1      F d d d:1      F 2 2 -1d
70:2      F d d d:2     -F 2uv 2vw
71        I m m m       -I 2 2
72        I b a m       -I 2 2c
72:cab    I m c b       -I 2a 2
72:bca    I c m a       -I 2b 2b
73        I b c a       -I 2b 2c
73:ba-c   I c a b       -I 2a 2b
74        I m m a       -I 2b 2
74:ba-c   I m m b       -I 2a 2a
74:cab    I b m m       -I 2c 2c
74:-cba   I c m m       -I 2 2b
74:bca    I m c m       -I 2 2a
74:a-cb   I m a m       -I 2c 2
75        P 4            P 4
76        P 41           P 4w
77        P 42           P 4c
78        P 43           P 4cw
79        I 4            I 4
80        I 41           I 4bw
81        P -4           P -4
82        I -4           I -4
83        P 4/m         -P 4
84        P 42/m        -P 4c
85:1      P 4/n:1        P 4ab -1ab
85:2      P 4/n:2       -P 4a
86:1      P 42/n:1       P 4n -1n
86:2      P 42/n:2      -P 4bc
87        I 4/m         -I 4
88:1      I 41/a:1       I 4bw -1bw
88:2      I 41/a:2      -I 4ad
89        P 4 2 2        P 4 2
90        P 42 1 2       P 4ab 2ab
91        P 41 2 2       P 4w 2c
92        P 41 21 2      P 4abw 2nw
93        P 42 2 2       P 4c 2
94        P 42 21 2      P 4n 2n
95        P 43 2 2       P 4cw 2c
96        P 43 21 2      P 4nw 2abw
97        I 4 2 2        I 4 2
98        I 41 2 2       I 4bw 2bw
99        P 4 m m        P 4 -2
100        P 4 b m        P 4 -2ab
101        P 42 c m       P 4c -2c
102        P 42 n m       P 4n -2n
103        P 4 c c        P 4 -2c
104        P 4 n c        P 4 -2n
105        P 42 m c       P 4c -2
106        P 42 b c       P 4c -2ab
107        I 4 m m        I 4 -2
108        I 4 c m        I 4 -2c
109        I 41 m d       I 4bw -2
110        I 41 c d       I 4bw -2c
111        P -4 2 m       P -4 2
112        P -4 2 c       P -4 2c
113        P -4 21 m      P -4 2ab
114        P -4 21 c      P -4 2n
115        P -4 m 2       P -4 -2
116        P -4 c 2       P -4 -2c
117        P -4 b 2       P -4 -2ab
118        P -4 n 2       P -4 -2n
119        I -4 m 2       I -4 -2
120        I -4 c 2       I -4 -2c
121        I -4 2 m       I -4 2
122        I -4 2 d       I -4 2bw
123        P 4/m m m     -P 4 2
124        P 4/m c c     -P 4 2c
125:1      P 4/n b m:1    P 4 2 -1ab
125:2      P 4/n b m:2   -P 4a 2b
126:1      P 4/n n c:1    P 4 2 -1n
126:2      P 4/n n c:2   -P 4a 2bc
127        P 4/m b m     -P 4 2ab
128        P 4/m n c     -P 4 2n
129:1      P 4/n m m:1    P 4ab 2ab -1ab
129:2      P 4/n m m:2   -P 4a 2a
130:1      P 4/n c c:1    P 4ab 2n -1ab
130:2      P 4/n c c:2   -P 4a 2ac
131        P 42/m m c    -P 4c 2
132        P 42/m c m    -P 4c 2c
133:1      P 42/n b c:1   P 4n 2c -1n
133:2      P 42/n b c:2  -P 4ac 2b
134:1      P 42/n n m:1   P 4n 2 -1n
134:2      P 42/n n m:2  -P 4ac 2bc
135        P 42/m b c    -P 4c 2ab
136        P 42/m n m    -P 4n 2n
137:1      P 42/n m c:1   P 4n 2n -1n
137:2      P 42/n m c:2  -P 4ac 2a
138:1      P 42/n c m:1   P 4n 2ab -1n
138:2      P 42/n c m:2  -P 4ac 2ac
139        I 4/m m m     -I 4 2
140        I 4/m c m     -I 4 2c
141:1      I 41/a m d:1   I 4bw 2bw -1bw
141:2      I 41/a m d:2  -I 4bd 2
142:1      I 41/a c d:1   I 4bw 2aw -1bw
142:2      I 41/a c d:2  -I 4bd 2c
143        P 3            P 3
144        P 31           P 31
145        P 32           P 32
146:H      R 3:H          R 3
146:R      R 3:R          P 3*
147        P -3          -P 3
148:H      R -3:H        -R 3
148:R      R -3:R        -P 3*
149        P 3 1 2        P 3 2
150        P 3 2 1        P 3 2"
151        P 31 1 2       P 31 2c (0 0 1)
152        P 31 2 1       P 31 2"
153        P 32 1 2       P 32 2c (0 0 -1)
154        P 32 2 1       P 32 2"
155:H      R 32:H         R 3 2"
155:R      R 32:R         P 3* 2
156        P 3 m 1        P 3 -2"
157        P 3 1 m        P 3 -2
158        P 3 c 1        P 3 -2"c
159        P 3 1 c        P 3 -2c
160:H      R 3 m:H        R 3 -2"
160:R      R 3 m:R        P 3* -2
161:H      R 3 c:H        R 3 -2"c
161:R      R 3 c:R        P 3* -2n
162        P -3 1 m      -P 3 2
163        P -3 1 c      -P 3 2c
164        P -3 m 1      -P 3 2"
165        P -3 c 1      -P 3 2"c
166:H      R -3 m:H      -R 3 2"
166:R      R -3 m:R      -P 3* 2
167:H      R -3 c:H      -R 3 2"c
167:R      R -3 c:R      -P 3* 2n
168        P 6            P 6
169        P 61           P 61
170        P 65           P 65
171        P 62           P 62
172        P 64           P 64
173        P 63           P 6c
174        P -6           P -6
175        P 6/m         -P 6
176        P 63/m        -P 6c
177        P 6 2 2        P 6 2
178        P 61 2 2       P 61 2 (0 0 -1)
179        P 65 2 2       P 65 2 (0 0 1)
180        P 62 2 2       P 62 2c (0 0 1)
181        P 64 2 2       P 64 2c (0 0 -1)
182        P 63 2 2       P 6c 2c
183        P 6 m m        P 6 -2
184        P 6 c c        P 6 -2c
185        P 63 c m       P 6c -2
186        P 63 m c       P 6c -2c
187        P -6 m 2       P -6 2
188        P -6 c 2       P -6c 2
189        P -6 2 m       P -6 -2
190        P -6 2 c       P -6c -2c
191        P 6/m m m     -P 6 2
192        P 6/m c c     -P 6 2c
193        P 63/m c m    -P 6c 2
194        P 63/m m c    -P 6c 2c
195        P 2 3          P 2 2 3
196        F 2 3          F 2 2 3
197        I 2 3          I 2 2 3
198        P 21 3         P 2ac 2ab 3
199        I 21 3         I 2b 2c 3
200        P m -3        -P 2 2 3
201:1      P n -3:1       P 2 2 3 -1n
201:2      P n -3:2      -P 2ab 2bc 3
202        F m -3        -F 2 2 3
203:1      F d -3:1       F 2 2 3 -1d
203:2      F d -3:2      -F 2uv 2vw 3
204        I m -3        -I 2 2 3
205        P a -3        -P 2ac 2ab 3
206        I a -3        -I 2b 2c 3
207        P 4 3 2        P 4 2 3
208        P 42 3 2       P 4n 2 3
209        F 4 3 2        F 4 2 3
210        F 41 3 2       F 4d 2 3
211        I 4 3 2        I 4 2 3
212        P 43 3 2       P 4acd 2ab 3
213        P 41 3 2       P 4bd 2ab 3
214        I 41 3 2       I 4bd 2c 3
215        P -4 3 m       P -4 2 3
216        F -4 3 m       F -4 2 3
217        I -4 3 m       I -4 2 3
218        P -4 3 n       P -4n 2 3
219        F -4 3 c       F -4c 2 3
220        I -4 3 d       I -4bd 2c 3
221        P m -3 m      -P 4 2 3
222:1      P n -3 n:1     P 4 2 3 -1n
222:2      P n -3 n:2    -P 4a 2bc 3
223        P m -3 n      -P 4n 2 3
224:1      P n -3 m:1     P 4n 2 3 -1n
224:2      P n -3 m:2    -P 4bc 2bc 3
225        F m -3 m      -F 4 2 3
226        F m -3 c      -F 4c 2 3
227:1      F d -3 m:1     F 4d 2 3 -1d
227:2      F d -3 m:2    -F 4vw 2vw 3
228:1      F d -3 c:1     F 4d 2 3 -1cd
228:2      F d -3 c:2    -F 4cvw 2vw 3
229        I m -3 m      -I 4 2 3
230        I a -3 d      -I 4bd 2c 3
"""
