"""
Module for handling layer groups for generation of 2d crystals. A layer group is
the combination of a crystallographic point group with 2 directions of periodic
translational symmetry. There are 80 layer groups, each of which corresponds to
a 3d space group, with a possible permutation of the axes.
"""
from optparse import OptionParser

class Layergroup:
    """
    Class for storing layer groups. Used for the generation of 2d crystals.
    Each layer group corresponds to the point group of a 3d space group. Thus, we
    can access the symmetry information using code for standard space groups. The
    primary class, Layergroup, allows access to the layer group number, space group
    number, and space group symbol.

    Args:
        input_value: The layer group number or symbol
    """
    def __init__(self, input_value):

        # list of layer group number, symbol, space group number and permutation
        # international b-axis is the speical axis
        self.group_list = [
            (1,     'P1',       1,      [1,2,3,1,2,3,3]),
            (2,     'P-1',      2,      [1,2,3,1,2,3,3]),
            (3,     'P112',     3,      [1,3,2,1,3,2,2]),
            (4,     'P11m',     6,      [1,3,2,1,3,2,2]),
            (5,     'P11a',     7,      [1,3,2,1,3,2,2]),
            (6,     'P112/m',   10,     [1,3,2,1,3,2,2]),
            (7,     'P112/a',   13,     [1,3,2,1,3,2,2]),
            (8,     'P211',     3,      [3,1,2,2,3,1,1]),
            (9,     'P2111',    4,      [3,1,2,2,3,1,1]),
            (10,    'C211',     5,      [2,1,3,2,1,3,3]),
            (11,    'Pm11',     6,      [3,1,2,2,3,1,1]),
            (12,    'Pb11',     7,      [3,1,2,2,3,1,1]),
            (13,    'Cm11',     8,      [2,1,3,2,1,3,3]),
            (14,    'P2/m11',   10,     [3,1,2,2,3,1,1]),
            (15,    'P21/m11',  11,     [3,1,2,2,3,1,1]),
            (16,    'P2/b11',   13,     [3,1,2,2,3,1,1]),
            (17,    'P21/b11',  14,     [3,1,2,2,3,1,1]),
            (18,    'C2/m11',   12,     [2,1,3,2,1,3,3]),
            (19,    'P222',     16,     [1,2,3,1,2,3,3]),
            (20,    'P2122',    17,     [3,2,1,3,2,1,1]),
            (21,    'P21211',   18,     [1,2,3,1,2,3,3]),
            (22,    'C222',     21,     [1,2,3,1,2,3,3]),
            (23,    'Pmm2',     25,     [1,2,3,1,2,3,3]),
            (24,    'Pma2',     28,     [1,2,3,1,2,3,3]),
            (25,    'Pba2',     32,     [1,2,3,1,2,3,3]),
            (26,    'Cmm2',     35,     [1,2,3,1,2,3,3]),
            (27,    'Pm2m',     25,     [1,3,2,1,3,2,2]),
            (28,    'Pm21b',    26,     [1,3,2,1,3,2,2]),
            (29,    'Pb21m',    26,     [3,1,2,2,3,1,1]),
            (30,    'Pb2b',     27,     [1,3,2,1,3,2,2]),
            (31,    'Pm2a',     28,     [1,3,2,1,3,2,2]),
            (32,    'Pm21n',    31,     [1,3,2,1,3,2,2]),
            (33,    'Pb21a',    29,     [1,3,2,1,3,2,2]),
            (34,    'Pb2n',     30,     [3,1,2,2,3,1,1]),
            (35,    'Cm2m',     38,     [1,2,3,1,2,3,3]),
            (36,    'Cm2e',     39,     [1,2,3,1,2,3,3]),
            (37,    'Pmmm',     47,     [1,2,3,1,2,3,3]),
            (38,    'Pmaa',     49,     [3,2,1,3,2,1,1]),
            (39,    'Pban',     50,     [1,2,3,1,2,3,3]),
            (40,    'Pmam',     51,     [1,3,2,1,3,2,2]),
            (41,    'Pmma',     51,     [1,2,3,1,2,3,3]),
            (42,    'Pman',     53,     [1,3,2,1,3,2,2]),
            (43,    'Pbaa',     54,     [1,3,2,1,3,2,2]),
            (44,    'Pbam',     55,     [1,2,3,1,2,3,3]),
            (45,    'Pbma',     57,     [3,2,1,3,2,1,1]),
            (46,    'Pmmn',     59,     [1,2,3,1,2,3,3]),
            (47,    'Cmmm',     65,     [1,2,3,1,2,3,3]),
            (48,    'Cmme',     67,     [1,2,3,1,2,3,3]),
            (49,    'P4',       75,     [1,2,3,1,2,3,3]),
            (50,    'P-4',      81,     [1,2,3,1,2,3,3]),
            (51,    'P4/m',     83,     [1,2,3,1,2,3,3]),
            (52,    'P4/n',     85,     [1,2,3,1,2,3,3]),
            (53,    'P422',     89,     [1,2,3,1,2,3,3]),
            (54,    'P4212',    90,     [1,2,3,1,2,3,3]),
            (55,    'P4mm',     99,     [1,2,3,1,2,3,3]),
            (56,    'P4bm',     100,    [1,2,3,1,2,3,3]),
            (57,    'P-42m',    111,    [1,2,3,1,2,3,3]),
            (58,    'P-421m',   113,    [1,2,3,1,2,3,3]),
            (59,    'P-4m2',    115,    [1,2,3,1,2,3,3]),
            (60,    'P-4b2',    117,    [1,2,3,1,2,3,3]),
            (61,    'P4/mmmm',  123,    [1,2,3,1,2,3,3]),
            (62,    'P4/nbm',   125,    [1,2,3,1,2,3,3]),
            (63,    'P4/mbm',   127,    [1,2,3,1,2,3,3]),
            (64,    'P4/nmm',   129,    [1,2,3,1,2,3,3]),
            (65,    'P3',       143,    [1,2,3,1,2,3,3]),
            (66,    'P-3',      147,    [1,2,3,1,2,3,3]),
            (67,    'P312',     149,    [1,2,3,1,2,3,3]),
            (68,    'P321',     150,    [1,2,3,1,2,3,3]),
            (69,    'P3m1',     156,    [1,2,3,1,2,3,3]),
            (70,    'P31m',     157,    [1,2,3,1,2,3,3]),
            (71,    'P-31m',    162,    [1,2,3,1,2,3,3]),
            (72,    'P-3m1',    164,    [1,2,3,1,2,3,3]),
            (73,    'P6',       168,    [1,2,3,1,2,3,3]),
            (74,    'P-6',      174,    [1,2,3,1,2,3,3]),
            (75,    'P6/m',     175,    [1,2,3,1,2,3,3]),
            (76,    'P622',     177,    [1,2,3,1,2,3,3]),
            (77,    'P6mm',     183,    [1,2,3,1,2,3,3]),
            (78,    'P-6m2',    187,    [1,2,3,1,2,3,3]),
            (79,    'P-62m',    189,    [1,2,3,1,2,3,3]),
            (80,    'P6/mmm',   191,    [1,2,3,1,2,3,3]),
            ]

        self.input = str(input_value)
        self.error = False
        if self.input.isdigit():
            self.lg = int(self.input)
            """The layer group number (between 1 and 80)"""
        else:
            for i, e1 in enumerate(self.group_list):
                if e1[1] == self.input:
                    self.lg = e1[0]
                    break
        if (self.lg is not None) and (0<self.lg<81):
            self.symbol=self.group_list[self.lg-1][1]
            """The Hermann-Mauguin symbol of the layer group"""
            self.sgnumber = self.group_list[self.lg-1][2]
            """The international space group number (between 1 and 230)"""
            self.permutation=self.group_list[self.lg-1][3]
        else:
            self.error = True
            print('Error:   unable to find the layer group, check your input: ', self.input)

    def print_all(self):
        if self.error is False:
            print('Layer group number: ', self.lg)
            print('Layer group symbol: ', self.symbol)
            print('Space group number: ', self.sgnumber)
            print('Permutation:        ', self.permutation)
#LG.Pointer = [ ...
#              1,   2,   3,   6,   7,  10,  13,   3,   4,   5, ...
#              6,   7,   8,  10,  11,  13,  14,  12,  16,  17, ...
#             18,  21,  25,  28,  32,  35,  25,  26,  26,  27, ...
#             28,  31,  29,  30,  38,  39,  47,  49,  50,  51, ...
#             51,  53,  54,  55,  57,  59,  65,  67,  75,  81, ...
#             83,  85,  89,  90,  99, 100, 111, 113, 115, 117, ...
#            123, 125, 127, 129, 143, 147, 149, 150, 156, 157, ...
#            162, 164, 168, 174, 175, 177, 183, 187, 189, 191, ...

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", metavar='sg', 
            help="input number or symbol")

    (options, args) = parser.parse_args()    
    Layergroup(options.input).print_all()









