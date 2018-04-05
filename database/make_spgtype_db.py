#!/usr/bin/env python

import sys

def read_spg_csv(filename="spg.csv"):
    spg_db = []
    for line in open(filename):
        data = line.split(',')
        spg_db.append(data)
    return spg_db

def extract_essense(spg_db):
    essense = []
    for x in spg_db:
        vals = []
        # IT number
        vals.append(int(x[4]))
        # SF
        vals.append(x[5])
        # Hall
        vals.append(x[6])
        # HM short all
        vals.append(x[7])
        # HM long
        vals.append(x[8])
        # HM short
        vals.append(x[7].split('=')[0].replace(' ',''))
        # Setting
        vals.append(x[2])
        # Centering
        vals.append(get_centering(x[8][0], ("%-9s" % x[2])[0]))
        essense.append(vals)
    return essense

def get_bravais(n, s):
    if 1 <= n and n <= 2:
        return 'triclinic'
    if 3 <= n and n <= 15:
        return 'monoclinic'
    if 16 <= n and n <= 74:
        return 'orthorhombic'
    if 75 <= n and n <= 142:
        return 'tetragonal'
    if 143 <= n and n <= 167:
        if is_rhombo( s ):
            return 'rhombohedral'
        else:
            return 'trigonal'
    if 168 <= n and n <= 194:
        return 'hexagonal'
    if 195 <= n and n <= 230:
        return 'cubic'

def get_pointgroup_number(spg_num):
    spg_num_head = [1, 2, 3, 6, 10, 16, 25, 47, 75, 81,
                    83, 89, 99, 111, 123, 143, 147, 149, 156, 162,
                    168, 174, 175, 177, 183, 187, 191, 195, 200, 207,
                    215, 221]
    pt_num = 0
    for n in spg_num_head:
        if spg_num >= n:
            pt_num += 1
            
    return pt_num

def get_centering(centering, setting):
    if centering == 'I':
        return "BODY"
    elif centering == 'F':
        return "FACE"
    elif centering == 'A':
        return "A_FACE"
    elif centering == 'C':
        return "C_FACE"
    elif centering == 'R':
        if setting == 'H':
            return "R_CENTER"
        else:
            return "PRIMITIVE"
    else:
        return "PRIMITIVE"

def is_rhombo( s ):
    if s=='R':
        return True

#maxlen = 0
#
#bravais = {'triclinic': 'TRICLI',
#           'monoclinic': 'MONOCLI',
#           'orthorhombic': 'ORTHO',
#           'tetragonal': 'TETRA',
#           'trigonal': 'TRIGO',
#           'rhombohedral': 'RHOMB',
#           'hexagonal': 'HEXA',
#           'cubic': 'CUBIC'}
#
#print "static const SpacegroupType spacegroup_types[] = {"
#
#print "  {%3d, \"%-6s\", \"%-16s\", \"%-31s\", \"%-19s\", \"%-10s\", \"%-5s\", %15s, %2d }, /* %3d */" % (0, "", "", "", "", "", "", "CENTERING_ERROR", 0, 0)
#
#for i, x in enumerate(extract_essense(read_spg_csv(sys.argv[1]))):
#    if len(x[6]) > maxlen:
#        maxlen = len(x[6])
#
#    pt_num = get_pointgroup_number(x[0])
#
#    print "  {%3d, \"%-6s\", \"%-16s\", \"%-31s\", \"%-19s\", \"%-10s\", \"%-5s\", %15s, %2d }, /* %3d */" % (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], pt_num, i + 1)
#print "};"
