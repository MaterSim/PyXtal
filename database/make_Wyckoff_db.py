#!/usr/bin/env python
'''
Code from https://github.com/atztogo/spglib/tree/master/database
Originally created by Atsushi Togo.
Modified by Scott Fredericks for python 3 compatibility.
'''

import sys
import numpy as np

lattice_symbols = {
    'P': [ [ 0, 0, 0 ] ],
    'A': [ [ 0, 0, 0 ], [ 0, 12, 12 ] ],
    'B': [ [ 0, 0, 0 ], [ 12, 0, 12 ] ],
    'C': [ [ 0, 0, 0 ], [ 12, 12, 0 ] ],
    'I': [ [ 0, 0, 0 ], [ 12, 12, 12 ] ],
    'R': [ [ 0, 0, 0 ] ],
    'H': [ [ 0, 0, 0 ], [ 16, 8, 8 ], [ 8, 16, 16 ] ],
    'F': [ [ 0, 0, 0 ], [ 0, 12, 12 ], [ 12, 0, 12 ], [ 12, 12, 0 ] ]
    }

def position2operation( p ):
    r = np.zeros( (3,3), dtype=int )
    xyz = 'xyz'

    for i in range(3):
        for j, c in enumerate( xyz ):
            pos = p[i].find(c)
            # if exist xyz
            if pos > -1:
                if pos == 0:
                    r[i,j] = 1
                if pos > 0:
                    if p[i][pos-1] == '+': # +y, +z
                        r[i,j] = 1
                    elif p[i][pos-1] == '-': # -x, -y, -z
                        r[i,j] = -1
                    elif p[i][pos-1].isdigit(): # 2x
                        if pos > 1:
                            if p[i][pos-2] == '-':
                                r[i,j] = - int( p[i][pos-1] )
                            else:
                                r[i,j] = int( p[i][pos-1] )
                        else:
                            r[i,j] = int( p[i][pos-1] )
                    else:
                        print("error")
                        sys.exit(1)

                # print p, p[i], r[i,j]

    t = np.zeros( 3, dtype=int )

    # translation multiplied by 24
    for i in range(3):
        #trans = p[i].translate(None,'xyz+-') #updated for Python 3
        trans = p[i].translate({ord(c): None for c in 'xyz+-'})
        pos = trans.find('/')
        if pos > -1:
            frac = trans.split('/')
            t[i] = int( frac[0] ) * ( 24 / int( frac[1] ) )

    return r,t

def encode_rotation( r ):
    enc = []
    for i in range(3):
        enc.append( ( r[i][0] + 2 ) * 9 + ( r[i][1] + 1 ) * 3 + r[i][2] + 1 )
    return enc[0] * 45**2 + enc[1] * 45 + enc[2]
    
def decode_rotation( c ):
    r = np.zeros( ( 3, 3 ), dtype=int )
    dec = [ c // 45**2, ( c % 45**2 ) // 45, c % 45 ]
    for i in range( 3 ):
        r[i][0] = dec[i] // 9 - 2
        r[i][1] = ( dec[i] % 9 ) // 3 - 1
        r[i][2] = dec[i] % 3 - 1
    return r

def decode_trans( c ):
    return c // 576, ( c % 576 ) // 24, ( c % 24 )

def run_test( len_wyckoffs, len_site, o_flat ):
    count = 1
    sum_elem = 1
    for i in range( 1, 531 ):
        sum_site = 0
        for j in range( len_wyckoffs[i], len_wyckoffs[i+1] ):
            print(len_site[j+1] - len_site[j],)
            
            sum_site += len_site[j+1] - len_site[j]
    
            for x in o_flat[len_site[j]:len_site[j+1]]:
                t = decode_trans( x // 3**9 )
                r = decode_rotation( x % 3**9 )
                print(r, t)
            
        print("  %d" % sum_site)
        sum_elem += sum_site
    return sum_elem

def damp_array_wyckoffs( len_wyckoff ):
    print("static const int position_wyckoff[] =")
    print("  { %4d," % 0,)
    for i, x in enumerate( len_wyckoffs[1:] ):
        if i%10 == 0:
            print()
            print("   ",)
        print("%4d," % x,)
    print("};")

def damp_array_numsites( len_site ):
    print("  { %5d," % 0,)
    for i, x in enumerate( len_site[1:] ):
        if i%10 == 0:
            print()
            print("   ",)
        print("%5d," % x,)
    print("};")
    
def damp_array_positions( orbits ):
    print("  { %9d," % 0,)
    for i, x in enumerate( orbits[1:] ):
        if i%5 == 0:
            print()
            print("   ",)
        print("%9d," % x,)
    print("};")

def damp_array_positions_short( orbits, len_site ):
    count = 0
    print("static const int coordinates_first[] =")
    print("  { %9d," % 0,)
    for i, x in enumerate( len_site[1:-1] ):
        if count%5 == 0:
            print()
            print("   ",)
        print("%9d," % orbits[x],)
        count += 1
    print("};")

    count = 0
    print("static const int num_sitesym[] =")
    print("  { %3d," % 0,)
    for i, x in enumerate( len_site[1:-1] ):
        if count%10 == 0:
            print()
            print("   ",)
        print("%3d," % ( len_site[i+2] - len_site[i+1] ),)
        count += 1
    print("};")
    
def read_wyckoff_csv( filename ):
    wyckoff_file = open( filename )
    rowdata = []
    points = []
    hP_nums = [433, 436, 444, 450, 452, 458, 460]
    for i, line in enumerate( wyckoff_file ):
        if line.strip() == 'end of data':
            break
        rowdata.append( line.strip().split(':') )

        # Remember the point where the first element is a number
        if rowdata[-1][0].isdigit():
            points.append(i)
    points.append(i)

    wyckoff = []
    for i in range( len( points ) - 1 ):
        symbol = rowdata[points[i]][1]
        if i + 1 in hP_nums:
            symbol = symbol.replace('R', 'H', 1)

        wyckoff.append(dict({'symbol': symbol}))
    
    # When the number of positions is larger than 4,
    # the positions are written in the next line.
    # So those positions are connected.
    for i in range( len( points ) - 1 ):
        count = 0
        wyckoff[i]['wyckoff'] = []
        for j in range( points[i]+1, points[i+1] ):
            # Hook if the third colum is a number (multiplicity)
            if rowdata[j][2].isdigit():
                pos = []
                w = { 'letter': rowdata[j][3],
                      'multiplicity': int( rowdata[j][2] ),
                      'site_symmetry': rowdata[j][4],
                      'positions': pos }
                wyckoff[i]['wyckoff'].append( w )
    
                for k in range( 4 ):
                    if len(rowdata[j][k+5]) > 0:
                        count += 1
                        pos.append( rowdata[j][k+5] )
    
            else:
                for k in range( 4 ):
                    if len(rowdata[j][k+5]) > 0:
                        count += 1
                        pos.append( rowdata[j][k+5] )

    return wyckoff

def get_wyckoff_positions( wyckoff ):
    positions = []
    for i, w in enumerate( wyckoff ):
        # print "%d:%s" % (i+1, w['symbol'])
        positions_hall = []
        for p in w['wyckoff']:
            # print "  %d %s %s" % ( p['multiplicity'],
            #                        p['letter'],
            #                        p['site_symmetry'] ),
            # print
            # for o in p['positions']:
            #     print o,
            # print
    
            positions_site = []
            for pure_trans in lattice_symbols[ w['symbol'][0] ]:
                for o in p['positions']:
                    r,t = position2operation( o.strip('()').split(',') )
                    t += np.array( pure_trans )
                    t = t % 24
                    positions_site.append( (r,t) )
            positions_hall.append( positions_site )
        positions.append( positions_hall )
    return positions

def encode_wyckoff_positions( positions ):
    maxval = 0
    maxmulti = 0
    maxr = None
    maxt = None
    number = 0
    pos_encoded = []
    numerators = np.zeros( 24, dtype=int )

    for i, positions_hall in enumerate( positions ):
        pos_encoded_hall = []
        for positions_site in positions_hall:
            pos_encoded_site = []
            for r, t in positions_site:
                # for j in range(3):
                #     numerators[ t[j] ] += 1
                
                t_encode = t[0] * 576 + t[1] * 24 + t[2]
                r_encode = encode_rotation( r )
                total = t_encode * 45**3 + r_encode

                pos_encoded_site.append( total )

                if (r - decode_rotation( r_encode )).any():
                    print("Error in r %d" % ( i+1 ))
                    sys.exit(1)
                # else:
                #     if ( abs( decode_rotation( r_encode ) ) > 1 ).any():
                #         print i+1
                #         print decode_rotation( r_encode )
                #         print t
                #         print
                if (t - decode_trans( t_encode )).any():
                    print("Error in t")
                    sys.exit(1)
                    
                if total > maxval:
                    maxval = total
                    maxr = r.copy()
                    maxt = t.copy()
                    number = i+1

            pos_encoded_hall.append( pos_encoded_site )
        pos_encoded.append( pos_encoded_hall )

    # for i, n in enumerate( numerators ):
    #     print i, n
    
    # print "MAX: %d" % maxmulti
    # print "MAX: %d (%d)" % (maxval,number)
    # print maxr
    # print maxt
    
    return pos_encoded

def get_data_arrays( encoded_positions ):
    len_wyckoffs = [0]
    len_site = [0]
    sum_wyckoffs = 1
    sum_site = 1
    encpos_flat = [0]
    for h in encoded_positions:
        len_wyckoffs.append( sum_wyckoffs )
        sum_wyckoffs += len(h)
        for w in h:
            len_site.append( sum_site )
            sum_site += len(w)
            encpos_flat += w

    len_wyckoffs.append( sum_wyckoffs )
    len_site.append( sum_site )

    return len_wyckoffs, len_site, encpos_flat
        
if __name__ == '__main__':
    
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults( is_site = False,
                         is_position = False,
                         is_wyckoff = False,
                         is_site_short = False,
                         is_wyckoff_short = False )
                         
    parser.add_option("--site", dest="is_site",
                      action="store_true",
                      help="Index of site symmetries of each space group")
    parser.add_option("--position", dest="is_position",
                      action="store_true",
                      help="Index of positions of each site symmetries and space group")
    parser.add_option("--wyckoff", dest="is_wyckoff",
                      action="store_true",
                      help="Encoded Wyckoff position operators")
    parser.add_option("--wyckoff_short", dest="is_wyckoff_short",
                      action="store_true",
                      help="Encoded Wyckoff position operators")
    (options, args) = parser.parse_args()

    wyckoff = read_wyckoff_csv( args[0] )
    operations = get_wyckoff_positions( wyckoff )
    encoded_positions = encode_wyckoff_positions( operations )
    len_wyckoffs, len_site, encpos_flat = \
        get_data_arrays( encoded_positions )
    
    # Test
    # sum_elem = run_test( len_wyckoffs, len_site, o_flat )
    # print "%d = %d" % ( len( o_flat ), sum_elem )
    
    # Damp full arrays
    if options.is_position:
        damp_array_wyckoffs( len_wyckoffs )
    if options.is_site:
        damp_array_numsites( len_site )
    if options.is_wyckoff:
        damp_array_positions( encpos_flat )

    if options.is_wyckoff_short:
        damp_array_positions_short( encpos_flat, len_site )
