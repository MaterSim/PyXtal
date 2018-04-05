#!/usr/bin/env python

'''
Code from https://github.com/atztogo/spglib/tree/master/database
Originally created by Atsushi Togo.
Modified by Scott Fredericks for python 3 compatibility.
'''

import numpy as np
import database.hall2operations as h2o
import database.make_Wyckoff_db as wdb

characters = "abcdefghijklmnopqrstuvwxyz"

def get_Hall_operations( hall_file, hall_number ):
    hall_symbols = h2o.read_spg_csv( hall_file  )
    return h2o.HallSymbol( hall_symbols[ hall_number - 1 ][0] ).get_full_operations()

def get_wyckoff_position_operators( wyckoff_file, hall_number ):
    wyckoff = wdb.read_wyckoff_csv( wyckoff_file )
    operations = wdb.get_wyckoff_positions( wyckoff )
    return operations[ hall_number - 1 ]

def multiply_operations( spg_op, pos_op ):
    r, t = spg_op
    o_r, o_t = pos_op
    return ( np.dot( r, o_r ), np.dot( o_t, r.T ) + t )

def show_site_symmetries( rotations, translations, operators ):
    for i, o_site in enumerate( operators ):
        for j, o in enumerate( o_site ):
            o_r = o[0]
            o_t = o[1] / 24.0
            print("--------- %s, %d ----------" % (characters[len(operators) - i - 1], j+1))
            for r, t in zip( rotations, translations ):
                mo = multiply_operations( ( r, t ), ( o_r, o_t ) )
                diff = mo[1] - o_t
                if ( mo[0] - o_r == 0 ).all() and ( abs( diff - diff.round() ) < 0.01 ).all():
                    print(r, t)


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults( is_site = False )
    parser.add_option("--site", dest="is_site",
                      action="store_true",
                      help="Show site symmetries")
    (options, args) = parser.parse_args()
    
    rotations, translations = get_Hall_operations( args[0], int( args[2] ) )
    operators = get_wyckoff_position_operators( args[1], int( args[2] ) )
    

    if options.is_site:
        show_site_symmetries( rotations, translations, operators )
