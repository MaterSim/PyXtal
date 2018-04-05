#!/usr/bin/env python

# This will (hopefully) be the code to extract symmetry operations
# from Hall symbols

'''
Code from https://github.com/atztogo/spglib/tree/master/database
Originally created by Atsushi Togo.
Modified by Scott Fredericks for python 3 compatibility.
'''

import numpy as np

lattice_symbols = {
    'P': [ [ 0, 0, 0 ] ],
    'A': [ [ 0, 0, 0 ], [ 0, 1./2, 1./2 ] ],
    'B': [ [ 0, 0, 0 ], [ 1./2, 0, 1./2 ] ],
    'C': [ [ 0, 0, 0 ], [ 1./2, 1./2, 0 ] ],
    'I': [ [ 0, 0, 0 ], [ 1./2, 1./2, 1./2 ] ],
    'R': [ [ 0, 0, 0 ], [ 2./3, 1./3, 1./3 ], [ 1./3, 2./3, 2./3 ] ],
    'H': [ [ 0, 0, 0 ], [ 2./3, 1./3, 0 ], [ 1./3, 2./3, 0 ] ],
    'F': [ [ 0, 0, 0 ], [ 0, 1./2, 1./2 ], [ 1./2, 0, 1./2 ], [ 1./2, 1./2, 0 ] ]
    }

rotation_matrices = {
    '1x': [ [ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ] ],
    '1y': [ [ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ] ],
    '1z': [ [ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ] ],
    '2x': [ [ 1, 0, 0 ],
            [ 0,-1, 0 ],
            [ 0, 0,-1 ] ],
    '2y': [ [-1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0,-1 ] ],
    '2z': [ [-1, 0, 0 ],
            [ 0,-1, 0 ],
            [ 0, 0, 1 ] ],
    '3x': [ [ 1, 0, 0 ],
            [ 0, 0,-1 ],
            [ 0, 1,-1 ] ],
    '3y': [ [-1, 0, 1 ],
            [ 0, 1, 0 ],
            [-1, 0, 0 ] ],
    '3z': [ [ 0,-1, 0 ],
            [ 1,-1, 0 ],
            [ 0, 0, 1 ] ],
    '4x': [ [ 1, 0, 0 ],
            [ 0, 0,-1 ],
            [ 0, 1, 0 ] ],
    '4y': [ [ 0, 0, 1 ],
            [ 0, 1, 0 ],
            [-1, 0, 0 ] ],
    '4z': [ [ 0,-1, 0 ],
            [ 1, 0, 0 ],
            [ 0, 0, 1 ] ],
    '6x': [ [ 1, 0, 0 ],
            [ 0, 1,-1 ],
            [ 0, 1, 0 ] ],
    '6y': [ [ 0, 0, 1 ],
            [ 0, 1, 0 ],
            [-1, 0, 1 ] ],
    '6z': [ [ 1,-1, 0 ],
            [ 1, 0, 0 ],
            [ 0, 0, 1 ] ],
    '2px': [ [-1, 0, 0 ],     # b-c
             [ 0, 0,-1 ],
             [ 0,-1, 0 ] ],
    '2ppx': [ [-1, 0, 0 ],    # b+c
              [ 0, 0, 1 ],
              [ 0, 1, 0 ] ],
    '2py': [ [ 0, 0,-1 ],     # a-c
             [ 0,-1, 0 ],
             [-1, 0, 0 ] ],
    '2ppy': [ [ 0, 0, 1 ],    # a+c
              [ 0,-1, 0 ],
              [ 1, 0, 0 ] ],
    '2pz': [ [ 0,-1, 0 ],     # a-b
             [-1, 0, 0 ],
             [ 0, 0,-1 ] ],
    '2ppz': [ [ 0, 1, 0 ],    # a+b
              [ 1, 0, 0 ],
              [ 0, 0,-1 ] ],
    '3*': [ [ 0, 0, 1 ],     # a+b+c
            [ 1, 0, 0 ],
            [ 0, 1, 0 ] ]
    }

translations = {
    'a': [ 1./2, 0, 0 ],
    'b': [ 0, 1./2, 0 ],
    'c': [ 0, 0, 1./2 ],
    'n': [ 1./2, 1./2, 1./2 ],
    'u': [ 1./4, 0, 0 ],
    'v': [ 0, 1./4, 0 ],
    'w': [ 0, 0, 1./4 ],
    'd': [ 1./4, 1./4, 1./4 ]
    }


def read_spg_csv( filename="spg.csv" ):
    hall_symbols = []
    spg_num = []
    for line in open( filename ):
        data = line.split(',')
        hall_symbols.append( [data[6], data[4]] )
    return hall_symbols

class HallSymbol:
    def __init__( self, hall_symbol ):
        self.hall_symbol = hall_symbol.split()
        self.__decompose()
        self.__operations()
        self.__full_operations()

    def get_LNV( self ):
        return self.L, self.N, self.V

    def get_operations( self ):
        return self.R, self.T

    def get_full_operations( self ):
        return self.G_R, self.G_T

    def __full_operations( self ):
        E = np.array( rotation_matrices['1x'] )
        T0 = np.zeros( 3, dtype=float )
        if self.L[0] == '-':
            G_R = [ E, -E ]
            G_T = [ T0, T0 ]
        else:
            G_R = [ E ]
            G_T = [ T0 ]

        for r, t in zip( self.R, self.T ):
            G2_R, G2_T = self.__get_group( r, t )
            G_R, G_T = self.__multiply_groups( G_R, G_T, G2_R, G2_T )

        if self.V is not None:
            G_T = self.__change_of_basis( G_R, G_T )

        G_R_conventional = []
        G_T_conventional = []
        for t in lattice_symbols[ self.L[-1] ]:
            G_R_conventional, G_T_conventional = \
                self.__lattice_translation( G_R_conventional,
                                            G_T_conventional,
                                            G_R, G_T, np.array( t ) )

        self.G_R = G_R_conventional
        self.G_T = G_T_conventional

    def __change_of_basis( self, G_R, G_T ):
        G_T_new = []
        v = self.V.astype(float) / 12
        for r, t in zip( G_R, G_T ):
            G_T_new.append( -np.dot( r, v ) + t + v )
        return G_T_new

    def __lattice_translation( self, G_R, G_T, G_R0, G_T0, translation ):
        for r, t in zip( G_R0, G_T0 ):
            G_R.append( r.copy() )
            t_new = t + translation
            t_new -= np.round( t_new )
            for i, x in enumerate( t_new ):
                if x < -0.001:
                    t_new[i] += 1
            G_T.append( t_new )
        return G_R, G_T
            
    def __multiply_groups( self, G1_R, G1_T, G2_R, G2_T ): # G2xG1
        G_R = []
        G_T = []
        for r1, t1 in zip( G2_R, G2_T ):
            for r2, t2 in zip( G1_R, G1_T ):
                G_R.append( np.dot( r1, r2 ) )
                G_T.append( np.dot( r1, t2 ) + t1 )

        return G_R, G_T

    def __get_group( self, r, t ):
        G_R, G_T = self.__get_group_recursive( 
            [ np.array( r ) ], [ np.array( t ) ] )
        r = G_R.pop()
        t = G_T.pop()
        G_R.insert( 0, r )
        G_T.insert( 0, t )
        return G_R, G_T

    def __get_group_recursive( self, G_R, G_T ):
        if not ( G_R[-1] == rotation_matrices['1x'] ).all():
            r = np.dot( G_R[-1], G_R[0] )
            t = np.dot( G_R[-1], G_T[0] ) + G_T[-1]
            G_R.append( r )
            G_T.append( t )
            G_R, G_T = self.__get_group_recursive( G_R, G_T )
        return G_R, G_T
    
    def __operations( self ):
        R = []
        T = []
        for N in self.N:
            rot = np.array( rotation_matrices[ N[1]+N[2] ] )
            if N[0] == '-':
                rot = -rot
            R.append( rot )

            trans = np.zeros( 3, dtype=float )
            if N[3] is not None:
                for t in N[3]:
                    if t == '1' or t == '2' or t == '3' or t == '4' or t == '5':
                        trans_screw = float( t ) / int( N[1] )
                        if N[2] == 'x':
                            trans[0] += trans_screw
                        elif N[2] == 'y':
                            trans[1] += trans_screw
                        elif N[2] == 'z':
                            trans[2] += trans_screw
                        else:
                            raise
                    else:
                        trans += np.array( translations[ t ] )
            T.append( trans )
                    
        self.R = np.array( R, dtype=int )
        self.T = np.array( T, dtype=float )

    def __rotation_matrix( self, str ):
        pass

    # Decompose Hall symbol
    # The following methods are used by __decompose().
    def __decompose( self ):
        L = self.hall_symbol.pop(0)
        N = []
        V = None
        precededN = 0  
        for i, ms in enumerate( self.hall_symbol ):
            if ms[0] == '(':
                V = self.__change_of_basis_symbol( self.hall_symbol[i+2] )
                break
            else:
                N.append( self.__matrix_symbol( ms, i, precededN ) )
                precededN = int( N[-1][1][0] )

        self.L = L
        self.N = N
        self.V = V

    def __matrix_symbol( self, N, i, precededN ):
        if N[0] == '-':
            improper = '-'
            N = N[1:]
        else:
            improper = None

        N, R, A = self.__rotation( N, i, precededN )

        if len( N ) > 0:
            T = self.__translation( N )
        else:
            T = None

        return [ improper, R, A, T ]

    def __rotation( self, N, i, precededN ):
        A = None
        if N[0] == '2':
            if len( N ) > 1:  # 2"
                if N[1] == '=':
                    R = '2pp'
                    A = 'z'
                    N = N[2:]
            if i == 1 and A is None:
                if precededN == 2 or precededN == 4:   # 2x
                    R = '2'
                    A = 'x'
                    N = N[1:]
                elif precededN == 3 or precededN == 6: # 2'
                    R = '2p'
                    A = 'z'
                    N = N[1:]
        elif N[0] == '3':     # 3*
            if i == 2:
                R = '3'
                A = '*'
                N = N[1:]
            elif len( N ) > 1:    
                if N[1] == '*':
                    R = '3'
                    A = '*'
                    N = N[2:]

        if A is None:
            R = N[0]
            N = N[1:]
            if len( N ) > 0 and i == 0:
                N, A = self.__principal_axis( N )
            else:
                A = 'z'

        return N, R, A

    def __principal_axis( self, N ):
        if N[0] == 'x':
            return N[1:], 'x'
        if N[0] == 'y':
            return N[1:], 'y'
        return N, 'z'
            
    def __translation( self, N ):
        T = []
        for i in range( len( N ) ):
            T.append( N[i] )
        return T

    def __change_of_basis_symbol( self, V ):
        if V[0] == '-':
            return np.array([0,0,-1])
        else:
            return np.array([0,0,1])
            

def damp_operations( filename ):
    hall_symbols = read_spg_csv( filename )
    count = 0
    print("  0       ,  /* dummy */")
    for i in range( 530 ):
        hs = HallSymbol( hall_symbols[ i ][0] )
        G_R, G_T = hs.get_full_operations()
        for j, (r, t) in enumerate( zip( G_R, G_T ) ):
            count += 1
            r_encode = encode_rotation( r )
            x = (t*12+0.1).astype(int)
            t_encode = x[0] * 144 + x[1] * 12 + x[2]
            total = t_encode * 3 ** 9 + r_encode
            print("  %-8d," % ( total ),)
            print(" /* %4d (%3d) [" % ( count, i+1 ),)
            print("%2d," * 9 % tuple( decode_rotation( total % (3**9) ) ),)
            print("%2d,%2d,%2d] */" % tuple( decode_trans( total // (3**9) ) ))

def damp_operations_old( filename ):
    hall_symbols = read_spg_csv( filename )
    count = 0
    for i in range( 530 ):
        hs = HallSymbol( hall_symbols[ i ][0] )
        G_R, G_T = hs.get_full_operations()
        for j, (r, t) in enumerate( zip( G_R, G_T ) ):
            count += 1
            print("{%3d," % (i+1),)
            print("%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d,%2d," % tuple( r.flatten() ),)
            print("%2d,%2d,%2d" % tuple( (t*12+0.1).astype(int) ),)
            print("}, /* %4d */" % count)

# Ternary numerical system
def encode_rotation( r ):
    r_sum = 0
    for i, x in enumerate( r.flatten() ):
        r_sum += ( x + 1 ) * 3**( 8 - i )
    return r_sum

def decode_rotation( c ):
    r = []
    for i in range( 8, -1, -1 ):
        r.append( ( c % (3**(i+1)) ) // ( 3**i ) - 1 )
    return np.array( r )

def decode_trans( c ):
    return c // 144, ( c % 144 ) // 12, ( c % 12 )

def get_reference_to_operations():
    hall_symbols = read_spg_csv( args[0] )
    count = 0
    for i in range( 530 ):
        hs = HallSymbol( hall_symbols[ i ][0] )
        G_R, G_T = hs.get_full_operations()
        print("{%4d,%5d}," % ( len( G_R ), count+1 ))
        count += len( G_R )

def watch_hs( filename, number ):
    num = number - 1
    hall_symbols = read_spg_csv( filename )
    h = HallSymbol( hall_symbols[ num ][0] )
    print(h.get_LNV())
    G_R, G_T = h.get_full_operations()
    print(number, ":", hall_symbols[ num ][0], "(", len(G_R), ")")
    for i, (r, t) in enumerate( zip( G_R, G_T ) ):
        print("-----", i+1, "-----")
        print(r, t)

if __name__ == '__main__':
    # Parse options
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults( watch_hs = False,
                         damp_operations = False,
                         shift = None,
                         origin = None )
    parser.add_option("--hs", dest="watch_hs",
                      action="store_true",
                      help="spg.csv [spg NUM]")
    parser.add_option("--damp", dest="is_damp",
                      action="store_true" )
    parser.add_option("--reference", dest="is_reference",
                      action="store_true" )
    
    (options, args) = parser.parse_args()
    
    if options.is_damp:
        damp_operations( args[0] )
    if options.is_reference:
        get_reference_to_operations()
    elif options.watch_hs:
        watch_hs( args[0], int( args[1] ) )



