'''
Created on Aug 7, 2016
@author: Niels Lubbes
'''

from sage.all import *


class OrbRing:

    num_field = QQ

    vstr = ''
    vstr += 'x0,x1,x2,x3,x4,x5,x6,x7,x8,'
    vstr += 'v0,v1,v2,v3,v4,v5,v6,v7,v8,'
    vstr += 'c0,s0,c1,s1,'
    vstr += 't0,t1,t2,t3,t4,t5,t6,t7'

    R = PolynomialRing( num_field, var( vstr ), order = 'degrevlex' )


    @staticmethod
    def coerce( expr ):
        return sage_eval( str( expr ), OrbRing.R.gens_dict() )


    @staticmethod
    def random_int( val ):
        '''
        INPUT:
            - "val" -- An integer.
        OUTPUT:
            - A random element in the interval [-val,val]
        '''
        return int( ZZ.random_element( -val, val + 1 ) )

    @staticmethod
    def random_elt( lst ):
        '''
        INPUT:
            - "lst" -- A list.
        OUTPUT:
            - A random element in "lst".
        '''
        idx = int( ZZ.random_element( 0, len( lst ) ) )
        return lst[idx]

