'''
Created on Aug 7, 2016
@author: Niels Lubbes
'''
from sage.all import *

from class_orb_ring import *
from orb_matrices import *
from orb_verbose import *


class OrbInput:
    '''
    OrbInput represents a 
        * 1-parameter subgroup in Aut(S^7) 
        * a circle in S^7
    where S^7 is the projective 7-sphere.
    
    This object is input to: 
        "orb_product.orb_product(input)"
        
    Additional attributes indicate which 
    aspects should be computed from the
    orbital product.
    '''



    def __init__( self ):
        #
        # Matrix representing a 1-parameter subgroup in Aut(S^7).
        #
        self.omat = None

        #
        # Matrix over "OrbRing.R" representing a circle C in S^7.
        # The circle C is obtained by applying "self.vmat" to the circle:
        # ZeroSet(<-x0^2+x1^2+x2^2, x3, x4, x5, x6, x7, x8>).
        #
        self.vmat = None

        #
        # Matrix representing a projection from P^8 to P^3.
        #
        self.pmat = None


        #
        # A dictionary containing info how the matrices
        # were constructed.
        #
        # The value corresponding to the key 'pmat' is
        # a 3-tuple whose elements are arguments for
        # the method "orb_matrices.get_mat()", which
        # return "self.pmat".
        #
        # Similarly with keys 'omat' and 'vmat'.
        #
        self.info_dct = {}


        #
        # A dictionary with booleans
        # which indicate what attributes of the
        # surface S should be computed, where surface S
        # is  the orbits of points on a circle w.r.t. a
        # 1-parameter subgroup.
        #
        self.do = {}
        self.do['pmz'] = True  # compute parametrization of S
        self.do['bpt'] = True  # compute base points of parametrization of S
        self.do['imp'] = True  # compute implicit equation of S
        self.do['dim'] = True  # compute dimension of S
        self.do['deg'] = True  # compute degree of S
        self.do['prj'] = True  # compute projection of S
        self.do['fct'] = True  # compute components of projection of S
        self.do['gen'] = True  # compute geometric genus of S
        self.do['sng'] = True  # compute singular locus of projection of S
        self.do['tst'] = True  # test parametrization/implicitization


    def set_short_str( self, short_str ):
        '''
        INPUT:
            - "short_str" -- A string of a list of length two
                             The first element is a string and 
                             the second element a dictionary 
                             with specifications as "self.info_dct".
                             See also "OrbOutput.get_short_str()".
        OUTPUT:
            - Return "self" where 
                  "self.info_dct", 
                  "self.omat", 
                  "self.vmat", 
                  "self.pmat" 
              are set according to input "short_str". 
        '''
        dct = sage_eval( short_str )[1]
        self.set( dct['pmat'], dct['omat'], dct['vmat'] )

        return self


    def set( self, p_tup, o_tup, v_tup ):
        '''
        INPUT:
            - "p_tup" -- A 3-tuple of strings with format as in
                         the docs of "orb_matrices.get_mat(...)".
            
            - "o_tup" -- Same specs as "p_tup".
            
            - "v_tup" -- Same specs as "p_tup".
        OUTPUT:
            - Sets "self.pmat", "self.omat" and "self.vmat"
              according "p_tup", "o_tup" and "v_tup" respectively. 
              The matrices are obtained with 
                  
                  "orb_matrices.get_mat()". 
              
              Set "self.info_dct" with info about the decomposition
              of the matrices:
              
                      "self.info_dct['pmat'] = p_tup"
                      "self.info_dct['omat'] = o_tup"
                      "self.info_dct['vmat'] = v_tup"
            
            - Return "self".            
        '''

        self.pmat = get_mat( p_tup[0], p_tup[1], p_tup[2] )
        self.omat = get_mat( o_tup[0], o_tup[1], o_tup[2] )
        self.vmat = get_mat( v_tup[0], v_tup[1], v_tup[2] )

        self.info_dct['pmat'] = p_tup
        self.info_dct['omat'] = o_tup
        self.info_dct['vmat'] = v_tup

        return self


    def random( self, coef_bnd = 3, random_pmat = True ):
        '''
        INPUT:
            - "coef_bnd"    -- Positive integer.
            - "random_pmat" -- A Boolean. 
             
        OUTPUT:
            - Sets (restricted) random values for "self.omat", "self.vmat" and "self.pmat".
              The rational coefficients of the matrices are in the interval:
              [-coef_bnd, +coef_bnd].
              
              If "random_pmat== True" then "self.pmat" is set random (P1),
              and to a standard projection otherwise (P0).  
              
              
            - Return "self".
        '''

        ch_lst = ['r', 's', 'p', 'm', 'a']

        #
        # All the coefficients are bound by coef_size-n
        # for some integer n.
        #
        coef_size = OrbRing.random_elt( range( 1, coef_bnd + 1 ) )

        #
        # random self.pmat
        #
        Pn = {True:'P1', False:'P0'}[random_pmat]
        self.pmat = get_mat( Pn, 'I', 'I' )
        self.info_dct['pmat'] = ( Pn, 'I', 'I' )


        #
        # random self.omat
        #
        rnd = OrbRing.random_elt( [0, 1] )
        if rnd == 0:
            ch_str = ''.join( [OrbRing.random_elt( ch_lst ) for i in range( 4 ) ] )
            B_str = 'O' + ch_str

            rnd2 = OrbRing.random_elt( [0, 1] )
            if rnd2 == 0:
                A_str = 'I'
                C_str = 'I'
            else:
                p = Permutations( range( 1, 8 + 1 ) )
                rp = p.random_element()
                rpI = Permutation( rp ).inverse()
                A_str = 'E' + str( list( rpI ) )
                C_str = 'E' + str( list( rp ) )

        elif rnd == 1:
            A_str = 'I'
            B_str = 'tT'
            C_str = 'I'

        self.omat = get_mat( A_str, B_str, C_str )
        self.info_dct['omat'] = ( A_str, B_str, C_str )

        #
        # random self.vmat
        #

        #     A -- self.vmat
        t_lst = [OrbRing.random_int( coef_size ) for i in range( 7 )]
        A_str = 'T' + str( t_lst )

        #     B -- self.vmat"
        #          note that small angles means small coefficients
        angle_lst = [ OrbRing.random_elt( range( 0, 360 ) ) for i in range( 4 ) ]
        ch_str = ''.join( [OrbRing.random_elt( ch_lst ) for i in range( 4 ) ] )
        B_str = 'R' + ch_str + str( angle_lst )

        #     C -- self.vmat
        rnd = OrbRing.random_elt( [0, 1, 2] )
        if rnd == 0:
            C_str = 'T' + str( [-t_lst[i] for i in range( 7 )] )  # inverse
        elif rnd == 1:
            C_str = 'T' + str( [OrbRing.random_int( coef_size ) for i in range( 7 )] )
        elif rnd == 2:
            C_str = 'I'

        #     A*B*C -- self.vmat
        self.vmat = get_mat( A_str, B_str, C_str )
        self.info_dct['vmat'] = ( A_str, B_str, C_str )

        return self


    # human readable string representation of object
    def __str__( self ):

        s = ''
        s += 15 * '.' + '\n'
        for key in self.info_dct.keys():
            s += key + '          = ' + str( self.info_dct[key][0] )
            s += ' ~~~ ' + str( self.info_dct[key][1] )
            s += ' ~~~ ' + str( self.info_dct[key][2] )
            s += '\n'

        s += 'do            = ' + str( self.do ) + '\n'
        s += 15 * '.'

        return s






