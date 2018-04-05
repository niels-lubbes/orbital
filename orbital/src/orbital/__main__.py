'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''

import os
import sys

from orbital.class_orb_tools import OrbTools

from orbital.poly_maps import image_map
from orbital.poly_maps import invert_birational_map
from orbital.poly_maps import compose_maps
from orbital.poly_maps import euclidean_type_form
from orbital.poly_maps import preimage_map
from orbital.poly_maps import ring_dict
from orbital.poly_maps import hilbert_poly

from orbital.sage_interface import sage_set_verbose
from orbital.sage_interface import sage_PolynomialRing
from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_NumberField
from orbital.sage_interface import sage_FractionField
from orbital.sage_interface import sage_matrix

from orbital.pov.pov_blum_cyclide import blum_cyclide
from orbital.pov.pov_ring_cyclide import ring_cyclide
from orbital.pov.pov_dp6_smooth import dp6_smooth
from orbital.pov.pov_dp6_sing import dp6_sing
from orbital.pov.pov_quadric_smooth import quadric_smooth
from orbital.pov.pov_perseus_cyclide import perseus_cyclide
from orbital.pov.pov_CH1_cyclide import CH1_cyclide
from orbital.pov.pov_spindle_cyclide import spindle_cyclide
from orbital.pov.pov_horn_cyclide import horn_cyclide
from orbital.pov.pov_veronese import veronese
from orbital.pov.pov_dp8_clifford import dp8_clifford

from linear_series.class_linear_series import LinearSeries
from linear_series.class_base_points import BasePointTree
from linear_series.class_poly_ring import PolyRing
from linear_series.get_linear_series import get_mon_lst

from orbital.surface_in_quadric import get_surf


def usecase__two_sphere_cyclide():
    '''
    Construct parametrization of families of circles
    on the sphere cyclide. The sphere cyclide is 
    a quartic del Pezzo surface embedded into
    projective 3-sphere. This surface has two families
    of circles, but is not rational. 
    '''

    R = sage_PolynomialRing( sage_QQ, 'a' )
    a = R.gens()[0]
    B = sage_NumberField( [a ** 2 - 1 / sage_QQ( 2 )], 'a0' )
    B = sage_FractionField( sage_PolynomialRing( B, 'k' ) )
    OrbTools.p( ring_dict( B ) )

    S = '[-y0^2+y1^2+y2^2+y3^2+y4^2]'
    X = '[(x1^2+x2^2+x3^2)^2-x0^2*(8*x1^2-2*x2^2-2*x3^2-2*x0^2)]'
    OrbTools.p( euclidean_type_form( X[1:-1] ) )

    p = '[y0-y4,y1,y2,y3]'
    m = '[y4,y1+a0*y0]'
    P = '[x0-k*x1]'

    f = invert_birational_map( p, S )
    OrbTools.p( 'f =', f )

    Y = image_map( f, X )
    OrbTools.p( 'Y =', len( Y ), Y )
    OrbTools.p( 'Y =', hilbert_poly( Y ) )

    F = preimage_map( m, Y, P, B )
    OrbTools.p( 'F =', len( F ), F )

    C = image_map( p, F, B )
    OrbTools.p( 'C =', len( C ), C )

    for i in range( len( C ) ):
        s = str( C[i] )
        for ( a, b ) in [( 'a0', 'sqrt(1/2)' ), ( 'k', 'p1' ), ( 'x0', '1' ), ( 'x1', 'x' ), ( 'x2', 'y' ), ( 'x3', 'z' )]:
            s = s.replace( a, b )
        print( s )


def usecase_povray():
    '''
    Create povray images. This may take a long time.
    '''
    blum_cyclide()
    ring_cyclide()
    dp6_smooth()
    dp6_sing()
    quadric_smooth()
    perseus_cyclide()
    CH1_cyclide()
    spindle_cyclide()
    horn_cyclide()
    veronese()
    dp8_clifford()


def usecase_celestial_types():
    '''
    Let surface X in P^m be the blowups of P^1xP^1 in either 0 or 2 complex conjugate  
    points so that m=8 and m=6 respectively. 
    Here P^1xP^1 denotes the fiber product of the projective line with itself. 
    We verify for n in [3,4,5,6,7] whether X can be linearly projected into the 
    projective n-sphere S^n.  
    '''

    #
    # P^1xP^1
    #
    # We construct a parametrization of the double Segre surface dP8 in
    # projective 8-space. This surface contains 2 conics through each point
#     ring = PolyRing( 'x,y,v,w', True )
#     ls_dP8 = LinearSeries( get_mon_lst( [2, 2], ring.gens() ), ring )
#     OrbTools.p( 'ls_dP8 =', ls_dP8 )
#     get_surf( ls_dP8, ( 7 + 1, 1 ), [-8, -8, -10, -6, -7, 6, 5, 6, 0, -8, -2, -7, -7, 7, -1, 0, -9, 7, 1, -9]  )
#     get_surf( ls_dP8, ( 6 + 1, 1 ) )
#     get_surf( ls_dP8, ( 5 + 1, 1 ) )
#     get_surf( ls_dP8, ( 4 + 1, 1 ) )
#     get_surf( ls_dP8, ( 3 + 1, 1 ) )

    #
    # P^1xP^1 blown up in two general complex conjugate points
    #
    # We construct a parametrization of a sextic del Pezzo surface dP6
    # in projective 6-space, that contains 3 conics through each point.
    # We show that dP6 can be projected into S^5 and S^4.
    #
    a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]
    bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bp_tree.add( 'xv', ( -a0, a0 ), 1 )
    bp_tree.add( 'xv', ( a0, -a0 ), 1 )
    ls_dP6 = LinearSeries.get( [2, 2], bp_tree )
    get_surf( ls_dP6, ( 5 + 1, 1 ), [-9, -6, 1, 4, -1, -8, -5, -5, -4, 8, 1] )
    prv_Q = sage_matrix( [( 0, 0, 0, 1, 0, 1, 1 ), ( 1, 1, 0, 1, 0, 0, 0 ), ( 0, 1, 1, 0, 1, 0, 1 ), ( 1, 1, 0, 0, 1, 0, 0 ), ( 1, 1, 1, 1, 1, 1, 0 ), ( 1, 0, 0, 1, 0, 1, 1 )] )
    get_surf( ls_dP6, ( 4 + 1, 1 ), [-1, -9, -10, -7, -10, -8, 0], prv_Q )

    #
    # P^1xP^1 blown up in two complex conjugate points that lie in the same fiber
    #
    # We construct a parametrization of a sextic weak del Pezzo surface wdP6
    # in projective 6-space, that contains 2 conics through each point.
    # We show that wdP6 can be projected into S^5 and S^4.
    #
    a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]
    bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bp_tree.add( 'xv', ( a0, 0 ), 1 )  # the complex conjugate base points lie in the same fiber
    bp_tree.add( 'xv', ( -a0, 0 ), 1 )
    ls_wdP6 = LinearSeries.get( [2, 2], bp_tree )
    get_surf( ls_wdP6, ( 5 + 1, 1 ), [-6, 8, -7, -8, 0, -8, 2, -5, -8] )
    prv_Q = sage_matrix( [( 1, 0, 0, 1, 1, 1, 0 ), ( 1, 0, 0, 1, 0, 1, 1 ), ( 0, 1, 1, 1, 0, 1, 0 ), ( 0, 0, 0, 0, 1, 0, 0 ), ( 0, 1, 0, 1, 1, 1, 0 ), ( 0, 0, 0, 1, 1, 1, 0 )] )
    get_surf( ls_wdP6, ( 4 + 1, 1 ), [-2, -7, -6, -10, -2, -4, 4] , prv_Q )


if __name__ == '__main__':

    # Debug output settings
    #
    sage_set_verbose( -1 )
    mod_lst = []
    mod_lst += ['__main__.py']
    OrbTools.filter( mod_lst )  # output only from specified modules
    OrbTools.filter( None )  # print all verbose output, comment to disable.
    OrbTools.start_timer()

    # set environment variables for output, Maple and Magma
    # Maple and Magma are optionally used in orbital.prod.orb_product.
    if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = '/home/niels/Desktop/n/src/output/orb/povray/'
    os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'
    os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'

    #########################################
    #                                       #
    # (Un)comment one or more use cases     #
    #                                       #
    #########################################

    # usecase__two_sphere_cyclide()
    # usecase_povray()  # takes a long time
    # usecase_celestial_types()

    dp6_sing()
    veronese()

    #########################################
    #                                       #
    # End of list of use case methods.      #
    #                                       #
    #########################################

    # end timing
    OrbTools.end_timer()

    print( '\nThe End' )



