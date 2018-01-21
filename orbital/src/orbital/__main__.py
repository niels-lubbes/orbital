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

from orbital.pov.pov_blum_cyclide import blum_cyclide
from orbital.pov.pov_ring_cyclide import ring_cyclide
from orbital.pov.pov_dp6_smooth import dp6_smooth
from orbital.pov.pov_quadric_smooth import quadric_smooth
from orbital.pov.pov_perseus_cyclide import perseus_cyclide
from orbital.pov.pov_CH1_cyclide import CH1_cyclide
from orbital.pov.pov_spindle_cyclide import spindle_cyclide
from orbital.pov.pov_horn_cyclide import horn_cyclide
from orbital.pov.pov_veronese import veronese


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
    TODO
    '''
    # blum_cyclide()
    # ring_cyclide()
    # dp6_smooth()
    # quadric_smooth()
    # perseus_cyclide()
    # CH1_cyclide()
    # spindle_cyclide()
    # horn_cyclide()
    veronese()

if __name__ == '__main__':

    #  Debug output settings
    #
    sage_set_verbose( -1 )
    OrbTools.filter( '__main__.py' )  # only print if verbose output by module <file_name>
    OrbTools.filter( None )  # print all verbose output
    OrbTools.start_timer()

    if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = '/home/niels/Desktop/n/src/output/orb/povray/'

    #########################################
    #                                       #
    # (Un)comment one or more use cases     #
    #                                       #
    #########################################

    # usecase__two_sphere_cyclide()
    usecase_povray()

    #########################################
    #                                       #
    # End of list of use case methods.      #
    #                                       #
    #########################################

    # end timing
    OrbTools.end_timer()

    print( '\nThe End' )



