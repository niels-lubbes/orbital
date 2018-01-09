'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes

'''

from sage_interface import sage_QQ
from sage_interface import sage_NumberField
from sage_interface import sage_FractionField
from sage_interface import sage_PolynomialRing
from sage_interface import sage__eval
from sage_interface import sage_ideal
from sage_interface import sage_factor
from sage_interface import sage_var
from sage_interface import sage_set_verbose

from class_orb_tools import OrbTools

from poly_maps import image_map
from poly_maps import invert_birational_map
from poly_maps import compose_maps
from poly_maps import euclidean_type_form
from poly_maps import preimage_map
from poly_maps import ring_dict
from poly_maps import hilbert_poly

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

import sys


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


if __name__ == '__main__':

    #  Debug output settings
    #
    sage_set_verbose( -1 )
    OrbTools.filter( '__main__.py' )  # only print if verbose output by module <file_name>
    # OrbTools.filter( None )  # print all verbose output
    OrbTools.start_timer()

    if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = '~/'


    #########################################
    #                                       #
    # (Un)comment one or more use cases     #
    #                                       #
    #########################################

    usecase__two_sphere_cyclide()


    #########################################
    #                                       #
    # End of list of use case methods.      #
    #                                       #
    #########################################

    # end timing
    OrbTools.end_timer()

    print( '\nThe End' )



