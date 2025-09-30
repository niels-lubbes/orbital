'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''

from orbital.poly_maps import ring_dict
from orbital.poly_maps import invert_map
from orbital.poly_maps import invert_birational_map
from orbital.poly_maps import image_map
from orbital.poly_maps import preimage_map
from orbital.poly_maps import compose_maps
from orbital.poly_maps import euclidean_type_form
from orbital.poly_maps import hilbert_poly

from orbital.sage_interface import sage_PolynomialRing
from orbital.sage_interface import sage_NumberField
from orbital.sage_interface import sage_FractionField
from orbital.sage_interface import sage_QQ


class TestPolyMaps( object ):

    def test__ring_dict( self ):

        R = sage_PolynomialRing( sage_QQ, 'a' )
        a = R.gens()[0]
        B = sage_NumberField( [a ** 2 - sage_QQ( 1 ) / 2], 'a0' )
        a0 = B.gens()[0]
        B = sage_FractionField( sage_PolynomialRing( B, 'k' ) )

        chk = "{'k': k, 'a0': a0, '1': 1}"

        out = ring_dict( B )
        print( out )
        assert str( out ) == chk

    def test__invert_map__stereographic_projection( self ):
        X = '[x1^2 + x2^2 + x3^2 + x4^2 - x0^2]'
        f = '[x0-x4, x1, x2, x3]'
        chk = '[x0 + (-y0^2 - y1^2 - y2^2 - y3^2)/(2*y0), x1 + (-y1), x2 + (-y2), x3 + (-y3), x4 + (y0^2 - y1^2 - y2^2 - y3^2)/(2*y0)]'
        out = invert_map( f, X, sage_QQ )
        print( out )
        assert str( out ) == chk

    def test__invert_birational_map__stereographic_projection( self ):
        f = '[x0-x4, x1, x2, x3]'
        X = '[x1^2 + x2^2 + x3^2 + x4^2 - x0^2]'
        chk = '[8*y0^2 + 8*y1^2 + 8*y2^2 + 8*y3^2, 16*y0*y1, 16*y0*y2, 16*y0*y3, -8*y0^2 + 8*y1^2 + 8*y2^2 + 8*y3^2]'
        out = invert_birational_map( f, X, sage_QQ )

        print( out )
        assert str( out ) == chk

    def test__image_map__two_sphere_cyclide( self ):
        f = '[ x0^2+x1^2+x2^2+x3^2,  2*x0*x1,  2*x0*x2,  2*x0*x3,-x0^2+x1^2+x2^2+x3^2]'
        X = '[(x1^2+x2^2+x3^2)^2-x0^2*(8*x1^2-2*x2^2-2*x3^2-2*x0^2)]'
        chk = '[5*y1^2 - 5*y2^2 - 5*y3^2 + 2*y0*y4 - 6*y4^2, 5*y0^2 - 10*y2^2 - 10*y3^2 + 2*y0*y4 - 11*y4^2]'

        out = image_map( f, X, sage_QQ )

        print( out )
        assert str( out ) == chk

    def test__preimage_map__central_projection( self ):

        # we compute the preimage of elliptic hyperboloid
        # of one sheet, with respect to the central
        # projection of the projective 3-sphere.
        # Afterward we compute its image in order to
        # verify that we obtain the same ideal of the
        # variety.

        X = '[x1^2 + x2^2 + x3^2 + x4^2 - x0^2]'
        f = '[x0, x1, x2, x3]'
        Y = '[ y1^2+2*y2^2-y3^2-y0^2 ]'  # EH1
        chk = '[[x2^2 - 2*x3^2 - x4^2, x0^2 - x1^2 - 3*x3^2 - 2*x4^2]]'

        out = preimage_map( f, X, Y, sage_QQ )

        print( out )
        assert str( out ) == chk

        out2 = image_map( f, str( out[0] ), sage_QQ )
        chk2 = '[y0^2 - y1^2 - 2*y2^2 + y3^2]'

        print( out2 )
        assert str( out2 ) == chk2

    def test__preimage_map__stereographic_projection( self ):

        # we compute the preimage of elliptic hyperboloid
        # of one sheet, with respect to the stereographic
        # projection of the projective 3-sphere.
        # Afterward we compute its image in order to
        # verify that we obtain the same ideal of the
        # variety.

        f = '[x0-x4, x1, x2, x3]'
        X = '[x1^2 + x2^2 + x3^2 + x4^2 - x0^2]'
        Y = '[ y1^2+2*y2^2-y3^2-y0^2 ]'  # EH1
        chk = '[[x2^2 - 2*x3^2 + 2*x0*x4 - 2*x4^2, x0^2 - x1^2 - 3*x3^2 + 2*x0*x4 - 3*x4^2]]'

        out = preimage_map( f, X, Y, sage_QQ )

        print( out )
        assert str( out ) == chk

        out2 = image_map( f, str( out[0] ), sage_QQ )
        chk2 = '[y0^2 - y1^2 - 2*y2^2 + y3^2]'

        print( out2 )
        assert str( out2 ) == chk2

    def test__compose_maps__P3( self ):

        f = '[x0,x1+x0,x2,x3+x0]'
        g = '[x0^2+x1^2+x2^2+x3^2, 2*x0*x1, 2*x0*x2, 2*x0*x3, -x0^2+x1^2+x2^2+x3^2]'
        chk = '[3*x0^2 + 2*x0*x1 + x1^2 + x2^2 + 2*x0*x3 + x3^2, 2*x0^2 + 2*x0*x1, 2*x0*x2, 2*x0^2 + 2*x0*x3, x0^2 + 2*x0*x1 + x1^2 + x2^2 + 2*x0*x3 + x3^2]'

        out = compose_maps( g, f )

        print( out )
        assert str( out ) == chk

    def test__euclidean_type_form__31( self ):
        X = 'x0^3 - 4*x0^2*x1 + 8*x1^3 + 6*x0*x2^2 + 8*x1*x2^2 - 4*x0^2*x3 + 8*x0*x1*x3 - 12*x1^2*x3 - 12*x2^2*x3 + 10*x0*x3^2 + 8*x1*x3^2 - 12*x3^3'
        chk1 = '(4) * (2*x1 - 3*x3) * (x1^2 + x2^2 + x3^2) + x0 * (x0^2 - 4*x0*x1 + 6*x2^2 - 4*x0*x3 + 8*x1*x3 + 10*x3^2)'
        chk2 = '4*(x^2 + y^2 + z^2)*(2*x - 3*z) + 6*y^2 + 8*x*z + 10*z^2 - 4*x - 4*z + 1'
        chk = chk1 + '\n' + chk2

        out = euclidean_type_form( X, sage_QQ )
        print( out )
        assert out == chk

    def test__hilbert_poly__P4( self ):
        Y = '[y1^4 - 10*y1^2*y2^2 - 6*y1^2*y3^2 - 16*y1^2*y4^2 + 25*y2^4 + 30*y2^2*y3^2 + 56*y2^2*y4^2 + 9*y3^4 + 32*y3^2*y4^2 + 32*y4^4, 2*y0*y4 + y1^2 - 5*y2^2 - 3*y3^2 - 6*y4^2, y0*y1^2 - 5*y0*y2^2 - 3*y0*y3^2 + 5*y1^2*y4 - 13*y2^2*y4 - 7*y3^2*y4 - 16*y4^3, y0^2 - y1^2 - y2^2 - y3^2 - y4^2]'
        chk = '2*t^2 + 2*t + 1'

        out = hilbert_poly( Y, sage_QQ )
        print( out )
        assert str( out ) == chk


if __name__ == '__main__':

    TestPolyMaps().test__ring_dict()
    TestPolyMaps().test__invert_map__stereographic_projection()
    TestPolyMaps().test__invert_birational_map__stereographic_projection()
    TestPolyMaps().test__image_map__two_sphere_cyclide()
    TestPolyMaps().test__preimage_map__central_projection()
    TestPolyMaps().test__preimage_map__stereographic_projection()
    TestPolyMaps().test__compose_maps__P3()
    TestPolyMaps().test__euclidean_type_form__31()
    TestPolyMaps().test__hilbert_poly__P4()

    pass
