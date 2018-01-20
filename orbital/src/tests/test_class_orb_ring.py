'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_n

from orbital.class_orb_ring import OrbRing

from linear_series.class_poly_ring import PolyRing



class TestOrbRing( object ):

    def test__coerce( self ):
        assert OrbRing.coerce( 'x8+v8+s1+c1+t7' ) in OrbRing.R

    def test__random_int( self ):
        val = 10
        rnd = OrbRing.random_int( val )
        print( rnd )
        assert rnd in range( -val, val + 1 )

    def test__random_elt( self ):

        lst = range( 10 )
        elt = OrbRing.random_elt( lst )
        print( elt )
        assert elt in lst

    def test__approx_QQ_coef__1( self ):

        ring = PolyRing( 'x,y,z', True ).ext_num_field( 't^2 - 3' )
        a0 = ring.root_gens()[0]
        ci_idx = 1
        for lbl, elt in [( 'a0', a0 ), ( '-a0', -a0 )]:
            print( lbl + ' = ' + str( elt.abs() ) )
            for ci in elt.complex_embeddings():
                print( '\t\t' + str( ci ) )
            s = str( elt.complex_embeddings()[ci_idx] )
            print( lbl + '.complex_embedding()[' + str( ci_idx ) + '] = ' + s )

        print( '--- a0 ---' )
        out = OrbRing.approx_QQ_coef( a0, ci_idx )
        chk = sage_QQ( 3900231685776981 ) / 2251799813685248
        print( out )
        print( sage_n( out ) )
        assert out == chk

        print( '--- (-a0) ---' )
        out = OrbRing.approx_QQ_coef( -a0, ci_idx )
        print( out )
        print( sage_n( out ) )
        assert out == -chk


    def test__approx_QQ_coef__2( self ):

        ring = PolyRing( 'x,y,v,w', True )
        ring.ext_num_field( 't^2 - 3' )
        ring.ext_num_field( 't^2 + 1' )
        ring.ext_num_field( 't^2 + 2/7*a0*a1 + 3/7' )
        ring.ext_num_field( 't^2 - 2/7*a0*a1 + 3/7' )
        a0, a1, a2, a3 = ring.root_gens()

        q2 = sage_QQ( 1 ) / 2
        a = 2 * a0 / 3
        b = ( -a0 * a1 / 3 - q2 ) * a3
        c = ( a0 * a1 / 3 - q2 ) * a2
        d = ( a1 / 2 - a0 / 3 ) * a3
        e = ( -a1 / 2 - a0 / 3 ) * a2

        ci_idx = 5  # index for complex embedding

        C1 = -( b + c - d - e )
        C2 = -( b + c + d + e )

        for lbl, elt in [( 'a0', a0 ), ( 'C1', C1 ), ( 'C2', C2 )]:
            print( lbl + ' = ' + str( elt.abs() ) )
            for ci in elt.complex_embeddings():
                print( '\t\t' + str( ci ) )
            s = str( elt.complex_embeddings()[ci_idx] )
            print( lbl + '.complex_embedding()[' + str( ci_idx ) + '] = ' + s )

        # C1
        print( '--- C1 ---' )
        out = OrbRing.approx_QQ_coef( C1, ci_idx )
        print( out )
        print( sage_n( out ) )
        chk = -sage_QQ( 3687885631267691 ) / 2251799813685248
        assert out == chk

        # C2
        print( '--- C2 ---' )
        out = OrbRing.approx_QQ_coef( C2, ci_idx )
        print( out )
        print( sage_n( out ) )
        chk = sage_QQ( 2749869658604311 ) / 4503599627370496
        assert out == chk



    def test__approx_QQ_pol_lst( self ):
        ring = PolyRing( 'x,y,z', True ).ext_num_field( 't^2 - 3' )
        x, y, z = ring.gens()
        a0 = ring.root_gens()[0]
        q2 = sage_QQ( 1 ) / 2
        pol_lst = [x * y * a0 + z ** 2 + q2 * y ** 2 - q2 - a0 * x,
                   a0,
                   - a0 * x  ]
        out_lst = OrbRing.approx_QQ_pol_lst( pol_lst, 1 )

        c = sage_QQ( 3900231685776981 ) / 2251799813685248

        print( out_lst )
        assert out_lst[0] == x * y * c + z ** 2 + q2 * y ** 2 - q2 - c * x
        assert out_lst[1] == c
        assert out_lst[2] == -c * x



if __name__ == '__main__':

    # TestOrbRing().test__coerce()
    # TestOrbRing().test__random_int()
    # TestOrbRing().test__random_elt()
    # TestOrbRing().test__approx_QQ_coef__1()
    # TestOrbRing().test__approx_QQ_coef__2()
    TestOrbRing().test__approx_QQ_pol_lst()

    pass
