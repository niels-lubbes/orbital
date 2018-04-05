'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''

from orbital.surface_in_quadric import get_S1xS1_pmz
from orbital.surface_in_quadric import get_deg_surf
from orbital.surface_in_quadric import get_prj_mat
from orbital.surface_in_quadric import approx_QQ
from orbital.surface_in_quadric import rand_surf_prj
from orbital.surface_in_quadric import get_surf
from orbital.surface_in_quadric import verify_get_surf
from orbital.surface_in_quadric import get_proj

from orbital.sage_interface import sage_PolynomialRing
from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage__eval
from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_diagonal_matrix

from linear_series.class_linear_series import LinearSeries
from linear_series.class_base_points import BasePointTree
from linear_series.class_poly_ring import PolyRing


class TestSurfaceInQuadric( object ):

    def test__get_S1xS1_pmz( self ):
        R = sage_PolynomialRing( sage_QQ, 'x,y,v,w' )
        x, y, v, w = R.gens()
        # xyvw_pmz_lst = [ x, y, x * w + y * v, x * v - y * w ]
        xyvw_pmz_lst = [ x, y, v, w ]

        out_lst = get_S1xS1_pmz( xyvw_pmz_lst )
        print( out_lst )
        assert str( out_lst ) == '[-s0 + 1, c0, -s1 + 1, c1]'


    def test__get_deg_surf( self ):
        R = sage_PolynomialRing( sage_QQ, 'x0,x1,x2,x3' )
        x0, x1, x2, x3 = R.gens()
        imp_lst = [ x0 * x1 - x2 ** 2 ]
        d = get_deg_surf( imp_lst, 3 )
        print( d )
        assert d == 2


    def test__get_prj_mat( self ):

        M = list( get_prj_mat( 3, 3 ) )
        print( str( M ) )
        assert str( M ) == '[(1, 0, 0), (0, 1, 0), (0, 0, 1)]'

        M = list( get_prj_mat( 2, 5, 0 ) )
        print( str( M ) )
        assert str( M ) == '[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0)]'

        M = list( get_prj_mat( 2, 5, 1 ) )
        print( str( M ) )
        assert str( M ) == '[(0, 1, 0, 0, 0), (0, 0, 1, 0, 0)]'

        M = list( get_prj_mat( 2, 5, 2 ) )
        print( str( M ) )
        assert str( M ) == '[(0, 0, 1, 0, 0), (0, 0, 0, 1, 0)]'


    def test__approx_QQ( self ):

        M = '''
        [
        (   0,   0, 1/2,    0,    0, -1/2,  1/2 ), 
        (   0,  -1,   0,    0,    0,    0,  1/2 ), 
        ( 1/2,   0,   0,  1/2,    0,    0,    1 ), 
        (   0,   0, 1/2,   -1,    0, -1/2,    0 ), 
        (   0,   0,   0,    0,   -1,    0, -1/2 ), 
        (-1/2,   0,   0, -1/2,    0,   -1, -1/2 ), 
        ( 1/2, 1/2,   1,    0, -1/2, -1/2,    0 )
        ]'''

        M = sage__eval( 'matrix( QQ,' + M + ')' )

        D, V = M.eigenmatrix_right()  # M*V==V*D
        assert M * V == V * D

        # M == V*D*~V == W.T*D*W == U.T*J*U
        # U == L*W and D == L.T*J*L
        W = sage_matrix( [col / col.norm() for col in V.columns()] )  # V.T normed
        L = sage_diagonal_matrix( [ d.abs().sqrt() for d in D.diagonal()] )
        J = sage_diagonal_matrix( [ d / abs( d ) for d in D.diagonal()] )
        U = L * W

        print( U.T * J * U )
        out = approx_QQ( U.T * J * U )
        print( out )

        assert M == out


    def test__rand_surf_prj( self ):

        a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]  # i
        bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
        bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )  # e1
        bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )  # e2

        ls = LinearSeries.get( [2, 2], bp_tree )  # |2(l1+l2)-e1-e2|
        prj_dim = 5
        prv_Q = sage_matrix( sage_QQ, [( 1, 1, 1, 1, 1, 0, 1 ), ( 0, 1, 0, 1, 0, 0, 0 ), ( 1, 0, 0, 1, 0, 1, 1 ), ( 1, 0, 1, 1, 0, 1, 0 ), ( 0, 0, 1, 0, 1, 0, 1 ), ( 0, 1, 1, 0, 1, 1, 1 )] )

        Q, pmz_lst, imp_lst = rand_surf_prj( ls, prj_dim, prv_Q )

        print( 'Q =\n' + str( list( Q ) ) )
        print( '\npmz_lst =\n' + str( pmz_lst ) )
        print( '\nimp_lst =\n' + str( imp_lst ) )

        assert Q == prv_Q
        assert str( pmz_lst ) == '[s0^2*c1^2 + c0*s0*c1*s1 - s0^2*c1*s1 + c0^2*s1^2 - c0*s0*s1^2 + s0^2*s1^2 - c0*s0*c1 + s0^2*c1 - 2*s0*c1^2 - 2*c0^2*s1 + 2*c0*s0*s1 - 2*s0^2*s1 - c0*c1*s1 + 2*s0*c1*s1 + c0*s1^2 - 2*s0*s1^2 + c0^2 - c0*s0 + s0^2 + c0*c1 - 2*s0*c1 + c1^2 - 2*c0*s1 + 4*s0*s1 - c1*s1 + s1^2 + c0 - 2*s0 + c1 - 2*s1 + 1, -s0^2*c1*s1 - c0*s0*s1^2 + s0^2*c1 + 2*c0*s0*s1 + 2*s0*c1*s1 + c0*s1^2 - c0*s0 - 2*s0*c1 - 2*c0*s1 - c1*s1 + c0 + c1, -c0*s0*c1^2 + c0^2*s1^2 - c0*s0*s1^2 + s0^2*s1^2 + c0*c1^2 - 2*c0^2*s1 + 2*c0*s0*s1 - 2*s0^2*s1 + c0*s1^2 - 2*s0*s1^2 + c0^2 - c0*s0 + s0^2 - 2*c0*s1 + 4*s0*s1 + s1^2 + c0 - 2*s0 - 2*s1 + 1, -c0*s0*c1^2 + s0^2*c1^2 - c0*s0*s1^2 + s0^2*s1^2 + c0*c1^2 - 2*s0*c1^2 + 2*c0*s0*s1 - 2*s0^2*s1 + c0*s1^2 - 2*s0*s1^2 - c0*s0 + s0^2 + c1^2 - 2*c0*s1 + 4*s0*s1 + s1^2 + c0 - 2*s0 - 2*s1 + 1, c0^2*c1^2 + s0^2*c1^2 + c0*s0*c1*s1 + c0^2*s1^2 - c0*s0*c1 - 2*s0*c1^2 - 2*c0^2*s1 - c0*c1*s1 + c0^2 + c0*c1 + c1^2, c0^2*c1^2 - c0*s0*c1^2 + s0^2*c1^2 - 2*c0^2*c1*s1 + c0*s0*c1*s1 - s0^2*c1*s1 + c0^2*s1^2 + 2*c0^2*c1 - c0*s0*c1 + s0^2*c1 + c0*c1^2 - 2*s0*c1^2 - 2*c0^2*s1 - c0*c1*s1 + 2*s0*c1*s1 + c0^2 + c0*c1 - 2*s0*c1 + c1^2 - c1*s1 + c1]'
        assert str( imp_lst ) == '[5*x0*x1 - 5*x1^2 + 12*x0*x2 - 13*x1*x2 + 4*x2^2 - 9*x0*x3 + 8*x1*x3 - 12*x2*x3 + 5*x3^2 - 8*x0*x4 - 5*x1*x4 - 8*x2*x4 + 13*x3*x4 + 4*x4^2 + 2*x0*x5 + 8*x1*x5 - 2*x2*x5 - 3*x3*x5 - 8*x4*x5 + 4*x5^2, 5*x0^2 - 5*x1^2 + 32*x0*x2 - 28*x1*x2 + 19*x2^2 - 34*x0*x3 + 28*x1*x3 - 47*x2*x3 + 25*x3^2 - 33*x0*x4 - 5*x1*x4 - 23*x2*x4 + 38*x3*x4 + 14*x4^2 + 7*x0*x5 + 28*x1*x5 - 7*x2*x5 - 8*x3*x5 - 28*x4*x5 + 14*x5^2, 8*x1^2*x4 - 46*x0*x2*x4 + 44*x1*x2*x4 - 38*x2^2*x4 + 42*x0*x3*x4 - 44*x1*x3*x4 + 83*x2*x3*x4 - 41*x3^2*x4 + 49*x0*x4^2 + 11*x1*x4^2 + 35*x2*x4^2 - 53*x3*x4^2 - 21*x4^3 - 8*x1^2*x5 + 20*x0*x2*x5 - 18*x1*x2*x5 + 12*x2^2*x5 - 16*x0*x3*x5 + 18*x1*x3*x5 - 31*x2*x3*x5 + 15*x3^2*x5 - 24*x0*x4*x5 - 58*x1*x4*x5 - 8*x2*x4*x5 + 31*x3*x4*x5 + 50*x4^2*x5 + x0*x5^2 + 21*x1*x5^2 - x2*x5^2 - 4*x3*x5^2 - 37*x4*x5^2 + 8*x5^3, 5*x0*x2*x3 - 5*x1*x2*x3 + 5*x2^2*x3 - 5*x0*x3^2 + 5*x1*x3^2 - 10*x2*x3^2 + 5*x3^3 + 6*x0*x2*x4 - 4*x1*x2*x4 + 2*x2^2*x4 - 12*x0*x3*x4 + 4*x1*x3*x4 - 11*x2*x3*x4 + 10*x3^2*x4 - 4*x0*x4^2 - 4*x2*x4^2 + 9*x3*x4^2 + 2*x4^3 - 6*x0*x2*x5 + 4*x1*x2*x5 - 2*x2^2*x5 + 7*x0*x3*x5 + x1*x3*x5 + 6*x2*x3*x5 - 5*x3^2*x5 + 5*x0*x4*x5 + 4*x1*x4*x5 + 3*x2*x4*x5 - 13*x3*x4*x5 - 6*x4^2*x5 - x0*x5^2 - 4*x1*x5^2 + x2*x5^2 + 4*x3*x5^2 + 6*x4*x5^2 - 2*x5^3, 40*x1^2*x3 - 10*x1*x2*x3 + 40*x2^2*x3 - 20*x0*x3^2 + 10*x1*x3^2 - 45*x2*x3^2 + 25*x3^3 + 276*x0*x2*x4 - 184*x1*x2*x4 + 92*x2^2*x4 - 307*x0*x3*x4 + 239*x1*x3*x4 - 331*x2*x3*x4 + 195*x3^2*x4 - 184*x0*x4^2 - 184*x2*x4^2 + 309*x3*x4^2 + 92*x4^3 - 120*x0*x2*x5 + 80*x1*x2*x5 - 40*x2^2*x5 + 135*x0*x3*x5 - 85*x1*x3*x5 + 125*x2*x3*x5 - 80*x3^2*x5 + 126*x0*x4*x5 + 184*x1*x4*x5 + 34*x2*x4*x5 - 219*x3*x4*x5 - 224*x4^2*x5 - 20*x0*x5^2 - 80*x1*x5^2 + 20*x2*x5^2 + 40*x3*x5^2 + 172*x4*x5^2 - 40*x5^3, 400*x2^3 - 80*x1*x2*x3 - 640*x2^2*x3 - 60*x0*x3^2 + 80*x1*x3^2 + 335*x2*x3^2 - 35*x3^3 - 792*x0*x2*x4 + 228*x1*x2*x4 - 264*x2^2*x4 + 419*x0*x3*x4 - 403*x1*x3*x4 + 637*x2*x3*x4 - 255*x3^2*x4 + 328*x0*x4^2 + 428*x2*x4^2 - 493*x3*x4^2 - 164*x4^3 + 352*x0*x2*x5 + 32*x1*x2*x5 - 16*x2^2*x5 - 159*x0*x3*x5 + 223*x1*x3*x5 - 187*x2*x3*x5 + 70*x3^2*x5 - 150*x0*x4*x5 - 428*x1*x4*x5 - 386*x2*x4*x5 + 441*x3*x4*x5 + 412*x4^2*x5 - 8*x0*x5^2 + 168*x1*x5^2 + 108*x2*x5^2 - 98*x3*x5^2 - 332*x4*x5^2 + 84*x5^3, 100*x1*x2^2 - 140*x1*x2*x3 - 70*x2^2*x3 + 20*x0*x3^2 + 40*x1*x3^2 + 105*x2*x3^2 - 55*x3^3 - 236*x0*x2*x4 + 224*x1*x2*x4 - 312*x2^2*x4 + 277*x0*x3*x4 - 249*x1*x3*x4 + 596*x2*x3*x4 - 290*x3^2*x4 + 349*x0*x4^2 + 75*x1*x4^2 + 199*x2*x4^2 - 294*x3*x4^2 - 137*x4^3 + 166*x0*x2*x5 - 144*x1*x2*x5 + 122*x2^2*x5 - 147*x0*x3*x5 + 109*x1*x3*x5 - 296*x2*x3*x5 + 135*x3^2*x5 - 200*x0*x4*x5 - 474*x1*x4*x5 - 88*x2*x4*x5 + 228*x3*x4*x5 + 346*x4^2*x5 + 11*x0*x5^2 + 219*x1*x5^2 - 11*x2*x5^2 - 34*x3*x5^2 - 281*x4*x5^2 + 72*x5^3, 400*x0*x2^2 - 480*x1*x2*x3 - 40*x2^2*x3 - 260*x0*x3^2 + 480*x1*x3^2 - 315*x2*x3^2 + 215*x3^3 + 288*x0*x2*x4 - 92*x1*x2*x4 - 704*x2^2*x4 - 791*x0*x3*x4 + 167*x1*x3*x4 + 407*x2*x3*x4 + 395*x3^2*x4 + 208*x0*x4^2 + 200*x1*x4^2 - 92*x2*x4^2 + 577*x3*x4^2 - 4*x4^3 - 328*x0*x2*x5 + 352*x1*x2*x5 + 24*x2^2*x5 + 451*x0*x3*x5 + 53*x1*x3*x5 - 57*x2*x3*x5 - 230*x3^2*x5 + 150*x0*x4*x5 - 508*x1*x4*x5 - 46*x2*x4*x5 - 949*x3*x4*x5 - 68*x4^2*x5 - 88*x0*x5^2 + 48*x1*x5^2 + 188*x2*x5^2 + 322*x3*x5^2 + 148*x4*x5^2 - 76*x5^3, 100*x1^2*x2 - 40*x1*x2*x3 + 80*x2^2*x3 - 30*x0*x3^2 + 40*x1*x3^2 - 95*x2*x3^2 + 45*x3^3 + 754*x0*x2*x4 - 536*x1*x2*x4 + 618*x2^2*x4 - 903*x0*x3*x4 + 711*x1*x3*x4 - 1444*x2*x3*x4 + 810*x3^2*x4 - 811*x0*x4^2 - 125*x1*x4^2 - 561*x2*x4^2 + 916*x3*x4^2 + 343*x4^3 + 50*x1^2*x5 - 344*x0*x2*x5 + 246*x1*x2*x5 - 248*x2^2*x5 + 373*x0*x3*x5 - 281*x1*x3*x5 + 564*x2*x3*x5 - 315*x3^2*x5 + 480*x0*x4*x5 + 936*x1*x4*x5 + 112*x2*x4*x5 - 622*x3*x4*x5 - 834*x4^2*x5 - 49*x0*x5^2 - 371*x1*x5^2 + 49*x2*x5^2 + 106*x3*x5^2 + 639*x4*x5^2 - 148*x5^3, 800*x1^3 + 10*x1*x2*x3 + 280*x2^2*x3 - 80*x0*x3^2 - 10*x1*x3^2 - 295*x2*x3^2 + 95*x3^3 + 2444*x0*x2*x4 - 2296*x1*x2*x4 - 1052*x2^2*x4 - 1733*x0*x3*x4 + 1221*x1*x3*x4 + 91*x2*x3*x4 - 15*x3^2*x4 + 104*x0*x4^2 + 600*x1*x4^2 - 1496*x2*x4^2 + 2351*x3*x4^2 + 348*x4^3 - 100*x1^2*x5 - 1024*x0*x2*x5 + 816*x1*x2*x5 + 192*x2^2*x5 + 833*x0*x3*x5 - 551*x1*x3*x5 + 119*x2*x3*x5 - 40*x3^2*x5 + 190*x0*x4*x5 + 196*x1*x4*x5 + 642*x2*x4*x5 - 1977*x3*x4*x5 - 904*x4^2*x5 - 104*x0*x5^2 - 216*x1*x5^2 + 104*x2*x5^2 + 376*x3*x5^2 + 764*x4*x5^2 - 208*x5^3]'


    def test__get_surf( self ):

        a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]  # i
        bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
        bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )  # e1
        bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )  # e2

        ls = LinearSeries.get( [2, 2], bp_tree )  # |2(l1+l2)-e1-e2|
        sig = ( 5, 1 )
        coef_lst = [-1, 1, -1, 1, 0, 1, 1, -1]
        prv_Q = [( 1, 0, 0, 1, 0, 1, 1 ), ( 1, 0, 1, 0, 0, 1, 0 ), ( 1, 0, 0, 1, 0, 0, 0 ), ( 1, 0, 0, 0, 0, 1, 1 ), ( 1, 0, 0, 0, 0, 0, 1 ), ( 0, 1, 1, 0, 1, 1, 0 )]

        dct = get_surf( ls, sig, coef_lst, prv_Q )

        print( list( dct['Q'] ) )
        print( dct['pmz_lst'] )
        print( dct['imp_lst'] )
        print( list( dct['M'] ) )
        print( list( dct['UJ'][0] ) )
        print( list( dct['UJ'][1] ) )

        assert str( list( dct['Q'] ) ) == str( prv_Q )
        assert str( dct['pmz_lst'] ) == '[-c0*s0*c1^2 + c0^2*s1^2 - c0*s0*s1^2 + s0^2*s1^2 + c0*c1^2 - 2*c0^2*s1 + 2*c0*s0*s1 - 2*s0^2*s1 + c0*s1^2 - 2*s0*s1^2 + c0^2 - c0*s0 + s0^2 - 2*c0*s1 + 4*s0*s1 + s1^2 + c0 - 2*s0 - 2*s1 + 1, -c0*s0*c1^2 + s0^2*c1^2 - c0^2*c1*s1 + s0^2*s1^2 + c0^2*c1 + c0*c1^2 - 2*s0*c1^2 - 2*s0^2*s1 - 2*s0*s1^2 + s0^2 + c1^2 + 4*s0*s1 + s1^2 - 2*s0 - 2*s1 + 1, -c0^2*c1^2 + c0^2*c1*s1 - c0*s0*s1^2 + s0^2*s1^2 - c0^2*c1 + 2*c0*s0*s1 - 2*s0^2*s1 + c0*s1^2 - 2*s0*s1^2 - c0*s0 + s0^2 - 2*c0*s1 + 4*s0*s1 + s1^2 + c0 - 2*s0 - 2*s1 + 1, -c0*s0*c1^2 - c0^2*c1*s1 + c0^2*s1^2 + s0^2*s1^2 + c0^2*c1 + c0*c1^2 - 2*c0^2*s1 - 2*s0^2*s1 - 2*s0*s1^2 + c0^2 + s0^2 + 4*s0*s1 + s1^2 - 2*s0 - 2*s1 + 1, c0^2*s1^2 + s0^2*s1^2 - 2*c0^2*s1 - 2*s0^2*s1 - 2*s0*s1^2 + c0^2 + s0^2 + 4*s0*s1 + s1^2 - 2*s0 - 2*s1 + 1, -c0*s0*c1^2 + s0^2*c1^2 - 2*c0^2*c1*s1 + c0*s0*c1*s1 - s0^2*c1*s1 + 2*c0^2*c1 - c0*s0*c1 + s0^2*c1 + c0*c1^2 - 2*s0*c1^2 - c0*c1*s1 + 2*s0*c1*s1 + c0*c1 - 2*s0*c1 + c1^2 - c1*s1 + c1]'
        assert str( dct['imp_lst'] ) == '[2*x1^2 + x0*x2 - x1*x2 - x2*x3 - 4*x0*x4 - 5*x1*x4 + 4*x2*x4 + 3*x3*x4 + x4^2 - x0*x5 - 4*x1*x5 + x2*x5 + x3*x5 + 3*x4*x5 + 2*x5^2, 2*x0*x1 - x0*x2 - x1*x2 - 6*x0*x3 - 2*x1*x3 + 3*x2*x3 + 4*x3^2 + 6*x0*x4 + x1*x4 - 2*x2*x4 - 5*x3*x4 + x4^2 + x0*x5 - x2*x5 - x3*x5 + x4*x5, 2*x0^2 - x0*x2 + x1*x2 - 4*x0*x3 + x2*x3 + 2*x3^2 - x1*x4 - x3*x4 + x4^2 + x0*x5 - x2*x5 - x3*x5 + x4*x5, 4*x1*x3^2 - x2*x3^2 + x3^3 + x0*x2*x4 + x1*x2*x4 - x2^2*x4 + x0*x3*x4 - 11*x1*x3*x4 + x2*x3*x4 - 10*x3^2*x4 - 2*x0*x4^2 + 6*x1*x4^2 + 22*x3*x4^2 - 12*x4^3 - x0*x3*x5 + 2*x2*x3*x5 - 4*x3^2*x5 - x0*x4*x5 + 15*x3*x4*x5 - 11*x4^2*x5 - 2*x0*x5^2 - x2*x5^2 + 2*x3*x5^2 + x4*x5^2, 4*x0*x3^2 - x2*x3^2 - 3*x3^3 + x0*x2*x4 + x1*x2*x4 - x2^2*x4 - 11*x0*x3*x4 + x1*x3*x4 + x2*x3*x4 + 10*x3^2*x4 + 6*x0*x4^2 - 2*x1*x4^2 - 10*x3*x4^2 + 4*x4^3 - 5*x0*x3*x5 + 2*x2*x3*x5 + 4*x3^2*x5 + 7*x0*x4*x5 - 4*x2*x4*x5 - 9*x3*x4*x5 + 5*x4^2*x5 + 2*x0*x5^2 - x2*x5^2 - 2*x3*x5^2 + x4*x5^2, 4*x1*x2*x3 - x2*x3^2 + x3^3 - 3*x0*x2*x4 - 7*x1*x2*x4 + 3*x2^2*x4 - 3*x0*x3*x4 - 3*x1*x3*x4 + 5*x2*x3*x4 - 2*x3^2*x4 + 6*x0*x4^2 + 6*x1*x4^2 - 4*x2*x4^2 + 2*x3*x4^2 - 4*x4^3 - 2*x0*x2*x5 - 2*x1*x2*x5 - x0*x3*x5 + 3*x0*x4*x5 + 2*x1*x4*x5 + 4*x2*x4*x5 + x3*x4*x5 - 5*x4^2*x5 + x2*x5^2 - x4*x5^2, 4*x0*x2*x3 - 5*x2*x3^2 + x3^3 - 7*x0*x2*x4 - 3*x1*x2*x4 + 3*x2^2*x4 - 3*x0*x3*x4 - 3*x1*x3*x4 + 17*x2*x3*x4 - 2*x3^2*x4 + 6*x0*x4^2 + 6*x1*x4^2 - 12*x2*x4^2 + 2*x3*x4^2 - 4*x4^3 - 2*x0*x2*x5 + 2*x1*x2*x5 + 7*x0*x3*x5 - 4*x3^2*x5 - 9*x0*x4*x5 - 2*x1*x4*x5 + 4*x2*x4*x5 + 9*x3*x4*x5 - 5*x4^2*x5 - 4*x0*x5^2 + x2*x5^2 + 4*x3*x5^2 - x4*x5^2, 2*x2^2*x3^2 - 4*x2*x3^3 + 2*x3^4 - 2*x0*x2^2*x4 - 2*x1*x2^2*x4 + 2*x2^3*x4 + 4*x2^2*x3*x4 + 19*x2*x3^2*x4 - 7*x3^3*x4 - x0*x2*x4^2 - x1*x2*x4^2 - 3*x2^2*x4^2 - 3*x0*x3*x4^2 - 3*x1*x3*x4^2 - 25*x2*x3*x4^2 + 4*x3^2*x4^2 + 6*x0*x4^3 + 6*x1*x4^3 + 8*x2*x4^3 + 10*x3*x4^3 - 12*x4^4 - 4*x2^2*x3*x5 + 4*x2*x3^2*x5 - 4*x0*x2*x4*x5 + 4*x2^2*x4*x5 + 3*x0*x3*x4*x5 - 2*x2*x3*x4*x5 - 4*x3^2*x4*x5 - 3*x0*x4^2*x5 + 6*x2*x4^2*x5 + 13*x3*x4^2*x5 - 13*x4^3*x5 + x0*x2*x5^2 - x1*x2*x5^2 + 2*x2^2*x5^2 - x2*x3*x5^2 - 6*x0*x4*x5^2 + x1*x4*x5^2 - x2*x4*x5^2 + 7*x3*x4*x5^2 - 2*x4^2*x5^2 - x0*x5^3 + x2*x5^3 + x3*x5^3 - x4*x5^3]'
        assert str( list( dct['M'] ) ) == '[(-2, 1, -1/2, -1, 5, 1/2), (1, -2, -1/2, -1, 7/2, 2), (-1/2, -1/2, 0, 3/2, -3, -1/2), (-1, -1, 3/2, 2, -7/2, -1/2), (5, 7/2, -3, -7/2, -1, -3/2), (1/2, 2, -1/2, -1/2, -3/2, -2)]'
        assert str( list( dct['UJ'][0] ) ) == '[(1.35786975221142?, 1.07346104931619?, -0.55476539365996?, -0.40827309297366?, -2.0166161907848?, -0.97578089340802?), (1.164220993205174?, -1.238241275959885?, -0.069749665334111?, -0.0755220593640689?, -0.253358176151464?, 0.852764409706635?), (0.2374970415389016?, 0.06861265576467416?, 0.8393349903238033?, -0.2501991469577176?, 0.08951234160953201?, -0.1515231571929268?), (0.05097640435508013?, -0.3562824973464569?, -0.0953003146422521?, 0.01449254245924731?, 0.1356081690977796?, -0.5531500949878551?), (0.0804899878314170?, 0.04919673465402113?, 0.00930980680657234?, 0.1455320285330222?, 0.0527745414065528?, -0.00912281285573814?), (1.1245900860635?, 0.9053316050139?, -1.0130556284196?, -1.5021250674132?, 1.7776721895235?, 0.0915041576522?)]'
        assert str( list( dct['UJ'][1] ) ) == '[(-1.000000000000000?, 0, 0, 0, 0, 0), (0, -1.000000000000000?, 0, 0, 0, 0), (0, 0, -1.000000000000000?, 0, 0, 0), (0, 0, 0, -1.000000000000000?, 0, 0), (0, 0, 0, 0, -1.000000000000000?, 0), (0, 0, 0, 0, 0, 1.000000000000000?)]'


    def test__verify_get_surf( self ):

        a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]  # i
        bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
        bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )  # e1
        bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )  # e2

        ls = LinearSeries.get( [2, 2], bp_tree )  # |2(l1+l2)-e1-e2|
        sig = ( 6, 1 )
        coef_lst = [-1, -1, 0, 0, 0, -1, 1, 0, -1, -1, -1]
        prv_Q = None
        dct = get_surf( ls, sig, coef_lst, prv_Q )

        test_dct = verify_get_surf( dct )
        print( test_dct )

        assert test_dct['all']


    def test__get_proj( self ):

        # construct dct
        a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]  # i
        bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
        bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )  # e1
        bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )  # e2
        ls = LinearSeries.get( [2, 2], bp_tree )  # |2(l1+l2)-e1-e2|
        sig = ( 6, 1 )
        coef_lst = [-1, -1, 0, 0, 0, -1, 1, 0, -1, -1, -1]
        prv_Q = None
        dct = get_surf( ls, sig, coef_lst, prv_Q )

        # construct projection matrix P
        U, J = dct['UJ']
        U.swap_rows( 0, 6 )
        J.swap_columns( 0, 6 )
        J.swap_rows( 0, 6 )
        assert dct['M'] == approx_QQ( U.T * J * U )
        approxU = approx_QQ( U )
        P = get_prj_mat( 4, 7, 0 )
        P[0, 6] = -1;P[3, 3] = 0;P[3, 4] = 1
        P = P * approxU

        # call get_proj
        f_xyz, pmz_AB_lst = get_proj( dct['imp_lst'], dct['pmz_lst'], P )

        print( f_xyz )
        print( pmz_AB_lst )

        assert str( f_xyz ).startswith( '232136158126691416007724114435018582183550370830165237' )
        assert str( pmz_AB_lst ).startswith( '[3687965761301911367110366337978449028075575/' )


if __name__ == '__main__':

    # TestSurfaceInQuadric().test__get_S1xS1_pmz()
    # TestSurfaceInQuadric().test__get_deg_surf()
    # TestSurfaceInQuadric().test__get_prj_mat()
    # TestSurfaceInQuadric().test__approx_QQ()
    # TestSurfaceInQuadric().test__rand_surf_prj()
    # TestSurfaceInQuadric().test__get_surf()
    # TestSurfaceInQuadric().test__verify_get_surf()
    TestSurfaceInQuadric().test__get_proj()


