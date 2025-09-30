'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Aug 30, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_identity_matrix
from orbital.sage_interface import sage_PolynomialRing
from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_var

from orbital.cossin.cos_sin import get_cs

from orbital.sphere.sphere_transform import get_rot_mat
from orbital.sphere.sphere_transform import get_rot_S3
from orbital.sphere.sphere_transform import get_trn_S3
from orbital.sphere.sphere_transform import get_scale_S3
from orbital.sphere.sphere_transform import get_hp_P4


class TestSphereTransform( object ):

    def test__get_rot_mat( self ):

        angle = 45
        c, s = get_cs( angle )
        out = list( get_rot_mat( 3, 0, 1, angle ) )
        chk = [( c, -s, 0 ), ( s, c, 0 ), ( 0, 0, 1 )]
        print( out )
        print( chk )
        assert str( out ) == str( chk )

    def test__get_rot_S3( self ):

        a12, a13, a23, a14, a24, a34 = 6 * [0]
        a23 = 45
        c, s = get_cs( a23 )

        out = get_rot_S3( [a12, a13, a23, a14, a24, a34] )
        print( 'out =' )
        print( out )

        chk = sage_identity_matrix( sage_QQ, 5 )
        chk[2, 2] = c;chk[2, 3] = -s
        chk[3, 2] = s;chk[3, 3] = c;
        print( 'chk =' )
        print( chk )

        assert str( list( out ) ) == str( list( chk ) )

    def test__get_trn_S3( self ):

        R = sage_PolynomialRing( sage_QQ, 't0,t1,t2,t3,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4' )
        t0, t1, t2, t3, x0, x1, x2, x3, x4, y0, y1, y2, y3, y4 = R.gens()

        #
        # We obtain p:S^3 ---> S^3 by composing
        # stereographic projection S^3--->P^3, x |---> (x0-x4:x1:x2:x3)
        # translation P^3--->P^3, and
        # inverse stereographic projection P^3--->S^3.
        #
        z0 = t0 * ( x0 - x4 )
        z1 = t0 * x1 + t1 * ( x0 - x4 )
        z2 = t0 * x2 + t2 * ( x0 - x4 )
        z3 = t0 * x3 + t3 * ( x0 - x4 )
        ZZ = z1 ** 2 + z2 ** 2 + z3 ** 2
        p = [ z0 ** 2 + ZZ, 2 * z0 * z1, 2 * z0 * z2, 2 * z0 * z3, -z0 ** 2 + ZZ]
        #
        # We look at the coefficients of p with the following code
        #
        # for i in [0, 1, 2, 3, 4]:
        #     print '\np[', i, ']'
        #     for cf in [t0 ** 2, t1 ** 2, t2 ** 2, t3 ** 2, t0 * t1, t0 * t2, t0 * t3, t1 * t2, t1 * t3, t2 * t3 ]:
        #         pcf = p[i].coefficient( cf )
        #         if pcf == 0: continue
        #         pcf = sage_factor( pcf )
        #         print( cf, ':', pcf )
        #
        # this leads to the following output
        #
        # p[ 0 ]
        # (t0^2, ':', x0^2 + x1^2 + x2^2 + x3^2 - 2*x0*x4 + x4^2)
        # (t1^2, ':', (x0 - x4)^2)
        # (t2^2, ':', (x0 - x4)^2)
        # (t3^2, ':', (x0 - x4)^2)
        # (t0*t1, ':', (2) * x1 * (x0 - x4))
        # (t0*t2, ':', (2) * x2 * (x0 - x4))
        # (t0*t3, ':', (2) * x3 * (x0 - x4))
        #
        # p[ 1 ]
        # (t0^2, ':', (2) * x1 * (x0 - x4))
        # (t0*t1, ':', (2) * (x0 - x4)^2)
        #
        # p[ 2 ]
        # (t0^2, ':', (2) * x2 * (x0 - x4))
        # (t0*t2, ':', (2) * (x0 - x4)^2)
        #
        # p[ 3 ]
        # (t0^2, ':', (2) * x3 * (x0 - x4))
        # (t0*t3, ':', (2) * (x0 - x4)^2)
        #
        # p[ 4 ]
        # (t0^2, ':', -x0^2 + x1^2 + x2^2 + x3^2 + 2*x0*x4 - x4^2)
        # (t1^2, ':', (x0 - x4)^2)
        # (t2^2, ':', (x0 - x4)^2)
        # (t3^2, ':', (x0 - x4)^2)
        # (t0*t1, ':', (2) * x1 * (x0 - x4))
        # (t0*t2, ':', (2) * x2 * (x0 - x4))
        # (t0*t3, ':', (2) * x3 * (x0 - x4))
        #

        #
        # We use -x0^2+x1^2+x2^2+x3^2+x4^2==0
        # to simplify the coefficients so that we
        # obtain the map q:S^3--->S^3.
        #
        TT = ( sage_QQ( 1 ) / 2 ) * ( t1 ** 2 + t2 ** 2 + t3 ** 2 )
        XT = x1 * t1 + x2 * t2 + x3 * t3
        q = [
            x0 + ( x0 - x4 ) * TT + XT,
            x1 + ( x0 - x4 ) * t1,
            x2 + ( x0 - x4 ) * t2,
            x3 + ( x0 - x4 ) * t3,
            x4 + ( x0 - x4 ) * TT + XT
            ]

        #
        # Obtain 5x5 matrix Mq of linear transformation q.
        #
        mat = []
        for i in [0, 1, 2, 3, 4]:
            row = []
            for cf in [x0, x1, x2, x3, x4]:
                qcf = q[i].coefficient( cf )
                row += [ qcf ]
            mat += [row]
        Mq = sage_matrix( mat )
        T = get_trn_S3( [t1, t2, t3] )
        print( Mq )
        assert Mq == T

        #
        # Check how the map acts on P^3
        #
        YY = y1 ** 2 + y2 ** 2 + y3 ** 2
        u = Mq * sage_vector( [y0 ** 2 + YY, 2 * y0 * y1, 2 * y0 * y2, 2 * y0 * y3, -y0 ** 2 + YY] )
        u = [ u[0] - u[4], u[1], u[2], u[3]]
        assert '[2*y0^2, 2*t1*y0^2 + 2*y0*y1, 2*t2*y0^2 + 2*y0*y2, 2*t3*y0^2 + 2*y0*y3]' == str( u )

        #
        # Some additional sanity checks
        #
        Q = T.subs( {t1:1, t2:2, t3:5} )
        v = Q * sage_vector( [1, 0, 0, 0, -1] )
        assert -v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + v[3] ** 2 + v[4] ** 2 == 0
        assert Q * sage_vector( [1, 0, 0, 0, 1] ) == sage_vector( [1, 0, 0, 0, 1] )

    def test__get_scale_S3( self ):

        R = sage_PolynomialRing( sage_QQ, 's,t,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4' )
        s, t, x0, x1, x2, x3, x4, y0, y1, y2, y3, y4 = R.gens()

        z0 = t * ( x0 - x4 )
        z1 = s * x1
        z2 = s * x2
        z3 = s * x3
        ZZ = z1 ** 2 + z2 ** 2 + z3 ** 2
        p = [ z0 ** 2 + ZZ, 2 * z0 * z1, 2 * z0 * z2, 2 * z0 * z3, -z0 ** 2 + ZZ]
        #
        # We look at the coefficients of p with the following code
        #
        # for i in [0, 1, 2, 3, 4]:
        #     print '\np[', i, ']'
        #     for cf in [s ** 2, s * t, t ** 2]:
        #         pcf = p[i].coefficient( cf )
        #         if pcf == 0: continue
        #         pcf = sage_factor( pcf )
        #         print( cf, ':', pcf )
        #
        # this leads to the following output
        #
        # p[ 0 ]
        # (s^2, ':', x1^2 + x2^2 + x3^2)
        # (t^2, ':', (x0 - x4)^2)
        #
        # p[ 1 ]
        # (s*t, ':', (2) * x1 * (x0 - x4))
        #
        # p[ 2 ]
        # (s*t, ':', (2) * x2 * (x0 - x4))
        #
        # p[ 3 ]
        # (s*t, ':', (2) * x3 * (x0 - x4))
        #
        # p[ 4 ]
        # (s^2, ':', x1^2 + x2^2 + x3^2)
        # (t^2, ':', (-1) * (x0 - x4)^2)

        #
        # We use -x0^2+x1^2+x2^2+x3^2+x4^2==0
        # to simplify the coefficients so that we
        # obtain the map q:S^3--->S^3. We set t=1.
        #
        q = [ ( x0 + x4 ) * s ** 2 + ( x0 - x4 ),
              2 * x1 * s,
              2 * x2 * s,
              2 * x3 * s,
              ( x0 + x4 ) * s ** 2 - ( x0 - x4 )]

        #
        # Obtain 5x5 matrix Mq of linear transformation q.
        #
        Mq = []
        for i in [0, 1, 2, 3, 4]:
            row = []
            for cf in [x0, x1, x2, x3, x4]:
                qcf = q[i].coefficient( cf )
                row += [ qcf ]
            Mq += [row]
        Mq = sage_matrix( Mq )
        S = get_scale_S3( s )
        print( Mq )
        assert Mq == S

        #
        # Check how the map acts on P^3
        #
        YY = y1 ** 2 + y2 ** 2 + y3 ** 2
        u = Mq * sage_vector( [y0 ** 2 + YY, 2 * y0 * y1, 2 * y0 * y2, 2 * y0 * y3, -y0 ** 2 + YY] )
        u = [ u[0] - u[4], u[1], u[2], u[3]]
        assert '[4*y0^2, 4*s*y0*y1, 4*s*y0*y2, 4*s*y0*y3]' == str( u )

        #
        # Some additional sanity checks
        #
        Q = S.subs( {s:2} )
        v = Q * sage_vector( [1, 0, 1, 0, 0] )
        w = Q * sage_vector( [1, 0, 0, 0, 1] )
        assert -v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + v[3] ** 2 + v[4] ** 2 == 0
        assert list( w ) == [8, 0, 0, 0, 8]

    def test__get_hp_P4( self ):

        c0, s0, c1, s1 = sage_var( 'c0,s0,c1,s1' )
        v = sage_vector( [1, c0, s0, 0, 0] )
        w = sage_vector( [1, c1, s1, 0, 0] )

        a12, a13, a23, a14, a24, a34 = [90, 0, 0, 0, 0, 0]
        M = get_rot_S3( [a12, a13, a23, a14, a24, a34] )
        out = get_hp_P4( v, M * w )
        print( out )
        assert str( out ) == '[1, -c1*s0 - c0*s1, c0*c1 - s0*s1, 0, 0]'


if __name__ == '__main__':

    TestSphereTransform().test__get_rot_mat()
    TestSphereTransform().test__get_rot_S3()
    TestSphereTransform().test__get_trn_S3()
    TestSphereTransform().test__get_scale_S3()
    TestSphereTransform().test__get_hp_P4()
    pass

