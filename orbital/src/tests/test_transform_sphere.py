'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 22, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_identity_matrix
from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_var

from orbital.cossin.cos_sin import get_cs

from orbital.transform_sphere import get_rot_mat
from orbital.transform_sphere import get_rot_S3
from orbital.transform_sphere import get_trn_S3
from orbital.transform_sphere import get_xfer_S3
from orbital.transform_sphere import get_hp_S3
from orbital.transform_sphere import get_prj_S3


class TestTransformSphere( object ):

    def test__get_rot_mat( self ):

        angle = 2
        c, s = get_cs( angle )
        out = list( get_rot_mat( 3, 0, 1, angle ) )
        chk = [( c, -s, 0 ), ( s, c, 0 ), ( 0, 0, 1 )]
        print( out )
        print( chk )
        assert str( out ) == str( chk )


    def test__get_rot_S3( self ):

        a01, a02, a03, a12, a13, a23 = 6 * [0]
        a23 = 2
        c, s = get_cs( a23 )

        out = get_rot_S3( a01, a02, a03, a12, a13, a23 )
        print( 'out =' )
        print( out )

        chk = sage_identity_matrix( sage_QQ, 5 )
        chk[2, 2] = c;chk[2, 3] = -s
        chk[3, 2] = s;chk[3, 3] = c;
        print( 'chk =' )
        print( chk )

        assert str( list( out ) ) == str( list( chk ) )


    def test__get_trn_S3( self ):

        x0, y0, z0, s = 1, 2, 3, 1
        out = get_trn_S3( x0, y0, z0, s )
        print( out )

        chk = [( 2, 0, 0, -2, 2 ), ( 0, 2, 0, -4, 4 ), ( 0, 0, 2, -6, 6 ), ( 2, 4, 6, -12, 14 ), ( 2, 4, 6, -14, 16 )]
        assert str( list( out ) ) == str( chk )


    def test__get_xfer_S3( self ):

        x0, y0, z0, s = 1, 2, 3, 1
        a01, a02, a03, a12, a13, a23 = 5 * [0] + [2]

        out = get_xfer_S3( a01, a02, a03, a12, a13, a23, x0, y0, z0, s )
        chk = get_trn_S3( x0, y0, z0, s ) * get_rot_S3( a01, a02, a03, a12, a13, a23 )

        assert out == chk


    def test__get_hp_S3( self ):
        a01, a02, a03, a12, a13, a23 = 5 * [0] + [2]
        M = get_rot_S3( a01, a02, a03, a12, a13, a23 )

        c0, s0, c1, s1 = sage_var( 'c0,s0,c1,s1' )
        v = sage_vector( [c0, s0, 0, 0, 1] )
        w = sage_vector( [c1, s1, 0, 0, 1] )

        out = get_hp_S3( v, M * w )

        print( out )

        assert str( out ) == '(c0*c1 - s0*s1, c1*s0 + c0*s1, 0, 0, 1)'


    def test__get_prj_S3( self ):

        x0, x1, x2, x3, x4 = sage_var( 'x0,x1,x2,x3,x4' )
        v = sage_vector( [x1, x2, x3, x4, x0] )
        out = get_prj_S3( v )
        print( out )
        assert out == [-x1 / ( x0 - x4 ), -x2 / ( x0 - x4 ), -x3 / ( x0 - x4 )]




if __name__ == '__main__':

    TestTransformSphere().test__get_rot_mat()
    TestTransformSphere().test__get_rot_S3()
    TestTransformSphere().test__get_trn_S3()
    TestTransformSphere().test__get_xfer_S3()
    TestTransformSphere().test__get_hp_S3()
    TestTransformSphere().test__get_prj_S3()
    pass
