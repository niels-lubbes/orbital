'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Apr 4, 2018
@author: Niels Lubbes
'''

from orbital.prod.orb_matrices import get_tmat
from orbital.prod.orb_matrices import get_omat
from orbital.prod.orb_matrices import get_rmat
from orbital.prod.orb_matrices import get_pmat
from orbital.prod.orb_matrices import get_emat
from orbital.prod.orb_matrices import get_xmat
from orbital.prod.orb_matrices import get_mat

from orbital.class_orb_ring import OrbRing

from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_identity_matrix
from orbital.sage_interface import sage_Permutation


class TestOrbMatrices( object ):


    def test__get_tmat( self ):

        #
        # Setup inverse stereographic projection:
        # S: P^7 ---> S^7
        #
        v = OrbRing.coerce( '[v0,v1,v2,v3,v4,v5,v6,v7,v8]' )
        d = v[1] ** 2 + v[2] ** 2 + v[3] ** 2
        v0_2 = v[0] ** 2
        vec_lst = []
        vec_lst += [v0_2 + d]
        vec_lst += [2 * v[0] * v[1]]
        vec_lst += [2 * v[0] * v[2]]
        vec_lst += [2 * v[0] * v[3]]
        vec_lst += [2 * v[0] * v[4]]
        vec_lst += [2 * v[0] * v[5]]
        vec_lst += [2 * v[0] * v[6]]
        vec_lst += [2 * v[0] * v[7]]
        vec_lst += [-v0_2 + d]
        x = sage_vector( OrbRing.R, vec_lst )
        print( x )

        #
        # Setup Euclidean translation in S^7
        #    T:S^7--->S^7
        #
        tmat = get_tmat( None )
        print( tmat )

        #
        # Setup stereographic projection
        #     P:S^7--->P^7
        #
        pmat = []
        pmat += [[1, 0, 0, 0] + [ 0, 0, 0, 0, -1]]
        pmat += [[0, 1, 0, 0] + [ 0, 0, 0, 0, 0]]
        pmat += [[0, 0, 1, 0] + [ 0, 0, 0, 0, 0]]
        pmat += [[0, 0, 0, 1] + [ 0, 0, 0, 0, 0]]
        pmat += [[0, 0, 0, 0] + [ 1, 0, 0, 0, 0]]
        pmat += [[0, 0, 0, 0] + [ 0, 1, 0, 0, 0]]
        pmat += [[0, 0, 0, 0] + [ 0, 0, 1, 0, 0]]
        pmat += [[0, 0, 0, 0] + [ 0, 0, 0, 1, 0]]
        pmat = sage_matrix( pmat )
        print( pmat )

        # Check whether composition is an Euclidean
        # translation in P^7:
        #    PoToS
        #
        tv = pmat * tmat * x
        tv = tv / ( 2 * v[0] )
        print( tv )

        assert str( tv ) == '(v0, v0*t1 + v1, v0*t2 + v2, v0*t3 + v3, v0*t4 + v4, v0*t5 + v5, v0*t6 + v6, v0*t7 + v7)'


    def test__get_omat_1( self ):
        out = get_omat( 'Orsma' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, c0, -s0, 0, 0, 0, 0, 0, 0), (0, s0, c0, 0, 0, 0, 0, 0, 0), (0, 0, 0, -c0, s0, 0, 0, 0, 0), (0, 0, 0, -s0, -c0, 0, 0, 0, 0), (0, 0, 0, 0, 0, -1, 0, 0, 0), (0, 0, 0, 0, 0, 0, -1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, -1)]'


    def test__get_omat_2( self ):
        out = get_omat( 'Oprpp' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0), (0, 0, 0, c0, -s0, 0, 0, 0, 0), (0, 0, 0, s0, c0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_rmat_1( self ):
        out = get_rmat( 'Rprpp[0,90,0,0]' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, -1, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_rmat_2( self ):
        out = get_rmat( 'Rprpp[0,45,0,0]' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 119/169, -120/169, 0, 0, 0, 0), (0, 0, 0, 120/169, 119/169, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_pmat__True( self ):
        out = get_pmat( True )
        print( out )
        for elt in out.list():
            assert elt in [0, 1]

    def test__get_pmat__False( self ):
        out = get_pmat( False )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, -1), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0, 0, 0, 0)]'


    def test__get_emat__E_21435678( self ):
        out = get_emat( 'E[2,1,4,3,5,6,7,8]' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_emat__E_random( self ):
        out = get_emat( 'E' )
        print( out )
        print( list( out ) )
        for row in out:
            assert sum( row ) == 1

    def test__get_xmat__90_0_0( self ):
        out = get_xmat( 'X[90,0,0]' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, -1, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_xmat__0_90_0( self ):
        out = get_xmat( 'X[0,90,0]' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, -1, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_xmat__0_0_90( self ):
        out = get_xmat( 'X[0,0,90]' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, -1, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_mat__I_Orppp_I( self ):
        out = get_mat( 'I', 'Orppp', 'I' )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, c0, -s0, 0, 0, 0, 0, 0, 0), (0, s0, c0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'


    def test__get_mat_1( self ):
        #
        # Permute a block matrices to obtain rotation along x-axis.
        #
        print( get_mat( 'I', 'Oprpp', 'I' ) )


        e_lst = [1, 4, 2, 3, 5, 6, 7, 8]
        i_lst = sage_Permutation( e_lst ).inverse()

        A_str = 'E' + str( i_lst )
        B_str = 'Oprpp'
        C_str = 'E' + str( e_lst )
        tup = ( A_str, B_str, C_str )
        print( tup )

        mati = get_emat( A_str )
        mate = get_emat( C_str )

        assert mate * mati == sage_identity_matrix( 9 )

        out = get_mat( *tup )
        print( out )
        print( list( out ) )
        assert str( list( out ) ) == '[(1, 0, 0, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, c0, -s0, 0, 0, 0, 0, 0), (0, 0, s0, c0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0, 1)]'




if __name__ == '__main__':

#     TestOrbMatrices().test__get_tmat()
#     TestOrbMatrices().test__get_omat_1()
#     TestOrbMatrices().test__get_omat_2()
#     TestOrbMatrices().test__get_rmat_1()
#     TestOrbMatrices().test__get_rmat_2()
#     TestOrbMatrices().test__get_pmat__True()
#     TestOrbMatrices().test__get_pmat__False()
#     TestOrbMatrices().test__get_emat__E_21435678()
#     TestOrbMatrices().test__get_emat__E_random()
#     TestOrbMatrices().test__get_xmat__90_0_0()
#     TestOrbMatrices().test__get_xmat__0_90_0()
#     TestOrbMatrices().test__get_xmat__0_0_90()
#     TestOrbMatrices().test__get_mat__I_Orppp_I()
    TestOrbMatrices().test__get_mat_1()
    pass
