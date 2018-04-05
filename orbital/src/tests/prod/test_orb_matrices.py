'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Apr 4, 2018
@author: Niels Lubbes
'''

from orbital.prod.orb_matrices import get_xmat


class TestOrbMatrices( object ):

    def test__get_xmat( self ):
        print( get_xmat( 'X[0,90,0]' ) )

    def test__get_emat( self ):
        print( get_emat( 'E[2,1,4,3,5,6,7,8]' ) )

    def test__get_omat( self ):
        print( get_omat( 'Oprpp' ) )

    def test__get_rmat_1( self ):
        print( get_rmat( 'Rprpp[0,90,0,0]' ) )

    def test__get_rmat_2( self ):
        print( get_rmat( 'Rprpp[0,45,0,0]' ) )

    def test__get_mat( self ):
        print( get_mat( 'I', 'Orppp', 'I' ) )



if __name__ == '__main__':

    TestOrbMatrices().test__get_xmat()

    pass
