'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Aug 30, 2018
@author: Niels Lubbes
'''

from orbital.sphere.class_sphere_input import SphereInput


class TestSphereInput( object ):

    def test__random( self ):

        sinp = SphereInput().random( 10 )

        assert [sinp.rota[i] in range( 360 ) for i in range( 6 )]
        assert [sinp.rotb[i] in range( 360 ) for i in range( 6 )]

        assert [sinp.trna[i] in range( 360 ) for i in range( 3 )]
        assert [sinp.trnb[i] in range( 360 ) for i in range( 3 )]

        assert -10 < sinp.sa < 10
        assert -10 < sinp.sb < 10


    def test__set( self ):

        lsta = [( 0, 0, 45 ), ( 0, 0, 0 ), ( 0, 0, 0 ), 1]
        lstb = [( 0, 0, 0 ), ( 0, 0, 0 ), ( 0, 0, 0 ), 1]
        sinp = SphereInput().set( [lsta] + [lstb] )
        print( sinp )



if __name__ == '__main__':

    # TestSphereInput().test__random()
    TestSphereInput().test__set()

