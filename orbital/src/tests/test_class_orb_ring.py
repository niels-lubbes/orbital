'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''

from orbital.class_orb_ring import OrbRing


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



if __name__ == '__main__':

    TestOrbRing().test__coerce()
    TestOrbRing().test__random_int()
    TestOrbRing().test__random_elt()
    pass
