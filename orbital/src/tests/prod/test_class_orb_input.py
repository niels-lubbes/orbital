'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Apr 5, 2018
@author: Niels Lubbes
'''

from orbital.prod.class_orb_input import OrbInput
from orbital.sage_interface import sage__eval


class TestOrbInput( object ):

    def test__random( self ):
        coef_bnd = 3
        in1 = OrbInput().random( 3, False )
        print( in1 )

        print( in1.info_dct['vmat'][0] )
        print( in1.info_dct['vmat'][2] )

        s = ''
        s += in1.info_dct['vmat'][0][1:]
        if in1.info_dct['vmat'][2] != 'I':
            s += '+'
            s += in1.info_dct['vmat'][2][1:]

        for elt in sage__eval( s ):
            assert elt >= -coef_bnd and elt <= coef_bnd

        in2 = OrbInput().set( in1.info_dct['pmat'], in1.info_dct['omat'], in1.info_dct['vmat'] )

        assert in1.pmat == in2.pmat
        assert in1.omat == in2.omat
        assert in1.vmat == in2.vmat


if __name__ == '__main__':

    TestOrbInput().test__random()

    pass
