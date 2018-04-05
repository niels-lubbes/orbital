'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 22, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_pi
from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_n
from orbital.sage_interface import sage_cos
from orbital.sage_interface import sage_sin


from orbital.cossin.cos_sin import get_cs


class TestCosSin( object ):

    def test__get_cs( self ):

        for angle in range( 0, 90 ):
            rad = sage_QQ( angle ) / 180 * sage_pi

            cs = ( c, s ) = get_cs( angle )
            csn = ( sage_n( c ), sage_n( s ) )
            csr = ( sage_n( sage_cos( rad ) ), sage_n( sage_sin( rad ) ) )

            dc = abs( cs[0] - csr[0] )
            ds = abs( cs[1] - csr[1] )

            print( cs, csn, csr, ( dc, ds ) )

            assert dc < 0.0086 and ds < 0.0086





if __name__ == '__main__':

    TestCosSin().test__get_cs()
    pass
