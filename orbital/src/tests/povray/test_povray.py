'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''
import os

from orbital.sage_interface import sage_var
from orbital.sage_interface import sage_pi
from orbital.sage_interface import sage_QQ

from orbital.povray.class_pov_input import PovInput

from orbital.povray.povray import create_pov

from orbital.povray.povray_aux import get_time_str

from orbital.povray.povray_aux import rgbt2pov


class TestPovray( object ):

    def test__rgbt2_pov( self ):
        # blue
        col = rgbt2pov( ( 28, 125, 154, 0 ) )
        print( col )
        assert str( col ) == '(0.00775102739766061, 0.20835965596076741, 0.3297290329675149, 0.0)'

        # beige
        col = rgbt2pov( ( 254, 242, 190, 0 ) )
        print( col )
        assert str( col ) == '(0.9913928435929399, 0.8912620368134188, 0.5234431552143247, 0.0)'


    def test__povray( self ):

        c0, s0, c1, s1 = sage_var( 'c0,s0,c1,s1' )
        x, y, z = sage_var( 'x,y,z' )
        r = 1
        R = 2
        pmz_AB_lst = [1, ( c0 * r + R ) * c1, ( c0 * r + R ) * s1, r * s0]
        f_xyz = ( x ** 2 + y ** 2 + z ** 2 + R ** 2 - r ** 2 ) ** 2 - 4 * R ** 2 * ( x ** 2 + y ** 2 )

        v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 150 )]
        v1_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 300 )]

        # set PovInput as container
        # put very low quality for testing purposes
        pin = PovInput()
        pin.path = './' + get_time_str() + '_TEST_POVRAY_REMOVE_ME/'
        pin.fname = 'orb'
        pin.scale = 1
        pin.cam_dct['location'] = ( 0, -7, 0 )
        pin.cam_dct['lookat'] = ( 0, 0, 0 )
        pin.cam_dct['rotate'] = ( 45, 0, 0 )
        pin.light_radius = 5
        pin.axes_dct['show'] = True
        pin.axes_dct['len'] = 3
        pin.width = 2
        pin.height = 2
        pin.quality = 1
        pin.ani_delay = 1

        pin.impl = f_xyz

        pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
        pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )

        pin.curve_dct['A'] = {'step0':v0_lst, 'step1': v1_lst, 'prec':1, 'width':0.08}
        pin.curve_dct['B'] = {'step0':v0_lst, 'step1': v1_lst, 'prec':1, 'width':0.08}

        pin.text_dct['A'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
        pin.text_dct['B'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
        pin.text_dct['SURF'] = [True, ( 0.2, 0.7, 0.3, 0.0 ), 'F_Glass10']

        # raytrace image/animation
        show_surf = True
        ani = False
        ft_lst = []

        lst = create_pov( pin, ['A', 'B'], show_surf, ani, ft_lst )
        lst = create_pov( pin, ['A'], False, True, [] )


if __name__ == '__main__':

    TestPovray().test__rgbt2_pov()
    # TestPovray().test__povray()
    pass
