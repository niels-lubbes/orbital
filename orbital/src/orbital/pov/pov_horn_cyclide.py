'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 21, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_var
from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_factor
from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_pi

from orbital.class_orb_tools import OrbTools

from orbital.class_pov_input import PovInput

from orbital.povray.povray import create_pov
from orbital.povray.povray_aux import get_time_str


def horn_cyclide():

    # Construct linear series for the horn cyclide, which is
    # an inversion of a circular cone.
    #
    # We construct a trigonometric parametrization
    # by rotating a circle.
    #
    r = 1; R = sage_QQ( 1 ) / 2;
    x, y, v, w = sage_var( 'x,y,v,w' )
    c0, s0, c1, s1 = sage_var( 'c0,s0,c1,s1' )
    V = sage_vector( [r * c0 + R, 0, r * s0] )
    M = sage_matrix( [( c1, -s1, 0 ), ( s1, c1, 0 ), ( 0, 0, 1 )] )
    pmz_AB_lst = [1] + list( M * V )
    OrbTools.p( 'pmz_AB_lst =', pmz_AB_lst )
    for pmz in pmz_AB_lst:
        OrbTools.p( '\t\t', sage_factor( pmz ) )

    # PovInput horn cyclide
    #
    pin = PovInput()

    pin.path = './' + get_time_str() + '_horn_cyclide/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, -4, 0 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 20, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 1
    pin.ani_delay = 10

    pin.impl = None

    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )

    pin.pmz_dct['FA'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['FB'] = ( pmz_AB_lst, 1 )

    v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 10 )]

    v1_lst_A = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 270, 15 )]
    v1_lst_B = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 180, 15 )]

    v1_lst_F = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 1 )]

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst_A, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst_B, 'prec':10, 'width':0.05}
    pin.curve_dct['FA'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.02}
    pin.curve_dct['FB'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.02}

    col_F = ( 0.1, 0.1, 0.1, 0.0 )
    pin.text_dct['A'] = [True, ( 0.4, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.0, 0.2, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FA'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FB'] = [True, col_F, 'phong 0.2 phong_size 5' ]

    # raytrace image/animation
    create_pov( pin, ['A', 'B'] )
    return
    create_pov( pin, ['A', 'B', 'FA', 'FB'] )
