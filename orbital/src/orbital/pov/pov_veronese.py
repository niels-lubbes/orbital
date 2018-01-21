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
from orbital.sage_interface import sage_n

from orbital.class_orb_tools import OrbTools

from orbital.class_orb_ring import OrbRing

from orbital.povray.class_pov_input import PovInput
from orbital.povray.povray import create_pov
from orbital.povray.povray_aux import get_time_str

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries


def veronese():

    # Construct projection of Veronese surface.
    #
    c0, s0, c1, s1 = sage_var( 'c0,s0,c1,s1' )
    x, y = sage_var( 'x,y' )

    pmz_A_lst = [ 1, c0 * s0 * s1, c0 * s0 * c1, c0 * c0 * c1 * s1 ]

    P1 = c0 / ( s0 - 1 )
    P2 = c1 / ( s1 - 1 )
    P3 = ( s0 / c0 ) * ( c1 / ( s1 - 1 ) )

    dct_CD = {x:P1, y:P2 }
    den_CD = ( s0 - 1 ) ** 2 * ( s1 - 1 ) ** 2

    dct_ED = {x:P3, y:P2 }
    den_ED = c0 ** 2 * ( s1 - 1 ) ** 2

    pmz_lst = [ x ** 2 + y ** 2 + 1, -x, -x * y, y ]
    pmz_B_lst = [ ( pmz.subs( dct_CD ) * den_CD ).expand() for pmz in pmz_lst  ]
    pmz_C_lst = [ ( pmz.subs( dct_ED ) * den_ED ).expand() for pmz in pmz_lst  ]

    for A, pmz_lst in [( 'A', pmz_A_lst ), ( 'B', pmz_B_lst ), ( 'C', pmz_C_lst )]:
        OrbTools.p( 'pmz_' + A + '_lst =', pmz_lst )
        for pmz in pmz_lst:
            OrbTools.p( '\t\t', sage_factor( pmz ) )

    # PovInput Veronese surface
    #
    pin = PovInput()

    pin.path = './' + get_time_str() + '_veronese/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, -1.5, 0 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 35, 0, 45 )
    pin.light_radius = 2
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 0.5
    pin.width = 400
    pin.height = 200
    pin.quality = 1
    pin.ani_delay = 10

    pin.impl = None

    pin.pmz_dct['A'] = ( pmz_A_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_B_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_C_lst, 1 )

    pin.pmz_dct['FA'] = ( pmz_A_lst, 0 )
    pin.pmz_dct['FB'] = ( pmz_B_lst, 1 )
    pin.pmz_dct['FC'] = ( pmz_C_lst, 1 )

    v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 5 )]
    v1_A_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 9 )]
    v1_B_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 18 )]
    v1_C_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 9 )]

    v1_lst_F = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 4 )]

    prec = 50

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_A_lst, 'prec':prec, 'width':0.01}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_B_lst, 'prec':prec, 'width':0.01}
    pin.curve_dct['C'] = {'step0':v0_lst, 'step1':v1_C_lst, 'prec':prec, 'width':0.01}

    pin.curve_dct['FA'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.001}
    pin.curve_dct['FB'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.001}
    pin.curve_dct['FC'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.001}

    pin.text_dct['A'] = [True, ( 0.4, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['C'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]

    col_F = ( 0.1, 0.1, 0.1, 0.0 )
    pin.text_dct['FA'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FB'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FC'] = [True, col_F, 'phong 0.2 phong_size 5' ]


    # raytrace image/animation
    create_pov( pin, ['A', 'B', 'C'] )
    create_pov( pin, ['A', 'B', 'C', 'FA', 'FB', 'FC'] )
    create_pov( pin, ['A', 'B', 'C', 'FA', 'FB'] )

