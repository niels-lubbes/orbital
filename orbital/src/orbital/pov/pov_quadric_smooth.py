'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 16, 2018
@author: Niels Lubbes
'''

from sage_interface import sage_QQ
from sage_interface import sage_PolynomialRing
from sage_interface import sage__eval
from sage_interface import sage_var
from sage_interface import sage_matrix
from sage_interface import sage_vector
from sage_interface import sage_factor
from sage_interface import sage_n
from sage_interface import sage_factor
from sage_interface import sage_sqrt
from sage_interface import sage_pi

from class_orb_tools import OrbTools

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

from class_pov_input import PovInput
from povray import create_pov
from povray_aux import get_time_str

from class_orb_ring import OrbRing


def quadric_smooth():

    c0, s0, c1, s1, t0, t1 = OrbRing.coerce( 'c0,s0,c1,s1,t0,t1' )
    P0 = sage_vector( [-2, -1, -0.5] )
    Q0 = sage_vector( [2, -1, -0.5] )
    P1 = sage_vector( [2, -1, -0.5] )
    Q1 = sage_vector( [-2, -1, -0.5] )
    V0 = t0 * P0 + ( t0 - 1 ) * Q0;
    V1 = t0 * P1 + ( t0 - 1 ) * Q1;

    M0 = sage_matrix( [( c1, s1, 0 ), ( -s1, c1, 0 ), ( 0, 0, 1 )] )
    M1 = sage_matrix( [( c1, s1, 0 ), ( -s1, c1, 0 ), ( 0, 0, 1 )] )

    pmz_A_lst = [1] + list( M0 * V0 )
    pmz_B_lst = [1] + list( M1 * V1 )

    OrbTools.p( 'pmz_A_lst =', pmz_A_lst )
    for pmz in pmz_A_lst:
        OrbTools.p( '\t\t', pmz )

    OrbTools.p( 'pmz_B_lst =', pmz_B_lst )
    for pmz in pmz_B_lst:
        OrbTools.p( '\t\t', pmz )

    # PovInput ring cyclide
    #
    pin = PovInput()

    pin.path = './' + get_time_str() + '_quadric_smooth/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, -10, 0 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 0, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 400
    pin.height = 200
    pin.quality = 1
    pin.ani_delay = 10

    pin.impl = None

    pin.pmz_dct['A'] = ( pmz_A_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_B_lst, 0 )
    pin.pmz_dct['FA'] = ( pmz_A_lst, 0 )
    pin.pmz_dct['FB'] = ( pmz_B_lst, 0 )

    v0_lst = [sage_QQ( i ) / 10 for i in range( -15, 30, 5 )]  # -15, 35
    v1_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 36 )]
    v1_lst_F = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 2 )]

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.1}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.1}
    pin.curve_dct['FA'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}
    pin.curve_dct['FB'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}


    pin.text_dct['A'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FA'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FB'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]


    # raytrace image/animation
    create_pov( pin, ['A', 'B', 'FA', 'FB'] )


