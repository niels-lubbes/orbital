'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 17, 2018
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
from sage_interface import sage_cos
from sage_interface import sage_sin

from class_orb_tools import OrbTools

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

from class_orb_ring import OrbRing

from class_pov_input import PovInput
from povray import create_pov
from povray_aux import get_time_str

import sys



def perseus_cyclide():

    # Construct linear series for Perseus cyclide.
    #
    # We first construct a trigonometric parametrization
    # by rotating a circle.
    #
    # We convert via the following formulas, the trigonometric
    # parametrization to a rational parametrization.
    #
    #     cos(s) = (y^2-x^2) / (y^2+x^2)
    #     sin(s) = 2*x*y / (y^2+x^2)
    #     y=1; x = arctan( s/2 )
    #
    #     cos(pi/3)=1/2
    #     sin(pi/3)=sqrt(3)/2
    #
    # After that we do a basepoint analysis on the obtained rational
    # parametrization.
    #
    r, R = 1, 2
    c0, s0, c1, s1 = sage_var( 'c0,s0,c1,s1' )
    x, y, v, w, a0 = sage_var( 'x,y,v,w,a0' )
    q2 = sage_QQ( 1 ) / 2

    MZ = sage_matrix( [( c1, s1, 0 ), ( -s1, c1, 0 ), ( 0, 0, 1 )] )
    MZc = MZ.subs( {c1:q2, s1:q2 * a0} )

    V = sage_vector( [r * c0, 0, r * s0] )
    V = MZc * V
    V[0] = V[0] + R
    OrbTools.p( 'V =', V )

    pmz_AB_lst = list( MZ * V )

    OrbTools.p( 'pmz_AB_lst =', pmz_AB_lst )
    for pmz in pmz_AB_lst:
        OrbTools.p( '\t\t', sage_factor( pmz ) )

    C0 = ( y ** 2 - x ** 2 ) / ( y ** 2 + x ** 2 )
    S0 = 2 * x * y / ( y ** 2 + x ** 2 )
    C1 = ( w ** 2 - v ** 2 ) / ( w ** 2 + v ** 2 )
    S1 = 2 * v * w / ( w ** 2 + v ** 2 )
    den = ( y ** 2 + x ** 2 ) * ( w ** 2 + v ** 2 )
    dct = {c0:C0, s0:S0, c1:C1, s1:S1 }

    pmz_lst = [den] + [ ( elt.subs( dct ) * den ).simplify_full() for elt in list( MZ * V ) ]
    OrbTools.p( 'pmz_lst =', pmz_lst )
    for pmz in pmz_lst:
        OrbTools.p( '\t\t', sage_factor( pmz ) )

    # do a basepoint analysis on the rational parametrization
    # The True argument is for resetting the number field to QQ!
    ring = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2-3' )
    ls = LinearSeries( [str( pmz ) for pmz in pmz_lst], ring )
    OrbTools.p( ls.get_bp_tree() )

    # construct linear series for families of conics
    ring = PolyRing( 'x,y,v,w' )  # construct polynomial ring over new ground field
    OrbTools.p( ring )
    x, y, v, w = ring.gens()
    a0, a1, a2, a3 = ring.root_gens()

    p1 = [ 'xv', ( -a3, a1 ) ]
    p2 = [ 'xv', ( -a2, -a1 ) ]
    p3 = [ 'xv', ( a3, a1 ) ]
    p4 = [ 'xv', ( a2, -a1 ) ]

    bpt_1234 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bpt_1234.add( p1[0], p1[1], 1 )
    bpt_1234.add( p2[0], p2[1], 1 )
    bpt_1234.add( p3[0], p3[1], 1 )
    bpt_1234.add( p4[0], p4[1], 1 )

    bpt_12 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bpt_12.add( p1[0], p1[1], 1 )
    bpt_12.add( p2[0], p2[1], 1 )

    bpt_34 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bpt_34.add( p3[0], p3[1], 1 )
    bpt_34.add( p4[0], p4[1], 1 )

    if True:
        ls_22 = LinearSeries.get( [2, 2], bpt_1234 )  # |2(l1+l2)-e1-e2-e3-e4|
        ls_21 = LinearSeries.get( [2, 1], bpt_1234 )
        ls_12 = LinearSeries.get( [1, 2], bpt_1234 )
        ls_11a = LinearSeries.get( [1, 1], bpt_12 )
        ls_11b = LinearSeries.get( [1, 1], bpt_34 )

        OrbTools.p( 'linear series 22 =\n', ls_22 )
        OrbTools.p( 'linear series 21 =\n', ls_21 )
        OrbTools.p( 'linear series 12 =\n', ls_12 )
        OrbTools.p( 'linear series 11a =\n', ls_11a )
        OrbTools.p( 'linear series 11b =\n', ls_11b )


    # compute reparametrization from the linear series of families
    ring = PolyRing( 'x,y,v,w,c0,s0,c1,s1' )  # construct polynomial ring with new generators
    OrbTools.p( ring )
    x, y, v, w, c0, s0, c1, s1 = ring.gens()
    a0, a1, a2, a3 = ring.root_gens()
    pmz_AB_lst = [1] + sage__eval( str( pmz_AB_lst ), ring.ring_dct )
    pmz_lst = sage__eval( str( pmz_lst ), ring.ring_dct )

    q2 = sage_QQ( 1 ) / 2
    a = 2 * a0 / 3
    b = ( -a0 * a1 / 3 - q2 ) * a3
    c = ( a0 * a1 / 3 - q2 ) * a2
    d = ( a1 / 2 - a0 / 3 ) * a3
    e = ( -a1 / 2 - a0 / 3 ) * a2
    bc = b + c
    de = d + e

    X = 1 - s0; Y = c0;
    V = 1 - s1; W = c1;
    CB_dct = {x:X, y:Y, v:W * X + bc * W * Y - de * V * Y, w:V * X + bc * V * Y + de * W * Y};
    DB_dct = {x:X, y:Y, v:W * X - bc * W * Y + de * V * Y, w:V * X - bc * V * Y - de * W * Y};
    EB_dct = {x:X, y:Y, v:W * X ** 2 + W * Y ** 2 - a * V * Y ** 2, w:V * X ** 2 + V * Y ** 2 + a * W * Y ** 2};
    pmz_CB_lst = [ pmz.subs( CB_dct ) for pmz in pmz_lst ]  # CB  11a
    pmz_DB_lst = [ pmz.subs( DB_dct ) for pmz in pmz_lst ]  # CB  11b
    pmz_EB_lst = [ pmz.subs( EB_dct ) for pmz in pmz_lst ]  # CB  21

    # output
    OrbTools.p( 'pmz_AB_lst =\n', pmz_AB_lst )
    OrbTools.p( 'pmz_CB_lst =\n', pmz_CB_lst )
    OrbTools.p( 'pmz_DB_lst =\n', pmz_DB_lst )
    OrbTools.p( 'pmz_EB_lst =\n', pmz_EB_lst )

    # approximate by map defined over rational numbers
    ci_idx = 5  # index defining the complex embedding
    pmz_AB_lst = OrbRing.approx_QQ_pol_lst( pmz_AB_lst, ci_idx )
    pmz_CB_lst = OrbRing.approx_QQ_pol_lst( pmz_CB_lst, ci_idx )
    pmz_DB_lst = OrbRing.approx_QQ_pol_lst( pmz_DB_lst, ci_idx )
    pmz_EB_lst = OrbRing.approx_QQ_pol_lst( pmz_EB_lst, ci_idx )

    # mathematica input
    pmz_lst = [ ( pmz_lst, 'ZZ' ),
                ( pmz_AB_lst, 'AB' ),
                ( pmz_CB_lst, 'CB' ),
                ( pmz_DB_lst, 'DB' ),
                ( pmz_EB_lst, 'EB' )
                ]
    ms = ''
    for pmz, AB in pmz_lst:
        s = 'pmz' + AB + '=' + str( pmz ) + ';'
        s = s.replace( '[', '{' ).replace( ']', '}' )
        ms += '\n' + s
    OrbTools.p( 'Mathematica input =', ms )

    # PovInput ring cyclide
    #
    pin = PovInput()

    pin.path = './' + get_time_str() + '_perseus_cyclide/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, -7, 0 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 45, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 400
    pin.height = 200
    pin.quality = 1
    pin.ani_delay = 10

    pin.impl = None

    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['D'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['E'] = ( pmz_EB_lst, 0 )

    pin.pmz_dct['FA'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['FB'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['FC'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['FD'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['FE'] = ( pmz_EB_lst, 0 )

    v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 10 )]  # 10
    v1_lst_A = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 15 )]  # 15
    v1_lst_B = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 36 )]  # 15
    v1_lst_C = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 36 )]  # 15
    v1_lst_D = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 36 )]  # 15
    v1_lst_E = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 15 )]  # 15


    v1_lst_F = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 2 )]

    prec = 50

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst_A, 'prec':prec, 'width':0.05}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst_B, 'prec':prec, 'width':0.05}
    pin.curve_dct['C'] = {'step0':v0_lst, 'step1':v1_lst_C, 'prec':prec, 'width':0.05}
    pin.curve_dct['D'] = {'step0':v0_lst, 'step1':v1_lst_D, 'prec':prec, 'width':0.05}
    pin.curve_dct['E'] = {'step0':v0_lst, 'step1':v1_lst_E, 'prec':prec, 'width':0.05}

    pin.curve_dct['FA'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}
    pin.curve_dct['FB'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}
    pin.curve_dct['FC'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}
    pin.curve_dct['FD'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}
    pin.curve_dct['FE'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}


    pin.text_dct['A'] = [True, ( 0.4, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.0, 0.0, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['C'] = [True, ( 0.0, 0.2, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['D'] = [True, ( 0.2, 0.0, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['E'] = [True, ( 0.1, 0.3, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]

    col_F = ( 0.1, 0.1, 0.1, 0.0 )
    pin.text_dct['FA'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FB'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FC'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FD'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FE'] = [True, col_F, 'phong 0.2 phong_size 5' ]


    # raytrace image/animation
    create_pov( pin, ['A', 'B', 'C', 'D', 'E', 'FA', 'FB', 'FC', 'FD', 'FE'] )
    create_pov( pin, ['C', 'D', 'FC', 'FD'] )
    create_pov( pin, ['A', 'E', 'FC', 'FD'] )
    create_pov( pin, ['A', 'B', 'E', 'FC', 'FD'] )