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


def ring_cyclide():

    # cos(a) = (1-m^2) / (1+m^2)
    # sin(a) = 2*m / (1+m^2)
    # m = arctan( a/2 )
    #
    x, y, v, w, c0, s0, c1, s1, R, r = sage_var( 'x,y,v,w,c0,s0,c1,s1,R,r' )
    V = sage_vector( [r * c0 + R, 0, r * s0] )
    M = sage_matrix( [( c1, -s1, 0 ), ( s1, c1, 0 ), ( 0, 0, 1 )] )
    pmz_AB_Rr_lst = [1] + list( M * V )
    OrbTools.p( 'pmz_AB_Rr_lst =', pmz_AB_Rr_lst )
    for pmz in pmz_AB_Rr_lst:
        OrbTools.p( '\t\t', sage_factor( pmz ) )

    pmz_AB_lst = [ pmz.subs( {r:1, R:2} ) if type( pmz ) != int else pmz for pmz in pmz_AB_Rr_lst ]
    OrbTools.p( 'pmz_AB_lst =', pmz_AB_lst )
    for pmz in pmz_AB_lst:
        OrbTools.p( '\t\t', sage_factor( pmz ) )

    C0 = ( y ** 2 - x ** 2 ) / ( y ** 2 + x ** 2 )
    S0 = 2 * x * y / ( y ** 2 + x ** 2 )
    C1 = ( w ** 2 - v ** 2 ) / ( w ** 2 + v ** 2 )
    S1 = 2 * v * w / ( w ** 2 + v ** 2 )
    den = ( y ** 2 + x ** 2 ) * ( w ** 2 + v ** 2 )
    dct = {c0:C0, s0:S0, c1:C1, s1:S1 }

    pmz_lst = [den] + [ ( elt.subs( dct ) * den ).simplify_full() for elt in list( M * V ) ]
    OrbTools.p( 'pmz_lst =', pmz_lst )

    pmz_rR_lst = [ pmz.subs( {r:1, R:2} ) for pmz in pmz_lst ]
    OrbTools.p( 'pmz_rR_lst =', pmz_rR_lst )
    for pmz in pmz_rR_lst:
        OrbTools.p( '\t\t', sage_factor( pmz ) )

    ls = LinearSeries( [str( pmz ) for pmz in pmz_rR_lst], PolyRing( 'x,y,v,w' ) )
    OrbTools.p( ls.get_bp_tree() )

    # construct dct
    a0, a1 = PolyRing( 'x,y,v,w' ).ext_num_field( 't^2+1/3' ).ext_num_field( 't^2+1' ).root_gens()

    p1 = [ 'xv', ( -a0, a1 ) ]
    p2 = [ 'xv', ( a0, -a1 ) ]
    p3 = [ 'xv', ( -a0, -a1 ) ]
    p4 = [ 'xv', ( a0, a1 ) ]

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

    ls_22 = LinearSeries.get( [2, 2], bpt_1234 )  # |2(l1+l2)-e1-e2-e3-e4|
    OrbTools.p( 'linear series 22 =\n', ls_22 )

    ls_21 = LinearSeries.get( [2, 1], bpt_1234 )
    OrbTools.p( 'linear series 21 =\n', ls_21 )

    ls_12 = LinearSeries.get( [1, 2], bpt_1234 )
    OrbTools.p( 'linear series 12 =\n', ls_12 )

    ls_11a = LinearSeries.get( [1, 1], bpt_12 )
    OrbTools.p( 'linear series 11a =\n', ls_11a )

    ls_11b = LinearSeries.get( [1, 1], bpt_34 )
    OrbTools.p( 'linear series 11b =\n', ls_11b )

    # compute reparametrization
    R_xyvw = sage_PolynomialRing( sage_QQ, 'x,y,v,w' )
    x, y, v, w = R_xyvw.gens()
    X, Y, V, W, q = sage_var( 'x,y,v,w,q' )
    c0, s0, c1, s1 = sage_var( 'c0,s0,c1,s1' )
    xyvw_dct = { X:x, Y:y, V:v, W:w }
    trig_dct = {X:1 - s0, Y:c0, V:1 - s1, W:c1}
    pol_lst = sage__eval( str( pmz_rR_lst ), R_xyvw.gens_dict() )
    q = sage_n( sage_sqrt( 3 ) ).exact_rational()  # approximation of sqrt(3)

    # CB
    CB_dct = { x:X, y:Y, v: W * X + q * V * Y, w: V * X - q * W * Y }
    pmz_CB_lst = [ pol.subs( CB_dct ).subs( trig_dct ) for pol in pol_lst ]

    # DB
    DB_dct = { x:X, y:Y, v: W * X - q * V * Y, w: V * X + q * W * Y  }
    pmz_DB_lst = [ pol.subs( DB_dct ).subs( trig_dct ) for pol in pol_lst ]

    # output
    OrbTools.p( 'pmz_AB_lst =\n', pmz_AB_lst )
    OrbTools.p( 'pmz_CB_lst =\n', pmz_CB_lst )
    OrbTools.p( 'pmz_DB_lst =\n', pmz_DB_lst )

    # mathematica
    pmz_lst = [ ( pmz_AB_lst, 'AB' ),
                ( pmz_CB_lst, 'CB' ),
                ( pmz_DB_lst, 'DB' )]

    for pmz, AB in pmz_lst:
        s = 'pmz' + AB + '=' + str( pmz ) + ';'
        s = s.replace( '[', '{' ).replace( ']', '}' )
        print( s )

    # PovInput ring cyclide
    #
    pin = PovInput()

    pin.path = './' + get_time_str() + '_ring_cyclide/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, -7, 0 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 45, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10

    pin.impl = None

    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['D'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['FA'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['FB'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['FC'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['FD'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['WA'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['WB'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['WD'] = ( pmz_DB_lst, 0 )

    v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 10 )]
    v1_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 15 )]

    v1_lst_A = [ sage_pi / 2 + ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 30 )]
    v1_lst_F = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 1 )]

    v1_lst_WA = [0.1, 0.52, 0.94, 1.36, 1.78, 2.2, 2.61, 3.04, 3.45, 3.88, 4.3, 4.712, 5.13, 5.55, 5.965]
    v1_lst_WB = [0, 0.7, 1.31, 1.8, 2.18, 2.5, 2.77, 3.015, 3.26, 3.51, 3.78, 4.099, 4.49, 4.97, 5.579];
    v1_lst_WD = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 15 )]

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst_A, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['C'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['D'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['FA'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.02}
    pin.curve_dct['FB'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.02}
    pin.curve_dct['FC'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.02}
    pin.curve_dct['FD'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.02}
    pin.curve_dct['WA'] = {'step0':v0_lst, 'step1':v1_lst_WA, 'prec':10, 'width':0.05}
    pin.curve_dct['WB'] = {'step0':v0_lst, 'step1':v1_lst_WB, 'prec':10, 'width':0.05}
    pin.curve_dct['WD'] = {'step0':v0_lst, 'step1':v1_lst_WD, 'prec':10, 'width':0.05}

    pin.text_dct['A'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['C'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['D'] = [True, ( 0.0, 0.2, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FA'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FB'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FC'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FD'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['WA'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['WB'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['WD'] = [True, ( 0.0, 0.2, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]

    # raytrace image/animation
    create_pov( pin, ['A', 'B', 'C', 'D', 'FA', 'FB', 'FC', 'FD', 'WA', 'WB', 'WD'] )

    F_lst = ['FA', 'FC', 'FD']
    create_pov( pin, ['A', 'C', 'D'] )
    create_pov( pin, ['A', 'C', 'D'] + F_lst )
    create_pov( pin, ['WA', 'WB', 'WD'] )
    create_pov( pin, ['WA', 'WB', 'WD'] + F_lst )


