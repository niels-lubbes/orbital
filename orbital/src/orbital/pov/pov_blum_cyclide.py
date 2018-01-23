'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 16, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_PolynomialRing
from orbital.sage_interface import sage__eval
from orbital.sage_interface import sage_var
from orbital.sage_interface import sage_identity_matrix
from orbital.sage_interface import sage_invariant_theory
from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_pi

from orbital.class_orb_tools import OrbTools

from orbital.surface_in_quadric import get_S1xS1_pmz
from orbital.surface_in_quadric import approx_QQ
from orbital.surface_in_quadric import get_surf
from orbital.surface_in_quadric import get_proj

from orbital.povray.class_pov_input import PovInput
from orbital.povray.povray import create_pov
from orbital.povray.povray_aux import get_time_str

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries


def blum_cyclide():

    # construct dct
    a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]  # i

    bpt_1234 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bpt_1234.add( 'xv', ( -1 * a0, 1 * a0 ), 1 )  # e1
    bpt_1234.add( 'xv', ( 1 * a0, -1 * a0 ), 1 )  # e2
    bpt_1234.add( 'xw', ( -2 * a0, 2 * a0 ), 1 )  # e3
    bpt_1234.add( 'xw', ( 2 * a0, -2 * a0 ), 1 )  # e4

    bpt_12 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bpt_12.add( 'xv', ( -1 * a0, 1 * a0 ), 1 )  # e1
    bpt_12.add( 'xv', ( 1 * a0, -1 * a0 ), 1 )  # e2

    bpt_34 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bpt_34.add( 'xw', ( -2 * a0, 2 * a0 ), 1 )  # e3
    bpt_34.add( 'xw', ( 2 * a0, -2 * a0 ), 1 )  # e4

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

    sig = ( 4, 1 )
    pol_lst = ls_22.get_implicit_image()

    # determine signature
    x_lst = sage_PolynomialRing( sage_QQ, [ 'x' + str( i ) for i in range( sum( sig ) )] ).gens()
    for pol in pol_lst:

        if pol.degree() == 2:
            M = sage_invariant_theory.quadratic_form( pol, x_lst ).as_QuadraticForm().matrix()
            D, V = sage_matrix( sage_QQ, M ).eigenmatrix_right()  # D has first all negative values on diagonal
            cur_sig = ( len( [ d for d in D.diagonal() if d < 0 ] ), len( [ d for d in D.diagonal() if d > 0 ] ) )
        else:
            cur_sig = '[no signature]'
        OrbTools.p( '\t\t', pol, cur_sig )

    # obtain surface in sphere
    coef_lst = [0, -1, -1]
    dct = get_surf( ls_22, sig, coef_lst )

    # construct projection matrix P
    U, J = dct['UJ']
    U.swap_rows( 0, 4 )
    J.swap_columns( 0, 4 )
    J.swap_rows( 0, 4 )
    assert dct['M'] == approx_QQ( U.T * J * U )
    approxU = approx_QQ( U )
    P = sage_identity_matrix( 5 ).submatrix( 0, 0, 4, 5 )
    P[0, 4] = -1;
    P = P * approxU
    OrbTools.p( ' approx_QQ( U ) =', list( approx_QQ( U ) ) )
    OrbTools.p( ' approx_QQ( J ) =', list( approx_QQ( J ) ) )
    OrbTools.p( ' P              =', list( P ) )

    # call get_proj
    f_xyz, pmz_AB_lst = get_proj( dct['imp_lst'], dct['pmz_lst'], P )
    f_xyz_deg_lst = [f_xyz.degree( sage_var( v ) ) for v in ['x', 'y', 'z']]

    # compute reparametrization
    R_xyvw = sage_PolynomialRing( sage_QQ, 'x,y,v,w' )
    x, y, v, w = R_xyvw.gens()
    X, Y, V, W = sage_var( 'x,y,v,w' )
    xyvw_dct = { X:x, Y:y, V:v, W:w }
    pol_lst = sage__eval( str( ls_22.pol_lst ), R_xyvw.gens_dict() )

    # CB
    CB_dct = { x:X, y:Y, v:X * W + Y * V, w: X * V - Y * W }
    pmz_CB_lst = get_S1xS1_pmz( [ pol.subs( CB_dct ).subs( xyvw_dct ) for pol in pol_lst] )
    pmz_CB_lst = list( P * sage_vector( pmz_CB_lst ) )

    # DB
    DB_dct = { x:X, y:Y, v:4 * X * W - Y * V, w: X * V + Y * W }
    pmz_DB_lst = get_S1xS1_pmz( [ pol.subs( DB_dct ).subs( xyvw_dct ) for pol in pol_lst] )
    pmz_DB_lst = list( P * sage_vector( pmz_DB_lst ) )

    # EB
    EB_dct = { x:X, y:Y, v:40 * W * X ** 2 + 25 * W * Y ** 2 + 24 * V * X * Y, w:40 * V * X ** 2 + 16 * V * Y ** 2 - 15 * W * X * Y  }
    pmz_EB_lst = get_S1xS1_pmz( [ pol.subs( EB_dct ).subs( xyvw_dct ) for pol in pol_lst] )
    pmz_EB_lst = list( P * sage_vector( pmz_EB_lst ) )

    # AF
    AF_dct = { x:-10 * Y * V ** 2 - 25 * Y * W ** 2 + 9 * X * V * W,
               y:15 * X * V ** 2 + 24 * X * W ** 2 - 15 * Y * V * W,
               v:V, w:W  }
    pmz_AF_lst = get_S1xS1_pmz( [ pol.subs( AF_dct ).subs( xyvw_dct ) for pol in pol_lst] )
    pmz_AF_lst = list( P * sage_vector( pmz_AF_lst ) )

    # output
    OrbTools.p( 'f_xyz =', f_xyz_deg_lst, '\n', f_xyz )
    OrbTools.p( 'pmz_AB_lst =\n', pmz_AB_lst )
    OrbTools.p( 'pmz_CB_lst =\n', pmz_CB_lst )
    OrbTools.p( 'pmz_DB_lst =\n', pmz_DB_lst )
    OrbTools.p( 'pmz_EB_lst =\n', pmz_EB_lst )
    OrbTools.p( 'pmz_AF_lst =\n', pmz_AF_lst )

    # mathematica
    pmz_lst = [ ( pmz_AB_lst, 'AB' ),
                ( pmz_CB_lst, 'CB' ),
                ( pmz_DB_lst, 'DB' ),
                ( pmz_EB_lst, 'EB' ),
                ( pmz_AF_lst, 'AF' )]

    OrbTools.p( 'Mathematica input for ParametricPlot3D:' )
    for pmz, AB in pmz_lst:
        s = 'pmz' + AB + '=' + str( pmz )
        s = s.replace( '[', '{' ).replace( ']', '}' )
        print( s )

    # PovInput for Blum cyclide
    #
    pin = PovInput()
    pin.path = './' + get_time_str() + '_blum_cyclide/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, -7, 0 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 20, 180, 20 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10
    pin.impl = f_xyz

    start0 = sage_QQ( 1 ) / 10  # step0=10 step1=15
    v0_lst = [ start0 + ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 10 )]
    v1_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 15 )]
    v1_lst_F = [ start0 + ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 1 )]

    v1_lst_WE = [1.8, 2.3, 2.7, 3.1, 3.5, 3.8, 4.134, 4.31, 4.532, 4.7, 4.9, 5.08, 5.25, 5.405, 5.553, 5.7, 5.84]
    v1_lst_WF = [1.69, 1.87, 2.07, 2.26, 2.5, 2.72, 2.96, 3.2, 3.42, 3.65, 3.81]
    v1_lst_WD = [5.01, 5.12, 5.22, 5.32, 5.44, 5.56, 5.68, 5.81, 5.95, 6.1, 6.27, 6.474]

    v1_lst_SA = [6.5]; v1_lst_SE = [5.4];
    v1_lst_SB = [5.95]; v1_lst_SF = [2.28];
    v1_lst_SC = [4.83]; v1_lst_SD = [5.55];

    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['D'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['E'] = ( pmz_EB_lst, 0 )
    pin.pmz_dct['F'] = ( pmz_AF_lst, 1 )
    pin.pmz_dct['WD'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['WE'] = ( pmz_EB_lst, 0 )
    pin.pmz_dct['WF'] = ( pmz_AF_lst, 1 )
    pin.pmz_dct['SA'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['SB'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['SC'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['SD'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['SE'] = ( pmz_EB_lst, 0 )
    pin.pmz_dct['SF'] = ( pmz_AF_lst, 1 )
    pin.pmz_dct['FA'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['FB'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['FC'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['FD'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['FE'] = ( pmz_EB_lst, 0 )
    pin.pmz_dct['FF'] = ( pmz_AF_lst, 1 )

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['C'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['D'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['E'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['F'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}

    pin.curve_dct['WD'] = {'step0':v0_lst, 'step1':v1_lst_WD, 'prec':10, 'width':0.05}
    pin.curve_dct['WE'] = {'step0':v0_lst, 'step1':v1_lst_WE, 'prec':10, 'width':0.05}
    pin.curve_dct['WF'] = {'step0':v0_lst, 'step1':v1_lst_WF, 'prec':10, 'width':0.05}

    pin.curve_dct['SA'] = {'step0':v0_lst, 'step1':v1_lst_SA, 'prec':10, 'width':0.05}
    pin.curve_dct['SB'] = {'step0':v0_lst, 'step1':v1_lst_SB, 'prec':10, 'width':0.05}
    pin.curve_dct['SC'] = {'step0':v0_lst, 'step1':v1_lst_SC, 'prec':10, 'width':0.05}
    pin.curve_dct['SD'] = {'step0':v0_lst, 'step1':v1_lst_SD, 'prec':10, 'width':0.06}
    pin.curve_dct['SE'] = {'step0':v0_lst, 'step1':v1_lst_SE, 'prec':10, 'width':0.05}
    pin.curve_dct['SF'] = {'step0':v0_lst, 'step1':v1_lst_SF, 'prec':10, 'width':0.05}

    pin.curve_dct['FA'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}
    pin.curve_dct['FB'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}
    pin.curve_dct['FC'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}
    pin.curve_dct['FD'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}
    pin.curve_dct['FE'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}
    pin.curve_dct['FF'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':10, 'width':0.01}

    pin.text_dct['A'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['C'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['E'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['F'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['D'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['WE'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['WF'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['WD'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['SA'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['SB'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['SC'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['SE'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['SF'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['SD'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FA'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FB'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FC'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FE'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FF'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['FD'] = [True, ( 0.1, 0.1, 0.1, 0.0 ), 'phong 0.2 phong_size 5' ]

    # raytrace image/animation
    fam_lst = []
    fam_lst += ['A', 'B', 'C', 'D', 'E', 'F']
    fam_lst += ['WD', 'WE', 'WF']
    fam_lst += ['SA', 'SB', 'SC', 'SD', 'SE', 'SF']
    fam_lst += ['FA', 'FB', 'FC', 'FD', 'FE', 'FF']
    create_pov( pin, fam_lst )

    F_lst = ['FA', 'FB', 'FC']
    S_lst = ['SA', 'SB', 'SC', 'SD', 'SE', 'SF']
    create_pov( pin, ['A', 'B', 'C'] )
    create_pov( pin, ['A', 'B', 'C'] + F_lst )
    create_pov( pin, ['WD', 'WE', 'WF'] )
    create_pov( pin, ['WD', 'WE', 'WF'] + F_lst )
    create_pov( pin, S_lst + F_lst )

