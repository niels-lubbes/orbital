'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes

'''

import os
import sys

from sage_interface import sage_QQ
from sage_interface import sage_NumberField
from sage_interface import sage_FractionField
from sage_interface import sage_PolynomialRing
from sage_interface import sage__eval
from sage_interface import sage_ideal
from sage_interface import sage_factor
from sage_interface import sage_var
from sage_interface import sage_set_verbose
from sage_interface import sage_invariant_theory
from sage_interface import sage_matrix
from sage_interface import sage_vector
from sage_interface import sage_identity_matrix
from sage_interface import sage_n
from sage_interface import sage_factor

from class_orb_tools import OrbTools

from poly_maps import image_map
from poly_maps import invert_birational_map
from poly_maps import compose_maps
from poly_maps import euclidean_type_form
from poly_maps import preimage_map
from poly_maps import ring_dict
from poly_maps import hilbert_poly

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

from surface_in_quadric import get_S1xS1_pmz
from surface_in_quadric import get_prj_mat
from surface_in_quadric import approx_QQ
from surface_in_quadric import get_surf
from surface_in_quadric import get_proj

from class_pov_input import PovInput
from povray import create_pov
from povray_aux import get_time_str


def usecase__two_sphere_cyclide():
    '''
    Construct parametrization of families of circles
    on the sphere cyclide. The sphere cyclide is 
    a quartic del Pezzo surface embedded into
    projective 3-sphere. This surface has two families
    of circles, but is not rational. 
    '''

    R = sage_PolynomialRing( sage_QQ, 'a' )
    a = R.gens()[0]
    B = sage_NumberField( [a ** 2 - 1 / sage_QQ( 2 )], 'a0' )
    B = sage_FractionField( sage_PolynomialRing( B, 'k' ) )
    OrbTools.p( ring_dict( B ) )

    S = '[-y0^2+y1^2+y2^2+y3^2+y4^2]'
    X = '[(x1^2+x2^2+x3^2)^2-x0^2*(8*x1^2-2*x2^2-2*x3^2-2*x0^2)]'
    OrbTools.p( euclidean_type_form( X[1:-1] ) )

    p = '[y0-y4,y1,y2,y3]'
    m = '[y4,y1+a0*y0]'
    P = '[x0-k*x1]'

    f = invert_birational_map( p, S )
    OrbTools.p( 'f =', f )

    Y = image_map( f, X )
    OrbTools.p( 'Y =', len( Y ), Y )
    OrbTools.p( 'Y =', hilbert_poly( Y ) )

    F = preimage_map( m, Y, P, B )
    OrbTools.p( 'F =', len( F ), F )

    C = image_map( p, F, B )
    OrbTools.p( 'C =', len( C ), C )

    for i in range( len( C ) ):
        s = str( C[i] )
        for ( a, b ) in [( 'a0', 'sqrt(1/2)' ), ( 'k', 'p1' ), ( 'x0', '1' ), ( 'x1', 'x' ), ( 'x2', 'y' ), ( 'x3', 'z' )]:
            s = s.replace( a, b )
        print( s )


def usecase__blum_cyclide():

    # construct dct
    a0 = PolyRing( 'x,y,v,w' ).ext_num_field( 't^2 + 1' ).root_gens()[0]  # i

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
    OrbTools.p( 'linear series 22 =\n', ls_22 )

    ls_21 = LinearSeries.get( [2, 1], bpt_1234 )
    OrbTools.p( 'linear series 21 =\n', ls_21 )

    ls_12 = LinearSeries.get( [1, 2], bpt_1234 )
    OrbTools.p( 'linear series 12 =\n', ls_12 )

    ls_11a = LinearSeries.get( [1, 1], bpt_12 )
    OrbTools.p( 'linear series 11a =\n', ls_11a )

    ls_11b = LinearSeries.get( [1, 1], bpt_34 )
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
    print( '---' )
    print( approx_QQ( U ) )
    print( '---' )
    print( approx_QQ( J ) )
    print( '---' )
    assert dct['M'] == approx_QQ( U.T * J * U )
    approxU = approx_QQ( U )
    P = sage_identity_matrix( 5 ).submatrix( 0, 0, 4, 5 )
    P[0, 4] = -1;
    print( P )
    P = P * approxU

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
    pin.cam_dct['location'] = ( 0, 0, -7 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 0, 40, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10

    pin.impl = f_xyz

    v0_lstE = [1.8, 2.3, 2.7, 3.1, 3.5, 3.8, 4.134, 4.31, 4.532, 4.7, 4.9, 5.08, 5.25, 5.405, 5.553, 5.7, 5.84]
    v0_lstF = [1.69, 1.87, 2.07, 2.26, 2.5, 2.72, 2.96, 3.2, 3.42, 3.65, 3.81]
    v0_lstD = [5.01, 5.12, 5.22, 5.32, 5.44, 5.56, 5.68, 5.81, 5.95, 6.1, 6.27, 6.474]

    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['D'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['E'] = ( pmz_EB_lst, 0 )
    pin.pmz_dct['F'] = ( pmz_AF_lst, 1 )

    pin.curve_dct['A'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.08}
    pin.curve_dct['B'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.08}
    pin.curve_dct['C'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.08}
    pin.curve_dct['D'] = {'step0':v0_lstD, 'step1':15, 'prec':10, 'width':0.08}
    pin.curve_dct['E'] = {'step0':v0_lstE, 'step1':15, 'prec':10, 'width':0.08}
    pin.curve_dct['F'] = {'step0':v0_lstF, 'step1':15, 'prec':10, 'width':0.08}

    pin.text_dct['A'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['C'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['E'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['F'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['D'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]


    # raytrace image/animation
    create_pov( pin, ['A'] )
    return
    create_pov( pin, ['A', 'B', 'C', 'D', 'E', 'F'] )


def usecase__ring_cyclide():

    # cos(a) = (1-m^2) / (1+m^2)
    # sin(a) = 2*m / (1+m^2)
    # m = arctan( a/2 )
    #
    x, y, v, w, c0, s0, c1, s1, R, r = sage_var( 'x,y,v,w,c0,s0,c1,s1,R,r' )
    V = sage_vector( [r * c0 + R, 0, r * s0] )
    M = sage_matrix( [( c1, -s1, 0 ), ( s1, c1, 0 ), ( 0, 0, 1 )] )
    print( M * V )
    C0 = ( y ** 2 - x ** 2 ) / ( y ** 2 + x ** 2 )
    S0 = 2 * x * y / ( y ** 2 + x ** 2 )
    C1 = ( w ** 2 - v ** 2 ) / ( w ** 2 + v ** 2 )
    S1 = 2 * v * w / ( w ** 2 + v ** 2 )
    den = ( y ** 2 + x ** 2 ) * ( w ** 2 + v ** 2 )

    dct = {c0:C0, s0:S0, c1:C1, s1:S1 }
    pmz_lst = [ ( elt.subs( dct ) * den ).simplify_full() for elt in list( M * V ) ]
    print( '---' )
    for pmz in pmz_lst:
        print( sage_factor( pmz ) )
    print( '---' )

    pmz_lst = [ pmz.subs( {r:1, R:2} ) for pmz in pmz_lst ]
    print( '---' )
    for pmz in pmz_lst:
        print( sage_factor( pmz ) )
    print( '---' )

    ls = LinearSeries( [str( pmz ) for pmz in pmz_lst], PolyRing( 'x,y,v,w' ) )
    OrbTools.p( ls.get_bp_tree() )
    return

    # construct dct
    a0, a1 = PolyRing( 'x,y,v,w' ).ext_num_field( 't^2+1/3' ).ext_num_field( 't^2+1' ).root_gens()

    p1 = [ 'xv', ( a0, 0 ) ]
    p2 = [ 'xv', ( -a0, 0 ) ]
    p3 = [ 'xw', ( a0, 0 ) ]
    p4 = [ 'xw', ( -a0, 0 ) ]

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

    sig = ( 1, 4 )  # ring cyclide
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
    coef_lst = [2, 1]
    dct = get_surf( ls_22, sig, coef_lst )

    # construct projection matrix P
    U, J = dct['UJ']
    U.swap_rows( 0, 2 )
    J.swap_columns( 0, 2 )
    J.swap_rows( 0, 2 )
    print( '---' )
    print( approx_QQ( sage_n( U ) ) )
    print( '---' )
    print( approx_QQ( J ) )
    print( '---' )
    assert dct['M'] == approx_QQ( U.T * J * U )
    approxU = approx_QQ( sage_n( U ) )
    P = sage_identity_matrix( 5 ).submatrix( 0, 0, 4, 5 )
    P[0, 4] = -1;
    print( P )
    P = P * approxU

    return
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

    for pmz, AB in pmz_lst:
        s = 'pmz' + AB + '=' + str( pmz )
        s = s.replace( '[', '{' ).replace( ']', '}' )
        print( s )

    # PovInput ring cyclide
    #
    pin = PovInput()

    pin.path = './' + get_time_str() + '_ring_cyclide/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, 0, sage_QQ( -21 ) / 10 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 310, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10

    pin.impl = f_xyz

    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
    pin.pmz_dct['D'] = ( pmz_DB_lst, 0 )
    pin.pmz_dct['E'] = ( pmz_EB_lst, 0 )
    pin.pmz_dct['F'] = ( pmz_AF_lst, 1 )

    pin.curve_dct['A'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['C'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['D'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['E'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['F'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}

    # raytrace image/animation
    create_pov( pin, ['A', 'B', 'C', 'D', 'E', 'F'] )


if __name__ == '__main__':

    #  Debug output settings
    #
    sage_set_verbose( -1 )
    OrbTools.filter( '__main__.py' )  # only print if verbose output by module <file_name>
    OrbTools.filter( None )  # print all verbose output
    OrbTools.start_timer()

    if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = '/home/niels/Desktop/n/src/output/orb/povray/'

    #########################################
    #                                       #
    # (Un)comment one or more use cases     #
    #                                       #
    #########################################

    # usecase__two_sphere_cyclide()
    usecase__blum_cyclide()
    # usecase__ring_cyclide()


    #########################################
    #                                       #
    # End of list of use case methods.      #
    #                                       #
    #########################################

    # end timing
    OrbTools.end_timer()

    print( '\nThe End' )



