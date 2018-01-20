'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 16, 2018
@author: Niels Lubbes
'''

from sage_interface import sage_PolynomialRing
from sage_interface import sage_QQ
from sage_interface import sage_var
from sage_interface import sage__eval
from sage_interface import sage_vector
from sage_interface import sage_pi

from linear_series.class_linear_series import LinearSeries
from linear_series.class_base_points import BasePointTree
from linear_series.class_poly_ring import PolyRing

from surface_in_quadric import get_surf
from surface_in_quadric import approx_QQ
from surface_in_quadric import get_prj_mat
from surface_in_quadric import get_proj
from surface_in_quadric import get_S1xS1_pmz

from class_pov_input import PovInput

from povray import create_pov

from povray_aux import get_time_str



def dp6_smooth():

    # compute parametrizations of canonical model
    a0 = PolyRing( 'x,y,v,w' ).ext_num_field( 't^2 + 1' ).root_gens()[0]
    bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )
    bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )
    ls_AB = LinearSeries.get( [2, 2], bp_tree )
    ls_CB = LinearSeries.get( [1, 1], bp_tree )

    # compute surface in quadric of signature (6,1)
    c_lst = [-1, -1, 0, 0, 0, -1, 1, 0, -1, -1, -1]
    dct = get_surf( ls_AB, ( 6, 1 ), c_lst )

    # compute projection to R^3
    approxU = approx_QQ( dct['UJ'][0] )
    P = get_prj_mat( 4, 7, 0 )
    P[0, 6] = -1;P[3, 3] = 0;P[3, 4] = 1
    P = P * approxU
    f_xyz, pmz_AB_lst = get_proj( dct['imp_lst'], dct['pmz_lst'], P )

    # compute reparametrization
    R_xyvw = sage_PolynomialRing( sage_QQ, 'x,y,v,w' )
    x, y, v, w = R_xyvw.gens()
    X, Y, V, W = sage_var( 'x,y,v,w' )
    xyvw_dct = { X:x, Y:y, V:v, W:w }
    XYZW_dct = { x:X, y:Y, v:V, w:W }
    CB_dct = { x:X, y:Y, v:X * W + Y * V, w: X * V - Y * W }
    pol_lst = sage__eval( str( ls_AB.pol_lst ), R_xyvw.gens_dict() )
    pmz_CB_lst = [ pol.subs( CB_dct ).subs( xyvw_dct ) for pol in pol_lst]
    pmz_CB_lst = get_S1xS1_pmz( pmz_CB_lst )
    pmz_CB_lst = list( P * dct['Q'] * sage_vector( pmz_CB_lst ) )

    # set PovInput as container
    # put very low quality for testing purposes
    pin = PovInput()

    pin.path = './' + get_time_str() + '_dp6_smooth/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, 0, sage_QQ( -21 ) / 10 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 310, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 400
    pin.height = 800
    pin.quality = 11
    pin.ani_delay = 1

    pin.impl = f_xyz

    v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 10 )]
    v1_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 15 )]

    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}
    pin.curve_dct['C'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.05}

    # raytrace image/animation
    create_pov( pin, ['A', 'B'] )


