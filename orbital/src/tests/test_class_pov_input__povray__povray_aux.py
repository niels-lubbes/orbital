'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''
import os

from orbital.sage_interface import sage_PolynomialRing
from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_var
from orbital.sage_interface import sage__eval
from orbital.sage_interface import sage_vector

from linear_series.class_linear_series import LinearSeries
from linear_series.class_base_points import BasePointTree
from linear_series.class_poly_ring import PolyRing

from orbital.surface_in_quadric import get_surf
from orbital.surface_in_quadric import approx_QQ
from orbital.surface_in_quadric import get_prj_mat
from orbital.surface_in_quadric import get_proj
from orbital.surface_in_quadric import get_S1xS1_pmz

from orbital.class_pov_input import PovInput

from orbital.povray import create_pov

from orbital.povray_aux import get_time_str


class TestPovray( object ):

    def test__povray( self ):

        # compute parametrizations of canonical model
        ring = PolyRing( 'x,y,v,w' )
        ring.ext_num_field( 't^2 + 1' )
        a0 = ring.root_gens()[0]
        bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
        bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )
        bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )
        ls_AB = LinearSeries.get( 2, bp_tree )
        ls_CB = LinearSeries.get( 1, bp_tree )

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


        # if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = './'

        # set PovInput as container
        # put very low quality for testing purposes
        pin = PovInput()
        pin.impl = f_xyz
        pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
        pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
        pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
        pin.path = os.environ['OUTPUT_PATH'] + get_time_str() + '_TEST_POVRAY_REMOVE_ME/'
        pin.fname = 'orb'
        pin.scale = 1
        pin.cam_dct['location'] = ( 0, 0, -10 )
        pin.cam_dct['lookat'] = ( 0, 0, 0 )
        pin.cam_dct['rotate'] = ( 310, 0, 0 )
        pin.light_radius = 5
        pin.axes_dct['show'] = False
        pin.axes_dct['len'] = 1.2
        pin.width = 5
        pin.height = 5
        pin.quality = 1
        pin.ani_delay = 1
        pin.curve_dct['A'] = {'step0':120, 'step1': 120, 'prec':1, 'width':0.02}
        pin.curve_dct['B'] = {'step0':120, 'step1': 360, 'prec':1, 'width':0.02}
        pin.curve_dct['C'] = {'step0':2 * 36, 'step1':2 * 36, 'prec':5, 'width':0.02}

        # raytrace image/animation
        show_surf = False
        ani = False
        ft_lst = []

        lst = create_pov( pin, ['A'], show_surf, ani, ft_lst )
        lst = create_pov( pin, ['B'], show_surf, True, ft_lst )

        print( pmz_AB_lst )
        print( pmz_CB_lst )

if __name__ == '__main__':

    TestPovray().test__povray()
    pass
