'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 23, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage__eval
from orbital.sage_interface import sage_var
from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_factor
from orbital.sage_interface import sage_n
from orbital.sage_interface import sage_pi

from orbital.sage_interface import sage_gcd
from orbital.sage_interface import sage_Combinations

from orbital.class_orb_tools import OrbTools

from orbital.class_orb_ring import OrbRing

from orbital.povray.class_pov_input import PovInput
from orbital.povray.povray import create_pov
from orbital.povray.povray_aux import get_time_str

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

from orbital.poly_maps import compose_maps


def trace_bp_tree( bp_lst ):
    out_lst = []
    for bp in bp_lst:
        if bp.mult > 0:
            out_lst += [ bp.mult ]
        out_lst += trace_bp_tree( bp.bp_lst_t )
        out_lst += trace_bp_tree( bp.bp_lst_s )

    return out_lst


def veronese_blowup():
    '''
    Construct a povray image of two minimal families on the blowup of the 
    Veronese surface in 12 complex conjugate points.
    '''

    # construct linear series
    p1 = [ 'z', ( 0, 0 ) ]  # 1
    p2 = [ 'x', ( 0, 0 ) ]  # 2
    p3 = [ 'y', ( 0, 0 ) ]  # 3

    bp_123 = BasePointTree()
    bp_123.add( p1[0], p1[1], 1 )
    bp_123.add( p2[0], p2[1], 1 )
    bp_123.add( p3[0], p3[1], 1 )

    bp_1 = BasePointTree()
    bp_1.add( p1[0], p1[1], 1 )

    bp_2 = BasePointTree()
    bp_2.add( p2[0], p2[1], 1 )

    bp_3 = BasePointTree()
    bp_3.add( p3[0], p3[1], 1 )

    ls = LinearSeries.get( [4], bp_123 )
    ls_1 = LinearSeries.get( [1], bp_1 )
    ls_2 = LinearSeries.get( [1], bp_2 )
    ls_3 = LinearSeries.get( [1], bp_3 )

    OrbTools.p( 'ls =\n', ls )
    OrbTools.p( 'ls_1 =\n', ls_1 )
    OrbTools.p( 'ls_2 =\n', ls_2 )
    OrbTools.p( 'ls_3 =\n', ls_3 )


    # compute projection
    p_lst = []
    p_lst += [ls.pol_lst[0]]
    p_lst += [ls.pol_lst[1]]
    p_lst += [ls.pol_lst[5]]
    p_lst += [ls.pol_lst[11]]
    OrbTools.p( 'p_lst =' )
    for p in p_lst: OrbTools.p( '\t\t', p )
    ring = PolyRing( 'x,y,z' )
    lst = LinearSeries( p_lst, ring )
    bpt = lst.get_bp_tree()
    OrbTools.p( bpt )

    # compose with inverse stereographic projection with center (1:0:0:0:1)
    ring = PolyRing( 'x,y,z,x0,x1,x2,x3' )
    x, y, z, x0, x1, x2, x3 = ring.gens()
    f_lst = '[x0^2+x1^2+x2^2+x3^2, 2*x0*x1, 2*x0*x2, 2*x0*x3, -x0^2+x1^2+x2^2+x3^2]'
    g_lst = [ p.subs( {x:x0, y:x1, z:x2} ) for p in ring.coerce( p_lst ) ]
    fg_lst = ring.coerce( compose_maps( f_lst, g_lst ) )
    OrbTools.p( 'fg_lst =' )
    for fg in fg_lst: OrbTools.p( '\t\t', sage_factor( fg ) )

    # stereographic project from center (1:0:0:0:-1)
    q_lst = [ fg_lst[0] + fg_lst[4], fg_lst[1], fg_lst[2], fg_lst[3] ]
    q_lst = [ q / 2 for q in q_lst ]
    q_lst = [ q.subs( {x0:x, x1:y, x2:z} ) for q in q_lst ]
    OrbTools.p( 'q_lst =' )
    for q in q_lst: OrbTools.p( '\t\t', sage_factor( q ) )
    OrbTools.p( 'gcd(q_lst) =', sage_gcd( q_lst ) )

    # check base locus
    ring = PolyRing( 'x,y,z' )
    ls = LinearSeries( q_lst, ring )
    OrbTools.p( ls.get_bp_tree() )
    return


    #
    # Lines through (0:0:1) in projective plane.
    #
    #           P^1xP^1 ---> P^2
    # ( s : t ; v : w ) |--> ( s/t : (v/w)*(s/t) : 1 ) = ( w*s : v*s : w*t )
    #
    #         S^1 ---> P^1
    #     (c0,s0) |--> (1-s0:c0 )
    #

    # compute reparametrization from the linear series of families
    ring = PolyRing( 'x,y,z,c0,s0,c1,s1', True )  # construct polynomial ring with new generators
    OrbTools.p( ring )
    pmz_lst = ring.coerce( q_lst )

    x, y, z, c0, s0, c1, s1 = ring.gens()
    X = 1 - s0; Y = c0;
    V = 1 - s1; W = c1;

    dct_1 = {x:W * X, y:V * X, z:W * Y };
    pmz_1_lst = [ pmz.subs( dct_1 ) for pmz in pmz_lst ]

    dct_2 = {z:W * X, y:V * X, x:W * Y };
    pmz_2_lst = [ pmz.subs( dct_2 ) for pmz in pmz_lst ]

    dct_3 = {x:W * X, z:V * X, y:W * Y };
    pmz_3_lst = [ pmz.subs( dct_3 ) for pmz in pmz_lst ]

    # output
    OrbTools.p( 'pmz_1_lst =\n', pmz_1_lst )
    OrbTools.p( 'pmz_2_lst =\n', pmz_2_lst )
    OrbTools.p( 'pmz_3_lst =\n', pmz_3_lst )


    # mathematica input
    ms = ''
    for pmz, AB in [[pmz_1_lst, '1'], [pmz_2_lst, '2'], [pmz_3_lst, '3']]:
        s = 'pmz' + AB + '=' + str( pmz ) + ';'
        s = s.replace( '[', '{' ).replace( ']', '}' )
        ms += '\n' + s
    OrbTools.p( 'Mathematica input =', ms )

    return

    # PovInput ring cyclide
    #
    pin = PovInput()

    pin.path = './' + get_time_str() + '_veronese_blowup/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, -7, 0 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 45, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.height = 200
    pin.width = 400
    pin.quality = 1
    pin.ani_delay = 10

    pin.impl = None

    pin.pmz_dct['A'] = ( pmz_1_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_2_lst, 0 )
    pin.pmz_dct['C'] = ( pmz_3_lst, 0 )

    pin.pmz_dct['FA'] = ( pmz_1_lst, 0 )
    pin.pmz_dct['FB'] = ( pmz_2_lst, 0 )
    pin.pmz_dct['FC'] = ( pmz_3_lst, 0 )


    v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 5 )]  # 10
    v1_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 10 )]  # 15

    v1_lst_F = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 2 )]

    prec = 50

    pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':prec, 'width':0.05}
    pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':prec, 'width':0.05}
    pin.curve_dct['C'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':prec, 'width':0.05}

    pin.curve_dct['FA'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}
    pin.curve_dct['FB'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}
    pin.curve_dct['FC'] = {'step0':v0_lst, 'step1':v1_lst_F, 'prec':prec, 'width':0.02}

    col_a = ( 0.6, 0.4, 0.1, 0.0 )
    col_b = ( 0.1, 0.15, 0.0, 0.0 )
    col_c = ( 0.1, 0.0, 0.0, 0.0 )
    pin.text_dct['A'] = [True, col_a, 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, col_b, 'phong 0.2 phong_size 5' ]
    pin.text_dct['C'] = [True, col_c, 'phong 0.2 phong_size 5' ]

    col_F = ( 0.1, 0.1, 0.1, 0.0 )
    pin.text_dct['FA'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FB'] = [True, col_F, 'phong 0.2 phong_size 5' ]
    pin.text_dct['FC'] = [True, col_F, 'phong 0.2 phong_size 5' ]


    # raytrace image/animation
    create_pov( pin, ['A', 'B', 'C'] )
    return
    create_pov( pin, ['A', 'B', 'C', 'FA', 'FB', 'FC'] )





    # p_lst += [ls.pol_lst[0]]
    # p_lst += [ls.pol_lst[4] + ls.pol_lst[1]]
    # p_lst += [ls.pol_lst[8] + ls.pol_lst[9]]
    # p_lst += [ls.pol_lst[11]]


#     for c_lst in sage_Combinations( range( 12 ), 4 ):
#         OrbTools.p( 'c_lst =', c_lst )
#         OrbTools.p( 'ls.pol_lst =', len( ls.pol_lst ) , ls.pol_lst )
#
#         p_lst = [ ls.pol_lst[c] for c in c_lst ]
#         OrbTools.p( 'p_lst =' )
#         for p in p_lst: OrbTools.p( '\t\t', p )
#
#         # check base locus
#         ring = PolyRing( 'x,y,z' )
#         lst = LinearSeries( p_lst, ring )
#
#         try:
#             bpt = lst.get_bp_tree()
#             OrbTools.p( bpt )
#         except Exception as e:
#             continue
#
#         out_lst = []
#         for chart in bpt.chart_lst:
#             out_lst += trace_bp_tree( bpt[chart] )
#         OrbTools.p( 'out_lst =', out_lst )
#         if sum( out_lst ) == len( out_lst ) and len( out_lst ) > 7:
#             return
#
#     return



