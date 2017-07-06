'''
Created on Sep 15, 2016

@author: Niels Lubbes
'''
from sage.all import *

from class_orb_ring import *
from orb_verbose import *
from class_pov_input import *
from povray_aux import *
from surf_in_quad import *

from linear_series.class_linear_series import *
from linear_series.class_base_points import *
from linear_series.verbose import *



def con_s5d6_smooth_0():
    '''
    We use the "surf_in_quad" module.
    
    We consider 3 families of conics of a degree 6 
    del Pezzo surface in S^5 which is a priori not 
    an orbital product.
    
    The NS lattice of deg 6 del Pezzo is
    
        ZZ< h1, h2, e1, e2 >
    
    with h1*h2=1 and other intersections are 0.
    
    Here e1 and e2 are classes of the blow up of 
    the points:
        (1:a0;1:-a0) and (1:-a0;1:a0) 
    in P^1xP^1.
    
    Note that if we fix eg. a coordinate in the 
    right P^1 then we parametrize a conic in family.
    Thus a map P^1xP^1--->P^6 of bidegree (d1,d1)
    represents 2 families of degree d1 and d2 curves
    respectively. 
    
    Let "pmz_AB_lst" represent a map P^1xP^1--->P^6
    with image a degree 6 del Pezzo contained in 
    a quadric of signature (1,6). Let "M" be the 
    Gramm matrix of this quadric. We compute
    the orthogonal diagonalization of "M", such
    that we obtain a projective transformation to
    a degree 6 del Pezzo surface in the 5-sphere.
    We compose this transformation to a projection
    into projective 3-space.
    
    The class of "pmz_AB_lst" is 2(h1+h2)-e1-e2
    and thus of bidegree (2,2).
    
    The classes of the families are
        family A: h1
        family B: h2
        family C: h1+h2-e1-e2.
    
    The linear series corresponding to h1+h2-e1-e2 is 
        [x*v-y*w, y*v+x*w ]
    thus a conic in P^1xP^1 whose parametric image 
    is a conic in family C has equation:
    
        a*(x*v-y*w) == b*(y*v+x*w).
    
    We can parametrize a conic---that is the preimage of a
    conic in family C---in P^1xP^1 via the map
    
        P^1xP^1  ---->  P^1xP^1
       (x:y;a:b) |--->  (x:y;x*b+y*a:x*a-y*b)
    
    We obtain the parametrization of the conic in the 
    the family C by composing with the parametrization of 
    the surface. 
        
    See also <linear_series.ls_main.py> for computing 
    the linear series corresponding to  
    2(h1+h2)-e1-e2 and h1+h2-e1-e2.    
    
    See below for some of the values of the parameters that
    are computed:
    
    pmz_AB_lst = [
         'x^2*v^2 - y^2*w^2',
         'x^2*v*w + y^2*v*w',
         'x^2*w^2 + y^2*w^2',
         'x*y*v^2 - y^2*v*w',
         'x*y*v*w - y^2*w^2',
         'y^2*v*w + x*y*w^2',
         'y^2*v^2 + y^2*w^2'
         ]
    
    imp_lst = [ 'x3*x5 + x5^2 - x2*x6 - x4*x6',
                'x4^2 + x5^2 - x2*x6',
                'x3*x4 + x4*x5 - x1*x6 + x5*x6',
                'x2*x4 - x1*x5 + x2*x6',
                'x1*x4 - x0*x5 + x1*x6 - x5*x6',
                'x3^2 - x5^2 - x0*x6 + x2*x6 + 2*x4*x6',
                'x2*x3 - x0*x5 + x1*x6 - x5*x6',
                'x1*x3 - x0*x4 - x4*x6',
                'x1^2 - x0*x2 - x2*x6',
                'x0*x5^2 + x2*x5^2 - x2^2*x6 - 2*x1*x5*x6 + x5^2*x6 + x2*x6^2',
                'x0*x4*x5 + x1*x5^2 - x1*x2*x6 - x0*x5*x6 + x4*x5*x6 + x1*x6^2 - x5*x6^2'
                ]
    
    M = [
        (   0,   0, 1/2,    0,    0, -1/2,  1/2 ), 
        (   0,  -1,   0,    0,    0,    0,  1/2 ), 
        ( 1/2,   0,   0,  1/2,    0,    0,    1 ), 
        (   0,   0, 1/2,   -1,    0, -1/2,    0 ), 
        (   0,   0,   0,    0,   -1,    0, -1/2 ), 
        (-1/2,   0,   0, -1/2,    0,   -1, -1/2 ), 
        ( 1/2, 1/2,   1,    0, -1/2, -1/2,    0 )
        ]
    '''

    # Uncomment following line for disabling output of debug info:
    #
    # op( False, 'orb_main.py' )
    dprint( False, 'no output' )

    # QQ[i]/<i^2+1>
    #
    ring = PolyRing( 'x,y,v,w' )
    ring.ext_num_field( 't^2 + 1' )
    a0 = ring.root_gens()[0]

    # e1, e2
    #
    bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )
    bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )
    op( 'bp_tree =', bp_tree )

    # |2(h1+h2)-e1-e2|
    #
    ls_AB = LinearSeries.get( 2, bp_tree )
    op( 'ls_AB =', ls_AB.get_bp_tree() )

    # |h1+h2-e1-e2|
    #
    ls_CB = LinearSeries.get( 1, bp_tree )
    op( 'ls_CB =', ls_CB.get_bp_tree() )

    # implicit equation and parametrization in P^7
    #
    c_lst = [-1, -1, 0, 0, 0, -1, 1, 0, -1, -1, -1]
    dct = get_surf( ls_AB, ( 6, 1 ), c_lst )
    U, J = dct['UJ']
    U.swap_rows( 0, 6 )
    J.swap_columns( 0, 6 )
    J.swap_rows( 0, 6 )
    op( 'test_get_surf =', test_get_surf( dct ) )

    # define projection P: P^7--->P^3
    # such that image is compact
    #
    approxU = approx_QQ( dct['UJ'][0] )
    P = get_prj_mat( 4, 7, 0 )
    P[0, 6] = -1;P[3, 3] = 0;P[3, 4] = 1
    P = P * approxU
    op( 'P             =', list( P ) )

    # families A, B in P^3
    #
    f_xyz, pmz_AB_lst = get_proj( dct['imp_lst'], dct['pmz_lst'], P )

    # families C, B in P^3
    #
    R_xyvw = PolynomialRing( QQ, var( 'x,y,v,w' ) )
    x, y, v, w = R_xyvw.gens()
    X, Y, V, W = var( 'x,y,v,w' )
    xyvw_dct = { X:x, Y:y, V:v, W:w }
    XYZW_dct = { x:X, y:Y, v:V, w:W }
    CB_dct = { x:X, y:Y, v:X * W + Y * V, w: X * V - Y * W }
    pol_lst = sage_eval( str( ls_AB.pol_lst ), R_xyvw.gens_dict() )
    pmz_CB_lst = [ pol.subs( CB_dct ).subs( xyvw_dct ) for pol in pol_lst]
    pmz_CB_lst = get_S1xS1_pmz( pmz_CB_lst )
    pmz_CB_lst = list( P * dct['Q'] * vector( pmz_CB_lst ) )

    op( 'f_xyz         =', f_xyz )
    op( 'pmz_AB_lst    =', pmz_AB_lst )
    op( 'pmz_CB_lst    =', pmz_CB_lst )

    # set PovInput as container
    #
    pin = PovInput()
    pin.impl = f_xyz
    pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
    pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
    pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
    pin.path = os.environ['OUTPUT_PATH'] + get_time_str() + '_fam_ABC/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, 0, QQ( -21 ) / 10 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 310, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10
    pin.curve_dct['A'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['C'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}

    # save intermediate results
    #
    dct = {}
    dct['f_xyz'] = f_xyz
    dct['pmz_AB_lst'] = pmz_AB_lst
    dct['pmz_CB_lst'] = pmz_CB_lst
    dct['pin'] = pin
    get_orb_dct()['con_s5d6_smooth_0'] = dct
    save_orb_dct()

    return dct













