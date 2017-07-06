'''
Created on Aug 7, 2016
@author: Niels Lubbes

Here functionality of the "orbital" package is tested.
Our methods use the Sage computer algebra libraries.  
The method names in this module are of the form: 
    "test_[method name to be tested]_[index]()"
For output we use "op()" in "verbose.py".
'''

from sage.all import *

import os
import inspect

from class_orb_input import *
from class_orb_output import *
from class_orb_ring import *
from orb_matrices import *
from orb_product import *
from orb_verbose import *
from orb_examples import *
from class_pov_input import *
from povray_aux import *
from povray import *
from surf_in_quad import *
from surf_examples import *

from linear_series.class_linear_series import *
from linear_series.class_base_points import *
from linear_series.verbose import *


def test_pov_coef_lst_0():
    '''
    Test "povray_aux.pov_exp_lst()"
    '''

    d = 2
    v = 3 + 1
    lst = []
    for elt in pov_exp_lst( d, v ):
        lst += [var( 'x' ) ** elt[0] * var( 'y' ) ** elt[1] * var( 'z' ) ** elt[2]]

    op( len( lst ), lst )
    op( binomial( d + v - 1, v - 1 ) )
    op( pov_coef_lst( 'x^2+y^2+z^2-1' ) )


def test_get_cs_0():
    '''
    Show case how rational points on "x^2+y^2-1=0" are obtained
    using the module "cos_sin".
    '''
    pt_dct = get_pt_dct();

    for angle in pt_dct.keys():
        op( angle, '=', pt_dct[angle] )


def test_orb_matrices_0():
    '''
    Show case the module "orb_matrices".
    The 9x9 matrices over QQ(c0,s0) represent automorphisms of 
    the projective 7-sphere S^7. These automorphisms are 
    called Moebius transformations.  
    '''

    #
    # Compute examples of orthogonal matrices.
    #
    op( 10 * '=' )
    op( 'Compute examples of orthogonal matrices.' )
    op( 10 * '=' )
    op( 'get_xmat( "X[0,90,0]")\n' + str( get_xmat( 'X[0,90,0]' ) ) )
    op( 'get_emat( "E[2,1,4,3,5,6,7,8]" )\n' + str( get_emat( 'E[2,1,4,3,5,6,7,8]' ) ) )
    op( 'get_omat( "Oprpp" )             \n' + str( get_omat( 'Oprpp' ) ) )
    op( 'get_rmat( "Rprpp[0,90,0,0]" )   \n' + str( get_rmat( 'Rprpp[0,90,0,0]' ) ) )
    op( 'get_rmat( "Rprpp[0,45,0,0]" )   \n' + str( get_rmat( 'Rprpp[0,45,0,0]' ) ) )
    op( 'get_mat( "I", "Orppp", "I" )    \n' + str( get_mat( 'I', 'Orppp', 'I' ) ) )

    #
    # Permute a block matrices to obtain rotation along x-axis.
    #
    op( 10 * '=' )
    op( 'Permute a block matrices to obtain rotation along x-axis.' )
    op( 10 * '=' )
    e_lst = [1, 4, 2, 3, 5, 6, 7, 8]
    eI_lst = Permutation( e_lst ).inverse()
    o_tup = ( 'E' + str( eI_lst ), 'Oprpp', 'E' + str( e_lst ) )
    o_mat = get_mat( *o_tup )
    op( 'Permuted omat:\n' + str( o_mat ) )
    op( list( o_mat ) )

    #
    # Compute random permutation of omat.
    #
    op( 10 * '=' )
    op( 'Compute random permutation of omat.' )
    op( 10 * '=' )
    p = Permutations( range( 1, 8 + 1 ) )
    rp = p.random_element()
    rpI = Permutation( rp ).inverse()
    A_str = 'E' + str( list( rpI ) )
    B_str = 'Orrrr'
    C_str = 'E' + str( list( rp ) )
    emat = get_emat( C_str )
    ematI = get_emat( A_str )
    pomat = get_mat( A_str, B_str, C_str )
    op( 'A_str = ', A_str )
    op( 'B_str = ', B_str )
    op( 'C_str = ', C_str )
    op( 'emat*ematI==1: ', emat * ematI == identity_matrix( 9 ) )
    op( 'Randomly permuted orbital matrix:\n' + str( pomat ) )

    #
    # Verify Euclidean translations of S^7.
    #
    op( 10 * '=' )
    op( 'Verify Euclidean translations of S^7.' )
    op( 10 * '=' )

    #
    # Setup inverse stereographic projection:
    # S: P^7 ---> S^7
    #
    v = OrbRing.coerce( '[v0,v1,v2,v3,v4,v5,v6,v7,v8]' )
    d = v[1] ** 2 + v[2] ** 2 + v[3] ** 2
    v0_2 = v[0] ** 2
    vec_lst = []
    vec_lst += [v0_2 + d]
    vec_lst += [2 * v[0] * v[1]]
    vec_lst += [2 * v[0] * v[2]]
    vec_lst += [2 * v[0] * v[3]]
    vec_lst += [2 * v[0] * v[4]]
    vec_lst += [2 * v[0] * v[5]]
    vec_lst += [2 * v[0] * v[6]]
    vec_lst += [2 * v[0] * v[7]]
    vec_lst += [-v0_2 + d]
    x = vector( OrbRing.R, vec_lst )
    op( 'Inverse stereographic projection S:P^7--->S^7,\n' + str( x ) )

    # Setup Euclidean translation in S^7
    #    T:S^7--->S^7
    #
    tmat = get_tmat( None )
    op( 'Euclidean rotation T:S^7--->S^7,\n' + str( tmat ) )
    op( list( tmat ) )

    # Setup stereographic projection
    #     P:S^7--->P^7
    #
    pmat = []
    pmat += [[1, 0, 0, 0] + [ 0, 0, 0, 0, -1]]
    pmat += [[0, 1, 0, 0] + [ 0, 0, 0, 0, 0]]
    pmat += [[0, 0, 1, 0] + [ 0, 0, 0, 0, 0]]
    pmat += [[0, 0, 0, 1] + [ 0, 0, 0, 0, 0]]
    pmat += [[0, 0, 0, 0] + [ 1, 0, 0, 0, 0]]
    pmat += [[0, 0, 0, 0] + [ 0, 1, 0, 0, 0]]
    pmat += [[0, 0, 0, 0] + [ 0, 0, 1, 0, 0]]
    pmat += [[0, 0, 0, 0] + [ 0, 0, 0, 1, 0]]
    pmat = matrix( pmat )
    op( 'Projection P:S^7--->P^7,\n' + str( pmat ) )


    # Check whether composition is an Euclidean
    # translation in P^7:
    #    PoToS
    #
    tv = pmat * tmat * x
    op( 'Composition PoToS is Euclidean translation of P^7: P^7--->P^7,\n', tv / 2 )


def test_orbital_product_0():
    '''
    Test "orb_product.orb_product()" with 
    examples from "orb_examples" module.
    '''

    #
    # Show long string representation of "OrbOutput" object.
    #
    OrbOutput.short_str = False

    #
    # compute orbital product for examples in "orb_examples"
    #
    exm_name_lst = get_exm_name_lst()  # all example names
    for exm_name in ['s3_perseus']:
        oin, pin = get_exm( exm_name )
        o = orb_product( oin )
        op( o )


def test_orbital_product_1( return_if_found = True, analyze_deg_4 = False ):
    '''
    INPUT:
        - "return_if_found" -- Boolean. 
        - "analyze_deg_4"   -- Boolean.
    
    OUTPUT:
        - Outputs via "orb_verbose.op(OrbOutput)" random orbital product 
          surfaces in S^7, in an infinite while-loop.

        - If "analyze_deg4" is True then the quartic product surfaces in 
          the projective 3-sphere S^3 are analyzed in more detail.
        
        - If "return_if_found" is True then the methods returns if 
          the randomly computed surface is one of the following:
          
              * A sextic surface in S^5. 
                
              * The quartic Veronese surface in S^4. 
              
              * An Euclidean translational surface in S^n for n>3.
              
              * A quartic surface in S^3 whose projection into 
                projective 3-space has only 1 singular component.
                Only applies if "analyze_deg4==True".
           

        
    '''
    OrbOutput.short_str = True
    halt = False
    while True:
        try:
            # random input
            input = OrbInput().random( 3, False )

            # only compute dimension and degree
            for key in input.do.keys(): input.do[key] = False
            input.do['imp'] = True
            input.do['dim'] = True
            input.do['deg'] = True

            # compute orbital product
            o = orb_product( input )
            op( o )

            # handle interesting surfaces
            #
            if  ( o.dim, o.deg ) in [( 5, 6 ), ( 4, 4 )]:
                halt = True

            if 'tT' in input.info_dct['omat'] and o.dim > 3:
                halt = True

            if ( o.dim, o.deg ) == ( 3, 4 ) and analyze_deg_4:

                for key in input.do.keys():
                    input.do[key] = True
                input.do['bpt'] = False

                o = orb_product( input )

                bln = o.prj_pol.total_degree() > 2
                bln = bln and o.sng_lst != None
                if bln and len( o.sng_lst ) == 1:
                    halt = True

            # halt if interesting surface is found
            if return_if_found and halt:
                OrbOutput.short_str = False
                op( o )
                op( 'Returning because interesting surface was found...' )
                return

        except BaseException as e:

            op( 'EXCEPTION:', e )
            op( input )


def test_orbital_product_2():
    '''
    Test "orb_product.orb_product()" 
    using "OrbOutput.get_short_str()" values.    
    '''

    #
    # Show long string representation of "OrbOutput" object.
    #
    OrbOutput.short_str = False

    ######################################################
    # Start specification of "s_lst" which is a list of  #
    # "OrbOutput.get_short_str()" outputs.               #
    ######################################################

    s_lst = []

    #
    # A dictionary with elements in s_lst as keys and
    # "OrbitInput.do_dct" dictionaries as values.
    #
    dct = {}
    true_do_dct = OrbInput().do  # default values are all true
    false_do_dct = { key:False for key in true_do_dct.keys() }

    #
    # --- degree 6 in S^5
    #
    s = "['@(5,6)=(dim,deg)', {'pmat': ('P1', 'I', 'I'), 'omat': ('T[1, -1, 1, 1, -1, 0, 0]', 'Opsms', 'I'), 'vmat': ('T[1, 0, -1, 0, 0, 0, 0]', 'Rramr[340, 225, 264, 320]', 'T[-1, 0, 1, 0, 0, 0, 0]')}]"
    dct[s] = copy( false_do_dct )
    dct[s]['pmz'] = True
    dct[s]['bpt'] = False
    dct[s]['imp'] = True
    dct[s]['dim'] = True
    dct[s]['deg'] = True
    dct[s]['prj'] = True
    dct[s]['fct'] = True
    # s_lst += [s]

    s = "['@(5,6)=(dim,deg)', {'pmat': ('P1', 'I', 'I'), 'omat': ('T[0, 0, 0, 1, 1, -1, -1]', 'Oprps', 'I'), 'vmat': ('T[0, -1, 0, -1, 0, 0, 0]', 'Rspps[148, 344, 284, 304]', 'T[0, 1, 0, 1, 0, 0, 0]')}]"
    dct[s] = copy( false_do_dct )
    dct[s]['pmz'] = True
    dct[s]['bpt'] = True  # singular
    # s_lst += [s]

    #
    # --- deg 4 in S^3, stereoprojection has no isolated singularities
    #
    s = "['@(3,4)=(dim,deg)', {'pmat': ('P0', 'I', 'I'), 'omat': ('T[1, 0, 0, 0, 0, 0, 0]', 'Orppp', 'T[-1, 0, 0, 0, 0, 0, 0]'), 'vmat': ('T[0, 1, 1, 0, 0, 0, 0]', 'Rrrrs[37, 0, 0, 0]', 'T[0, -1, -1, 0, 0, 0, 0]')}]"
    dct[s] = copy( true_do_dct )
    dct[s]['bpt'] = False  # True terminates
    dct[s]['sng'] = True
    dct[s]['tst'] = True
    # s_lst += [s]


    # s="['@(3,4)=(dim,deg)', {'pmat': ('P0', 'I', 'I'), 'omat': ('T[1, 0, -1, 0, 1, 0, -1]', 'Osaap', 'I'), 'vmat': ('T[1, 0, -1, -1, -1, 1, -1]', 'Rmrrs[169, 23, 153, 164]', 'T[0, 1, 0, -1, -1, 0, 0]')}]"
    # s="['@(4,4)=(dim,deg)', {'pmat': ('P0', 'I', 'I'), 'omat': ('I', 'Orsrr', 'I'), 'vmat': ('M[(0, 0, 0, 0, 1, 1, 0, 1, 0), (0, 0, 0, 1, 1, 1, 0, 0, 0), (0, 1, 0, 0, 0, 1, 1, 1, 1), (0, 0, 0, 1, 0, 0, 0, 1, 1), (1, 1, 0, 0, 1, 1, 0, 0, 1), (0, 1, 0, 0, 1, 1, 0, 1, 1), (1, 1, 1, 0, 1, 1, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 1, 1, 1, 1, 0, 1)]', 'I', 'I')}]"
    # s="['@(4,4)=(dim,deg)', {'pmat': ('P0', 'I', 'I'), 'omat': ('I', 'Orrss', 'I'), 'vmat': ('M[(0, 0, 0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 1, 1, 0, 0, 1, 1), (0, 0, 0, 0, 0, 1, 1, 0, 1), (1, 1, 1, 0, 1, 0, 1, 0, 0), (1, 0, 1, 0, 0, 0, 1, 1, 1), (1, 1, 1, 1, 0, 0, 1, 0, 0), (0, 0, 1, 0, 1, 0, 1, 1, 1), (0, 0, 0, 1, 0, 0, 1, 0, 0), (1, 1, 1, 1, 1, 1, 0, 1, 1)]', 'I', 'I')}]"


    ######################################################
    # End specifiation of "s_lst".                       #
    ######################################################

    #
    # Compute orbital product of entries in s_lst.
    #
    for s in s_lst:
        op( 'short_str =', s )
        input = OrbInput().set_short_str( s )
        input.do = dct[s]
        o = orb_product( input )
        op( o )


def test_orbital_product_3():
    '''
    Create random examples where the circle which
    is rotated is obtained as a random hyperplane section.
    The hyperplane section is obtained by applying 
    a random matrix in PGL(8+1) to the hyperplane 
    section { x | x3=..=x8=0=-x0^2+x1^2+...+x8^2=0 }.  
    Since the random matrix is not orthogonal
    we cannot use the matrix to compute parametric 
    presentations for these examples. 
    See "test_orbital_product_3()" for computing
    random orthogonal matrices.
    '''
    # enable short output string for
    # OrbOutput object
    OrbOutput.short_str = True

    # continue computing random surfaces
    # until it has the desired embedding
    # dimension and degree
    dim, deg = -1, -1
    while ( dim, deg ) not in [( 4, 6 ), ( 5, 6 ), ( 4, 4 )]:
        try:

            # compute random input
            input = OrbInput()
            pmat = list( MatrixSpace( ZZ, 9, 9 ).random_element() )
            ch_lst = ['r', 's', 'p', 'm', 'a']
            ch_str = ''.join( [OrbRing.random_elt( ch_lst ) for i in range( 4 ) ] )
            p_tup = ( 'P1', 'I', 'I' )
            o_tup = ( 'I', 'O' + ch_str, 'I' )
            v_tup = ( 'M' + str( pmat ), 'I', 'I' )
            input.set( p_tup, o_tup, v_tup )

            # only compute dimension and degree
            for key in input.do.keys(): input.do[key] = False
            input.do['imp'] = True
            input.do['dim'] = True
            input.do['deg'] = True
            o = orb_product( input )
            op( o )

            # update dim an deg for while loop
            dim, deg = o.dim, o.deg

        except BaseException as e:

            # continue while loop regardless of exceptions
            op( 'EXCEPTION:', e )

    #
    # Compute with previous input in more detail
    #
    OrbOutput.short_str = False
    for key in input.do.keys(): input.do[key] = True
    input.do['pmz'] = False
    input.do['bpt'] = False
    input.do['tst'] = False

    o = orb_product( input )
    op( o )
    op( list( pmat ) )


def test_create_pov_0( fast_low_quality = True ):
    '''
    Test "povray.create_pov()" for creating images 
    and animations of "OrbOutput" objects of
    examples in the "orb_examples" module.
    '''

    exm_name_lst = get_exm_name_lst()
    for exm_name in ['s3_ring']:

        oin, pin = get_exm( exm_name )

        o = orb_product( oin )
        op( o )
        pin.pmz_dct['A'] = ( o.prj_pmz_lst, 0 )
        pin.pmz_dct['B'] = ( o.prj_pmz_lst, 1 )

        # execute faster with much too low quality,
        # for debugging purposes.
        #
        if fast_low_quality:
            pin.width = 800
            pin.height = 400
            pin.quality = 4
            pin.ani_delay = 1
            pin.curve_dct['A'] = {'step0':36, 'step1':2 * 36, 'prec':5, 'width':0.02}
            pin.curve_dct['B'] = {'step0':36, 'step1':2 * 36, 'prec':5, 'width':0.02}
        op( 'pin.path = ', pin.path )

        # raytrace image/animation
        #
        lst = create_pov( pin, ['A', 'B'], False, False, [] )
        create_pov( pin, ['A'], False, False, [] )
        create_pov( pin, ['B'], False, False, [] )
        if not fast_low_quality:
            create_pov( pin, ['A', 'B'], False, True, [] )
        # create_pov( pin, ['A', 'B'], False, True, [( 'A', 0.5 ), ( 'B', 0.5 )] )

        # save output
        #
        dct = {}
        dct['pin'] = pin
        dct['pov_str_A'] = lst[0]
        dct['pov_str_B'] = lst[1]
        dct['pov_str_C'] = ' '
        get_orb_dct()['test_create_pov_0_' + exm_name] = dct
        save_orb_dct()


def test_create_pov_1( fast_low_quality = False ):
    '''
    Test "povray.create_pov()" for creating images 
    and animations.
    We load examples from "surf_examples"
    '''
    # obtain output
    #
    dct = con_s5d6_smooth_0()
    f_xyz = dct['f_xyz']
    pmz_AB_lst = dct['pmz_AB_lst']
    pmz_CB_lst = dct['pmz_CB_lst']
    pin = dct['pin']

    # execute faster with much too low quality,
    # for debugging purposes.
    #
    if fast_low_quality:
        pin.width = 100
        pin.height = 75
        pin.quality = 4
        pin.ani_delay = 1
        pin.curve_dct['A'] = {'step0':2 * 36, 'step1':2 * 36, 'prec':5, 'width':0.02}
        pin.curve_dct['B'] = {'step0':2 * 36, 'step1':2 * 36, 'prec':5, 'width':0.02}
        pin.curve_dct['C'] = {'step0':2 * 36, 'step1':2 * 36, 'prec':5, 'width':0.02}
    op( 'pin.path = ', pin.path )

    # raytrace image/animation
    #
    lst = create_pov( pin, ['A', 'B', 'C'], False, False, [] )
    create_pov( pin, ['A'], False, False, [] )
    create_pov( pin, ['B'], False, False, [] )
    create_pov( pin, ['C'], False, False, [] )
    if not fast_low_quality:
        create_pov( pin, ['A', 'B', 'C'], False, True, [] )
        create_pov( pin, ['A', 'B', 'C'], False, True, [( 'A', 0.5 ), ( 'B', 0.5 ), ( 'C', 0.5 )] )

    # save output
    #
    dct['pin'] = pin
    dct['pov_str_A'] = lst[0]
    dct['pov_str_B'] = lst[1]
    dct['pov_str_C'] = lst[2]
    get_orb_dct()['test_create_pov_1'] = dct
    save_orb_dct()


def test_surf_in_quad():
    '''
    Test "surf_in_quad" module.
    
    We construct a parametrization of a surface via
    the map associated to a linear series.
    We compute its implicit equations and find a quadratic 
    form with given signature in its ideal.
    We find a projective transformation of the surface
    into the hyperquadric which is defined by a diagonal form
    with given signature.        
    '''

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

    # implicit equation and parametrization in P^7.
    #
    sig = ( 1, 5 )
    coef_lst = None
    prv_Q = None
    # coef_lst = [1, 1, -1, -1, 0, 1, 0, 0, 1, 0, -1] # Uncomment for using precomputed quadric in ideal with signature (1,6).
    # coef_lst = [1, 0, 1, -1, 0, 0, 0]  # Uncomment for using precomputed quadric in ideal with signature (1,5) with Q as below:
    # prv_Q = [( 0, 1, 1, 1, 0, 1, 0 ), ( 1, 0, 0, 1, 1, 1, 0 ), ( 0, 1, 1, 1, 0, 1, 1 ), ( 1, 0, 1, 1, 0, 1, 1 ), ( 0, 0, 0, 0, 1, 0, 0 ), ( 0, 0, 0, 1, 1, 1, 0 )]

    dct = get_surf( ls_AB, sig, coef_lst, prv_Q )
    op( 'test_get_surf =', test_get_surf( dct ) )

    # define projection P: P^5--->P^3
    # such that image is compact.
    #
    approxU = approx_QQ( dct['UJ'][0] )
    P = get_prj_mat( 4, 6, 0 )
    P[0, 5] = -1
    P = P * approxU
    op( 'P             =', list( P ) )

    # families A, B in P^3.
    #
    f_xyz, pmz_AB_lst = get_proj( dct['imp_lst'], dct['pmz_lst'], P )

    # families C, B in P^3.
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

    # save intermediate results.
    #
    dct = {}
    dct['f_xyz'] = f_xyz
    dct['pmz_AB_lst'] = pmz_AB_lst
    dct['pmz_CB_lst'] = pmz_CB_lst
    get_orb_dct()['test_surf_in_quad'] = dct
    save_orb_dct()


if __name__ == '__main__':

    # for some functionality we use maple and/or magma
    #
    os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'
    os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'
    if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = '/home/niels/Desktop/n/src/output/orb/povray/'

    #  Debug output settings
    #
    ot( True )  # show timing
    op( True )  # show all output
    dprint( False, 'no output' )  # surpress output from linear series package

    #
    # Uncomment below for showing only output when called
    # from 'orb_main.py' module.
    #
    # op( False, 'orb_main.py' )


    #########################################
    #                                       #
    # Uncomment one or more test methods    #
    #                                       #
    #########################################

    # test_get_cs_0()  # create list of rational points on circle
    # test_pov_coef_lst_0() # create input string for polynomial for povray's POLY.
    # test_orb_matrices_0()  # showcase matrix representations of Aut(S^7)
    # test_orbital_product_0()  # investigate examples in "orb_examples" module
    # test_orbital_product_1()  # compute random examples (stop/analyze deg4)
    # test_orbital_product_2()  # investigate examples from "short_str" examples
    # test_orbital_product_3()  # compute random examples without parametrization
    # test_create_pov_0()  # create images and animations using "orb_examples".
    # test_create_pov_1()  # raytrace examples from "surf_examples"
    test_surf_in_quad()  # construct deg 6 del Pezzo surface in S^5

    #########################################
    #                                       #
    # End of list of test methods.          #
    #                                       #
    #########################################

    # end timing
    ot()

    print
    print 'The End'
