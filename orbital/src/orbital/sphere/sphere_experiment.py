'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Aug 30, 2018
@author: Niels Lubbes

In order to interactively use the functionality of this module,
copy paste the following code into a Sage notebook:

------------------------------------------------------------------------------
reset()
from orbital.class_orb_tools import OrbTools
from orbital.sphere.class_sphere_input import SphereInput
from orbital.sphere.sphere_experiment import clifford

os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'
OrbTools.filter( [] ) 

lay = []
lay+=[['pmz','bas','mrk','frm']]
lay+=[['fam','famA','famB']]
lay+=[['stp'],['famt']]
lay+=[['ppt','rng','opa']]
lay+=[['prj']]
lay+=[['imp','snp','sng']]
lay+=[['a12','a13','a23']]
lay+=[['a14','a24','a34']]
lay+=[['u1','u2','u3']]
lay+=[['sa']]
lay+=[['b12','b13','b23']]
lay+=[['b14','b24','b34']]
lay+=[['v1','v2','v3']]
lay+=[['sb']]

@interact( layout = lay )
def orbit( 
           a12 = slider( 0, 360, 1, default = 0 ),
           a13 = slider( 0, 360, 1, default = 0 ),
           a23 = slider( 0, 360, 1, default = 0 ),
           a14 = slider( 0, 360, 1, default = 0 ),           
           a24 = slider( 0, 360, 1, default = 0 ),
           a34 = slider( 0, 360, 1, default = 0 ),
           u1 = slider( -10, 10, 1/4, default = 0 ),
           u2 = slider( -10, 10, 1/4, default = 0 ),
           u3 = slider( -10, 10, 1/4, default = 0 ),
           sa = slider( 0, 10, 1/4, default = 1 ),
           b12 = slider( 0, 360, 1, default = 0 ),
           b13 = slider( 0, 360, 1, default = 0 ),
           b23 = slider( 0, 360, 1, default = 45 ),
           b14 = slider( 0, 360, 1, default = 0 ),           
           b24 = slider( 0, 360, 1, default = 0 ),
           b34 = slider( 0, 360, 1, default = 0 ),
           v1 = slider( -10, 10, 1/4, default = 0 ),
           v2 = slider( -10, 10, 1/4, default = 0 ),
           v3 = slider( -10, 10, 1/4, default = 0 ),
           sb = slider( 0, 10, 1/4, default = 1 ),
           prj = slider(0,3,1, default = 0 ),           
           ppt = slider( 10, 500, 10, default = 30 ),
           rng = slider( 1, 10, 1, default = 2 ),
           opa = slider( 0, 1, 1/10, default = 1 ),           
           frm = checkbox(False),
           pmz = checkbox(False),           
           bas = checkbox(True),
           mrk = checkbox(True),
           fam = checkbox(False),
           famA = checkbox(True),
           famB = checkbox(True),
           famt = slider( 1, 10, 1, default = 2 ),     
           stp = slider( 1, 36, 1, default = 6 ),
           imp = checkbox(False),
           sng = checkbox(False),
           snp = checkbox(True)):

    # construct SphereInput
    lsta = [( a12, a13, a23 ), ( a14, a24, a34 ), ( u1, u2, u3 ), sa ]
    lstb = [( b12, b13, b23 ), ( b14, b24, b34 ), ( v1, v2, v3 ), sb ]
    sinp = SphereInput().set( [lsta] + [lstb] )
    sinp.prj = prj
    sinp.imp = imp
    sinp.rng = rng
    sinp.sng = sng
    sinp.snp = snp
    sinp.pmz = pmz
    sinp.opa = opa
    sinp.bas = bas
    sinp.mrk = mrk
    sinp.fam = fam
    sinp.famA = famA
    sinp.famB = famB
    sinp.famt = famt
    sinp.stp = stp
    sinp.col_pmz = Color( "orange" )
    sinp.col_imp = Color( "#ffefb0" )
    sinp.col_famA = Color( "red" )
    sinp.col_famB = Color( "blue" )    

    # compute product of circles
    plt, out = clifford( sinp )

    # show output
    print out
    show( plt, frame = frm )
------------------------------------------------------------------------------
'''

from orbital.sage_interface import sage_var
from orbital.sage_interface import sage_SR
from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_pi
from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_cos
from orbital.sage_interface import sage_sin
from orbital.sage_interface import sage_PolynomialRing
from orbital.sage_interface import sage_factor
from orbital.sage_interface import sage_Graphics
from orbital.sage_interface import sage_parametric_plot3d
from orbital.sage_interface import sage_implicit_plot3d
from orbital.sage_interface import sage_point3d
from orbital.sage_interface import sage_Color

from orbital.sphere.sphere_transform import get_rot_S3
from orbital.sphere.sphere_transform import get_trn_S3
from orbital.sphere.sphere_transform import get_scale_S3
from orbital.sphere.sphere_transform import get_hp_P4

from orbital.prod.orb_product import get_sing_lst
from orbital.class_orb_ring import OrbRing


def clifford( sinp ):
    '''
    Parameters
    ----------
    sinp : SphereInput    
    
    Returns
    -------
    tuple
        A tuple (plt,out) where 
        * plt : sage_Graphics object with plots of surface defined 
                by "sinp".
        * out : Information string.
        
    Notes
    -----
    See documentation at top of this file for possible usage of 
    this method.
    '''
    #
    # compute matrices A and B
    #
    T = get_trn_S3( sinp.trna )
    R = get_rot_S3( sinp.rota )
    S = get_scale_S3( sinp.sa )
    A = S * R * T

    T = get_trn_S3( sinp.trnb )
    R = get_rot_S3( sinp.rotb )
    S = get_scale_S3( sinp.sb )
    B = S * R * T

    #
    # Compute implicit and parametric form of a stereographic
    # projection of the Hamiltonian product of A and B.
    #
    baseA, baseB, pmzAB = get_pmz( A, B, sinp.prj )
    if sinp.imp:
        dct = get_imp( A, B, sinp.prj, sinp.sng, sinp.snp )
        key_lst = ['Agreat', 'Bgreat', 'eqn_x', 'eqn_str', 'eqn_xyz', 'sng_lst']
        Agreat, Bgreat, eqn_x, eqn_str, eqn_xyz, sng_lst = [ dct[key] for key in key_lst ]

    #
    # create graphics object
    #
    plt = sage_Graphics()
    a, b = sage_var( 'a,b' )

    if sinp.pmz:
        plt += sage_parametric_plot3d( pmzAB, ( a, 0, 2 * sage_pi ), ( b, 0, 2 * sage_pi ),
                                       color = sinp.col_pmz,
                                       aspect_ratio = 1, plot_points = sinp.ppt, opacity = sinp.opa )

    if sinp.bas:
        plt += sage_parametric_plot3d( baseA, ( a, 0, 2 * sage_pi ), color = sinp.col_famA, thickness = 10, aspect_ratio = 1 )
        plt += sage_parametric_plot3d( baseB, ( b, 0, 2 * sage_pi ), color = sinp.col_famB, thickness = 10, aspect_ratio = 1 )

    if sinp.mrk:
        plt += sage_point3d( ( 0, 0, 0 ), size = 30, color = sage_Color( "magenta" ) )
        plt += sage_point3d( ( 1, 0, 0 ), size = 30, color = sage_Color( "green" ) )

    if sinp.fam:
        if sinp.famA:
            for b1 in range( 0, 360, sinp.stp ):
                ps = [ pmzAB[i].subs( {b:b1 * sage_pi / 180} ) for i in [0, 1, 2]]
                plt += sage_parametric_plot3d( ps, ( a, 0, 2 * sage_pi ), color = sage_Color( "red" ),
                                               thickness = sinp.famt, aspect_ratio = 1, plot_points = sinp.ppt )
        if sinp.famB:
            for a1 in range( 0, 360, sinp.stp ):
                ps = [ pmzAB[i].subs( {a:a1 * sage_pi / 180} ) for i in [0, 1, 2]]
                plt += sage_parametric_plot3d( ps, ( b, 0, 2 * sage_pi ), color = sage_Color( "blue" ),
                                               thickness = sinp.famt, aspect_ratio = 1, plot_points = sinp.ppt )

    if sinp.imp:
        rng = sinp.rng
        x, y, z = sage_var( 'x,y,z' )
        plt += sage_implicit_plot3d( eqn_xyz, ( x, -rng, rng ), ( y, -rng, rng ), ( z, -rng, rng ),
                                     color = sinp.col_imp,
                                     plot_points = sinp.ppt, opacity = sinp.opa )

    #
    # create output string
    #
    out = ''
    out += str( sinp )
    out += '\n'
    if sinp.imp:
        out += '\neqn_str = ' + eqn_str
        out += '\nAgreat  = ' + str( Agreat )
        out += '\nBgreat  = ' + str( Bgreat )
        out += '\nA       = ' + str( list( A ) )
        out += '\nB       = ' + str( list( B ) )
        out += '\npmzAB   = ' + str( pmzAB )
        if sinp.sng:
            out += '\nsng_lst (long)  ='
            out += '\n-----'
            for sng in sng_lst:
                out += '\n' + str( sng )
            out += '\n-----'
            out += '\nsng_lst (short) ='
            out += '\n-----'
            for sng in sng_lst:
                out += '\n' + str( sng[1] )
            out += '\n-----'

    return plt, out


def get_pmz( A, B, prj ):
    '''
    Computes parametrization of a stereographic projection of 
    the pointwise Hamiltonian product of circles in the sphere S^3 
    defined by transformations of the standard circle [1,cos(a),sin(a),0,0]
    by A and B respectively.  

        
    Parameters
    ----------
    A : sage_Matrix<sage_QQ>
        Represents a linear transformation S^3--->S^3
    B : sage_Matrix<sage_QQ>
        Represents a linear transformation S^3--->S^3
    prj : int 
        Choice for stereographic projection S^3--->P^3: 
        0: (x0:x1:x2:x3:x4) |--> (x0-x4:x1:x2:x3)
        1: (x0:x1:x2:x3:x4) |--> (x0-x1:x4:x2:x3)
        2: (x0:x1:x2:x3:x4) |--> (x0-x2:x1:x4:x3)
        3: (x0:x1:x2:x3:x4) |--> (x0-x3:x1:x2:x4)  
    
    Returns
    -------
    tuple
        Returns tuple (baseA, baseB, pmzAB) where
        * baseA: Parametrization of projection of A in cos(a) and sin(a).
        * baseB: Parametrization of projection of B in cos(b) and sin(b).
        * pmzAB: Parametrization of projection of A*B.                     
    '''
    dct = {}

    # Hamiltonian product of circles
    a, b = sage_var( 'a,b' )
    u = list( A * sage_vector( [1, sage_cos( a ), sage_sin( a ), 0, 0] ) )
    v = list( B * sage_vector( [1, sage_cos( b ), sage_sin( b ), 0, 0] ) )
    p = get_hp_P4( u, v )

    # stereographic projection
    if prj == 0: j, i_lst = 4, [1, 2, 3]
    if prj == 1: j, i_lst = 1, [4, 2, 3]
    if prj == 2: j, i_lst = 2, [1, 4, 3]
    if prj == 3: j, i_lst = 3, [1, 2, 4]
    p = [ p[i] / ( p[0] - p[j] ) for i in i_lst ]

    # put in dictionary
    pmzAB = [ elt.full_simplify() for elt in p ]
    baseA = [ u[i] / ( u[0] - u[j] ) for i in i_lst ]
    baseB = [ v[i] / ( v[0] - v[j] ) for i in i_lst ]

    return baseA, baseB, pmzAB


def get_imp( A, B, prj, sng, snp ):
    '''
    Computes implicit equation of a stereographic projection S of 
    the pointwise Hamiltonian product of circles in the sphere S^3 
    defined by transformations of the standard circle [1,cos(a),sin(a),0,0]
    by A and B respectively.     
        
    Parameters
    ----------
    A : sage_Matrix<sage_QQ>
        Represents a linear transformation S^3--->S^3
    
    B : sage_Matrix<sage_QQ>
        Represents a linear transformation S^3--->S^3
    
    prj : int
        Choice for stereographic projection S^3--->P^3: 
        0: (x0:x1:x2:x3:x4) |--> (x0-x4:x1:x2:x3)
        1: (x0:x1:x2:x3:x4) |--> (x0-x1:x4:x2:x3)
        2: (x0:x1:x2:x3:x4) |--> (x0-x2:x1:x4:x3)
        3: (x0:x1:x2:x3:x4) |--> (x0-x3:x1:x2:x4)
    
    sng : boolean
        If true computes singular locus of S. Needs Magma path set
        in os.environ['PATH']. Otherwise the empty-list is returned.

    snp : boolean
        If true and if sng is True, then the singular locus is 
        computed with a probablistic method, which is faster but
        the correctness of the output is not guaranteed.
        
    Returns
    -------
    dict
        {   
            'Agreat' : boolean
                If True, then circle A is great.
                
            'Bgreat' : boolean
                If True, then circle B is great.
            
            'eqn_x'  : sage_PolynomialRing
                Equation of S in x0,...,x3

            'eqn_str':
                Formatted equation in x0,...,x3 of the 
                form f(x1:x2:x3)+x0*g(x0:x1:x2:x3).
            
            'eqn_xyz':
                Equation of S in x,y,z
                        
            'sng_lst':
                Empty-list if Magma is not installed or 
                list of singularities of S otherwise.
        }
        
    '''

    dct = {}  # output

    # create polynomial ring
    #
    R = sage_PolynomialRing( sage_QQ, 'x0,x1,x2,x3,a0,a1,a2,a3,a4,b0,b1,b2,b3,b4' )
    x0, x1, x2, x3, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4 = R.gens()

    # construct ideal for A
    #
    sv = [0, 0, 0, 1, 0]
    tv = [0, 0, 0, 0, 1]
    u0, u1, u2, u3, u4 = list( A * sage_vector( sv ) )
    v0, v1, v2, v3, v4 = list( A * sage_vector( tv ) )
    eqA = [-a0 ** 2 + a1 ** 2 + a2 ** 2 + a3 ** 2 + a4 ** 2]
    eqA += [u0 * a0 + u1 * a1 + u2 * a2 + u3 * a3 + u4 * a4]
    eqA += [v0 * a0 + v1 * a1 + v2 * a2 + v3 * a3 + v4 * a4]
    dct['Agreat'] = u0 == v0 == 0

    # construct ideal for B
    #
    u0, u1, u2, u3, u4 = list( B * sage_vector( sv ) )
    v0, v1, v2, v3, v4 = list( B * sage_vector( tv ) )
    eqB = [-b0 ** 2 + b1 ** 2 + b2 ** 2 + b3 ** 2 + b4 ** 2]
    eqB += [u0 * b0 + u1 * b1 + u2 * b2 + u3 * b3 + u4 * b4]
    eqB += [v0 * b0 + v1 * b1 + v2 * b2 + v3 * b3 + v4 * b4]
    dct['Bgreat'] = u0 == v0 == 0

    # stereographic projection
    #
    if prj == 0: j, i_lst = 4, [1, 2, 3]
    if prj == 1: j, i_lst = 1, [4, 2, 3]
    if prj == 2: j, i_lst = 2, [1, 4, 3]
    if prj == 3: j, i_lst = 3, [1, 2, 4]

    # construct equation of for projection of A*B
    #
    c = c0, c1, c2, c3, c4 = get_hp_P4( [a0, a1, a2, a3, a4], [b0, b1, b2, b3, b4] )
    x = [x0, x1, x2, x3]
    i1, i2, i3 = i_lst
    id = [ x[0] - ( c[0] - c[j] ), x[1] - c[i1], x2 - c[i2], x3 - c[i3] ] + eqA + eqB
    dct['eqn_x'] = eqn_x = R.ideal( id ).elimination_ideal( [a0, a1, a2, a3, a4, b0, b1, b2, b3, b4] ).gens()[0]

    # get equation in string form
    #
    f = eqn_x.subs( {x0:0} )
    dct['eqn_str'] = str( sage_factor( f ) ) + '+' + str( sage_factor( eqn_x - f ) )
    xs, ys, zs = sage_var( 'x,y,z' )
    dct['eqn_xyz'] = sage_SR( eqn_x.subs( {x0:1, x1:xs, x2:ys, x3:zs} ) )

    # compute singular locus
    #
    dct['sng_lst'] = []
    if sng:
        dct['sng_lst'] = get_sing_lst( OrbRing.coerce( eqn_x ), snp )

    return dct

