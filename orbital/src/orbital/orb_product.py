'''
Created on Jul 8, 2016
@author: Niels Lubbes
'''
from sage.all import *

from class_orb_input import *
from class_orb_output import *
from class_orb_ring import *
from orb_matrices import *
from orb_verbose import *

from linear_series.verbose import *
from linear_series.class_linear_series import *
from linear_series.class_base_points import *

import warnings


def get_dim( imp_lst ):
    '''
    INPUT:
        - "imp_lst" -- A list of homogenous polynomials in QQ[x0,...,x8]
                       representing a variety S in the projective 7-sphere.
    OUTPUT:
        - Return dimension of the variety S.
    '''
    dim = 7
    for imp in imp_lst:
        if imp.total_degree() == 1:
            dim = dim - 1
    return dim

def get_degree( pol_lst ):
    '''
    INPUT:
        - "pol_lst" -- A list of polynomials in x0,...,x8. The 
                       ideal of "pol_lst" should define a 
                       projective surface S in P^8.
    OUTPUT:
        - The degree of S.
    '''
    # consider ideal in ring of the right dimension.
    Ry = PolynomialRing( QQ, var( 'y0,y1,y2,y3,y4,y5,y6,y7,y8' ), order = 'degrevlex' )
    I = Ry.ideal( sage_eval( str( pol_lst ).replace( 'x', 'y' ), Ry.gens_dict() ) )

    # compute Hilbert polynomial: deg*t^dim + ...
    hpol = I.hilbert_polynomial()
    hdeg = hpol.diff().diff()

    op( hdeg, hpol )

    return hdeg

def get_project( pol_lst, pmat, retry = True ):
    '''
    INPUT:
        - "pol_lst" -- A list of homogeneous polynomials in QQ[x0,...,x8].
        - "pmat"    -- A matrix defined over the rationals QQ.    
    OUTPUT:
        - Returns a 2-tuple of polynomials:
        
            * a homogeneous polynomial F in QQ[x0,x1,x2,x3]. 
            
            * F(1,x,y,z) in QQ[x,y,z] (affine polynomial)
    '''

    Ry = PolynomialRing( GF( 2 ), var( 'y0,y1,y2,y3,y4,y5,y6,y7,y8' ), order = 'degrevlex' )
    v = OrbRing.coerce( '[v0,v1,v2,v3,v4,v5,v6,v7,v8]' )
    x = OrbRing.coerce( '[x0,x1,x2,x3,x4,x5,x6,x7,x8]' )
    vx_dct = {v[i]:x[i] for i in range( 9 )}

    op( "\n" + str( pmat ) )

    tries = 0
    projected = False
    while not projected:

        # obtain the linear equations of the projection map
        pmat = matrix( OrbRing.R, list( pmat ) )
        leq_lst = list( pmat * vector( x ) )

        # compute the image of this projection map
        proj_lst = [ v[i] - leq_lst[i] for i in range( len( leq_lst ) ) ]
        p_lst = ideal( pol_lst + proj_lst ).elimination_ideal( x ).gens()

        # obtain a polynomial in x0,...,x8
        p_lst = [p.subs( vx_dct ) for p in p_lst]
        fx = p_lst[0]

        tries += 1
        if len( fx.variables() ) < 4 or len( p_lst ) != 1:

            pmat = get_pmat( True )

            if tries % 100 == 0:
                op( 'tries =', tries, p_lst )

            if tries == 1000:
                return -1
        else:
            projected = True

    w0, w1, w2, w3 = fx.variables()
    fx = fx.subs( {w0:x[0], w1:x[1], w2:x[2], w3:x[3]} )

    x0, x1, x2, x3 = OrbRing.coerce( 'x0,x1,x2,x3' )
    x, y, z = var( 'x,y,z' )
    fxyz = fx.subs( {x0:1, x1:x, x2:y, x3:z} )

    op( fx )
    op( fxyz )

    return fx, fxyz

def get_factor_lst( pol ):
    '''
    INPUT:
        - "pol" -- A polynomial.
    OUTPUT:
        - Return a list of factors of "pol": 
          [ (<factor>,<multiplicity>),... ] 
    '''
    maple.eval( 'fct := evala(AFactors(' + str( pol ) + '));' )
    fct_lst = maple.eval( 'lprint(fct);' )
    fct_lst = str( fct_lst ).split( ',' )

    op( fct_lst )

    cf = fct_lst[0].replace( '[', '' ).replace( ']', '' ).strip()
    new_lst = []
    for i in range( 1, len( fct_lst ), 2 ):
        fact = fct_lst[i].replace( '[', '' ).replace( ']', '' ).strip()
        mult = fct_lst[i + 1].replace( '[', '' ).replace( ']', '' ).strip()
        new_lst += [( fact, mult )]

    op( cf )
    op( len( new_lst ), new_lst )

    return new_lst

def get_genus( pol, plane = 'x1+2*x2+17*x3+11*x0' ):
    '''
    INPUT:
        - "pol"   -- A polynomial in QQ[x0,x1,x2,x3].
        - "plane" -- A String of a linear polynomial in QQ[x0,x1,x2,x3].
    OUTPUT:
        - An integer denoting the geometric genus of
          the curve given by the intersection of the 
          surface defined by "pol", with the plane defined 
          by "plane". If this geometric genus is not defined
          then -2 is returned.
    '''
    op( pol )

    # obtain an equation for the curve defined
    # intersecting a plane with the zero-set of pol.
    plane = OrbRing.coerce( plane )
    K = ideal( pol, plane ).groebner_basis()
    P = K[0]
    if P.total_degree() <= 1:
        P = K[1]
    P = P.subs( {OrbRing.coerce( 'x3' ):1} )

    # compute geometric genus with Maple.
    maple.eval( 'with(algcurves);' )
    maple.eval( 'P := ' + str( P ) + ';' )
    gen = maple.eval( 'genus(P,x1,x2);' )

    op( gen )

    try:
        return int( gen )
    except ValueError:
        return -2


def get_sing_lst( pol, probable = True ):
    '''
    INPUT:
        - "pol"      -- A homogeneous polynomial in QQ[x0,x1,x2,x3].
        - "probable" -- If True, performs a non-deterministic version of a 
                        radical decomposition algorithm, which is faster.
                        The correctness of output is not guaranteed in this case.
    OUTPUT:
        - Suppose that X is the surface defined by the zero-set of "pol".
          The output is a list
           
            [ (<I>, <H>), ... ]
          
          where <I> is the ideal of a component in singular locus of X, 
          and <H> is the Hilbert polynomial of this ideal.
    '''
    x0, x1, x2, x3 = OrbRing.coerce( 'x0,x1,x2,x3' )
    df_str = str( [diff( pol, x0 ), diff( pol, x1 ), diff( pol, x2 ), diff( pol, x3 )] )[1:-1]

    op( df_str )

    mi = ''
    mi += 'P<x0,x1,x2,x3> := PolynomialRing(RationalField(), 4);\n'
    mi += 'MI := ideal< P |' + df_str + '>;\n'

    if probable:
        mi += 'MD := ProbableRadicalDecomposition( MI );\n'
    else:
        mi += 'MD := RadicalDecomposition( MI );\n'
    mi += '#MD;\n'

    mo1 = int( magma.eval( mi ) )

    sing_lst = []
    Ry = PolynomialRing( OrbRing.num_field, var( 'y0,y1,y2,y3' ), order = 'degrevlex' )
    for idx in range( mo1 ):

        comp = str( magma.eval( 'Basis(MD[' + str( idx + 1 ) + ']);\n' ) )
        comp = comp.replace( '\n', '' )

        # compute hilbert polynomial of component in singular locus
        compy = sage_eval( comp.replace( 'x', 'y' ), Ry.gens_dict() )
        idy = Ry.ideal( compy )
        hpol = idy.hilbert_polynomial()

        sing_lst += [( comp, hpol )]
        op( idx, sing_lst[-1] )

    op( sing_lst )

    return sing_lst


def get_pmz( pmat, omat, vmat ):
    '''
    INPUT:
        - "pmat" -- A 4x9 invertible matrix with entries in QQ. 
                    A matrix represents a projection from 
                    a 7-sphere in projective 8-space  
                    to projective 3-space.
        
        - "omat" -- A 9x9 invertible matrix with entries in QQ[c1,s1]. 
                    This matrix represents a projective curve 
                    in the automorphism group of the projective 7-sphere.
                    Thus "omat" represents a 1-parameter subgroup in Aut(S^7).
                    
        - "vmat" -- A 9x9 invertible matrix with entries in QQ. 
                    This matrix represents an element in Aut(S^7),
                    which transforms a standard circle (see below).
    OUTPUT:
        - Returns two lists of elements in QQ[c0,s0,c1,s1].
          The 1st list has length 9. The 2nd list has length 4.
          
          The 1st list represent a parametrization 
              S^1xS^1--->S^7
          of a surface in S^7 and the 2nd list the 
          parametrization 
              S^1xS^1--->P^3
          of the surface projected into
          P^3 (projective 3-space) using "pmat". 
          Here (c0,s0) and (c1,s1) are points on S^1.
                
          The surface in S^7 is obtained by applying 
          a 1-parameter subgroup to a circle C. 
          Here C is the "vmat"-transform of 
          the standard circle B in S^7 where:
          
                B = { x | -x0^2+x1^2+x2^2==0 } 
              
              S^7 = { x | -x0^2+x1^2+...+x8^2==0 }. 
          
        
    '''
    c1, s1 = OrbRing.coerce( 'c1,s1' )
    pmz_lst = list( omat * vmat * vector( [1, c1, s1, 0, 0, 0, 0, 0, 0] ) )
    prj_pmz_lst = list( pmat * omat * vmat * vector( [1, c1, s1, 0, 0, 0, 0, 0, 0] ) )

    return pmz_lst, prj_pmz_lst


def get_orb_bp_tree( pmz_lst ):
    '''
    INPUT:
        - "pmz_lst" -- A list of 9 elements p0,...,p8 in QQ[c0,s0,c1,s1]
                       such that -p0^2+p1^2+...+p8^2==0.
                       Some of the polynomials can be equal to zero. 
                       The list should represent a paramatrization:
                       
                           S^1xS^1--->S^7.
                       
                       Here (c0,s0) is a points on S^1 such that
                       thus c0^2+s0^2-1==0. Similarly for (c1,s1).
                                                                                 
    OUTPUT:
    
        - Returns "linear_series.BasePointTree" object, 
          which represents the base points of a parametrizing 
          map given by the composition:
          
                P^1xP^1---->S^1xS^1--->S^7--->S^n
          
          with 2<=n<=7. The composition of the latter two maps 
          are defined by omitting the zero polynomials from "pmz_lst". 
    '''

    # setup dictionary for reparametrization'
    #
    c0, s0, c1, s1 = OrbRing.coerce( 'c0,s0,c1,s1' )
    dct1 = {}
    dct1[c0] = '2*t0*t1/(t0^2+t1^2)'
    dct1[s0] = '(-t0^2+t1^2)/(t0^2+t1^2)'
    dct1[c1] = '2*v0*v1/(v0^2+v1^2)'
    dct1[s1] = '(-v0^2+v1^2)/(v0^2+v1^2)'
    for key in dct1: dct1[key] = OrbRing.coerce( dct1[key] )

    # apply reparametrization and multiply out denominators
    # where we only consider non-zero entries
    #
    ps_lst = [ pmz for pmz in pmz_lst if pmz != 0 ]
    gcm1 = OrbRing.coerce( '(t0^2+t1^2)*(v0^2+v1^2)' )
    ps_lst = [ OrbRing.coerce( ps.subs( dct1 ) * gcm1 ) for ps in ps_lst ]

    # ensure that polynomials are co-prime
    #
    gcd1 = gcd( ps_lst )
    ps_lst = [ OrbRing.coerce( ps / gcd1 ) for ps in ps_lst ]
    op( 'gcd =', gcd1 )
    op( 'ps_lst =', ps_lst )

    # Verify whether "ps_lst" represents a map P^1xP^1--->S^n
    # where "n==len(ps_lst)".
    #
    sum1 = sum( [-ps_lst[0] ** 2] + [ ps ** 2 for ps in ps_lst[1:] ] )
    op( 'sum =', sum1 )
    if sum1 != 0:
        warnings.warn( 'Warning: Not parametrization of surface in S^7: ' + str( sum1 ), )

    # set coordinates x,y,v,w
    #
    t0, t1, v0, v1 = OrbRing.coerce( 't0,t1,v0,v1' )
    dct2 = {}
    dct2[t0] = var( 'x' )
    dct2[t1] = var( 'y' )
    dct2[v0] = var( 'v' )
    dct2[v1] = var( 'w' )
    xyvw_lst = [ str( ps.subs( dct2 ) ) for ps in ps_lst ]

    #
    # verbose output of linear_series package
    #
    if op.input_file_name != None:
        dprint( False, 'no output' )
    else:
        dprint( True )

    #
    # Compute base point tree using "linear_series" package
    #
    ls = LinearSeries( xyvw_lst, PolyRing( 'x,y,v,w' ) )
    bp_tree = ls.get_bp_tree()
    op( ls )
    op( bp_tree )

    return bp_tree


def get_imp( omat, vmat ):
    '''
    INPUT:
        - "omat" -- A 9x9 invertible matrix with entries in QQ[c0,s0]. 
                    This matrix represents a projective curve 
                    in the automorphism group of the projective 7-sphere.
                    Thus "omat" represents a 1-parameter subgroup in Aut(S^7).
                    
        - "vmat" -- A 9x9 invertible matrix with entries in QQ. 
                    This matrix represents an element in Aut(S^7),
                    which transforms a standard circle (see OUTPUT below).    

    OUTPUT:
        - Returns a list of elements in QQ[x0,...,x8].
        
          This list represent the generators of the ideal corresponding to 
          the variety, which is obtained by applying a 1-parameter subgroup to
          a circle C. Here C is the "vmat"-transform of 
          the standard circle B in S^7 where:
          
                B = { x | -x0^2+x1^2+x2^2==0 } 
              
              S^7 = { x | -x0^2+x1^2+...+x8^2==0 }.   
    '''

    # declare list of coordinate variables
    v_lst = OrbRing.coerce( '[v0,v1,v2,v3,v4,v5,v6,v7,v8]' )
    x_lst = OrbRing.coerce( '[x0,x1,x2,x3,x4,x5,x6,x7,x8]' )
    c_lst = OrbRing.coerce( '[c0,s0]' )

    # construct generators of ideal of orbital product
    g_lst = []

    # v[i] coordinates of points on a circle C in S^7
    vmatI = vmat.inverse()
    g_lst += list( vmatI * vector( v_lst ) )[3:]
    g_lst += OrbRing.coerce( '[-v0^2+v1^2+v2^2+v3^2+v4^2+v5^2+v6^2+v7^2+v8^2]' )

    # consider all point that are an orbit of C
    # the 1-parameter subgroup of Aut(S^7)
    # represented by the matrix omat
    e_lst = list( omat * vector( v_lst ) )
    g_lst += [ x_lst[i] - e_lst[i] for i in range( 0, 8 + 1 ) ]
    g_lst += OrbRing.coerce( '[-x0^2+x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2]' )
    g_lst += OrbRing.coerce( '[s0^2+c0^2-1]' )

    # compute the resulting variety by elimination
    g_ideal = OrbRing.R.ideal( g_lst )
    imp_lst = list( g_ideal.elimination_ideal( v_lst + c_lst ).gens() )

    return imp_lst


def get_pmz_test( o ):
    '''
    INPUT:
    
        - "o" -- OrbOutput object.
    
    OUTPUT:
        
        - Return:
         
            * "None" if the test could not be performed because
              either parametrization or implicit equation was not available.
            
            * Return "True" if implicit equation of projected surface
              agrees with parametrization composed with projection.
          
            * Return "False" if both parametrization and implicit equation were 
              available but do not agree.               
    '''

    if not o.input.do['pmz'] or not o.input.do['imp']:
        return None

    if o.prj_pol == -1 or o.gen == -2:
        return False

    op( 'Testing parametrization...' )
    c0, s0, c1, s1 = OrbRing.coerce( 'c0,s0,c1,s1' )
    x0, x1, x2, x3 = OrbRing.coerce( 'x0,x1,x2,x3' )
    f = o.prj_pol
    p = o.prj_pmz_lst
    fp = o.prj_pol.subs( {x0:p[0], x1:p[1], x2:p[2], x3:p[3]} )
    test = fp.reduce( [c0 * c0 + s0 * s0 - 1, c1 * c1 + s1 * s1 - 1] )
    op( test )
    if test != 0:
        return False

    return True


def orb_product( input ):
    '''
    For additional info see the documentation of 
    "OrbInput" and "OrbOutput".
    
    INPUT:
        - "input" -- "OrbInput"
    OUTPUT:
        - Returns "OrbOutput" object. 
          
    '''

    o = OrbOutput( input )

    if input.do['pmz']:
        o.pmz_lst, o.prj_pmz_lst = get_pmz( input.pmat, input.omat, input.vmat )

    if input.do['bpt']:
        o.bp_tree = get_orb_bp_tree( o.pmz_lst )

    if input.do['imp']:
        o.imp_lst = get_imp( input.omat, input.vmat )
    else:
        return o  # cannot obtain remaining attributes without "o.imp_lst"

    # Compute remaining attributes
    #
    if input.do['dim']:
        o.dim = get_dim( o.imp_lst )

    if input.do['deg']:
        o.deg = get_degree( o.imp_lst )

    if input.do['prj']:
        o.prj_pol, o.xyz_pol = get_project( o.imp_lst, input.pmat )

    if input.do['fct']:
        o.fct_lst = get_factor_lst( o.prj_pol )

    if input.do['gen']:
        o.gen = get_genus( o.prj_pol )

    if input.do['sng']:
        o.sng_lst = get_sing_lst( o.prj_pol )

    # Test whether parametrization agrees with implicitization.
    #
    if input.do['tst']:
        o.pmz_test = get_pmz_test( o )

    return o






