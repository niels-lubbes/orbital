'''
Created on Aug 8, 2016
@author: Niels Lubbes
'''
from sage.all import *

from class_orb_ring import *
from orb_product import *
from orb_verbose import *

from linear_series.verbose import *
from linear_series.class_linear_series import *
from linear_series.class_base_points import *
from linear_series.class_poly_ring import *
from linear_series.get_implicit import *


def get_S1xS1_pmz( xyvw_pmz_lst ):
    '''
    INPUT:
        - "xyvw_pmz_lst" -- A list of bi-homogeneous polynomials in 
                            QQ[x,y;v,w] that represent a parametric map
                                P^1xP^1 ---> P^n.
    OUTPUT:
        - Returns a list of polynomials in 
            
            QQ[c0,s0,s1,c1] / <c0^2+s0^2-1,c1^2+s1^2-1>
          
          that represents the composition of the input parametric map
          with a birational map:
          
              S^1xS^1       --->       P^1xP^1
        (1:c0:s0)x(1:c1:s1) |--> (1-s0:c0)x(1-s1:c1)            
    '''

    R = PolynomialRing( QQ, 'x,y,v,w' )
    x, y, v, w = R.gens()
    xyvw_pmz_lst = sage_eval( str( xyvw_pmz_lst ), R.gens_dict() )
    c0, s0, c1, s1 = OrbRing.coerce( 'c0,s0,c1,s1' )
    sub_dct = {x:1 - s0, y:c0, v:1 - s1, w:c1 }
    pmz_lst = [ OrbRing.coerce( xyvw_pmz.subs( sub_dct ) ) for xyvw_pmz in xyvw_pmz_lst]

    return pmz_lst


def get_deg_surf( imp_lst, emb_dim ):
    '''
    INPUT:
        - "imp_lst" -- A list of polynomials in QQ[x0,...,xn]
                       defining a surface.
                       where n==emb_dim.
        
        - "emb_dim"  -- A non-negative integer, representing the embedding dimension. 
    OUTPUT:
        
        - Returns  degree of surface defined by "imp_lst".
    '''

    ring = PolynomialRing( QQ, [ 'x' + str( i ) for i in range( emb_dim + 1 )] )
    imp_lst = sage_eval( str( imp_lst ), ring.gens_dict() )
    hpol = ring.ideal( imp_lst ).hilbert_polynomial()

    return hpol.diff().diff()


def get_prj_mat( nrows, ncols, num_tries = 0 ):
    '''
    INPUT:
       
        - "nrows"     -- A non-negative integer.
       
        - "ncols"     -- A non-negative integer st. nrows<=ncols.  
       
        - "num_tries" -- A positive integer.
    
    OUTPUT:
    
        - Returns full rank nrow*ncol matrix over QQ
          with coefficients in [0,1].
          
          If nrows==ncols then the identity matrix 
          is returned.
          
          If "num_tries" is small then some simple 
          predefined matrices are returned.
          
          If "num_tries" is large enough then 
          a the returned matrix is random.          
    '''

    if nrows > ncols:
        raise ValueError( 'Expect nrows<=ncols. (nrows,ncols) =', ( nrows, ncols ) )

    if nrows == ncols:
        return identity_matrix( QQ, nrows )

    predef_lst = []
    for i in range( ncols - nrows ):
        predef_lst += [ identity_matrix( ncols ).submatrix( i, 0, nrows, ncols ) ]

    if num_tries < len( predef_lst ):
        return predef_lst[ num_tries ]

    MS = MatrixSpace( GF( 2 ), nrows, ncols )
    Q = matrix( QQ, list( MS.random_element() ) )
    while Q.rank() != nrows:
        Q = matrix( QQ, list( MS.random_element() ) )

    return Q


def approx_QQ( mat ):
    '''
    INPUT:
        - "mat"         -- A matrix defined over QQbar.
        - "close_exact" -- A boolean. 
    OUTPUT:
        - An approximation of the matrix over QQ.
    '''
    # compute rational approximation of U
    #
    A = []
    for row in list( mat ):
        for col in row:

            tmp = sage_eval( str( col ).replace( '?', '' ) )
            if tmp in QQ:
                A += [ tmp ]
            else:
                A += [ tmp.simplest_rational() ]

    return matrix( QQ, mat.nrows(), mat.ncols(), A )


def rand_surf_prj( ls, prj_dim, prv_Q = None ):
    '''
    Helper method for "get_surf".
    
    INPUT:
        - "ls"      -- LinearSeries on P^1xP^1 such that it 
                       associated map parametrizes a surface S 
                       in projective n-space P^n.
        
        - "prj_dim" -- An integer denoting the projected dimension.
        
        - "prv_Q"   -- Either None or a full rank (m+1)*(n+1) 
                       matrix over the rationals QQ.   
    
    OUTPUT:
        - Returns a 3-tuple
            
            ( <Q>, <pmz_lst>, <imp_lst> )
        
          where
          
            * <Q>: 
              An (m+1)*(n+1) matrix Q over the rationals QQ 
              corresponding to a linear projection of S in P^n to P^m 
              where m=="prj_dim". 
              If prv_Q is not equal None then Q is set to "prv_Q".

            * <pmz_lst>:
              A list of m+1 polynomials in 
                  QQ[c0,s0,c1,s1]/<c0^2+s0^1-1,c1^2+s1^2-1>
              which represent the birational/parametric map 
                  S^1xS^1 ---> Q(S) in P^m.
                                     
            * <imp_lst>:
              A list of polynomials in QQ[x0,x1,...,xm]
              which generate the ideal of the surface 
              Q(S) in P^m.        
    '''

    # compute implicit equations and degree of surface in P^n
    #
    pre_pmz_lst = ls.pol_lst
    pre_imp_lst = ls.get_implicit_image()
    pre_dim = len( ls.pol_lst ) - 1
    pre_deg = get_deg_surf( pre_imp_lst, pre_dim )

    # compute projection X of surface S in P^m where m==sum(sig)-1
    #
    if pre_dim < prj_dim:

        raise ValueError( 'Expect pre_dim < prj_dim: ', pre_dim, prj_dim )

    elif pre_dim == prj_dim:

        # no need for projection to lower dimension
        #
        Q = identity_matrix( QQ, prj_dim + 1, prj_dim + 1 )
        p1p1_pmz_lst = pre_pmz_lst
        imp_lst = pre_imp_lst

    else:  # pre_dim > prj_dim:

        # project to lower dimension while preserving degree
        #
        num_tries = 0
        deg = 0
        while deg != pre_deg:

            if prv_Q == None:
                Q = get_prj_mat( prj_dim + 1, pre_dim + 1, num_tries )
            else:
                Q = matrix( QQ, prv_Q )
            p1p1_pmz_lst = list( Q * vector( pre_pmz_lst ) )
            imp_lst = LinearSeries( p1p1_pmz_lst, ls.ring ).get_implicit_image()
            deg = get_deg_surf( imp_lst, len( p1p1_pmz_lst ) - 1 )

            num_tries += 1

    # set pmz_lst and imp_lst
    #
    pmz_lst = get_S1xS1_pmz( p1p1_pmz_lst )
    imp_lst = OrbRing.coerce( imp_lst )

    return Q, pmz_lst, imp_lst


def get_surf( ls, sig, coef_lst = None, prv_Q = None ):
    '''
    Constructs a surface S in a quadric with given signature. This 
    surface is a projection of the surface that is parametrized
    by a map associated to a given linear series.
    
    INPUT:
        
        - "ls"       -- LinearSeries on P^1xP^1 such that it 
                        associated map parametrizes a surface S 
                        in projective n-space P^n.
                    
        - "sig"      -- A 2-tuple of positive integers (a,b). 
        
        - "coef_lst" -- Either None or a list of coefficients 
                        over QQ of length equal to the number 
                        of generators for the ideal of S.                        
                                
        - "prv_Q"    -- Either None or an full rank (m+1)*(n+1) 
                        matrix over the rationals QQ.                                            
        
        
    OUTPUT:
        - Returns a dictionary "dct" where:
          
            * dct['Q']: 
              A (m+1)*(n+1) matrix over the rationals QQ corresponding to 
              a linear projection of S in P^n to P^m where m==sum(sig)-1. 
              Thus a quadric with signature "sig" is a hyper-quadric 
              in projective m-space. If "prv_Q" is not equal None then
              Q is set to "prv_Q".

            * dct['pmz_lst']
              A list of m+1 polynomials in 
              QQ[c0,s0,c1,s1]/<c0^2+s0^1-1,c1^2+s1^2-1>
              which represent the birational/parametric map 
                  S^1xS^1 ---> Q(S) < P^m.
                                     
            * dct['imp_lst']
              A list of polynomials in QQ[x0,x1,...,xm]
              which represent the ideal of the surface Q(S) in P^m.
                  
            * dct['M']
              The symmetric matrix M associated to an hyperquadric 
              that contains the surface image Q(S) under the linear
              transformation Q, or None if no such hyperquadric exists.
              If "coef_lst" is None then M is computed
              randomly, which is the default. 
              The ideal for Q(S) is internally 
              represented by an ordered list of polynomials 
              in Q[x0,...,xn]. The sum of generators with 
              given coefficients in "coef_lst" should be 
              the equation of a quadric with signature "sig".              
            
            * dct['UJ']
              A tuple (U,J) such that M == U.T*J*U and J is a diagonal 
              matrix with entries in [0,1,-1] and signature "sig" 
              (ie. an orthogonal diagonalization of M).
              The matrices U and J are defined over QQbar. 
              The diagonal of J contains first all -1 and then 
              all +1 values (thus the order of the 2-tuple "sig" 
              does matter).
    '''

    sig_set = set( [] )
    cur_sig = ()
    op( 'Computing random quadrics in ideal...' )

    while cur_sig != sig:

        # compute random projection
        #
        Q, pmz_lst, imp_lst = rand_surf_prj( ls, sum( sig ) - 1, prv_Q )

        # set coefficient list
        #
        if coef_lst != None:
            c_lst = coef_lst
        else:
            c_lst = []
            lst = [-1, 0, 1]
            for imp in imp_lst:
                idx = int( ZZ.random_element( 0, len( lst ) ) )
                c_lst += [ lst[idx] ]

        # obtain quadric in ideal from "c_lst"
        #
        i_lst = range( len( imp_lst ) )
        M_pol = [ c_lst[i] * imp_lst[i] for i in i_lst if imp_lst[i].total_degree() == 2 ]
        M_pol = sum( M_pol )

        # eigendecomposition (M*V == V*D)
        #
        ring = PolynomialRing( QQ, [ 'x' + str( i ) for i in range( sum( sig ) )] )
        x_lst = ring.gens()
        M = invariant_theory.quadratic_form( M_pol, x_lst ).as_QuadraticForm().matrix()
        M = matrix( QQ, M )
        D, V = M.eigenmatrix_right()  # D has first all negative values on diagonal

        # determine signature of quadric
        #
        num_neg = len( [ d for d in D.diagonal() if d < 0 ] )
        num_pos = len( [ d for d in D.diagonal() if d > 0 ] )
        cur_sig = ( num_neg, num_pos )
        if cur_sig not in sig_set:
            sig_set.add( cur_sig )
            op( '\tsig =', sig, ', sig_set =', sig_set )

    op( 'Q        =', list( Q ) )
    op( 'pmz_lst  =', pmz_lst )
    op( 'imp_lst  =', imp_lst )
    op( 'c_lst    =', c_lst )
    op( 'M_pol    =', M_pol )

    # diagonal orthonormalization
    #
    # M == ~V*D*V == W.T*D*W == U.T*J*U
    # where U == L*W and D == L.T*J*L
    #
    D, V = M.eigenmatrix_right()  # M*V==V*D
    W = matrix( [col / col.norm() for col in V.columns()] )
    L = diagonal_matrix( [ d.abs().sqrt() for d in D.diagonal()] )
    J = diagonal_matrix( [ d / abs( d ) for d in D.diagonal()] )
    U = L * W

    op( 'M        =', list( M ) )
    op( 'U        =', list( U ) )
    op( 'J diag.  =', J.diagonal() )

    # Return output
    #
    dct = {}
    dct['Q'] = Q
    dct['pmz_lst'] = pmz_lst
    dct['imp_lst'] = imp_lst
    dct['M'] = M
    dct['UJ'] = ( U, J )

    return dct


def test_get_surf( dct ):
    '''
    INPUT:
        - "dct" -- The dictionary which is the output of "get_surf()".
    OUTPUT:
        - Returns a dictionary "test_dct" with the keys "dct.keys()" 
          and the key 'all'.
          The values are True if it passed the test and False otherwise.
          If test_dct['all'] is True then all test were passed successfully. 
    '''
    x_var = ','.join( ['x' + str( i ) for i in range( dct['M'].nrows() )] )
    x_lst = OrbRing.coerce( x_var )
    x_vec = vector( OrbRing.R, x_lst )

    Q = dct['Q']
    M = dct['M']
    U, J = dct['UJ']
    imp_lst = dct['imp_lst']
    pmz_lst = dct['pmz_lst']


    pmz_dct = { x_lst[i]:pmz_lst[i] for i in range( len( pmz_lst ) )  }
    sub_lst = [ imp.subs( pmz_dct ) for imp in imp_lst ]

    test_dct = {}
    test_dct['Q'] = Q.rank() == Q.nrows()
    test_dct['M'] = x_vec * M * x_vec in ideal( imp_lst )
    test_dct['UJ'] = M == approx_QQ( U.T * J * U )
    test_dct['imp_lst'] = len( imp_lst ) > 0
    test_dct['pmz_lst'] = set( sub_lst ) == set( [0] )

    test_all = True
    for key in test_dct:
        test_all = test_all and test_dct[key]
    test_dct['all'] = test_all

    return test_dct


def get_proj( imp_lst, pmz_lst, mat ):
    '''
    INPUT:        
        - "imp_lst" -- A list of polynomials in QQ[x0,...,xn]
                       for n<9.
        
        - "pmz_lst" -- A list n+1 of polynomials in
                           QQ[c0,s0,c1,s1]/<c0^2+s0^2-1,c1^2+s1^2-1>.

        - "mat"     -- A 4*n matrix defined over QQ and of full rank.
    
    OUTPUT:
        - Returns ( f_xyz, prj_pmz_lst )
          where 
              
              * "f_xyz" is a polynomial in QQ[x,y,z,w] whose zeroset 
                is the linear projection of the zeroset of "imp_lst".
                The linear projection is defined by "mat". 
              
              * "prj_pmz_lst" is the linear projection of the vector
                 defined by "pmz_lst". The linear projection is 
                 defined by "mat".
    '''
    # init lists of polynomial ring generators
    #
    v = OrbRing.coerce( '[v0,v1,v2,v3,v4,v5,v6,v7,v8]' )
    x = OrbRing.coerce( '[x0,x1,x2,x3,x4,x5,x6,x7,x8]' )

    # obtain the linear equations of the projection map
    #
    leq_lst = list( mat * vector( x[:mat.ncols()] ) )

    # compute the image of this projection map
    #
    proj_lst = [ v[i] - leq_lst[i] for i in range( len( leq_lst ) ) ]
    p_lst = ideal( imp_lst + proj_lst ).elimination_ideal( x ).gens()

    # subs (v0:v1:v2:v3) => (1:x:y:z)
    #
    v_xyz_dct = { v[0]:1, v[1]:var( 'x' ), v[2]:var( 'y' ), v[3]:var( 'z' ) }
    p_lst = [p.subs( v_xyz_dct ) for p in p_lst]
    if len( p_lst ) != 1:
        raise ValueError( 'Expect projection to be a surface: p_lst =', p_lst )
    f_xyz = p_lst[0]

    # project parametrization
    #
    prj_pmz_lst = list( mat * vector( pmz_lst ) )

    return f_xyz, prj_pmz_lst





