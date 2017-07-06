'''
Created on Aug 8, 2016
@author: Niels Lubbes
'''

from class_orb_ring import *
from orb_verbose import *
from cos_sin import *



def get_mat( A_str, B_str, C_str ):
    '''
    INPUT:
        - "A_str" -- A string of either one of the 
                     following 6 forms:
                            
                        'O****'
                        'tT'                  
                        'T[#,#,#,#,#,#,#]'                   
                        'R****[%,%,%,%]'                   
                        'I'                   
                        'P1'                  
                        'P0'
                        'M<list of a 9x9 matrix>'
                        'E[$,$,$,$,$,$,$,$]'
                        'E'
                        'X[%,%,%]'

                     where: 
                         
                         *-symbol is a place holder for
                           one of the following characters: 
                           r,s,m,p,a                           
                         
                         #-symbol denotes an element of "OrbRing.num_field".                     
                         
                         %-symbol denotes an integer in [0,360].
                         
                         $-symbol denotes an integer in [1,8] such that each
                           integer occurs only once. For example 'E[2,1,3,4,5,6,7,8]'
                           is correct, but 'E[1,1,3,4,5,6,7,8]' is incorrect.
                     
        - "B_str" -- Same specs as "A_str".
        - "C_str" -- Same specs as "A_str".  
    OUTPUT:
        - Let matrix A over QQ[c0,s0,c1,s1] be defined as follows:  
              
             "A_str[0] == 'O' " : "A = get_omat( A_str )".
             "A_str[0] == 'T' " : "A = get_tmat( A_str )".
             "A_str[0] == 'E' " : "A = get_emat( A_str )".
             "A_str[0] == 'X' " : "A = get_xmat( A_str )".  
             "A_str    == 'tT'" : "A = get_tmat('tT')".
             "A_str[0] == 'R' " : "A = get_rmat( A_str )".
             "A_str    == 'I' " : "A = identity_matrix(9,9)".      
             "A_str    == 'P1'" : "A = get_pmat( True )".
             "A_str    == 'P0'" : "A = get_pmat( False )".
             "A_str[0] == 'M' " : "A = matrix(<list of a matrix>)".
          
          Similarly, we obtain matrices B and C.
          
          The output is the matrix A*B*C.
    '''

    M_lst = []
    for M_str in [A_str, B_str, C_str]:

        if M_str[0] == 'O': M_lst += [get_omat( M_str )]
        elif M_str[0] == 'T': M_lst += [get_tmat( M_str )]
        elif M_str[0] == 'E': M_lst += [get_emat( M_str )]
        elif M_str[0] == 'X': M_lst += [get_xmat( M_str )]
        elif M_str == 'tT': M_lst += [get_tmat( 'tT' )]
        elif M_str[0] == 'R': M_lst += [get_rmat( M_str )]
        elif M_str == 'I': M_lst += [identity_matrix( OrbRing.num_field, 9, 9 )]
        elif M_str == 'P1': M_lst += [get_pmat( True )]
        elif M_str == 'P0': M_lst += [get_pmat( False )]
        elif M_str[0] == 'M':
            mat_lst = OrbRing.coerce( M_str[1:] )
            M_lst += [matrix( mat_lst )]


    return M_lst[0] * M_lst[1] * M_lst[2]


def get_tmat( t_str = None ):
    '''
    INPUT:
        - "t_str" -- A String with either one of the 
                     following 3 formats:
                     
                     * A string with format:
                         'T[#,#,#,#,#,#,#]'
                       where # are in "OrbRing.num_field".
                     
                     * 'tT'.
                     
                     * "None".
    OUTPUT:
        - A 9x9 matrix defined over "OrbRing.num_field" 
          representing an Euclidean translation of S^7. 
          This map is obtained as the composition 
          of a stereographic projection with center 
          (1:0:0:0:0:0:0:0:1), 
          an Euclidean translation in R^7, and
          the inverse stereographic projection.
          
          If "t_str==None" then the entries of the 
          translation matrix are indeterminates  
          t1,...,t7 in "OrbRing.R".
          
          If "t_str=='tT'" then the indeterminates
          are set to [c0,s0,0,0,0,0,0]. 
          Thus the translations along a circle.
    '''

    t = OrbRing.coerce( '[t1,t2,t3,t4,t5,t6,t7]' )

    # construct a translation matrix with undetermined
    # translations in t1,...,t7
    a = ( QQ( 1 ) / 2 ) * sum( [ ti ** 2 for ti in t] )
    mat = []
    mat += [[ 1 + a ] + list( t ) + [-a]]
    for i in range( 0, 7 ):
        mat += [[t[i]] + identity_matrix( OrbRing.R, 7 ).row( i ).list() + [-t[i]]]
    mat += [[ a ] + list( t ) + [1 - a]]

    if t_str == None:
        # return matrix with indeterminates
        return  matrix( OrbRing.R, mat )

    elif t_str == 'tT':
        # translations along a circle
        c0, s0 = OrbRing.coerce( 'c0,s0' )
        q = [c0, s0, 0, 0, 0, 0, 0]

    else:
        # Substitute coordinates [#,#,#,#,#,#,#]
        # for the t0,...,t7 variables.
        q = sage_eval( t_str[1:] )
        if len( q ) != 7:
            raise ValueError( 'Expect 7 translation coordinates: ', t_str )

    # substitute q for t
    smat = []
    for row in mat:
        srow = []
        for col in row:
            srow += [ col.subs( {t[i]:q[i] for i in range( 7 )} )]
        smat += [srow]

    if t_str == 'tT':
        return matrix( OrbRing.R, smat )
    else:
        return matrix( OrbRing.num_field, smat )


def get_omat( o_str ):
    '''
    INPUT:
        - "o_str" -- A string with format:
                        
                        'O****'
                    
                    where: 
                         
                         *-symbol is a place holder for
                           one of the following characters: 
                           r,s,m,p,a  
                    
    OUTPUT:
        - Return a 9x9 matrix over "OrbRing.num_field"
          of the shape
          
              1  0  0  0  0  0  0  0  0 
              0  *  *  0  0  0  0  0  0 
              0  *  *  0  0  0  0  0  0  
              0  0  0  *  *  0  0  0  0  
              0  0  0  *  *  0  0  0  0  
              0  0  0  0  0  *  *  0  0  
              0  0  0  0  0  *  *  0  0  
              0  0  0  0  0  0  0  *  *  
              0  0  0  0  0  0  0  *  *
          
          where the each 2x2 matrix on the diagonal is of eihter 
          one of the following shapes:
          
            'r': c0 -s0    's': -c0  s0      
                 s0  c0         -s0 -c0

            'p': 1  0    'm': -1  0    'a': 1  0
                 0  1          0 -1         0 -1
    
    '''
    # parse input
    if o_str[0] != 'O' or len( o_str ) != 5:
        raise ValueError( 'Incorrect input string: ', o_str )


    br = [[c0, -s0], [s0, c0]]
    bs = [[-c0, s0], [-s0, -c0]]
    bp = [[1, 0], [0, 1]]
    bm = [[-1, 0], [0, -1]]
    ba = [[1, 0], [0, -1]]
    b_dct = { 'r': br, 's': bs, 'p': bp, 'm': bm, 'a':ba }

    bmat_lst = []
    for ch in o_str[1:]:
        bmat_lst += [b_dct.get( ch, 'error' )]

    omat = identity_matrix( OrbRing.R , 9 )
    idx = 1
    for bmat in bmat_lst:
        omat.set_block( idx, idx, matrix( OrbRing.R, bmat ) )
        idx = idx + 2

    return omat

def get_rmat( r_str ):
    '''
    INPUT:
        - "r_str" -- A string with the following format
                                             
                        'R****[%,%,%,%]'                   

                     where: 
                         
                         *-symbol is a place holder for
                           one of the following characters: 
                           r,s,m,p,a                                                                   
                         %-symbol denotes an integer in [0,360]
        
    OUTPUT:
        - Return a 9x9 matrix over "OrbRing.num_field"
          of the same shape as 
              
              "get_omat( 'O****' )"
          
          The "c0" and "s0" in each block are substituted:
          
              c0:=cos( get_cs(%)[0] )
              s0:=sin( get_cs(%)[1] )
                  
          where get_cs(%) returns a rational approximation
          of the sin and cos of the angle %.
    '''

    omat = get_omat( 'O' + r_str[1:5] )
    cs_lst = sage_eval( r_str[5:] )
    c0, s0 = OrbRing.coerce( 'c0,s0' )

    rmat = []
    idx = -1
    dct = {}
    for row in list( omat ):

        if idx % 2 == 0:
            cv, sv = get_cs( cs_lst[ idx / 2 ] )
            dct = { c0:cv, s0:sv}
            if cv ** 2 + sv ** 2 != 1:
               raise Exception( 'Expect cos(%)^2+sin(%)^2=1.' )

        idx = idx + 1

        rrow = []
        for col in row:
            rrow += [col.subs( dct )]
        rmat += [rrow]

    return matrix( rmat )


def get_pmat( random = False ):
    '''
    INPUT:
        - "random" -- A boolean.
    OUTPUT:
        - A full rank 9*4 matrix over QQ. 
          
          If "random" is True then a random 9x4 matrix over GF(2)
          is computed, which is then coerced to a matrix over QQ.
          
          If "random" is False then the returned matrix corresponds 
          to the linear transformation:
          (x0:...:x8) |-> (x0-x8:x1:x2:x3)  
    '''

    if random:
        MS = MatrixSpace( GF( 2 ), 4, 9 )
        pmat = MS.random_element()
        while pmat.rank() != 4:
            pmat = MS.random_element()
        pmat = matrix( QQ, list( pmat ) )
    else:
        pmat = []
        pmat += [[1, 0, 0, 0] + 4 * [0] + [ -1 ] ]
        pmat += [[0, 1, 0, 0] + 5 * [0]]
        pmat += [[0, 0, 1, 0] + 5 * [0]]
        pmat += [[0, 0, 0, 1] + 5 * [0]]
        pmat = matrix( QQ, pmat )

    return pmat


def get_emat( perm_str ):
    '''
    INPUT:
        - Either 'E' or a string with format:
            'E[$,$,$,$,$,$,$,$]'
          where each integer $ in range(1,8+1) occurs
          exactly once. For example 'E[2,1,3,4,5,6,7,8]'
          is correct, but 'E[1,1,3,4,5,6,7,8]' is incorrect.
          
    OUTPUT:
        - A 9*9 permutation matrix over QQ 
          such that the first row an column are: 
          (1,0,0,0,0,0,0,0,0).
          
          If "perm_str=='E'" then a random permutation
          matrix is computed (the first coordinate is not permuted).
          
          Otherwise the returned matrix corresponds 
          to the permutation:
              [0,1,2,3,4,5,6,7,8] ---> [0]+[$,$,$,$,$,$,$,$]
    '''


    #
    # Permutations only permute list of the form [1,2,...]
    # so we add one to the elements of the input list
    #
    if perm_str[0] != 'E':
        raise ValueError( 'Incorrect perm_str: ', perm_str )
    elif perm_str == 'E':
        p = Permutations( range( 1, 8 + 1 ) )
        perm_lst = list( p.random_element() )
    else:
        perm_lst = sage_eval( perm_str[1:] )
        if len( perm_lst ) != 8:
            raise ValueError( 'Incorrect perm_str (8 integers expected): ', perm_str )

    perm_lst = [1] + [perm + 1 for perm in perm_lst ]

    return Permutation( perm_lst ).to_matrix().T


def get_xmat( r_str ):
    '''
    INPUT:
        - r_str -- A string of the following format:
        
                        'X[%,%,%]'
                    
                    where % denotes an integer in [0,360].
                    
    OUTPUT:
        - A matrix corresponding to the map S^7--->S^7 
          that factors via the following rotations matrices 
          subsequently:
          
            Angle 1:
                [  1   0   0   0   0   0   0   0   0]
                [  0  c0 -s0   0   0   0   0   0   0]
                [  0  s0  c0   0   0   0   0   0   0]
                [  0   0   0   1   0   0   0   0   0]
                [  0   0   0   0   1   0   0   0   0]
                [  0   0   0   0   0   1   0   0   0]
                [  0   0   0   0   0   0   1   0   0]
                [  0   0   0   0   0   0   0   1   0]
                [  0   0   0   0   0   0   0   0   1] 

            Angle 2:
                [  1   0   0   0   0   0   0   0   0]
                [  0  c0   0 -s0   0   0   0   0   0]
                [  0   0   1   0   0   0   0   0   0]
                [  0  s0   0  c0   0   0   0   0   0]
                [  0   0   0   0   1   0   0   0   0]
                [  0   0   0   0   0   1   0   0   0]
                [  0   0   0   0   0   0   1   0   0]
                [  0   0   0   0   0   0   0   1   0]
                [  0   0   0   0   0   0   0   0   1] 
            
            Angle 3:
                [  1   0   0   0   0   0   0   0   0]
                [  0   1   0   0   0   0   0   0   0]
                [  0   0  c0 -s0   0   0   0   0   0]
                [  0   0  s0  c0   0   0   0   0   0]
                [  0   0   0   0   1   0   0   0   0]
                [  0   0   0   0   0   1   0   0   0]
                [  0   0   0   0   0   0   1   0   0]
                [  0   0   0   0   0   0   0   1   0]
                [  0   0   0   0   0   0   0   0   1]                 
    '''

    # parse input string
    if r_str[0] != 'X':
        raise ValueError( 'Incorrect format of input: ', r_str )
    r_lst = sage_eval( r_str[1:] )
    if len( r_lst ) != 3:
        raise ValueError( 'Incorrect format of input: ', r_str )

    # for debug purposes
    DBG = False

    # mat1
    a_lst = [r_lst[0], 0, 0, 0]
    mat1 = get_rmat( 'Rrppp' + str( a_lst ) )
    if DBG:
        op( 'mat1 =\n' + str( get_omat( 'Orppp' ) ) )

    # mat2
    a_lst = [r_lst[1], 0, 0, 0]
    e_lst = [1, 3, 2, 4, 5, 6, 7, 8]
    eI_lst = Permutation( e_lst ).inverse()
    m_tup = ( 'E' + str( eI_lst ), 'Rrppp' + str( a_lst ), 'E' + str( e_lst ) )
    mat2 = get_mat( *m_tup )
    if DBG:
        m_tup = ( 'E' + str( eI_lst ), 'Orppp', 'E' + str( e_lst ) )
        op( 'mat2 =\n' + str( get_mat( *m_tup ) ) )

    # mat3
    a_lst = [0, r_lst[2], 0, 0]
    e_lst = [1, 4, 2, 3, 5, 6, 7, 8]
    eI_lst = Permutation( e_lst ).inverse()
    m_tup = ( 'E' + str( eI_lst ), 'Rprpp' + str( a_lst ), 'E' + str( e_lst ) )
    mat3 = get_mat( *m_tup )
    if DBG:
        m_tup = ( 'E' + str( eI_lst ), 'Oprpp', 'E' + str( e_lst ) )
        op( 'mat3 =\n' + str( get_mat( *m_tup ) ) )

    return mat3 * mat2 * mat1


