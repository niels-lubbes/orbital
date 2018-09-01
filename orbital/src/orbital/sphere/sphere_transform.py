'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 22, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_identity_matrix
from orbital.sage_interface import sage_QQ

from orbital.cossin.cos_sin import get_cs


def get_rot_mat( dim, i, j, angle ):
    '''    
    Constructs a higher dimensional rotation matrix.
        
    Parameters
    ----------
    dim: int
        Width and height of square matrix.        
    i: int
    j: int
    angle: int
        Angle in [0,360].
        
    Returns
    -------
    sage_matrix
        A matrix which looks like a square identity matrix of
        width "dim", except the entries 
            (i,i), (i,j), 
            (j,i), (j,j)
        define the rotation matrix:
            [ cos(angle) -sin(angle) ]
            [ sin(angle)  cos(angle) ]
    '''

    c, s = get_cs( angle )


    mat = []
    for a in range( dim ):
        row = []
        for b in range( dim ):
            e = 0
            if a == b: e = 1
            if a != b: e = 0
            if ( a, b ) == ( i, i ): e = c
            if ( a, b ) == ( i, j ): e = -s
            if ( a, b ) == ( j, i ): e = s
            if ( a, b ) == ( j, j ): e = c
            row += [e]
        mat += [row]

    return sage_matrix( mat )


def get_rot_S3( rot ):
    '''
    Constructs a rotation matrix for S^3.
        
    Parameters
    ----------
    rot : list<int>
        List of 6 rotation angles. Each angle is an integer in [0,360].
    
    Returns
    -------
    sage_matrix
        A 5x5 matrix that rotates the projective closure of S^3 in 
        projective 4-space along the angle parameters.         
    '''

    a12, a13, a23, a14, a24, a34 = rot

    M = sage_identity_matrix( 5 )
    M *= get_rot_mat( 5, 1, 2, a12 )
    M *= get_rot_mat( 5, 1, 3, a13 )
    M *= get_rot_mat( 5, 2, 3, a23 )
    M *= get_rot_mat( 5, 1, 4, a14 )
    M *= get_rot_mat( 5, 2, 4, a24 )
    M *= get_rot_mat( 5, 3, 4, a34 )

    return sage_matrix( M )


def get_trn_S3( trn ):
    '''
    Constructs a translation matrix for S3 that preserves (1:0:0:0:1).
    Via a stereographic projection the transformation 
    correspond to a translation in Euclidean 3-space.
    
    Parameters
    ----------
    trn: list<sage_QQ>
        List of three rationals tx, ty, tz representing 
        translation parameter in x-, y- and z- direction.
      
    Returns
    -------
    sage_matrix
        Transformation matrix that preserves the projective closure
        of S^3 in P^4 and defines a translation in Euclidean 3-space.    
    '''

    tx, ty, tz = trn

    T = ( sage_QQ( 1 ) / 2 ) * ( tx ** 2 + ty ** 2 + tz ** 2 )

    M = []
    M += [[T + 1, tx, ty, tz, -T]]
    M += [[tx, 1, 0, 0, -tx]]
    M += [[ty, 0, 1, 0, -ty]]
    M += [[tz, 0, 0, 1, -tz]]
    M += [[T, tx, ty, tz, -T + 1]]

    return sage_matrix( M )


def get_scale_S3( s ):
    '''
    Constructs a scaling matrix for S3
    that preserves (1:0:0:0:1).
    Via the stereographic projection the transformation 
    correspond to a scalings in Euclidean 3-space.
    
    Parameters
    ----------
    s: sage_QQ
        Scaling factor.
      
    Returns
    -------
    sage_matrix
        Transformation matrix that preserves the projective closure
        of S^3 in P^4 and defines a scaling in Euclidean 3-space.    
    '''

    M = []
    M += [[s ** 2 + 1, 0, 0, 0, s ** 2 - 1]]
    M += [[0, 2 * s, 0, 0, 0]]
    M += [[0, 0, 2 * s, 0, 0]]
    M += [[0, 0, 0, 2 * s, 0]]
    M += [[s ** 2 - 1, 0, 0, 0, s ** 2 + 1]]

    return sage_matrix( M )



def get_hp_P4( a, b ):
    '''
    Computes the Hamiltonian product of two vectors 
    in the projective closure P^4 of the quaternions. 
    
    Parameters
    ----------
    a: list  
        List of length 5. 
    b: list
        List of length 5.
        
    Returns
    -------
    list
        List of length 5.
    '''

    a0, a1, a2, a3, a4 = a
    b0, b1, b2, b3, b4 = b

    lst = [ a0 * b0,
            a1 * b1 - a2 * b2 - a3 * b3 - a4 * b4,
            a1 * b2 + a2 * b1 + a3 * b4 - a4 * b3,
            a1 * b3 + a3 * b1 + a4 * b2 - a2 * b4,
            a1 * b4 + a4 * b1 + a2 * b3 - a3 * b2 ]

    return lst


