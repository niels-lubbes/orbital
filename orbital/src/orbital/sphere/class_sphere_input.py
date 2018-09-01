'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Aug 30, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage__eval
from orbital.sage_interface import sage_flatten
from orbital.sage_interface import sage_Color

from orbital.class_orb_ring import OrbRing


class SphereInput:
    '''
    This value object represents two circles CA and CB in the 3-sphere
    and is input for sphere_experiment.clifford().

    Notice that the projective 3-sphere S^3 is a projective compactification of 
    Euclidean 3-space R^3. Translations, rotations and scalings of S^3 act on R^3
    as expected. We remark that S^3 is also a projective compactification of the 
    unit quaternions.
    
    The input defines two circles in S^3 that are obtained via
    rotations, scalings and translations of the standard circle
    [1,cos(t),sin(t),0,0].
    We denote the linear transformations for the two circles CA and CB 
    by A and B resp.
    
    We denote by CA@CB the surface (or curve) resulting from taking the pointwise
    Hamiltonian product of the circles CA and CB.   
    We denote by S is the stereographic projection of the Hamiltonian product CA@CB. 
    
    
    Attributes
    ----------
    rota : list<sage_QQ>
        List of of the form [a12, a13, a23, a14, a24, a34] defining
        rotations by A. 
        
    trna : list<sage_QQ>
        List of of the form [tx, ty, tz] defining translations by A. 
    
    sa : int
        Scaling factor by A.

    rotb : list<sage_QQ>
        List of of the form [a12, a13, a23, a14, a24, a34] defining
        rotations by B. 
        
    trnb : list<sage_QQ>
        List of of the form [tx, ty, tz] defining translations by B. 
    
    sb : int
        Scaling factor by B.
    
    prj : int
        Choice for stereographic projection S^3--->P^3: 
        0: (x0:x1:x2:x3:x4) |--> (x0-x4:x1:x2:x3)
        1: (x0:x1:x2:x3:x4) |--> (x0-x1:x4:x2:x3)
        2: (x0:x1:x2:x3:x4) |--> (x0-x2:x1:x4:x3)
        3: (x0:x1:x2:x3:x4) |--> (x0-x3:x1:x2:x4)
                        
    imp : boolean 
        If True, compute and plot implicit equation of S.                        
    
    sng : boolean 
        If True, compute singular locus of S (uses Magma).
    
    snp : boolean 
        If True, compute (if sng is True) the singular locus of S probablistic.

    pmz : boolean.
        If True, plot parametric surface.
            
    bas : boolean.
        If True, plot base circles.

    mrk : boolean.
        If True, plot points (0,0,0) and (1,0,0).
    
    fam : boolean
        If True, plot families.
    
    famA : boolean
        If True and fam==True, then plot family A.

    famB : boolean
        If True and fam==True, then plot family B.

    famt : int
        Thickness of the circles in family A and B.

    ppt : int
        Number of plot_points in plots. Higher values may lead to higher quality.

    opa : sage_QQ
        Rational number between 0 and 1, denoting the opacity of the surface plot.
        
    rng : float
        Range for implicit plot.
             
    stp : int
        A number in [0,360] that denotes the stepsize between subsequent curves in family
        
    col_pmz : string
        Color for parametric surface plot.
        
    col_imp : string
        Color for implicit surface plot.
        
    col_famA : string
        Color for circles in family A.
        
    col_famB : string
        Color for circles in family B.   

    Notes
    -----
    Interesting examples:
    
    perseus = '[[(0, 0, 0), (0, 0, 0), (0, 0, -2), 1], [(0, 0, 0), (0, 0, 0), (0, 0, 1), 1]]'
    CH1 = '[[(0, 0, 0), (0, 0, 0), (-1, 0, 0), 1], [(0, 0, 0), (0, 0, 0), (1, 0, 0), 1]]'
    great8a = '[[(0, 0, 0), (0, 0, 0), (0, 0, 0), 1], [(0, 0, 0), (0, 0, 0), (3/2, 0, 0), 1]]'
    great8b = '[[(0, 0, 0), (0, 0, 0), (0, 0, 0), 1], [(55, 30, 0), (0,65, 0), (3/2, 0, 0), 1]]'

    Use the following command to set corresponding attributes:

    sinp = SphereInput().set(great8)
    '''

    def __init__( self ):

        # circle A
        self.rota = 0
        self.trna = 0
        self.sa = 0

        # circle B
        self.rotb = 0
        self.trnb = 0
        self.sb = 0

        # choice for stereographic projection
        self.prj = 0

        # settings for computing implicit equation
        self.imp = False  # compute implicit equation of S
        self.sng = False  # compute singular locus of S
        self.snp = False  # compute singular locus probablistic

        # parameters for plots
        self.pmz = True  # plot parametric surface
        self.bas = True  # plot base circles
        self.mrk = True  # plot base points
        self.fam = False  # plot families
        self.famA = True
        self.famB = True
        self.famt = 2  # thickness
        self.ppt = 50  # plot_points
        self.opa = 1  # opacity
        self.rng = 2  # range
        self.stp = 6  # stepsize

        # colors
        self.col_pmz = sage_Color( "orange" )
        self.col_imp = sage_Color( "#ffefb0" )
        self.col_famA = sage_Color( "red" )
        self.col_famB = sage_Color( "blue" )


    def random( self, bnd = 10 ):
        '''
        Sets attributes to random values. 
                
        Parameters
        ----------
        coef_bnd : int    
            Positive integer.
                             
        Returns
        -------
        self
            Sets self.rota, self.trna, self.sa and
            self.rotb, self.trnb, self.sb with random
            values.
            
        '''
        # OrbRing.random_int( coef_size )

        q4 = sage_QQ( 1 ) / 4

        s_lst = [ i * q4 for i in range( 4 * bnd ) ]
        s_lst += [ -s for s in s_lst ]

        # circle A
        self.rota = [ OrbRing.random_elt( range( 360 ) ) for i in range( 6 )]
        self.trna = [ OrbRing.random_elt( range( 360 ) ) for i in range( 3 )]
        self.sa = OrbRing.random_elt( s_lst )

        # circle B
        self.rotb = [ OrbRing.random_elt( range( 360 ) ) for i in range( 6 )]
        self.trnb = [ OrbRing.random_elt( range( 360 ) ) for i in range( 3 )]
        self.sb = OrbRing.random_elt( s_lst )

        return self


    def set( self, short_input ):
        '''
        Convenience method for setting attributes of self.
            
        Parameters
        ----------
        short_input : object
            A list or string of a list that has the following form
            [
                [( a12, a13, a23 ), ( a14, a24, a34 ), ( u1, u2, u3 ), sa ],
                 [( b12, b13, b23 ), ( b14, b24, b34 ), ( v1, v2, v3 ), sb ]
            ]
            See also SphereInput.__str__().
        
        Returns
        -------
        self
            The following attributes of self are set: 
            rota, trna, sa, 
            rotb, trnb, sb.
        
        '''

        lsta, lstb = sage__eval( str( short_input ) )
        a12, a13, a23, a14, a24, a34, u1, u2, u3, sa = sage_flatten( lsta )
        b12, b13, b23, b14, b24, b34, v1, v2, v3, sb = sage_flatten( lstb )

        # circle A
        self.rota = [a12, a13, a23, a14, a24, a34]
        self.trna = [u1, u2, u3]
        self.sa = sa

        # circle B
        self.rotb = [b12, b13, b23, b14, b24, b34]
        self.trnb = [v1, v2, v3]
        self.sb = sb

        return self


    # human readable string representation of object
    def __str__( self ):

        a12, a13, a23, a14, a24, a34 = self.rota
        b12, b13, b23, b14, b24, b34 = self.rotb
        u1, u2, u3 = self.trna
        v1, v2, v3 = self.trnb
        sa = self.sa
        sb = self.sb

        lsta = [( a12, a13, a23 ), ( a14, a24, a34 ), ( u1, u2, u3 ), sa ]
        lstb = [( b12, b13, b23 ), ( b14, b24, b34 ), ( v1, v2, v3 ), sb ]

        s = ''
        s += '\n--- SphereInput ---'
        s += '\n ' + str( lsta )
        s += '\n ' + str( lstb )
        s += '\n short_input = ' + str( [lsta] + [lstb] )
        s += '\n-------------------'

        return s
