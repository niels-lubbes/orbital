'''
Created on Aug 8, 2016
@author: Niels Lubbes
'''
from sage.all import *

import subprocess
from time import gmtime, strftime

from orb_verbose import *
from class_orb_ring import *
from class_pov_input import *


def pov_exp_lst( d, v, tbl = [] ):
    '''
    Called by "pov_coef_lst".
    
    INPUT:
        - "d"   -- A positive integer.
        - "v"   -- A positive integer.
        - "tbl" -- A list of lists of length "v", should always be []
                   unless called recursively on its output. 
    OUTPUT:
        - An ordered list of exponents in "v" variables.
          This is used for the povray poly object:
          
          http://www.povray.org/documentation/3.7.0/r3_4.html#r3_4_5_3_2
    '''

    if tbl == []:
        tbl = [ v * [0] ]

    if d == 0:
        return tbl

    ntbl = []
    for row in tbl:
        for idx in range( len( row ) ):
            nrow = copy( row )
            nrow[idx] += 1
            if nrow not in tbl + ntbl:
                ntbl += pov_exp_lst( d - 1, v, [nrow] )

    out_lst = []
    add_lst = []
    for val in ntbl:
        added = False
        for add in add_lst:
            added = added or ( add == val )
        if not added:
            out_lst += [val]
            add_lst += [val]

    return out_lst

def pov_coef_lst( poly ):
    '''
    INPUT:
        - "poly" -- A (String of a) polynomial in QQ[x,y,z].
    OUTPUT:
        - A 2-tuple (<degree poly>, <coefficient list of poly>)
          where the coefficient list is ordered according to 
          povray's poly object:
          
         http://www.povray.org/documentation/3.7.0/r3_4.html#r3_4_5_3_2
    '''

    R = PolynomialRing( QQ, var( 'x,y,z' ), order = 'degrevlex' )  # lower degree equations first
    x, y, z = R.gens()
    poly = sage_eval( str( poly ), R.gens_dict() )

    d = poly.total_degree()
    v = len( poly.variables() )

    exp_lst = pov_exp_lst( d, v + 1 )

    coef_lst = []
    for exp in exp_lst:
        coef_lst += [ poly.coefficient( {x:exp[0], y:exp[1], z:exp[2]} ) ]

    return d, coef_lst

def pov_nopow( poly ):
    '''
    INPUT:
        - "poly" -- A (String of a) polynomial in QQ[x,y,z].
    OUTPUT:
        - A string representation of the polynomial "poly" 
          without '^'-power symbols.
    EXAMPLE:
        - "x^3*y^2*z+x^2*y" --> "x*x*x*y*y*z+x*x*y"
    '''

    R = PolynomialRing( QQ, var( 'x,y,z' ), order = 'degrevlex' )  # lower degree equations first
    x, y, z = R.gens()
    poly = sage_eval( str( poly ), R.gens_dict() )

    dct = {}
    for i in range( 2, poly.total_degree() + 1 ):
        dct.update( {'x^' + str( i ): ( i * 'x*' )[:-1] } )
        dct.update( {'y^' + str( i ): ( i * 'y*' )[:-1] } )
        dct.update( {'z^' + str( i ): ( i * 'z*' )[:-1] } )

    return reduce( lambda x, y: x.replace( y, dct[y] ), dct, str( poly ) )


def get_pmz_value( pmz_lst, a, b, prec = 50 ):
    '''
    INPUT:
        - "pmz_lst" -- A list of 4 polynomials in 
                       QQ[c0,s0,c1,s1].
        - "a"       -- An integer in [0,360)
        - "b"       -- An integer in [0,360)
        - "prec"    -- Number of digits.
        - "scale"   -- A positive integer.
    OUTPUT:
        - Returns a point in R^3 
                F(a,b)=[ x, y, z  ] 
          with precision "prec".
          
          Here F(a,b) is represented by "pmz_lst" and
          has the following domain and range:  
          
              F: [0,360)X[0,360) ---> R^3.
          
          The parametrization is a map terms of cosine and sine.   
          
          If F is not defined at the point then we 
          return None.       
    '''

    ra = ( Rational( a ) / 180 ) * pi
    rb = ( Rational( b ) / 180 ) * pi

    c0, s0, c1, s1 = OrbRing.coerce( 'c0,s0,c1,s1' )

    W = pmz_lst[0].subs( {c0:cos( ra ), s0:sin( ra ), c1:cos( rb ), s1:sin( rb )} )
    X = pmz_lst[1].subs( {c0:cos( ra ), s0:sin( ra ), c1:cos( rb ), s1:sin( rb )} )
    Y = pmz_lst[2].subs( {c0:cos( ra ), s0:sin( ra ), c1:cos( rb ), s1:sin( rb )} )
    Z = pmz_lst[3].subs( {c0:cos( ra ), s0:sin( ra ), c1:cos( rb ), s1:sin( rb )} )

    if W == 0:
        return None

    XW = round( X / W, prec )
    YW = round( Y / W, prec )
    ZW = round( Z / W, prec )

    return [ XW, YW, ZW ]


def get_curve_lst( pin, fam ):
    '''
    INPUT:
        - "pin"  -- PovInput object where the following attributes 
                    are used:
                        * "pin.pmz_dct"
                        * "pin.scale"
                        * "pin.curve_dct"
                        * "pin.curve_lst_dct"
                                         
        - "fam"  -- A key string for a family id (eg. 'A')
        
    OUTPUT: 
    
        - Returns "pin.curve_lst_dct[<fam>]" 
          if 
              "pin.curve_lst_dct[<fam>]!=None", 
          else let us assume w.l.o.g. that "fam=='A'".
          Let 
        
            pin.curve_lst_dct['A'] = [ <curveA<0>>, <curveA<1>>, ... ]
            
          where 
          
              * <curveA<n>> is a list of points 
                    [x,y,z] 
                on a curve in family with id 'A', which are 
                scaled with
                    "pin.scale". 
          
              * The space between points in <curveA<n>> is 
                determined by "pin.curve_dct['A']['step0']".

              * The space between <curve<n>> and <curve<n+1>> is  
                determined by "pin.curve_dct['A']['step1']".

              * The precision of points in <curveA<n>> is 
                determined by "pin.curve_dct['A']['prec']".
          
          Returns pin.curve_lst_dct['A'].
    '''

    op( fam )

    if fam in pin.curve_lst_dct and pin.curve_lst_dct[fam] != None:
        op( 'Already computed ', fam )
        return pin.curve_lst_dct[fam]

    pin.curve_lst_dct[fam] = []
    for i1 in range( 0, 360, pin.curve_dct[fam]['step1'] ):
        curve = []
        for i0 in range( 0, 360, pin.curve_dct[fam]['step0'] ):
            pmz_lst, fam_id = pin.pmz_dct[fam]

            if fam_id == 0:
                point = get_pmz_value( pmz_lst, i0, i1, pin.curve_dct[fam]['prec'] )
            elif fam_id == 1:
                point = get_pmz_value( pmz_lst, i1, i0, pin.curve_dct[fam]['prec'] )
            else:
                raise ValueError( 'Expect pin.pmz_dct[fam][1] in [0,1]: ', fam_id )

            # add points to curve
            if point != None:  # map not defined
                point = [ coord * pin.scale for coord in point  ]
                curve += [point]

        # need at least 3 points for cubic interpolation
        if len( curve ) >= 3:

            # add curve to family
            pin.curve_lst_dct[fam] += [curve]

    return pin.curve_lst_dct[fam]


def get_time_str():
    '''
    OUTPUT:
        - A string of the current local time.
    EXAMPLE:
        - '2016-08-09__15-01-31'
    '''
    return strftime( "%Y-%m-%d__%H-%M-%S" )


def create_dir( file_name ):
    '''
    INPUT:
        - "file_name" -- An absolute path to a file.

    OUTPUT:
        - Creates the directory in which the file name 
          resides, in case this directory does not 
          exists.
              
        - Returns True if the directory was created and
          False otherwise.
    '''
    if not os.path.exists( os.path.dirname( file_name ) ):

        try:

            os.makedirs( os.path.dirname( file_name ) )

        except OSError as exc:
            # race condition ?
            if exc.errno != errno.EEXIST:
                raise

        return True
    else:
        return False




def convert_pngs_gif( path, fname, num_curves, ani_delay ):
    '''
    INPUT: 
        - "path"       -- Location (ending with '/')
                          of "ani/" directory containing
                          png-files with name 
                             "[fname]-#.png"
                          where # is a number in [0,num_curves]
        - "fname"      -- A string. 
        - "num_curves" -- A positive integer.
        - "ani_delay"  -- A positive integer.
    OUTPUT:
        An animated gif file with name "[path]+[fname].gif".
    '''
    # We now convert the ray traced images into an animated
    # gif using the following linux command:
    #
    #    convert -resize 768x576 -delay 20 -loop 0 `ls -v orb?*.png` orb-ani.gif
    #

    file_name_prefix = path + 'ani/' + fname
    op( file_name_prefix )
    cmd = ['convert']
    cmd += ['-delay', str( ani_delay )]
    cmd += [ '-loop', '0']
    for idx in range( 0, num_curves ):
        cmd += [file_name_prefix + '-' + str( idx ) + '.png']
    cmd += [path + fname + '.gif']

    p = subprocess.Popen( cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE )

    out, err = p.communicate()
    op( 'out =', out )
    op( 'err =', err )




