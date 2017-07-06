'''
Created on Aug 10, 2016

@author: Niels Lubbes
'''
from sage.all import *
import os
from orb_verbose import *

pt_dct = None  # global variable used by get_pt_lst


def get_cs( angle ):
    '''
    INPUT:
        - "angle" -- An integer in [0,360).
    OUTPUT:
        - A pair of rational numbers (a,b)
          such that a^2+b^2=1
          and ( a, b ) is a close approximation 
          of ( cos(angle/180*pi), sin(angle/180*pi) ).
    '''
    angle = angle % 360

    radi = ( Rational( angle ) / 180 ) * pi

    if sin( radi ) in QQ and cos( radi ) in QQ:
        return ( QQ( cos( radi ) ), QQ( sin( radi ) ) )

    p0, p1, p2 = get_pt_dct()[angle % 90]

    c = QQ( p0 ) / p2
    s = QQ( p1 ) / p2

    #  90 = pi/2
    #  sin( a + 90 ) =  cos(a)
    #  cos( a + 90 ) = -sin(a)
    dct = {0:( c, s ), 1:( -s, c ), 2:( -c, -s ), 3:( c, s )}

    rc, rs = dct[int( angle ) / int( 90 )]

    return rc, rs




def get_pt_dct( fname = 'cos_sin' ):
    '''
    INPUT:
        - "fname" -- Name of file without extention
                     The file should contain 3 integers on each line 
                     separated by spaces. We expect them to be 
                     Pythagorian triples.
    OUTPUT:
        - Sets "global pt_lst".
        - Returns a dictionary 
            { <angle>:[<a>,<b>,<c>],... }
          where <a>^2+<b>^2=<c>^2
          and <angle> corresponds to round(arctan( <b>/<a> )*180/pi)
          and <angle> runs from 1 to 89 degrees.
           
          The list of triples was obtained from:
          http://www.tsm-resources.com/alists/PythagTriples.txt
    '''

    path = os.path.dirname( os.path.abspath( __file__ ) ) + '/'
    file_name = path + fname

    global pt_dct
    if pt_dct != None:
        return pt_dct

    try:
        pt_dct = load( file_name )

    except:

        op( 'Calculating Pythagorian triples and angles from:', file_name )

        angle_lst = []
        pt_dct = {}
        with open( file_name + '.txt', 'r' ) as f:
            for line in f:
                ps_lst = line.replace( '\r\n', '' ).split()
                pt0 = [ QQ( ps ) for ps in ps_lst ]  # need to be divisable!

                # Try all combinations
                # while a triple is still small.
                # We assume that the pythagorian triples are
                # ordered on coeffitient size in the input file.
                #
                pt1a = [ pt0[0], pt0[1], pt0[2]]
                pt1b = [ -pt0[0], pt0[1], pt0[2]]
                pt1c = [ pt0[0], -pt0[1], pt0[2]]
                pt1d = [ -pt0[0], -pt0[1], pt0[2]]

                pt2a = [ pt0[1], pt0[0], pt0[2]]
                pt2b = [ -pt0[1], pt0[0], pt0[2]]
                pt2c = [ pt0[1], -pt0[0], pt0[2]]
                pt2d = [ -pt0[1], -pt0[0], pt0[2]]

                for pt in [pt1a, pt1b, pt1c, pt1d, pt2a, pt2b, pt2c, pt2d]:

                    if pt[0] ** 2 + pt[1] ** 2 != pt[2] ** 2:
                        raise ValueError( 'Expect a file containing Pythagorian triples:', pt )

                    # cos = pt[0]/pt[2], sin = pt[1]/pt[2], tan=sin/cos
                    angle = round( arctan( pt[1] / pt[0] ) * 180 / pi )

                    if angle not in angle_lst and angle > 0:
                        angle_lst += [angle]
                        pt_dct.update( {angle:pt} )

        op( len( pt_dct.keys() ) )
        save( pt_dct, file_name )

    return pt_dct
