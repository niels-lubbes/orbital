'''
Created on Aug 8, 2016
@author: Niels Lubbes
'''

class PovInput:

    def __init__( self ):
        '''
        Constructor
        '''

        #
        # Implicit equation of a surface in QQ[x,y,z].
        #
        self.impl = None

        #
        # A dictionary of parametrizations.
        #
        # For the keys we use characters in [A-Z].
        #
        # The values are 2-tuples:
        #
        #    (<pmz>,<fam12>)
        #
        # Where
        #
        #   * <pmz> a list of 4 polynomials in
        #         QQ[c0,s0,c1,s1]/<c0^2+s0^2-1,c1^2+s1^2-1>
        #     which define a parametrization of surface
        #     in projective 3-space:
        #         S^1xS^1 ---> P^3.
        #
        #   * <fam12> is an integer either 0 or 1.
        #     If 0 then the 1st S^1 parametrizes a curve in
        #     a family. If 1 then the 2nd S^1.
        #
        self.pmz_dct = {}

        #
        # Path and filename
        #
        # Path should end with '/'-character
        # fname should be without extension
        #
        self.path = '/home/niels/Desktop/n/src/output/orb/povray/'
        self.fname = 'orb'

        #
        # Scaling
        #
        self.scale = 1

        #
        # Show axis
        #
        self.axes_dct = {}
        self.axes_dct['show'] = True
        self.axes_dct['len'] = 3

        #
        # Camera is first positioned at "self.cam_location"
        # while pointing at "self.cam_lookat".
        # Then the camera is rotated with angles
        # "self.cam_rotate".
        #
        self.cam_dct = {}
        self.cam_dct['location'] = ( 0, 0, -3 )
        self.cam_dct['lookat'] = ( 0, 0, 0 )
        self.cam_dct['rotate'] = ( 0, 0, 0 )

        #
        # Lights are positioned on a sphere of radius
        # " self.light_radius"
        #
        self.light_radius = 10

        #
        # The width, height and quality of image.
        #
        # For quality see:
        # http://www.povray.org/documentation/view/3.6.1/223/
        #
        # 0, 1      Just show quick colors.
        #           Use full ambient lighting only.
        #           Quick colors are used only at 5 or below.
        # 2, 3      Show specified diffuse and ambient light.
        # 4         Render shadows, but no extended lights.
        # 5         Render shadows, including extended lights.
        # 6, 7      Compute texture patterns, compute photons
        # 8         Compute reflected, refracted, and transmitted rays.
        # 9, 10, 11 Compute media and radiosity.
        #
        self.width = 100
        self.height = 75
        self.quality = 3  # quality is between 0-11

        #
        # Delay in ms between frames of animation
        #
        self.ani_delay = 10

        #
        # The OrbOutput object "self.o"
        # represents to families of curves A and B.
        #
        # Family A is represented as a list:
        #
        #   self.curve_lst_dct['A']=[ <curveA(0)>, <curveA(1)>, ... ]
        #
        # where
        #     * <curveA(#)> is a list of points on
        #       a curve in family A.
        #
        #     * The space between points in <curveA(#)>
        #       is determined by "step0".
        #
        #     * The space between <curve(n)> and <curve(n+1)>
        #       is determined by "step1".
        #
        #     * The precision of points in <curveA(#)> is
        #       determined by "prec".
        #
        # The curves in family A are sweeped out by spheres
        # with radius "width".
        #
        # See also "povray_aux.get_curve_lst()"
        # for how the following attributes are used.
        #
        self.curve_dct = {}
        self.curve_dct['A'] = {'step0':2 * 36, 'step1':36, 'prec':5, 'width':0.02}
        self.curve_dct['B'] = {'step0':2 * 36, 'step1':36, 'prec':5, 'width':0.02}
        self.curve_dct['C'] = {'step0':2 * 36, 'step1':36, 'prec':5, 'width':0.02}

        #
        # Place holders for the curve lists if computed
        # by method "povray_aux.get_curve_lst()".
        #
        self.curve_lst_dct = {}

        #
        # The textures where 'SURF' is the identifier
        # for surface and A-Z are id's for families of curves.
        #     [<show>, <rgbt>, <finish>]
        # where
        #    <show>   : A boolean. In .pov file: #if <show>...#end
        #    <rgbt>   : red, green, blue, transparency
        #    <finish> : see povray documentation for texture{ finish{...} }
        #
        self.text_dct = {}
        self.text_dct['SURF'] = [True, ( 0.2, 0.7, 0.3, 0.0 ), 'F_Glass10']
        self.text_dct['A'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
        self.text_dct['B'] = [True, ( 0.2, 0.3, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]
        self.text_dct['C'] = [True, ( 0.8, 0.6, 0.2, 0.0 ), 'phong 0.2 phong_size 5' ]



    # human readable string representation of object
    def __str__( self ):
        return 'PovInput<' + self.path + self.fname + '>'



