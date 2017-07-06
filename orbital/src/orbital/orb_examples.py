'''
Created on Aug 7, 2016
@author: Niels Lubbes
'''
from sage.all import *

from class_orb_input import *
from class_pov_input import *
from povray_aux import *


def get_exm_name_lst():
    return ['s3_perseus', 's3_ring', 's3_spindle', 's3_horn', 's3_sphere',
            's5d6', 's5d6_zoom', 's3d8']

def get_exm( name ):
    '''
    Returns output of method in the module by name.
    '''

    # The following inspection method does not work well and lists also
    # methods from imported modules.
    #    import inspect
    #    all_functions = inspect.getmembers(module, inspect.isfunction)

    if name == 's3_perseus':
        return exm_s3( ( 0, 0, 90 ), ( 1, 1, 0 ) )
    elif name == 's3_ring':
        return exm_s3( ( 0, 0, 90 ), ( 2, 0, 0 ) )
    elif name == 's3_spindle':
        return exm_s3( ( 0, 0, 90 ), ( 1, 0, 0 ) )
    elif name == 's3_horn':
        return exm_s3( ( 0, 0, 90 ), ( QQ( 1 ) / 2, 0, 0 ) )
    elif name == 's3_sphere':
        return exm_s3( ( 0, 0, 90 ), ( 0, 0, 0 ) )
    elif name == 's3d8':
        return exm_s3d8()
    elif name == 's5d6':
        return exm_s5d6( False )
    elif name == 's5d6_zoom':
        return exm_s5d6( True )
    else:
        raise ValueError( 'Unknown example name: ' + name )


def exm_s3( r_tup, t_tup ):
    '''
    Setup OrbInput and PovInput for 
    surfaces of revolution in S^3 of degree 4. 
    
    INPUT:
        - r_tup -- A tuple of 3 integers in [0,360]. 
        - t_tup -- A tuple of 3 integers.
    OUTPUT:
        - OrbInput and PovInput objects.
          
          * OrbInput specs:
            First rotate the unit circle in the xy-plane
                "r_tup[0]" degrees along z-axis, then
                "r_tup[1]" degrees along y-axis, and finally,
                "r_tup[2]" degrees along x-axis
            After that we translate along the vector
            represented by t_tup.
            The surface of revolution is obtained 
            by rotating the resulting circle along 
            the z-axis.
          
          * Perseus cyclide: [r_tup, t_tup] = [ (0,0,90), (1,1,0) ]
            It's stereographic projection has no isolated singularities, 
            since these are projected onto the Euclidean absolute conic 
            at infinity.
          * Ring cyclide:     [r_tup, t_tup] = [ (0,0,90), (2  ,0,0) ]
          * Spindle cyclide:  [r_tup, t_tup] = [ (0,0,90), (1  ,0,0) ]  
          * Horn cyclide:     [r_tup, t_tup] = [ (0,0,90), (1/2,0,0) ]
          * Sphere:           [r_tup, t_tup] = [ (0,0,90), (0  ,0,0) ]
    '''
    #
    # Set OrbInput object.
    #
    oin = OrbInput()
    p_tup = ( 'P0', 'I', 'I' )
    o_tup = ( 'I', 'Orppp', 'I' )  # rotations along z-axis
    v_tup = ( 'I', 'T' + str( list( t_tup ) + 4 * [0] ), 'X' + str( list( r_tup ) ) )
    oin.set( p_tup, o_tup, v_tup )
    for key in oin.do.keys(): oin.do[key] = True
    oin.do['bpt'] = False

    #
    # Init PovInput.
    #
    pin = PovInput()
    rt_str = str( [ r_tup, t_tup ] ).replace( ',', '_' ).replace( '(', '' ).replace( ')', '' ).replace( ' ', '' )[1:-1]
    pin.path = os.environ['OUTPUT_PATH'] + get_time_str() + '_exm_s3_' + rt_str + '/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, 0, -7 )
    pin.cam_dct['lookat'] = ( 0, -0.5, 0 )
    pin.cam_dct['rotate'] = ( 140, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10
    pin.text_dct['A'] = [True, ( 0.4, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.0, 0.3, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.curve_dct['A'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['C'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}

    return oin, pin


def exm_s5d6( zoom_sing ):
    '''
    Degree 6 in S^5. 
    
    INPUT:
        - zoom_sing -- A Boolean.
    OUTPUT:
        - OrbInput and PovInput objects for 
          a surface of degree 6 in S^5 with 
          an isolated singularity.
          If "zoom_sing" then the raytracing is 
          setup such that we zoom into the singularity.
    '''

    #
    # Set OrbInput object
    #
    oin = OrbInput()
    pmat = []
    pmat += [[1, 0, 0, 0] + [ 0, 0, 0, -1, -1]]
    pmat += [[0, 1, 0, 0] + [ 0, 0, 0, 0, 0]]
    pmat += [[0, 0, 1, 0] + [ 0, 0, 0, 0, 0]]
    pmat += [[0, 0, 0, 1] + [ -1, 0, 0, 0, 0]]
    p_tup = ( 'M' + str( pmat ), 'I', 'I' )
    o_tup = ( 'I', 'Oprpr', 'I' )
    v_tup = ( 'T[0, 1, 1, 0, 0, 0, 0]', 'Rrppr[37,0,0,90]', 'T[0, -1, -1, 0, 0, 0, 0]' )
    oin.set( p_tup, o_tup, v_tup )
    oin.do['bpt'] = False
    oin.do['sng'] = False

    # Init PovInput
    #
    pin = PovInput()
    pin.path = os.environ['OUTPUT_PATH'] + get_time_str() + '_exm_s5d6_a/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, 0, -3 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 250, 330, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10
    pin.text_dct['A'] = [True, ( 0.5, 0.0, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.text_dct['B'] = [True, ( 0.0, 0.3, 0.0, 0.0 ), 'phong 0.2 phong_size 5' ]
    pin.curve_dct['A'] = {'step0':5, 'step1':10, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':5, 'step1':10, 'prec':10, 'width':0.05}

    # set camera location such that
    # we get a close up of the singularity
    if zoom_sing:
        pin.cam_dct['location'] = ( 0, 0, -3 )
        pin.cam_dct['lookat'] = ( 0.5, 1, 0 )
        pin.cam_dct['rotate'] = ( 20, 60, 320 )
        pin.light_radius = 2

    return oin, pin


def exm_s3d8():
    '''
    Degree 8 in S^3. 
    '''

    #
    # Set OrbInput object
    #
    oin = OrbInput()
    p_tup = ( 'P0', 'I', 'I' )
    o_tup = ( 'I', 'Orrrp', 'I' )
    v_tup = ( 'T[1,0,1,0,1,0,0]', 'Rspmr[90,0,0,90]', 'T[-1,0,-1,0,-1,0,0]' )
    oin.set( p_tup, o_tup, v_tup )
    oin.do['bpt'] = False
    oin.do['sng'] = False

    #
    # Init PovInput
    #
    pin = PovInput()
    pin.path = os.environ['OUTPUT_PATH'] + get_time_str() + '_exm_s7d8/'
    pin.fname = 'orb'
    pin.scale = 1
    pin.cam_dct['location'] = ( 0, 0, QQ( -21 ) / 10 )
    pin.cam_dct['lookat'] = ( 0, 0, 0 )
    pin.cam_dct['rotate'] = ( 0, 0, 0 )
    pin.light_radius = 5
    pin.axes_dct['show'] = False
    pin.axes_dct['len'] = 1.2
    pin.width = 800
    pin.height = 400
    pin.quality = 11
    pin.ani_delay = 10
    pin.curve_dct['A'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['B'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}
    pin.curve_dct['C'] = {'step0':10, 'step1':15, 'prec':10, 'width':0.05}

    return oin, pin

