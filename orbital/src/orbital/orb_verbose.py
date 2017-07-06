'''
Created on Aug 4, 2016
@author: Niels Lubbes
'''
from sage.all import *

import inspect
import time
import sys


# global variable used by methods in this module.
orb_dct = None


def get_orb_dct( fname = "orb_dct" ):
    '''
    INPUT:
        - "fname" -- Name of file without extension.
    OUTPUT:
        - Loads global variable "orb_dct" 
          in memory from file "<local path>/<fname>.sobj"
          if called for the first time.
          
        - Returns "orb_dct".
    '''
    path = os.path.dirname( os.path.abspath( __file__ ) ) + '/'
    file_name = path + fname
    global orb_dct
    if orb_dct == None:
        try:
            # raise
            op( 'Loading from:', file_name )
            orb_dct = load( file_name )
        except Exception as e:
            op( 'Cannot load "orb_dct": ', e )
            orb_dct = {}

    return orb_dct


def save_orb_dct( fname = "orb_dct" ):
    '''
    INPUT:
        - "fname" -- Name of file without extension.        
    OUTPUT:
        - Saves global "orb_dct" to  "fname".
    '''
    path = os.path.dirname( os.path.abspath( __file__ ) ) + '/'
    file_name = path + fname
    global orb_dct

    op( 'Saving to:', file_name )
    save( orb_dct, file_name )



def op( *arg_lst ):
    '''
    INPUT:
        - "*arg_lst" -- List of arguments.
    OUTPUT:
        - * op(True): from now on all "op" calls are handled.
          * op(False, [file_name] ): only calls from [file_name] are handled.
                                     [file_name] should not be "" or "None".
          * Otherwise prints arguments to "sys.stdout" together with  
            reflection info from "inspect.stack()".
        - Returns output string or "None" if "arg_lst[0]" is "True" or "False".
    '''

    # check arguments
    if type( arg_lst[0] ) == type( True ):  # note that 0==False
        if arg_lst[0] == False:
            if len( arg_lst ) != 2:
                raise ValueError( 'If the 1st argument equals "False" then ' +
                                  'the 2nd argument is ' +
                                  'expected to be a file name. ', arg_lst )
            op.input_file_name = arg_lst[1]
            op( op.input_file_name )
            return None
        elif arg_lst[0] == True:
            op.input_file_name = None
            return None

    # collect relevant info from stack trace
    sk_lst_lst = inspect.stack()
    file_name = str( sk_lst_lst[1][1] )
    line = str( sk_lst_lst[1][2] )
    method_name = str( sk_lst_lst[1][3] )

    # only output when op is called from "op.input_file_name"
    if op.input_file_name != None:
        if not file_name.endswith( op.input_file_name ):
            return

    # construct output string
    s = method_name + '(' + line + ')' + ': '
    for arg in arg_lst:
        s += str( arg ) + ' '

    # print output
    print s
    sys.stdout.flush()

    return s

def ot( start = False ):
    '''
    INPUT:
        - "start" -- A boolean.
    OUTPUT:
        - If "start==True" then 0 is printed to "sys.stdout". 
          If "start==False" then outputs the following to "sys.stdout": 
          seconds passed since the last "ot(True)" call.       
          The outputs also contain reflection info from "inspect.stack()".
        - Returns output string.
    '''
    # get time
    if start:
        ot.t0 = time.clock()  # set static variable.
        dt = 0
    else:
        ct = time.clock()
        dt = ct - ot.t0
        ot.t0 = ct

    # collect relevant info from stack trace
    sk = inspect.stack()
    line = str( sk[1][2] )
    method_name = str( sk[1][3] )

    # construct output string
    s = method_name + '(' + line + ')[' + str( dt ) + ']'
    print s
    sys.stdout.flush()

    return s

