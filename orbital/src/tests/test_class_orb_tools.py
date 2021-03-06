'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 23, 2017
@author: Niels Lubbes
'''
from orbital.class_orb_tools import OrbTools

class TestClassOrbTools:


    def test__p( self ):

        OrbTools.filter( None )
        assert OrbTools.p( 'Hello world!' ) != None

        OrbTools.filter( ['another_class.py'] )
        assert OrbTools.p( 'No output since called from another class.' ) == None

        OrbTools.filter_unset()
        assert OrbTools.p( 'Filter is disabled so output this string.' ) != None

        OrbTools.filter_reset()
        assert OrbTools.p( 'Filter is enabled again so do not output.' ) == None

        OrbTools.filter( ['test_class_orb_tools.py'] )
        assert OrbTools.p( 'Only output if called from this class' ) != None


    def test__tool_dct( self ):

        orb = OrbTools()
        orb2 = OrbTools()

        # watch out to not use the default file name
        # otherwise it might take long to load the data
        test_fname = 'test_tools'
        key = 'test__tool_dct'

        dct = orb.get_tool_dct( fname = test_fname )
        dct[key] = True
        orb.save_tool_dct( fname = test_fname )

        assert key in orb.get_tool_dct( fname = test_fname )
        assert key in orb2.get_tool_dct( fname = test_fname )

        orb.set_enable_tool_dct( False )
        assert key not in orb.get_tool_dct( fname = test_fname )
        assert key not in orb2.get_tool_dct( fname = test_fname )

        orb.set_enable_tool_dct( True )
        assert key in orb.get_tool_dct( fname = test_fname )
        assert key in orb2.get_tool_dct( fname = test_fname )


