'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Aug 30, 2018
@author: Niels Lubbes
'''

from orbital.sage_interface import sage_matrix
from orbital.sage_interface import sage_identity_matrix
from orbital.sage_interface import sage__eval
from orbital.sage_interface import sage_var

from orbital.sphere.class_sphere_input import SphereInput

from orbital.sphere.sphere_experiment import get_pmz
from orbital.sphere.sphere_experiment import get_imp
from orbital.sphere.sphere_experiment import clifford

import os


class TestSphereExperiment( object ):

    def test__get_pmz( self ):

        A = sage__eval( '[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 119/169, -120/169, 0), (0, 0, 120/169, 119/169, 0), (0, 0, 0, 0, 1)]' )
        A = sage_matrix( A )
        B = sage_identity_matrix( 5 )

        baseA, baseB, pmzAB = get_pmz( A, B, 0 )

        print( baseA )
        print( baseB )
        print( pmzAB )

        assert len( baseA ) == 3
        assert len( baseB ) == 3
        assert len( pmzAB ) == 3


    def test__get_imp( self ):

        A = sage__eval( '[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 119/169, -120/169, 0), (0, 0, 120/169, 119/169, 0), (0, 0, 0, 0, 1)]' )
        A = sage_matrix( A )
        B = sage_identity_matrix( 5 )

        dct = get_imp( A, B, 0, False, False )

        key_lst = ['Agreat', 'Bgreat', 'eqn_x', 'eqn_str', 'eqn_xyz', 'sng_lst']
        Agreat, Bgreat, eqn_x, eqn_str, eqn_xyz, sng_lst = [ dct[key] for key in key_lst ]

        assert Agreat and Bgreat
        assert eqn_x.total_degree() == 4
        assert sng_lst == []


    def test__imp_pmz( self ):

        A = sage__eval( '[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 119/169, -120/169, 0), (0, 0, 120/169, 119/169, 0), (0, 0, 0, 0, 1)]' )
        A = sage_matrix( A )
        B = sage_identity_matrix( 5 )

        baseA, baseB, pmzAB = get_pmz( A, B, 0 )

        dct = get_imp( A, B, 0, False, False )
        key_lst = ['Agreat', 'Bgreat', 'eqn_x', 'eqn_str', 'eqn_xyz', 'sng_lst']
        Agreat, Bgreat, eqn_x, eqn_str, eqn_xyz, sng_lst = [ dct[key] for key in key_lst ]

        x, y, z = sage_var( 'x,y,z' )
        seqn = eqn_xyz.subs( {x:pmzAB[0], y:pmzAB[1], z:pmzAB[2]} ).simplify_trig()

        print( seqn )
        assert seqn == 0


    def test__clifford( self ):

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'

        lsta = [( 0, 0, 45 ), ( 0, 0, 0 ), ( 0, 0, 0 ), 1]
        lstb = [( 0, 0, 0 ), ( 0, 0, 0 ), ( 0, 0, 0 ), 1]

        sinp = SphereInput().set( [lsta] + [lstb] )
        sinp.imp = True
        sinp.sng = True
        sinp.snp = True
        sinp.srf = True
        sinp.bas = True
        sinp.fam = True
        sinp.famA = True
        sinp.famB = True

        plt, out = clifford( sinp )

        print( out )

        print( type( plt ) )
        print( type( out ) )

        assert str( type( plt ) ) == "<class 'sage.plot.plot3d.base.Graphics3dGroup'>"
        assert type( out ) == str


if __name__ == '__main__':

    # TestSphereExperiment().test__get_pmz()
    # TestSphereExperiment().test__get_imp()
    # TestSphereExperiment().test__imp_pmz()
    TestSphereExperiment().test__clifford()
