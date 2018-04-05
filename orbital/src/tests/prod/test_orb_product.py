'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Apr 5, 2018
@author: Niels Lubbes
'''
import os

from orbital.prod.orb_product import get_emb_dim
from orbital.prod.orb_product import get_deg_dim
from orbital.prod.orb_product import get_project
from orbital.prod.orb_product import get_factor_lst
from orbital.prod.orb_product import get_genus
from orbital.prod.orb_product import get_sing_lst
from orbital.prod.orb_product import get_pmz
from orbital.prod.orb_product import get_orb_bp_tree
from orbital.prod.orb_product import get_imp
from orbital.prod.orb_product import get_pmz_verify
from orbital.prod.orb_product import orb_product

from orbital.prod.class_orb_input import OrbInput
from orbital.prod.class_orb_output import OrbOutput

from orbital.class_orb_tools import OrbTools

from orbital.class_orb_ring import OrbRing

from orbital.prod.orb_matrices import get_pmat
from orbital.prod.orb_matrices import get_mat

from orbital.sage_interface import sage__eval
from orbital.sage_interface import sage_maple
from orbital.sage_interface import sage_magma
from orbital.sage_interface import sage_identity_matrix


class TestOrbInput( object ):


    def test__get_emb_dim( self ):
        imp_lst = OrbRing.coerce( '[-x0^2+x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2, x8, x7, x6, x5, x1*x2-x0*x4]' )
        emb = get_emb_dim( imp_lst )
        print( emb )
        assert emb == 3


    def test__get_deg_dim_1( self ):

        imp_lst = OrbRing.coerce( '[-x0^2+x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2, x8, x7]' )
        deg_dim = get_deg_dim( imp_lst )
        print( deg_dim )
        assert deg_dim == ( 2, 5 )


    def test__get_deg_dim_2( self ):

        imp_lst = OrbRing.coerce( '[-x0^2+x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2, x8, x7, x6, x5, x4, x3]' )
        deg_dim = get_deg_dim( imp_lst )
        print( deg_dim )
        assert deg_dim == ( 2, 1 )


    def test__get_deg_dim_3( self ):

        imp_lst = OrbRing.coerce( '[-x0^2+x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2, x8, x7, x6, x5, x1*x2-x0*x4]' )
        deg_dim = get_deg_dim( imp_lst )
        print( deg_dim )
        assert deg_dim == ( 4, 2 )


    def test__get_project( self ):

        pol_lst = OrbRing.coerce( '[-x0^2+x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2, x8, x7, x6, x5, x1*x2-x0*x4]' )
        pmat = get_pmat( False )
        print( pmat )

        out = get_project( pol_lst, pmat )
        print( out )
        assert str( out ) == '(x0^4 - x0^2*x1^2 - x0^2*x2^2 - x1^2*x2^2 - x0^2*x3^2, -x^2*y^2 - x^2 - y^2 - z^2 + 1)'


    def test__get_factor_lst( self ):

        pol = OrbRing.coerce( '(x1-x0)^3*(x1-3*x0)^2*(x2^2+x3^2)' )
        print( pol )

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'
        try:
            sage_maple.eval( 'with(algcurves);' )
            print( 'Maple is installed.' )
        except:
            print( 'Maple is NOT installed.' )
            fct_lst = get_factor_lst( pol )
            print( fct_lst )
            assert fct_lst == []
            return

        fct_lst = get_factor_lst( pol )
        print( fct_lst )
        assert str( fct_lst ) == "[('x2-RootOf(_Z^2+1)*x3', '1'), ('x2+RootOf(_Z^2+1)*x3', '1'), ('x0-1/3*x1', '2'), ('x0-x1', '3')]"


    def test__get_genus( self ):

        pol = OrbRing.coerce( '(x1^2+x2^2+x3^2+3*x0^2)^2-16*(x1^2+x2^2)' )
        print( pol )

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'
        try:
            sage_maple.eval( '1 + 1' )
            print( 'Maple is installed.' )
        except:
            print( 'Maple is NOT installed.' )
            gen = get_genus( pol )
            print( gen )
            assert gen == -3
            return

        gen = get_genus( pol )
        print( gen )
        assert gen == 1


    def test__get_sing_lst( self ):
        pol = OrbRing.coerce( '(x1^2+x2^2+x3^2+3*x0^2)^2-16*x0^2*(x1^2+x2^2)' )
        print( pol )

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'
        try:
            sage_magma.eval( 1 + 1 )
            print( 'Magma is installed.' )
        except:
            print( 'Magma is NOT installed.' )
            sng_lst = get_sing_lst( pol )
            print( sng_lst )
            assert sng_lst == []
            return

        sng_lst = get_sing_lst( pol )
        print( sng_lst )
        assert str( sng_lst ) == "[('[x0^2 + 1/3*x3^2,x1,x2]', 2), ('[x0,x1^2 + x2^2 + x3^2]', 2*t + 1)]"


    def test__get_pmz( self ):
        pmat = get_pmat( False )
        omat = sage_identity_matrix( 9 )
        vmat = sage_identity_matrix( 9 )
        pmz_lst, prj_pmz_lst = get_pmz( pmat, omat, vmat )
        print( pmz_lst )
        print( prj_pmz_lst )
        assert str( pmz_lst ) == '[1, c1, s1, 0, 0, 0, 0, 0, 0]'
        assert str( prj_pmz_lst ) == '[1, c1, s1, 0]'


    def test__get_orb_bp_tree( self ):
        pmz_lst = OrbRing.coerce( '[2,c0,s0,c1,s1,c0,s0,c1,s1]' )
        bp_tree = get_orb_bp_tree( pmz_lst )
        print( bp_tree )

        out = str( bp_tree )

        assert 'chart=xv, depth=0, mult=1, sol=(a0, (a0))' in out
        assert 'chart=xv, depth=0, mult=1, sol=(a0, (-a0))' in out
        assert 'chart=xv, depth=0, mult=1, sol=(-a0, (-a0))' in out
        assert 'chart=xv, depth=0, mult=1, sol=(-a0, (a0))' in out


    def test__get_imp( self ):

        omat = get_mat( 'T[1, 0, 0, 0, 0, 0, 0]', 'Orppp', 'T[-1, 0, 0, 0, 0, 0, 0]' )
        vmat = get_mat( 'T[0, 1, 1, 0, 0, 0, 0]', 'Rrrrs[37, 0, 0, 0]', 'T[0, -1, -1, 0, 0, 0, 0]' )

        imp_lst = get_imp( omat, vmat )  # ideal of Perseus cyclide
        print( imp_lst )
        assert str( imp_lst ) == '[x7, x6, x5, x4, 100*x1^2 - 4*x0*x3 - 40*x1*x3 + 9*x3^2 - 200*x1*x8 + 44*x3*x8 + 100*x8^2, 100*x0^2 - 100*x2^2 - 4*x0*x3 - 40*x1*x3 - 91*x3^2 - 200*x1*x8 + 44*x3*x8]'


    def test__get_pmz_verify__perseus( self ):

        o = OrbOutput( OrbInput() )

        o.input.do['pmz'] = True
        o.input.do['imp'] = True

        o.prj_pmz_lst = '[-2*s1 + 3, 4/5*c0*c1 - 3/5*s0*c1 + 7/5*c0*s1 + 6/5*s0*s1 - 12/5*c0 - 11/5*s0 - 2*s1 + 3, 3/5*c0*c1 + 4/5*s0*c1 - 6/5*c0*s1 + 7/5*s0*s1 + 11/5*c0 - 12/5*s0, -2*s1 + 2]'
        o.prj_pmz_lst = OrbRing.coerce( o.prj_pmz_lst )

        o.prj_pol = '25*x0^4 + 100*x0^3*x1 + 50*x0^2*x1^2 - 100*x0*x1^3 + 25*x1^4 - 50*x0^2*x2^2 - 100*x0*x1*x2^2 + 50*x1^2*x2^2 + 25*x2^4 - 24*x0^3*x3 - 40*x0^2*x1*x3 + 20*x0*x1^2*x3 + 20*x0*x2^2*x3 - 41*x0^2*x3^2 - 100*x0*x1*x3^2 + 50*x1^2*x3^2 + 50*x2^2*x3^2 + 20*x0*x3^3 + 25*x3^4'
        o.prj_pol = OrbRing.coerce( o.prj_pol )

        o.gen = 1

        print( o )
        tst = get_pmz_verify( o )
        print( tst )
        assert tst == True


    def test__orb_product__65_smooth( self ):

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'
        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'

        s = "['@(6,5)=(deg,emb)', {'pmat': ('P1', 'I', 'I'), 'omat': ('T[1, -1, 1, 1, -1, 0, 0]', 'Opsms', 'I'), 'vmat': ('T[1, 0, -1, 0, 0, 0, 0]', 'Rramr[340, 225, 264, 320]', 'T[-1, 0, 1, 0, 0, 0, 0]')}]"

        input = OrbInput().set_short_str( s )
        input.do['pmz'] = True
        input.do['bpt'] = False
        input.do['imp'] = True
        input.do['dde'] = True
        input.do['prj'] = True
        input.do['fct'] = False
        input.do['gen'] = False
        input.do['sng'] = False
        input.do['tst'] = False

        o = orb_product( input )
        print( o )

        assert ( o.deg, o.emb, o.dim ) == ( 6, 5, 2 )


    def test__orb_product__65_sing( self ):

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'
        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'

        s = "['@(6,5)=(deg,emb)', {'pmat': ('P1', 'I', 'I'), 'omat': ('T[0, 0, 0, 1, 1, -1, -1]', 'Oprps', 'I'), 'vmat': ('T[0, -1, 0, -1, 0, 0, 0]', 'Rspps[148, 344, 284, 304]', 'T[0, 1, 0, 1, 0, 0, 0]')}]"

        input = OrbInput().set_short_str( s )
        input.do['pmz'] = True
        input.do['bpt'] = True
        input.do['imp'] = True
        input.do['dde'] = True
        input.do['prj'] = True
        input.do['fct'] = False
        input.do['gen'] = False
        input.do['sng'] = False
        input.do['tst'] = False

        print( input )

        o = orb_product( input )
        print( o )

        assert ( o.deg, o.emb, o.dim ) == ( 6, 5, 2 )

        assert 'chart=xv, depth=0, mult=1, sol=(a1, 0)' in str( o.bp_tree )
        assert 'chart=xv, depth=0, mult=1, sol=(-a1, 0)' in str( o.bp_tree )


    def test__orb_product__43_perseus( self ):

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'
        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/magma/link'

        # perseus cyclide
        s = "['@(4,3)=(deg,emb)', {'pmat': ('P0', 'I', 'I'), 'omat': ('T[1, 0, 0, 0, 0, 0, 0]', 'Orppp', 'T[-1, 0, 0, 0, 0, 0, 0]'), 'vmat': ('T[0, 1, 1, 0, 0, 0, 0]', 'Rrrrs[37, 0, 0, 0]', 'T[0, -1, -1, 0, 0, 0, 0]')}]"

        input = OrbInput().set_short_str( s )
        input.do['pmz'] = True
        input.do['bpt'] = False
        input.do['imp'] = True
        input.do['dde'] = True
        input.do['prj'] = True
        input.do['fct'] = True
        input.do['gen'] = True
        input.do['sng'] = True
        input.do['tst'] = True

        print( input )

        o = orb_product( input )
        print( o )

        assert ( o.deg, o.emb, o.dim ) == ( 4, 3, 2 )

        # ------------------------------
        # ...............
        # pmat          = P0 ~~~ I ~~~ I
        # omat          = T[1, 0, 0, 0, 0, 0, 0] ~~~ Orppp ~~~ T[-1, 0, 0, 0, 0, 0, 0]
        # vmat          = T[0, 1, 1, 0, 0, 0, 0] ~~~ Rrrrs[37, 0, 0, 0] ~~~ T[0, -1, -1, 0, 0, 0, 0]
        # do            = {'sng': True, 'fct': True, 'pmz': True, 'prj': True, 'imp': True, 'bpt': False, 'tst': True, 'gen': True, 'dde': True, 'deg': True}
        # ...............
        # pmz_lst       = [4/5*c0*c1 - 3/5*s0*c1 + 7/5*c0*s1 + 6/5*s0*s1 - 12/5*c0 - 11/5*s0 - 1/5*c1 - 18/5*s1 + 28/5, 4/5*c0*c1 - 3/5*s0*c1 + 7/5*c0*s1 + 6/5*s0*s1 - 12/5*c0 - 11/5*s0 - 2*s1 + 3, 3/5*c0*c1 + 4/5*s0*c1 - 6/5*c0*s1 + 7/5*s0*s1 + 11/5*c0 - 12/5*s0, -2*s1 + 2, 0, 0, 0, 0, 4/5*c0*c1 - 3/5*s0*c1 + 7/5*c0*s1 + 6/5*s0*s1 - 12/5*c0 - 11/5*s0 - 1/5*c1 - 8/5*s1 + 13/5]
        # prj_pmz_lst   = [-2*s1 + 3, 4/5*c0*c1 - 3/5*s0*c1 + 7/5*c0*s1 + 6/5*s0*s1 - 12/5*c0 - 11/5*s0 - 2*s1 + 3, 3/5*c0*c1 + 4/5*s0*c1 - 6/5*c0*s1 + 7/5*s0*s1 + 11/5*c0 - 12/5*s0, -2*s1 + 2]
        # imp_lst       = [x7, x6, x5, x4, 100*x1^2 - 4*x0*x3 - 40*x1*x3 + 9*x3^2 - 200*x1*x8 + 44*x3*x8 + 100*x8^2, 100*x0^2 - 100*x2^2 - 4*x0*x3 - 40*x1*x3 - 91*x3^2 - 200*x1*x8 + 44*x3*x8]
        # emb           = 3
        # dim           = 2
        # deg           = 4
        # gen           = 1
        # prj_pol       = 25*x0^4 + 100*x0^3*x1 + 50*x0^2*x1^2 - 100*x0*x1^3 + 25*x1^4 - 50*x0^2*x2^2 - 100*x0*x1*x2^2 + 50*x1^2*x2^2 + 25*x2^4 - 24*x0^3*x3 - 40*x0^2*x1*x3 + 20*x0*x1^2*x3 + 20*x0*x2^2*x3 - 41*x0^2*x3^2 - 100*x0*x1*x3^2 + 50*x1^2*x3^2 + 50*x2^2*x3^2 + 20*x0*x3^3 + 25*x3^4
        # prj_pol{x0:0} = (25) * (x1^2 + x2^2 + x3^2)^2
        # xyz_pol       = 25*x^4 + 50*x^2*y^2 + 25*y^4 + 50*x^2*z^2 + 50*y^2*z^2 + 25*z^4 - 100*x^3 - 100*x*y^2 + 20*x^2*z + 20*y^2*z - 100*x*z^2 + 20*z^3 + 50*x^2 - 50*y^2 - 40*x*z - 41*z^2 + 100*x - 24*z + 25
        # pmz_test      = True
        # fct_lst       = 1 factors
        # sng_lst       = 1 components
        #               ~ ('[x0,x1^2 + x2^2 + x3^2]', 2*t + 1)
        # short_str     =
        #     s="['@(4,3)=(deg,emb)', {'pmat': ('P0', 'I', 'I'), 'omat': ('T[1, 0, 0, 0, 0, 0, 0]', 'Orppp', 'T[-1, 0, 0, 0, 0, 0, 0]'), 'vmat': ('T[0, 1, 1, 0, 0, 0, 0]', 'Rrrrs[37, 0, 0, 0]', 'T[0, -1, -1, 0, 0, 0, 0]')}]"
        # ------------------------------




if __name__ == '__main__':

    OrbTools.filter( None )

#     TestOrbInput().test__get_emb_dim()
#     TestOrbInput().test__get_deg_dim_1()
#     TestOrbInput().test__get_deg_dim_2()
#    TestOrbInput().test__get_deg_dim_3()
#    TestOrbInput().test__get_project()
#     TestOrbInput().test__get_factor_lst()
#     TestOrbInput().test__get_genus()
#     TestOrbInput().test__get_sing_lst()
#     TestOrbInput().test__get_pmz()
#     TestOrbInput().test__get_orb_bp_tree()
#     TestOrbInput().test__get_imp()
#     TestOrbInput().test__get_pmz_verify__perseus()
    TestOrbInput().test__orb_product__65_smooth()
    TestOrbInput().test__orb_product__65_sing()
    TestOrbInput().test__orb_product__43_perseus()

    pass
