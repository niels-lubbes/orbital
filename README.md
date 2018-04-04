# Orbital 


## Introduction

Orbital is a Python library for constructing and rendering curves on surfaces.   

For the underlying theory of this package we refer to the article
"Surfaces that are covered by two families of circles" [pdf](https://arxiv.org/abs/1302.6710).

This library depends on [SageMath](https://SageMath.org) and [Povray](http://povray.org) libraries.

## Installation

* Install Sage from [SageMath](https://SageMath.org).
We assume that `sage` is accessible from your commandline interface.

* Install [Povray](http://povray.org).
We assume that `povray` is accessible from your commandline interface.

* Install the `orbital` package: 
```    
sage -pip install orbital
```    
If you do not have root access use the following command instead:
```    
sage -pip install --user orbital
```    

* We advice to upgrade the `orbital` package regularly:
```
sage -pip install --upgrade orbital
```
 
* To execute some [usecases](https://github.com/niels-lubbes/moebius_aut/blob/master/orbital/src/orbital/__main__.py) type:
```    
sage -python -m orbital
```

* For showing which files were installed 
or for uninstalling the `orbital` package, 
use one of the following commands:
```
sage -pip show --files orbital
sage -pip uninstall orbital
```


## Examples

See also [this file](https://github.com/niels-lubbes/orbital/blob/master/orbital/src/orbital/__main__.py) 
for example usecases. 
See the [source code](https://github.com/niels-lubbes/orbital/blob/master/orbital/src/orbital)
the io-specification of each function.
The [test functions](https://github.com/niels-lubbes/orbital/blob/master/orbital/src/tests)
might be informative for how to call each function.

For running the examples below, either copy paste the code into the Sage interface or run them as a Python module:

    sage -python -m my_module_name.py

### Example 1: Constructing celestial surfaces

A "celestial surface" is a surface that contains at least two circles through almost each point. 
Such a surface can be embedded in the projective n-sphere S^n for some n>=2. 
The n-sphere S^n is a hyperquadric of signature (n+1,1).
A surface in S^n of degree d that contains l circles through
almost each point has "celestial type": (l,d,n).
We have [shown](https://arxiv.org/abs/1302.6710) that n<=7 and l is either
infinity, or at most 6.

We denote the fiber product of the projective line with itself by P^1xP^1. 
Celestial surfaces that contain finitely many circles through each point
are the blowup of P^1xP^1 in either 0,2 or 4 complex conjugate points. 
We can parametize such blowups by constructing [linear series](https://github.com/niels-lubbes/linear_series) 
of forms of bidegree (2,2) that pass through the base points. 

We use  
[surface_in_quadric.get_surf()](https://github.com/niels-lubbes/orbital/blob/master/orbital/src/orbital/surface_in_quadric.py)
to compute from the parametrization the ideal of the surface.
We look in the ideal for a quadratic form of signature (n+1,1).  

```python
# We explicitly import the required modules
#
from linear_series.class_linear_series import LinearSeries
from linear_series.class_base_points import BasePointTree
from linear_series.class_poly_ring import PolyRing

from orbital.surface_in_quadric import get_surf

from orbital.class_orb_tools import OrbTools

# We enable verbose output (nothing is filtered out)
#
OrbTools.filter( None )

# We construct a parametrization of a sextic del Pezzo surface dP6
# in projective 6-space, that contains 3 conics through each point 
#
a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]
bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )
bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )
ls_dP6 = LinearSeries.get( [2, 2], bp_tree )
print( ls_dP6 )

# We show that dP6 can be projected into S^5 and S^4
#
dct61 = get_surf( ls_dP6, ( 5+1, 1 ) )
dct51 = get_surf( ls_dP6, ( 4+1, 1 ) )
```
Output:
    
    { 7, <<x^2*v^2 - y^2*w^2, x^2*v*w + y^2*v*w, x^2*w^2 + y^2*w^2, x*y*v^2 - y^2*v*w, x*y*v*w - y^2*w^2, y^2*v*w + x*y*w^2, y^2*v^2 + y^2*w^2>>, QQ( <a0|t^2 + 1> )[x, y, v, w] }
    
    get_surf(316): Computing random quadrics in ideal... 
    get_surf(356):           sig = (6, 1) , sig_set = set([(3, 4)]) 
    get_surf(356):           sig = (6, 1) , sig_set = set([(3, 4), (4, 3)]) 
    get_surf(356):           sig = (6, 1) , sig_set = set([(2, 5), (3, 4), (4, 3)]) 
    get_surf(356):           sig = (6, 1) , sig_set = set([(1, 6), (2, 5), (3, 4), (4, 3)]) 
    get_surf(356):           sig = (6, 1) , sig_set = set([(1, 6), (2, 5), (3, 4), (6, 1), (4, 3)]) 
    get_surf(358): Q        = [(1, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 1)] 
    get_surf(359): pmz_lst  = [-c0^2*c1^2 + s0^2*s1^2 - 2*s0^2*s1 - 2*s0*s1^2 + s0^2 + 4*s0*s1 + s1^2 - 2*s0 - 2*s1 + 1, -c0^2*c1*s1 - s0^2*c1*s1 + c0^2*c1 + s0^2*c1 + 2*s0*c1*s1 - 2*s0*c1 - c1*s1 + c1, c0^2*c1^2 + s0^2*c1^2 - 2*s0*c1^2 + c1^2, c0^2*c1*s1 - c0*s0*s1^2 - c0^2*c1 + 2*c0*s0*s1 + c0*s1^2 - c0*s0 - 2*c0*s1 + c0, -c0^2*c1^2 + c0*s0*c1*s1 - c0*s0*c1 - c0*c1*s1 + c0*c1, -c0*s0*c1^2 - c0^2*c1*s1 + c0^2*c1 + c0*c1^2, c0^2*c1^2 + c0^2*s1^2 - 2*c0^2*s1 + c0^2] 
    get_surf(360): imp_lst  = [x3*x5 + x5^2 - x2*x6 - x4*x6, x4^2 + x5^2 - x2*x6, x3*x4 + x4*x5 - x1*x6 + x5*x6, x2*x4 - x1*x5 + x2*x6, x1*x4 - x0*x5 + x1*x6 - x5*x6, x3^2 - x5^2 - x0*x6 + x2*x6 + 2*x4*x6, x2*x3 - x0*x5 + x1*x6 - x5*x6, x1*x3 - x0*x4 - x4*x6, x1^2 - x0*x2 - x2*x6, x0*x5^2 + x2*x5^2 - x2^2*x6 - 2*x1*x5*x6 + x5^2*x6 + x2*x6^2, x0*x4*x5 + x1*x5^2 - x1*x2*x6 - x0*x5*x6 + x4*x5*x6 + x1*x6^2 - x5*x6^2] 
    get_surf(361): c_lst    = [-7, -3, -2, 5, 5, -5, -2, 3, -8, -9, -4] 
    get_surf(362): M_pol    = -8*x1^2 + 8*x0*x2 + 3*x1*x3 - 2*x2*x3 - 5*x3^2 - 3*x0*x4 + 5*x1*x4 + 5*x2*x4 - 2*x3*x4 - 3*x4^2 - 3*x0*x5 - 5*x1*x5 - 7*x3*x5 - 2*x4*x5 - 5*x5^2 + 5*x0*x6 + 5*x1*x6 + 18*x2*x6 - 6*x4*x6 - 5*x5*x6 
    get_surf(375): M        = [(0, 0, 4, 0, -3/2, -3/2, 5/2), (0, -8, 0, 3/2, 5/2, -5/2, 5/2), (4, 0, 0, -1, 5/2, 0, 9), (0, 3/2, -1, -5, -1, -7/2, 0), (-3/2, 5/2, 5/2, -1, -3, -1, -3), (-3/2, -5/2, 0, -7/2, -1, -5, -5/2), (5/2, 5/2, 9, 0, -3, -5/2, 0)] 
    get_surf(376): U        = [(0.4222851601889938?, -1.642986325709679?, -1.741082149831869?, 0.566552639363350?, 1.613758669096323?, 0.622664009006899?, 1.936454943842840?), (0.4183110304040526?, 1.987338683834588?, -0.363622947184111?, 0.894221829882205?, 0.0615319212670577?, 2.095276628554762?, 0.2813654745496045?), (0.2607984896339884?, 1.200467513228684?, -1.194884742729225?, -1.963345463013132?, -0.03107664706610797?, -0.652905646874609?, 0.697593039586024?), (1.127122586930573?, -0.00467979593646299?, 0.01033916593830126?, -0.0703669096326930?, 0.591100268511707?, -0.1154273187487600?, -0.6753623868208163?), (0.1260788700673402?, -0.1311875787915946?, -0.03643156712249499?, -0.0834441446065686?, -0.2265073812454545?, 0.1354604011123243?, -0.001938362197660932?), (0.0487572393598491?, 0.0320855964509422?, -0.0779939768440050?, 0.1136406792879351?, -0.0816933670163349?, -0.1016817237269874?, 0.01399274430607368?), (1.30767090764666?, 0.3290941018473?, 2.14449422159360?, 0.0129850903388?, -0.12785892545682?, -0.4961779530816?, 2.18448489039277?)] 
    get_surf(377): J diag.  = [-1.000000000000000?, -1.000000000000000?, -1.000000000000000?, -1.000000000000000?, -1.000000000000000?, -1.00000000000000?, 1.000000000000000?] 

    get_surf(316): Computing random quadrics in ideal... 
    get_surf(356):           sig = (5, 1) , sig_set = set([(3, 3)]) 
    get_surf(356):           sig = (5, 1) , sig_set = set([(3, 3), (2, 2)]) 
    get_surf(356):           sig = (5, 1) , sig_set = set([(2, 4), (3, 3), (2, 2)]) 
    get_surf(356):           sig = (5, 1) , sig_set = set([(2, 4), (3, 3), (4, 2), (2, 2)]) 
    get_surf(356):           sig = (5, 1) , sig_set = set([(2, 4), (2, 3), (3, 3), (4, 2), (2, 2)]) 
    get_surf(356):           sig = (5, 1) , sig_set = set([(3, 2), (3, 3), (2, 3), (2, 2), (4, 2), (2, 4)]) 
    get_surf(356):           sig = (5, 1) , sig_set = set([(3, 2), (3, 3), (1, 5), (2, 3), (2, 2), (4, 2), (2, 4)]) 
    get_surf(356):           sig = (5, 1) , sig_set = set([(3, 2), (3, 3), (1, 5), (2, 3), (2, 2), (5, 1), (4, 2), (2, 4)]) 
    get_surf(358): Q        = [(0, 0, 1, 1, 1, 1, 1), (1, 1, 0, 1, 1, 0, 0), (0, 0, 1, 0, 1, 0, 0), (1, 0, 1, 0, 0, 0, 1), (0, 1, 0, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1, 1)] 
    get_surf(359): pmz_lst  = [c0^2*c1^2 - c0*s0*c1^2 + s0^2*c1^2 + c0*s0*c1*s1 + c0^2*s1^2 - c0*s0*s1^2 - c0*s0*c1 + c0*c1^2 - 2*s0*c1^2 - 2*c0^2*s1 + 2*c0*s0*s1 - c0*c1*s1 + c0*s1^2 + c0^2 - c0*s0 + c0*c1 + c1^2 - 2*c0*s1 + c0, -2*c0^2*c1^2 + c0*s0*c1*s1 - s0^2*c1*s1 - c0*s0*s1^2 + s0^2*s1^2 - c0*s0*c1 + s0^2*c1 + 2*c0*s0*s1 - 2*s0^2*s1 - c0*c1*s1 + 2*s0*c1*s1 + c0*s1^2 - 2*s0*s1^2 - c0*s0 + s0^2 + c0*c1 - 2*s0*c1 - 2*c0*s1 + 4*s0*s1 - c1*s1 + s1^2 + c0 - 2*s0 + c1 - 2*s1 + 1, s0^2*c1^2 + c0*s0*c1*s1 - c0*s0*c1 - 2*s0*c1^2 - c0*c1*s1 + c0*c1 + c1^2, c0^2*c1^2 + s0^2*c1^2 + c0^2*s1^2 + s0^2*s1^2 - 2*s0*c1^2 - 2*c0^2*s1 - 2*s0^2*s1 - 2*s0*s1^2 + c0^2 + s0^2 + c1^2 + 4*s0*s1 + s1^2 - 2*s0 - 2*s1 + 1, -c0*s0*c1^2 - c0^2*c1*s1 + c0*s0*c1*s1 - s0^2*c1*s1 + c0^2*s1^2 - c0*s0*s1^2 + c0^2*c1 - c0*s0*c1 + s0^2*c1 + c0*c1^2 - 2*c0^2*s1 + 2*c0*s0*s1 - c0*c1*s1 + 2*s0*c1*s1 + c0*s1^2 + c0^2 - c0*s0 + c0*c1 - 2*s0*c1 - 2*c0*s1 - c1*s1 + c0 + c1, -c0*s0*c1^2 + s0^2*c1^2 - c0^2*c1*s1 + c0*s0*c1*s1 - s0^2*c1*s1 + c0^2*s1^2 - c0*s0*s1^2 + s0^2*s1^2 + c0^2*c1 - c0*s0*c1 + s0^2*c1 + c0*c1^2 - 2*s0*c1^2 - 2*c0^2*s1 + 2*c0*s0*s1 - 2*s0^2*s1 - c0*c1*s1 + 2*s0*c1*s1 + c0*s1^2 - 2*s0*s1^2 + c0^2 - c0*s0 + s0^2 + c0*c1 - 2*s0*c1 + c1^2 - 2*c0*s1 + 4*s0*s1 - c1*s1 + s1^2 + c0 - 2*s0 + c1 - 2*s1 + 1] 
    get_surf(360): imp_lst  = [x1^2 - 2*x0*x2 + 2*x2^2 + 3*x1*x3 + 2*x3^2 - x0*x4 + 3*x1*x4 + 2*x2*x4 + 4*x3*x4 + 3*x4^2 - 5*x1*x5 - x2*x5 - 7*x3*x5 - 7*x4*x5 + 6*x5^2, x0*x1 + 2*x0*x2 - x1*x2 - 2*x2^2 + 2*x0*x3 - x1*x3 - 2*x2*x3 - 2*x3^2 + x0*x4 - 2*x1*x4 - 2*x2*x4 - 5*x3*x4 - 3*x4^2 - 3*x0*x5 + 2*x1*x5 + 4*x2*x5 + 7*x3*x5 + 9*x4*x5 - 6*x5^2, x0^2 - 2*x0*x2 + x2^2 - 2*x0*x3 + 2*x2*x3 + x3^2 - 2*x0*x4 + 2*x2*x4 + 3*x3*x4 + 2*x4^2 + 2*x0*x5 - 2*x2*x5 - 3*x3*x5 - 4*x4*x5 + 2*x5^2, x0*x2*x4 - x1*x2*x4 - x2^2*x4 - 2*x2*x3*x4 - x0*x4^2 - x1*x4^2 - x2*x4^2 - x3*x4^2 + x0*x2*x5 + x1*x2*x5 - x2^2*x5 + x0*x4*x5 + x1*x4*x5 + x2*x4*x5 + x3*x4*x5 + 2*x4^2*x5 - 2*x4*x5^2, 2*x1*x2*x3 + 2*x2*x3^2 + 3*x1*x2*x4 + 2*x0*x3*x4 + x1*x3*x4 + 3*x2*x3*x4 + 2*x0*x4^2 + x1*x4^2 - x3*x4^2 - 2*x4^3 - 3*x1*x2*x5 - x1*x3*x5 - 5*x2*x3*x5 - 2*x3^2*x5 - 4*x0*x4*x5 - 3*x1*x4*x5 - 2*x2*x4*x5 - 4*x3*x4*x5 + 2*x1*x5^2 + 4*x2*x5^2 + 7*x3*x5^2 + 8*x4*x5^2 - 6*x5^3, 2*x0*x2*x3 - 2*x2^2*x3 - 2*x2*x3^2 + 3*x1*x2*x4 - x1*x3*x4 + x2*x3*x4 - 2*x3^2*x4 + 2*x0*x4^2 + x1*x4^2 - x3*x4^2 - 2*x4^3 - 6*x0*x2*x5 - 3*x1*x2*x5 + 6*x2^2*x5 + x1*x3*x5 + 5*x2*x3*x5 + 2*x3^2*x5 - 2*x0*x4*x5 + x1*x4*x5 + 4*x2*x4*x5 + 8*x3*x4*x5 + 4*x4^2*x5 - 2*x1*x5^2 - 4*x2*x5^2 - 7*x3*x5^2 - 8*x4*x5^2 + 6*x5^3, 4*x0*x2^2 - 4*x2^3 - x1*x2*x4 - 2*x2^2*x4 + x1*x3*x4 - x2*x3*x4 + 2*x3^2*x4 + 2*x0*x4^2 + x1*x4^2 - 4*x2*x4^2 + 3*x3*x4^2 - 2*x0*x2*x5 + x1*x2*x5 + 4*x2^2*x5 - x1*x3*x5 + 3*x2*x3*x5 - 2*x3^2*x5 - 2*x0*x4*x5 - 3*x1*x4*x5 + 8*x2*x4*x5 - 10*x3*x4*x5 - 6*x4^2*x5 + 2*x1*x5^2 - 4*x2*x5^2 + 7*x3*x5^2 + 12*x4*x5^2 - 6*x5^3, 4*x1*x2^2*x4 + 8*x2^2*x3*x4 + 7*x1*x2*x4^2 + 6*x2^2*x4^2 + x1*x3*x4^2 + 11*x2*x3*x4^2 + 2*x3^2*x4^2 + 6*x0*x4^3 + 5*x1*x4^3 + 7*x3*x4^3 - 4*x1*x2^2*x5 - 18*x1*x2*x4*x5 - 8*x2^2*x4*x5 - 22*x2*x3*x4*x5 - 14*x0*x4^2*x5 - 16*x1*x4^2*x5 - 18*x2*x4^2*x5 - 21*x3*x4^2*x5 - 14*x4^3*x5 + 8*x0*x2*x5^2 + 11*x1*x2*x5^2 - 6*x2^2*x5^2 - x1*x3*x5^2 + 3*x2*x3*x5^2 - 2*x3^2*x5^2 + 8*x0*x4*x5^2 + 9*x1*x4*x5^2 + 22*x2*x4*x5^2 + 7*x3*x4*x5^2 + 34*x4^2*x5^2 + 2*x1*x5^3 - 4*x2*x5^3 + 7*x3*x5^3 - 14*x4*x5^3 - 6*x5^4, 4*x2^2*x3^2 + 8*x2^2*x3*x4 + x1*x3^2*x4 + 2*x2*x3^2*x4 + 2*x3^3*x4 + 3*x1*x2*x4^2 + 5*x2^2*x4^2 + x0*x3*x4^2 + 2*x1*x3*x4^2 + 5*x2*x3*x4^2 + 6*x3^2*x4^2 + 4*x0*x4^3 + 3*x1*x4^3 + 7*x3*x4^3 + x4^4 - 8*x2^2*x3*x5 - x1*x3^2*x5 - 2*x2*x3^2*x5 - 2*x3^3*x5 - 8*x1*x2*x4*x5 - 10*x2^2*x4*x5 - x0*x3*x4*x5 - 6*x1*x3*x4*x5 - 17*x2*x3*x4*x5 - 17*x3^2*x4*x5 - 10*x0*x4^2*x5 - 12*x1*x4^2*x5 - 10*x2*x4^2*x5 - 34*x3*x4^2*x5 - 16*x4^3*x5 + 4*x0*x2*x5^2 + 5*x1*x2*x5^2 + x2^2*x5^2 + 4*x1*x3*x5^2 + 8*x2*x3*x5^2 + 11*x3^2*x5^2 + 6*x0*x4*x5^2 + 13*x1*x4*x5^2 + 18*x2*x4*x5^2 + 47*x3*x4*x5^2 + 45*x4^2*x5^2 - 4*x1*x5^3 - 8*x2*x5^3 - 20*x3*x5^3 - 42*x4*x5^3 + 12*x5^4] 
    get_surf(361): c_lst    = [-10, 5, -5, -9, 1, -4, -4, 5, -2] 
    get_surf(362): M_pol    = -5*x0^2 + 5*x0*x1 - 10*x1^2 + 40*x0*x2 - 5*x1*x2 - 35*x2^2 + 20*x0*x3 - 35*x1*x3 - 20*x2*x3 - 35*x3^2 + 25*x0*x4 - 40*x1*x4 - 40*x2*x4 - 80*x3*x4 - 55*x4^2 - 25*x0*x5 + 60*x1*x5 + 40*x2*x5 + 120*x3*x5 + 135*x4*x5 - 100*x5^2 
    get_surf(375): M        = [(-5, 5/2, 20, 10, 25/2, -25/2), (5/2, -10, -5/2, -35/2, -20, 30), (20, -5/2, -35, -10, -20, 20), (10, -35/2, -10, -35, -40, 60), (25/2, -20, -20, -40, -55, 135/2), (-25/2, 30, 20, 60, 135/2, -100)] 
    get_surf(376): U        = [(1.6933941857178?, -2.89541050102?, -2.653367105320?, -5.86311521468?, -7.099211920161?, 9.871055997465?), (2.498117209683007?, 0.981499026628145?, -5.36367372770992?, 0.882766141511829?, -0.376882420455794?, -1.329144377493666?), (0.2322210338038690?, -0.06231997155305?, -0.397636002501175?, -0.621489333941490?, 2.11065326066461?, 0.983821111227560?), (0.188616770955087?, 0.74593096227414?, 0.172743575353930?, -0.552555715764062?, -0.063248859624552?, -0.140814616405223?), (0.181539062149018?, -0.31031003745072?, 0.045379811107120?, -0.27834291356717?, 0.019303561944726?, -0.261410666733980?), (2.0568502102608?, 0.0577183450160?, 0.9996834066278?, 0.9614768171519?, 0.00897099245090?, 0.5103327866067?)] 
    get_surf(377): J diag.  = [-1.000000000000000?, -1.000000000000000?, -1.000000000000000?, -1.000000000000000?, -1.00000000000000?, 1.000000000000000?] 


### Example 2: Constructing the Clifford translation of a circle along a circle

### Example 3: Computing and rendering a hexagonal web of conics on a surface. 

Creates povray image of a linear projection to 3-space 
of a smooth sextic del Pezzo surface in the projective 5-sphere. 
This surface contains 3 families of conics that form a hexagonal web. 
We render the families of conics by using Povray.

```python
# We explicitly import the required libraries.
#
from linear_series.class_linear_series import LinearSeries
from linear_series.class_base_points import BasePointTree
from linear_series.class_poly_ring import PolyRing

from orbital.surface_in_quadric import get_surf
from orbital.surface_in_quadric import approx_QQ
from orbital.surface_in_quadric import get_prj_mat
from orbital.surface_in_quadric import get_proj

from orbital.povray.class_pov_input import PovInput
from orbital.povray.povray import create_pov
from orbital.povray.povray_aux import get_time_str

from orbital.sage_interface import sage_QQ
from orbital.sage_interface import sage_var
from orbital.sage_interface import sage_vector
from orbital.sage_interface import sage_pi
from orbital.sage_interface import sage_set_verbose

from orbital.class_orb_tools import OrbTools

# Disable verbose output
#
sage_set_verbose( -1 )
OrbTools.filter( 'no output' )

# Compute linear series on P^1xP^1 of bi-degree (2,2) and (1,1)
# passing through two complex conjugate base points. 
# The linear series ls_AB of bi-degree (2,2) defines a map 
# that parametrizes a surface isomorphic to the blow up of P^1xP^1 
# in two points. This surface is a del Pezzo surface of degree 6. 
# Its ideal is generated by quadratic forms.
#
a0 = PolyRing( 'x,y,v,w', True ).ext_num_field( 't^2 + 1' ).root_gens()[0]
bp_tree = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
bp = bp_tree.add( 'xv', ( -a0, a0 ), 1 )
bp = bp_tree.add( 'xv', ( a0, -a0 ), 1 )
ls_AB = LinearSeries.get( [2, 2], bp_tree )
ls_CB = LinearSeries.get( [1, 1], bp_tree )

# We compute a surface that is contained in a hyperquadric of 
# signature (6,1) and that is the projection of the surface
# parametrized by ls_AB. The quadratic form of signature (6,1)
# has associated matrix M with orthogonal diagonalization
# M=U.T*J*U, where J is a diagonal matrix.
#
c_lst = [-1, -1, 0, 0, 0, -1, 1, 0, -1, -1, -1] # precomputed; set to None for computing new example
dct = get_surf( ls_AB, ( 6, 1 ), c_lst )
U, J = dct['UJ']
U.swap_rows( 0, 6 );J.swap_columns( 0, 6 );J.swap_rows( 0, 6 )
assert dct['M'] == approx_QQ( U.T * J * U )

# In order to visualize the surface we project it to 3-space
# with a linear map defined by the matrix P. The list 
# pmz_AB_lst defines a map that parametrizes the projected surface.
# If we fix the 1st and 2nd parameter we obtain a family of conics    
# called A and B respectively.
#
approxU = approx_QQ( U )
P = get_prj_mat( 4, 7, 0 )
P[0, 6] = -1;P[3, 3] = 0;P[3, 4] = 1
P = P * approxU
f_xyz, pmz_AB_lst = get_proj( dct['imp_lst'], dct['pmz_lst'], P )

# We compute a reparametrization pmz_CB_lst of the projected surface.
# If we fix the 1st and 2nd parameter we obtain a family of conics    
# called C and B respectively. 
# See also orbital.surface_in_quadric.get_S1xS1_pmz().
#
ring = PolyRing( 'x,y,v,w,c0,s0,c1,s1' )  # construct polynomial ring with new generators
x, y, v, w, c0, s0, c1, s1 = ring.gens()
X = 1 - s0; Y = c0; V = 1 - s1; W = c1;
CB_dct = { x:X, y:Y, v:X * W + Y * V, w: X * V - Y * W }
pmz_CB_lst = [ p.subs( CB_dct ) for p in ring.coerce( ls_AB.pol_lst )]
pmz_CB_lst = list( P * dct['Q'] * sage_vector( pmz_CB_lst ) )

# In order to render the projected surface we create a PovInput object.
#
pin = PovInput()
pin.path = '/home/niels/Desktop/' + get_time_str() + '_dp6_smooth/'
pin.fname = 'orb'
pin.scale = 1
pin.cam_dct['location'] = ( 0, 0, sage_QQ( -21 ) / 10 )
pin.cam_dct['lookat'] = ( 0, 0, 0 )
pin.cam_dct['rotate'] = ( 310, 0, 0 )
pin.shadow = True
pin.light_lst = [( 0, 0, -5 ), ( 0, -5, 0 ), ( -5, 0, 0 ), ( 0, 0, 5 ), ( 0, 5, 0 ), ( 5, 0, 0 ) ]
pin.axes_dct['show'] = False
pin.axes_dct['len'] = 1.2
pin.height = 400
pin.width = 800
pin.quality = 11
pin.ani_delay = 1
pin.impl = None
pin.pmz_dct['A'] = ( pmz_AB_lst, 0 )
pin.pmz_dct['B'] = ( pmz_AB_lst, 1 )
pin.pmz_dct['C'] = ( pmz_CB_lst, 0 )
v0_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 10 )]
v1_lst = [ ( sage_QQ( i ) / 180 ) * sage_pi for i in range( 0, 360, 15 )]
pin.curve_dct['A'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.02}
pin.curve_dct['B'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.02}
pin.curve_dct['C'] = {'step0':v0_lst, 'step1':v1_lst, 'prec':10, 'width':0.02}
col_A = ( 0.4, 0.0, 0.0, 0.0 )
col_B = ( 0.2, 0.3, 0.2, 0.0 )
col_C = ( 0.8, 0.6, 0.2, 0.0 )
pin.text_dct['A'] = [True, col_A, 'phong 0.2 phong_size 5' ]
pin.text_dct['B'] = [True, col_B, 'phong 0.2 phong_size 5' ]
pin.text_dct['C'] = [True, col_C, 'phong 0.2 phong_size 5' ]
print('pin.path =', pin.path)

# raytrace all families of conics on the projected surface using Povray.
# This takes a long time.
#
lst = create_pov( pin, ['A', 'B', 'C'] )

```
Output:

    pin.path = /home/niels/Desktop/2018-04-04__12-42-31_dp6_smooth/


![output image](https://raw.githubusercontent.com/niels-lubbes/orbital/master/orbital/img/deg6-dp6.png?token=AX_Io9IlHHRjrwsxr_I03zGKcN8FF31Eks5azdscwA%3D%3D "Hexagonal web on sextic del Pezzo surface")    
    

