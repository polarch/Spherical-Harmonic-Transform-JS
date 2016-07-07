# Spherical-Harmonic-Transform-JS
Javascript routines for real spherical harmonic transform and interpolation of spherical data.

---
>
> Archontis Politis, 2016
>
> archontis.politis@aalto.fi
>
---

## Description

This small library contains the basic routines for computation of the spherical harmonic transform of spherical functions or data. Real spherical harmonics ```R_{nm}``` of order ```n``` and degree ```m``` are used, of the form

```
R_{nm}(\theta,\phi) = 
\sqrt{(2n+1)\frac{(n-|m|)!}{(n+|m|)!}} P_l^{|m|}(cos\theta) N_m(\phi)
```

where

```
N_m(\phi) = \sqrt{2} cos(m\phi},    m>0
N_m(\phi) = 1,                      m=0
N_m(\phi) = \sqrt{2} sin(|m|\phi},  m<0
```

and with ```P_l^{|m|}``` being the unnormalized associated Legendre functions. The functions in this convention are orthogonal, with an integrated product of ```4\pi```.

The library contains routines for the forward transform, based on either a weighted summation formula, suitable for certain sampling grids, or by a least-squares solution, suitable for more general and unstructured smapling points. The inverse transform can be used for reconstruction and spherical interpolation of the underlying distribution of the sampled data, at any direction.

Additional to the transform, rotation of spherical functions directly in the spherical harmonic domain is implemented, using the fast recursive algorithm of [Ivanic & Ruedenberg](http://pubs.acs.org/doi/pdf/10.1021/jp953350u). The library uses some matrix algebra routines from the [numeric.js](http://www.numericjs.com/) library.

Some uses of the library may be spectral analysis of spherical functions, analysis of directional data, interpolation of data on a sphere, and whatever involves spherical harmonics.

## Installation

If you use Node.js and NPM, you can install the library by

```
npm install spherical-harmonic-transform
```

in your ```node_modules``` directory, and then load it in your code by

```javascript
var sht = require('spherical-harmonic-transform');
```

If you intend to use it in the browser, you can directly include the generated bundle file inside the ```bundle``` directory in a ```<script/>``` tag

```html
<script src="spherical-harmonic-transform.js"></script>
```

The library is exposed globally as the variable ```sht```.

## Usage

The forward transform accepts data defined either as Cartesian vectors 
```javascript 
data = [ [x_1,y_1,y_1], ..., [x_K, y_K, z_K] ]
```
or in spherical coordinates as 
```javascript
data = [ [azi_1,elev_1,r_1], ..., [azi_K, elev_K, r_K] ]
```
with the syntax
```javascript
var sphSpectrum = sht.forwardSHT(maxOrder, data, CART_OR_SPH, DIRECT_OR_PINV)
```

The ```maxOrder``` should be a positive integer or zero, and defines the maximum order up to which all spectral coefficients are going to be computed. For a ```maxOrder``` of ```N```, ```(N+1)^2``` coeffs are returned. The choice of ```maxOrder``` should depend on the variability of the data and the number and arrangement of the sampling points. The flag ```CART_OR_SPH``` should be ```0``` for Cartesian vectors or ```1``` for spherical vectors. The flag ```DIRECT_OR_PINV``` should be ```0``` for an SHT by weighted summation, or ```1``` for a least-squares SHT. Spherical directions are azimuth and elevation, defined in rads.

The inverse transform is applied as
```javascript
var data = sht.inverseSHT(sphSpectrum, aziElev)
```
and reconstructs the function from its spherical spectrum, at directions defined in ```aziElev = [ [azi1,elev1], ..., [azi_M, elev_M] ]```. The ```data``` vector of length ```M``` returns the reconstructed or interpolated values.

Spherical harmonics can be computed for multiple directions ```aziElev = [ [azi1,elev1], ..., [azi_M, elev_M] ]``` by
```javascript
var sphHarm_mtx = sht.computeRealSH(maxOrder, aziElev)
```
with all harmonics up to order ```N``` returned as 
```javascript
sphHarm_mtx = [ [R_00(azi1,elev1), ..., R_00(azi_M,elev_M)], 
                ..., 
                [R_nm(azi1,elev1), ..., R_nm(azi_M,elev_M)], 
                ..., 
                [[R_NN(azi1,elev1), ..., R_NN(azi_M,elev_M)] ]
```

A spherical distribution can be rotated directly by applying a rotation matrix on its spherical spectrum. Such a matrix can be computed by
```javascript
var sphRot_mtx = sht.getSHrotMtx(Rxyz, maxOrder)
```
where ```maxOrder``` defines the maximum order of coefficients that are to be rotated, and ```Rxyz``` is a standard rotation matrix that describes the intended rotation of the coordinate system. If ```maxOrder = N```, the resulting matrix is of size ```(N+1)^2 x (N+1)^2```.

Associated legendre functions are computed internally by ```sht.computeRealSH```. If you need to compute them directly, a fast recursive function is included that computes all ```n+1``` polynomials for a certain degree ```n```, based on the polynomials of the two preceding degrees
```javascript
var P_n = sht.recurseLegendrePoly(n, x, Pn_minus1, Pn_minus2)
```
For ```n<4``` you don't need to implement the recursion, and the values are returned directly.

---

Apart from the above there are a few utility functions included.

The two routines
```javascript
var aziElevR = sht.convertCart2Sph(xyz, OMIT_MAG)
var xyz = sht.convertSph2Cart(aziElevR)
```
convert Cartesian vectors to spherical vectors, and inversely. If ```OMIT_MAG=1``` then only azimuth and elevation is returned (e.g. in case of unit vectors). Similarly if ```aziElevR``` does not include radius, unit vectors are returned.

A recursive factorial function is included as
```javascript
var fact = sht.factorial(N)
```

Two pseudo-inverse routines are included, borrowed from examples found in ```numeric.js```:
```javascript
var pinvA = sht.pinv_direct(A)
var pinvB = sht.pinv_svd(B)
```
where the first implements the left pseudo-inverse directly, while the second is based on the singular value decomposition.

A routine for generation of a standard rotation matrix is included, based on a Z-Y'-X'' convention
```javascript
var Rzyx = yawPitchRoll2Rzyx(angleZ, angleY, angleX)
```

## License

The library is released under the [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause).
