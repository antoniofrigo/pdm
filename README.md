# pdm
Personal Phase Dispersion Minimization Calculator.

* Written in Fortran 2008.
* Uses [Stellingwerf 1978](https://ui.adsabs.harvard.edu/abs/1978ApJ...224..953S) phase dispersion algorithm.

### Basic Method

It's very unlikely that a variable star (such as a Cepheid) would be observed continuously, and instead repeated measurements are taken over the course of months or weeks at fairly large intervals. In order to determine the period, we have to fold the light curve into itself and somehow determine whether the resulting curve is correct or not. Stellingwerf's method takes the light curve which has been phased by a trial period t, bins the data by phase, and then compares the pooled variance of each of the bins to the overall variance of the data (F-test).

Stellingwerf also showed that the spectral line width in the frequency domain of the true period is 1/(2T) where T is the size of the interval in time. In order to accelerate our computations, we segment the data beforehand, Stellingwerf's method over these individual segments, and then use all of those variances to compute the pooled variance overall. Variances are calculated in one pass with Welford's algorithm.

### Compilation and Test

Compile with `gfortran` using `make` and then run `./pdm`. Make sure to create a directory `MODULES/` beforehand in the same directory.

Sample output for OGLE-BLG-CEP-027: 
```
 Total cases:       29640
 
#       Theta  Period (d)
  1    0.159180    0.296519
  2    0.227178    0.296578
  3    0.350247    0.296459
 Total Time:   3.88123894    
STOP 0
```
This compares nicely to the published period of `P_1 = 0.2965282`.
