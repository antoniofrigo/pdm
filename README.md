# pdm
Personal Phase Dispersion Minimization Calculator.

* Written in Fortran 2008.
* Uses [Stellingwerf 1978](https://ui.adsabs.harvard.edu/abs/1978ApJ...224..953S) phase dispersion algorithm.

### Basic Method

It's very unlikely that a variable star (such as a Cepheid) would be observed continuously, and instead repeated measurements are taken over the course of months or weeks at fairly large intervals. In order to determine the period, we have to fold the light curve into itself and somehow determine whether the resulting curve is correct or not. Stellingwerf's method takes the light curve which has been phased by a trial period t, bins the data by phase, and then compares the pooled variance of each of the bins to the overall variance of the data (F-test).

Stellingwerf also showed that the spectral line width in the frequency domain of the true period is 1/(2T) where T is the size of the interval in time. In order to accelerate our computations, we segment the data beforehand, Stellingwerf's method over these individual segments, and then use all of those variances to compute the pooled variance overall. Variances are calculated in one pass with Welford's algorithm.

### Compilation and Test

Simply run `python3 pdm_control.py`. This'll automatically run `make` and supply the relevant command-line arguments.

If one wishes to run it standalone, `pdm` takes four arguments:
```
./pdm data_file num_bins stdev_sep max_freq
```
- `data_file`: path to file containing two columns of data in ASCII format with first column of time and second magnitude
- `num_bins`: number of bins used
- `stdev_sep`: threshold for segmentation in standard deviations from the mean difference
- `max_freq`: the maximum frequency searched

For instance, to run the program for OGLE-BLG-CEP-027.dat with 10 bins, a 10 standard deviation threshold for segmentation, and a maximum frequency of 20/day, we'd run

```
make
./pdm OGLE-BLG-CEP-027.dat 10 10 20
```

Sample output for OGLE-BLG-CEP-027: 
```
Trial          14480 of           14480
Trial           2000 of            2000
  #     Period (d)          Theta
  1   0.2965302682   0.0453676414
 Time elapsed:    3.31755495   
STOP 0
```
This compares nicely to the published period of `P_1 = 0.2965282`.
