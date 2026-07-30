# Signal to Noise Estimators for the Squared Multiple Correlation Coefficient

Computes five estimators of the population signal to noise ratio
\\\phi^2 = \rho^2 / (1 - \rho^2)\\ associated with the squared multiple
correlation coefficient \\\rho^2\\. Two are functions of the unadjusted
and the Wherry-adjusted sample \\R^2\\; the other three are the Muirhead
(1985) unique minimum variance unbiased estimators that improve
substantially on the plug-in estimator at small *N* and modest numbers
of predictors.

## Usage

``` r
signal_to_noise_R2(R2, N, p)
```

## Arguments

- R2:

  The usual sample estimate of the squared multiple correlation
  coefficient (no degrees of freedom adjustment). Numeric scalar in
  \\(0, 1)\\.

- N:

  Sample size.

- p:

  Number of predictor variables.

## Value

A `data.frame` with columns `term` and `value` and one row per
estimator: `phi2_hat`, `phi2_adj_hat`, `phi2_umvue`, `phi2_umvue_l`, and
`phi2_umvue_nl`.

## Details

The signal to noise ratio \\\phi^2 = \rho^2 / (1 - \rho^2)\\ is a
natural reparameterization of \\\rho^2\\ that is bounded only below (at
zero) and so behaves more like a variance ratio than a proportion. It is
also the noncentrality parameter (up to a factor of *N*) for the omnibus
*F*-test of \\\rho^2 = 0\\ under fixed predictors; see
[`convert_R2_f`](https://yelleknek.github.io/DMAR/reference/convert_R2.md).

The five estimators returned, in increasing order of bias-correction
machinery, are:

- `phi2_hat`: the plug-in estimator \\\hat\phi^2 = R^2 / (1 - R^2)\\.
  Biased upward in small samples because the sample \\R^2\\ is itself
  biased upward.

- `phi2_adj_hat`: the plug-in estimator applied to the Wherry-adjusted
  \\R^2\\. Removes the leading-order bias in \\R^2\\ but is not itself
  unbiased for \\\phi^2\\.

- `phi2_umvue`: Muirhead's (1985) unique minimum variance unbiased
  estimator (his \\\theta_U\\); equivalent to Stuart, Ord, and
  Arnold's (1999) equation 28.97. Requires \\N \ge p + 6\\ (the gate on
  all three Muirhead estimators); for smaller \\N\\ the value is `NA`.

- `phi2_umvue_l`: Muirhead's (1985) linearly-improved unique minimum
  variance unbiased estimator (his \\\theta_L\\); equivalent to Stuart
  et al.\\ (1999) equation 28.98. Dominates `phi2_umvue` in mean squared
  error.

- `phi2_umvue_nl`: Muirhead's (1985) nonlinearly-improved estimator (his
  \\\theta\_{NL}\\). Dominates the linear improvement in MSE but
  requires \\p \ge 5\\; for smaller \\p\\ the value is `NA`.

Muirhead's nonlinear estimator is the recommended choice when \\p \ge
5\\; otherwise the linear estimator should be preferred over the plug-in
and adjusted-\\R^2\\ forms. All three Muirhead estimators are truncated
at zero, so the returned value is on the same scale as \\\phi^2\\.

As \\N\\ grows with \\p\\ fixed, the five estimators converge to a
common value (the population \\\phi^2\\); the difference between them is
the small-sample bias machinery in operation. The `@examples` block
illustrates that convergence.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Muirhead, R. J. (1985). Estimating a particular function of the multiple
correlation coefficient. *Journal of the American Statistical
Association, 80*, 923–925.

Stuart, A., Ord, J. K., & Arnold, S. (1999). *Kendall's advanced theory
of statistics, volume 2A: Classical inference and the linear model* (6th
ed.). Arnold.

## See also

[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`convert_R2_f`](https://yelleknek.github.io/DMAR/reference/convert_R2.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Fixed R^2 = 0.5 and p = 2, growing N: the five estimators agree
#    to within a small fraction once N is moderate.
signal_to_noise_R2(R2 = .5, N = 50,   p = 2)
#>  term          value
#>  phi2_hat      1    
#>  phi2_adj_hat  0.918
#>  phi2_umvue    0.878
#>  phi2_umvue_l  0.806
#>  phi2_umvue_nl <NA> 
signal_to_noise_R2(R2 = .5, N = 100,  p = 2)
#>  term          value
#>  phi2_hat      1    
#>  phi2_adj_hat  0.96 
#>  phi2_umvue    0.939
#>  phi2_umvue_l  0.901
#>  phi2_umvue_nl <NA> 
signal_to_noise_R2(R2 = .5, N = 500,  p = 2)
#>  term          value
#>  phi2_hat      1    
#>  phi2_adj_hat  0.992
#>  phi2_umvue    0.988
#>  phi2_umvue_l  0.98 
#>  phi2_umvue_nl <NA> 

# 2. With p = 5 the nonlinear estimator is available; it differs
#    most from the plug-in at small N.
signal_to_noise_R2(R2 = .5, N = 50,  p = 5)
#>  term          value
#>  phi2_hat      1    
#>  phi2_adj_hat  0.796
#>  phi2_umvue    0.755
#>  phi2_umvue_l  0.691
#>  phi2_umvue_nl 0.692
signal_to_noise_R2(R2 = .5, N = 500, p = 5)
#>  term          value
#>  phi2_hat      1    
#>  phi2_adj_hat  0.98 
#>  phi2_umvue    0.976
#>  phi2_umvue_l  0.968
#>  phi2_umvue_nl 0.968
```
