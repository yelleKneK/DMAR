# A Priori Monte Carlo Simulation for Sample Size Planning for RMSEA in SEM

Conduct a priori Monte Carlo simulation to empirically study the effects
of (mis)specifications of input information on the calculated sample
size. The sample size is planned so that the expected width of a
confidence interval for the population RMSEA is no larger than desired.
Random data are generated from the true covariance matrix but fit to the
proposed model, whereas the sample size is calculated based on the input
covariance matrix and proposed model.

## Usage

``` r
ss_aipe_rmsea_sensitivity(
  width,
  model,
  Sigma,
  N = NULL,
  conf_level = 0.95,
  G = 200,
  filename = NULL,
  ...
)
```

## Arguments

- width:

  desired confidence interval width for the population RMSEA.

- model:

  the model the researcher proposes, which may or may not be the true
  model, written in lavaan model syntax (see
  [`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)).
  The observed variable names in the model must match the row and column
  names of `Sigma`.

- Sigma:

  the true population covariance matrix, which is used to generate
  random data for the simulation study. The row and column names of
  `Sigma` must match the observed variables in `model`.

- N:

  if `N` is specified, random samples of the specified size are
  generated. Otherwise the sample size is calculated with the sample
  size planning method so that the expected width of a confidence
  interval for the population RMSEA is no larger than `width`.

- conf_level:

  confidence level (i.e., 1 - the Type I error rate).

- G:

  number of replications in the Monte Carlo simulation.

- filename:

  an optional path for a comma separated file recording every converged
  replication (its index, the RMSEA estimate, the two confidence limits,
  and the interval width): nothing is written when `filename` is `NULL`
  (the default), a new file with a header row is created otherwise, an
  existing file at that path is appended to, and a throwaway run should
  point it at `tempfile(fileext = ".csv")`.

- ...:

  additional arguments passed to
  [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) when fitting the
  model (for example `estimator` or `missing`).

## Value

A `data.frame` with columns `term` and `value` summarizing the a priori
Monte Carlo study. The `term` entries are: `"mean_rmsea"`,
`"median_rmsea"`, `"sd_rmsea"` (summaries of the realized RMSEA
estimates across the converged replications); `"mean_ci_width"`,
`"median_ci_width"`, `"sd_ci_width"` (summaries of the realized interval
widths); `"pct_ci_less_w"` (proportion of intervals narrower than the
target width); `"pct_ci_miss_low"` and `"pct_ci_miss_high"`
(tail-specific empirical non-coverage of the population RMSEA);
`"total_type_I_error"` (overall empirical non-coverage, the sum of the
two tails); and the echoes `"suc_rep"` (number of converged
replications), `"total_N"` (the *N* evaluated), `"df"` (model degrees of
freedom), `"true_rmsea"` (the population RMSEA recovered from fitting
`model` to `Sigma`), `"width"`, and `"conf_level"`. The proportion rows
are on the 0 to 1 scale, not percentages.

## Details

This function implements the sample size planning method proposed in
Kelley and Lai (2011). It uses
[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) to fit the proposed
model to the population covariance matrix, which recovers the population
RMSEA (the model misspecification) and the model degrees of freedom, and
to fit the model to each simulated sample, and it uses
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md) to
construct the confidence interval for the population RMSEA in each
replication. The model is specified in lavaan syntax, so lavaan must be
installed.

Earlier versions of this function used the sem package to fit the model.
The fit is now carried out with lavaan, the structural equation modeling
backend used throughout DMAR. The population RMSEA is read from
[`lavaan::fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html),
which is computed reliably for the large-sample population fit.

## Note

Replications in which lavaan fails to converge, or for which the RMSEA
is undefined, are skipped; the number of converged replications is
reported as `suc_rep`. Increase `G` if many replications fail to
converge.

## References

Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix
that yields a specified minimizer and a specified minimum discrepancy
function value. *Psychometrika, 57*, 357–369.
[doi:10.1007/BF02295424](https://doi.org/10.1007/BF02295424)

Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
root mean square error of approximation: Sample size planning for narrow
confidence intervals. *Multivariate Behavioral Research, 46*, 1–32.
[doi:10.1080/00273171.2011.543027](https://doi.org/10.1080/00273171.2011.543027)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

## See also

[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md),
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# True data generating model: two correlated factors, each measured by
# three standardized indicators. The factor correlation is 0.5 and every
# loading is 0.7. The implied population covariance matrix is assembled
# from the loading matrix, the factor correlation matrix, and the
# residual variances.
Lambda <- matrix(0, 6, 2)
Lambda[1:3, 1] <- 0.7
Lambda[4:6, 2] <- 0.7
Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
Sigma <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
dimnames(Sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))

# Proposed (misspecified) model: a single common factor.
proposed <- "g =~ x1 + x2 + x3 + x4 + x5 + x6"

# The proposed model is fit once at a very large N to recover the
# population RMSEA, the sample size is planned so that the expected width
# of the 95 percent interval is 0.05, and a fresh sample of that size is
# drawn and fit on every replication. Notice that true_rmsea is about
# 0.20, since a single factor is a poor description of two-factor data,
# that the realized widths sit close to the target, and that
# pct_ci_less_w is near one half, which is what planning for the expected
# width delivers. G = 20 keeps the example quick; a reported sensitivity
# study deserves the default G = 200 or more.
set.seed(113)
ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed, Sigma = Sigma,
                          G = 20)
#>  term               value   
#>  mean_rmsea         0.198   
#>  median_rmsea       0.2     
#>  sd_rmsea           0.0176  
#>  mean_ci_width      0.05    
#>  median_ci_width    0.05    
#>  sd_ci_width        9.65e-05
#>  pct_ci_less_w      0.6     
#>  pct_ci_miss_low    0       
#>  pct_ci_miss_high   0.05    
#>  total_type_I_error 0.05    
#>  suc_rep            20      
#>  total_N            695     
#>  df                 9       
#>  true_rmsea         0.204   
#>  width              0.05    
#>  conf_level         0.95    
#> 
#> Confidence level: 95%
```
