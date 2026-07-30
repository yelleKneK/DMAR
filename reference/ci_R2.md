# Confidence Interval for the Population Squared Multiple Correlation Coefficient

Constructs a confidence interval for the population squared multiple
correlation coefficient \\\rho^2\\ by inverting the sampling
distribution of the sample \\R^2\\. The confidence interval is for the
population value \\\rho^2\\; the required input is the corresponding
sample value, the observed sample squared multiple correlation
coefficient \\R^2\\ (or, equivalently, the observed *F*-statistic and
degrees of freedom). The right choice of sampling distribution, and so
the right interval, depends on whether the predictors are treated as
random draws from a joint multivariate normal distribution (the default
and the typical case in the behavioral, educational, and social
sciences) or as fixed by design (planned dosing levels, factorial
covariates, etc.). The function selects the sampling distribution via
`random_predictors` and inverts the corresponding noncentral
distribution; the construction is the regression analogue of a
noncentral distribution based CI on a standardized effect size (Steiger
& Fouladi, 1992; Kelley, 2007).

## Usage

``` r
ci_R2(
  R2 = NULL,
  df_1 = NULL,
  df_2 = NULL,
  conf_level = 0.95,
  random_predictors = TRUE,
  F_value = NULL,
  N = NULL,
  p = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  tol = 1e-09
)
```

## Arguments

- R2:

  Observed value of the sample squared multiple correlation coefficient

- df_1:

  Numerator degrees of freedom

- df_2:

  Denominator degrees of freedom

- conf_level:

  Confidence interval coverage; 1-Type I error rate

- random_predictors:

  Whether or not the predictor variables are random or fixed (random is
  default)

- F_value:

  Obtained *F*-value

- N:

  Sample size

- p:

  Number of predictors

- alpha_lower:

  Type I error for the lower confidence limit

- alpha_upper:

  Type I error for the upper confidence limit

- tol:

  Tolerance for iterative convergence

## Value

A 3-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The rows are ordered `"lower_limit"` (lower confidence
limit on the population \\\rho^2\\), `"R2"` (the sample squared multiple
correlation coefficient supplied by the user, the point estimate), and
`"upper_limit"` (upper confidence limit on the population \\\rho^2\\),
so the point estimate sits between its confidence limits. The
`prob_less` and `prob_greater` columns report the achieved lower-tail
and upper-tail error probabilities at each limit (they are `NA` for the
`"R2"` estimate row). For random-predictor mode
(`random_predictors = TRUE`) the limits are computed via the Lee (1971)
bisection over the multiple-correlation sampling distribution; for
fixed-predictor mode they are computed by inversion of the noncentral
*F* distribution.

## Details

**Fixed vs.\\ random predictors.** The two regression models give
*different* sampling distributions for the sample \\R^2\\, and so
different confidence intervals. Under fixed predictors the design matrix
is treated as constant in hypothetical replications of the study, and
the omnibus \\F\\-statistic \\F = (R^2 / p) / ((1 - R^2) / (N - p -
1))\\ follows a noncentral *F* with \\p\\ and \\N - p - 1\\ degrees of
freedom and noncentrality \\\lambda = N \rho^2 / (1 - \rho^2)\\ (Cohen,
1988); the CI is obtained by inverting that distribution at the supplied
confidence level (see
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)).
Under random predictors the design matrix is itself a draw from a joint
multivariate normal distribution and the unconditional sampling
distribution of the sample \\R^2\\ is given by Lee (1971); `ci_R2` uses
the Lee (1971) bisection (the same construction Algina and Olejnik 2000
implemented in SAS) to invert that distribution. Gatsonis and Sampson
(1989) document the comparison and show that treating random predictors
as fixed tends to over-state precision (and so under-state the CI
width); the discrepancy is modest at moderate to large \\N\\ but
non-trivial at small \\N\\ with moderate-to-large effects. In the
behavioral, educational, and social sciences predictor variables are
almost always random, so the default is `random_predictors = TRUE`; pass
`random_predictors = FALSE` for designs in which the predictor variables
are fixed by design.

## References

Algina, J. & Olejnik, S. (2000). Determining sample size for accurate
estimation of the squared multiple correlation coefficient.
*Multivariate Behavioral Research, 35*, 119–137.
[doi:10.1207/s15327906mbr3501_5](https://doi.org/10.1207/s15327906mbr3501_5)

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106*(3), 516–524.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Lee, Y. S. (1971). Some results on the sampling distribution of the
multiple correlation coefficient. *Journal of the Royal Statistical
Society, Series B, 33*(1), 117–130.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Smithson, M. (2003). *Confidence intervals*. New York, NY: Sage
Publications.

Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
interval estimation, power calculations, sample size estimation, and
hypothesis testing in multiple regression. *Behavior Research Methods,
Instruments, & Computers, 24*(4), 581–582.
[doi:10.3758/BF03203611](https://doi.org/10.3758/BF03203611)

## See also

[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

Other confidence intervals for effect sizes:
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# For random predictor variables.
ci_R2(R2 = .25, N = 100, p = 5, conf_level = .95, random_predictors = TRUE)
#>  term        value  prob_less prob_greater
#>  lower_limit 0.0801 0.025     0.975       
#>  R2          0.25   <NA>      <NA>        
#>  upper_limit 0.372  0.975     0.025       
#> 
#> Confidence level: 95%

ci_R2(F_value = 6.266667, N = 100, p = 5, conf_level = .95, random_predictors = TRUE)
#>  term        value  prob_less prob_greater
#>  lower_limit 0.0801 0.025     0.975       
#>  R2          0.25   <NA>      <NA>        
#>  upper_limit 0.372  0.975     0.025       
#> 
#> Confidence level: 95%

# For fixed predictor variables.
ci_R2(R2 = .25, N = 100, p = 5, conf_level = .95, random_predictors = FALSE)
#>  term        value  prob_less prob_greater
#>  lower_limit 0.0808 0.025     0.975       
#>  R2          0.25   <NA>      <NA>        
#>  upper_limit 0.354  0.975     0.025       
#> 
#> Confidence level: 95%

ci_R2(F_value = 6.266667, N = 100, p = 5, conf_level = .95, random_predictors = FALSE)
#>  term        value  prob_less prob_greater
#>  lower_limit 0.0808 0.025     0.975       
#>  R2          0.25   <NA>      <NA>        
#>  upper_limit 0.354  0.975     0.025       
#> 
#> Confidence level: 95%

# One sided confidence intervals when predictors are random.
ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = .05, alpha_upper = 0,
      conf_level = NULL, random_predictors = TRUE)
#>  term        value  prob_less prob_greater
#>  lower_limit 0.0992 0.05      0.95        
#>  R2          0.25   <NA>      <NA>        
#>  upper_limit 1      1         0           
#> 
#> Confidence level: 95%

ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = 0, alpha_upper = .05,
      conf_level = NULL, random_predictors = TRUE)
#>  term        value prob_less prob_greater
#>  lower_limit 0     0         1           
#>  R2          0.25  <NA>      <NA>        
#>  upper_limit 0.347 0.95      0.05        
#> 
#> Confidence level: 95%

# One sided confidence intervals when predictors are fixed.
ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = .05, alpha_upper = 0,
      conf_level = NULL, random_predictors = FALSE)
#>  term        value prob_less prob_greater
#>  lower_limit 0.1   0.05      0.95        
#>  R2          0.25  <NA>      <NA>        
#>  upper_limit 1     1         0           
#> 
#> Confidence level: 95%

ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = 0, alpha_upper = .05,
      conf_level = NULL, random_predictors = FALSE)
#>  term        value prob_less prob_greater
#>  lower_limit 0     0         1           
#>  R2          0.25  <NA>      <NA>        
#>  upper_limit 0.332 0.95      0.05        
#> 
#> Confidence level: 95%
```
