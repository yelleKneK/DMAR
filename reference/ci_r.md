# Confidence Interval for the Population Multiple Correlation Coefficient

Constructs a confidence interval for the population multiple correlation
coefficient \\\rho = \sqrt{\rho^2}\\ from the sample multiple
correlation coefficient. The confidence interval is for the population
value \\\rho\\; the required input is the corresponding sample value,
the observed sample multiple correlation coefficient \\R\\ (or,
equivalently, the observed *F*-statistic and degrees of freedom). The
interval is obtained by inverting the sampling distribution of the
sample \\R^2\\ and propagating the limits through the monotone (square
root) transform. The right choice of sampling distribution depends on
whether the predictors are treated as random draws from a joint
multivariate normal distribution (the default and the typical case in
the behavioral, educational, and social sciences) or as fixed by design.

## Usage

``` r
ci_r(
  R = NULL,
  df_1 = NULL,
  df_2 = NULL,
  conf_level = 0.95,
  random_predictors = TRUE,
  F_value = NULL,
  N = NULL,
  p = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)

ci_R(
  R = NULL,
  df_1 = NULL,
  df_2 = NULL,
  conf_level = 0.95,
  random_predictors = TRUE,
  F_value = NULL,
  N = NULL,
  p = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- R:

  Observed value of the sample multiple correlation coefficient

- df_1:

  Numerator degrees of freedom

- df_2:

  Denominator degrees of freedom

- conf_level:

  Confidence interval coverage (i.e., 1- Type I error rate); default is
  .95

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

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 3-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The rows are ordered `"lower_limit"`, `"R"` (the sample
multiple correlation coefficient supplied by the user, the point
estimate), and `"upper_limit"`, so the point estimate sits between its
confidence limits. The lower and upper limits are the confidence limits
on the population multiple correlation coefficient \\\rho\\ (square
roots of the corresponding limits on \\\rho^2\\). The `prob_less` and
`prob_greater` columns report the achieved lower-tail and upper-tail
error probabilities at each limit (they are `NA` for the `"R"` estimate
row).

## Details

**Fixed vs.\\ random predictors.** The two regression models give
*different* sampling distributions for the sample \\R^2\\, and so
different confidence intervals on \\\rho\\. Under fixed predictors the
design matrix is treated as constant in hypothetical replications of the
study, and the omnibus \\F\\-statistic follows a noncentral *F* with
\\p\\ and \\N - p - 1\\ degrees of freedom and noncentrality \\\lambda =
N \rho^2 / (1 - \rho^2)\\ (Cohen, 1988); the CI on \\\rho^2\\ is
obtained by inverting that distribution and then taking the square root
(see
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)).
Under random predictors the design matrix is itself a draw from a joint
multivariate normal distribution and the unconditional sampling
distribution of the sample \\R^2\\ is given by Lee (1971); the same Lee
bisection that
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) uses for
the random-predictor CI on \\\rho^2\\ is applied here and the limits are
mapped to \\\rho\\. Gatsonis and Sampson (1989) document the comparison;
in the behavioral, educational, and social sciences predictor variables
are almost always random, so the default is `random_predictors = TRUE`.
Pass `random_predictors = FALSE` for designs in which the predictor
variables are fixed by design.

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

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Lee, Y. S. (1971). Some results on the sampling distribution of the
multiple correlation coefficient. *Journal of the Royal Statistical
Society, Series B, 33*(1), 117–130.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Smithson, M. (2003). *Confidence intervals*. New York, NY: Sage
Publications.

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
interval estimation, power calculations, sample size estimation, and
hypothesis testing in multiple regression. *Behavior Research Methods,
Instruments, & Computers, 24*(4), 581–582.
[doi:10.3758/BF03203611](https://doi.org/10.3758/BF03203611)

## See also

[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
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
ci_r(R = .7071, df_1 =5, df_2 = 50, conf_level = .95, random_predictors=TRUE)
#>  term        value prob_less prob_greater
#>  lower_limit 0.489 0.025     0.975       
#>  R           0.707 <NA>      <NA>        
#>  upper_limit 0.799 0.975     0.025       
#> 
#> Confidence level: 95%
```
