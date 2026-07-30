# Confidence Interval for a Regression Coefficient, Raw or Standardized

The general engine behind
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md)
(unstandardized) and
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)
(standardized): computes a confidence interval for a population
regression coefficient by the standard *t*-based approach or the
noncentral *t* approach, in whichever metric the inputs are supplied.

## Usage

``` r
ci_reg_coef(
  b_j,
  SE_b_j = NULL,
  s_Y = NULL,
  s_X = NULL,
  N,
  p,
  R2_Y_X = NULL,
  R2_j_X_without_j = NULL,
  conf_level = 0.95,
  R2_Y_X_without_j = NULL,
  t_value = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  noncentral = FALSE,
  ...
)
```

## Arguments

- b_j:

  Value of the regression coefficient for the *j*th predictor variable

- SE_b_j:

  Standard error for the *j*th predictor variable

- s_Y:

  Standard deviation of *Y*, the dependent variable

- s_X:

  Standard deviation of \\X_j\\, the predictor variable of interest

- N:

  Sample size

- p:

  The number of predictors

- R2_Y_X:

  The squared multiple correlation coefficient predicting `Y` from the
  `p` predictor variables

- R2_j_X_without_j:

  The squared multiple correlation coefficient predicting the *`j`*th
  predictor variable (i.e., the predictor of interest) from the
  remaining `p`-1 predictor variables

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate)

- R2_Y_X_without_j:

  The squared multiple correlation coefficient predicting `Y` from the
  *`p`*-1 predictor variable with the `j`th predictor of interest
  excluded

- t_value:

  The *t*-value evaluating the null hypothesis that the population
  regression coefficient for the `j`th predictor equals zero

- alpha_lower:

  The Type I error rate for the lower confidence interval limit

- alpha_upper:

  The Type I error rate for the upper confidence interval limit

- noncentral:

  `TRUE` or `FALSE`, specifying whether or not the noncentral approach
  to confidence intervals should be used

- ...:

  Optional additional specifications for nested functions

## Value

A 3-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The rows are ordered `"lower_limit"`, `"reg_coef"` (the
regression coefficient point estimate), and `"upper_limit"`, so the
point estimate sits between its confidence limits. The lower and upper
rows give the confidence limits on the regression coefficient. The
`prob_less` and `prob_greater` columns report the achieved tail
probabilities at each limit when the noncentral t method is used (they
are `NA` for the `"reg_coef"` estimate row).

## Details

For standardized variables, do not specify the standard deviation of the
variables and input the standardized regression coefficient for `b_j`.

When `b_j` is reconstructed from squared multiple correlations (that is,
from `R2_Y_X`, `R2_Y_X_without_j`, and `R2_j_X_without_j` rather than a
supplied `b_j`, `SE_b_j`, or `t_value`), only the magnitude of the
coefficient is identifiable; its sign is not. The positive root is
returned and a warning is issued. If the coefficient is negative, negate
the point estimate and swap and negate the confidence limits, or supply
`b_j` directly.

## Note

Not all of the values need to be specified, only those that contain all
of the necessary information in order to compute the confidence interval
(options are thus given for the values that need to be specified).

The function `ci_rc` in DMAR also calculates the confidence interval for
the population (unstandardized) regression coefficient. The function
`ci_src` also calculates the confidence interval for the population
(standardized) regression coefficient. These two functions perform the
same tasks as `ci_reg_coef` does and are preferred to it because of
simpler arguments.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

Smithson, M. (2003). *Confidence intervals*. New York, NY: Sage
Publications.

## See also

[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)

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
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
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
ci_reg_coef(b_j = 0.61319, SE_b_j = 0.16098, N = 30, p = 6)
#>  term        value prob_less prob_greater
#>  lower_limit 0.28  0.025     0.975       
#>  reg_coef    0.613 <NA>      <NA>        
#>  upper_limit 0.946 0.975     0.025       
#> 
#> Confidence level: 95%
```
