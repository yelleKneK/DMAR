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
  reg_coef_j = NULL,
  se_reg_coef_j = NULL,
  sd_Y = NULL,
  sd_X = NULL,
  N,
  p,
  R2_Y_X = NULL,
  R2_j_X_without_j = NULL,
  conf_level = 0.95,
  R2_Y_X_without_j = NULL,
  t_value = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  noncentral = FALSE
)
```

## Arguments

- reg_coef_j:

  Value of the regression coefficient for the *j*th predictor variable,
  or `NULL` (the default) when the other inputs determine it: `t_value`
  with `R2_Y_X` and `R2_j_X_without_j`, or the three squared multiple
  correlations, which fix its magnitude but not its sign (see Details)

- se_reg_coef_j:

  Standard error for the *j*th predictor variable

- sd_Y:

  Standard deviation of *Y*, the dependent variable

- sd_X:

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
  Type I error rate). A level above .999999998, whose tail areas are
  below 1e-9, is refused

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

  The Type I error rate for the upper confidence interval limit. The
  table records the coverage, `1 - alpha_lower - alpha_upper`, as its
  level and, when the two tail areas differ, the areas themselves, which
  the line beneath the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound. A
  nonzero tail area below 1e-9 is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- noncentral:

  `TRUE` or `FALSE`, specifying whether or not the noncentral approach
  to confidence intervals should be used

## Value

A 3-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The rows are ordered `"ci_lower"`, `"reg_coef"` (the
regression coefficient point estimate), and `"ci_upper"`, so the point
estimate sits between its confidence limits. The lower and upper rows
give the confidence limits on the regression coefficient. The
`prob_less` and `prob_greater` columns report the achieved tail
probabilities at each limit when the noncentral t method is used (they
are `NA` for the `"reg_coef"` estimate row).

## Details

For standardized variables, do not specify the standard deviation of the
variables and input the standardized regression coefficient for
`reg_coef_j`.

**Exactness.** For a coefficient in the units of the data, the central
interval (`noncentral = FALSE`, the default) is exact under normal
errors given the values of the predictors, and so for random predictors
as well. The noncentral interval multiplies the confidence limits for
the noncentrality parameter by the estimated standard error, while the
noncentrality parameter is the coefficient divided by its population
standard error, so that interval is approximate on either scale. In a
simulation with one predictor held fixed and normal errors (4,000
samples per condition, *N* of 15, 40, and 200, Monte Carlo standard
errors of .001 to .004), the 95 percent noncentral interval for a
coefficient in the units of the data covered .944 to .955 when the
noncentrality parameter was 1, and from .961 (*N* = 200) to .993 (*N* =
15) when it was 7; the central interval covered .943 to .953 in every
condition.

**Standardized data.** For a coefficient of standardized data, or one
rescaled by the sample standard deviations, the standard error treats
those standard deviations as fixed constants (Kelley & Maxwell, 2003,
pp. 315–316). The standard deviation of *Y* varies from sample to sample
in every design, and those of the predictors vary as well when the
predictors are random, so this standard error is approximate either way
and the interval does not have exactly its nominal coverage. In the
simulation of multivariate normal data described in
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md), the
central 95 percent interval (the default here) covered the population
standardized coefficient .957 to .993 of the time with a single random
predictor correlated .3 to .7 with *Y*, .959 to .967 for coefficients of
.34 to .52 among two to four random predictors correlated positively
with *Y* and with one another and explaining less than 60 percent of its
variance, and .952 to .956 in the conditions whose coefficient was 0 or
near it, but only .882 to .952 for suppressors, and .938 and .881 with
two uncorrelated predictors that explained 72 and 90 percent of the
variance of *Y*. With the values of the predictors held fixed it covered
.955 to .989 whenever the coefficient was not 0. The noncentral
interval, the one
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)
reports, covered more in every condition except those whose coefficient
was 0 or near it, where the two agreed within Monte Carlo error.

When `reg_coef_j` is reconstructed from squared multiple correlations
(that is, from `R2_Y_X`, `R2_Y_X_without_j`, and `R2_j_X_without_j`
rather than a supplied `reg_coef_j`, `se_reg_coef_j`, or `t_value`),
only the magnitude of the coefficient is identifiable; its sign is not.
The positive root is returned and a warning is issued. If the
coefficient is negative, negate the point estimate and swap and negate
the confidence limits, or supply `reg_coef_j` directly.

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

Smithson, M. (2003). *Confidence intervals*. Thousand Oaks, CA: Sage
Publications.

## See also

[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_ratio_of_means()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md),
[`ci_ratio_of_means_paired()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md),
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
ci_reg_coef(reg_coef_j = 0.61319, se_reg_coef_j = 0.16098, N = 30, p = 6)
#>  term     value prob_less prob_greater
#>  ci_lower 0.280 0.025     0.975       
#>  reg_coef 0.613 <NA>      <NA>        
#>  ci_upper 0.946 0.975     0.025       
#> 
#> Confidence level: 95%
```
