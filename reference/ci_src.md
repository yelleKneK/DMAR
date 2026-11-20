# Confidence Interval for a Standardized Regression Coefficient

Computes a confidence interval for a population standardized regression
coefficient, by the standard *t*-based approach or the noncentral *t*
approach. A thin convenience wrapper around
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
which is the general engine; for the coefficient in its raw metric use
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md).

## Usage

``` r
ci_src(
  std_reg_coef_j = NULL,
  se_std_reg_coef_j = NULL,
  N = NULL,
  p = NULL,
  R2_Y_X = NULL,
  R2_j_X_without_j = NULL,
  conf_level = 0.95,
  R2_Y_X_without_j = NULL,
  t_value = NULL,
  reg_coef_j = NULL,
  se_reg_coef_j = NULL,
  sd_Y = NULL,
  sd_X = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL
)
```

## Arguments

- std_reg_coef_j:

  The standardized regression coefficient, or `NULL` (the default) when
  it is given in the units of the data as `reg_coef_j` or when the other
  inputs determine it (see Details)

- se_std_reg_coef_j:

  The standard error of the standardized regression coefficient

- N:

  Sample size

- p:

  The number of predictors

- R2_Y_X:

  The squared multiple correlation coefficient predicting *Y* from the
  *p* predictor variables

- R2_j_X_without_j:

  The squared multiple correlation coefficient predicting the *j*th
  predictor variable (i.e., the predictor of interest) from the
  remaining *p*-1 predictor variables

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate). A level above .999999998, whose tail areas are
  below 1e-9, is refused

- R2_Y_X_without_j:

  The squared multiple correlation coefficient predicting *Y* from the
  *p*-1 predictor variable with the *j*th predictor of interest excluded

- t_value:

  The *t*-value evaluating the null hypothesis that the population
  regression coefficient for the *j*th predictor equals zero

- reg_coef_j:

  The unstandardized regression coefficient

- se_reg_coef_j:

  The standard error of the unstandardized regression coefficient

- sd_Y:

  Standard deviation of *Y*, the dependent variable

- sd_X:

  Standard deviation of *X*, the predictor variable of interest

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

## Value

A 3-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The `term` values are `"ci_lower"`, `"src"` (the
standardized regression coefficient point estimate), and `"ci_upper"`,
so the estimate sits between its confidence limits. The `prob_less` and
`prob_greater` columns report the achieved tail probabilities at each
limit when the noncentral t method is used (`NA` for the estimate row).

## Details

For standardized variables, do not specify the standard deviation of the
variables and input the standardized regression coefficient for
`std_reg_coef_j`.

The interval needs `N`, `p`, and one of these: the coefficient
(`std_reg_coef_j`, or `reg_coef_j` with `sd_Y` and `sd_X`) with its
standard error (`se_std_reg_coef_j` or `se_reg_coef_j`), with `t_value`,
or with `R2_Y_X` and `R2_j_X_without_j`; `t_value` with
`se_std_reg_coef_j`, or with `R2_Y_X` and `R2_j_X_without_j`; or the
three squared multiple correlations `R2_Y_X`, `R2_Y_X_without_j`, and
`R2_j_X_without_j`. The squared multiple correlations fix the magnitude
of the coefficient but not its sign: the positive root is returned, with
a warning, as in
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md).
A call that gives none of these is refused with a condition of class
`dmar_argument_invalid`. The interval is computed by
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
and its conditions name the arguments of this function
(`std_reg_coef_j`, not `reg_coef_j`). A standardized coefficient larger
than 1.1 in absolute value is reported with a warning of class
`dmar_result_qualified`, since a coefficient that large is more often
one in the units of the data.

**Fixed and random predictors.** The standard error of the standardized
coefficient, whether given as `se_std_reg_coef_j`, rescaled from
`se_reg_coef_j` by `sd_X / sd_Y`, or computed from the squared multiple
correlations, treats the sample standard deviations of *Y* and of
\\X_j\\ as fixed constants (Kelley & Maxwell, 2003, pp. 315–316). The
standard deviation of *Y* varies from sample to sample in every design,
and those of the predictors vary as well when the predictors are random,
as they are in most observational studies, so this standard error is
approximate either way and the interval does not have exactly its
nominal coverage. When the values of the predictors are set by the
design, the standard error is, in large samples, too large whenever the
population coefficient is not 0, so the interval covers above its
nominal level. When the predictors are random, the standard error can be
too large or too small, and which it is depends on the correlations
among *Y* and the predictors, not on \\R^2\\ alone. In 2,193
two-predictor structures drawn with each of the three correlations
uniform between -.9 and .9 (structures whose correlation matrix had an
eigenvalue below .02 were set aside), the large-sample standard error of
the first coefficient was too small in a quarter of them: in none of the
743 in which the predictors explained at most 30 percent of the variance
of *Y*, in 8 percent of those explaining 30 to 60 percent, 56 percent of
those explaining 60 to 80, and 79 percent of those explaining 80 to 90,
and in all 70 that explained more than 90 percent. A structure outside
that range can reverse the last: with one predictor correlated .96 with
*Y*, a second correlated .1 with *Y*, and the two uncorrelated (\\R^2 =
.93\\), the large-sample standard error of the first coefficient is 2.2
times its sampling standard deviation, and at *N* = 200 both intervals
covered it in every one of 4,000 samples.

In a simulation of multivariate normal data with random predictors (20
conditions with 1 to 4 predictors and *N* from 50 to 1,000; 4,000 or
20,000 samples per condition, with Monte Carlo standard errors of .001
to .005), the 95 percent interval of this function covered the
population standardized coefficient .961 to .998 of the time with a
single predictor correlated .3 to .7 with *Y*; .964 to .979 for
coefficients of .34 to .52 among two to four predictors correlated
positively with *Y* and with one another and explaining less than 60
percent of its variance; .952 to .956 in the five conditions whose
coefficient was 0 or near it (at most .17 in absolute value); and .955
to .965 for suppressors (a predictor correlated 0 to .1 with *Y* and .6
to .8 with another). With two uncorrelated predictors that explained 72
and 90 percent of the variance of *Y*, it covered .982 and .937. The
central *t* interval with the same standard error (that of
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)
with `noncentral = FALSE`) covered less: .882 to .952 for the
suppressors, and .938 and .881 in the last two conditions. With the
values of the predictors held fixed from sample to sample instead (six
conditions), this function's interval covered .971 to .998 whenever the
coefficient was not 0.

## Note

This function calls upon `ci_reg_coef` in DMAR, but has a different
naming scheme. See
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)
for more details.

To form a confidence interval for the unstandardized regression
coefficient, use `ci_rc`. This function is used to form a confidence
interval for the standardized regression coefficient.

Not all of the values need to be specified, only those that contain all
of the necessary information in order to compute the confidence interval
(options are thus given for the values that need to be specified).

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

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md)

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
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ci_src(std_reg_coef_j = .6707, .1761, N = 30, p = 6, conf_level = .95)
#>  term     value prob_less prob_greater
#>  ci_lower 0.270 0.025     0.975       
#>  src      0.671 <NA>      <NA>        
#>  ci_upper 1.060 0.975     0.025       
#> 
#> Confidence level: 95%
```
