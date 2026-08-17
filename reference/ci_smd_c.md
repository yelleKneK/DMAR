# Confidence Limits for the Standardized Mean Difference Using the Control Group Standard Deviation as the Divisor

Computes the exact noncentral *t*-based confidence limits for the
standardized mean difference that uses the control group standard
deviation as the divisor (Glass's *g*). Standardizing by the control
group alone keeps the scale of the effect anchored in the untreated
population, which matters when the treatment may alter variability as
well as the mean.

## Usage

``` r
ci_smd_c(
  ncp = NULL,
  smd_c = NULL,
  n_C = NULL,
  n_E = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  tol = 1e-09,
  ...
)
```

## Arguments

- ncp:

  The estimated noncentrality parameter, this is generally the observed
  *t*-statistic from comparing the control and experimental group
  (assuming homogeneity of variance)

- smd_c:

  The standardized mean difference (using the control group standard
  deviation in the denominator)

- n_C:

  The sample size for the control group

- n_E:

  The sample size for experimental group

- conf_level:

  The confidence level (1-Type I error rate)

- alpha_lower:

  The Type I error rate for the lower tail

- alpha_upper:

  The Type I error rate for the upper tail

- tol:

  The tolerance of the iterative method for determining the critical
  values

- ...:

  Potentially include parameter for inner functions

## Value

A 3-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` (the lower bound of the confidence interval),
`"smd_c"` (the standardized mean difference standardized by the control
group standard deviation), and `"upper_limit"` (the upper bound).

## Warning

This function uses `conf_limits_nct`, which has as one of its arguments
`tol` (and can be modified with `tol` of the present function). If the
present function fails to converge (i.e., if it runs but does not report
a solution), it is likely that the `tol` value is too restrictive and
should be increased by a factor of 10, but probably by no more than 100.
Running the function `conf_limits_nct` directly will report the actual
probability values of the limits found. This should be done if any
modification to `tol` is necessary in order to ensure acceptable
confidence limits for the noncentral *t* parameter have been achieved.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574.
[doi:10.1177/0013164401614002](https://doi.org/10.1177/0013164401614002)

Glass, G. V. (1976). Primary, secondary, and meta-analysis of research.
*Educational Researcher, 5*, 3–8.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md),
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
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
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ci_smd_c(smd_c = .5, n_C = 100, n_E = 100, conf_level = .95)
#>  term        value
#>  lower_limit 0.213
#>  smd_c       0.5  
#>  upper_limit 0.785
#> 
#> Confidence level: 95%
```
