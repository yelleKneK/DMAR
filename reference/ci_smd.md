# Confidence Interval for the Standardized Mean Difference (Two Independent Groups)

Constructs an exact-coverage confidence interval for the population
standardized mean difference \\\delta = (\mu_1 - \mu_2)/\sigma\\
(Cohen's *d* when expressed as a sample quantity) for two independent
groups under bivariate normality with equal variances. The interval is
obtained by inverting the noncentral *t* sampling distribution of the
rescaled statistic \\t = \hat d \sqrt{n_1 n_2 / (n_1 + n_2)}\\, which
under the independent groups equal variances model is exactly noncentral
*t* with \\n_1 + n_2 - 2\\ degrees of freedom and noncentrality
parameter \\\lambda = \delta \sqrt{n_1 n_2 / (n_1 + n_2)}\\ (Hedges,
1981). The confidence limits are then rescaled back to the \\\delta\\
metric. This is the same construction Steiger and Fouladi (1997) and
Kelley (2007) describe for noncentral effect size CIs.

## Usage

``` r
ci_smd(
  ncp = NULL,
  smd = NULL,
  n_1 = NULL,
  n_2 = NULL,
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
  *t*-statistic from comparing the two groups and assumes homogeneity of
  variance

- smd:

  The standardized mean difference (using the pooled standard deviation
  in the denominator)

- n_1:

  The sample size for Group 1

- n_2:

  The sample size for Group 2

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

  Allows one to potentially include parameter values for inner functions

## Value

A 3-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` (the lower bound of the confidence interval on the
standardized mean difference), `"smd"` (the point estimate), and
`"upper_limit"` (the upper bound).

## Details

**ncp-input vs. smd-input paths.** The function accepts the effect size
in either of two equivalent metrics: the observed *t*-statistic (via
`ncp`) or the sample standardized mean difference (via `smd`). The two
paths are mathematically equivalent under the equal variances assumption
(since \\t = \hat d \sqrt{n_1 n_2 / (n_1 + n_2)}\\); pick whichever is
easier to obtain. Supply exactly one. Both paths internally call
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)
to invert the noncentral *t* distribution at the specified two-tailed
(or asymmetric, via `alpha_lower` / `alpha_upper`) confidence level.

**Independent vs.\\ paired comparison.** `ci_smd` assumes two
*independent* groups with a common variance. DMAR does not currently
provide a confidence interval for the standardized mean difference in a
paired or within-subject design, whose sampling distribution depends on
the correlation between the paired measurements; applying the
independent groups interval to paired data gives the wrong coverage.
([`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md) is
not a paired interval either; it is the interval for Glass's estimator,
which standardizes the difference between two independent groups by the
control group standard deviation.)

**Bias correction (Hedges' g).** `ci_smd` reports the CI on *d*; if the
bias-corrected *g* is desired, multiply the bounds by the Hedges and
Olkin (1985) correction factor \\J(\nu) = 1 - 3/(4 \nu - 1)\\ (with
\\\nu = n_1 + n_2 - 2\\). Because \\J(\nu)\\ is a constant, the
rescaling preserves coverage.

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

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Hedges, L. V., & Olkin, I. (1985). *Statistical methods for
meta-analysis*. Academic Press.

Kelley, K. (2005). The effects of nonnormal distributions on confidence
intervals around the standardized mean difference: Bootstrap and
parametric confidence intervals. *Educational and Psychological
Measurement, 65*(1), 51–69.
[doi:10.1177/0013164404264850](https://doi.org/10.1177/0013164404264850)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample size planning.
*Evaluation and the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`plot_smd`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
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
# Steiger and Fouladi (1997) example values.
ci_smd(ncp = 2.6, n_1 = 10, n_2 = 10, conf_level = 1 - .05)
#>  term        value
#>  lower_limit 0.195
#>  smd         1.16 
#>  upper_limit 2.1  
#> 
#> Confidence level: 95%
ci_smd(ncp = 2.4, n_1 = 300, n_2 = 300, conf_level = 1 - .05)
#>  term        value 
#>  lower_limit 0.0355
#>  smd         0.196 
#>  upper_limit 0.356 
#> 
#> Confidence level: 95%
```
