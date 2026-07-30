# Confidence Interval for a Standardized Contrast in ANCOVA With One Covariate

Calculate the confidence interval for a standardized contrast in ANCOVA
with one covariate. The standardizer (i.e., the divisor) can be either
the error standard deviation of the ANOVA model (i.e., the model
excluding the covariate) or of the ANCOVA model.

## Usage

``` r
ci_sc_ancova(
  psi = NULL,
  adj_means = NULL,
  s_anova = NULL,
  s_ancova = NULL,
  standardizer = "s_ancova",
  c_weights,
  n,
  cov_means,
  SSwithin_x,
  conf_level = 0.95
)
```

## Arguments

- psi:

  Unstandardized contrast of adjusted means

- adj_means:

  The vector that contains the adjusted mean of each group on the
  dependent variable

- s_anova:

  The standard deviation of the errors from the ANOVA model (i.e., the
  square root of the mean square error from ANOVA)

- s_ancova:

  The standard deviation of the errors from the ANCOVA model (i.e., the
  square root of the mean square error from ANCOVA)

- standardizer:

  Which error standard deviation the user wants to use, the value of
  which can be either `"s_ancova"` or `"s_anova"`

- c_weights:

  The contrast weights (chose weights so that the positive *c*-weights
  sum to 1 and the negative *c*-weights sum to -1; i.e., use fractional
  values not integers).

- n:

  Either a single number that indicates the sample size per group, or a
  vector that contains the sample size of each group

- cov_means:

  A vector that contains the group means of the covariate

- SSwithin_x:

  The sum of squares within groups obtained from the summary table for
  ANOVA on the covariate

- conf_level:

  The desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

## Value

A 3-row `data.frame` with columns `term` and `value` (numeric). The
`term` values are `"lower_limit"` (the lower confidence limit on the
standardized ANCOVA contrast), `"psi"` (the standardized contrast), and
`"upper_limit"` (the upper limit). The divisor used in standardization
(either `"s_anova"` or `"s_ancova"`) is attached as the `"standardizer"`
attribute of the returned data.frame.

## Details

The argument `SSwithin_x` is the sum of squares within groups for the
covariate, taken from the ANOVA source table in which the covariate (not
the outcome) is the dependent variable. Published reports do not always
print this quantity directly. When a report gives the covariate group
means, the group sample sizes, and the *F* statistic from the one-way
ANOVA on the covariate, `SSwithin_x` can be recovered algebraically. The
worked example below follows Lai and Kelley (2012): three groups of
sizes 19, 18, and 19 (so \\N = 56\\) have covariate means 60.08, 57.08,
and 57.97, and the covariate ANOVA reports \\F = 0.756\\ with 2 and 53
degrees of freedom. The sum of squares between groups for the covariate,
computed from the group means and sample sizes, is approximately 88.5,
so the mean square between groups is approximately \\88.5 / 2 = 44.3\\.
Because *F* is the ratio of the mean square between groups to the mean
square within groups, the mean square within groups is approximately
\\44.3 / 0.756 = 58.6\\, and the sum of squares within groups is that
mean square times its degrees of freedom, approximately \\58.6 \times 53
= 3103\\. That recovered value is what you would pass to `SSwithin_x`.
The “Examples” section reproduces this computation in code.

## Note

Be sure to use the standard deviations and not the error variances for
`s_anova` and `s_ancova`, not the squares of these values which would
come from the source tables (i.e., do not use the variance of the errors
but rather use its square root, the standard deviation).

If `n` receives a single number, that number is considered as the sample
size per group. If `n` is assigned to a vector, the vector is considered
as the sample size of each group.

Be sure to use fractional *c*-weights when doing complex contrasts (not
integers) to specify `c_weights`. For example, in an ANCOVA of four
groups, if the user wants to compare the mean of group 1 and 2 with the
mean of group 3 and 4, `c_weights` should be specified as c(0.5, 0.5,
-0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the
contrast weights are zero.

The argument to be assigned to `standardizer` must be either
`"s_ancova"` or `"s_anova"`.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*, 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md)

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
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
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
# Maxwell, Delaney, & Kelley (2027) offer an example that 30 depressive
# individuals are randomly assigned to three groups, 10 in each, and ANCOVA
# is performed on the posttest scores using the participants' pretest
# scores as the covariate. The means of pretest scores of group 1, 2, and 3 are
# 17, 17.7, and 17.4, respectively, whereas the adjusted means of groups 1, 2, and 3
# are 7.5, 12, and 14, respectively. The error variance in ANCOVA is 29 and thus
# 5.385165 is the error standard deviation, with the sum of squares within groups
# from an ANOVA on the covariate is 752.5.

# To obtained the confidence interval for the standardized adjusted mean difference
# between group 1 and 2, using the ANCOVA error standard deviation:
ci_sc_ancova(adj_means = c(7.5, 12, 14), s_ancova = 5.385165, c_weights = c(1, -1, 0),
             n = 10, cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
#>  term        value 
#>  lower_limit -1.73 
#>  psi         -0.836
#>  upper_limit 0.0785
#> 
#> Confidence level: 95%

# Or, with less error in rounding:
ci_sc_ancova(adj_means = c(7.54, 11.98, 13.98), s_ancova = 5.393, c_weights = c(-1, 0, 1),
             n = 10, cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
#>  term        value
#>  lower_limit 0.249
#>  psi         1.19 
#>  upper_limit 2.12 
#> 
#> Confidence level: 95%

# Now, using the standard deviation from ANOVA (and not ANCOVA as above), we have:
ci_sc_ancova(adj_means = c(7.54, 11.98, 13.98), s_anova = 6.294, s_ancova = 5.393,
             c_weights = c(-1, 0, 1),n = 10, cov_means = c(17, 17.7, 17.4),
             SSwithin_x = 752.5, standardizer = "s_anova", conf_level = .95)
#>  term        value
#>  lower_limit 0.214
#>  psi         1.02 
#>  upper_limit 1.82 
#> 
#> Confidence level: 95%

# Recovering SSwithin_x from a covariate ANOVA F when a report does not print
# it directly (see the Details section). This example follows Lai and Kelley
# (2012): three groups of sizes 19, 18, and 19 have covariate means 60.08,
# 57.08, and 57.97, and the one-way ANOVA on the covariate reports F = 0.756.
cov_means_ex  <- c(60.08, 57.08, 57.97)
n_ex          <- c(19, 18, 19)
grand_x       <- sum(n_ex * cov_means_ex) / sum(n_ex)
ss_between_x  <- sum(n_ex * (cov_means_ex - grand_x)^2)
ms_between_x  <- ss_between_x / (length(n_ex) - 1)
ms_within_x   <- ms_between_x / 0.756
SSwithin_x_ex <- ms_within_x * (sum(n_ex) - length(n_ex))

ci_sc_ancova(adj_means = c(63.88, 62.39, 56.48), s_ancova = 20.267,
             c_weights = c(0.5, 0.5, -1), n = n_ex, cov_means = cov_means_ex,
             SSwithin_x = SSwithin_x_ex)
#>  term        value
#>  lower_limit -0.23
#>  psi         0.328
#>  upper_limit 0.884
#> 
#> Confidence level: 95%
```
