# Confidence Interval for an (Unstandardized) Contrast in ANCOVA With One Covariate

Calculates the confidence interval for an unstandardized contrast in the
one-covariate ANCOVA. Two procedures are available through the
`procedure` argument. The default (`"t"`) returns a single
*per-comparison* interval based on the *t* distribution: it gives the
correct \\(1 -\\ `conf_level`\\)\\ coverage for one contrast chosen in
advance, and its standard error includes the \\(\sum c_i \bar X_i)^2 /
SS\_{\mathrm{within}(x)}\\ term that accounts for the covariate
separation between the groups in that one contrast. The
`"bryant_paulson"` procedure instead returns Bryant–Paulson
*simultaneous* (familywise) intervals over a whole family of contrasts
of adjusted means; when selected, `ci_c_ancova` simply forwards its
arguments to
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
and returns that result. See Details for which to use when.

## Usage

``` r
ci_c_ancova(
  psi = NULL,
  adj_means = NULL,
  s_ancova = NULL,
  c_weights,
  n,
  cov_means,
  SSwithin_x,
  conf_level = 0.95,
  procedure = c("t", "bryant_paulson"),
  ...
)
```

## Arguments

- psi:

  The unstandardized contrast of adjusted means

- adj_means:

  The vector that contains the adjusted mean of each group on the
  dependent variable

- s_ancova:

  The standard deviation of the errors from the ANCOVA model (i.e., the
  square root of the mean square error from ANCOVA)

- c_weights:

  The contrast weights

- n:

  Either a single number that indicates the sample size *per group* or a
  vector that contains the sample size of each group

- cov_means:

  A vector that contains the group means of the covariate

- SSwithin_x:

  The sum of squares within groups obtained from the summary table for
  ANOVA on the covariate

- conf_level:

  The desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- procedure:

  The interval procedure, one of `"t"` (the default, a single
  per-comparison interval based on the *t* distribution) or
  `"bryant_paulson"` (Bryant–Paulson simultaneous intervals for a family
  of contrasts of adjusted means). When `"bryant_paulson"` is chosen the
  call is forwarded to
  [`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md);
  see Details.

- ...:

  Allows one to potentially include parameter values for inner
  functions. When `procedure = "bryant_paulson"`, these are passed on to
  [`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
  (for example `num_covariates`, `df`, or `contrast_type`).

## Value

A 3-row `data.frame` with columns `term` and `value` (numeric). The
`term` values are `"lower_limit"` (the lower confidence limit on the
unstandardized ANCOVA contrast), `"psi"` (the unstandardized contrast
point estimate), and `"upper_limit"` (the upper limit).

When `procedure = "bryant_paulson"`, the return value is whatever
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
returns (a table with one row per contrast and columns `contrast`,
`estimate`, `lower_limit`, and `upper_limit`).

## Details

**Per-comparison versus simultaneous.** The two procedures answer
different questions and are not interchangeable. Use the default
`procedure = "t"` when a single contrast was planned in advance: the
interval has exact per-comparison coverage and its width reflects the
covariate adjustment for that specific contrast through the \\(\sum c_i
\bar X_i)^2 / SS\_{\mathrm{within}(x)}\\ term. Use
`procedure = "bryant_paulson"` when several contrasts (for example all
pairwise comparisons of adjusted means) are examined together and the
coverage statement must hold simultaneously across the family: the
Bryant–Paulson generalized studentized range supplies a larger critical
value that controls the familywise error rate and, because the
covariates are random, correctly absorbs the extra sampling uncertainty
from estimating the covariate adjustment (which holds on average over
the covariate distribution, so the per-contrast separation term is not
added again). A per-comparison interval used as if it were simultaneous
understates the family error rate; a simultaneous interval used for one
planned contrast is wider than necessary. See
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
for the full description of the simultaneous procedure and its
arguments.

## Note

Be sure to use the standard deviation and not the error variance for
`s_ancova`, not the square of this value which would come from the
source table (i.e., do not use the variance of the error but rather use
the square root).

If `n` receives a single number, that number is considered as the sample
size *per group*. If `n` receives a vector, the vector is considered as
the sample size of each group.

Be sure to use fractions not the integers to specify `c_weights`. For
example, in an ANCOVA of four groups, if the user wants to compare the
mean of group 1 and 2 with the mean of group 3 and 4, `c_weights` should
be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1).
Make sure the sum of the contrast weights are zero.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md),[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
for the Bryant–Paulson simultaneous procedure

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
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
# scores as the covariate. The means of pretest scores of group 1 to 3 are
# 17, 17.7, and 17.4, respectively, and the adjusted means of groups 1 to 3
# are 7.5, 12, and 14, respectively. The error variance in ANCOVA is 29,
# and the sum of squares within groups from ANOVA on the covariate is 752.5.

# To obtain the confidence interval for adjusted mean of group 1 versus group 2:
ci_c_ancova(adj_means = c(7.5, 12, 14), s_ancova = sqrt(29),
            c_weights = c(1, -1, 0), n = 10,
            cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
#>  term        value
#>  lower_limit -9.46
#>  psi         -4.5 
#>  upper_limit 0.458
#> 
#> Confidence level: 95%

# That interval is the right one for a single contrast planned in advance.
# For the family of all three pairwise comparisons of the adjusted means,
# with coverage that holds simultaneously across the family, select
# procedure = "bryant_paulson"; the call is forwarded to ci_c_ancova_bp().
# That route inverts the Bryant-Paulson distribution numerically to get its
# multiplier, which takes about half a second, so it is shown here rather
# than run.
# ci_c_ancova(adj_means = c(7.5, 12, 14), s_ancova = sqrt(29), n = 10,
#             procedure = "bryant_paulson")
# The simultaneous limits are wider, which is the price of the family
# statement: group 1 against group 2 runs from -10.61 to 1.61, where the
# single planned contrast above ran from -9.47 to 0.47.
```
