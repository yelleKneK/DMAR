# Confidence Interval for the Population Root Mean Square Error of Approximation

Constructs a confidence interval for the population root mean square
error of approximation (RMSEA), a population badness-of-fit index for
structural equation models, from a reported RMSEA, its degrees of
freedom, and the sample size. The interval is obtained by inverting the
noncentral chi square distribution of the test statistic \\T =
a\hat{F}\\, where \\\hat{F}\\ is the minimized fit function and \\a\\
the sample size multiplier, and mapping the resulting noncentrality
limits to the RMSEA metric (Steiger & Lind, 1980; Browne & Cudeck,
1993).

## Usage

``` r
ci_rmsea(
  rmsea,
  df,
  N,
  conf_level = 0.9,
  alpha_lower = NULL,
  alpha_upper = NULL,
  sample_size_convention = c("N-1", "N")
)
```

## Arguments

- rmsea:

  Observed root mean square error of approximation

- df:

  Degrees of freedom of the model, greater than 0 (a model with zero
  degrees of freedom is saturated, and its RMSEA is undefined)

- N:

  Sample size: greater than 1 under the default
  `sample_size_convention`, and greater than 0 under `"N"`

- conf_level:

  Desired confidence level, between 0 and 1. The default, 0.90, gives
  the conventional interval for the RMSEA, whose two limits decide the
  tests of close and not close fit at .05 (see Details).

- alpha_lower:

  The Type I error rate for the lower tail

- alpha_upper:

  The Type I error rate for the upper tail. The table records the
  coverage, `1 - alpha_lower - alpha_upper`, as its level and, when the
  two tail areas differ, the areas themselves, which the line beneath
  the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound.

- sample_size_convention:

  The multiplier \\a\\ with which the reported RMSEA was computed:
  `"N-1"` (the default), so \\a = N - 1\\, the convention of Browne and
  Cudeck (1993), of Kelley and Lai (2011), and of
  [`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md);
  or `"N"`, so \\a = N\\, the convention of lavaan's default normal
  likelihood and of Mplus.

## Value

A 3-row `data.frame` with columns `term` and `value`. The `term` values
are `"ci_lower"` (the lower bound of the confidence interval on the
population RMSEA, truncated at zero by definition), `"rmsea"` (the
observed point estimate), and `"ci_upper"` (the upper bound).

## Details

The RMSEA expresses the badness of model fit per degree of freedom on
the noncentrality scale. Under the noncentral chi square model for the
sample fit statistic, \\T = a\hat{F}\\ has approximate noncentral chi
square distribution with \\df\\ degrees of freedom and noncentrality
parameter \\\lambda = a \cdot df \cdot \mathrm{RMSEA}^2\\. The function
rebuilds \\T = a \cdot df \cdot \mathrm{RMSEA}^2 + df\\ from the
reported RMSEA, obtains the interval for \\\lambda\\ by inverting the
noncentral chi square distribution at the requested confidence level
([`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md)
does the inversion), and maps its limits back to the RMSEA scale as
\\\sqrt{\lambda / (a \cdot df)}\\. When the lower noncentrality limit
hits zero (*i.e.*, the data are compatible with a well-fitting model),
the lower RMSEA limit is truncated at zero because RMSEA is nonnegative
by construction.

The multiplier must be the one with which the RMSEA was computed. An
RMSEA from lavaan under its default likelihood was computed with \\N\\;
given to the default `"N-1"`, it gives an interval that matches neither
convention, off by an amount of order \\1/N\\. Two further cases need
the model's chi square statistic rather than its RMSEA, and
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
takes it: an RMSEA truncated at zero, which cannot distinguish the
statistics below \\df\\ that produce it, and the RMSEA of a multiple
group model, which carries the factor \\\sqrt{G}\\ for \\G\\ groups
(Steiger, 1998); this function assumes a single group.

The 90 percent CI (rather than the usual 95 percent) is the conventional
reporting choice for the RMSEA, from the .05 and .95 percentile points
of Steiger and Lind (1980) through Browne and Cudeck (1993), because a
\\100(1 - 2\alpha)\\ percent interval decides the two one-sided tests at
level \\\alpha\\, so that the limits of the 90 percent interval settle
them at level .05 (MacCallum, Browne, & Sugawara, 1996, p. 136 and Table
1): the hypothesis of close fit, \\H_0: \mathrm{RMSEA} \le \epsilon_0\\,
is rejected when the lower limit exceeds \\\epsilon_0\\, and the
hypothesis of not close fit, \\H_0: \mathrm{RMSEA} \ge \epsilon_0\\, is
rejected, so that fit within \\\epsilon_0\\ is supported, when the upper
limit falls below \\\epsilon_0\\. A lower limit of zero likewise says
that the test of exact fit does not reject at .05. `ci_rmsea` therefore
defaults to `conf_level = 0.90`, as every RMSEA interval in DMAR does.
`conf_level = 0.95`, the level of the package's other intervals, gives
the 95 percent interval, whose limits decide the one-sided tests at
.025. The printed table states the level beneath it. See
[`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md)
for the test of not close fit with its *p*-value and verdict.

## References

Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model
fit. In K. A. Bollen & J. S. Long (Eds.), *Testing structural equation
models* (pp. 136–162). Sage.

Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
root mean square error of approximation: Sample size planning for narrow
confidence intervals. *Multivariate Behavioral Research, 46*, 1–32.
[doi:10.1080/00273171.2011.543027](https://doi.org/10.1080/00273171.2011.543027)

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA
fit index. *Structural Equation Modeling, 5*(4), 411–419.
[doi:10.1080/10705519809540115](https://doi.org/10.1080/10705519809540115)

Steiger, J. H., & Lind, J. C. (1980). *Statistically-based tests for the
number of common factors*. Paper presented at the annual Spring meeting
of the Psychometric Society, Iowa City, IA.

## See also

[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
[`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md)

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
# 1. The conventional 90 percent CI on RMSEA, the default. Its upper
#    limit, 0.052, lies just above 0.05, so the test of not close fit at
#    .05 does not reject an RMSEA of 0.05 or more, and fit within 0.05 is
#    not established even though the point estimate sits well below it.
ci_rmsea(rmsea = .035, df = 40, N = 425)
#>  term     value
#>  ci_lower 0.015
#>  rmsea    0.035
#>  ci_upper 0.052
#> 
#> Confidence level: 90%

# 2. The 95 percent CI is wider.
ci_rmsea(rmsea = .035, df = 40, N = 425, conf_level = .95)
#>  term     value
#>  ci_lower 0.008
#>  rmsea    0.035
#>  ci_upper 0.055
#> 
#> Confidence level: 95%

# 3. An RMSEA reported by lavaan under its default likelihood, which
#    forms the statistic with N.
ci_rmsea(rmsea = .0921, df = 24, N = 301, sample_size_convention = "N")
#>  term     value
#>  ci_lower 0.071
#>  rmsea    0.092
#>  ci_upper 0.114
#> 
#> Confidence level: 90%

# 4. Fewer degrees of freedom and a smaller N: a wider CI.
ci_rmsea(rmsea = .055, df = 10, N = 100)
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#>  term     value
#>  ci_lower 0.000
#>  rmsea    0.055
#>  ci_upper 0.129
#> 
#> Confidence level: 90%
```
