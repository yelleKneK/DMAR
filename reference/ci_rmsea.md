# Confidence Interval for the Population Root Mean Square Error of Approximation

Constructs a confidence interval for the population root mean square
error of approximation (RMSEA), a population badness-of-fit index for
structural equation models. The interval is obtained by inverting the
noncentral chi square distribution of the sample fit function \\T =
(N - 1) \hat F\_{ML}\\ under the model implied covariance structure,
mapping the resulting noncentrality limits to the RMSEA metric (Steiger
& Lind, 1980; Browne & Cudeck, 1993).

## Usage

``` r
ci_rmsea(
  rmsea,
  df,
  N,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL
)
```

## Arguments

- rmsea:

  Observed root mean square error of approximation

- df:

  Degrees of freedom of the model

- N:

  Sample size

- conf_level:

  Desired confidence level (e.g., .90, .95, .99)

- alpha_lower:

  The Type I error rate for the lower tail

- alpha_upper:

  The Type I error rate for the upper tail

## Value

A 3-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` (the lower bound of the confidence interval on the
population RMSEA, truncated at zero by definition), `"rmsea"` (the
observed point estimate), and `"upper_limit"` (the upper bound).

## Details

The RMSEA expresses the badness of model fit per degree of freedom on
the noncentrality scale. Under the noncentral chi square model for the
sample fit statistic, the sample \\T = (N - 1) \hat F\_{ML}\\ has
approximate noncentral chi square distribution with \\df\\ degrees of
freedom and noncentrality parameter \\\lambda = (N - 1) df \cdot
\mathrm{RMSEA}^2\\. The CI on \\\mathrm{RMSEA}^2\\ is obtained by
inverting the noncentral chi square distribution at the requested
confidence level
([`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md)
does the inversion); the bounds are then mapped back to the RMSEA scale
via the square root. When the lower noncentrality limit hits zero
(*i.e.*, the data are compatible with a well-fitting model), the lower
RMSEA limit is truncated at zero because RMSEA is non-negative by
construction.

The 90 percent CI (rather than the usual 95 percent) is the conventional
reporting choice for RMSEA (Browne & Cudeck, 1993) because the upper
limit of the 90 percent CI plays a one-sided role in the test of close
fit (\\H_0: \mathrm{RMSEA} \le 0.05\\). `ci_rmsea` defaults to
`conf_level = 0.95` in line with the rest of the package; pass
`conf_level = 0.90` when the close fit test is the intended use.

## References

Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model
fit. In K. A. Bollen & J. S. Long (Eds.), *Testing structural equation
models* (pp. 136–162). Sage.

Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
root mean square error of approximation: Sample size planning for narrow
confidence intervals. *Multivariate Behavioral Research, 46*, 1–32.
[doi:10.1080/00273171.2011.543027](https://doi.org/10.1080/00273171.2011.543027)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Steiger, J. H., & Lind, J. C. (1980). *Statistically-based tests for the
number of common factors*. Paper presented at the annual Spring meeting
of the Psychometric Society, Iowa City, IA.

## See also

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
# 1. A typical 95 percent CI on RMSEA.
ci_rmsea(rmsea = .055, df = 40, N = 425, conf_level = .95)
#>  term        value 
#>  lower_limit 0.037 
#>  rmsea       0.055 
#>  upper_limit 0.0727
#> 
#> Confidence level: 95%

# 2. The 90 percent CI is the conventional choice when interpretation
#    will follow the Browne and Cudeck (1993) close fit decision rule
#    (the test of H_0: RMSEA <= 0.05 vs. the upper CI limit). Here the
#    upper limit lands below 0.05, which is the close fit threshold
#    Browne and Cudeck recommend.
ci_rmsea(rmsea = .035, df = 40, N = 425, conf_level = .90)
#>  term        value 
#>  lower_limit 0.0147
#>  rmsea       0.035 
#>  upper_limit 0.052 
#> 
#> Confidence level: 90%

# 3. Wider model with smaller N: more uncertainty, wider CI.
ci_rmsea(rmsea = .055, df = 10, N = 100, conf_level = .90)
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#>  term        value
#>  lower_limit 0    
#>  rmsea       0.055
#>  upper_limit 0.129
#> 
#> Confidence level: 90%
```
