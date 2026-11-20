# Asymptotic Variance of the Pearson Correlation Coefficient

Computes the large-sample variance of the sample Pearson correlation
under bivariate normality and the variance of its Fisher *Z* transform.
This standalone function is useful for AIPE planning, meta-analytic
weighting, and Wald-style inference on *r*.

## Usage

``` r
var_r(true_cor, N)
```

## Arguments

- true_cor:

  Population correlation coefficient \\\rho\\. Numeric scalar or vector
  in \\(-1, 1)\\.

- N:

  Total sample size on which the Pearson \\r\\ would be computed. Scalar
  or vector.

## Value

A `data.frame` with the rows

- `var_r_normal`, the normal-theory asymptotic variance \\(1 - \rho^2)^2
  / (N - 1)\\ (Hotelling, 1953, Section 7),

- `var_fisher_z`, the variance of the Fisher \\Z\\ transform, \\1/(N -
  3)\\.

## Details

The variance of *r* is evaluated at the population correlation
`true_cor` and the sample size `N`. Substituting a sample correlation
for `true_cor` gives a plug-in estimate of that variance, which carries
the sampling error of the correlation it was computed from. The variance
of the Fisher *Z* transform, \\1/(N - 3)\\, depends on `N` alone.

**Normal-theory variance (default).** Under bivariate normality the
large-sample variance is \$\$\mathrm{Var}(\hat r) \\\approx\\ (1 -
\rho^2)^2 / (N - 1),\$\$ the leading term of the exact moment expansion
(Hotelling, 1953, Section 7; the exact density of \\r\\ is Fisher's,
1915). As \\N\\ grows, \\N\\ times this variance approaches \\(1 -
\rho^2)^2\\, the variance of the limiting distribution of
\\\sqrt{N}(\hat r - \rho)\\. Carried to the Fisher \\Z\\ scale by the
delta method it becomes \\1/(N - 1)\\; the interval of
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
uses Fisher's \\1/(N - 3)\\ on that scale instead (see below).

**Non-normal margins.** The normal-theory variance does not hold when
the margins are not normal: marginal kurtosis can change the
large-sample variance of \\r\\ considerably when \\\rho \ne 0\\ (Bonett
& Wright, 2000, p. 24), and in general the variance depends on every
fourth-order cumulant of the joint distribution (Duncan & Layard, 1973,
Equation 1.1).

**Connection to Fisher's Z transform.** On the variance- stabilized
scale \\Z = \tanh^{-1}(r)\\ the asymptotic variance is \\1/(N-3)\\
regardless of \\\rho\\ (Fisher, 1921). This is reported alongside the
raw-scale variance because it is the natural working scale for CI
construction
([`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)).

## References

Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
estimating Pearson, Kendall and Spearman correlations. *Psychometrika,
65*(1), 23–28.
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183)

Duncan, G. T., & Layard, M. W. J. (1973). A Monte-Carlo study of
asymptotically robust tests for correlation coefficients. *Biometrika,
60*(3), 551–558.
[doi:10.1093/biomet/60.3.551](https://doi.org/10.1093/biomet/60.3.551)

Fisher, R. A. (1915). Frequency distribution of the values of the
correlation coefficient in samples from an indefinitely large
population. *Biometrika, 10*(4), 507–521.

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Hotelling, H. (1953). New light on the correlation coefficient and its
transforms. *Journal of the Royal Statistical Society, Series B, 15*(2),
193–232. (Section 7 gives the exact moments of \\r\\.)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin, 118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

## See also

[`expected_r`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)

Other variance utilities:
[`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
[`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`var_ete()`](https://yelleknek.github.io/DMAR/reference/var_ete.md),
[`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md),
[`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
[`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
[`var_smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Normal-theory variance:
var_r(true_cor = 0.30, N = 50)
#>  term         value
#>  var_r_normal 0.017
#>  var_fisher_z 0.021
```
