# Asymptotic Variance of the Pearson Correlation Coefficient

Computes the asymptotic (large-sample) variance of the sample Pearson
product-moment correlation \\r\\ under bivariate normality (Fisher,
1915) and, optionally, the Bonett-Wright (2000) kurtosis-corrected
variance for non-normal margins. Stand-alone variance utility:
surprisingly absent from CRAN despite being a building block for AIPE
planning, meta-analytic weighting, and Wald-style inference on \\r\\.

## Usage

``` r
var_r(rho, n, kurtosis_x = NULL, kurtosis_y = NULL)
```

## Arguments

- rho:

  Population correlation coefficient. Numeric scalar or vector in \\(-1,
  1)\\.

- n:

  Total sample size on which the Pearson \\r\\ would be computed. Scalar
  or vector.

- kurtosis_x:

  Optional excess kurtosis of the marginal distribution of \\X\\. When
  `kurtosis_x` and `kurtosis_y` are both supplied (along with
  `rho_xy_2x2y` if available), the Bonett-Wright (2000) corrected
  variance is returned in addition to the normal-theory variance.

- kurtosis_y:

  Optional excess kurtosis of \\Y\\; see `kurtosis_x`.

## Value

A `data.frame` with the rows

- `var_r_normal`, the normal-theory asymptotic variance \\(1 - \rho^2)^2
  / (n - 1)\\ (Hotelling, 1953, Section 7),

- `var_r_bonett_wright` (when kurtoses are supplied) the Bonett &
  Wright (2000) kurtosis-corrected variance,

- `var_fisher_z`, the variance of the Fisher \\Z\\ transform, \\1/(n -
  3)\\.

## Details

**Normal-theory variance (default).** Under bivariate normality the
large-sample variance is \$\$\mathrm{Var}(\hat r) \\\approx\\ (1 -
\rho^2)^2 / (n - 1),\$\$ the leading term of the exact moment expansion
(Hotelling, 1953, Section 7; the exact density of \\r\\ is Fisher's,
1915). This is the workhorse variance and is exact in the limit; it is
also what Fisher's \\Z\\ CI
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
uses on the transformed scale.

**Bonett-Wright kurtosis correction.** When the marginals are *not*
normal, the asymptotic variance picks up a kurtosis- dependent
correction (Bonett & Wright, 2000): \$\$\mathrm{Var}(\hat r) \\\approx\\
(1 - \rho^2)^2 / (n - 1) \cdot \bigl(1 + \rho^2 (\gamma_2^{(X)} +
\gamma_2^{(Y)}) / 4\bigr),\$\$ where \\\gamma_2^{(X)}\\,
\\\gamma_2^{(Y)}\\ are the excess kurtoses of the two marginals. The
correction is exact when the joint distribution is elliptical; for
non-elliptical joints it is a first-order approximation. When the
kurtosis arguments are `NULL`, only the normal-theory variance is
returned.

**Connection to Fisher's Z transform.** On the variance- stabilized
scale \\Z = \tanh^{-1}(r)\\ the asymptotic variance is \\1/(n-3)\\
regardless of \\\rho\\ (Fisher, 1921). This is reported alongside the
raw-scale variance because it is the natural working scale for CI
construction
([`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)).

## References

Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
estimating Pearson, Kendall and Spearman correlations. *Psychometrika,
65*(1), 23–28.
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183)

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
var_r(rho = 0.30, n = 50)
#>  term         value 
#>  var_r_normal 0.0169
#>  var_fisher_z 0.0213

# 2. With Bonett-Wright correction for leptokurtic margins
#        (excess kurtosis = 3 for each variable):
var_r(rho = 0.30, n = 50, kurtosis_x = 3, kurtosis_y = 3)
#>  term                value 
#>  var_r_normal        0.0169
#>  var_fisher_z        0.0213
#>  var_r_bonett_wright 0.0192
```
