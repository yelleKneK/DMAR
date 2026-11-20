# Exact Expected Value of the Sample Partial Correlation

Computes \\\mathrm{E}\[r\_{XY \cdot Z_1 \cdots Z\_{p-1}} \mid \rho, N,
p\]\\, the exact expected value of the sample partial Pearson
correlation coefficient under multivariate normality, applying the
Olkin-Pratt (1958) / Hotelling (1953) result for the simple Pearson
correlation to the partial-correlation setting, an extension Olkin and
Pratt (1958, Section 2.3) give themselves. Like the simple *r*, the
partial *r* is downward-biased as an estimator of the population partial
correlation \\\rho\_{XY \cdot Z}\\, with the magnitude of the bias
growing in the number of predictors \\p\\ and shrinking in the sample
size \\N\\.

## Usage

``` r
expected_partial_r(true_cor, N, p)
```

## Arguments

- true_cor:

  Population partial correlation \\\rho\_{XY \cdot Z_1 \cdots
  Z\_{p-1}}\\. Scalar or vector in \\\[-1, 1\]\\.

- N:

  Total sample size; \\N - p \ge 3\\ is required for well-conditioned
  computation. Scalar or vector.

- p:

  Number of predictors in the regression of *Y* on *X* and the variables
  partialed out, so `p - 1` variables are partialed out and `p` is at
  least 2. A categorical variable with *a* levels counts as *a* - 1
  predictors. It is the `p` of
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
  [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md), and
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md):
  the correlation concerns one predictor of a regression with `p`
  predictors. In the notation of a model comparison, \\q\\ is the number
  of predictors removed from the full model to give the restricted model
  (\\q = p\_{full} - p\_{restricted}\\, the predictors being tested),
  and the number of variables partialed out is \\p - q\\; this
  correlation concerns one predictor, so \\q = 1\\ and \\p - 1\\
  variables are partialed out. Scalar or vector.

## Value

A `data.frame` with one row per (`true_cor`, `N`, `p`) input and the
columns `true_cor`, `N`, `p`, `expected_partial_r`, `bias`, and
`relative_bias`.

## Details

The formula is the same as for the simple *r* but with the degrees of
freedom reduced from \\N\\ to \\N - (p - 1)\\: \$\$\mathrm{E}\[r\_{XY
\cdot Z} \mid \rho\_{XY \cdot Z}, N, p\] \\=\\ \rho\_{XY \cdot Z}
\\\cdot\\ {}\_2F_1\\\left(\tfrac{1}{2},\\ \tfrac{1}{2};\\ \tfrac{N - p +
2}{2};\\ \rho\_{XY \cdot Z}^{\\2}\right) \\\cdot\\
\frac{\Gamma\\\left(\tfrac{N - p + 1}{2}\right)^2}
{\Gamma\\\left(\tfrac{N - p}{2}\right)\\ \Gamma\\\left(\tfrac{N - p +
2}{2}\right)}.\$\$

`expected_partial_r()` is especially useful at the *design stage*:
sample size plans that assume \\\rho\_{XY \cdot Z}\\ as the value to be
observed will under-deliver on the expected width or power of a CI on
the partial correlation by an amount that grows in \\p\\.

**Generalization of Olkin-Pratt.** Under multivariate normality, the
partial correlation \\r\_{XY \cdot Z}\\ computed from a sample of size
\\N\\ has the same sampling distribution as a simple Pearson correlation
from a sample of size \\N - (p - 1)\\ (Anderson, 2003, Theorem 4.3.5,
Section 4.3.2, p. 143, who credits the derivation to Fisher). The bias
formula for the simple *r* (Hotelling, 1953; Olkin & Pratt, 1958)
therefore applies directly to the partial *r* with the substitution \\N
\to N - (p - 1)\\. The Olkin-Pratt (1958) unbiased estimator extends in
the same way: given an observed \\r\_{XY \cdot Z}\\, the unbiased
estimator of \\\rho\_{XY \cdot Z}\\ is \\r \cdot {}\_2F_1(1/2, 1/2; (N -
p - 1)/2; 1 - r^2)\\ (Olkin & Pratt, 1958, Section 2.3, where the
partial correlation is shown to have the simple correlation's density at
the reduced sample size, so their estimator applies with that
substitution).

**Magnitude of the bias.** To leading order, from Hotelling's (1953)
expansion of the simple-correlation expectation applied at the reduced
sample size: \$\$\rho\_{XY \cdot Z} - \mathrm{E}\[r\_{XY \cdot Z}\]
\\\approx\\ \rho\_{XY \cdot Z} (1 - \rho\_{XY \cdot Z}^2) / \[2 (N -
p)\],\$\$ matching the sign convention of the returned `bias` column,
\\\rho\_{XY \cdot Z} - \mathrm{E}\[r\_{XY \cdot Z}\]\\, which is
positive for positive \\\rho\_{XY \cdot Z}\\. For \\\rho\_{XY \cdot Z} =
0.4\\, \\N = 30\\, \\p = 3\\, the exact bias is about \\+0.00624\\; for
\\N = 30\\, \\p = 11\\, about \\+0.00888\\.

**Tuning the series convergence.** The underlying \\{}\_2F_1\\ series is
summed by forward recurrence and stops when the relative contribution
falls below `getOption("dmar.expected_r.tol", 1e-15)`, with a maximum of
`getOption("dmar.expected_r.max_iter", 5000L)` terms. These are the same
options as for
[`expected_r`](https://yelleknek.github.io/DMAR/reference/expected_r.md);
users almost never need to change them.

## References

Anderson, T. W. (2003). *An introduction to multivariate statistical
analysis* (3rd ed.), Sections 4.2 and 4.3. Wiley. (Theorem 4.3.5,
Section 4.3.2, p. 143: the sample partial correlation based on \\N\\
observations, with \\k\\ variables partialed out, has cdf \\F\[r \mid
N - k, \rho\]\\; Anderson writes the count of partialed variables as
\\p - q\\ in his notation, and here it is \\p - 1\\.)

Hotelling, H. (1953). New light on the correlation coefficient and its
transforms. *Journal of the Royal Statistical Society, Series B, 15*(2),
193–232.

Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain
correlation coefficients. *The Annals of Mathematical Statistics,
29*(1), 201–211.

## See also

[`expected_r`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. A single value: true_cor = 0.4, N = 30, three variables partialed out
#    (p = 4 predictors: X and the three variables).
expected_partial_r(true_cor = 0.4, N = 30, p = 4)
#>  true_cor N  p expected_partial_r bias  relative_bias
#>  0.400    30 4 0.394              0.006 0.016        

# 2. Bias grows with p at fixed (true_cor, N):
expected_partial_r(true_cor = 0.4, N = 30, p = c(2, 3, 6, 11, 21))
#>  true_cor N  p  expected_partial_r bias  relative_bias
#>  0.400    30 2  0.394              0.006 0.015        
#>  0.400    30 3  0.394              0.006 0.016        
#>  0.400    30 6  0.393              0.007 0.018        
#>  0.400    30 11 0.391              0.009 0.022        
#>  0.400    30 21 0.381              0.019 0.047        

# 3. Partialing one variable (p = 2) at N = 11 is distributed exactly as
#    a simple correlation at N = 10 (Anderson, 2003, Theorem 4.3.5):
expected_partial_r(true_cor = 0.5, N = 11, p = 2)$expected_partial_r
#> [1] 0.4786588
expected_r(true_cor = 0.5, N = 10)$expected_r
#> [1] 0.4786588
```
