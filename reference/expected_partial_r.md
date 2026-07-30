# Exact Expected Value of the Sample Partial Correlation

Computes \\\mathrm{E}\[r\_{XY \cdot Z_1 \cdots Z_J} \mid \rho, n, J\]\\,
the exact expected value of the sample partial Pearson correlation
coefficient under multivariate normality, generalizing the Olkin-Pratt
(1958) / Hotelling (1953) result for the simple Pearson correlation to
the partial-correlation setting (Olkin & Finn, 1995). Like the simple
*r*, the partial *r* is downward-biased as an estimator of the
population partial correlation \\\rho\_{XY \cdot Z}\\, with the
magnitude of the bias growing in the number of controls \\J\\ and
shrinking in the sample size \\n\\.

## Usage

``` r
expected_partial_r(rho, n, J)
```

## Arguments

- rho:

  Population partial correlation \\\rho\_{XY \cdot Z_1 \cdots Z_J}\\.
  Scalar or vector in \\\[-1, 1\]\\.

- n:

  Total sample size; \\n - J - 1 \ge 3\\ is required for
  well-conditioned computation.

- J:

  Number of variables partialled out (the count of \\Z_1, \ldots,
  Z_J\\); at least 1.

## Value

A `data.frame` with one row per (`rho`, `n`, `J`) input and the columns
`rho`, `n`, `J`, `expected_partial_r`, `bias`, and `relative_bias`.

## Details

The formula is the same as for the simple *r* but with the degrees of
freedom reduced from \\n\\ to \\n - J\\: \$\$\mathrm{E}\[r\_{XY \cdot Z}
\mid \rho\_{XY \cdot Z}, n, J\] \\=\\ \rho\_{XY \cdot Z} \\\cdot\\
{}\_2F_1\\\left(\tfrac{1}{2},\\ \tfrac{1}{2};\\ \tfrac{n - J - 1}{2};\\
\rho\_{XY \cdot Z}^{\\2}\right) \\\cdot\\ \frac{\Gamma\\\left(\tfrac{n -
J - 1}{2}\right)^2} {\Gamma\\\left(\tfrac{n - J - 2}{2}\right)\\
\Gamma\\\left(\tfrac{n - J}{2}\right)}.\$\$

`expected_partial_r()` is especially useful at the *design stage*:
sample size plans that assume \\\rho\_{XY \cdot Z}\\ as the value to be
observed will under-deliver on the expected width or power of a CI on
the partial correlation by an amount that grows in \\J\\.

**Generalization of Olkin-Pratt.** Under multivariate normality, the
partial correlation \\r\_{XY \cdot Z}\\ computed from a sample of size
\\n\\ has the same sampling distribution as a simple Pearson correlation
from a sample of size \\n - J\\ (Anderson, 2003, Section 4.3). The bias
formula for the simple *r* (Hotelling, 1953; Olkin & Pratt, 1958)
therefore applies directly to the partial *r* with the substitution \\n
\to n - J\\. The Olkin-Pratt (1958) unbiased estimator extends in the
same way: given an observed \\r\_{XY \cdot Z}\\, the unbiased estimator
of \\\rho\_{XY \cdot Z}\\ is \\r \cdot {}\_2F_1(1/2, 1/2; (n - J - 2)/2;
1 - r^2)\\ (Olkin & Finn, 1995).

**Magnitude of the bias.** To leading order (Olkin & Finn, 1995,
equation 11): \$\$\mathrm{E}\[r\_{XY \cdot Z}\] - \rho\_{XY \cdot Z}
\\\approx\\ -\rho\_{XY \cdot Z} (1 - \rho\_{XY \cdot Z}^2) / \[2 (n -
J - 1)\].\$\$ For \\\rho\_{XY \cdot Z} = 0.4\\, \\n = 30\\, \\J = 2\\,
the bias is about \\-0.006\\; for \\n = 30\\, \\J = 10\\, about
\\-0.009\\.

**Tuning the series convergence.** The underlying \\{}\_2F_1\\ series is
summed by forward recurrence and stops when the relative contribution
falls below `getOption("DMAR.expected_r.tol", 1e-15)`, with a maximum of
`getOption("DMAR.expected_r.max_iter", 5000L)` terms. These are the same
options as for
[`expected_r`](https://yelleknek.github.io/DMAR/reference/expected_r.md);
users almost never need to change them.

## References

Anderson, T. W. (2003). *An introduction to multivariate statistical
analysis* (4th ed.), Sections 4.2 and 4.3. Wiley.

Hotelling, H. (1953). New light on the correlation coefficient and its
transforms. *Journal of the Royal Statistical Society, Series B, 15*(2),
193–232.

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin, 118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

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
[`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. A single value: rho = 0.4, n = 30, J = 3 controls.
expected_partial_r(rho = 0.4, n = 30, J = 3)
#>  rho n  J expected_partial_r bias    relative_bias
#>  0.4 30 3 0.394              0.00648 0.0162       

# 2. Bias grows with J at fixed (rho, n):
expected_partial_r(rho = 0.4, n = 30, J = c(1, 2, 5, 10, 20))
#>  rho n  J  expected_partial_r bias    relative_bias
#>  0.4 30 1  0.394              0.00602 0.015        
#>  0.4 30 2  0.394              0.00624 0.0156       
#>  0.4 30 5  0.393              0.00702 0.0176       
#>  0.4 30 10 0.391              0.00888 0.0222       
#>  0.4 30 20 0.381              0.0187  0.0469       

# 3. Reduces to expected_r() when J = 1 ... 0:
#        ("partialing zero variables" is the simple correlation case
#        with one fewer df by convention; see Anderson (2003).)
expected_partial_r(rho = 0.5, n = 11, J = 1)$expected_partial_r
#> [1] 0.4786588
expected_r(rho = 0.5, n = 10)$expected_r
#> [1] 0.4786588
```
