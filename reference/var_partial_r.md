# Asymptotic Variance of the Partial Correlation Coefficient

Computes the large-sample variance of the sample partial correlation
coefficient \\r\_{XY \cdot Z_1 \cdots Z\_{p-1}}\\ under multivariate
normality. Two formulas are available: the first-order variance on the
raw \\r\\ scale (the default), and the variance \\1/(N - p - 2)\\ of the
Fisher's *Z* transformation (Fisher, 1921, 1924), which is the
appropriate quantity for constructing a confidence interval by
transformation and back-transformation.

## Usage

``` r
var_partial_r(r, N, p = 2, fisher_z = FALSE)
```

## Arguments

- r:

  The sample partial correlation coefficient \\r\_{XY \cdot Z_1 \cdots
  Z\_{p-1}}\\. Must be in \\\[-1, 1\]\\.

- N:

  Total sample size.

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
  variables are partialed out. Defaults to 2 (one variable partialed
  out).

- fisher_z:

  Logical. If `FALSE` (the default), the function returns the raw-scale
  first-order variance. If `TRUE`, it returns the variance of the
  Fisher's *Z* transformation \\Z = \mathrm{arctanh}(r)\\ under the
  Fisher (1924) reduction.

## Value

A one-row `data.frame` with columns `term` (either `"var_partial_r"` or
`"var_fisher_z_partial_r"`) and `value` (the requested variance).

## Details

**Raw scale, first-order variance.** The first-order variance of a
simple correlation is \\(1 - \rho^2)^2 / N\\ (Olkin & Finn, 1995,
Equation 3). The sample partial correlation given \\p - 1\\ variables is
distributed as a simple correlation from \\N - (p - 1)\\ observations,
and the functions here use \$\$\mathrm{Var}(\hat r\_{XY \cdot Z})
\\\approx\\ \frac{(1 - \rho^2\_{XY \cdot Z})^2}{N - p},\$\$ which
differs from that at order \\1/N^2\\. The function evaluates this with
\\\hat r\\ substituted for \\\rho\\. Olkin and Finn (1995) give the
variance of a simple correlation; the divisor \\N - p\\ is this
package's choice and is not taken from that article.

**Fisher's *Z* Transformation.** Fisher (1921) showed that for a Pearson
correlation, the transformation \\Z = \tfrac{1}{2}\log\\(1+r)/(1-r)\\ =
\mathrm{arctanh}(r)\\ is approximately normal with variance \\1/(N -
3)\\. Fisher (1924, p. 330) showed that the partial correlation based on
\\N\\ observations with \\p - 1\\ variables partialed out is distributed
as a simple correlation from a sample reduced in size by \\p - 1\\;
combined with the Fisher (1921) variance of \\Z\\, the Fisher's *Z*
transformation of \\\hat r\_{XY \cdot Z_1 \cdots Z\_{p-1}}\\ is
therefore approximately normal with variance \\1/(N - p - 2)\\. This is
the standard ingredient for constructing a confidence interval on
\\\rho\_{XY \cdot Z}\\ by transforming, building a Wald interval on
\\Z\\, and back-transforming with \\\tanh\\.

## References

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum.

Fisher, R. A. (1915). Frequency distribution of the values of the
correlation coefficient in samples from an indefinitely large
population. *Biometrika, 10*(4), 507–521.

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Fisher, R. A. (1924). The distribution of the partial correlation
coefficient. *Metron, 3*, 329–332.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on contrasts, Chapter 5 on multiple
comparisons, and Chapter 9 on ANCOVA.)

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin, 118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

## See also

[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
[`var_R2`](https://yelleknek.github.io/DMAR/reference/var_R2.md),
[`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# First-order variance for r_p = .35, N = 80, and two variables
#     partialed out (p = 3 predictors: X and the two variables).
var_partial_r(r = 0.35, N = 80, p = 3)
#>  term          value
#>  var_partial_r 0.010

# Variance of Fisher's Z transformation (Fisher, 1921, 1924) for the
#     same setting, useful for building a CI on rho_p.
var_partial_r(r = 0.35, N = 80, p = 3, fisher_z = TRUE)
#>  term                   value
#>  var_fisher_z_partial_r 0.013
```
