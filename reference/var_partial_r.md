# Asymptotic Variance of the Partial Correlation Coefficient

Computes the large-sample variance of the sample partial correlation
coefficient \\r\_{XY \cdot Z_1 \cdots Z_J}\\ under multivariate
normality. Two formulas are available: the classical asymptotic variance
on the raw \\r\\ scale (the default), and the variance \\1/(n - J - 3)\\
of the Fisher's *Z* transformation (Fisher, 1921, 1924), which is the
appropriate quantity for constructing a confidence interval by
transformation and back-transformation.

## Usage

``` r
var_partial_r(r, n, J = 1, fisher_z = FALSE)
```

## Arguments

- r:

  The sample partial correlation coefficient \\r\_{XY \cdot Z_1 \cdots
  Z_J}\\. Must be in \\\[-1, 1\]\\.

- n:

  Total sample size.

- J:

  Number of variables partialled out (i.e., the count of \\Z_1, \ldots,
  Z_J\\); must be at least 1. Defaults to 1.

- fisher_z:

  Logical. If `FALSE` (the default), the function returns the raw-scale
  asymptotic variance. If `TRUE`, it returns the variance of the
  Fisher's *Z* transformation \\Z = \mathrm{arctanh}(r)\\ under the
  Fisher (1924) reduction.

## Value

A one-row `data.frame` with columns `term` (either `"var_partial_r"` or
`"var_fisher_z_partial_r"`) and `value` (the requested variance).

## Details

**Raw-scale asymptotic variance.** Under multivariate normality the
partial correlation \\\hat r\_{XY \cdot Z}\\ has the large-sample
variance (Fisher, 1924, for the reduction; the simple-correlation
building block is, e.g., Olkin & Finn, 1995, their Equation 3):
\$\$\mathrm{Var}(\hat r\_{XY \cdot Z}) \\\approx\\ \frac{(1 -
\rho^2\_{XY \cdot Z})^2}{n - J - 1},\$\$ a direct generalization of the
classical asymptotic variance \\(1 - \rho^2)^2 / (n - 1)\\ of the simple
Pearson correlation (Fisher, 1915) with the degrees of freedom reduced
by the number of partialled variables. The function evaluates this with
\\\hat r\\ substituted for \\\rho\\.

**Fisher's *Z* Transformation.** Fisher (1921) showed that for a Pearson
correlation, the transformation \\Z = \tfrac{1}{2}\log\\(1+r)/(1-r)\\ =
\mathrm{arctanh}(r)\\ is approximately normal with variance \\1/(n -
3)\\. Fisher (1924) showed that the partial correlation based on \\n\\
observations with \\J\\ variables partialled out is distributed as a
simple correlation from a sample reduced in size by \\J\\; combined with
the Fisher (1921) variance of \\Z\\, the Fisher's *Z* transformation of
\\\hat r\_{XY \cdot Z_1 \cdots Z_J}\\ is therefore approximately normal
with variance \\1/(n - J - 3)\\. This is the standard ingredient for
constructing a confidence interval on \\\rho\_{XY \cdot Z}\\ by
transforming, building a Wald interval on \\Z\\, and back-transforming
with \\\tanh\\.

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
# Olkin-Finn (1995) asymptotic variance for r_p = .35, n = 80,
#     J = 2 control variables.
var_partial_r(r = 0.35, n = 80, J = 2)
#>  term          value
#>  var_partial_r 0.01 

# Variance of Fisher's Z transformation (Fisher, 1921, 1924) for the
#     same setting, useful for building a CI on rho_p.
var_partial_r(r = 0.35, n = 80, J = 2, fisher_z = TRUE)
#>  term                   value 
#>  var_fisher_z_partial_r 0.0133
```
