# Confidence Interval on the Largest Eigenvalue of a Sample Covariance Matrix

Computes an asymptotic confidence interval on the largest population
eigenvalue \\\lambda_1\\ of a population covariance matrix, given a
sample covariance matrix from \\n\\ observations on \\p\\ variables
under multivariate normality. Useful in principal-components analysis
and dimension-reduction settings to gauge whether the largest eigenvalue
is well-separated from the second.

## Usage

``` r
ci_eigenvalue(cov_matrix, n = NULL, conf_level = 0.95, k = 1)
```

## Arguments

- cov_matrix:

  Sample covariance matrix (a symmetric, positive- semidefinite numeric
  matrix), or a data frame whose columns are the variables (in which
  case `cov_matrix` is computed internally).

- n:

  Sample size (number of rows of the original data). Required when
  `cov_matrix` is supplied as a matrix; ignored when `cov_matrix` is a
  data frame (`nrow(cov_matrix)` is used instead).

- conf_level:

  Confidence level. Default `0.95`.

- k:

  Which eigenvalue (1 = largest, 2 = next, ...) to bracket. Default `1`.

## Value

A 3-row `data.frame` with rows ordered `"lower_limit"`, `"eigenvalue"`
(the sample eigenvalue point estimate), and `"upper_limit"`, so the
point estimate sits between its confidence limits.

## Details

**Asymptotic distribution.** Under multivariate normality with
eigenvalues \\\lambda_1 \> \lambda_2 \ge \cdots \ge \lambda_p\\, the
sample eigenvalues \\\hat\lambda_j\\ are asymptotically independent and
approximately normal with mean \\\lambda_j\\ and variance \\2
\lambda_j^2 / (n - 1)\\ when the eigenvalues are simple (well-separated)
(Anderson, 2003, Theorem 13.3.1; Muirhead, 1982, Section 9.7). The
asymptotic CI is therefore \$\$\hat\lambda_j \cdot \exp\\\left(\pm
z\_{1 - \alpha/2} \sqrt{\frac{2}{n - 1}}\right),\$\$ on the
multiplicative scale (equivalently, a Wald CI on \\\log \lambda_j\\ with
variance \\2/(n - 1)\\). The log scale is the natural
variance-stabilizing transformation for an eigenvalue.

**Caveats.** The asymptotic CI assumes well-separated population
eigenvalues. When the largest two eigenvalues are close, the sample
eigenvalue exhibits a "repulsion" phenomenon and the CI is biased
(typically too narrow). Diagnostic: if \\\hat\lambda_1 / \hat\lambda_2\\
is close to 1, the asymptotic CI should not be relied upon; a bootstrap
is preferable.

## References

Anderson, T. W. (2003). *An introduction to multivariate statistical
analysis* (4th ed.). Wiley. (See Chapter 13.)

Muirhead, R. J. (1982). *Aspects of multivariate statistical theory*.
Wiley. (See Section 9.7.)

## See also

[`prcomp`](https://rdrr.io/r/stats/prcomp.html),
[`eigen`](https://rdrr.io/r/base/eigen.html)

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. From a data frame:
set.seed(113)
X <- data.frame(matrix(rnorm(200), nrow = 50))
ci_eigenvalue(X, k = 1)
#>  term        value
#>  lower_limit 0.925
#>  eigenvalue  1.37 
#>  upper_limit 2.04 
#> 
#> Confidence level: 95%

# 2. From an explicit covariance matrix:
S <- cov(X)
ci_eigenvalue(S, n = nrow(X), k = 1)
#>  term        value
#>  lower_limit 0.925
#>  eigenvalue  1.37 
#>  upper_limit 2.04 
#> 
#> Confidence level: 95%
```
