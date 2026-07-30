# Generate a Population Covariance Matrix From a One-Factor Confirmatory Factor Model

Given a vector of factor loadings (\\\lambda\\) and the corresponding
vector of unique (error) variances (\\\psi^2\\), this function computes
the implied population covariance matrix under a single-factor
confirmatory factor model: \$\$\Sigma = \Lambda \Lambda^\top + \Psi.\$\$
This builds the model implied covariance for a confirmatory factor model
with uncorrelated errors. The unique-variances matrix \\\Psi\\ is
strictly diagonal, so no residual covariances (correlated uniquenesses)
are allowed; a model with correlated errors is outside the scope of this
function. Because the model is a single common factor, the loadings
matrix \\\Lambda\\ reduces to a column vector \\\lambda\\ of length
\\p\\ (one per indicator), and the unique-variances matrix \\\Psi\\ is
the diagonal \\\mathrm{diag}(\psi^2)\\ of a length-\\p\\ vector. The
formals are named in the lowercase vector form (`lambda` and
`psi_squared`) to reflect this; the matrix-form symbols (\\\Lambda\\,
\\\Psi\\) remain in the mathematical exposition above.

## Usage

``` r
covmat_from_cfa(lambda, psi_squared, ...)

covmat_from_cfm(lambda, psi_squared, ...)
```

## Arguments

- lambda:

  A numeric vector of factor loadings (one per indicator). Can also be
  supplied as a single-row or single-column matrix and will be coerced
  to a vector.

- psi_squared:

  A numeric vector of unique (error) variances, one per indicator.
  Recycled to match the length of `lambda` when a single value is
  supplied (equal-error-variance model).

- ...:

  Optional advanced controls. Currently the only recognized passthrough
  is `tol_det`: the tolerance below which the determinant of the implied
  covariance matrix triggers a positive- definite warning (default
  `1e-05`). The argument is hidden in `...` because the default is
  appropriate for almost every application; users who pass an
  unrecognized name through `...` (for example, a misspelling) are
  notified with a warning so the typo is not silently ignored.

## Value

A list with the single element `population_cov`: the implied population
covariance matrix of the manifest indicators (\\p \times p\\,
symmetric).

## Details

Under the single-factor common-factor model each indicator score is
\\x_i = \lambda_i \xi + \delta_i\\, where \\\xi\\ is the (standardized)
latent factor and \\\delta_i\\ is the indicator-specific residual with
variance \\\psi_i^2\\. The population covariance among the manifest
indicators is therefore \$\$\Sigma = \Lambda \Lambda^\top + \Psi,\$\$
where \\\Lambda\\ is the \\p \times 1\\ column of loadings \\(\lambda_1,
\ldots, \lambda_p)^\top\\ and \\\Psi = \mathrm{diag}(\psi_1^2, \ldots,
\psi_p^2)\\. In code we work with the vectors `lambda` and `psi_squared`
directly. Because \\\Psi\\ is built with
[`diag()`](https://rdrr.io/r/base/diag.html) from a length-\\p\\ vector,
the errors are uncorrelated by construction: there is no way to specify
a residual covariance between two indicators. A confirmatory factor
model with correlated errors (correlated uniquenesses) requires a more
general formulation than this function provides.

`covmat_from_cfm()` is a backward-compatible alias for
`covmat_from_cfa()`; the two are identical, and the `cfa` spelling
matches the
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) naming.

## See also

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five indicators with equal loadings and equal error variances
covmat_from_cfa(lambda = rep(0.7, 5), psi_squared = rep(0.51, 5))
#> $population_cov
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 1.00 0.49 0.49 0.49 0.49
#> [2,] 0.49 1.00 0.49 0.49 0.49
#> [3,] 0.49 0.49 1.00 0.49 0.49
#> [4,] 0.49 0.49 0.49 1.00 0.49
#> [5,] 0.49 0.49 0.49 0.49 1.00
#> 

# Unequal loadings
covmat_from_cfa(lambda      = c(0.5, 0.6, 0.7, 0.8),
                psi_squared = c(0.75, 0.64, 0.51, 0.36))
#> $population_cov
#>      [,1] [,2] [,3] [,4]
#> [1,] 1.00 0.30 0.35 0.40
#> [2,] 0.30 1.00 0.42 0.48
#> [3,] 0.35 0.42 1.00 0.56
#> [4,] 0.40 0.48 0.56 1.00
#> 
```
