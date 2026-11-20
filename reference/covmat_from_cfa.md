# Model Implied Covariance Matrix of a One-Factor Confirmatory Factor Model

Computes the covariance matrix of the indicators that a single-factor
confirmatory factor model implies, from the factor loadings and the
residual variances of the indicators: \$\$\Sigma = \lambda
\lambda^\top + \Theta, \qquad \Theta = \mathrm{diag}(\theta_1, \ldots,
\theta_p).\$\$ The loadings and residual variances may be population
values, for a covariance matrix to plan a study or to simulate from, or
estimates from a fitted model, for the covariance matrix the fitted
model implies. The calculation is the same in both cases, and the result
is the model implied covariance matrix: the population covariance matrix
of the model when the inputs are population values, and the fitted
model's covariance matrix when they are estimates.

## Usage

``` r
covmat_from_cfa(factor_loading, residual_var, tol_det = 1e-05)
```

## Arguments

- factor_loading:

  A numeric vector of factor loadings, one per indicator (at least two),
  on the metric of a factor with variance 1. A single-row or
  single-column matrix is accepted and read as a vector. Population
  loadings or estimates from a fitted model.

- residual_var:

  A numeric vector of residual variances (unique or error variances),
  one per indicator, or a single value used for every indicator (the
  equal residual variance model). Population values or estimates from a
  fitted model.

- tol_det:

  The value below which the determinant of the implied covariance matrix
  triggers a warning that the matrix may not be positive definite
  (default `1e-05`).

## Value

A list with the single element `model_implied_cov`: the model implied
covariance matrix of the indicators, a symmetric \\p \times p\\ matrix.

## Details

Under the single-factor model, each indicator score is \\x_i = \lambda_i
\xi + \delta_i\\, where \\\xi\\ is the common factor and \\\delta_i\\ is
the residual of indicator *i*, with variance \\\theta_i\\, uncorrelated
with the factor and with the other residuals. The general form of the
model implied covariance matrix, \\\Lambda \Phi \Lambda^\top + \Theta\\,
reduces to \\\lambda \lambda^\top + \Theta\\ here: \\\Lambda\\ is the
\\p \times 1\\ column of loadings \\\lambda = (\lambda_1, \ldots,
\lambda_p)^\top\\ (`factor_loading`), the factor variance \\\Phi\\ is
fixed at 1, and \\\Theta\\ is the diagonal matrix of the residual
variances (`residual_var`). Because \\\Theta\\ is diagonal, the
residuals are uncorrelated by construction; a model with correlated
residuals needs a more general formulation than this function provides.

Estimates from a fit that fixed one loading at 1 and estimated the
factor variance \\\phi\\ are put on this metric by multiplying each
loading by \\\sqrt{\phi}\\; the residual variances are unchanged, and so
is the implied matrix.
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) fixes the
factor variance at 1, so its loadings and residual variances can be
given as they are.

Estimated residual variances at or below 0 (a Heywood case) are used as
given; the function does not check their sign. The determinant check
warns only when the determinant of the implied matrix is below
`tol_det`, which a Heywood case need not produce: with loadings .9, .8,
and 1.1 and residual variances .2, .3, and -.01, the implied matrix is
positive definite and no warning is given. A determinant above `tol_det`
does not show that the matrix is positive definite either, because an
even number of negative eigenvalues gives a positive determinant. Check
the residual variances themselves for a Heywood case, and the
eigenvalues of the result for positive definiteness.

The `cfa` spelling matches the
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) naming.

## See also

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five indicators with equal loadings and equal residual variances
covmat_from_cfa(factor_loading = rep(0.7, 5), residual_var = rep(0.51, 5))
#> $model_implied_cov
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 1.00 0.49 0.49 0.49 0.49
#> [2,] 0.49 1.00 0.49 0.49 0.49
#> [3,] 0.49 0.49 1.00 0.49 0.49
#> [4,] 0.49 0.49 0.49 1.00 0.49
#> [5,] 0.49 0.49 0.49 0.49 1.00
#> 

# Unequal loadings
covmat_from_cfa(factor_loading = c(0.5, 0.6, 0.7, 0.8),
                residual_var   = c(0.75, 0.64, 0.51, 0.36))
#> $model_implied_cov
#>      [,1] [,2] [,3] [,4]
#> [1,] 1.00 0.30 0.35 0.40
#> [2,] 0.30 1.00 0.42 0.48
#> [3,] 0.35 0.42 1.00 0.56
#> [4,] 0.40 0.48 0.56 1.00
#> 

# From the estimates of a fitted one-factor model, whose factor variance
# cfa_1() fixes at 1: the covariance matrix the fitted model implies for
# four verbal tests of holzinger_swineford. The tests have no missing
# scores, so the model is fit to their covariance matrix, which gives
# the estimates the raw scores give.
items <- c("t5_general_information", "t6_paragraph_comprehension",
           "t7_sentence", "t9_word_meaning")
fit <- cfa_1(S = cov(holzinger_swineford[, items]), N = 301)
covmat_from_cfa(
  factor_loading = fit$estimate[grepl("=~", fit$syntax)],
  residual_var = fit$estimate[fit$syntax %in% paste(items, "~~", items)])
#> $model_implied_cov
#>           [,1]     [,2]     [,3]     [,4]
#> [1,] 152.77322 29.55780 45.69592 67.70201
#> [2,]  29.55780 12.15598 12.72080 18.84684
#> [3,]  45.69592 12.72080 26.55657 29.13693
#> [4,]  67.70201 18.84684 29.13693 58.62156
#> 
```
