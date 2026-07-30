# A Reliability Coefficient Estimate

Returns a one-row `data.frame` in the column convention used by the
broom ecosystem (`term`, `estimate`, `se`, `ci_lower`, `ci_upper`). The
`term` is the coefficient name (`"alpha"`, `"omega"`, etc.).

## Usage

``` r
# S3 method for class 'dmar_reliability'
tidy(x, ...)
```

## Arguments

- x:

  A `dmar_reliability` object returned by any of
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  [`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
  or
  [`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md).

- ...:

  Unused.

## Value

A one-row `data.frame`.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
S <- cov(mtcars[, 1:5])
res <- reliability_alpha(S = S, N = 32, ci_method = "feldt")
generics::tidy(res)
#>    term  estimate se  ci_lower  ci_upper
#> 1 alpha 0.4671591 NA 0.1066876 0.7112825
generics::glance(res)
#>   coefficient  estimate se  ci_lower  ci_upper conf_level nobs n_items
#> 1       alpha 0.4671591 NA 0.1066876 0.7112825       0.95   32       5
#>   ci_method
#> 1     feldt
```
