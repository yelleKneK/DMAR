# Tidy / Glance Methods for ci_R2 Output

Returns the standard broom-style one-row summary of the \\R^2\\
confidence interval: `term` (always `"R2"`), `estimate`, `ci_lower`,
`ci_upper`, `conf_level`.

## Usage

``` r
# S3 method for class 'dmar_ci_R2'
tidy(x, ...)

# S3 method for class 'dmar_ci_R2'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_ci_R2` object returned by `ci_R2`.

- ...:

  Unused.

## Value

A one-row `data.frame` in broom convention.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
res <- ci_R2(R2 = 0.25, N = 100, p = 5)
generics::tidy(res)
#>   term estimate   ci_lower  ci_upper conf_level
#> 1   R2     0.25 0.08007896 0.3720553       0.95
generics::glance(res)
#>   term estimate   ci_lower  ci_upper conf_level
#> 1   R2     0.25 0.08007896 0.3720553       0.95
```
