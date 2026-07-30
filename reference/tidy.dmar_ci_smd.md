# Tidy / Glance Methods for ci_smd Output

Returns the standard broom-style one-row summary of the SMD confidence
interval: `term` (always `"smd"`), `estimate`, `ci_lower`, `ci_upper`,
`conf_level`.

## Usage

``` r
# S3 method for class 'dmar_ci_smd'
tidy(x, ...)

# S3 method for class 'dmar_ci_smd'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_ci_smd` object returned by `ci_smd`.

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `term` (always `"smd"`), `estimate`
(the point estimate of *d*), `ci_lower` (the lower confidence limit on
\\\delta\\), `ci_upper` (the upper confidence limit on \\\delta\\), and
`conf_level` (the confidence level used in construction; `NA` when
asymmetric `alpha_lower` / `alpha_upper` were supplied instead). Column
names follow the broom convention so the result composes with the broom
ecosystem.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
res <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
generics::tidy(res)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95
generics::glance(res)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95
```
