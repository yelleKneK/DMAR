# Glance at a one-factor CFA fit

Returns a one-row `data.frame` of model-level summaries in the column
convention used by the broom ecosystem (`cfi`, `tli`, `rmsea`,
`rmsea_low`, `rmsea_high`, `AIC`, `BIC`, `logLik`, `comp_rel`, `npar`).

## Usage

``` r
# S3 method for class 'dmar_cfa_1'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_cfa_1` object returned by
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md).

- ...:

  Unused.

## Value

A one-row `data.frame`.

## Author

Ken Kelley <kkelley@nd.edu>
