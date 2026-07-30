# Glance at a Multiple-Factor CFA Fit

Returns a one-row `data.frame` of model-level summaries from a
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) table, in
the column convention used by the broom ecosystem.

## Usage

``` r
# S3 method for class 'dmar_cfa_k'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_cfa_k` object returned by
  [`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md).

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `chi_square`, `df`, `p_value`,
`cfi`, `tli`, `rmsea`, `rmsea_low`, `rmsea_high`, `srmr`, `AIC`, `BIC`,
`logLik`.

## Author

Ken Kelley <kkelley@nd.edu>
