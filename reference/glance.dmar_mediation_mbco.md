# Glance at an MBCO Mediation Fit

Returns a one-row `data.frame` of full-model summaries from a
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
table, in the column convention used by the broom ecosystem.

## Usage

``` r
# S3 method for class 'dmar_mediation_mbco'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_mediation_mbco` object returned by
  [`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md).

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `nobs`, `npar`, `deviance`, `AIC`,
`BIC`, `logLik`.

## Author

Ken Kelley <kkelley@nd.edu>
