# Tidy an MBCO Mediation Table

Returns the effect rows of a
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
table in the column convention used by the broom ecosystem. The
`statistic` column is the MBCO likelihood ratio statistic.

## Usage

``` r
# S3 method for class 'dmar_mediation_mbco'
tidy(x, ...)
```

## Arguments

- x:

  A `dmar_mediation_mbco` object returned by
  [`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md).

- ...:

  Unused.

## Value

A `data.frame` with columns `term`, `estimate`, `se`, `statistic`,
`p_value`, `ci_lower`, `ci_upper`.

## Author

Ken Kelley <kkelley@nd.edu>
