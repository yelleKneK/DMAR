# Glance at a Reliability Coefficient Estimate

Returns a one-row `data.frame` of model-level summaries in the column
convention used by the broom ecosystem (`coefficient`, `estimate`, `se`,
`ci_lower`, `ci_upper`, `conf_level`, `nobs`, `n_items`, `ci_method`,
`wald_df`, the degrees of freedom of the *t* reference of a Wald
interval, `Inf` for the normal reference or another interval).

## Usage

``` r
# S3 method for class 'dmar_reliability'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_reliability` object.

- ...:

  Unused.

## Value

A one-row `data.frame`.

## Author

Ken Kelley <kkelley@nd.edu>
