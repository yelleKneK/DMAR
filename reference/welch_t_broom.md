# Broom-Style Tidy / Glance Methods for `welch_t()`

`tidy()` returns the single mean-difference estimate and its confidence
interval in the broom convention; `glance()` coincides with it, since a
two-sample *t* test reports one estimand and there are no extra
model-level statistics to add.

## Usage

``` r
# S3 method for class 'dmar_welch_t'
tidy(x, ...)

# S3 method for class 'dmar_welch_t'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_welch_t` object returned by
  [`welch_t`](https://yelleknek.github.io/DMAR/reference/welch_t.md).

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `term`, `estimate`, `ci_lower`,
`ci_upper`, `statistic`, `df`, `p_value`, and `conf_level`.

## Author

Ken Kelley <kkelley@nd.edu>
