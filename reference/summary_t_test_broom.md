# Broom-Style Tidy / Glance Methods for `summary_t_test()`

`tidy()` returns the single mean-difference estimate and its confidence
interval in the broom convention; `glance()` coincides with it, since a
two-sample *t* test reports one estimand and there are no extra
model-level statistics to add.

## Usage

``` r
# S3 method for class 'dmar_summary_t_test'
tidy(x, ...)

# S3 method for class 'dmar_summary_t_test'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_summary_t_test` object returned by
  [`summary_t_test`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md).

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `term`, `estimate`, `ci_lower`,
`ci_upper`, `statistic`, `df`, `p_value`, and `conf_level`.

## Author

Ken Kelley <kkelley@nd.edu>
