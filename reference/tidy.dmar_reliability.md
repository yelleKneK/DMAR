# A Reliability Coefficient Estimate

Returns a one-row `data.frame` in the column convention used by the
broom ecosystem (`term`, `estimate`, `se`, `ci_lower`, `ci_upper`,
`conf_level`, `wald_df`). The `term` is the coefficient name (`"alpha"`,
`"omega"`, etc.), `conf_level` is the level of the interval, so that
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
on the row states it, and `wald_df` the degrees of freedom of the *t*
reference of a Wald interval (the `wald_df` argument of
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
and
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)),
`Inf` for the normal reference or an interval that is not a Wald
interval, so that rows on the two references can be told apart after
they are stacked.

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
# Coefficient alpha for the three verbal tests of the Holzinger and
# Swineford battery, from their covariance matrix.
S <- cov(holzinger_swineford[, c("t6_paragraph_comprehension",
                                 "t7_sentence", "t9_word_meaning")])
res <- reliability_alpha(S = S, N = 301, ci_method = "feldt")
generics::tidy(res)
#>    term  estimate se  ci_lower ci_upper conf_level wald_df h1_diagnostic
#> 1 alpha 0.8305929 NA 0.7945197 0.861239       0.95     Inf          <NA>
#>   wald_calibration_note
#> 1                  <NA>
generics::glance(res)
#>   coefficient  estimate se  ci_lower ci_upper conf_level nobs n_items ci_method
#> 1       alpha 0.8305929 NA 0.7945197 0.861239       0.95  301       3     feldt
#>   wald_df h1_diagnostic wald_calibration_note
#> 1     Inf          <NA>                  <NA>
```
