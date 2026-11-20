# Glance at a Multiple-Factor CFA Fit

Returns a one-row `data.frame` of model-level summaries from a
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
table, in the column convention used by the broom ecosystem.

## Usage

``` r
# S3 method for class 'dmar_fa_confirmatory'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_fa_confirmatory` object returned by
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).

- ...:

  Unused.

## Value

A one-row `data.frame` with columns `chi_square`, `df`, `p_value`,
`cfi`, `tli`, `rmsea`, `rmsea_low`, `rmsea_high`, `srmr`, `AIC`, `BIC`,
`logLik`, and the columns that say what produced them, so that rows of
different fits can be stacked and still read: `estimator`; `missing`,
the missing data method of the fit (`"listwise"`, `"pairwise"`,
`"fiml"`, `"two_stage"`, `"robust_two_stage"`, or `"none"` for a
covariance matrix); `test`, the statistic in `chi_square` (`"standard"`,
`"scaled"`, or a two-stage statistic); `fit_indices`, the version of
`cfi`, `tli`, and `rmsea` (`"robust"`, `"scaled"`, or `"standard"`); and
`rmsea_ci_level`, the level of the RMSEA interval (`rmsea_conf_level`,
0.90 by default, where the `conf_level` column of `tidy()` is the level
of the parameter intervals). After two-stage estimation `logLik` is the
complete data normal log likelihood at the EM estimates and `AIC` and
`BIC` are the two-stage criteria of Lai (2021), comparable only with
two-stage fits of the same variables, not with a FIML or listwise row;
and the difference of two rows' `chi_square` or `logLik` is not a chi
square test of nested models under a scaled or two-stage statistic
([`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html)
on `output = "fit"` gives the right one).

## References

Lai, K. (2021). Using information criteria under missing data: Full
information maximum likelihood versus two-stage estimation. *Structural
Equation Modeling, 28*(2), 278–291.
[doi:10.1080/10705511.2020.1780925](https://doi.org/10.1080/10705511.2020.1780925)

## Author

Ken Kelley <kkelley@nd.edu>
