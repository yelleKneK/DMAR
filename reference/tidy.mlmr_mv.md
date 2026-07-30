# A Multivariate FIML Regression Fit

Broom-style `tidy()` and `glance()` for
[`mlmr_mv`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md) fits.
`tidy()` returns one row per coefficient per outcome with the broom
dotted columns plus a leading `response` column identifying the outcome;
`glance()` returns a one-row model-level summary with the per-outcome
\\R^2\\ averaged and the number of responses reported.

## Usage

``` r
# S3 method for class 'mlmr_mv'
tidy(x, conf.int = FALSE, conf_level = NULL, standardized = FALSE, ...)

# S3 method for class 'mlmr_mv'
glance(x, ...)
```

## Arguments

- x:

  An `mlmr_mv` fit.

- conf.int:

  Logical: include `ci_lower` / `ci_upper`?

- conf_level:

  Ignored (the interval level is fixed at fit time and stored on the
  object); present for broom signature compatibility.

- standardized:

  Logical: include `std_estimate`?

- ...:

  Unused.

## Value

For `tidy.mlmr_mv`, a `data.frame` with columns `response`, `term`,
`estimate`, `se`, `statistic`, `p_value`, and optionally `ci_lower`,
`ci_upper`, `std_estimate`. For `glance.mlmr_mv`, a one-row `data.frame`
with `R2` (mean across outcomes), `df`, `logLik`, `AIC`, `BIC`,
`deviance`, `nobs`, and `n.responses`.

## See also

[`mlmr_mv`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md);
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) for the
univariate methods these mirror.

## Author

Ken Kelley <kkelley@nd.edu>
