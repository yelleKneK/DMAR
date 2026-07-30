# broom-style tidy / glance methods for the CI family

These methods dispatch on two shared S3 classes that the `ci_*` family
of functions tags onto its returned `data.frame`s, so the broom
ecosystem can consume DMAR confidence-interval output uniformly.

## Usage

``` r
# S3 method for class 'dmar_ci_long'
tidy(x, ...)

# S3 method for class 'dmar_ci_long'
glance(x, ...)

# S3 method for class 'dmar_ci_anova'
tidy(x, ...)

# S3 method for class 'dmar_ci_anova'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_ci_long` or `dmar_ci_anova` object returned by the
  corresponding `ci_*` function.

- ...:

  Unused.

## Value

A one-row `data.frame` in broom convention.

## Details

The class is determined by the output shape:

- `dmar_ci_long`:

  Long-format CI tables (rows for `lower_limit`, `upper_limit`, and
  optionally an estimate row whose `term` is the parameter name). Used
  by [`ci_cc`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
  [`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
  [`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
  and
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md).

- `dmar_ci_anova`:

  Wide-format ANOVA effect-size CI tables (one row, columns for the
  effect name, point estimate, bounds, and metadata). Used by
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
  [`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
  [`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
  and
  [`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).

Both methods produce a one-row `data.frame` in the broom convention with
`term`, `estimate`, `conf.low`, `conf.high`, and (when available)
`conf.level`. `glance()` on these objects calls `tidy()`; the output is
already model-level.

## Author

Ken Kelley <kkelley@nd.edu>
