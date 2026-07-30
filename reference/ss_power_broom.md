# broom-style tidy / glance methods for the power-based sample size planners

Functions in the `ss_power_*` family return a tidy long- format
`data.frame` with `term` and `value` columns: a row for the recommended
sample size, a row for the realized statistical power, and rows that
echo the planning inputs. When those returns are tagged with the class
`dmar_ss_power` (which the `ss_power_*` functions do),
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
produce a broom-style one-row summary.

## Usage

``` r
# S3 method for class 'dmar_ss_power'
tidy(x, ...)

# S3 method for class 'dmar_ss_power'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_ss_power` object returned by one of
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
  [`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
  [`ss_power_r`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
  [`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
  [`ss_power_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
  and the rest of the `ss_power_*` family.

- ...:

  Unused.

## Value

A one-row `data.frame` in broom convention.

## Details

`tidy()` returns the recommended sample size and the achieved power as a
single row with broom-conventional column names. `glance()` returns the
same row extended with the input planning parameters that were echoed in
the long output.

## Author

Ken Kelley <kkelley@nd.edu>
