# Tidy a one-factor CFA fit

Returns a one-row-per-parameter `data.frame` in the column convention
used by the broom ecosystem (`term`, `estimate`, `std.error`,
`statistic`, `p.value`). Use
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) with the
default `output = "verbose"` (which attaches a `dmar_cfa_1` class to the
returned data.frame) and call
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) on
it.

## Usage

``` r
# S3 method for class 'dmar_cfa_1'
tidy(x, ...)
```

## Arguments

- x:

  A `dmar_cfa_1` object returned by
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md).

- ...:

  Unused.

## Value

A `data.frame` with columns `term`, `estimate`, `std.error`,
`statistic`, `p.value`.

## Details

Fit-index rows (CFI, TLI, RMSEA, AIC, BIC, composite reliability, etc.)
are excluded from the tidy output; they belong in
[`glance.dmar_cfa_1`](https://yelleknek.github.io/DMAR/reference/glance.dmar_cfa_1.md).

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
cov_mat <- matrix(
  c(1.384, 1.484, 1.988, 2.429, 3.031,
    1.484, 2.756, 2.874, 3.588, 4.390,
    1.988, 2.874, 4.845, 4.894, 6.080,
    2.429, 3.588, 4.894, 6.951, 7.476,
    3.031, 4.390, 6.080, 7.476, 10.313),
  nrow = 5)
fit <- cfa_1(N = 300, S = cov_mat)
generics::tidy(fit)
#>           term  estimate   std.error  statistic      p.value
#> 1     lambda_1 0.9965424 0.054801415  18.184610 0.000000e+00
#> 2     lambda_2 1.4510476 0.075847938  19.131009 0.000000e+00
#> 3     lambda_3 1.9897431 0.098274764  20.246734 0.000000e+00
#> 4     lambda_4 2.4524282 0.115205053  21.287505 0.000000e+00
#> 5     lambda_5 3.0347148 0.138559250  21.901928 0.000000e+00
#> 6          phi 1.0000000 0.000000000         NA           NA
#> 7        psi_1 0.3862898 0.035138804  10.993255 0.000000e+00
#> 8        psi_2 0.6412741 0.060291955  10.636147 0.000000e+00
#> 9        psi_3 0.8697722 0.087067507   9.989630 0.000000e+00
#> 10       psi_4 0.9134254 0.102258184   8.932541 0.000000e+00
#> 11       psi_5 1.0691287 0.134709305   7.936562 1.998401e-15
#> 12 loading_sum 9.9244760 0.421216261  23.561474 0.000000e+00
#> 13   error_sum 3.8798902 0.175431675  22.116247 0.000000e+00
#> 14       omega 0.9621012 0.003536945 272.014789 0.000000e+00
generics::glance(fit)
#>   cfi   tli rmsea rmsea_low rmsea_high      AIC      BIC    logLik comp_rel
#> 1   1 1.003     0         0      0.058 4793.531 4830.568 -2386.765    0.962
```
