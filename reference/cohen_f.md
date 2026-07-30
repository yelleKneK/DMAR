# Cohen's *f* Effect Size

Computes Cohen's *f* = \\\sigma_m / \sigma\\, the population standard
deviation of means relative to the within-group standard deviation, by
any of three equivalent specifications:

1.  raw population means and within-group variance,

2.  the population proportion of variance accounted for, \\\eta^2\\,

3.  \\\sigma_m\\ and \\\sigma\\ directly.

*Cohen's f is a population quantity*; supplied with population
parameters it returns the population value, supplied with sample
estimates it returns the corresponding sample value.

## Usage

``` r
cohen_f(
  mu = NULL,
  sigma_squared = NULL,
  n = NULL,
  eta_squared = NULL,
  sigma_m = NULL,
  sigma = NULL
)
```

## Arguments

- mu:

  Numeric vector of population means (one per group). Use together with
  `sigma_squared`.

- sigma_squared:

  The within-group variance. Use together with `mu`.

- n:

  Optional. Per-group sample sizes (a single number for equal group
  sizes, or a vector of length `length(mu)` for unequal). When `NULL`,
  equal weighting across groups is used (i.e., the population variance
  of `mu` is computed with the \\1/k\\ divisor).

- eta_squared:

  The population proportion of variance accounted for. Use this argument
  alone.

- sigma_m:

  The population standard deviation of the means (\\\sigma_m\\). Use
  together with `sigma`.

- sigma:

  The within-group standard deviation (\\\sigma\\). Use together with
  `sigma_m`.

## Value

A 1-row `data.frame` with columns `term` and `value`; `term` is
`"cohen_f"` and `value` is the computed value.

## Details

All three calling modes return the same value when applied to compatible
inputs (Cohen 1988, eq. 8.2.1):

- Raw form: \\f = \sqrt{\sum n_j (\mu_j - \bar\mu)^2 / N \cdot
  1/\sigma^2}\\.

- From \\\eta^2\\: \\f = \sqrt{\eta^2 / (1 - \eta^2)}\\.

- From the variance ratio: \\f = \sigma_m / \sigma\\.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

## See also

[`ci_srsnr`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`ci_snr`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# (1) From raw means and within-group variance:
cohen_f(mu = c(94, 91, 92, 83), sigma_squared = 67.375)
#>  term    value
#>  cohen_f 0.51 

# Equal n weights are the default; equivalent with explicit equal n:
cohen_f(mu = c(94, 91, 92, 83), sigma_squared = 67.375, n = 6)
#>  term    value
#>  cohen_f 0.51 

# Unequal n:
cohen_f(mu = c(94, 91, 92, 83), sigma_squared = 67.375, n = c(4, 6, 5, 5))
#>  term    value
#>  cohen_f 0.498

# (2) From eta_squared:
cohen_f(eta_squared = 0.10)
#>  term    value
#>  cohen_f 0.333

# (3) From sigma_m and sigma directly:
cohen_f(sigma_m = 4, sigma = 8)
#>  term    value
#>  cohen_f 0.5  
```
