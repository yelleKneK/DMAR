# Bias-Corrected Sample Skewness

Computes the sample skewness of a numeric vector using the
bias-corrected (SAS/SPSS Type 2) formula. Skewness measures asymmetry of
the distribution: zero is symmetric, positive values indicate a
right-tail heavier than the left, negative values the reverse.

## Usage

``` r
skewness(x, na_rm = TRUE)
```

## Arguments

- x:

  A numeric vector.

- na_rm:

  Logical. If `TRUE` (the default), missing values are removed before
  computation. If `FALSE`, the result is `NA` when `x` contains any
  `NA`.

## Value

A single numeric value: the bias-corrected sample skewness, or
`NA_real_` when fewer than three non-missing observations are available
or when the sample standard deviation is zero.

## Details

The reported value is \$\$\hat\gamma_1^{(2)} =
\frac{n}{(n-1)(n-2)}\sum\_{i=1}^{n}\left(\frac{x_i -
\bar{x}}{s}\right)^3,\$\$ where \\s\\ is the (divisor-\\n-1\\) sample
standard deviation. This is sometimes called the “Type 2” or
SAS/SPSS-default form; it is approximately unbiased under normality.

**Why isn't this in base R?** R Core has historically deferred
higher-order moment statistics to contributed packages, in part because
three popular formulas exist (biased Type 1, bias-corrected Type 2, and
Minitab Type 3) and choosing a default would be opinionated. DMAR adopts
Type 2, which is the form most often used in psychometric reporting and
the one already used internally by
[`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md).

**Diagnostic interpretation.** As a rough rule of thumb,
\\\|\mathrm{skewness}\| \> 2\\ is sometimes flagged as indicative of
departures from normality large enough to threaten normal-theory
inference (e.g., maximum likelihood estimation in factor analysis or
structural equation modeling).

## References

Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample
skewness and kurtosis. *The Statistician, 47*(1), 183–189.
[doi:10.1111/1467-9884.00122](https://doi.org/10.1111/1467-9884.00122)

## See also

[`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md),
[`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md)

Other descriptive statistics:
[`descriptives()`](https://yelleknek.github.io/DMAR/reference/descriptives.md),
[`kurtosis()`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Symmetric data: skewness near zero.
set.seed(113)
skewness(rnorm(1000))
#> [1] -0.03196212

# Right-skewed data: positive value.
skewness(rexp(1000, rate = 1))
#> [1] 1.735017

# The classic 1:5 example: exactly symmetric (returns 0).
skewness(1:5)
#> [1] 0
```
