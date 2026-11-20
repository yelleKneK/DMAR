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

**Three formulas.** Three measures of sample skewness are in use. With
\\m_r = \sum\_{i=1}^{n} (x_i - \bar{x})^r / n\\, they are \\g_1 = m_3 /
m_2^{3/2}\\ (Type 1), \\G_1 = g_1 \sqrt{n(n-1)} / (n-2)\\ (Type 2, the
value above), and \\b_1 = g_1 ((n-1)/n)^{3/2}\\ (Type 3). The three
differ by factors that depend only on \\n\\ and approach 1 as \\n\\
grows. `skewness()` reports \\G_1\\, which is also the skewness that
[`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md)
reports.
[`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)
describes the three corresponding measures of kurtosis.

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
