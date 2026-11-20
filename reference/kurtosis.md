# Bias-Corrected Sample Excess Kurtosis

Computes the sample excess kurtosis of a numeric vector using the
bias-corrected (SAS/SPSS Type 2) formula. Excess kurtosis measures
tailedness relative to the normal distribution: zero matches a normal,
positive values indicate heavier tails (“leptokurtic”), negative values
indicate lighter tails (“platykurtic”).

## Usage

``` r
kurtosis(x, na_rm = TRUE)
```

## Arguments

- x:

  A numeric vector.

- na_rm:

  Logical. If `TRUE` (the default), missing values are removed before
  computation. If `FALSE`, the result is `NA` when `x` contains any
  `NA`.

## Value

A single numeric value: the bias-corrected sample excess kurtosis, or
`NA_real_` when fewer than four non-missing observations are available
or when the sample standard deviation is zero.

## Details

The reported value is \$\$\hat\gamma_2^{(2)} =
\frac{n(n+1)}{(n-1)(n-2)(n-3)}\sum\_{i=1}^{n}\left(\frac{x_i -
\bar{x}}{s}\right)^4 - \frac{3(n-1)^2}{(n-2)(n-3)},\$\$ where \\s\\ is
the (divisor-\\n-1\\) sample standard deviation. Subtracting
\\3(n-1)^2/((n-2)(n-3))\\ makes the expected value of the statistic 0
for a sample from a normal distribution; that is, *excess* kurtosis is
reported (rather than “raw” kurtosis, whose value for a normal
distribution is 3).

**Three formulas.** Three measures of sample kurtosis are in use, as
there are three of sample skewness (see *Three formulas* in
[`skewness`](https://yelleknek.github.io/DMAR/reference/skewness.md)).
With \\m_r = \sum\_{i=1}^{n} (x_i - \bar{x})^r / n\\, they are \\g_2 =
m_4 / m_2^2 - 3\\ (Type 1), \\G_2 = ((n+1) g_2 + 6)(n-1) /
((n-2)(n-3))\\ (Type 2, the value above), and \\b_2 = (g_2 + 3)(1 -
1/n)^2 - 3\\ (Type 3). \\G_2\\ and \\b_2\\ are each \\g_2\\ multiplied
by a factor and shifted by a constant, both of which depend only on
\\n\\, and the three agree as \\n\\ grows. `kurtosis()` reports \\G_2\\,
which is also the kurtosis that
[`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md)
reports.

## References

Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample
skewness and kurtosis. *The Statistician, 47*(1), 183–189.
[doi:10.1111/1467-9884.00122](https://doi.org/10.1111/1467-9884.00122)

## See also

[`skewness`](https://yelleknek.github.io/DMAR/reference/skewness.md),
[`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md)

Other descriptive statistics:
[`descriptives()`](https://yelleknek.github.io/DMAR/reference/descriptives.md),
[`skewness()`](https://yelleknek.github.io/DMAR/reference/skewness.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Normal data: excess kurtosis near zero.
set.seed(113)
kurtosis(rnorm(1000))
#> [1] -0.1595807

# Heavy-tailed data: positive excess kurtosis.
kurtosis(rt(1000, df = 4))
#> [1] 5.076089

# The classic 1:5 example: bias-corrected excess kurtosis = -1.2.
kurtosis(1:5)
#> [1] -1.2
```
