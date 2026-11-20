# Unbiased Estimate of the Population Standard Deviation Under Normality

Returns the unbiased estimate of the population standard deviation under
normality. Although the sample variance \\s^2\\ is unbiased for
\\\sigma^2\\ when computed with \\N - 1\\ in the denominator, its square
root \\s\\ is biased downward for \\\sigma\\ because the square root
function is concave (Jensen's inequality). `sd_unbiased` applies the
classical Holtzman (1950) correction factor that exactly removes that
bias under normality.

## Usage

``` r
sd_unbiased(s = NULL, N = NULL, X = NULL)
```

## Arguments

- s:

  The usual estimate of the standard deviation (the square root of the
  unbiased variance \\s^2\\). Must not be negative.

- N:

  The sample size on which `s` is based. Must exceed 1; it need not be a
  whole number.

- X:

  Optional vector of raw scores from which `s` and `N` are inferred
  (`s = sd(X)`, `N = length(X)`). Mutually exclusive with `s` and `N`.

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` value
is `"sd"` and `value` is the unbiased estimate of \\\sigma\\.

## Details

The sample variance computed with \\N - 1\\ in the denominator is
unbiased for the population variance \\\sigma^2\\, but its square root
\\s\\ is biased downward for \\\sigma\\. Under normality, the
multiplicative bias is \$\$E\[s\] \\=\\ \sigma \cdot c_N^{-1}, \qquad
c_N \\=\\ \sqrt{(N - 1)/2}\\ \cdot \Gamma((N -
1)/2)\\/\\\Gamma(N/2),\$\$ (Holtzman, 1950). Multiplying \\s\\ by
\\c_N\\ therefore yields an unbiased estimator of \\\sigma\\. The
correction is non-trivial in small samples: \\c_N\\ is about 1.064 at
\\N = 5\\, 1.028 at \\N = 10\\, 1.009 at \\N = 30\\, and is essentially
1 by \\N = 100\\ (about 1.003). For most applied work the bias of \\s\\
is small enough to ignore, but it matters when \\s\\ feeds into
downstream quantities (variance components, standardizers, planning
calculations) at small \\N\\.

Implementation note: \\c_N\\ is computed as the reciprocal of \\c_4(N) =
E\[s\]/\sigma\\, and \\c_4(N)\\ as the exponential of its logarithm. For
\\N \ge 101\\ that logarithm is its asymptotic expansion, \\-1/\[4(N -
1)\] + 1/\[24(N - 1)^3\] - \cdots\\; for a smaller \\N\\ it is the
expansion at the first of \\N + 2, N + 4, \ldots\\ that reaches 101,
carried back down by the exact recurrence \\\log c_4(N) = \log c_4(N +
2) - \log(1 + 1/\[(N - 1)(N + 1)\])/2\\, whose terms all have one sign.
Nothing cancels at any \\N\\, so the computation avoids both the
overflow of `gamma` above \\N \approx 340\\ and the loss of digits in a
difference of two `lgamma` values, which grows with \\N\\. \\c_N\\ is
accurate to two units in the last place at any \\N\\, and it is never
below 1.

## References

Holtzman, W. H. (1950). The unbiased estimate of the population variance
and standard deviation. *American Journal of Psychology*, *63*, 615–617.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`sd`](https://rdrr.io/r/stats/sd.html),
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
X <- rnorm(10, 100, 15)

# The plug-in estimate sqrt(s^2) is biased downward for sigma.
var(X)^.5
#> [1] 14.31461

# Holtzman (1950) bias-corrected estimate, supplying s and N.
sd_unbiased(s = var(X)^.5, N = length(X))
#>  term value 
#>  sd   14.717

# Equivalent call from the raw vector.
sd_unbiased(X = X)
#>  term value 
#>  sd   14.717
```
