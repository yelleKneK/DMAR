# Coefficient of Variation (Biased or Unbiased Estimator)

Computes the sample coefficient of variation \\\hat\kappa = s / \bar Y\\
or, optionally, its first-order bias-corrected counterpart under
normality. Either supply a precomputed `cv` or supply the raw `mean` and
`sd`; with `unbiased = TRUE` the value is multiplied by the small-sample
correction \\(1 + 1/(4 N))\\. The (biased) sample coefficient of
variation (the default, `unbiased = FALSE`) is the form usually
reported. To accompany it with a confidence interval, use
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md).

## Usage

``` r
cv(cv = NULL, mean = NULL, sd = NULL, N = NULL, unbiased = FALSE)
```

## Arguments

- cv:

  The sample coefficient of variation, \\s/\bar Y\\. Optional; if
  supplied, `mean` and `sd` must not be.

- mean:

  Sample mean. Numeric scalar.

- sd:

  Sample standard deviation, using \\N - 1\\ in the variance
  denominator. Numeric scalar.

- N:

  Sample size. Required when `unbiased = TRUE`.

- unbiased:

  Logical. If `TRUE`, applies the first-order small-sample bias
  correction \\(1 + 1/(4 N))\\ to the plug-in estimator. Default
  `FALSE`.

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` value
is `"cv"` and `value` is either the plug-in estimator (default) or the
first-order bias-corrected estimator (when `unbiased = TRUE`).

## Details

The plug-in estimator \\\hat\kappa = s / \bar Y\\ is the workhorse
coefficient of variation in applied work and is the form usually
reported. It is what this function returns by default
(`unbiased = FALSE`). A point estimate is most informative when paired
with a confidence interval;
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) computes
one for \\\kappa\\. Under normality, however, the plug-in estimator is
biased downward, and the leading-order expansion is \\E\[\hat\kappa\] =
\kappa (1 - 1/(4 N)) + O(N^{-2})\\ (Sokal & Rohlf, 1995). The bias is
negligible at \\N\\ above about 100 but is non-trivial in small samples,
where multiplying the plug-in value by \\(1 + 1/(4 N))\\ removes the
leading-order term. The `unbiased = TRUE` option applies that
correction.

For confidence intervals on \\\kappa\\, the McKay (1932) noncentral *t*
based interval is implemented in
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md); the
corresponding asymptotic variances are in
[`var_cv`](https://yelleknek.github.io/DMAR/reference/var_cv.md). The
Vangel (1996) small-sample refinement of McKay's interval is a further
option described in the literature; it matters more than the bias
correction in this function when \\\kappa\\ is larger than about 0.3
(Kelley, 2007).

## References

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3.)

McKay, A. T. (1932). Distribution of the coefficient of variation and
the extended *t* distribution. *Journal of the Royal Statistical
Society, 95*(4), 695–698.

Sokal, R. R., & Rohlf, F. J. (1995). *Biometry: The principles and
practice of statistics in biological research* (3rd ed.). W. H. Freeman.

Vangel, M. G. (1996). Confidence intervals for a normal coefficient of
variation. *The American Statistician, 50*(1), 21–26.
[doi:10.1080/00031305.1996.10473537](https://doi.org/10.1080/00031305.1996.10473537)

## See also

[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`var_cv`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Point estimate from raw mean and SD.
cv(mean = 100, sd = 15)
#>  term value
#>  cv   0.15 

# 2. Bias-corrected estimate at N = 50; the correction is small but
#    non-negligible at this sample size.
cv(mean = 100, sd = 15, N = 50, unbiased = TRUE)
#>  term value
#>  cv   0.151

# 3. Bias correction at N = 10; the correction is larger here.
cv(cv = .15, N = 10, unbiased = TRUE)
#>  term value
#>  cv   0.154
```
