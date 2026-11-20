# Coefficient of Variation (Sample Value or Nearly Unbiased Estimate)

Computes the sample coefficient of variation \\k = s / M\\, the sample
standard deviation divided by the sample mean, or, optionally, the
nearly unbiased estimate of Kelley (2007, Equation 4) under normality.
Either supply a precomputed `cv` or supply the raw `mean` and `sd`; with
`unbiased = TRUE` the value is multiplied by \\(1 + 1/(4 N))\\. The
(biased) sample coefficient of variation (the default,
`unbiased = FALSE`) is the form usually reported. To accompany it with a
confidence interval, use
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md).

## Usage

``` r
cv(cv = NULL, mean = NULL, sd = NULL, N = NULL, unbiased = FALSE)
```

## Arguments

- cv:

  The sample coefficient of variation, \\k = s / M\\, a single finite
  number. Optional; if supplied, `mean` and `sd` must not be.

- mean:

  Sample mean. A single finite number other than 0. A negative mean
  gives a negative coefficient;
  [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)
  requires a positive one.

- sd:

  Sample standard deviation, using \\N - 1\\ in the variance
  denominator. A single finite number that is not negative.

- N:

  Sample size. Required when `unbiased = TRUE`, where it must exceed 1.

- unbiased:

  Logical. If `TRUE`, multiplies the sample value by \\(1 + 1/(4 N))\\.
  The result is nearly, not exactly, unbiased (see Details). Default
  `FALSE`.

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` value
is `"cv"` and `value` is either the sample coefficient (default) or the
nearly unbiased estimate (when `unbiased = TRUE`).

## Details

The sample coefficient of variation \\k = s / M\\ is the form usually
reported in applied work. It is what this function returns by default
(`unbiased = FALSE`). A point estimate is most informative when paired
with a confidence interval;
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) computes
one for the population value \\\kappa\\.

Under normality the sample coefficient is biased. Its expected value is
approximately \\\kappa \[1 - 1/(4 (N - 1)) + \kappa^2/N\]\\, the first
and second order terms of Expression 3 of Sokal and Braumann (1980, p.
51), which Kelley (2007) reproduces as his Equation 3. The term \\-1/(4
(N - 1))\\ comes from the sample standard deviation and the term
\\\kappa^2/N\\ from the sample mean in the denominator, so \\k\\
underestimates \\\kappa\\ on average when \\\kappa\\ is below about one
half and overestimates it when \\\kappa\\ is larger. Multiplying \\k\\
by \\(1 + 1/(4 N))\\, a correction that Sokal and Braumann (1980) and
Kelley (2007) credit to Haldane (1955), offsets the first term and not
the second. The relative bias that remains is about \\\kappa^2/N - 9/(32
N^2)\\; the second of these terms is derived from the expected value of
\\s\\ under normality and is not taken from those sources. Kelley (2007,
Equation 4) calls the result a nearly unbiased estimate and recommends
it when a point estimate is reported; the sample value \\k\\ is the one
used to form the confidence interval.

At \\\kappa = 0.1\\ and \\N = 10\\ the relative bias is about \\-2.6\\
percent for \\k\\ and \\-0.2\\ percent for the corrected value; at
\\\kappa = 0.3\\ and \\N = 10\\ it is about \\-1.8\\ percent and
\\+0.6\\ percent. With 10 cases the two are about equally biased, in
opposite directions, at a coefficient near 0.38 (near 0.35 in large
samples, and above 0.4 with 5 cases or fewer), and for a larger
\\\kappa\\ the corrected value is the more biased of the two. These
statements describe a mean that stays away from 0: under an unrestricted
normal distribution the expectation of \\1/M\\ does not exist as an
ordinary integral.

For confidence intervals on \\\kappa\\,
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) inverts
the noncentral *t* distribution of \\\sqrt{N}/k\\ (Johnson & Welch,
1940; McKay, 1932; Kelley, 2007). McKay (1932) and Vangel (1996) also
give approximations based on the chi square distribution, which `ci_cv`
does not use. The first-order variance of \\k\\ under normality,
\\\kappa^2 (1/2 + \kappa^2) / (N - 1)\\, is in
[`var_cv`](https://yelleknek.github.io/DMAR/reference/var_cv.md), whose
page relates it to Expression 5 of Sokal and Braumann (1980, p. 52).

## References

Haldane, J. B. S. (1955). The measurement of variation. *Evolution,
9*(4), 484.
[doi:10.1111/j.1558-5646.1955.tb01559.x](https://doi.org/10.1111/j.1558-5646.1955.tb01559.x)

Johnson, N. L., & Welch, B. L. (1940). Applications of the non-central
*t*-distribution. *Biometrika, 31*(3–4), 362–389.
[doi:10.1093/biomet/31.3-4.362](https://doi.org/10.1093/biomet/31.3-4.362)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

McKay, A. T. (1932). Distribution of the coefficient of variation and
the extended *t* distribution. *Journal of the Royal Statistical
Society, 95*(4), 695–698.

Sokal, R. R., & Braumann, C. A. (1980). Significance tests for
coefficients of variation and variability profiles. *Systematic Zoology,
29*(1), 50–66.
[doi:10.1093/sysbio/29.1.50](https://doi.org/10.1093/sysbio/29.1.50)

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
#>  cv   0.150

# 2. Nearly unbiased estimate at N = 50; the multiplier is close to 1
#    at this sample size.
cv(mean = 100, sd = 15, N = 50, unbiased = TRUE)
#>  term value
#>  cv   0.151

# 3. The same at N = 10, where the multiplier is larger.
cv(cv = .15, N = 10, unbiased = TRUE)
#>  term value
#>  cv   0.154
```
