# Provides the Critical Value(s) for a *t*-distribution

Provides the Critical Value(s) for a *t*-distribution

## Usage

``` r
cv_t(
  alpha_level,
  df,
  alternative = "not_equal",
  alpha_lower,
  alpha_upper,
  ncp = 0,
  verbose = TRUE
)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the false positive rate).

- df:

  The number of degrees of freedom (a positive number)

- alternative:

  The type of alternative hypothesis of interest.

- alpha_lower:

  The error rate on the lower (negative) side of the distribution.

- alpha_upper:

  The error rate on the upper (positive) side of the distribution.

- ncp:

  The noncentral parameter (if zero, the default, it is the central
  *t*-distribution).

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value(s), based on the input specifications, in a
output style.

## Details

Though a noncentral parameter can be included, that would not be done
for a standard null hypothesis significance test.

With a nonzero `ncp` the critical values are quantiles of the noncentral
*t* distribution, computed exactly. R's
[`qt()`](https://rdrr.io/r/stats/TDist.html) and
[`pt()`](https://rdrr.io/r/stats/TDist.html) with `ncp` switch to a
normal approximation above a noncentrality of 37.62 (and above 33 with
more than 2,000 degrees of freedom), where they can be far off: the .975
quantile with `ncp = 50` and `df = 3` is 186.49, and
[`qt()`](https://rdrr.io/r/stats/TDist.html) gives 429.18. There the
quantiles and the tail areas are computed by integrating over the
sampling distribution of the standard deviation; inside that range they
are R's. With fewer than one degree of freedom (`df` may be any positive
number) [`pt()`](https://rdrr.io/r/stats/TDist.html) is not accurate
inside that range either. At the exact quantiles of probability .0005 to
.9995, with noncentralities up to 37.5 in absolute value, the tail
probability it gave was off by as much as .14 with `df = 0.05`, .07 with
0.1, .01 with 0.2, .0013 with 0.3, \\3.8 \times 10^{-5}\\ with 0.5, and
\\2 \times 10^{-11}\\ with 0.8, so with `df < 1` every noncentral
quantile and tail area is integrated. The critical values then differ
from [`qt()`](https://rdrr.io/r/stats/TDist.html)'s by orders of
magnitude with `df` of 0.3 or less (and by as much as 15 percent with
0.5): with `df = 0.1`, `ncp = 0.5`, and `alpha_level = .05` they are
\\-1.0428 \times 10^{10}\\ and \\5.7095 \times 10^{13}\\, where
[`qt()`](https://rdrr.io/r/stats/TDist.html) gives -33554432 and
305284064. With so few degrees of freedom the tails are heavy enough
that a critical value can lie beyond the largest number R represents
(about \\1.8 \times 10^{308}\\); it is then returned as `-Inf` or `Inf`,
as [`qt()`](https://rdrr.io/r/stats/TDist.html) returns it with
`ncp = 0` (with `df = 0.001` one or both of the critical values are
infinite at every `alpha_level` from .000001 to .4 and noncentrality
from -100 to 100 tried).

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A basic call for finding critical values with equal area in the two tails.
cv_t(alpha_level = .05, df = 13)
#>  term     value  area_less area_greater
#>  lower_cv -2.160 0.025     0.975       
#>  upper_cv 2.160  0.975     0.025       

# A basic call for a single-sided confidence interval (for "a greater than" alternative hypothesis)
cv_t(alpha_level = .05, df = 13, alternative = "greater")
#>  term     value area_less area_greater
#>  lower_cv -Inf  0.000     1.000       
#>  upper_cv 1.771 0.950     0.050       

# A single-sided confidence interval (for "a greater than" alternative hypothesis); simple output.
cv_t(alpha_lower = 0, alpha_upper = .05, df = 13, verbose = FALSE)
#>  term     value
#>  lower_cv -Inf 
#>  upper_cv 1.771

# For a nonsymmetric 95% confidence interval.
cv_t(alpha_lower = .01, alpha_upper = .04, df = 13)
#>  term     value  area_less area_greater
#>  lower_cv -2.650 0.010     0.990       
#>  upper_cv 1.899  0.960     0.040       
```
