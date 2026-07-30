# Provides the Critical Value for the Tukey Honestly Significant Difference (HSD) Test

Provides the Critical Value for the Tukey Honestly Significant
Difference (HSD) Test

## Usage

``` r
cv_tukey_hsd(alpha_level, df, groups, verbose = TRUE)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the false positive rate). For the Tukey HSD
  test, the full `alpha_level` applies to the upper tail of the
  Studentized range distribution (see Details).

- df:

  The error degrees of freedom from the ANOVA (a positive number;
  typically `N - groups`).

- groups:

  The number of groups whose means are being compared (an integer of at
  least 2).

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value in a output style (a `data.frame` with one
row per critical value, following the format used by
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md) and
[`cv_z`](https://yelleknek.github.io/DMAR/reference/cv_z.md)).

## Details

The Tukey HSD test compares all pairs of group means using the
Studentized range distribution
([`qtukey`](https://rdrr.io/r/stats/Tukey.html)). Because the
Studentized range is the absolute difference between the largest and
smallest sample means (scaled by a standard error), it is non-negative
and the associated distribution has support on \\\[0, \infty)\\. As a
consequence, the Type I error rate `alpha_level` is *not* split between
two tails in the way it is for the (symmetric) *t*- and
*z*-distributions in
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md) and
[`cv_z`](https://yelleknek.github.io/DMAR/reference/cv_z.md); rather,
the full `alpha_level` applies to the upper tail.

The reported critical value is on the scale used for pairwise
comparisons of group means, i.e., \\q\_{1-\alpha, k, df} / \sqrt{2}\\,
where \\k\\ is the number of groups. A pair of means is declared
significantly different when the absolute standardized difference
between them exceeds this critical value.

The Tukey HSD critical value is a quantile of the Studentized range
distribution, which base R supplies through
[`qtukey`](https://rdrr.io/r/stats/Tukey.html), so unlike
[`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
and [`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md)
this function needs no multivariate distribution machinery and does not
require the mvtnorm package. Maxwell, Delaney, and Kelley (2027, Chapter
5) develop the Tukey method as the procedure for all-pairwise
comparisons within the multiple-comparisons problem.

## References

Tukey, J. W. (1953). *The problem of multiple comparisons*. Unpublished
manuscript, Princeton University.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem,
where the Tukey method for all-pairwise comparisons is developed.)

## See also

[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_z`](https://yelleknek.github.io/DMAR/reference/cv_z.md),
[`TukeyHSD`](https://rdrr.io/r/stats/TukeyHSD.html),
[`qtukey`](https://rdrr.io/r/stats/Tukey.html)

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Following the all-pairwise comparisons setting of Maxwell, Delaney, and
# Kelley (2027, Chapter 5): every pair of k = 3 group means is compared,
# with 27 error degrees of freedom, holding the family-wise error rate at
# .05.
cv_tukey_hsd(alpha_level = .05, df = 27, groups = 3)
#>  term     value area_less area_greater
#>  upper_cv 2.48  0.95      0.05        

# Using the built-in PlantGrowth dataset (3 groups, N = 30).
fit <- aov(weight ~ group, data = PlantGrowth)
cv_tukey_hsd(
  alpha_level  = .05,
  df     = df.residual(fit),
  groups = nlevels(PlantGrowth$group)
)
#>  term     value area_less area_greater
#>  upper_cv 2.48  0.95      0.05        

# A more stringent alpha with a simple (non-verbose) result.
cv_tukey_hsd(alpha_level = .01, df = 27, groups = 3, verbose = FALSE)
#>  term     value
#>  upper_cv 3.18 
```
