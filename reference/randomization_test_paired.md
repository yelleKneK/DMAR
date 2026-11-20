# Paired-Samples Randomization (Sign-Flip) Test

Computes an exact (or Monte Carlo) sign-flip randomization test for
paired observations \\(x_i, y_i)\\, treating the within-pair sign of
\\d_i = y_i - x_i\\ as the randomization mechanism (Fisher, 1971;
Edgington & Onghena, 2007). Under the null hypothesis of exchangeability
(\\H_0\\: the labeling of \\x\\ and \\y\\ within each pair is
arbitrary), each of the \\2^n\\ sign patterns is equally likely.

## Usage

``` r
randomization_test_paired(
  x,
  y,
  statistic = c("mean", "t"),
  alternative = c("two_sided", "less", "greater"),
  exact = NULL,
  n_resamples = 10000L,
  seed = NULL
)
```

## Arguments

- x, y:

  Paired numeric vectors of equal length. `NA`s are removed pairwise;
  the remaining values must be finite.

- statistic:

  One of `"mean"` (default) or `"t"`. `"mean"` reports the mean
  difference and `"t"` the paired *t* statistic. Over the sign patterns
  the paired *t* is a strictly increasing function of the mean
  difference, so the two give the same *p*-value and differ only in the
  statistic reported. `"t"` is refused when every difference has the
  same nonzero value (exactly, or to within the tolerance of the
  comparison), because the statistic then has a zero standard error and
  no value (see *When the studentized statistic is undefined*); when
  every difference is zero it is reported as 0 with a *p*-value of 1.

- alternative:

  One of `"two_sided"` (default; the base-R spelling `"two.sided"` is
  accepted as an alias), `"less"`, or `"greater"`.

- exact:

  Logical. If `NULL` (default), uses exact enumeration when \\n \le 20\\
  and Monte Carlo otherwise. If `TRUE`, forces exact enumeration (caps
  at \\n \le 25\\; above that the \\2^n\\ space is too large). If
  `FALSE`, forces Monte Carlo.

- n_resamples:

  A single whole number of at least 1: the number of Monte Carlo
  resamples when exact enumeration is not used. Default `10000L`.

- seed:

  Optional integer seed for reproducibility of the Monte Carlo branch.
  Default `NULL`, under which the randomizations are drawn from, and
  advance, the current random number stream; supply an integer for
  reproducibility, which is set for the call only, the caller's random
  number state being restored when it returns.

## Value

A `data.frame` (class `dmar_tbl`) with columns `term` and `value` and
five rows: `statistic` (the mean difference \\\bar d\\, or the paired
*t* under `statistic = "t"`, always finite), `p_value`, `n_pairs` (the
pairs with both scores present), `n_evaluated` (\\2^n\\ under exact
enumeration, `n_resamples` otherwise), and `exact` (1 when every sign
pattern was enumerated, 0 when patterns were sampled). No table is
returned when `statistic = "t"` is refused; the condition, of class
`dmar_studentized_statistic_undefined`, carries `reason`
(`"constant_differences"`), `alternative`, `n` (the number of pairs),
and `p_value_mean`, the exact *p*-value `statistic = "mean"` gives for
the same call.

## Details

For small \\n\\ (default \\n \le 20\\) the test enumerates all \\2^n\\
sign patterns exactly; for larger \\n\\ a Monte Carlo approximation is
used (default `n_resamples = 10000L`).

**Why randomization.** The randomization test makes no distributional
assumption on \\d_i\\; it only assumes that under the null, the sign of
each \\d_i\\ is arbitrary. This is exactly the inference that
pre-experimental random assignment licenses. Because the reference
distribution is built from the observed differences, the Type I error
rate under that null does not exceed the nominal level whatever the
distribution of the differences, and for any continuous distribution,
heavy-tailed, a mixture, or with outliers, it is the same attainable
level: under exact enumeration, the largest multiple of \\1/2^{n-1}\\
(two-sided) or \\1/2^n\\ (one-sided) that does not exceed the nominal
level (6/128 = .047 two-sided at .05 with 8 pairs). When different sign
patterns can give the same statistic, as with discrete scores, it can be
lower.

**Exact enumeration.** All \\2^n\\ sign patterns are enumerated for
\\n\\ of 20 or fewer by default, and for \\n\\ up to 25 when
`exact = TRUE`. The observed test statistic is compared with the full
reference distribution. The exact two-sided *p*-value is the proportion
of patterns yielding a test statistic at least as extreme (in absolute
value) as the observed. A pattern whose statistic equals the observed
one counts as at least as extreme; the comparison allows for the
rounding of the arithmetic, with a bound proportional to the size of the
scores, so that patterns tied in exact arithmetic are counted as tied.

**The two statistics order the patterns alike.** Write \\S\\ for the sum
of the signed differences and \\Q = \sum d_i^2\\, which no sign pattern
changes. The paired *t* is \\S \sqrt{n - 1} / \sqrt{n Q - S^2}\\, which
increases with \\S\\ because \\Q\\ is fixed. Both statistics therefore
rank the patterns as \\S\\ does, and the function counts on \\S\\ for
both.

**When the studentized statistic is undefined.** The studentized
statistic is a ratio, the mean difference over its standard error, and
it has a value only when that standard error is positive. If every
paired difference is the same nonzero number (exactly, or to within the
tolerance of the comparison), the standard error of the observed data is
0 and `statistic = "t"` is refused with a condition of class
`dmar_studentized_statistic_undefined`. The observed data are then at
the boundary of the reference set: among the \\2^n\\ sign patterns only
the observed one and its mirror image have a zero standard error. The
refusal concerns the statistic and not the evidence.
`statistic = "mean"` is defined for every arrangement and gives the
smallest *p*-value the design allows in the direction of the data,
\\2/2^n\\ two-sided and \\1/2^n\\ one-sided; the message states the
value, and the condition carries it as `p_value_mean`. The refusal is
made before any sign pattern is enumerated or drawn, so a refused Monte
Carlo call leaves the random number stream where it was.

When the observed statistic is defined, some arrangements in the
reference set can still have a zero standard error: with paired
differences of equal size and mixed signs, the two patterns that give
them all one sign. Those patterns are the two most extreme, which is
where the mean difference also places them, and since the function
counts on the signed sum \\S\\ for both statistics they are counted
without a standard error being formed. Refusing them would refuse
ordinary coarsely measured data. For differences of 2, 2, 2, \\-2\\, and
2 the paired *t* is 1.5 and the *p*-value, 12/32 = 0.375, is the one the
mean difference gives. When every difference is zero, every sign pattern
is the same data set; the statistic is reported as 0 and the *p*-value
is 1.

**Monte Carlo branch.** For larger \\n\\, `n_resamples` random sign
patterns are drawn uniformly from \\\\-1, +1\\^n\\; the Monte Carlo
*p*-value uses the standard (\\1 + \mathrm{count}\\) / (\\1 + B\\)
plug-in to avoid *p* = 0.

## References

Edgington, E. S., & Onghena, P. (2007). *Randomization tests* (4th ed.).
Chapman & Hall/CRC.

Fisher, R. A. (1971). *The design of experiments* (9th ed., reprint).
Hafner.

Pitman, E. J. G. (1937). Significance tests which may be applied to
samples from any populations. *Supplement to the Journal of the Royal
Statistical Society, 4*(1), 119–130.

## See also

[`t.test`](https://rdrr.io/r/stats/t.test.html) (parametric paired
test),
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_dunnett_t3()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
[`ci_games_howell()`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Small-n exact: Bayley scores on twin pairs.
control <- c(95, 102,  98, 107, 105)
treat   <- c(102, 108, 100, 112, 109)
randomization_test_paired(control, treat)
#>  term        value 
#>  statistic   4.800 
#>  p_value     0.0625
#>  n_pairs     5     
#>  n_evaluated 32    
#>  exact       1     

# 2. Larger n: Monte Carlo branch.
set.seed(113)
x <- rnorm(50, 100, 15)
y <- x + rnorm(50,   5, 12)
randomization_test_paired(x, y, n_resamples = 10000L)
#>  term        value 
#>  statistic   3.713 
#>  p_value     0.0566
#>  n_pairs     50    
#>  n_evaluated 10000 
#>  exact       0     

# 3. Every pair differs by the same amount. The mean difference is
#    defined, and gives the smallest p-value 5 pairs allow. With
#    statistic = "t" the same call is refused, with a message that says
#    why: the observed standard error is 0, so the studentized statistic
#    is undefined.
before <- c(10, 12, 9, 15, 11)
randomization_test_paired(before, before + 2)
#>  term        value 
#>  statistic   2.000 
#>  p_value     0.0625
#>  n_pairs     5     
#>  n_evaluated 32    
#>  exact       1     
```
