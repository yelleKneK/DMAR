# Paired-Samples Randomization (Sign-Flip) Test

Computes an exact (or Monte-Carlo) sign-flip randomization test for
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
  alternative = c("two.sided", "less", "greater"),
  exact = NULL,
  n_resamples = 10000L,
  seed = NULL
)
```

## Arguments

- x, y:

  Paired numeric vectors of equal length. `NA`s are removed pairwise.

- statistic:

  One of `"mean"` (default) or `"t"`. `"mean"` uses the mean difference
  \\\bar d\\ as the test statistic. `"t"` uses the paired *t*-statistic
  \\\bar d \sqrt{n} / s_d\\, which is more robust to pair-to-pair
  variability in \\\|d_i\|\\.

- alternative:

  One of `"two.sided"` (default), `"less"`, or `"greater"`.

- exact:

  Logical. If `NULL` (default), uses exact enumeration when \\n \le 20\\
  and Monte-Carlo otherwise. If `TRUE`, forces exact enumeration (caps
  at \\n \le 25\\; above that the \\2^n\\ space is too large). If
  `FALSE`, forces Monte-Carlo.

- n_resamples:

  Number of Monte-Carlo resamples when exact enumeration is not used.
  Default `10000L`.

- seed:

  Optional integer seed for reproducibility of the Monte- Carlo branch.
  Default `NULL`, which leaves the user's current RNG state intact;
  supply an integer for reproducibility.

## Value

A `data.frame` with rows for the observed test statistic, the *p*-value,
the number of pairs, the number of randomizations evaluated, and a flag
indicating whether the test was exact or Monte-Carlo.

## Details

For small \\n\\ (default \\n \le 20\\) the test enumerates all \\2^n\\
sign patterns exactly; for larger \\n\\ a Monte-Carlo approximation is
used (default `n_resamples = 10000L`).

**Why randomization.** The randomization test makes no distributional
assumption on \\d_i\\; it only assumes that under the null, the sign of
each \\d_i\\ is arbitrary. This is exactly the inference that
pre-experimental random assignment licenses, and it is robust to
heavy-tailed differences, mixtures, and outliers.

**Exact enumeration.** For \\n \le 25\\, all \\2^n\\ sign patterns are
enumerated. The observed test statistic is compared with the full
reference distribution. The exact two-sided *p*-value is the proportion
of patterns yielding a test statistic at least as extreme (in absolute
value) as the observed.

**Monte-Carlo branch.** For larger \\n\\, `n_resamples` random sign
patterns are drawn uniformly from \\\\-1, +1\\^n\\; the Monte-Carlo
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
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`dunnett_ci()`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`scheffe_ci()`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`tukey_kramer_ci()`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
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
#>  statistic   4.8   
#>  p_value     0.0625
#>  n_pairs     5     
#>  n_evaluated 32    
#>  exact       1     

# 2. Larger n: Monte-Carlo branch.
set.seed(113)
x <- rnorm(50, 100, 15)
y <- x + rnorm(50,   5, 12)
randomization_test_paired(x, y, n_resamples = 10000L)
#>  term        value 
#>  statistic   3.71  
#>  p_value     0.0566
#>  n_pairs     50    
#>  n_evaluated 10000 
#>  exact       0     
```
