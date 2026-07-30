# Paired Pairwise Comparisons With Multiple-Comparison Adjustment

Computes all pairwise paired-*t* comparisons among the levels of a
within-subjects factor and returns the mean difference, paired SD,
paired *t*-statistic, degrees of freedom, raw and adjusted *p*-values,
and a confidence interval on the mean difference, all in tidy long form.
*p*-values are adjusted across comparisons by the user-specified method
(Bonferroni, Holm, Hochberg, Hommel, BH, BY, or none).

## Usage

``` r
pairwise_within(
  data,
  subject = NULL,
  condition = NULL,
  outcome = NULL,
  adjust = c("holm", "bonferroni", "hochberg", "hommel", "BH", "BY", "none"),
  conf_level = 0.95,
  bonferroni_ci = FALSE
)
```

## Arguments

- data:

  Either an \\n \times k\\ wide numeric matrix / data.frame (one row per
  subject, one column per condition), or a long-format `data.frame`
  together with `subject`, `condition`, and `outcome` column names.

- subject:

  Long-format only: character name of the subject-id column.

- condition:

  Long-format only: character name of the within- subjects factor
  column.

- outcome:

  Long-format only: character name of the response column.

- adjust:

  Multiple-comparison adjustment method. One of `"holm"` (default),
  `"bonferroni"`, `"hochberg"`, `"hommel"`, `"BH"`, `"BY"`, or `"none"`.
  Passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- conf_level:

  Family-wise confidence level for the per-pair CIs. Default `0.95`. CIs
  are computed at the per-pair nominal level (\\1 - \alpha/m\\ under
  Bonferroni; `conf_level` otherwise).

- bonferroni_ci:

  Logical. If `TRUE`, the CIs use a per-pair confidence level of \\1 -
  (1 - \mathrm{conf\\level}) / m\\ to give simultaneous `conf_level`
  coverage across the \\m\\ comparisons (Bonferroni-corrected CIs).
  Default `FALSE`.

## Value

A `data.frame` with one row per pair. Columns: `contrast` (the labeled
difference, e.g. `"B - A"`), `mean_difference`, `sd_difference`,
`t_statistic`, `df`, `p_value` (raw), `p_adjusted`, `lower_limit`,
`upper_limit`, `n_pairs`.

## Details

**Why a paired pairwise.**
[`stats::pairwise.t.test()`](https://rdrr.io/r/stats/pairwise.t.test.html)
returns a square matrix of *p*-values, which doesn't compose with the
rest of the DMAR pipeline. This function returns one row per comparison,
matching the `data.frame(term, value)` style used elsewhere.

**CI scale.** CIs are on the mean-difference scale (unstandardized).
When `bonferroni_ci = TRUE`, the per-pair confidence level is \\1 - (1 -
\mathrm{conf\\level}) / m\\, giving Bonferroni-style simultaneous
coverage. The CI is built from the paired-*t* distribution with \\n -
1\\ degrees of freedom.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 11.)

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.

## See also

[`pairwise.t.test`](https://rdrr.io/r/stats/pairwise.t.test.html),
[`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md)

Other within-subjects analysis:
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
[`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)

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
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
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
# 1. Wide-format input: 4 timepoints x 10 subjects.
set.seed(113)
n <- 10; k <- 4
Y <- matrix(rnorm(n * k, 0, 1), n, k) +
     matrix(rep(seq(0, 0.9, length.out = k), n), n, k, byrow = TRUE) +
     rnorm(n, 0, 1.5)
colnames(Y) <- paste0("T", 1:k)
pairwise_within(Y)
#>   contrast mean_difference sd_difference t_statistic df    p_value p_adjusted
#> 1  T2 - T1       0.9655205      1.156502   2.6400670  9 0.02691082  0.1146525
#> 2  T3 - T1       1.4098513      1.432896   3.1114203  9 0.01249022  0.0749413
#> 3  T4 - T1       1.2894819      1.489415   2.7377861  9 0.02293049  0.1146525
#> 4  T3 - T2       0.4443308      1.669023   0.8418683  9 0.42166771  1.0000000
#> 5  T4 - T2       0.3239614      1.820780   0.5626468  9 0.58741461  1.0000000
#> 6  T4 - T3      -0.1203694      1.458101  -0.2610530  9 0.79992698  1.0000000
#>   lower_limit upper_limit n_pairs
#> 1   0.1382085   1.7928325      10
#> 2   0.3848194   2.4348832      10
#> 3   0.2240186   2.3549452      10
#> 4  -0.7496162   1.6382778      10
#> 5  -0.9785459   1.6264686      10
#> 6  -1.1634317   0.9226928      10

# 2. Long-format input:
long <- data.frame(
  subject = factor(rep(1:n, times = k)),
  time    = factor(rep(paste0("T", 1:k), each = n)),
  y       = as.vector(Y)
)
pairwise_within(long, subject = "subject", condition = "time", outcome = "y")
#>   contrast mean_difference sd_difference t_statistic df    p_value p_adjusted
#> 1  T2 - T1       0.9655205      1.156502   2.6400670  9 0.02691082  0.1146525
#> 2  T3 - T1       1.4098513      1.432896   3.1114203  9 0.01249022  0.0749413
#> 3  T4 - T1       1.2894819      1.489415   2.7377861  9 0.02293049  0.1146525
#> 4  T3 - T2       0.4443308      1.669023   0.8418683  9 0.42166771  1.0000000
#> 5  T4 - T2       0.3239614      1.820780   0.5626468  9 0.58741461  1.0000000
#> 6  T4 - T3      -0.1203694      1.458101  -0.2610530  9 0.79992698  1.0000000
#>   lower_limit upper_limit n_pairs
#> 1   0.1382085   1.7928325      10
#> 2   0.3848194   2.4348832      10
#> 3   0.2240186   2.3549452      10
#> 4  -0.7496162   1.6382778      10
#> 5  -0.9785459   1.6264686      10
#> 6  -1.1634317   0.9226928      10
```
