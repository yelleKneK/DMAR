# Paired Pairwise Comparisons With Multiple-Comparison Adjustment

Computes all pairwise paired-*t* comparisons among the levels of a
within-subjects factor and returns the mean difference, paired SD,
paired *t*-statistic, degrees of freedom, raw and adjusted *p*-values,
and a confidence interval on the mean difference, one row per pair.
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

A `data.frame` (class `dmar_tbl`) with one row per pair. Columns:
`contrast` (the labeled difference, e.g. `"B - A"`), `mean_difference`,
`sd_difference`, `t_value`, `df`, `p_value` (raw), `p_adjusted`,
`ci_lower`, `ci_upper`, `n_pairs`. The table records `conf_level` in its
`conf_level` attribute, which the console prints beneath the table and
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
reports. With `bonferroni_ci = TRUE` that is the simultaneous level of
the \\m\\ intervals, the convention of
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
with `adjust = "bonferroni"` and of the post hoc intervals
([`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
and the rest), and the line beneath the table names the Bonferroni
adjustment and the level of each interval, \\1 - (1 -
\mathrm{conf\\level}) / m\\.

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
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
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
#>  contrast mean_difference sd_difference t_value df p_value p_adjusted ci_lower
#>  T2 - T1  0.966           1.157         2.640   9  0.0269  0.1147     0.138   
#>  T3 - T1  1.410           1.433         3.111   9  0.0125  0.0749     0.385   
#>  T4 - T1  1.289           1.489         2.738   9  0.0229  0.1147     0.224   
#>  T3 - T2  0.444           1.669         0.842   9  0.4217  1.0000     -0.750  
#>  T4 - T2  0.324           1.821         0.563   9  0.5874  1.0000     -0.979  
#>  T4 - T3  -0.120          1.458         -0.261  9  0.7999  1.0000     -1.163  
#>  ci_upper n_pairs
#>  1.793    10     
#>  2.435    10     
#>  2.355    10     
#>  1.638    10     
#>  1.626    10     
#>  0.923    10     
#> 
#> Confidence level: 95%

# 2. Long-format input:
long <- data.frame(
  subject = factor(rep(1:n, times = k)),
  time    = factor(rep(paste0("T", 1:k), each = n)),
  y       = as.vector(Y)
)
pairwise_within(long, subject = "subject", condition = "time", outcome = "y")
#>  contrast mean_difference sd_difference t_value df p_value p_adjusted ci_lower
#>  T2 - T1  0.966           1.157         2.640   9  0.0269  0.1147     0.138   
#>  T3 - T1  1.410           1.433         3.111   9  0.0125  0.0749     0.385   
#>  T4 - T1  1.289           1.489         2.738   9  0.0229  0.1147     0.224   
#>  T3 - T2  0.444           1.669         0.842   9  0.4217  1.0000     -0.750  
#>  T4 - T2  0.324           1.821         0.563   9  0.5874  1.0000     -0.979  
#>  T4 - T3  -0.120          1.458         -0.261  9  0.7999  1.0000     -1.163  
#>  ci_upper n_pairs
#>  1.793    10     
#>  2.435    10     
#>  2.355    10     
#>  1.638    10     
#>  1.626    10     
#>  0.923    10     
#> 
#> Confidence level: 95%
```
