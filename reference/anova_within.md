# One Way Within-Subjects ANOVA With Sphericity Diagnostics and Corrections

Performs the univariate one-way within-subjects *F* test together with
Mauchly's test of sphericity and the three standard
\\\varepsilon\\-corrected *p*-values (Greenhouse-Geisser, Huynh-Feldt,
and lower-bound). Returns everything in a single tidy `data.frame` so
the user can decide which adjustment to report.

## Usage

``` r
anova_within(x, id = NULL, time = NULL, outcome = NULL)
```

## Arguments

- x:

  Either an \\n \times k\\ numeric matrix or `data.frame` (rows =
  subjects, columns = repeated measurements); *or* a long-format
  `data.frame` together with `id`, `time`, and `outcome` column names.

- id:

  Column name in `x` identifying the subject when `x` is in long format
  (`NULL` otherwise).

- time:

  Column name in `x` identifying the within-subjects factor level when
  `x` is in long format (`NULL` otherwise).

- outcome:

  Column name in `x` identifying the dependent variable when `x` is in
  long format (`NULL` otherwise).

## Value

A `data.frame` with one row per reported *F* test: `adjustment`
(`"none"`, `"Greenhouse-Geisser"`, `"Huynh-Feldt"`, `"lower_bound"`),
`F_value`, `df_1`, `df_2`, `p_value`, and `epsilon` (the correction
factor used; `NA` for the unadjusted row). `attr(<output>, "mauchly")`
contains the row from
[`mauchly_test`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
and the partial \\\eta^2\\ is attached as
`attr(<output>, "partial_eta_squared")`.

## Details

The unadjusted within-subjects *F* statistic is the same regardless of
sphericity; corrections shrink the numerator and denominator degrees of
freedom by a factor of \\\hat\varepsilon \in \[1/(k - 1),\\ 1\]\\, and
the *p*-value is recomputed against the adjusted reference *F*
distribution. When Mauchly's test rejects, prefer the
Huynh-Feldt-corrected *p*-value (less conservative than
Greenhouse-Geisser).

For multi-factor within-subjects designs or mixed designs, fit the model
with `stats::aov(... + Error(id/within))` or with
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) directly.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 11.)

## See also

[`mauchly_test`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`epsilon_corrections`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
[`aov`](https://rdrr.io/r/stats/aov.html)

Other within-subjects analysis:
[`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
[`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)

Other hypothesis tests:
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
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
# Simulated within-subjects data with no real effect.
set.seed(113)
Y <- matrix(rnorm(20 * 4), nrow = 20)
anova_within(Y)
#>  adjustment         F_value df_1 df_2 p_value epsilon
#>  none               0.618   3    57   0.6062  <NA>   
#>  Greenhouse-Geisser 0.618   2.36 44.9 0.5696  0.787  
#>  Huynh-Feldt        0.618   2.72 51.7 0.5911  0.906  
#>  lower_bound        0.618   1    19   0.4415  0.333  

# Built-in within-subjects example: nlme::Orthodont (distance ~ age).
res <- anova_within(nlme::Orthodont,
                    id = "Subject", time = "age", outcome = "distance")
res
#>  adjustment         F_value df_1 df_2 p_value  epsilon
#>  none               38      3    78   < 0.0001 <NA>   
#>  Greenhouse-Geisser 38      2.63 68.4 < 0.0001 0.877  
#>  Huynh-Feldt        38      2.95 76.8 < 0.0001 0.984  
#>  lower_bound        38      1    26   < 0.0001 0.333  
attr(res, "mauchly")
#>  W     statistic df p_value n_subjects n_levels method                      
#>  0.758 6.85      5  0.2326  27         4        Mauchly's test of sphericity
attr(res, "partial_eta_squared")
#> [1] 0.5940013
```
