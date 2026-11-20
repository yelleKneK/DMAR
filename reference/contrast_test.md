# Tests One or More Contrasts of Group Means in a One-Way Design

Given a fitted one-way [`aov`](https://rdrr.io/r/stats/aov.html) or
[`lm`](https://rdrr.io/r/stats/lm.html) object and a set of contrast
weights, computes for every contrast the estimate \\\hat{\psi} = \sum_i
c_i \bar{Y}\_i\\, its standard error, *t*-statistic, degrees of freedom,
two-sided *p*-value, and confidence interval. Supports several common
multiple-comparison adjustments and either equal-variance (pooled) or
Welch-style unequal-variance inference.

## Usage

``` r
contrast_test(
  object,
  contrasts = "pairwise",
  adjust = "none",
  conf_level = 0.95,
  var_equal = TRUE
)
```

## Arguments

- object:

  A fitted [`aov`](https://rdrr.io/r/stats/aov.html) or
  [`lm`](https://rdrr.io/r/stats/lm.html) object for a one-way design (a
  single grouping factor on the right-hand side of the formula).

- contrasts:

  Specification of one or more contrasts. Any of:

  `"pairwise"` (default)

  :   All pairwise comparisons among the group means.

  a named list of numeric vectors

  :   Each vector is one contrast and its name is used as the row label.

  a numeric matrix

  :   Each row is one contrast; `rownames`, if present, are used as
      labels.

  a numeric vector

  :   Treated as a single contrast.

  Each contrast vector must have length equal to the number of groups,
  and the weights are typically chosen to sum to zero.

- adjust:

  Multiple-comparison adjustment. One of `"none"` (default),
  `"bonferroni"`, `"scheffe"`, `"tukey"` (pairwise contrasts only), or
  any of the sequential methods supported by
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) (`"holm"`,
  `"hochberg"`, `"BH"`, `"BY"`).

- conf_level:

  Confidence level for the interval (default `0.95`).

- var_equal:

  Logical. If `TRUE` (default), uses the pooled error variance
  \\\mathit{MS}\_{\text{error}}\\ and the residual degrees of freedom
  from `object`. If `FALSE`, uses each group's own sample variance and a
  Welch-Satterthwaite approximate \\df\\ per contrast.

## Value

A `data.frame` with one row per contrast and columns `contrast`,
`estimate`, `se`, `t_value`, `df`, `p_value`, `p_adjusted`, `ci_lower`,
and `ci_upper`. The adjustment, confidence level, and variance
assumption are stored as `attr(*, "adjust")`, `attr(*, "conf_level")`,
and `attr(*, "var_equal")`. With `adjust` `"bonferroni"`, `"scheffe"`,
or `"tukey"`, the recorded level is the simultaneous level of the family
of intervals, the convention of the post hoc intervals and of
[`pairwise_within`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md)
with `bonferroni_ci = TRUE`; under Bonferroni the line beneath the table
also names the level of each interval, \\1 - \alpha / m\\. The table
prints through the
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
display layer and works with
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

**Test statistic.** For a contrast with weights \\c_1, \ldots, c_k\\
(\\k\\ = number of groups), the estimate is \\\hat{\psi} = \sum_i c_i
\bar{Y}\_i\\. Under equal variances, the standard error is
\\\sqrt{\mathit{MS}\_{\text{error}} \sum_i c_i^2 / n_i}\\ with \\df =
N - k\\; under unequal variances, the standard error is \\\sqrt{\sum_i
c_i^2 s_i^2 / n_i}\\ with the Welch-Satterthwaite df,
\$\$df\_{\text{Welch}} = \frac{\left(\sum_i c_i^2 s_i^2 /
n_i\right)^2}{\sum_i (c_i^2 s_i^2 / n_i)^2 / (n_i - 1)}.\$\$ The
unadjusted *p*-value is two-sided based on the *t* reference
distribution.

**Adjustments.** The `p_adjusted` and confidence interval critical value
are computed as follows.

- `"none"`: no adjustment; the CI uses \\t\_{1-\alpha/2,df}\\.

- `"bonferroni"`: \\p\_{\text{adj}} = \min(1, m\\ p)\\ for \\m\\
  contrasts, with CI based on \\t\_{1-\alpha/(2m),df}\\.

- `"scheffe"`: appropriate for any contrast (or family of contrasts).
  \\p\_{\text{adj}}\\ comes from the upper tail of an *F* reference
  distribution applied to \\t^2 / (k-1)\\, and the CI uses
  \\\sqrt{(k-1)\\ F\_{1-\alpha,\\k-1,df}}\\.

- `"tukey"`: requires every contrast to be pairwise. Uses the
  studentized range distribution \\Q\_{k, df}\\, so that
  \\p\_{\text{adj}} = P(Q\_{k, df} \> \|t\|\sqrt{2})\\ and the CI uses
  \\q\_{1-\alpha,\\k,df} / \sqrt{2}\\. Both come from base R's
  [`ptukey`](https://rdrr.io/r/stats/Tukey.html) and
  [`qtukey`](https://rdrr.io/r/stats/Tukey.html) where those are
  accurate, and from the exact distribution where they are not (at small
  df, for example fewer than 9 df with up to 6 groups; from 10000 to
  \\10^8\\ df; at confidence levels above .999; with more than 50
  groups, below the .75 point of the distribution (.80 with more than
  1000 groups); and wherever `qtukey` stops short of its root; with 20
  groups and 2 df `qtukey(.99, 20, 2)` is 44.526 where the exact value
  is 37.943);
  [`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
  describes the computation and the regions.

- `"holm"`, `"hochberg"`, `"BH"`, `"BY"`:
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) is applied to the
  unadjusted *p*-values; the CI uses the unadjusted \\t\\-critical value
  because these methods do not give simultaneous CIs in closed form.

**Variance assumption with adjustments.** The Tukey and Scheffé
procedures assume equal variances; combining them with
`var_equal = FALSE` is at the user's risk (the resulting Type I error
rate is no longer guaranteed). For unequal variances, common
alternatives are Games-Howell (Tukey-style) and Brown-Forsythe
(Scheffé-style); these are not currently supported here.

**Scope.** Only one-way designs are supported in v1 (one outcome, one
grouping factor). Multi-way designs throw an informative error.

## References

Hsu, J. C. (1996). *Multiple comparisons: Theory and methods*. Chapman &
Hall.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Scheffe, H. (1953). A method for judging all contrasts in the analysis
of variance. *Biometrika, 40*, 87–104.

Smith, J. E., Meyers, R. J., & Delaney, H. D. (1998). The community
reinforcement approach with homeless alcohol-dependent individuals.
*Journal of Consulting and Clinical Psychology, 66*(3), 541–548.

Tukey, J. W. (1953). The problem of multiple comparisons. Unpublished
manuscript, Princeton University.

## See also

[`TukeyHSD`](https://rdrr.io/r/stats/TukeyHSD.html),
[`pairwise.t.test`](https://rdrr.io/r/stats/pairwise.t.test.html),
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html),
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
and
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)
for the tidy methods

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
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# All pairwise comparisons among the three arms of the depression_bdi
# treatment study.
fit <- aov(bdi_post ~ condition, data = depression_bdi)
contrast_test(fit, contrasts = "pairwise")
#>  contrast            estimate se    t_value df p_value p_adjusted ci_lower
#>  placebo - ssri      4.900    2.815 1.741   27 0.0931  0.0931     -0.876  
#>  wait_list - ssri    6.700    2.815 2.380   27 0.0246  0.0246     0.924   
#>  wait_list - placebo 1.800    2.815 0.639   27 0.5279  0.5279     -3.976  
#>  ci_upper
#>  10.676  
#>  12.476  
#>  7.576   
#> 
#> Confidence level: 95%

# Custom contrasts with a Tukey-protected family-wise error rate.
contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
#>  contrast            estimate se    t_value df p_value p_adjusted ci_lower
#>  placebo - ssri      4.900    2.815 1.741   27 0.0931  0.2088     -2.079  
#>  wait_list - ssri    6.700    2.815 2.380   27 0.0246  0.0617     -0.279  
#>  wait_list - placebo 1.800    2.815 0.639   27 0.5279  0.7998     -5.179  
#>  ci_upper
#>  11.879  
#>  13.679  
#>  8.779   
#> 
#> Confidence level: 95%, simultaneous interval (Tukey)

# A user-defined contrast: the SSRI arm vs. the average of the placebo
# and wait list arms. With levels ordered ssri, placebo, wait_list, the
# weights c(1, -0.5, -0.5) estimate that difference.
contrast_test(
  fit,
  contrasts = list("ssri vs non-drug arms" = c(1, -0.5, -0.5)),
  adjust = "scheffe"
)
#>  contrast              estimate se    t_value df p_value p_adjusted ci_lower
#>  ssri vs non-drug arms -5.800   2.438 -2.379  27 0.0247  0.0766     -12.114 
#>  ci_upper
#>  0.514   
#> 
#> Confidence level: 95%, simultaneous interval (Scheffe)

# Welch-style inference: the wait list variance is about twice the
# SSRI variance, so the pooled error term is worth questioning.
contrast_test(fit, contrasts = "pairwise", var_equal = FALSE)
#>  contrast            estimate se    t_value df     p_value p_adjusted ci_lower
#>  placebo - ssri      4.900    2.480 1.976   17.918 0.0638  0.0638     -0.313  
#>  wait_list - ssri    6.700    2.933 2.285   16.202 0.0361  0.0361     0.489   
#>  wait_list - placebo 1.800    3.003 0.599   16.764 0.5569  0.5569     -4.543  
#>  ci_upper
#>  10.113  
#>  12.911  
#>  8.143   
#> 
#> Confidence level: 95%

# Pairwise treatment comparisons in the Smith, Meyers, and Delaney
# (1998) drinking trial, on the normalizing log scale. Each row is
# one pairwise contrast of the three treatment means.
fit_drinks <- aov(log_drinks ~ treatment, data = drinks_trial)
contrast_test(fit_drinks, contrasts = "pairwise")
#>  contrast                    estimate se    t_value df p_value p_adjusted
#>  CRA - Standard              -0.454   0.198 -2.294  85 0.0242  0.0242    
#>  CRA + Disulfiram - Standard -0.594   0.232 -2.567  85 0.0120  0.0120    
#>  CRA + Disulfiram - CRA      -0.140   0.238 -0.589  85 0.5575  0.5575    
#>  ci_lower ci_upper
#>  -0.848   -0.061  
#>  -1.055   -0.134  
#>  -0.612   0.332   
#> 
#> Confidence level: 95%

# An a priori contrast: the two active CRA arms (averaged) versus
# standard care. With levels ordered Standard, CRA, CRA + Disulfiram,
# the weights c(-1, 0.5, 0.5) compare the active arms against Standard.
contrast_test(
  fit_drinks,
  contrasts = list("CRA arms vs Standard" = c(-1, 0.5, 0.5))
)
#>  contrast             estimate se    t_value df p_value p_adjusted ci_lower
#>  CRA arms vs Standard -0.524   0.180 -2.917  85 0.0045  0.0045     -0.882  
#>  ci_upper
#>  -0.167  
#> 
#> Confidence level: 95%
```
