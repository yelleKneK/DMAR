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
`estimate`, `se`, `t`, `df`, `p`, `p_adj`, `conf_lower`, and
`conf_upper`. The adjustment, confidence level, and variance assumption
are stored as `attr(*, "adjust")`, `attr(*, "conf_level")`, and
`attr(*, "var_equal")`. The table prints through the
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

**Adjustments.** The `p_adj` and confidence interval critical value are
computed as follows.

- `"none"`: no adjustment; the CI uses \\t\_{1-\alpha/2,df}\\.

- `"bonferroni"`: \\p\_{\text{adj}} = \min(1, m\\ p)\\ for \\m\\
  contrasts, with CI based on \\t\_{1-\alpha/(2m),df}\\.

- `"scheffe"`: appropriate for any contrast (or family of contrasts).
  \\p\_{\text{adj}}\\ comes from the upper tail of an *F* reference
  distribution applied to \\t^2 / (k-1)\\, and the CI uses
  \\\sqrt{(k-1)\\ F\_{1-\alpha,\\k-1,df}}\\.

- `"tukey"`: requires every contrast to be pairwise. Uses the
  studentized range distribution (`ptukey`/`qtukey`) so that
  \\p\_{\text{adj}} = 1 - \mathrm{ptukey}(\|t\|\sqrt{2}; k, df)\\ and
  the CI uses \\q\_{1-\alpha,\\k,df} / \sqrt{2}\\.

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
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
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
# All pairwise comparisons among the three PlantGrowth groups.
fit <- aov(weight ~ group, data = PlantGrowth)
contrast_test(fit, contrasts = "pairwise")
#>  contrast    estimate se    t     df p       p_adj   conf_lower conf_upper
#>  trt1 - ctrl -0.371   0.279 -1.33 27 0.194   0.194   -0.943     0.201     
#>  trt2 - ctrl 0.494    0.279 1.77  27 0.0877  0.0877  -0.078     1.07      
#>  trt2 - trt1 0.865    0.279 3.1   27 0.00446 0.00446 0.293      1.44      
#> 
#> Confidence level: 95%

# Custom contrasts with a Tukey-protected family-wise error rate.
contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
#>  contrast    estimate se    t     df p       p_adj conf_lower conf_upper
#>  trt1 - ctrl -0.371   0.279 -1.33 27 0.194   0.391 -1.06      0.32      
#>  trt2 - ctrl 0.494    0.279 1.77  27 0.0877  0.198 -0.197     1.19      
#>  trt2 - trt1 0.865    0.279 3.1   27 0.00446 0.012 0.174      1.56      
#> 
#> Confidence level: 95%

# A user-defined contrast: ctrl vs. average of the two treatments.
contrast_test(
  fit,
  contrasts = list("ctrl vs trts" = c(1, -0.5, -0.5)),
  adjust = "scheffe"
)
#>  contrast     estimate se    t      df p     p_adj conf_lower conf_upper
#>  ctrl vs trts -0.0615  0.241 -0.255 27 0.801 0.968 -0.687     0.564     
#> 
#> Confidence level: 95%

# Welch-style inference (unequal group variances).
contrast_test(fit, contrasts = "pairwise", var_equal = FALSE)
#>  contrast    estimate se    t     df   p      p_adj  conf_lower conf_upper
#>  trt1 - ctrl -0.371   0.311 -1.19 16.5 0.25   0.25   -1.03      0.288     
#>  trt2 - ctrl 0.494    0.231 2.13  16.8 0.0479 0.0479 0.00513    0.983     
#>  trt2 - trt1 0.865    0.287 3.01  14.1 0.0093 0.0093 0.249      1.48      
#> 
#> Confidence level: 95%

# Pairwise treatment comparisons in the Smith, Meyers, and Delaney
# (1998) drinking trial, on the normalizing log scale. Each row is
# one pairwise contrast of the three treatment means.
fit_drinks <- aov(log_drinks ~ treatment, data = drinks_trial)
contrast_test(fit_drinks, contrasts = "pairwise")
#>  contrast                    estimate se    t      df p      p_adj  conf_lower
#>  CRA - Standard              -0.454   0.198 -2.29  85 0.0242 0.0242 -0.848    
#>  CRA + Disulfiram - Standard -0.594   0.232 -2.57  85 0.012  0.012  -1.05     
#>  CRA + Disulfiram - CRA      -0.14    0.238 -0.589 85 0.557  0.557  -0.612    
#>  conf_upper
#>  -0.0606   
#>  -0.134    
#>  0.332     
#> 
#> Confidence level: 95%

# An a priori contrast: the two active CRA arms (averaged) versus
# standard care. With levels ordered Standard, CRA, CRA + Disulfiram,
# the weights c(-1, 0.5, 0.5) compare the active arms against Standard.
contrast_test(
  fit_drinks,
  contrasts = list("CRA arms vs Standard" = c(-1, 0.5, 0.5))
)
#>  contrast             estimate se   t     df p       p_adj   conf_lower
#>  CRA arms vs Standard -0.524   0.18 -2.92 85 0.00451 0.00451 -0.882    
#>  conf_upper
#>  -0.167    
#> 
#> Confidence level: 95%
```
