# Scheffe-Adjusted Simultaneous Confidence Intervals for Contrasts

Computes the Scheffe (1953, 1959) simultaneous confidence intervals on
user-specified contrasts among the means of a one-way design. The
Scheffe procedure controls the family-wise error rate for *any* set of
contrasts, however many and however post-hoc, which makes it more
conservative than Tukey-Kramer or Bonferroni for the specific case of
all-pairwise comparisons but optimal for arbitrary post-hoc contrasts.

## Usage

``` r
scheffe_ci(x, group = NULL, contrasts = NULL, conf_level = 0.95)
```

## Arguments

- x:

  A fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object with a single one-way
  factor predictor, or a numeric vector of observations with `group`
  supplied.

- group:

  Optional factor of group labels when `x` is a numeric vector.

- contrasts:

  An \\a \times m\\ matrix or vector of contrast coefficients (rows =
  levels, columns = contrasts). Each column must sum to zero. If `NULL`
  (default), the function returns intervals for all \\a - 1\\ pairwise
  contrasts.

- conf_level:

  Family-wise confidence level. Default `0.95`.

## Value

A `data.frame` with one row per contrast. Columns: `contrast` (a printed
label), `contrast_value`, `se`, `F_statistic`, `lower_limit`,
`upper_limit`, `p_adjusted`.

## Details

**Critical value.** For \\a\\ groups with \\\nu\\ error degrees of
freedom, the Scheffe critical value is \$\$S \\=\\ \sqrt{(a - 1) F\_{1 -
\alpha, a - 1, \nu}},\$\$ where \\F\_{1 - \alpha, a - 1, \nu}\\ is the
upper \\\alpha\\ quantile of the central *F* distribution. The Scheffe
simultaneous CI on a contrast \\\psi = \sum_i c_i \mu_i\\ is
\$\$\hat\psi \\\pm\\ S \cdot \mathit{SE}(\hat\psi),\$\$ where
\\\mathit{SE}(\hat\psi) = \sqrt{\mathit{MS}\_E \sum_i c_i^2 / n_i}\\.

**Scope.** The Scheffe family-wise coverage holds for *any* number of
contrasts, pairwise, complex, or chosen after looking at the data. The
trade-off is conservativeness: for all-pairwise comparisons,
Tukey-Kramer is uniformly more powerful.

## References

Scheffe, H. (1953). A method for judging all contrasts in the analysis
of variance. *Biometrika, 40*(1/2), 87–104.

Scheffe, H. (1959). *The analysis of variance*. Wiley.

## See also

[`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`tukey_kramer_ci`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
[`dunnett_ci`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md)

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
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
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
# 1. All pairwise contrasts via the default:
fit <- lm(weight ~ group, data = PlantGrowth)
scheffe_ci(fit)
#>  contrast    contrast_value se    F_statistic lower_limit upper_limit
#>  trt1 - ctrl -0.371         0.279 0.886       -1.09       0.351      
#>  trt2 - ctrl 0.494          0.279 1.57        -0.228      1.22       
#>  trt2 - trt1 0.865          0.279 4.81        0.143       1.59       
#>  p_adjusted
#>  0.4241    
#>  0.2265    
#>  0.0163    
#> 
#> Confidence level: 95%

# 2. A complex contrast: average of treatments vs control.
#        Coefficients (ctrl, trt1, trt2) = (-1, 0.5, 0.5).
cmat <- matrix(c(-1, 0.5, 0.5), nrow = 3,
               dimnames = list(c("ctrl", "trt1", "trt2"), "avg(trt) - ctrl"))
scheffe_ci(fit, contrasts = cmat)
#>  contrast        contrast_value se    F_statistic lower_limit upper_limit
#>  avg(trt) - ctrl 0.0615         0.241 0.0324      -0.564      0.687      
#>  p_adjusted
#>  0.9681    
#> 
#> Confidence level: 95%
```
