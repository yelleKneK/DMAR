# Scheffe-Adjusted Simultaneous Confidence Intervals for Contrasts

Computes the Scheffe (1953, 1959) simultaneous confidence intervals on
user-specified contrasts among the means of a one-way design. The
Scheffe procedure controls the family-wise error rate for *any* set of
contrasts, however many and however post-hoc, which makes it more
conservative than Tukey-Kramer or Bonferroni for the specific case of
all-pairwise comparisons but optimal for arbitrary post-hoc contrasts.

## Usage

``` r
ci_scheffe(x, group = NULL, contrasts = NULL, conf_level = 0.95)
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
  (default), the function returns intervals for all \\a (a - 1) / 2\\
  pairwise contrasts.

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
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
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
# 1. All pairwise contrasts among the six marketing panels of the
#    test_market data via the default:
fit <- lm(brand_movement ~ panel, data = test_market)
ci_scheffe(fit)
#>  contrast contrast_value se   F_statistic lower_limit upper_limit p_adjusted
#>  2 - 1    0.11           0.38 0.0168      -1.3        1.52        0.9999    
#>  3 - 1    0.685          0.38 0.65        -0.73       2.1         0.6652    
#>  4 - 1    0.87           0.38 1.05        -0.545      2.28        0.4201    
#>  5 - 1    1.25           0.38 2.15        -0.17       2.66        0.1061    
#>  6 - 1    1.3            0.38 2.32        -0.12       2.71        0.0855    
#>  3 - 2    0.575          0.38 0.458       -0.84       1.99        0.8021    
#>  4 - 2    0.76           0.38 0.8         -0.655      2.17        0.5639    
#>  5 - 2    1.14           0.38 1.78        -0.28       2.55        0.1668    
#>  6 - 2    1.19           0.38 1.95        -0.23       2.6         0.1364    
#>  4 - 3    0.185          0.38 0.0474      -1.23       1.6         0.9984    
#>  5 - 3    0.56           0.38 0.434       -0.855      1.97        0.8186    
#>  6 - 3    0.61           0.38 0.516       -0.805      2.02        0.7611    
#>  5 - 4    0.375          0.38 0.195       -1.04       1.79        0.9605    
#>  6 - 4    0.425          0.38 0.25        -0.99       1.84        0.9342    
#>  6 - 5    0.05           0.38 0.00346     -1.36       1.46        1.0000    
#> 
#> Confidence level: 95%

# 2. A contrast chosen after inspecting the means: the two panels with
#    the highest brand movement (5 and 6) against the two with the
#    lowest (1 and 2). The Scheffe coverage holds for a contrast picked
#    this way, and the interval still excludes zero even though none of
#    the pairwise intervals above does.
cmat <- matrix(c(-0.5, -0.5, 0, 0, 0.5, 0.5), nrow = 6,
               dimnames = list(levels(test_market$panel),
                               "panels 5,6 - panels 1,2"))
ci_scheffe(fit, contrasts = cmat)
#>  contrast                contrast_value se    F_statistic lower_limit
#>  panels 5,6 - panels 1,2 1.22           0.269 4.09        0.215      
#>  upper_limit p_adjusted
#>  2.22        0.0117    
#> 
#> Confidence level: 95%
```
