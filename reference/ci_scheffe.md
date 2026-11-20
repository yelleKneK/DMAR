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
label), `estimate`, `se`, `F_value`, `ci_lower`, `ci_upper`,
`p_adjusted`.

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
[`ci_dunnett_t3()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
[`ci_games_howell()`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
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
#>  contrast estimate se    F_value ci_lower ci_upper p_adjusted
#>  2 - 1    0.110    0.380 0.017   -1.305   1.525    0.9999    
#>  3 - 1    0.685    0.380 0.650   -0.730   2.100    0.6652    
#>  4 - 1    0.870    0.380 1.049   -0.545   2.285    0.4201    
#>  5 - 1    1.245    0.380 2.147   -0.170   2.660    0.1061    
#>  6 - 1    1.295    0.380 2.323   -0.120   2.710    0.0855    
#>  3 - 2    0.575    0.380 0.458   -0.840   1.990    0.8021    
#>  4 - 2    0.760    0.380 0.800   -0.655   2.175    0.5639    
#>  5 - 2    1.135    0.380 1.785   -0.280   2.550    0.1668    
#>  6 - 2    1.185    0.380 1.945   -0.230   2.600    0.1364    
#>  4 - 3    0.185    0.380 0.047   -1.230   1.600    0.9984    
#>  5 - 3    0.560    0.380 0.434   -0.855   1.975    0.8186    
#>  6 - 3    0.610    0.380 0.516   -0.805   2.025    0.7611    
#>  5 - 4    0.375    0.380 0.195   -1.040   1.790    0.9605    
#>  6 - 4    0.425    0.380 0.250   -0.990   1.840    0.9342    
#>  6 - 5    0.050    0.380 0.003   -1.365   1.465    1.0000    
#> 
#> Confidence level: 95%, simultaneous interval (Scheffe)

# 2. A contrast chosen after inspecting the means: the two panels with
#    the highest brand movement (5 and 6) against the two with the
#    lowest (1 and 2). The Scheffe coverage holds for a contrast picked
#    this way, and the interval still excludes zero even though none of
#    the pairwise intervals above does.
cmat <- matrix(c(-0.5, -0.5, 0, 0, 0.5, 0.5), nrow = 6,
               dimnames = list(levels(test_market$panel),
                               "panels 5,6 - panels 1,2"))
ci_scheffe(fit, contrasts = cmat)
#>  contrast                estimate se    F_value ci_lower ci_upper p_adjusted
#>  panels 5,6 - panels 1,2 1.215    0.269 4.090   0.215    2.215    0.0117    
#> 
#> Confidence level: 95%, simultaneous interval (Scheffe)
```
