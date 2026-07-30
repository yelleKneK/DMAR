# Tukey-Kramer Simultaneous Confidence Intervals for Pairwise Contrasts

Computes the Tukey-Kramer simultaneous confidence intervals for all \\a
(a - 1) / 2\\ pairwise contrasts among \\a\\ group means in a one-way
design with possibly unequal group sample sizes (Tukey, 1953; Kramer,
1956; Hayter, 1984), and returns the result in tidy long form. Each
interval has individual coverage at least the specified `conf_level` and
the family-wise coverage is at least `conf_level`.

## Usage

``` r
tukey_kramer_ci(x, group = NULL, conf_level = 0.95)
```

## Arguments

- x:

  Either (a) a fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object with a single one-way
  factor predictor, or (b) a numeric vector of observations, in which
  case `group` must also be supplied.

- group:

  When `x` is a vector, a factor (or coercible to factor) of group
  labels, same length as `x`.

- conf_level:

  Family-wise confidence level. Default `0.95`.

## Value

A `data.frame` with one row per pairwise contrast. Columns: `contrast`
(e.g., `"B - A"`), `mean_difference`, `se`, `q_statistic` (the
studentized-range \\q\\), `lower_limit`, `upper_limit`, `p_adjusted`.

## Details

**Formula.** For groups \\i, j\\ with means \\\bar y_i, \bar y_j\\ and
sample sizes \\n_i, n_j\\, the Tukey-Kramer simultaneous CI is \$\$\bar
y_i - \bar y_j \\\pm\\ q\_{\alpha, a, \nu}
\sqrt{\frac{\mathit{MS}\_E}{2} \left(\frac{1}{n_i} +
\frac{1}{n_j}\right)},\$\$ where \\q\_{\alpha, a, \nu}\\ is the upper
\\\alpha\\ quantile of the studentized-range distribution with \\a\\
groups and \\\nu\\ error degrees of freedom
([`stats::qtukey()`](https://rdrr.io/r/stats/Tukey.html)).

**Why Tukey-Kramer.** Hayter (1984) proved that the Tukey-Kramer
procedure is conservative for unbalanced designs (the coverage
probability is at least `conf_level`). For balanced designs it reduces
to Tukey's HSD and the coverage is exactly `conf_level`.

**Adjusted *p*-values.** Each pairwise \\p\\-value is computed from the
studentized-range distribution: \\p = 1 - \mathrm{ptukey}(\|q\|, a,
\nu)\\.

## References

Hayter, A. J. (1984). A proof of the conjecture that the Tukey-Kramer
multiple comparisons procedure is conservative. *Annals of Statistics,
12*(1), 61–75.

Kramer, C. Y. (1956). Extension of multiple range tests to group means
with unequal numbers of replications. *Biometrics, 12*(3), 307–310.

Tukey, J. W. (1953). *The problem of multiple comparisons*. Unpublished
manuscript, Princeton University.

## See also

[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`dunnett_ci`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md),
[`scheffe_ci`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md),
[`TukeyHSD`](https://rdrr.io/r/stats/TukeyHSD.html)

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
[`scheffe_ci()`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Balanced one-way: built-in PlantGrowth dataset.
fit <- lm(weight ~ group, data = PlantGrowth)
tukey_kramer_ci(fit)
#>  contrast    mean_difference se    q_statistic lower_limit upper_limit
#>  trt1 - ctrl -0.371          0.197 -1.88       -1.06       0.32       
#>  trt2 - ctrl 0.494           0.197 2.51        -0.197      1.19       
#>  trt2 - trt1 0.865           0.197 4.39        0.174       1.56       
#>  p_adjusted
#>  0.3909    
#>  0.1980    
#>  0.0120    
#> 
#> Confidence level: 95%

# 2. Same data via vector / group interface:
tukey_kramer_ci(PlantGrowth$weight, group = PlantGrowth$group)
#>  contrast    mean_difference se    q_statistic lower_limit upper_limit
#>  trt1 - ctrl -0.371          0.197 -1.88       -1.06       0.32       
#>  trt2 - ctrl 0.494           0.197 2.51        -0.197      1.19       
#>  trt2 - trt1 0.865           0.197 4.39        0.174       1.56       
#>  p_adjusted
#>  0.3909    
#>  0.1980    
#>  0.0120    
#> 
#> Confidence level: 95%
```
