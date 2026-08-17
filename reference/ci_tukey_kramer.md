# Tukey-Kramer Simultaneous Confidence Intervals for Pairwise Contrasts

Computes the Tukey-Kramer simultaneous confidence intervals for all \\a
(a - 1) / 2\\ pairwise contrasts among \\a\\ group means in a one-way
design with possibly unequal group sample sizes (Tukey, 1953; Kramer,
1956; Hayter, 1984), and returns the result in tidy long form. Each
interval has individual coverage at least the specified `conf_level` and
the family-wise coverage is at least `conf_level`.

## Usage

``` r
ci_tukey_kramer(x, group = NULL, conf_level = 0.95)
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
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_scheffe`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`TukeyHSD`](https://rdrr.io/r/stats/TukeyHSD.html)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
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
# 1. Balanced one-way: the six marketing panels of the test_market
#    data, four outlets per panel, so the procedure is exactly Tukey's
#    HSD. Panels 5 and 6 separate from panel 1.
fit <- lm(brand_movement ~ panel, data = test_market)
ci_tukey_kramer(fit)
#>  contrast mean_difference se    q_statistic lower_limit upper_limit p_adjusted
#>  2 - 1    0.11            0.269 0.409       -1.1        1.32        0.9997    
#>  3 - 1    0.685           0.269 2.55        -0.522      1.89        0.4882    
#>  4 - 1    0.87            0.269 3.24        -0.337      2.08        0.2481    
#>  5 - 1    1.25            0.269 4.63        0.0375      2.45        0.0411    
#>  6 - 1    1.3             0.269 4.82        0.0875      2.5         0.0315    
#>  3 - 2    0.575           0.269 2.14        -0.632      1.78        0.6607    
#>  4 - 2    0.76            0.269 2.83        -0.447      1.97        0.3797    
#>  5 - 2    1.14            0.269 4.22        -0.0725     2.34        0.0725    
#>  6 - 2    1.19            0.269 4.41        -0.0225     2.39        0.0562    
#>  4 - 3    0.185           0.269 0.689       -1.02       1.39        0.9961    
#>  5 - 3    0.56            0.269 2.08        -0.647      1.77        0.6840    
#>  6 - 3    0.61            0.269 2.27        -0.597      1.82        0.6055    
#>  5 - 4    0.375           0.269 1.4         -0.832      1.58        0.9162    
#>  6 - 4    0.425           0.269 1.58        -0.782      1.63        0.8674    
#>  6 - 5    0.05            0.269 0.186       -1.16       1.26        1.0000    
#> 
#> Confidence level: 95%

# 2. Same data via vector / group interface:
ci_tukey_kramer(test_market$brand_movement, group = test_market$panel)
#>  contrast mean_difference se    q_statistic lower_limit upper_limit p_adjusted
#>  2 - 1    0.11            0.269 0.409       -1.1        1.32        0.9997    
#>  3 - 1    0.685           0.269 2.55        -0.522      1.89        0.4882    
#>  4 - 1    0.87            0.269 3.24        -0.337      2.08        0.2481    
#>  5 - 1    1.25            0.269 4.63        0.0375      2.45        0.0411    
#>  6 - 1    1.3             0.269 4.82        0.0875      2.5         0.0315    
#>  3 - 2    0.575           0.269 2.14        -0.632      1.78        0.6607    
#>  4 - 2    0.76            0.269 2.83        -0.447      1.97        0.3797    
#>  5 - 2    1.14            0.269 4.22        -0.0725     2.34        0.0725    
#>  6 - 2    1.19            0.269 4.41        -0.0225     2.39        0.0562    
#>  4 - 3    0.185           0.269 0.689       -1.02       1.39        0.9961    
#>  5 - 3    0.56            0.269 2.08        -0.647      1.77        0.6840    
#>  6 - 3    0.61            0.269 2.27        -0.597      1.82        0.6055    
#>  5 - 4    0.375           0.269 1.4         -0.832      1.58        0.9162    
#>  6 - 4    0.425           0.269 1.58        -0.782      1.63        0.8674    
#>  6 - 5    0.05            0.269 0.186       -1.16       1.26        1.0000    
#> 
#> Confidence level: 95%
```
