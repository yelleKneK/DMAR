# Dunnett's Simultaneous Confidence Intervals Against a Control

Computes Dunnett's (1955, 1964) simultaneous confidence intervals for
the \\a - 1\\ comparisons of \\a - 1\\ treatment means against a single
control mean, with family-wise coverage at the specified `conf_level`.
Returns the result in tidy long form.

## Usage

``` r
ci_dunnett(
  x,
  group = NULL,
  control = NULL,
  alternative = c("two_sided", "less", "greater"),
  conf_level = 0.95
)
```

## Arguments

- x:

  Either (a) a fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object with a one-way factor
  predictor, or (b) a numeric vector of observations, in which case
  `group` must also be supplied.

- group:

  When `x` is a vector, a factor of group labels of the same length.

- control:

  Character name of the control level (must be one of the factor
  levels). If `NULL` (default), the first level is used.

- alternative:

  One of `"two_sided"` (default; the base-R spelling `"two.sided"` is
  accepted as an alias), `"less"`, or `"greater"`.

- conf_level:

  Family-wise confidence level. Default `0.95`.

## Value

A `data.frame` with one row per non-control level. Columns: `contrast`,
`mean_difference`, `se`, `t_statistic`, `lower_limit`, `upper_limit`,
`p_adjusted`.

## Details

**Critical value.** The two-sided Dunnett critical value \\d\_{\alpha,
a - 1, \nu}^{(2)}\\ is obtained from the multivariate *t* distribution
with \\a - 1\\ dimensions, common correlation \\0.5\\ (the Dunnett
correlation under balanced *n*; the function does not adjust for unequal
*n*), and \\\nu\\ error degrees of freedom. The function uses the
existing
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
critical value.

**Adjusted *p*-values.** Computed exactly from the same equicorrelated
multivariate *t* distribution. The one common correlation \\1/2\\ admits
a one-factor representation, so the probability that all comparisons
fall inside (or below) the observed statistic collapses to two nested
one-dimensional integrals, evaluated by quadrature. The adjusted
*p*-value is one minus that probability. The computation is
deterministic (no Monte Carlo) and needs no additional package.

## References

Dunnett, C. W. (1955). A multiple comparison procedure for comparing
several treatments with a control. *Journal of the American Statistical
Association, 50*(272), 1096–1121.

Dunnett, C. W. (1964). New tables for multiple comparisons with a
control. *Biometrics, 20*(3), 482–491.

Hsu, J. C. (1996). *Multiple comparisons: Theory and methods*. Chapman &
Hall.

## See also

[`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`ci_scheffe`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
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
# 1. Compare the SSRI and placebo arms of the depression_bdi study
#    against the wait list control:
fit <- lm(bdi_post ~ condition, data = depression_bdi)
ci_dunnett(fit, control = "wait_list")
#>  contrast            mean_difference se   t_statistic lower_limit upper_limit
#>  ssri - wait_list    -6.7            2.81 -2.38       -13.3       -0.132     
#>  placebo - wait_list -1.8            2.81 -0.639      -8.37       4.77       
#>  p_adjusted
#>  0.0452    
#>  0.7494    
#> 
#> Confidence level: 95%

# 2. One-sided: a treatment that works pulls the posttest BDI down,
#    so the directional alternative is "less":
ci_dunnett(fit, control = "wait_list", alternative = "less")
#>  contrast            mean_difference se   t_statistic lower_limit upper_limit
#>  ssri - wait_list    -6.7            2.81 -2.38       -Inf        -1.08      
#>  placebo - wait_list -1.8            2.81 -0.639      -Inf        3.82       
#>  p_adjusted
#>  0.0226    
#>  0.3973    
#> 
#> Confidence level: 95%
```
