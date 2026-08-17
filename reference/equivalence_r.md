# Equivalence Test for the Pearson Correlation via Two One-Sided Tests (TOST)

Performs a two one-sided tests procedure for equivalence of a Pearson
correlation \\\rho\\ to zero against user-specified equivalence bounds
\\\[-\rho_L, \rho_U\]\\ (Counsell & Cribbie, 2015; Goertzen & Cribbie,
2010). Uses the Fisher's \\Z\\ transformation throughout, a large-sample
approximation whose accuracy under bivariate normality improves quickly
with *n*. Equivalence is declared when the 100(1 - 2\\\alpha\\)%
Fisher's \\Z\\ CI on \\\rho\\ lies entirely inside the equivalence
region.

## Usage

``` r
equivalence_r(
  r = NULL,
  n = NULL,
  x = NULL,
  y = NULL,
  rho_lower = NULL,
  rho_upper = NULL,
  alpha_level = 0.05
)
```

## Arguments

- r, n:

  Observed sample correlation *r* and sample size. Alternatively supply
  `x` and `y` to compute *r* from raw data.

- x, y:

  Numeric vectors of paired observations. If supplied, `r` and `n` are
  computed from the data and the `r`/`n` arguments are ignored.

- rho_lower, rho_upper:

  Equivalence bounds on the correlation scale, both positive. The
  equivalence region is \\\[-\rho_L, +\rho_U\]\\. If only `rho_upper` is
  supplied, the bounds are symmetric.

- alpha_level:

  One-sided significance level. Default `0.05`.

## Value

A `data.frame` with rows for the observed *r*, the two one-sided test
statistics on the Fisher's \\Z\\ scale, their *p*-values, the joint TOST
*p*-value, the 100(1 - 2\\\alpha\\)% CI on \\\rho\\, the equivalence
bounds, a binary equivalence flag, and the sample size (`n`).

## Details

**Fisher's \\Z\\ transformation.** \$\$Z = \tfrac{1}{2}
\log\left(\frac{1 + r}{1 - r}\right), \quad \mathrm{Var}(Z) =
\frac{1}{n - 3}.\$\$ The TOST is run on the Fisher's \\Z\\ scale:

- Lower test: \\(Z - Z\_{-\rho_L}) \sqrt{n - 3}\\ compared against the
  upper \\\alpha\\ of \\N(0, 1)\\.

- Upper test: \\(Z - Z\_{\rho_U}) \sqrt{n - 3}\\ compared against the
  lower \\\alpha\\ of \\N(0, 1)\\.

The CI bounds are back-transformed from the \\Z\\ scale via \\r =
\tanh(Z)\\ so they remain in \\\[-1, 1\]\\.

**Choosing \\\rho_L\\ and \\\rho_U\\.** Common choices in psychology are
\\0.1\\, or domain-specific meaningfulness thresholds (e.g., \\0.2\\ for
cognitive task correlations). The bounds must be set *before* data
collection.

## References

Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for comparing
correlation and regression coefficients. *British Journal of
Mathematical and Statistical Psychology, 68*(2), 292–309.
[doi:10.1111/bmsp.12045](https://doi.org/10.1111/bmsp.12045)

Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of
association: An equivalence testing approach. *British Journal of
Mathematical and Statistical Psychology, 63*(3), 527–537.
[doi:10.1348/000711009X475853](https://doi.org/10.1348/000711009X475853)

Lakens, D. (2017). Equivalence tests: A practical primer for *t* tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*(4), 355–362.
[doi:10.1177/1948550617697177](https://doi.org/10.1177/1948550617697177)

## See also

[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
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

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Equivalence test that |rho| < 0.10 with n = 200 and r = 0.05:
equivalence_r(r = 0.05, n = 200, rho_upper = 0.10)
#>  term         value 
#>  r            0.05  
#>  z_lower_test 2.11  
#>  z_upper_test -0.706
#>  p_lower      0.0174
#>  p_upper      0.2401
#>  p_tost       0.2401
#>  lower_limit  -0.067
#>  upper_limit  0.166 
#>  rho_lower    -0.1  
#>  rho_upper    0.1   
#>  equivalent   0     
#>  n            200   
#> 
#> Confidence level: 90%

# 2. From raw data:
set.seed(113)
x <- rnorm(150); y <- 0.04 * x + rnorm(150)
equivalence_r(x = x, y = y, rho_upper = 0.15)
#>  term         value  
#>  r            -0.0329
#>  z_lower_test 1.43   
#>  z_upper_test -2.23  
#>  p_lower      0.0759 
#>  p_upper      0.0128 
#>  p_tost       0.0759 
#>  lower_limit  -0.167 
#>  upper_limit  0.102  
#>  rho_lower    -0.15  
#>  rho_upper    0.15   
#>  equivalent   0      
#>  n            150    
#> 
#> Confidence level: 90%
```
