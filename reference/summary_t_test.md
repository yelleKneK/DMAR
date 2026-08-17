# Two-Sample *t* Test From Summary Statistics

Computes a two-sample *t* test (pooled or Welch) directly from the
per-group means, standard deviations, and sample sizes, without
requiring access to the raw observations. Returns the test statistic,
degrees of freedom, *p*-value, and a CI on the mean difference in a
`data.frame`. Useful for re-analyses from published papers that report
only the summary numbers.

## Usage

``` r
summary_t_test(
  mean_1,
  sd_1,
  n_1,
  mean_2,
  sd_2,
  n_2,
  mu = 0,
  var_equal = TRUE,
  alternative = c("two_sided", "less", "greater"),
  conf_level = 0.95
)
```

## Arguments

- mean_1, mean_2:

  Group sample means.

- sd_1, sd_2:

  Group sample standard deviations.

- n_1, n_2:

  Group sample sizes.

- mu:

  Null value of the mean difference \\\mu_1 - \mu_2\\. Default `0`.

- var_equal:

  Logical. If `TRUE` (default), uses Student's pooled-variance *t*. If
  `FALSE`, uses Welch's separate- variance *t* with Satterthwaite
  degrees of freedom.

- alternative:

  One of `"two_sided"` (default; the base-R spelling `"two.sided"` is
  accepted as an alias), `"less"`, or `"greater"`.

- conf_level:

  Confidence level for the CI on the mean difference. Default `0.95`.

## Value

A `data.frame` with rows for the mean difference, the *t* statistic,
degrees of freedom, *p*-value, and the CI lower and upper limits on the
mean difference.

## Details

**Pooled-variance *t* (Student, 1908).** Under \\\sigma_1 = \sigma_2\\,
the pooled SD is \\s_p = \sqrt{((n_1 - 1) s_1^2 + (n_2 - 1) s_2^2) /
(n_1 + n_2 - 2)}\\, the test statistic is \\t = (\bar x_1 - \bar x_2 -
\mu_0) / (s_p \sqrt{1 / n_1 + 1 / n_2})\\, and \\df = n_1 + n_2 - 2\\.

**Welch's *t* (Welch, 1947).** Under unequal variances, \\t = (\bar
x_1 - \bar x_2 - \mu_0) / \sqrt{s_1^2 / n_1 + s_2^2 / n_2}\\ with
Satterthwaite degrees of freedom (see
[`welch_t`](https://yelleknek.github.io/DMAR/reference/welch_t.md)).

**Choosing pooled vs Welch.** Methodological reviews now recommend Welch
as the default (Delacre, Lakens, & Leys, 2017; Ruxton, 2006).
Pooled-variance *t* is preserved here primarily for reproducing analyses
from older sources that used it.

## References

Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should by
default use Welch's *t*-test instead of Student's *t*-test.
*International Review of Social Psychology, 30*(1), 92–101.
[doi:10.5334/irsp.82](https://doi.org/10.5334/irsp.82)

Ruxton, G. D. (2006). The unequal variance *t*-test is an underused
alternative to Student's *t*-test and the Mann-Whitney *U* test.
*Behavioral Ecology, 17*(4), 688–690.
[doi:10.1093/beheco/ark016](https://doi.org/10.1093/beheco/ark016)

Snedecor, G. W., & Cochran, W. G. (1989). *Statistical methods* (8th
ed.). Iowa State University Press.

Student. (1908). The probable error of a mean. *Biometrika, 6*(1), 1–25.
[doi:10.2307/2331554](https://doi.org/10.2307/2331554)

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1/2), 28–35.

## See also

[`welch_t`](https://yelleknek.github.io/DMAR/reference/welch_t.md),
[`t.test`](https://rdrr.io/r/stats/t.test.html),
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md)

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
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Re-analysis from published summary statistics:
#        Group A: M = 100, SD = 15, n = 30
#        Group B: M = 108, SD = 18, n = 25
summary_t_test(mean_1 = 100, sd_1 = 15, n_1 = 30,
               mean_2 = 108, sd_2 = 18, n_2 = 25)
#>  term            value 
#>  mean_difference -8    
#>  t_statistic     -1.8  
#>  df              53    
#>  p_value         0.0778
#>  lower_limit     -16.9 
#>  upper_limit     0.922 
#> 
#> Confidence level: 95%

# 2. Welch version for the same data:
summary_t_test(mean_1 = 100, sd_1 = 15, n_1 = 30,
               mean_2 = 108, sd_2 = 18, n_2 = 25,
               var_equal = FALSE)
#>  term            value 
#>  mean_difference -8    
#>  t_statistic     -1.77 
#>  df              46.8  
#>  p_value         0.0835
#>  lower_limit     -17.1 
#>  upper_limit     1.1   
#> 
#> Confidence level: 95%
```
