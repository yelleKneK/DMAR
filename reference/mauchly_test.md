# Mauchly's Test of Sphericity for a One-Way Within-Subjects Design

Tests the null hypothesis that the covariance matrix of the orthonormal
contrasts among the \\k\\ repeated measurements is proportional to the
identity (the *sphericity* assumption underlying univariate repeated
measures *F*-tests).

## Usage

``` r
mauchly_test(x, id = NULL, time = NULL, outcome = NULL)
```

## Arguments

- x:

  Either an \\n \times k\\ numeric matrix or `data.frame` (rows =
  subjects, columns = repeated measurements); *or* a long-format
  `data.frame` together with `id`, `time`, and `outcome` column names.

- id:

  Column name in `x` identifying the subject when `x` is in long format
  (`NULL` otherwise).

- time:

  Column name in `x` identifying the within-subjects factor level when
  `x` is in long format (`NULL` otherwise).

- outcome:

  Column name in `x` identifying the dependent variable when `x` is in
  long format (`NULL` otherwise).

## Value

A one-row `data.frame` with columns `W` (Mauchly's statistic),
`statistic` (the chi square approximation), `df`, `p_value`,
`n_subjects`, `n_levels`, and `method`.

## Details

Sphericity is the assumption that the variances of all pairwise
differences among the \\k\\ levels are equal, equivalently, that the
covariance matrix \\\Sigma_C\\ of any orthonormal set of \\k - 1\\
contrasts among the levels is proportional to the identity. Mauchly's
(1940) test statistic is \$\$W =
\frac{\det(\hat\Sigma_C)}{\bigl(\mathrm{tr}(\hat\Sigma_C) / (k -
1)\bigr)^{k - 1}},\$\$ and the chi square approximation \$\$X^2 = -\\m
\\\log W \quad \mathrm{with}\\ m = (n - 1) - \frac{2(k - 1)^2 +
(k - 1) + 2}{6\\(k - 1)}\$\$ has approximately \\(k - 1)k/2 - 1\\
degrees of freedom under \\H_0\\. The reported *p*-value uses Box's
(1949) second-order correction, a weighted combination of the chi square
tails on \\(k - 1)k/2 - 1\\ and \\(k - 1)k/2 + 3\\ degrees of freedom,
which improves the first-order approximation in small samples; this
matches [`mauchly.test`](https://rdrr.io/r/stats/mauchly.test.html).

When sphericity is rejected, the univariate *F* test is liberal; correct
using the Greenhouse-Geisser, Huynh-Feldt, or lower-bound epsilon
adjustments via
[`epsilon_corrections`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md)
or directly via
[`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md).

The test is only defined for \\k \ge 3\\; with \\k = 2\\, sphericity is
trivially true and the function returns `W = 1, p = 1`.

## References

Mauchly, J. W. (1940). Significance test for sphericity of a normal
\\n\\-variate distribution. *Annals of Mathematical Statistics, 11*(2),
204–209.

Box, G. E. P. (1949). A general distribution theory for a class of
likelihood criteria. *Biometrika, 36*(3/4), 317–346.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 11 for sphericity in within-subjects
designs.)

## See also

[`epsilon_corrections`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
[`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md)

Other within-subjects analysis:
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
[`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)

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
# Wide-format example: simulated within-subjects data with 4 levels.
set.seed(113)
Y <- matrix(rnorm(20 * 4), nrow = 20)
mauchly_test(Y)
#>  W     statistic df p_value n_subjects n_levels method                      
#>  0.597 9.13      5  0.1045  20         4        Mauchly's test of sphericity

# Long-format example using built-in nlme::Orthodont (4 ages per subject).
mauchly_test(nlme::Orthodont, id = "Subject", time = "age",
             outcome = "distance")
#>  W     statistic df p_value n_subjects n_levels method                      
#>  0.758 6.85      5  0.2326  27         4        Mauchly's test of sphericity
```
