# Mixed-Design Multivariate ANOVA With All Four Test Statistics

Computes the multivariate analysis of variance for a mixed-design (one
between-subjects factor and one within-subjects factor), returning
Wilks's \\\Lambda\\, Pillai's trace, Hotelling-Lawley trace, and Roy's
largest root, each with the associated *F*-approximation, degrees of
freedom, and *p*-value. The three effects, between-subjects (*A*),
within-subjects (*B*), and the interaction (\\A \times B\\), are tested
separately. Wraps [`Anova`](https://rdrr.io/pkg/car/man/Anova.html) and
returns the result in the tidy DMAR style.

## Usage

``` r
manova_split_plot(data, within, between, ss_type = 3L)
```

## Arguments

- data:

  A `data.frame` in wide format, one row per subject, with the
  within-subjects measurements in separate columns and the
  between-subjects factor as one additional column.

- within:

  Character vector of column names holding the repeated measures values
  (one column per level of the within- subjects factor). Must be in the
  canonical level order.

- between:

  Character name of the between-subjects factor column in `data`.

- ss_type:

  The sum-of-squares type for the between-subjects effects, passed
  through to [`Anova`](https://rdrr.io/pkg/car/man/Anova.html). Accepts
  the integers `1`, `2`, or `3` or the equivalent Roman- numeral strings
  `"I"`, `"II"`, or `"III"`, and defaults to `3L` (Type III). The chosen
  type is recorded in the returned table (see **Value**). Type I (`1` or
  `"I"`) is not available for this mixed design, because
  [`Anova`](https://rdrr.io/pkg/car/man/Anova.html) computes only Type
  II and Type III sums of squares for the multivariate repeated measures
  path; requesting it raises an error.

## Value

A `data.frame` with rows for each of the three effects crossed with each
of the four multivariate statistics, plus one trailing row recording the
sum-of-squares type. Columns: `effect`, `statistic_name`,
`statistic_value`, `F_approx`, `df_1`, `df_2`, `p_value`. The final row
has `effect == "sum_of_squares_type"` and carries the chosen type (`1`,
`2`, or `3`) in its numeric `statistic_value`; its remaining numeric
columns are `NA`, so the `statistic_value` column stays numeric.

## Details

**The four statistics.** For an effect with \\H\\ and \\E\\ hypothesis-
and error-cross-products matrices:

- Wilks's \\\Lambda = \det(E) / \det(E + H)\\

- Pillai's trace \\V = \mathrm{tr}(H (E + H)^{-1})\\

- Hotelling-Lawley trace \\T_0^2 = \mathrm{tr}(H E^{-1})\\

- Roy's largest root \\\theta = \lambda_1(H E^{-1})\\

**When to use which.** Pillai's trace is the most robust to departures
from the multivariate normal / homogeneous-covariance assumptions.
Wilks's \\\Lambda\\ is the most widely reported. Roy's largest root is
the most powerful when the alternative concentrates on a single
dimension. The four statistics agree exactly when the effect has 1
numerator degree of freedom.

**Sum-of-squares type.** The between-subjects effects are computed by
[`Anova`](https://rdrr.io/pkg/car/man/Anova.html) using the
sum-of-squares type selected through `ss_type` (Type III by default).
Type II conditions each effect on the others that do not contain it, and
Type III conditions each effect on every other effect in the model; for
a single between-subjects factor the two coincide, and they can differ
once additional between-subjects terms are present. Type I, the
sequential decomposition, is not available here:
[`Anova`](https://rdrr.io/pkg/car/man/Anova.html) computes only Type II
and Type III for the multivariate repeated measures path. The type in
force is reported in the returned table so the analysis is
self-documenting.

**Dependency.** Requires the car package on CRAN.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 7 on higher-order designs and Chapter 14.)

Rencher, A. C., & Christensen, W. F. (2012). *Methods of multivariate
analysis* (4th ed.). Wiley.

## See also

[`Anova`](https://rdrr.io/pkg/car/man/Anova.html),
[`manova`](https://rdrr.io/r/stats/manova.html),
[`anova_within_two_way`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md)

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

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two groups of ten measured at three times. Both groups start at the
# same place and rise over time, and group B rises twice as fast, so the
# data generating means differ in slope as well as in level.
set.seed(113)
n_per <- 10
d <- data.frame(
  subject = factor(1:(2 * n_per)),
  group   = factor(rep(c("A", "B"), each = n_per)),
  t1      = c(rnorm(n_per, 0,   1), rnorm(n_per, 0,   1)),
  t2      = c(rnorm(n_per, 0.4, 1), rnorm(n_per, 0.8, 1)),
  t3      = c(rnorm(n_per, 0.8, 1), rnorm(n_per, 1.6, 1))
)

# Rows are the between-subjects effect (labeled A), the within-subjects
# effect (labeled B), and their interaction, each with all four
# multivariate criteria, followed by a row recording the sum-of-squares
# type. The multivariate tests make no sphericity assumption, which is
# what recommends them over the univariate repeated measures F and its
# epsilon corrections.
manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")
#>                 effect   statistic_name statistic_value  F_approx df_1 df_2
#> 1                    A           Pillai      0.08540425 1.6808262    1   18
#> 2                    A            Wilks      0.91459575 1.6808262    1   18
#> 3                    A Hotelling-Lawley      0.09337924 1.6808262    1   18
#> 4                    A              Roy      0.09337924 1.6808262    1   18
#> 5                    B           Pillai      0.37096011 5.0126566    2   17
#> 6                    B            Wilks      0.62903989 5.0126566    2   17
#> 7                    B Hotelling-Lawley      0.58972431 5.0126566    2   17
#> 8                    B              Roy      0.58972431 5.0126566    2   17
#> 9          interaction           Pillai      0.08217611 0.7610359    2   17
#> 10         interaction            Wilks      0.91782389 0.7610359    2   17
#> 11         interaction Hotelling-Lawley      0.08953364 0.7610359    2   17
#> 12         interaction              Roy      0.08953364 0.7610359    2   17
#> 13 sum_of_squares_type             <NA>      3.00000000        NA   NA   NA
#>       p_value
#> 1  0.21119035
#> 2  0.21119035
#> 3  0.21119035
#> 4  0.21119035
#> 5  0.01944306
#> 6  0.01944306
#> 7  0.01944306
#> 8  0.01944306
#> 9  0.48245247
#> 10 0.48245247
#> 11 0.48245247
#> 12 0.48245247
#> 13         NA
```
