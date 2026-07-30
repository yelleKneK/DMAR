# Mixed-Model ANOVA F-Ratios for One- and Two-Way Designs

Computes the classical *F*-ratios for one- and two-way ANOVA designs
when one or more factors are random rather than fixed, using the
expected-mean-square (EMS) rules that determine the correct denominator
for each test (Searle, Casella, & McCulloch, 1992). The function is the
closed-form alternative to fitting via
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) and is the standard
treatment in classical psychometrics and design-of-experiments texts
(Maxwell, Delaney, & Kelley, 2027, Ch. 10).

## Usage

``` r
mixed_anova(
  data,
  outcome,
  factor_A,
  factor_B = NULL,
  A_type = c("fixed", "random"),
  B_type = c("random", "fixed")
)
```

## Arguments

- data:

  A `data.frame` containing the response, the factor(s), and (when both
  factors are crossed) the subject identifier.

- outcome:

  Character name of the response column.

- factor_A:

  Character name of factor *A*.

- factor_B:

  Character name of factor *B*, or `NULL` for one-way designs. Default
  `NULL`.

- A_type:

  One of `"fixed"` (default) or `"random"`: the EMS classification of
  factor *A*.

- B_type:

  One of `"fixed"` or `"random"` (default): the EMS classification of
  factor *B*. Ignored when `factor_B` is `NULL`.

## Value

A `data.frame` with one row per testable effect. Columns: `effect`,
`ss`, `df`, `ms`, `denominator`, `F_value`, `p_value`.

## Details

**One way design (only `factor_A`).**

- Both *A* fixed and *A* random use the same observed *F*-ratio \\MS_A /
  MS\_{\mathrm{within}}\\; the test of "is there an effect of *A*?" is
  identical. The interpretation differs: the random-effects test asks
  whether the variance component \\\sigma^2_A\\ is zero.

**Two-way design, both fixed (Model I).**

- \\F_A = MS_A / MS\_{AB}\\ (when interaction is present in the model
  and treated as error) or \\F_A = MS_A / MS\_{\mathrm{within}}\\ (when
  interaction is pooled into error). The function uses
  \\MS\_{\mathrm{within}}\\ as the denominator throughout for Model I.

**Two-way design, both random (Model II).**

- \\F_A = MS_A / MS\_{AB}\\, \\F_B = MS_B / MS\_{AB}\\, \\F\_{AB} =
  MS\_{AB} / MS\_{\mathrm{within}}\\.

**Two-way mixed design (Model III, e.g., *A* fixed, *B* random).**

- \\F_A = MS_A / MS\_{AB}\\ (fixed factor against the interaction with
  the random factor)

- \\F_B = MS_B / MS\_{\mathrm{within}}\\ (random factor against the
  within-cell residual)

- \\F\_{AB} = MS\_{AB} / MS\_{\mathrm{within}}\\.

**Balanced data assumed.** The classical EMS rules require equal cell
sizes. The function errors out on unbalanced data and recommends a
mixed-effects fit via [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

**Sums of squares.** For the balanced designs this function targets, the
Type I, Type II, and Type III sums of squares for each effect coincide,
so the decomposition is unambiguous and no sums-of-squares type needs to
be chosen (Maxwell, Delaney, & Kelley, 2027, Ch. 7). The returned object
carries a numeric `sum_of_squares_type` attribute equal to 3, with the
understanding that it equals Types I and II here; it records the
convention without implying a choice that would matter for these
designs.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 7 on the sums of squares for balanced
designs and Chapter 10 on random and mixed effects.)

Searle, S. R., Casella, G., & McCulloch, C. E. (1992). *Variance
components*. Wiley.

## See also

[`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`anova_within_two_way`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html)

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
[`tukey_kramer_ci()`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two-way mixed design: A fixed, B random.
set.seed(113)
grid <- expand.grid(A = factor(1:3), B = factor(1:5), rep = 1:4)
grid$y <- with(grid, 0.8 * as.integer(A) + rep(rnorm(5, 0, 1.5), each = 1)[as.integer(B)] +
                      rnorm(nrow(grid), 0, 1))
mixed_anova(grid, outcome = "y", factor_A = "A", factor_B = "B",
            A_type = "fixed", B_type = "random")
#>   effect         ss df         ms denominator    F_value      p_value
#> 1      A  42.287564  2 21.1437821       MS_AB 54.7505466 2.148782e-05
#> 2      B 153.417534  4 38.3543835   MS_within 31.2357357 1.817787e-12
#> 3    A:B   3.089472  8  0.3861839   MS_within  0.3145075 9.565106e-01
```
