# Between-Subjects Factorial ANOVA for Unbalanced Designs

Fits a between-subjects factorial ANOVA (two or more crossed factors)
and returns each effect's sum of squares, *F* test, and partial
\\\eta^2\\ and partial \\\omega^2\\ with confidence intervals, using a
sum-of-squares type the user chooses. It is built for the case that
makes the choice matter, an *unbalanced* design (unequal cell sizes),
where the Type I, II, and III sums of squares differ. For a balanced
design the three types coincide and the choice is immaterial.

## Usage

``` r
factorial_anova(formula, data, ss_type = 3L, conf_level = 0.95)
```

## Arguments

- formula:

  A two-sided [`formula`](https://rdrr.io/r/stats/formula.html) naming
  the numeric response and the crossed factors, for example `y ~ A * B`
  or `y ~ A * B * C`. Predictors are coerced to factors. Write the full
  crossing you want tested; `A * B` expands to `A + B + A:B`.

- data:

  A `data.frame` containing the response and the factors.

- ss_type:

  The sum-of-squares type: `1`, `2`, or `3` (equivalently `"I"`, `"II"`,
  or `"III"`). Default `3` (Type III), the common reporting default. See
  Details for what each type tests.

- conf_level:

  Confidence level for the effect size confidence intervals. Default
  `0.95`.

## Value

A `data.frame` (class `dmar_tbl`) with one row per effect plus a
`Residuals` row. Columns are the effect label (`effect`), the sum of
squares (`SS`), degrees of freedom (`df`), the *F* statistic (`F_value`)
and its *p*-value (`p_value`), and partial \\\eta^2\\ and partial
\\\omega^2\\ with their lower and upper confidence limits. The chosen
sum-of-squares type is recorded on the object (`attr(x, "ss_type")`) and
printed beneath the table. The residual row carries only `SS` and `df`.
Stored values keep full precision; the display rounds (see
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)).

## Details

Every sum of squares is computed as a model comparison, the increase in
error sum of squares when an effect's parameters are removed from a
model, following the model comparison development of Maxwell, Delaney,
and Kelley (2027, Chapter 7). The computation uses only base R
(`stats`); it does not depend on the car package, though it agrees with
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) to numerical
precision.

**The types as model comparisons.** A sum of squares for an effect is
the increase in the error sum of squares when the effect's parameters
are dropped from the model, `SS = E(restricted) - E(full)`. The three
conventional types differ only in which other effects the full and
restricted models hold in common (Maxwell, Delaney, and Kelley, 2027,
Chapter 7; Overall and Spiegel, 1969):

- **Type I (sequential).** Each effect is adjusted only for the effects
  listed before it in `formula`: `SS(A)`, then `SS(B | A)`, then
  `SS(A:B | A, B)`. The parts sum to the model sum of squares, but the
  answer depends on the order of the terms.

- **Type II.** Each effect is adjusted for every other effect that does
  not contain it (it respects marginality): a main effect is adjusted
  for the other main effects but not for the interactions that contain
  it. Type II is the most powerful choice when the interaction is null
  and does not depend on how the factors are coded (Overall and
  Spiegel's Method 2; Appelbaum and Cramer, 1974).

- **Type III.** Each effect is adjusted for all other effects, including
  the higher order interactions that contain it. This tests each main
  effect as a contrast on the unweighted marginal means, so it is
  computed here with sum-to-zero contrasts
  ([`contr.sum`](https://rdrr.io/r/stats/contrast.html)), which is what
  makes the Type III main-effect test the intended one. It is the
  default reported by many programs.

For a balanced design the effects are orthogonal and the three types are
identical; the distinction is a property of unbalanced (nonorthogonal)
data.

**Reading the main effects when an interaction is present.** When an
interaction is real, the marginal main-effect tests, of any type, are
usually not the question of interest; examine the interaction and the
simple effects instead (Maxwell, Delaney, and Kelley, 2027). See the
vignette *Sums of Squares in Nonorthogonal Designs* for a worked
comparison.

**Effect sizes.** Partial \\\eta^2\\ and partial \\\omega^2\\ are formed
for each effect from its *F*, its degrees of freedom, and the error
degrees of freedom, with noncentral *F* confidence intervals
([`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)).
The confidence limits are those of the interval for the population
proportion of variance the effect accounts for (Kelley, 2007). Partial
\\\eta^2\\ and partial \\\omega^2\\ are two point estimators of that
same population quantity, partial \\\omega^2\\ correcting the upward
bias of partial \\\eta^2\\, so the two estimators differ but share the
interval.

**Estimability.** All cells must be filled. If a factor combination is
empty the design is rank deficient and the factorial effects are not all
estimable; the function stops with a message rather than return a value
that depends on an arbitrary choice.

## References

Appelbaum, M. I., & Cramer, E. M. (1974). Some problems in the
nonorthogonal analysis of variance. *Psychological Bulletin, 81*(6),
335–343.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 7 on higher order between-subjects
designs.)

Overall, J. E., & Spiegel, D. K. (1969). Concerning least squares
analysis of experimental data. *Psychological Bulletin, 72*(5), 311–322.

## See also

[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md) for a
covariate-adjusted one-way design,
[`mixed_anova`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md)
for mixed-model (fixed and random) *F* ratios,
[`manova_split_plot`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md)
for the multivariate mixed design, and
[`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md)
/
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
for the effect size intervals.

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
# An unbalanced two-factor design: IQ gain in the pygmalion expectancy
# experiment, by treatment and by lower (grades 1 and 2) versus upper
# (grades 3 through 6) grades, where the expectancy effect concentrated
# in the lower grades. The cell sizes are unequal, so the Types differ.
pyg <- pygmalion
pyg$grade_band <- factor(ifelse(pyg$grade <= 2, "lower", "upper"))
factorial_anova(iq_gain ~ treatment * grade_band, data = pyg)
#>  effect               SS    df  F_value p_value eta_squared_partial
#>  treatment            1510  1   8.24    0.0044  0.0262             
#>  grade_band           1830  1   10      0.0017  0.0317             
#>  treatment:grade_band 1120  1   6.11    0.0140  0.0196             
#>  Residuals            55900 306 <NA>    <NA>    <NA>               
#>  eta_squared_partial_lower eta_squared_partial_upper omega_squared_partial
#>  0.00257                   0.0703                    0.0228               
#>  0.00453                   0.0785                    0.0283               
#>  0.00071                   0.0598                    0.0162               
#>  <NA>                      <NA>                      <NA>                 
#>  omega_squared_partial_lower omega_squared_partial_upper
#>  0.00257                     0.0703                     
#>  0.00453                     0.0785                     
#>  0.00071                     0.0598                     
#>  <NA>                        <NA>                       
#> 
#> Sum of squares: Type III
#> 
#> Confidence level: 95%

# The same design under Type II (adjusts each main effect for the other
# main effect, but not for the interaction). Here the main effect F
# statistics drop, and a warning reports that the affected noncentral F
# lower limits are clamped to 0.
factorial_anova(iq_gain ~ treatment * grade_band, data = pyg, ss_type = 2)
#> Warning: The noncentral F lower-limit clamp in conf_limits_ncf() fired for 4 of the effect size confidence intervals; the affected lower limits were clamped to 0. See ?conf_limits_ncf for the meaning of the clamp.
#>  effect               SS    df  F_value p_value eta_squared_partial
#>  treatment            798   1   4.37    0.0374  0.0141             
#>  grade_band           790   1   4.33    0.0384  0.0139             
#>  treatment:grade_band 1120  1   6.11    0.0140  0.0196             
#>  Residuals            55900 306 <NA>    <NA>    <NA>               
#>  eta_squared_partial_lower eta_squared_partial_upper omega_squared_partial
#>  0                         0.0504                    0.0108               
#>  0                         0.0501                    0.0106               
#>  0.00071                   0.0598                    0.0162               
#>  <NA>                      <NA>                      <NA>                 
#>  omega_squared_partial_lower omega_squared_partial_upper
#>  0                           0.0504                     
#>  0                           0.0501                     
#>  0.00071                     0.0598                     
#>  <NA>                        <NA>                       
#> 
#> Sum of squares: Type II
#> 
#> Confidence level: 95%
```
