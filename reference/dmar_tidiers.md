# Tidy and Glance Methods for DMAR Result Tables

A DMAR function returns a tidy `data.frame` built to be read: one row
per quantity, a numeric `value` column, and a display layer that rounds
the printed values, not the stored ones (see
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)).
The tidy verbs [`tidy`](https://generics.r-lib.org/reference/tidy.html)
and [`glance`](https://generics.r-lib.org/reference/glance.html) give
the same numbers in the two shapes a programmer usually wants instead:
one row per term with a typed column for each quantity, and a one-row
summary of the result as a whole. This page states that contract once,
for the confidence interval family, the post hoc family, the contrast
tests, and the power-based sample size planners.

## Usage

``` r
# S3 method for class 'dmar_cliff_delta'
tidy(x, ...)

# S3 method for class 'dmar_cliff_delta'
glance(x, ...)

# S3 method for class 'dmar_contrast_test'
tidy(x, ...)

# S3 method for class 'dmar_contrast_test'
glance(x, ...)

# S3 method for class 'dmar_df_barnard_rubin'
tidy(x, ...)

# S3 method for class 'dmar_df_barnard_rubin'
glance(x, ...)

# S3 method for class 'dmar_df_missing'
tidy(x, ...)

# S3 method for class 'dmar_df_missing'
glance(x, ...)

# S3 method for class 'dmar_ci_long'
tidy(x, ...)

# S3 method for class 'dmar_ci_long'
glance(x, ...)

# S3 method for class 'dmar_ci_anova'
tidy(x, ...)

# S3 method for class 'dmar_ci_anova'
glance(x, ...)

# S3 method for class 'dmar_anova_es'
tidy(x, ...)

# S3 method for class 'dmar_anova_es'
glance(x, ...)

# S3 method for class 'dmar_post_hoc_ci'
tidy(x, ...)

# S3 method for class 'dmar_post_hoc_ci'
glance(x, ...)

# S3 method for class 'dmar_ss_power'
tidy(x, ...)

# S3 method for class 'dmar_ss_power'
glance(x, ...)

# S3 method for class 'dmar_ss_aipe'
tidy(x, ...)

# S3 method for class 'dmar_ss_aipe'
glance(x, ...)

# S3 method for class 'dmar_tbl'
tidy(x, ...)

# S3 method for class 'dmar_tbl'
glance(x, ...)

# S3 method for class 'dmar_content_validity'
tidy(x, ...)

# S3 method for class 'dmar_content_validity'
glance(x, ...)

# S3 method for class 'dmar_dmacs'
tidy(x, ...)

# S3 method for class 'dmar_dmacs'
glance(x, ...)

# S3 method for class 'dmar_measurement_invariance'
tidy(x, ...)

# S3 method for class 'dmar_measurement_invariance'
glance(x, ...)

# S3 method for class 'dmar_measurement_alignment'
tidy(x, ...)

# S3 method for class 'dmar_measurement_alignment'
glance(x, ...)

# S3 method for class 'dmar_factor_congruence'
tidy(x, ...)

# S3 method for class 'dmar_factor_congruence'
glance(x, ...)

# S3 method for class 'dmar_ss_power_sensitivity'
tidy(x, ...)

# S3 method for class 'dmar_ss_power_sensitivity'
glance(x, ...)

# S3 method for class 'dmar_factor_analysis'
tidy(x, ...)

# S3 method for class 'dmar_factor_analysis'
glance(x, ...)

# S3 method for class 'dmar_factor_count'
tidy(x, ...)

# S3 method for class 'dmar_parallel_analysis'
glance(x, ...)

# S3 method for class 'dmar_probability_of_superiority_paired'
tidy(x, ...)

# S3 method for class 'dmar_probability_of_superiority_paired'
glance(x, ...)

# S3 method for class 'dmar_ci_ratio_of_means'
tidy(x, ...)

# S3 method for class 'dmar_ci_ratio_of_means'
glance(x, ...)

# S3 method for class 'dmar_vargha_delaney_A'
tidy(x, ...)

# S3 method for class 'dmar_vargha_delaney_A'
glance(x, ...)
```

## Arguments

- x:

  A DMAR result object carrying one of the classes listed above.

- ...:

  Unused, present for consistency with the generics.

## Value

`tidy()` returns a `data.frame` with one row per term and
broom-convention column names. `glance()` returns a one-row `data.frame`
summarizing the result as a whole. Both return values at full precision.

## Details

**What the verbs return.** `tidy(x)` returns a `data.frame` with one row
per term, where a term is whatever the family produces one of: a
parameter estimate, a contrast, a planned design. Its columns follow the
naming convention the broom ecosystem uses, which separates words with
dots rather than the underscores DMAR uses everywhere else: `term`,
`estimate`, `se`, `statistic`, `p_value`, `ci_lower`, `ci_upper`, and
`conf_level`. A method reports the subset of those columns its family
can fill, plus any column the family genuinely adds, such as
`p_adjusted` for a multiplicity-adjusted set of comparisons or `power`
for a sample size planner.

`glance(x)` returns a one-row `data.frame` summarizing the result as a
whole, in the same dotted convention: how many comparisons were made, at
what confidence level, with which planning inputs. When a result has a
single estimand and nothing further to say at the model level, as for a
lone effect size and its confidence interval, `glance()` coincides with
`tidy()`. That is expected rather than a defect, since there is no
model-level quantity that the single row does not already carry.

Neither verb rounds. The `dmar_tbl` layer formats what is printed, while
`tidy()` and `glance()` return full precision, which is what makes them
the right input to a downstream calculation or plot.

**Why broom is not a dependency.** The `tidy()` and `glance()` generics
live in generics, a small package that holds the generics and little
else. broom imports them from there, and so does DMAR, which registers
its methods against
[`generics::tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`generics::glance`](https://generics.r-lib.org/reference/glance.html)
rather than against broom itself. A user with broom or the tidymodels
stack loaded gets DMAR methods on the generic they already call; a user
with neither installed can still call
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
directly. DMAR never loads broom, and does not need it installed.

**The families and the classes they carry.** Each family tags its return
with a leading S3 class, ahead of `dmar_tbl` and `data.frame`, so the
verbs dispatch while printing and data-frame behavior are untouched.

|  |  |  |
|----|----|----|
| **S3 class** | **Family** | **One `tidy()` row is** |
| `dmar_ci_long` | confidence intervals, long form | an estimate and its limits |
| `dmar_ci_anova` | ANOVA effect size intervals | an effect and estimand, its estimate and limits |
| `dmar_anova_es` | ANOVA effect size estimates | an effect and estimand, its estimate |
| `dmar_post_hoc_ci` | simultaneous intervals | one pairwise or one contrast comparison |
| `dmar_contrast_test` | contrast tests | one contrast, with its test and its interval |
| `dmar_ss_power` | sample size planners | a planned size and the power it buys |
| `dmar_ss_power_sensitivity` | planner sensitivity studies | a planned size and two powers |

**The confidence interval family.** Two classes cover the two output
shapes.

- `dmar_ci_long`:

  Long-format interval tables, with rows for `ci_lower` and `ci_upper`
  and, when the function reports one, an estimate row whose `term` is
  the name of the parameter. Carried by
  [`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
  [`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
  [`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
  [`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md),
  and
  [`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md).

- `dmar_ci_anova`:

  Wide-format ANOVA effect size interval tables, with one row per effect
  and columns for the effect name, the point estimate of each estimand
  with its limits (`eta_squared_partial`,
  `eta_squared_partial_ci_lower`, `eta_squared_partial_ci_upper`, and
  the same for the total), and the design metadata. Carried by
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
  and
  [`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
  and by
  [`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
  whose one estimate sits beside `ci_lower` and `ci_upper`.

Both produce a `data.frame` with `term`, `estimate`, `ci_lower`,
`ci_upper`, and `conf_level` (`NA` when the object records no level).
For `dmar_ci_long` it is one row, and `glance()` calls `tidy()`, since
the row is already the whole result. For `dmar_ci_anova` it is one row
per effect and estimand, the rows of one estimand together: `term` names
the estimand (`"eta_squared_total"`, `"eta_squared_partial"`), an
`effect` column names the effect, a `stratum` column follows it for a
fit with an `Error()` term, and the limits of the estimand are
`ci_lower` and `ci_upper` (`NA` where no interval is formed, as for the
total of a factorial effect). Its `glance()` coincides with `tidy()`
when that is one row, and otherwise is one row giving the estimands
(`term`, separated by commas when there are two), the number of effects
(`n_effects`), `N`, and `conf_level`. The point estimates of
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)
and
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md)
carry `dmar_anova_es`, whose `tidy()` has the same rows with `term`,
`effect`, `stratum` (for an `aovlist` fit), and `estimate`, and whose
`glance()` gives `term`, `n_effects`, and `N` (`NA` for
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
which reports none). The two ratio of means functions carry
`dmar_ci_ratio_of_means` ahead of `dmar_ci_long`, because a Fieller
confidence set need not be an interval: their one row adds, after the
limits, `excluded_lower` and `excluded_upper` (the ends of the gap of a
set that is the complement of an interval, `NA` otherwise) and
`confidence_set`, the form of the set, so that a complement is not read
from its limits of `-Inf` and `Inf` as the whole line (`NA` when a
bootstrap interval was withheld and the limits are `NA`).
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)
carries `dmar_vargha_delaney_A`, whose one `tidy()` row, named `"A"`,
adds the test of stochastic equality (`se`, `statistic`, `p_value`) and
the two sample sizes; its `glance()` coincides with it. Its verbs have
`df` and `test` after `p_value` and a logical column
`interval_degenerate`.
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
and
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md)
carry `dmar_cliff_delta` and `dmar_probability_of_superiority_paired`,
whose verbs are the default ones with the same logical column, `TRUE`
when the estimated variance is 0 and the limits are not a confidence
interval.

A table with no method of its own falls to the default methods for
`dmar_tbl`. Its `tidy()` keeps the confidence level the table records in
a `conf_level` column: beside the interval of a wide table (`NA` when no
interval was formed), and, for a long table of `term` and `value` rows
with `ci_lower` and `ci_upper` rows, on the rows of the interval, its
estimate and its two limits, with `NA` on every other row (an interval
held in rows named for its estimate, as
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
holds the RMSEA's in `rmsea_ci_lower` and `rmsea_ci_upper`, is marked
the same way, at its own level, except the interval for the standardized
mean difference that
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
and
[`proportion_of_superiority`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md)
echo in `smd_ci_lower` and `smd_ci_upper`, which their own intervals are
computed from and whose rows keep `NA`). The rows carry the level and
say which estimate the limits belong to, so that
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
on a `tidy()` result, or on rows of several `tidy()` results stacked
with [`rbind()`](https://rdrr.io/r/base/cbind.html), states the level
and the estimate of its own interval, as it does on the table itself. In
a wide table that carries an interval, the `estimate` is the column the
interval belongs to (the table's own `estimate` column, or its first
numeric column that is not a limit, standard error, test statistic,
degrees of freedom, *p*-value, confidence level, or count), and the
other numeric columns follow under their own names. When that column has
a name of its own (`kappa` of
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
`D2` of
[`ci_mahalanobis`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
`eta_squared_partial` of
[`simple_effects_AB`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md)),
an `estimand` column after `estimate` keeps the name, since `term` names
the row (from the table's label columns, or the row number when it has
none) and not the quantity;
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
labels the sentence of the row with it, as it labels the sentence of the
table's own row. The `tidy()` of
[`content_validity_index`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
whose estimate is the item's `i_cvi`, carries the same column, and so
does that of
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
whose estimate is the parameter's FMI (`"fmi"`) while `term` names the
parameter. The `tidy()` of
[`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md),
whose RMSEA intervals sit in its `rmsea_ci_lower` and `rmsea_ci_upper`
columns beside the other fit measures of each count, keeps the count as
`estimate` and adds a `conf_level` column holding the level of those
intervals, 0.90 unless the call set another. The default `glance()` of a
long table is one row with a column for each term; for a table with an
interval it adds an `estimand` column, the name of the column that holds
the estimate the interval belongs to (the corrected correlation of
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
listed after the observed one, and `NA` for a table that reports its
limits alone, as
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md) and
[`ci_snr`](https://yelleknek.github.io/DMAR/reference/ci_snr.md) do),
and a `conf_level` column, and names the column of the *p*-value the
table's sentence reports, so that
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
on the row gives the table's own sentence.

**Unequal tail areas.** An interval formed from `alpha_lower` and
`alpha_upper` that differ records its level as the coverage \\1 -
\alpha\_{lower} - \alpha\_{upper}\\ and the two tail areas beside it.
Every verb above that returns a `conf_level` for such an interval then
adds `alpha_lower` and `alpha_upper` columns after it, holding the two
areas on the rows of the interval (`NA` on the other rows of a long
`tidy()`), so that
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
on the rows states the split, or the one-sided bound, as it does on the
table. An equal-tailed interval, the usual one, adds neither column, so
its columns are those it has without the record of the split; a stack of
an equal-tailed and an unequal-tailed result therefore has different
columns, which [`rbind()`](https://rdrr.io/r/base/cbind.html) refuses
and
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
fills with `NA`, read as equal tails.

**The post hoc family.**
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
[`ci_dunnett_t3`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
[`ci_scheffe`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
and
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)
all carry `dmar_post_hoc_ci`. Their source table is wide, with one row
per comparison: a `contrast` label, a point estimate (`mean_difference`
for the pairwise and many-to-one procedures, `estimate` for Scheffe), a
standard error, a test statistic, the `ci_lower` and `ci_upper` of the
simultaneous interval, and the multiplicity-adjusted `p_adjusted`.
`tidy()` maps that to `term`, `estimate`, `ci_lower`, `ci_upper`,
`p_adjusted`, and `conf_level`, one row per comparison. `glance()`
describes the family of comparisons as a whole: how many there were, and
the simultaneous confidence level they hold jointly.

**The contrast tests.**
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
carries `dmar_contrast_test`. Its source table is wide, with one row per
contrast: a `contrast` label, the estimate \\\hat{\psi} = \sum_i c_i
\bar{Y}\_i\\, its standard error, the *t*-statistic and the degrees of
freedom it is referred to, the unadjusted *p*-value, the
multiplicity-adjusted `p_adjusted`, and the `ci_lower` and `ci_upper`
limits. `tidy()` renames those to `term`, `estimate`, `ci_lower`,
`ci_upper`, `statistic`, `df`, `p_value`, `p_adjusted`, and
`conf_level`, one row per contrast. Both *p*-values are kept, because
the pair is what a contrast table is read for: what the contrast would
show on its own, and what it shows once the family it belongs to is
accounted for.

Where a post hoc procedure fixes its adjustment as part of the method, a
contrast test chooses one, and the same weights tested under
`adjust = "none"` and under `adjust = "tukey"` are two different
inferences. `glance()` therefore records the choice alongside the
family-level numbers: `n_contrasts`, `adjust`, `var_equal`, the smallest
adjusted *p*-value `p_adjusted_min`, and `conf_level`. `adjust` and
`var_equal` name a procedure rather than measure a quantity, so this
one-row summary, unlike a DMAR result table, is not numeric throughout.
The post hoc and contrast methods, and wide result tables that carry
adjustment metadata, retain `ci_adjustment` and `p_adjustment` in both
`tidy()` and `glance()`. These are distinct because an interval and a
test can use different multiplicity procedures.
[`results_sentence()`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
reports an adjusted *p*-value beside a simultaneous interval only when
the recorded procedures match.

**The power-based sample size planners.** A planner in the `ss_power_*`
family returns a long table with a row for the recommended sample size,
a row for the realized power, and rows echoing the planning inputs. A
planner that reports one size and one power for one design tags its
return `dmar_ss_power`. This covers the closed-form effect size planners
([`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_r`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md)),
the contrast and ANCOVA planners
([`ss_power_c`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_sc`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_contrast`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)),
the ANOVA and cluster designs
([`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_factorial_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_split_plot_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md),
[`ss_power_rm_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md)),
and the mediation planner
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
whose reported power is the joint power to detect the indirect effect
(its row `power_ab`) and whose component path powers `glance()` carries
as extra columns.

The size `tidy()` reports is the design's planning unit: per group, per
cell, per subject, or per cluster. The one-way ANOVA planner, whose
natural unit is the total, is summarized by its total \\N\\. A design
that reports two group sizes reports one of them beside the realized
power, falling through to the total \\N\\ when the per-group sizes are
unequal, and `glance()` keeps every group size as a column so none is
lost. A planner whose result spans several effects, with no single
size-and-power summary to give, returns a plain `dmar_tbl` and does not
gain these verbs at all.

The Monte Carlo sensitivity siblings
[`ss_power_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md)
and
[`ss_power_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md)
report two powers at one planned sample size, the empirical (simulated)
power and the analytic power, and comparing the two is the object of the
study. They carry `dmar_ss_power_sensitivity` instead: `tidy()` places
both powers beside the planned `sample_size`, and `glance()` adds the
simulated distribution of the estimator.

**Adding a planner to the family.** A planner opts in by setting
`dmar_ss_power` as a leading class before routing its return through
`.as_dmar_tbl()`. The rows the verbs read are named in the internal
vectors `.SS_POWER_SIZE_TERMS` and `.SS_POWER_POWER_TERMS`. A planner
whose size or power row is not named there reports `NA` rather than
failing, so a new row name has to be added to those vectors when a
planner introduces one. `glance()` leaves out only the two rows `tidy()`
reads, the size and the power, and keeps every other row as a column,
other powers included: the marginal powers of the composite planners
(`power_a`, `power_b`, and `power_ab` for parameters labeled `a`, `b`,
and `ab` in
[`ss_power_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md))
and the path powers of
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md).

## See also

[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md) for
the printing layer these tables share, and the "Reading DMAR result
tables" vignette for the wider output convention.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A single interval: tidy() and glance() coincide, because there is
# nothing at the model level the one row does not already carry.
res <- ci_r(r = 0.5, N = 50)
generics::tidy(res)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1    r      0.5 0.2574879 0.6832563       0.95
generics::glance(res)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1    r      0.5 0.2574879 0.6832563       0.95

# A family of simultaneous intervals: one tidy() row per comparison,
# one glance() row describing the family.
set.seed(113)
y <- c(rnorm(10, 0), rnorm(10, 1), rnorm(10, 2))
g <- factor(rep(c("a", "b", "c"), each = 10))
gh <- ci_games_howell(y, group = g)
generics::tidy(gh)
#>    term estimate    ci_lower ci_upper   p_adjusted conf_level ci_adjustment
#> 1 b - a 1.665520  0.62193264 2.709108 1.976750e-03       0.95  games_howell
#> 2 c - a 2.809851  1.61418391 4.005519 3.511181e-05       0.95  games_howell
#> 3 c - b 1.144331 -0.01274098 2.301403 5.283057e-02       0.95  games_howell
#>   p_adjustment
#> 1 games_howell
#> 2 games_howell
#> 3 games_howell
generics::glance(gh)
#>   n_contrasts conf_level ci_adjustment p_adjustment
#> 1           3       0.95  games_howell games_howell

# A set of contrasts: tidy() keeps both the unadjusted and the
# adjusted p-value, and glance() names the adjustment that produced
# the second of them.
fit <- aov(bdi_post ~ condition, data = depression_bdi)
ct <- contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
generics::tidy(ct)
#>                  term estimate   ci_lower  ci_upper statistic df    p_value
#> 1      placebo - ssri      4.9 -2.0793312 11.879331 1.7407322 27 0.09311706
#> 2    wait_list - ssri      6.7 -0.2793312 13.679331 2.3801849 27 0.02462568
#> 3 wait_list - placebo      1.8 -5.1793312  8.779331 0.6394526 27 0.52791718
#>   p_adjusted conf_level ci_adjustment p_adjustment
#> 1  0.2088175       0.95         tukey        tukey
#> 2  0.0617362       0.95         tukey        tukey
#> 3  0.7998146       0.95         tukey        tukey
generics::glance(ct)
#>   n_contrasts adjust var_equal p_adjusted_min conf_level ci_adjustment
#> 1           3  tukey      TRUE      0.0617362       0.95         tukey
#>   p_adjustment
#> 1        tukey

# A sample size planner: tidy() gives the size and the power it buys,
# glance() adds the planning inputs that produced them.
plan <- ss_power_smd(true_smd = 0.5, desired_power = 0.80)
generics::tidy(plan)
#>          term estimate     power
#> 1 sample_size       64 0.8014596
generics::glance(plan)
#>          term estimate     power noncentral_t_parm true_smd desired_power
#> 1 sample_size       64 0.8014596          2.828427      0.5           0.8
#>   alpha_level tails
#> 1        0.05     2
```
