# Tidy and Glance Methods for DMAR Result Tables

A DMAR function returns a tidy `data.frame` built to be read: one row
per quantity, a numeric `value` column, and a display layer that rounds
sensibly on the way to the console (see
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
# S3 method for class 'dmar_contrast_test'
tidy(x, ...)

# S3 method for class 'dmar_contrast_test'
glance(x, ...)

# S3 method for class 'dmar_ci_long'
tidy(x, ...)

# S3 method for class 'dmar_ci_long'
glance(x, ...)

# S3 method for class 'dmar_ci_anova'
tidy(x, ...)

# S3 method for class 'dmar_ci_anova'
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

# S3 method for class 'dmar_ss_power_sensitivity'
tidy(x, ...)

# S3 method for class 'dmar_ss_power_sensitivity'
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
| `dmar_ci_anova` | ANOVA effect size intervals | an effect size and its limits |
| `dmar_post_hoc_ci` | simultaneous intervals | one pairwise or one contrast comparison |
| `dmar_contrast_test` | contrast tests | one contrast, with its test and its interval |
| `dmar_ss_power` | sample size planners | a planned size and the power it buys |
| `dmar_ss_power_sensitivity` | planner sensitivity studies | a planned size and two powers |

**The confidence interval family.** Two classes cover the two output
shapes.

- `dmar_ci_long`:

  Long-format interval tables, with rows for `lower_limit` and
  `upper_limit` and, when the function reports one, an estimate row
  whose `term` is the name of the parameter. Carried by
  [`ci_cc`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
  [`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
  [`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
  and
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md).

- `dmar_ci_anova`:

  Wide-format ANOVA effect size interval tables, with one row and
  columns for the effect name, the point estimate, the limits, and the
  design metadata. Carried by
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
  [`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
  [`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
  and
  [`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).

Both produce a one-row `data.frame` with `term`, `estimate`, `ci_lower`,
`ci_upper`, and, when the object records it, `conf_level`. `glance()` on
either class calls `tidy()`, since the row is already the whole result.

**The post hoc family.**
[`tukey_kramer_ci`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
[`games_howell_ci`](https://yelleknek.github.io/DMAR/reference/games_howell_ci.md),
[`scheffe_ci`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md),
and
[`dunnett_ci`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md)
all carry `dmar_post_hoc_ci`. Their source table is wide, with one row
per comparison: a `contrast` label, a point estimate (`mean_difference`
for the pairwise and many-to-one procedures, `contrast_value` for
Scheffe), a standard error, a test statistic, the `lower_limit` and
`upper_limit` of the simultaneous interval, and the
multiplicity-adjusted `p_adjusted`. `tidy()` maps that to `term`,
`estimate`, `ci_lower`, `ci_upper`, `p_adjusted`, and `conf_level`, one
row per comparison. `glance()` describes the family of comparisons as a
whole: how many there were, and the simultaneous confidence level they
hold jointly.

**The contrast tests.**
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
carries `dmar_contrast_test`. Its source table is wide, with one row per
contrast: a `contrast` label, the estimate \\\hat{\psi} = \sum_i c_i
\bar{Y}\_i\\, its standard error, the *t*-statistic and the degrees of
freedom it is referred to, the unadjusted *p*-value, the
multiplicity-adjusted `p_adj`, and the `conf_lower` and `conf_upper`
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
and whose component path powers `glance()` carries as extra columns.

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
planner introduces one.

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
res <- ci_cc(r = 0.5, n = 50)
generics::tidy(res)
#>      term estimate  ci_lower  ci_upper conf_level
#> 1 est_cor      0.5 0.2574879 0.6832563       0.95
generics::glance(res)
#>      term estimate  ci_lower  ci_upper conf_level
#> 1 est_cor      0.5 0.2574879 0.6832563       0.95

# A family of simultaneous intervals: one tidy() row per comparison,
# one glance() row describing the family.
set.seed(113)
y <- c(rnorm(10, 0), rnorm(10, 1), rnorm(10, 2))
g <- factor(rep(c("a", "b", "c"), each = 10))
gh <- games_howell_ci(y, group = g)
generics::tidy(gh)
#>    term estimate    ci_lower ci_upper   p_adjusted conf_level
#> 1 b - a 1.665520  0.62193264 2.709108 1.976750e-03       0.95
#> 2 c - a 2.809851  1.61418391 4.005519 3.511181e-05       0.95
#> 3 c - b 1.144331 -0.01274098 2.301403 5.283057e-02       0.95
generics::glance(gh)
#>   n_contrasts conf_level
#> 1           3       0.95

# A set of contrasts: tidy() keeps both the unadjusted and the
# adjusted p-value, and glance() names the adjustment that produced
# the second of them.
fit <- aov(weight ~ group, data = PlantGrowth)
ct <- contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
generics::tidy(ct)
#>          term estimate   ci_lower  ci_upper statistic df     p_value p_adjusted
#> 1 trt1 - ctrl   -0.371 -1.0622161 0.3202161 -1.330791 27 0.194387880 0.39087114
#> 2 trt2 - ctrl    0.494 -0.1972161 1.1852161  1.771996 27 0.087681675 0.19799599
#> 3 trt2 - trt1    0.865  0.1737839 1.5562161  3.102787 27 0.004459236 0.01200642
#>   conf_level
#> 1       0.95
#> 2       0.95
#> 3       0.95
generics::glance(ct)
#>   n_contrasts adjust var_equal p_adjusted_min conf_level
#> 1           3  tukey      TRUE     0.01200642       0.95

# A sample size planner: tidy() gives the size and the power it buys,
# glance() adds the planning inputs that produced them.
plan <- ss_power_smd(smd = 0.5, desired_power = 0.80)
generics::tidy(plan)
#>          term estimate     power
#> 1 sample_size       64 0.8014596
generics::glance(plan)
#>          term estimate     power noncentral_t_parm supposed_smd desired_power
#> 1 sample_size       64 0.8014596          2.828427          0.5           0.8
#>   alpha_level tails
#> 1        0.05     2
```
