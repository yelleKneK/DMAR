# Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for an Intraclass Correlation Coefficient

Quantifies how much misspecification of the population ICC can distort
an AIPE-based sample size plan. Given a true (population) ICC and the
value used in planning, the function simulates draws of size \\n \times
k\\ from the relevant variance-components model, computes the ICC and
its *F*-distribution confidence interval on each replication via
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md), and
summarizes how often the realized interval width is below the desired
target and how often the interval covers the population value. This is
the standard sensitivity-analysis workflow described in Kelley (2007)
and Maxwell, Delaney, and Kelley (2027, Section 3.11 on sample size
planning).

## Usage

``` r
ss_aipe_icc_sensitivity(
  true_icc = NULL,
  supposed_icc = NULL,
  n_raters,
  width,
  assurance = NULL,
  specified_N = NULL,
  conf_level = 0.95,
  type = c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)", "ICC(1,k)", "ICC(2,k)", "ICC(3,k)"),
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_icc:

  Population intraclass correlation coefficient (the data generating
  value), at the level matching `type`. Must lie in \\\[0, 1)\\.

- supposed_icc:

  The supposed ICC, the planning value from which
  [`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)
  plans the study (a value the researcher supposes, not an estimate
  computed from data). Supply this or `specified_N` but not both.

- n_raters:

  Number of raters (or repeated measurements) per subject; must be at
  least 2.

- width:

  Desired full width of the (back-transformed) confidence interval on
  the ICC.

- assurance:

  Probability with which the realized interval should be no wider than
  `width` (must be `NULL` or strictly between 0 and 1; `NULL` means plan
  to the *expected* width without an assurance constraint).

- specified_N:

  Pre-specified number of subjects to evaluate (use this when you want
  the sensitivity results at a fixed *n* rather than at the *n* that
  [`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)
  would recommend); incompatible with `supposed_icc`.

- conf_level:

  Desired confidence level (i.e., 1 minus the Type I error rate);
  default 0.95.

- type:

  Which Shrout-Fleiss (1979) ICC form is being planned. One of
  `"ICC(1,1)"` (default; one-way random model), `"ICC(2,1)"` (two-way
  random model, absolute agreement), `"ICC(3,1)"` (two-way mixed model,
  consistency), or the average-of-\\k\\ versions `"ICC(1,k)"`,
  `"ICC(2,k)"`, `"ICC(3,k)"`. The simulation generates data from the
  variance- components model matching the requested type, so the
  realized ICC on each replication is comparable to `true_icc`.

- G:

  Number of Monte Carlo replications; defaults to 1000. Increase (e.g.,
  5000 or 10000) for stable empirical-coverage estimates.

- print_iter:

  Logical. If `TRUE` the simulation prints the iteration index after
  each replication (helpful for long runs); default `FALSE`.

- filename:

  An optional path for a comma separated file recording every
  replication (the realized ICC, the two confidence limits, the interval
  width, and the two tail-specific non-coverage indicators): nothing is
  written when `filename` is `NULL` (the default), a new file with a
  header row is created otherwise, an existing file at that path is
  appended to, and a throwaway run should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo results across the `G` replications. The `term` entries are:
`"mean_icc"`, `"median_icc"`, `"sd_icc"` (mean / median / SD of the *G*
observed ICC estimates); `"mean_ci_width"`, `"median_ci_width"`,
`"sd_ci_width"` (corresponding summaries of the realized interval
widths); `"pct_ci_less_width"` (proportion of intervals at or below the
planning width `width`); `"pct_ci_miss_low"` and `"pct_ci_miss_high"`
(tail-specific non-coverage of `true_icc`); `"total_type_I_error"`
(overall empirical non-coverage of `true_icc`); plus the input echoes
`"total_N"`, `"n_raters"`, `"true_icc"`, `"supposed_icc"`, `"width"`,
`"conf_level"`, and `"assurance"` (present only when an assurance was
supplied). The ICC type is not a row; it is stored as the `"icc_type"`
attribute on the returned object so the `value` column stays numeric.

## Details

Sample size planning for the intraclass correlation coefficient under
the Accuracy in Parameter Estimation framework chooses *n* so that the
expected (or, with assurance, the high-probability) confidence interval
width is no larger than `width` (Bonett, 2002; Kelley & Maxwell, 2003).
Because the procedure assumes the planning value `supposed_icc` matches
the population value `true_icc`, in practice the realized width will
deviate from the planned width whenever the planning value is wrong.
This sensitivity analysis quantifies the deviation by Monte Carlo
simulation: the planned *n* is obtained from
[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)
with `supposed_icc`, then samples are drawn from the *true* population
(with population ICC `true_icc`) and the realized confidence interval
widths are summarized.

**Data generating model.** For `type` starting with `"ICC(1,"` the
simulation uses the one-way random model: each row (subject) gets a
subject random effect, every cell adds independent Gaussian noise, and
the population ICC equals \\\sigma^2\_{\mathrm{subj}} /
(\sigma^2\_{\mathrm{subj}} + \sigma^2\_{\mathrm{err}})\\. For `type`
starting with `"ICC(2,"` the simulation also adds a rater random effect,
so the population `ICC(2,1)` equals \\\sigma^2\_{\mathrm{subj}} /
(\sigma^2\_{\mathrm{subj}} + \sigma^2\_{\mathrm{rater}} +
\sigma^2\_{\mathrm{err}})\\. For `type` starting with `"ICC(3,"` the
rater effect is fixed (centered constants), so the population `ICC(3,1)`
equals \\\sigma^2\_{\mathrm{subj}} / (\sigma^2\_{\mathrm{subj}} +
\sigma^2\_{\mathrm{err}})\\ but the two-way decomposition is used in the
estimator. Within each cell, the total variance is unity by
construction. For the single-rater forms `true_icc` therefore maps
directly to the subject-variance share; for the average-of-\\k\\ forms
`true_icc` is first mapped to the single-rater scale through the inverse
Spearman-Brown relation \\\rho = \rho_k / \[k - (k - 1)\rho_k\]\\, so
that the population ICC at the average-of-\\k\\ level equals `true_icc`.
Coverage is checked against `true_icc` on its own scale, matching the
scale on which
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md) reports each
form.

## The record of the simulation

Every sensitivity function of the package summarizes its replications
over one denominator: the replications that produced a confidence
interval (for a power sensitivity function, a test). The proportion of
intervals that met the width target (`pct_ci_less_width`), the
proportions of intervals that missed the true value, and the summaries
of the widths are all taken over those replications. A replication
without an interval is in none of these rows. It is counted, and the
returned table carries the counts in its attribute `simulation`, a list
read with `attr(x, "simulation")`:

- `attempted`:

  the number of replications asked for (`G`).

- `returned`:

  the number that produced an interval (a test).

- `failed`:

  the number that did not (`attempted - returned`).

- `denominator`:

  the number of replications the rows named above are taken over, which
  is `returned`.

- `width_target_probability`:

  the proportion of the replications with an interval whose interval met
  the width target, the value of the `pct_ci_less_width` row; `NA` for a
  function or a call without a width target, and when no replication
  produced an interval.

- `width_target_mc_ci`:

  the exact binomial (Clopper-Pearson) limits for that proportion, at
  the confidence level of the table. They describe the error of the
  simulation, which shrinks as `G` grows, and nothing else.

A function that records more about its run adds elements to the list,
and its help page names them. The line "Replications" beneath the
printed table gives `G`, and the three counts when a replication was
without an interval. With many replications without an interval the rows
describe only the studies in which the interval could be formed, so
`failed` is to be reported beside them.

## References

Bonett, D. G. (2002). Sample size requirements for estimating intraclass
correlations with desired precision. *Statistics in Medicine, 21*(9),
1331–1335. [doi:10.1002/sim.1108](https://doi.org/10.1002/sim.1108)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 10 and 11 on intraclass correlation and
reliability.)

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428.

## See also

[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
[`ss_aipe_crossed_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md),
[`ss_aipe_mean_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Reduced G and a wide target width keep this Monte Carlo example fast;
# raise G (e.g., 1000 or more) for stable empirical-coverage estimates.

# Well-specified case: plan with the true population ICC.
set.seed(113)
ss_aipe_icc_sensitivity(
  true_icc      = 0.70,
  supposed_icc  = 0.70,
  n_raters      = 3,
  width         = 0.40,
  conf_level    = 0.95,
  G             = 25,
  print_iter    = FALSE
)
#>  term               value
#>  mean_icc           0.675
#>  median_icc         0.670
#>  sd_icc             0.105
#>  mean_ci_width      0.389
#>  median_ci_width    0.404
#>  sd_ci_width        0.090
#>  pct_ci_less_width  0.440
#>  pct_ci_miss_low    0.040
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.040
#>  total_N            19   
#>  n_raters           3    
#>  true_icc           0.700
#>  supposed_icc       0.700
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 25
#> 
#> Confidence level: 95%

# Misspecified case: the planner used .70 but the truth is .50.
# The realized interval widths will tend to be wider than the target.
set.seed(113)
ss_aipe_icc_sensitivity(
  true_icc      = 0.50,
  supposed_icc  = 0.70,
  n_raters      = 3,
  width         = 0.40,
  conf_level    = 0.95,
  G             = 25,
  print_iter    = FALSE
)
#>  term               value
#>  mean_icc           0.472
#>  median_icc         0.445
#>  sd_icc             0.141
#>  mean_ci_width      0.509
#>  median_ci_width    0.537
#>  sd_ci_width        0.072
#>  pct_ci_less_width  0.120
#>  pct_ci_miss_low    0.040
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.040
#>  total_N            19   
#>  n_raters           3    
#>  true_icc           0.500
#>  supposed_icc       0.700
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 25
#> 
#> Confidence level: 95%

# Fixed-n mode: skip the planner and evaluate at a chosen sample size.
set.seed(113)
ss_aipe_icc_sensitivity(
  true_icc    = 0.50,
  specified_N = 40,
  n_raters    = 3,
  width       = 0.40,
  G           = 25,
  print_iter  = FALSE
)
#>  term               value
#>  mean_icc           0.512
#>  median_icc         0.502
#>  sd_icc             0.082
#>  mean_ci_width      0.348
#>  median_ci_width    0.356
#>  sd_ci_width        0.030
#>  pct_ci_less_width  1.000
#>  pct_ci_miss_low    0.080
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.080
#>  total_N            40   
#>  n_raters           3    
#>  true_icc           0.500
#>  supposed_icc       <NA> 
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 25
#> 
#> Confidence level: 95%
```
