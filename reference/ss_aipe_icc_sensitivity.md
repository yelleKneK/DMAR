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
  true_rho = NULL,
  estimated_rho = NULL,
  k,
  width,
  assurance = NULL,
  specified_N = NULL,
  conf_level = 0.95,
  type = c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)", "ICC(1,k)", "ICC(2,k)", "ICC(3,k)"),
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_icc_sensitivity_result.csv"
)
```

## Arguments

- true_rho:

  Population intraclass correlation coefficient (the data generating
  value), at the level matching `type`. Must lie in \\\[0, 1)\\.

- estimated_rho:

  ICC used to plan the study (the value the researcher guessed when
  invoking
  [`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)).
  Supply this or `specified_N` but not both.

- k:

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
  would recommend); incompatible with `estimated_rho`.

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
  realized ICC on each replication is comparable to `true_rho`.

- G:

  Number of Monte Carlo replications; defaults to 1000. Increase (e.g.,
  5000 or 10000) for stable empirical-coverage estimates.

- print_iter:

  Logical. If `TRUE` the simulation prints the iteration index after
  each replication (helpful for long runs); default `FALSE`.

- save:

  Logical. If `TRUE` the per-replication results are appended to a CSV
  file at `filename`; default `FALSE`.

- filename:

  Path used when `save = TRUE`; default
  `"ss_aipe_icc_sensitivity_result.csv"` in the current working
  directory.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo results across the `G` replications. The `term` entries are:
`"mean_icc"`, `"median_icc"`, `"sd_icc"` (mean / median / SD of the *G*
observed ICC estimates); `"mean_ci_width"`, `"median_ci_width"`,
`"sd_ci_width"` (corresponding summaries of the realized interval
widths); `"pct_ci_less_w"` (proportion of intervals at or below the
planning width `width`); `"pct_ci_miss_low"` and `"pct_ci_miss_high"`
(tail-specific non-coverage of `true_rho`); `"total_type_I_error"`
(overall empirical non-coverage of `true_rho`); plus the input echoes
`"total_N"`, `"k"`, `"true_rho"`, `"estimated_rho"`, `"width"`,
`"conf_level"`, and `"assurance"` (present only when an assurance was
supplied). The ICC type is not a row; it is stored as the `"icc_type"`
attribute on the returned object so the `value` column stays numeric.

## Details

Sample size planning for the intraclass correlation coefficient under
the Accuracy in Parameter Estimation framework chooses *n* so that the
expected (or, with assurance, the high-probability) confidence interval
width is no larger than `width` (Bonett, 2002; Kelley & Maxwell, 2003).
Because the procedure assumes the planning value `estimated_rho` matches
the population value `true_rho`, in practice the realized width will
deviate from the planned width whenever the planning value is wrong.
This sensitivity analysis quantifies the deviation by Monte Carlo
simulation: the planned *n* is obtained from
[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)
with `estimated_rho`, then samples are drawn from the *true* population
(with population ICC `true_rho`) and the realized confidence interval
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
construction. For the single-rater forms `true_rho` therefore maps
directly to the subject-variance share; for the average-of-\\k\\ forms
`true_rho` is first mapped to the single-rater scale through the inverse
Spearman-Brown relation \\\rho = \rho_k / \[k - (k - 1)\rho_k\]\\, so
that the population ICC at the average-of-\\k\\ level equals `true_rho`.
Coverage is checked against `true_rho` on its own scale, matching the
scale on which
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md) reports each
form.

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
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
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
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Reduced G and a wide target width keep this Monte Carlo example fast;
# raise G (e.g., 1000 or more) for stable empirical-coverage estimates.

# Well-specified case: plan with the true population ICC.
set.seed(113)
ss_aipe_icc_sensitivity(
  true_rho      = 0.70,
  estimated_rho = 0.70,
  k             = 3,
  width         = 0.40,
  conf_level    = 0.95,
  G             = 25,
  print_iter    = FALSE
)
#>  term               value 
#>  mean_icc           0.675 
#>  median_icc         0.67  
#>  sd_icc             0.105 
#>  mean_ci_width      0.389 
#>  median_ci_width    0.404 
#>  sd_ci_width        0.0897
#>  pct_ci_less_w      0.44  
#>  pct_ci_miss_low    0.04  
#>  pct_ci_miss_high   0     
#>  total_type_I_error 0.04  
#>  total_N            19    
#>  k                  3     
#>  true_rho           0.7   
#>  estimated_rho      0.7   
#>  width              0.4   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%

# Misspecified case: the planner used .70 but the truth is .50.
# The realized interval widths will tend to be wider than the target.
set.seed(113)
ss_aipe_icc_sensitivity(
  true_rho      = 0.50,
  estimated_rho = 0.70,
  k             = 3,
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
#>  sd_ci_width        0.0717
#>  pct_ci_less_w      0.12  
#>  pct_ci_miss_low    0.04  
#>  pct_ci_miss_high   0     
#>  total_type_I_error 0.04  
#>  total_N            19    
#>  k                  3     
#>  true_rho           0.5   
#>  estimated_rho      0.7   
#>  width              0.4   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%

# Fixed-n mode: skip the planner and evaluate at a chosen sample size.
set.seed(113)
ss_aipe_icc_sensitivity(
  true_rho    = 0.50,
  specified_N = 40,
  k           = 3,
  width       = 0.40,
  G           = 25,
  print_iter  = FALSE
)
#>  term               value 
#>  mean_icc           0.512 
#>  median_icc         0.502 
#>  sd_icc             0.0824
#>  mean_ci_width      0.348 
#>  median_ci_width    0.356 
#>  sd_ci_width        0.0302
#>  pct_ci_less_w      1     
#>  pct_ci_miss_low    0.08  
#>  pct_ci_miss_high   0     
#>  total_type_I_error 0.08  
#>  total_N            40    
#>  k                  3     
#>  true_rho           0.5   
#>  estimated_rho      <NA>  
#>  width              0.4   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%
```
