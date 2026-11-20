# Sample Size Planning From the AIPE Perspective for the Least Costly Two-Phase Design

Takes planning values for a two-phase design (the standard deviation of
a costly target measurement, its correlation with an inexpensive proxy,
and the cost of a proxy relative to a target measurement) and a desired
confidence interval width for the target mean. Returns the number of
labeled units \\n\\, which carry both measurements, and the number of
pool units \\M\\, which carry the proxy, of the least costly design in
whole numbers whose interval from
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
has an expected width no larger than desired. The table also gives the
cost relative to measuring the target alone and the correlation the
proxy must have before it can save anything.

## Usage

``` r
ss_aipe_two_phase_cost(
  true_sd_Y,
  true_cor,
  width,
  cost_ratio,
  conf_level = 0.95
)
```

## Arguments

- true_sd_Y:

  Planning value of the standard deviation of the target measurement,
  \\\sigma_Y\\.

- true_cor:

  Planning value of the correlation between the target and the proxy,
  \\\rho\\. Only its absolute value matters.

- width:

  The desired full width of the confidence interval (the upper limit
  minus the lower limit), in the units of the target.

- cost_ratio:

  The cost of one proxy measurement relative to one target measurement,
  \\r = c_P / c_T\\, a number in (0, 1).

- conf_level:

  Confidence level (default `0.95`).

## Value

A `data.frame` (classes `dmar_ss_aipe` and `dmar_tbl`) with rows

- necessary_n_labeled:

  the number of labeled units, \\n\\; when `necessary_n_pool` is 0, the
  number of units for the target alone.

- necessary_n_pool:

  the number of units in the pool, \\M\\, every one of which carries the
  proxy; 0 when the design measures the target alone.

- relative_cost:

  the cost of the returned design, \\n + rM\\, over the cost of the
  target alone, `necessary_n_target_only`, both in whole numbers under
  the criteria above; 1 when the returned design is the target alone.

- relative_cost_continuous:

  the same ratio for fractional sizes under the normal quantile,
  \\(\sqrt{1 - \rho^2} + \|\rho\|\sqrt{r})^2\\, or 1 at or below
  `threshold_cor`; reported for comparison.

- expected_width:

  the planned full width at the returned sizes and the planning values:
  the left side of the two-phase inequality above, or of the
  target-alone inequality when the pool is 0. It is at most `width`.

- necessary_n_target_only:

  the number of units the target alone needs for the same width.

- threshold_cor:

  the absolute correlation above which the proxy lowers the cost of the
  fractional solution, \\2\sqrt{r}/(1 + r)\\.

- width, true_sd_Y, true_cor, cost_ratio, conf_level:

  the inputs.

## Details

**Estimand and interval.** The quantity to be estimated is \\\mu_Y\\,
the population mean of the target measurement. The sampling unit is one
unit of the pool: \\M\\ units are drawn at random, the proxy is observed
on all of them, and the target on a simple random subsample of \\n\\ of
them. The interval planned for is the one
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
reports for the regression estimator of double sampling, a Wald-type
interval whose standard error estimates \\\sigma_Y^2(1 - \rho^2)/n +
\rho^2\sigma_Y^2/M\\ and whose reference distribution is *t* on \\n -
2\\ degrees of freedom. This is a plan for the expected width, with no
assurance: about half of all realized intervals are wider than `width`.
[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
adds an assurance for a given pool.

**The cost model.** Every unit of the pool, the labeled ones included,
costs \\c_P\\ for its proxy, and each labeled unit costs \\c_T\\ more
for its target, so a design costs \\c_T n + c_P M\\, which is \\n + rM\\
in units of one target measurement. In an evaluation of a system whose
outputs an automated judge labels, for example, \\c_T\\ is the cost of a
human assessor's label and \\c_P\\ that of the judge's. The target alone
on \\n\\ units costs \\n\\.

**Which designs are admissible.** A two-phase design has \\n \ge 3\\ and
\\M \ge n + 1\\, the least
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
accepts, and meets the width when, with \\q = 1 - (1 -
\mathrm{conf\\level})/2\\ and \\\omega\\ the desired width (`width`),
\$\$2\\ t\_{q,\\ n-2}\\ \sigma_Y \left\\\frac{1 - \rho^2}{n} +
\frac{\rho^2}{M}\right\\^{1/2} \le \omega,\$\$ the criterion
[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
uses for a given pool, so that planner, given
`n_pool = necessary_n_pool`, returns `necessary_n_labeled`. The left
side sets the two estimates in the standard error at their planning
values, which overstates the expected width slightly (see
[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)).
The alternative without a proxy is the *t* interval for one mean on the
smallest \\n\\ with \\2\\ t\_{q,\\ n-1}\\ c_4(n)\\ \sigma_Y / \sqrt{n}
\le \omega\\, where \\c_4(n)\sigma_Y\\ is the expected sample standard
deviation under normality; this is the `necessary_n_target_only` row
here and in
[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md).
Because the *t* quantile is part of the requirement, a design with very
few labeled units is admissible only when it is wide enough to pay for
its quantile (12.7 at a 95% level with 3 labeled units), and no separate
minimum is imposed.

**The search.** For each \\n\\ the smallest admissible pool follows from
the inequality above, so the search is over \\n\\ alone. Since
\\t\_{q,\\ n-2}\\ is at least the normal quantile \\z_q\\, no design
with \\n\\ labeled units can cost less than \\n + r \max\\n + 1,\\
b/(1 - a/n)\\\\ for \\n \> a\\, and none exists for \\n \le a\\, where
\\a = (1 - \rho^2) s\\ and \\b = \rho^2 s\\ with \\s = (2 z_q
\sigma_Y/\omega)^2\\. That bound is convex in \\n\\, so the function
starts at its minimizer with the target alone as the design to beat,
examines every \\n\\ above and then below the minimizer for as long as
the bound is less than the cost of the best design found, and returns
that design. No design in whole numbers costs less: every \\n\\ not
examined has a bound at least as large as the cost returned. A two-phase
design is returned only when it costs less than the target alone, and
among designs of equal cost the first found is kept, so a tie with the
target alone returns the target alone. The search and its bound were
derived for DMAR; in 2,200 planning sets drawn at random (correlations
of either sign from .02 to .999, cost ratios from .001 to .98, four
confidence levels) the result equaled an exhaustive search over pairs
\\(n, M)\\ every time.

**The fractional solution, for comparison.** If the sizes could be
fractional and the normal quantile were used, minimizing \\n + rM\\
subject to \\a/n + b/M = 1\\ (a Lagrange multiplier argument) gives \\n
= a + \sqrt{rab}\\ and \\M = b + \sqrt{ab/r}\\, at a cost of
\\(\sqrt{1 - \rho^2} + \|\rho\|\sqrt{r})^2\\ times the \\s\\ units the
target alone would need. That factor is the `relative_cost_continuous`
row. It is below 1 exactly when \\\|\rho\| \> 2\sqrt{r}/(1 + r)\\, the
`threshold_cor` row; at or below that correlation the row is 1. The two
sizes are in the ratio \\M/n = \sqrt{b/(ra)}\\: each is proportional to
its share of the variance and inversely proportional to the square root
of its cost. Neyman (1938, Equations 49, 52, and 53, pp. 110–111)
derived that square root rule for double sampling in which the
inexpensive variable forms strata and the variance is minimized for a
fixed total cost; the form here, with a regression on the proxy and the
cost minimized for a fixed width, is the same allocation problem stated
the other way round. For the regression estimator itself, Tamhane (1978,
Equation 2.2, p. 421) gives the same ratio of the two sizes for a fixed
budget, and the correlation above which double sampling is no worse than
the target alone, \\\rho^2 \ge 4 c_1 c_2/(c_1 + c_2)^2\\ with \\c_1\\
and \\c_2\\ the costs of the proxy and of the target, which is
`threshold_cor`; the result is for a bivariate normal population with a
known covariance matrix. The fractional solution is a lower bound: no
design in whole numbers that meets the width costs less than
`relative_cost_continuous` times \\s\\ target measurements. The design
returned costs more than that bound because sizes are whole numbers and
because the interval uses the *t* quantile (a median of 7% more, and far
more when only a few units are needed, in the 891 random planning sets
above in which a two-phase design was returned).

**What can happen at the edges.**

- A correlation at or below `threshold_cor`: the target alone is
  returned (`necessary_n_pool` is 0 and `relative_cost` is 1). In the
  random planning sets above no two-phase design was ever the least
  costly there.

- A correlation above `threshold_cor`: the threshold is exact for
  fractional sizes and the normal quantile only. The target alone can
  still be the least costly design in whole numbers, most often when few
  units are needed or the correlation is a little above the threshold
  (376 of the 1,267 random planning sets above the threshold), and it is
  then returned.

- `true_cor = 0`: the proxy carries no information and the target alone
  is returned.

- A proxy that costs as much as the target or more (`cost_ratio` of 1 or
  more) never saves anything and is refused with an error of class
  `dmar_cost_ratio_invalid`.

- A very wide `width`: the target alone needs at least 2 units, the
  fewest from which a standard deviation can be estimated.

- Sizes that cannot be represented as whole numbers (above \\2^{52}\\,
  about 4.5e15), or a search that would examine more than \\10^8\\
  values of \\n\\: refused with an error of class
  `dmar_two_phase_no_design`. Neither arises for a width within a factor
  of a million of `true_sd_Y`.

## References

Neyman, J. (1938). Contribution to the theory of sampling human
populations. *Journal of the American Statistical Association, 33*(201),
101–116.
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)

Tamhane, A. C. (1978). Inference based on regression estimator in double
sampling. *Biometrika, 65*(2), 419–427.
[doi:10.1093/biomet/65.2.419](https://doi.org/10.1093/biomet/65.2.419)

## See also

[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
for the labeled units given a pool, with or without an assurance,
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
for the interval, and
[`ss_aipe_two_phase_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)
to evaluate a design by simulation.

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
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
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
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A proxy measurement costs 1% of a target measurement. Planning values
# of 0.43 for the target's standard deviation and 0.55 for its
# correlation with the proxy; a 90% interval no wider than 0.10. The
# least costly design labels 151 units of a pool of 957 and costs 79%
# of what the 202 units of the target alone would.
ss_aipe_two_phase_cost(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
                       cost_ratio = 0.01, conf_level = 0.90)
#>  term                     value
#>  necessary_n_labeled      151  
#>  necessary_n_pool         957  
#>  relative_cost            0.795
#>  relative_cost_continuous 0.792
#>  expected_width           0.100
#>  necessary_n_target_only  202  
#>  threshold_cor            0.198
#>  width                    0.100
#>  true_sd_Y                0.430
#>  true_cor                 0.550
#>  cost_ratio               0.010
#>  conf_level               0.900
#> 
#> Confidence level: 90%

# At 10% of the cost the same proxy no longer pays: a correlation of 0.55
# is below the threshold, and the design measures the target alone.
ss_aipe_two_phase_cost(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
                       cost_ratio = 0.10, conf_level = 0.90)
#>  term                     value
#>  necessary_n_labeled      202  
#>  necessary_n_pool         0    
#>  relative_cost            1    
#>  relative_cost_continuous 1    
#>  expected_width           0.100
#>  necessary_n_target_only  202  
#>  threshold_cor            0.575
#>  width                    0.100
#>  true_sd_Y                0.430
#>  true_cor                 0.550
#>  cost_ratio               0.100
#>  conf_level               0.900
#> 
#> Confidence level: 90%

# A worked example in standard deviation units: a proxy that correlates
# 0.90 with the target at 5% of its cost, and a 95% interval of width
# 0.20. The design is 109 labeled units in a pool of 1011, a cost of
# 109 + 0.05 * 1011 = 159.55 target measurements against 387.
plan <- ss_aipe_two_phase_cost(true_sd_Y = 1, true_cor = 0.90, width = 0.20,
                               cost_ratio = 0.05)
plan
#>  term                     value
#>  necessary_n_labeled      109  
#>  necessary_n_pool         1011 
#>  relative_cost            0.412
#>  relative_cost_continuous 0.406
#>  expected_width           0.200
#>  necessary_n_target_only  387  
#>  threshold_cor            0.426
#>  width                    0.200
#>  true_sd_Y                1    
#>  true_cor                 0.900
#>  cost_ratio               0.050
#>  conf_level               0.950
#> 
#> Confidence level: 95%

# The companion planner, given that pool, returns the same number of
# labeled units.
ss_aipe_two_phase(true_sd_Y = 1, true_cor = 0.90, width = 0.20, n_pool = 1011)
#>  term                    value
#>  necessary_n_labeled     109  
#>  necessary_n_target_only 387  
#>  ratio                   0.282
#>  expected_width          0.200
#>  minimum_width           0.123
#>  width                   0.200
#>  true_sd_Y               1    
#>  true_cor                0.900
#>  n_pool                  1011 
#>  conf_level              0.950
#> 
#> Confidence level: 95%
```
