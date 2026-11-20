# Sensitivity Analysis for Sample Size Planning for the Standardized ANOVA Contrast From the Accuracy in Parameter Estimation (AIPE) Perspective

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation (AIPE) Perspective for the standardized
ANOVA contrast.

## Usage

``` r
ss_aipe_sc_sensitivity(
  true_std_contrast = NULL,
  supposed_std_contrast = NULL,
  c_weights,
  width = NULL,
  specified_n_per_group = NULL,
  assurance = NULL,
  conf_level = 0.95,
  G = 10000,
  print_iter = TRUE,
  filename = NULL,
  assurance_method = "analytic",
  seed = NULL
)
```

## Arguments

- true_std_contrast:

  population standardized contrast

- supposed_std_contrast:

  supposed standardized contrast, the planning value given to
  [`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md)
  to choose the sample size

- c_weights:

  the contrast weights, fractions whose positive values sum to 1 and
  whose sum is zero, the rule of
  [`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md) applied
  to the values as
  [`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md)
  applies it: weights
  [`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md)
  refuses are refused here, with a condition of class
  `dmar_argument_invalid`, and integer weights such as `-1:1` are
  accepted

- width:

  the desired full width of the obtained confidence interval

- specified_n_per_group:

  selected sample size to use in order to determine distributional
  properties of at a given value of sample size

- assurance:

  parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be NULL or between zero and unity)

- conf_level:

  the desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- G:

  number of generations (i.e., replications) of the simulation

- print_iter:

  to print the current value of the iterations

- filename:

  an optional path for a comma separated file recording every
  replication (the realized standardized contrast, `std_contrast`, the
  full and the two one-sided interval widths, the three non-coverage
  indicators, and the two confidence limits): nothing is written when
  `filename` is `NULL` (the default), a new file with a header row is
  created otherwise, an existing file at that path is appended to, and a
  throwaway run should point it at `tempfile(fileext = ".csv")`

- assurance_method:

  passed to
  [`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md)
  when the sample size is planned from `supposed_std_contrast` with
  `assurance`. Only `"analytic"` (the default) is offered: the plan is
  the smallest sample size per group, at or above the plan without
  `assurance`, at which the exact probability of an interval no wider
  than `width` is at least `assurance`. `"simulation"` is refused with a
  condition of class `dmar_argument_invalid`, as it is by
  [`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
  since that probability is exact, whether or not `assurance` or
  `specified_n_per_group` is given.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across the `G` replications. The `term`
entries are: `mean_std_contrast`, `median_std_contrast`,
`sd_std_contrast` (summaries of the realized standardized contrast);
`mean_ci_width`, `median_ci_width`, `sd_ci_width` (summaries of the full
interval widths); `mean_ci_width_lower` and `mean_ci_width_upper` (mean
one-sided widths, measured from the observed contrast to each limit);
`pct_ci_less_width` (proportion of intervals at or below the target
width); `pct_ci_miss_low` and `pct_ci_miss_high` (tail-specific
empirical non-coverage of `true_std_contrast`); `total_type_I_error`
(overall empirical non-coverage, the sum of the two tails); and the
input echoes `n_per_group` (the size simulated at, planned or
specified), `total_N`, `true_std_contrast`, `supposed_std_contrast` (NA
when `specified_n_per_group` was supplied instead), `width`,
`conf_level`, and `assurance` (present only when an assurance was
supplied). The proportion and Type I error rows are proportions on the 0
to 1 scale, not percentages.

## Details

Each replication draws *J* groups of *n* normal observations with
standard deviation 1 and forms the
[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md) interval
from the group means and the pooled standard deviation, the square root
of the mean square error that
[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md) takes as
`s_anova`.

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

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574.
[doi:10.1177/0013164401614002](https://doi.org/10.1177/0013164401614002)

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in Parameter Estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Sensitivity analysis for a standardized three-group ANOVA contrast
# (-1, 0, 1) with a standardized contrast of 0.5 and target full width 0.40. G is kept small
# here so the example runs quickly; raise it for a stable sweep.
set.seed(113)
ss_aipe_sc_sensitivity(
  true_std_contrast = 0.5, supposed_std_contrast = 0.5,
  c_weights = c(-1, 0, 1),
  width = 0.40,
  conf_level = 0.95, G = 50, print_iter = FALSE
)
#>  term                  value
#>  mean_std_contrast     0.495
#>  median_std_contrast   0.493
#>  sd_std_contrast       0.103
#>  mean_ci_width         0.399
#>  median_ci_width       0.399
#>  sd_ci_width           0.002
#>  mean_ci_width_lower   0.200
#>  mean_ci_width_upper   0.199
#>  pct_ci_less_width     0.720
#>  pct_ci_miss_low       0.020
#>  pct_ci_miss_high      0.040
#>  total_type_I_error    0.060
#>  n_per_group           197  
#>  total_N               591  
#>  true_std_contrast     0.500
#>  supposed_std_contrast 0.500
#>  width                 0.400
#>  conf_level            0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
