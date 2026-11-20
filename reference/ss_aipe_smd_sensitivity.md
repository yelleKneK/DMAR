# Sensitivity Analysis for Sample Size Given the Accuracy in Parameter Estimation Approach for the Standardized Mean Difference

Performs sensitivity analysis for sample size determination for the
standardized mean difference given a population and a standardized mean
difference. Allows one to determine the effect of being wrong when
estimating the population standardized mean difference in terms of the
width of the obtained (two-sided) confidence intervals.

## Usage

``` r
ss_aipe_smd_sensitivity(
  true_smd = NULL,
  supposed_smd = NULL,
  width = NULL,
  specified_n_per_group = NULL,
  assurance = NULL,
  assurance_method = "analytic",
  conf_level = 0.95,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_smd:

  population standardized mean difference

- supposed_smd:

  the supposed standardized mean difference: the planning value from
  which
  [`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
  plans the per-group sample size (not an estimate computed from data);
  it can equal `true_smd` to simulate a correctly specified plan

- width:

  describe full width for the confidence interval around the population
  standardized mean difference

- specified_n_per_group:

  selected sample size to use in order to determine distributional
  properties of at a given value of sample size

- assurance:

  parameter to ensure confidence interval width with a specified degree
  of certainty (must be `NULL` or between zero and unity). When the
  sample size is planned from `supposed_smd`, an assurance above
  .999999998 is refused, as in
  [`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md).

- assurance_method:

  passed to
  [`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
  with `assurance` when the sample size is planned from `supposed_smd`.
  Only `"analytic"` (the default) is offered: the plan is the smallest
  per-group sample size, at or above the plan without `assurance`, at
  which the exact probability that the interval is no wider than `width`
  is at least `assurance`. `"simulation"` is refused with a condition of
  class `dmar_argument_invalid`, as it is by
  [`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
  since that probability is exact, whether or not `assurance` or
  `specified_n_per_group` is given.

- conf_level:

  the desired degree of confidence (i.e., 1-Type I error rate)

- G:

  number of generations (i.e., replications) of the simulation

- print_iter:

  to print the current value of the iterations

- filename:

  an optional path for a comma separated file recording every
  replication (the realized standardized mean difference, the full and
  the two one-sided interval widths, the three non-coverage indicators,
  and the two confidence limits): nothing is written when `filename` is
  `NULL` (the default), a new file with a header row is created
  otherwise, an existing file at that path is appended to, and a
  throwaway run should point it at `tempfile(fileext = ".csv")`

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across the `G` replications. The `term`
entries are: `mean_smd`, `median_smd`, `sd_smd` (summaries of the
realized standardized mean difference); `mean_ci_width`,
`median_ci_width`, `sd_ci_width` (summaries of the full interval
widths); `mean_ci_width_lower` and `mean_ci_width_upper` (mean one-sided
widths, measured from the observed standardized mean difference to each
limit); `pct_ci_less_width` (proportion of intervals at or below the
target width); `pct_ci_miss_low` and `pct_ci_miss_high` (tail-specific
empirical non-coverage of `true_smd`); `total_type_I_error` (overall
empirical non-coverage, the sum of the two tails); and the input echoes
`n_per_group` (the size simulated at, planned or specified), `total_N`,
`true_smd`, `supposed_smd` (NA when `specified_n_per_group` was supplied
instead), `width`, `conf_level`, and `assurance` (present only when an
assurance was supplied). The proportion and Type I error rows are
proportions on the 0 to 1 scale, not percentages.

## Details

For sensitivity analysis when planning sample size given the desire to
obtain narrow confidence intervals for the population standardized mean
difference. Given a population value and a supposed (planning) value,
one can determine the effects of incorrectly specifying the population
standardized mean difference (`true_smd`) on the obtained widths of the
confidence intervals. Also, one can evaluate the percent of the
confidence intervals that are less than the desired width (especially
when modifying the `assurance` parameter); see `ss_aipe_smd`)
Alternatively, one can specify `specified_n_per_group` to determine the
results at a particular sample size (when doing this `supposed_smd`
cannot be specified).

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

Kelley, K. (2005). The effects of nonnormal distributions on confidence
intervals around the standardized mean difference: Bootstrap and
parametric confidence intervals, *Educational and Psychological
Measurement, 65*, 51–69.
[doi:10.1177/0013164404264850](https://doi.org/10.1177/0013164404264850)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Since 'true_smd' equals 'supposed_smd', this usage
# returns the results of a correctly specified situation.
# Note that 'G' should be large (50 is used to make the example run easily)
set.seed(113)
Res.1 <- ss_aipe_smd_sensitivity(true_smd=.5, supposed_smd=.5, width=.30,
                                 assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)

# Objects contained in the 'summary'.
Res.1$term
#>  [1] "mean_smd"            "median_smd"          "sd_smd"             
#>  [4] "mean_ci_width"       "median_ci_width"     "sd_ci_width"        
#>  [7] "mean_ci_width_lower" "mean_ci_width_upper" "pct_ci_less_width"  
#> [10] "pct_ci_miss_low"     "pct_ci_miss_high"    "total_type_I_error" 
#> [13] "n_per_group"         "total_N"             "true_smd"           
#> [16] "supposed_smd"        "width"               "conf_level"         

# True standardized mean difference is .4, but specified at .5.
# Change 'G' to some large number (e.g., G=5,000)
Res.2 <- ss_aipe_smd_sensitivity(true_smd=.4, supposed_smd=.5, width=.30,
                                 assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)

# The effect of the misspecification on mean confidence intervals is:
Res.2[1,]
#>  term     value
#>  mean_smd 0.372
#> 
#> Replications: 50
#> 
#> Confidence level: 95%

# True standardized mean difference is .5, but specified at .4.
Res.3 <- ss_aipe_smd_sensitivity(true_smd=.5, supposed_smd=.4, width=.30,
                                 assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)

# The effect of the misspecification on mean confidence intervals is:
Res.3[1,]
#>  term     value
#>  mean_smd 0.484
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
