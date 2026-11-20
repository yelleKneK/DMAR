# Sensitivity Analysis for Sample Size Planning for the Standardized Mean From the Accuracy in Parameter Estimation (AIPE) Perspective

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation (AIPE) Perspective for the standardized
mean.

## Usage

``` r
ss_aipe_sm_sensitivity(
  true_sm = NULL,
  supposed_sm = NULL,
  width = NULL,
  specified_N = NULL,
  assurance = NULL,
  conf_level = 0.95,
  G = 10000,
  print_iter = TRUE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_sm:

  population standardized mean

- supposed_sm:

  supposed standardized mean, the planning value from which
  [`ss_aipe_sm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.md)
  plans the sample size (not an estimate computed from data)

- width:

  desired full width of the confidence interval for the population
  standardized mean

- specified_N:

  selected sample size to use in order to determine distributional
  properties of a given value of sample size

- assurance:

  parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be `NULL` or between zero and unity)

- conf_level:

  the desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- G:

  number of generations (i.e., replications) of the simulation

- print_iter:

  to print the current value of the iterations

- filename:

  Optional path to a CSV file; when supplied, the per-replication
  results (the observed standardized mean, the full and one-sided
  widths, the tail misses, and the interval limits) are written there,
  appended when the file already exists, and a throwaway run should
  point it at `tempfile(fileext = ".csv")`; the default `NULL` writes
  nothing.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across the `G` replications. The `term`
entries are: `mean_sm`, `median_sm`, `sd_sm` (summaries of the realized
standardized mean); `mean_ci_width`, `median_ci_width`, `sd_ci_width`
(summaries of the full interval widths); `mean_ci_width_lower` and
`mean_ci_width_upper` (mean one-sided widths, measured from the observed
standardized mean to each limit); `pct_ci_less_width` (proportion of
intervals at or below the target width); `pct_ci_miss_low` and
`pct_ci_miss_high` (tail-specific empirical non-coverage of `true_sm`);
`total_type_I_error` (overall empirical non-coverage, the sum of the two
tails); and the input echoes `total_N`, `true_sm`, `supposed_sm` (NA
when `specified_N` was supplied instead), `width`, `conf_level`, and
`assurance` (present only when an assurance was supplied). The
proportion and Type I error rows are proportions on the 0 to 1 scale,
not percentages.

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

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ss_aipe_sm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# With true_sm equal to supposed_sm the planning value is correct, and
# the sweep reports what a correctly specified plan delivers: the share of
# intervals no wider than the target should sit near the assurance. G = 10
# keeps the example quick; a reported sensitivity study deserves the
# default G = 10000.
set.seed(113)
res_correct <- ss_aipe_sm_sensitivity(true_sm = 1.5, supposed_sm = 1.5,
  width = 0.5, assurance = 0.95, conf_level = 0.95, G = 10,
  print_iter = FALSE)
res_correct
#>  term                value
#>  mean_sm             1.562
#>  median_sm           1.541
#>  sd_sm               0.112
#>  mean_ci_width       0.471
#>  median_ci_width     0.468
#>  sd_ci_width         0.019
#>  mean_ci_width_lower 0.237
#>  mean_ci_width_upper 0.235
#>  pct_ci_less_width   0.900
#>  pct_ci_miss_low     0.100
#>  pct_ci_miss_high    0.000
#>  total_type_I_error  0.100
#>  total_N             154  
#>  true_sm             1.500
#>  supposed_sm         1.500
#>  width               0.500
#>  conf_level          0.950
#>  assurance           0.950
#> 
#> Replications: 10
#> 
#> Confidence level: 95%

# The terms the summary reports.
res_correct$term
#>  [1] "mean_sm"             "median_sm"           "sd_sm"              
#>  [4] "mean_ci_width"       "median_ci_width"     "sd_ci_width"        
#>  [7] "mean_ci_width_lower" "mean_ci_width_upper" "pct_ci_less_width"  
#> [10] "pct_ci_miss_low"     "pct_ci_miss_high"    "total_type_I_error" 
#> [13] "total_N"             "true_sm"             "supposed_sm"        
#> [16] "width"               "conf_level"          "assurance"          

# The proportion of realized full widths no wider than the target.
res_correct$value[res_correct$term == "pct_ci_less_width"]
#> [1] 0.9

# The population standardized mean is 1.5 but the plan assumed 2, so the
# planner sizes the study for a wider sampling distribution than the data
# will show, and the realized intervals come in narrower than the target.
set.seed(113)
res_misspecified <- ss_aipe_sm_sensitivity(true_sm = 1.5, supposed_sm = 2,
  width = 0.5, G = 20, print_iter = FALSE)

# The effect of the misspecification on the mean interval width.
res_misspecified$value[res_misspecified$term == "mean_ci_width"]
#> [1] 0.4293829
```
