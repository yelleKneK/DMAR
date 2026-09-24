# Sensitivity Analysis for Sample Size Planning for the Standardized Mean From the Accuracy in Parameter Estimation (AIPE) Perspective

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation (AIPE) Perspective for the standardized
mean.

## Usage

``` r
ss_aipe_sm_sensitivity(
  true_sm = NULL,
  estimated_sm = NULL,
  desired_width = NULL,
  specified_N = NULL,
  assurance = NULL,
  conf_level = 0.95,
  G = 10000,
  print_iter = TRUE,
  filename = NULL
)
```

## Arguments

- true_sm:

  population standardized mean

- estimated_sm:

  estimated standardized mean

- desired_width:

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

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across the `G` replications. The `term`
entries are: `mean_sm`, `median_sm`, `sd_sm` (summaries of the realized
standardized mean); `mean_ci_width`, `median_ci_width`, `sd_ci_width`
(summaries of the full interval widths); `mean_ci_width_lower` and
`mean_ci_width_upper` (mean one-sided widths, measured from the observed
standardized mean to each limit); `pct_ci_less_w` (proportion of
intervals at or below the target width); `pct_ci_miss_low` and
`pct_ci_miss_high` (tail-specific empirical non-coverage of `true_sm`);
`total_type_I_error` (overall empirical non-coverage, the sum of the two
tails); and the input echoes `total_N`, `true_sm`, `estimated_sm` (NA
when `specified_N` was supplied instead), `width`, `conf_level`, and
`assurance` (present only when an assurance was supplied). The
proportion and Type I error rows are proportions on the 0 to 1 scale,
not percentages.

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
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

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
# With true_sm equal to estimated_sm the planning value is correct, and
# the sweep reports what a correctly specified plan delivers: the share of
# intervals no wider than the target should sit near the assurance. G = 10
# keeps the example quick; a reported sensitivity study deserves the
# default G = 10000.
set.seed(113)
res_correct <- ss_aipe_sm_sensitivity(true_sm = 1.5, estimated_sm = 1.5,
  desired_width = 0.5, assurance = 0.95, conf_level = 0.95, G = 10,
  print_iter = FALSE)
res_correct
#>  term                value 
#>  mean_sm             1.56  
#>  median_sm           1.54  
#>  sd_sm               0.112 
#>  mean_ci_width       0.471 
#>  median_ci_width     0.468 
#>  sd_ci_width         0.0188
#>  mean_ci_width_lower 0.237 
#>  mean_ci_width_upper 0.235 
#>  pct_ci_less_w       0.9   
#>  pct_ci_miss_low     0.1   
#>  pct_ci_miss_high    0     
#>  total_type_I_error  0.1   
#>  total_N             154   
#>  true_sm             1.5   
#>  estimated_sm        1.5   
#>  width               0.5   
#>  conf_level          0.95  
#>  assurance           0.95  
#> 
#> Confidence level: 95%

# The terms the summary reports.
res_correct$term
#>  [1] "mean_sm"             "median_sm"           "sd_sm"              
#>  [4] "mean_ci_width"       "median_ci_width"     "sd_ci_width"        
#>  [7] "mean_ci_width_lower" "mean_ci_width_upper" "pct_ci_less_w"      
#> [10] "pct_ci_miss_low"     "pct_ci_miss_high"    "total_type_I_error" 
#> [13] "total_N"             "true_sm"             "estimated_sm"       
#> [16] "width"               "conf_level"          "assurance"          

# The proportion of realized full widths no wider than the target.
res_correct$value[res_correct$term == "pct_ci_less_w"]
#> [1] 0.9

# The population standardized mean is 1.5 but the plan assumed 2, so the
# planner sizes the study for a wider sampling distribution than the data
# will show, and the realized intervals come in narrower than the target.
set.seed(113)
res_misspecified <- ss_aipe_sm_sensitivity(true_sm = 1.5, estimated_sm = 2,
  desired_width = 0.5, G = 20, print_iter = FALSE)

# The effect of the misspecification on the mean interval width.
res_misspecified$value[res_misspecified$term == "mean_ci_width"]
#> [1] 0.4293829
```
