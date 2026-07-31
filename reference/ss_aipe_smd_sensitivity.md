# Sensitivity Analysis for Sample Size Given the Accuracy in Parameter Estimation Approach for the Standardized Mean Difference

Performs sensitivity analysis for sample size determination for the
standardized mean difference given a population and a standardized mean
difference. Allows one to determine the effect of being wrong when
estimating the population standardized mean difference in terms of the
width of the obtained (two-sided) confidence intervals.

## Usage

``` r
ss_aipe_smd_sensitivity(
  true_delta = NULL,
  estimated_delta = NULL,
  desired_width = NULL,
  n_per_group = NULL,
  assurance = NULL,
  conf_level = 0.95,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_smd_sensitivity_result.csv",
  ...
)
```

## Arguments

- true_delta:

  population standardized mean difference

- estimated_delta:

  estimated standardized mean difference; can be `true_delta` to perform
  standard simulations

- desired_width:

  describe full width for the confidence interval around the population
  standardized mean difference

- n_per_group:

  selected sample size to use in order to determine distributional
  properties of at a given value of sample size

- assurance:

  parameter to ensure confidence interval width with a specified degree
  of certainty (must be `NULL` or between zero and unity)

- conf_level:

  the desired degree of confidence (i.e., 1-Type I error rate)

- G:

  number of generations (i.e., replications) of the simulation

- print_iter:

  to print the current value of the iterations

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

- ...:

  for modifying parameters of functions this function calls

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis. Summary rows include the mean / median / SD
of the realized full confidence interval width (`mean_full_width`,
`median_full_width`, `sd_full_width`), the proportion of intervals
narrower than the planning target (`pct_less_desired`), the mean lower-
and upper-tail widths (`mean_width_lower`, `mean_width_upper`), and the
empirical Type I error rates on each tail (`type_I_error_upper`,
`type_I_error_lower`).

## Details

For sensitivity analysis when planning sample size given the desire to
obtain narrow confidence intervals for the population standardized mean
difference. Given a population value and an estimated value, one can
determine the effects of incorrectly specifying the population
standardized mean difference (`true_delta`) on the obtained widths of
the confidence intervals. Also, one can evaluate the percent of the
confidence intervals that are less than the desired width (especially
when modifying the `assurance` parameter); see `ss_aipe_smd`)
Alternatively, one can specify `n_per_group` to determine the results at
a particular sample size (when doing this `estimated_delta` cannot be
specified).

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
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

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
# Since 'true_delta' equals 'estimated_delta', this usage
# returns the results of a correctly specified situation.
# Note that 'G' should be large (50 is used to make the example run easily)
Res.1 <- ss_aipe_smd_sensitivity(true_delta=.5, estimated_delta=.5, desired_width=.30,
                                 assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)

# Objects contained in the 'summary'.
Res.1$term
#> [1] "mean_full_width"    "median_full_width"  "sd_full_width"     
#> [4] "pct_less_desired"   "mean_width_lower"   "mean_width_upper"  
#> [7] "type_I_error_upper" "type_I_error_lower"

# True standardized mean difference is .4, but specified at .5.
# Change 'G' to some large number (e.g., G=5,000)
Res.2 <- ss_aipe_smd_sensitivity(true_delta=.4, estimated_delta=.5, desired_width=.30,
                                 assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)

# The effect of the misspecification on mean confidence intervals is:
Res.2[1,]
#>  term            value
#>  mean_full_width 0.298
#> 
#> Confidence level: 95%

# True standardized mean difference is .5, but specified at .4.
Res.3 <- ss_aipe_smd_sensitivity(true_delta=.5, estimated_delta=.4, desired_width=.30,
                                 assurance=NULL, conf_level=.95, G=50, print_iter=FALSE)

# The effect of the misspecification on mean confidence intervals is:
Res.3[1,]
#>  term            value
#>  mean_full_width 0.302
#> 
#> Confidence level: 95%
```
