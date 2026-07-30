# Sensitivity Analysis for Sample Size Planning for the Standardized ANOVA Contrast From the Accuracy in Parameter Estimation (AIPE) Perspective

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation (AIPE) Perspective for the standardized
ANOVA contrast.

## Usage

``` r
ss_aipe_sc_sensitivity(
  true_psi = NULL,
  estimated_psi = NULL,
  c_weights,
  desired_width = NULL,
  n_per_group = NULL,
  assurance = NULL,
  conf_level = 0.95,
  G = 10000,
  print_iter = TRUE,
  save = FALSE,
  filename = "ss_aipe_sc_sensitivity_result.csv",
  ...
)
```

## Arguments

- true_psi:

  population standardized contrast

- estimated_psi:

  estimated standardized contrast

- c_weights:

  the contrast weights

- desired_width:

  the desired full width of the obtained confidence interval

- n_per_group:

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

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

- ...:

  allows one to potentially include parameter values for inner functions

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `rep` replications. Inputs echoed back
as rows include `true_psi`, `est_psi`, `desired_width`, `assurance`,
`sample_size`, and `groups`. Realized-distribution rows include the mean
/ median / SD of the full interval width (`mean_full_width`,
`median_full_width`, `sd_full_width`), the mean lower- and upper-tail
widths (`mean_width_lower`, `mean_width_upper`), the proportion of
intervals at or below the target width (`pct_Width_obs_narrower`), and
the empirical Type I error rates on each tail (`type_I_error_upper`,
`type_I_error_lower`).

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
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Sensitivity analysis for a standardized three-group ANOVA contrast
# (-1, 0, 1) at psi = 0.5 and target full width 0.40. G is kept small
# here so the example runs quickly; raise it for a stable sweep.
set.seed(113)
ss_aipe_sc_sensitivity(
  true_psi = 0.5, estimated_psi = 0.5,
  c_weights = c(-1, 0, 1),
  desired_width = 0.40,
  conf_level = 0.95, G = 50, print_iter = FALSE
)
#>  term                   value  
#>  mean_full_width        0.399  
#>  median_full_width      0.399  
#>  sd_full_width          0.00164
#>  pct_Width_obs_narrower 0.72   
#>  mean_width_lower       0.2    
#>  mean_width_upper       0.199  
#>  type_I_error_upper     4      
#>  type_I_error_lower     2      
#> 
#> Confidence level: 95%
# }
```
