# Sensitivity Analysis for Sample Size Planning With the Goal of Accuracy in Parameter Estimation (I.e., a Narrow Observed Confidence Interval)

Given `estimated_R2` and `true_R2`, one can perform a sensitivity
analysis to determine the effect of a misspecified population squared
multiple correlation coefficient using the Accuracy in Parameter
Estimation (AIPE) approach to sample size planning. The function
evaluates the effect of a misspecified `true_R2` on the width of
obtained confidence intervals.

## Usage

``` r
ss_aipe_R2_sensitivity(
  true_R2 = NULL,
  estimated_R2 = NULL,
  w = NULL,
  p = NULL,
  random_predictors = TRUE,
  specified_N = NULL,
  assurance = NULL,
  conf_level = 0.95,
  generate_random_predictors = TRUE,
  rho_yx = 0.3,
  rho_xx = 0.3,
  G = 10000,
  print_iter = TRUE,
  save = FALSE,
  filename = "ss_aipe_r2_sensitivity_result.csv",
  ...
)
```

## Arguments

- true_R2:

  Value of the population squared multiple correlation coefficient

- estimated_R2:

  Value of the estimated (for sample size planning) squared multiple
  correlation coefficient

- w:

  Full confidence interval width of interest

- p:

  Number of predictors

- random_predictors:

  Whether or not the sample size procedure and the simulation itself
  should be based on random (set to `TRUE`) or fixed predictors (set to
  `FALSE`)

- specified_N:

  Selected sample size to use in order to determine distributional
  properties at a given value of sample size

- assurance:

  Parameter to ensure confidence interval width with a specified degree
  of certainty

- conf_level:

  Confidence interval coverage (symmetric coverage)

- generate_random_predictors:

  Specify whether the simulation should be based on random (default) or
  fixed regressors.

- rho_yx:

  Value of the correlation between *y* (dependent variable) and each of
  the *x* variables (independent variables)

- rho_xx:

  Value of the correlation among the *x* variables (independent
  variables)

- G:

  Number of generations (i.e., replications) of the simulation

- print_iter:

  Should the iteration number (between 1 and `G`) during the run of the
  function

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

- ...:

  for modifying parameters of functions this function calls upon

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. Each row corresponds
to one summary statistic; the `term` entries include the mean / median /
SD of the lower and upper confidence limits on \\R^2\\, the mean /
median / SD of the observed \\R^2\\, the mean / median / SD of the full
and one-sided interval widths, the percentage of intervals with width at
or below the planning target (`pct_less_w`), the tail-specific empirical
non-coverage of `true_R2` (`pct_ci_miss_low`, `pct_ci_miss_high`) and
the overall empirical type I error rate (`total_type_I_error`), and the
number of replications on which a confidence interval could not be
obtained (`num_probs_with_cis`).

## Details

When `estimated_R2`=`true_R2`, the results are that of a simulation
study when all assumptions are satisfied. Rather than specifying
`estimated_R2`, one can specify `specified_N` to determine the results
of a particular sample size (when doing this `estimated_R2` cannot be
specified).

The sample size estimation procedure technically assumes multivariate
normal variables (`p`+1) with fixed predictors (`x`/independent
variables), yet the function assumes random multivariate normal
predictors (having a `p`+1 multivariate distribution). As Gatsonis and
Sampson (1989) note in the context of statistical power analysis (recall
this function is used in the context of precision), there is little
difference in the outcome.

In the behavioral, educational, and social sciences, predictor variables
are almost always random, and thus `random_predictors` should generally
be used. `random_predictors=TRUE` specifies how both the sample size
planning procedure and the confidence intervals are calculated based on
the random predictors/regressors. The internal simulation generates
random or fixed predictors/regressors based on whether variables
predictor variables are random or fixed. However, when
`random_predictors=FALSE`, only the sample size planning procedure and
the confidence intervals are calculated based on the parameter. The
parameter `generate_random_predictors` (where the default is `TRUE` so
that random predictors/regressors are generated) allows random or fixed
predictor variables to be generated. Because the sample size planning
procedure and the internal simulation are both specified, for purposes
of sensitivity analysis random/fixed can be crossed to examine the
effects of specifying sample size based on one but using it on data
based on the other.

## References

Algina, J. & Olejnik, S. (2000). Determining sample size for accurate
estimation of the squared multiple correlation coefficient.
*Multivariate Behavioral Research, 35*, 119–137.
[doi:10.1207/s15327906mbr3501_5](https://doi.org/10.1207/s15327906mbr3501_5)

Gatsonis, C. & Sampson, A. R. (1989). Multiple Correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106*(3), 516–524.

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals, *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
interval estimation, power calculations, sample size estimation, and
hypothesis testing in multiple regression. *Behavior Research Methods,
Instruments, & Computers, 24*(4), 581–582.
[doi:10.3758/BF03203611](https://doi.org/10.3758/BF03203611)

## See also

[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Change 'G' to some large number (e.g., G=10,000)
ss_aipe_R2_sensitivity(true_R2 = .5, estimated_R2 = .4, w = .10, p = 5, conf_level = 0.95, G = 25)
#> 1 
#> 2 
#> 3 
#> 4 
#> 5 
#> 6 
#> 7 
#> 8 
#> 9 
#> 10 
#> 11 
#> 12 
#> 13 
#> 14 
#> 15 
#> 16 
#> 17 
#> 18 
#> 19 
#> 20 
#> 21 
#> 22 
#> 23 
#> 24 
#> 25 
#>  term                     value  
#>  mean_low_lim_r2          0.444  
#>  median_low_lim_r2        0.444  
#>  sd_low_lim_r2            0.0252 
#>  mean_up_lim_r2           0.538  
#>  median_up_lim_r2         0.538  
#>  sd_up_lim_r2             0.0231 
#>  mean_r2                  0.495  
#>  median_r2                0.494  
#>  sd_r2                    0.0242 
#>  mean_lower_ci_width_r2   0.0505 
#>  median_lower_ci_width_r2 0.0506 
#>  sd_lower_ci_width_r2     0.00109
#>  mean_upper_ci_width_r2   0.0433 
#>  median_upper_ci_width_r2 0.0433 
#>  sd_upper_ci_width_r2     0.00109
#>  mean_ci_width_r2         0.0937 
#>  median_ci_width_r2       0.0939 
#>  sd_ci_width_r2           0.00217
#>  pct_less_w               1      
#>  pct_ci_miss_low          0      
#>  pct_ci_miss_high         0.04   
#>  total_type_I_error       0.04   
#>  num_probs_with_cis       0      
#> 
#> Confidence level: 95%
```
