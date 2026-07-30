# Sensitivity Analysis for Sample Size Planing From the Accuracy in Parameter Estimation Perspective for the Unstandardized Regression Coefficient

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation Perspective for the unstandardized
regression coefficient.

## Usage

``` r
ss_aipe_rc_sensitivity(
  true_var_Y = NULL,
  true_cov_YX = NULL,
  true_cov_XX = NULL,
  estimated_var_Y = NULL,
  estimated_cov_YX = NULL,
  estimated_cov_XX = NULL,
  specified_N = NULL,
  which_predictor = 1,
  w = NULL,
  noncentral = FALSE,
  standardize = FALSE,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = TRUE,
  save = FALSE,
  filename = "ss_aipe_rc_sensitivity_result.csv"
)
```

## Arguments

- true_var_Y:

  Population variance of the dependent variable (*Y*)

- true_cov_YX:

  Population covariances vector between the *p* predictor variables and
  the dependent variable (*Y*)

- true_cov_XX:

  Population covariance matrix of the *p* predictor variables

- estimated_var_Y:

  Estimated variance of the dependent variable (*Y*)

- estimated_cov_YX:

  Estimated covariances vector between the *p* predictor variables and
  the dependent variable (*Y*)

- estimated_cov_XX:

  Estimated Population covariance matrix of the *p* predictor variables

- specified_N:

  Directly specified sample size (instead of using `Estimated.Rho.YX`
  and `Estimated.RHO.XX`)

- which_predictor:

  identifies which of the *p* predictors is of interest

- w:

  desired confidence interval width for the regression coefficient of
  interest

- noncentral:

  specify with a `TRUE` or `FALSE` statement whether or not the
  noncentral approach to sample size planning should be used

- standardize:

  specify with a `TRUE` or `FALSE` statement whether or not the
  regression coefficient will be standardized; default is `TRUE`

- conf_level:

  desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate)

- assurance:

  degree of certainty that the obtained confidence interval will be
  sufficiently narrow (i.e., the probability that the observed interval
  will be no larger than desired)

- G:

  the number of generations/replication of the simulation student within
  the function

- print_iter:

  specify with a `TRUE/FALSE` statement if the iteration number should
  be printed as the simulation within the function runs

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis. This function delegates to
[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
and inherits its return structure: each row reports one summary
statistic from the realized distribution across replications, including
means, medians, and standard deviations of the standardized regression
coefficient and its lower and upper confidence limits, the realized
interval widths, the percentage of intervals at or below the planning
target, and the empirical Type I error rates.

## Details

Direct specification of `True.Rho.YX` and `True.RHO.XX` is necessary,
even if one is interested in a single regression coefficient, so that
the covariance/correlation structure can be specified when the
simulation student within the function runs.

## Note

Note that when `True.Rho.YX=Estimated.Rho.YX` and
`True.RHO.XX=Estimated.RHO.XX`, the results are not literally from a
sensitivity analysis, rather the function performs a standard simulation
study. A simulation study can be helpful in order to determine if the
sample size procedure under or overestimates necessary sample size. See
`ss_aipe_reg_coef_sensitivity` in DMAR for more details.

## References

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

## See also

[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md),
[`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md),
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Sensitivity analysis for an unstandardized regression coefficient
# with two correlated predictors. G is kept small here so the example
# runs quickly; raise G for a stable Monte Carlo summary.
set.seed(113)
Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
rho_YX  <- c(0.4, 0.3)
cov_YX  <- rho_YX
ss_aipe_rc_sensitivity(
  true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
  estimated_var_Y = 1, estimated_cov_YX = cov_YX, estimated_cov_XX = Sigma_X,
  which_predictor = 1, w = 0.20, conf_level = 0.95,
  G = 20, print_iter = FALSE
)
#>  term               value  
#>  mean_b_j           0.337  
#>  median_b_j         0.332  
#>  sd_b_j             0.0428 
#>  mean_ci_width      0.195  
#>  median_ci_width    0.194  
#>  sd_ci_width        0.00901
#>  pct_ci_less_w      65     
#>  pct_ci_miss_low    0      
#>  pct_ci_miss_high   0      
#>  total_type_I_error 0      
#>  mean_r2            0.199  
#>  median_r2          0.195  
#>  sd_r2              0.0361 
#> 
#> Confidence level: 95%
# }
```
