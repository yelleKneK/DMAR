# Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the Standardized Regression Coefficient

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation Perspective for the standardized
regression coefficient.

## Usage

``` r
ss_aipe_src_sensitivity(
  true_var_Y = NULL,
  true_cov_YX = NULL,
  true_cov_XX = NULL,
  estimated_var_Y = NULL,
  estimated_cov_YX = NULL,
  estimated_cov_XX = NULL,
  specified_N = NULL,
  which_predictor = 1,
  w = NULL,
  noncentral = TRUE,
  standardize = TRUE,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = TRUE,
  save = FALSE,
  filename = "ss_aipe_src_sensitivity_result.csv"
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

  Directly specified sample size (instead of planning one from the
  estimated covariance structure)

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
  sufficiently narrow

- G:

  the number of generations/replication of the simulation study within
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
and inherits its return structure: mean / median / SD summaries of the
realized standardized regression coefficient, the realized interval
widths, and the realized squared multiple correlation coefficient; the
proportion of intervals at or below the planning target
(`pct_ci_less_w`); the tail-specific and overall empirical non-coverage
rates (`pct_ci_miss_low`, `pct_ci_miss_high`, `total_type_I_error`), all
proportions on the 0 to 1 scale; and the input echoes (`total_N`, `p`,
`which_predictor`, `true_b_j`, `estimated_b_j`, `width`, `conf_level`,
and, when one was supplied, `assurance`). See
[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
for the full row list.

## Details

Direct specification of `true_cov_YX` and `true_cov_XX` is necessary,
even if one is interested in a single regression coefficient, so that
the covariance/correlation structure can be specified when the
simulation study within the function runs.

## Note

Note that when the true and estimated covariance structures agree
(`true_cov_YX` equals `estimated_cov_YX` and `true_cov_XX` equals
`estimated_cov_XX`), the results are not literally from a sensitivity
analysis, rather the function performs a standard simulation study. A
simulation study can be helpful in order to determine if the sample size
procedure under or overestimates necessary sample size. See
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
[`ss_aipe_rc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc_sensitivity.md),
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
# Sensitivity analysis for a standardized regression coefficient
# with two correlated predictors. A production run uses many more
# generations (G = 1000 is typical); G is reduced here so the
# example runs quickly.
set.seed(113)
Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
cov_YX  <- c(0.4, 0.3)
ss_aipe_src_sensitivity(
  true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
  estimated_var_Y = 1, estimated_cov_YX = cov_YX, estimated_cov_XX = Sigma_X,
  which_predictor = 1, w = 0.20, conf_level = 0.95,
  G = 50, print_iter = FALSE
)
#>  term               value  
#>  mean_b_j           0.349  
#>  median_b_j         0.357  
#>  sd_b_j             0.043  
#>  mean_ci_width      0.2    
#>  median_ci_width    0.2    
#>  sd_ci_width        0.00376
#>  pct_ci_less_w      0.56   
#>  pct_ci_miss_low    0      
#>  pct_ci_miss_high   0      
#>  total_type_I_error 0      
#>  mean_R2            0.201  
#>  median_R2          0.198  
#>  sd_R2              0.0366 
#>  total_N            366    
#>  p                  2      
#>  which_predictor    1      
#>  true_b_j           0.341  
#>  estimated_b_j      0.341  
#>  width              0.2    
#>  conf_level         0.95   
#> 
#> Confidence level: 95%
```
