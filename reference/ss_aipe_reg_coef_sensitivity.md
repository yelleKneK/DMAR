# Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the (Standardized and Unstandardized) Regression Coefficient

This function performs a sensitivity analysis when planning sample size
from the Accuracy in Parameter Estimation Perspective for the
standardized or unstandardized regression coefficient.

## Usage

``` r
ss_aipe_reg_coef_sensitivity(
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
  filename = "ss_aipe_reg_coef_sensitivity_result.csv"
)
```

## Arguments

- true_var_Y:

  Population variance of the dependent variable (*Y*)

- true_cov_YX:

  Population covariances vector between the `p` predictor variables and
  the dependent variable (*Y*)

- true_cov_XX:

  Population covariance matrix of the `p` predictor variables

- estimated_var_Y:

  Estimated variance of the dependent variable (*Y*)

- estimated_cov_YX:

  Estimated covariances vector between the `p` predictor variables and
  the dependent variable (`Y`)

- estimated_cov_XX:

  Estimated Population covariance matrix of the `p` predictor variables

- specified_N:

  Directly specified sample size (instead of using `Estimated.Rho.YX`
  and `Estimated.RHO.XX`)

- which_predictor:

  Identifies which of the *p* predictors is of interest

- w:

  desired Confidence interval width for the regression coefficient of
  interest

- noncentral:

  Specify with a `TRUE` or `FALSE` statement whether or not the
  noncentral approach to sample size planning should be used

- standardize:

  Specify with a `TRUE` or `FALSE` statement whether or not the
  regression coefficient will be standardized

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate)

- assurance:

  Degree of certainty that the obtained confidence interval will be
  sufficiently narrow

- G:

  The number of generations/replication of the simulation student within
  the function

- print_iter:

  Specify with a `TRUE`/`FALSE` statement if the iteration number should
  be printed as the simulation within the function runts

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis. Each row reports one summary statistic from
the realized distribution across replications, including means, medians,
and standard deviations of the regression-coefficient point estimates
and their lower and upper confidence limits, the mean / median / SD of
the realized interval widths, the percentage of intervals at or below
the planning target (`pct_less_w`), and the empirical Type I error
rates.

## Details

Direct specification of `True.Rho.YX` and `True.RHO.XX` is necessary,
even if one is interested in a single regression coefficient, so that
the covariance/correlation structure can be specified when the
simulation student within the function runs.

## Note

Note that when `True.Rho.YX`=`Estimated.Rho.YX` and
`True.RHO.XX`=`Estimated.RHO.XX`, the results are not literally from a
sensitivity analysis, rather the function performs a standard simulation
study. A simulation study can be helpful in order to determine if the
sample size procedure under or overestimates necessary sample size.

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
# with two predictors at a modest R squared. The Monte Carlo loop is
# run with a small number of generations (G) here so the example is
# fast; use a larger G (for example G = 1000) in real applications.
set.seed(113)
Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
cov_YX <- c(0.4, 0.3)
ss_aipe_reg_coef_sensitivity(
  true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = Sigma_X,
  estimated_var_Y = 1, estimated_cov_YX = cov_YX, estimated_cov_XX = Sigma_X,
  which_predictor = 1, w = 0.20, conf_level = 0.95,
  G = 100, print_iter = FALSE
)
#>  term               value 
#>  mean_b_j           0.338 
#>  median_b_j         0.339 
#>  sd_b_j             0.0528
#>  mean_ci_width      0.199 
#>  median_ci_width    0.199 
#>  sd_ci_width        0.0109
#>  pct_ci_less_w      52    
#>  pct_ci_miss_low    2     
#>  pct_ci_miss_high   4     
#>  total_type_I_error 6     
#>  mean_r2            0.198 
#>  median_r2          0.194 
#>  sd_r2              0.0381
#> 
#> Confidence level: 95%
# }
```
