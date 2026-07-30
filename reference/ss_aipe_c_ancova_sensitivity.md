# Sensitivity Analysis for Sample Size Planning for the (Unstandardized) Contrast in Randomized ANCOVA From the Accuracy in Parameter Estimation (AIPE) Perspective

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation (AIPE) Perspective for the
(unstandardized) contrast in randomized ANCOVA design.

## Usage

``` r
ss_aipe_c_ancova_sensitivity(
  true_error_var_ancova = NULL,
  est_error_var_ancova = NULL,
  true_error_var_anova = NULL,
  est_error_var_anova = NULL,
  rho,
  est_rho = NULL,
  G = 10000,
  mu_y,
  sigma_y,
  mu_x,
  sigma_x,
  c_weights,
  width,
  conf_level = 0.95,
  assurance = NULL,
  save = FALSE,
  filename = "ss_aipe_c_ancova_sensitivity_result.csv"
)
```

## Arguments

- true_error_var_ancova:

  population error variance of the ANCOVA model

- est_error_var_ancova:

  estimated error variance of the ANCOVA model

- true_error_var_anova:

  population error variance of the ANOVA model (i.e., excluding the
  covariate)

- est_error_var_anova:

  estimated error variance of the ANOVA model (i.e., excluding the
  covariate)

- rho:

  population correlation coefficient of the response and the covariate

- est_rho:

  estimated correlation coefficient of the response and the covariate

- G:

  number of generations (i.e., replications) of the simulation

- mu_y:

  vector that contains the response's population mean of each group

- sigma_y:

  the population standard deviation of the response

- mu_x:

  the population mean of the covariate

- sigma_x:

  the population standard deviation of the covariate

- c_weights:

  the contrast weights

- width:

  the desired full width of the obtained confidence interval

- conf_level:

  the desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- assurance:

  parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be NULL or between zero and unity)

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis. Summary rows include the mean and median
realized full-CI width (`mean_width_obs`, `median_width_obs`), the
proportion of intervals narrower than the planning target
(`width_narrower`), the mean ratio of the ANCOVA standard error to the
ANOVA-only (`mean_ratio`), and tail-specific and overall empirical Type
I error rates (`type_I_error_upper`, `type_I_error_lower`,
`type_I_error`). The per-replication vectors (`psi_obs`, `se_psi`,
`se_psi_restricted`, `width_obs`) are not returned; they are written to
the CSV named by `filename` when `save = TRUE`.

## Details

The arguments `mu_y`, `mu_x`, `sigma_y`, and `sigma_x` are used to
generate random data in the simulations for the sensitivity analysis.
The value of `sigma_y` should be the same as the square root of
`true_error_var_anova`.

So far this function is based on one-covariate randomized ANCOVA design
only. The argument `mu_x` should be a single number, because it is
assumed that the population mean of the covariate is equal across groups
in randomized ANCOVA.

## References

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Monte Carlo sensitivity sweep; G is small here so the example runs quickly.
# Raise G (e.g., G = 1000 or more) for a stable sensitivity analysis.
set.seed(113)
ss_aipe_c_ancova_sensitivity(true_error_var_ancova=30,
                             est_error_var_ancova=30, rho=.2, mu_y=c(10,12,15,13), mu_x=2,
                             G=50, sigma_x=1.3, sigma_y=2, c_weights=c(1,0,-1,0), width=3)
#>  term               value
#>  type_I_error       6    
#>  type_I_error_upper 2    
#>  type_I_error_lower 4    
#>  width_narrower     100  
#>  mean_width_obs     1.07 
#>  median_width_obs   1.07 
#>  mean_ratio         0.999
#> 
#> Confidence level: 95%

ss_aipe_c_ancova_sensitivity(true_error_var_anova=36,
                             est_error_var_anova=36, rho=.2, est_rho=.2, G=50,
                             mu_y=c(10,12,15,13), mu_x=2, sigma_x=1.3, sigma_y=6,
                             c_weights=c(1,0,-1,0), width=3, assurance=NULL)
#>  term               value
#>  type_I_error       2    
#>  type_I_error_upper 2    
#>  type_I_error_lower 0    
#>  width_narrower     50   
#>  mean_width_obs     3.02 
#>  median_width_obs   3    
#>  mean_ratio         0.999
#> 
#> Confidence level: 95%
# }
```
