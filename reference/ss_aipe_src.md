# Sample Size Necessary for the Accuracy in Parameter Estimation Approach for a Standardized Regression Coefficient of Interest

A function used to plan sample size from the accuracy in parameter
estimation approach for a standardized regression coefficient of
interest given the input specification.

## Usage

``` r
ss_aipe_src(
  rho2_Y_X = NULL,
  Rho2_j_X_without_j = NULL,
  p = NULL,
  beta_j = NULL,
  width,
  which_width = "Full",
  sigma_Y = 1,
  sigma_X_j = 1,
  rho_XX = NULL,
  rho_YX = NULL,
  which_predictor = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  conf_level = 0.95,
  assurance = NULL
)
```

## Arguments

- rho2_Y_X:

  Population value of the squared multiple correlation coefficient

- Rho2_j_X_without_j:

  Population value of the squared multiple correlation coefficient
  predicting the *j*th predictor variable from the remaining *p*-1
  predictor variables

- p:

  The number of predictor variables

- beta_j:

  The regression coefficient for the *j*th predictor variable (i.e., the
  predictor of interest)

- width:

  The desired width of the confidence interval

- which_width:

  Which width (`"Full"`, `"Lower"`, or `"Upper"`) the width refers to
  (at present, only `"Full"` can be specified)

- sigma_Y:

  The population standard deviation of *Y* (i.e., the dependent
  variables)

- sigma_X_j:

  The population standard deviation of the *j*th *X* variable (i.e., the
  predictor variable of interest)

- rho_XX:

  Population correlation matrix for the *p* predictor variables

- rho_YX:

  Population *p* length vector of correlation between the dependent
  variable (*Y*) and the *p* independent variables

- which_predictor:

  Identifies which of the *p* predictors is of interest

- alpha_lower:

  Type I error rate for the lower confidence interval limit

- alpha_upper:

  Type I error rate for the upper confidence interval limit

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate)

- assurance:

  Degree of certainty that the obtained confidence interval will be
  sufficiently narrow, which yields an approximate sample size to be
  verified with function `ss_aipe_reg_coef_sensitivity` to determine if
  it is appropriate

## Value

Returns the necessary sample size in order for the goals of accuracy in
parameter estimation to be satisfied for the confidence interval for a
particular regression coefficient given the input specifications.

## Details

Not all of the arguments need to be specified, only those that provide
all of the necessary information so that the sample size can be
determined for the conditions specified.

## Note

This function calls upon
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
in DMAR but has a different naming scheme. See
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
for more details.

## Warning

As discussed in Kelley and Maxwell (2008), the sample size planning
approach from the AIPE perspective used in this function is only an
approximation.

## References

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

## See also

[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Exchangable correlation structure
rho_YX <- c(.3, .3, .3, .3, .3)
rho_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5), c(.5, .5, 1, .5, .5),
                c(.5, .5, .5, 1, .5), c(.5, .5, .5, .5, 1))

ss_aipe_src(width = .1, which_width = "Full", sigma_Y = 1, sigma_X_j = 1, rho_XX = rho_XX,
            rho_YX = rho_YX, which_predictor = 1, conf_level = 1 - .05)
#>  term        value
#>  necessary_N 2191 
#> 
#> Confidence level: 95%

ss_aipe_src(width = .1, which_width = "Full", sigma_Y = 1, sigma_X_j = 1, rho_XX = rho_XX,
            rho_YX = rho_YX, which_predictor = 1, conf_level = 1 - .05,
            assurance = .85)
#>  term        value
#>  necessary_N 2241 
#> 
#> Confidence level: 95%
```
