# Sample Size Necessary for the Accuracy in Parameter Estimation Approach for an Unstandardized Regression Coefficient of Interest

A function used to plan sample size from the accuracy in parameter
estimation perspective for an unstandardized regression coefficient of
interest given the input specification.

## Usage

``` r
ss_aipe_rc(
  true_R2_Y_X = NULL,
  true_R2_j_X_without_j = NULL,
  p = NULL,
  true_reg_coef_j = NULL,
  width,
  which_width = "Full",
  true_sd_Y = 1,
  true_sd_X = 1,
  true_cor_XX = NULL,
  true_cor_YX = NULL,
  which_predictor = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  conf_level = 0.95,
  assurance = NULL,
  assurance_method = "analytic"
)
```

## Arguments

- true_R2_Y_X:

  Population value of the squared multiple correlation coefficient

- true_R2_j_X_without_j:

  Population value of the squared multiple correlation coefficient
  predicting the *j*th predictor variable from the remaining *p*-1
  predictor variables

- p:

  The number of predictor variables

- true_reg_coef_j:

  The population regression coefficient for the *j*th predictor variable
  (i.e., the predictor of interest), in the units `true_sd_Y` and
  `true_sd_X` describe: the standardized coefficient times
  `true_sd_Y / true_sd_X`. On the covariance structure path it may be
  left out; when it is given it must equal that value (see
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md))

- width:

  The desired width of the confidence interval

- which_width:

  Which Width (`"Full"`, `"Lower"`, or `"Upper"`) the width refers to
  (at present, only `"Full"` can be specified)

- true_sd_Y:

  The population standard deviation of *Y* (i.e., the dependent
  variable)

- true_sd_X:

  The population standard deviation of the *j*th *X* variable (i.e., the
  predictor variable of interest)

- true_cor_XX:

  Population correlation matrix for the *p* predictor variables

- true_cor_YX:

  Population *p* length vector of correlation between the dependent
  variable (*Y*) and the *p* independent variables

- which_predictor:

  Identifies which of the *p* predictors is of interest

- alpha_lower:

  Type I error rate for the lower confidence interval limit, at least
  1e-9 and at most .5; a tail area below 1e-9 is refused, as in
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)

- alpha_upper:

  Type I error rate for the upper confidence interval limit, at least
  1e-9 and at most .5

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate); a level above .999999998, whose tail areas are
  below 1e-9, is refused

- assurance:

  Optional probability with which the realized confidence interval is to
  be no wider than `width`; when `NULL` (the default) the plan is for
  the width of the interval at the expected values of the sample squared
  multiple correlations (see
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)),
  not for the mean width over samples

- assurance_method:

  How the probability that the realized interval is no wider than
  `width` is found, used only with `assurance`. Only `"analytic"` (the
  default) is offered here: the probability is computed exactly (an *F*
  probability) and reported as `assurance_attained`, and the plan is the
  smallest sample size, at or above the plan without `assurance`, at
  which it is at least `assurance`. `"simulation"`, which
  [`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md),
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
  with `noncentral = TRUE`, and
  [`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)
  offer, where the probability has no exact form, is refused with a
  condition of class `dmar_argument_invalid`. See "The assurance step"
  in
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md).

## Value

A `data.frame` with columns `term` and `value`: `necessary_N`, the
necessary total sample size; with `assurance`, the rows `assurance` and
`assurance_attained`, the exact probability at `necessary_N` that the
realized interval is no wider than `width`. With `assurance` the printed
footer names the assurance method, which the attribute
`assurance_method` records, with `assurance_probability` `"exact"`.

## Details

Not all of the arguments need to be specified, only those that provide
all of the necessary information so that the sample size can be
determined for the conditions specified.

The interval planned for is the central *t* interval for the coefficient
in the units of the data, the interval
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md) reports,
for random predictors from a multivariate normal population, with the
*t* quantiles of its two tails. The plan without `assurance` is the
smallest sample size, of at least \\p + 2\\, at which that interval, at
the expected squared multiple correlations, is no wider than `width`; it
can be a case or more smaller than that of MBESS (see
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)).
With `assurance` the plan is the smallest sample size at which that
interval is no wider than `width` with at least that probability,
computed exactly (`assurance_method = "analytic"`). The published chi
square approximation of Kelley and Maxwell (2003) can give a smaller or
a larger sample size than the exact plan.

## Note

This function calls upon `ss_aipe_reg_coef` in DMAR but has a different
naming scheme. See `ss_aipe_reg_coef` for more details.

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
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Exchangeable correlation structure
true_cor_YX <- c(.3, .3, .3, .3, .3)
true_cor_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5), c(.5, .5, 1, .5, .5),
c(.5, .5, .5, 1, .5), c(.5, .5, .5, .5, 1))

ss_aipe_rc(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
           true_cor_XX = true_cor_XX,
           true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05)
#>  term        value
#>  necessary_N 2184 
#> 
#> Confidence level: 95%

# The smallest N at which the interval is no wider than .1 with
# probability at least .85; the last row is that probability.
ss_aipe_rc(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
           true_cor_XX = true_cor_XX,
           true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05, assurance = .85)
#>  term               value
#>  necessary_N        2281 
#>  assurance          0.850
#>  assurance_attained 0.850
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; assurance_attained is the exact probability that the interval is no wider than the width

# In the units of the data (Y with standard deviation 10, the predictor
# with standard deviation 2), a width of .1 * 10 / 2 gives the same plan.
ss_aipe_rc(width = .5, true_sd_Y = 10, true_sd_X = 2, true_cor_XX = true_cor_XX,
           true_cor_YX = true_cor_YX, which_predictor = 1, assurance = .85)
#>  term               value
#>  necessary_N        2281 
#>  assurance          0.850
#>  assurance_attained 0.850
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; assurance_attained is the exact probability that the interval is no wider than the width
```
