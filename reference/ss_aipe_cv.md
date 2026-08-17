# Sample Size Planning for the Coefficient of Variation Given the Goal of Accuracy in Parameter Estimation Approach to Sample Size Planning

Determines the necessary sample size so that the expected confidence
interval width for the coefficient of variation will be sufficiently
narrow, optionally with a desired degree of certainty that the interval
will not be wider than desired. The population coefficient of variation
may be given directly as `C_of_V` or through `mu` and `sigma`, in which
case `C_of_V` is taken as `sigma / mu`. The value of `C_of_V` should be
positive.

## Usage

``` r
ss_aipe_cv(
  C_of_V = NULL,
  width = NULL,
  conf_level = 0.95,
  assurance = NULL,
  mu = NULL,
  sigma = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- C_of_V:

  Population coefficient of variation on which the sample size procedure
  is based

- width:

  Desired (full) width of the confidence interval

- conf_level:

  Confidence interval coverage; 1-Type I error rate

- assurance:

  Value with which confidence can be placed that describes the
  likelihood of obtaining a confidence interval less than the value
  specified (e.g., .80, .90, .95)

- mu:

  Population mean (specified with `sigma` when `C_of_V` is not
  specified)

- sigma:

  Population standard deviation (specified with `mu` when `C_of_V` is
  not specified)

- alpha_lower:

  Type I error for the lower confidence limit

- alpha_upper:

  Type I error for the upper confidence limit

- ...:

  For modifying parameters of functions this function calls

## Value

Returns the necessary sample size given the input specifications.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3.)

## See also

[`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md),
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose one wishes to have a confidence interval with an expected width of .10
# for a 99% confidence interval when the population coefficient of variation is .10.
ss_aipe_cv(C_of_V = .1, width = .1, conf_level = .99)
#> Warning: During the iterative sample size search, the noncentrality parameter exceeded 37.62 in magnitude (the limit of R's noncentral t accuracy) in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_nct.
#>  term        value
#>  necessary_N 20   
#> 
#> Confidence level: 99%

# The same planning problem parameterized by the population mean and standard
# deviation: mu = 10 and sigma = 1 imply the same coefficient of variation, .10.
ss_aipe_cv(mu = 10, sigma = 1, width = .1, conf_level = .99)
#> Warning: During the iterative sample size search, the noncentrality parameter exceeded 37.62 in magnitude (the limit of R's noncentral t accuracy) in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_nct.
#>  term        value
#>  necessary_N 20   
#> 
#> Confidence level: 99%

# Ensuring that the confidence interval will be sufficiently narrow with a 99\%
# certainty for the situation above.
ss_aipe_cv(C_of_V = .1, width = .1, conf_level = .99, assurance = .99)
#> Warning: During the iterative sample size search, the noncentrality parameter exceeded 37.62 in magnitude (the limit of R's noncentral t accuracy) in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_nct.
#> Warning: During the iterative sample size search, the noncentrality parameter exceeded 37.62 in magnitude (the limit of R's noncentral t accuracy) in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_nct.
#>  term        value
#>  necessary_N 33   
#> 
#> Confidence level: 99%
```
