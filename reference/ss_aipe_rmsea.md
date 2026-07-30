# Sample Size Planning for RMSEA in SEM

Sample size planning for the population root mean square error of
approximation (RMSEA) from the accuracy in parameter estimation (AIPE)
perspective. The sample size is planned so that the expected width of a
confidence interval for the population RMSEA is no larger than desired.

## Usage

``` r
ss_aipe_rmsea(RMSEA, df, width, conf_level = 0.95)
```

## Arguments

- RMSEA:

  The input RMSEA value

- df:

  Degrees of freedom of the model

- width:

  Desired confidence interval width

- conf_level:

  Desired confidence level (e.g., .90, .95, .99, etc.)

## Value

Returns the necessary total sample size in order to achieve the desired
degree of accuracy (i.e., the sufficiently narrow confidence interval).

## References

Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
root mean square error of approximation: Sample size planning for narrow
confidence intervals. *Multivariate Behavioral Research, 46*, 1–32.
[doi:10.1080/00273171.2011.543027](https://doi.org/10.1080/00273171.2011.543027)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ss_aipe_rmsea(RMSEA = .035, df = 50, width = .05, conf_level = .95)
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#>  term        value
#>  necessary_N 361  
#> 
#> Confidence level: 95%
```
