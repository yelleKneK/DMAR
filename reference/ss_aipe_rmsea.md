# Sample Size Planning for RMSEA in SEM

Sample size planning for the population root mean square error of
approximation (RMSEA) from the accuracy in parameter estimation (AIPE)
perspective. The sample size is planned so that the expected width of a
confidence interval for the population RMSEA is no larger than desired.

## Usage

``` r
ss_aipe_rmsea(true_rmsea, df, width, conf_level = 0.9)
```

## Arguments

- true_rmsea:

  The population RMSEA planned for

- df:

  Degrees of freedom of the model

- width:

  Desired confidence interval width

- conf_level:

  Desired confidence level of the planned interval. The default, 0.90,
  plans the conventional interval for the RMSEA, the one
  [`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)
  forms by default (see Details).

## Value

Returns the necessary total sample size in order to achieve the desired
degree of accuracy (i.e., the sufficiently narrow confidence interval).

## Details

The width planned is that of the interval
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)
forms when the sample RMSEA equals `true_rmsea`, with the sample size
multiplier \\N - 1\\ of Kelley and Lai (2011). The RMSEA is
conventionally reported with a 90 percent interval, because a \\100(1 -
2\alpha)\\ percent interval decides the one-sided tests of close and not
close fit at level \\\alpha\\ (MacCallum, Browne, & Sugawara, 1996, p.
136), and every RMSEA interval in the package defaults to it, the
interval planned here included. Kelley and Lai (2011) planned 95 percent
intervals in their examples; a higher level needs a larger sample for
the same width. For an RMSEA of .045 on 5 degrees of freedom and a width
of .060, the plan is 802 at the default and 1068 with
`conf_level = 0.95`. The printed table states the level beneath it.

## References

Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
root mean square error of approximation: Sample size planning for narrow
confidence intervals. *Multivariate Behavioral Research, 46*, 1–32.
[doi:10.1080/00273171.2011.543027](https://doi.org/10.1080/00273171.2011.543027)

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

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
# The expected width of the 90 percent interval, the default, is .05.
ss_aipe_rmsea(true_rmsea = .035, df = 50, width = .05)
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#>  term        value
#>  necessary_N 294  
#> 
#> Confidence level: 90%

# The same width for a 95 percent interval needs a larger sample.
ss_aipe_rmsea(true_rmsea = .035, df = 50, width = .05, conf_level = .95)
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#> Note: The lower confidence limit of the noncentrality parameter is at its lower bound, so the lower RMSEA limit is set to 0 based on RMSEA's definition.
#>  term        value
#>  necessary_N 361  
#> 
#> Confidence level: 95%
```
