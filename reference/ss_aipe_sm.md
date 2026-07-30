# Sample Size Planning for Accuracy in Parameter Estimation (AIPE) of the Standardized Mean

Plans the sample size needed for a sufficiently narrow confidence
interval for the population standardized mean, the mean divided by the
standard deviation, from the accuracy in parameter estimation (AIPE)
perspective.

## Usage

``` r
ss_aipe_sm(sm, width, conf_level = 0.95, assurance = NULL, ...)
```

## Arguments

- sm:

  The population standardized mean

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- assurance:

  Parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be `NULL` or between zero and unity)

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

- sample_size:

  The necessary sample size in order to achieve the desired degree of
  accuracy (i.e., the sufficiently narrow confidence interval)

## References

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574.
[doi:10.1177/0013164401614002](https://doi.org/10.1177/0013164401614002)

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Kelley, K. (2005). The effects of nonnormal distributions on confidence
intervals around the standardized mean difference: Bootstrap and
parametric confidence intervals, *Educational and Psychological
Measurement, 65*, 51–69.
[doi:10.1177/0013164404264850](https://doi.org/10.1177/0013164404264850)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`ci_sm`](https://yelleknek.github.io/DMAR/reference/ci_sm.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose the population mean is believed to be 20, and the population
# standard deviation is believed to be 2; thus the population standardized
# mean is believed to be 10. To determine the necessary sample size for a
# study so that the full width of the 95 percent confidence interval
# obtained in the study will be, with 90% assurance, no wider than 2.5,
# the function should be specified as follows.

ss_aipe_sm(sm = 10, width = 2.5, conf_level = .95, assurance = .90)
#> Warning: During the iterative sample size search, conf_limits_nct() reported a noncentrality parameter exceeding 37.62 in magnitude in 254 intermediate evaluations, the limit at which R's pt()/qt() can return accurate noncentral t probabilities. The returned sample size may be affected; see ?conf_limits_nct.
#>  term        value
#>  necessary_N 150  
#> 
#> Confidence level: 95%
```
