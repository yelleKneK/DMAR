# Sample Size Planning for Accuracy in Parameter Estimation (AIPE) of the Standardized Contrast in ANOVA

Plans the sample size per group so that the confidence interval for a
standardized contrast of means in a fixed effects analysis of variance,
the interval computed by
[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md), is
sufficiently narrow, an application of the accuracy in parameter
estimation (AIPE) approach to the standardized contrast.

## Usage

``` r
ss_aipe_sc(
  psi_standardized,
  c_weights,
  width,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  assurance = NULL,
  ...
)
```

## Arguments

- psi_standardized:

  Population standardized contrast

- c_weights:

  The contrast weights

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage (i.e., 1 - Type I error
  rate). Default is `.95`, which gives a symmetric two-sided interval.
  Specify either `conf_level` or both of `alpha_lower` and
  `alpha_upper`, not both.

- alpha_lower:

  Lower-tail Type I error rate, used to plan an asymmetric confidence
  interval. When supplied together with `alpha_upper`, the planned
  interval has lower-tail probability `alpha_lower` and upper-tail
  probability `alpha_upper`. Set `conf_level = NULL` when supplying
  these.

- alpha_upper:

  Upper-tail Type I error rate, used together with `alpha_lower` to plan
  an asymmetric confidence interval.

- assurance:

  Parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be NULL or between zero and unity)

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

- necessary_n_per_group:

  Necessary sample size *per group*

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

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose the population standardized contrast is believed to be .6
# in some 5-group ANOVA model. The researcher is interested in comparing
# the average of means of group 1 and 2 with the average of group 3 and 4.

# To calculate the necessary sample size per group such that the width
# of 95 percent confidence interval of the standardized
# contrast is, with 90 percent assurance, no wider than .4:

ss_aipe_sc(psi_standardized=.6, c_weights=c(.5, .5, -.5, -.5, 0), width=.4, assurance=.90)
#>  term                  value
#>  necessary_n_per_group 102  
#> 
#> Confidence level: 95%

# Asymmetric confidence interval: most of the alpha goes in the upper tail
# (e.g., when a one-sided concern dominates). Pass alpha_lower and
# alpha_upper instead of conf_level.
ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
           conf_level = NULL, alpha_lower = .01, alpha_upper = .04)
#>  term                  value
#>  necessary_n_per_group 108  
```
