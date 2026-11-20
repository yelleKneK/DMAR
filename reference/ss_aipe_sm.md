# Sample Size Planning for Accuracy in Parameter Estimation (AIPE) of the Standardized Mean

Plans the sample size needed for a sufficiently narrow confidence
interval for the population standardized mean, the mean divided by the
standard deviation, from the accuracy in parameter estimation (AIPE)
perspective.

## Usage

``` r
ss_aipe_sm(true_sm, width, conf_level = 0.95, assurance = NULL)
```

## Arguments

- true_sm:

  The planning value of the population standardized mean (the population
  mean divided by the population standard deviation)

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage, (i.e., 1 - Type I error
  rate). A level above .999999998, whose tail areas are below 1e-9, is
  refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md))

- assurance:

  Optional. The probability that the confidence interval of the planned
  study is no wider than `width`: `NULL` (the default), which plans for
  the expected width, or a single number strictly between 0 and 1. Any
  such value may be used. Values such as `.80`, `.90`, and `.99` are the
  usual ones, and a value at or below `.50` is accepted as well (see
  *Details*). Anything else is refused with a condition of class
  `dmar_assurance_invalid`. An assurance above .999999998 is refused
  with a condition of class `dmar_argument_invalid`: the assurance step
  finds noncentral *t* limits with tail areas of `(1 - assurance) / 2`,
  which must be at least 1e-9.

## Value

A 1-row `data.frame` with columns `term` and `value`:

- necessary_N:

  The necessary total sample size in order to achieve the desired degree
  of accuracy (i.e., the sufficiently narrow confidence interval)

## Details

**The plan for the expected width.** Without `assurance` the function
finds the smallest sample size at which the noncentral *t* confidence
interval for the standardized mean
([`ci_sm`](https://yelleknek.github.io/DMAR/reference/ci_sm.md)),
computed at `true_sm`, is no wider than `width`. The interval of the
study itself is computed at the sample standardized mean, which varies
from sample to sample, and the interval is wider the larger that value
is in magnitude. A plan for the expected width therefore leaves a
probability near one half that the interval obtained is wider than
`width`.

**The assurance step.** With `assurance` the function follows the
one-sample analog of the procedure Kelley and Rausch (2006, Equations 18
to 20 and 25 to 31) give for the standardized mean difference. At the
sample size of the plan for the expected width it finds the value that
the sample standardized mean exceeds in magnitude with probability
`1 - assurance`, by a search between the one-sided and the two-sided
confidence limits their Equations 25 and 28 name as its bounds, and then
plans for the expected width with that value in place of `true_sm`. That
value is defined for every probability strictly between 0 and 1, so
every such `assurance` is accepted. A value above about one half gives a
larger sample size than the plan for the expected width, and a value
below it gives a smaller one, which is the sample size at which the
interval is no wider than `width` with only that probability: for
`true_sm = 1`, `width = 0.5`, and 95 percent confidence the function
returns 93 without `assurance`, 93 at `.50`, 87 at `.20`, 100 at `.80`,
and 105 at `.90`. Because the value planned for is found at the first
sample size and the returned sample size differs from it, the
probability at the returned sample size is near `assurance` and not
exactly it; in that example it is .519, .228, .812, and .925 at the four
values of `assurance`, computed from the noncentral *t* distribution.

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
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
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

ss_aipe_sm(true_sm = 10, width = 2.5, conf_level = .95, assurance = .90)
#>  term        value
#>  necessary_N 150  
#> 
#> Confidence level: 95%

# The plan for the expected width, and the sample size at which the
# interval is no wider than 0.5 with probability .50, for a standardized
# mean of 1.
ss_aipe_sm(true_sm = 1, width = 0.5)
#>  term        value
#>  necessary_N 93   
#> 
#> Confidence level: 95%
ss_aipe_sm(true_sm = 1, width = 0.5, assurance = .50)
#>  term        value
#>  necessary_N 93   
#> 
#> Confidence level: 95%
```
