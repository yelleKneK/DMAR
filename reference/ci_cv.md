# Confidence Interval for the Coefficient of Variation

Computes the confidence interval for the population coefficient of
variation, the standard deviation relative to the mean, from the
noncentral *t* distribution, so that variability can be reported on a
scale that is free of the measurement units. The interval is exact under
normality and is not robust to departures from it (see Details).

## Usage

``` r
ci_cv(
  cv = NULL,
  mean = NULL,
  sd = NULL,
  N = NULL,
  data = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL
)
```

## Arguments

- cv:

  Sample coefficient of variation, the sample standard deviation divided
  by the sample mean, given with `N`. Must be positive.

- mean:

  Sample mean; must be positive. It is given with `sd` and `N`.

- sd:

  Sample standard deviation (square root of the unbiased estimate of the
  variance). Must be positive.

- N:

  Total sample size, given with `cv` or with `mean` and `sd`. Must
  exceed 1; it need not be a whole number.

- data:

  Numeric vector of at least two finite scores, given in place of the
  summary statistics: the coefficient of variation and the sample size
  are taken from it.

- conf_level:

  Desired confidence level (1-Type I error rate). A level above
  .999999998, whose tail areas are below 1e-9, is refused

- alpha_lower:

  The lower tail area. An omitted confidence level yields to supplied
  tail areas; an explicit level must agree with them.

- alpha_upper:

  The upper tail area, supplied together with `alpha_lower`. The table
  records the coverage, `1 - alpha_lower - alpha_upper`, as its level
  and, when the two tail areas differ, the areas themselves, which the
  line beneath the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound. A
  nonzero tail area below 1e-9 is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

## Value

A 4-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The rows are ordered so the two point estimates sit
between the confidence limits: `"ci_lower"` (lower confidence limit on
the coefficient of variation), `"cv"` (the sample coefficient of
variation), `"cv_unbiased"` (the nearly unbiased estimate, \\(1 + 1/(4
N))\\ times the sample value; see
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)), and
`"ci_upper"` (upper confidence limit). The `prob_less` and
`prob_greater` columns hold the requested tail areas: `alpha_lower` and
`1 - alpha_lower` at the lower limit, and `1 - alpha_upper` and
`alpha_upper` at the upper limit, also when the upper limit is infinite.
They are `NA` for the point-estimate rows.

## Details

**The interval.** For a random sample of \\N\\ scores from a normal
population with a positive mean, \\\sqrt{N} / k\\, where \\k = s / M\\
is the sample coefficient of variation, has a noncentral *t*
distribution on \\N - 1\\ degrees of freedom whose noncentrality
parameter is \\\sqrt{N} / \kappa\\, with \\\kappa = \sigma / \mu\\ the
population coefficient of variation (Johnson & Welch, 1940; Kelley,
2007, Equations 5 to 8; McKay, 1932). The function finds the confidence
limits for the noncentrality parameter
([`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)) and
divides \\\sqrt{N}\\ by each: the lower limit for \\\kappa\\ comes from
the upper limit for the noncentrality parameter, and the upper limit
from the lower.

**Normality.** The interval is exact when the scores are a random sample
from a normal distribution with a positive mean: its coverage is then
the stated level at every sample size. It is not robust to departures
from normality, and a larger sample does not repair it, because the
distribution that is inverted is that of a normal sample whatever its
size. In simulations of a nominal 95 percent interval, with 10,000
samples in each condition (population coefficients of variation of .3
and .5, sample sizes of 30 and 200), coverage was .81 to .85 for Laplace
data, .90 to .94 for lognormal data, and .96 to .97 for gamma data,
where it was .945 to .951 for normal data, the stated level within the
error of the simulation. For the Laplace and the lognormal data it was
lower at 200 than at 30. A large-sample interval that does not assume
normality follows from the asymptotic normal distribution of \\k\\
(Kelley, Darku, & Chattopadhyay, 2018, Equations 29 to 33). It needs a
large sample: with 400,000 samples in each of the same conditions it
covered .806 to .905 at \\N = 30\\ and .894 to .942 at \\N = 200\\,
below the stated level in every condition, normal data included. The
function does not compute it.

**An infinite upper limit.** The upper limit is infinite when \\\sqrt{N}
/ k\\ is at most the `1 - alpha_upper` quantile of the central *t*
distribution on \\N - 1\\ degrees of freedom, which happens when the
sample coefficient is large for the sample size. The lower confidence
limit of the noncentrality parameter is then not positive: at that tail
area the data do not rule out a population mean of zero, and no finite
coefficient of variation bounds the interval above. The infinite limit
is part of the exact interval. Over repeated samples from a normal
population with a positive mean, the interval covers the population
coefficient with the stated probability, counting the intervals that are
bounded above and those that are not. In 20,000 simulated samples with a
population coefficient of .8 and \\N = 6\\, the upper limit of the 95
percent interval was infinite in 31 percent and coverage was .952; with
a coefficient of .5 and \\N = 5\\ it was infinite in 9 percent and
coverage was .948. The stated level describes all of the intervals
together and neither kind alone (computed from the noncentral *t*
distribution, coverage is .92 and .73 among the intervals with an
infinite limit and .96 and .97 among the others). Nothing is signaled
for an infinite limit, and the `prob_less` and `prob_greater` columns
hold the requested tail areas there as at a finite limit. A sample whose
mean is not positive has no interval and is refused; with a population
coefficient of .5 or less and five or more scores the probability of
such a sample is below .00001. Where such samples are not rare, coverage
among the samples that have an interval is somewhat above the stated
level: it is .962 for a 95 percent interval with a coefficient of 1 and
\\N = 5\\, where 1.3 percent of samples have a mean that is not
positive. Of the 20,000 samples with a coefficient of .8 above, 25 had
such a mean, and the figures given describe the other 19,975.

## References

Johnson, N. L., & Welch, B. L. (1940). Applications of the non-central
*t*-distribution. *Biometrika, 31*(3–4), 362–389.
[doi:10.1093/biomet/31.3-4.362](https://doi.org/10.1093/biomet/31.3-4.362)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods, 23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on confidence intervals for effect
sizes.)

McKay, A. T. (1932). Distribution of the coefficient of variation and
the extended *t* distribution. *Journal of the Royal Statistical
Society, 95*(4), 695–698.

## See also

[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_ratio_of_means()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md),
[`ci_ratio_of_means_paired()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Fifteen scores from a normal population with a mean of 5 and a
# standard deviation of 1, whose coefficient of variation is .20.
set.seed(113)
x <- rnorm(15, mean = 5, sd = 1)

# The interval from the scores, from their mean and standard deviation,
# and from the sample coefficient itself: three routes to one interval.
ci_cv(data = x)
#>  term        value prob_less prob_greater
#>  ci_lower    0.150 0.025     0.975       
#>  cv          0.207 <NA>      <NA>        
#>  cv_unbiased 0.210 <NA>      <NA>        
#>  ci_upper    0.335 0.975     0.025       
#> 
#> Confidence level: 95%
ci_cv(mean = mean(x), sd = sd(x), N = 15)
#>  term        value prob_less prob_greater
#>  ci_lower    0.150 0.025     0.975       
#>  cv          0.207 <NA>      <NA>        
#>  cv_unbiased 0.210 <NA>      <NA>        
#>  ci_upper    0.335 0.975     0.025       
#> 
#> Confidence level: 95%
ci_cv(cv = sd(x) / mean(x), N = 15)
#>  term        value prob_less prob_greater
#>  ci_lower    0.150 0.025     0.975       
#>  cv          0.207 <NA>      <NA>        
#>  cv_unbiased 0.210 <NA>      <NA>        
#>  ci_upper    0.335 0.975     0.025       
#> 
#> Confidence level: 95%

# Unequal tail areas, with .01 below the interval and .04 above it.
ci_cv(data = x, alpha_lower = .01, alpha_upper = .04)
#>  term        value prob_less prob_greater
#>  ci_lower    0.142 0.010     0.990       
#>  cv          0.207 <NA>      <NA>        
#>  cv_unbiased 0.210 <NA>      <NA>        
#>  ci_upper    0.316 0.960     0.040       
#> 
#> Confidence level: 95% (unequal tail areas: 0.01 below the lower limit, 0.04 above the upper limit)
```
