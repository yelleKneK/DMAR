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
  true_std_contrast,
  c_weights,
  width,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  assurance = NULL,
  assurance_method = "analytic"
)
```

## Arguments

- true_std_contrast:

  Population standardized contrast

- c_weights:

  The contrast weights, as
  [`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md) takes
  them: fractions whose positive values sum to 1 and whose sum is zero
  (for example `c(1/2, 1/2, -1)`, not `c(1, 1, -2)`), so that the
  standardized contrast is a difference of averages of means in units of
  the error standard deviation. Weights that break the rule are refused
  with the message of
  [`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md) and a
  condition of class `dmar_argument_invalid`. The rule is applied to the
  values, so integer weights such as `-1:1` are accepted, as they are by
  [`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md).

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage (i.e., 1 - Type I error
  rate). Default is `.95`, which gives a symmetric two-sided interval.
  Specify either `conf_level` or both of `alpha_lower` and
  `alpha_upper`, not both. A level above .999999998, whose tail areas
  are below 1e-9, is refused.

- alpha_lower:

  Lower-tail Type I error rate, used to plan an asymmetric confidence
  interval. When supplied together with `alpha_upper`, the planned
  interval has lower-tail probability `alpha_lower` and upper-tail
  probability `alpha_upper`. Set `conf_level = NULL` when supplying
  these.

- alpha_upper:

  Upper-tail Type I error rate, used together with `alpha_lower` to plan
  an asymmetric confidence interval. Each tail area must be at least
  1e-9: below it the noncentral *t* limits lose accuracy (see "Small
  tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- assurance:

  Parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be NULL or greater than .5 and less than 1). An assurance above
  .999999998 is refused: the floor of 1e-9 on a tail area (see "Small
  tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md))
  applies to `(1 - assurance) / 2`.

- assurance_method:

  How the probability that the
  [`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md)
  interval is no wider than `width` is found, used only with
  `assurance`. Only `"analytic"` (the default) is offered here: the
  probability is computed exactly, and the plan is the smallest sample
  size per group, not below the plan without `assurance`, at which it
  reaches `assurance`. `"simulation"`, which
  [`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md),
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
  with `noncentral = TRUE`, and
  [`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)
  offer, where the probability has no exact form, is refused with a
  condition of class `dmar_argument_invalid`, whether or not `assurance`
  is given, since a simulated estimate of an exact probability would
  only add Monte Carlo error.

## Value

- necessary_n_per_group:

  Necessary sample size *per group*

- assurance, assurance_attained:

  With `assurance` only: the request, and the exact probability at
  `necessary_n_per_group` that the
  [`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md)
  interval is no wider than `width` (at least `assurance`). The printed
  footer names the assurance method, which the attribute
  `assurance_method` records, with `assurance_probability` `"exact"`.

## Details

The plan is for the interval of
[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md), the
limits for the noncentrality parameter of a noncentral *t* on \\J(n -
1)\\ degrees of freedom at the observed *t* times \\\sqrt{\sum_j c_j^2 /
n}\\. Without `assurance` it is the smallest *n* per group, and at least
4, at which that interval, evaluated at the planning value, is no wider
than `width`; a `width` that would need more than \\2^{53}\\ (about
\\9.0 \times 10^{15}\\, the largest number up to which R holds every
whole number exactly) per group is refused with a condition of class
`dmar_argument_invalid`.

With `assurance` the plan is the smallest sample size per group, not
below the plan without assurance, at which the probability that the
interval is no wider than `width` is at least `assurance`
(`assurance_method = "analytic"`), and that probability, computed
exactly, is `assurance_attained`. The observed *t* is noncentral *t* on
\\J(n - 1)\\ degrees of freedom with noncentrality \\\delta /
\sqrt{\sum_j c_j^2 / n}\\, the width is smallest at \\t = 0\\ and rises
on either side, so the event is \\t_L \le T \le t_U\\ for the two values
at which the width equals the target, and its probability is that of the
noncentral *t*. The plan is found by doubling a step upward from the
plan without assurance and then bisecting, which assumes that the
probability does not decrease as *n* grows from there. In 96 planning
sets (contrasts of 0, .5, -1, and 2; weights for two, three, and five
groups; widths of .5 and 1; assurances of .80 and .99; and tail areas of
.025 and .025 or .01 and .04) it did not decrease at any *n* from the
plan without assurance to ten past the plan, and the plan was the
smallest *n* found by raising *n* one at a time.

The step of Kelley and Rausch (2006), as Lai and Kelley (2012, pp.
358–359 and 361) apply it to a standardized contrast, plans without
assurance for the value \\t^\*\\ that \\\|T\|\\ does not exceed with
probability `assurance` at the plan without assurance, carried back to
the scale of the standardized contrast; it does not compute the
probability. That step can fall short of the request: at
`true_std_contrast = 0`, `c_weights = c(1, -1)`, `width = .8`, and
`assurance = .95` it gives 49 per group, which attains .9499, where the
plan here is 50, which attains .9946. It can also be larger than the
smallest plan that attains the request: in the 96 sets above it attained
the request in every set and was one case or more larger in 17.

Past about \\10^{11}\\ per group the plan and `assurance_attained` are
limited by the precision of doubles. The width changes by about
\\1/(2n)\\ of itself from one *n* to the next, and for a nonzero
contrast it is the difference of two limits near the noncentrality
\\\delta / \sqrt{\sum_j c_j^2 / n}\\, which a double holds to about
\\1.1 \times 10^{-16}\\ of its size. Against a 512-bit computation of
the interval, with `c_weights = c(1, -1)` and `conf_level = .95`, and
with `width` and the tail areas as R holds them (each tail area, (1 -
.95) / 2, is 0.025000000000000022 in doubles): the plans checked up to
4.6e9 per group were the smallest *n*; at `true_std_contrast = 2` and
`width` of 1e-5, 1e-6, 3e-7, and 1e-7 the plan was 2, 1,850, 487,182,
and 6,485,787 below the smallest *n* (about 4.6e11, 4.6e13, 5.1e14, and
4.6e15), and at .5 it was 4 above, 1,956 below, 129 above, and 4,854,160
above it (about 3.2e11, 3.2e13, 3.5e14, and 3.2e15). At a contrast of 0
the width is computed from the normal quantiles, and the plan is below
the smallest *n* only when the *n* at which the width equals `width` is
within about \\10^{-15} n\\ above a whole number, and then by at most
the larger of 1 and about \\10^{-15} n\\. At `width = 1e-6` that *n* is
30,731,670,565,552.9985, and the plan, 30,731,670,565,553 per group, is
the smallest *n*; at 3e-7 that *n* is 341,463,006,283,922.21, and the
plan, 341,463,006,283,922, is 1 below the smallest *n*; and the plan was
2 below at 1e-7 (about 3.1e15 per group), 4 below at 7e-8 (6.3e15), and
8 below at 6e-8 (8.5e15). With `assurance = .9`, `assurance_attained`
was within 1e-6 of the exact probability up to 4.6e9 per group.

With unequal tails the width at an observed *t* is not the width at
\\-t\\: the interval at \\-t\\ with tail areas (`alpha_lower`,
`alpha_upper`) is the reflection of the interval at *t* with the two
exchanged. A negative `true_std_contrast` therefore gets the plan of its
magnitude with `alpha_lower` and `alpha_upper` exchanged, so that
reversing the contrast and the tails together leaves the plan unchanged.
Without `assurance` the signed planning value already gives that plan.
With equal tails the plan does not depend on the sign.

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
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
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
# contrast is, with 90 percent assurance, no wider than .4 (the last row
# is the exact probability at that sample size):

ss_aipe_sc(true_std_contrast=.6, c_weights=c(.5, .5, -.5, -.5, 0), width=.4, assurance=.90)
#>  term                  value
#>  necessary_n_per_group 102  
#>  assurance             0.900
#>  assurance_attained    0.964
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; assurance_attained is the exact probability that the interval is no wider than the width

# Asymmetric confidence interval: most of the alpha goes in the upper tail
# (e.g., when a one-sided concern dominates). Pass alpha_lower and
# alpha_upper instead of conf_level.
ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
           conf_level = NULL, alpha_lower = .01, alpha_upper = .04)
#>  term                  value
#>  necessary_n_per_group 108  
```
