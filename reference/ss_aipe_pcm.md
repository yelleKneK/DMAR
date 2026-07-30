# Sample Size Planning for Polynomial Change Models in Longitudinal Study

This function plans sample size with respect to the group-by-time
interaction in the context of a longitudinal design with two groups. It
plans sample size from the accuracy in parameter estimation (AIPE)
perspective, where the goal is to obtain a sufficiently narrow
confidence interval for the fixed effect polynomial change coefficient
parameter (e.g., linear, quadratic, etc.). The sample size returned can
be one such that (a) the expected confidence interval width is
sufficiently narrow, or (b) the observed confidence interval will be
sufficiently narrow with a specified high degree of assurance (e.g.,
.99, .95, .90, etc.). This function accompanies Kelley and Rausch
(2011).

## Usage

``` r
ss_aipe_pcm(
  variance_trend,
  error_variance = NULL,
  variance_true_minus_estimated_trend = NULL,
  duration,
  frequency,
  width,
  conf_level = 0.95,
  trend = "linear",
  assurance = NULL
)
```

## Arguments

- variance_trend:

  The variance of the individuals' true change coefficients (i.e.,
  \\\sigma^2\_{\upsilon_m}\\ in Kelley & Rausch, 2011) for the
  polynomial trend (e.g., linear, quadratic, etc.) of interest

- error_variance:

  The true level one error variance (i.e., \\\sigma^2\_{\epsilon}\\ in
  Kelley & Rausch, 2011). Either `error_variance` or
  `variance_true_minus_estimated_trend` must be supplied; if
  `variance_true_minus_estimated_trend` is given directly,
  `error_variance` may be omitted.

- variance_true_minus_estimated_trend:

  The variance of the difference between the \\m\\th true change
  coefficient minus the \\m\\th estimated change coefficient (i.e.,
  \\\sigma^2\_{\hat{\pi}\_{m} - \pi\_{m}}\\ from Equation 19 in Kelley &
  Rausch, 2011). When derived from `error_variance` this equals
  \\\sigma^2\_{\epsilon} f^{2p} / \sum_t c\_{mt}^2\\, where \\f\\ is the
  frequency, \\p\\ the polynomial order, and \\\sum_t c\_{mt}^2\\ the
  sum of squared orthogonal polynomial contrast weights over the
  measurement occasions. A user who already has this variance may supply
  it directly and omit `error_variance`.

- duration:

  The duration of the study

- frequency:

  The number of times measurement occurs within each unit of time

- width:

  Width of the confidence interval

- conf_level:

  The desired level of confidence for the confidence interval that will
  be computed at the completion of the study

- trend:

  The polynomial trend (1st-3rd) of interest specified as "linear",
  "quadratic", or "cubic"

- assurance:

  Value with which confidence can be placed that describes the
  likelihood of obtaining a confidence interval less than the value
  specified (e.g, .80, .90, .95)

## Value

A `data.frame` (class `dmar_tbl`) with a single row,
`necessary_n_per_group`, giving the necessary number of subjects *per
group* (the total study size is twice this value) for the combination of
the desired confidence interval width, confidence level, optional
assurance, and the population parameters at the specified design.

## Note

Like in all formal sample size planning methods that require the value
of one or more population parameter(s), if the population parameters are
incorrectly specified, there is no guarantee that the sample size this
function returns will be accurate. Of course, the further away from the
true values, the further away the true sample size will tend to be.

The number of timepoints in a study (say \\M\\) is defined by \\f \times
D + 1\\, where \\f\\ is the frequency and \\D\\ is the duration.

## References

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Kelley, K., & Rausch, J. R. (2011). Sample size planning for
longitudinal models: Accuracy in parameter estimation for polynomial
change parameters. *Psychological Methods, 16*(4), 391–405.
[doi:10.1037/a0023352](https://doi.org/10.1037/a0023352)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 11, 15.)

## See also

[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)
for the power analytic analog (planning to detect the group-by-time
change difference rather than to estimate it precisely) on the same
model, and
[`ss_aipe_pcm_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md)
for a Monte Carlo check of how parameter misspecification affects the
plan.

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The examples reproduce the tolerance-of-antisocial-thinking illustration
# of Kelley and Rausch (2011, Tables 1 and 2), which draws on the National
# Youth Survey data also used by Raudenbush and Liu (2001). The level-one
# error variance is 0.0262 and the between-subject slope variance is 0.003.
# The planner finds the sample size needed for a confidence interval on the
# group-by-time slope difference that is no wider than `width`. The returned
# necessary_n_per_group is per group, so the total study size is twice that. Unlike
# power analysis, the value of the slope is not needed here: the confidence
# interval width does not depend on it.

# (1) Expected-width planning. With five measurement occasions
#     (M = frequency * duration + 1 = 1 * 4 + 1) and a target width of
#     0.025, the expected 95% confidence interval is sufficiently narrow at
#     278 subjects per group (Kelley & Rausch, 2011, Table 1, T = 5).
ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
            duration = 4, frequency = 1, width = 0.025, conf_level = .95)
#>  term                  value
#>  necessary_n_per_group 278  
#> 
#> Confidence level: 95%

# (2) More measurement occasions sharpen the estimate. Extending the study
#     so that M = 10 (duration = 9, frequency = 1) cuts the expected-width
#     requirement from 278 to 165 per group (Kelley & Rausch, 2011, Table 1,
#     T = 10).
ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
            duration = 9, frequency = 1, width = 0.025, conf_level = .95)
#>  term                  value
#>  necessary_n_per_group 165  
#> 
#> Confidence level: 95%

# (3) A wider tolerated interval costs less. Relaxing the target width from
#     0.025 to 0.05 at M = 5 drops the requirement from 278 to 71 per group
#     (Kelley & Rausch, 2011, Table 1, T = 5).
ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
            duration = 4, frequency = 1, width = 0.05, conf_level = .95)
#>  term                  value
#>  necessary_n_per_group 71   
#> 
#> Confidence level: 95%

# (4) Adding an assurance parameter. Requiring 85% assurance that the
#     realized confidence interval will be no wider than 0.025 raises the
#     M = 5 requirement from 278 to 295 per group (Kelley & Rausch, 2011,
#     Table 2, T = 5). Assurance guards against the expected-width plan being
#     too small for the particular sample obtained.
ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
            duration = 4, frequency = 1, width = 0.025, conf_level = .95,
            assurance = .85)
#>  term                  value
#>  necessary_n_per_group 295  
#> 
#> Confidence level: 95%

# (5) A higher assurance costs more. Demanding 99% assurance rather than 85%
#     raises the per-group requirement further, from 295 to 316.
ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
            duration = 4, frequency = 1, width = 0.025, conf_level = .95,
            assurance = .99)
#>  term                  value
#>  necessary_n_per_group 316  
#> 
#> Confidence level: 95%
```
