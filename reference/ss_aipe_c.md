# Sample Size Planning for an ANOVA Contrast From the Accuracy in Parameter Estimation (AIPE) Perspective

Plans the sample size *per group* so that the confidence interval for an
unstandardized contrast of means in a fixed effects analysis of variance
is sufficiently narrow, following the accuracy in parameter estimation
(AIPE) approach: the design goal is a contrast estimated with the
precision the research question requires, not merely one detected as
nonzero. AIPE sample size planning for ANOVA and ANCOVA contrasts is
developed in Lai and Kelley (2012).

## Usage

``` r
ss_aipe_c(
  true_var_error,
  c_weights,
  width,
  conf_level = 0.95,
  assurance = NULL,
  true_kurtosis = 0
)
```

## Arguments

- true_var_error:

  The population error variance: the variance of the errors common to
  the groups, which the mean square error (the mean square within
  groups) of an analysis of variance estimates. It is a variance, not a
  standard deviation: a planning value stated as an error standard
  deviation is squared before it is given. Must be a single finite
  positive number. This is the one way to state the error variability,
  as in
  [`ss_aipe_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md)
  and
  [`ss_power_contrast`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md).

- c_weights:

  The contrast weights

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage, (i.e., 1 - Type I error
  rate). A level above .999999998, whose tail areas are below 1e-9, is
  refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- assurance:

  Parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be NULL or between zero and unity)

- true_kurtosis:

  The population excess kurtosis of the errors (their fourth
  standardized moment less 3), on the scale
  [`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)
  reports, so the value it gives for the residuals of a pilot study can
  be passed directly. It is used only with `assurance`, through the
  effective degrees of freedom described in Details: a negative value
  (lighter tails than the normal, as bounded outcomes such as
  proportions and rating scale scores often have) makes the sample
  variance less variable than normal theory implies, and a positive
  value (heavier tails) makes it more variable. The default `0` is the
  normal theory value. Must be a single finite number greater than -2;
  no distribution has a smaller excess kurtosis, and only a symmetric
  two-point distribution attains -2.

## Value

A 1-row `data.frame` with columns `term` and `value`:

- necessary_n_per_group:

  the necessary sample size *per group*

A plan made with `assurance` and a nonzero `true_kurtosis` keeps that
one row, carries the value as the attribute `true_kurtosis`, prints a
line under the table saying that the assurance allows for that excess
kurtosis of the errors, and reports it as a `true_kurtosis` column of
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html).
With the default `true_kurtosis = 0`, or without `assurance`, no
attribute is added.

## Details

Without `assurance`, the planned *n* per group is the smallest integer
(and at least 2) satisfying \$\$n \ge \frac{4\\ t^2\_{1-\alpha/2,\\
J(n-1)}\\ \sigma^2 \sum_j c_j^2}{\omega^2},\$\$ where \\\sigma^2\\ is
the error variance (`true_var_error`), \\c_j\\ are the contrast weights
of the \\J\\ groups, and \\\omega\\ is the desired full width. That is
the width the interval would have if the pooled sample variance \\S^2\\
equaled \\\sigma^2\\; because \\S^2\\ varies from sample to sample, the
realized width exceeds \\\omega\\ in a large proportion of studies, a
proportion that approaches one half as *n* grows. With `assurance` equal
to \\\gamma\\, the right side is multiplied by \\\chi^2\_{\gamma,\\ \nu}
/ \nu\\, the \\\gamma\\ quantile of a chi square on \\\nu\\ degrees of
freedom divided by \\\nu\\, which is the \\\gamma\\ quantile of
\\S^2/\sigma^2\\ when \\\nu S^2/\sigma^2\\ is chi square on \\\nu\\
degrees of freedom; the realized width is then at most \\\omega\\ with
probability \\\gamma\\. For normal errors \\\nu = J(n - 1)\\ and the chi
square law is exact.

For errors that are not normal, the variability of \\S^2\\ depends on
their excess kurtosis \\\kappa\\. With \\J\\ groups of \\n\\, the pooled
variance has \\\mathrm{var}(S^2) = \sigma^4 \[2/\\J(n - 1)\\ +
\kappa/(nJ)\]\\, so, for an assurance above one half, a plan made under
normality delivers more assurance than requested when \\\kappa \< 0\\
and less when \\\kappa \> 0\\. The `true_kurtosis` argument replaces
\\J(n - 1)\\ in the chi square quantile with the effective degrees of
freedom \$\$\nu = \frac{2}{2/\\J(n - 1)\\ + \kappa/(nJ)},\$\$ the
degrees of freedom for which \\\chi^2\_\nu/\nu\\ has the variance of
\\S^2/\sigma^2\\; \\\kappa = 0\\ gives \\\nu = J(n - 1)\\ and the normal
theory plan, and the degrees of freedom of the *t* quantile are not
changed. This is the approximation of Box (1953, p. 319, Eq. 3): a
variance estimate on \\\phi\\ degrees of freedom from a population with
excess kurtosis \\\gamma_2\\ has variance \\2\sigma^4/(d\phi)\\, with
\\d = \[1 + \phi\gamma_2/\\2(\phi + 1)\\\]^{-1}\\, and is distributed
approximately as a normal theory estimate on \\d\phi\\ degrees of
freedom; with \\\phi = n - 1\\ in each group, \\\nu\\ is \\d\\ times
\\J(n - 1)\\. The approximation matches two moments and leaves the shape
to the chi square. It requires errors with a finite fourth moment, and
in Box's sampling experiments (p. 325) the large sample result was
approached more slowly when the kurtosis was large. The kurtosis of a
pilot study's residuals is itself an imprecise estimate in samples of
the size pilots have, so the value supplied is a planning value, as the
error variance is. The sensitivity companion,
[`ss_aipe_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
simulates normal errors and has no `true_kurtosis` argument, so it
evaluates a plan under the normal theory law, not under the errors a
nonzero `true_kurtosis` allows for.

A Monte Carlo check of the approximation used two groups,
`c_weights = c(1, -1)`, an error variance of 1, a width of .8, and an
assurance of .80, with 10,000 replications per condition. With uniform
errors (excess kurtosis -1.2), the plan made with the kurtosis supplied
(53 per group) met the width target in a proportion .811 of the
replications, where the plan made under normality (55 per group) met it
in .915; with *t* errors on 10 degrees of freedom (excess kurtosis 1),
the plan with the kurtosis supplied (56 per group) met it in .811, where
the normal theory plan (55 per group) met it in .776, short of the .80
requested. The study's script and results are kept with the package's
development records and are not distributed with the package.

## Note

Be sure to use the error variance and not its square root (i.e., the
standard deviation of the errors).

## References

Box, G. E. P. (1953). Non-normality and tests on variances. *Biometrika,
40*(3–4), 318–335.
[doi:10.1093/biomet/40.3-4.318](https://doi.org/10.1093/biomet/40.3-4.318)

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample-size planning.
*Evaluation & the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
[`ss_aipe_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova.md),
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose the population error variance of some three-group ANOVA model
# is believed to be 40. The researcher is interested in the difference
# between the mean of group 1 and the average of means of group 2 and 3.
# To plan the sample size so that, with 90 percent certainty, the
# obtained 95 percent full confidence interval width is no wider than 3:

ss_aipe_c(true_var_error = 40, c_weights = c(1, -0.5, -0.5),
          width = 3, assurance = .90)
#>  term                  value
#>  necessary_n_per_group 114  
#> 
#> Confidence level: 95%

# The planning values can come from a pilot study. The posttest scores in
# depression_bdi have an error variance near 40, and their residuals are
# lighter tailed than the normal, as scores on a bounded inventory often
# are. Thirty residuals estimate the kurtosis imprecisely, so it is a
# planning value to weigh, as the error variance is.
pilot <- lm(bdi_post ~ condition, data = depression_bdi)
sigma(pilot)^2
#> [1] 39.61852
kurtosis(residuals(pilot))
#> [1] -0.9045304
ss_aipe_c(true_var_error = sigma(pilot)^2, c_weights = c(1, -0.5, -0.5),
          width = 3, assurance = .90, true_kurtosis = kurtosis(residuals(pilot)))
#>  term                  value
#>  necessary_n_per_group 110  
#> 
#> Confidence level: 95%
#> 
#> The assurance allows for errors with an excess kurtosis of -0.905 (normal errors have 0).
```
