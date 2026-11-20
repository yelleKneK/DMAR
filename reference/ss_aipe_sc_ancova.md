# Sample Size Planning From the AIPE Perspective for Standardized ANCOVA Contrasts

Sample size planning from the accuracy in parameter estimation (AIPE)
perspective for standardized ANCOVA contrasts.

## Usage

``` r
ss_aipe_sc_ancova(
  true_contrast = NULL,
  true_sd_anova = NULL,
  true_sd_ancova = NULL,
  true_std_contrast = NULL,
  true_sd_ratio = NULL,
  true_cor = NULL,
  divisor = "s_ancova",
  c_weights,
  width,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  assurance = NULL,
  assurance_method = c("analytic", "simulation"),
  G = 10000,
  seed = NULL
)
```

## Arguments

- true_contrast:

  The population unstandardized ANCOVA (adjusted) contrast

- true_sd_anova:

  The population error standard deviation of the ANOVA model

- true_sd_ancova:

  The population error standard deviation of the ANCOVA model

- true_std_contrast:

  The population standardized ANCOVA (adjusted) contrast

- true_sd_ratio:

  The ratio of `true_sd_ancova` over `true_sd_anova`, greater than 0 and
  at most 1

- true_cor:

  The population correlation coefficient between the response and the
  covariate

- divisor:

  Which error standard deviation to be used in standardizing the
  contrast; the value can be either `"s_ancova"` or `"s_anova"`

- c_weights:

  Contrast weights

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage (i.e., 1 - Type I error
  rate). Default is `.95`, which gives a symmetric two-sided interval,
  the interval
  [`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
  reports. A level above .999999998, whose tail areas are below 1e-9, is
  refused.

- alpha_lower, alpha_upper:

  The lower-tail and upper-tail Type I error rates, an alternative to
  `conf_level` (give `conf_level = NULL` with them).
  [`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
  reports intervals with equal tail areas only, so the two must be
  equal, and `alpha_lower = alpha_upper = a` plans the interval of
  `conf_level = 1 - 2a`. Unequal values are refused with a condition of
  class `dmar_argument_invalid`: they would plan for an interval no DMAR
  function reports. Each must be at least 1e-9: below it the noncentral
  *t* limits lose accuracy (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- assurance:

  Parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be `NULL` or greater than .5 and less than 1). An assurance
  above .999999998 is refused, with either `assurance_method`: the
  published assurance step finds noncentral *t* limits with tail areas
  of `(1 - assurance) / 2`, which must be at least 1e-9.

- assurance_method:

  How the plan with `assurance` is made, used only with `assurance`.
  `"analytic"` (the default) is the published plan of Lai and Kelley
  (2012), computed without simulation, which does not compute the
  probability it attains and can attain less than `assurance`;
  `assurance_attained` is then `NA`. `"simulation"` returns the smallest
  sample size per group, at or above the plan without `assurance`, at
  which the probability that the
  [`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
  interval is no wider than `width` is at least `assurance`, with that
  probability estimated from `G` simulated samples and reported as
  `assurance_attained` (see Details).

- G:

  Number of simulated samples at each sample size the assurance step
  tries, and in its estimate of the attained probability; a whole number
  of at least 100 (default 10000, which gives a Monte Carlo standard
  error of .004 for a probability of .80). Used and checked only with
  `assurance` and `assurance_method = "simulation"`.

- seed:

  Optional single number. When given, the assurance step is
  reproducible, and the random number stream of the calling session is
  left as it was. With `NULL` (the default) the draws come from the
  session's current stream. Used only where `G` is.

## Value

A `data.frame` with columns `term` and `value`: `necessary_n_per_group`,
the per-group sample size needed for the planned ANCOVA contrast; with
`assurance`, the rows `assurance` (the request) and
`assurance_attained`: `NA` for the published plan of
`assurance_method = "analytic"`, which does not compute it, and for
`assurance_method = "simulation"` the simulated probability at
`necessary_n_per_group` that the
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
interval is no wider than `width`, with its Monte Carlo standard error
in the attribute `assurance_mc_se` and the row `G`. With `assurance` the
printed footer names the assurance method and says how
`assurance_attained` was found, or that the published plan does not
compute it; the attributes `assurance_method` and
`assurance_probability` (`"simulation"` or `"not computed"`) record the
same.

## Details

The sample size planning method this function is based on is developed
in the context of simple (i.e., one-response-one-covariate) ANCOVA model
and randomized design (i.e., same population covariate mean across
groups).

An ANCOVA contrast can be standardized in at least two ways: (a) divided
by the error standard deviation of the ANOVA model, (b) divided by the
error standard deviation of the ANCOVA model. This function can be used
to analyze both types of standardized ANCOVA contrasts.

Not all of the effect size arguments need to be specified. When
`divisor="s_ancova"` the input is either (a) `true_std_contrast`, or (b)
`true_contrast` (the unstandardized ANCOVA contrast) and
`true_sd_ancova`. When `divisor="s_anova"`, the valid input combinations
are (a) `true_std_contrast` and `true_sd_ratio`; (b) `true_std_contrast`
and `true_cor`; or (c) `true_contrast`, `true_sd_anova`, and
`true_sd_ancova`. With `divisor = "s_anova"` the ratio of the two
standard deviations, \\\sqrt{1 - \rho^2}\\, must be greater than 0 and
at most 1, as it is in every population (a covariate cannot raise the
error variance), so `true_cor` of 1 or -1 is refused there.

**The interval planned for.** The plan is for the interval of
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md):
the limits for the noncentrality parameter of a noncentral *t* on \\N -
J - 1\\ degrees of freedom at the observed *t* of the adjusted contrast,
carried back to the scale of the standardized contrast by the standard
error \\\sqrt{C/n + (\sum_j c_j \bar{X}\_j)^2 / SS\_{W,X}}\\, with \\C =
\sum_j c_j^2\\, \\\bar{X}\_j\\ the covariate means, and \\SS\_{W,X}\\
the within-group sum of squares of the covariate (and, for the divisor
`"s_anova"`, by the ratio \\s\_{ANCOVA}/s\_{ANOVA}\\). The plan without
`assurance` is for the criterion that this interval, evaluated at the
planning value, is no wider than `width` (Lai & Kelley, 2012, Equations
21 and 28); it uses \\\sigma\_{ANCOVA}\sqrt{C/n}\\ as the standard error
of the contrast, leaving out the covariate means term, which in a
randomized design is small (Lai & Kelley, 2012, pp. 354–355). The plan
is the smallest *n* per group that meets the criterion, searched from 2
per group, the smallest sample size with positive error degrees of
freedom. The search doubles a step and then bisects, which finds the
smallest *n* because the interval at the planning value narrows as *n*
grows from 2 (checked at every *n* up to 80, and in steps of 10 up to
600, for 420 combinations of the contrast, the ratio
\\\sigma\_{ANCOVA}/\sigma\_{ANOVA}\\, weights for two to five groups,
and the confidence level). At `true_std_contrast = 0`,
`c_weights = c(1, -1)`, and `width = .5` the plan is 123 per group,
where the interval is .49985 wide (.50190 at 122), and at
`true_std_contrast = .5`, `true_cor = .5`, `divisor = "s_anova"`,
`c_weights = c(1, -1)`, and `width = 3.6` it is 2 per group, where the
interval is 3.596 wide. A `width` that would need more than \\2^{53}\\
(about \\9.0 \times 10^{15}\\, the largest number up to which R holds
every whole number exactly) per group is refused with a condition of
class `dmar_argument_invalid`. The intervals have equal tail areas, as
those of
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
do and as Lai and Kelley (2012, p. 354) assume.

Past about \\10^{11}\\ per group the plan and `assurance_attained` are
limited by the precision of doubles, as for
[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md):
the width changes by about \\1/(2n)\\ of itself from one *n* to the
next, and for a nonzero contrast it is the difference of two limits near
the noncentrality. Against a 512-bit computation of the interval, with
`c_weights = c(1, -1)`, `conf_level = .95`, and the divisor
`"s_ancova"`, and with `width` and the tail areas as R holds them (each
tail area, (1 - .95) / 2, is 0.025000000000000022 in doubles): at
`true_std_contrast = .5` and `width` of 1e-5, 1e-6, and 1e-7 the plan
was 4, 1,035, and 525,379 above the smallest *n* that meets the
criterion (about 3.2e11, 3.2e13, and 3.2e15 per group). At a contrast of
0, where the width is computed from the normal quantiles, the plan was
the smallest *n* at `width` of 1e-3, 1e-5, and 1e-6 (about 3.1e7,
3.1e11, and 3.1e13 per group) and 2 below it at 1e-7 (about 3.1e15).
Past about \\10^{15}\\ per group the plan at a contrast of 0 can be
below the smallest *n*, by up to about \\10^{-15} n\\: at `width = 6e-8`
it is 8,536,575,157,098,048 per group, 8 below the smallest *n*,
8,536,575,157,098,056. With `assurance = .9`,
`assurance_method = "simulation"`, `G = 200`, and `seed = 1`,
`assurance_attained` was the share of the same simulated statistics
whose interval, computed in 512 bits, is no wider than `width`, up to
3.1e11 per group; at about 3.1e13 per group it was .950 where that share
is .940 (contrast 0) and .900 where it is .895 (contrast .5), and at
about 3.1e15 it was 1.000 where it is 0 (the plan is below the smallest
*n*) and .905 where it is .930.

**The assurance step.** With `assurance` the plan is by default
(`assurance_method = "analytic"`) the published one (Lai & Kelley, 2012,
pp. 358–359, a generalization of Kelley & Rausch, 2006): the value
\\t^\*\\ that \\\|T\|\\ does not exceed with probability `assurance`
under the planner's model, carried back to the scale of the standardized
contrast, is planned for without assurance. It is computed without
simulation and does not compute the probability it attains, so
`assurance_attained` is `NA`, and the printed footer says so. That model
leaves out the covariate means term and treats
\\s\_{ANCOVA}/s\_{ANOVA}\\ as known, and the widths vary little, so the
plan attains less than the request: at `true_std_contrast = .5`,
`true_cor = .5`, `c_weights = c(1, -1)`, `width = .5`, and
`assurance = .80`, the published plans (129 per group with
`divisor = "s_ancova"`, 98 with `"s_anova"`) give an interval no wider
than .5 with probability .77 and .55.

With `assurance_method = "simulation"` the plan is the smallest sample
size per group, not below the plan without `assurance`, at which the
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
interval is no wider than `width` with at least that probability, for
normal within-group errors and a normal covariate with the same
distribution in every group (the randomized design). The width depends
on the sample through five independent statistics: the contrast of the
covariate means, the within-group sum of squares of the covariate, the
contrast of the group means of the errors, the error of the pooled
within-group slope, and the ANCOVA error sum of squares. The probability
is estimated from `G` draws of these statistics from their exact
sampling distributions, the same draws at every sample size the search
tries, so that the estimate changes smoothly as the sample size grows,
and both the covariate means term and, with `divisor = "s_anova"`, the
sampling of \\s\_{ANCOVA}/s\_{ANOVA}\\ enter it. `assurance_attained` is
the estimate at the returned sample size, at least `assurance`, with its
Monte Carlo standard error in the attribute `assurance_mc_se`. Give
`seed` for a plan that can be reproduced. For the divisor `"s_anova"`,
Lai and Kelley (2012, pp. 360–361) wrote that their assurance method is
not appropriate, because the width varies with
\\s\_{ANCOVA}/s\_{ANOVA}\\ as well as with *t*, and gave as a practical
solution an a priori Monte Carlo simulation that raises the sample size
until the proportion of intervals no wider than the target reaches the
assurance; this step is that simulation, made from the sufficient
statistics rather than from raw data, and applied to both divisors. In
the example above the simulated plans (130 and 104 per group,
`seed = 1`) give an interval no wider than .5 with probability .87 and
.84. Because the probability is estimated, a plan can differ by one from
the smallest sample size at which the probability itself reaches
`assurance`: with `"s_anova"` above, 103 gives .803.
[`ss_aipe_sc_ancova_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md)
simulates raw data at a plan, and the attained probability of any plan
can be checked there.

## References

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

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
[`ss_aipe_sc_ancova_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ss_aipe_sc_ancova(true_std_contrast = .8, width = .5, c_weights = c(.5, .5, 0, -1))
#>  term                  value
#>  necessary_n_per_group 98   
#> 
#> Confidence level: 95%

ss_aipe_sc_ancova(true_std_contrast = .8, true_sd_ratio = .6, width = .5,
                  c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#>  term                  value
#>  necessary_n_per_group 39   
#> 
#> Confidence level: 95%

ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .4, width = .3,
               c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#>  term                  value
#>  necessary_n_per_group 221  
#> 
#> Confidence level: 95%

# With assurance, the published plan of Lai and Kelley (2012), which
# does not compute the probability it attains.
ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, width = .5,
                  c_weights = c(1, -1), assurance = .90)
#>  term                  value
#>  necessary_n_per_group 130  
#>  assurance             0.900
#>  assurance_attained    <NA> 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; the plan is the published approximation, which does not compute the probability that the interval is no wider than the width (assurance_attained is NA) and can attain less than the assurance. assurance_method = "simulation" computes it and returns the smallest plan that attains it

# By simulation: the smallest sample size per group at which the
# interval is no wider than .5 with probability .90. G = 2000 keeps the
# example quick; a reported plan deserves the default G = 10000.
plan <- ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, width = .5,
                          c_weights = c(1, -1), assurance = .90,
                          assurance_method = "simulation",
                          G = 2000, seed = 113)
plan
#>  term                  value
#>  necessary_n_per_group 131  
#>  assurance             0.900
#>  assurance_attained    0.930
#>  G                     2000 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: simulation; assurance_attained is estimated from G = 2000 simulated samples (Monte Carlo SE 0.006)
attr(plan, "assurance_mc_se")
#> [1] 0.005705261
```
