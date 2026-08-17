# Sample Size for AIPE on a Mediated (Indirect) Effect \\ab\\

Determines the sample size needed for the confidence interval on a
mediated effect \\ab\\ (the product of the \\X \to M\\ and \\M \to Y\\
coefficients in a simple three-variable mediation model) to have a
desired full width. Two methods are available. `"closed_form"` (the
default) plans for the symmetric Wald interval built on the delta method
standard error of the product (Sobel, 1982) and answers instantly.
`"monte_carlo"` plans for the Monte Carlo confidence interval
(MacKinnon, Lockwood, & Williams, 2004; Tofighi & MacKinnon, 2011), the
interval
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
reports under `ci_method = "monte_carlo"`, and it plans by a priori
Monte Carlo simulation: at each candidate sample size the mediation
model is fit to `G` simulated data sets and the realized interval widths
are recorded. Because the Monte Carlo interval respects the skewness of
the sampling distribution of a product, and because the simulation
measures the widths that fitted models actually deliver,
`method = "monte_carlo"` is the recommended way to settle on the final
sample size; the closed form is its fast first approximation.

## Usage

``` r
ss_aipe_indirect_effect(
  a,
  b,
  width,
  method = c("closed_form", "monte_carlo"),
  conf_level = 0.95,
  n_max = 10000L,
  B = 5000L,
  G = 1000L,
  seed = NULL
)
```

## Arguments

- a:

  Anticipated population coefficient for \\X \to M\\, on the
  standardized scale. Numeric scalar in \\(-1, 1)\\.

- b:

  Anticipated population coefficient for \\M \to Y\\ controlling for
  \\X\\, standardized. Numeric scalar in \\(-1, 1)\\.

- width:

  Desired full width of the confidence interval on \\ab\\.

- method:

  One of `"closed_form"` (default) or `"monte_carlo"`; see Details.

- conf_level:

  Desired confidence level (default `0.95`).

- n_max:

  Upper bound on the search; default `10000`.

- B:

  Number of Monte Carlo draws forming the interval within each simulated
  study when `method = "monte_carlo"`; default `5000`. This is the same
  `B` the analysis-stage interval uses (see
  [`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)).

- G:

  Number of simulated studies per candidate sample size when
  `method = "monte_carlo"`; default `1000`. The mean simulated width
  carries a simulation error of about its standard deviation over
  \\\sqrt{G}\\; raise `G` for a sharper answer.

- seed:

  Optional integer seed for the Monte Carlo method, used locally (the
  caller's random number generator state is restored on exit). Default
  `NULL` leaves the random number generator state alone.

## Value

A `data.frame` with rows for the recommended sample size, the expected
CI width at that size, and the inputs echoed back. Under
`method = "closed_form"` the expected width is the delta method width
evaluated at the returned sample size; under `method = "monte_carlo"` it
is the mean simulated width there. The method is carried on the returned
object as the `ci_method` attribute.

## Details

**The mediation model.** The simple mediator model is \$\$M = \alpha_1 +
a X + \varepsilon_M,\$\$ \$\$Y = \alpha_2 + c' X + b M +
\varepsilon_Y,\$\$ with the indirect (mediated) effect of \\X\\ on \\Y\\
through \\M\\ equal to \\ab\\ (MacKinnon, Lockwood, Hoffman, West, &
Sheets, 2002). Both methods plan on the standardized scale with no
direct effect: the planning population takes \\X\\, \\M\\, and \\Y\\
with unit variances and \\c' = 0\\, so `a` and `b` are the standardized
paths.

**The closed form.** Under the planning population the sampling variance
of \\\hat a\\ is \\(1 - a^2)/(n - 2)\\. In the equation for \\Y\\ the
mediator is regressed alongside \\X\\, with which it is correlated at
\\a\\, so the sampling variance of \\\hat b\\ carries the variance
inflation factor \\1/(1 - a^2)\\: \$\$\mathrm{Var}(\hat b) \\=\\
\frac{1 - b^2}{(n - 3)(1 - a^2)}.\$\$ The estimators come from two
separate equations and are uncorrelated, so the delta method (Sobel,
1982) standard error of the product is \$\$\mathrm{SE}(\hat a \hat b)
\\=\\ \sqrt{\\a^2 \mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat
a)\\},\$\$ and the closed form returns the smallest \\n\\ at which the
Wald width \\2 z\_{1 - \alpha/2}\\ \mathrm{SE}(\hat a \hat b)\\ is at or
below `width`. Two approximations remain. The Wald interval is symmetric
while the sampling distribution of a product is skewed, so the Wald
interval is not the interval an indirect effect should be reported with
(Tofighi & Kelley, 2020). And the closed form evaluates the standard
error at the planning values, while a fitted model evaluates it at the
estimates, which leaves a discrepancy of a percent or two in realized
width at moderate sample sizes. Both are reasons to treat the closed
form as the first approximation and to verify the final plan with
`method = "monte_carlo"`, which measures the realized widths directly.

**Planning for the Monte Carlo interval.** With
`method = "monte_carlo"`, each candidate \\n\\ is evaluated by a priori
Monte Carlo simulation (Muthén & Muthén, 2002; Schoemann, Boulton, &
Short, 2017): `G` data sets of size \\n\\ are drawn from the planning
population, the two mediation regressions are fit to each, and the Monte
Carlo interval is formed by drawing `B` pairs \\(\tilde a, \tilde b)\\
from normal distributions centered at the estimates with the estimated
standard errors, multiplying, and reading off the empirical \\(\alpha/2,
1 - \alpha/2)\\ quantiles (MacKinnon, Lockwood, & Williams, 2004). Since
\\\hat a\\ and \\\hat b\\ are uncorrelated here, the independent draws
realize the joint normal approximation of the estimates, the same
construction
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
uses for its Monte Carlo interval. The necessary sample size is the
smallest \\n\\ whose mean simulated width is at or below `width`; the
search starts from the closed-form answer, brackets the crossing
geometrically, and bisects. A planning call fits the mediation model
several thousand times and takes a few seconds, which is why the Monte
Carlo example below is shown rather than run. The necessary sample size
inherits the simulation error of the mean widths; raising `G` narrows
it, and supplying `seed` makes a plan reproducible.

**Relation to the MBCO procedure.** The model-based constrained
optimization (MBCO) likelihood ratio test of Tofighi and Kelley (2020),
implemented in
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md),
is the recommended test of a mediation effect, and the intervals that
suit an indirect effect are the profile likelihood interval and the
Monte Carlo interval, both of which accommodate the skewness of the
product. This planner targets the Monte Carlo interval. Planning for the
profile likelihood interval would require inverting a pair of
constrained optimizations in every simulated study (two constrained
OpenMx fits per interval, times `G`, times every candidate sample size),
while the Monte Carlo interval costs `B` products of normal draws per
study and is the inexpensive interval that also accommodates the
skewness, the one Tofighi and Kelley (2020) report for their memory
example. A study planned with `method = "monte_carlo"` and analyzed with
`mediation_mbco(ci_method = "monte_carlo")` is therefore planned and
analyzed on the same interval.

**Beyond the simple model.** The planning population assumes
standardized observed variables, one mediator, no covariates, and no
direct effect. With a nonzero direct effect the residual variance of
\\Y\\ is \\1 - b^2 - c'^2 - 2abc'\\ rather than \\1 - b^2\\, so assuming
\\c' = 0\\ errs toward a larger sample whenever \\c'(c' + 2ab) \> 0\\
(consistent mediation) and toward a smaller one otherwise. When the
direct effect, covariates, several mediators, or latent variables matter
to the design, plan by simulation from the full model with
[`ss_aipe_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
labeling the paths and defining the indirect effect via `ab := a*b`; its
intervals are the Wald intervals of the fitted model, the same target as
the closed form here.
[`ss_aipe_indirect_effect_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md)
quantifies what a plan from this page delivers when the population paths
differ from the planning values.

## References

Fritz, M. S., & MacKinnon, D. P. (2007). Required sample size to detect
the mediated effect. *Psychological Science, 18*(3), 233–239.
[doi:10.1111/j.1467-9280.2007.01882.x](https://doi.org/10.1111/j.1467-9280.2007.01882.x)

Lachowicz, M. J., Preacher, K. J., & Kelley, K. (2018). A novel measure
of effect size for mediation analysis. *Psychological Methods, 23*,
244–261. [doi:10.1037/met0000165](https://doi.org/10.1037/met0000165)

MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
Sheets, V. (2002). A comparison of methods to test mediation and other
intervening variable effects. *Psychological Methods, 7*(1), 83–104.
[doi:10.1037/1082-989X.7.1.83](https://doi.org/10.1037/1082-989X.7.1.83)

MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
limits for the indirect effect: Distribution of the product and
resampling methods. *Multivariate Behavioral Research, 39*(1), 99–128.
[doi:10.1207/s15327906mbr3901_4](https://doi.org/10.1207/s15327906mbr3901_4)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Muthén, L. K., & Muthén, B. O. (2002). How to use a Monte Carlo study to
decide on sample size and determine power. *Structural Equation
Modeling, 9*(4), 599–620.
[doi:10.1207/S15328007SEM0904_8](https://doi.org/10.1207/S15328007SEM0904_8)

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

Schoemann, A. M., Boulton, A. J., & Short, S. D. (2017). Determining
power and sample size for simple and complex mediation models. *Social
Psychological and Personality Science, 8*(4), 379–386.
[doi:10.1177/1948550617715068](https://doi.org/10.1177/1948550617715068)

Sobel, M. E. (1982). Asymptotic confidence intervals for indirect
effects in structural equation models. *Sociological Methodology, 13*,
290–312.

Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
analysis: Introducing the model-based constrained optimization
procedure. *Psychological Methods, 25*, 496–515.
[doi:10.1037/met0000259](https://doi.org/10.1037/met0000259)

Tofighi, D., & MacKinnon, D. P. (2011). RMediation: An R package for
mediation analysis confidence intervals. *Behavior Research Methods,
43*(3), 692–700.
[doi:10.3758/s13428-011-0076-x](https://doi.org/10.3758/s13428-011-0076-x)

## See also

[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
for the analysis the plan feeds;
[`ss_aipe_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md)
for AIPE planning of an indirect effect in an arbitrary lavaan model;
[`ss_aipe_indirect_effect_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md);
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Plan n so the 95% CI on ab has full width <= 0.20, with
#        anticipated standardized a = 0.40 and b = 0.40. The closed
#        form answers instantly:
ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20)
#>  term           value
#>  necessary_N    116  
#>  expected_width 0.2  
#>  a              0.4  
#>  b              0.4  
#>  ab             0.16 
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# 2. The recommended plan targets the Monte Carlo interval directly:
#        every candidate sample size fits the mediation model to G
#        simulated data sets and measures the realized widths. The
#        call takes a few seconds, so it is shown here rather than
#        run. It returns a slightly larger sample size than the
#        closed form because the interval it plans for is a little
#        wider than the Wald interval:
# ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20,
#                         method = "monte_carlo", seed = 113)
```
