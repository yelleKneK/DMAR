# Sample Size Planning for the Standardized Mean Difference (AIPE)

Determines the per-group sample size needed for a two-independent-groups
design so that the confidence interval for the population standardized
mean difference, denoted \\\delta\\, is no wider than a user-specified
value when the sample standardized mean difference (Cohen's *d*) equals
\\\delta\\. This is the accuracy in parameter estimation (AIPE) approach
of Kelley and Rausch (2006), the standardized mean difference companion
to power-based planning via
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md).
Optionally, supplying `assurance` returns the smallest sample size, not
below the plan without it, at which the realized interval is no wider
than the target with at least that probability, computed exactly.

## Usage

``` r
ss_aipe_smd(
  true_smd,
  conf_level = 0.95,
  width,
  assurance = NULL,
  assurance_method = "analytic"
)
```

## Arguments

- true_smd:

  The supposed value of the population standardized mean difference
  \\\delta\\ the sample size is planned against: a value the researcher
  posits, either a minimally important effect or a value believed to be
  true in the population, never a sample estimate. Echoed in the
  returned table as the `true_smd` row.

- conf_level:

  Desired confidence level (i.e., \\1-\alpha\\, where \\\alpha\\ is the
  Type I error rate). Default `0.95`. A level above .999999998, whose
  tail areas are below 1e-9, is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- width:

  Desired (full) width of the two-sided confidence interval on
  \\\delta\\.

- assurance:

  Optional probability with which the realized confidence interval is to
  be no wider than `width`. When `NULL` (the default), the plan is for
  the width of the interval at \\d = \delta\\; when supplied (e.g.,
  0.80, 0.90, 0.99), the function returns the smallest per-group sample
  size, not below that plan, at which the probability that the realized
  interval is no wider than `width` is at least `assurance` (see
  Details). Must be `NULL` or strictly between 0.50 and 1. An assurance
  above .999999998 is refused with a condition of class
  `dmar_argument_invalid`: the floor of 1e-9 on a tail area (see "Small
  tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md))
  applies to `(1 - assurance) / 2`.

- assurance_method:

  How the probability that the realized interval is no wider than
  `width` is found, used only with `assurance`. Only `"analytic"` (the
  default) is offered here: the probability is computed exactly, and the
  plan is the smallest sample size, not below the plan without
  `assurance`, at which it reaches `assurance`. `"simulation"`, which
  [`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md),
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
  with `noncentral = TRUE`, and
  [`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)
  offer, where the probability has no exact form, is refused with a
  condition of class `dmar_argument_invalid`, whether or not `assurance`
  is given, since a simulated estimate of an exact probability would
  only add Monte Carlo error.

## Value

A `data.frame` with columns `term` and `value`:

- `necessary_n_per_group`:

  the necessary per-group sample size.

- `true_smd`, `width`:

  the planning inputs, so the assumptions the sample size was planned
  under travel with the result. The `true_smd` row is the supposed
  effect the plan is built on, as given (with its sign): a value the
  researcher posits, either a minimally important effect or a value
  believed to be true in the population, never a sample estimate.

- `assurance`:

  the requested probability; present only when `assurance` was given.

- `assurance_attained`:

  the exact probability, at `necessary_n_per_group` and `true_smd`, that
  the realized interval is no wider than `width`, at least `assurance`;
  present only when `assurance` was given.

The confidence level is reported in the printed footer, and with
`assurance` the footer names the assurance method; the attribute
`assurance_method` records it, and `assurance_probability` records that
`assurance_attained` is `"exact"`.

## Details

**The plan without assurance.** The interval is the one
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)
reports: the confidence interval for the noncentrality parameter of the
noncentral *t* distribution, with \\2n - 2\\ degrees of freedom, at the
observed *t* statistic, rescaled to \\\delta\\ by \\\sqrt{2/n}\\, with
*n* the per-group sample size. The function uses \\\delta\\ as if it
were *d* and returns the smallest *n* at which that interval is no wider
than `width` (Kelley & Rausch, 2006, pp. 371–372). The width depends on
\\\delta\\ only through \\\|\delta\|\\.

**The assurance step.** At a given *n* the width of the interval depends
on the sample only through *d* and grows as \\\|d\|\\ grows, so the
realized interval is no wider than the target exactly when \\\|d\| \le
d^\*(n)\\, the value of \\\|d\|\\ at which the width equals the target.
With \\\lambda = \|\delta\| \sqrt{n/2}\\, \\\nu = 2n - 2\\, and \\t^\* =
d^\*(n) \sqrt{n/2}\\, the probability of that event under normality with
equal variances is \$\$P_n \\=\\ \int_0^\infty \Big\\ \Phi\big(t^\*
\sqrt{v/\nu} - \lambda\big) - \Phi\big(-t^\* \sqrt{v/\nu} - \lambda\big)
\Big\\\\ f\_\nu(v)\\ dv,\$\$ where \\\Phi\\ is the standard normal
distribution function and \\f\_\nu\\ the density of a chi square on
\\\nu\\ degrees of freedom: the *t* statistic is \\(Z + \lambda) /
\sqrt{V/\nu}\\ with \\Z\\ standard normal and \\V\\ chi square,
independent. The integral involves only the normal distribution and the
chi square, so it does not depend on R's noncentral *t*, which is a
normal approximation above a noncentrality of 37.62. It depends on
\\\delta\\ only through \\\|\delta\|\\, and the plan is the same for
`true_smd` and `-true_smd`.

The returned sample size is the smallest *n* at or above the plan
without assurance, \\n_0\\, at which \\P_n\\ is at least `assurance`
(`assurance_method = "analytic"`), found by doubling a step upward from
\\n_0\\ and then bisecting, and `assurance_attained` is \\P_n\\ there.
The search assumes that \\P_n\\ does not decrease as *n* grows from
\\n_0\\. In 45 planning sets (\\\|\delta\|\\ of 0, .2, .5, 1, and 2,
widths of .25, .5, and 1, and assurances of .80, .90, and .99) it did
not decrease at any *n* from \\n_0\\ to ten past the plan, and the plan
was the smallest *n* found by raising *n* one at a time from \\n_0\\.
Because the width changes little with *d*, \\P_n\\ rises steeply with
*n*: at `true_smd = .5` and a width of .5 it is .838 at 129 and .916 at
130.

The probability was checked against the noncentral *t* distribution
computed by R's [`pt()`](https://rdrr.io/r/stats/TDist.html) where it is
exact, against a second integral that conditions on the normal variate
instead of the chi square, and against the proportion of intervals no
wider than the target among samples of two groups of normal data.

**Large noncentrality.** R's
[`pt()`](https://rdrr.io/r/stats/TDist.html) is a normal approximation
to the noncentral *t* distribution above a noncentrality of 37.62, and
its series loses accuracy above 33 with more than 2,000 degrees of
freedom, ranges that narrow targets with a large \\\|\delta\|\\ reach
(`true_smd = 2` with a width of .25, or `true_smd = 1` with a width of
.1). The interval
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) reports
is computed exactly there (see
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)), so
the width of the interval increases with \\\|d\|\\ at every sample size,
the event above is \\\|d\| \le d^\*(n)\\, and \\P_n\\ is exact.

**The published procedure.** Kelley and Rausch (2006, pp. 372–374) plan
with `assurance` as follows: at \\n_0\\, find the value
\\\delta\_\gamma\\ that \\\|d\|\\ exceeds with probability \\1 -
\gamma\\, where \\\gamma\\ is `assurance`, by searching between the
upper limit of a one-sided interval and the upper limit of a two-sided
interval computed at \\d = \delta\\ (their Equations 25 to 31), and plan
for \\\delta\_\gamma\\ as for \\\delta\\. The procedure does not compute
\\P_n\\. Since \\\delta\_\gamma\\ is found at \\n_0\\ and the plan for
it is larger, that plan can be larger than the smallest sample size that
attains the request: in the 45 sets above it attained the request in
every set and was one case larger in 7, all but one at an assurance of
.99; at `true_smd = 1`, a width of .5, and an assurance of .99 it is
150, where \\P_n\\ is .9943, and this function's plan is 149, where it
is .9903. It can also fall short: at `true_smd = 0`, a width of .8, and
an assurance of .95 it is 49 per group, where \\P_n\\ is .9499, and this
function's plan is 50, where it is .9946.

## Warning

The returned value is the sample size *per group*.

## References

Anderson, S. F., & Kelley, K. (2024). Sample size planning for
replication studies: The devil is in the design. *Psychological Methods,
29*(5), 844–867.
[doi:10.1037/met0000520](https://doi.org/10.1037/met0000520)

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*(11), 1547–1562.
[doi:10.1177/0956797617723724](https://doi.org/10.1177/0956797617723724)

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

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

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample-size planning.
*Evaluation & the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`stats::power.t.test()`](https://rdrr.io/r/stats/power.t.test.html)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ss_aipe_smd(true_smd = .5, conf_level = .95, width = .30)
#>  term                  value
#>  necessary_n_per_group 353  
#>  true_smd              0.500
#>  width                 0.300
#> 
#> Confidence level: 95%

# With assurance the last row is the exact probability, at the planned
# sample size, that the interval is no wider than the target.
ss_aipe_smd(true_smd = .5, conf_level = .95, width = .30, assurance = .8)
#>  term                  value
#>  necessary_n_per_group 356  
#>  true_smd              0.500
#>  width                 0.300
#>  assurance             0.800
#>  assurance_attained    0.860
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; assurance_attained is the exact probability that the interval is no wider than the width
ss_aipe_smd(true_smd = .5, conf_level = .95, width = .30, assurance = .95)
#>  term                  value
#>  necessary_n_per_group 359  
#>  true_smd              0.500
#>  width                 0.300
#>  assurance             0.950
#>  assurance_attained    0.966
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; assurance_attained is the exact probability that the interval is no wider than the width

# Reversing the groups reverses the sign of the effect and leaves the
# plan unchanged.
ss_aipe_smd(true_smd = -.5, conf_level = .95, width = .30, assurance = .95)
#>  term                  value 
#>  necessary_n_per_group 359   
#>  true_smd              -0.500
#>  width                 0.300 
#>  assurance             0.950 
#>  assurance_attained    0.966 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; assurance_attained is the exact probability that the interval is no wider than the width
```
