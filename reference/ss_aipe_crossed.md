# Sample Size Planning From the AIPE Perspective for the Mean of a Crossed Design of Units and Raters

Plans the numbers of units and raters for a fully crossed design, in
which every rater scores every unit once, so that the confidence
interval for the mean score is sufficiently narrow at the least cost.
The units are the objects of measurement (persons, essays, recordings,
or the topics of an information retrieval evaluation) and the raters are
whoever or whatever produces the scores (judges, coders, assessors, or
automated scoring models). When the raters are a sample from a
population of raters, as the units are a sample of units, the interval
for the population mean carries the variability of both samples, and no
number of units can make up for too few raters. The interval planned for
is the one
[`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)
computes from the data.

## Usage

``` r
ss_aipe_crossed(
  true_var_unit,
  true_var_rater,
  true_var_residual,
  width,
  conf_level = 0.95,
  assurance = NULL,
  n_raters = NULL,
  unit_cost = 0,
  rater_cost = 0,
  judgment_cost = 1,
  quantile = c("t", "normal"),
  max_units = 1e+05,
  max_raters = 500
)
```

## Arguments

- true_var_unit, true_var_rater, true_var_residual:

  Planning values of the three variance components: of the unit effects
  (\\\sigma^2_u\\, how much units differ in their mean score over the
  population of raters), of the rater effects (\\\sigma^2_r\\, how much
  raters differ in leniency), and of the residuals (\\\sigma^2_e\\, the
  unit by rater interaction together with error, which one score per
  cell cannot separate). The first two must be nonnegative and the third
  positive. The `var_unit`, `var_rater`, and `var_residual` rows of
  [`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)
  applied to a pilot array estimate them.

- width:

  Desired full width of the confidence interval for the mean.

- conf_level:

  Confidence level of the interval (default `0.95`).

- assurance:

  Optional probability, strictly between 0 and 1, that the realized
  interval is no wider than `width`. The default `NULL` plans for the
  width at the planning values. A plan with an assurance is never
  smaller than the plan without one (see Details), so an assurance at or
  below one half returns the plan without one. An assurance requires
  `quantile = "t"`, since it concerns the sampling variability of the
  estimated variance that the *t* quantile accounts for.

- n_raters:

  Optional fixed number of raters (a whole number of at least 2), for a
  study whose panel of raters is already set; only the number of units
  is then planned. The default `NULL` plans both.

- unit_cost, rater_cost, judgment_cost:

  Costs of adding a unit, of recruiting a rater, and of one judgment
  (one rater's score of one unit), in any common currency. A design with
  \\n_u\\ units and \\n_r\\ raters costs \\c_u n_u + c_r n_r + c_j n_u
  n_r\\. The defaults `0`, `0`, and `1` make the cost the number of
  judgments, so the plan is the design that collects the fewest scores.
  The costs must be nonnegative and not all zero.

- quantile:

  `"t"` (the default) or `"normal"`: the quantile the planned interval
  is formed with (see Details).

- max_units, max_raters:

  The largest numbers of units and of raters the search considers.
  Defaults `100000` and `500`.

## Value

A `data.frame` (class `dmar_ss_aipe`) with rows `necessary_n_units` and
`necessary_n_raters` (the design; the second is `specified_n_raters`
when `n_raters` was supplied), `total_judgments` (\\n_u n_r\\, the
number of scores to collect), `cost`, `expected_width` (the full width
of the interval at the planning values, at most `width`; a plan with an
assurance guards against an unlucky sample, so its expected width is no
wider, and narrower once the assurance asks for more than the plan
without one delivers), `df` (the Satterthwaite degrees of freedom the
interval will have at the planned design, computed from the planning
values), `rater_floor` (see Details), the echoes `true_var_unit`,
`true_var_rater`, `true_var_residual`, `width`, and `conf_level`, and
`assurance` when one was supplied. The quantile the plan used is the
`"quantile"` attribute. `tidy()` reports the number of units as the
sample size, and `glance()` keeps every row.

## Details

**The model and the variance of the mean.** The score of unit \\u\\ by
rater \\r\\ is \\y\_{ur} = \mu + a_u + b_r + e\_{ur}\\, with independent
unit effects, rater effects, and residuals of mean zero and variances
\\\sigma^2_u\\, \\\sigma^2_r\\, and \\\sigma^2_e\\. Both units and
raters are random, so the target is \\\mu\\, the mean over the
population of units and the population of raters, and the grand mean of
a complete array of \\n_u\\ units by \\n_r\\ raters has variance \$\$V =
\frac{\sigma^2_u}{n_u} + \frac{\sigma^2_r}{n_r} + \frac{\sigma^2_e}{n_u
n_r}.\$\$ In the vocabulary of generalizability theory (Brennan, 2001)
the units are the objects of measurement and the raters a random facet
crossed with them.

**The rater floor.** Whatever the number of units, \\V\\ is at least
\\\sigma^2_r / n_r\\, so with \\n_r\\ raters the half-width cannot fall
below \\z\\\sigma_r / \sqrt{n_r}\\, with \\z\\ the normal quantile. A
design meeting a target half-width \\h\\ (half of `width`) under the
normal quantile therefore needs more than \\z^2 \sigma^2_r / h^2\\
raters; that number is the `rater_floor` row. The *t* quantile raises
the requirement further. When raters vary, adding units cannot make up
for a panel that is too small, and a plan with `n_raters` at or below
the floor is refused with that explanation.

**The interval planned for.**
[`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)
estimates \\V\\ by the ANOVA estimates of the three components, a
negative estimate truncated at zero, and refers the grand mean to a *t*
distribution whose degrees of freedom are Satterthwaite's (1946) for
that combination of mean squares. With `quantile = "t"`, the default,
the plan forms the interval the same way at the planning values: the
planned half-width is \\t\_{\nu} \sqrt{V}\\, where the degrees of
freedom \\\nu\\ are those of the plug-in estimate \\(MS_u + MS_r - MS_e)
/ (n_u n_r)\\ computed from the expected mean squares. With
`quantile = "normal"` the planned half-width is \\z \sqrt{V}\\, the
half-width of an interval whose variance is known; `quantile = "normal"`
reproduces planning that uses the normal quantile. Its designs are
smaller, and because the analysis interval uses the *t* quantile, its
realized intervals exceed the target width more often than those of a
plan made with the *t* quantile. The gap is largest when the rater
component is at or near zero. The normal plan then uses two raters,
since raters add nothing to \\V\\, but the analysis still estimates the
rater component from those two. With no rater variance, \\MS_r\\ exceeds
\\MS_e\\ by chance in about a third of samples and the estimate is kept;
the Satterthwaite degrees of freedom then fall, to one or two when the
unit component is small beside the residual one, and the interval can be
several times wider than the target. The *t* plan sees the few degrees
of freedom two raters give at the planning values and adds raters or
units to make up for them.
[`ss_aipe_crossed_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md)
shows how a given plan fares.

**Assurance.** The plan without an assurance gives an interval no wider
than the target at the planning values, and the realized interval is
narrower or wider as the estimated variance falls below or above \\V\\.
An assurance \\\gamma\\ (Kelley & Maxwell, 2003) asks instead for a
realized width at most `width` with probability at least \\\gamma\\. The
plan approximates the distribution of the plug-in variance by a scaled
chi square on the Satterthwaite degrees of freedom, as Satterthwaite
(1946) does, and requires \$\$t\_{\nu} \sqrt{V \max(1, q\_{\gamma} /
\nu)} \le h,\$\$ with \\q\_{\gamma}\\ the \\\gamma\\ quantile of the chi
square distribution on \\\nu\\ degrees of freedom. Under the same
approximation the plan without an assurance already meets the target
with probability at least \\P(\chi^2\_{\nu} \le \nu)\\, which exceeds
one half because the median of a chi square lies below its mean: about
.58 at \\\nu = 5\\, .53 at \\\nu = 30\\, and nearer one half as \\\nu\\
grows. An assurance at or below that probability asks for nothing more,
and \\q\_{\gamma} / \nu\\ is then at most 1; taking the larger of 1 and
\\q\_{\gamma} / \nu\\ keeps such an assurance from shrinking the design.
The rule is that a plan with an assurance is never smaller than the plan
without one: it costs at least as much, and with a fixed panel it has at
least as many units. An assurance at or below one half always returns
the plan without one.

**The search.** For every number of raters from 2 to `max_raters` (or
for `n_raters` alone) the planner finds the smallest number of units
that meets the target, prices each such design, and returns the
cheapest; designs of equal cost are resolved toward fewer raters. With
the *t* quantile the planned width is not monotone in the number of
units: as units are added, \\V\\ approaches its rater floor while the
degrees of freedom eventually fall toward \\n_r - 1\\, and past some
point the larger quantile outweighs the smaller variance, so more units
lengthen the planned interval. The search therefore looks for the
smallest feasible number of units on the stretch where the planned width
is still falling, and it counts a number of raters as infeasible when
the shortest interval that number allows is too wide. A warning says so
when the cheapest design found uses `max_raters` raters, since a larger
limit might find a cheaper one. When no design within the limits meets
the target, the error names the limit that binds: if more units would
do, it gives the design they would allow and points to `max_units`;
otherwise more raters are needed, and it says whether the reason is the
rater floor or the loss of degrees of freedom under the *t* quantile.

**A closed form for the allocation.** The search needs no closed form,
but one special case has one, which serves as a check. When only
judgments cost (`unit_cost = 0`, `rater_cost = 0`), the problem of
minimizing \\V\\ for a fixed number of judgments has the form
Marcoulides and Goldstein (1990) solved for a two-facet generalizability
study, whose error variance is \\\sigma^2_1 / n_1 + \sigma^2_2 / n_2 +
\sigma^2\_{12} / (n_1 n_2)\\ under a budget on \\n_1 n_2\\: the optimum
has \\n_1 / n_2 = \sigma^2_1 / \sigma^2_2\\ (their Equations 5 and 6).
Here that is \\n_u / n_r = \sigma^2_u / \sigma^2_r\\. The cost is flat
near the optimum, so the least costly design in whole numbers can stray
from that ratio while collecting barely more judgments than the
continuous optimum.

**An information retrieval evaluation.** In a comparison of two
retrieval systems over a set of topics, each judged by a set of
relevance assessors (people, or automated judges such as language
models), the units are the topics, the raters are the assessors, and
each cell holds the difference between the two systems' effectiveness
scores on the topic as computed from that assessor's judgments. The
grand mean is then the mean difference between the systems over the
populations of topics and assessors, and `true_var_unit` and
`true_var_rater` describe how much the difference between the systems
varies over topics and over assessors (the topic by system and assessor
by system interactions).

**When the planning values are uncertain.** The plan is only as good as
the planning values, and a rater component that is underestimated is
costly, since it governs the floor.
[`ss_aipe_crossed_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md)
shows by Monte Carlo simulation what a plan delivers when the true
components differ from the planning values. When no trustworthy planning
values exist, a two-stage sequential procedure, in which a pilot array
determines the final numbers of units and raters, can take their place.

## References

Brennan, R. L. (2001). *Generalizability theory*. Springer.
[doi:10.1007/978-1-4757-3456-0](https://doi.org/10.1007/978-1-4757-3456-0)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Marcoulides, G. A., & Goldstein, Z. (1990). The optimization of
generalizability studies with resource constraints. *Educational and
Psychological Measurement, 50*(4), 761–768.
[doi:10.1177/0013164490504004](https://doi.org/10.1177/0013164490504004)

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
[doi:10.2307/3002019](https://doi.org/10.2307/3002019)

## See also

[`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)
for the interval the plan is for,
[`ss_aipe_crossed_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md)
for what a plan delivers when the planning values are wrong, and
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md) for the
reliability of a single rater's scores in the same design.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md),
[`ss_aipe_mean_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md),
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
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Essays scored by raters. Planning values: essays differ with
# variance 0.50, raters differ in leniency with variance 0.10, and the
# residual variance is 0.40. The goal is a 95% interval no wider than
# 0.40 for the mean score; by default the plan collects the fewest
# scores.
ss_aipe_crossed(true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
                width = 0.40)
#>  term               value 
#>  necessary_n_units  94    
#>  necessary_n_raters 22    
#>  total_judgments    2068  
#>  cost               2068  
#>  expected_width     0.400 
#>  df                 72.459
#>  rater_floor        9.604 
#>  true_var_unit      0.500 
#>  true_var_rater     0.100 
#>  true_var_residual  0.400 
#>  width              0.400 
#>  conf_level         0.950 
#> 
#> Confidence level: 95%

# Recruiting and training a rater costs as much as 40 scores, so the
# cheapest design uses fewer raters and more essays.
ss_aipe_crossed(true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
                width = 0.40, rater_cost = 40, judgment_cost = 1)
#>  term               value 
#>  necessary_n_units  112   
#>  necessary_n_raters 19    
#>  total_judgments    2128  
#>  cost               2888  
#>  expected_width     0.399 
#>  df                 53.264
#>  rater_floor        9.604 
#>  true_var_unit      0.500 
#>  true_var_rater     0.100 
#>  true_var_residual  0.400 
#>  width              0.400 
#>  conf_level         0.950 
#> 
#> Confidence level: 95%

# A panel of 30 raters is already set: plan the essays only, with 80%
# assurance that the interval is no wider than 0.40.
ss_aipe_crossed(true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
                width = 0.40, n_raters = 30, assurance = 0.80)
#>  term               value  
#>  necessary_n_units  89     
#>  specified_n_raters 30     
#>  total_judgments    2670   
#>  cost               2670   
#>  expected_width     0.378  
#>  df                 104.006
#>  rater_floor        9.604  
#>  true_var_unit      0.500  
#>  true_var_rater     0.100  
#>  true_var_residual  0.400  
#>  width              0.400  
#>  conf_level         0.950  
#>  assurance          0.800  
#> 
#> Confidence level: 95%
```
