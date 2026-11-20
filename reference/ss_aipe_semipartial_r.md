# Sample Size for AIPE on a Semipartial (Part) Correlation

Takes a planning value of the population semipartial correlation
\\\rho\_{Y(X \cdot Z_1 \cdots Z\_{p-1})}\\ (the correlation of \\Y\\
with the part of \\X\\ that \\Z_1, \ldots, Z\_{p-1}\\ do not explain,
with \\Y\\ not residualized), the number of predictors (*X* and the
variables partialed out), and a desired confidence interval width, and
returns the smallest sample size at which the interval is expected to be
that narrow (accuracy in parameter estimation, AIPE; Kelley & Maxwell,
2003). With `assurance` it returns, by simulation, the smallest sample
size at which the realized interval is that narrow with about the stated
probability (the search and the probability reported are both Monte
Carlo estimates; see Details).

## Usage

``` r
ss_aipe_semipartial_r(
  true_r_sp,
  p,
  width,
  which_width = c("Full", "Lower", "Upper"),
  conf_level = 0.95,
  assurance = NULL,
  true_R2_full = NULL,
  G = 20000,
  seed = NULL
)
```

## Arguments

- true_r_sp:

  Population semipartial correlation the sample size is planned for, in
  \\(-1, 1)\\. The planner treats this value as the truth;
  [`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)
  shows what the plan delivers when the truth (`true_r_sp` there)
  differs from the value planned with (`supposed_r_sp` there).

- p:

  Number of predictors in the regression of *Y* on *X* and the variables
  partialed out, so `p - 1` variables are partialed out and `p` is at
  least 2. A categorical variable with *a* levels counts as *a* - 1
  predictors. It is the `p` of
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
  [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md), and
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md):
  the correlation concerns one predictor of a regression with `p`
  predictors. In the notation of a model comparison, \\q\\ is the number
  of predictors removed from the full model to give the restricted model
  (\\q = p\_{full} - p\_{restricted}\\, the predictors being tested),
  and the number of variables partialed out is \\p - q\\; this
  correlation concerns one predictor, so \\q = 1\\ and \\p - 1\\
  variables are partialed out.

- width:

  Desired full width of the confidence interval on the semipartial
  correlation.

- which_width:

  Whether `width` refers to the `"Full"` width (default) or the
  `"Lower"`/`"Upper"` half-width.

- conf_level:

  Desired confidence level (default `0.95`).

- assurance:

  Optional. The probability that the realized interval is no wider than
  `width`. The function returns the smallest sample size, at or above
  the sample size planned without it, at which the simulated probability
  is at least `assurance` (see "The assurance step" in Details). A
  single number above 0.5 and below 1; when `NULL` (the default) the
  plan is for the expected width only and nothing is simulated.

- true_R2_full:

  Population squared multiple correlation of the full model the sample
  size is planned for, the model with \\X\\ and every variable partialed
  out (`true_R2_Y_X` in
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)),
  treated as the truth as `true_r_sp` is, in \\\[0, 1\]\\ and not
  smaller than `true_r_sp` squared (a smaller value is refused with an
  error of class `dmar_semipartial_r_incoherent`, unless the pair could
  have come from rounding to two decimals, when the plan proceeds as if
  `true_R2_full` were `true_r_sp` squared and gives a warning of class
  `dmar_semipartial_r_rounded`; see
  [`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)).
  When `NULL` (the default), the plan uses the largest variance over
  every full-model \\R^2\\ that is possible with `true_r_sp`, so the
  sample size is sufficient whatever that \\R^2\\ is.

- G:

  Number of simulated samples at each sample size the assurance step
  tries, and in the final estimate of the attained probability; a whole
  number of at least 100 (default 20000). A smaller value is refused,
  because the step estimates a probability from these samples and with
  fewer than 100 its standard error can exceed .05. Used and checked
  only when `assurance` is given.

- seed:

  Optional single number. When given, the assurance step is
  reproducible, and the random number stream of the calling session is
  left as it was. With `NULL` (the default) the draws come from the
  session's current stream. Used only when `assurance` is given.

## Value

A `data.frame` with the columns `term` and `value` and one row per term:

- `necessary_N`:

  the planned total sample size.

- `expected_width`:

  the full width of the interval at `necessary_N` when the two sample
  values equal the planning values.

- `true_r_sp`, `p`, `width`, `conf_level`:

  the inputs, with `width` the `width` as given.

- `true_R2_full`:

  the planning value of the full-model \\R^2\\; `NA` when none was
  supplied.

- `assurance`:

  the requested probability; present only when `assurance` was given, as
  are the next two rows.

- `assurance_attained`:

  the proportion of `G` fresh simulated samples of size `necessary_N`
  whose interval was no wider than the target. Its Monte Carlo standard
  error is the attribute `assurance_mc_se`.

- `G`:

  the number of simulated samples.

## Details

**Estimand and interval.** The quantity planned for is the population
semipartial correlation of \\Y\\ with \\X\\ given \\Z_1, \ldots,
Z\_{p-1}\\, for independent observations from a multivariate normal
population. The interval is the two-sided Wald interval at `conf_level`,
the estimate plus and minus \\z\_{1 - \alpha/2} \sqrt{\hat V / (N -
p)}\\, where \\\hat V\\ is the first-order variance below evaluated at
the sample semipartial correlation and the sample full-model \\R^2\\
([`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)
with `under_null = FALSE`). Its limits are not truncated at \\-1\\ and
\\1\\.

**Asymptotic variance of the semipartial.** Under multivariate
normality, to first order, the sample semipartial correlation has
variance \\V / N\\, with \$\$V \\=\\ (1 - \rho\_{sp}^2)^2 \\-\\
P_r^2\\(1 - 2 \rho\_{sp}^2),\$\$ where \\\rho\_{sp}\\ is the population
semipartial correlation (`true_r_sp`) and \\P_r^2\\ is the full-model
\\R^2\\ (`true_R2_full`) minus \\\rho\_{sp}^2\\, that is, the squared
multiple correlation of \\Y\\ with the variables partialed out. This is
the first-order variance of Aloe and Becker (2012, Equation 4), which
they obtain from Alf and Graf (1999), with \\N - p\\ in place of their
\\n\\: in their notation the numerator is \\P^4 - 2 P^2 + P\_{(f)}^2 +
1 - P\_{(f)}^4\\, with \\P^2\\ the full-model and \\P\_{(f)}^2\\ the
restricted-model squared multiple correlation, which equals \\V\\ with
\\P_r^2 = P\_{(f)}^2\\. It follows from the delta method variance of the
difference between two nested squared multiple correlations (Model A of
Olkin & Finn, 1995, for one variable partialed out; Alf & Graf, 1999,
for any number; see
[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)).
The planner uses the divisor \\N - p\\, as
[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)
does. Inverting for the sample size needed to achieve a target
half-width \\w\_{1/2}\\ at confidence level \\1 - \alpha\\: \$\$N \\=\\
p + \Big\lceil z\_{1 - \alpha/2}^{2} \cdot V / w\_{1/2}^{2}
\Big\rceil.\$\$

**Planning without `true_R2_full`.** \\V\\ decreases in \\P_r^2\\ when
\\\rho\_{sp}^2 \< 1/2\\ and increases in \\P_r^2\\ when \\\rho\_{sp}^2
\> 1/2\\, so its largest value over every possible full-model \\R^2\\ is
\$\$\max\\(1 - \rho\_{sp}^2)^2,\\ \rho\_{sp}^2 (1 - \rho\_{sp}^2)\\,\$\$
and the planner uses that value when `true_R2_full` is `NULL`. For
\\\rho\_{sp}^2 \le 1/2\\ it is \\(1 - \rho\_{sp}^2)^2\\, the variance
when the variables partialed out explain none of the outcome; Algina and
Olejnik (2003, p. 317) note, for the squared semipartial correlation,
that the required sample size is largest when the squared multiple
correlation of the restricted model is zero. A planning value of
`true_R2_full` then gives a smaller sample size. For \\\rho\_{sp}^2 \>
1/2\\ the largest value is \\\rho\_{sp}^2 (1 - \rho\_{sp}^2)\\,
approached as the full-model \\R^2\\ approaches 1.

**Comparison with partial-r planning.** The partial correlation \\r\_{XY
\cdot Z}\\ divides the covariance after residualizing both \\X\\ and
\\Y\\ on \\Z\\; the semipartial divides after residualizing only \\X\\.
The semipartial is the natural effect size companion to a standardized
regression coefficient: its square equals the \\\Delta R^2\\ contributed
by \\X\\ above and beyond the variables partialed out. See
[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)
for the asymptotic variance, and
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md)
for the partial-correlation analog of this function.

**Limits of the plan.** The plan is a first-order one. When the returned
sample size is small, which happens when \\\|\rho\_{sp}\|\\ is large and
`true_R2_full` is near `true_r_sp` squared, the first-order variance
understates the variance of the estimate and the realized intervals are
wider than planned; check such a plan with
[`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md).
The planned sample size is also small when `true_R2_full` is near 1, for
any `true_r_sp`, and the same check applies: the plan for
`true_r_sp = .1`, `p = 4`, a width of .20, and `true_R2_full = .99` is
12, at which the intervals missed the population value in .13 of 200,000
simulated samples. When \\\|\rho\_{sp}\|\\ is above about .9 and the
variables partialed out explain some of the outcome, the interval misses
more often on one side than the other even with several hundred
observations (in one simulated condition, \\\rho\_{sp} = .95\\, a
full-model \\R^2\\ of .97, \\p = 2\\, and \\N = 395\\, it missed in .064
of samples, .051 of them on one side). A plan with `true_r_sp = 0` and
`true_R2_full = 1` is refused, because the first-order variance is zero
there and no sample size follows from it.

**The assurance step.** A plan without `assurance` makes the interval as
narrow as the target when the two sample values equal their planning
values; in roughly half of all samples it is wider. With `assurance` the
function looks for the sample size at which the realized interval is no
wider than the target with that probability. The width of the interval
depends on two statistics, the sample semipartial correlation and the
sample full-model \\R^2\\, and no closed form for their joint
distribution is used here, so the probability is simulated:

1.  *What is simulated.* Samples of size \\N\\ from a multivariate
    normal population with semipartial correlation `true_r_sp` and
    full-model \\R^2\\ `true_R2_full`: \\X\\, \\Z_1, \ldots, Z\_{p-1}\\,
    and an error \\E\\ are independent standard normal, and \\Y =
    \rho\_{sp} X + P_r\\ Z_1 + \sqrt{1 - \rho\_{sp}^2 - P_r^2}\\ E\\.
    Under normality the distribution of the two statistics depends on
    the population only through \\\rho\_{sp}\\, the full-model \\R^2\\,
    \\p\\, and \\N\\ (see
    [`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)),
    so this population stands for all of them. The sums of squares and
    products of the \\p + 1\\ independent variables are Wishart with
    \\N - 1\\ degrees of freedom, and the function draws the few entries
    of that matrix's triangular factor that the two statistics depend on
    (the Bartlett decomposition) in place of the data.

2.  *The interval.* From each draw, the sample semipartial correlation
    and sample full-model \\R^2\\ give the half-width \\z\_{1 -
    \alpha/2} \sqrt{\hat V / (N - p)}\\, the interval that is planned
    and that
    [`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)
    forms from data.

3.  *The search.* The simulated probability at \\N\\ is the share of `G`
    half-widths no larger than the target half-width. Starting from the
    sample size planned without `assurance`, the step upward doubles
    until that share reaches `assurance`, and a bisection finds the
    first sample size at which it does, with fresh draws at every sample
    size tried. The answer is never below the plan without `assurance`.

4.  *The attained probability.* The share that chose \\N\\ is biased
    upward by the choice, so `assurance_attained` is estimated again at
    the returned \\N\\ from `G` draws that had no part in choosing it.
    Its standard error, \\\sqrt{A (1 - A) / G}\\ with \\A\\ that
    proportion, is the attribute `assurance_mc_se`.

Plans made this way were checked from data with
[`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)
in 11 conditions (5,000 simulated data sets each; a Monte Carlo standard
error of .003 to .006): the proportion of intervals no wider than the
target was .798 to .819 for a requested .80, .854 to .858 for .85, .906
to .908 for .90, and .953 for .95.

The step has properties a closed-form plan does not.

- The answer has Monte Carlo error. Without `seed`, two calls can return
  sample sizes that differ by a few cases (eight unseeded plans at the
  default `G` for one condition gave 252 to 254, and 40 seeded plans of
  each of two other conditions gave 470 to 473 and 306 to 310). The
  probability at the smallest of those sizes was below the request (.897
  and .895 for a requested .90), so a plan is the size at which the
  probability is about, and not at least, the request. Give `seed` for a
  plan that can be reproduced, and report it with `G`.

- `assurance_attained` can fall below `assurance` by a few of its
  standard errors, since it is estimated from draws that did not choose
  \\N\\. Its standard error is .003 at the default `G` and an assurance
  of .80; in 200 seeded plans the lowest value was .009 below a request
  of .90. It is an unbiased estimate of the probability at the sample
  size returned. A larger `G` narrows both this and the spread of the
  sample sizes.

- Without `true_R2_full` the simulated population is the one the plan
  without `assurance` is made for: variables partialed out that explain
  none of the outcome (\\P_r^2 = 0\\) when `true_r_sp` squared is at
  most 1/2, and a full-model \\R^2\\ of 1 (\\P_r^2 = 1 - \rho\_{sp}^2\\,
  no error term) when it is above 1/2. The second is the limit that
  gives the largest first-order variance;
  [`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)
  refuses a population \\R^2\\ of exactly 1, so check such a plan there
  with a value just below 1 or give `true_R2_full`.

- `which_width = "Lower"` or `"Upper"` treats `width` as a half-width,
  in the assurance step as in the plan.

- `G` and `seed` are checked and used only when `assurance` is given;
  without it nothing is drawn and the three rows are absent.

## References

Alf, E. F., Jr., & Graf, R. G. (1999). Asymptotic confidence limits for
the difference between two squared multiple correlations: A simplified
approach. *Psychological Methods, 4*(1), 70–75.
[doi:10.1037/1082-989X.4.1.70](https://doi.org/10.1037/1082-989X.4.1.70)

Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
analysis with applications in partial correlation and multiple
regression analysis. *Multivariate Behavioral Research, 38*(3), 309–323.
[doi:10.1207/s15327906mbr3803_02](https://doi.org/10.1207/s15327906mbr3803_02)

Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
predictors in meta-analysis. *Journal of Educational and Behavioral
Statistics, 37*(2), 278–297.
[doi:10.3102/1076998610396901](https://doi.org/10.3102/1076998610396901)

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on the one-way ANOVA and Chapter 4 on
contrasts.)

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin, 118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

## See also

[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
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
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Plan N so the 95% CI on r_sp (three variables partialed out, so
#        p = 4 predictors) has full width <= 0.15 when the planning
#        value of the semipartial is 0.25, whatever the full-model R^2.
ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15)
#>  term           value
#>  necessary_N    605  
#>  expected_width 0.150
#>  true_r_sp      0.250
#>  p              4    
#>  width          0.150
#>  conf_level     0.950
#>  true_R2_full   <NA> 
#> 
#> Confidence level: 95%

# 2. With a planning value of 0.42 for the full-model R^2, fewer are
#        needed.
ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15,
                      true_R2_full = 0.42)
#>  term           value
#>  necessary_N    391  
#>  expected_width 0.150
#>  true_r_sp      0.250
#>  p              4    
#>  width          0.150
#>  conf_level     0.950
#>  true_R2_full   0.420
#> 
#> Confidence level: 95%

# 3. With an assurance of .80, found by simulation. G = 2000 keeps the
#        example quick; a reported plan deserves the default G = 20000.
#        The last three rows are the request, the proportion of fresh
#        simulated samples that met the width at the returned N, and G.
plan <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15,
                              assurance = 0.80, G = 2000, seed = 113)
plan
#>  term               value
#>  necessary_N        621  
#>  expected_width     0.148
#>  true_r_sp          0.250
#>  p                  4    
#>  width              0.150
#>  conf_level         0.950
#>  true_R2_full       <NA> 
#>  assurance          0.800
#>  assurance_attained 0.804
#>  G                  2000 
#> 
#> Confidence level: 95%
attr(plan, "assurance_mc_se")
#> [1] 0.008876486
```
