# AIPE Sample Size Planning for an Equivalence Test on the Standardized Mean Difference

Computes the smallest per-group sample size at which the confidence
interval that
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
reports, the 100(1 - 2\\\alpha\\)% interval on the standardized mean
difference whose position against the equivalence bounds decides the two
one-sided tests (Lakens, 2017), is no wider than `width` when the sample
standardized mean difference equals the planning value `true_smd`
(accuracy in parameter estimation, AIPE; Kelley & Rausch, 2006). With
the conventional \\\alpha = .05\\ the interval is the 90% interval. With
`assurance` the function returns the smallest per-group sample size at
which the realized interval is no wider than `width` with at least that
probability.

## Usage

``` r
ss_aipe_equivalence_smd(
  true_smd = 0,
  width,
  alpha_level = 0.05,
  assurance = NULL,
  balanced = TRUE
)
```

## Arguments

- true_smd:

  Planning value of the population standardized mean difference
  \\\delta\\, a single finite number. Default `0`. The interval is
  narrowest when the sample standardized mean difference is 0, so for a
  given width the plan at 0 is the smallest over planning values.

- width:

  Target full width of the interval on the standardized mean difference
  scale (e.g., `0.20` for a 90% interval of width 0.20).

- alpha_level:

  One-sided significance level of each of the two tests. The interval
  planned for is at confidence level \\1 - 2\alpha\\, as in
  [`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md).
  Default `0.05` (90% interval).

- assurance:

  Optional probability in \\(0, 1)\\ that the realized interval is no
  wider than `width`. The function returns the smallest per-group sample
  size, at or above the one planned without it, at which that
  probability is at least `assurance` (Kelley, Maxwell, & Rausch, 2003;
  Kelley & Rausch, 2006). Default `NULL` (the plan is for the width at
  `true_smd` only).

- balanced:

  Logical; `TRUE` (default) plans equal-\\n\\ groups. Unequal-\\n\\
  planning is not yet supported.

## Value

A `data.frame` with columns `term` and `value` and one row per term:

- `necessary_n_per_group`:

  the planned sample size of each group.

- `total_N`:

  twice `necessary_n_per_group`.

- `width`:

  the target width, as given.

- `true_smd`:

  the planning value.

- `expected_width`:

  the full width of the interval
  [`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
  reports with `necessary_n_per_group` in each group when the sample
  standardized mean difference equals `true_smd`, that is, `ci_upper`
  minus `ci_lower` of
  `equivalence_smd(smd = true_smd, n_1 = n, n_2 = n, alpha_level = alpha_level)`
  with `n` the planned size. It is no larger than `width`. Without
  `assurance`, the same interval with one fewer case in each group is
  wider than `width`, unless the plan is the smallest size, 2 per group.

- `assurance`:

  the requested probability; present only when `assurance` was given.

- `assurance_attained`:

  the probability, with `necessary_n_per_group` in each group and
  population standardized mean difference `true_smd`, that the realized
  interval is no wider than `width`; present only when `assurance` was
  given. It is computed, not simulated (see Details), and it is at least
  `assurance`.

The `conf_level` attribute records the confidence level of the interval,
\\1 - 2\alpha\\.

## Details

**The interval planned for.**
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
reports the interval \$\$\hat d \pm t\_{1 - \alpha,\\ n_1 + n_2 - 2}
\sqrt{\frac{n_1 + n_2}{n_1 n_2} + \frac{\hat d^{\\2}}{2 (n_1 +
n_2)}},\$\$ a large-sample standard error of \\\hat d\\ with a critical
value from Student's *t* distribution on \\n_1 + n_2 - 2\\ degrees of
freedom, and it declares equivalence when that interval lies inside the
bounds. This function plans for the same interval, because it is the
interval the analysis reports; a plan for any other interval would
describe a width the analysis does not show. With \\n\\ in each group
the width is \$\$w(\hat d, n) \\=\\ t\_{1 - \alpha,\\ 2n - 2} \sqrt{(8 +
\hat d^{\\2}) / n}.\$\$ The plan is the smallest \\n\\ of at least 2
with \\w(\delta, n) \le \omega\\, where \\\delta\\ is `true_smd` and
\\\omega\\ is `width`. Because \\w\\ decreases in \\n\\, the interval
with \\n - 1\\ in each group is wider than \\\omega\\. The search starts
at \\z\_{1 - \alpha}^2 (8 + \delta^2) / \omega^2\\, the answer with the
standard normal quantile in place of the *t* quantile, which is not
above the answer because the *t* quantile is the larger.

**Assurance.** Without `assurance` the plan makes the interval no wider
than the target when the sample standardized mean difference equals
`true_smd`; a sample with a larger \\\|\hat d\|\\ gives a wider
interval. With `assurance` the function finds the sample size at which
the realized interval is no wider than the target with at least that
probability, from the sampling distribution of \\\hat d\\, as Kelley and
Rausch (2006, pp. 372–373) do for the noncentral *t* interval:

1.  At a given \\n\\ the width depends on the sample only through \\\hat
    d\\ and increases with \\\|\hat d\|\\, so the interval is no wider
    than \\\omega\\ exactly when \\\|\hat d\| \le c(n)\\, where \\c(n)^2
    = n \omega^2 / t\_{1 - \alpha,\\ 2n - 2}^2 - 8\\. When \\c(n)^2 \le
    0\\, even the interval at \\\hat d = 0\\ is wider than \\\omega\\
    and the probability is 0.

2.  For two independent samples from normal populations with a common
    variance, \\T = \hat d \sqrt{n / 2}\\ is the two-group *t*
    statistic, distributed as noncentral *t* on \\\nu = 2n - 2\\ degrees
    of freedom with noncentrality \\\lambda = \delta \sqrt{n / 2}\\. The
    probability is \\P_n = P(\|T\| \le c(n) \sqrt{n / 2})\\, the form of
    Kelley and Rausch (2006, Equation 20). With \\T = (Z + \lambda) /
    \sqrt{V / \nu}\\, where \\Z\\ is standard normal and \\V\\ is chi
    square on \\\nu\\ degrees of freedom, independent of \\Z\\,
    \$\$P(\|T\| \le k) \\=\\ E_V\Big\[\Phi\big(k \sqrt{V / \nu} -
    \lambda\big) - \Phi\big(-k \sqrt{V / \nu} - \lambda\big)\Big\],\$\$
    which is computed by numerical integration over \\V\\. Only the
    normal distribution function and the chi square density enter,
    because R's noncentral *t* distribution function is an approximation
    at noncentrality values above 37.62, which plans for narrow
    intervals reach.

3.  The function returns the smallest \\n\\, at or above the sample size
    planned without `assurance`, with \\P_n \ge\\ `assurance`, found by
    doubling a step upward and then bisecting. `assurance_attained` is
    \\P_n\\ at that \\n\\.

The probability is exact for the normal model, up to numerical
integration. Nothing is simulated, so the function draws no random
numbers and returns the same plan on every call. The probability agreed
with R's noncentral *t* distribution function to within \\10^{-10}\\
wherever that function is accurate, and with simulation from normal
scores: in 24 planning conditions of 20,000 simulated studies each, the
proportion of intervals no wider than the target at the planned size
differed from \\P_n\\ by at most 2.0 Monte Carlo standard errors.

Three properties of the step are worth knowing.

- The answer is never below the sample size planned without `assurance`,
  and an `assurance` at or below the probability at that size returns it
  unchanged. At or above that size \\P_n\\ did not decrease in any of 48
  planning conditions scanned over every whole \\n\\ (6,639 sizes in
  all).

- At `true_smd = 0` the plan without `assurance` often meets a moderate
  request already. The probability at that size depends on how far
  rounding up to a whole number narrows the interval: over widths from
  0.05 to 3 at \\\alpha = .05\\ it ranged from .11 to .97, with a median
  of .85, and for a width of 0.50 it is .839 with 88 per group.

- As \\\|\delta\|\\ grows, the probability at the plan without
  `assurance` falls toward one half, because \\\|\hat d\|\\ exceeds
  \\\|\delta\|\\ about half the time, and the step adds more cases, as
  it does when the target narrows: for `true_smd = 0.5`, a width of
  0.30, and `assurance = 0.99` the plan moves from 249 to 257 per group,
  and for `true_smd = 1`, a width of 0.20, and `assurance = 0.90` from
  610 to 621.

**The noncentral t interval.** Kelley and Rausch (2006) plan for the
interval from the noncentral *t* distribution, which
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)
computes (Kelley, 2007) and
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
plans for; they use the interval with the standard normal quantile only
as the starting value of their search (pp. 371–372). The noncentral *t*
interval has exact coverage for the normal model, where the interval of
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
is a large-sample approximation. At the same \\\hat d\\, \\n\\, and
level the interval of
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
is usually the wider of the two: for 90% intervals with \\\|\hat d\|\\
up to 3, by 18 to 30 percent with 3 per group and by less than 1 percent
with 100 per group. At the 80 percent level it can be the narrower, by
less than 1 percent, when \\\|\hat d\|\\ is 5 or more
([`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
says where the comparison was computed). A study that will report
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) at
level \\1 - 2\alpha\\ is planned with
`ss_aipe_smd(conf_level = 1 - 2 * alpha_level)`.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample-size planning.
*Evaluation & the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Lakens, D. (2017). Equivalence tests: A practical primer for *t* tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*(4), 355–362.
[doi:10.1177/1948550617697177](https://doi.org/10.1177/1948550617697177)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
the analysis planned for;
[`ss_aipe_equivalence_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
which simulates the plan;
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
and [`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
for the noncentral *t* interval.

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
# 1. Plan for a 90% interval on the standardized mean difference of
#    width 0.20, with a planning value of 0.
ss_aipe_equivalence_smd(true_smd = 0, width = 0.20)
#>  term                  value
#>  necessary_n_per_group 543  
#>  total_N               1086 
#>  width                 0.200
#>  true_smd              0    
#>  expected_width        0.200
#> 
#> Confidence level: 90%

# 2. The analysis at the planned size, when the sample standardized mean
#    difference equals the planning value, reports an interval no wider
#    than the target (expected_width is its width).
plan <- ss_aipe_equivalence_smd(true_smd = 0.05, width = 0.20)
n <- plan$value[plan$term == "necessary_n_per_group"]
fit <- equivalence_smd(smd = 0.05, n_1 = n, n_2 = n, smd_upper = 0.3)
fit$value[fit$term == "ci_upper"] - fit$value[fit$term == "ci_lower"]
#> [1] 0.1998534

# 3. With 80% assurance that the realized interval is no wider than 0.20.
#    The last row is the probability at the returned size.
ss_aipe_equivalence_smd(true_smd = 0.05, width = 0.20, assurance = 0.80)
#>  term                  value
#>  necessary_n_per_group 543  
#>  total_N               1086 
#>  width                 0.200
#>  true_smd              0.050
#>  expected_width        0.200
#>  assurance             0.800
#>  assurance_attained    0.870
#> 
#> Confidence level: 90%
```
