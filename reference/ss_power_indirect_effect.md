# Sample Size Planning for Power for the Indirect (Mediation) Effect

Power, or the sample size required for a desired power, for the test of
the indirect effect \\a b\\ in the simple mediation model, with paths
specified in standardized metric (unit-variance \\X\\, \\M\\, and
\\Y\\). The default test is joint significance (the indirect effect is
declared when both \\\hat a\\ and \\\hat b\\ are individually
significant), which tracks the resampling tests' power closely and far
exceeds the Sobel test in small samples (Fritz & MacKinnon, 2007); the
Sobel normal-theory test is available for comparison. Power comes from
the exact sampling distribution of the two path *t* statistics, with
\\X\\ random by default. This is the power counterpart of the accuracy
in parameter estimation (AIPE) planner
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
and the planning complement of the analysis function
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md).

## Usage

``` r
ss_power_indirect_effect(
  true_std_reg_coef_a,
  true_std_reg_coef_b,
  true_std_reg_coef_direct = 0,
  desired_power = NULL,
  specified_N = NULL,
  alpha_level = 0.05,
  method = c("joint_significance", "sobel"),
  random_predictors = TRUE
)
```

## Arguments

- true_std_reg_coef_a:

  Population standardized regression coefficient of the \\X \to M\\ path
  \\a\\.

- true_std_reg_coef_b:

  Population standardized regression coefficient of the \\M \to Y\\ path
  \\b\\, holding \\X\\.

- true_std_reg_coef_direct:

  Population standardized regression coefficient of the direct effect
  \\c'\\ of \\X\\ on \\Y\\ holding \\M\\. Defaults to 0; it enters only
  through the residual variance of \\Y\\.

- desired_power:

  Desired power; supply this to solve for \\N\\.

- specified_N:

  Total sample size \\N\\; supply this to evaluate the realized power.
  Specify exactly one of `desired_power` and `specified_N`.

- alpha_level:

  Two-sided Type I error rate for each component test. Defaults to 0.05.

- method:

  `"joint_significance"` (default) or `"sobel"`.

- random_predictors:

  `TRUE` (default) when \\X\\ is random, as when it is measured, and
  `FALSE` when it is fixed by design, as when groups are assigned. The
  mediator is random under either setting.

## Value

A tidy `data.frame` with `necessary_N` (or `specified_N`), `power_ab`
(the power to detect the indirect effect \\a b\\, the quantity the
sample size is planned against), the powers of the tests of the two
paths (`power_a`, `power_b`; `NA` for the Sobel method), the paths
(`true_std_reg_coef_a`, `true_std_reg_coef_b`,
`true_std_reg_coef_direct`), their implied standardized indirect effect
`true_std_indirect_effect`, and `alpha_level`. The method is recorded in
the `"method"` attribute and the predictor model in the
`"random_predictors"` attribute. The result carries the `dmar_ss_power`
class, so [`tidy`](https://generics.r-lib.org/reference/tidy.html)
reports the sample size and the power to detect the indirect effect (the
row `power_ab`), and
[`glance`](https://generics.r-lib.org/reference/glance.html) adds the
component powers and the planning inputs.

## Details

With unit-variance variables, \\M = a X + e_M\\ and \\Y = c' X + b M +
e_Y\\, with normal errors of variances \\1 - a^2\\ and \\\sigma^2\_{e_Y}
= 1 - (b^2 + c'^2 + 2 a b c')\\, which must be positive, or the function
stops. As in
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md),
\\t_a\\ tests the slope of \\M\\ on \\X\\ (\\N - 2\\ degrees of freedom)
and \\t_b\\ the coefficient of \\M\\ in the regression of \\Y\\ on \\X\\
and \\M\\ (\\N - 3\\). Given \\X\\, \\t_a\\ is noncentral *t* with
noncentrality \\a \sqrt{SS_X / (1 - a^2)}\\; given \\X\\ and \\M\\,
\\t_b\\ is noncentral *t* with noncentrality \\b \sqrt{SS\_{M \cdot X}}
/ \sigma\_{e_Y}\\, where \\SS\_{M \cdot X}\\ is the residual sum of
squares of \\M\\ on \\X\\. \\SS\_{M \cdot X} / (1 - a^2)\\ is chi square
on \\N - 2\\ degrees of freedom in every design, because \\M\\ is random
given \\X\\; \\SS_X\\ is chi square on \\N - 1\\ with \\X\\ random, and
\\N\\ with \\X\\ fixed at variance 1 (divisor \\N\\; two groups of equal
size, for example). Each power averages the noncentral *t* power over
these chi squares, by Gauss-Hermite sums on their normal scores.

\\SS\_{M \cdot X}\\ is the error sum of squares of one test and the
predictor sum of squares of the other, so the tests are dependent:
`power_ab` is the probability that both reject, which is below `power_a`
times `power_b`. The first-order Sobel statistic of
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md)
equals \\t_a t_b / \sqrt{t_a^2 + t_b^2}\\, and its power is the
probability, integrated numerically, that it exceeds the normal critical
value; `power_a` and `power_b` are then `NA`. Because it refers *t*
statistics to the normal, its power need not rise with \\N\\: with paths
of .9 and .05 it falls from .0729 at \\N = 10\\ to .0562 at 59, and with
paths of .8 and .02 and `alpha_level = .5` it is .5072 at 12 and .5035
at 59, and does not reach .507 again until 216. The joint significance
power rose with \\N\\ in each of 1,000 settings checked at every \\N\\
from 10 to 1,500 (paths from .001 to .99 in absolute value,
`alpha_level` from .0001 to .997, \\X\\ random or fixed), which is not a
proof. So the search assumes neither: `necessary_N` is the smallest
\\N\\ of at least 10 whose power reaches `desired_power`. It scans
upward, passing over a range of \\N\\ only where a bound shows that none
reaches it: written as its quantile at a fixed normal score, each chi
square rises with \\N\\, so the power at every \\N\\ of the range is at
most the power with the noncentralities at the top of the range and each
test's threshold (critical value times scale) at its least over the
range. Where a path is near zero, the power is also at most the rate at
which the test of that path rejects when the path is zero plus
\\\phi(1)\\ (the standard normal density at 1) times the expected square
of its noncentrality (for the Sobel method, the test of \\\|t\|\\
against the normal critical value, since the Sobel statistic is at most
the smaller of \\\|t_a\|\\ and \\\|t_b\|\\), which lets the search pass
at once over a long range in which the power stays just below
`desired_power`. It stops at \\N\\ = 10,000,000 with an error of class
`dmar_computation_failed`. The tests check `power_ab` against raw-data
simulation within Monte Carlo error at the plan for power .80 and at one
case fewer.

Fixing each noncentrality at its population standard error and
multiplying the two powers overstates the joint significance power: with
paths of .39 that approximation gives .8017 at \\N = 65\\, where the
power is .7675, and the plan for .80 is 69.

## References

Fritz, M. S., & MacKinnon, D. P. (2007). Required sample size to detect
the mediated effect. *Psychological Science, 18*(3), 233–239.
[doi:10.1111/j.1467-9280.2007.01882.x](https://doi.org/10.1111/j.1467-9280.2007.01882.x)

MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
Sheets, V. (2002). A comparison of methods to test mediation and other
intervening variable effects. *Psychological Methods, 7*(1), 83–104.
[doi:10.1037/1082-989X.7.1.83](https://doi.org/10.1037/1082-989X.7.1.83)

## See also

[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md) to
analyze the study this plans;
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
to plan for confidence interval width instead of detection;
[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what the chosen design delivers.

Other mediation:
[`mediate()`](https://yelleknek.github.io/DMAR/reference/mediate.md),
[`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_change()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md),
[`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Fritz and MacKinnon's (2007) running scenario (a = b = .39), X random.
ss_power_indirect_effect(true_std_reg_coef_a = .39,
                         true_std_reg_coef_b = .39, desired_power = .80)
#>  term                     value
#>  necessary_N              69   
#>  power_ab                 0.801
#>  power_a                  0.921
#>  power_b                  0.871
#>  true_std_reg_coef_a      0.390
#>  true_std_reg_coef_b      0.390
#>  true_std_reg_coef_direct 0    
#>  true_std_indirect_effect 0.152
#>  alpha_level              0.050

# X fixed by design (two groups of equal size); the mediator is random.
ss_power_indirect_effect(true_std_reg_coef_a = .39,
                         true_std_reg_coef_b = .39, desired_power = .80,
                         random_predictors = FALSE)
#>  term                     value
#>  necessary_N              68   
#>  power_ab                 0.805
#>  power_a                  0.931
#>  power_b                  0.866
#>  true_std_reg_coef_a      0.390
#>  true_std_reg_coef_b      0.390
#>  true_std_reg_coef_direct 0    
#>  true_std_indirect_effect 0.152
#>  alpha_level              0.050

# A near-zero a path against a larger b: the weak link drives the requirement.
ss_power_indirect_effect(true_std_reg_coef_a = .14,
                         true_std_reg_coef_b = .39, desired_power = .80)
#>  term                     value
#>  necessary_N              398  
#>  power_ab                 0.801
#>  power_a                  0.801
#>  power_b                  1.000
#>  true_std_reg_coef_a      0.140
#>  true_std_reg_coef_b      0.390
#>  true_std_reg_coef_direct 0    
#>  true_std_indirect_effect 0.055
#>  alpha_level              0.050

# Realized power at a given N, and the Sobel test at the same N.
ss_power_indirect_effect(true_std_reg_coef_a = .39,
                         true_std_reg_coef_b = .39, specified_N = 75)
#>  term                     value
#>  specified_N              75   
#>  power_ab                 0.844
#>  power_a                  0.940
#>  power_b                  0.898
#>  true_std_reg_coef_a      0.390
#>  true_std_reg_coef_b      0.390
#>  true_std_reg_coef_direct 0    
#>  true_std_indirect_effect 0.152
#>  alpha_level              0.050
ss_power_indirect_effect(true_std_reg_coef_a = .39,
                         true_std_reg_coef_b = .39, specified_N = 75,
                         method = "sobel")
#>  term                     value
#>  specified_N              75   
#>  power_ab                 0.727
#>  power_a                  <NA> 
#>  power_b                  <NA> 
#>  true_std_reg_coef_a      0.390
#>  true_std_reg_coef_b      0.390
#>  true_std_reg_coef_direct 0    
#>  true_std_indirect_effect 0.152
#>  alpha_level              0.050
```
