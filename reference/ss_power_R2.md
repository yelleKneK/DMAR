# Plan Sample Size to Make the Test of the Squared Multiple Correlation Coefficient Sufficiently Powerful

Determine the necessary sample size for the omnibus test of the squared
multiple correlation coefficient (\\R^2\\), or the realized statistical
power given a specified sample size, under either fixed or random
predictors. The fixed-predictors path uses Cohen's (1988) noncentral *F*
formulation; the random-predictors path uses the three-moment noncentral
*F* approximation of Lee (1971, Section 5.1) to the sampling
distribution of the sample \\R^2\\ under joint multivariate normality.

## Usage

``` r
ss_power_R2(
  true_R2 = NULL,
  alpha_level = 0.05,
  desired_power = 0.85,
  p,
  specified_N = NULL,
  true_f2 = NULL,
  null_R2 = 0,
  random_predictors = TRUE,
  print_iter = FALSE
)
```

## Arguments

- true_R2:

  Population squared multiple correlation coefficient the sample size or
  power is planned for. The planner treats this value as the truth;
  [`ss_power_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md)
  shows what the plan delivers when the truth (`true_R2` there) differs
  from the value planned with (`supposed_R2` there).

- alpha_level:

  Type I error rate

- desired_power:

  Desired degree of statistical power, a single number strictly between
  0 and 1. The values 0 and 1 are refused, with a condition of class
  `dmar_desired_power_invalid`.

- p:

  The number of predictor variables

- specified_N:

  The sample size used to calculate power (rather than determine
  necessary sample size). This is the *total* sample size across all
  groups or observations.

- true_f2:

  Cohen's (1988) effect size for multiple regression: `true_R2`/(1 -
  `true_R2`).

- null_R2:

  The null value of the omnibus test, which must be 0. Any other value
  stops with an error of class `dmar_null_R2_unsupported` that explains
  why and gives the call to
  [`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md);
  see Details.

- random_predictors:

  Whether the predictor variables are treated as random (`TRUE`, the
  default) or fixed (`FALSE`). See Details.

- print_iter:

  If the progress of the iterative procedure is printed to the screen as
  the iterations are occurring

## Value

A `data.frame` with columns `term` and `value`. For an \\N\\ search the
rows are `necessary_N`, `actual_power`, `noncentral_f_parm` (only
meaningful for `random_predictors = FALSE`; `NA` otherwise), and
`true_f2` (Cohen's \\f^2 = \rho^2 / (1 - \rho^2)\\ at `true_R2`). For a
power-at-specified-`N` computation the first row is `specified_N`
instead of `necessary_N`.

## Details

Determine the necessary sample size given a particular `true_R2`,
`alpha_level`, `p`, and `desired_power`. Alternatively, given `true_R2`,
`alpha_level`, `p`, and `specified_N`, the function can be used to
determine the statistical power.

**Fixed vs.\\ random predictors.** The two regression models give
*different* sampling distributions for the omnibus \\F\\-statistic, and
so different power. Under fixed predictors the design matrix is treated
as constant in hypothetical replications of the study, and \\F\\ follows
a noncentral *F* with \\p\\ and \\N-p-1\\ degrees of freedom and
noncentrality \\\lambda = N \cdot f^2\\, where \\f^2 = \rho^2 / (1 -
\rho^2)\\ and \\\rho^2\\ is `true_R2` (Cohen, 1988). Cohen's
noncentrality is \\\lambda = f^2 (u + v + 1)\\, with \\u\\ and \\v\\ the
numerator and denominator degrees of freedom (as given by Maxwell, 2000,
Equation 1); for the omnibus test \\u = p\\ and \\v = N - p - 1\\, so
\\u + v + 1 = N\\ and the two forms are the same. Under random
predictors the design matrix is itself a draw from a joint multivariate
normal distribution, and the unconditional distribution of the sample
\\R^2\\ is due to Fisher (1928); `ss_power_R2()` uses the three-moment
noncentral *F* approximation of Lee (1971, Section 5.1) to that
distribution, the same approximation
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) uses
for random-predictor confidence intervals. Gatsonis and Sampson (1989)
give the exact calculation for random predictors and compare it with
Cohen's tables (the 1977 edition of the book cited below), whose sample
sizes were slightly smaller than the exact ones for most designs with up
to 10 predictors. The sample size this function returns under random
predictors is the smallest at which the power by Lee's approximation
attains `desired_power`. Against the exact distribution, the power by
Lee's approximation differed by at most .0043 over the combinations
examined (squared correlations up to .95, `alpha_level` from .001 to
.10, and one or more error degrees of freedom) and by at most .0004 with
25 or more error degrees of freedom; the necessary sample size differed
from the one the exact distribution gives, by one case, in 20 of 1,500
searches, and in those the power by the exact distribution at the size
returned was within .0011 of `desired_power`. In the behavioral,
educational, and social sciences predictor variables are almost always
random, so the default is `random_predictors = TRUE`; pass
`random_predictors = FALSE` for designs in which the predictor variables
are fixed by design (for example, planned dosing levels).

**Tests other than the omnibus test.** `ss_power_R2()` plans the test
that the population squared multiple correlation coefficient of all `p`
predictors is zero. The test that predictors added to a model contribute
nothing, a test of the change in \\R^2\\, is planned with
[`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md),
which takes the number of predictors in both models.

`null_R2` exists so that a call written for MBESS with a nonzero null
value is answered with an explanation; it must be 0. A nonzero value
stops with an error of class `dmar_null_R2_unsupported` (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)),
for these reasons. In MBESS a nonzero null value \\\rho_0^2\\ forms the
effect size \\(\rho^2 - \rho_0^2) / (1 - \rho^2)\\, which is Cohen's
(1988) \\f^2\\ for the increment in \\R^2\\ when predictors are added to
a model that already explains \\\rho_0^2\\, and refers it to an *F*
distribution with \\p\\ and \\N - p - 1\\ degrees of freedom and
noncentrality \\N f^2\\. What that calculation is depends on how `p` is
read. With `p` read as the number of predictors added, it is Cohen's
calculation, \\\lambda = f^2 (u + v + 1)\\, for a study whose total
sample size is \\N\\ plus the number of predictors in the restricted
model: the sample size it reports is Cohen's total sample size less that
number. With `p` read as all of the predictors, as this function defines
it, the calculation is not the test of the increment, whose numerator
degrees of freedom are the number of predictors added and not `p`. The
argument alone cannot say which calculation is meant, because the test
of an increment needs the number of predictors in the restricted model,
so the value is refused and the message writes the call to
[`ss_power_R2_change()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md).

The number that MBESS returns for a population value \\a\\ and a null
value \\b\\ is
`ss_power_R2(true_f2 = (a - b) / (1 - a), p = p, random_predictors = FALSE)`,
the omnibus test of `p` predictors with an effect of that size; the last
example shows it. Whether the population squared multiple correlation
coefficient of all `p` predictors exceeds a stated value above zero is a
third question. DMAR does not plan sample size for that test; the test
itself rejects when the one-sided lower confidence limit from
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) (with
`alpha_lower = alpha_level` and `alpha_upper = 0`) exceeds the stated
value.

**Limit of the calculation.** R computes the noncentral *F* distribution
accurately for a noncentrality parameter up to about one million, and
warns beyond it. `ss_power_R2()` returns no power from an evaluation at
which R gave that warning: a power requested at such a `specified_N` is
refused, and a search whose answer lies there stops and says so. A
search that only passes through such sizes on its way to a smaller
answer returns that answer, which is evaluated without the warning. The
search does not go beyond a sample size of ten million, and stops with
an error when `desired_power` is not attained by then. With `true_R2` of
0 the null hypothesis is true: the power at a `specified_N` is
`alpha_level`, and a search stops with an error of class
`dmar_effect_zero` because no sample size attains `desired_power`.

## Note

When determining sample size for a desired degree of power, there will
always be a slightly larger degree of actual power. This is the case
because the algorithm employed determines sample size until the actual
power is no less than the desired power (given sample size is a whole
number power will almost certainly not be exactly the specified value).
This is the same as other statistical power procedures that return whole
numbers for necessary sample size.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Fisher, R. A. (1928). The general sampling distribution of the multiple
correlation coefficient. *Proceedings of the Royal Society of London,
Series A, 121*(788), 654–673.
[doi:10.1098/rspa.1928.0224](https://doi.org/10.1098/rspa.1928.0224)

Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106*(3), 516–524.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Lee, Y. S. (1971). Some results on the sampling distribution of the
multiple correlation coefficient. *Journal of the Royal Statistical
Society, Series B, 33*(1), 117–130.

Maxwell, S. E. (2000). Sample size and multiple regression analysis.
*Psychological Methods, 5*(4), 434–458.
[doi:10.1037/1082-989X.5.4.434](https://doi.org/10.1037/1082-989X.5.4.434)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*(11), 1547–1562.
[doi:10.1177/0956797617723724](https://doi.org/10.1177/0956797617723724)

## See also

[`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md)
for the test of predictors added to a model,
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`ss_power_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
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
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
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
# Random predictors (default; appropriate for most behavioral / social
# science applications).
ss_power_R2(true_R2 = .5, alpha_level = .05, desired_power = .85, p = 5)
#>  term              value
#>  necessary_N       23   
#>  actual_power      0.852
#>  noncentral_f_parm <NA> 
#>  true_f2           1    

# Fixed predictors (Cohen 1988): predictor variables fixed by design.
ss_power_R2(true_R2 = .5, alpha_level = .05, desired_power = .85,
            p = 5, random_predictors = FALSE)
#>  term              value 
#>  necessary_N       21    
#>  actual_power      0.858 
#>  noncentral_f_parm 21.000
#>  true_f2           1     

# Effect size input (Cohen's f^2).
ss_power_R2(true_f2 = 1, alpha_level = .05, desired_power = .85, p = 5)
#>  term              value
#>  necessary_N       23   
#>  actual_power      0.852
#>  noncentral_f_parm <NA> 
#>  true_f2           1    

# Realized power at a specified N.
ss_power_R2(true_R2 = .5, specified_N = 15, alpha_level = .05,
            desired_power = .85, p = 5)
#>  term              value
#>  specified_N       15   
#>  actual_power      0.543
#>  noncentral_f_parm <NA> 
#>  true_f2           1    

# The number MBESS reported for a population value of .30 against a null
# value of .25 with p = 4: the omnibus test of four fixed predictors with
# Cohen's f squared for that increment. The test of the increment itself
# is planned with ss_power_R2_change().
ss_power_R2(true_f2 = (.30 - .25) / (1 - .30), p = 4, specified_N = 85,
            alpha_level = .01, random_predictors = FALSE)
#>  term              value
#>  specified_N       85   
#>  actual_power      0.223
#>  noncentral_f_parm 6.071
#>  true_f2           0.071
```
