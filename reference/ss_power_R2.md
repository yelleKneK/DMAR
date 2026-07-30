# Plan Sample Size to Make the Test of the Squared Multiple Correlation Coefficient Sufficiently Powerful

Determine the necessary sample size for the omnibus test of the squared
multiple correlation coefficient (\\R^2\\), or the realized statistical
power given a specified sample size, under either fixed or random
predictors. The fixed-predictors path uses Cohen's (1988) noncentral *F*
formulation; the random-predictors path uses the Lee (1971) two-moment
approximation to the sampling distribution of the sample \\R^2\\ under
joint multivariate normality.

## Usage

``` r
ss_power_R2(
  population_R2 = NULL,
  alpha_level = 0.05,
  desired_power = 0.85,
  p,
  specified_N = NULL,
  cohen_f2 = NULL,
  null_R2 = 0,
  random_predictors = TRUE,
  print_progress = FALSE,
  ...
)
```

## Arguments

- population_R2:

  Population squared multiple correlation coefficient

- alpha_level:

  Type I error rate

- desired_power:

  Desired degree of statistical power

- p:

  The number of predictor variables

- specified_N:

  The sample size used to calculate power (rather than determine
  necessary sample size). This is the *total* sample size across all
  groups or observations.

- cohen_f2:

  Cohen's (1988) effect size for multiple regression:
  `population_R2`/(1-`population_R2`)

- null_R2:

  Value of the null hypothesis that the squared multiple correlation
  will be evaluated against (this will typically be zero)

- random_predictors:

  Whether the predictor variables are treated as random (`TRUE`, the
  default) or fixed (`FALSE`). See Details.

- print_progress:

  If the progress of the iterative procedure is printed to the screen as
  the iterations are occurring

- ...:

  Possible additional parameters for internal functions

## Value

A `data.frame` with columns `term` and `value`. For an \\N\\ search the
rows are `necessary_N`, `actual_power`, `noncentral_f_parm` (only
meaningful for `random_predictors = FALSE`; `NA` otherwise), and
`effect_size` (Cohen's \\f^2\\). For a power-at-specified-`N`
computation the first row is `specified_N` instead of `necessary_N`.

## Details

Determine the necessary sample size given a particular `population_R2`,
`alpha_level`, `p`, and `desired_power`. Alternatively, given
`population_R2`, `alpha_level`, `p`, and `specified_N`, the function can
be used to determine the statistical power.

**Fixed vs.\\ random predictors.** The two regression models give
*different* sampling distributions for the omnibus \\F\\-statistic, and
so different power. Under fixed predictors the design matrix is treated
as constant in hypothetical replications of the study, and \\F\\ follows
a noncentral *F* with \\p\\ and \\N-p-1\\ degrees of freedom and
noncentrality \\\lambda = N \cdot f^2\\, where \\f^2 = \rho^2 / (1 -
\rho^2)\\ (Cohen, 1988). Under random predictors the design matrix is
itself a draw from a joint multivariate normal distribution, and the
unconditional distribution of the sample \\R^2\\ is given by Lee (1971);
`ss_power_R2()` uses Lee's two-moment Patnaik (1949) approximation to
that distribution, the same approximation
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) uses
for random-predictor confidence intervals. Gatsonis and Sampson (1989)
document the comparison and show that Cohen's fixed-predictor formula
tends to over-state power (and so under-state required \\N\\) relative
to the random model; the discrepancy is modest at moderate to large
\\N\\ but non-trivial for small \\N\\ with moderate-to-large effects. In
the behavioral, educational, and social sciences predictor variables are
almost always random, so the default is `random_predictors = TRUE`; pass
`random_predictors = FALSE` for designs in which the predictor variables
are fixed by design (for example, planned dosing levels).

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

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Patnaik, P. B. (1949). The non-central \\\chi^2\\- and *F*-distributions
and their applications. *Biometrika, 36*(1–2), 202–232.
[doi:10.1093/biomet/36.1-2.202](https://doi.org/10.1093/biomet/36.1-2.202)

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*(11), 1547–1562.
[doi:10.1177/0956797617723724](https://doi.org/10.1177/0956797617723724)

## See also

[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`ss_power_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_tost_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md),
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
ss_power_R2(population_R2 = .5, alpha_level = .05, desired_power = .85, p = 5)
#>  term              value
#>  necessary_N       23   
#>  actual_power      0.852
#>  noncentral_f_parm <NA> 
#>  effect_size       1    

# Fixed predictors (Cohen 1988): predictor variables fixed by design.
ss_power_R2(population_R2 = .5, alpha_level = .05, desired_power = .85,
            p = 5, random_predictors = FALSE)
#>  term              value
#>  necessary_N       21   
#>  actual_power      0.858
#>  noncentral_f_parm 21   
#>  effect_size       1    

# Effect size input (Cohen's f^2).
ss_power_R2(cohen_f2 = 1, alpha_level = .05, desired_power = .85, p = 5)
#>  term              value
#>  necessary_N       23   
#>  actual_power      0.852
#>  noncentral_f_parm <NA> 
#>  effect_size       1    

# Realized power at a specified N.
ss_power_R2(population_R2 = .5, specified_N = 15, alpha_level = .05,
            desired_power = .85, p = 5)
#>  term              value
#>  specified_N       15   
#>  actual_power      0.543
#>  noncentral_f_parm <NA> 
#>  effect_size       1    
```
