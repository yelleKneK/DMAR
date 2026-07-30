# Sample Size Planning for Power in Factorial ANCOVA

Power and sample size for any effect (a main effect or any interaction)
in a between-subjects factorial design with covariates: the analysis of
covariance generalization of
[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md).
Covariates earn their keep by absorbing error variance: with a joint
squared multiple correlation \\R^2\\ between the covariates and the
outcome within cells, the error variance falls by the factor \\1 -
R^2\\, so an effect size \\f\\ defined on the original (ANOVA) metric
grows to \\f / \sqrt{1 - R^2}\\ in the covariate-adjusted analysis, at
the price of one error degree of freedom per covariate.

## Usage

``` r
ss_power_factorial_ancova(
  factor_levels,
  effect_indices,
  f = NULL,
  partial_eta_squared = NULL,
  covariate_R2 = 0,
  n_covariates = 0,
  desired_power = 0.85,
  alpha_level = 0.05,
  n_per_cell = NULL
)
```

## Arguments

- factor_levels:

  Integer vector giving the number of levels of each factor, for example
  `c(2, 4, 3)` for a 2 x 4 x 3 design.

- effect_indices:

  Integer vector identifying the factors that define the effect of
  interest: `1` for the first factor's main effect, `c(2, 3)` for the B
  x C interaction, and so on.

- f:

  Cohen's \\f\\ for the chosen effect *on the unadjusted (ANOVA)
  metric*, that is, with the within-cell standard deviation of the
  outcome in its denominator. Supply this or `partial_eta_squared`, not
  both.

- partial_eta_squared:

  Partial eta squared for the chosen effect on the unadjusted metric.

- covariate_R2:

  Joint squared multiple correlation between the covariates and the
  outcome within cells, in \\\[0, 1)\\. `0` reproduces the ANOVA
  analysis (with the covariate degrees of freedom still spent if
  `n_covariates > 0`).

- n_covariates:

  Number of covariates, a non-negative integer.

- desired_power:

  Desired power; the per-cell sample size is solved when `n_per_cell` is
  `NULL`. Defaults to 0.85 to match
  [`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md).

- alpha_level:

  Type I error rate.

- n_per_cell:

  Per-cell sample size; when supplied, the realized power at that size
  is returned instead of solving for size.

## Value

A `data.frame` with the per-cell and total sample sizes (or the supplied
ones), the realized `actual_power`, the numerator and error degrees of
freedom, the unadjusted and covariate-adjusted effect sizes (`f`,
`f_adjusted`), `covariate_R2`, `n_covariates`, the noncentrality
parameter, and `alpha_level`. The result carries the `dmar_ss_power`
class, so [`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention (the reported size is the per-cell count).

## Details

The test of an effect with numerator degrees of freedom
\\\mathit{df}\_h\\ (the product of the involved factors' levels each
minus one) is a noncentral *F* with noncentrality \\\lambda = N
f\_{\mathrm{adj}}^2\\, where \\N\\ is the total sample size,
\\f\_{\mathrm{adj}} = f / \sqrt{1 - R^2}\\, and error degrees of freedom
\\N - (\prod \mathrm{levels}) - q\\ for \\q\\ covariates (the standard
one-line ANCOVA adjustment; Maxwell, Delaney, & Kelley, 2027, Chapter
9). The covariate slopes are assumed homogeneous across cells and the
covariates measured at baseline, so that adjusting does not bias the
treatment effects in a randomized design.

The complete worked example for this function, a 2 x 4 x 3 ANCOVA with
two baseline covariates, planned effect by effect and then analyzed with
Type III sums of squares, interaction plots, and focused follow-up
contrasts, is the “Power for factorial ANCOVA” vignette:
[`vignette("ancova_2x4x3_power", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/ancova_2x4x3_power.md).

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 on designs with covariates.)

## See also

[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md)
for the no-covariate case this wraps;
[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md) and
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
for the analysis side;
[`vignette("ancova_2x4x3_power", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/ancova_2x4x3_power.md)
for the full worked design.

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_tost_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
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
# A 2 x 4 x 3 design, planning the f = .10 main effect of the
# two-level factor A. Two baseline covariates with a modest joint
# R^2 = .25 cut the required total N by roughly a quarter:
ss_power_factorial_ancova(factor_levels = c(2, 4, 3), effect_indices = 1,
                          f = 0.10, covariate_R2 = 0,    n_covariates = 0,
                          desired_power = 0.80)
#>  term                 value
#>  necessary_n_per_cell 33   
#>  total_N              792  
#>  actual_power         0.803
#>  df_effect            1    
#>  df_error             768  
#>  f                    0.1  
#>  f_adjusted           0.1  
#>  covariate_R2         0    
#>  n_covariates         0    
#>  noncentrality        7.92 
#>  alpha_level          0.05 
ss_power_factorial_ancova(factor_levels = c(2, 4, 3), effect_indices = 1,
                          f = 0.10, covariate_R2 = 0.25, n_covariates = 2,
                          desired_power = 0.80)
#>  term                 value
#>  necessary_n_per_cell 25   
#>  total_N              600  
#>  actual_power         0.806
#>  df_effect            1    
#>  df_error             574  
#>  f                    0.1  
#>  f_adjusted           0.115
#>  covariate_R2         0.25 
#>  n_covariates         2    
#>  noncentrality        8    
#>  alpha_level          0.05 

# Realized power for the three-way interaction at 6 per cell.
ss_power_factorial_ancova(factor_levels = c(2, 4, 3),
                          effect_indices = c(1, 2, 3), f = 0.15,
                          covariate_R2 = 0.25, n_covariates = 2,
                          n_per_cell = 6)
#>  term                 value
#>  specified_n_per_cell 6    
#>  total_N              144  
#>  actual_power         0.276
#>  df_effect            6    
#>  df_error             118  
#>  f                    0.15 
#>  f_adjusted           0.173
#>  covariate_R2         0.25 
#>  n_covariates         2    
#>  noncentrality        4.32 
#>  alpha_level          0.05 
```
