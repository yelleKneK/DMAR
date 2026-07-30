# Sample Size or Composite Power for a One-Way or Factorial ANOVA

Determine the necessary per-cell sample size to achieve a desired level
of composite statistical power in a balanced analysis of variance with
\\a\\ groups (a one-way design) or a factorial arrangement of factors,
or, given a per-cell sample size, return the realized composite power.
Composite power is the probability that every effect named in `effects`
is statistically significant in the same study, the quantity a design
must be planned against when its conclusion requires more than one
result to hold at once. Each effect is a main effect or an interaction,
tested by its own *F* test, and any subset of them can make up the
composite.

## Usage

``` r
ss_power_composite_anova(
  factor_levels,
  effects,
  means = NULL,
  sigma = NULL,
  desired_power = 0.85,
  alpha_level = 0.05,
  n_per_cell = NULL
)
```

## Arguments

- factor_levels:

  Integer vector of the number of levels of each factor, one entry per
  factor (each at least 2). A single value `a` is a one-way design with
  \\a\\ groups; `c(2, 3)` is a 2 by 3 factorial.

- effects:

  A non-empty list naming the effects in the composite. Each element is
  a list with `factors` (a vector of factor indices into
  `factor_levels`: one index for a main effect, several for an
  interaction) and, unless `means` is supplied, exactly one of `f`
  (Cohen's *f*) or `partial_eta_squared`. An optional `label` names the
  effect in the output and the figure; the default label is the factor
  indices joined by `x`. The purported effect sizes are population
  values the researcher posits, never sample estimates.

- means:

  Optional array of population cell means whose dimensions are
  `factor_levels` (a matrix for two factors), or a numeric vector of
  length `prod(factor_levels)` in array order. When supplied, each named
  effect's Cohen's *f* is computed from the means and `sigma`, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  mean pattern.

- sigma:

  The common within-cell population standard deviation of the outcome,
  required with `means` and used only there.

- desired_power:

  Desired composite statistical power (default 0.85). Used only when
  `n_per_cell` is `NULL`.

- alpha_level:

  Type I error rate for each individual *F* test (default 0.05), the
  per-test rate, not a rate for the composite event.

- n_per_cell:

  Per-cell sample size (balanced); if supplied, the realized composite
  power is returned rather than a sample size planned.

## Value

A `data.frame` with `term` and `value` columns: the recommended (or
supplied) `n_per_cell` and total `N`, the `composite_power`, the
`residual_df`, then for each effect its marginal `power_<label>`,
purported `f_<label>`, numerator `df_<label>`, and
`noncentral_parm_<label>`, followed by `cells` and `alpha_level`. The
result carries the `dmar_ss_power` class for
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html), and a
`dmar_composite_power_factorial` class so
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
figure.

## Details

The population effects can be stated two ways: as effect sizes, a
Cohen's *f* or partial eta squared per effect, or as a full array of
population cell means together with a common within-cell standard
deviation, from which each named effect's *f* is read off the analysis
of variance decomposition of the means. Supplying means lets the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method draw the
mean pattern itself, so the size of the population effects is on the
page.

If the design includes one or more covariates, use
[`ss_power_composite_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md).

## References

Maxwell, S. E. (2004). The persistence of underpowered studies in
psychological research: Causes, consequences, and remedies.
*Psychological Methods, 9*, 147–163.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 7 on factorial designs and Chapter 3 on
statistical power.)

## See also

[`ss_power_composite_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md)
for the design with one or more covariates;
[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md)
and
[`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md)
for a single effect

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

Other composite power:
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A 2 by 3 factorial ANOVA whose conclusion needs both main effects, so the
# design is planned against their composite.
ss_power_composite_anova(
  factor_levels = c(2, 3),
  effects = list(list(factors = 1, f = 0.25),
                 list(factors = 2, f = 0.20)),
  desired_power = 0.80)
#>  term                 value
#>  necessary_n_per_cell 43   
#>  necessary_N          258  
#>  composite_power      0.806
#>  residual_df          252  
#>  power_1              0.979
#>  power_2              0.823
#>  f_1                  0.25 
#>  f_2                  0.2  
#>  df_1                 1    
#>  df_2                 2    
#>  noncentral_parm_1    16.1 
#>  noncentral_parm_2    10.3 
#>  cells                6    
#>  alpha_level          0.05 
#>  desired_power        0.8  

# Realized composite power at 40 per cell for a main effect and the
# interaction of a 2 by 2 design, sizes given as partial eta squared.
ss_power_composite_anova(
  factor_levels = c(2, 2),
  effects = list(list(factors = 1,       partial_eta_squared = 0.06),
                 list(factors = c(1, 2), partial_eta_squared = 0.04)),
  n_per_cell = 40)
#>  term                 value
#>  specified_n_per_cell 40   
#>  specified_N          160  
#>  composite_power      0.647
#>  residual_df          156  
#>  power_1              0.888
#>  power_1x2            0.728
#>  f_1                  0.253
#>  f_1x2                0.204
#>  df_1                 1    
#>  df_1x2               1    
#>  noncentral_parm_1    10.2 
#>  noncentral_parm_1x2  6.67 
#>  cells                4    
#>  alpha_level          0.05 

# The effects can instead be a full pattern of population cell means with a
# common within-cell SD; plot() then draws the mean pattern itself.
cell_means <- matrix(c(10, 12, 11,
                       13, 12, 16), nrow = 2, byrow = TRUE)
fit <- ss_power_composite_anova(
  factor_levels = c(2, 3), means = cell_means, sigma = 4,
  effects = list(list(factors = 1, label = "A"),
                 list(factors = 2, label = "B")),
  n_per_cell = 30)
fit
#>  term                 value
#>  specified_n_per_cell 30   
#>  specified_N          180  
#>  composite_power      0.712
#>  residual_df          174  
#>  power_A              0.994
#>  power_B              0.717
#>  f_A                  0.333
#>  f_B                  0.212
#>  df_A                 1    
#>  df_B                 2    
#>  noncentral_parm_A    20   
#>  noncentral_parm_B    8.12 
#>  cells                6    
#>  alpha_level          0.05 
plot(fit)

```
