# Sample Size or Power for a Standardized Contrast in a One-Way Between-Subjects ANOVA

Determine the necessary per-group sample size to achieve a desired level
of statistical power for the test of a single planned standardized
contrast in a one-way between-subjects analysis of variance, or, given a
per-group sample size, return the realized statistical power.

## Usage

``` r
ss_power_sc(
  psi_standardized,
  c_weights,
  desired_power = 0.85,
  alpha_level = 0.05,
  n = NULL,
  directional = FALSE
)
```

## Arguments

- psi_standardized:

  The population standardized contrast effect, \\\psi / \sigma\\, where
  \\\sigma\\ is the within-group standard deviation

- c_weights:

  Vector of contrast weights (must sum to zero); use fractional weights
  so that the positive weights sum to 1 (e.g.,
  `c(0.5, 0.5, -0.5, -0.5)`)

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- n:

  Per-group sample size (assumed balanced); if specified, returns the
  realized power

- directional:

  Logical: `TRUE` for a one-sided test (in the same sign as
  `psi_standardized`), `FALSE` (default) for a two-sided test

## Value

A `data.frame` with rows for `necessary_n_per_group` (or
`specified_n_per_group`), `actual_power`, and `noncentral_t_parm`. The
result carries the `dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention.

## Details

Under the alternative hypothesis the contrast *t*-statistic follows a
noncentral *t*-distribution with degrees of freedom \\N - J\\ (\\J =
\\`length(c_weights)`) and noncentrality parameter \\\lambda = \psi^\* /
\sqrt{\sum c_j^2 / n}\\, where \\\psi^\*\\ is the standardized contrast.

The function searches over per-group sample sizes \\n\\ until power
first reaches `desired_power`; when `n` is supplied it returns the
realized power.

## References

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_power_c`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md)

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
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Power for a standardized contrast of 0.5 across 4 groups,
# contrast (G1 + G2)/2 vs (G3 + G4)/2, desired power = .80
ss_power_sc(psi_standardized = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
            desired_power = 0.80)
#>  term                  value
#>  necessary_n_per_group 32   
#>  actual_power          0.801
#>  noncentral_t_parm     2.83 

# Realized power at n = 30 per group
ss_power_sc(psi_standardized = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5), n = 30)
#>  term                  value
#>  specified_n_per_group 30   
#>  actual_power          0.775
#>  noncentral_t_parm     2.74 
```
