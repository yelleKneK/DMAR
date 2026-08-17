# Sample Size for Equivalence or Noninferiority of a Linear Contrast

Computes the smallest per-group sample size at which the two one-sided
tests procedure (Schuirmann, 1987) for a linear contrast \\\psi = \sum_j
c_j \mu_j\\, or the companion one-sided noninferiority test, attains a
desired power. Power is computed exactly through
[`power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md).
This is the declaration-probability route to planning; the accuracy in
parameter estimation (AIPE) route, which targets the confidence interval
width directly, is
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
and the two answer the same question whenever `true_psi = 0` and the
width target is calibrated to the bounds.

## Usage

``` r
ss_power_equivalence_c(
  c_weights,
  sigma,
  delta_lower = NULL,
  delta_upper = NULL,
  true_psi = 0,
  desired_power = 0.85,
  alpha_level = 0.05,
  side = c("equivalence", "noninferiority")
)
```

## Arguments

- c_weights:

  The contrast weights. The weights must sum to zero with the positive
  weights summing to 1 and the negative weights to -1, so that the
  bounds are on the raw scale of the response.

- sigma:

  The anticipated error standard deviation (the square root of the mean
  square error).

- delta_lower, delta_upper:

  Equivalence bounds on the raw scale of the response. Both must be
  positive; the equivalence region is \\(-\delta_L, +\delta_U)\\. If
  only `delta_upper` is supplied, the bounds are symmetric.
  Noninferiority uses \\-\delta_L\\ alone.

- true_psi:

  The population contrast the design should be able to detect as
  equivalent (or noninferior). Default `0`. For equivalence it must lie
  strictly inside the bounds; for noninferiority, strictly above
  \\-\delta_L\\. The farther `true_psi` sits from the center of the
  region, the larger the required sample size.

- desired_power:

  The target probability of declaring equivalence (or noninferiority) at
  `true_psi`. Default `0.85`, matching
  [`ss_power_contrast`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md).

- alpha_level:

  One-sided significance level for each test. Default `0.05`.

- side:

  `"equivalence"` (default) or `"noninferiority"`.

## Value

A `data.frame` with rows `necessary_n_per_group` (the recommended sample
size for each of the \\J\\ groups named by `c_weights`), `total_N` (the
implied total, \\J \times n\\), and `actual_power` (the exact power
achieved at the recommendation). The result carries the `dmar_ss_power`
class, so [`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention.

## Details

**Design.** Planning assumes equal allocation across the \\J\\ groups
named by `c_weights` and a pooled error term on \\N - J\\ degrees of
freedom. Groups with zero weight still contribute error degrees of
freedom, which is why they belong in `c_weights` when the fitted model
will include them.

**The search.** The function starts from the normal-theory approximation
and moves to the smallest integer \\n\\ whose exact power reaches
`desired_power`. Power is monotone in \\n\\ once the design is feasible,
so the search is a short walk.

**Relation to the half-width rule.** With symmetric bounds
\\\pm\delta\\, `true_psi = 0`, and \\\alpha = .05\\, targeting a 90% CI
half-width of \\\delta/2\\ yields a declaration probability of about
.90; this function makes the probability the target directly rather than
through the width.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`ss_power_contrast`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
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
# 1. Two-group equivalence with bounds of 5 raw-scale points and an
#    anticipated error SD of 15.67: n per group for 90% power at
#    true equivalence.
ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                       delta_upper = 5, desired_power = 0.90)
#>  term                  value
#>  necessary_n_per_group 214  
#>  total_N               428  
#>  actual_power          0.901

# 2. Noninferiority is cheaper than equivalence at the same bound.
ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                       delta_upper = 5, desired_power = 0.90,
                       side = "noninferiority")
#>  term                  value
#>  necessary_n_per_group 169  
#>  total_N               338  
#>  actual_power          0.9  

# 3. A true contrast off center raises the requirement.
ss_power_equivalence_c(c_weights = c(1, -1), sigma = 15.67,
                       delta_upper = 5, true_psi = 2,
                       desired_power = 0.90)
#>  term                  value
#>  necessary_n_per_group 468  
#>  total_N               936  
#>  actual_power          0.9  
```
