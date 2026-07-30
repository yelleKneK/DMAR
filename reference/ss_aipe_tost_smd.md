# AIPE Sample Size Planning for TOST on the Standardized Mean Difference

Computes the minimum per-group sample size needed so that the
equivalence CI (the 100(1 - 2\\\alpha\\)% CI on Cohen's *d*) has
expected half-width \\\le \omega\\ (Kelley, 2007; Lakens, 2017). With
the standard \\\alpha = 0.05\\ TOST level, the equivalence CI is the 90%
CI. The function inverts the approximate variance of *d* on the
noncentral *t* scale.

## Usage

``` r
ss_aipe_tost_smd(
  population_smd = 0,
  width,
  alpha_level = 0.05,
  assurance = NULL,
  balanced = TRUE
)
```

## Arguments

- population_smd:

  Anticipated population standardized mean difference \\\delta\\ used to
  plan the variance. Default `0`: plans under the most-conservative
  null-effect case.

- width:

  Target full CI width on the *d* scale (e.g., `0.20` for a 90% CI of
  width 0.20).

- alpha_level:

  One-sided TOST significance level. The CI used in planning is at
  confidence level \\1 - 2\alpha\\. Default `0.05` (90% CI).

- assurance:

  Optional assurance probability in \\(0, 1)\\. When supplied, the
  function chooses *n* so that the probability of achieving `width` or
  less is at least `assurance` (Kelley, Maxwell, & Rausch, 2003).
  Default `NULL` (no assurance correction).

- balanced:

  Logical; `TRUE` (default) plans equal-\\n\\ groups. Unequal-\\n\\
  planning is not yet supported.

## Value

A 1-row `data.frame` with the per-group recommended sample size
`necessary_n_per_group`, the implied total `N`, the target `width`, the
planning `delta`, and the resulting `ci_half_width_expected` at the
chosen *n*.

## Details

**Approximate-variance plan.** The large-sample variance of *d* is
\$\$\mathrm{Var}(\hat d) \\\approx\\ (n_1 + n_2) / (n_1 n_2) + d^2 / (2
(n_1 + n_2)).\$\$ For a balanced design with per-group size \\n\\, the
half-width of the equivalence CI at level \\1 - 2\alpha\\ is
approximately \\z\_{1-\alpha} \sqrt{\mathrm{Var}(\hat d)}\\. The
function solves for the smallest integer \\n\\ giving expected
half-width \\\le \omega / 2\\.

**Assurance.** Under `assurance = q`, the function increments \\n\\
until the bootstrap-style probability that the realized half-width is
\\\le \omega / 2\\ is at least \\q\\. (Implemented as a thin Monte-Carlo
overlay; for the most common case `assurance = 0.8` the upward shift is
roughly +1 to +5 per group.)

**Note on conservatism of the assurance plan.** The empirical simulation
reported in
[`vignette("aipe_simulation_study", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/aipe_simulation_study.md)
finds that `ss_aipe_tost_smd()` is tight at \\\gamma = 0.80\\ but
operates on the boundary of its valid range at \\\gamma = 0.99\\: the
realized assurance at the recommended sample size is within Monte-Carlo
error of the target, typically a few tenths of a percentage point below
0.99. The mechanism is that the planner inverts a normal approximation
to \\\Pr(\widehat W \> \omega)\\, and at the 99% level the upper tail of
\\\widehat W\\ is heavier than the approximation accounts for. Adding a
small safety margin (5 to 10 subjects per group) restores the desired
probability statement when planning at high assurance. See the
simulation vignette for the per- condition results.

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
obtaining precision: Delineating methods of sample size planning.
*Evaluation and the Health Professions, 26*(3), 258–287.
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

[`tost_smd`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

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
# 1. Plan for a 90% CI on d of width <= 0.20, under d_planning = 0:
ss_aipe_tost_smd(population_smd = 0, width = 0.20)
#>  term                   value 
#>  necessary_n_per_group  542   
#>  total_N                1084  
#>  width                  0.2   
#>  population_smd         0     
#>  ci_half_width_expected 0.0999

# 2. Plan for the same width assuming a true d = 0.05:
ss_aipe_tost_smd(population_smd = 0.05, width = 0.20)
#>  term                   value 
#>  necessary_n_per_group  542   
#>  total_N                1084  
#>  width                  0.2   
#>  population_smd         0.05  
#>  ci_half_width_expected 0.0999

# 3. With 80% assurance (the assurance path is Monte Carlo, so seed for
#    a reproducible result):
set.seed(113)
ss_aipe_tost_smd(population_smd = 0.05, width = 0.20, assurance = 0.80)
#>  term                   value 
#>  necessary_n_per_group  542   
#>  total_N                1084  
#>  width                  0.2   
#>  population_smd         0.05  
#>  ci_half_width_expected 0.0999
```
