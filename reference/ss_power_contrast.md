# Sample Size and Statistical Power for a Contrast in a Fixed-Effects ANOVA

Determine the necessary per-group sample size for a contrast in a
one-way fixed-effects ANOVA so as to achieve a desired level of
statistical power, or, alternatively, compute the achieved power for a
given sample size. The contrast is specified by a vector of weights and
the population means (or a population contrast value) along with the
within-group variance.

## Usage

``` r
ss_power_contrast(
  c_weights,
  mu = NULL,
  sigma_squared = NULL,
  psi = NULL,
  desired_power = 0.85,
  alpha_level = 0.05,
  directional = FALSE,
  n_per_group = NULL,
  print_progress = FALSE
)
```

## Arguments

- c_weights:

  Vector of contrast weights. Required to satisfy `sum(c_weights) == 0`,
  `sum(c_weights[c_weights > 0]) == 1`, and
  `sum(c_weights[c_weights < 0]) == -1`; that is, write the contrast in
  normalized fractional form so that the positive coefficients sum to 1
  and the negative coefficients sum to -1. The contrast estimate is then
  directly interpretable as a (weighted) mean difference.

- mu:

  Vector of population group means, of length `length(c_weights)`.
  Either `mu` or `psi` must be supplied (but not both).

- sigma_squared:

  Population within-group variance (\\\sigma^2\\); must be a positive
  number.

- psi:

  Optional. Directly specify the population contrast value \\\psi =
  \sum_j c_j \mu_j\\ as a single number. Use this when the means
  themselves are not of interest, only the contrast value.

- desired_power:

  Target statistical power for the test of the contrast (default
  `0.85`). Used only when `n_per_group` is `NULL`.

- alpha_level:

  Type I error rate (default `0.05`).

- directional:

  Logical. `FALSE` (default) yields a two-sided test. `TRUE` yields a
  one-sided test in the direction of the population contrast.

- n_per_group:

  Optional per-group sample size at which to evaluate power. May be a
  scalar (the common per-group size, equal across groups) or a numeric
  vector of length `length(c_weights)` giving each group's size. When
  `NULL` (default), the function instead solves for the smallest
  per-group *n* achieving `desired_power`.

- print_progress:

  If `TRUE`, print the trial \\n\\ and the corresponding power as the
  iterative search proceeds.

## Value

A `data.frame` with two columns, `term` and `value`:

- When solving for *n* (`n_per_group = NULL`):

  Rows are `necessary_n_per_group`, `total_N`, `actual_power`,
  `noncentral_t_parm`, and `effect_size_f` (Cohen's \\f\\ for a
  one-degree-of-freedom contrast).

- When evaluating power (`n_per_group` supplied):

  Rows are `specified_n_per_group` (or `NA` for unequal *n*), `total_N`,
  `actual_power`, `noncentral_t_parm`, and `effect_size_f`.

The result carries the `dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention; the reported size is the per-group *n* (or `NA`
when unequal group sizes are supplied).

## Details

Let \\\psi = \sum_j c_j \mu_j\\ be the population contrast. Under the
usual fixed-effects ANOVA model with common within-group variance
\\\sigma^2\\ and per-group sizes \\n_j\\, the standard error of the
contrast estimate is \$\$\mathrm{SE}\_{\hat\psi} = \sqrt{\\\sigma^2
\sum_j c_j^2 / n_j\\},\$\$ the test statistic \\t = \hat\psi /
\mathrm{SE}\_{\hat\psi}\\ follows a central *t* distribution with \\N -
a\\ degrees of freedom under \\H_0\\: \psi = 0\\, and a noncentral *t*
distribution with the same df and noncentrality parameter \\\lambda =
\psi / \mathrm{SE}\_{\hat\psi}\\ under the alternative.

For a two-sided test at level \\\alpha\\, \$\$\text{Power} = \Pr(t \>
t\_{1-\alpha/2,\\df}) + \Pr(t \< -t\_{1-\alpha/2,\\df}),\$\$ computed
exactly from the noncentral \\t\\ distribution; for a one-sided test,
only the appropriate tail contributes.

Cohen's \\f\\ for a one-df contrast is reported as \\f = \|\lambda\| /
\sqrt{N}\\, equivalently \\f^2 = F\_{\text{pop}} / N\\ (Cohen, 1988, Ch.
8); this matches the value that `pwr.f2.test` expects when used with
`u = 1` and `v = N - a`.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md)

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
# Four-group example from Maxwell & Delaney's textbook tradition: contrast
# the average of three treatment means with a fourth (control), under
# population means (90, 92, 88, 81), within-group variance 144.
#
# 1. Power achieved at n = 20 per group (total N = 80). Should be ~ .80.
ss_power_contrast(
  c_weights     = c(1/3, 1/3, 1/3, -1),
  mu            = c(90, 92, 88, 81),
  sigma_squared = 144,
  n_per_group   = 20
)
#>  term                  value
#>  specified_n_per_group 20   
#>  total_N               80   
#>  actual_power          0.818
#>  noncentral_t_parm     2.9  
#>  effect_size_f         0.325

# 2. Per-group sample size needed for power = .90.
ss_power_contrast(
  c_weights     = c(1/3, 1/3, 1/3, -1),
  mu            = c(90, 92, 88, 81),
  sigma_squared = 144,
  desired_power = 0.90
)
#>  term                  value
#>  necessary_n_per_group 26   
#>  total_N               104  
#>  actual_power          0.907
#>  noncentral_t_parm     3.31 
#>  effect_size_f         0.325

# 3. Same effect size specification using a directly given psi.
ss_power_contrast(
  c_weights     = c(1/3, 1/3, 1/3, -1),
  psi           = 9,
  sigma_squared = 144,
  n_per_group   = 20
)
#>  term                  value
#>  specified_n_per_group 20   
#>  total_N               80   
#>  actual_power          0.818
#>  noncentral_t_parm     2.9  
#>  effect_size_f         0.325

# 4. Unequal per-group sample sizes.
ss_power_contrast(
  c_weights     = c(0.5, 0.5, -0.5, -0.5),
  mu            = c(90, 92, 88, 81),
  sigma_squared = 144,
  n_per_group   = c(15, 25, 25, 15)
)
#>  term                  value
#>  specified_n_per_group <NA> 
#>  total_N               80   
#>  actual_power          0.639
#>  noncentral_t_parm     2.35 
#>  effect_size_f         0.262
```
