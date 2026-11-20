# Sample Size or Power for a Pearson Correlation Coefficient (Fisher Z Transformation)

Determine the necessary sample size to achieve a desired level of
statistical power for the test of a Pearson correlation against a null
value (typically zero), or, given a sample size, return the realized
statistical power. The computation uses the Fisher's *Z* transformation,
which has a near-normal sampling distribution.

## Usage

``` r
ss_power_r(
  true_cor,
  null_cor = 0,
  desired_power = 0.85,
  alpha_level = 0.05,
  specified_N = NULL,
  directional = FALSE
)
```

## Arguments

- true_cor:

  Population correlation coefficient the sample size or power is planned
  for (the value under the alternative hypothesis), in \\(-1, 1)\\. The
  planner treats this value as the truth.

- null_cor:

  The null hypothesis value of the correlation (default 0)

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- specified_N:

  Total sample size (the *number of pairs*); if specified, the function
  returns the power at that sample size.

- directional:

  Logical: `TRUE` for a one-sided test (in the same direction as the
  difference `true_cor - null_cor`), `FALSE` (default) for a two-sided
  test

## Value

A `data.frame` with rows for `necessary_N` (or `specified_N`),
`actual_power`, `true_cor`, and `null_cor`.

## Details

Under the alternative the Fisher-transformed correlation \\Z_r =
\tanh^{-1}(r)\\ is approximately normal with mean \\Z\_\rho =
\tanh^{-1}(\rho)\\ and variance \\1 / (N - 3)\\. Power is computed from
this normal approximation.

For sample size, a closed-form expression is used as the starting point,
\$\$N = ((z\_{\alpha} + z\_{\beta}) / (Z\_\rho - Z\_{\rho_0}))^2 +
3,\$\$ where \\\rho\\ is `true_cor`, \\\rho_0\\ is `null_cor`,
\\z\_{\alpha}\\ is the standard normal critical value at `alpha_level`
(one- or two-sided, as `directional` says), and \\z\_{\beta}\\ is the
standard normal quantile at `desired_power`, which is then verified
iteratively to ensure power exactly meets or exceeds `desired_power`.
The search is bounded at \\N = 10^7\\. When `true_cor` and `null_cor`
are so close that `desired_power` is unreachable within that bound, the
function stops with an error rather than searching indefinitely.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on the one-way ANOVA and Chapter 4 on
contrasts.)

## See also

[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_Z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

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
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
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
# true_cor = 0.30, null_cor = 0, desired power = .80, two-sided
ss_power_r(true_cor = 0.30, desired_power = 0.80)
#>  term         value
#>  necessary_N  85   
#>  actual_power 0.800
#>  true_cor     0.300
#>  null_cor     0    

# Same with a directional alternative
ss_power_r(true_cor = 0.30, desired_power = 0.80, directional = TRUE)
#>  term         value
#>  necessary_N  68   
#>  actual_power 0.802
#>  true_cor     0.300
#>  null_cor     0    

# Realized power for N = 100 pairs
ss_power_r(true_cor = 0.30, specified_N = 100)
#>  term         value
#>  specified_N  100  
#>  actual_power 0.862
#>  true_cor     0.300
#>  null_cor     0    

# Test against a non-zero null (null_cor = 0.20), looking for evidence that rho > 0.20
ss_power_r(true_cor = 0.40, null_cor = 0.20, desired_power = 0.80, directional = TRUE)
#>  term         value
#>  necessary_N  130  
#>  actual_power 0.801
#>  true_cor     0.400
#>  null_cor     0.200
```
