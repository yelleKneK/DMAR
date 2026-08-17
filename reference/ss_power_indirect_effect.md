# Sample Size Planning for Power for the Indirect (Mediation) Effect

Power, or the sample size required for a desired power, for the test of
the indirect effect \\a b\\ in the simple mediation model, with paths
specified in standardized metric (unit-variance \\X\\, \\M\\, and
\\Y\\). The default test is joint significance (the indirect effect is
declared when both \\\hat a\\ and \\\hat b\\ are individually
significant), which tracks the resampling tests' power closely and far
exceeds the Sobel test in small samples (Fritz & MacKinnon, 2007); the
Sobel normal-theory test is available for comparison. This is the power
counterpart of the accuracy in parameter estimation (AIPE) planner
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
and the planning complement of the analysis function
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md).

## Usage

``` r
ss_power_indirect_effect(
  a,
  b,
  c_prime = 0,
  desired_power = NULL,
  N = NULL,
  alpha_level = 0.05,
  method = c("joint_significance", "sobel")
)
```

## Arguments

- a:

  Standardized \\X \to M\\ path.

- b:

  Standardized \\M \to Y\\ path, holding \\X\\.

- c_prime:

  Standardized direct effect of \\X\\ on \\Y\\ holding \\M\\. Defaults
  to 0; it enters only through the residual variance of \\Y\\.

- desired_power:

  Desired power; supply this to solve for \\N\\.

- N:

  Total sample size; supply this to evaluate the realized power. Specify
  exactly one of `desired_power` and `N`. This is the total sample size.

- alpha_level:

  Two-sided Type I error rate for each component test. Defaults to 0.05.

- method:

  `"joint_significance"` (default) or `"sobel"`.

## Value

A tidy `data.frame` with `necessary_N` (or `specified_N`),
`actual_power` (the power to detect the indirect effect, the quantity
the sample size is planned against), the component powers (`power_a`,
`power_b`; `NA` for the Sobel method), the paths (`a`, `b`, `c_prime`),
the implied `indirect_effect`, and `alpha_level`. The method is recorded
in the `"method"` attribute. The result carries the `dmar_ss_power`
class, so [`tidy`](https://generics.r-lib.org/reference/tidy.html)
reports the sample size and the power to detect the indirect effect, and
[`glance`](https://generics.r-lib.org/reference/glance.html) adds the
component powers and the planning inputs.

## Details

With unit-variance variables, the large-sample standard errors are
\\\mathrm{se}\_a = \sqrt{(1 - a^2)/N}\\ and \\\mathrm{se}\_b =
\sqrt{\sigma^2\_{e_Y} / \[N (1 - a^2)\]}\\ with \\\sigma^2\_{e_Y} = 1 -
(b^2 + c'^2 + 2 a b c')\\. Because \\\hat a\\ and \\\hat b\\ are
asymptotically independent in this model, the joint significance power
is the product of the two component powers; the Sobel power refers \\ab
/ \mathrm{se}\_{ab}\\ (first-order delta method, via the same variance
as
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md))
to the normal. The joint-significance component powers use the exact
noncentral *t* (with \\n - 2\\ and \\n - 3\\ degrees of freedom), so its
only approximations are the population standard errors and component
independence; the tests validate the result against raw-data simulation.
A specified parameter combination must be admissible (positive residual
variances), or the function stops.

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
# Fritz and MacKinnon's (2007) running scenario (a = b = .39). The
# joint significance approximation returns necessary_N = 65; raw-data
# simulation puts the power at N = 65 nearer .77 and reaches .80 near
# N = 70, which is why Fritz and MacKinnon's simulation-based table
# reports a somewhat larger requirement.
ss_power_indirect_effect(a = .39, b = .39, desired_power = .80)
#>  term            value
#>  necessary_N     65   
#>  actual_power    0.802
#>  power_a         0.92 
#>  power_b         0.872
#>  a               0.39 
#>  b               0.39 
#>  c_prime         0    
#>  indirect_effect 0.152
#>  alpha_level     0.05 

# A near-zero a path against a larger b: the weak link drives the requirement.
ss_power_indirect_effect(a = .14, b = .39, desired_power = .80)
#>  term            value 
#>  necessary_N     395   
#>  actual_power    0.8   
#>  power_a         0.8   
#>  power_b         1     
#>  a               0.14  
#>  b               0.39  
#>  c_prime         0     
#>  indirect_effect 0.0546
#>  alpha_level     0.05  

# Realized power at a given N, and the Sobel comparison (always lower).
ss_power_indirect_effect(a = .39, b = .39, N = 75)
#>  term            value
#>  specified_N     75   
#>  actual_power    0.871
#>  power_a         0.951
#>  power_b         0.915
#>  a               0.39 
#>  b               0.39 
#>  c_prime         0    
#>  indirect_effect 0.152
#>  alpha_level     0.05 
ss_power_indirect_effect(a = .39, b = .39, N = 75, method = "sobel")
#>  term            value
#>  specified_N     75   
#>  actual_power    0.7  
#>  power_a         <NA> 
#>  power_b         <NA> 
#>  a               0.39 
#>  b               0.39 
#>  c_prime         0    
#>  indirect_effect 0.152
#>  alpha_level     0.05 
```
