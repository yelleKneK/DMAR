# AIPE Sample Size Planning for an Equivalence Test on the Standardized Mean Difference

Computes the minimum per-group sample size needed so that the
equivalence CI (the 100(1 - 2\\\alpha\\)% CI on the standardized mean
difference) has expected full width \\\le \omega\\, that is, expected
half-width \\\le \omega / 2\\ (Kelley, 2007; Lakens, 2017). With the
standard \\\alpha = 0.05\\ TOST level, the equivalence CI is the 90% CI.
The function inverts the large-sample variance of *d*; the noncentral
*t* distribution of *d* enters through the optional assurance step.

## Usage

``` r
ss_aipe_equivalence_smd(
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

A 5-row `data.frame` with columns `term` and `value`: the per-group
recommended sample size `necessary_n_per_group`, the implied total
`total_N`, the target `width`, the planning value `population_smd`, and
`ci_width_expected`, the expected full CI width at the chosen *n*.

## Details

**Approximate-variance plan.** The large-sample variance of *d* is
\$\$\mathrm{Var}(\hat d) \\\approx\\ (n_1 + n_2) / (n_1 n_2) + d^2 / (2
(n_1 + n_2)).\$\$ For a balanced design with per-group size \\n\\, the
half-width of the equivalence CI at level \\1 - 2\alpha\\ is
approximately \\z\_{1-\alpha} \sqrt{\mathrm{Var}(\hat d)}\\. The
function solves for the smallest integer \\n\\ giving expected
half-width \\\le \omega / 2\\.

**Assurance.** Under `assurance = q`, the function increments \\n\\
until the simulated probability that the realized half-width is \\\le
\omega / 2\\ is at least \\q\\. (Implemented as a thin Monte Carlo
overlay. At the default planning value `population_smd = 0` the shift is
typically zero; it grows with the planning value, reaching several per
group by `population_smd = 0.5` with a narrow target width.)

**Note on conservatism of the assurance plan.** The empirical simulation
study of the AIPE planner family finds that `ss_aipe_equivalence_smd()`
is tight at \\\gamma = 0.80\\ but operates on the boundary of its valid
range at \\\gamma = 0.99\\: the realized assurance at the recommended
sample size is within Monte Carlo error of the target, typically a few
tenths of a percentage point below 0.99. The mechanism is that the
planner inverts a normal approximation to \\\Pr(\widehat W \> \omega)\\,
and at the 99% level the upper tail of \\\widehat W\\ is heavier than
the approximation accounts for. Adding a small safety margin (5 to 10
subjects per group) restores the desired probability statement when
planning at high assurance;
[`ss_aipe_equivalence_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md)
reproduces the check for any one condition.

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

[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Plan for a 90% CI on d of width <= 0.20, under d_planning = 0:
ss_aipe_equivalence_smd(population_smd = 0, width = 0.20)
#>  term                  value
#>  necessary_n_per_group 542  
#>  total_N               1084 
#>  width                 0.2  
#>  population_smd        0    
#>  ci_width_expected     0.2  

# 2. Plan for the same width assuming a true d = 0.05:
ss_aipe_equivalence_smd(population_smd = 0.05, width = 0.20)
#>  term                  value
#>  necessary_n_per_group 542  
#>  total_N               1084 
#>  width                 0.2  
#>  population_smd        0.05 
#>  ci_width_expected     0.2  

# 3. With 80% assurance (the assurance path is Monte Carlo, so seed for
#    a reproducible result):
set.seed(113)
ss_aipe_equivalence_smd(population_smd = 0.05, width = 0.20, assurance = 0.80)
#>  term                  value
#>  necessary_n_per_group 542  
#>  total_N               1084 
#>  width                 0.2  
#>  population_smd        0.05 
#>  ci_width_expected     0.2  
```
