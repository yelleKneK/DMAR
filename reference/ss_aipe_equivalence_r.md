# AIPE Sample Size Planning for an Equivalence Test on the Pearson Correlation

Computes the minimum sample size needed so that the equivalence CI (the
100(1 - 2\\\alpha\\)% CI on the Pearson correlation \\\rho\\) has
expected full width \\\le \omega\\ (Kelley, 2007; Lakens, 2017). With
the standard \\\alpha = 0.05\\ TOST level, the equivalence CI is the 90%
CI. The interval is the Fisher's \\Z\\ construction that
[`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md)
and
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
use, so the plan and the analysis invert the same interval.

## Usage

``` r
ss_aipe_equivalence_r(
  population_r = 0,
  width,
  alpha_level = 0.05,
  assurance = NULL
)
```

## Arguments

- population_r:

  Anticipated population correlation \\\rho\\ used to plan the width.
  Default `0`: the confidence interval for a correlation is widest at
  \\\rho = 0\\, so the default plans under the widest-interval case and
  is conservative for any other value.

- width:

  Target full CI width on the correlation scale (e.g., `0.20` for a 90%
  CI of width 0.20). Must be in \\(0, 2)\\, the width of the correlation
  scale itself.

- alpha_level:

  One-sided TOST significance level. The CI used in planning is at
  confidence level \\1 - 2\alpha\\. Default `0.05` (90% CI).

- assurance:

  Optional assurance probability in \\(0, 1)\\. When supplied, the
  function chooses *N* so that the probability of achieving `width` or
  less is at least `assurance` (Kelley, Maxwell, & Rausch, 2003).
  Default `NULL` (no assurance correction).

## Value

A 4-row `data.frame` with columns `term` and `value`: the recommended
sample size `necessary_N`, the target `width`, the planning value
`population_r`, and the resulting `ci_width_expected` at the chosen *N*.

## Details

**Closed form on the Fisher's \\Z\\ scale.** The equivalence CI has
half-width \\h = z\_{1-\alpha} / \sqrt{N - 3}\\ on the Fisher's \\Z\\
scale, and its width on the correlation scale is \$\$w(N) \\=\\
\tanh(Z\_\rho + h) - \tanh(Z\_\rho - h),\$\$ where \\Z\_\rho =
\tanh^{-1}(\rho)\\. The function returns the smallest integer \\N \ge
4\\ with \\w(N) \le \omega\\. At \\\rho = 0\\ this is available in
closed form, \\N = \lceil 3 + (z\_{1-\alpha} / \tanh^{-1}(\omega / 2))^2
\rceil, \\ and away from zero the back-transform shortens the interval,
so the required *N* can only decrease as \\\|\rho\|\\ grows.

**Choosing the width from equivalence bounds.** To leave room for an
equivalence verdict inside bounds \\(-b, b)\\, the interval must at
minimum fit inside the bounds when centered at the anticipated \\\rho\\,
so a width somewhat below \\2 b\\ (for \\\rho\\ near 0) is the natural
target; the Monte Carlo sensitivity sibling
[`ss_aipe_equivalence_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md)
reports the realized proportion of equivalence verdicts at the planned
*N*.

**Assurance.** Under `assurance = q`, the function increments \\N\\
until the Monte Carlo probability that the realized width is \\\le
\omega\\ is at least \\q\\, drawing the sampling distribution of
\\\widehat Z\\ as normal with mean \\Z\_\rho\\ and variance \\1 / (N -
3)\\.

## References

Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for comparing
correlation and regression coefficients. *British Journal of
Mathematical and Statistical Psychology, 68*(2), 292–309.
[doi:10.1111/bmsp.12045](https://doi.org/10.1111/bmsp.12045)

Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of
association: An equivalence testing approach. *British Journal of
Mathematical and Statistical Psychology, 63*(3), 527–537.
[doi:10.1348/000711009X475853](https://doi.org/10.1348/000711009X475853)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample size planning.
*Evaluation and the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Lakens, D. (2017). Equivalence tests: A practical primer for *t* tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*(4), 355–362.
[doi:10.1177/1948550617697177](https://doi.org/10.1177/1948550617697177)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md)

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
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
# 1. Plan for a 90% CI on the correlation of width <= 0.20, under
#    the widest-interval planning value rho = 0:
ss_aipe_equivalence_r(population_r = 0, width = 0.20)
#>  term              value
#>  necessary_N       272  
#>  width             0.2  
#>  population_r      0    
#>  ci_width_expected 0.2  

# 2. The same width assuming a true correlation of 0.30 requires
#    fewer participants, since the interval narrows away from zero:
ss_aipe_equivalence_r(population_r = 0.30, width = 0.20)
#>  term              value
#>  necessary_N       226  
#>  width             0.2  
#>  population_r      0.3  
#>  ci_width_expected 0.2  

# 3. With 80% assurance (the assurance path is Monte Carlo, so seed
#    for a reproducible result):
set.seed(113)
ss_aipe_equivalence_r(population_r = 0.30, width = 0.20, assurance = 0.80)
#>  term              value
#>  necessary_N       240  
#>  width             0.2  
#>  population_r      0.3  
#>  ci_width_expected 0.194
```
