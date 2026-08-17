# Sample Size for AIPE on an Intraclass Correlation Coefficient

Determines the sample size needed for a confidence interval on a
population intraclass correlation coefficient (ICC) to have a desired
width, using Bonett's (2002) Fisher-style variance-stabilizing
transformation. The function inverts the asymptotic variance on the
transformed scale (where the CI is symmetric and approximately normal),
solves for the smallest \\n\\ that achieves the target half-width on the
back-transformed (raw-ICC) scale, and optionally inflates the result by
a chi squared assurance correction (Kelley & Maxwell, 2003).

## Usage

``` r
ss_aipe_icc(
  rho,
  k,
  width,
  which_width = c("Full", "Lower", "Upper"),
  conf_level = 0.95,
  type = c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)", "ICC(1,k)", "ICC(2,k)", "ICC(3,k)"),
  assurance = NULL
)
```

## Arguments

- rho:

  Anticipated population ICC at the level matching `type`, in \\\[0,
  1)\\. This is the value the researcher expects the truth to be near,
  motivated by prior literature, a pilot, or substantive theory. Because
  the planning formula inverts the asymptotic variance *at* this value,
  a wrong guess inflates or deflates the realized confidence interval
  width relative to `width`; the function
  [`ss_aipe_icc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md)
  quantifies the impact of misspecification by Monte Carlo.

- k:

  Number of raters (or measurements per subject); must be at least 2.

- width:

  Desired full width of the back-transformed CI on the ICC.

- which_width:

  Whether `width` is the `"Full"` width of the interval (default) or a
  half-width: `"Lower"` and `"Upper"` both interpret `width` as half the
  full width, so they plan for a full width of twice `width` and return
  the same sample size. Because the interval is not generally symmetric
  about the estimate, its realized lower and upper half-widths can
  differ from each other and from half the full width; the planner does
  not target them separately. A genuinely one-sided width target is not
  currently offered.

- conf_level:

  Desired confidence level (default `0.95`).

- type:

  Which Shrout-Fleiss (1979) ICC form is being planned. One of the
  single-rater forms `"ICC(1,1)"`, `"ICC(2,1)"`, `"ICC(3,1)"` (which
  share the planning variance) or the average-of-\\k\\ forms
  `"ICC(1,k)"`, `"ICC(2,k)"`, `"ICC(3,k)"`; default `"ICC(1,1)"`. Both
  `rho` and `width` are interpreted on the scale of the requested form:
  for an average-of-\\k\\ form the planning value is mapped to the
  single-rater scale through the inverse Spearman-Brown relation and
  each candidate confidence limit is mapped back, so the returned sample
  size targets the width of the interval on the average-of-\\k\\ ICC
  itself (see Details).

- assurance:

  Optional. Probability that the realized CI is no wider than `width`;
  when supplied, the sample size is inflated by the standard chi squared
  correction.

## Value

A `data.frame` with rows for the recommended sample size (*number of
subjects*), the expected back-transformed CI width, and the inputs
echoed back. The Shrout-Fleiss form the plan targets is stored as the
`"icc_type"` attribute so the `value` column stays numeric.

## Details

**Bonett's (2002) Fisher-style transform.** Bonett (2002) showed that
the transformation \$\$L(\rho) \\=\\ \frac{1}{2} \log\\\left( \frac{1 +
(k - 1)\rho}{1 - \rho}\right)\$\$ approximately variance-stabilizes the
single-rater ICC, with \$\$\mathrm{Var}(L(\hat\rho)) \\\approx\\
\frac{k}{2\\(k - 1)\\(n - 2)}.\$\$ A confidence interval is constructed
by adding \\\pm z\_{1-\alpha/2}\\ standard errors on the \\L\\ scale and
back-transforming to the raw-ICC scale via \\\rho = (e^{2L} - 1) /
(e^{2L} - 1 + k)\\. The minimum sample size is found by searching for
the smallest \\n\\ whose back-transformed CI width is below the target.

**Single-rater vs.\\ average-of-\\k\\ ICC.** The Bonett (2002) variance
applies directly to the single-rater forms (`ICC(1,1)`, `ICC(2,1)`,
`ICC(3,1)`). For the average-of-\\k\\ forms the planning value
\\\rho_k\\ is first mapped to the single-rater scale through the inverse
Spearman-Brown relation \\\rho = \rho_k / \[k - (k - 1)\rho_k\]\\, the
interval is formed on the \\L\\ scale as above, and each candidate limit
is mapped back to the average-of-\\k\\ scale (composing the inverse
\\L\\ transform with the Spearman-Brown formula reduces to \\\rho_k =
1 - e^{-2L}\\). The width criterion therefore applies to the confidence
interval on the average-of-\\k\\ ICC itself, following the convention
used by
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md).
Because the two scales differ, an average-of-\\k\\ plan generally
recommends a different sample size than a single-rater plan at the same
numeric `rho` and `width`; at `rho = 0.7`, `k = 3`, and `width = 0.20`,
planning for `ICC(1,1)` recommends \\n = 69\\ subjects while planning
for `ICC(1,k)` recommends \\n = 110\\.

**The assurance correction can fall slightly short.** Monte Carlo
evaluation with
[`ss_aipe_icc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md)
shows that the chi squared inflation tends to deliver a little less
assurance than requested. At the condition of the second example below
(`rho = 0.7`, `k = 3`, `width = 0.20`, `assurance = 0.80`), the
recommended \\n = 79\\ yields an empirical assurance of about .77
against the requested .80 (10,000 replications of the *F*-based interval
computed by [`icc`](https://yelleknek.github.io/DMAR/reference/icc.md)),
and the smallest sample size whose empirical assurance reaches .80 is
\\n = 81\\. The mechanism is that the realized interval widths on the
raw-ICC scale have a heavier upper tail than the chi squared inflation
on the transformed scale accounts for, so the buffer the correction adds
is slightly too small. When meeting the assurance target matters, check
the recommended sample size with
[`ss_aipe_icc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md)
and increase \\n\\ until the empirical assurance reaches the target.

## References

Bonett, D. G. (2002). Sample size requirements for estimating intraclass
correlations with desired precision. *Statistics in Medicine, 21*(9),
1331–1335. [doi:10.1002/sim.1108](https://doi.org/10.1002/sim.1108)

Donner, A. (1986). A review of inference procedures for the intraclass
correlation coefficient in the one-way random effects model.
*International Statistical Review, 54*(1), 67–82.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428.

Smith, C. A. B. (1956). On the estimation of intraclass correlation.
*Annals of Human Genetics, 21*(4), 363–373.

## See also

[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md)

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
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
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
# 1. Plan n so the 95% CI on a single-rater ICC has full width <= 0.20
#        with k = 3 raters and an anticipated ICC of 0.7.
ss_aipe_icc(rho = 0.7, k = 3, width = 0.20)
#>  term           value
#>  necessary_N    69   
#>  expected_width 0.199
#>  rho            0.7  
#>  k              3    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# 2. With 80% assurance:
ss_aipe_icc(rho = 0.7, k = 3, width = 0.20, assurance = 0.80)
#>  term           value
#>  necessary_N    79   
#>  expected_width 0.186
#>  rho            0.7  
#>  k              3    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%
```
