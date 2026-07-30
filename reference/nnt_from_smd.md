# Number Needed to Treat (NNT) From Cohen's *d*

Converts a standardized mean difference (Cohen's *d*) into the number
needed to treat (NNT) using the Kraemer-Kupfer (2006) framework, which
connects *d* to the success-rate difference (SRD) under the assumption
of a continuous, normally distributed outcome with equal variances and a
hypothetical median-cut criterion for treatment success. When a
confidence interval on *d* is supplied (or noncentrality parameters /
sample sizes are provided so that one can be constructed), the bounds
are propagated through the same conversion to give a CI on the NNT.

## Usage

``` r
nnt_from_smd(
  smd,
  n_1 = NULL,
  n_2 = NULL,
  conf_level = 0.95,
  smd_lower = NULL,
  smd_upper = NULL
)
```

## Arguments

- smd:

  Sample standardized mean difference (Cohen's *d*); a numeric scalar.
  Positive values correspond to the treatment group exceeding the
  control group.

- n_1, n_2:

  Sample sizes in the treatment and control groups; both required when a
  noncentral *t*-based CI on the NNT is desired.

- conf_level:

  Confidence level for the CI on the NNT (when `n_1` and `n_2` are
  supplied). Default `0.95`.

- smd_lower, smd_upper:

  Optional pre-computed confidence limits on *d*. If supplied, these are
  used directly to propagate the interval through the SRD-to-NNT map and
  the noncentral *t* computation is skipped.

## Value

A `data.frame` with rows for the success-rate difference (`srd`), the
point estimate of NNT (`nnt`), and (when an interval is constructable)
the lower and upper NNT limits. When the lower CI on *d* is exactly
zero, the corresponding NNT bound is reported as `Inf`; when it is
negative, that bound is a finite negative value (the NNT to harm), so a
CI on *d* that spans zero yields an NNT interval passing through the
infinite point that separates benefit from harm.

## Details

**The conversion.** Under bivariate normality with equal variances,
Kraemer & Kupfer (2006) showed that the proportion of times a randomly
drawn treatment-group observation exceeds a randomly drawn control-group
observation is \$\$p \\=\\ \Pr(Y_T \> Y_C) \\=\\ \Phi\\\bigl(d /
\sqrt{2}\bigr),\$\$ from which the success-rate difference (their effect
size) is \$\$\mathrm{SRD} \\=\\ 2 p - 1 \\=\\ 2 \Phi\\\bigl(d /
\sqrt{2}\bigr) - 1,\$\$ and the number needed to treat is its
reciprocal, \$\$\mathrm{NNT} \\=\\ 1 / \mathrm{SRD}.\$\$ Larger *d*
produces smaller NNT; \\d = 0\\ produces \\\mathrm{NNT} = \infty\\ (no
advantage). The conversion is monotone, so the SRD/NNT confidence
interval is obtained by applying the same transformation to the
endpoints of the CI on *d*; the lower NNT limit comes from the *upper*
*d* limit and vice versa (Furukawa & Leucht, 2011).

**When NNT becomes infinite or negative.** If the lower CI on *d* is
exactly zero, the corresponding upper NNT bound is `Inf`: the data do
not exclude the possibility that the treatment produces no advantage (or
even harm). Negative values of *d* are allowed; the function returns
negative NNT values which are conventionally read as the NNT to *harm*.

**Assumption check.** The Kraemer-Kupfer conversion assumes a
continuous, normally distributed outcome with equal variances across
groups. For skewed outcomes, ordinal outcomes, or unequal variances, the
empirical common-language effect size
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) or the
Vargha-Delaney *A* statistic is more defensible. Furukawa & Leucht
(2011) compare four methods and recommend the Kraemer-Kupfer formula as
the most accurate under normality.

## References

Furukawa, T. A., & Leucht, S. (2011). How to obtain NNT from Cohen's
*d*: Comparison of four methods. *PLoS ONE, 6*(4), e19070.
[doi:10.1371/journal.pone.0019070](https://doi.org/10.1371/journal.pone.0019070)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. (The noncentral *t* interval on the standardized
mean difference that is mapped to the NNT here.)
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kraemer, H. C., & Kupfer, D. J. (2006). Size of treatment effects and
their importance to clinical research and practice. *Biological
Psychiatry, 59*(11), 990–996.
[doi:10.1016/j.biopsych.2005.09.014](https://doi.org/10.1016/j.biopsych.2005.09.014)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Point estimate only.
nnt_from_smd(smd = 0.5)
#>  term value
#>  smd  0.5  
#>  srd  0.276
#>  nnt  3.62 

# 2. With a noncentral t CI from sample sizes:
nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
#>  term      value
#>  smd       0.5  
#>  srd       0.276
#>  nnt       3.62 
#>  smd_lower 0.101
#>  smd_upper 0.897
#>  nnt_lower 2.11 
#>  nnt_upper 17.6 
#> 
#> Confidence level: 95%

# 3. With pre-computed CI on d:
nnt_from_smd(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80)
#>  term      value
#>  smd       0.5  
#>  srd       0.276
#>  nnt       3.62 
#>  smd_lower 0.2  
#>  smd_upper 0.8  
#>  nnt_lower 2.33 
#>  nnt_upper 8.89 
#> 
#> Confidence level: 95%

# 4. Lower d below zero: upper NNT bound is a finite negative value (NNT to harm).
nnt_from_smd(smd = 0.3, smd_lower = -0.10, smd_upper = 0.70)
#>  term      value
#>  smd       0.3  
#>  srd       0.168
#>  nnt       5.95 
#>  smd_lower -0.1 
#>  smd_upper 0.7  
#>  nnt_lower 2.64 
#>  nnt_upper -17.7
#> 
#> Confidence level: 95%
```
