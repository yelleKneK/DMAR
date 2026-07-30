# Proportion of Superiority (Sometimes Called Cohen's \\U_3\\)

Computes the proportion of the treatment-group population that exceeds
the *control-group mean* under bivariate normality with equal variances.
This quantity is sometimes called Cohen's \\U_3\\ (Cohen, 1988). Under
those assumptions it equals \\\Phi(\delta)\\, where \\\delta\\ is the
population standardized mean difference. When sample sizes are supplied,
the CI on the proportion of superiority is constructed by transforming
the noncentral *t* CI on Cohen's *d* via \\\Phi\\, which is monotone and
therefore preserves coverage exactly.

## Usage

``` r
proportion_of_superiority(
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

  Sample standardized mean difference (Cohen's *d*). Numeric scalar.

- n_1, n_2:

  Group sample sizes (required if a CI is wanted).

- conf_level:

  Confidence level for the CI. Default `0.95`.

- smd_lower, smd_upper:

  Optional pre-computed CI limits on *d*; when supplied directly, the
  function skips the noncentral *t* step and just transforms these
  limits.

## Value

A `data.frame` with rows for *d*, the proportion of superiority, and
(when a CI is constructable) the lower / upper limits on *d* and on the
proportion of superiority.

## Details

The proportion of superiority is one of three "U" indices Cohen (1988)
defined; the other two (\\U_1\\, the proportion of non-overlap, and
\\U_2\\, the proportion of either population that exceeds the same
percentile in the other) can be derived from it directly: with \\U_3 =
\Phi(\delta)\\ for the proportion of superiority, Cohen's \\U_2 =
\Phi(\delta/2)\\ and \\U_1 = (2 \cdot U_2 - 1) / U_2\\ (Cohen, 1988,
Table 2.2.1).

**Why this rather than `cles`.** The proportion of superiority answers
the question "what fraction of the treatment population exceeds the
*control-group mean*,” whereas
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) answers
"what fraction of randomly drawn pairs favor the treatment over the
control.” Both are unitless probability-scale summaries of a Cohen's-*d*
difference, but the proportion of superiority is marginal while CLES is
paired. Specifically, \\\Phi(d)\\ versus \\\Phi(d/\sqrt{2})\\; for \\d =
0.5\\, the proportion of superiority is 0.69 and CLES is 0.64.

**CI construction.** Because \\\Phi(\cdot)\\ is monotone, the CI on the
proportion of superiority is just \\\[\Phi(d_L),\\ \Phi(d_U)\]\\ where
\\\[d_L,\\ d_U\]\\ is the noncentral *t* CI on *d* from
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md).

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum. (See Section 2.2
for the \\U_1\\, \\U_2\\, and \\U_3\\ indices.)

Hedges, L. V., & Olkin, I. (1985). *Statistical methods for
meta-analysis*. Academic Press.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. (The noncentral *t* interval on the standardized
mean difference that is transformed here.)
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md)

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
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Proportion of superiority at three reference d values:
proportion_of_superiority(smd = 0.2)
#>  term                      value
#>  smd                       0.2  
#>  proportion_of_superiority 0.579
proportion_of_superiority(smd = 0.5)
#>  term                      value
#>  smd                       0.5  
#>  proportion_of_superiority 0.691
proportion_of_superiority(smd = 0.8)
#>  term                      value
#>  smd                       0.8  
#>  proportion_of_superiority 0.788

# 2. With a noncentral t CI from sample sizes:
proportion_of_superiority(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
#>  term                      value
#>  smd                       0.5  
#>  proportion_of_superiority 0.691
#>  smd_lower                 0.101
#>  smd_upper                 0.897
#>  lower_limit               0.54 
#>  upper_limit               0.815
#> 
#> Confidence level: 95%
```
