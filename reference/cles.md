# Common-Language Effect Size (McGraw & Wong, 1992)

Computes the common-language (CL) effect size for two independent
groups, defined as the probability that a randomly drawn observation
from group 1 exceeds a randomly drawn observation from group 2 under
bivariate normality with equal variances: \$\$\mathrm{CL} \\=\\ \Pr(Y_1
\> Y_2) \\=\\ \Phi\\\bigl(\delta / \sqrt{2}\bigr),\$\$ where \\\delta\\
is the population standardized mean difference and \\\Phi\\ is the
standard normal cumulative distribution function. When sample sizes are
supplied, the confidence interval on CL is constructed by transforming
the noncentral *t*-based CI on Cohen's *d* (Steiger & Fouladi, 1997;
Kelley, 2007) through \\\Phi(\cdot / \sqrt{2})\\, which is
monotone-increasing so the coverage probability is preserved exactly.
This is preferred over the normal-approximation CI on CL commonly seen
in applied work (Brooks, Dalal, & Nolan, 2014).

## Usage

``` r
cles(
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
  Positive means group 1 exceeds group 2.

- n_1, n_2:

  Sample sizes in the two groups; both required when a confidence
  interval on CL is desired.

- conf_level:

  Confidence level for the CI. Default `0.95`.

- smd_lower, smd_upper:

  Optional pre-computed confidence limits on *d*. If supplied, these are
  used directly and the noncentral computation is skipped.

## Value

A `data.frame` with rows for the point estimate (`cl`) and, when sample
sizes are supplied, the lower and upper CI limits. The *d*-equivalent of
each row is also reported for transparency.

## Details

The common-language idea extends to other effect sizes; the common
language effect size for correlations is developed by Liu, Carlson, and
Kelley (2019).

**Background.** McGraw & Wong (1992) introduced the CL effect size to
make Cohen's *d* more interpretable: instead of "the means differ by 0.5
SD," one can say "in 64 treated person scores higher than the control
person." Under bivariate normality with equal variances, the population
probability \\\Pr(Y_1 \> Y_2)\\ equals \\\Phi(\delta/\sqrt{2})\\, where
\\\delta = (\mu_1 - \mu_2)/\sigma\\ (McGraw & Wong, 1992).

**Connection to other measures.** CL is identical to the AUC (Area Under
the Curve) interpretation of *d* in receiver-operating analysis. Vargha
& Delaney (2000) generalized CL to the nonparametric setting (their *A*
measure) by replacing the population *p* with its empirical Mann-Whitney
estimate; under bivariate normality the two coincide. The
success-rate-difference and number-needed-to- treat scales (Kraemer &
Kupfer, 2006; see
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md))
are linear transformations of CL: \\\mathrm{SRD} = 2 \mathrm{CL} - 1\\.

**Confidence interval construction.** Because \\\Phi(\cdot/\sqrt{2})\\
is monotone-increasing, the CI on CL is obtained by transforming the CI
on *d*: \\\[\Phi(d_L/\sqrt 2),\\ \Phi(d_U/\sqrt 2)\]\\. This is an
exact-coverage interval (under the noncentral *t* sampling model) and is
more accurate than the normal-approximation CI on CL that uses a
Wald-style variance for \\\hat p\\ (Brooks, Dalal, & Nolan, 2014).

## References

Brooks, M. E., Dalal, D. K., & Nolan, K. P. (2014). Are common language
effect sizes easier to understand than traditional effect sizes?
*Journal of Applied Psychology, 99*(2), 332–340.
[doi:10.1037/a0034745](https://doi.org/10.1037/a0034745)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kraemer, H. C., & Kupfer, D. J. (2006). Size of treatment effects and
their importance to clinical research and practice. *Biological
Psychiatry, 59*(11), 990–996.
[doi:10.1016/j.biopsych.2005.09.014](https://doi.org/10.1016/j.biopsych.2005.09.014)

Liu, X. S., Carlson, R., & Kelley, K. (2019). Common language effect
size for correlations. *The Journal of General Psychology, 146*(3),
325–338.
[doi:10.1080/00221309.2019.1585321](https://doi.org/10.1080/00221309.2019.1585321)

McGraw, K. O., & Wong, S. P. (1992). A common language effect size
statistic. *Psychological Bulletin, 111*(2), 361–365.
[doi:10.1037/0033-2909.111.2.361](https://doi.org/10.1037/0033-2909.111.2.361)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md)

Other effect size estimates:
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
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Point estimate only.
cles(smd = 0.5)
#>  term value
#>  smd  0.5  
#>  cl   0.638

# 2. With a noncentral t CI from sample sizes (preferred):
cles(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
#>  term      value
#>  smd       0.5  
#>  cl        0.638
#>  smd_lower 0.101
#>  smd_upper 0.897
#>  cl_lower  0.528
#>  cl_upper  0.737
#> 
#> Confidence level: 95%

# 3. With a pre-computed CI on d:
cles(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80)
#>  term      value
#>  smd       0.5  
#>  cl        0.638
#>  smd_lower 0.2  
#>  smd_upper 0.8  
#>  cl_lower  0.556
#>  cl_upper  0.714
#> 
#> Confidence level: 95%

# 4. CL at three reference d values:
cles(smd = 0.2)
#>  term value
#>  smd  0.2  
#>  cl   0.556
cles(smd = 0.5)
#>  term value
#>  smd  0.5  
#>  cl   0.638
cles(smd = 0.8)
#>  term value
#>  smd  0.8  
#>  cl   0.714
```
