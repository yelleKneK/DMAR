# Sample Size for AIPE on a Pearson Correlation

Determines the sample size needed for a confidence interval on a
population Pearson correlation \\\rho\\ to have a desired width
(accuracy in parameter estimation; Kelley & Maxwell, 2003). The interval
planned for is the Fisher's *Z* interval that
[`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md)
reports for `method = "pearson"`: the correlation is transformed as
\\z(r) = \mathrm{atanh}(r)\\, an interval with standard error
\\1/\sqrt{n - 3}\\ is formed on the *z* scale, and the limits are
back-transformed through \\\tanh(\cdot)\\ (Fisher, 1921; Bonett &
Wright, 2000, Equation 2). Because the plan targets the same interval
the analysis will report, the planned width and the analyzed width
agree.

## Usage

``` r
ss_aipe_r(rho, width, conf_level = 0.95, assurance = NULL)
```

## Arguments

- rho:

  Anticipated population Pearson correlation, in \\(-1, 1)\\.

- width:

  Desired full width of the confidence interval on the correlation.

- conf_level:

  Desired confidence level (default `0.95`).

- assurance:

  Optional. Probability that the realized CI is no wider than `width`
  (\\1 - \gamma\\). When supplied, the sample size is inflated using the
  standard chi squared correction (Kelley, 2008); when `NULL`, the
  assurance is fixed at 0.5.

## Value

A `data.frame` with the rows `necessary_N` (the recommended total sample
size, rounded up), `expected_width` at that sample size, and the inputs
echoed back.

## Details

**Closed-form first pass.** On the Fisher's *Z* scale the interval has
half-width \\z\_{1 - \alpha/2} / \sqrt{n - 3}\\, and the delta method
maps it back to the correlation scale as approximately \$\$w \\\approx\\
2\\ z\_{1 - \alpha/2} \\ \frac{1 - \rho^2}{\sqrt{n - 3}},\$\$ which
solves to the first-stage approximation of Bonett and Wright (2000),
\$\$n_0 \\=\\ 3 + \Big\lceil 4\\ (z\_{1 - \alpha/2})^2 (1 - \rho^2)^2 /
w^2 \Big\rceil.\$\$

**Exact iteration.** The back-transformed width depends on \\\rho\\
through \\\tanh(\cdot)\\, so the delta method approximation can land a
few observations off in either direction. Starting from \\n_0\\, the
function evaluates the exact back-transformed width \\\tanh(z\_\rho +
z\_{1 - \alpha/2}/\sqrt{n - 3}) - \tanh(z\_\rho - z\_{1 -
\alpha/2}/\sqrt{n - 3})\\ and steps the integer \\n\\ until it is the
smallest sample size whose width is at or below `width`. Where Bonett
and Wright (2000) stop after a single second-stage adjustment, this
search is exact.

**The planning value matters least near zero.** At a fixed sample size
the back-transformed width is largest at \\\rho = 0\\ and shrinks as
\\\|\rho\|\\ grows, so a planning value closer to zero yields a larger,
more conservative sample size. When little is known about the population
correlation, `rho = 0` gives the sample size that suffices for any
population value.

**When to use simple vs. partial correlation planning.** Use this
function when the inferential target is the correlation between two
variables with nothing partialed out. When the target is the correlation
after statistically controlling for other variables, see
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md).

The Monte Carlo companion
[`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md)
evaluates how the plan behaves when the population correlation differs
from the planning value.

## References

Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
estimating Pearson, Kendall and Spearman correlations. *Psychometrika,
65*(1), 23–28.
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183)

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 on correlations.)

## See also

[`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_power_r`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)

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
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Plan n so the 95% CI on the Pearson correlation has full width
# at most 0.20, when the anticipated correlation is 0.30.
ss_aipe_r(rho = 0.30, width = 0.20)
#>  term           value
#>  necessary_N    320  
#>  expected_width 0.2  
#>  rho            0.3  
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# A narrower target width requires a larger sample size.
ss_aipe_r(rho = 0.30, width = 0.10)
#>  term           value
#>  necessary_N    1274 
#>  expected_width 0.1  
#>  rho            0.3  
#>  width_target   0.1  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# With 80% assurance that the realized interval is no wider than
# the target (Kelley, 2008):
ss_aipe_r(rho = 0.30, width = 0.20, assurance = 0.80)
#>  term           value
#>  necessary_N    342  
#>  expected_width 0.193
#>  rho            0.3  
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# Planning at rho = 0 gives the sample size that suffices for any
# population correlation, since the interval is widest there.
ss_aipe_r(rho = 0, width = 0.20)
#>  term           value
#>  necessary_N    385  
#>  expected_width 0.2  
#>  rho            0    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%
```
