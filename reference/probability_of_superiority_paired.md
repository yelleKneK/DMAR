# Probability of Superiority for a Paired-Samples Design

Computes the probability-of-superiority effect size for paired
observations (Grissom & Kim, 2005, 2012), \\P_S = \Pr(Y_1 \> Y_2)\\,
along with an analytic confidence interval based on the Brunner-Munzel
(2000) U-statistic standard error and a Fisher- \\\mathrm{arctanh}\\
transformation to keep the bounds inside \\\[0, 1\]\\. The paired
counterpart of the Vargha-Delaney (2000) \\A\\ statistic /
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
for two independent groups.

## Usage

``` r
probability_of_superiority_paired(x, y, conf_level = 0.95)
```

## Arguments

- x, y:

  Paired numeric vectors of equal length. \\x\\ and \\y\\ are
  interpreted as repeated measurements on the same units (e.g.,
  pre/post, condition 1 / condition 2, sibling-1 / sibling-2).

- conf_level:

  Confidence level. Default `0.95`.

## Value

A `data.frame` with rows for the point estimate of \\P_S\\, the lower /
upper CI bounds, the variance, and the counts of within-pair wins / ties
/ losses for \\y_1\\.

## Details

**Definition.** For paired observations \\(x_i, y_i)\\, \$\$P_S \\=\\
\Pr(Y \> X) + 0.5 \cdot \Pr(Y = X),\$\$ where ties are split. The sample
estimator is the proportion of pairs with \\y_i \> x_i\\, plus half the
proportion of ties. This is the natural paired-data analog of
Vargha-Delaney's \\A\\ statistic and is unbiased under exchangeability
of paired observations.

**Why paired-specific.** The independent-groups `cles` and `cliff_delta`
estimators are biased when the two samples are paired, because their
variance formulas assume independence of the two groups. For paired data
the within-pair correlation reduces the effective sampling variance,
which is captured by the Brunner-Munzel (2000) variance used here.

**Confidence interval.** The standard error is built from the
within-pair sign indicators (Brunner-Munzel, 2000):
\$\$\mathrm{Var}(\hat P_S) \\=\\ \frac{1}{n^2}\sum\_{i=1}^{n} (s_i -
\bar s)^2,\$\$ where \\s_i = \mathrm{I}(y_i \> x_i) + 0.5 \cdot
\mathrm{I}(y_i = x_i)\\. The CI is built on the \\\mathrm{arctanh}(2
P_S - 1)\\ scale (mapping \\P_S \in \[0, 1\]\\ to the real line) and
back-transformed to keep the limits inside the unit interval, exactly
mirroring
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md).

## References

Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher
problem: Asymptotic theory and a small-sample approximation.
*Biometrical Journal, 42*(1), 17–25.
[doi:10.1002/(SICI)1521-4036(200001)42:1\<17::AID-BIMJ17\>3.0.CO;2-U](https://doi.org/10.1002/%28SICI%291521-4036%28200001%2942%3A1%3C17%3A%3AAID-BIMJ17%3E3.0.CO%3B2-U)

Grissom, R. J., & Kim, J. J. (2005). *Effect sizes for research: A broad
practical approach*. Lawrence Erlbaum.

Grissom, R. J., & Kim, J. J. (2012). *Effect sizes for research:
Univariate and multivariate applications* (2nd ed.). Routledge.

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`proportion_of_superiority`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md)

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
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Paired pre/post data:
set.seed(113)
pre  <- rnorm(30, mean = 100, sd = 15)
post <- pre + rnorm(30, mean =  5, sd = 10)
probability_of_superiority_paired(x = pre, y = post)
#>  term                       value  
#>  probability_of_superiority 0.667  
#>  lower_limit                0.484  
#>  upper_limit                0.81   
#>  var_ps                     0.00741
#>  wins_y_over_x              20     
#>  losses_y_under_x           10     
#>  ties                       0      
#> 
#> Confidence level: 95%
```
