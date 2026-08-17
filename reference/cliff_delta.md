# Cliff's \\\delta\\ Ordinal Effect Size

Computes Cliff's (1993) \\\delta\\ statistic for two independent groups,
the difference between the probability that a randomly drawn observation
from group 1 exceeds one from group 2 and the reverse probability,
together with an analytic confidence interval built from the U-statistic
variance (Cliff, 1996). Most R implementations of Cliff's \\\delta\\
fall back to a bootstrap CI; the analytic CI here is faster,
deterministic, and exact in the large-sample limit.

## Usage

``` r
cliff_delta(group_1, group_2, conf_level = 0.95)
```

## Arguments

- group_1, group_2:

  Numeric vectors of observations in the two groups. Ordinal data are
  fine; the statistic uses only ranks.

- conf_level:

  Confidence level for the CI. Default `0.95`.

## Value

A `data.frame` with rows for the point estimate `cliff_delta` and the
lower/upper CI bounds. The output also reports the proportion of pairs
with \\y_1 \> y_2\\, the proportion with \\y_1 \< y_2\\, and the
proportion of ties.

## Details

**Definition.** Cliff's \\\delta\\ is \$\$\delta \\=\\ \Pr(Y_1 \> Y_2) -
\Pr(Y_1 \< Y_2) \\=\\ 2 \cdot A - 1,\$\$ where \\A\\ is the
Vargha-Delaney (2000) statistic. The sample estimator is \\\hat\delta =
(\\\\(i,j): y\_{1i} \> y\_{2j}\\ - \\\\(i,j): y\_{1i} \< y\_{2j}\\) /
(n_1 n_2)\\. Ties contribute zero to both counts. \\\delta\\ ranges over
\\\[-1, 1\]\\, with 0 indicating no stochastic dominance.

**Analytic CI.** The asymptotic variance of \\\hat\delta\\ is (Cliff,
1993; restated as Feng & Cliff, 2004, Equation 2, p. 323)
\$\$\mathrm{Var}(\hat\delta) \\=\\ \frac{(n_2 - 1) \sigma^2\_{d_1} +
(n_1 - 1) \sigma^2\_{d_2} + \sigma^2_d}{n_1 n_2},\$\$ where
\\\sigma^2\_{d_i}\\ is the variance of the per-observation dominance
scores within each group. (Feng & Cliff's printed equation transposes
the \\(n_1 - 1)\\ and \\(n_2 - 1)\\ coefficients, which matters only for
unequal group sizes; the pairing above is the correct one, checked by
simulation against the empirical variance of \\\hat\delta\\.) The CI is
constructed on the Fisher-style \\\mathrm{arctanh}\\-transformed scale
and back-transformed to respect the bounded range of \\\delta\\
(analogous to Fisher's *Z* CI for Pearson \\r\\). Feng & Cliff (2004,
Equation 5, p. 324) recommend an alternative asymmetric interval that
models the dependence of the variance on \\\delta\\; the two
constructions agree to first order.

**Connection to other measures.** Cliff's \\\delta\\ is a linear
transformation of the Vargha-Delaney (2000) \\A\\ statistic (\\\delta =
2A - 1\\) and of the Mann-Whitney \\U\\ statistic (\\U / (n_1 n_2) =
A\\). It is the ordinal analog of the common-language effect size
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) and is
preferable when bivariate normality is implausible (skewed outcomes,
ordinal scales).

## References

Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
ordinal questions. *Psychological Bulletin, 114*(3), 494–509.
[doi:10.1037/0033-2909.114.3.494](https://doi.org/10.1037/0033-2909.114.3.494)

Cliff, N. (1996). *Ordinal methods for behavioral data analysis*.
Lawrence Erlbaum.

Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
improved confidence interval. *Journal of Modern Applied Statistical
Methods, 3*(2), 322–332.
[doi:10.22237/jmasm/1099267560](https://doi.org/10.22237/jmasm/1099267560)

Long, J. D., Feng, D., & Cliff, N. (2003). Ordinal analysis of
behavioral data. In I. B. Weiner (Ed.), *Handbook of psychology, Vol. 2:
Research methods* (pp.\\ 635–661). Wiley.

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`ss_aipe_cliff_delta`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
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
# 1. Two groups of different sizes, no ties:
set.seed(113)
a <- rnorm(30, mean = 0, sd = 1)
b <- rnorm(40, mean = 0.5, sd = 1)
cliff_delta(a, b)
#>  term            value 
#>  cliff_delta     -0.223
#>  lower_limit     -0.47 
#>  upper_limit     0.0557
#>  var_cliff_delta 0.0188
#>  p_y1_greater    0.388 
#>  p_y1_less       0.612 
#>  p_tied          0     
#> 
#> Confidence level: 95%

# 2. With ties (ordinal data):
o1 <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)
o2 <- c(2, 3, 3, 4, 4, 5, 5, 5)
cliff_delta(o1, o2)
#>  term            value 
#>  cliff_delta     -0.403
#>  lower_limit     -0.776
#>  upper_limit     0.179 
#>  var_cliff_delta 0.0675
#>  p_y1_greater    0.194 
#>  p_y1_less       0.597 
#>  p_tied          0.208 
#> 
#> Confidence level: 95%

# 3. Robust to right skew. Cliff's delta on the raw, untransformed
#    drinking outcome from the Smith, Meyers, and Delaney (1998)
#    trial, comparing the Community Reinforcement Approach (CRA)
#    against standard care. Because the statistic uses only ranks it
#    needs no normalizing transformation of the heavily skewed
#    outcome, unlike the standardized mean difference.
data(drinks_trial)
cra <- drinks_trial$drinks_per_week[drinks_trial$treatment == "CRA"]
std <- drinks_trial$drinks_per_week[drinks_trial$treatment == "Standard"]
cliff_delta(cra, std)
#>  term            value  
#>  cliff_delta     -0.3   
#>  lower_limit     -0.535 
#>  upper_limit     -0.0213
#>  var_cliff_delta 0.0179 
#>  p_y1_greater    0.293  
#>  p_y1_less       0.593  
#>  p_tied          0.114  
#> 
#> Confidence level: 95%
```
