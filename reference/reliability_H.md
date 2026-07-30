# Maximal Reliability Coefficient *H* (Hancock & Mueller, 2001)

Computes the Hancock-Mueller (2001) maximal-reliability coefficient *H*
from a vector of standardized factor loadings; *H* is the reliability of
the optimally-weighted composite of a set of indicators of a single
latent construct. It is uniformly greater than or equal to coefficient
alpha and McDonald's omega for the same data, so it sets a useful upper
bound on what reliability can plausibly be for that indicator set. A
delta method confidence interval is reported when standard errors of the
standardized loadings are supplied.

## Usage

``` r
reliability_H(loadings, se_loadings = NULL, conf_level = 0.95)
```

## Arguments

- loadings:

  Numeric vector of standardized factor loadings, each in \\(-1, 1)\\.
  At least 2 loadings are required.

- se_loadings:

  Optional vector of standard errors of the standardized loadings (same
  length as `loadings`). When supplied, a delta method CI on *H* is
  reported.

- conf_level:

  Confidence level for the CI. Default `0.95`.

## Value

A `data.frame` with rows for the point estimate `reliability_H` and
(when SEs are supplied) the lower / upper CI bounds and the delta method
variance.

## Details

**Definition.** For \\p\\ indicators of a single latent factor with
standardized loadings \\\lambda_1, \ldots, \lambda_p\\, Hancock &
Mueller (2001) showed that the maximum reliability achievable by any
linear composite of the indicators is \$\$H \\=\\ \frac{\sum\_{i=1}^{p}
\lambda_i^2 / (1 - \lambda_i^2)} {1 + \sum\_{i=1}^{p} \lambda_i^2 / (1 -
\lambda_i^2)}.\$\$ Equivalently, defining \\\theta_i = \lambda_i^2 /
(1 - \lambda_i^2)\\ (the signal-to-noise ratio for indicator \\i\\), \\H
= \sum \theta_i / (1 + \sum \theta_i)\\. As \\p\\ grows or as the
individual loadings grow toward 1, \\H \to 1\\.

**Relationship to coefficient alpha and omega.** Coefficient alpha
([`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md))
is the reliability of the *equally-weighted* sum of indicators; *H* is
the reliability of the *optimally-weighted* composite. Hancock & Mueller
(2001) prove \\H \ge \omega \ge \alpha\\ for a unidimensional indicator
set, with equality only when all loadings are equal. *H* is therefore
most useful for diagnostics: if *H* is much higher than alpha, the
standard composite is leaving reliability on the table.

**Confidence interval via the delta method.** Conditional on standard
errors \\\mathrm{SE}(\hat \lambda_i)\\, the delta method variance of *H*
is \$\$\mathrm{Var}(\hat H) \\\approx\\ \sum\_{i=1}^{p}
\left(\frac{\partial H}{\partial \lambda_i}\right)^2 \mathrm{SE}(\hat
\lambda_i)^2,\$\$ with \\\partial H / \partial \lambda_i = 2 \lambda_i /
\[(1 - \lambda_i^2)^2 (1 + \sum_j \theta_j)^2\]\\. The CI is built on
the \\\mathrm{logit}(H)\\ scale (mapping \\\[0, 1\]\\ to the real line)
and back-transformed, as recommended by Browne (1968) for bounded
reliability coefficients.

## References

Browne, M. W. (1968). A comparison of factor analytic techniques.
*Psychometrika, 33*(3), 267–334.

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Raykov, T. (1997). Estimation of composite reliability for congeneric
measures. *Applied Psychological Measurement, 21*(2), 173–184.
[doi:10.1177/01466216970212006](https://doi.org/10.1177/01466216970212006)

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

## See also

[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Five indicators with standardized loadings 0.6, 0.7, ..., 0.8:
reliability_H(loadings = c(0.6, 0.65, 0.70, 0.75, 0.80))
#>  term          value
#>  reliability_H 0.842
#> 
#> Confidence level: 95%

# 2. With per-loading standard errors from a CFA output:
reliability_H(loadings    = c(0.6, 0.7, 0.8),
               se_loadings = c(0.05, 0.04, 0.03))
#>  term          value   
#>  reliability_H 0.767   
#>  lower_limit   0.716   
#>  upper_limit   0.812   
#>  var_H         0.000599
#> 
#> Confidence level: 95%
```
