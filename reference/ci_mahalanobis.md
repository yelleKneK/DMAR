# Confidence Interval for the Squared Mahalanobis Distance

Computes the squared Mahalanobis distance \\D^2\\ together with an exact
confidence interval for the population squared distance \\\Delta^2\\,
obtained by inverting Hotelling's \\T^2\\ statistic through its
(noncentral) *F*-distribution as in Reiser (2001). Both the one-sample
setting (mean vector against a hypothesized population mean) and the
two-sample setting (between-groups distance from discriminant analysis
or multivariate group comparison) are supported, and either raw data or
a pre-computed \\D^2\\ with sample sizes can be supplied.

## Usage

``` r
ci_mahalanobis(
  D2 = NULL,
  group_1 = NULL,
  group_2 = NULL,
  mu_0 = NULL,
  n_1 = NULL,
  n_2 = NULL,
  p = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- D2:

  Optional pre-computed squared Mahalanobis distance. Ignored if
  `group_1` is supplied.

- group_1:

  Optional numeric matrix or data frame for the first sample (\\n_1
  \times p\\). Rows are observations and columns are variables.

- group_2:

  Optional numeric matrix or data frame for the second sample (\\n_2
  \times p\\). When `NULL`, the function operates in one-sample mode
  against `mu_0`.

- mu_0:

  Optional hypothesized population mean for the one-sample case
  (length-\\p\\ numeric vector). Defaults to a vector of zeros.

- n_1:

  Sample size for group 1 (required when supplying `D2` directly).

- n_2:

  Sample size for group 2 (required for two-sample mode when supplying
  `D2` directly; leave `NULL` for one-sample mode).

- p:

  Dimensionality (number of variables) when supplying `D2` directly.

- conf_level:

  Confidence coverage for a symmetric interval (default `0.95`). Used
  only when both `alpha_lower` and `alpha_upper` are `NULL`.

- alpha_lower, alpha_upper:

  Optional Type I error rates for the lower and upper tail. If both are
  `NULL`, a symmetric interval at `conf_level` is used. If both are
  supplied, `conf_level` is recomputed as
  `1 - alpha_lower - alpha_upper`.

- ...:

  Additional arguments passed to
  [`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
  (for example `tol`).

## Value

A one-row `data.frame` with columns `sample_type` (`"one-sample"` or
`"two-sample"`), `D2` (point estimate of the squared distance),
`lower_limit` and `upper_limit` (the confidence limits on the population
squared distance \\\Delta^2\\), `F_value`, `df_1`, `df_2`, `n_1`, `n_2`
(`NA` in one-sample mode), and `p`.

## Details

**Definition.** For a \\p\\-vector \\\mathbf{x}\\ drawn from a
multivariate normal with mean \\\boldsymbol{\mu}\\ and covariance
\\\boldsymbol{\Sigma}\\, Mahalanobis's (1936) squared distance from a
reference vector \\\boldsymbol{\mu}\_0\\ is \$\$\Delta^2 =
(\boldsymbol{\mu} - \boldsymbol{\mu}\_0)^\top \boldsymbol{\Sigma}^{-1}
(\boldsymbol{\mu} - \boldsymbol{\mu}\_0).\$\$ In the two-sample case the
population distance between groups is \\\Delta^2 =
(\boldsymbol{\mu}\_1 - \boldsymbol{\mu}\_2)^\top
\boldsymbol{\Sigma}^{-1} (\boldsymbol{\mu}\_1 - \boldsymbol{\mu}\_2)\\,
assuming a common covariance. The corresponding sample estimates plug
the sample means and the sample (or pooled) covariance into the same
quadratic form.

**Link to Hotelling's \\T^2\\.** Hotelling's (1931) \\T^2\\ statistic is
\\T^2 = n D^2\\ (one sample) or \\T^2 = \\n_1 n_2 / (n_1 + n_2)\\ D^2\\
(two samples). Under multivariate normality \$\$\frac{n_1 + n_2 - p -
1}{(n_1 + n_2 - 2)\\p}\\T^2 \sim F'\\\left(p,\\ n_1 + n_2 - p - 1,\\
\lambda = \frac{n_1 n_2}{n_1 + n_2}\\\Delta^2\right)\$\$ in the
two-sample case, and analogously \\\\(n - p)/\[(n-1)p\]\\\\T^2 \sim
F'(p, n-p, n\Delta^2)\\ in the one-sample case (see Anderson, 2003,
Section 5.2).

**Confidence interval.** The CI on \\\Delta^2\\ is obtained by inverting
these distributional results (Reiser, 2001): a CI on the noncentrality
parameter \\\lambda\\ is constructed via
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
and then mapped back to \\\Delta^2\\ by \\\Delta^2 = \lambda\\(n_1 +
n_2)/(n_1 n_2)\\ (two sample) or \\\Delta^2 = \lambda / n\\ (one
sample). When the observed \\F\\ is below the lower-tail critical value
of the central *F*-distribution at the requested confidence level, the
lower CI on \\\lambda\\ (and hence on \\\Delta^2\\) is clamped to zero,
in keeping with the
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
convention.

**Bias.** The plug-in estimator \\D^2\\ is upward biased for
\\\Delta^2\\; the CI from this function is exact for \\\Delta^2\\ under
multivariate normality and reflects the bias structure correctly, but
the point estimate reported is the standard plug-in \\D^2\\.

## References

Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
Analysis* (4th ed.). Wiley.

Hotelling, H. (1931). The generalization of Student's ratio. *The Annals
of Mathematical Statistics, 2*(3), 360–378.

Mahalanobis, P. C. (1936). On the generalized distance in statistics.
*Proceedings of the National Institute of Sciences of India, 2*(1),
49–55.

Reiser, B. (2001). Confidence intervals for the Mahalanobis distance.
*Communications in Statistics–Simulation and Computation, 30*(1), 37–45.
[doi:10.1081/SAC-100001856](https://doi.org/10.1081/SAC-100001856)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

## See also

[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two-sample distance between the two iris species 'setosa' and
#     'versicolor' on the four sepal/petal measurements.
g1 <- as.matrix(iris[iris$Species == "setosa",     1:4])
g2 <- as.matrix(iris[iris$Species == "versicolor", 1:4])
ci_mahalanobis(group_1 = g1, group_2 = g2)
#>  sample_type D2  lower_limit upper_limit F_value df_1 df_2 n_1 n_2 p
#>  two-sample  103 72.6        132         625     4    95   50  50  4
#> 
#> Confidence level: 95%

# One-sample distance: how far is the centroid of versicolor from a
#     reference vector of (5.5, 3.0, 4.0, 1.5)?
ci_mahalanobis(group_1 = g2, mu_0 = c(5.5, 3.0, 4.0, 1.5))
#>  sample_type D2   lower_limit upper_limit F_value df_1 df_2 n_1 n_2  p
#>  one-sample  5.26 2.79        7.49        61.7    4    46   50  <NA> 4
#> 
#> Confidence level: 95%

# Pre-computed D^2 (no raw data needed).
ci_mahalanobis(D2 = 103.2, n_1 = 50, n_2 = 50, p = 4)
#>  sample_type D2  lower_limit upper_limit F_value df_1 df_2 n_1 n_2 p
#>  two-sample  103 72.6        131         625     4    95   50  50  4
#> 
#> Confidence level: 95%
```
