# Lin's Concordance Correlation Coefficient

Computes Lin's (1989) concordance correlation coefficient (CCC) for a
pair of vectors of paired observations, together with a confidence
interval based on either Lin's (1989) original *z*-transformed standard
error or the King-Chinchilli (2001) variant that is corrected for
skewness when agreement is high. The CCC measures agreement (not merely
correlation) between two methods of measurement: a CCC of 1 means
perfect agreement (\\y_i = x_i\\ for all \\i\\), while Pearson's \\r\\
would still be 1 for any straight-line relationship, even one with
non-unit slope.

## Usage

``` r
lin_ccc(x, y, conf_level = 0.95, method = c("king_chinchilli", "lin"))
```

## Arguments

- x, y:

  Paired numeric vectors of equal length (e.g., the two measurement
  methods).

- conf_level:

  Confidence level for the CI. Default `0.95`.

- method:

  One of `"lin"` (the original Lin 1989 SE) or `"king_chinchilli"` (the
  King & Chinchilli 2001 SE corrected for skewness). The default is
  `"king_chinchilli"` because the original SE is biased when the
  underlying agreement is high; the two coincide for low-to-moderate
  CCC.

## Value

A `data.frame` with rows for the CCC point estimate, the lower and upper
CI limits, and decomposition components (Pearson \\r\\, accuracy
\\C_b\\, location-shift \\u\\, scale-shift \\v\\).

## Details

**Definition.** Lin (1989) defined the CCC as \$\$\rho_c \\=\\ \frac{2
\rho\\ \sigma_x \sigma_y} {\sigma_x^2 + \sigma_y^2 + (\mu_x -
\mu_y)^2},\$\$ where \\\rho = \mathrm{Cor}(X, Y)\\ is the Pearson
correlation and the denominator subtracts the disagreement attributable
to mean shift and scale shift. \\\rho_c\\ factors as \\\rho_c = \rho
\cdot C_b\\, where \\C_b \in \[0, 1\]\\ is the "bias correction factor"
that captures location and scale agreement, and \\C_b = 1\\ iff \\\mu_x
= \mu_y\\ and \\\sigma_x = \sigma_y\\.

**Lin (1989) SE.** The Fisher-style *z*-transform of the CCC, \\z =
\frac{1}{2} \log\\(1 + \rho_c)/(1 - \rho_c)\\\\, has approximate
variance \$\$\mathrm{Var}(z) \\\approx\\ \frac{1}{n - 2} \left\[
\frac{(1 - \rho^2) \rho_c^2}{(1 - \rho_c^2) \rho^2} + \frac{2 \rho_c^3
(1 - \rho_c) u^2}{\rho (1 - \rho_c^2)^2} - \frac{\rho_c^4 u^4}{2 \rho^2
(1 - \rho_c^2)^2}\right\],\$\$ where \\u = (\mu_x - \mu_y) /
\sqrt{\sigma_x \sigma_y}\\. The CI is built on the *z*-scale and
back-transformed via \\\tanh\\.

**King-Chinchilli (2001) SE.** King and Chinchilli pointed out that
Lin's variance is downward-biased when the agreement is high, which
produces over-narrow CIs. They derive a corrected variance that accounts
for the third and fourth moments of the joint distribution. We use their
general-form estimator (their equation 10), which converges to Lin's for
low CCC and is recommended for routine use.

## References

King, T. S., & Chinchilli, V. M. (2001). A generalized concordance
correlation coefficient for continuous and categorical data. *Statistics
in Medicine, 20*(14), 2131–2147.
[doi:10.1002/sim.845](https://doi.org/10.1002/sim.845)

Lin, L. I.-K. (1989). A concordance correlation coefficient to evaluate
reproducibility. *Biometrics, 45*(1), 255–268.

Lin, L. I.-K. (2000). A note on the concordance correlation coefficient.
*Biometrics, 56*(1), 324–325.
[doi:10.1111/j.0006-341X.2000.00324.x](https://doi.org/10.1111/j.0006-341X.2000.00324.x)

## See also

[`loa`](https://yelleknek.github.io/DMAR/reference/loa.md)

Other agreement and measurement:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
[`loa()`](https://yelleknek.github.io/DMAR/reference/loa.md),
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two methods of measuring the same quantity:
set.seed(113)
method_a <- rnorm(40, mean = 100, sd = 15)
method_b <- method_a + rnorm(40, mean = 2, sd = 5)
lin_ccc(method_a, method_b)
#>  term        value 
#>  ccc         0.928 
#>  lower_limit -0.977
#>  upper_limit 1     
#>  pearson_r   0.941 
#>  C_b         0.986 
#>  u           -0.154
#>  v           0.933 
#> 
#> Confidence level: 95%

# 2. Compare CCC with Pearson r when there is a systematic offset:
lin_ccc(method_a, method_a + 5)$value[1:2]   # CCC < r
#> [1]  0.9492315 -0.9957660
cor(method_a, method_a + 5)                  # Pearson r = 1
#> [1] 1
```
