# Correct a Correlation for Attenuation Due to Measurement Error

Applies the Spearman (1904) correction for attenuation: the observed
correlation between two fallible measures understates the correlation
between the constructs they measure, and dividing by the square root of
the product of the two reliabilities recovers it, \$\$r_c \\=\\
\frac{r\_{XY}}{\sqrt{\rho\_{XX'}\\\rho\_{YY'}}}.\$\$ Within classical
test theory (Lord & Novick, 1968), the disattenuated correlation
estimates the correlation between the true scores, that is, how strongly
the two constructs would correlate if each were measured without error.
When `N` is supplied, a confidence interval for the corrected
correlation is formed by disattenuating the endpoints of the Fisher *z*
interval for the observed correlation, the standard practice when the
reliabilities are treated as known.

## Usage

``` r
correction_for_attenuation(
  r,
  reliability_x,
  reliability_y,
  N = NULL,
  conf_level = 0.95
)
```

## Arguments

- r:

  The observed correlation between the two measures, in \\\[-1, 1\]\\.

- reliability_x:

  Reliability of the first measure, in \\(0, 1\]\\. Any reliability
  estimate appropriate to the use may be supplied (for example,
  coefficient alpha or omega from
  [`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)).

- reliability_y:

  Reliability of the second measure, in \\(0, 1\]\\. For a correlation
  between a measure and an error-free criterion, supply 1 for that side
  (correcting for criterion unreliability only yields what the validity
  generalization literature calls the operational validity).

- N:

  Optional sample size on which `r` is based. When supplied, a
  `conf_level` confidence interval for the corrected correlation is
  reported by correcting the endpoints of the Fisher *z* interval for
  `r`.

- conf_level:

  Confidence level for the interval when `N` is supplied. Defaults to
  0.95.

## Value

A `data.frame` (class `dmar_tbl`) in `term` / `value` layout with the
observed correlation (`correlation_observed`), the corrected correlation
(`correlation_corrected`), the corrected interval (`lower_limit`,
`upper_limit`; present only when `N` is supplied), and the
`reliability_x`, `reliability_y`, and `N` inputs.

## Details

The correction treats the two reliabilities as known constants, which is
the conventional assumption; uncertainty in the reliabilities themselves
would widen the interval further. Because the observed correlation can
exceed what the supplied reliabilities allow (sampling error, or
reliabilities that understate the truth), the corrected value can exceed
1 in magnitude; when that happens the value is reported as computed,
with a warning, rather than silently truncated, since a corrected
correlation beyond 1 is itself diagnostic information about the inputs.

**Prefer the factor model when you have the items.** The Spearman
formula is the summary-statistics route: it is exactly right when all
you have are the observed correlation and reliability estimates. When
the item-level data are available, the better practice is to estimate
the construct-level correlation directly as the factor correlation in a
two-factor model (each scale loading on its own factor, factors free to
correlate): the latent correlation is then estimated jointly with the
measurement model rather than assembled from plug-in reliabilities, and
it comes with a standard error that propagates the sampling variability
of all the moving parts. The example below shows both routes on the same
data, the formula route using
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
the model route using lavaan; with congeneric items the two agree
closely, and when they disagree the factor model is the one to trust.

## References

Lord, F. M., & Novick, M. R. (1968). *Statistical theories of mental
test scores*. Addison-Wesley.

Spearman, C. (1904). The proof and measurement of association between
two things. *The American Journal of Psychology, 15*(1), 72–101.

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
and its family for estimating the reliabilities supplied here;
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) and
lavaan for the latent variable route the Details recommend when items
are available;
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_r.md) for
inference on the observed correlation itself;
[`convert_r_z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)
and
[`convert_z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md)
for the Fisher transformation the interval uses.

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
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
# An observed correlation of .30 between measures with reliabilities .80
# and .70 corresponds to a construct-level correlation of about .40.
correction_for_attenuation(r = 0.30, reliability_x = 0.80, reliability_y = 0.70)
#>  term                  value
#>  correlation_observed  0.3  
#>  correlation_corrected 0.401
#>  reliability_x         0.8  
#>  reliability_y         0.7  

# With the sample size, the corrected interval comes along.
correction_for_attenuation(r = 0.30, reliability_x = 0.80, reliability_y = 0.70,
                    N = 120)
#>  term                  value
#>  correlation_observed  0.3  
#>  correlation_corrected 0.401
#>  lower_limit           0.171
#>  upper_limit           0.608
#>  reliability_x         0.8  
#>  reliability_y         0.7  
#>  N                     120  
#> 
#> Confidence level: 95%

# Correct one side only (error-free criterion).
correction_for_attenuation(r = 0.30, reliability_x = 0.80, reliability_y = 1)
#>  term                  value
#>  correlation_observed  0.3  
#>  correlation_corrected 0.335
#>  reliability_x         0.8  
#>  reliability_y         1    

# The two routes to the construct-level correlation, on the same data
# (requires lavaan). Two congeneric scales of three items each whose
# latent variables correlate .50:
set.seed(113)
n <- 400
fx <- rnorm(n); fy <- 0.5 * fx + sqrt(1 - 0.25) * rnorm(n)
lam <- c(.8, .7, .6)
items <- data.frame(
  x1 = lam[1] * fx + rnorm(n, 0, sqrt(1 - lam[1]^2)),
  x2 = lam[2] * fx + rnorm(n, 0, sqrt(1 - lam[2]^2)),
  x3 = lam[3] * fx + rnorm(n, 0, sqrt(1 - lam[3]^2)),
  y1 = lam[1] * fy + rnorm(n, 0, sqrt(1 - lam[1]^2)),
  y2 = lam[2] * fy + rnorm(n, 0, sqrt(1 - lam[2]^2)),
  y3 = lam[3] * fy + rnorm(n, 0, sqrt(1 - lam[3]^2)))

# Route 1, summary statistics: omega reliabilities into the formula.
x_score <- rowMeans(items[, 1:3]); y_score <- rowMeans(items[, 4:6])
om_x <- reliability_omega(data = items[, 1:3])$value[1]
#> Robust omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "percentile" (or "bca"); B = 10000 replications is the default when you do.
om_y <- reliability_omega(data = items[, 4:6])$value[1]
#> Robust omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "percentile" (or "bca"); B = 10000 replications is the default when you do.
correction_for_attenuation(r = cor(x_score, y_score),
                           reliability_x = om_x, reliability_y = om_y,
                           N = n)
#>  term                  value
#>  correlation_observed  0.307
#>  correlation_corrected 0.413
#>  lower_limit           0.29 
#>  upper_limit           0.529
#>  reliability_x         0.76 
#>  reliability_y         0.725
#>  N                     400  
#> 
#> Confidence level: 95%

# Route 2, the factor model: the latent correlation estimated directly.
fit <- lavaan::cfa("X =~ x1 + x2 + x3\nY =~ y1 + y2 + y3",
                   data = items, std.lv = TRUE)
lavaan::parameterEstimates(fit)[
  lavaan::parameterEstimates(fit)$op == "~~" &
  lavaan::parameterEstimates(fit)$lhs == "X" &
  lavaan::parameterEstimates(fit)$rhs == "Y", ]
#>    lhs op rhs  est    se     z pvalue ci.lower ci.upper
#> 15   X ~~   Y 0.43 0.056 7.661      0     0.32    0.539
```
