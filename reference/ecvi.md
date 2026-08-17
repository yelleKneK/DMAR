# Expected Cross-Validation Index (ECVI) for a Covariance-Structure Model

The ECVI of Browne and Cudeck (1989) estimates how well a fitted model's
implied covariance matrix would fit an independent sample of the same
size from the same population. It is the single-sample estimate of the
cross-validation discrepancy, so a smaller ECVI indicates a model
expected to generalize better; ECVI is most useful for comparing
competing models fit to the same data. A confidence interval, derived
from the noncentral chi square distribution, accompanies the point
estimate.

## Usage

``` r
ecvi(
  fit = NULL,
  chisq = NULL,
  df = NULL,
  npar = NULL,
  n = NULL,
  conf_level = 0.95
)
```

## Arguments

- fit:

  A fitted lavaan model. Supply this, or the summary statistics below.

- chisq, df, npar, n:

  The model chi square, its degrees of freedom, the number of free
  parameters, and the total sample size. Used when `fit` is not
  supplied, so an ECVI can be obtained from a published fit table.

- conf_level:

  Confidence level for the interval. Defaults to 0.95.

## Value

A `data.frame` (class `dmar_tbl`) with rows `ecvi`, `lower_limit`, and
`upper_limit` in the `value` column.

## Details

With \\q\\ free parameters and total sample size \\N\\, \\\mathrm{ECVI}
= (\chi^2 + 2q)/N\\, the value the *Journal of Statistical Software*
reference implementation in lavaan reports. Writing \\\hat\lambda =
\chi^2 - df\\ for the estimated noncentrality, this is \\(\hat\lambda +
df + 2q)/N\\; the confidence interval replaces \\\hat\lambda\\ by the
lower and upper noncentrality limits from
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md),
the same inversion used for the RMSEA interval (see
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)).
ECVI differs from the AIC only by the constant factor \\N\\, so the two
rank models identically; ECVI is reported because its metric (a
discrepancy per observation) and its confidence interval are
interpretable on their own.

## References

Browne, M. W., & Cudeck, R. (1989). Single sample cross-validation
indices for covariance structures. *Multivariate Behavioral Research,
24*(4), 445–455.

## See also

[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md).

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# From a published fit table (no model object needed).
ecvi(chisq = 24.361, df = 8, npar = 13, n = 301)
#>  term        value
#>  ecvi        0.167
#>  lower_limit 0.125
#>  upper_limit 0.243
#> 
#> Confidence level: 95%

fit <- lavaan::cfa(
  "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
   verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
  data = holzinger_swineford, std.lv = TRUE)
ecvi(fit)
#>  term        value
#>  ecvi        0.167
#>  lower_limit 0.125
#>  upper_limit 0.243
#> 
#> Confidence level: 95%
```
