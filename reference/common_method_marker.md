# Marker-Variable Adjustment for Common Method Variance

The marker-variable technique of Lindell and Whitney (2001) estimates
common method variance from the correlation of a *marker* variable that
is theoretically unrelated to the substantive constructs: any non-zero
correlation it shows is attributed to shared method, and that amount is
partialled out of the substantive correlations. When no a priori marker
is available, the smallest positive correlation among the substantive
items is used as a proxy, the common marker-free variant of the method.
Correlations that survive the adjustment are evidence that a
relationship is not an artifact of method variance.

## Usage

``` r
common_method_marker(R, marker_r = NULL)
```

## Arguments

- R:

  A correlation matrix among the substantive items.

- marker_r:

  The marker variable's (CMV) correlation. When `NULL` (default) the
  smallest positive off-diagonal correlation in `R` is used as the proxy
  marker.

## Value

A `data.frame` (class `dmar_tbl`) with rows `marker_correlation`,
`mean_abs_r_unadjusted`, and `mean_abs_r_adjusted` in the `value`
column. The full adjusted correlation matrix is the `"adjusted"`
attribute.

## Details

Writing \\r_M\\ for the marker (or proxy) correlation, each substantive
correlation is adjusted as \\r^{A}\_{ij} = (r\_{ij} - r_M) / (1 - r_M)\\
(Lindell & Whitney, 2001, Equation 3). The CMV-adjusted correlation
matrix is returned as the `"adjusted"` attribute; the reported table
summarizes the marker correlation and the average absolute correlation
before and after adjustment.

## References

Lindell, M. K., & Whitney, D. J. (2001). Accounting for common method
variance in cross-sectional research designs. *Journal of Applied
Psychology, 86*(1), 114–121.
[doi:10.1037/0021-9010.86.1.114](https://doi.org/10.1037/0021-9010.86.1.114)

## See also

[`common_method_single_factor`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md)
for the single-factor screen.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
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
R <- matrix(c(1, .5, .4, .5, 1, .45, .4, .45, 1), 3, 3,
            dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
res <- common_method_marker(R, marker_r = 0.10)
res
#>  term                  value
#>  marker_correlation    0.1  
#>  mean_abs_r_unadjusted 0.45 
#>  mean_abs_r_adjusted   0.389
attr(res, "adjusted")
#>           a         b         c
#> a 1.0000000 0.4444444 0.3333333
#> b 0.4444444 1.0000000 0.3888889
#> c 0.3333333 0.3888889 1.0000000
```
