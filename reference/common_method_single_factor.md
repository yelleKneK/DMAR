# Single-Common-Factor Screen for Common Method Variance

This function implements Harman's single-factor test, the most widely
used (and weakest) screen for common method variance. Extract a single
unrotated factor from all of the items and inspect how much of their
variance it accounts for. The rationale is that if a single method
factor dominated the responses, one unrotated factor would capture a
large share of the common variance. A first factor accounting for more
than half of the variance is the customary red flag (Podsakoff,
MacKenzie, Lee, & Podsakoff, 2003). The screen is coarse and cannot by
itself rule method variance in or out; the marker-variable and
latent-method-factor approaches are stronger (see
[`common_method_marker`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md)).

## Usage

``` r
common_method_single_factor(data = NULL, S = NULL, R = NULL)
```

## Arguments

- data:

  A `data.frame` or numeric matrix of item responses. Supply this, a
  covariance matrix `S`, or a correlation matrix `R` (exactly one).

- S:

  A symmetric covariance matrix among the items, when raw data are not
  available but the summary statistics a paper reports are. It is
  converted to a correlation matrix internally, so the test acts on the
  same scale-free quantity regardless of which input is supplied.

- R:

  A correlation matrix among the items, when raw data are not available.

## Value

A `data.frame` (class `dmar_tbl`) with rows `variance_explained` (the
first-factor proportion) and `n_items` in the `value` column.

## Details

Harman's single-factor test (the proportion of variance explained by one
common factor) is related to but distinct from the marker-variable
technique. A *marker variable* (or common-method marker) is a variable
chosen to be theoretically unrelated to the substantive constructs under
study, so that any observed correlation between it and the substantive
items is attributable to shared method rather than to a true
relationship; it is used to estimate or partial out common method
variance (Lindell & Whitney, 2001). The single-factor test uses no such
marker, it asks only whether a single dimension dominates the item set,
so it can flag a strong common factor but cannot identify whether that
factor is method or substance.

The statistic is the proportion of total variance explained by the first
unrotated principal component, the largest eigenvalue of the item
correlation matrix divided by the number of items. Correlations from raw
data use pairwise-complete observations. A supplied covariance matrix is
first standardized to a correlation matrix with
[`cov2cor`](https://rdrr.io/r/stats/cor.html).

## References

Lindell, M. K., & Whitney, D. J. (2001). Accounting for common method
variance in cross-sectional research designs. *Journal of Applied
Psychology, 86*(1), 114–121.
[doi:10.1037/0021-9010.86.1.114](https://doi.org/10.1037/0021-9010.86.1.114)

Podsakoff, P. M., MacKenzie, S. B., Lee, J.-Y., & Podsakoff, N. P.
(2003). Common method biases in behavioral research: A critical review
of the literature and recommended remedies. *Journal of Applied
Psychology, 88*(5), 879–903.
[doi:10.1037/0021-9010.88.5.879](https://doi.org/10.1037/0021-9010.88.5.879)

## See also

[`common_method_marker`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md)
for the marker-variable adjustment.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
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
set.seed(113)
f <- rnorm(200)
d <- data.frame(
  x1 = f + rnorm(200), x2 = f + rnorm(200), x3 = f + rnorm(200),
  x4 = rnorm(200),     x5 = rnorm(200),     x6 = rnorm(200))
common_method_single_factor(d)
#>  term               value
#>  variance_explained 0.327
#>  n_items            6    

# The same screen from the summary statistics a paper reports.
common_method_single_factor(S = cov(d))
#>  term               value
#>  variance_explained 0.327
#>  n_items            6    
```
