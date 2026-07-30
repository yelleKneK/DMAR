# Quantify Simple Structure in a Factor Loading Matrix

Summarizes how closely a rotated loading matrix approaches Thurstone's
simple structure, in which each item loads on as few factors as possible
so that the factors are interpretable. Three complementary quantities
are reported: the mean item complexity (the average number of factors an
item effectively loads on, one for a perfectly simple item), the
hyperplane proportion (the share of loadings near zero, which Thurstone
sought to maximize), and the counts of pure versus complex items at a
salience cutoff. Together they turn a visual impression of a loading
matrix into numbers.

## Usage

``` r
simple_structure(Lambda, salient = 0.3, hyperplane = 0.1)
```

## Arguments

- Lambda:

  A numeric matrix of factor loadings, items in rows and factors in
  columns (for example `unclass(psych::fa(...)$loadings)` or a lavaan
  standardized loading matrix). Row names, if present, label the items.

- salient:

  Absolute loading at or above which an item is counted as loading
  saliently on a factor. Defaults to 0.30 (about ten percent of an
  item's variance), a common floor for a meaningful loading.

- hyperplane:

  Absolute loading below which a loading is treated as lying in the
  hyperplane (effectively zero). Defaults to 0.10.

## Value

A `data.frame` (class `dmar_tbl`) with one row per summary quantity
(`term`, `value`): the number of items and factors, the mean and median
item complexity, the hyperplane proportion, and the counts and
proportion of pure items. The per-item complexities are attached as the
`"complexity"` attribute (a named numeric vector), and the salience and
hyperplane cutoffs as the `"salient"` and `"hyperplane"` attributes.

## Details

Item complexity is the index of Hoffman (1978), \\c_i = (\sum_j
\lambda\_{ij}^2)^2 / \sum_j \lambda\_{ij}^4\\, which equals one when an
item loads on a single factor and rises toward the number of factors as
the loadings spread out; it is the same complexity that
[`psych::fa`](https://rdrr.io/pkg/psych/man/fa.html) reports. An item is
*pure* when exactly one of its loadings is salient and *complex* when
more than one is. The hyperplane proportion is the fraction of all
loadings whose absolute value is below `hyperplane`; a clean simple
structure is mostly such near-zero loadings.

## References

Hoffman, P. J. (1978). Hierarchical factoring schema by item complexity.
*Multivariate Behavioral Research, 13*(1), 91–104.

Thurstone, L. L. (1947). *Multiple-factor analysis*. University of
Chicago Press.

## See also

[`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
and [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) for the
convergent and discriminant sides of an exploratory solution.

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
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A nearly simple two-factor structure: six items, three per factor.
Lambda <- rbind(
  i1 = c(0.80, 0.05), i2 = c(0.75, 0.10), i3 = c(0.70, -0.05),
  i4 = c(0.08, 0.78), i5 = c(-0.04, 0.72), i6 = c(0.30, 0.60))
simple_structure(Lambda)
#>  term                  value
#>  items                 6    
#>  factors               2    
#>  mean_complexity       1.09 
#>  median_complexity     1.02 
#>  hyperplane_proportion 0.333
#>  n_pure                5    
#>  n_complex             1    
#>  proportion_pure       0.833
attr(simple_structure(Lambda), "complexity")
#>       i1       i2       i3       i4       i5       i6 
#> 1.007812 1.035544 1.010204 1.021036 1.006173 1.470588 
```
