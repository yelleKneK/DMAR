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
simple_structure(factor_loading_matrix, salient = 0.3, hyperplane = 0.1)
```

## Arguments

- factor_loading_matrix:

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

Item complexity is Hofmann's complexity index, proposed in Hofmann
(1977) and given as Equation 1 of Hofmann (1978), \\c_i = (\sum_j
\lambda\_{ij}^2)^2 / \sum_j \lambda\_{ij}^4\\, with \\\lambda\_{ij}\\
the loading of item \\i\\ on factor \\j\\ in `factor_loading_matrix`,
which equals one when an item loads on a single factor and rises toward
the number of factors as the loadings spread out; it is the same
complexity that `psych::fa` reports. The `"complexity"` attribute holds
these per-item values, and the `mean_complexity` row is their arithmetic
average, what Hofmann (1978) calls the total matrix complexity. These
are complexities, not Kaiser's (1974) simplicity index; Hofmann (1978)
shows that either can be derived from the other at the item level, with
the simplicity of item \\i\\ in an \\m\\-factor solution given by his
Equation 3, \\s_i = \[1/(m - 1)\]\[(m / c_i) - 1\]\\. An item is *pure*
when exactly one of its loadings is salient and *complex* when more than
one is. The hyperplane proportion is the fraction of all loadings whose
absolute value is below `hyperplane`; a clean simple structure is mostly
such near-zero loadings.

## References

Hofmann, R. J. (1977). Indices descriptive of factor complexity. *The
Journal of General Psychology, 96*(1), 103–110.
[doi:10.1080/00221309.1977.9920803](https://doi.org/10.1080/00221309.1977.9920803)

Hofmann, R. J. (1978). Complexity and simplicity as objective indices
descriptive of factor solutions. *Multivariate Behavioral Research,
13*(2), 247–250.
[doi:10.1207/s15327906mbr1302_9](https://doi.org/10.1207/s15327906mbr1302_9)

Kaiser, H. F. (1974). An index of factorial simplicity. *Psychometrika,
39*(1), 31–36.

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
[`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_adjusted_r()`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`fa_bifactor()`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
[`fa_second_order()`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
[`factor_congruence()`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
[`factor_scores()`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`sem_exact_fit_bootstrap()`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
[`sem_fit_measures()`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
[`sem_information_criteria()`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A nearly simple two-factor structure: six items, three per factor.
loading_matrix <- rbind(
  i1 = c(0.80, 0.05), i2 = c(0.75, 0.10), i3 = c(0.70, -0.05),
  i4 = c(0.08, 0.78), i5 = c(-0.04, 0.72), i6 = c(0.30, 0.60))
simple_structure(factor_loading_matrix = loading_matrix)
#>  term                  value
#>  items                 6    
#>  factors               2    
#>  mean_complexity       1.092
#>  median_complexity     1.016
#>  hyperplane_proportion 0.333
#>  n_pure                5    
#>  n_complex             1    
#>  proportion_pure       0.833
attr(simple_structure(loading_matrix), "complexity")
#>       i1       i2       i3       i4       i5       i6 
#> 1.007812 1.035544 1.010204 1.021036 1.006173 1.470588 
```
