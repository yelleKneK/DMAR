# Tucker's Congruence Coefficient \\\phi\\ (Factor Similarity)

Computes Tucker's (1951) congruence coefficient \\\phi\\, a measure of
similarity between two factor-loading patterns (typically the
standardized loadings of the same factor estimated on two different
samples or with different methods), together with a permutation-based
*p*-value testing the null hypothesis of unrelated loading patterns.
\\\phi\\ is the standard tool for factor-replication studies
(Lorenzo-Seva & ten Berge, 2006).

## Usage

``` r
procrustes_phi(loadings_1, loadings_2, n_perm = 10000L)
```

## Arguments

- loadings_1, loadings_2:

  Numeric vectors of factor loadings on the same indicator set, of equal
  length. Either standardized or raw loadings work; the coefficient is
  scale-invariant.

- n_perm:

  Number of permutations for the significance test. Default `10000`. Set
  to 0 to skip the test.

## Value

A `data.frame` (class `dmar_tbl`) in `term` / `value` layout with the
row `tucker_phi`, the point estimate of \\\phi\\. When `n_perm > 0` the
table also carries `p_value_perm`, the two-sided permutation *p*-value,
and `n_perm`, the number of permutations requested.

## Details

**Definition.** For two vectors of loadings \\\bm\lambda_1,
\bm\lambda_2\\ on a shared set of \\p\\ indicators, Tucker's congruence
coefficient is \$\$\phi(\bm\lambda_1, \bm\lambda_2) \\=\\
\frac{\sum\_{i=1}^{p} \lambda\_{1i} \lambda\_{2i}}
{\sqrt{\sum\_{i=1}^{p} \lambda\_{1i}^2 \cdot \sum\_{i=1}^{p}
\lambda\_{2i}^2}}.\$\$ \\\phi\\ is the cosine of the angle between the
two loading vectors and ranges over \\\[-1, 1\]\\; values near \\\pm 1\\
indicate high (anti-)congruence, values near 0 indicate orthogonality.

**Permutation test.** Under the null hypothesis that the two loading
patterns are unrelated, randomly permuting one of the loading vectors
and recomputing \\\phi\\ produces a sampling distribution against which
the observed \\\phi\\ can be evaluated. The two-sided *p*-value is
\\(r + 1) / (m + 1)\\, where *r* counts the permuted \\\|\phi\|\\ values
at least as large as the observed \\\|\phi\|\\ and *m* is `n_perm`.
Adding one to each part counts the observed arrangement, which is itself
a legitimate permutation; without it a *p*-value of exactly zero could
be reported, a value a sampled permutation test cannot support (Phipson
& Smyth, 2010). The smallest reportable *p*-value is therefore \\1 /
(m + 1)\\.

**Interpretation.** Benchmark values for \\\phi\\ have been proposed in
the literature (Lorenzo-Seva & ten Berge, 2006), but context always
matters; this package reports the coefficient with its uncertainty and
leaves interpretation to the context of the application.

## References

Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
coefficient as a meaningful index of factor similarity. *Methodology,
2*(2), 57–64.
[doi:10.1027/1614-2241.2.2.57](https://doi.org/10.1027/1614-2241.2.2.57)

Phipson, B., & Smyth, G. K. (2010). Permutation *p*-values should never
be zero: Calculating exact *p*-values when permutations are randomly
drawn. *Statistical Applications in Genetics and Molecular Biology,
9*(1), Article 39.
[doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)

Tucker, L. R. (1951). *A method for synthesis of factor analysis
studies* (Personnel Research Section Report No. 984). Department of the
Army.

## See also

[`cor`](https://rdrr.io/r/stats/cor.html),
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)

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
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# 1. Two highly similar loading patterns:
l1 <- c(0.72, 0.65, 0.81, 0.55, 0.69)
l2 <- c(0.70, 0.62, 0.83, 0.58, 0.66)
procrustes_phi(l1, l2)
#>  term         value 
#>  tucker_phi   0.999 
#>  p_value_perm 0.0080
#>  n_perm       10000 

# 2. Loadings on different factors should show low congruence:
l3 <- c(0.10, 0.05, 0.20, 0.85, 0.78)
procrustes_phi(l1, l3, n_perm = 5000)
#>  term         value 
#>  tucker_phi   0.702 
#>  p_value_perm 0.8244
#>  n_perm       5000  
```
