# Bifactor Model Dimensionality and Reliability Indices

Computes the indices used to judge whether a multidimensional scale is
nonetheless unidimensional enough to score as a single total (Rodriguez,
Reise, and Haviland, 2016): the explained common variance (ECV),
coefficient omega and omega hierarchical (omega_H) for the general
factor, omega hierarchical subscale (omega_HS) for each group factor,
the percentage of uncontaminated correlations (PUC), and coefficient H,
the construct reliability (maximal reliability) of Hancock and Mueller
(2001). The input is a fitted bifactor model in which one general factor
loads on every item and each item loads on exactly one orthogonal group
factor.

## Usage

``` r
bifactor_indices(fit, general = NULL)
```

## Arguments

- fit:

  A fitted bifactor lavaan model: one general factor on all items plus
  orthogonal group factors, each item on one group factor. The factors
  must be orthogonal (the bifactor specification).

- general:

  Optional name of the general factor. When `NULL` (default) the general
  factor is detected as the one that loads on every item.

## Value

A `data.frame` (class `dmar_tbl`) with one row per factor (the general
factor first, then each group factor) and columns `factor`, `ECV`,
`omega`, `omega_H`, `omega_HS`, `PUC`, and `H`. Quantities that do not
apply to a row are `NA` (for example omega_H and PUC on a group factor).
The `"improper"` attribute flags a Heywood solution.

## Details

Let \\\lambda^g_i\\ be item \\i\\'s standardized loading on the general
factor, \\\lambda^s_i\\ its loading on its group factor, and
\\\theta_i\\ its standardized residual variance. With orthogonal factors
the overall ECV is \\\sum_i (\lambda^g_i)^2 / \sum_i \[(\lambda^g_i)^2 +
(\lambda^s_i)^2\]\\; omega and omega_H share the total-score variance
\\(\sum_i \lambda^g_i)^2 + \sum_g (\sum\_{i \in g} \lambda^s_i)^2 +
\sum_i \theta_i\\ as denominator, with the general part \\(\sum_i
\lambda^g_i)^2\\ in the numerator of omega_H. Each group factor's
omega_HS uses the analogous numerator and that subscale's own total
variance. PUC is one minus the share of item pairs that fall within the
same group factor. Coefficient H is computed from \\\sum \lambda^2 /
(1 - \lambda^2)\\ on the relevant loadings (see
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)).
High ECV and PUC with a high omega_H support scoring a single total;
substantial subscale omega_HS argues for reporting subscales as well.

A bifactor model that is over-parameterized for the data frequently
yields an improper (Heywood) solution – a standardized loading outside
\\\[-1, 1\]\\ or a negative residual variance – whose indices are not
trustworthy. When that happens the indices are still returned but a
warning is issued and the `"improper"` attribute is set to `TRUE`.

## References

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Reise, S. P. (2012). The rediscovery of bifactor measurement models.
*Multivariate Behavioral Research, 47*(5), 667–696.
[doi:10.1080/00273171.2012.715555](https://doi.org/10.1080/00273171.2012.715555)

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
bifactor models: Calculating and interpreting statistical indices.
*Psychological Methods, 21*(2), 137–150.

## See also

[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md).

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
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
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Nine items: one general factor and three orthogonal group factors.
set.seed(113)
n <- 600
g <- rnorm(n); grp <- list(rnorm(n), rnorm(n), rnorm(n))
X <- vapply(1:9, function(i)
  0.5 * g + 0.5 * grp[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n), numeric(n))
colnames(X) <- paste0("x", 1:9)
model <- "g  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
          f1 =~ x1 + x2 + x3
          f2 =~ x4 + x5 + x6
          f3 =~ x7 + x8 + x9"
fit <- lavaan::cfa(model, data = as.data.frame(X),
                   orthogonal = TRUE, std.lv = TRUE)
bifactor_indices(fit)
#>  factor ECV   omega omega_H omega_HS PUC  H    
#>  g      0.547 0.858 0.673   <NA>     0.75 0.769
#>  f1     0.609 0.727 <NA>    0.28     <NA> 0.412
#>  f2     0.543 0.756 <NA>    0.345    <NA> 0.477
#>  f3     0.49  0.743 <NA>    0.377    <NA> 0.506
```
