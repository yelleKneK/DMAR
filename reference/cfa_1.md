# One Factor Confirmatory Factor Analysis Model

Fits a single factor (congeneric by default) confirmatory factor
analysis model to raw item data or a sample covariance matrix. This is
the one factor special case of
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md): the
function is a convenience wrapper that only requires the data (and, when
the data hold more than the items, a vector of item names), builds the
one factor specification, and forwards everything else to
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md). The
factor is named `f1`, so the rows of the returned table are
`lambda_f1_1`, `lambda_f1_2`, ..., `psi_f1_1`, ..., and `omega_f1`,
exactly as a one factor
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) call
would report them (the `syntax` column names the item behind each
number).

## Usage

``` r
cfa_1(
  data = NULL,
  items = NULL,
  S = NULL,
  N = NULL,
  equal_loading = FALSE,
  equal_error = FALSE,
  estimator = "ML",
  missing = "listwise",
  se = "standard",
  conf_level = 0.95,
  output = c("verbose", "measurement", "summary", "standardized", "fit"),
  ...
)
```

## Arguments

- data:

  A raw data matrix or data frame, rows are respondents and columns
  include the items. A matrix without column names is given the names
  `y1`, `y2`, ... Supply exactly one of `data` or `S`.

- items:

  Character vector naming the items of the factor (three or more; two
  are accepted with `equal_loading = TRUE`, the just identified
  tau-equivalent case). The default `NULL` uses every column of `data`
  (or of `S`), so a data set that holds only the items needs no `items`
  at all.

- S:

  A symmetric covariance matrix of the items; `N` is then required.
  Dimnames are optional here: a matrix without them is given the item
  names `y1`, `y2`, ... Supply exactly one of `data` or `S`.

- N:

  Total sample size. Required with `S`; ignored (inferred from the rows)
  with `data`.

- equal_loading:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the loadings within a factor to a single value (across
  factors nothing is equated). Defaults to `FALSE`.

- equal_error:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the error variances within a factor to a single value.
  Defaults to `FALSE`.

- estimator:

  Character; estimator passed to lavaan. Must be one of `"ML"` (default;
  maximum likelihood, fully efficient under multivariate normality),
  `"MLR"` (robust maximum likelihood: maximum likelihood estimates with
  standard errors and test statistic corrected for nonnormality; Satorra
  & Bentler, 1994), `"WLS"` (the asymptotic distribution free estimator
  of Browne, 1984; raw data and a large sample required), `"WLSMV"`
  (diagonally weighted least squares with mean- and variance-adjusted
  test statistic; Muthén, 1984; Muthén, du Toit, & Spisic, 1997; the
  standard choice for ordered categorical items, and what `ordered`
  switches to), or `"GLS"` (generalized least squares; Browne, 1974).
  With a robust estimator the reported fit indices are the robust
  versions.

- missing:

  Character; missing-data handling passed to lavaan when raw data are
  supplied. Common values are `"listwise"` (default, listwise deletion)
  and `"ml"` (full information maximum likelihood). Ignored with `S`.
  With `ordered` items, `"ml"`/`"fiml"` are not available; use
  `"pairwise"` or `"listwise"`.

- se:

  Standard error type passed to lavaan; see
  [`cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html). Common values are
  `"standard"` (default), `"robust.sem"` (with `estimator = "MLR"`), and
  `"none"` (point estimates only; fastest).

- conf_level:

  Confidence level for the parameter confidence intervals, including the
  delta method intervals for omega, AVE, and *H*. Defaults to 0.95. The
  RMSEA interval is a separate convention (see Details).

- output:

  Format of the returned object:

  `"verbose"`

  :   (default) Parameter estimates with confidence intervals, the
      per-factor defined parameters (`loading_sum`, `error_sum`,
      `omega`, `ave`, `H`), and fit information.

  `"measurement"`

  :   The measurement-property rows only: per factor `omega`, `ave`, and
      `H` (with delta method standard errors and confidence intervals),
      the latent correlation `phi` for every factor pair (with its
      confidence interval), and, for raw data, the heterotrait-monotrait
      ratio `htmt` for every factor pair via
      [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md).

  `"summary"`

  :   The raw [`summary()`](https://rdrr.io/r/base/summary.html) output
      from lavaan (not a data frame).

  `"standardized"`

  :   The standardized parameter estimates from
      [`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html).

  `"fit"`

  :   The raw lavaan fit object. Escape hatch for direct lavaan access,
      including likelihood ratio tests between two
      [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)
      fits via
      [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html).

- ...:

  Additional arguments forwarded to
  [`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) and,
  through it, to [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

## Value

The value of the corresponding
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) call: a
`data.frame` (classes `dmar_cfa_k`, `dmar_tbl`) with one row per
parameter (`estimate`, `se`, `z_value`, `p_value`, `ci_lower`,
`ci_upper`) followed by the fit rows, or the alternative shapes selected
by `output` (see
[`?cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)).

## Details

The model is identified by fixing the factor variance to 1 and
estimating every loading. `equal_loading` and `equal_error` impose the
classical measurement structures on the single factor, and the header of
the printed table names the structure implied by the constraints. For
the composite reliability coefficient with the observed total variance
in the denominator, use
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
with `denominator = "observed"`; for the model implied omega, the
`omega_f1` row of this function's output and
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
agree.

Ordered categorical items are not supported here; use
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md), whose
`ordered` argument fits WLSMV with the theta parameterization and
reports the Green and Yang (2009) categorical sum score omega, or
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md).

## References

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*(1), 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*(1), 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

McDonald, R. P. (1999). *Test theory: A unified treatment*. Lawrence
Erlbaum Associates.

## See also

[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) for the
general function this wraps (factor analysis with any number of factors,
intercept constraints, ordered categorical items, and the measurement
output); [`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md)
for the two factor wrapper;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for coefficient omega with confidence intervals.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
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
set.seed(113)
f <- rnorm(200)
loadings <- c(0.5, 0.6, 0.65, 0.7, 0.8)
X <- sapply(loadings, function(l) l * f + rnorm(200, sd = sqrt(1 - l^2)))
colnames(X) <- paste0("y", 1:5)

# All columns are items, so the data are all that is needed.
cfa_1(X)
#> Measurement structure, per factor:
#>   f1: congeneric (no equality constraints)
#> 
#>  syntax   term           estimate  se     z_value p_value  ci_lower ci_upper
#>  f1 =~ y1 lambda_f1_1    0.454     0.0708 6.4     < 0.0001 0.315    0.592   
#>  f1 =~ y2 lambda_f1_2    0.553     0.0721 7.67    < 0.0001 0.412    0.695   
#>  f1 =~ y3 lambda_f1_3    0.644     0.071  9.07    < 0.0001 0.505    0.783   
#>  f1 =~ y4 lambda_f1_4    0.732     0.0677 10.8    < 0.0001 0.599    0.865   
#>  f1 =~ y5 lambda_f1_5    0.796     0.0649 12.3    < 0.0001 0.669    0.924   
#>  f1 ~~ f1 phi_f1         1         0      <NA>    <NA>     1        1       
#>  y1 ~~ y1 psi_f1_1       0.721     0.0771 9.36    < 0.0001 0.57     0.872   
#>  y2 ~~ y2 psi_f1_2       0.7       0.0777 9       < 0.0001 0.547    0.852   
#>  y3 ~~ y3 psi_f1_3       0.614     0.0729 8.43    < 0.0001 0.472    0.757   
#>  y4 ~~ y4 psi_f1_4       0.46      0.0643 7.15    < 0.0001 0.334    0.586   
#>  y5 ~~ y5 psi_f1_5       0.324     0.0599 5.4     < 0.0001 0.206    0.441   
#>           loading_sum_f1 3.18      0.204  15.6    < 0.0001 2.78     3.58    
#>           error_sum_f1   2.82      0.142  19.8    < 0.0001 2.54     3.1     
#>           omega_f1       0.782     0.0242 32.3    < 0.0001 0.735    0.829   
#>           ave_f1         0.426     0.0317 13.4    < 0.0001 0.364    0.488   
#>           H_f1           0.819     0.023  35.6    < 0.0001 0.774    0.864   
#>           chi_square     3.61      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           df             5         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           p_chi_square   0.6061    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           cfi            1         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           tli            1.01      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           nnfi           1.01      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea          0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_lower 0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_upper 0.083     <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_level 0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>           srmr           0.0182    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           AIC            2584.466  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           BIC            2617.449  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H0             -1282.233 <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H1             -1280.426 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# Equal loadings (essentially tau-equivalent), named in the header.
cfa_1(X, equal_loading = TRUE)
#> Measurement structure, per factor:
#>   f1: essentially tau-equivalent (equal loadings)
#> 
#>  syntax   term           estimate  se     z_value p_value  ci_lower ci_upper
#>  f1 =~ y1 lambda_f1      0.654     0.0414 15.8    < 0.0001 0.573    0.735   
#>  f1 =~ y2 lambda_f1      0.654     0.0414 15.8    < 0.0001 0.573    0.735   
#>  f1 =~ y3 lambda_f1      0.654     0.0414 15.8    < 0.0001 0.573    0.735   
#>  f1 =~ y4 lambda_f1      0.654     0.0414 15.8    < 0.0001 0.573    0.735   
#>  f1 =~ y5 lambda_f1      0.654     0.0414 15.8    < 0.0001 0.573    0.735   
#>  f1 ~~ f1 phi_f1         1         0      <NA>    <NA>     1        1       
#>  y1 ~~ y1 psi_f1_1       0.699     0.0801 8.72    < 0.0001 0.542    0.856   
#>  y2 ~~ y2 psi_f1_2       0.662     0.0765 8.65    < 0.0001 0.512    0.812   
#>  y3 ~~ y3 psi_f1_3       0.609     0.0714 8.53    < 0.0001 0.469    0.749   
#>  y4 ~~ y4 psi_f1_4       0.491     0.0601 8.17    < 0.0001 0.373    0.609   
#>  y5 ~~ y5 psi_f1_5       0.417     0.0532 7.84    < 0.0001 0.313    0.521   
#>           loading_sum_f1 3.27      0.207  15.8    < 0.0001 2.87     3.68    
#>           error_sum_f1   2.88      0.145  19.8    < 0.0001 2.59     3.16    
#>           omega_f1       0.788     0.0235 33.5    < 0.0001 0.742    0.834   
#>           ave_f1         0.432     0.0343 12.6    < 0.0001 0.364    0.499   
#>           H_f1           0.794     0.0232 34.3    < 0.0001 0.749    0.84    
#>           chi_square     22.8      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           df             9         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           p_chi_square   0.0068    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           cfi            0.945     <NA>   <NA>    <NA>     <NA>     <NA>    
#>           tli            0.939     <NA>   <NA>    <NA>     <NA>     <NA>    
#>           nnfi           0.939     <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea          0.0874    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_lower 0.0432    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_upper 0.133     <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_level 0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>           srmr           0.111     <NA>   <NA>    <NA>     <NA>     <NA>    
#>           AIC            2595.614  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           BIC            2615.404  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H0             -1291.807 <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H1             -1280.426 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# From a covariance matrix and sample size, as a paper reports them.
cfa_1(S = cov(X), N = 200)
#> Measurement structure, per factor:
#>   f1: congeneric (no equality constraints)
#> 
#>  syntax   term           estimate  se     z_value p_value  ci_lower ci_upper
#>  f1 =~ y1 lambda_f1_1    0.454     0.0708 6.4     < 0.0001 0.315    0.592   
#>  f1 =~ y2 lambda_f1_2    0.553     0.0721 7.67    < 0.0001 0.412    0.695   
#>  f1 =~ y3 lambda_f1_3    0.644     0.071  9.07    < 0.0001 0.505    0.783   
#>  f1 =~ y4 lambda_f1_4    0.732     0.0677 10.8    < 0.0001 0.599    0.865   
#>  f1 =~ y5 lambda_f1_5    0.796     0.0649 12.3    < 0.0001 0.669    0.924   
#>  f1 ~~ f1 phi_f1         1         0      <NA>    <NA>     1        1       
#>  y1 ~~ y1 psi_f1_1       0.721     0.0771 9.36    < 0.0001 0.57     0.872   
#>  y2 ~~ y2 psi_f1_2       0.7       0.0777 9       < 0.0001 0.547    0.852   
#>  y3 ~~ y3 psi_f1_3       0.614     0.0729 8.43    < 0.0001 0.472    0.757   
#>  y4 ~~ y4 psi_f1_4       0.46      0.0643 7.15    < 0.0001 0.334    0.586   
#>  y5 ~~ y5 psi_f1_5       0.324     0.0599 5.4     < 0.0001 0.206    0.441   
#>           loading_sum_f1 3.18      0.204  15.6    < 0.0001 2.78     3.58    
#>           error_sum_f1   2.82      0.142  19.8    < 0.0001 2.54     3.1     
#>           omega_f1       0.782     0.0242 32.3    < 0.0001 0.735    0.829   
#>           ave_f1         0.426     0.0317 13.4    < 0.0001 0.364    0.488   
#>           H_f1           0.819     0.023  35.6    < 0.0001 0.774    0.864   
#>           chi_square     3.61      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           df             5         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           p_chi_square   0.6061    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           cfi            1         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           tli            1.01      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           nnfi           1.01      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea          0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_lower 0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_upper 0.083     <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_level 0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>           srmr           0.0182    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           AIC            2584.466  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           BIC            2617.449  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H0             -1282.233 <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H1             -1280.426 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# A subset of columns via \'items\'.
cfa_1(X, items = c("y1", "y2", "y3", "y4"))
#> Measurement structure, per factor:
#>   f1: congeneric (no equality constraints)
#> 
#>  syntax   term           estimate  se     z_value p_value  ci_lower ci_upper
#>  f1 =~ y1 lambda_f1_1    0.447     0.0768 5.81    < 0.0001 0.296    0.597   
#>  f1 =~ y2 lambda_f1_2    0.606     0.0789 7.67    < 0.0001 0.451    0.761   
#>  f1 =~ y3 lambda_f1_3    0.605     0.0798 7.58    < 0.0001 0.448    0.761   
#>  f1 =~ y4 lambda_f1_4    0.732     0.0798 9.18    < 0.0001 0.576    0.889   
#>  f1 ~~ f1 phi_f1         1         0      <NA>    <NA>     1        1       
#>  y1 ~~ y1 psi_f1_1       0.727     0.0818 8.89    < 0.0001 0.567    0.888   
#>  y2 ~~ y2 psi_f1_2       0.639     0.0848 7.53    < 0.0001 0.473    0.805   
#>  y3 ~~ y3 psi_f1_3       0.663     0.0868 7.64    < 0.0001 0.493    0.833   
#>  y4 ~~ y4 psi_f1_4       0.46      0.0893 5.15    < 0.0001 0.285    0.635   
#>           loading_sum_f1 2.39      0.174  13.8    < 0.0001 2.05     2.73    
#>           error_sum_f1   2.49      0.144  17.3    < 0.0001 2.21     2.77    
#>           omega_f1       0.696     0.0348 20      < 0.0001 0.628    0.765   
#>           ave_f1         0.369     0.0366 10.1    < 0.0001 0.297    0.44    
#>           H_f1           0.72      0.0367 19.6    < 0.0001 0.648    0.792   
#>           chi_square     0.0393    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           df             2         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           p_chi_square   0.9805    <NA>   <NA>    <NA>     <NA>     <NA>    
#>           cfi            1         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           tli            1.05      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           nnfi           1.05      <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea          0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_lower 0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_upper 0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>           rmsea_ci_level 0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>           srmr           0.00293   <NA>   <NA>    <NA>     <NA>     <NA>    
#>           AIC            2149.746  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           BIC            2176.133  <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H0             -1066.873 <NA>   <NA>    <NA>     <NA>     <NA>    
#>           H1             -1066.853 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%
```
