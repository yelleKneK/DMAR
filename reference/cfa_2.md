# Two Factor Confirmatory Factor Analysis Model

Fits a two factor confirmatory factor analysis model to raw item data or
a sample covariance matrix. This is the two factor special case of
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md): the
function is a convenience wrapper that only requires the items of each
factor, builds the two factor specification, and forwards everything
else to
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md). The
factors are named `f1` and `f2`, so the rows of the returned table are
`lambda_f1_1`, `lambda_f2_1`, `phi_f1_f2` (the factor correlation),
`omega_f1`, `omega_f2`, and so on, exactly as a two factor
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) call
would report them (the `syntax` column names the item behind each
number). To name the factors substantively, or for three or more
factors, call
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) directly.

## Usage

``` r
cfa_2(
  data = NULL,
  factor_1,
  factor_2,
  S = NULL,
  N = NULL,
  equal_loading = FALSE,
  equal_error = FALSE,
  correlated_factors = TRUE,
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
  include the items named in `factor_1` and `factor_2`. Supply exactly
  one of `data` or `S`.

- factor_1:

  Character vector naming the items of the first factor (two or more).

- factor_2:

  Character vector naming the items of the second factor (two or more).
  No item may appear in both factors.

- S:

  A symmetric covariance matrix of the items, with dimnames naming the
  items; `N` is then required. Supply exactly one of `data` or `S`.

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

- correlated_factors:

  Logical. If `TRUE` (default) the factors covary freely; because each
  factor variance is fixed to 1, the `phi` terms for factor pairs are
  the latent correlations. If `FALSE` the factor covariances are fixed
  to zero.

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
  through it, to [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html)
  (for example `ordered`, `equal_intercept`, or `M`).

## Value

The value of the corresponding
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) call: a
`data.frame` (classes `dmar_cfa_k`, `dmar_tbl`) with one row per
parameter (`estimate`, `se`, `z_value`, `p_value`, `ci_lower`,
`ci_upper`) followed by the fit rows, or the alternative shapes selected
by `output` (see
[`?cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)).

## Details

Each factor is identified by fixing its variance to 1 and estimating
every loading; each item loads on exactly one factor (simple structure).
The factor correlation is estimated by default and
`correlated_factors = FALSE` fixes it to zero. Per-factor constraint
vectors use the factor names, for example
`equal_loading = c(f1 = TRUE, f2 = FALSE)`.

## See also

[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) for the
general function this wraps;
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) for the
one factor wrapper;
[`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) and
[`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
for the discriminant and convergent validity summaries the
`output = "measurement"` table reports alongside omega.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
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
data(holzinger_swineford)

# Two factors, each named by its items.
cfa_2(holzinger_swineford,
      factor_1 = c("t6_paragraph_comprehension", "t7_sentence",
                   "t9_word_meaning"),
      factor_2 = c("t20_deduction", "t22_problem_reasoning",
                   "t23_series_completion"))
#> Measurement structure, per factor:
#>   f1: congeneric (no equality constraints)
#>   f2: congeneric (no equality constraints)
#> 
#>  syntax                                                   term          
#>  f1 =~ t6_paragraph_comprehension                         lambda_f1_1   
#>  f1 =~ t7_sentence                                        lambda_f1_2   
#>  f1 =~ t9_word_meaning                                    lambda_f1_3   
#>  f1 ~~ f1                                                 phi_f1        
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_f1_1      
#>  t7_sentence ~~ t7_sentence                               psi_f1_2      
#>  t9_word_meaning ~~ t9_word_meaning                       psi_f1_3      
#>  f2 =~ t20_deduction                                      lambda_f2_1   
#>  f2 =~ t22_problem_reasoning                              lambda_f2_2   
#>  f2 =~ t23_series_completion                              lambda_f2_3   
#>  f2 ~~ f2                                                 phi_f2        
#>  t20_deduction ~~ t20_deduction                           psi_f2_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_f2_2      
#>  t23_series_completion ~~ t23_series_completion           psi_f2_3      
#>  f1 ~~ f2                                                 phi_f1_f2     
#>                                                           loading_sum_f1
#>                                                           error_sum_f1  
#>                                                           omega_f1      
#>                                                           ave_f1        
#>                                                           H_f1          
#>                                                           loading_sum_f2
#>                                                           error_sum_f2  
#>                                                           omega_f2      
#>                                                           ave_f2        
#>                                                           H_f2          
#>                                                           chi_square    
#>                                                           df            
#>                                                           p_chi_square  
#>                                                           cfi           
#>                                                           tli           
#>                                                           nnfi          
#>                                                           rmsea         
#>                                                           rmsea_ci_lower
#>                                                           rmsea_ci_upper
#>                                                           rmsea_ci_level
#>                                                           srmr          
#>                                                           AIC           
#>                                                           BIC           
#>                                                           H0            
#>                                                           H1            
#>  estimate  se     z_value p_value  ci_lower ci_upper
#>  2.95      0.169  17.4    < 0.0001 2.61     3.28    
#>  4.39      0.249  17.6    < 0.0001 3.9      4.88    
#>  6.49      0.371  17.5    < 0.0001 5.76     7.22    
#>  1         0      <NA>    <NA>     1        1       
#>  3.47      0.418  8.32    < 0.0001 2.65     4.29    
#>  7.29      0.902  8.09    < 0.0001 5.53     9.06    
#>  16.5      2.01   8.23    < 0.0001 12.6     20.4    
#>  11.1      1.15   9.71    < 0.0001 8.9      13.4    
#>  6.74      0.524  12.9    < 0.0001 5.71     7.77    
#>  6.68      0.521  12.8    < 0.0001 5.66     7.7     
#>  1         0      <NA>    <NA>     1        1       
#>  248       23.5   10.6    < 0.0001 202      294     
#>  38.9      4.77   8.16    < 0.0001 29.6     48.2    
#>  38.6      4.71   8.2     < 0.0001 29.4     47.9    
#>  0.73      0.0428 17      < 0.0001 0.646    0.814   
#>  13.8      0.647  21.4    < 0.0001 12.6     15.1    
#>  27.3      2      13.6    < 0.0001 23.4     31.2    
#>  0.875     0.0136 64.3    < 0.0001 0.848    0.902   
#>  0.719     0.0228 31.6    < 0.0001 0.675    0.764   
#>  0.885     0.0115 77      < 0.0001 0.862    0.907   
#>  24.6      1.54   15.9    < 0.0001 21.5     27.6    
#>  326       23.8   13.7    < 0.0001 279      372     
#>  0.649     0.0354 18.3    < 0.0001 0.58     0.719   
#>  0.469     0.0332 14.1    < 0.0001 0.404    0.535   
#>  0.738     0.0273 27.1    < 0.0001 0.685    0.792   
#>  13.8      <NA>   <NA>    <NA>     <NA>     <NA>    
#>  8         <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.0880    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.993     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.987     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.987     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.0489    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.0915    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.0234    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  11755.495 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  11803.688 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5864.748 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# The measurement properties: omega, ave, and H per factor, the
# factor correlation, and the htmt ratio.
cfa_2(holzinger_swineford,
      factor_1 = c("t6_paragraph_comprehension", "t7_sentence",
                   "t9_word_meaning"),
      factor_2 = c("t20_deduction", "t22_problem_reasoning",
                   "t23_series_completion"),
      output = "measurement")
#> Measurement structure, per factor:
#>   f1: congeneric (no equality constraints)
#>   f2: congeneric (no equality constraints)
#> 
#>  syntax   term       estimate se     z_value p_value  ci_lower ci_upper
#>  f1 ~~ f2 phi_f1_f2  0.73     0.0428 17      < 0.0001 0.646    0.814   
#>           omega_f1   0.875    0.0136 64.3    < 0.0001 0.848    0.902   
#>           ave_f1     0.719    0.0228 31.6    < 0.0001 0.675    0.764   
#>           H_f1       0.885    0.0115 77      < 0.0001 0.862    0.907   
#>           omega_f2   0.649    0.0354 18.3    < 0.0001 0.58     0.719   
#>           ave_f2     0.469    0.0332 14.1    < 0.0001 0.404    0.535   
#>           H_f2       0.738    0.0273 27.1    < 0.0001 0.685    0.792   
#>  f1 ~~ f2 htmt_f1_f2 0.731    <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%
```
