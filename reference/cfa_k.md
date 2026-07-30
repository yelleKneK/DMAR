# Multiple-Factor Confirmatory Factor Analysis Model

Fits a confirmatory factor analysis model with one or more factors,
where each factor is specified by naming its indicator variables and the
measurement structure is specified by describing what is constrained
(`equal_loading`, `equal_intercept`, `equal_error`) rather than by more
technical terms. The function then reports which classical measurement
structure the description implies (congeneric, essentially
tau-equivalent, tau-equivalent, essentially parallel, or parallel), the
parameter estimates with confidence intervals, fit information, and, per
factor, coefficient omega, the average variance extracted (AVE), and
coefficient *H*, each with a delta method standard error and confidence
interval computed by lavaan from defined parameters (no additional
packages are involved).

## Usage

``` r
cfa_k(
  data = NULL,
  factors,
  S = NULL,
  N = NULL,
  M = NULL,
  equal_loading = FALSE,
  equal_intercept = FALSE,
  equal_error = FALSE,
  correlated_factors = TRUE,
  meanstructure = NULL,
  estimator = "ML",
  missing = "listwise",
  ordered = NULL,
  se = "standard",
  conf_level = 0.95,
  output = c("verbose", "measurement", "summary", "standardized", "fit"),
  ...
)
```

## Arguments

- data:

  A raw data matrix or data frame, rows are respondents and columns
  include the items named in `factors`. Supply exactly one of `data` or
  `S`.

- factors:

  Named list. Each element names a factor and gives the character vector
  of its indicator columns (two or more per factor; three or more when
  only one factor is specified). Each item loads on exactly one factor
  (simple structure).

- S:

  A symmetric covariance matrix of the items, with dimnames naming the
  items; `N` is then required. Supply exactly one of `data` or `S`.

- N:

  Total sample size. Required with `S`; ignored (inferred from the rows)
  with `data`.

- M:

  Optional named numeric vector of item means, used with a covariance
  matrix to model the mean structure (required there when
  `equal_intercept` is used). Ignored for raw data.

- equal_loading:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the loadings within a factor to a single value (across
  factors nothing is equated). Defaults to `FALSE`.

- equal_intercept:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the item intercepts within a factor to a single value; this
  requires the mean structure (raw data, or `M` with a covariance
  matrix). Defaults to `FALSE`.

- equal_error:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the error variances within a factor to a single value.
  Defaults to `FALSE`.

- correlated_factors:

  Logical. If `TRUE` (default) the factors covary freely; because each
  factor variance is fixed to 1, the `phi` terms for factor pairs are
  the latent correlations. If `FALSE` the factor covariances are fixed
  to zero.

- meanstructure:

  Logical or `NULL`. `NULL` (default) models the mean structure exactly
  when it is needed: when any `equal_intercept` is `TRUE` or when `M` is
  supplied. Set `TRUE` to model intercepts regardless (raw data or `M`
  required), or `FALSE` to suppress them (an error if `equal_intercept`
  is used).

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

- ordered:

  Ordered-categorical items: `NULL` (none, the default), `TRUE` (every
  item), or a character vector of item names. Requires raw data. Each
  factor must be all ordered or all continuous. Declaring ordered items
  switches the estimator to `"WLSMV"` (with a message, unless a
  categorical estimator was requested), fits thresholds in place of
  intercepts, and reports each ordered factor's omega on the categorical
  sum score metric via the Green and Yang (2009) computation; see
  *Details*.

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
      including likelihood ratio tests between two `cfa_k()` fits via
      [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html).

- ...:

  Additional arguments forwarded to
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html) (e.g., `group`,
  `cluster`, `bootstrap`).

## Value

For `output = "verbose"` (default) and `output = "measurement"`, a
`data.frame` (classes `dmar_cfa_k`, `dmar_tbl`) with columns `syntax`,
`term`, `estimate`, `se`, `z_value`, `p_value`, `ci_lower`, `ci_upper`.
The `"model"` attribute is a named character vector giving, per factor,
the implied classical structure; the printed header displays it. For
`output = "summary"`, the lavaan summary object; for `"standardized"`,
the standardized solution; for `"fit"`, the lavaan fit object.

## Details

With `ordered` items, the model is fit by WLSMV to polychoric
correlations with thresholds. A sum score of ordered items lives on the
metric of the observed categories, not on the latent response metric of
the polychoric loadings, so for a factor whose items are ordered the
reported omega is the Green and Yang (2009) categorical sum score omega
computed from the same fit; the substitution is announced in a message
and recorded in the `omega_metric` attribute, and the delta method
interval columns are `NA` for those rows because the delta method
interval describes the latent response metric. The categorical sum score
omega is the coefficient Kelley and Pornprasertmanit (2016) call
categorical omega. AVE and coefficient *H* concern the latent response
correlations themselves and are reported unchanged on that metric. For a
bootstrap confidence interval on a categorical omega, use
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
on the factor's items.

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) and
[`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) are
convenience wrappers around this function for the one and two factor
cases: [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
takes a vector of items and fits one factor over them, and
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) takes
the items of each of two factors. Both forward every argument here, so
their results are this function's results, with the factors named `f1`
(and `f2`).

**Describing the model instead of naming it.** The classical measurement
structures are nested patterns of within-factor equality constraints
(Lord & Novick, 1968; Graham, 2006):

- congeneric:

  loadings, intercepts, and error variances all free.

- essentially tau-equivalent:

  equal loadings; intercepts and error variances free.

- tau-equivalent:

  equal loadings and equal intercepts; error variances free.

- essentially parallel:

  equal loadings and equal error variances; intercepts free.

- parallel:

  equal loadings, equal intercepts, and equal error variances.

The caller states the constraints; the function reports the implied
name, per factor, in the printed header and in the `"model"` attribute
of the returned table. The distinction between tau-equivalent and
essentially tau-equivalent (and between parallel and essentially
parallel) lives entirely in the mean structure: the covariance structure
of the two members of each pair is identical, so without intercepts in
the model only the "essentially" form can be claimed. That is why
`equal_intercept` requires raw data or `M`: covariances alone cannot
speak to it. A constraint pattern outside the classical list (for
example equal error variances with free loadings) is fit as requested
and labeled descriptively, since it has no conventional name.

Identification fixes each factor variance to 1 (and each factor mean to
0 when the mean structure is modeled), so all loadings are estimated and
within-factor equality constraints are meaningful. Two `cfa_k()` fits
that differ only in descriptor settings are nested, so `output = "fit"`
feeds
[`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html)
directly (with `estimator = "MLR"`, lavaan applies the scaled difference
test).

**Measurement properties.** For factor \\f\\ with unstandardized
loadings \\\lambda_j\\ and error variances \\\psi_j\\ (factor variance
1): coefficient omega \\= (\sum_j \lambda_j)^2 / ((\sum_j \lambda_j)^2 +
\sum_j \psi_j)\\ (McDonald, 1999), the reliability of the unit-weighted
composite; the average variance extracted \\= J^{-1} \sum_j \lambda_j^2
/ (\lambda_j^2 + \psi_j)\\ (Fornell & Larcker, 1981), the mean
proportion of item variance the factor accounts for; and coefficient \\H
= (1 + (\sum_j \lambda_j^2/\psi_j)^{-1})^{-1}\\ (Hancock & Mueller,
2001), the reliability of the optimally weighted composite, which no
single item can drag below its value for any subset. All three are
computed as lavaan defined parameters, so each carries a delta method
standard error and a `conf_level` confidence interval in the same table
as the model parameters.

**Discriminant validity.** Three complementary readings come from
`output = "measurement"`: (a) the latent correlation `phi` for a factor
pair, with a confidence interval whose upper limit near 1 means the data
cannot distinguish the two factors; (b) the Fornell and Larcker (1981)
comparison, which asks whether each factor's `ave` exceeds the squared
`phi` of its pairs (the factor should share more variance with its own
items than with the other factor); and (c) for raw data, the model-free
`htmt` ratio (Henseler, Ringle, & Sarstedt, 2015). The rows report the
numbers and their uncertainty; the judgment is the researcher's.

**Confidence interval conventions.** Parameter rows (including omega,
`ave`, *H*, and `phi`) use `conf_level`. The RMSEA interval follows its
own convention: the `rmsea_ci_level` row records the level actually used
(0.90, the conventional level for RMSEA, as in lavaan), and
`rmsea_ci_lower` / `rmsea_ci_upper` are that interval.

Common row names under `term`: `lambda_<factor>_<j>` (loadings;
`lambda_<factor>` when equated), `psi_<factor>_<j>` (error variances;
`psi_<factor>` when equated), `nu_<factor>_<j>` (intercepts, when the
mean structure is modeled; `nu_<factor>` when equated), `phi_<factor>`
(factor variance, fixed to 1), `phi_<factor1>_<factor2>` (latent
correlation), the per-factor defined parameters (`loading_sum_<factor>`,
`error_sum_<factor>`, `omega_<factor>`, `ave_<factor>`, `H_<factor>`),
and the fit rows `chi_square`, `df`, `p_chi_square`, `cfi`, `tli`,
`nnfi`, `rmsea`, `rmsea_ci_lower`, `rmsea_ci_upper`, `rmsea_ci_level`,
`srmr`, `AIC`, `BIC`, `H0`, `H1`.

## References

Browne, M. W. (1974). Generalized least squares estimators in the
analysis of covariance structures. *South African Statistical Journal,
8*, 1–24.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*, 62–83.

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Muthén, B., du Toit, S. H. C., & Spisic, D. (1997). *Robust inference
using weighted least squares and quadratic estimating equations in
latent variable modeling with categorical and continuous outcomes*.
Unpublished technical report.

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
standard errors in covariance structure analysis. In A. von Eye & C. C.
Clogg (Eds.), *Latent variables analysis: Applications for developmental
research* (pp. 399–419). Thousand Oaks, CA: Sage.

Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation
models with unobservable variables and measurement error. *Journal of
Marketing Research, 18*(1), 39–50.

Graham, J. M. (2006). Congeneric and (essentially) tau-equivalent
estimates of score reliability: What they are and how to use them.
*Educational and Psychological Measurement, 66*(6), 930–944.
[doi:10.1177/0013164406288165](https://doi.org/10.1177/0013164406288165)

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*(1), 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for
assessing discriminant validity in variance-based structural equation
modeling. *Journal of the Academy of Marketing Science, 43*(1), 115–135.
[doi:10.1007/s11747-014-0403-8](https://doi.org/10.1007/s11747-014-0403-8)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Lord, F. M., & Novick, M. R. (1968). *Statistical theories of mental
test scores*. Addison-Wesley.

McDonald, R. P. (1999). *Test theory: A unified treatment*. Erlbaum.

## See also

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) and
[`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) for the
one and two factor convenience wrappers;
[`plot_cfa_k`](https://yelleknek.github.io/DMAR/reference/plot_cfa_k.md)
to display the estimates and the equality question visually;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
and [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) for the
measurement properties as standalone functions;
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
for the across-group analog of these within-factor constraints;
[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html),
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html).

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
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
hs_factors <- list(
  verbal    = c("t6_paragraph_comprehension", "t7_sentence",
                "t9_word_meaning"),
  deduction = c("t20_deduction", "t22_problem_reasoning",
                "t23_series_completion"))

# Congeneric measurement model for both factors (the default:
# nothing is constrained, and the header names the structure).
cfa_k(holzinger_swineford, hs_factors)
#> Measurement structure, per factor:
#>   verbal: congeneric (no equality constraints)
#>   deduction: congeneric (no equality constraints)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal_1      
#>  verbal =~ t7_sentence                                    lambda_verbal_2      
#>  verbal =~ t9_word_meaning                                lambda_verbal_3      
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal_1         
#>  t7_sentence ~~ t7_sentence                               psi_verbal_2         
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal_3         
#>  deduction =~ t20_deduction                               lambda_deduction_1   
#>  deduction =~ t22_problem_reasoning                       lambda_deduction_2   
#>  deduction =~ t23_series_completion                       lambda_deduction_3   
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction_2      
#>  t23_series_completion ~~ t23_series_completion           psi_deduction_3      
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
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

# Equal loadings within every factor. Because only the covariance
# structure identifies this constraint, the implied structure is
# essentially tau-equivalent, and the header says so.
cfa_k(holzinger_swineford, hs_factors, equal_loading = TRUE)
#> Measurement structure, per factor:
#>   verbal: essentially tau-equivalent (equal loadings)
#>   deduction: essentially tau-equivalent (equal loadings)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal        
#>  verbal =~ t7_sentence                                    lambda_verbal        
#>  verbal =~ t9_word_meaning                                lambda_verbal        
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal_1         
#>  t7_sentence ~~ t7_sentence                               psi_verbal_2         
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal_3         
#>  deduction =~ t20_deduction                               lambda_deduction     
#>  deduction =~ t22_problem_reasoning                       lambda_deduction     
#>  deduction =~ t23_series_completion                       lambda_deduction     
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction_2      
#>  t23_series_completion ~~ t23_series_completion           psi_deduction_3      
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
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
#>  3.46      0.163  21.2    < 0.0001 3.14     3.78    
#>  3.46      0.163  21.2    < 0.0001 3.14     3.78    
#>  3.46      0.163  21.2    < 0.0001 3.14     3.78    
#>  1         0      <NA>    <NA>     1        1       
#>  1.22      0.508  2.41    0.0159   0.229    2.22    
#>  10.8      1.03   10.5    < 0.0001 8.81     12.8    
#>  31.3      2.65   11.8    < 0.0001 26.1     36.5    
#>  7.02      0.399  17.6    < 0.0001 6.24     7.8     
#>  7.02      0.399  17.6    < 0.0001 6.24     7.8     
#>  7.02      0.399  17.6    < 0.0001 6.24     7.8     
#>  1         0      <NA>    <NA>     1        1       
#>  281       23.8   11.8    < 0.0001 234      328     
#>  38.4      4.53   8.46    < 0.0001 29.5     47.2    
#>  36.2      4.41   8.22    < 0.0001 27.6     44.9    
#>  0.667     0.0473 14.1    < 0.0001 0.574    0.759   
#>  10.4      0.489  21.2    < 0.0001 9.42     11.3    
#>  43.3      2.8    15.5    < 0.0001 37.9     48.8    
#>  0.713     0.0233 30.7    < 0.0001 0.667    0.758   
#>  0.569     0.0231 24.7    < 0.0001 0.524    0.615   
#>  0.918     0.0299 30.8    < 0.0001 0.86     0.977   
#>  21.1      1.2    17.6    < 0.0001 18.7     23.4    
#>  356       24.5   14.5    < 0.0001 308      404     
#>  0.555     0.0338 16.4    < 0.0001 0.489    0.621   
#>  0.429     0.0295 14.6    < 0.0001 0.372    0.487   
#>  0.738     0.029  25.5    < 0.0001 0.682    0.795   
#>  149       <NA>   <NA>    <NA>     <NA>     <NA>    
#>  12        <NA>   <NA>    <NA>     <NA>     <NA>    
#>  < 0.0001  <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.83      <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.787     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.787     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.195     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.168     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.223     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.179     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  11882.876 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  11916.240 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5932.438 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# Descriptors can differ by factor.
cfa_k(holzinger_swineford, hs_factors,
      equal_loading = c(verbal = TRUE, deduction = FALSE))
#> Measurement structure, per factor:
#>   verbal: essentially tau-equivalent (equal loadings)
#>   deduction: congeneric (no equality constraints)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal        
#>  verbal =~ t7_sentence                                    lambda_verbal        
#>  verbal =~ t9_word_meaning                                lambda_verbal        
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal_1         
#>  t7_sentence ~~ t7_sentence                               psi_verbal_2         
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal_3         
#>  deduction =~ t20_deduction                               lambda_deduction_1   
#>  deduction =~ t22_problem_reasoning                       lambda_deduction_2   
#>  deduction =~ t23_series_completion                       lambda_deduction_3   
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction_2      
#>  t23_series_completion ~~ t23_series_completion           psi_deduction_3      
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
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
#>  3.46      0.163  21.2    < 0.0001 3.14     3.78    
#>  3.46      0.163  21.2    < 0.0001 3.14     3.78    
#>  3.46      0.163  21.2    < 0.0001 3.14     3.78    
#>  1         0      <NA>    <NA>     1        1       
#>  1.22      0.508  2.4     0.0163   0.225    2.21    
#>  10.8      1.03   10.5    < 0.0001 8.82     12.8    
#>  31.3      2.65   11.8    < 0.0001 26.1     36.5    
#>  11.2      1.16   9.67    < 0.0001 8.91     13.4    
#>  6.57      0.533  12.3    < 0.0001 5.53     7.61    
#>  6.86      0.527  13      < 0.0001 5.82     7.89    
#>  1         0      <NA>    <NA>     1        1       
#>  248       23.7   10.5    < 0.0001 201      294     
#>  41.1      4.93   8.34    < 0.0001 31.5     50.8    
#>  36.2      4.86   7.46    < 0.0001 26.7     45.8    
#>  0.667     0.047  14.2    < 0.0001 0.575    0.759   
#>  10.4      0.489  21.2    < 0.0001 9.42     11.3    
#>  43.4      2.8    15.5    < 0.0001 37.9     48.8    
#>  0.713     0.0233 30.7    < 0.0001 0.667    0.758   
#>  0.57      0.0231 24.7    < 0.0001 0.524    0.615   
#>  0.919     0.0299 30.8    < 0.0001 0.86     0.977   
#>  24.6      1.55   15.9    < 0.0001 21.6     27.6    
#>  325       23.8   13.6    < 0.0001 278      372     
#>  0.651     0.0355 18.3    < 0.0001 0.581    0.72    
#>  0.471     0.0332 14.2    < 0.0001 0.406    0.536   
#>  0.74      0.0274 27      < 0.0001 0.687    0.794   
#>  134       <NA>   <NA>    <NA>     <NA>     <NA>    
#>  10        <NA>   <NA>    <NA>     <NA>     <NA>    
#>  < 0.0001  <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.846     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.769     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.769     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.203     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.173     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.234     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.164     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  11871.888 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  11912.666 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5924.944 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# Equal loadings and intercepts (tau-equivalent), then also equal
# error variances (parallel). The mean structure is added because
# equal_intercept asks about it.
cfa_k(holzinger_swineford, hs_factors, equal_loading = TRUE,
      equal_intercept = TRUE)
#> Measurement structure, per factor:
#>   verbal: tau-equivalent (equal loadings and intercepts)
#>   deduction: tau-equivalent (equal loadings and intercepts)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal        
#>  verbal =~ t7_sentence                                    lambda_verbal        
#>  verbal =~ t9_word_meaning                                lambda_verbal        
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal_1         
#>  t7_sentence ~~ t7_sentence                               psi_verbal_2         
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal_3         
#>  t6_paragraph_comprehension ~1                            nu_verbal            
#>  t7_sentence ~1                                           nu_verbal            
#>  t9_word_meaning ~1                                       nu_verbal            
#>  deduction =~ t20_deduction                               lambda_deduction     
#>  deduction =~ t22_problem_reasoning                       lambda_deduction     
#>  deduction =~ t23_series_completion                       lambda_deduction     
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction_2      
#>  t23_series_completion ~~ t23_series_completion           psi_deduction_3      
#>  t20_deduction ~1                                         nu_deduction         
#>  t22_problem_reasoning ~1                                 nu_deduction         
#>  t23_series_completion ~1                                 nu_deduction         
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
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
#>  4.37      0.258  17      < 0.0001 3.86     4.87    
#>  4.37      0.258  17      < 0.0001 3.86     4.87    
#>  4.37      0.258  17      < 0.0001 3.86     4.87    
#>  1         0      <NA>    <NA>     1        1       
#>  49.6      4.44   11.2    < 0.0001 40.9     58.3    
#>  13.9      1.75   7.94    < 0.0001 10.5     17.4    
#>  25.5      2.54   10.1    < 0.0001 20.6     30.5    
#>  15.5      0.298  52      < 0.0001 14.9     16.1    
#>  15.5      0.298  52      < 0.0001 14.9     16.1    
#>  15.5      0.298  52      < 0.0001 14.9     16.1    
#>  6.39      0.473  13.5    < 0.0001 5.47     7.32    
#>  6.39      0.473  13.5    < 0.0001 5.47     7.32    
#>  6.39      0.473  13.5    < 0.0001 5.47     7.32    
#>  1         0      <NA>    <NA>     1        1       
#>  304       25.7   11.8    < 0.0001 253      354     
#>  52.4      5.96   8.79    < 0.0001 40.7     64.1    
#>  79.2      7.76   10.2    < 0.0001 64       94.4    
#>  23.4      0.48   48.7    < 0.0001 22.4     24.3    
#>  23.4      0.48   48.7    < 0.0001 22.4     24.3    
#>  23.4      0.48   48.7    < 0.0001 22.4     24.3    
#>  0.901     0.0577 15.6    < 0.0001 0.788    1.01    
#>  13.1      0.773  17      < 0.0001 11.6     14.6    
#>  89        5.33   16.7    < 0.0001 78.6     99.5    
#>  0.659     0.0315 20.9    < 0.0001 0.597    0.72    
#>  0.428     0.0333 12.8    < 0.0001 0.363    0.493   
#>  0.715     0.0316 22.6    < 0.0001 0.653    0.777   
#>  19.2      1.42   13.5    < 0.0001 16.4     22      
#>  435       27.8   15.6    < 0.0001 381      490     
#>  0.458     0.0423 10.8    < 0.0001 0.375    0.541   
#>  0.299     0.0349 8.58    < 0.0001 0.231    0.368   
#>  0.589     0.0454 13      < 0.0001 0.5      0.678   
#>  1100      <NA>   <NA>    <NA>     <NA>     <NA>    
#>  16        <NA>   <NA>    <NA>     <NA>     <NA>    
#>  < 0.0001  <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -0.259    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -0.259    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.474     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.451     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.498     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.989     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  12837.111 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  12877.889 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -6407.556 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%
cfa_k(holzinger_swineford, hs_factors, equal_loading = TRUE,
      equal_intercept = TRUE, equal_error = TRUE)
#> Measurement structure, per factor:
#>   verbal: parallel (equal loadings, intercepts, and error variances)
#>   deduction: parallel (equal loadings, intercepts, and error variances)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal        
#>  verbal =~ t7_sentence                                    lambda_verbal        
#>  verbal =~ t9_word_meaning                                lambda_verbal        
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal           
#>  t7_sentence ~~ t7_sentence                               psi_verbal           
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal           
#>  t6_paragraph_comprehension ~1                            nu_verbal            
#>  t7_sentence ~1                                           nu_verbal            
#>  t9_word_meaning ~1                                       nu_verbal            
#>  deduction =~ t20_deduction                               lambda_deduction     
#>  deduction =~ t22_problem_reasoning                       lambda_deduction     
#>  deduction =~ t23_series_completion                       lambda_deduction     
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction        
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction        
#>  t23_series_completion ~~ t23_series_completion           psi_deduction        
#>  t20_deduction ~1                                         nu_deduction         
#>  t22_problem_reasoning ~1                                 nu_deduction         
#>  t23_series_completion ~1                                 nu_deduction         
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
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
#>  3.75      0.274  13.7    < 0.0001 3.22     4.29    
#>  3.75      0.274  13.7    < 0.0001 3.22     4.29    
#>  3.75      0.274  13.7    < 0.0001 3.22     4.29    
#>  1         0      <NA>    <NA>     1        1       
#>  30.4      1.75   17.3    < 0.0001 27       33.8    
#>  30.4      1.75   17.3    < 0.0001 27       33.8    
#>  30.4      1.75   17.3    < 0.0001 27       33.8    
#>  13.9      0.284  49.2    < 0.0001 13.4     14.5    
#>  13.9      0.284  49.2    < 0.0001 13.4     14.5    
#>  13.9      0.284  49.2    < 0.0001 13.4     14.5    
#>  7.56      0.585  12.9    < 0.0001 6.41     8.7     
#>  7.56      0.585  12.9    < 0.0001 6.41     8.7     
#>  7.56      0.585  12.9    < 0.0001 6.41     8.7     
#>  1         0      <NA>    <NA>     1        1       
#>  139       8      17.3    < 0.0001 123      154     
#>  139       8      17.3    < 0.0001 123      154     
#>  139       8      17.3    < 0.0001 123      154     
#>  23.8      0.586  40.5    < 0.0001 22.6     24.9    
#>  23.8      0.586  40.5    < 0.0001 22.6     24.9    
#>  23.8      0.586  40.5    < 0.0001 22.6     24.9    
#>  0.979     0.0704 13.9    < 0.0001 0.841    1.12    
#>  11.3      0.823  13.7    < 0.0001 9.65     12.9    
#>  91.2      5.26   17.3    < 0.0001 80.9     102     
#>  0.582     0.0418 13.9    < 0.0001 0.5      0.664   
#>  0.317     0.0371 8.53    < 0.0001 0.244    0.39    
#>  0.582     0.0418 13.9    < 0.0001 0.5      0.664   
#>  22.7      1.75   12.9    < 0.0001 19.2     26.1    
#>  416       24     17.3    < 0.0001 369      463     
#>  0.553     0.0447 12.4    < 0.0001 0.465    0.64    
#>  0.292     0.0373 7.81    < 0.0001 0.219    0.365   
#>  0.553     0.0447 12.4    < 0.0001 0.465    0.64    
#>  1340      <NA>   <NA>    <NA>     <NA>     <NA>    
#>  20        <NA>   <NA>    <NA>     <NA>     <NA>    
#>  < 0.0001  <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0         <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -0.228    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -0.228    <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.468     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.447     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.49      <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.9       <NA>   <NA>    <NA>     <NA>     <NA>    
#>  0.743     <NA>   <NA>    <NA>     <NA>     <NA>    
#>  13070.624 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  13096.574 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -6528.312 <NA>   <NA>    <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>    <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# Measurement properties: omega, ave, and H per factor (each with a
# delta method standard error and confidence interval), the latent
# correlations, and the htmt ratios.
cfa_k(holzinger_swineford, hs_factors, output = "measurement")
#> Measurement structure, per factor:
#>   verbal: congeneric (no equality constraints)
#>   deduction: congeneric (no equality constraints)
#> 
#>  syntax              term                  estimate se     z_value p_value 
#>  verbal ~~ deduction phi_verbal_deduction  0.73     0.0428 17      < 0.0001
#>                      omega_verbal          0.875    0.0136 64.3    < 0.0001
#>                      ave_verbal            0.719    0.0228 31.6    < 0.0001
#>                      H_verbal              0.885    0.0115 77      < 0.0001
#>                      omega_deduction       0.649    0.0354 18.3    < 0.0001
#>                      ave_deduction         0.469    0.0332 14.1    < 0.0001
#>                      H_deduction           0.738    0.0273 27.1    < 0.0001
#>  verbal ~~ deduction htmt_verbal_deduction 0.731    <NA>   <NA>    <NA>    
#>  ci_lower ci_upper
#>  0.646    0.814   
#>  0.848    0.902   
#>  0.675    0.764   
#>  0.862    0.907   
#>  0.58     0.719   
#>  0.404    0.535   
#>  0.685    0.792   
#>  <NA>     <NA>    
#> 
#> Confidence level: 95%

# \donttest{
# Ordered-categorical items: the model is fit by WLSMV to polychoric
# correlations, and each ordered factor's omega is reported on the
# categorical sum score metric (Green & Yang, 2009).
set.seed(113)
eta <- rnorm(200)
lat <- sweep(matrix(rep(eta, 6), 200, 6), 2,
             seq(0.5, 0.8, length.out = 6), `*`) +
  matrix(rnorm(200 * 6), 200, 6) %*%
  diag(sqrt(1 - seq(0.5, 0.8, length.out = 6)^2))
likert <- as.data.frame(apply(lat, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1, 0, 1, Inf)))))
names(likert) <- paste0("item_", 1:6)
cfa_k(likert,
      list(scale_a = paste0("item_", 1:3),
           scale_b = paste0("item_", 4:6)),
      ordered = TRUE, output = "measurement")
#> Ordered items declared: switching to estimator = "WLSMV" with robust standard errors.
#> Warning: lavaan->lav_object_post_check():  
#>    covariance matrix of latent variables is not positive definite ; use 
#>    lavInspect(fit, "cov.lv") to investigate.
#> omega_scale_a, omega_scale_b reported on the categorical sum score metric (Green & Yang, 2009).
#> Measurement structure, per factor:
#>   scale_a: congeneric (no equality constraints)
#>   scale_b: congeneric (no equality constraints)
#> 
#>  syntax             term                 estimate se     z_value p_value 
#>  scale_a ~~ scale_b phi_scale_a_scale_b  1.07     0.0882 12.1    < 0.0001
#>                     omega_scale_a        0.437    <NA>   <NA>    <NA>    
#>                     ave_scale_a          0.235    0.0537 4.38    < 0.0001
#>                     H_scale_a            0.482    0.0753 6.4     < 0.0001
#>                     omega_scale_b        0.774    <NA>   <NA>    <NA>    
#>                     ave_scale_b          0.601    0.0389 15.5    < 0.0001
#>                     H_scale_b            0.827    0.0252 32.8    < 0.0001
#>  scale_a ~~ scale_b htmt_scale_a_scale_b 1.06     <NA>   <NA>    <NA>    
#>  ci_lower ci_upper
#>  0.893    1.24    
#>  <NA>     <NA>    
#>  0.13     0.34    
#>  0.335    0.63    
#>  <NA>     <NA>    
#>  0.525    0.677   
#>  0.778    0.877   
#>  <NA>     <NA>    
#> 
#> Confidence level: 95%

# Does the equal-loadings description hold? Compare nested fits.
fit_free  <- cfa_k(holzinger_swineford, hs_factors, output = "fit")
fit_equal <- cfa_k(holzinger_swineford, hs_factors,
                   equal_loading = TRUE, output = "fit")
lavaan::lavTestLRT(fit_free, fit_equal)
#> 
#> Chi-Squared Difference Test
#> 
#>           Df   AIC   BIC   Chisq Chisq diff   RMSEA Df diff Pr(>Chisq)    
#> fit_free   8 11756 11804  13.769                                          
#> fit_equal 12 11883 11916 149.150     135.38 0.33033       4  < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
