# Two Factor Confirmatory Factor Analysis Model

Fits a two factor confirmatory factor analysis model to raw item data or
a sample covariance matrix. This is the two factor special case of
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md):
the function is a convenience wrapper that only requires the items of
each factor, builds the two factor specification, and forwards
everything else to
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).
The factors are named `f1` and `f2`, so the rows of the returned table
are `lambda_f1_1`, `lambda_f2_1`, `phi_f1_f2` (the factor correlation),
`omega_f1`, `omega_f2`, and so on, exactly as a two factor
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
call would report them (the `syntax` column names the item behind each
number). To name the factors substantively, or for three or more
factors, call
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
directly.

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
  missing = "fiml",
  se = NULL,
  conf_level = 0.95,
  output = c("verbose", "measurement", "summary", "standardized", "fit"),
  ...,
  wald_df = Inf,
  rmsea_conf_level = 0.9
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
  No item may appear in both factors. An item of either factor that is
  not a column of `data` or an item of `S` is refused with an error of
  class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

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
  `"MLR"` (robust maximum likelihood, raw data required: maximum
  likelihood estimates with the test statistic corrected for
  nonnormality, the scaled statistic of Satorra and Bentler, 1994, in
  the form Yuan and Bentler, 2000, gave it for incomplete data; its
  matching Huber-White sandwich standard errors are used when `se` is
  omitted), `"WLS"` (the asymptotic distribution free estimator of
  Browne, 1984; raw data and a large sample required), `"WLSMV"`
  (diagonally weighted least squares with mean- and variance-adjusted
  test statistic; Muthén, 1984; Muthén, du Toit, & Spisic, 1997; the
  standard choice for ordered categorical items, and what `ordered`
  switches to), or `"GLS"` (generalized least squares; Browne, 1974).
  With a robust estimator the reported fit indices are the robust
  versions.

- missing:

  Character; missing-data handling when raw data are supplied: `"fiml"`
  (default; full information maximum likelihood, which uses every case
  with an observed item; lavaan's name `"ml"` is accepted), `"listwise"`
  (listwise deletion), `"pairwise"` (passed to lavaan, for ordered
  items), `"two_stage"` (two-stage maximum likelihood with normal theory
  standard errors and test), or `"robust_two_stage"` (two-stage maximum
  likelihood with standard errors and test robust to nonnormality);
  lavaan's spellings `"two.stage"` and `"robust.two.stage"` are accepted
  for the last two. See *Details* for two-stage estimation. Ignored with
  `S`, except that a two-stage request with `S` is an error, since a
  covariance matrix has no incomplete cases. With `ordered` items,
  `"ml"`/`"fiml"` and the two-stage options are not available; use
  `"pairwise"` or `"listwise"`. FIML needs a maximum likelihood
  estimator of continuous items with normal theory or Huber-White
  standard errors. When the default meets ordered items, a least squares
  estimator (`"WLS"`, `"WLSMV"`, `"GLS"`), `se = "robust.sem"`, or a
  Satorra-Bentler kind of test passed through `...`, the items are
  fitted as they stand if no item value is missing (listwise deletion
  then deletes nothing), and otherwise the call stops and names the
  option to choose, rather than deleting cases unasked. With complete
  data FIML and listwise deletion give the same estimates and the same
  chi square, and for `estimator = "ML"` the same fit indices; the
  standard errors differ slightly, since FIML takes them from the
  observed information, as lavaan does and as Savalei (2010) recommends
  with incomplete data, and listwise deletion from the expected
  information.

- se:

  Standard error method passed to lavaan. The default, `NULL`, follows
  the estimator: model-based `"standard"` errors for ML, GLS and WLS;
  Huber-White sandwich `"robust.huber.white"` errors for MLR; and
  `"robust.sem"` errors for WLSMV. Explicit compatible choices are
  retained. The ML and MLR point estimates use the same likelihood;
  their standard errors and tests use different reference assumptions.
  This difference does not depend on whether optimization is iterative
  or an equivalent solution is computed analytically. `"robust.sem"`
  requires complete data and is refused with FIML; use
  `"robust.huber.white"` to retain incomplete respondents. No
  standard-error request silently changes the missing-data treatment.
  `"none"` requests point estimates only. Under two-stage estimation the
  matching corrected two-stage errors are used; `se` must be `NULL`,
  `"standard"`, or `"none"`. MLR and Huber-White errors require raw
  data. Robust SEM errors from covariance input require ML and a
  matching `NACOV` in `...`.

- conf_level:

  Confidence level for the parameter confidence intervals, including the
  delta method intervals for omega, AVE, and *H*. Defaults to 0.95. It
  does not reach the RMSEA interval, whose level is `rmsea_conf_level`.

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
      [`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html),
      whose defined rows describe the standardized items: `omega` there
      is the reliability of the unit-weighted sum of the standardized
      items, while `ave` and `H` do not depend on the metric of the
      items and equal the rows of the other tables. For a factor with
      `equal_loading` or `equal_error` these rows are computed from each
      item's own standardized loading and error variance (see *The
      standardized solution* under *Details* of
      [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)).
      Under `se = "none"` the `se`, `z`, `pvalue`, `ci.lower`, and
      `ci.upper` columns are `NA`.

  `"fit"`

  :   The raw lavaan fit object. Escape hatch for direct lavaan access,
      including likelihood ratio tests between two
      [`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
      fits via
      [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html).

- ...:

  Additional arguments forwarded to
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  and, through it, to
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html) (for example
  `ordered`, `equal_intercept`, or `M`).

- wald_df:

  The reference distribution of the Wald intervals and *p*-values of the
  parameter rows: `Inf` (the default), the normal distribution, which
  leaves lavaan's intervals and *p*-values as they are; a single
  positive number, the degrees of freedom of a *t* distribution for
  every parameter row, the intervals and *p*-values recomputed from the
  same estimates and standard errors; or `"missing_information"`, von
  Hippel's (2016) degrees of freedom \\\max\\3, N(1 - \gamma_j) - 1\\\\
  for each loading, from its fraction of missing information
  \\\gamma_j\\ (\\N - 1\\ with complete data), with every other row left
  on the normal reference. A number applies to the tables of
  `output = "verbose"`, `"measurement"`, and `"standardized"`, the
  keyword to that of `"verbose"` only (the `"measurement"` table has no
  loading), and either needs standard errors other than `"none"` and
  `"bootstrap"`; the keyword also needs the fits its calibration covers,
  maximum likelihood (`estimator = "ML"` or `"MLR"`) of continuous items
  by FIML or with complete data, with normal theory or Huber-White
  standard errors. Anything else is an error that says why. It follows
  `...`, so it must be named in full. See *Details*.

- rmsea_conf_level:

  Confidence level of the RMSEA interval of the fit rows, which the
  `rmsea_ci_level` row records. Defaults to 0.90, the conventional
  level, whose limits decide the tests of close and not close fit at .05
  (see
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md));
  `rmsea_conf_level = 0.95` gives a 95 percent interval. It sets the
  level in the table DMAR returns and in its `glance()`;
  `output = "summary"` returns the lavaan summary, which shows lavaan's
  90 percent interval whatever this argument is.

## Value

The value of the corresponding
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
call: a `data.frame` (classes `dmar_fa_confirmatory`, `dmar_tbl`) with
one row per parameter (`estimate`, `se`, `test_statistic`, `p_value`,
`ci_lower`, `ci_upper`) followed by the fit rows, or the alternative
shapes selected by `output` (see
[`?fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)).

## Details

Each factor is identified by fixing its variance to 1 and estimating
every loading; each item loads on exactly one factor (simple structure).
The factor correlation is estimated by default and
`correlated_factors = FALSE` fixes it to zero. Per-factor constraint
vectors use the factor names, for example
`equal_loading = c(f1 = TRUE, f2 = FALSE)`. A fit whose likelihood has
no maximum (an error variance that runs off below zero while the
likelihood is flat or still rising) stops with an error of class
`dmar_no_maximum` that prints the lavaan code of the fit, as
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
describes; there, `bounds = "pos.var"` (which `...` also passes on from
here) or a different model are the choices.

## References

Browne, M. W. (1974). Generalized least squares estimators in the
analysis of covariance structures. *South African Statistical Journal,
8*, 1–24.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*(1), 62–83.
[doi:10.1111/j.2044-8317.1984.tb00789.x](https://doi.org/10.1111/j.2044-8317.1984.tb00789.x)

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

Savalei, V. (2010). Expected versus observed information in SEM with
incomplete normal and nonnormal data. *Psychological Methods, 15*(4),
352–367. [doi:10.1037/a0020143](https://doi.org/10.1037/a0020143)

von Hippel, P. T. (2016). New confidence intervals and bias comparisons
show that maximum likelihood can beat multiple imputation in small
samples. *Structural Equation Modeling, 23*(3), 422–437.
[doi:10.1080/10705511.2015.1047931](https://doi.org/10.1080/10705511.2015.1047931)

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
*Sociological Methodology, 30*, 165–200.

## See also

[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
for the general function this wraps;
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
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
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
#>  t6_paragraph_comprehension ~1                            nu_f1_1       
#>  t7_sentence ~1                                           nu_f1_2       
#>  t9_word_meaning ~1                                       nu_f1_3       
#>  f2 =~ t20_deduction                                      lambda_f2_1   
#>  f2 =~ t22_problem_reasoning                              lambda_f2_2   
#>  f2 =~ t23_series_completion                              lambda_f2_3   
#>  f2 ~~ f2                                                 phi_f2        
#>  t20_deduction ~~ t20_deduction                           psi_f2_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_f2_2      
#>  t23_series_completion ~~ t23_series_completion           psi_f2_3      
#>  t20_deduction ~1                                         nu_f2_1       
#>  t22_problem_reasoning ~1                                 nu_f2_2       
#>  t23_series_completion ~1                                 nu_f2_3       
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
#>  estimate  se     test_statistic p_value  ci_lower ci_upper
#>  2.947     0.169  17.387         < 0.0001 2.614    3.279   
#>  4.389     0.249  17.608         < 0.0001 3.900    4.877   
#>  6.489     0.373  17.417         < 0.0001 5.759    7.220   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  3.473     0.419  8.280          < 0.0001 2.651    4.296   
#>  7.295     0.899  8.117          < 0.0001 5.533    9.056   
#>  16.509    2.043  8.082          < 0.0001 12.506   20.512  
#>  9.183     0.201  45.694         < 0.0001 8.789    9.577   
#>  17.362    0.297  58.452         < 0.0001 16.780   17.944  
#>  15.299    0.441  34.667         < 0.0001 14.434   16.164  
#>  11.149    1.152  9.678          < 0.0001 8.891    13.407  
#>  6.738     0.525  12.844         < 0.0001 5.710    7.767   
#>  6.681     0.524  12.751         < 0.0001 5.654    7.708   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  248.254   23.592 10.523         < 0.0001 202.015  294.493 
#>  38.903    4.782  8.135          < 0.0001 29.529   48.276  
#>  38.626    4.773  8.093          < 0.0001 29.272   47.981  
#>  26.890    1.113  24.170         < 0.0001 24.710   29.071  
#>  26.239    0.529  49.579         < 0.0001 25.202   27.276  
#>  18.136    0.526  34.483         < 0.0001 17.105   19.167  
#>  0.730     0.043  16.885         < 0.0001 0.645    0.815   
#>  13.825    0.647  21.376         < 0.0001 12.557   15.092  
#>  27.277    1.998  13.650         < 0.0001 23.360   31.194  
#>  0.875     0.014  64.416         < 0.0001 0.846    0.899   
#>  0.719     0.023  31.557         < 0.0001 0.673    0.762   
#>  0.885     0.012  76.912         < 0.0001 0.860    0.906   
#>  24.568    1.548  15.867         < 0.0001 21.534   27.603  
#>  325.783   23.890 13.637         < 0.0001 278.959  372.606 
#>  0.649     0.036  18.221         < 0.0001 0.577    0.716   
#>  0.469     0.033  14.095         < 0.0001 0.405    0.535   
#>  0.738     0.027  27.128         < 0.0001 0.682    0.788   
#>  13.769    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  8         <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.0880    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.993     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.987     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.987     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.049     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.000     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.092     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.900     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.023     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11767.495 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11837.931 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5864.748 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>           <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)

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
#>  syntax   term       estimate se    test_statistic p_value  ci_lower ci_upper
#>  f1 ~~ f2 phi_f1_f2  0.730    0.043 16.885         < 0.0001 0.645    0.815   
#>           omega_f1   0.875    0.014 64.416         < 0.0001 0.846    0.899   
#>           ave_f1     0.719    0.023 31.557         < 0.0001 0.673    0.762   
#>           H_f1       0.885    0.012 76.912         < 0.0001 0.860    0.906   
#>           omega_f2   0.649    0.036 18.221         < 0.0001 0.577    0.716   
#>           ave_f2     0.469    0.033 14.095         < 0.0001 0.405    0.535   
#>           H_f2       0.738    0.027 27.128         < 0.0001 0.682    0.788   
#>  f1 ~~ f2 htmt_f1_f2 0.731    <NA>  <NA>           <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%
```
