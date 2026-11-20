# One Factor Confirmatory Factor Analysis Model

Fits a single factor (congeneric by default) confirmatory factor
analysis model to raw item data or a sample covariance matrix. This is
the one factor special case of
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md):
the function is a convenience wrapper that only requires the data (and,
when the data hold more than the items, a vector of item names), builds
the one factor specification, and forwards everything else to
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).
The factor is named `f1`, so the rows of the returned table are
`lambda_f1_1`, `lambda_f1_2`, ..., `psi_f1_1`, ..., and `omega_f1`,
exactly as a one factor
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
call would report them (the `syntax` column names the item behind each
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
  include the items. A matrix without column names is given the names
  `y1`, `y2`, ... Supply exactly one of `data` or `S`.

- items:

  Character vector naming the items of the factor (three or more; two
  are accepted with `equal_loading = TRUE`, the just identified
  tau-equivalent case). The default `NULL` uses every column of `data`
  (or of `S`), so a data set that holds only the items needs no `items`
  at all. An item that is not a column of `data` or an item of `S` is
  refused with an error of class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

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
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

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
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
whose `ordered` argument fits WLSMV with the theta parameterization and
reports categorical omega on the sum score metric, or
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md).

A fit whose likelihood has no maximum (an error variance that runs off
below zero while the likelihood is flat or still rising, as it can on
three items whose correlations have a negative product) stops with an
error of class `dmar_no_maximum` that prints the lavaan code of the fit,
as
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
for the general function this wraps (factor analysis with any number of
factors, intercept constraints, ordered categorical items, and the
measurement output);
[`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) for the
two factor wrapper;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for coefficient omega with confidence intervals.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
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
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
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
#>  syntax   term           estimate  se    test_statistic p_value  ci_lower
#>  f1 =~ y1 lambda_f1_1    0.454     0.071 6.405          < 0.0001 0.315   
#>  f1 =~ y2 lambda_f1_2    0.553     0.073 7.604          < 0.0001 0.411   
#>  f1 =~ y3 lambda_f1_3    0.644     0.071 9.084          < 0.0001 0.505   
#>  f1 =~ y4 lambda_f1_4    0.732     0.068 10.792         < 0.0001 0.599   
#>  f1 =~ y5 lambda_f1_5    0.796     0.065 12.258         < 0.0001 0.669   
#>  f1 ~~ f1 phi_f1         1.000     0.000 <NA>           <NA>     1.000   
#>  y1 ~~ y1 psi_f1_1       0.721     0.077 9.358          < 0.0001 0.570   
#>  y2 ~~ y2 psi_f1_2       0.700     0.079 8.914          < 0.0001 0.546   
#>  y3 ~~ y3 psi_f1_3       0.614     0.073 8.447          < 0.0001 0.472   
#>  y4 ~~ y4 psi_f1_4       0.460     0.065 7.117          < 0.0001 0.333   
#>  y5 ~~ y5 psi_f1_5       0.324     0.060 5.370          < 0.0001 0.205   
#>  y1 ~1    nu_f1_1        -0.025    0.068 -0.361         0.7179   -0.158  
#>  y2 ~1    nu_f1_2        0.086     0.071 1.207          0.2275   -0.053  
#>  y3 ~1    nu_f1_3        0.088     0.072 1.227          0.2197   -0.053  
#>  y4 ~1    nu_f1_4        -0.087    0.071 -1.228         0.2194   -0.225  
#>  y5 ~1    nu_f1_5        -0.006    0.069 -0.081         0.9358   -0.141  
#>           loading_sum_f1 3.179     0.204 15.585         < 0.0001 2.780   
#>           error_sum_f1   2.819     0.142 19.828         < 0.0001 2.540   
#>           omega_f1       0.782     0.024 32.304         < 0.0001 0.731   
#>           ave_f1         0.426     0.032 13.426         < 0.0001 0.365   
#>           H_f1           0.819     0.023 35.503         < 0.0001 0.769   
#>           chi_square     3.615     <NA>  <NA>           <NA>     <NA>    
#>           df             5         <NA>  <NA>           <NA>     <NA>    
#>           p_chi_square   0.6061    <NA>  <NA>           <NA>     <NA>    
#>           cfi            1.000     <NA>  <NA>           <NA>     <NA>    
#>           tli            1.011     <NA>  <NA>           <NA>     <NA>    
#>           nnfi           1.011     <NA>  <NA>           <NA>     <NA>    
#>           rmsea          0.000     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_lower 0.000     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_upper 0.083     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_level 0.900     <NA>  <NA>           <NA>     <NA>    
#>           srmr           0.018     <NA>  <NA>           <NA>     <NA>    
#>           AIC            2594.466  <NA>  <NA>           <NA>     <NA>    
#>           BIC            2643.941  <NA>  <NA>           <NA>     <NA>    
#>           H0             -1282.233 <NA>  <NA>           <NA>     <NA>    
#>           H1             -1280.426 <NA>  <NA>           <NA>     <NA>    
#>  ci_upper
#>  0.592   
#>  0.696   
#>  0.783   
#>  0.865   
#>  0.924   
#>  1.000   
#>  0.872   
#>  0.854   
#>  0.757   
#>  0.587   
#>  0.442   
#>  0.109   
#>  0.225   
#>  0.229   
#>  0.052   
#>  0.130   
#>  3.579   
#>  3.097   
#>  0.826   
#>  0.489   
#>  0.860   
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)

# Equal loadings (essentially tau-equivalent), named in the header.
cfa_1(X, equal_loading = TRUE)
#> Measurement structure, per factor:
#>   f1: essentially tau-equivalent (equal loadings)
#> 
#>  syntax   term           estimate  se    test_statistic p_value  ci_lower
#>  f1 =~ y1 lambda_f1      0.654     0.042 15.543         < 0.0001 0.572   
#>  f1 =~ y2 lambda_f1      0.654     0.042 15.543         < 0.0001 0.572   
#>  f1 =~ y3 lambda_f1      0.654     0.042 15.543         < 0.0001 0.572   
#>  f1 =~ y4 lambda_f1      0.654     0.042 15.543         < 0.0001 0.572   
#>  f1 =~ y5 lambda_f1      0.654     0.042 15.543         < 0.0001 0.572   
#>  f1 ~~ f1 phi_f1         1.000     0.000 <NA>           <NA>     1.000   
#>  y1 ~~ y1 psi_f1_1       0.699     0.081 8.679          < 0.0001 0.541   
#>  y2 ~~ y2 psi_f1_2       0.662     0.077 8.635          < 0.0001 0.512   
#>  y3 ~~ y3 psi_f1_3       0.609     0.072 8.519          < 0.0001 0.469   
#>  y4 ~~ y4 psi_f1_4       0.491     0.060 8.127          < 0.0001 0.373   
#>  y5 ~~ y5 psi_f1_5       0.417     0.054 7.741          < 0.0001 0.311   
#>  y1 ~1    nu_f1_1        -0.025    0.075 -0.328         0.7432   -0.172  
#>  y2 ~1    nu_f1_2        0.086     0.074 1.159          0.2463   -0.059  
#>  y3 ~1    nu_f1_3        0.088     0.072 1.222          0.2216   -0.053  
#>  y4 ~1    nu_f1_4        -0.087    0.068 -1.278         0.2012   -0.220  
#>  y5 ~1    nu_f1_5        -0.006    0.065 -0.086         0.9316   -0.133  
#>           loading_sum_f1 3.272     0.210 15.543         < 0.0001 2.859   
#>           error_sum_f1   2.878     0.144 19.931         < 0.0001 2.595   
#>           omega_f1       0.788     0.024 33.216         < 0.0001 0.738   
#>           ave_f1         0.432     0.035 12.318         < 0.0001 0.365   
#>           H_f1           0.794     0.024 33.140         < 0.0001 0.743   
#>           chi_square     22.763    <NA>  <NA>           <NA>     <NA>    
#>           df             9         <NA>  <NA>           <NA>     <NA>    
#>           p_chi_square   0.0068    <NA>  <NA>           <NA>     <NA>    
#>           cfi            0.945     <NA>  <NA>           <NA>     <NA>    
#>           tli            0.939     <NA>  <NA>           <NA>     <NA>    
#>           nnfi           0.939     <NA>  <NA>           <NA>     <NA>    
#>           rmsea          0.087     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_lower 0.043     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_upper 0.133     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_level 0.900     <NA>  <NA>           <NA>     <NA>    
#>           srmr           0.111     <NA>  <NA>           <NA>     <NA>    
#>           AIC            2605.614  <NA>  <NA>           <NA>     <NA>    
#>           BIC            2641.896  <NA>  <NA>           <NA>     <NA>    
#>           H0             -1291.807 <NA>  <NA>           <NA>     <NA>    
#>           H1             -1280.426 <NA>  <NA>           <NA>     <NA>    
#>  ci_upper
#>  0.737   
#>  0.737   
#>  0.737   
#>  0.737   
#>  0.737   
#>  1.000   
#>  0.857   
#>  0.812   
#>  0.749   
#>  0.610   
#>  0.522   
#>  0.123   
#>  0.230   
#>  0.229   
#>  0.046   
#>  0.122   
#>  3.684   
#>  3.161   
#>  0.831   
#>  0.501   
#>  0.837   
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)

# From a covariance matrix and sample size, as a paper reports them,
# here with a subset of its variables named in \'items\' (which selects
# columns of raw data in the same way).
cfa_1(S = cov(X), N = 200, items = c("y1", "y2", "y3", "y4"))
#> Measurement structure, per factor:
#>   f1: congeneric (no equality constraints)
#> 
#>  syntax   term           estimate  se    test_statistic p_value  ci_lower
#>  f1 =~ y1 lambda_f1_1    0.447     0.077 5.814          < 0.0001 0.296   
#>  f1 =~ y2 lambda_f1_2    0.606     0.079 7.675          < 0.0001 0.451   
#>  f1 =~ y3 lambda_f1_3    0.605     0.080 7.578          < 0.0001 0.448   
#>  f1 =~ y4 lambda_f1_4    0.732     0.080 9.178          < 0.0001 0.576   
#>  f1 ~~ f1 phi_f1         1.000     0.000 <NA>           <NA>     1.000   
#>  y1 ~~ y1 psi_f1_1       0.727     0.082 8.887          < 0.0001 0.567   
#>  y2 ~~ y2 psi_f1_2       0.639     0.085 7.532          < 0.0001 0.473   
#>  y3 ~~ y3 psi_f1_3       0.663     0.087 7.642          < 0.0001 0.493   
#>  y4 ~~ y4 psi_f1_4       0.460     0.089 5.148          < 0.0001 0.285   
#>           loading_sum_f1 2.390     0.174 13.765         < 0.0001 2.049   
#>           error_sum_f1   2.489     0.144 17.318         < 0.0001 2.207   
#>           omega_f1       0.696     0.035 20.034         < 0.0001 0.624   
#>           ave_f1         0.369     0.037 10.076         < 0.0001 0.300   
#>           H_f1           0.720     0.037 19.629         < 0.0001 0.643   
#>           chi_square     0.039     <NA>  <NA>           <NA>     <NA>    
#>           df             2         <NA>  <NA>           <NA>     <NA>    
#>           p_chi_square   0.9805    <NA>  <NA>           <NA>     <NA>    
#>           cfi            1.000     <NA>  <NA>           <NA>     <NA>    
#>           tli            1.048     <NA>  <NA>           <NA>     <NA>    
#>           nnfi           1.048     <NA>  <NA>           <NA>     <NA>    
#>           rmsea          0.000     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_lower 0.000     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_upper 0.000     <NA>  <NA>           <NA>     <NA>    
#>           rmsea_ci_level 0.900     <NA>  <NA>           <NA>     <NA>    
#>           srmr           0.003     <NA>  <NA>           <NA>     <NA>    
#>           AIC            2149.746  <NA>  <NA>           <NA>     <NA>    
#>           BIC            2176.133  <NA>  <NA>           <NA>     <NA>    
#>           H0             -1066.873 <NA>  <NA>           <NA>     <NA>    
#>           H1             -1066.853 <NA>  <NA>           <NA>     <NA>    
#>  ci_upper
#>  0.597   
#>  0.761   
#>  0.761   
#>  0.889   
#>  1.000   
#>  0.888   
#>  0.805   
#>  0.833   
#>  0.635   
#>  2.730   
#>  2.771   
#>  0.760   
#>  0.443   
#>  0.786   
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#>  <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)
```
