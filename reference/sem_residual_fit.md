# The SRMR With Its Unbiased Estimate, Confidence Interval, and Tests of Fit

Summarizes how far a fitted model's implied covariances are from the
observed ones with the standardized root mean square residual (SRMR) and
gives it the inferential treatment the RMSEA has long had: the
asymptotically unbiased estimate of the population SRMR of
Maydeu-Olivares (2017), its standard error and confidence interval, a
test of exact fit, and an optional test of close fit, together with the
largest standardized residual. The SRMR is an effect size of misfit in
the metric of a correlation, and Maydeu-Olivares, Shi, and Rosseel
(2018) found its interval and test of close fit more accurate than those
of the RMSEA in models with many variables. The computations are
lavaan's
([`lavaan::lavResiduals()`](https://rdrr.io/pkg/lavaan/man/lavResiduals.html)),
and the function requires lavaan 0.7-2 or newer (in Suggests).

## Usage

``` r
sem_residual_fit(fit, conf_level = 0.9, close_fit = NULL)
```

## Arguments

- fit:

  A fitted lavaan model: a maximum likelihood fit (estimator ML, MLM,
  MLMV, or MLR) under lavaan's default normal likelihood, of one group
  of continuous variables at a single level, from complete raw data
  (fitted under listwise deletion, FIML, or pairwise deletion) or from a
  covariance matrix. It must have converged to an admissible solution,
  as
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
  judges it. See Details for what is refused and why.

- conf_level:

  Confidence level of the interval for the population SRMR. Defaults to
  0.90, as lavaan's does and as every RMSEA interval in DMAR does: a 90
  percent interval pairs with a test of close fit at .05 (see Details).
  `conf_level = 0.95` gives a 95 percent interval.

- close_fit:

  Optional bound \\c\\ for the test of close fit, \\H_0: \mathrm{SRMR}
  \le c\\. Without it, no bound is imposed and the `close_fit_bound` and
  `p_close_fit` rows are `NA`. There is no universal value:
  Maydeu-Olivares, Shi, and Rosseel (2018) treat the bound as the
  researcher's choice. Shi, Maydeu-Olivares, and DiStefano (2018)
  proposed population reference values that scale with the average
  \\R^2\\ (communality) of the observed variables, in a criterion with
  two parts that must both hold: a model fits closely when its SRMR is
  at most .05 times the average \\R^2\\ and its largest absolute
  standardized residual covariance is at most .10 (and acceptably at .10
  times the average \\R^2\\ and .15). They recommend applying these
  population values through an interval or a test of close fit rather
  than to the sample values, so `close_fit = 0.05` times the average
  \\R^2\\ tests the first part (`lavaan::lavInspect(fit, "rsquare")`
  gives the \\R^2\\ values). The second part is not tested here; the row
  `largest_absolute_residual` is its sample value, which they caution
  against comparing with .10 unless the sample is very large. Ximénez,
  Maydeu-Olivares, Shi, and Revuelta (2022) found the bound of the first
  part, applied to the unbiased SRMR, to work reasonably well in finite
  samples, except that with small loadings (an average \\R^2\\ of .09)
  close fit called for samples of 500 or more, and recommend it.

## Value

A `dmar_tbl` (class `dmar_sem_fit`) with columns `term` and `value` and,
always in this order, the rows `srmr` (the sample SRMR),
`srmr_se_exact_fit` (its standard error under exact fit only),
`p_exact_fit`, `unbiased_srmr`, `unbiased_srmr_se`,
`unbiased_srmr_ci_lower`, `unbiased_srmr_ci_upper`, `close_fit_bound`
and `p_close_fit` (`NA` unless `close_fit` is given), `srmr_mean` (`NA`
without a mean structure), `largest_absolute_residual`, and `N`.
Attributes: `standardized_residuals`, the matrix of standardized
covariance residuals; `largest_residual`, the pair of variables with the
largest absolute residual and its signed value; and `sem_fit`, which
records the residual covariance, the information setting, the cases
used, and the lines of the printed footer.

## Details

### The SRMR and Its Unbiased Estimate

With \\p\\ observed variables, sample covariances \\s\_{ij}\\, and model
implied covariances \\\hat{\sigma}\_{ij}\\, the standardized residuals
are \\e\_{ij} = (s\_{ij} - \hat{\sigma}\_{ij}) / \sqrt{s\_{ii}
s\_{jj}}\\, and the SRMR is \\\sqrt{\sum\_{i \le j} e\_{ij}^2 / t}\\
over the \\t = p(p + 1)/2\\ variances and covariances (row `srmr`;
lavaan's `srmr_bentler_nomean`). Its population counterpart is the same
quantity computed from the population covariances, a measure of the size
of misfit that does not depend on the sample size. The sample SRMR
overestimates it, since sampling error adds to every residual.
Maydeu-Olivares (2017) removes the bias on the squared scale: with
\\\Xi\\ the asymptotic covariance matrix of the standardized residuals
and \\e\\ their vector, the unbiased estimate (row `unbiased_srmr`) is
\\k^{-1}\sqrt{\max(e'e - \mathrm{tr}\\\Xi, 0)/t}\\, with \\k = 1 -
\[\mathrm{tr}(\Xi^2) + 2e'\Xi e\]/\[4(e'e)^2\]\\ (Maydeu-Olivares et
al., 2018, Eq. 16). It is truncated at zero when the estimated sampling
contribution exceeds \\e'e\\, and the square root and the truncation
leave some bias in small samples.

### The Interval and the Test of Close Fit

The interval for the population SRMR is the large sample normal
reference interval \\\hat{\rho}\_u \pm z\_{1 - \alpha/2}
\widehat{SE}(\hat{\rho}\_u)\\ for the unbiased estimate
\\\hat{\rho}\_u\\ (Maydeu-Olivares et al., 2018, Eqs. 20 and 22), where
\\\alpha\\ is 1 minus `conf_level`. Its limits are reported as lavaan
returns them, without truncation, so the lower limit can fall below zero
although the population SRMR cannot, and when the estimate is truncated
to zero the interval is centered at zero. In the simulations of
Maydeu-Olivares et al. (2018) the interval computed without the
normality assumption covered at close to the nominal rate in most
conditions, including models of 60 variables, and did less well than the
RMSEA interval only for very small misfit in small models, and then
mostly in small samples; computed under normality, it covered poorly
when the data were not normal. The test of close fit, \\H_0:
\mathrm{SRMR} \le c\\, refers \\z = (\hat{\rho}\_u -
c)/\widehat{SE}(\hat{\rho}\_u)\\ to the standard normal (Maydeu-Olivares
et al., 2018, Eq. 21; row `p_close_fit`); it rejects at level \\\alpha\\
exactly when the lower limit of the \\100(1 - 2\alpha)\\ percent
interval exceeds \\c\\, so the default `conf_level = 0.90` pairs with a
test at .05, as the RMSEA interval of
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
does with its tests, and `conf_level = 0.95` with a test at .025. The
printed table states the level beneath it. Failing to reject close fit
does not establish it.

### The Test of Exact Fit

The test of exact fit (row `p_exact_fit`) asks whether the population
SRMR is zero. It uses the sample (biased) SRMR, whose expected value and
standard error under exact fit are known (Maydeu-Olivares, 2017): \\z =
\max(\mathrm{SRMR} - E_0, 0) / SE_0\\, referred to the upper tail of the
standard normal. The row `srmr_se_exact_fit` is that standard error,
\\SE_0 = \sqrt{\mathrm{tr}(\Xi^2)/(2t\\\mathrm{tr}\\\Xi)}\\, which holds
only under exact fit: it is not the standard error of the sample SRMR of
a misspecified model and must not be used to form an interval, which is
what the unbiased estimate and its standard error are for. Because of
the truncation at zero the *p*-value is at most 0.5, and it equals 0.5
whenever the sample SRMR is below its expectation under exact fit, as it
is more than half the time when the model is correct.

### Which Fits, and the Residual Covariance

The asymptotic covariance \\\Xi\\ of the residuals, and so the unbiased
estimate, its standard error, the interval, and both tests, depend on
the fit, as lavaan computes them. For an ML fit with standard errors
under normality it is normal theory. When the fit computed the
asymptotic covariance matrix of the sample moments from their fourth
order moments (MLM, MLMV, or a Satorra-Bentler test), or uses sandwich
standard errors (MLR), it is asymptotically distribution free, which is
the version to use with nonnormal item data, for which normal theory
intervals perform poorly (Maydeu-Olivares et al., 2018). The projection
onto the residual space uses the fit's information setting, so observed
information (the default under MLR) gives slightly different values from
expected information; an MLR fit with `information = "expected"`
reproduces MLM. The standardizing standard deviations are treated as
fixed. The `sem_fit` attribute records which residual covariance and
information were used.

The function refuses what these computations do not cover, each with a
message that says why: estimators other than ML (least squares and
categorical estimators), `likelihood = "wishart"`, multiple groups,
multilevel models, ordered variables, sampling weights,
`conditional.x = TRUE` (whose residuals are conditional on the exogenous
covariates rather than the marginal covariances the SRMR summarizes),
fits with incomplete data under FIML or pairwise deletion, and two-stage
fits. A FIML or pairwise fit of data with no missing values is accepted:
its residuals, their standard errors, and the tests are those of a
listwise fit with the same mean structure and information setting
(FIML's default is the observed information, which the footer names). A
fit with auxiliary variables (lavaan's `aux`) is refused: lavaan's
residual summary averages over them, and their residuals of about zero
lower the SRMR, while with complete data they change nothing in the
model's fit, so the fit without them is the one to summarize. A fit
whose cases were removed by listwise deletion is evaluated on the
complete cases with a warning, since that assumes the data are missing
completely at random. A model fitted from a covariance matrix
(`sample.cov`) gives the same values as one fitted to the raw data under
normal theory, since these need only the moments. For a saturated model
(zero degrees of freedom) the residuals are zero up to the convergence
tolerance, and the standard errors, limits, and *p*-values are `NA`,
with a warning that says so.

The SRMR here summarizes the covariance residuals only. For a model with
a mean structure, row `srmr_mean` is the root mean square of the
standardized mean residuals, zero when the means are saturated; when the
model restricts the means, a warning points to it, since the covariance
SRMR does not reflect misfit in the means. Always look at the residuals
themselves as well as their average (the attribute
`standardized_residuals`): a good average can hide a large local misfit.

## References

Maydeu-Olivares, A. (2017). Assessing the size of model misfit in
structural equation models. *Psychometrika, 82*(3), 533–558.
[doi:10.1007/s11336-016-9552-7](https://doi.org/10.1007/s11336-016-9552-7)

Maydeu-Olivares, A., Shi, D., & Rosseel, Y. (2018). Assessing fit in
structural equation models: A Monte-Carlo evaluation of RMSEA versus
SRMR confidence intervals and tests of close fit. *Structural Equation
Modeling, 25*(3), 389–402.
[doi:10.1080/10705511.2017.1389611](https://doi.org/10.1080/10705511.2017.1389611)

Shi, D., Maydeu-Olivares, A., & DiStefano, C. (2018). The relationship
between the standardized root mean square residual and model
misspecification in factor analysis models. *Multivariate Behavioral
Research, 53*(5), 676–694.
[doi:10.1080/00273171.2018.1476221](https://doi.org/10.1080/00273171.2018.1476221)

Ximénez, C., Maydeu-Olivares, A., Shi, D., & Revuelta, J. (2022).
Assessing cutoff values of SEM fit indices: Advantages of the unbiased
SRMR index and its cutoff criterion based on communality. *Structural
Equation Modeling, 29*(3), 368–380.
[doi:10.1080/10705511.2021.1992596](https://doi.org/10.1080/10705511.2021.1992596)

## See also

[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
for the chi square, RMSEA, and incremental indices beside the sample
SRMR, and
[`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md)
for the test of not close fit of the RMSEA.

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
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The three-factor model of nine holzinger_swineford tests.
fit <- lavaan::cfa(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
   textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
   speed   =~ t10_addition + t12_counting_groups_of_dots +
              t13_straight_and_curved_capitals",
  data = holzinger_swineford, std.lv = TRUE)
sem_residual_fit(fit)
#>  term                      value   
#>  srmr                      0.065   
#>  srmr_se_exact_fit         0.006   
#>  p_exact_fit               < 0.0001
#>  unbiased_srmr             0.058   
#>  unbiased_srmr_se          0.010   
#>  unbiased_srmr_ci_lower    0.042   
#>  unbiased_srmr_ci_upper    0.074   
#>  close_fit_bound           <NA>    
#>  p_close_fit               <NA>    
#>  srmr_mean                 <NA>    
#>  largest_absolute_residual 0.189   
#>  N                         301     
#> 
#> Confidence level: 90%
#> 
#> Residual covariance: normal theory; expected information
#> Largest standardized residual: t2_cubes with t10_addition, -0.189

# The 90 percent interval, the default, pairs with a test of close fit
# at .05. The bound is the researcher's choice; here it is .05 times the
# average R squared of the nine tests, the population value for close
# fit of Shi, Maydeu-Olivares, and DiStefano (2018).
r2 <- lavaan::lavInspect(fit, "rsquare")
sem_residual_fit(fit, close_fit = 0.05 * mean(r2))
#>  term                      value   
#>  srmr                      0.065   
#>  srmr_se_exact_fit         0.006   
#>  p_exact_fit               < 0.0001
#>  unbiased_srmr             0.058   
#>  unbiased_srmr_se          0.010   
#>  unbiased_srmr_ci_lower    0.042   
#>  unbiased_srmr_ci_upper    0.074   
#>  close_fit_bound           0.025   
#>  p_close_fit               0.0003  
#>  srmr_mean                 <NA>    
#>  largest_absolute_residual 0.189   
#>  N                         301     
#> 
#> Confidence level: 90%
#> 
#> Residual covariance: normal theory; expected information
#> Largest standardized residual: t2_cubes with t10_addition, -0.189

# A 95 percent interval instead.
sem_residual_fit(fit, conf_level = 0.95)
#>  term                      value   
#>  srmr                      0.065   
#>  srmr_se_exact_fit         0.006   
#>  p_exact_fit               < 0.0001
#>  unbiased_srmr             0.058   
#>  unbiased_srmr_se          0.010   
#>  unbiased_srmr_ci_lower    0.039   
#>  unbiased_srmr_ci_upper    0.077   
#>  close_fit_bound           <NA>    
#>  p_close_fit               <NA>    
#>  srmr_mean                 <NA>    
#>  largest_absolute_residual 0.189   
#>  N                         301     
#> 
#> Confidence level: 95%
#> 
#> Residual covariance: normal theory; expected information
#> Largest standardized residual: t2_cubes with t10_addition, -0.189
```
