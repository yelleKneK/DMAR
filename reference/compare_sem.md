# Compare Structural Equation Models by an Information Criterion

Ranks a set of structural equation models fitted to the same data by one
information criterion of
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
(BIC by default) and reports, for each model, its difference from the
best, its relative weight within the set, and, for the criteria that
approximate a log marginal likelihood, the approximate log Bayes factor
of the best model over it. The models may be nested or not. A criterion
compares models only through their log likelihoods, so before anything
is computed the function checks that the log likelihoods are comparable:
the same cases and data values, the same observed variables and groups,
the same missing data method and mean structure, and the same exogenous
observed variables conditioned on. Requires lavaan (in Suggests).

## Usage

``` r
compare_sem(
  fits,
  criterion = c("BIC", "SABIC", "HBIC", "IBIC", "SPBIC", "AIC"),
  prior_means = NULL
)

# S3 method for class 'dmar_compare_sem'
tidy(x, ...)

# S3 method for class 'dmar_compare_sem'
glance(x, ...)
```

## Arguments

- fits:

  A named list of two or more fitted lavaan models: maximum likelihood
  fits of continuous variables (estimator ML or its robust variants MLM,
  MLMV, and MLR, which maximize the same likelihood) fitted to the same
  data. That is, the same cases, data values, observed variables (in any
  order), and groups; the same missing data method (`missing`) and mean
  structure (`meanstructure`); and, when `fixed.x` or `conditional.x` is
  `TRUE`, the same exogenous observed variables. Models fitted to
  summary statistics (`sample.cov`) are compared when they share the
  sample moments and *N*. Two-level models and models fitted with
  sampling weights are refused, and the criterion can impose further
  conditions (see
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)).

- criterion:

  The criterion: `"BIC"` (the default), `"SABIC"`, `"HBIC"`, `"IBIC"`,
  `"SPBIC"`, or `"AIC"`, in either case.

- prior_means:

  For SPBIC, an optional list of prior mean vectors, one per model in
  the order of `fits` (and with its names, if named), each in the order
  of `lavaan::coef()` of its model. `NULL` applies the default rule of
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md):
  free loadings and regression coefficients centered at zero and the
  other free parameters at their estimates. A saturated model needs an
  explicit prior mean.

- x:

  A `dmar_compare_sem` table returned by `compare_sem()`.

- ...:

  Unused.

## Value

A `dmar_tbl` (class `dmar_compare_sem`) with one row per model, sorted
by the criterion (ties keep the order of `fits`), and the columns `term`
(the model's name), `value` (its criterion), `delta`, `relative_weight`,
`log_evidence_best_over_model`, and `npar` (its number of free
parameters). The `criterion` attribute names the criterion, the
`comparison` attribute records what was checked (*N*, the number of
groups, whether the models were fitted to raw data, the `missing`
settings used, the mean structure, the exogenous observed variables
conditioned on, and for IBIC and SPBIC the identification), and the
`sem_information` attribute holds the details of each model's
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
result, named by model. `tidy()` returns every column with the criterion
beside it, and `glance()` one row with the criterion, the number of
models, *N*, the best model, its weight, and the difference of the
second best.

## Details

### Differences, Weights, and Evidence

For the criterion values \\C_k\\ of the models, \\\Delta_k = C_k -
\min_j C_j\\ (column `delta`), and the relative weight of model *k* is
\\\exp(-\Delta_k/2)/\sum_j \exp(-\Delta_j/2)\\. BIC, HBIC, IBIC, and
SPBIC each approximate \\-2\\ times the log marginal likelihood of a
model (Bollen, Ray, Zavisca, & Harden, 2012), so \\\exp(-\Delta_k/2)\\
approximates the Bayes factor of model *k* against the best model, the
weights approximate the posterior probabilities of the models when every
model in the set has the same prior probability (Burnham & Anderson,
2004), and \\\Delta_k/2\\ (column `log_evidence_best_over_model`)
approximates the log Bayes factor of the best model over model *k*. For
AIC the weights are Akaike weights, the weight of evidence that each
model is the one of the set closest to the data generating process in
Kullback-Leibler discrepancy (Burnham & Anderson, 2004); they are not
posterior probabilities under equal prior odds, and the log evidence
column is `NA`. SABIC rests on the minimum description length principle
rather than on a marginal likelihood, so its weights carry no further
interpretation and its log evidence column is `NA` too. So it is for
every criterion of models fitted by two-stage maximum likelihood
(`missing = "two.stage"` or `"robust.two.stage"`), whose criteria are
those of Lai (2021): built from the complete data normal log likelihood
at the EM estimates of the means and covariances, they select in large
samples the model complete data would select, but they are not a
likelihood of the observed data and approximate no marginal likelihood
(see
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)).

Every weight is conditional on the set of models compared: adding or
removing a model changes all of them, and none of them says how well a
model fits in absolute terms
([`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
speaks to that). Nor do the weights carry the sampling variability of
the criteria; the uncertainty of the selection itself is not estimated.

### What Is Checked

Each model must have converged to an admissible solution and have a log
likelihood the criteria are defined on (see
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
for what is refused and why). Then every model is compared with the
first, and the first difference found is reported with the two models it
concerns and what to do.

- Data:

  The same groups, the same observed variables, the same number of
  cases, the same cases, and the same data values (or, for models fitted
  to summary statistics, the same sample moments). The order of the
  variables does not matter; a set of different variables does, since a
  log likelihood is a density of the variables it models.

- Missing data:

  The same missing data method. Listwise deletion, FIML, and the
  two-stage method give different likelihoods, and two listwise fits of
  different variables keep different cases. The two two-stage settings,
  `missing = "two.stage"` and `"robust.two.stage"`, count as one method:
  they fit the same EM moments by the same likelihood and differ only in
  their standard errors and test statistic, which no criterion uses.
  Two-stage fits must also share the auxiliary variables of their first
  stage (lavaan's `aux`), since those change the EM estimates their log
  likelihoods are evaluated at. A two-stage fit whose first stage read
  cases observed on an auxiliary variable but on none of the model's
  variables, which lavaan leaves out of its *N*, is compared by the
  criterion
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
  computes for it with every case the first stage read, its log
  likelihood and the *N* of the penalty included (in several groups held
  together by a constraint, at the estimates with every case read, which
  move with the groups' counts), and a line under the table says so;
  such a fit and the fit
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)`(missing = "two_stage")`
  or
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  returns for the same data count the same cases.

- Mean structure:

  The same `meanstructure` setting: without a mean structure lavaan
  leaves the saturated means out of the parameter count, which moves
  every criterion by an amount that depends on the setting.

- Exogenous variables:

  With `fixed.x = TRUE` (the default of
  [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html)) or
  `conditional.x = TRUE`, lavaan's log likelihood is the density
  conditional on the exogenous observed variables, the joint density
  less the marginal density of those variables, and its parameter count
  leaves out their moments. Two models that treat different variables as
  exogenous therefore have log likelihoods of different conditional
  distributions, whose difference can be hundreds of points for two
  saturated models of the same three variables. Models with the same
  exogenous variables are compared, and the others are refused with the
  advice to refit every model with `fixed.x = FALSE`.

- Identification:

  For IBIC, and for SPBIC under its default prior, every latent variable
  of every model must be identified the same way (each by a marker
  indicator, or each by a fixed variance), since both criteria change
  with the identification even though the likelihood does not (see
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)).
  A set that uses both is refused, also when each model mixes them in
  the same way: one identification throughout is a convention that can
  be stated and reported, and it catches the models that differ, which
  may share no latent variable to compare. Within one data set the units
  of the observed variables are shared, but an IBIC ranking still
  depends on them.

## References

Bollen, K. A., Harden, J. J., Ray, S., & Zavisca, J. (2014). BIC and
alternative Bayesian information criteria in the selection of structural
equation models. *Structural Equation Modeling, 21*(1), 1–19.
[doi:10.1080/10705511.2014.856691](https://doi.org/10.1080/10705511.2014.856691)

Bollen, K. A., Ray, S., Zavisca, J., & Harden, J. J. (2012). A
comparison of Bayes factor approximation methods including two new
methods. *Sociological Methods & Research, 41*(2), 294–324.
[doi:10.1177/0049124112452393](https://doi.org/10.1177/0049124112452393)

Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
Understanding AIC and BIC in model selection. *Sociological Methods &
Research, 33*(2), 261–304.
[doi:10.1177/0049124104268644](https://doi.org/10.1177/0049124104268644)

Lai, K. (2021). Using information criteria under missing data: Full
information maximum likelihood versus two-stage estimation. *Structural
Equation Modeling, 28*(2), 278–291.
[doi:10.1080/10705511.2020.1780925](https://doi.org/10.1080/10705511.2020.1780925)

## See also

[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
for the criteria of one model and what each assumes,
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
for the fit of each model in absolute terms, and the vignette
[`vignette("sem_fit_measures", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/sem_fit_measures.md).

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
# One factor against two correlated factors for six of the
# holzinger_swineford tests.
one <- lavaan::cfa(
  "g =~ t1_visual_perception + t2_cubes + t4_lozenges +
        t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
  data = holzinger_swineford)
two <- lavaan::cfa(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
   textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
  data = holzinger_swineford)
x <- compare_sem(list(one_factor = one, two_factors = two))
x
#>  term        value     delta  relative_weight log_evidence_best_over_model npar
#>  two_factors 10947.060 0.000  1.000           0.000                        13  
#>  one_factor  11020.221 73.161 0.000           36.581                       12  
#> 
#> Criterion: BIC (smaller is preferred); N = 301; 2 models.
#> Checked: the same cases and data values, observed variables, missing data method (missing = "listwise"), and mean structure (meanstructure = FALSE); no exogenous variables conditioned on.
#> relative_weight: approximate posterior probability of each model within this set under equal prior odds; log_evidence_best_over_model: approximate log Bayes factor of the best model over each.
generics::glance(x)
#>   criterion n_models   N  best_model best_weight delta_second
#> 1       BIC        2 301 two_factors           1     73.16114

# The same comparison by SPBIC, under its default prior.
compare_sem(list(one_factor = one, two_factors = two), criterion = "SPBIC")
#>  term        value     delta  relative_weight log_evidence_best_over_model npar
#>  two_factors 10939.120 0.000  1.000           0.000                        13  
#>  one_factor  11008.596 69.475 0.000           34.738                       12  
#> 
#> Criterion: SPBIC (smaller is preferred); N = 301; 2 models.
#> Checked: the same cases and data values, observed variables, missing data method (missing = "listwise"), and mean structure (meanstructure = FALSE); no exogenous variables conditioned on; identification: marker indicator.
#> relative_weight: approximate posterior probability of each model within this set under equal prior odds; log_evidence_best_over_model: approximate log Bayes factor of the best model over each.
```
