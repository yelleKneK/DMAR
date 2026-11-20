# Expected Cross-Validation Index (ECVI) for a Covariance-Structure Model

The ECVI of Browne and Cudeck (1989) estimates how well a fitted model's
implied covariance matrix would fit an independent sample of the same
size from the same population. It is the single-sample estimate of the
cross-validation discrepancy, so a smaller ECVI indicates a model
expected to generalize better; ECVI is most useful for comparing
competing models fit to the same data. A confidence interval, derived
from the noncentral chi square distribution, accompanies the point
estimate.

## Usage

``` r
ecvi(
  fit = NULL,
  chisq = NULL,
  df = NULL,
  npar = NULL,
  multiplier = NULL,
  conf_level = 0.95
)
```

## Arguments

- fit:

  A fitted lavaan model. Supply this, or the summary statistics below.

- chisq, df, npar:

  The model chi square, its degrees of freedom, and the number of free
  parameters. Used when `fit` is not supplied, so an ECVI can be
  obtained from a published fit table.

- multiplier:

  The multiplier with which the chi square was formed from the minimized
  fit function, which is the divisor of the ECVI (\\n\\ in Details): the
  total sample size \\N\\ for a statistic formed as \\N\hat{F}\\, and
  \\N - G\\ for \\G\\ groups for a statistic formed as \\(N -
  G)\hat{F}\\. Given with `chisq`, `df`, and `npar` when `fit` is not
  supplied; with a fitted model it is taken from the fit.

- conf_level:

  Confidence level for the interval. Defaults to 0.95.

## Value

A `data.frame` (class `dmar_tbl`) with rows `ecvi`, `ci_lower`, and
`ci_upper` in the `value` column.

## Details

With \\q\\ free parameters, \\\mathrm{ECVI} = (\chi^2 + 2q)/n\\, where
\\n\\ is the multiplier with which the chi square was formed from the
minimized fit function (the argument `multiplier`), \\\chi^2 =
n\hat{F}\\, so that \\\mathrm{ECVI} = \hat{F} + 2q/n\\ (Browne & Cudeck,
1989). With a fitted model the divisor is the one lavaan uses, and the
value is the one
[`lavaan::fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
reports: \\n = N\\, the total sample size, for maximum likelihood under
lavaan's default normal likelihood (including MLR, whose unscaled
statistic is used, and FIML), and \\n = N - G\\ for \\G\\ groups under
`likelihood = "wishart"` and for the least squares estimators (GLS, WLS,
ULS, DWLS). With summary statistics, `multiplier` is used as given: pass
\\N - G\\ when the published chi square was formed as \\(N -
G)\hat{F}\\, as by programs that follow the Wishart convention. Writing
\\\hat\lambda = \chi^2 - df\\ for the estimated noncentrality, the ECVI
is \\(\hat\lambda + df + 2q)/n\\; the confidence interval replaces
\\\hat\lambda\\ by the lower and upper noncentrality limits from
[`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md),
the same inversion used for the RMSEA interval (see
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)).

The noncentrality limits are at least 0, while the point estimate is the
untruncated \\\hat{F} + 2q/n\\ of Browne and Cudeck (1989), so for a
model with \\\chi^2 \< df\\ the estimate lies below the lower limit,
\\(df + 2q)/n\\. When the chi square lies at or below the upper
\\\alpha/2\\ quantile of the central chi square distribution, for
\\\alpha\\ = 1 - `conf_level`, the lower noncentrality limit is 0 and
the lower limit is that same \\(df + 2q)/n\\. When it lies at or below
the lower \\\alpha/2\\ quantile, no noncentrality is consistent with it
at that level, and the upper limit is `NA`, with a warning.
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
reports the upper RMSEA limit as 0 in that case, so that it agrees with
the test of not close fit, but the ECVI has no such test, and an upper
limit of \\(df + 2q)/n\\ would make an interval of zero width that
excludes the estimate.

A fitted model must have converged to an admissible solution, as for
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md):
a solution with a negative variance estimate, or with a covariance
matrix of the latent variables, of the residuals of the observed
variables, or of the structural residuals that is not positive definite,
is refused (with an error of class `dmar_inadmissible`), since its fit
describes estimates outside the parameter space; with auxiliary
variables the check is made on the model's own parameters.

When the chi square is the likelihood ratio statistic \\2(\ell_s -
\ell_M)\\, as under lavaan's default normal likelihood, \\N \cdot
\mathrm{ECVI} = \mathrm{AIC} + 2\ell_s\\, where \\\ell_s\\ is the
maximized log likelihood of the saturated model and \\\ell_M\\ that of
the fitted model. ECVI and AIC then rank models fitted to the same data
identically, and their values differ only by the saturated log
likelihood constant before division by \\N\\. Under the Wishart
convention the statistic is \\(N - G)/N\\ times the likelihood ratio
statistic, and the identity holds only approximately. ECVI is reported
because its metric (a discrepancy per observation) and its confidence
interval are interpretable on their own.

## References

Browne, M. W., & Cudeck, R. (1989). Single sample cross-validation
indices for covariance structures. *Multivariate Behavioral Research,
24*(4), 445–455.

## See also

[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md).

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
# From a published fit table (no model object needed).
ecvi(chisq = 24.361, df = 8, npar = 13, multiplier = 301)
#>  term     value
#>  ecvi     0.167
#>  ci_lower 0.125
#>  ci_upper 0.243
#> 
#> Confidence level: 95%

fit <- lavaan::cfa(
  "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
   verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
  data = holzinger_swineford, std.lv = TRUE)
ecvi(fit)
#>  term     value
#>  ecvi     0.167
#>  ci_lower 0.125
#>  ci_upper 0.243
#> 
#> Confidence level: 95%
```
