# Coefficient Omega (McDonald) With a Confidence Interval

Estimates McDonald's (1999) coefficient \\\omega\\ for a homogeneous
composite score from a single-factor confirmatory factor analysis model
and returns a confidence interval for the population coefficient. The
`denominator` argument selects whether the total variance in the
denominator of \\\omega\\ is estimated directly from the data
(`"observed"`, the default: robust omega) or taken from the fitted model
(`"model_implied"`); see *Details* for the properties of each choice. A
bootstrap confidence interval is never run unless requested; see
`ci_method`.

## Usage

``` r
reliability_omega(
  data = NULL,
  S = NULL,
  N = NULL,
  ci_method = c("mlr", "ml", "mlr_logistic", "ml_logistic", "likelihood", "adf",
    "adf_logistic", "feldt", "fisher", "bonett", "hakstian_whalen", "bootstrap_se",
    "bootstrap_se_logistic", "percentile", "bca", "none"),
  denominator = c("observed", "model_implied"),
  missing = c("listwise", "fiml"),
  aux = NULL,
  conf_level = 0.95,
  B = 10000,
  seed = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of item scores (rows are respondents,
  columns are items). Either `data` or `S` must be supplied. How
  incomplete rows are handled is governed by `missing`: listwise
  deletion by default, or full information maximum likelihood with
  `missing = "fiml"`. When `aux` is supplied, `data` contains both the
  item columns and the auxiliary columns, and the columns not named in
  `aux` are the items.

- S:

  A symmetric covariance matrix among the items. If supplied, `N` must
  also be supplied; methods that require raw data (`"mlr*"`, `"adf*"`,
  `"bootstrap_*"`, `"percentile"`, `"bca"`) are then unavailable, as are
  `missing = "fiml"` and `aux`, since a covariance matrix has no
  incomplete cases.

- N:

  Total sample size; required when `S` is supplied.

- ci_method:

  Method for constructing the confidence interval. See *Details*. When
  not supplied: for `denominator = "model_implied"` the default is
  `"mlr"` with raw data (use `"ml"` with covariance input); for
  `denominator = "observed"` (robust omega) no interval is computed by
  default, because its interval is bootstrap based and a bootstrap is
  never run unless requested. Ask for `"percentile"` or `"bca"` to
  obtain the recommended interval.

- denominator:

  How the total variance in the denominator of \\\omega\\ is estimated:
  `"observed"` (default; robust omega) uses the variance of the
  composite estimated directly from the data; `"model_implied"` uses the
  total variance reproduced by the fitted single-factor model. See
  *Details*. With `"observed"`, `ci_method` must be a bootstrap method
  or `"none"`.

- missing:

  How incomplete rows of `data` are handled: `"listwise"` (the default;
  complete-case analysis, the historical behavior) or `"fiml"` (full
  information maximum likelihood, using every case with at least one
  observed item). See *Details*.

- aux:

  Optional character vector naming auxiliary variable columns of `data`,
  entered as saturated correlates under full information maximum
  likelihood. Supplying `aux` implies `missing = "fiml"`. See *Details*.

- conf_level:

  Confidence level for the interval. Defaults to `0.95`.

- B:

  Number of bootstrap replications when a bootstrap method is selected.
  Defaults to `10000`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, which leaves the user's current RNG state intact; supply an
  integer for reproducibility.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample coefficient \\\omega\\), `"se"` (standard error from the chosen
method, `NA` for closed-form transformation CIs), `"lower_limit"` and
`"upper_limit"` (clamped to \[0, 1\]), `"conf_level"`, `"N"` (the cases
the analysis used: the complete cases under listwise deletion, every
case with at least one observed item under `missing = "fiml"`),
`"N_complete"` (the complete cases; equal to `"N"` under listwise
deletion and for covariance input), and `"J"`. Attributes `coefficient`
(`"omega"`), `ci_method`, `denominator`, `missing`, and (when supplied)
`aux` record the computation; bootstrap calls also record `B`.

## Details

Coefficient \\\omega\\ is a population coefficient of determination.
Write the composite \\Y = \sum_j X_j\\ through its measurement
decomposition, \\Y = E\[Y \mid \eta\] + e\\, the regression of the
observed composite on the latent variable \\\eta\\ it is intended to
measure. By the law of total variance, the proportion of composite
variance attributable to \\\eta\\ is \\\mathrm{Var}(E\[Y \mid \eta\]) /
\mathrm{Var}(Y)\\ (McDonald, 1999, 2011). Under the congeneric
single-factor model, in which item \\j\\ has factor loading
\\\lambda_j\\ and error variance \\\psi_j^{2}\\, the numerator equals
\\\left(\sum_j \lambda_j\right)^2\\ and the population coefficient is
\$\$\omega = \frac{\left(\sum\_{j}
\lambda\_{j}\right)^{2}}{\sigma\_{Y}^{2}},\$\$ with \\\sigma\_{Y}^{2}\\
the variance of the composite. Coefficient \\\omega\\ relaxes the
equal-loadings assumption underlying coefficient \\\alpha\\ (Guttman,
1945; Cronbach, 1951), so on congeneric scales \\\omega\\ estimates
population reliability whereas \\\alpha\\ underestimates it.

**The denominator argument.** The sample estimate substitutes lavaan
estimates of the loadings into the numerator; the two `denominator`
settings differ in how \\\sigma\_{Y}^{2}\\ is estimated.

- `"observed"` (default):

  The variance of the composite estimated directly from the data (the
  sum of all elements of the unrestricted covariance matrix, on the same
  maximum likelihood divisor as the fitted loadings). This is the
  coefficient that Kelley and Pornprasertmanit (2016) call hierarchical
  omega (`MBESS::ci.reliability(type = "hierarchical")`); the DMAR
  documentation refers to it as *robust omega*, and it is the default
  because its denominator estimates the variance of the composite
  consistently whether or not the single-factor model is correctly
  specified.

- `"model_implied"`:

  The total variance reproduced by the fitted single-factor model,
  \\\left(\sum_j \hat\lambda_j\right)^2 + \sum_j \hat\psi_j^{2}\\. This
  is the textbook form of \\\hat\omega\\, correct exactly when the
  one-factor model reproduces the composite variance.

The choice matters only insofar as the single-factor model is
misspecified, and the properties of each setting can be stated exactly.

- If the single-factor model is correctly specified, the two definitions
  coincide in the population, and both estimators converge to the same
  value. Nothing is given up, in that sense, by either choice.

- The observed total variance is a consistent estimator of
  \\\mathrm{Var}(Y)\\ whether or not the single-factor model is
  correctly specified; the model implied total variance is consistent
  for \\\mathrm{Var}(Y)\\ only when the model is correct.

- Under misspecification (for example, a minor unmodeled factor or
  correlated errors), only `"observed"` retains the coefficient of
  determination interpretation: the proportion of the variance of the
  composite users actually compute that is attributable to the fitted
  common factor. With `"model_implied"` the estimate becomes a ratio of
  two model derived quantities whose denominator is no longer the
  variance of any composite a user scores.

- Under a correctly specified model with normal items, the model implied
  denominator uses the model structure and can be a slightly more
  efficient estimator of \\\sigma\_{Y}^{2}\\ in finite samples. In the
  simulations of Kelley and Pornprasertmanit (2016) the practical
  differences between the two coefficients were negligible when the
  model held, while interval coverage under modest model error favored
  the observed denominator paired with a bootstrap interval.

Those simulation results are why Kelley and Pornprasertmanit (2016)
recommend the observed denominator with a bootstrap confidence interval
when unidimensionality is only approximate, which is common with real
items.

Two cautions frame the choice. First, no denominator repairs a
misspecified measurement model: the fitted loadings absorb part of
whatever structure the single-factor model omits, so the numerator is
affected under either setting. Assess the single-factor model (for
example with
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)) before
interpreting any \\\omega\\ variant, and model real multidimensionality
directly rather than patching over it. Second, the naming history is
worth knowing. The observed denominator coefficient was introduced as
hierarchical omega (Kelley & Pornprasertmanit, 2016; MBESS type
`"hierarchical"`), a name motivated by a hierarchical factor logic:
model misfit is viewed as a set of minor common factors (visible as
residual correlations), the single factor is retained as an
approximation, and the coefficient isolates the variance attributable to
the general factor alone, expressed relative to the observed variance of
the unweighted composite. It is not the bifactor coefficient
\\\omega_H\\ of Zinbarg, Revelle, Yovel, and Li (2005), whose numerator
comes from the general factor loadings of an explicitly multidimensional
model. Upon reflection, the authors would have named the coefficient for
its behavior rather than for the hierarchical motivation, as observed
omega or robust omega; DMAR uses *robust omega*, with the qualifications
that word requires. The robustness is to misspecification of the *total
variance* only, since the numerator remains model based under either
denominator; it is not the outlier robustness of Zhang and Yuan (2016),
and it is separate from the robust maximum likelihood standard errors
available through `ci_method`. Robust omega also shares a design
principle with categorical omega: in both, the total variance in the
denominator is not taken from the fitted factor model. The reliability
vignette develops this framing.

Available confidence interval methods (set via `ci_method`):

- `"ml"`, `"ml_logistic"`:

  Wald interval using the maximum likelihood standard error from lavaan
  (Raykov, 2002). The `_logistic` variant applies Browne's (1982) logit
  transformation.

- `"mlr"`, `"mlr_logistic"`:

  Wald interval using the robust (Satorra & Bentler, 1994) standard
  error. Default; recommended among closed-form methods when item
  distributions deviate from normality (Kelley & Pornprasertmanit,
  2016).

- `"likelihood"`:

  Profile likelihood interval: the set of population values not rejected
  by the likelihood ratio test, located by refitting the model under the
  nonlinear constraint that the model implied reliability equals each
  candidate value. It respects the \[0, 1\] range and is not forced to
  be symmetric about the estimate. Maximum likelihood; available with
  raw data or covariance input, for `denominator = "model_implied"`
  only.

- `"adf"`, `"adf_logistic"`:

  Wald interval using weighted least squares (“ADF”) estimation (Browne,
  1984). Requires raw data and a relatively large sample size.

- `"feldt"`, `"fisher"`, `"bonett"`, `"hakstian_whalen"`:

  Closed-form intervals derived for coefficient \\\alpha\\. They apply
  mechanically to \\\omega\\ but their coverage performance for
  \\\omega\\ is generally inferior to the maximum likelihood and
  bootstrap intervals (Kelley & Pornprasertmanit, 2016).

- `"bootstrap_se"`, `"bootstrap_se_logistic"`, `"percentile"`, `"bca"`:

  Nonparametric bootstrap intervals (Efron & Tibshirani, 1993): the rows
  of `data` are resampled with replacement `B` times and \\\omega\\ is
  recomputed, with the factor model refit, on each replication.
  `"percentile"` takes the interval limits from the empirical quantiles
  of the bootstrap estimates; it respects \[0, 1\] and is not forced to
  be symmetric about the estimate. `"bca"` (bias-corrected and
  accelerated) adjusts the two quantile positions for median bias,
  estimated from the bootstrap distribution, and for the rate at which
  the estimator's variance changes with the parameter, the acceleration,
  estimated by the jackknife; those two adjustments make it second-order
  accurate where the percentile interval is first-order accurate
  (DiCiccio & Efron, 1996). `"bootstrap_se"` uses the standard deviation
  of the bootstrap estimates as a standard error in a normal-theory
  interval, built on the logit scale for the `_logistic` variant so the
  endpoints respect \[0, 1\]. Replications whose model refit does not
  converge are dropped, and the interval is computed from the
  replications that return a value. The default `B = 10000` is an
  accuracy choice: the BCa adjustment pushes the working quantiles
  farther into the tails of the bootstrap distribution than the
  percentile interval uses, and stabilizing them takes more replications
  than the customary 2000; reduce `B` for exploration, not for a
  reported analysis. Recommended when assumptions of parametric methods
  are questionable. Requires raw data and the boot package; supply
  `seed` for run-to-run reproducibility.

- `"none"`:

  Return only the point estimate.

With `denominator = "observed"`, only the bootstrap methods and `"none"`
are available: the delta method standard errors and the alpha-derived
closed forms are derived under the model implied ratio and do not
account for sampling of the observed denominator. This pairing is not a
limitation in practice, since the bootstrap is the interval Kelley and
Pornprasertmanit (2016) recommend for the observed denominator
coefficient in any case. Because no analysis in DMAR runs a bootstrap
unless the user requests one, the default for robust omega is the point
estimate with no interval, accompanied by a message naming the call that
produces the recommended interval; request `ci_method = "percentile"` or
`"bca"` to obtain it.

**Missing data and auxiliary variables.** By default incomplete rows are
listwise-deleted; `missing = "fiml"` keeps every case with at least one
observed item and fits the single-factor model by full information
maximum likelihood, which is consistent and efficient under the missing
at random (MAR) assumption. The `aux` argument names auxiliary
variables: columns of `data` that are not part of the composite but are
correlated with the items or with the reasons values are missing. They
enter as *saturated correlates* (Graham, 2003): correlated freely with
each other and with every item's residual, never loading on the factor,
so the measurement model is undisturbed while FIML uses their
information; a good auxiliary also makes MAR itself more plausible
(Collins, Schafer, & Kam, 2001). Supplying `aux` implies
`missing = "fiml"`; combining it with an explicit `missing = "listwise"`
is an error, and listwise deletion remains the default so no existing
result changes. Under `missing = "fiml"` the denominators are estimated
as follows: with `"model_implied"` the fitted FIML model supplies the
total variance directly, and with `"observed"` the total variance is the
sum of the FIML estimate of the item covariance matrix (from a saturated
model over the items and any auxiliaries), which is the estimate of
\\\mathrm{Var}(Y)\\ that uses the partially observed rows; it is already
on the maximum likelihood divisor, matching the fitted loadings. The
available intervals under `"fiml"` are `"ml"`, `"mlr"` (Yuan & Bentler,
2000, robust), their `_logistic` variants, and the bootstrap methods,
which resample rows including the partially observed ones and refit by
FIML; the complete-data closed forms, `"adf"`, and `"likelihood"` are
errors rather than silent fallbacks to listwise deletion. Multiple
imputation is a different feature with a different interface and is out
of scope here.

**Comparison with other packages.** The psych package provides
[`omega`](https://rdrr.io/pkg/psych/man/omega.html), which fits a
Schmid-Leiman hierarchical factor model and reports several variants of
\\\omega\\ (\\\omega_t\\, \\\omega_h\\) alongside extensive psychometric
diagnostics. `reliability_omega` in DMAR differs in emphasis: it
implements McDonald's \\\omega\\ from a single-factor (congeneric) model
and accompanies the point estimate with a confidence interval drawn from
the methods compared in Kelley and Pornprasertmanit (2016). The same
denominator distinction appears in semTools' `compRelSEM()` as its
`obs.var` argument.

## References

Bonett, D. G. (2002). Sample size requirements for testing and
estimating coefficient alpha. *Journal of Educational and Behavioral
Statistics, 27*(4), 335–340.
[doi:10.3102/10769986027004335](https://doi.org/10.3102/10769986027004335)

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
*Topics in applied multivariate analysis* (pp. 72–141). Cambridge, UK:
Cambridge University Press.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*, 62–83.

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
standard errors in covariance structure analysis. In A. von Eye & C. C.
Clogg (Eds.), *Latent variables analysis: Applications for developmental
research* (pp. 399–419). Thousand Oaks, CA: Sage.

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
*Sociological Methodology, 30*, 165–200.

Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
inclusive and restrictive strategies in modern missing data procedures.
*Psychological Methods, 6*, 330–351.
[doi:10.1037/1082-989X.6.4.330](https://doi.org/10.1037/1082-989X.6.4.330)

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika, 16*(3), 297–334.

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science, 11*(3), 189–228.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Feldt, L. S. (1965). The approximate sampling distribution of
Kuder-Richardson reliability coefficient twenty. *Psychometrika, 30*,
357–370.

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

Guttman, L. (1945). A basis for analyzing test-retest reliability.
*Psychometrika, 10*(4), 255–282.

Hakstian, A. R., & Whalen, T. E. (1976). A *k*-sample significance test
for independent alpha coefficients. *Psychometrika, 41*, 219–231.

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
[doi:10.3758/BF03192993](https://doi.org/10.3758/BF03192993)

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

McDonald, R. P. (1999). *Test theory: A unified treatment*. Mahwah, NJ:
Lawrence Erlbaum Associates.

McDonald, R. P. (2011). Measuring latent quantities. *Psychometrika,
76*, 511–536.
[doi:10.1007/s11336-011-9223-7](https://doi.org/10.1007/s11336-011-9223-7)

Raykov, T. (2002). Analytic estimation of standard error and confidence
interval for scale reliability. *Multivariate Behavioral Research, 37*,
89–103.
[doi:10.1207/S15327906MBR3701_04](https://doi.org/10.1207/S15327906MBR3701_04)

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. [doi:10.1007/BF02296192](https://doi.org/10.1007/BF02296192)

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

Zhang, Z., & Yuan, K.-H. (2016). Robust coefficients alpha and omega and
confidence intervals with outlying observations and missing data:
Methods and software. *Educational and Psychological Measurement, 76*,
387–411.
[doi:10.1177/0013164415594658](https://doi.org/10.1177/0013164415594658)

Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's
\\\alpha\\, Revelle's \\\beta\\, and McDonald's \\\omega_H\\: Their
relations with each other and two alternative conceptualizations of
reliability. *Psychometrika, 70*, 123–133.
[doi:10.1007/s11336-003-0974-7](https://doi.org/10.1007/s11336-003-0974-7)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
(categorical omega for ordered items),
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
(single-factor CFA used internally),
[`omega`](https://rdrr.io/pkg/psych/man/omega.html).

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
set.seed(113)
J <- 6
loadings <- seq(0.5, 0.8, length.out = J)
eta <- rnorm(200)
errors <- matrix(rnorm(200 * J), 200, J) %*% diag(sqrt(1 - loadings^2))
items <- sweep(matrix(rep(eta, J), 200, J), 2, loadings, `*`) + errors
colnames(items) <- paste0("y", seq_len(J))

# Default: robust omega, point estimate only (no bootstrap is run
# unless requested; a message names the call that produces the
# recommended interval).
reliability_omega(data = items)
#> Robust omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "percentile" (or "bca"); B = 10000 replications is the default when you do.
#>  term        value
#>  estimate    0.821
#>  se          <NA> 
#>  lower_limit <NA> 
#>  upper_limit <NA> 
#>  conf_level  0.95 
#>  N           200  
#>  N_complete  200  
#>  J           6    

# Robust omega with the recommended bootstrap interval
# (few replications for a quick example).
reliability_omega(data = items, ci_method = "percentile", B = 200)
#>  term        value 
#>  estimate    0.821 
#>  se          0.0204
#>  lower_limit 0.773 
#>  upper_limit 0.853 
#>  conf_level  0.95  
#>  N           200   
#>  N_complete  200   
#>  J           6     

# Model implied denominator with its closed-form robust ML interval.
reliability_omega(data = items, denominator = "model_implied")
#>  term        value 
#>  estimate    0.82  
#>  se          0.0191
#>  lower_limit 0.782 
#>  upper_limit 0.857 
#>  conf_level  0.95  
#>  N           200   
#>  N_complete  200   
#>  J           6     

# ML CI from a covariance matrix (model implied denominator).
reliability_omega(S = cov(items), N = 200,
                  denominator = "model_implied", ci_method = "ml")
#>  term        value 
#>  estimate    0.82  
#>  se          0.0197
#>  lower_limit 0.781 
#>  upper_limit 0.858 
#>  conf_level  0.95  
#>  N           200   
#>  N_complete  200   
#>  J           6     

# Full information maximum likelihood with an auxiliary variable:
# missingness on y2 depends on an auxiliary z (missing at random
# given z). Supplying aux implies missing = "fiml".
z <- eta + rnorm(200, sd = 0.5)
d <- data.frame(items, z = z)
d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
reliability_omega(data = d, aux = "z", denominator = "model_implied")
#>  term        value 
#>  estimate    0.819 
#>  se          0.0202
#>  lower_limit 0.779 
#>  upper_limit 0.859 
#>  conf_level  0.95  
#>  N           200   
#>  N_complete  124   
#>  J           6     
# }
```
