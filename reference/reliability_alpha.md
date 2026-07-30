# Coefficient Alpha With a Confidence Interval

Estimates coefficient \\\alpha\\ for a homogeneous composite score and
returns a confidence interval for the population coefficient using one
of several documented methods.

## Usage

``` r
reliability_alpha(
  data = NULL,
  S = NULL,
  N = NULL,
  estimator = c("analytic", "model_implied"),
  missing = c("listwise", "fiml"),
  aux = NULL,
  ci_method = NULL,
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

  A symmetric covariance matrix among the items. If `S` is supplied, `N`
  must also be supplied; raw-data-only confidence interval methods
  (`"adf"`, `"bootstrap_*"`, `"percentile"`, `"bca"`) are then
  unavailable, as are `missing = "fiml"` and `aux`, since a covariance
  matrix has no incomplete cases.

- N:

  Total sample size; required when `S` (rather than `data`) is supplied.

- estimator:

  How \\\alpha\\ is estimated: `"analytic"` (the default) applies the
  closed-form equation to the observed covariance matrix, and
  `"model_implied"` takes the reliability implied by the
  \\\tau\\-equivalent single-factor model fit by maximum likelihood. See
  *Details*.

- missing:

  How incomplete rows of `data` are handled: `"listwise"` (the default;
  complete-case analysis, the historical behavior) or `"fiml"` (full
  information maximum likelihood, using every case with at least one
  observed item). See *Details*.

- aux:

  Optional character vector naming auxiliary variable columns of `data`,
  entered as saturated correlates under full information maximum
  likelihood. Supplying `aux` implies `missing = "fiml"`. See *Details*.

- ci_method:

  Method for constructing the confidence interval. See *Details* for
  which methods each estimator supports. The default `NULL` resolves to
  `"bonett"` for the analytic estimator and to `"mlr"` for the model
  implied estimator (or `"ml"` when only a covariance matrix is
  supplied, since `"mlr"` needs raw data). With `missing = "fiml"` the
  analytic default is `"ml"`, whose standard error comes from the FIML
  information matrix.

- conf_level:

  Confidence level for the interval (1 - Type I error rate). Defaults to
  `0.95`.

- B:

  Number of bootstrap replications when a bootstrap method is selected.
  Defaults to `10000`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, which leaves the user's current RNG state intact; supply an
  integer for reproducibility. When set, the function saves and restores
  `.Random.seed` so the user's global RNG state is not polluted.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample coefficient \\\alpha\\), `"se"` (standard error, `NA` for
methods that do not produce one), `"lower_limit"` and `"upper_limit"`
(clamped to \[0, 1\]), `"conf_level"`, `"N"` (the cases the analysis
used: the complete cases under listwise deletion, every case with at
least one observed item under `missing = "fiml"`), `"N_complete"` (the
complete cases, so the cost of listwise deletion is visible at a glance;
equal to `"N"` under listwise deletion and for covariance input), and
`"J"` (number of items). The selected coefficient, CI method, and
missing-data treatment travel as the attributes `coefficient`,
`ci_method`, `missing`, and (when supplied) `aux`; bootstrap calls also
record `B`.

## Details

Coefficient \\\alpha\\ was first derived by Guttman (1945) as his
\\\lambda_3\\ and was subsequently popularized by Cronbach (1951), under
whose name the coefficient is often informally cited. The attribution to
Cronbach is historically incomplete; the modern literature increasingly
refers to the coefficient simply as “coefficient alpha” (see, e.g.,
Sijtsma, 2009; Revelle & Zinbarg, 2009). DMAR follows that convention.

For a \\J\\-item composite \\Y = \sum_j X_j\\ the population coefficient
is \$\$\alpha = \frac{J}{J-1}\left(1 - \frac{\sum\_{j}
\sigma\_{j}^{2}}{\sigma\_{Y}^{2}}\right),\$\$ where \\\sigma\_{j}^{2}\\
is the variance of item *j* and \\\sigma\_{Y}^{2}\\ is the variance of
the composite. The sample estimate substitutes sample variances. Under
classical test theory, \\\alpha\\ equals the population reliability of
the composite when the items are essentially \\\tau\\-equivalent (i.e.,
equal factor loadings); when loadings differ, \\\alpha\\ is a lower
bound on the population reliability. For well-behaved homogeneous
measurement instruments coefficient \\\alpha\\ and McDonald's (1999)
coefficient \\\omega\\
([`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md))
typically yield very similar values; \\\omega\\ extends to congeneric
items (heterogeneous loadings) without the lower-bound caveat. For
ordered-categorical items see
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md).

**Two estimators of the same coefficient.** The argument `estimator`
selects how \\\alpha\\ is estimated from the data. Both target the same
population quantity, and they agree in the population whenever the
\\\tau\\-equivalent model holds; they differ in a finite sample because
they take different routes to it.

- `"analytic"`:

  The default. The closed-form equation above, applied to the observed
  covariance matrix. This is the classical coefficient, the number a
  hand calculation produces, and it makes no assumption beyond those of
  classical test theory.

- `"model_implied"`:

  The reliability implied by the \\\tau\\-equivalent (equal loadings)
  single-factor model fit by maximum likelihood. With a shared loading
  \\\lambda\\ and error variances \\\psi_j^{2}\\ the model implied
  reliability of the \\J\\-item composite is \$\$\alpha = \frac{(J
  \lambda)^{2}}{(J \lambda)^{2} + \sum\_{j} \psi\_{j}^{2}},\$\$ with
  maximum likelihood estimates substituted. Estimating through the model
  brings inference the formula cannot provide: a delta method standard
  error, a robust (Satorra-Bentler) variant under nonnormality, the
  profile likelihood interval, and a fit assessment of the
  \\\tau\\-equivalence claim itself. When the equal-loadings claim is
  doubtful, the congeneric
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  is the appropriate coefficient rather than either \\\alpha\\.

Users of MBESS will recognize these as `ci.reliability(type = "alpha")`
and `ci.reliability(type = "alpha-cfa")` respectively.

**Missing data and auxiliary variables.** By default incomplete rows are
listwise-deleted, which is unbiased only when the data are missing
completely at random and is inefficient always; the `mlmr` vignette
develops the argument at length. Setting `missing = "fiml"` keeps every
case with at least one observed item and estimates by full information
maximum likelihood, which is consistent and efficient under the weaker
missing at random (MAR) assumption. The `aux` argument names auxiliary
variables: columns of `data` that are not part of the composite but are
correlated with the items or with the reasons values are missing. They
are entered as *saturated correlates* (Graham, 2003): correlated freely
with each other and with every item's residual, never loading on the
factor and never entering the composite, so the measurement model is
undisturbed while FIML uses their information. Beyond recovering
information, a good auxiliary makes the MAR assumption itself more
plausible, since missingness that depends on the auxiliary becomes MAR
once the auxiliary is conditioned on (Collins, Schafer, & Kam, 2001).
Supplying `aux` implies `missing = "fiml"`; combining it with an
explicit `missing = "listwise"` is an error. Listwise deletion remains
the default so no existing result changes and the missing-data treatment
is always a visible, deliberate choice. How each estimator uses FIML:
the model implied estimator simply fits its model with `missing = "ml"`;
the analytic estimator applies the classical formula, unchanged, to the
FIML estimate of the item covariance matrix (from a saturated model over
the items and any auxiliaries), so the estimand stays the classical
coefficient and only the covariance matrix it is computed from improves.
Under `missing = "fiml"` the available intervals are `"ml"` and
`"ml_logistic"` (both estimators; the standard error comes from the FIML
information matrix), `"mlr"` and `"mlr_logistic"` (model implied
estimator; the Yuan-Bentler robust standard error), and the bootstrap
methods, which resample rows (including the partially observed ones) and
refit by FIML on each replication. The complete-data closed forms
(`"feldt"`, `"fisher"`, `"bonett"`, `"hakstian_whalen"`), `"adf"`, and
`"likelihood"` are errors with `missing = "fiml"` rather than silently
reverting to listwise deletion. Multiple imputation is a different
feature with a different interface and is out of scope here.

Available confidence interval methods (set via `ci_method`). Some belong
to one estimator only, because a closed-form interval for the sample
coefficient and a model-based interval are not interchangeable;
requesting a method the chosen estimator cannot supply is an error that
names the estimator to use instead:

- `"feldt"`:

  *Analytic estimator only.* The *F*-distribution interval of Feldt
  (1965), exact under multivariate normality and parallel items.

- `"fisher"`:

  *Analytic estimator only.* Fisher's \\z'\\ transformation (Fisher,
  1950). Tends to overcover (Padilla, Divers, & Newton, 2012) and is
  generally not recommended.

- `"bonett"`:

  *Analytic estimator only.* Bonett's (2002) log transformation. The
  default for that estimator; well-behaved under normality.

- `"hakstian_whalen"`:

  *Analytic estimator only.* Cube-root transformation of Hakstian and
  Whalen (1976).

- `"mlr"`, `"mlr_logistic"`:

  *Model implied estimator only.* Wald interval using the robust
  (Satorra & Bentler, 1994) standard error from the fitted model. The
  default for that estimator; recommended among the closed forms when
  item distributions deviate from normality (Kelley & Pornprasertmanit,
  2016). Requires raw data.

- `"likelihood"`:

  *Model implied estimator only.* Profile likelihood interval: the set
  of population values not rejected by the likelihood ratio test under
  the \\\tau\\-equivalent model, located by refitting under the
  nonlinear constraint that the model implied reliability equals each
  candidate value. Respects \[0, 1\], is not forced to be symmetric, and
  works from raw data or covariance input. Requires lavaan. It is
  unavailable with the analytic estimator because the interval and the
  point estimate would then refer to different quantities: the interval
  profiles the model implied coefficient while the estimate is the
  sample coefficient, so under a misspecified model the interval can
  exclude the estimate it accompanies.

- `"ml"`, `"ml_logistic"`:

  Available to both estimators, by the route each affords. With
  `"analytic"` it is the closed-form ML standard error of van Zyl,
  Neudecker, and Nel (2000), computed directly from the covariance
  matrix; with `"model_implied"` it is the delta method standard error
  from the fitted model. The `_logistic` variant applies Browne's (1982)
  logit transformation.

- `"adf"`, `"adf_logistic"`:

  Available to both estimators. With `"analytic"` it is the asymptotic
  distribution-free standard error of Maydeu-Olivares, Coffman, and
  Hartmann (2007); with `"model_implied"` the model is fit by weighted
  least squares (Browne, 1984) and the delta method applied. Requires
  raw data and a relatively large sample size.

- `"bootstrap_se"`, `"bootstrap_se_logistic"`, `"percentile"`, `"bca"`:

  Available to both estimators. Nonparametric bootstrap intervals (Efron
  & Tibshirani, 1993): the rows of `data` are resampled with replacement
  `B` times and the chosen estimator is recomputed on each replication.
  `"percentile"` takes the interval limits from the empirical quantiles
  of the bootstrap estimates (for a 95 percent interval, the 2.5th and
  97.5th percentiles); it respects \[0, 1\] and is not forced to be
  symmetric about the estimate, but its coverage degrades when the
  estimator is biased or its variance changes with the parameter.
  `"bca"` (bias-corrected and accelerated) adjusts the two quantile
  positions for exactly those features, estimating the median bias from
  the bootstrap distribution and the acceleration from the jackknife,
  and is second-order accurate where the percentile interval is
  first-order accurate (DiCiccio & Efron, 1996). `"bootstrap_se"` uses
  the standard deviation of the bootstrap estimates as a standard error
  in an ordinary normal-theory interval; the `_logistic` variant builds
  that interval on the logit scale so the endpoints respect \[0, 1\].
  Replications on which the estimator cannot be computed (for example, a
  model refit that does not converge) are dropped, and the interval is
  computed from the replications that return a value. Requires raw data
  and the boot package. The default `B = 10000` is an accuracy choice:
  the BCa adjustment pushes the working quantiles farther into the tails
  of the bootstrap distribution than the percentile interval uses, and
  stabilizing them takes more replications than the customary 2000;
  reduce `B` for exploration, not for a reported analysis. Bootstrap
  results vary from run to run; supply `seed` for reproducibility.

- `"none"`:

  Return only the point estimate.

**Comparison with other packages.** The psych package provides
[`alpha`](https://rdrr.io/pkg/psych/man/alpha.html), which reports
coefficient \\\alpha\\ along with item-level diagnostics (item-total
correlations, alpha-if-item-deleted, and several alternative
coefficients). The emphasis in psych is broad exploratory psychometric
reporting. `reliability_alpha` in DMAR differs in emphasis: it returns a
single point estimate alongside a principled confidence interval drawn
from the methods compared in Kelley and Pornprasertmanit (2016).

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

Fisher, R. A. (1950). *Statistical methods for research workers*.
Edinburgh, UK: Oliver & Boyd.

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

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
standard errors in covariance structure analysis. In A. von Eye & C. C.
Clogg (Eds.), *Latent variables analysis: Applications for developmental
research* (pp. 399–419). Thousand Oaks, CA: Sage.

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
*Sociological Methodology, 30*, 165–200.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Maydeu-Olivares, A., Coffman, D. L., & Hartmann, W. M. (2007).
Asymptotically distribution-free (ADF) interval estimation of
coefficient alpha. *Psychological Methods, 12*, 157–176.
[doi:10.1037/1082-989X.12.2.157](https://doi.org/10.1037/1082-989X.12.2.157)

McDonald, R. P. (1999). *Test theory: A unified treatment*. Mahwah, NJ:
Lawrence Erlbaum Associates.

Padilla, M. A., Divers, J., & Newton, M. (2012). Coefficient alpha
bootstrap confidence interval under nonnormality. *Applied Psychological
Measurement, 36*, 331–348.
[doi:10.1177/0146621612445470](https://doi.org/10.1177/0146621612445470)

Revelle, W., & Zinbarg, R. E. (2009). Coefficients alpha, beta, omega,
and the GLB: Comments on Sijtsma. *Psychometrika, 74*, 145–154.
[doi:10.1007/s11336-008-9102-z](https://doi.org/10.1007/s11336-008-9102-z)

Sijtsma, K. (2009). On the use, the misuse, and the very limited
usefulness of Cronbach's alpha. *Psychometrika, 74*, 107–120.
[doi:10.1007/s11336-008-9101-0](https://doi.org/10.1007/s11336-008-9101-0)

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution
of the maximum likelihood estimator of Cronbach's alpha. *Psychometrika,
65*(3), 271–280.
[doi:10.1007/BF02296146](https://doi.org/10.1007/BF02296146)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper that dispatches by coefficient),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
(single-factor CFA used internally),
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
[`alpha`](https://rdrr.io/pkg/psych/man/alpha.html).

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Simulate six tau-equivalent items with population reliability ~ .8.
J <- 6
loadings <- rep(0.6, J)
eta <- rnorm(200)
errors <- matrix(rnorm(200 * J, sd = sqrt(1 - 0.6^2)), 200, J)
items <- outer(eta, loadings) + errors
colnames(items) <- paste0("y", seq_len(J))

# Default (Bonett's transformation) CI from raw data.
reliability_alpha(data = items)
#>  term        value
#>  estimate    0.777
#>  se          0.11 
#>  lower_limit 0.723
#>  upper_limit 0.82 
#>  conf_level  0.95 
#>  N           200  
#>  N_complete  200  
#>  J           6    

# Same point estimate from a covariance matrix; CI requires N.
S <- cov(items)
reliability_alpha(S = S, N = 200, ci_method = "feldt")
#>  term        value
#>  estimate    0.777
#>  se          <NA> 
#>  lower_limit 0.725
#>  upper_limit 0.821
#>  conf_level  0.95 
#>  N           200  
#>  N_complete  200  
#>  J           6    

# Percentile bootstrap with 200 reps for a quick example.
reliability_alpha(data = items, ci_method = "percentile", B = 200)
#>  term        value 
#>  estimate    0.777 
#>  se          0.0248
#>  lower_limit 0.72  
#>  upper_limit 0.817 
#>  conf_level  0.95  
#>  N           200   
#>  N_complete  200   
#>  J           6     

# Full information maximum likelihood with an auxiliary variable:
# missingness on y2 depends on an auxiliary z (missing at random
# given z), so listwise deletion is biased and FIML with z is not.
# Supplying aux implies missing = "fiml".
z <- eta + rnorm(200, sd = 0.5)
d <- data.frame(items, z = z)
d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
reliability_alpha(data = d, aux = "z")
#>  term        value 
#>  estimate    0.778 
#>  se          0.0256
#>  lower_limit 0.728 
#>  upper_limit 0.829 
#>  conf_level  0.95  
#>  N           200   
#>  N_complete  124   
#>  J           6     
```
