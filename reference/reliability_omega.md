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
  ci_method = c("mlr_logistic", "ml_logistic", "mlr", "ml", "likelihood", "adf",
    "adf_logistic", "feldt", "fisher", "bonett", "hakstian_whalen", "bootstrap_se",
    "bootstrap_se_logistic", "percentile", "bca", "none"),
  denominator = c("observed", "model_implied"),
  missing = c("fiml", "listwise", "two_stage"),
  aux = NULL,
  conf_level = 0.95,
  B = 10000,
  seed = NULL,
  ncpus = 1L,
  improper = c("keep", "drop"),
  wald_df = Inf,
  progress = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of item scores (rows are respondents,
  columns are items). Either `data` or `S` must be supplied. How
  incomplete rows are handled is governed by `missing`: full information
  maximum likelihood by default, listwise deletion with
  `missing = "listwise"`, or two-stage maximum likelihood with
  `missing = "two_stage"`. When `aux` is supplied, `data` contains both
  the item columns and the auxiliary columns, and the columns not named
  in `aux` are the items. An item column that is not numeric (or
  logical), and an item that takes a single value in the rows analyzed
  (the complete rows under listwise deletion, the rows that observe it
  otherwise), are refused by name before anything is computed, as is a
  zero variance on the diagonal of `S`: such an item adds a constant to
  the composite, whose reliability it leaves as it is, while the
  coefficient would count it among the \\J\\ items (alpha with it is
  \\J(J-2)/(J-1)^2\\ times alpha without it). A composite that does not
  vary is refused too, as when the items are ipsative (their scores sum
  to the same total in every row): its reliability is not defined. Under
  `missing = "fiml"` and `"two_stage"` it is checked in the rows that
  observe every item when they outnumber the items, where a composite
  that does not vary leaves the likelihood of the saturated model
  without a maximum, and otherwise on the EM estimate of the covariance
  matrix. In `holzinger_swineford`, `t1_visual_perception`, `t2_cubes`,
  and `t4_lozenges`, each centered on the pupil's mean of the three, are
  ipsative items of this kind. SPSS user-defined missing values, which
  `haven::read_sav(user_na = TRUE)` keeps as values of a
  `haven_labelled_spss` column, are missing, as SPSS treats them.

- S:

  A symmetric covariance matrix among the items. If supplied, `N` must
  also be supplied; methods that require raw data (`"mlr*"`, `"adf*"`,
  `"bootstrap_*"`, `"percentile"`, `"bca"`) are then unavailable, as are
  `missing = "fiml"` and `missing = "two_stage"` named explicitly, and
  `aux`, since a covariance matrix has no incomplete cases (the default
  of `missing` analyzes it as it stands).

- N:

  Total sample size; required when `S` is supplied.

- ci_method:

  Method for constructing the confidence interval. See *Details*. When
  not supplied (or passed as `NULL`): for
  `denominator = "model_implied"` the default is `"mlr_logistic"` with
  raw data, the robust maximum likelihood standard error with the limits
  formed on the logit scale and transformed back (Browne, 1982; Kelley,
  2007a), whatever the treatment of missing values, except under
  `wald_df = "missing_information"`, where it is `"mlr"`, the same
  standard error with the limits formed on the scale of omega, the
  interval that rule is calibrated for; with covariance input, which the
  robust standard error cannot use, it is `"ml_logistic"`, the maximum
  likelihood standard error with the limits formed the same way (`"ml"`
  under `wald_df = "missing_information"`); for
  `denominator = "observed"` (robust omega) no interval is computed by
  default, because its interval is bootstrap based and a bootstrap is
  never run unless requested. Ask for `"percentile"` or `"bca"` to
  obtain the recommended interval. The closed forms, the Wald intervals,
  and the profile likelihood apply to `denominator = "model_implied"`
  only; a combination that does not apply stops with a message naming
  the denominator that does (see *Details*).

- denominator:

  How the total variance in the denominator of \\\omega\\ is estimated:
  `"observed"` (default; robust omega) uses the variance of the
  composite estimated directly from the data; `"model_implied"` uses the
  total variance reproduced by the fitted single-factor model. See
  *Details*. With `"observed"`, `ci_method` must be a bootstrap method
  or `"none"`.

- missing:

  How incomplete rows of `data` are handled: `"fiml"` (the default; full
  information maximum likelihood, using every case with at least one
  observed item or observed auxiliary variable), `"listwise"`
  (complete-case analysis), or `"two_stage"` (two-stage maximum
  likelihood, using the same cases as FIML). The default analysis is
  FIML whenever an item value is missing and the complete data analysis
  otherwise, unlike
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  which fits complete data by FIML too: data in which no value is
  missing, data complete but for rows with every item missing (which
  both treatments drop), and `S` are analyzed as they stand, since the
  complete-data intervals need that analysis; incomplete items are
  analyzed by FIML, or, when `ci_method` names an interval FIML cannot
  give, the call stops and names the choices. See *Details*. With
  `missing = "listwise"`, when fewer than two rows are complete,
  listwise deletion leaves nothing to analyze, and the function stops
  before any fit, says so, and names the treatment that uses the
  incomplete rows. The single-factor model is fit by maximum likelihood,
  which needs the covariance matrix positive definite, and \\N\\
  complete rows give a matrix of rank at most \\N - 1\\, so with no more
  complete rows than items the function stops the same way;
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  which needs only the sums of the matrix, is computed from two or more
  complete rows. Full information maximum likelihood does not impose an
  ordinary complete-case count rule. Observed-denominator omega and
  two-stage estimation require unrestricted moments, checked for affine
  collapse or numerical instability as described in
  [`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md).
  EM convergence alone does not establish existence or uniqueness.
  Model-implied omega can retain a proper identified constrained FIML
  fit and H0-based uncertainty while withholding unavailable H1
  diagnostics and residual correlations. The result states this
  qualification. Bootstrap resamples retain a distinct-case estimability
  screen as a separate policy, not a universal likelihood-existence
  theorem. Listwise complete rows are counted as distinct rows. Two
  columns of one name are refused before any of this, under every
  treatment, since the items are read by name.

- aux:

  Optional character vector naming auxiliary variable columns of `data`:
  entered as saturated correlates under full information maximum
  likelihood, and into the first stage under `missing = "two_stage"`.
  Supplying `aux` without `missing` implies `missing = "fiml"`. See
  *Details*.

- conf_level:

  Confidence level for the interval. Defaults to `0.95`.

- B:

  Number of bootstrap replications when a bootstrap method is selected,
  a single whole number. Defaults to `10000`. The fewest accepted
  depends on the interval method and the level: at the default
  `conf_level`, 400 for `"percentile"` (10 replications beyond each
  limit), 2000 for `"bca"` (50 beyond each nominal limit, because the
  adjusted limits sit farther into the tails), and 100 for
  `"bootstrap_se"` and `"bootstrap_se_logistic"` (a standard deviation
  of the replications); at 0.99 the first two are 2000 and 10000. A
  smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`). The count is held to the number requested.
  Replications set aside afterward are counted and reported in a warning
  of class `dmar_bootstrap_replications_dropped`, which names the
  function and says to request more when those that return a value are
  fewer than the number required; when fewer than half of the `B`
  replications return a value, the interval is withheld (limits and
  bootstrap standard error `NA`) under a warning of class
  `dmar_bootstrap_insufficient`. See
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md),
  which also describes `options(dmar.bootstrap_exploratory = TRUE)` for
  a first look with a smaller `B`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, under which the resamples are drawn from, and advance, the
  current random number stream; supply an integer for a reproducible
  interval, which is set for the call only, the caller's random number
  state being restored when it returns.

- ncpus:

  The number of processes over which the bootstrap replications are
  spread. The default `1` runs them in this R session; a larger number
  uses forked workers (boot's `parallel = "multicore"`) or, on Windows,
  which cannot fork, a socket cluster (`"snow"`), whose workers load the
  installed DMAR. boot draws every resample before any replication runs,
  so a given `seed` returns the same interval, counts, and replicates
  whatever `ncpus` is. A session that forks workers and also starts
  background processes through callr or processx needs the environment
  variable `PROCESSX_NOTIFY_OLD_SIGCHLD` set before processx is loaded
  (see
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)).
  Ignored unless a bootstrap method is selected.

- improper:

  How a bootstrap replication whose single-factor solution is improper
  (a negative error variance, a Heywood case) is treated when the
  solution still yields a coefficient: `"keep"` (the default) uses its
  estimate like any other, and `"drop"` sets it aside. Either way such
  replications are counted in the `bootstrap` attribute and the printed
  line, and a replication whose solution yields no coefficient is set
  aside under both (see *Details*). With `"drop"` an improper point
  estimate stops the call, since an interval built from proper
  replications only would not describe it. Ignored unless a bootstrap
  method is selected.

- wald_df:

  The reference distribution of a Wald interval: `Inf` (the default),
  the normal distribution; a single positive number, the degrees of
  freedom of a *t* distribution; or `"missing_information"`, von
  Hippel's (2016) degrees of freedom \\\max\\3, N(1 - \gamma) - 1\\\\
  from omega's fraction of missing information \\\gamma\\ (\\N - 1\\
  with complete data). A number applies to the Wald intervals of the
  model implied coefficient (`"ml"`, `"mlr"`, `"adf"`, and their
  `_logistic` variants), the keyword to `"ml"` and `"mlr"` only, without
  two-stage estimation, the intervals its calibration covers; anything
  else is an error that says why. See *Details*.

- progress:

  `NULL` (the default), or a function of two arguments that the
  bootstrap calls as it runs, so that an application can show how far a
  long bootstrap has gone and let its user stop it (a Shiny app, for
  example). It is called as `progress(done, total)`, with the number of
  replications attempted so far and `B`: after each replication when
  `ncpus = 1`, and with more processes after each block of replications
  they complete together, about a fiftieth of `B`
  (`ncpus * ceiling(B / (50 * ncpus))` replications, the last block
  possibly fewer), so that `done` advances by the block. A BCa interval
  takes its acceleration from a regression on the resampling frequencies
  when at least as many replications return a value as there are cases,
  and from the jackknife, one leave-one-out refit per case, when fewer
  do or when that regression fails or gives influence values that are
  not all finite; the jackknife runs at most once in a call, and not at
  all when no interval is formed or every returned replicate has the
  same value. While it runs, the callback is called again after each
  refit (with more processes, after each block of
  `ncpus * ceiling(n / (50 * ncpus))` refits, `n` the number of cases
  resampled), with `done` and `total` both `B`, so it can be called more
  than `B` times, and `done` equal to `total` does not mean that the
  bootstrap has finished. If the callback returns `FALSE`, the bootstrap
  stops before the next replication or block with an error of class
  `dmar_bootstrap_cancelled`, whose elements `done` and `total` say how
  far it got (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md));
  any other value is ignored. An error inside the callback stops the
  call with an error of class `dmar_progress_error`, whose message names
  the callback and which carries `done`, `total`, and, as `parent`, the
  callback's own condition with its class intact, so that an application
  can tell a failing callback from a failing analysis. The callback runs
  outside the random number stream, so a given `seed` returns the same
  result with or without it, and a callback that draws random numbers
  leaves the caller's stream as it found it. With `NULL` nothing is
  called. Ignored unless a bootstrap method is selected.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample coefficient \\\omega\\), `"se"` (standard error on the
coefficient scale, `NA` for methods that do not produce one; for the
transformation-based intervals `"fisher"`, `"bonett"`, and
`"hakstian_whalen"` it is the delta method back-transform of the
transformation-scale standard error), `"se_transformed"` (only for those
transformation-based intervals: the standard error on the transformation
scale, with the scale named in the attribute `se_transform_scale`:
`"fisher_z"`, `"log(1-alpha)"`, or `"cube_root"`), `"ci_lower"` and
`"ci_upper"` (each clamped into \[0, 1\], so the limits stay in order),
`"conf_level"`, `"N"` (the cases the analysis used: the complete cases
under listwise deletion, and under `"fiml"` and `"two_stage"` every case
with at least one observed item or observed auxiliary variable),
`"N_complete"` (the complete cases; equal to `"N"` under listwise
deletion and for covariance input), and `"n_items"` (number of items).
Attributes `coefficient` (`"omega"`), `ci_method`, `denominator`,
`missing`, `wald_df` (the degrees of freedom of the Wald interval's *t*
reference, `Inf` for the normal one and for intervals that are not Wald
intervals; under `wald_df = "missing_information"` the attribute
`wald_df_table` holds the
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
table they came from), and (when supplied) `aux` record the computation;
bootstrap calls also record `B`. The attribute `residual_cor` carries
the local dependence summary of a maximum likelihood fit of the
single-factor model to the same cases, whatever the `ci_method` (the
weighted least squares fit behind the ADF intervals has residuals of its
own, so it is not used for the summary): a list with `matrix` (the
observed minus the model implied item correlations, named by item),
`largest` (the residual correlation largest in absolute value, with its
sign), `items` (the pair it belongs to), and `saturated`. With three
items the single-factor model is saturated (zero degrees of freedom) and
a proper solution reproduces every correlation, so its residuals are
zero and carry no information about local dependence: `saturated` is
then `TRUE`, `largest` is `NA`, and `items` is `NULL`, rather than
naming whichever pair carries the largest rounding error. The printed
table names the denominator and the largest residual correlation (or
says that the model is saturated) in two lines under the table. A
bootstrap call also carries the attribute `bootstrap`, a list recording
the replications the interval conditions on: `attempted` (`B`),
`returned` (the replications whose value the interval uses; the standard
error is their standard deviation and the percentile and BCa limits are
their quantiles), `failed` (the difference), `reasons` (a named count of
why replications were set aside: `error` when the refit raised an error,
`non_finite` when it was set aside without one, as a fit that did not
converge is, and as an improper one is under `improper = "drop"` even
though its estimate is finite), `error_messages` (the distinct error
messages with their counts, empty when none), `undetermined` (the
replications set aside because the resample left a polychoric
correlation undetermined, which only
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
counts; zero here), `improper` and `improper_kept` (the replications
whose single-factor solution was improper, with a negative error
variance, or with three items had no single-factor solution at all:
those set aside, which are among the `non_finite` ones, and those whose
estimate the interval uses, which are among the `returned` ones),
`no_solution` (the replications of three items set aside because their
correlations admit no single-factor solution, a subset of `improper`),
`improper_rule` (`"keep"` or `"drop"`; see the argument `improper`),
`not_converged` (zero here: a fit that did not converge is counted among
the `non_finite` ones), `too_few_cases` (the replications set aside
because the resample held too few distinct cases for the rule on the
number of cases, which are among the `error` ones; when they number more
than half of `B`, no interval is formed), `caution` (`TRUE` when the
replications set aside and the improper ones kept number at least \\(1 -
C)/2\\ of `B`, with \\C\\ the confidence level, so that how they were
treated can decide where a limit falls), `require_half` (`TRUE`: an
interval requires at least half of all attempted replications to return
defined estimates; this is an estimability safeguard, not a coverage
guarantee), `required` (the fewest replications a request for an
interval of this method and level may ask for, a count held to `B` and
not to the replications that return a value: 400 for `"percentile"`,
2000 for `"bca"`, 100 for the two bootstrap standard error methods at
the 95 percent level; see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)),
`bootstrap_status` (`"complete"`, `"dropped"`, `"exploratory"`,
`"insufficient_count"`, or `"insufficient_share"`; under the last two
the limits are `NA`), `B_used` (the smallest count any reported limit
rests on), `target_returned` and `target_withheld` (counts and
withholding flags for each target in a vector bootstrap), `acceleration`
and `acceleration_method` (the BCa acceleration used and whether it came
from the regression or the jackknife; `NA` unless a BCa interval was
formed), `jackknife` (the same `attempted`, `returned`, `failed`,
`reasons`, `error_messages`, `undetermined`, `improper`,
`improper_kept`, `no_solution`, `not_converged`, and `too_few_cases` for
the leave-one-out refits behind the acceleration, all zero when no
jackknife ran), `ncpus` (the number of processes used), `t` (the
estimate of every replication, in order, `NA` where there is none; an
improper estimate set aside under `improper = "drop"` is kept here, so
the interval under the other rule can be formed without refitting), and
`status` (a factor giving the fate of every replication, with levels
`"proper"`, `"improper"`, `"error"`, `"non_finite"`, `"undetermined"`,
and `"not_converged"`). The interval uses the replications whose status
is `"proper"` and, under `improper = "keep"`, those whose status is
`"improper"` and whose `t` is finite. The attribute `B_used` is the
number of replications the interval uses (`returned`). The attributes
`n_dropped` (`B` minus `B_used`), `B_required`, and `bootstrap_status`
repeat the counts in the names every bootstrap result of the package
carries. The warning has class `dmar_bootstrap_replications_dropped`
when replications were set aside and the interval was formed, and
`dmar_bootstrap_insufficient` when the interval was withheld. A BCa
interval is also withheld, with the second class, when fewer than 10
replications lie beyond a limit where the bias correction and the
acceleration place it (see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)).
The printed table shows these counts in one line, and a warning states
them whenever a replication was set aside or the caution applies.

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

Stated as one identity, the two denominators differ by
\\\mathbf{1}'(\mathbf{S} - \hat{\Sigma})\mathbf{1}\\, the sum of every
residual covariance of the single-factor fit. The model implied setting
treats covariance the model omits as absent from the composite; the
observed setting counts it as variance of the composite that the factor
does not explain. The first describes the composite as the fitted model
describes it, the second the composite users compute. Estimating the
factor model from its own implied moments and reporting the reliability
of a scored composite with the observed variance are therefore not in
tension: the model is needed for the numerator, which partitions the
variance of the composite into the part due to the factor and the rest,
and is not needed for the denominator, which is the variance of the
composite itself. The omega that
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) reports
beside the loadings, and that
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
reports for continuous items, is the model implied coefficient, a
property of the fitted model.

Two cautions frame the choice. First, no denominator repairs a
misspecified measurement model: the fitted loadings absorb part of
whatever structure the single-factor model omits, so the numerator is
affected under either setting. Assess the single-factor model (for
example with
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)) before
interpreting any \\\omega\\ variant, and model real multidimensionality
directly rather than patching over it. A gap between the two
denominators shows that the model does not reproduce the composite
variance, but a small gap does not show that the loadings are free of
such structure: a doublet or a method factor shared by a few same-signed
items inflates their loadings while leaving the two denominators close.
The residual correlations of the single-factor fit are where that
structure usually shows, so the result carries them (attribute
`residual_cor`; see *Value*) and the printed table names the largest.
One limit applies: when the shared structure sits on the items that
dominate the factor, the single factor absorbs it and the misfit appears
as residual correlations among the other items, so the largest residual
signals local dependence more reliably than it locates it. Second, the
naming history is worth knowing. The observed denominator coefficient
was introduced as hierarchical omega (Kelley & Pornprasertmanit, 2016;
MBESS type `"hierarchical"`), a name motivated by a hierarchical factor
logic: model misfit is viewed as a set of minor common factors (visible
as residual correlations), the single factor is retained as an
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
principle with categorical omega as Kelley and Pornprasertmanit (2016)
computed it: in both, the total variance in the denominator is not taken
from the fitted factor model. Green and Yang (2009) built the
denominator of their categorical coefficient from the model implied
polychoric correlations instead, and
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
offers both through its own `denominator` argument. The reliability
vignette develops this framing.

Available confidence interval methods (set via `ci_method`):

- `"ml"`, `"ml_logistic"`:

  Wald interval using the maximum likelihood standard error from lavaan
  (Raykov, 2002). The `_logistic` variant applies Browne's (1982) logit
  transformation (see below).

- `"mlr"`, `"mlr_logistic"`:

  Wald interval using the robust (Satorra & Bentler, 1994) standard
  error, recommended among closed-form methods when item distributions
  deviate from normality (Kelley & Pornprasertmanit, 2016).
  `"mlr_logistic"` is the default for `denominator = "model_implied"`.

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
  (DiCiccio & Efron, 1996). The acceleration is estimated as boot
  estimates it, by regression on the resampling frequencies when at
  least as many replications as cases returned a value and by the
  jackknife otherwise, or when that regression fails; a leave-one-out
  refit that fails is counted and set aside rather than stopping the
  call, and when the acceleration or the bias correction cannot be
  estimated at all the percentile interval is returned in its place,
  with a warning and the substitution recorded in the `ci_method`
  attribute; see
  [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  for the full account. `"bootstrap_se"` uses the standard deviation of
  the bootstrap estimates as a standard error in a normal-theory
  interval, built on the logit scale for the `_logistic` variant so the
  endpoints respect \[0, 1\] (the standard deviation is carried to the
  logit scale by the delta method, divided by \\\hat\omega(1 -
  \hat\omega)\\, as in Equation 32 of Kelley and Pornprasertmanit, 2016,
  and in `MBESS::ci.reliability()`). When every replication returns the
  same value, the bootstrap distribution is a single point and every
  method returns that point at both limits, with a warning (items that
  are exact linear functions of one variable, which give alpha such
  replications, do not reach the bootstrap here: their covariance matrix
  is singular, and the single-factor fit behind the estimate stops on
  it). Replications whose model refit does not converge are set aside,
  and so is a replication in which an item takes a single value (a
  resample that has left out every case on which a rare item differs
  from the rest), which is refused before the refit, by name, as the
  data would be, and so is one with no more distinct cases than the rule
  on the number of cases for the data asks for (see `missing`; a
  resample repeats rows, so near the rule's bound many are, and when
  more than half are no interval is reported, since those left hold more
  distinct cases than a resample typically does and their spread
  understates the sampling variability, and at the bound each holds
  every distinct case; without repeated rows that happens at fewer cases
  than about 1.65 times the count the rule compares them with); the
  interval is computed from the replications that return a value, whose
  quantiles the percentile and BCa limits are; the counts, and why the
  rest were set aside, travel in the `bootstrap` attribute (see
  *Value*), print under the table, and are stated in a warning whenever
  a replication was set aside. The default `B = 10000` is an accuracy
  choice: the BCa adjustment pushes the working quantiles farther into
  the tails of the bootstrap distribution than the percentile interval
  uses, and stabilizing them takes more replications than the customary
  2000; reduce `B` for exploration, not for a reported analysis.
  Recommended when assumptions of parametric methods are questionable.
  Each replication refits the single-factor model by maximum likelihood:
  robust maximum likelihood has the same estimates and differs only in
  the standard errors and the test statistic, which a replication does
  not use, so the robust fit is not repeated `B` times. Requires raw
  data and the boot package; supply `seed` for run-to-run
  reproducibility, and `ncpus` to spread the replications over several
  processes, which leaves the interval for a given seed unchanged. How
  replications with an improper solution are treated is described next.

- `"none"`:

  Return only the point estimate.

**Which intervals apply to which denominator.**

- `denominator = "model_implied"`:

  Every method on the menu: the Wald intervals `"ml"`, `"mlr"`, and
  `"adf"` with their `_logistic` variants, the profile `"likelihood"`
  interval, the alpha-derived closed forms `"feldt"`, `"fisher"`,
  `"bonett"`, and `"hakstian_whalen"`, the four bootstrap methods, and
  `"none"`.

- `denominator = "observed"` (robust omega, the default):

  The four bootstrap methods (`"percentile"`, `"bca"`, `"bootstrap_se"`,
  `"bootstrap_se_logistic"`) and `"none"` only. The delta method
  standard errors, the profile likelihood, and the alpha-derived closed
  forms are derived under the model implied ratio and do not account for
  sampling of the observed denominator, so `"mlr_logistic"`, the default
  interval of the model implied coefficient, is not an interval for
  robust omega.

A request for a combination that does not apply, such as
`ci_method = "mlr"` with the default denominator, stops with a message
naming the denominator that does; nothing falls back silently.
Covariance input (`S`) and the incomplete data methods
(`missing = "fiml"` and `"two_stage"`) narrow the menu further, as the
`S` argument and the paragraphs on missing data below describe, and each
stops the same way. The pairing of robust omega with the bootstrap is
not a limitation in practice, since the bootstrap is the interval Kelley
and Pornprasertmanit (2016) recommend for the observed denominator
coefficient in any case. Because no analysis in DMAR runs a bootstrap
unless the user requests one, the default for robust omega is the point
estimate with no interval, accompanied by a message naming the call that
produces the recommended interval; request `ci_method = "percentile"` or
`"bca"` to obtain it.

**Limits on the scale of omega or on the logit scale.** The Wald
intervals come in two forms, both built from the same standard error.
`"ml"`, `"mlr"`, and `"adf"` take the estimate plus or minus a normal
quantile times the standard error, an interval symmetric about the
estimate whose upper limit can pass one (the limits are then clamped
into \[0, 1\]; see *Value*). Their `_logistic` variants form the limits
on the logit scale, \\\log\\\hat\omega / (1 - \hat\omega)\\\\, whose
standard error the delta method gives as \\\mathrm{SE} / \\\hat\omega
(1 - \hat\omega)\\\\, and transform them back (Browne, 1982): the limits
of a monotone transformation of a parameter are the transformation of
its limits, the confidence interval transformation principle that Kelley
(2007a) sets out. The interval then lies within (0, 1) and is asymmetric
near one, as the sampling distribution of a high reliability is; the
`se` row is the standard error of \\\hat\omega\\ itself under either
form. The default is `"mlr_logistic"` on the evidence of a Monte Carlo
study (in the package's source repository,
`methodology/reliability_limits`) of factors of three to six items with
*N* from 50 to 500. The study formed both kinds of limits from the same
standard error, the normal theory maximum likelihood one for normal data
and the Huber-White one of MLR for nonnormal data, so what it compared
is the scale on which the limits are formed. With normal data, under
listwise deletion and FIML alike, the logit limits covered from .939 to
.956 in every cell (a mean absolute departure from .95 of .003), where
the limits on the scale of omega covered from .910 to .952 (.009 to
.010), their misses falling mostly on one side, with the population
value below the lower limit. With nonnormal data fitted by MLR both
forms covered below .95, the logit form by less (a mean absolute
departure of .013 against .022). No cell used the robust.sem standard
error that `"mlr_logistic"` takes with complete data; it estimates the
same asymptotic variance as the Huber-White one, which the default takes
under FIML.
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
forms the limits of its omega row the same way, and
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
those of *H* by default, so the three give limits of one form. With a
covariance matrix, which the robust standard error cannot use, the
default is `"ml_logistic"`, as in
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md).
Under `wald_df = "missing_information"` the default is `"mlr"` (`"ml"`
with a covariance matrix), since that rule is calibrated for the
intervals on the scale of omega only (see below).

**A *t* reference distribution for the Wald intervals.** The Wald
intervals (`"ml"`, `"mlr"`, `"adf"`, and their `_logistic` variants) are
the estimate plus or minus a normal quantile times its standard error.
The normal quantile treats the standard error as known, which it is only
in the limit; in small samples, and more so when values are missing,
such intervals are too short. Von Hippel (2016, pp. 431–433) showed this
for maximum likelihood estimates from incomplete data and proposed
keeping the estimate and its standard error while replacing the normal
quantile with a *t* quantile (his Eq. 29). The argument `wald_df` does
that: the interval becomes \\\hat{\omega} \pm t\_{1 - \alpha/2, \nu}
\mathrm{SE}\\ with \\\nu\\ the value of `wald_df`, and for the
`_logistic` variants the *t* quantile is applied on the logit scale
before the limits are transformed back. The default, `Inf`, is the
normal reference. A number is the analyst's own choice of degrees of
freedom, best used to see how much a plausible range of them would
lengthen the interval; with a few hundred cases the answer is very
little (on 200 degrees of freedom the 95% quantile is 1.972 rather than
1.960).

The calibration study described next selected its simulated data sets by
an eligibility rule that keeps some draws the affine-collapse diagnostic
sets aside, so its results are conditional on that rule, and coverage
under the diagnostic's selection has not been established by that study.
A raw-data eligibility replay alone cannot validate the FMI refits and
information diagnostics. The keyword is an approximate, opt-in rule, not
a coverage guarantee.

The keyword `wald_df = "missing_information"` takes the degrees of
freedom from omega's own fraction of missing information \\\gamma\\ by
von Hippel's (2016, Eqs. 31 and 34) effective sample size rule with \\k
= 1\\, \$\$\nu = \max\\3, N(1 - \gamma) - 1\\,\$\$ where \\N\\ is the
number of cases the fit used and \\\gamma\\ is estimated from the fit as
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
estimates it; with complete data (listwise deletion, or a covariance
matrix) \\\gamma = 0\\ and \\\nu = N - 1\\. The rule comes from a
calibration study of these intervals under full information maximum
likelihood (in the package's source repository,
`methodology/missing_data`;
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
describes its design), which simulated the model implied omega of one
factor models with six items (*N* from 30 to 200) or 12 (*N* = 50 and
100), normal data fit by ML and nonnormal data by MLR. Over the cells
with \\N \le 100\\ and missing values the coverage of the *t* interval
was .0180 from .95 on average, against .0242 for the normal one, and no
higher than .9550 in any cell, with normal theory standard errors
(`"ml"`) and Huber-White ones (`"mlr"`) alike; with nonnormal data both
intervals still covered below .95, the *t* one less so. On the logit
scale the rule did not pass: with normal data the `_logistic` interval
on the normal reference was already within .0044 of .95 on average, and
the *t* reference raised its coverage as high as .9645, above the .96
the study allowed. The keyword therefore applies to `ci_method = "ml"`
and `"mlr"` only (with `ci_method` not supplied, the interval is
`"mlr"`, or `"ml"` with a covariance matrix), under listwise deletion or
`missing = "fiml"`, and is an error with the `_logistic` variants, with
`"adf"` and `"adf_logistic"` (the study evaluated maximum likelihood
standard errors only), with two-stage estimation (for which no fraction
of missing information is computed), and with every interval that is not
a Wald interval. The study fit no auxiliary variables: with `aux` the
rule takes omega's fraction of missing information from the fit that
carries them, which is smaller when they predict the missing values, and
the study supports that case by analogy rather than directly. Under
listwise deletion the standard errors here come from the expected
information (`"ml"`, from a covariance matrix too) and the robust.sem
sandwich (`"mlr"`, which needs raw data), where the study's cells
without missing values used the observed information and the Huber-White
sandwich; for omega at *N* = 30 and 50 they are shorter by less than 1
and by 1 to 2 percent, so the complete data intervals are, if anything,
slightly shorter than those the study evaluated. The printed table names
the reference distribution in a line under it, under the keyword with
*N* and the FMI; the `wald_df` attribute records the degrees of freedom,
and under the keyword the `wald_df_table` attribute holds the
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
table they came from.

**Missing data and auxiliary variables.** By default
(`missing = "fiml"`) incomplete data are analyzed by full information
maximum likelihood (FIML), which keeps every case with at least one
observed item or observed auxiliary variable and fits the single-factor
model to all of them, consistent and efficient under the missing at
random (MAR) assumption; `missing = "listwise"` fits the complete rows
alone. Under the default, the model implied omega here and the omega
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
reports for the same items agree on incomplete data, the estimate and,
with `ci_method = "ml_logistic"`, the interval; the default interval,
`"mlr_logistic"`, forms its limits the same way from the robust standard
error. Data in which no value is missing, and a covariance matrix, are
analyzed as they stand: FIML would give the same estimate, and the
standard errors and intervals come from the expected information (and
for `"mlr"` the robust.sem sandwich) where FIML takes the observed
information (and the Huber-White sandwich). A row with every item
missing does not count as incomplete data: without auxiliary variables,
listwise deletion and FIML both drop it, so data that are complete but
for such rows (a blank line at the end of a spreadsheet, say) are
analyzed as they stand too. The default analysis is therefore FIML
whenever an item value is missing and the complete data analysis
otherwise.
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
differs: it fits complete data by FIML as well. The complete data
analysis is used here because the closed-form intervals (`"feldt"`,
`"fisher"`, `"bonett"`, `"hakstian_whalen"`), `"adf"`, `"adf_logistic"`,
and `"likelihood"` need it, FIML giving none of them. On complete data
the standard error of omega from
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
therefore differs slightly from the one here: by 0.4 percent for
`t1_visual_perception`, `t2_cubes`, `t3_paper_form_board`, and
`t4_lozenges` in `holzinger_swineford`. When the default meets
incomplete items and one of those intervals, the call stops rather than
deleting cases unasked, and names the two choices:
`missing = "listwise"`, which computes the interval from the complete
rows alone, or an interval FIML gives (`"ml"`, `"mlr"`, their
`_logistic` variants, or a bootstrap method). The `aux` argument names
auxiliary variables: columns of `data` that are not part of the
composite but are correlated with the items or with the reasons values
are missing. They enter as *saturated correlates* (Graham, 2003):
correlated freely with each other and with every item's residual, never
loading on the factor, so the measurement model is undisturbed while
FIML uses their information; a good auxiliary also makes MAR itself more
plausible (Collins, Schafer, & Kam, 2001). A case whose items are all
missing is kept when an auxiliary variable is observed on it: the model
estimates the items and the auxiliaries jointly, so such a case informs
the item moments through its auxiliaries, and setting it aside would
make the analysis conditional on having an observed item, which biases
the estimates when the items go missing together at random given the
auxiliaries. Supplying `aux` without `missing` implies
`missing = "fiml"`, whether or not a value is missing; combining it with
an explicit `missing = "listwise"` is an error. Under `missing = "fiml"`
the denominators are estimated as follows: with `"model_implied"` the
fitted FIML model supplies the total variance directly, and with
`"observed"` the total variance is the sum of the FIML estimate of the
item covariance matrix (the EM estimate under a saturated model over the
items and any auxiliaries, the matrix
[`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md)
returns), which is the estimate of \\\mathrm{Var}(Y)\\ that uses the
partially observed rows; it is already on the maximum likelihood
divisor, matching the fitted loadings. The available intervals under
`"fiml"` are `"ml"`, `"mlr"` (Yuan & Bentler, 2000, robust), their
`_logistic` variants, and the bootstrap methods, which resample rows
including the partially observed ones and refit by FIML; the
complete-data closed forms, `"adf"`, and `"likelihood"` are errors
rather than silent fallbacks to listwise deletion.

`missing = "two_stage"` estimates in two stages instead (Yuan & Bentler,
2000; Savalei & Bentler, 2009): the means and covariances of the items
are first estimated by maximum likelihood under the saturated normal
model (the EM estimates, the same saturated estimate FIML's observed
denominator uses), and the single-factor model is then fit to them by
the complete data maximum likelihood discrepancy function. The auxiliary
variables named in `aux` enter the first stage only, where they inform
the estimated moments, and the measurement model is fit as if they were
absent, which is simpler than the saturated correlates FIML needs.
Two-stage estimation uses the cases FIML uses, every case with an
observed item or an observed auxiliary variable: the first stage
estimates the items and the auxiliaries jointly, so a case observed on
an auxiliary alone informs the item moments through it, and setting such
a case aside would condition the analysis on an observed item and bias
the estimates when the items went missing together at random given the
auxiliaries. lavaan's two-stage estimator reads such a case in its first
stage but leaves it out of its sample size, which would make its
standard errors too large; the fits go through
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
which computes them again with every case the first stage read, and a
message counts such cases (when a fit cannot be computed again, a
warning says that the standard error is too large and that a bootstrap
interval counts every case). The estimates are consistent under MAR and
only marginally less efficient than FIML's (Savalei & Bentler, 2009),
and with incomplete nonnormal data Savalei and Falk (2014) found the
coverage of robust two-stage intervals better than that of robust FIML
intervals. With `"model_implied"` the `"ml"` interval uses the normal
theory two-stage standard error (Savalei & Bentler, 2009, Eq. 6) and
`"mlr"` the standard error robust to nonnormality (Yuan & Bentler, 2000;
lavaan's `missing = "robust.two.stage"`); with `"observed"` the
numerator comes from the two-stage loadings and the denominator from the
EM covariance matrix, so the coefficient differs from its FIML
counterpart only through the loadings, and its interval is a bootstrap
that refits both stages on every replication. A delta method interval
for that coefficient can be derived under two-stage estimation, but its
behavior has not been studied, so it is not offered. With three items
the single-factor model is saturated and two-stage estimation gives the
FIML estimate and standard error. The intervals that FIML refuses are
refused here too, and `"likelihood"` for a further reason: no likelihood
of the observed data is maximized. Multiple imputation is a different
feature with a different interface and is out of scope here.

Both methods need the covariance of every pair of items and auxiliary
variables, and the data carry no information about the covariance of a
pair that no case observes together: lavaan only warns, and its EM
returns an arbitrary value for it. Nor do they carry any about how two
variables vary together when the cases that observe both hold a single
value of one of them (or are one case): the estimate of their covariance
then rests on the form of the normal model rather than on the data, or
is not determined at all. Either pair is refused by name before any fit,
with a pointer to
[`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md),
and a refusal of listwise deletion offers `missing = "fiml"` only when
the data hold no such pair. In a bootstrap, a resample that leaves such
a pair (or, for the BCa acceleration, a jackknife deletion that does) is
set aside with the same message and counted, under either method.

**Improper solutions.** A single-factor fit can converge to an improper
solution, one with a negative error variance (a Heywood case). McDonald
(2004), following van Driel (1978), separates improper estimates that
reflect sampling error around a proper structure, typically a negative
error variance near zero in a small sample or for an item that nearly
duplicates the factor, from those that reflect a misspecified or
empirically underidentified model: the first kind needs no change of
model, while the second calls for a different model rather than an
estimate. One sample cannot always tell the two apart, so the function
draws the line at what the coefficient can mean. When the error
variances still sum to a nonnegative number, the model implies a
nonnegative error variance for the composite and a model implied
\\\omega\\ between zero and one, so the coefficient is returned with a
warning. When they sum to less than zero, the model implies a negative
error variance for the composite and a model implied \\\omega\\ above
one, which no reliability coefficient can take; sampling error near a
zero error variance does not produce this unless every item is nearly
error free, and it is the signature of a solution that has run away, one
error variance falling without bound as its loading grows. The function
then stops with a message naming the item. So it does for three items
whose correlations have a negative product: a single factor implies that
the product of the three correlations equals the square of the product
of the three loadings, so no real loadings reproduce them and the
maximum likelihood estimate does not exist. Whatever the sum, a fit
whose likelihood has no maximum, an error variance running off below
zero while the likelihood is flat or still rising (checked by fitting
the model again with the variance held further out, as
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
describes), stops the function with an error of class `dmar_no_maximum`
that names the item and prints the lavaan code of the single-factor fit,
where `bounds = "pos.var"` or a different model are the choices.

A bootstrap meets improper solutions more often than the sample does,
since a resample can push an error variance that is near zero in the
sample below it, which is a standardized loading past one (the fits are
on the items' own metric, where a loading above one is routine). The
argument `improper` governs the replications whose solution is improper
but still yields a coefficient: `"keep"`, the default, uses their
estimates like any other, and `"drop"` sets them aside. Keeping is the
default because setting them aside conditions the bootstrap distribution
on admissibility and so trims its upper tail, the replications in which
an error variance fell below zero; Yuan and Hayashi (2003) make the
parallel argument for test statistics, counting a replication that fails
because its resample lies far from the model rather than discarding it.
A replication whose solution yields no coefficient is set aside under
either setting. When its error variances sum below zero, the solution
has run away and implies an \\\omega\\ above one, beyond the top of the
distribution like the improper replications kept. When three items have
no single-factor solution, the side depends on the scale the resample is
read on: on the estimator's own scale it continues the top of the
distribution, since as a resample approaches that boundary one loading,
and \\\omega\\ with it, grows without bound, while its correlations, one
of which is negative or zero, point to low reliability. The argument of
Yuan and Hayashi does not settle it, since their extreme value belongs
to a statistic that measures misfit and so has a direction, so such a
replication is counted, on its own as well as among the improper ones,
rather than placed on either side. A replication whose likelihood has no
maximum is set aside as well, and counted among those that raised an
error, under a message that names the item whose error variance ran off;
the run goes on. Every such replication is counted in the `bootstrap`
attribute and the printed line, and when the replications set aside and
the improper ones kept together reach the tail area of the interval,
\\(1 - C)/2\\ of `B`, a caution says that the limits depend on how they
were treated. When some of them had no single-factor solution, the
caution adds that the lower limit is likely too high: with three weakly
related items the estimator itself runs high near that boundary, and the
resamples beyond it are those whose correlations point to low
reliability. When about half of the replications are improper, neither
setting gives an interval worth reporting, and the counts are what the
result can state.

**Comparison with other packages.** The psych package provides `omega`,
which fits a Schmid-Leiman hierarchical factor model and reports several
variants of \\\omega\\ (\\\omega_t\\, \\\omega_h\\) alongside extensive
psychometric diagnostics. `reliability_omega` in DMAR differs in
emphasis: it implements McDonald's \\\omega\\ from a single-factor
(congeneric) model and accompanies the point estimate with a confidence
interval drawn from the methods compared in Kelley and Pornprasertmanit
(2016). The same denominator distinction appears in semTools'
`compRelSEM()` as its `obs.var` argument.

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

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables. *Structural
Equation Modeling, 16*, 477–497.
[doi:10.1080/10705510903008238](https://doi.org/10.1080/10705510903008238)

Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
robust full information maximum likelihood with incomplete nonnormal
data. *Structural Equation Modeling, 21*(2), 280–302.
[doi:10.1080/10705511.2014.882692](https://doi.org/10.1080/10705511.2014.882692)

van Driel, O. P. (1978). On various causes of improper solutions in
maximum likelihood factor analysis. *Psychometrika, 43*, 225–243.

von Hippel, P. T. (2016). New confidence intervals and bias comparisons
show that maximum likelihood can beat multiple imputation in small
samples. *Structural Equation Modeling, 23*(3), 422–437.
[doi:10.1080/10705511.2015.1047931](https://doi.org/10.1080/10705511.2015.1047931)

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

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*, 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

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

McDonald, R. P. (2004). Respecifying improper structures. *Structural
Equation Modeling, 11*(2), 194–209.

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

Yuan, K.-H., & Hayashi, K. (2003). Bootstrap approach to inference and
power analysis based on three test statistics for covariance structure
models. *British Journal of Mathematical and Statistical Psychology,
56*, 93–110.

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
(single-factor CFA used internally), `omega`.

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
#>  term       value
#>  estimate   0.821
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Denominator: observed composite variance (robust omega)
#> Largest residual correlation under the single-factor model: -0.043 (y2 with y5)

# The closed-form standard errors are derived under the model implied
# ratio, so with the observed denominator the interval comes from a
# bootstrap, which refits the single factor model once per
# replication. B = 10 keeps the example quick and is below the 400 a
# reported interval needs; the option lets it run and the warning says
# so. A reported interval deserves the default B = 10000, an accuracy
# choice rather than a formality (see Details). The percentile
# interval is what Kelley and Pornprasertmanit (2016) recommend for
# this coefficient, with ci_method = "bca" as the alternative, and
# seed makes the interval reproducible.
old <- options(dmar.bootstrap_exploratory = TRUE)
reliability_omega(data = items, ci_method = "percentile", B = 10,
                  seed = 113)
#> Warning: reliability_omega(): B = 10 is below the 400 replications a 95% interval is reported from (10 beyond each limit). The limits are shown because options(dmar.bootstrap_exploratory = TRUE) is set. With 10 replications each limit rests on about 0.3 of them and will move from seed to seed; use the result to explore or to demonstrate, and report an interval only from B of at least 400. B, B_used, and B_required are recorded in the result.
#>  term       value
#>  estimate   0.821
#>  se         0.016
#>  ci_lower   0.794
#>  ci_upper   0.840
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Denominator: observed composite variance (robust omega)
#> Largest residual correlation under the single-factor model: -0.043 (y2 with y5)
#> 
#> Bootstrap: 10 replications attempted, 10 returned a value; at least 400 are required for a reported interval. Exploratory: below the number required for a reported interval
options(old)

# Model implied denominator with its closed-form robust ML interval,
# whose limits are formed on the logit scale; ci_method = "mlr" gives
# limits formed on the scale of omega.
reliability_omega(data = items, denominator = "model_implied")
#>  term       value
#>  estimate   0.820
#>  se         0.019
#>  ci_lower   0.779
#>  ci_upper   0.854
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Denominator: composite variance implied by the fitted model
#> Largest residual correlation under the single-factor model: -0.043 (y2 with y5)

# Two further routes into the same coefficient. The first works from
# the summary statistics a paper reports, a covariance matrix and its
# sample size, with the maximum likelihood interval that summary
# input supports.
reliability_omega(S = cov(items), N = 200,
                  denominator = "model_implied", ci_method = "ml")
#>  term       value
#>  estimate   0.820
#>  se         0.020
#>  ci_lower   0.781
#>  ci_upper   0.858
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Denominator: composite variance implied by the fitted model
#> Largest residual correlation under the single-factor model: -0.043 (y2 with y5)

# How much would a t reference lengthen that interval? On 30 degrees of
# freedom, an arbitrary small value chosen to show the effect, the
# limits move out a little; the line under the table says so.
reliability_omega(S = cov(items), N = 200,
                  denominator = "model_implied", ci_method = "ml",
                  wald_df = 30)
#>  term       value
#>  estimate   0.820
#>  se         0.020
#>  ci_lower   0.780
#>  ci_upper   0.860
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Wald interval: t reference on 30 degrees of freedom (wald_df), in place of the normal
#> 
#> Denominator: composite variance implied by the fitted model
#> Largest residual correlation under the single-factor model: -0.043 (y2 with y5)

# The second is full information maximum likelihood with an auxiliary
# variable, where missingness on y2 depends on an auxiliary z (missing
# at random given z). Supplying aux implies missing = "fiml"; the
# N_complete row shows how many rows listwise deletion would have
# kept. With values missing, the interval can also refer to a t
# distribution on degrees of freedom taken from omega's own fraction of
# missing information, here in the fit that carries z (the calibration
# study behind the rule had no auxiliary variable, so this case rests
# on analogy; see Details); the line under the table gives N, the FMI,
# and the degrees of freedom that result. The rule is calibrated for
# limits on the scale of omega, so the interval here is the "mlr" one.
z <- eta + rnorm(200, sd = 0.5)
d <- data.frame(items, z = z)
d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
reliability_omega(data = d, aux = "z", denominator = "model_implied",
                  wald_df = "missing_information")
#>  term       value
#>  estimate   0.819
#>  se         0.020
#>  ci_lower   0.780
#>  ci_upper   0.859
#>  conf_level 0.950
#>  N          200  
#>  N_complete 140  
#>  n_items    6    
#> 
#> Wald interval: t reference on 187.3 degrees of freedom, from von Hippel's (2016) max(3, N(1 - FMI) - 1) with N = 200 and FMI = 0.05845 (wald_df = "missing_information"), in place of the normal
#> 
#> The FMI degrees-of-freedom rule is approximate. Its calibration study admitted data sets that the affine-collapse screening sets aside, so coverage under that screening has not been established. See the method's help for its scope and conditional evidence.
#> 
#> Denominator: composite variance implied by the fitted model
#> Largest residual correlation under the single-factor model: -0.075 (y2 with y5)

# Two-stage estimation instead: z enters only the first stage, where the
# means and covariances are estimated, and the interval uses the robust
# two-stage standard error, its limits again formed on the logit scale.
reliability_omega(data = d, aux = "z", missing = "two_stage",
                  denominator = "model_implied")
#>  term       value
#>  estimate   0.819
#>  se         0.020
#>  ci_lower   0.776
#>  ci_upper   0.855
#>  conf_level 0.950
#>  N          200  
#>  N_complete 140  
#>  n_items    6    
#> 
#> Denominator: composite variance implied by the fitted model
#> Largest residual correlation under the single-factor model: -0.070 (y2 with y5)
```
