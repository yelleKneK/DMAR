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
  missing = c("fiml", "listwise", "two_stage"),
  aux = NULL,
  ci_method = NULL,
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
  matrix. SPSS user-defined missing values, which
  `haven::read_sav(user_na = TRUE)` keeps as values of a
  `haven_labelled_spss` column, are missing, as SPSS treats them.

- S:

  A symmetric covariance matrix among the items. If `S` is supplied, `N`
  must also be supplied; raw-data-only confidence interval methods
  (`"adf"`, `"bootstrap_*"`, `"percentile"`, `"bca"`) are then
  unavailable, as are `missing = "fiml"` and `missing = "two_stage"`
  named explicitly, and `aux`, since a covariance matrix has no
  incomplete cases (the default of `missing` analyzes it as it stands).

- N:

  Total sample size; required when `S` (rather than `data`) is supplied.

- estimator:

  How \\\alpha\\ is estimated: `"analytic"` (the default) applies the
  closed-form equation to the observed covariance matrix, and
  `"model_implied"` takes the reliability implied by the
  \\\tau\\-equivalent single-factor model fit by maximum likelihood. See
  *Details*.

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
  both treatments drop), and `S` are analyzed as they stand, with the
  complete data default interval, since the complete-data intervals need
  that analysis; incomplete items are analyzed by FIML, or, when
  `ci_method` names an interval FIML cannot give, the call stops and
  names the choices. See *Details*. With `missing = "listwise"`, when
  fewer than two rows are complete, listwise deletion leaves nothing to
  analyze, and the function stops before any fit, says so, and names the
  treatment that uses the incomplete rows. The analytic estimator needs
  only the trace and the sum of the covariance matrix, which any two
  rows give, so it is computed from two or more complete rows. The model
  implied estimator fits the single-factor model by maximum likelihood,
  which needs the covariance matrix positive definite, and \\N\\
  complete rows give a matrix of rank at most \\N - 1\\; with no more
  complete rows than items it stops before the fit, says why, and names
  the treatment that uses the incomplete rows. Ordinary FIML is
  evaluated by the observed values and missingness patterns, not a
  blanket count of complete cases. Analytic alpha and two-stage
  estimation need unrestricted moments. The likelihood diagnostic
  refuses a certified affine-collapse direction or a numerical affine
  instability, as described in
  [`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md);
  EM convergence alone does not establish existence or uniqueness.
  Model-implied alpha can retain a proper identified constrained FIML
  fit with H0-based errors when H1 diagnostics are unavailable; the
  returned result states this limit. The complete rows behind the
  model-implied listwise estimator are counted as distinct rows.
  Bootstrap resamples retain a separate distinct-case estimability
  screen, which is a policy, not a proof that every smaller constrained
  FIML likelihood is unbounded. The refusal message for the analytic
  estimator names listwise deletion when two or more rows are complete
  and listwise deletion would accept them, and otherwise says why it
  cannot help either. A refusal of listwise deletion names
  `missing = "fiml"` only when that method could use the data, so the
  two refusals never point at each other, and when two columns share a
  name (which listwise deletion without `aux` reads by position, and
  that method refuses) it says to name them apart first.

- aux:

  Optional character vector naming auxiliary variable columns of `data`:
  entered as saturated correlates under full information maximum
  likelihood, and into the first stage under `missing = "two_stage"`.
  Supplying `aux` without `missing` implies `missing = "fiml"`. See
  *Details*.

- ci_method:

  Method for constructing the confidence interval. See *Details* for
  which methods each estimator supports. The default `NULL` resolves to
  `"bonett"` for the analytic estimator and to `"mlr"` for the model
  implied estimator (or `"ml"` when only a covariance matrix is
  supplied, since `"mlr"` needs raw data). Under FIML (the default of
  `missing` with incomplete items) and two-stage estimation the analytic
  default is `"ml"`, whose standard error comes from the information
  matrix of the saturated model; with complete data the default of
  `missing` keeps `"bonett"`.

- conf_level:

  Confidence level for the interval (1 - Type I error rate). Defaults to
  `0.95`.

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
  interval. When set, the seed applies inside the function only, and the
  user's generator state is restored when the function returns.

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

  For `estimator = "model_implied"` with a bootstrap interval: how a
  replication whose tau-equivalent single-factor solution is improper (a
  negative error variance, a Heywood case) is treated when the solution
  still yields a coefficient. `"keep"` (the default) uses its estimate
  like any other, and `"drop"` sets it aside; either way such
  replications are counted, and one whose solution yields no coefficient
  is set aside under both. With `"drop"` an improper point estimate
  stops the call. See
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  for the reasoning. The analytic estimator fits no model, so no
  replication of it can be improper, and the argument has no effect
  there beyond being recorded as the `improper_rule` of the `bootstrap`
  attribute.

- wald_df:

  The degrees of freedom of the *t* distribution to which a Wald
  interval refers: `Inf` (the default) is the normal reference, and a
  single positive number replaces the normal quantile with the *t*
  quantile on that many degrees of freedom. It applies to `"ml"`,
  `"mlr"`, `"adf"`, and their `_logistic` variants; a finite value with
  any other interval is an error. The keyword `"missing_information"` of
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  is an error here: its calibration found that the *t* reference gives
  alpha's FIML interval no gain and too much coverage. See *Details*.

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
(sample coefficient \\\alpha\\), `"se"` (standard error on the
coefficient scale, `NA` for methods that do not produce one; for the
transformation-based intervals `"fisher"`, `"bonett"`, and
`"hakstian_whalen"` it is the delta method back-transform of the
transformation-scale standard error, evaluated at the estimate),
`"se_transformed"` (only for those transformation-based intervals: the
standard error on the transformation scale, the quantity the interval is
built from, with the scale named in the attribute `se_transform_scale`:
`"fisher_z"`, `"log(1-alpha)"`, or `"cube_root"`), `"ci_lower"` and
`"ci_upper"` (each clamped into \[0, 1\], so the limits stay in order),
`"conf_level"`, `"N"` (the cases the analysis used: the complete cases
under listwise deletion, and under `"fiml"` and `"two_stage"` every case
with at least one observed item or observed auxiliary variable),
`"N_complete"` (the complete cases, so the cost of listwise deletion is
visible at a glance; equal to `"N"` under listwise deletion and for
covariance input), and `"n_items"` (number of items). The selected
coefficient, CI method, and missing-data treatment travel as the
attributes `coefficient`, `ci_method`, `missing`, `wald_df` (the degrees
of freedom of the Wald interval's *t* reference, `Inf` for the normal
one and for intervals that are not Wald intervals), and (when supplied)
`aux`; bootstrap calls also record `B`. A bootstrap call also carries
the attribute `bootstrap`, a list recording the replications the
interval conditions on: `attempted` (`B`), `returned` (the replications
whose value the interval uses; the standard error is their standard
deviation and the percentile and BCa limits are their quantiles),
`failed` (the difference), `reasons` (a named count of why replications
were set aside: `error` when the refit raised an error, `non_finite`
when it was set aside without one, as a fit that did not converge is,
and as an improper one is under `improper = "drop"` even though its
estimate is finite), `error_messages` (the distinct error messages with
their counts, empty when none), `undetermined` (the replications set
aside because the resample left a polychoric correlation undetermined,
which only
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
counts; zero here), `improper` and `improper_kept` (the replications
whose tau-equivalent single-factor solution, under
`estimator = "model_implied"`, was improper, with a negative error
variance: those set aside, which are among the `non_finite` ones, and
those whose estimate the interval uses, which are among the `returned`
ones), `no_solution` (zero here: it counts the replications of three
items with no congeneric single-factor solution, which
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
meets), `improper_rule` (`"keep"` or `"drop"`; see the argument
`improper`), `not_converged` (zero here: a fit that did not converge is
counted among the `non_finite` ones), `too_few_cases` (the replications
set aside because the resample held too few distinct cases for the rule
on the number of cases, which are among the `error` ones; when they
number more than half of `B`, no interval is formed), `caution` (`TRUE`
when the replications set aside and the improper ones kept number at
least \\(1 - C)/2\\ of `B`, with \\C\\ the confidence level, so that how
they were treated can decide where a limit falls), `require_half`
(`TRUE`: an interval requires at least half of all attempted
replications to return defined estimates; this is an estimability
safeguard, not a coverage guarantee), `required` (the fewest
replications a request for an interval of this method and level may ask
for, a count held to `B` and not to the replications that return a
value: 400 for `"percentile"`, 2000 for `"bca"`, 100 for the two
bootstrap standard error methods at the 95 percent level; see
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
ordered categorical items see
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md).

Coefficient \\\alpha\\ is at most one but, as a function of the
covariance matrix, has no lower bound: it falls below zero when the item
covariances are negative on balance, which can arise when an item is
scored in the opposite direction to the others and needs reverse
scoring. A reliability coefficient cannot be negative, so a negative
estimate is returned as computed, with a warning that says so, and each
confidence limit is clamped into \[0, 1\]; the interval around an
estimate far enough below zero is then \[0, 0\], which says that the
data support no positive reliability for the items as keyed.

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

**Missing data and auxiliary variables.** Listwise deletion is unbiased
only when the data are missing completely at random and is inefficient
always; the `mlmr` vignette develops the argument at length. By default
(`missing = "fiml"`) incomplete data are therefore analyzed by full
information maximum likelihood (FIML), which keeps every case with at
least one observed item and is consistent and efficient under the weaker
missing at random (MAR) assumption; `missing = "listwise"` analyzes the
complete rows alone. Data in which no value is missing, and a covariance
matrix, are analyzed as they stand: FIML would give the same
coefficient, and the interval is the complete data one, Bonett's by
default for the analytic estimator. A row with every item missing does
not count as incomplete data: without auxiliary variables, listwise
deletion and FIML both drop it, so data that are complete but for such
rows (a blank line at the end of a spreadsheet, say) are analyzed as
they stand too. The default analysis is therefore FIML whenever an item
value is missing and the complete data analysis otherwise.
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
differs: it fits complete data by FIML as well. The complete data
analysis is used here because the closed-form intervals (`"bonett"`,
`"feldt"`, `"fisher"`, `"hakstian_whalen"`), `"adf"`, `"adf_logistic"`,
and `"likelihood"` need it, FIML giving none of them. When the default
meets incomplete items and one of those intervals, the call stops rather
than deleting cases unasked, and names the two choices:
`missing = "listwise"`, which computes the interval from the complete
rows alone, or an interval FIML gives for the estimator in use (`"ml"`,
`"ml_logistic"`, with the model implied estimator `"mlr"` and
`"mlr_logistic"`, or a bootstrap method); when `ci_method` is not
supplied, the analytic estimator's interval under FIML is `"ml"`. The
`aux` argument names auxiliary variables: columns of `data` that are not
part of the composite but are correlated with the items or with the
reasons values are missing. They are entered as *saturated correlates*
(Graham, 2003): correlated freely with each other and with every item's
residual, never loading on the factor and never entering the composite,
so the measurement model is undisturbed while FIML uses their
information. Beyond recovering information, a good auxiliary makes the
MAR assumption itself more plausible, since missingness that depends on
the auxiliary becomes MAR once the auxiliary is conditioned on (Collins,
Schafer, & Kam, 2001). A case whose items are all missing is kept when
an auxiliary variable is observed on it: the saturated model estimates
the items and the auxiliaries jointly, so such a case informs the item
moments through its auxiliaries, and setting it aside would make the
analysis conditional on having an observed item, which biases the
estimates when the items go missing together at random given the
auxiliaries. Supplying `aux` without `missing` implies
`missing = "fiml"`, whether or not a value is missing; combining it with
an explicit `missing = "listwise"` is an error. How each estimator uses
FIML: the model implied estimator simply fits its model with
`missing = "ml"`; the analytic estimator applies the classical formula,
unchanged, to the FIML estimate of the item covariance matrix (the EM
estimate under a saturated model over the items and any auxiliaries, the
matrix
[`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md)
returns), so the estimand stays the classical coefficient and only the
covariance matrix it is computed from improves. Under `missing = "fiml"`
the available intervals are `"ml"` and `"ml_logistic"` (both estimators;
the standard error comes from the FIML information matrix), `"mlr"` and
`"mlr_logistic"` (model implied estimator; the Yuan-Bentler robust
standard error), and the bootstrap methods, which resample rows
(including the partially observed ones) and refit by FIML on each
replication. The complete-data closed forms (`"feldt"`, `"fisher"`,
`"bonett"`, `"hakstian_whalen"`), `"adf"`, and `"likelihood"` are errors
with `missing = "fiml"` rather than silently reverting to listwise
deletion.

`missing = "two_stage"` estimates in two stages (Yuan & Bentler, 2000;
Savalei & Bentler, 2009): the means and covariances of the items are
first estimated by maximum likelihood under the saturated normal model
(the EM estimates, with any `aux` variables entering this stage only),
and the model is then fit to them by the complete data maximum
likelihood discrepancy function. The analytic coefficient is a function
of the first stage alone, which is the saturated maximum likelihood
estimate FIML uses too, so its value and its intervals are those of
`missing = "fiml"` computed from the same cases. Those are the cases
FIML uses: every case with an observed item or, with `aux`, an observed
auxiliary variable, since the first stage estimates the items and the
auxiliaries jointly and a case observed on an auxiliary alone informs
the item moments through it. Setting such a case aside would condition
the analysis on an observed item and bias the estimates when the items
went missing together at random given the auxiliaries. lavaan's
two-stage estimator reads such a case in its first stage but leaves it
out of its sample size, which would make its standard errors too large;
the model implied fits go through
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
which computes them again with every case the first stage read, and a
message counts such cases (when a fit cannot be computed again, a
warning says that the standard error is too large and that a bootstrap
interval counts every case). The model implied coefficient comes from
the tau-equivalent model fit in the second stage, with the normal theory
two-stage standard error (Savalei & Bentler, 2009, Eq. 6) for `"ml"` and
the one robust to nonnormality (Yuan & Bentler, 2000; Savalei & Falk,
2014) for `"mlr"`. The intervals FIML refuses are refused here too,
`"likelihood"` also because no likelihood of the observed data is
maximized. Multiple imputation is a different feature with a different
interface and is out of scope here.

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
  first-order accurate (DiCiccio & Efron, 1996). The acceleration is
  estimated as boot estimates it, by regression on the resampling
  frequencies when at least as many replications as cases returned a
  value and by the jackknife otherwise, or when that regression fails; a
  leave-one-out refit that fails is counted and set aside rather than
  stopping the call, and when the acceleration or the bias correction
  cannot be estimated at all the percentile interval is returned in its
  place, with a warning and the substitution recorded in the `ci_method`
  attribute; see
  [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  for the full account. `"bootstrap_se"` uses the standard deviation of
  the bootstrap estimates as a standard error in an ordinary
  normal-theory interval; the `_logistic` variant builds that interval
  on the logit scale so the endpoints respect \[0, 1\], with the
  standard deviation carried to that scale by the delta method, divided
  by \\\hat\alpha(1 - \hat\alpha)\\, as in Equation 32 of Kelley and
  Pornprasertmanit (2016) and in `MBESS::ci.reliability()`. When every
  replication returns the same value, as when every item is an exact
  linear function of the same variable (an item entered twice beside
  other items does not do it), the bootstrap distribution is a single
  point and every method returns that point at both limits, with a
  warning; so does the `_logistic` variant at a value of one (identical
  items), where the logit is infinite. Replications on which the
  estimator cannot be computed (for example, a model refit that does not
  converge) are set aside, and so is a replication in which an item
  takes a single value (a resample that has left out every case on which
  a rare item differs from the rest, or the jackknife deletion of such a
  case), which is refused as the data would be, by name. A replication
  follows the rule on the number of cases as the data do (see
  `missing`): a resample repeats rows, and one with no more distinct
  cases than the rule asks for is set aside, so near the rule's bound
  many are (at 12 cases and eight variables, about three resamples in
  four). When more than half are, no interval is reported (the standard
  error and the limits are `NA`): the replications left hold more
  distinct cases than a resample typically does, so their spread
  understates the sampling variability of the coefficient, and at the
  bound itself each holds every distinct case and, without repeated
  rows, gives the estimate. Without repeated rows that happens at fewer
  cases than about 1.65 times the count the rule compares them with
  (fewer than 14 for eight variables). The interval is computed from the
  replications that return a value, whose quantiles the percentile and
  BCa limits are; the counts, and why the rest were set aside, travel in
  the `bootstrap` attribute (see *Value*), print under the table, and
  are stated in a warning whenever a replication was set aside. Requires
  raw data and the boot package. The default `B = 10000` is an accuracy
  choice: the BCa adjustment pushes the working quantiles farther into
  the tails of the bootstrap distribution than the percentile interval
  uses, and stabilizing them takes more replications than the customary
  2000; reduce `B` for exploration, not for a reported analysis.
  Bootstrap results vary from run to run; supply `seed` for
  reproducibility, and `ncpus` to spread the replications over several
  processes, which leaves the interval for a given seed unchanged. Under
  `estimator = "model_implied"` each replication refits the
  tau-equivalent model by maximum likelihood, whose estimate robust
  maximum likelihood shares, and a replication whose solution is
  improper is kept or set aside according to `improper`, as
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  describes; an improper point estimate whose error variances sum below
  zero, so that the model implies a negative error variance for the
  composite, stops with a message rather than returning a coefficient
  above one. (The likelihood of the tau-equivalent model has a maximum:
  one item's error variance cannot run off below zero while its loading
  grows, since the loadings are equal and every covariance would grow
  with them, so the stop of class `dmar_no_maximum` that
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  describes does not arise here.)

- `"none"`:

  Return only the point estimate.

**Two complete rows.** Coefficient \\\alpha\\ needs only the trace and
the sum of the covariance matrix, so it is computed from any two rows,
but the covariance matrix of two rows has rank one: every pair of items
correlates +1 or -1, whatever the items, and the estimate is a function
of the difference between the two rows alone. Feldt's interval is exact
for any *N* of at least 2 under multivariate normality and parallel
items (for 0/1 items, as in
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
it is an approximation; Feldt, 1965), and it is the interval the
refusals below name. The Hakstian-Whalen interval is defined there too,
but at *N* = 2 its upper limit lies above one for every data set whose
alpha is below one at any confidence level above about .90, and is
clamped to one (at an alpha of one, which two rows can give, the
interval is \\\[1, 1\]\\). The rest are refused at *N* = 2, each for its
own reason, as `"fisher"` is below *N* = 4: Bonett's variance,
\\2J/((J-1)(N-2))\\, is infinite; the `"ml"` and `"adf"` standard errors
(and their `_logistic` variants), which treat that matrix as the
population's, are zero for any data (the asymptotic variance of van Zyl,
Neudecker, and Nel, 2000, vanishes at a matrix of rank one, and each
case's centered cross products are then proportional to the covariance
matrix, which alpha does not change with); and every bootstrap resample
is one row twice or both rows, so the bootstrap distribution has no
spread. The model implied estimator needs more complete rows than items
(see `missing`).

**A *t* reference distribution for the Wald intervals.** The intervals
`"ml"`, `"mlr"`, and `"adf"`, with their `_logistic` variants, are Wald
intervals: the estimate plus or minus a normal quantile times its
standard error. The argument `wald_df` replaces the normal quantile with
the *t* quantile on `wald_df` degrees of freedom and leaves the estimate
and the standard error as they are (von Hippel, 2016, Eq. 29); for the
`_logistic` variants the *t* quantile is applied on the logit scale. The
default, `Inf`, is the normal reference, and a finite value with any
other interval is an error. A number is the analyst's own choice of
degrees of freedom, a sensitivity check;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
discusses it further.

The keyword `wald_df = "missing_information"`, which
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
and
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
accept, is an error here. The calibration study behind it (in the
package's source repository, `methodology/missing_data`;
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
describes its design) evaluated alpha's FIML interval
(`missing = "fiml"`, `ci_method = "ml"`) with normal data, one factor
models of six items (*N* from 30 to 200) or 12 (*N* = 50 and 100), and
15 or 30 percent of the values of half the items missing. On the normal
reference that interval was already close to .95: over the cells with
\\N \le 100\\ and missing values its coverage was .0043 from .95 on
average (between .940 and .958 in every cell). The *t* reference on von
Hippel's \\\max\\3, N(1 - \gamma) - 1\\\\ degrees of freedom came barely
closer (.0042), and in one cell its coverage reached .966, above the .96
that the study's decision rule, fixed before the run, allowed, so no
rule passed. That cell, with 12 items, *N* = 50, and 30 percent missing
completely at random, is the one where the study's fit of the saturated
model converged in only 53 percent of the replications (in 100
replications of that design this function's fit, which starts at the EM
estimates, converged in 91 and the study's in 49), so the coverage there
conditions on the samples the study's fit could handle; without it the
*t* reference would have passed (.0035 against .0042, at most .9595).
The verdict is the rule's as fixed in advance, and alpha keeps the
normal reference unless a number is given; the *t* reference would
change its interval little in any case. On the logit scale none of von
Hippel's four rules passed either. The printed table names the reference
distribution in a line under it, and the `wald_df` attribute records it.

**Comparison with other packages.** The psych package provides `alpha`,
which reports coefficient \\\alpha\\ along with item-level diagnostics
(item-total correlations, alpha-if-item-deleted, and several alternative
coefficients). The emphasis in psych is broad exploratory psychometric
reporting. `reliability_alpha` in DMAR differs in emphasis: it returns a
single point estimate alongside a confidence interval drawn from the
methods compared in Kelley and Pornprasertmanit (2016).

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

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables. *Structural
Equation Modeling, 16*, 477–497.
[doi:10.1080/10705510903008238](https://doi.org/10.1080/10705510903008238)

Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
robust full information maximum likelihood with incomplete nonnormal
data. *Structural Equation Modeling, 21*(2), 280–302.
[doi:10.1080/10705511.2014.882692](https://doi.org/10.1080/10705511.2014.882692)

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

von Hippel, P. T. (2016). New confidence intervals and bias comparisons
show that maximum likelihood can beat multiple imputation in small
samples. *Structural Equation Modeling, 23*(3), 422–437.
[doi:10.1080/10705511.2015.1047931](https://doi.org/10.1080/10705511.2015.1047931)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper that dispatches by coefficient),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
(single-factor CFA used internally),
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
`alpha`.

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
#>  term           value
#>  estimate       0.777
#>  se             0.025
#>  se_transformed 0.110
#>  ci_lower       0.723
#>  ci_upper       0.820
#>  conf_level     0.950
#>  N              200  
#>  N_complete     200  
#>  n_items        6    

# Same point estimate from a covariance matrix; CI requires N.
S <- cov(items)
reliability_alpha(S = S, N = 200, ci_method = "feldt")
#>  term       value
#>  estimate   0.777
#>  se         <NA> 
#>  ci_lower   0.725
#>  ci_upper   0.821
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    

# The bootstrap intervals resample the rows and recompute the
# coefficient on each replication. The percentile interval reads its
# limits off the empirical quantiles of the bootstrap estimates.
# B = 500 keeps the example quick (a 95 percent percentile interval
# is reported from 400 or more); a reported interval deserves the
# default B = 10000, and seed makes the interval reproducible.
reliability_alpha(data = items, ci_method = "percentile", B = 500,
                  seed = 113)
#>  term       value
#>  estimate   0.777
#>  se         0.023
#>  ci_lower   0.730
#>  ci_upper   0.815
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Bootstrap: 500 replications attempted, 500 returned a value

# The bias-corrected and accelerated interval adjusts those two
# quantile positions for median bias and for acceleration, and is the
# better choice for a reported interval. The default B = 10000 is an
# accuracy choice, not a formality, since BCa works farther into the
# tails of the bootstrap distribution than the percentile interval
# does (see Details). For the same reason B = 2000 is the fewest
# replications a 95 percent BCa interval is run with.
reliability_alpha(data = items, ci_method = "bca", B = 2000, seed = 113)
#>  term       value
#>  estimate   0.777
#>  se         0.024
#>  ci_lower   0.725
#>  ci_upper   0.819
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Bootstrap: 2000 replications attempted, 2000 returned a value

# An application that runs a long bootstrap (a Shiny app, for example)
# can follow it through a progress callback, which is called with the
# replications attempted so far and B, and can stop it by returning
# FALSE. The callback changes no result: this is the percentile
# interval above.
reliability_alpha(data = items, ci_method = "percentile", B = 500,
                  seed = 113, progress = function(done, total) {
                    if (done %% 250 == 0) message(done, " of ", total)
                  })
#> 250 of 500
#> 500 of 500
#>  term       value
#>  estimate   0.777
#>  se         0.023
#>  ci_lower   0.730
#>  ci_upper   0.815
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Bootstrap: 500 replications attempted, 500 returned a value

# Full information maximum likelihood with an auxiliary variable.
# Missingness on y2 depends on an auxiliary z (missing at random given
# z), so listwise deletion is biased and FIML with z is not. Supplying
# aux implies missing = "fiml"; the N_complete row shows how many
# rows listwise deletion would have kept.
z <- eta + rnorm(200, sd = 0.5)
d <- data.frame(items, z = z)
d$y2[runif(200) < plogis(-1 + 1.5 * as.numeric(scale(z)))] <- NA
reliability_alpha(data = d, aux = "z")
#>  term       value
#>  estimate   0.776
#>  se         0.025
#>  ci_lower   0.727
#>  ci_upper   0.826
#>  conf_level 0.950
#>  N          200  
#>  N_complete 140  
#>  n_items    6    
```
