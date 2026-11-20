# Coefficient *H* (Maximal Reliability) With a Confidence Interval

Estimates coefficient *H* (Hancock & Mueller, 2001), the reliability of
the optimally weighted composite of the items that measure one factor,
and returns a confidence interval for the population coefficient. *H* is
the largest reliability that any weighted sum of the items can reach
under the factor model, whereas McDonald's \\\omega\\
([`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md))
is the reliability of the unit-weighted sum that people usually score.
The coefficient can be computed from a fitted lavaan model (`fit`, the
recommended input), from raw item scores (`data`), from a covariance
matrix and its sample size (`S` and `N`), or from the standardized
loadings a paper reports (`loadings`). The interval comes from the delta
method with the full covariance matrix of the parameter estimates, built
on the logit scale by default, or from a nonparametric bootstrap when
one is requested.

## Usage

``` r
reliability_H(
  loadings = NULL,
  se_loadings = NULL,
  conf_level = 0.95,
  fit = NULL,
  factor = NULL,
  group = NULL,
  data = NULL,
  S = NULL,
  N = NULL,
  vcov_loadings = NULL,
  ci_method = c("mlr_logistic", "mlr", "ml_logistic", "ml", "delta_logistic", "delta",
    "bootstrap_se", "bootstrap_se_logistic", "percentile", "bca", "none"),
  missing = c("fiml", "listwise"),
  aux = NULL,
  B = 10000,
  seed = NULL,
  ncpus = 1L,
  progress = NULL
)
```

## Arguments

- loadings:

  Numeric vector of the standardized loadings of one factor, each in
  \\(-1, 1)\\, at least two. A fitted lavaan model supplied in this
  position is taken as `fit`, so `reliability_H(fit)` works.

- se_loadings:

  Optional vector of the marginal standard errors of the standardized
  loadings, the same length as `loadings`. The interval it gives ignores
  the covariances among the loading estimates and warns that it does;
  prefer `vcov_loadings` or `fit`.

- conf_level:

  Confidence level for the interval. Defaults to `0.95`.

- fit:

  A fitted lavaan model (from
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html),
  [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html), or
  `fa_confirmatory(output = "fit")`), the recommended input; see
  *Details*.

- factor:

  With `fit`, the name of the factor whose *H* is wanted; required when
  the model has more than one factor measured by observed items.

- group:

  With a multiple group `fit`, the label of the group (as
  `lavaan::lavInspect(fit, "group.label")` gives it) whose *H* is
  wanted; required for such a fit.

- data:

  A numeric matrix or data frame of item scores of one scale (rows are
  respondents, columns are items, three or more), and, when `aux` is
  supplied, the auxiliary columns. An item column that is not numeric
  (or logical), and an item that takes a single value in the rows
  analyzed, are refused by name before any fit, as is a zero variance on
  the diagonal of `S`: the single-factor model cannot be fit with such
  an item (lavaan would stop on it with an error that names no
  variable). SPSS user-defined missing values, which
  `haven::read_sav(user_na = TRUE)` keeps as values of a
  `haven_labelled_spss` column, are missing, as SPSS treats them.

- S:

  A symmetric covariance matrix among the items (three or more); `N` is
  then required, and only `"ml_logistic"`, `"ml"`, and `"none"` apply.

- N:

  Sample size. Required with `S`; optional with `loadings`, where it is
  recorded in the table; not used with `fit` or `data`, which count
  their own cases.

- vcov_loadings:

  With `loadings`, the \\J \times J\\ covariance matrix of the
  standardized loading estimates, for the delta method interval with the
  full covariance matrix.

- ci_method:

  Method for the confidence interval. `"delta_logistic"` and `"delta"`
  use the delta method with the covariance matrix that comes with the
  input (the fit's, or `vcov_loadings`); `"mlr_logistic"`, `"mlr"`,
  `"ml_logistic"`, and `"ml"` name the standard error of the
  single-factor model fitted to `data` or `S` (robust maximum likelihood
  or maximum likelihood); `"percentile"`, `"bca"`, `"bootstrap_se"`, and
  `"bootstrap_se_logistic"` are nonparametric bootstraps (with `data`,
  or with a single-group `fit` that carries its raw data); `"none"`
  returns the point estimate. The `_logistic` methods build the interval
  on the logit scale. When not supplied, the default is `"mlr_logistic"`
  with `data`, `"ml_logistic"` with `S`, and `"delta_logistic"` with
  `fit` and with `loadings` that come with `vcov_loadings` or
  `se_loadings`; for `loadings` without either it is `"none"`, since
  nothing is then known about their sampling variability. A method that
  does not apply to the input stops with a message naming the ones that
  do.

- missing:

  With `data`, how incomplete rows are handled: `"fiml"` (the default;
  full information maximum likelihood) or `"listwise"`. The default
  analysis is FIML whenever an item value is missing and the complete
  data analysis otherwise (data in which no value is missing, or none
  but in rows with every item missing, which both treatments drop, are
  fitted as they stand), as in
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  and unlike
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md);
  see *Details*. Under listwise deletion the single-factor model needs
  more complete rows than items (with no more, their covariance matrix
  is singular and the maximum likelihood fit is not defined), and the
  function stops before the fit when there are not. Ordinary FIML
  instead evaluates the constrained measurement model. A failure of the
  unrestricted H1 diagnostic does not alone invalidate its proper
  identified estimate or H0-based uncertainty; the result carries an
  explicit qualification. Supplied two-stage fits require valid
  unrestricted first-stage moments even for a point estimate. See
  [`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md)
  for the sufficient affine-collapse diagnostic and its limits.
  Bootstrap resamples retain the family's separate distinct-case screen.

- aux:

  With `data`, an optional character vector naming auxiliary variable
  columns of `data`, entered as saturated correlates under full
  information maximum likelihood; supplying it implies
  `missing = "fiml"`.

- B:

  Number of bootstrap replications when a bootstrap method is selected,
  a single whole number. Defaults to `10000`, the reliability family's
  default; every replication refits a model (see *Details* for the
  cost). The fewest accepted depends on the interval method and the
  level: at the default `conf_level`, 400 for `"percentile"` (10
  replications beyond each limit), 2000 for `"bca"` (50 beyond each
  nominal limit, because the adjusted limits sit farther into the
  tails), and 100 for `"bootstrap_se"` and `"bootstrap_se_logistic"` (a
  standard deviation of the replications); at 0.99 the first two are
  2000 and 10000. A smaller `B` is refused before anything is computed
  (class `dmar_count_invalid`). The count is held to the number
  requested. Replications set aside afterward are counted and reported
  in a warning of class `dmar_bootstrap_replications_dropped`, which
  names the function and says to request more when those that return a
  value are fewer than the number required; when fewer than half of the
  `B` replications return a value, the interval is withheld (limits and
  bootstrap standard error `NA`) under a warning of class
  `dmar_bootstrap_insufficient`. See
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md),
  which also describes `options(dmar.bootstrap_exploratory = TRUE)` for
  a first look with a smaller `B`.

- seed:

  Random number seed for the bootstrap. Defaults to `NULL`, under which
  the resamples are drawn from, and advance, the current random number
  stream; a supplied seed is set for the call, and the caller's RNG
  state is restored when it returns.

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
  A replication whose solution is improper is always set aside, since
  *H* is not defined for it. Ignored unless a bootstrap method is
  selected.

  Bootstrap intervals from supplied two-stage fits currently refuse
  separately stored auxiliary variables: each refit must resample and
  carry those values with the items. Use an eligible delta interval
  where available, or explicitly resample and refit them together.

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

The reliability family's table (class
`c("dmar_reliability", "dmar_tbl", "data.frame")`), with columns `term`
and `value` and rows `"estimate"` (*H*), `"se"` (its standard error:
delta method, or the standard deviation of the bootstrap replicates;
`NA` for `ci_method = "none"`), `"ci_lower"` and `"ci_upper"`,
`"conf_level"`, `"N"` (the cases the fit used, in the named group of a
multiple group fit; `NA` for `loadings` without `N`), `"N_complete"`
(the cases complete on the factor's items), and `"n_items"` (the number
of items). Attributes record the computation: `coefficient` (`"H"`),
`ci_method` (the method used), `missing` (`"listwise"`, `"fiml"`, the
fit's own label such as `"pairwise"`, or `"unknown"` for `loadings`),
`vcov_basis` (`"full"` when the delta method used the full covariance
matrix of the estimates, `"diagonal"` when only marginal standard errors
were supplied, `"none"` when no delta method standard error was formed),
`se_type` (the lavaan standard error type of the fit; `NA` for
`loadings`), `factor` and `group` (with `fit`; `NA` otherwise), `metric`
(`"observed"`, `"latent_response"` for ordered items, or `"unknown"` for
`loadings`), `residual_covariances` (the number of residual covariances
inside the block that entered the weights), and `loadings` (the
standardized loadings of the block, named by item). With `loadings`,
`vcov_supplied` records what came with them, whether or not an interval
was formed: `"full"` for `vcov_loadings`, `"diagonal"` for
`se_loadings`, and `"none"` for neither, in which case the printed table
says what an interval would need. With `data`, `aux` is recorded when
supplied. A bootstrap call also records `B` and the family's `bootstrap`
attribute (the replications attempted, returned, and set aside, with the
reasons; see
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)),
which here also carries `bias`, the mean of the returned replicates
minus the estimate (`NA` when no interval is formed). Its `no_solution`
counts the resamples of three items whose correlations admit no
single-factor solution, which are among the `improper` ones (on the data
path; zero with `fit`, whose refits are not checked for it), and
`improper_kept` is zero, since *H* is not defined for an improper
solution.

## Details

**The coefficient.** For a factor \\f\\ measured by items \\b\\ that
load on it alone, with unstandardized loadings \\\lambda\\, residual
covariance matrix \\\Theta\_{bb}\\, and model implied factor variance
\\\phi_f\\, \$\$q = \phi_f \\ \lambda' \Theta\_{bb}^{-1} \lambda, \qquad
H = \frac{q}{1 + q}.\$\$ With uncorrelated residuals and standardized
loadings \\\ell_j\\ this is the form Hancock and Mueller (2001) gave,
\$\$H = \frac{\sum_j \ell_j^2 / (1 - \ell_j^2)} {1 + \sum_j \ell_j^2 /
(1 - \ell_j^2)},\$\$ in which each item contributes its ratio of common
to unique variance, \\\ell_j^2 / (1 - \ell_j^2)\\. The weights that
reach *H* are proportional to \\\Theta\_{bb}^{-1} \lambda\\, so an item
counts in proportion to its loading divided by its residual variance,
and adding an item can never lower *H*. Gabler and Raykov (2017, Eq. 9)
showed that the same expression with a nondiagonal \\\Theta\_{bb}\\ is
the maximal reliability when residuals inside the block are correlated;
a fitted model with residual covariances among the factor's items gets
that coefficient here, and the printed table says how many entered the
weights. Rescaling the factor rescales \\\lambda\\ and \\\phi_f\\
inversely, so *H* does not depend on whether the factor was identified
by a marker item or by fixing its variance.

Two further readings of the number are useful. *H* is the squared
correlation between the factor and its best linear predictor from the
block's items, the regression factor score computed from those items
alone, so it is the reliability of that score: for a one-factor model
the determinacy that
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md)
reports is the square root of *H*. And because *H* is the largest
reliability any weighting of the items can reach under the fitted model,
it is at least the model implied \\\omega\\ of the unit-weighted sum
from the same fit (Raykov, 1997), with equality when the weights
\\\Theta\_{bb}^{-1} \lambda\\ are all equal, which with uncorrelated
residuals means equal standardized loadings. Robust omega, whose
denominator is the observed variance of the composite, need not lie
below it.

**Reading *H* with care.** The sample *H* is a positively biased
estimator of the population coefficient, and a composite formed with
weights estimated from the sample is on average less reliable than the
population *H*. Aguirre-Urreta, Rönkkö, and McIntosh (2019) found both
effects pronounced in small samples (around *N* = 100), diminishing as
the sample grows and barely noticeable beyond about 1,000, although with
three items of modest loadings even 1,000 cases left more than 1% bias;
the reliability estimate of the unit-weighted composite was nearly
unbiased in the same conditions. Report \\\omega\\, with its interval,
for a scale that people will score by summing or averaging its items,
and *H* as the ceiling that optimal weights could reach, or as the
reliability of scores computed from the factor model. A bootstrap call
estimates the bias of the sample *H* (the `bias` element of the
`bootstrap` attribute), and the printed table shows it.

**The four inputs.** Exactly one of `fit`, `data`, `S`, and `loadings`
is supplied.

- `fit`:

  A fitted lavaan model, from
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html),
  [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html), or
  `fa_confirmatory(output = "fit")`; the recommended input, because it
  carries everything the coefficient and its interval need, and the
  model can be the one the analysis reports. Any factor measured by
  observed items can be named with `factor` (a factor of a confirmatory
  model, or an exogenous or endogenous factor of a structural model),
  and any group of a multiple group model with `group`. *H* is read from
  the fitted matrices, and its standard error comes from the fit's own
  covariance matrix of the estimates, so it is the standard, robust, or
  FIML standard error the model was fitted with. With ordered items
  fitted by WLSMV (either parameterization) *H* is the maximal
  reliability of the latent responses behind the categories, not of the
  observed categories; the `metric` attribute and the printed table say
  so. A loading or a regression the model holds at zero is not one,
  whether it is fixed (`0*item`) or constrained (`a*item` with `a == 0`,
  directly or through a chain of equality constraints between labels, as
  in `a == b` with `b == 0`): the model written with that term is the
  model written without it and gets the same *H*, so an item listed
  under the factor with its loading held at zero is not one of the
  factor's items, and a cross-loading held at zero refuses nothing. A
  fit is refused, with a message naming the remedy, when an item of the
  factor also loads on another factor or is regressed on another
  variable (for the factors of a bifactor model,
  [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
  reports the index that Rodriguez, Reise, and Haviland, 2016, call *H*,
  which counts the other factors' common variance as uncorrelated error
  and is not this maximal reliability, and, in its
  `"maximal_reliability"` attribute, the maximal reliability the fitted
  model implies for each factor, and
  [`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
  gives both a confidence interval), when the factor is measured by
  other factors, when an observed covariate predicts the factor, when
  the solution is improper (an item's residual variance that is not
  positive, for which *H* is not defined), and when a factor or a group
  must be named and is not.

- `data`:

  Raw scores on the items of one scale. The single-factor congeneric
  model is fitted exactly as
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  fits it, so *H* and \\\omega\\ rest on the same fit, the same cases,
  and the same kind of standard error, and `missing` and `aux` work as
  they do there. Three or more items are required, since the model is
  not identified with two.

- `S`:

  A covariance matrix among the items, with its sample size `N`; the
  same single-factor model, fitted by maximum likelihood.

- `loadings`:

  The standardized loadings of one factor, as a paper reports them, with
  the closed form above (which assumes uncorrelated residuals). For an
  interval, supply `vcov_loadings`, the covariance matrix of the
  standardized loading estimates (for example the matching block of
  `lavaan::lavInspect(fit, "vcov.std.all")`); `se_loadings`, their
  marginal standard errors, is accepted, but it ignores the covariances
  and warns that it does. Loadings alone carry nothing an interval can
  be formed from, so without either the estimate is returned with `se`
  and the limits `NA`, and the printed table says what an interval would
  need.

**The interval.** The delta method variance of \\\hat H\\ is \\g' V g\\,
where \\V\\ is the covariance matrix of all the free parameter estimates
of the fit and \\g\\ is the gradient of *H* with respect to them,
computed by central differences of the closed form above. \\V\\ is used
only when it can be, by the check the factor analysis tables
([`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md))
make before they report a standard error: when lavaan could not compute
it, or computed it from an information matrix, observed or expected,
that is singular in the directions the fit's equality constraints leave
free (scaled to a unit diagonal, its smallest eigenvalue is below
\\10^{-6}\\), a parameter of the fit is not identified at the estimates
and \\g' V g\\ means nothing. The estimate is then returned with `se`
and the limits `NA`, and with the warning, of class
`dmar_result_qualified`, that those tables give for the same fit: the
same words for a singular matrix, and for one lavaan could not compute a
sentence that ends by naming the standard error and interval of *H*
rather than the tables' columns. The check is of the whole fit, as the
tables make it, so a factor whose own parameters are identified loses
its interval beside one whose parameters are not. An interval from
`se_loadings` uses only the marginal standard errors of the standardized
loadings, the diagonal of \\V\\. The loading estimates of a fitted
factor model are correlated, and the diagonal can understate or
overstate the variance: in seven blocks of three to six tests from
[`holzinger_swineford`](https://yelleknek.github.io/DMAR/reference/holzinger_swineford.md),
each fitted as one factor, the standard error from the full matrix
ranged from about 0.7 to 1.3 times the one from the diagonal. By default
the Wald interval is built on the logit scale and transformed back
(Browne, 1982): the limits of a monotone transformation of a parameter
are the transformation of its limits, the confidence interval
transformation principle that Kelley (2007) sets out, so the interval
stays inside (0, 1) and is asymmetric near the upper bound, as the
sampling distribution of a high reliability is; `"delta"`, `"ml"`, and
`"mlr"` give the untransformed Wald interval. The default standard error
depends on the input: with `data`, the robust maximum likelihood
standard error (Satorra-Bentler, or Yuan-Bentler under FIML), as for
\\\omega\\, because Kelley and Pornprasertmanit (2016) found the normal
theory Wald interval wanting when items depart from normality; with `S`,
the maximum likelihood standard error; and with `fit`, the fit's own.
With `fit`, `"ml"` and `"mlr"` (and their `_logistic` forms) are
accepted when they name the standard error the fit already carries, and
otherwise stop with a message saying which one it carries.

The delta method interval is the fast default. For a reported analysis,
the percentile or BCa bootstrap with the model refitted on every
replication (by FIML when data are missing) is the safer choice, as it
is for \\\omega\\ (Kelley & Pornprasertmanit, 2016), and it also
estimates the bias of the sample *H*. With `data` the rows are resampled
and the single-factor model refitted; with `fit` the rows the model was
fitted to are resampled and the model itself refitted through
[`lavaan::lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html), with
the fitted estimates as starting values (for a single-group fit that
carries its raw data). Before any resample is drawn, the model is
refitted this way to the rows of `fit` itself, and the call stops unless
that refit reproduces the estimate, so the replications are refits of
the model the estimate describes. A fit with sampling weights is
refused, since a resample of its rows does not carry them and every
replication would be an unweighted fit, and so is a fit with a `cluster`
variable, whose cases are not independent and cannot be resampled one at
a time; the delta method interval, built from the fit's own covariance
matrix of the estimates, serves both. The bootstrap engine is the one
the whole reliability family uses (Efron & Tibshirani, 1993):
replications that fail to converge, give an improper solution (or, with
three items, a resample whose correlations admit no single-factor
solution), or raise an error are set aside and counted (a resample in
which an item takes a single value is refused before the refit, by name,
as the data would be; lavaan would stop on the item after printing its
table of the variables; so is a resample with too few distinct cases, by
the rule of `missing` with `data`, and with `fit` by the rule of the
fit's own treatment of incomplete rows: under listwise deletion, or with
complete data, one with no more distinct rows than the model's
continuous variables, whose sample covariance matrix is singular, and
under FIML or two-stage estimation one with no more distinct cases that
observe a variable, while a fit by pairwise deletion is held to
neither), and the counts travel in the `bootstrap` attribute, print
under the table, and are stated in a warning. A replication whose
likelihood has no maximum (an error variance running off below zero
while the likelihood is flat or still rising; see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md))
is counted among those that raised an error, and a fit with no maximum
behind the estimate itself stops with an error of class
`dmar_no_maximum`, which with `data` prints the lavaan code of the
single-factor fit, where `bounds = "pos.var"` or a different model are
the choices. *H* is not defined for an improper solution, so the
attribute records `improper_rule = "drop"`. Each replication is a model
fit, a few hundredths of a second for a one-factor model on complete
data and longer under FIML, so the default `B = 10000` (the family's
default, chosen for the stability of the BCa quantiles) takes minutes.
The examples therefore stop at the delta method; the bootstrap call for
their five verbal tests is
`reliability_H(data = holzinger_swineford[, verbal], ci_method = "bca", B = 2000, seed = 113)`.

**Missing data.** With `data`, incomplete items are fitted by default by
full information maximum likelihood (FIML), which keeps every case with
at least one observed item or observed auxiliary variable, with
auxiliary variables (`aux`) entered as saturated correlates (Graham,
2003), and `missing = "listwise"` fits the complete rows alone; data in
which no value is missing, and `S`, are fitted as they stand. A row with
every item missing does not count as incomplete data, since without
auxiliary variables both treatments drop it. The default analysis is
therefore FIML whenever an item value is missing and the complete data
analysis otherwise, unlike
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
which fits complete data by FIML too. That is the rule of
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
whose complete-data intervals need the complete data analysis, and it is
followed here so that the two fit the same single-factor model the same
way; on incomplete data it agrees with
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).
With `fit`, the fit's own treatment applies: a model fitted with
`missing = "ml"` gives *H* and its standard error by FIML. Either way
`N` counts the cases the fit used and `N_complete` those complete on the
factor's items, so the cost of listwise deletion is visible. Under
`missing = "fiml"` a pair of items or auxiliary variables that no case
observes together is refused by name before any fit, since the data
carry no information about its covariance (lavaan only warns, and its EM
returns an arbitrary value for it), and so is a pair the cases observe
together only where one of them takes a single value, since those cases
carry none about how the two vary together; in a bootstrap, a resample
that leaves such a pair is set aside with the same message and counted.

**From standardized loadings.** With `loadings` alone the function
returns the estimate, with `se` and the limits `NA`, because
standardized loadings carry nothing an interval can be formed from (the
default `ci_method` is then `"none"`). With `vcov_loadings` the interval
is the delta method interval from the full covariance matrix of the
loading estimates. With `se_loadings` it is the delta method interval
from their marginal standard errors, and the call warns every time, with
class `dmar_H_marginal_se`, that the interval ignores the covariances
among the loading estimates. Every input returns the reliability family
table described under *Value*; the sampling variance of *H* is the
square of the row `se`. `reliability(type = "H")` dispatches here.

The coefficient that
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
reports as `H_<factor>` is the same quantity, with lavaan's delta method
standard error and an untransformed interval;
`reliability_H(fit, factor = "speed", ci_method = "delta")` reproduces
its estimate and standard error for the factor `speed`.

## References

Aguirre-Urreta, M. I., Rönkkö, M., & McIntosh, C. N. (2019). A
cautionary note on the finite sample behavior of maximal reliability.
*Psychological Methods, 24*(2), 236–252.
[doi:10.1037/met0000176](https://doi.org/10.1037/met0000176)

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
*Topics in applied multivariate analysis* (pp. 72–141). Cambridge, UK:
Cambridge University Press.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Gabler, S., & Raykov, T. (2017). Evaluation of maximal reliability for
unidimensional measuring instruments with correlated errors. *Structural
Equation Modeling, 24*(1), 104–111.
[doi:10.1080/10705511.2016.1159916](https://doi.org/10.1080/10705511.2016.1159916)

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Raykov, T. (1997). Estimation of composite reliability for congeneric
measures. *Applied Psychological Measurement, 21*(2), 173–184.
[doi:10.1177/01466216970212006](https://doi.org/10.1177/01466216970212006)

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
bifactor models: Calculating and interpreting statistical indices.
*Psychological Methods, 21*(2), 137–150.
[doi:10.1037/met0000045](https://doi.org/10.1037/met0000045)

## See also

[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for the reliability of the unit-weighted composite,
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
for the family's general entry point,
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md)
for the scores whose reliability *H* is,
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
for the fit, and
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
for the bifactor index of Rodriguez, Reise, and Haviland (2016) and the
maximal reliability of each factor of a bifactor model.

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
verbal <- c("t5_general_information", "t6_paragraph_comprehension",
            "t7_sentence", "t8_word_classification", "t9_word_meaning")

# From raw item scores: the single-factor model is fitted here, and the
# interval uses its robust maximum likelihood standard error on the
# logit scale.
reliability_H(data = holzinger_swineford[, verbal])
#>  term       value
#>  estimate   0.918
#>  se         0.007
#>  ci_lower   0.903
#>  ci_upper   0.931
#>  conf_level 0.950
#>  N          301  
#>  N_complete 301  
#>  n_items    5    

# The omega of the unit-weighted sum from the same model, for
# comparison; H is the ceiling that optimal weights could reach.
reliability_omega(data = holzinger_swineford[, verbal],
                  denominator = "model_implied")
#>  term       value
#>  estimate   0.902
#>  se         0.010
#>  ci_lower   0.881
#>  ci_upper   0.919
#>  conf_level 0.950
#>  N          301  
#>  N_complete 301  
#>  n_items    5    
#> 
#> Denominator: composite variance implied by the fitted model
#> Largest residual correlation under the single-factor model: -0.040 (t8_word_classification with t9_word_meaning)

# From a fitted model with several factors, name the factor. The
# standard error is the fit's own, from the full covariance matrix of
# its estimates (requires lavaan).
fit <- lavaan::cfa(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
   textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
   speed   =~ t10_addition + t12_counting_groups_of_dots +
              t13_straight_and_curved_capitals",
  data = holzinger_swineford, std.lv = TRUE)
h <- reliability_H(fit, factor = "textual")
h
#>  term       value
#>  estimate   0.885
#>  se         0.011
#>  ci_lower   0.861
#>  ci_upper   0.906
#>  conf_level 0.950
#>  N          301  
#>  N_complete 301  
#>  n_items    3    
#> 
#> Factor: textual
attr(h, "loadings")
#> t6_paragraph_comprehension                t7_sentence 
#>                  0.8515822                  0.8550655 
#>            t9_word_meaning 
#>                  0.8380103 

# From the standardized loadings a paper reports, without an interval.
reliability_H(loadings = c(0.6, 0.7, 0.8), ci_method = "none")
#>  term       value
#>  estimate   0.767
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          <NA> 
#>  N_complete <NA> 
#>  n_items    3    
#> 
#> No interval (ci_method = "none"): an interval from standardized loadings needs 'vcov_loadings', the covariance matrix of their estimates

# The broom verbs.
generics::tidy(h)
#>   term  estimate        se  ci_lower  ci_upper conf_level wald_df h1_diagnostic
#> 1    H 0.8852807 0.0114807 0.8607882 0.9059349       0.95     Inf          <NA>
#>   wald_calibration_note
#> 1                  <NA>
generics::glance(h)
#>   coefficient  estimate        se  ci_lower  ci_upper conf_level nobs n_items
#> 1           H 0.8852807 0.0114807 0.8607882 0.9059349       0.95  301       3
#>        ci_method wald_df h1_diagnostic wald_calibration_note
#> 1 delta_logistic     Inf          <NA>                  <NA>
```
