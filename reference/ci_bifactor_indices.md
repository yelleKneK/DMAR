# Bootstrap Confidence Intervals for the Bifactor Indices

Forms nonparametric bootstrap confidence intervals for the indices
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
reports for a fitted bifactor model: the explained common variance
(ECV), coefficient omega, omega hierarchical (omega_H) for the general
factor, omega hierarchical subscale (omega_HS) for each group factor,
the index Rodriguez, Reise, and Haviland (2016) call coefficient H, the
closed form of Hancock and Mueller (2001) applied to one factor's
loadings, which counts the other factors' common variance as
uncorrelated error, and the maximal reliability of each factor under the
fitted model, which the index is not (see
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)).
Each replication resamples the cases the fit used, refits the same model
with the same estimator, and recomputes every index, so the intervals
rest on the same replications (those of the maximal reliability, under
`improper = "keep"`, on the proper ones among them). Replications whose
solution is improper are set aside by default and counted, and the table
says how many.

## Usage

``` r
ci_bifactor_indices(
  fit,
  general = NULL,
  denominator = c("model_implied", "observed"),
  ci_method = c("percentile", "bca"),
  conf_level = 0.95,
  B = 1000L,
  seed = NULL,
  ncpus = 1L,
  improper = c("drop", "keep"),
  progress = NULL
)

# S3 method for class 'dmar_ci_bifactor'
tidy(x, ...)

# S3 method for class 'dmar_ci_bifactor'
glance(x, ...)
```

## Arguments

- fit:

  A fitted bifactor lavaan model to raw data, as for
  [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md):
  one general factor on all items plus two or more orthogonal group
  factors, each item on one group factor at most (an item on none is on
  the general factor alone). A fit that
  [`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
  refuses is refused here, and so are a multiple-group fit, a fit with a
  cluster variable or sampling weights, a fit to a covariance matrix,
  and a fit that did not converge (see *Details*).

- general:

  Optional name of the general factor, as for
  [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md);
  by default the factor that loads on every item.

- denominator:

  How the total-score variance in the denominators of omega, omega_H,
  and omega_HS is estimated, as for
  [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md):
  `"model_implied"` (the default) or `"observed"`. Each replication uses
  the same setting.

- ci_method:

  The bootstrap interval: `"percentile"` (the default) or `"bca"`, the
  bias-corrected and accelerated interval. See *Details*.

- conf_level:

  Confidence level of the intervals. Defaults to `0.95`.

- B:

  Number of bootstrap replications, a single whole number. Defaults to
  `1000`, the package's count for a bootstrap whose every replication
  refits a model; see *Details* for the cost. The fewest accepted at the
  default `conf_level` is 400 for `"percentile"` (10 replications beyond
  each limit; 200 at 0.90, and 2000 at 0.99, which the default does not
  meet) and 2000 for `"bca"` (50 beyond each nominal limit), which the
  default does not meet either: pass `B = 2000` or more with
  `ci_method = "bca"`. A smaller `B` is refused before anything is
  computed (class `dmar_count_invalid`). The count is held to the number
  requested: replications set aside afterward are counted and reported
  in a warning of class `dmar_bootstrap_replications_dropped`, and an
  index for which fewer than half of the `B` replications return a value
  has its interval withheld (`NA`) under a warning of class
  `dmar_bootstrap_insufficient`. Where refits can fail, request more
  than the minimum. See
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md),
  which also describes `options(dmar.bootstrap_exploratory = TRUE)` for
  a first look with a smaller `B`.

- seed:

  Random number seed for the bootstrap. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  internally and the caller's state is restored on exit.

- ncpus:

  The number of processes over which the replications are spread. The
  default `1` runs them in this R session; a larger number uses forked
  workers (boot's `parallel = "multicore"`) or, on Windows, a socket
  cluster (`"snow"`), whose workers load the installed DMAR. boot draws
  every resample before any replication runs, so a given `seed` returns
  the same intervals, counts, and replicates whatever `ncpus` is. A
  session that forks workers and also starts background processes
  through callr or processx needs the environment variable
  `PROCESSX_NOTIFY_OLD_SIGCHLD` set before processx is loaded (see
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)).

- improper:

  How a replication whose solution is improper (a standardized loading
  beyond one or a negative residual variance, a Heywood case) is
  treated: `"drop"` (the default) sets it aside, and `"keep"` uses its
  indices when every one of them but the maximal reliability, which is
  not defined for an improper solution, is a proportion. Either way such
  replications are counted. With `"drop"` an improper point estimate
  stops the call. See *Details*.

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

- x:

  A `dmar_ci_bifactor` table returned by `ci_bifactor_indices()`.

- ...:

  Unused.

## Value

A `data.frame` (class `dmar_ci_bifactor`, `dmar_tbl`) with one row for
each index that applies to a factor, the general factor first (ECV,
omega, omega_H, PUC, H, maximal_reliability) and then each group factor
(ECV, omega, omega_HS, H, maximal_reliability; its ECV is the general
factor's share of the common variance of that factor's items, and its
omega the reliability of their total score, as defined under *Details*
of
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
which the printed table repeats), and columns `factor`, `index`,
`estimate` (the value
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
reports, in its table or, for `maximal_reliability`, in its
`"maximal_reliability"` attribute), `se` (the standard deviation of the
replications used), `ci_lower`, and `ci_upper`. PUC has `NA` in `se` and
the limits, and so does every column of the maximal_reliability rows
when the point estimate is improper, and so do the `se` and limits of
any index whose interval is withheld (see `require_half` under the
`bootstrap` attribute below). The rows depend on the model alone, and
the columns never change. Attributes:

- `bootstrap`:

  The counts and draws: `attempted` (`B`), `returned` (the replications
  used), `failed` (those set aside), `reasons`, `error_messages`,
  `undetermined` (always zero here: it counts the replications that left
  a polychoric correlation undetermined, which only
  [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
  meets), `improper` (improper replications set aside), `improper_kept`,
  `no_solution` (always zero here: it counts the resamples of three
  items whose correlations admit no single-factor solution, which only
  the single-factor coefficients of
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  and
  [`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
  meet), `not_converged`, `too_few_cases` (always zero here: it counts
  the replications the single-factor coefficients set aside for too few
  distinct cases), `improper_rule`, `caution` (`TRUE` when the
  replications set aside and the improper ones kept number at least
  \\(1 - C)/2\\ of `B`, with \\C\\ the confidence level, so that how
  they were treated can decide where a limit falls), `require_half`
  (`TRUE`: an index's interval requires at least half of all attempted
  replications to return a value for it; otherwise its `se`, `ci_lower`,
  and `ci_upper` are `NA`, and the warning and the printed line say so;
  this is an estimability safeguard, not a coverage guarantee),
  `required` (the fewest replications a request may ask for, 400 for the
  percentile interval and 2000 for BCa at the 95 percent level; the
  count is held to `B`, not to the replications that return a value) and
  `bootstrap_status` (`"complete"`, `"dropped"`, `"exploratory"`,
  `"insufficient_count"`, or `"insufficient_share"`, the last two when
  some index has its interval withheld; see
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)),
  `target_returned` and `target_withheld` (for each bootstrapped index,
  that is, each column of `t` and named as those columns, the
  replications that returned a value for it and whether its interval was
  withheld), `acceleration` and `acceleration_method` (one per index),
  `jackknife`, `ncpus`, `t` (a `B` by index matrix of the estimates of
  every replication, one column for each row with an interval, named
  `"factor:index"`, and `NA` where a replication gave none; the indices
  of an improper replication set aside under `"drop"` are kept here when
  they are proportions, and its maximal reliabilities are `NA`),
  `status` (the fate of every replication, a factor with levels
  `"proper"`, `"improper"`, `"error"`, `"non_finite"`, `"undetermined"`,
  and `"not_converged"`), and `bias` (one per column of `t`, named as
  its columns: the mean of the replicates the index's interval uses
  minus its estimate, the bootstrap estimate of its bias; see
  *Details*).

- `B`, `B_used`:

  The replications requested and the smallest number any reported
  interval rests on (under `improper = "keep"`, the intervals of the
  maximal reliability rest on the proper replications, the count of
  finite values in a maximal reliability column of `t`, which can be
  fewer than `bootstrap$returned`).

- `n_dropped`, `B_required`, `bootstrap_status`:

  `B` minus `B_used`, the fewest replications a request may ask for, and
  the verdict, in the names every bootstrap result of the package
  carries.

- `ci_method`:

  The interval of each index, named `"factor:index"`: the method
  requested, or the percentile label with the reason when a BCa interval
  could not be formed.

- `ci_method_requested`, `conf_level`, `improper_rule`, `denominator`:

  The settings of the call.

- `metric`:

  As for
  [`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md):
  `"observed"` for continuous items and `"latent_response"` for ordered
  items.

- `general`:

  The name of the general factor.

- `improper`:

  Whether the point estimate is improper (possible only with
  `improper = "keep"`).

- `N`, `N_incomplete`:

  The rows resampled and how many of them have a missing value.

[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
returns one row per index with `term` (`"factor:index"`), `estimate`,
`se`, `ci_lower`, `ci_upper`, and `conf_level`;
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns one row with the counts and the settings.

## Details

**Why intervals.** The bifactor indices are estimates, and the decisions
they are read for (whether a total score is close enough to
unidimensional to report alone, whether a subscale carries reliable
variance beyond the general factor) are decisions about the population
values. A confidence interval says how precisely the sample has located
each of them. In the bifactor model of 15 Holzinger and Swineford tests
described under *Improper replications* below, omega_H is 0.577 with a
95 percent interval from 0.468 to 0.672 (`B = 1000`, `seed = 113`): 301
cases locate the share of total score variance due to the general factor
only to within an interval about 0.2 wide.
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
reports the point estimates alone. The intervals are a separate
function, rather than an argument of
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
because they need columns that its one-row-per-factor table cannot gain
without its shape depending on an argument.

**What is resampled.** The cases are the rows of data lavaan stored with
`fit`, which are the rows the fit used: under listwise deletion (the
default, and lavaan's default for ordered items) the complete rows, and
under `missing = "ml"` (full information maximum likelihood) every row
with at least one observed item, each with its own pattern of observed
and missing values, so that a resample has the missing data the sample
has and each replication is again a full information fit. A row with
nothing observed, which lavaan keeps in its stored data but leaves out
of the sample, is left out here too, and the number of rows resampled is
checked against the number of cases lavaan reports. Rows are drawn with
replacement, `B` times (Efron & Tibshirani, 1993).

**How each replication is fit.** Each resample is fit with the parameter
table and the options of `fit`: the same model, the same estimator, the
same missing data treatment, the same ordered items, and default
starting values, as for `fit` itself. Standard errors and test
statistics are not computed in a replication, since no index uses them;
robust maximum likelihood (MLR) and WLSMV have the estimates of maximum
likelihood and diagonally weighted least squares and differ from them
only in those, so a replication is the cheaper fit with the same
estimates. Before any resample is drawn, the model is refit this way to
the rows of `fit` itself, and the call stops unless that refit
reproduces the point estimates, so the replications are refits of the
model the point estimates describe. The indices of each replication are
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
of its refit, with the general factor named, so every index is computed
as the point estimate is.

**Improper replications.** A bifactor model is prone to improper
solutions, a standardized loading past one or a residual variance below
zero, which can come from sampling error around a small residual
variance, from a model that does not describe the items, or from a group
factor the data barely separate from the general factor, and a resample
can give such a solution where the sample does not. A replication is
improper when
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
flags its solution or lavaan's own check after estimation fails, and the
point estimate is judged the same way. For a fit with auxiliary
variables (lavaan's `aux` under full information maximum likelihood)
that check is made on the model's own parameters, as
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
explains: lavaan's own check also tests the residual covariances of the
auxiliary variables, which are not positive definite once one of them is
related to the general factor strongly enough, however proper the
bifactor solution. A bifactor model of 15 of the Holzinger and Swineford
tests (`holzinger_swineford`: a general factor on all of them, with
spatial, verbal, and memory group factors on tests 1–4, 5–9, and 14–19)
gave an improper solution in 212 of 1,000 resamples (`seed = 113`) and
did not converge in 17 more. The omega of the improper solutions ranged
from 0.84 to 26.5, 132 of them above one; kept, all of them put the
upper limit of the 95 percent interval for omega at 13.7 and the lower
limit for ECV at 0.004, where the proper replications give \[0.855,
0.902\] and \[0.309, 0.509\]. Such values are not sampling variability
around a proper structure: a solution in which a residual variance has
run away, or in which the general factor has all but vanished (an ECV
near zero), is not a perturbation of the sample's solution but a
different one, and its indices do not estimate the quantities the point
estimate does. The default, `improper = "drop"`, therefore sets them
aside. This differs from the default of the reliability family
([`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)),
where an improper single-factor replication that still yields a
coefficient is typically a loading that ran just past one, and keeping
it preserves the upper tail of the bootstrap distribution. Setting
replications aside conditions the intervals on admissible solutions:
they describe the sampling distribution of the indices among resamples
that the model fits properly. Yuan and Hayashi (2003, Section 6) argue
that a replication that fails because its resample lies far from the
model is informative and should be counted rather than discarded; for a
test statistic they count it as an extreme value, but an index has no
side on which such a replication belongs, so here it is counted and
reported rather than placed. Every replication set aside, for any
reason, is counted in the `bootstrap` attribute, named in a warning, and
stated under the printed table, and when the replications set aside (or
the improper ones kept) reach the tail area of an interval, \\(1 -
C)/2\\ of `B` for a confidence level \\C\\, the warning and the printed
line add a caution that the limits depend on how they were treated. When
fewer than half of the `B` replications, or fewer than the number the
interval method and the level require (400 for the percentile interval
and 2000 for BCa at the 95 percent level; see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)),
return a value for an index, its interval is withheld, and the warning
(class `dmar_bootstrap_insufficient`) and the printed line say so, as in
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md).
A run at exactly the required `B` therefore loses its intervals to a
single replication set aside; the default `B = 1000` leaves room for
them at the 95 percent level. With `improper = "keep"` an improper
replication is used when every index it gives but the maximal
reliability, which is not defined for an improper solution (see *The
maximal reliability* below), is a proportion; one whose solution gives
an index outside \[0, 1\] (an omega above one belongs to a composite
whose implied error variance is negative) or an H that is not defined (a
standardized loading at or past one leaves its item no error variance in
H's formula) yields no coefficient and is set aside under either
setting, as in the reliability family. So is a replication whose
likelihood has no maximum, a residual variance running off below zero
while the likelihood is flat or still rising (see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)):
it is counted among the replications that raised an error, under a
message that names the item, and the run goes on (of the 22 improper
replications in 100 of the 15-test model above, `seed = 113`, 8 were of
this kind), while a refit of the data stored with `fit` that has no
maximum stops the call with an error of class `dmar_no_maximum`. Keeping
the admissible indices of improper replications can change their
interval limits; the maximal reliability uses only proper replications
under either setting. With `"drop"` an improper point estimate stops the
call, since an interval built from proper replications would not
describe it; with `"keep"` it is reported with a warning, unless one of
its indices other than the maximal reliability is not defined or lies
outside \[0, 1\]; its maximal reliability, which is not defined, is
reported as `NA`.

**The maximal reliability.** The rows `maximal_reliability` give each
factor's maximal reliability under the fitted model, \\\lambda_f' R^{-1}
\lambda_f\\ with \\R\\ the item correlations the model reproduces (the
`"maximal_reliability"` attribute of
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)):
the reliability the optimally weighted items reach as a measure of the
factor, and the square of the determinacy of its regression factor
score. Like every other index, it is read from
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
on each refit, so its interval costs nothing beyond the replications the
other indices need. It is a reliability estimate in its own right, as
much in need of an interval as omega (Raykov & Penev, 2006, give one for
a composite of congeneric measures), and it is the quantity to read
beside the index H, which can lie above or below it. It is not defined
for an improper solution. With `improper = "drop"` every replication
used is proper, and its intervals rest on the same replications as the
others; with `"keep"` the improper replications kept for the other
indices have no maximal reliability, so its intervals rest on the proper
replications among them, which are the replications `"drop"` uses, and
the printed table says how many. Its interval is therefore the same
under either rule. An improper point estimate, reported under `"keep"`,
has no maximal reliability, and these rows then carry no estimate and no
interval.

**The bias of H and the maximal reliability.** The sample H, the closed
form the index applies to a factor's loadings, is a positively biased
estimator of its population value, and a composite weighted with sample
estimates is on average less reliable than the population value:
Aguirre-Urreta, Rönkkö, and McIntosh (2019) found both effects
pronounced in samples of about 100 and slight beyond about 1,000, for H
of a single factor. The maximal reliability is the reliability of
weights that are optimal for the sample's estimates and is biased upward
in the same way, most of all for a group factor that the data determine
weakly. Run with `B = 200` in place of its 10, the example below gives a
mean of the replicates that exceeds the estimate by 0.031 for the first
group factor's maximal reliability (whose replicates have a standard
deviation of 0.078) and by 0.019 for its H, against 0.005 for its omega
and 0.003 for its omega_HS. The percentile interval takes its limits
from the replicates as they fall, and so carries the bias; the BCa
interval (`ci_method = "bca"`) moves its limits to correct for the
median bias of the replicates. The bias of each index with an interval,
the mean of the replicates its interval uses minus the estimate, is the
`bias` element of the `bootstrap` attribute (PUC, which has no interval,
has none, and neither has the maximal reliability of an improper point
estimate under `improper = "keep"`), and the printed table gives it for
H and the maximal reliability. Read the maximal reliability as a ceiling
that weights estimated from these data will not quite reach.

**Other replications set aside.** A refit that does not converge is set
aside and counted as such, and so is one that raises an error; the
warning quotes the most frequent error. For ordered items, a resample
that draws no response in some category of an item cannot refit that
item's thresholds, and is set aside as an error that says so. The
intervals are formed from the replications that remain, and the
effective count is the `B_used` attribute.

**The two intervals.** The percentile interval takes its limits from the
empirical quantiles of the bootstrap estimates; it respects the \[0, 1\]
range of the indices and is not forced to be symmetric. Kelley and
Pornprasertmanit (2016) found it the most consistently accurate of the
intervals they compared for the single-factor omega with an observed
denominator (their hierarchical omega, which is not the omega_H here),
which is the reason it is the default; the sources cited here do not
evaluate the coverage of bootstrap intervals for the bifactor indices
themselves. The bias-corrected and accelerated (BCa) interval moves the
two quantile positions to correct for median bias and for the rate at
which the variance of the estimator changes with the parameter, which
makes it second-order accurate where the percentile interval is
first-order accurate (DiCiccio & Efron, 1996). Its acceleration is
estimated as boot estimates it, by regression on the resampling
frequencies when at least as many replications as cases were used and
otherwise, or when that regression fails, by the jackknife, one
leave-one-out refit per case, which is run once for all of the indices;
an index whose BCa interval cannot be formed gets the percentile
interval, with a warning, and the `ci_method` attribute and the printed
line say which. PUC, the proportion of item pairs whose items lie on
different group factors, depends only on the model's structure, so it is
the same in every replication and has no interval; its row carries the
point value with `NA` limits.

**Cost.** Each replication is a model fit, so the call costs `B` fits
and, for `"bca"` with fewer replications used than cases, as many more
as there are cases. A proper fit of the 15-test Holzinger and Swineford
model takes a few hundredths of a second, but a solution that runs away
took seconds before the optimizer stopped, so a model that gives many
improper replications is slow to bootstrap as well as hard to interpret.
`ncpus` spreads the replications over several processes, and `seed`
makes a call reproducible. Requires the lavaan and boot packages.

**What is refused.** A fit that
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
refuses (correlated factors, no factor that loads on every item, fewer
than two group factors, an item on more than one) is refused, with the
same message. So is a multiple-group fit, since resampling the pooled
rows would let the group sizes vary from replication to replication (fit
each group separately); a fit with a `cluster` variable, whose cases are
not independent and cannot be resampled one at a time; a fit with
sampling weights, which a resample of rows does not carry; a fit to a
covariance matrix, which has no cases to resample; and a fit that did
not converge.

Bootstrap intervals from supplied two-stage fits currently refuse
separately stored auxiliary variables: each refit must resample and
carry those values with the items. Use an eligible delta interval where
available, or explicitly resample and refit them together.

## References

Aguirre-Urreta, M. I., Rönkkö, M., & McIntosh, C. N. (2019). A
cautionary note on the finite sample behavior of maximal reliability.
*Psychological Methods, 24*(2), 236–252.
[doi:10.1037/met0000176](https://doi.org/10.1037/met0000176)

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science, 11*(3), 189–228.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*(1), 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Raykov, T., & Penev, S. (2006). A direct method for obtaining
approximate standard error and confidence interval of maximal
reliability for composites with congeneric measures. *Multivariate
Behavioral Research, 41*(1), 15–28.
[doi:10.1207/s15327906mbr4101_2](https://doi.org/10.1207/s15327906mbr4101_2)

Reise, S. P. (2012). The rediscovery of bifactor measurement models.
*Multivariate Behavioral Research, 47*(5), 667–696.
[doi:10.1080/00273171.2012.715555](https://doi.org/10.1080/00273171.2012.715555)

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
bifactor models: Calculating and interpreting statistical indices.
*Psychological Methods, 21*(2), 137–150.
[doi:10.1037/met0000045](https://doi.org/10.1037/met0000045)

Yuan, K.-H., & Hayashi, K. (2003). Bootstrap approach to inference and
power analysis based on three test statistics for covariance structure
models. *British Journal of Mathematical and Statistical Psychology,
56*, 93–110.

## See also

[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
for the point estimates;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for a single-factor omega and its intervals;
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
for coefficient H.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
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
# The nine-item bifactor model of ?bifactor_indices: one general factor
# and three orthogonal group factors, simulated so that resamples give
# proper solutions and the example runs quickly. (The 15-test
# Holzinger and Swineford model of Details gives an improper solution
# in a fifth of its resamples, each of which took seconds to fit.)
set.seed(113)
n <- 600
g <- rnorm(n); grp <- list(rnorm(n), rnorm(n), rnorm(n))
X <- vapply(1:9, function(i)
  0.5 * g + 0.5 * grp[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n), numeric(n))
colnames(X) <- paste0("x", 1:9)
model <- "g  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
          f1 =~ x1 + x2 + x3
          f2 =~ x4 + x5 + x6
          f3 =~ x7 + x8 + x9"
fit <- lavaan::cfa(model, data = as.data.frame(X),
                   orthogonal = TRUE, std.lv = TRUE)

# Every replication refits the bifactor model. B = 10 keeps the
# example quick and is below the 400 a reported interval needs; the
# option lets it run and the warning says so. A reported interval
# deserves the default B = 1000. A resample whose solution is improper
# is set aside, and the warning and the printed lines count it.
old <- options(dmar.bootstrap_exploratory = TRUE)
ci <- ci_bifactor_indices(fit, B = 10, seed = 113)
#> Warning: ci_bifactor_indices(): B = 10 is below the 400 replications a 95% interval is reported from (10 beyond each limit). The limits are shown because options(dmar.bootstrap_exploratory = TRUE) is set. With 10 replications each limit rests on about 0.3 of them and will move from seed to seed; use the result to explore or to demonstrate, and report an interval only from B of at least 400. B, B_used, and B_required are recorded in the result.
options(old)
ci
#>  factor index               estimate se    ci_lower ci_upper
#>  g      ECV                 0.547    0.044 0.484    0.637   
#>  g      omega               0.858    0.007 0.851    0.869   
#>  g      omega_H             0.673    0.026 0.636    0.731   
#>  g      PUC                 0.750    <NA>  <NA>     <NA>    
#>  g      H                   0.769    0.019 0.747    0.814   
#>  g      maximal_reliability 0.678    0.030 0.647    0.752   
#>  f1     ECV                 0.609    0.085 0.487    0.813   
#>  f1     omega               0.727    0.017 0.707    0.760   
#>  f1     omega_HS            0.280    0.058 0.135    0.362   
#>  f1     H                   0.412    0.083 0.226    0.516   
#>  f1     maximal_reliability 0.396    0.089 0.220    0.537   
#>  f2     ECV                 0.543    0.077 0.456    0.708   
#>  f2     omega               0.756    0.017 0.736    0.787   
#>  f2     omega_HS            0.345    0.057 0.226    0.424   
#>  f2     H                   0.477    0.055 0.362    0.557   
#>  f2     maximal_reliability 0.456    0.057 0.335    0.545   
#>  f3     ECV                 0.490    0.079 0.425    0.646   
#>  f3     omega               0.743    0.017 0.707    0.761   
#>  f3     omega_HS            0.377    0.059 0.257    0.428   
#>  f3     H                   0.506    0.058 0.377    0.558   
#>  f3     maximal_reliability 0.494    0.064 0.358    0.556   
#> 
#> Confidence level: 95%
#> 
#> Denominator: composite variance implied by the fitted model
#> 
#> Intervals: percentile bootstrap over the 600 rows the fit used, with improper replications set aside; PUC depends only on which items share a group factor, so it has no interval
#> 
#> Bootstrap estimates of bias (mean of the replicates used minus the estimate): H g 0.00261, f1 0.0272, f2 0.014, f3 -0.0198; maximal_reliability g 0.00693, f1 0.0441, f2 0.0185, f3 -0.0185
#> 
#> For a group factor, ECV is the general factor's share of the common variance of that factor's items (one minus it is the group factor's share), and omega is the reliability of their total score
#> H is the index of Rodriguez, Reise, and Haviland (2016), which counts the other factors' common variance as uncorrelated error; maximal_reliability is the reliability the optimally weighted items reach as a measure of the factor under the fitted model
#> 
#> Bootstrap: 10 replications attempted, 10 returned a value; at least 400 are required for a reported interval. Exploratory: below the number required for a reported interval

# Every replication's indices, and what happened to it.
head(attr(ci, "bootstrap")$t[, c("g:omega", "g:omega_H", "g:H",
                                 "g:maximal_reliability")])
#>        g:omega g:omega_H       g:H g:maximal_reliability
#> [1,] 0.8521588 0.6695526 0.7664138             0.6787335
#> [2,] 0.8507987 0.6567387 0.7573378             0.6661037
#> [3,] 0.8586255 0.6597966 0.7656124             0.6704680
#> [4,] 0.8660473 0.6778651 0.7758604             0.6839321
#> [5,] 0.8519372 0.6848538 0.7755194             0.6945883
#> [6,] 0.8529592 0.6590697 0.7579440             0.6655351
table(attr(ci, "bootstrap")$status)
#> 
#>        proper      improper         error    non_finite  undetermined 
#>            10             0             0             0             0 
#> not_converged 
#>             0 

# The broom verbs: one row per index, and the counts in one row.
generics::tidy(ci)
#>                      term  estimate          se  ci_lower  ci_upper conf_level
#> 1                   g:ECV 0.5465883 0.044195296 0.4839474 0.6371550       0.95
#> 2                 g:omega 0.8579575 0.006517243 0.8507987 0.8694641       0.95
#> 3               g:omega_H 0.6734780 0.026377176 0.6357102 0.7307881       0.95
#> 4                   g:PUC 0.7500000          NA        NA        NA       0.95
#> 5                     g:H 0.7688763 0.019215757 0.7470697 0.8136444       0.95
#> 6   g:maximal_reliability 0.6775262 0.030165817 0.6471463 0.7521177       0.95
#> 7                  f1:ECV 0.6092285 0.084855419 0.4872606 0.8126619       0.95
#> 8                f1:omega 0.7273485 0.016972482 0.7066682 0.7595878       0.95
#> 9             f1:omega_HS 0.2797758 0.058145520 0.1345353 0.3617250       0.95
#> 10                   f1:H 0.4117808 0.082550447 0.2262811 0.5164175       0.95
#> 11 f1:maximal_reliability 0.3960647 0.089122416 0.2202457 0.5370388       0.95
#> 12                 f2:ECV 0.5427861 0.077494879 0.4560744 0.7081832       0.95
#> 13               f2:omega 0.7558981 0.017082830 0.7357629 0.7870244       0.95
#> 14            f2:omega_HS 0.3454961 0.057309359 0.2255118 0.4242169       0.95
#> 15                   f2:H 0.4771691 0.055102244 0.3624992 0.5570588       0.95
#> 16 f2:maximal_reliability 0.4562862 0.057094770 0.3346749 0.5453874       0.95
#> 17                 f3:ECV 0.4902972 0.078748894 0.4253930 0.6464365       0.95
#> 18               f3:omega 0.7432433 0.016653932 0.7068112 0.7607319       0.95
#> 19            f3:omega_HS 0.3770962 0.058573288 0.2572496 0.4276934       0.95
#> 20                   f3:H 0.5055844 0.058065932 0.3767877 0.5582657       0.95
#> 21 f3:maximal_reliability 0.4939142 0.063884581 0.3578673 0.5557906       0.95
generics::glance(ci)
#>   n_factors   N  B B_used set_aside improper_set_aside improper_kept
#> 1         4 600 10     10         0                  0             0
#>   not_converged caution conf_level  ci_method improper_rule   denominator
#> 1             0   FALSE       0.95 percentile          drop model_implied
```
