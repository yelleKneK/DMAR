# How DMAR Counts and Requires Bootstrap Replications

One rule governs every DMAR function that reads an interval, a bound, a
reference quantile, a standard error, or a Monte Carlo *p*-value from
`B` bootstrap or Monte Carlo replications: how many replications the
request must ask for, what happens to a replication that returns no
value, what the result records, and when no interval is reported. This
page states the rule and its reasons; each function's page says how the
rule applies to it.

## How many replications an interval needs

A percentile limit is a quantile of the `B` replications, so it is
located by the replications that lie beyond it. With `B` = 100 and a 95
percent interval there are two or three of them on each side, and the
limit moves with every new seed. DMAR runs an interval read from
quantiles only when the number of replications requested puts at least
10 beyond each limit, and never from a request of fewer than 100: 200,
400, and 2,000 replications at the 90, 95, and 99 percent levels. The
count is held to the number requested, `B`; what happens when some of
them return no value is decided afterward and is described under
*Replications that return no value*. In symbols, with \\t\\ the smallest
tail area a reported limit cuts off (\\(1 - C)/2\\ for a two-sided
interval at level \\C\\), \$\$B\_{\mathrm{required}} = \max\\100, \lceil
10 / t \rceil\\.\$\$ With 10, the share of the bootstrap distribution a
limit actually cuts off has a relative standard deviation of about 30
percent whatever the shape of that distribution (a nominal 2.5 percent
is between about 1.2 and 4.5 percent in 95 of 100 runs, depending
slightly on how the quantile is interpolated). How far the limit itself
moves depends on the shape: with 400 replications at the 95 percent
level its standard deviation from run to run is about 0.13 bootstrap
standard errors for a bell-shaped distribution, about 0.20 for one with
tails as heavy as a *t* distribution on 5 degrees of freedom, and about
0.29 on the long side of a skewed one. Ten is a floor and not a
recommendation. The defaults (1,000 where each replication refits a
model, 2,000 to 10,000 elsewhere) put 25 to 250 replications beyond each
limit at the 95 percent level, and at 1,000 the standard deviation of a
limit of a bell-shaped distribution is about 0.08 standard errors. A BCa
interval needs more, 1,000, 2,000, and 10,000 at the three levels, for
the reason given under *BCa intervals*. An interval built from the
standard deviation of the replications (`"bootstrap_se"`) is reported
from 100, at which that standard deviation has a coefficient of
variation of about 7 percent for a bell-shaped distribution (12 to 14
percent for the heavy-tailed and the skewed shapes above), and a Monte
Carlo *p*-value from 100, at which its standard error is at most 0.05. A
request below the number required is refused before anything is
computed, with a condition of class `dmar_count_invalid` that states the
number. None of these numbers concerns the coverage of the bootstrap for
a given statistic, which no number of replications repairs; they limit
only the error added by drawing a finite number of resamples. The
numbers in this paragraph were computed for DMAR by simulation (20,000
runs each); the statement about the share cut off follows from the Beta
distribution of an order statistic (the tail area beyond the \\k\\-th of
\\B\\ ordered draws from a continuous distribution is distributed as
Beta(\\k\\, \\B - k + 1\\)). No published source is cited for the
counts.

A limit at quantile position \\p\\ has the tail area \\\min(p, 1 - p)\\.
For an interval with unequal tail areas (`alpha_lower` and `alpha_upper`
of the ratio functions and of
[`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md))
the smaller area decides; for a one-sided bound
([`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md)) and for a
reference quantile
([`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md),
[`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md))
it is the area beyond that one quantile, so a bound at the 95th
percentile needs 200 replications and one at the median needs the floor
of 100.

The functions that compute their limits directly read them with
[`quantile`](https://rdrr.io/r/stats/quantile.html) (its default, type
7); the reliability family and
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
read them with
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html), which
places a limit at the \\(B + 1)\alpha\\ position among the ordered
replications. The two differ by a fraction of the distance between
adjacent replications.

## Replications that return no value

A resample can leave a statistic undefined (a constant variable, a model
that does not converge, a denominator of zero). Such a replication is
dropped and counted, never replaced. Every result records `B`
(requested), `B_used` (the number the interval rests on), `n_dropped`,
and `B_required` (the fewest replications a request may ask for) as
attributes, which
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)
returns, with `bootstrap_status`, and one warning of class
`dmar_bootstrap_replications_dropped` states the counts. An interval
from the replications that remain describes the resamples in which the
statistic is defined. When the dropped replications number at least the
share of the distribution beyond one limit (2.5 percent of `B` for a 95
percent interval), they could by themselves decide where that limit
falls, and the warning says so (`caution` is `TRUE` in the condition).
When more than half were dropped no interval is reported: the limits are
`NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`. The call does not stop for this, so an
estimate that took minutes to compute is kept. With more than half
dropped a larger `B` does not help: the replications that remain are a
minority of the resamples, selected by whether the statistic could be
computed, and their spread describes those resamples and not the
sampling distribution of the estimate.

The number required is counted on the request and not on the
replications that remain. A request of exactly 400 at the 95 percent
level from which 3 replications are dropped reports its interval from
the 397 that remain, with the warning. The two rules about shares, in
the paragraph above, then decide: more than half dropped withholds the
interval, and dropped replications that reach the share beyond one limit
bring the caution. The warning states how many were dropped, and when
those that remain are fewer than the number required it says so and says
to request more. Where replications can fail (a model refit in every
resample, sparse categories, a small sample), ask for more than the
minimum, so that the limits rest on at least the number the rule
intends; the defaults do.

`bootstrap_status` is one of `"complete"` (nothing dropped), `"dropped"`
(some dropped, the interval formed), `"exploratory"` (formed from a
request below the required number, under the option described below),
`"insufficient_count"` (no interval: a BCa interval whose adjusted
positions take more replications than were requested, or fewer than two
usable replications), and `"insufficient_share"` (more than half
dropped, no interval). With several statistics in one call (the pairs of
[`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md), the
effects of
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md),
the coefficients of
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)), `B_used`
is the smallest count any reported limit rests on and `n_dropped` is
`B - B_used`. Where a table already has rows for the counts it keeps
them ([`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md)
has rows `B` and `B_used`,
[`krippendorff_alpha`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md)
the row `B_used`); the reliability family keeps its fuller `bootstrap`
attribute, which gains `required` and `bootstrap_status`.

## BCa intervals

The BCa interval reads its limits from quantiles that the bias
correction and the acceleration move away from the nominal ones. For a
95 percent interval the positions might be 0.9 and 99.1 percent in place
of 2.5 and 97.5, and the limit at 0.9 percent has few replications
beyond it. The count of 10 is applied where each limit finally sits,
again to the number of replications requested, so a BCa interval can be
withheld after the replications are run (warning of class
`dmar_bootstrap_insufficient`, with the positions and the count that
would have served). To make that uncommon, a BCa interval is run only
with 50 replications beyond each nominal limit: 1,000, 2,000, and 10,000
at the 90, 95, and 99 percent levels. In simulations of an indirect
effect and of a correlation at the 95 percent level, the adjusted
positions asked for more than 2,000 replications in at most 14 percent
of data sets (an indirect effect with 50 cases and small paths) and in 0
to 3 percent of the others, and for more than 10,000 in about 1 percent
of the data sets of that hardest condition and in none of the others
(400 data sets per condition).

When the bias correction or the acceleration is itself undefined (a
leave-one-out estimate that cannot be computed, every replicate on one
side of the estimate, influence values without spread), the functions
differ by family.
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md),
[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
and
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
refuse the BCa interval with an error of class `dmar_bca_undefined`
whose message names the cause;
[`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md)
and
[`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md)
stop with an error that names `ci_method = "percentile"`; the
reliability family and
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
report the percentile interval in its place, with a warning, and record
the substitution in the `ci_method` attribute.

## Fewer replications, for exploration

A first look at an analysis, a help page example, and a simulation that
averages over many data sets do not need a reportable interval from each
run. `options(dmar.bootstrap_exploratory = TRUE)` lets a function run
with any `B` of 2 or more. Nothing else changes: the result still
records `B_required`, its footer says that the interval is below the
requirement, and every call gives a warning of class
`dmar_bootstrap_exploratory`. A planner that runs `B` replications
inside each of many simulated studies
([`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md))
gives the warning once for the call and records the request (`B_used` is
then `B`, the number each simulated study ran), and
[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md),
whose result is a plot, gives the warning and records nothing. The share
rules still apply. An interval obtained this way is not one to report.
Set the option for the calls that need it and restore it afterward:

    old <- options(dmar.bootstrap_exploratory = TRUE)
    mediate(d, x = "x", m = "m", y = "y", B = 50, seed = 113)
    options(old)

## Monte Carlo intervals and reference quantiles

The same rule holds for limits read from the quantiles of Monte Carlo
draws (`ci_method = "monte_carlo"` of
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md) and
of
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md),
the band of
[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md))
and for the reference quantile of
[`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
and
[`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md),
with "draw" or "reference data set" in place of "replication". The
Bollen-Stine *p*-value of
[`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)
is a proportion and not a quantile: it needs a request of 100 samples,
is withheld when more than half were set aside, and is reported with the
Wilson interval of its own Monte Carlo error.

## What is not a bootstrap

The rule does not govern `n_resamples` of
[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
and
[`procrustes_phi`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
or the number of simulated studies `G` of the sample size planning and
sensitivity functions. Those counts are checked as counts and their
functions report their own Monte Carlo error.

## Reproducibility

Every function that draws replications has a `seed` argument
(`boot_seed` in
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) and
[`mlmr_mv`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md)) that
defaults to `NULL`: the replications are then drawn from, and advance,
the current random number stream. A supplied seed is set for the call
only and the caller's generator state is restored when the function
returns, so a reproducible call does not disturb the random numbers of
the script around it. Two runs with the same seed and the same `B` give
the same limits; the same seed with a different `B` does not.

## Conditions

- `dmar_count_invalid`:

  Error (a `dmar_refusal`). `B` is not a single whole number of at least
  the number required. Fields `arg`, `min` (the number required), `max`,
  and `reason` (the sentence that explains the number).

- `dmar_bootstrap_replications_dropped`:

  Warning. Replications were dropped and the interval was formed from
  the rest. The message states how many, and says to request more when
  those that remain are fewer than `B_required`. Fields `B`, `B_used`,
  `n_dropped`, `B_required`, `share_dropped`, `caution`, and `status`.

- `dmar_bootstrap_insufficient`:

  Warning. No interval is reported: more than half of the replications
  were dropped (`status` `"insufficient_share"`), or the adjusted
  positions of a BCa interval take more replications than were
  requested, or fewer than two replications were usable
  (`"insufficient_count"`). The same fields; for a BCa interval withheld
  at its adjusted positions also `positions` (the two positions, as
  proportions) and `B_suggest` (one and a half times the count those
  positions take, rounded up to a hundred).

- `dmar_bootstrap_exploratory`:

  Warning. The result was computed from fewer replications than a
  reported result needs, because
  `options(dmar.bootstrap_exploratory = TRUE)` is set. The same fields
  as the dropped warning.

- `dmar_bootstrap_note`:

  Warning. Something else about the run that the reader of the interval
  should know, with every replication used and the interval formed: in
  the reliability family and
  [`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
  a BCa interval replaced by the percentile interval (the `ci_method`
  attribute records the substitution), or replicates that all have the
  same value, so that the interval has no width; in
  [`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
  and
  [`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md),
  the function's account of a run in which nothing was dropped
  (replications from an improper solution that were kept, reference data
  sets whose fits gave warnings). The same fields as the dropped warning
  where the function records them.

- `dmar_bca_undefined`:

  Error (a `dmar_refusal`). The bias correction or the acceleration of a
  BCa interval is undefined
  ([`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md),
  [`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
  [`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)).

The warnings of the reliability family and of
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
which share one bootstrap engine, begin with the name of the function
that was called (`"reliability_omega(): "`) and carry that name in the
field `fun`;
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
passes its call on, so its warning names the function that ran the
bootstrap. A program can muffle or read any of the warnings by class:

    withCallingHandlers(
      mediate(d, x = "x", m = "m", y = "y", seed = 113),
      dmar_bootstrap_replications_dropped = function(w) {
        message(w$n_dropped, " of ", w$B, " dropped")
        invokeRestart("muffleWarning")
      })

## The functions the rule governs

"Required" is the fewest replications a request may ask for at the 95
percent level (for the reference quantiles, at the 95th percentile).
Kind Q reads quantiles, S takes a standard deviation of the
replications, and P is a Monte Carlo *p*-value.

|  |  |  |  |  |
|----|----|----|----|----|
| **Function** | **Kind** | **Default `B`** | **Required** | **Why that default** |
| [`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md), percentile and Monte Carlo | Q | 2000 | 400 | each replication refits two regressions |
| [`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md), BCa | Q, BCa | 2000 | 2000 | the same |
| [`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md), [`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md) (Monte Carlo) | Q | 10000 | 400 | a draw evaluates the effects, no refit |
| [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md), [`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md), [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md), [`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md), [`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md): percentile | Q | 10000 | 400 | a closed form or a small model per replication |
| the same functions: BCa | Q, BCa | 10000 | 2000 | the same |
| the same functions: `"bootstrap_se"`, `"bootstrap_se_logistic"` | S | 10000 | 100 | the same |
| [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md) | Q, BCa, or S | 10000 (1000 under `missing = "mi"`) | 400, 2000, or 100 | under multiple imputation each replication is imputed several times |
| [`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md) | Q or BCa | 1000 | 400 or 2000 | each replication refits the bifactor model |
| [`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md) | Q | 1000 | 400 | each replication refits the model |
| [`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md) | Q | 1000 | 400 | each replication refits the model three times |
| [`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md), [`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md) | Q or BCa | 10000 | 400 or 2000 | kappa is cheap to recompute |
| [`krippendorff_alpha`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md) | Q | 1000 | 400 | each resample recomputes all coincidences |
| [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) | Q, one-sided | 0 (no bootstrap) | 200 | the bound is opt-in |
| [`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md), [`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md) | Q or BCa | 10000 | 400 or 2000 | a ratio of means is cheap to recompute |
| [`R2_mixed_effects`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md) | Q | 1000 | 400 | each replication refits the mixed model |
| [`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md) | Q | 10000 | 1000 | its own floor of 1000 |
| [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md), [`mlmr_mv`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md) | Q | 1000 | 400 | each resample refits the model |
| [`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md), [`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md) | Q, reference quantile | 1000 | 200 | each reference data set is generated and fit |
| [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md) | P | 1000 | 100 | each sample refits the model |

The planning functions
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
and
[`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md)
run an interval inside each simulated study; their pages state the inner
`B` they use.

## References

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC. (For the percentile and
BCa intervals themselves, as the function pages cite it. The counts on
this page are not taken from it.)

## See also

[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)
for reading the counts from a result, and the pages of the functions in
the table.

## Author

Ken Kelley <kkelley@nd.edu>
