# Find Target Sample Sizes for the Accuracy in Standardized Conditions Means Estimation in CRD

Finds the number of clusters, the number of individuals per cluster, or
both, for a two-condition cluster randomized design (CRD), so that the
confidence interval for the standardized difference between the
condition means has a desired width or is as narrow as a budget allows.
Four of the functions fix one of the two sizes and take a desired width
or a budget; the two `both` functions take a budget and return the sizes
that give the narrowest interval, or take a desired width and return the
least costly sizes. Each returns a table with the size or sizes found,
the cost, and the expected width of the interval where it was computed.

The interval is the likelihood-based confidence interval (Cheung, 2009)
computed by the `OpenMx` package (Boker et al., 2011). See
Pornprasertmanit and Schneider (2014) for the method.

## Usage

``` r
ss_aipe_crd_es_n_clusters_fixed_width(
  width,
  cluster_size,
  true_es,
  es_type = 1,
  true_icc,
  pr_treat,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  G = 1000,
  true_icc_covariate = NULL,
  seed = NULL,
  ncpus = 1,
  cluster_cost = NULL,
  individual_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_es_n_individuals_fixed_width(
  width,
  n_clusters,
  true_es,
  es_type = 1,
  true_icc,
  pr_treat,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  G = 1000,
  true_icc_covariate = NULL,
  seed = NULL,
  ncpus = 1,
  cluster_cost = NULL,
  individual_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_es_n_clusters_fixed_budget(
  budget,
  cluster_size,
  cluster_cost,
  individual_cost,
  G = NULL,
  pr_treat = NULL,
  true_icc = NULL,
  true_es = NULL,
  es_type = 1,
  n_covariates = 0,
  true_icc_covariate = NULL,
  true_R2_within = 0,
  true_R2_between = 0,
  assurance = NULL,
  seed = NULL,
  ncpus = 1,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_es_n_individuals_fixed_budget(
  budget,
  n_clusters,
  cluster_cost,
  individual_cost,
  G = NULL,
  pr_treat = NULL,
  true_icc = NULL,
  true_es = NULL,
  es_type = 1,
  n_covariates = 0,
  true_icc_covariate = NULL,
  true_R2_within = 0,
  true_R2_between = 0,
  assurance = NULL,
  seed = NULL,
  ncpus = 1,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_es_both_fixed_budget(
  budget,
  cluster_cost = 0,
  individual_cost = 1,
  true_es,
  es_type = 1,
  true_icc,
  pr_treat,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  G = 1000,
  true_icc_covariate = NULL,
  seed = NULL,
  ncpus = 1,
  diff_size = NULL
)

ss_aipe_crd_es_both_fixed_width(
  width,
  cluster_cost = 0,
  individual_cost = 1,
  true_es,
  es_type = 1,
  true_icc,
  pr_treat,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  G = 1000,
  true_icc_covariate = NULL,
  seed = NULL,
  ncpus = 1,
  diff_size = NULL
)
```

## Arguments

- width:

  The desired width of the confidence interval of the standardized
  difference between the condition means

- cluster_size:

  The number of individuals in each cluster: a single whole number of at
  least 2

- true_es:

  The population effect size the interval is planned for (on the scale
  chosen by `es_type`): a single finite number, such as `0.5`. It is
  required wherever a width is computed, because the width is read from
  data simulated at that value: in the four planners that take or return
  a width, and in the two planners that fix a size and a budget when `G`
  is supplied. A call that asks for a width without it is refused with
  an error of class `dmar_argument_invalid` that names what is missing

- es_type:

  The type of effect size. There are only three possible options: 0 =
  the effect size using total standard deviation, 1 = the effect size
  using the individual-level standard deviation (level 1), 2 = the
  effect size using the cluster-level standard deviation (level 2)

- true_icc:

  The population intraclass correlation of the dependent variable (a
  planning value): a single number in \[0, 1), and above 0 with
  `es_type = 2`. At 1 no variance lies within clusters, and the model
  fitted to each simulated data set has no within-cluster variance to
  estimate; with `es_type = 2` the effect is in units of the
  between-cluster standard deviation, which an intraclass correlation of
  0 makes 0

- pr_treat:

  The proportion of treatment clusters: a single number strictly between
  0 and 1

- true_R2_between:

  The population proportion of the between-level variance that the
  covariate explains (a planning value in \[0, 1), used when
  `n_covariates = 1`; the default is 0 in all six functions). It must be
  0 when `n_covariates = 0`

- true_R2_within:

  The population proportion of the within-level variance that the
  covariate explains (a planning value in \[0, 1), used when
  `n_covariates = 1`; the default is 0 in all six functions)

- n_covariates:

  The number of covariates at the cluster (between) level, not counting
  the treatment indicator: 0 (the default, the treatment indicator
  alone) or 1. A larger value is refused. It is not `p`, the total
  number of predictors, which here is `p = n_covariates + 1`.

- assurance:

  The assurance: the probability, over repeated studies of the design,
  that the study gives a confidence interval and that the interval is no
  wider than the width planned (e.g., .80, .90, .95). With `NULL` (the
  default) the plan is for the expected width instead. A simulated data
  set whose interval cannot be formed counts as one whose interval is
  wider than any width (see the section on simulated data sets without a
  usable interval)

- conf_level:

  The desired level of confidence for the confidence interval

- G:

  The number of replications used in the a priori Monte Carlo
  simulation: a single whole number from 1 to 999,999 (default 1000 in
  the four planners that need a simulation). Each simulated data set
  takes its own seed, drawn without replacement from the whole numbers 1
  to 999,999, so no more can be simulated. Anything else (a fractional
  number, 0, a negative number, a number above 999,999, `NA`, a
  character string, or more than one value) is refused with an error of
  class `dmar_count_invalid` before any simulation starts. In the two
  planners that fix a size and a budget the default is `NULL`, which
  skips the simulation: the size the budget buys is returned without an
  expected width. Supplying `G` there asks for the expected width, which
  also needs `pr_treat`, `true_icc`, and `true_es`; a call with `G` and
  without any of them is refused, with a message that names the ones
  missing, in place of a table without the width that was asked for.

- true_icc_covariate:

  The population intraclass correlation of the covariate (a planning
  value, used when `n_covariates = 1`; when it is `NULL` the value of
  `true_icc` is used). If `true_icc_covariate = 0`, the covariate varies
  within clusters only. If `true_icc_covariate = 1`, the covariate
  varies between clusters only

- seed:

  An optional integer seed for the a priori Monte Carlo simulation. The
  default `NULL` draws the seeds of the simulated data sets from the
  current random number stream, which advances that stream, so repeated
  calls reflect the sampling variability of the simulation. An integer
  makes the result reproducible: the seed is set for the call only, and
  the caller's random number state is restored when the function returns

- ncpus:

  The number of processes over which the `G` simulated data sets are
  spread, the name this count has in the bootstrap functions of the
  package
  ([`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  and its relatives). The default `1` runs them in this R session; a
  larger whole number uses forked workers or, on Windows, which cannot
  fork, a socket cluster. Each simulated data set has its own seed,
  drawn before any is run, so a given `seed` returns the same result
  whatever `ncpus` is. A socket cluster is stopped when the function
  exits, whether it returns, stops with an error, or is interrupted.

- cluster_cost:

  The cost of adding one cluster, whatever the number of individuals in
  it: a single finite number of 0 or more. Where the costs choose the
  design (the two budget planners that fix a size and the two `both`
  planners), `cluster_cost` and `individual_cost` cannot both be 0

- individual_cost:

  The cost of adding one individual to a cluster: a single finite number
  of 0 or more, and positive where the costs choose the cluster size
  (`ss_aipe_crd_es_n_individuals_fixed_budget` and the two `both`
  planners)

- diff_size:

  Difference cluster size specification. The differences in cluster
  sizes can be specified in two ways, and the specified vector is
  recycled across the clusters. First, users may specify differences as
  integers, which can be negative or positive; the resulting cluster
  sizes add the specified values to the estimated cluster size. For
  example, if the cluster size is 25, the number of clusters is 10, and
  `diff_size = c(-1, 0, 1)`, the cluster sizes will be 24, 25, 26, 24,
  25, 26, 24, 25, 26, and 24. Second, users may specify multipliers of
  the cluster size as positive decimals; at least one value must be
  non-integer, which is what selects the multiplicative form. The
  resulting cluster sizes multiply the estimated cluster size by the
  specified values and round to the nearest integer. For example, if the
  cluster size is 25, the number of clusters is 10, and
  `diff_size = c(0.8, 1, 1.2)`, the cluster sizes will be 20, 25, 30,
  20, 25, 30, 20, 25, 30, and 20. In either form a resulting cluster
  size below 1 is set to 1. If `NULL`, the cluster size is equal across
  clusters

- n_clusters:

  The number of clusters: a single whole number of at least
  `3 + n_covariates`, which with `pr_treat` puts a cluster in each
  condition

- budget:

  The amount available for the study, in the units of `cluster_cost` and
  `individual_cost`: a single finite number of 0 or more

## Value

The `ss_aipe_crd_es_n_clusters_fixed_width` and
`ss_aipe_crd_es_n_clusters_fixed_budget` functions provide the number of
clusters. The `ss_aipe_crd_es_n_individuals_fixed_width` and
`ss_aipe_crd_es_n_individuals_fixed_budget` functions provide the
cluster size. The `ss_aipe_crd_es_both_fixed_budget` and
`ss_aipe_crd_es_both_fixed_width` provide the number of clusters and the
cluster size, respectively.

Each returns a table with the columns `term` and `value`. The rows,
where they apply, are `necessary_n_clusters` (the number of clusters,
where it is solved for), `necessary_cluster_size` (the number of
individuals in each cluster, where it is solved for; with `diff_size`
the cluster sizes and their frequencies are printed ahead of the table
instead), `width` (the desired width, in the functions that plan from
one), the width of the interval for the standardized effect size chosen
by `es_type`, and `budget` (the cost of the design, when the two costs
are supplied). Without `assurance` the width row is `expected_width`,
the mean width of the interval over the `G` simulated data sets. With
`assurance` it is `assured_width`, the width that the interval does not
exceed with the requested probability, followed by the row `assurance`
that echoes that probability.

A size that was given is an argument of the call and is not repeated in
the table: only a size the function solved for has a row, and its name
begins with `necessary_`.
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) report the
planned size as `estimate` (the number of clusters where it was solved
for, otherwise the cluster size) and the desired width as `width`, which
is `NA` for a plan made from a budget; `glance()` carries
`expected_width` or `assured_width` as a column of that name.

A table whose width was simulated carries the record of the simulation
of the design reported in its attribute `simulation`, a list read with
`attr(x, "simulation")`: `attempted` (`G`), `returned` (the data sets
with a usable interval), `failed` (those without one), `denominator`
(the number of data sets the width row is taken over: `returned` for
`expected_width`, `attempted` for `assured_width`), `failure_reasons`
(the data sets without a usable interval by reason, a named integer
vector: `"fit_error"`, `"optimizer_status_<code>"`, `"limit_not_found"`,
`"limit_status_<code>"`, `"limit_not_finite"`, `"limits_not_ordered"`),
`failure_messages` (for a fit that stopped with an error, the first
error message, named by reason), `width` (the width reported),
`width_mc_ci` (the Monte Carlo interval of that width at `conf_level`,
which describes the error of the simulation only: the mean plus and
minus the normal quantile times its standard error, or for the assured
width the interval of two order statistics of the `G` widths),
`n_clusters` and `cluster_size` (the design), and `search`, the totals
over every design the call simulated (`n_designs`, `attempted`,
`returned`, `failed`, `failure_reasons`, `failure_messages`) and the
warnings the fits gave (`warnings`, one element per kind, with its
`class`, first `message`, all its `messages`, and `n`, the number of
data sets that gave it). The line "Replications" beneath the printed
table gives `G`, and the three counts when a data set gave no usable
interval.

## Details

Here are the functions' descriptions:

- `ss_aipe_crd_es_n_clusters_fixed_width`:

  Find the number of clusters given a specified width of the confidence
  interval and the cluster size

- `ss_aipe_crd_es_n_individuals_fixed_width`:

  Find the cluster size given a specified width of the confidence
  interval and the number of clusters

- `ss_aipe_crd_es_n_clusters_fixed_budget`:

  Find the number of clusters given a budget and the cluster size

- `ss_aipe_crd_es_n_individuals_fixed_budget`:

  Find the cluster size given a budget and the number of clusters

- `ss_aipe_crd_es_both_fixed_budget`:

  Find the sample size combinations (the number of clusters and that
  cluster size) providing the narrowest confidence interval given the
  fixed budget

- `ss_aipe_crd_es_both_fixed_width`:

  Find the sample size combinations (the number of clusters and that
  cluster size) providing the lowest cost given the specified width of
  the confidence interval

Every answer that targets a width, and the expected width a budget
planner reports when `G` and the population values are supplied, rests
on an a priori Monte Carlo simulation: a candidate design is evaluated
by generating `G` data sets and reading the likelihood-based confidence
interval on the standardized effect size from OpenMx. The planners whose
answer needs one such evaluation, or a handful, run in about a second at
a small `G` and are the ones shown in the examples. The two that search
over the number of clusters, `ss_aipe_crd_es_n_clusters_fixed_width` at
a fixed cluster size and `ss_aipe_crd_es_both_fixed_width`, which
repeats that search across candidate cluster sizes and keeps the least
costly combination that reaches the width, evaluate many candidates in
turn and run for minutes at the default `G = 1000`, so they are not
among the examples. A call such as
`ss_aipe_crd_es_n_clusters_fixed_width(width = 0.3, cluster_size = 20, true_es = 0.5, es_type = 1, true_icc = 0.25, pr_treat = 0.5, G = 1000, seed = 113)`
returns the number of clusters of 20 individuals that brings the
expected width of the interval on the individual-level standardized
effect size to 0.3, and
`ss_aipe_crd_es_both_fixed_width(width = 0.5, cluster_cost = 5, individual_cost = 1, true_es = 0.5, es_type = 1, true_icc = 0.25, pr_treat = 0.5, G = 1000, seed = 113)`
returns the least costly pairing of clusters and cluster size at those
costs that reaches a width of 0.5. Both accept `diff_size` for unequal
cluster sizes, in which case the simulated data sets carry the unequal
sizes and the confidence interval comes from a multiple-group model.

The width of a simulated interval is its upper limit minus its lower
limit.

## Simulated data sets without a usable interval

Each simulated data set is fitted with OpenMx, and its interval is used
only when the fit and the search for each limit succeeded: the optimizer
status code of the fit and of each limit's search is 0 or 1 (the codes
that OpenMx's `as.statusCode()` describes as `"OK"` and `"OK/green"`,
optimization succeeded; with a code such as 5, a Hessian that does not
appear to be convex, or 6, a solution that does not meet the first-order
optimality conditions, the fit may not be at the maximum the interval is
measured from, and that help page advises against trusting it), OpenMx's
diagnostic for each limit is `"success"`, and the two limits are finite
and bracket the estimate. A data set that fails any of these, or whose
fit stops with an error, gives no usable interval. How often that
happens depends on the design. With OpenMx 2.22.11 and a standardized
effect of 0.5 (the intraclass correlation in parentheses), the share of
1000 simulated data sets without a usable interval, with its Monte Carlo
standard error, was 0 with 20 clusters of 10 (0.2) and `es_type = 1`,
0.002 (0.0014) with six clusters of 5 (0.1), 0 with ten clusters of 2
(0.01), 0 with eight clusters of 4 (0.1) and `es_type = 0`, 0.415
(0.016) with six clusters of two to four individuals (0.05), and 0.208
(0.013) with 20 clusters of 10 (0.2) and `es_type = 2`; a second run of
400 data sets gave 0.378 (0.024) and 0.237 (0.021) for the last two.
These are single runs, and other versions of OpenMx can give other
shares. With the unequal clusters 412 of the 415 failed fits ended with
status 5 or 6, and in 200 further data sets each of the 86 such fits had
estimated the between-cluster variance below 0. The effect of
`es_type = 2` is divided by the between-cluster standard deviation,
which can be near 0: for some data sets the likelihood does not fall far
enough on one side for a finite limit to exist, for others OpenMx's
search does not reach a limit that exists, and some of the intervals it
does find are very wide (with 20 clusters of 10, the usable intervals of
300 data sets had a median width of 2.1, 10 of 242 were wider than 10,
and 2 wider than 100), so that their mean rests on a few data sets.

Without `assurance`, the expected width is the mean width over the data
sets that gave a usable interval, so it describes studies in which the
interval could be formed. With `assurance`, a data set without a usable
interval counts as one whose interval is wider than any width, and the
assured width is the `assurance` quantile of the widths of all `G` data
sets; the assurance is then the probability that a study gives an
interval and that the interval is no wider than the assured width. When
more than about `1 - assurance` of the data sets give no usable
interval, no width is reached with that probability.

When some of the data sets simulated at the design reported gave no
usable interval, a warning of class `dmar_replications_without_interval`
(with `scope` `"design"`) gives their number, the reasons, and the
number of data sets the width is taken over. When every data set at the
design reported gave an interval but some at other designs the search
simulated did not, a warning of the same class (with `scope` `"search"`)
gives the totals of the search, since those data sets bear on where it
stopped. When none gave one, or the assurance cannot be reached, the
function stops with an error of class `dmar_computation_failed` that
carries the counts (`attempted`, `returned`, `failed`,
`failure_reasons`), and no table is returned. Inside a search such a
design counts as one whose width does not meet the target. For a fit
that stopped with an error, the warning and the error give the first
error message, which names the cause when every data set fails for one
reason.

A warning that a fit gives is not discarded: each kind of warning is
given once when the function exits, as a warning of class
`dmar_result_qualified` that says how many simulated data sets gave it
and carries its own class as `warning_class`. Warnings of one class
whose messages differ only in their numbers count as one kind; the first
message is shown, and all of them are in the field `messages`. The
warning OpenMx gives for an optimizer status other than 0 is not passed
on: a fit of status 1 is used, as OpenMx's help describes it, and any
other status is counted as the reason the data set has no usable
interval.

## The searches

The planners that search for a size start from the size the
unstandardized difference between the condition means needs (the cluster
size search with `es_type = 1` can start lower; see below), and each
size they try is a simulation of `G` data sets. A width meets the target
when it is no larger than the target. From a start whose width meets the
target, the search steps down one cluster or one individual at a time to
the smallest size that still meets it. From a start whose width is above
the target, the search steps up 1, 3, 7, 15, and so on, past the start
until the width meets the target, and then halves the gap between the
last size that did not meet it and the first that did until the two are
adjacent. When the simulated width falls as the size grows, this gives
the smallest size that meets the target, the size that stepping up one
at a time gives, with far fewer sizes simulated. With a small `G` the
simulated width can rise and fall from one size to the next, and then
the two can stop at different sizes near the target. The same holds when
designs without a finite width (see below) lie among designs with one,
as can happen near the share of failed data sets an assurance allows:
the search returns a size that meets the target next to one that does
not, which need not be the smallest size that meets it.

The search for the number of clusters has no upper bound, since the
width falls toward 0 as clusters are added wherever it is finite. A
design at which no simulated data set gives a usable interval, or at
which too few do for the assurance, has no finite width; after such a
design the next size tried is at least twice as large, and three such
designs in a row stop the search with an error of class
`dmar_computation_failed`, with `reason` `"no_finite_width"`, that gives
the sizes, the counts, and the reasons.

With a fixed number of clusters the width levels off as the clusters
grow, because the clusters themselves vary, so a target can be out of
reach at any cluster size. For `es_type = 1` with no covariate at the
cluster level, the value it levels off at is known: as every cluster
grows, the within-cluster standard deviation becomes known and the
cluster means exact, and the likelihood-based interval for the
difference between *K* cluster means with a common unknown variance has
a width whose mean (or, with `assurance`, whose quantile) has a closed
form. A target at or below that value is refused before anything is
simulated, with an error of class `dmar_argument_invalid`, `reason`
`"width_not_reachable"`, that gives it as `width_limit`; a target above
it is reached by some cluster size. There the search starts from the
smaller of two sizes: the size the unstandardized difference needs, and
the size at which the expected width (or its quantile) with the
within-cluster standard deviation known falls below the target. That
width is the value the width levels off at, with the variance of a
cluster mean at the size in place of the between-cluster variance. The
unstandardized difference uses the *t* quantile on *K* - 2 degrees of
freedom, and with few clusters it needs far more: with 10 clusters and
an intraclass correlation of .2 it needs 1,875 individuals per cluster
for a target of 1.46, where about 10 reach it with `true_es = 0.5` (the
size the unstandardized difference needs does not depend on the effect,
but the size that reaches the target grows with it). With equal cluster
sizes, no covariate, and no `assurance`, the start is a lower bound on
the cluster size the target needs (the second of the two sizes uses
`true_R2_within` only with a covariate, since the simulated data use it
only then). The likelihood-based interval contains the interval computed
with the within-cluster standard deviation fixed at its maximum
likelihood estimate, whose expected width is the width with the
deviation known times the mean of the ratio of the deviation to its
estimate, which exceeds 1; so at every smaller cluster size the expected
width is above the target. With `es_type = 0` or `2`, or with a
covariate at the cluster level, the search starts from the size the
unstandardized difference needs, and only when that is at most 100,000
individuals per cluster. The cluster size search simulates cluster sizes
up to 100 above the size the unstandardized difference needs (above its
own start when that size is beyond 100,000) and never above 500, since
the model has one variable per individual of a cluster and a fit takes
far longer as the clusters grow; a search that would start above 500
individuals per cluster is refused before anything is simulated, with an
error of class `dmar_argument_invalid`, `reason`
`"cluster_size_beyond_reach"`, that gives the start as `start` and, as
`lower_bound`, whether it is a lower bound in the sense just given; the
message calls it one only then. A width still above the target there
stops it with an error of class `dmar_computation_failed`, `reason`
`"width_not_reached"`, that gives the sizes searched and the narrowest
width found with its Monte Carlo interval. The message says whether a
larger cluster size reaches the target (where the value the width levels
off at is known) or that this was not determined, and suggests a larger
`G` only when the target lies within that Monte Carlo interval.

The number of clusters the unstandardized difference needs, from which
the search over the number of clusters starts, is at least
`ceiling(2 / min(pr_treat, 1 - pr_treat)) + n_covariates`, which puts
two clusters in each condition (see
[`ss_aipe_crd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md)),
and the search steps down no further than the last design with two
clusters in each condition. A start beyond 1,000,000 clusters, the most
the unstandardized search considers, is refused before anything is
simulated, with an error of class `dmar_argument_invalid`, `reason`
`"n_clusters_beyond_reach"`. The planners that search over both sizes
start from the plan of the unstandardized difference, which keeps the
same floor, and do not simulate a number of clusters at which the budget
does not buy clusters of 2.

A call simulates at most 500 designs; the two planners that search over
both sizes share that bound among the searches they run, and three
neighboring designs with no finite width stop
`ss_aipe_crd_es_both_fixed_budget`. A call that reaches the bound stops
with an error of class `dmar_computation_failed`, `reason`
`"search_budget"`, that gives the totals. These bounds limit the
computation; none of them is a statistical rule.

## Refused calls

The design arguments, the costs, and the budget are checked before
anything is simulated, as in
[`ss_aipe_crd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md)
(see its section on refused calls): `n_covariates` a whole number of 0
or more (and here at most 1), `cluster_size` a whole number of at least
2, `n_clusters` a whole number of at least `3 + n_covariates` that puts
a cluster in each condition, `pr_treat` strictly between 0 and 1,
`true_R2_between` and `true_R2_within` in \[0, 1), `true_R2_between` 0
when `n_covariates = 0` (`reason` `"R2_between_without_covariates"`),
`diff_size` `NULL` or a vector of one or more finite numbers, and each
cost and the budget a single finite number of 0 or more. Wherever a
width is simulated (always in the four planners that take a width or
search over both sizes, and with `G` in the other two), `es_type` must
be 0, 1, or 2, `true_icc` as described under its argument, and, with a
covariate, `true_icc_covariate` in \[0, 1\]. Costs that are both 0 are
refused where they choose the design (`reason` `"costs_all_zero"`), and
`individual_cost = 0` where they choose the cluster size (`reason`
`"individual_cost_zero"`). A budget that does not buy the smallest
design a planner can return is refused with `reason`
`"budget_infeasible"` and the cost of that design as `min_cost` (with
`G`, the smallest design with a cluster in each condition and a degree
of freedom), and a budget whose design leaves a condition without a
cluster, when its width is asked for, with `reason`
`"arm_without_clusters"`. A fixed `n_clusters` that leaves a condition
without a cluster is refused with the same `reason` where a width is
simulated: always in `ss_aipe_crd_es_n_individuals_fixed_width`, and in
`ss_aipe_crd_es_n_individuals_fixed_budget` only with `G`, since without
it `pr_treat` is not used. Each of these has class
`dmar_argument_invalid`. `assurance` must be `NULL` or a single number
strictly between 0 and 1 in all six planners (class
`dmar_assurance_invalid`), and the three planners that take a width
refuse one that is not a single positive number.
`ss_aipe_crd_es_n_individuals_fixed_width` refuses, with `reason`
`"cluster_size_beyond_reach"` and `start = Inf`, a search whose start,
the cluster size the unstandardized difference needs, is beyond 100,000.

## References

Boker, S. M., Neale, M. C., Maes, H. H., Wilde, M., Spiegel, M., Brick,
T. R., ... Fox, J. (2011). OpenMx: An open source extended structural
equation modeling framework. *Psychometrika, 76*(2), 306–317.
[doi:10.1007/s11336-010-9200-6](https://doi.org/10.1007/s11336-010-9200-6)

Cheung, M. W.-L. (2009). Constructing approximate confidence intervals
for parameters with structural equation models. *Structural Equation
Modeling, 16*(2), 267–294.
[doi:10.1080/10705510902751291](https://doi.org/10.1080/10705510902751291)

Pornprasertmanit, S., & Schneider, W. J. (2010). *Efficient sample size
for power and desired accuracy in Cohen's d estimation in two-group
cluster randomized design* (Master Thesis). Illinois State University,
Normal, IL.

Pornprasertmanit, S., & Schneider, W. J. (2014). Accuracy in parameter
estimation in cluster randomized designs. *Psychological Methods,
19*(3), 356–379.
[doi:10.1037/a0037036](https://doi.org/10.1037/a0037036)

## See also

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two of these planners answer a question the budget alone settles. Given
# what it costs to open a cluster and what it costs to collect one more
# individual, the first reports how many clusters a budget buys at a fixed
# cluster size and the second reports how large each cluster can be at a
# fixed number of clusters. Clusters cost nothing to open here and each
# individual costs 1, so the budget buys 1000 individuals and the only
# question is how to arrange them.
ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
  cluster_cost = 0, individual_cost = 1)
#>  term                 value
#>  necessary_n_clusters 50   
#>  budget               1000 
#> 
#> Confidence level: 95%

ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 200,
  cluster_cost = 0, individual_cost = 1)
#>  term                   value
#>  necessary_cluster_size 5    
#>  budget                 1000 
#> 
#> Confidence level: 95%

# The interval width these planners can report, and every answer that
# targets a width, rests on an a priori Monte Carlo simulation: a candidate
# design is evaluated by generating G data sets and reading the
# likelihood-based confidence interval on the standardized effect size from
# OpenMx. The calls below use G = 2 so that the page runs quickly; a
# reported plan deserves the default G = 1000. Each call describes a
# population standardized effect size of 0.5, with es_type = 1 putting that
# effect size in individual-level standard deviation units and a quarter of
# the outcome variance lying between clusters.

# Supplying G and the population values to a budget planner adds the
# expected width of the interval the affordable design buys.
ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
  cluster_cost = 0, individual_cost = 1, true_es = 0.5, es_type = 1, true_icc = 0.25,
  pr_treat = 0.5, G = 2, seed = 113)
#>  term                 value
#>  necessary_n_clusters 50   
#>  expected_width       0.673
#>  budget               1000 
#> 
#> Replications: 2
#> 
#> Confidence level: 95%

# Cluster size needed for a target width, given the number of clusters.
# With 250 clusters the planner settles on the smallest cluster size it
# will consider, two individuals per cluster, and the expected width still
# comes in under the target: for a contrast between conditions that are
# assigned at the cluster level, precision is bought with clusters rather
# than with what happens inside them.
ss_aipe_crd_es_n_individuals_fixed_width(width = 0.5, n_clusters = 250,
  true_es = 0.5, es_type = 1, true_icc = 0.25, pr_treat = 0.5, G = 2,
  seed = 113)
#>  term                   value
#>  necessary_cluster_size 2    
#>  width                  0.500
#>  expected_width         0.473
#> 
#> Replications: 2
#> 
#> Confidence level: 95%

# Once recruiting a cluster costs 5, the number of clusters and the cluster
# size trade off against each other, and this planner searches the
# combinations the budget allows for the narrowest expected interval.
ss_aipe_crd_es_both_fixed_budget(budget = 1000, cluster_cost = 5,
  individual_cost = 1, true_es = 0.5, es_type = 1, true_icc = 0.25, pr_treat = 0.5,
  G = 2, seed = 113)
#>  term                   value
#>  necessary_n_clusters   111  
#>  necessary_cluster_size 4    
#>  expected_width         0.593
#>  budget                 999  
#> 
#> Replications: 2
#> 
#> Confidence level: 95%

# Unequal cluster sizes. Every planner accepts diff_size, which gives each
# cluster's deviation from cluster_size as an integer offset or as a
# multiplicative factor, recycled across the clusters; the planner prints
# the resulting cluster sizes and their frequencies ahead of its table.
# The budget planner shows both forms here.
ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
  cluster_cost = 0, individual_cost = 1, diff_size = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))
#>   cluster_size freq
#> 1           17    5
#> 2           18    6
#> 3           19    6
#> 4           20   16
#> 5           21    6
#> 6           22    6
#> 7           23    5
#>  term                 value
#>  necessary_n_clusters 50   
#>  budget               1000 
#> 
#> Confidence level: 95%

ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
  cluster_cost = 0, individual_cost = 1,
  diff_size = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
#>   cluster_size freq
#> 1           12    7
#> 2           16    6
#> 3           18    6
#> 4           20   12
#> 5           22    6
#> 6           24    7
#> 7           28    6
#>  term                 value
#>  necessary_n_clusters 50   
#>  budget               996  
#> 
#> Confidence level: 95%
```
