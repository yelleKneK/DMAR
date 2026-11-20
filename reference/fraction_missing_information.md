# Fraction of Missing Information for Each Parameter of a FIML Fit

Estimates the fraction of missing information (FMI) for every free and
every defined parameter of a model fit by full information maximum
likelihood (FIML), and states each in the sample size and interval width
terms of accuracy in parameter estimation (AIPE). The FMI of a parameter
is the proportion of the sampling variance of its estimate that is due
to the missing values (Savalei & Rhemtulla, 2012): zero when nothing is
missing, and closer to one the less the observed data say about the
parameter. Unlike the percentage of missing values, it depends on the
parameter, the model, the correlations among the variables, and the
missing data mechanism, so it is the number that says how much the
missing values cost the estimate a reader cares about, and how many more
cases a study with the same missing data would need to estimate it as
precisely as complete data would.

## Usage

``` r
fraction_missing_information(
  fit,
  information = c("unstructured", "structured", "hessian"),
  ci_method = c("none", "percentile"),
  conf_level = 0.95,
  B = 1000L,
  seed = NULL
)
```

## Arguments

- fit:

  A converged full information maximum likelihood fit: a lavaan fit with
  `missing = "ml"` (or `"ml.x"`), for example from
  `fa_confirmatory(..., missing = "fiml", output = "fit")` or
  `lavaan::cfa(..., missing = "ml")`, or a fit returned by
  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) or
  [`mlmr_mv`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md),
  whose lavaan fit is used. Requires the lavaan package.

- information:

  The estimate of the information matrix behind both standard errors:
  `"unstructured"` (the default), `"structured"`, or `"hessian"`, the
  three estimators of Chen and Savalei (2021). `"hessian"` is lavaan's
  default: for a fit made with lavaan's default information settings it
  agrees with `lavaan::parameterEstimates(fit, fmi = TRUE)` to about
  \\10^{-7}\\ when the fit runs lavaan's EM algorithm as far as the
  computation here does (`em_h1_args = list(tol = 1e-10)`), and
  otherwise to the precision of lavaan's EM, which stops earlier; at
  extreme scales lavaan's own values are the ones that are off. See
  *Details*.

- ci_method:

  `"none"` (the default) or `"percentile"`, a percentile bootstrap
  interval for each FMI. See *Details*.

- conf_level:

  Confidence level of the percentile interval (default `0.95`); used
  when `ci_method = "percentile"`.

- B:

  Number of bootstrap replications when `ci_method = "percentile"`, a
  single whole number. Default `1000`, smaller than the package's usual
  `10000` because every replication refits the model and then makes the
  two refits behind each FMI. The fewest accepted is 400 at the default
  `conf_level` (200 at 0.90, and 2000 at 0.99, which the default does
  not meet); a smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`). See *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- seed:

  Optional integer seed for the bootstrap. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  internally and the prior state restored on exit.

## Value

A `data.frame` (class `dmar_fmi`, `dmar_tbl`) with one row per free
parameter and per defined (`:=`) parameter, in the order of the lavaan
parameter table, and the same columns for every call: `term` (lavaan's
`lhs op rhs`, or the name of a defined parameter); `op` (the lavaan
operator, for filtering); `group` (the group label; `NA` for a single
group fit and for defined parameters); `fmi`; `fmi_ci_lower` and
`fmi_ci_upper` (the percentile limits; `NA` when `ci_method = "none"`);
`width_inflation`; `N`; `N_effective`; `N_target`; `estimate` (the
parameter estimate of `fit`); and `se_observed` and `se_complete`
(\\\mathrm{SE}\_O\\ and \\\mathrm{SE}\_C\\, normal theory standard
errors under `information`). The attributes record `information` and
`information_description`; `source_estimator`, `source_se`, and
`source_missing`, the settings of `fit`; `basis`; `N_used` (per group);
`N_excluded`; `n_patterns` (per group); `fmi_max`; `ci_method`; and,
with the percentile interval, `conf_level`, `B`, `B_used`, `n_dropped`,
`B_required`, and `bootstrap_status` (see *Bootstrap replications* under
Details). The print method leaves out the columns that carry nothing for
the call at hand (`op`; `group` for a single group; the interval limits
without an interval), shows an FMI within \\10^{-10}\\ of zero (a fully
observed variable's intercept, say) as 0, and states the information,
the sample size, and the imprecision of the estimate under the table.
The object keeps every column at full precision.
[`tidy.dmar_fmi`](https://yelleknek.github.io/DMAR/reference/tidy.dmar_fmi.md)
and
[`glance.dmar_fmi`](https://yelleknek.github.io/DMAR/reference/tidy.dmar_fmi.md)
give the one row per parameter and the one row summary views.

## Details

**What the fraction of missing information is.** The missing information
principle of Orchard and Woodbury (1972) splits the information that
complete data would carry about the parameters, \\J_X\\, into the
information the observed data carry, \\J_Y\\, and the information lost
with the missing values, \\J\_{X\|Y}\\: \\J_X = J_Y + J\_{X\|Y}\\. In
large samples the inverse of each information matrix is the covariance
matrix of the corresponding estimates, so for parameter \\j\\
\$\$\mathrm{FMI}\_j = 1 - \frac{\\J_X^{-1}\\\_{jj}}{\\J_Y^{-1}\\\_{jj}}
= 1 - \frac{\mathrm{SE}\_{C,j}^2}{\mathrm{SE}\_{O,j}^2},\$\$ where
\\\mathrm{SE}\_{O,j}\\ is the standard error of the FIML estimate and
\\\mathrm{SE}\_{C,j}\\ is the standard error it would have had if the
same cases had been observed completely (Savalei & Rhemtulla, 2012, Eq.
4). One minus the FMI is the asymptotic efficiency of the FIML estimate
relative to the complete data estimate. The FMI is not the proportion of
missing values. When the values are missing completely at random, the
FMI of a parameter is bounded by the largest univariate and pairwise
rates of missing values, and it falls below them when the observed
variables carry information about the missing ones; when the values are
missing at random it can be well above the missing rate (Savalei &
Rhemtulla, 2012, pp. 489–490). Savalei and Rhemtulla (2012, p. 489)
recommend reading the FMI as an effect size for the impact of the
missing values on each parameter rather than against a cutoff.

**The FMI in planning terms.** A Wald interval's width is proportional
to the standard error, and the standard error to \\1/\sqrt{N}\\, so the
precision the missing values cost a parameter can be stated as sample
sizes and widths (Savalei & Rhemtulla, 2012, p. 488). With *N* the
number of cases behind the parameter, the table reports

- `N_effective` \\= N(1 - \mathrm{FMI})\\, the number of complete cases
  that would estimate the parameter as precisely as the *N* incomplete
  cases did;

- `N_target` \\= N / (1 - \mathrm{FMI})\\, the number of cases with the
  same pattern and mechanism of missing values that would match the
  precision of *N* complete cases; and

- `width_inflation` \\= 1 / \sqrt{1 - \mathrm{FMI}} = \mathrm{SE}\_O /
  \mathrm{SE}\_C\\, how many times wider the confidence interval is than
  it would have been with complete data at the same *N*. An FMI of .5
  widens the interval by a factor of 1.41, and an FMI of .75 doubles it.

These are the statements of accuracy in parameter estimation, which
plans a sample size for the width of a confidence interval rather than
for power (Maxwell, Kelley, & Rausch, 2008), applied to the cost of the
missing values. A sample size planned for a confidence interval width
under complete data (by
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md),
or another `ss_aipe_*` planner) is, to a first approximation, divided by
one minus the FMI when the study is expected to have the same missing
data as the pilot or comparable study whose FMI is used. The
approximation is asymptotic: it treats the FMI as known and the interval
as a Wald interval.

**How the two standard errors are computed.** Both come from refits made
at the estimates of `fit` with the optimizer switched off, under normal
theory maximum likelihood standard errors and one type of information
for both, which is how lavaan computes the FMI since version 0.6-7 (Chen
& Savalei, 2021, p. 3). \\\mathrm{SE}\_O\\ comes from the fit's own data
under its missing data option. \\\mathrm{SE}\_C\\ comes from a fit to
the unrestricted FIML estimates of the means and covariances (the
saturated model's estimates, found by the EM algorithm), with the same
*N* and no rescaling of the covariance matrix. Savalei and Rhemtulla
(2012) evaluated the complete data standard errors at the model implied
moments, with expected information for the complete data and observed
information for the incomplete data; Chen and Savalei (2021, p. 3)
report that the mismatch produced negative estimates, particularly in
small samples, and the published recipe is not zero for complete data
when the model misfits. The matched computation used here gives an FMI
of zero, to rounding error, for every parameter when no value is
missing. The unrestricted moments are the FIML refit's own EM estimates,
computed in its own order of the variables, so a fit made with lavaan's
`ov_order = "data"`, which orders its variables as the data do, gets the
FMIs of the same fit made in the model's order (but see *Estimates that
are not at the maximum* for such a fit that lavaan rescales). Both
refits are made with the variables in units in which their standard
deviations are near 1, and their standard errors are carried back to the
original units exactly (see *Units of measurement*).

**Units of measurement.** The FMI is a ratio of two variances of the
same estimate, so it does not depend on the units in which the variables
are measured, and neither should its estimate. Three steps of the
computation do depend on the units when they are carried out in the
variables' own metric: lavaan's EM estimate of the unrestricted moments,
which stops when no mean or covariance changes by more than a fixed
amount, too coarse a rule for variables with small variances; the
inversion of the information matrix, which fails in floating point when
the variances span many orders of magnitude; and the Hessian estimator,
which differentiates the log likelihood numerically with a step suited
to standard deviations near 1. Both refits are therefore made with every
observed variable divided by the power of two nearest its standard
deviation (over all groups), the factors lavaan itself divides by when
it rescales a badly scaled fit, and with the parameter table carried
through the same map: every estimate and fixed value is divided by the
factor its units imply, every defined parameter and constraint is
rewritten in terms of the parameters in their original units, and the
standard errors are mapped back. Division by a power of two is exact in
floating point, so the refits hold the fit's estimates exactly, and with
no missing values their standard errors equal lavaan's own to about ten
digits. The EM is run until no mean or covariance changes by more than
\\10^{-10}\\ times the smallest variance, about \\10^{-10}\\ in that
metric, where lavaan's default is an absolute \\10^{-5}\\: EM converges
slowly where much information is missing, and stopped at lavaan's
tolerance the FMIs of a factor model of lavaan's nine Holzinger and
Swineford tests, with values missing at random on two of them, moved by
up to \\2 \times 10^{-5}\\ when the EM was run further. It is run for at
most 10,000 iterations (or the `max_iter` of the fit's `em_h1_args`, if
that is larger), and the computation warns when it stops there short of
its tolerance, which happens rarely and only where nearly all the
information about some combination of the means and covariances is
missing: in a one factor model of 12 items at *N* = 50 with 30 percent
of the values of half of them missing, one data set needed about 12,500
iterations, and its FMIs at the 10,000th were within \\3 \times
10^{-7}\\ of the converged ones (and within .002 at lavaan's tolerance).
Where the estimate does not exist the computation stops (see *Fits that
are refused*). With the fit's estimates held fixed, the three
estimators' FMIs then agreed to within \\2 \times 10^{-9}\\ when every
variable was multiplied by a constant from .001 to 1,000,000, in that
factor model and in a regression of one variable on three with values
missing at random. Computed from fits to the rescaled data, the
regression's FMIs agreed with those at a standard deviation of 1 to
within \\1.3 \times 10^{-7}\\ when one variable or every variable was
multiplied by .001, 10,000, or 100,000, and when the outcome was
multiplied by 100,000,000. lavaan's own
`parameterEstimates(fit, fmi = TRUE)`, which refits in the variables'
own metric, was off by up to .93 (by .04 to .07 for the slopes) with
every standard deviation near 100,000, and returned `NA` with every
standard deviation near 300,000, or with one near 10,000,000.

**Estimates that are not at the maximum.** The FMI is computed at the
fit's estimates, so it is the FMI of the maximum likelihood estimates
only when they are at the maximum. lavaan rescales its own fit only when
a standard deviation is .001 or less or 1000 or more, or two differ by a
factor of 1000, and a fit to variables measured in small or large units
between those limits can stop short of the maximum, although lavaan
reports it as converged and its check of the gradient, made in the
variables' own metric, passes it. With every standard deviation of that
regression near .003 the slopes were off by .14 percent and the FMIs
computed at them by up to .001 (by .00006 near .01 and .03). The factor
model of the nine tests, identified by a marker indicator for each
factor (the default of
[`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html)), with every
variable multiplied by 100 had a log likelihood .11 below the maximum
and FMIs off by up to .011, and multiplied by 500 a log likelihood 287
below it and FMIs off by .12; identified by fixing each factor's
variance (`std.lv = TRUE`) it reached the maximum at every scale from
.01 to 500. In lavaan 0.7-3 a fit made with `ov_order = "data"` and
rescaled by lavaan is not the maximum likelihood fit when the columns of
the data are not in the model's order: in such a regression the log
likelihood was 50 below the maximum, and the slopes' FMIs were .28, .22,
and .70 in place of .06, .19, and .70. The computation therefore takes
one Newton step from the fit's estimates in the metric of the refits.
Where the step moves an estimate by more than \\10^{-4}\\ of its size
(or of 1, for an estimate smaller than 1; fits at the maximum took steps
below \\6 \times 10^{-6}\\ in the examples studied here), it computes
the FMIs again one step on, at the maximum to working precision, and
warns when they differ from those at the fit's estimates by more than
\\10^{-4}\\, naming the largest change of an estimate, the gain in the
log likelihood, and the largest change of an FMI. The FMIs returned are
those at the fit's estimates. For variables in small or large units, or
in an order other than the model's, the remedy is to fit the model again
with lavaan's `rescale_data = TRUE`, which rescales whatever the units,
and with the columns of the data in the model's order: the regression
near .003 and the factor model multiplied by 100 then gave the FMIs of
the variables in their own units to within \\10^{-7}\\ and \\2 \times
10^{-6}\\. The warning names units as the cause, with that remedy, only
for a fit lavaan did not rescale with a standard deviation of .01 or
less or 100 or more (with every variable multiplied by one constant, the
factor model and the regression reached the maximum with standard
deviations from .01 to 80), and the order of the variables only for a
fit lavaan rescaled; a fit made with `rescale_data = TRUE`, or one
lavaan rescaled itself, is not told to take it again. A fit can also
stop short where the likelihood is flat in some direction, as in a small
sample with most of the outcome missing (60 cases, 16 of them observing
the outcome of a regression on four predictors), where lavaan's
optimizer stopped about \\5 \times 10^{-4}\\ standard errors short of
the maximum and neither rescaling, a tighter tolerance, nor a restart
from the estimates moved it; the FMIs were then within \\3 \times
10^{-5}\\ of the maximum's, which passes, or in another such data set
\\6 \times 10^{-4}\\ off, which the warning reports without naming units
as the cause. The check is not made for a model with an inequality
constraint or a parameter at a bound, where the gradient need not be
zero at the maximum.

The refits stay in the variables' own metric for a model with composites
(`<~`) or product terms, which lavaan does not rescale either, and for
one that shares a label between parameters of variables in different
units under `ceq.simple = TRUE` (one common factor for every variable is
tried first). For such a model the Hessian estimator warns when a
standard deviation falls outside .01 to 1000, an information matrix that
cannot be inverted stops the computation with a message that says so,
and the check for a maximum is not made; dividing the variables by
constants that bring their standard deviations nearer 1, which leaves
the FMI unchanged, is the remedy for all three.

**The three estimators and the default.** Chen and Savalei (2021)
studied three versions of the observed information, which coincide for a
saturated model such as a regression and estimate the same quantity when
the model is correct.

- `"hessian"` (their first estimator) uses the numeric Hessian of the
  log likelihood, lavaan's default. `information = "hessian"` agrees
  with `lavaan::parameterEstimates(fit, fmi = TRUE)` to about
  \\10^{-7}\\ for a fit made with `em_h1_args = list(tol = 1e-10)`, so
  that lavaan runs its EM as far, and otherwise to the precision of
  lavaan's EM, about \\10^{-5}\\ in the examples studied here (see
  *Units of measurement* above).

- `"structured"` (their second) keeps the first, dominant term of the
  observed information and evaluates it at the model implied means and
  covariances.

- `"unstructured"` (their third) has the same analytic form but uses the
  unrestricted estimates of the means and covariances wherever the
  computation needs them.

In their simulations the unstructured estimator was much less likely
than the others to break down and more precise than the structured one,
while the Hessian estimator was the most likely to break down when the
sample was small or the missing rate high (Chen & Savalei, 2021, p. 15).
It is therefore the default here, which differs from lavaan's. It also
cannot be negative, apart from rounding error: both of its informations
are evaluated at the saturated maximum likelihood estimates, where the
observed data information of the saturated model is bounded above by its
complete data information, and the bound carries through the model
Jacobian to every parameter, defined parameters included, so \\0 \le
\mathrm{FMI} \< 1\\. The structured and Hessian estimators evaluate the
information away from that point and can fall below zero, far below when
the model is misspecified: with the three-factor model for lavaan's nine
Holzinger and Swineford tests and values missing at random on six of
them, the Hessian estimates for the speed loadings are near \\-2\\ and
the unstructured ones between .17 and .34. When the model is wrong only
the Hessian estimator is consistent for its population value, an
advantage that Chen and Savalei (2021, p. 15) note has yet to be shown
to matter in practice. Under misspecification the three estimate
different quantities; the unstructured one is the information about the
parameters carried by the unrestricted moments.

**MLR and other standard errors.** For a maximum likelihood fit with
other standard errors (estimator `"MLR"`, whose standard errors are
robust to nonnormality; `"MLF"`, whose standard errors use the outer
product of the casewise gradients; or a fit with standard errors robust
to clustering and no multilevel structure) the estimates are the maximum
likelihood estimates, and only the standard errors differ. The FMI is
defined through the information of the likelihood, and the information
of the working normal likelihood is the bread of the sandwich, so the
FMI reported for such a fit is the normal theory FMI at its estimates,
the same number the ML fit gives. The source fit's standard errors are
named in the `source_se` attribute and under the printed table. A robust
counterpart would need the complete data meat of the sandwich, which
involves fourth moments of the values that were not observed and is not
identified from incomplete data without a distributional assumption. The
ratio of a normal theory complete data variance to a robust FIML
variance is not an FMI, because it attributes the inflation due to
nonnormality and misfit to the missing values. The robust and the normal
theory fractions coincide when nonnormality inflates the variances by
the same factor with and without the missing values, which is an
assumption, not a result.

**Defined parameters.** A parameter defined with `:=`, such as the
`omega_`, `H_`, and `ave_` rows that
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
defines for each factor or an indirect effect, has a delta method
standard error under each refit, so its FMI is the same ratio. For a
reliability coefficient this is the FMI that a sample size plan for its
confidence interval needs.

**Multiple groups and the *N* column.** Multiple group fits are
supported. The `N` of each row is the number of cases behind that
parameter: the group's *N* for a parameter of one group; the summed *N*
of the groups involved for a parameter constrained equal across groups
(by a shared label, an `==` constraint written with labels or with
defined parameters, or `group.equal`), whose FMI is then the same in
every group; and, for a defined parameter, the summed *N* of the groups
its expression reaches. *N* counts the cases the fit used. lavaan sets
aside a case with no observed value on any analyzed variable, as Savalei
and Rhemtulla (2012, p. 487) did in their example, and the number set
aside is the `N_excluded` attribute. In a regression fit directly with
lavaan, `missing = "ml"` with `fixed.x = TRUE` also drops every case
missing a predictor before anything is estimated; `missing = "ml.x"`, or
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md), which
models the predictors jointly by default, keeps those cases.

**Auxiliary variables.** A fit with auxiliary variables entered as
saturated correlates (Graham, 2003), as lavaan's `aux` argument enters
them (and `fa_confirmatory(..., aux = )` passes it on), is supported,
and the table includes the rows of the auxiliaries' variances, means,
and covariances. An auxiliary that shares a factor with the indicators
covaries with them through their residuals only, so the residual
covariance matrix that includes it is not positive definite, and the
covariances can imply correlations above 1 with the residuals, although
every variance is positive and the model implied covariance matrix is
positive definite. lavaan flags both, and a refit that took its checks
at face value would replace the estimates; the refits here keep them.

**Precision of the FMI itself.** A sample FMI is a noisy estimate of the
population value. In the simulations of Chen and Savalei (2021, p. 14)
it was nearly unbiased but varied greatly from sample to sample at small
*N*; with 10 to 20 percent of all values missing, samples of several
hundred gave reasonably precise estimates, and with 30 percent missing
it took more than 1,000. `ci_method = "percentile"` resamples the cases
the fit used (within groups for a multiple group fit), refits the model
with its own options, recomputes every FMI, and takes as the limits the
quantiles of the replications at `(1 - conf_level) / 2` and
`(1 + conf_level) / 2` (Efron & Tibshirani, 1993). The properties of
this interval have not been studied (Chen and Savalei, 2021, p. 15, list
bootstrap intervals for the sample FMI as future work); it shows the
sampling variability of the estimate, and its upper limit is the
conservative value to carry into a sample size plan.

**Bootstrap replications.** The two limits of each parameter's interval
are quantiles of the `B` replicated FMI values, so the interval is run
only when `B` puts at least 10 replications beyond each limit (400 at
the 95 percent level). A replication in which the FMI of some parameter
is undefined (the refit fails, does not converge, or returns an FMI that
is not a finite number) is dropped and counted, never replaced, and the
interval is computed from the rest. The result records `B` (requested),
`B_used` (the number the interval rests on), `n_dropped`, and
`B_required` (the fewest replications a request may ask for; that count
is held to `B`, not to `B_used`) as attributes, with `bootstrap_status`
(`"complete"`, `"dropped"`, `"exploratory"`, `"insufficient_count"`, or
`"insufficient_share"`). One warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the 2.5 percent of `B` that lies beyond one limit at the 95
percent level, because they could then decide by themselves where that
limit falls. When more than half were dropped, the limits of every
parameter are `NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`.
`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration; every such call gives a
warning of class `dmar_bootstrap_exploratory` and the result is not one
to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

Each replication costs a model fit and two refits, so the interval takes
about a minute at the default `B` for a model the size of the one in the
examples; the call is
`fraction_missing_information(fit, ci_method = "percentile", seed = 113)`.

**Fits that are refused.** The computation needs a converged and
unweighted maximum likelihood fit with one level, to continuous
variables, that used every observed value. A fit with ordered indicators
(and so WLSMV or another least squares estimator), with an estimator
outside the maximum likelihood family, with listwise deletion (lavaan's
default, which leaves no missing information to measure), with the
pairwise, available cases, or two-stage missing data options, with
`conditional.x = TRUE`, with more than one level, with sampling weights,
or that did not converge is refused with a message that says why, as are
a `lavaanList` or pooled multiple imputation object and a rotated
exploratory (`efa`) solution, whose loadings lavaan rotates after
estimation. As a last guard, the computation stops if the refit at the
estimates moves them. Besides a rotated solution, this happens with an
improper one: a Heywood case whose negative residual variance exceeds
the variable's observed variance in size is a solution lavaan (0.7-3)
reads as having run away, so it will not take those estimates as the
starting point of a refit and starts again from its own simple starting
values. The standard errors of such a refit belong to those values, not
to the fit, and `lavaan::parameterEstimates(fit, fmi = TRUE)`, which
makes the same refit, then reports FMIs computed at them without saying
so. Every other solution lavaan returns is refit at its estimates,
whatever its checks of a solution say of them, and its FMIs are computed
there: a Heywood case whose residual variance is less far below zero,
two factors correlated above 1 (a covariance matrix of the latent
variables that is not positive definite), and a covariance between
residuals that implies a correlation above 1. Before it refits,
`lavaan::parameterEstimates(fit, fmi = TRUE)` shrinks such a covariance
to a correlation below 1 (its check of starting values), so it reports
the FMIs of other estimates: for two factors correlated 1.08 in a sample
of 60, it gave the covariance an FMI of -1.62, where the FMI at the
estimates is .25. These are the FMIs of an improper solution, and they
deserve the same caution as the solution itself, of which lavaan warns
when it fits the model. Finally, the computation stops when the
information matrix of a refit cannot be inverted, with a message that
names the likeliest cause: variables too badly scaled for a model the
refits keep in their own metric (see *Units of measurement*); a
predictor (an observed variable regressed on) that lies so far from
zero, relative to its standard deviation, that an intercept and its
coefficients are nearly collinear, which no division by a power of two
changes (a predictor of a regression with a standard deviation of 1 did
so from about 8000 from zero, and gave the FMIs of the unshifted data at
7000; an indicator's or an outcome's location enters only its own
intercept and does not); a model that is not identified; or a structured
or Hessian information that is singular at the estimates. A predictor
1000 or more standard deviations from zero is named only when the same
model fit to the data with it centered has an information that can be
inverted, since centering, which changes no FMI but those of intercepts,
cannot repair a model that is not identified; where centering would
change a fixed or constrained intercept, so that this cannot be checked,
both causes are given.

The computation also stops when the data define no unrestricted estimate
of the means and covariances at which to evaluate the fraction (all
three estimators need it: the unstructured information is evaluated
there, and the complete data side of each is fitted to it). A pair of
variables that no case observes together, as in a planned missing
design, leaves their covariance undetermined: lavaan's EM leaves it at
its starting value and only warns, and in a two-form design on seven of
the Holzinger and Swineford tests the FMIs moved by up to .12 across
values of it that the data could not tell apart. The pairs are named, as
[`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md) and
the reliability family name them. The computation also checks for a
supported affine-collapse direction: all rows observing a subset must
satisfy the same affine relation, with a nonzero coefficient on every
member of the subset. Every partially overlapping row is included when a
candidate support is reduced. A subset with no more observed rows than
variables and a verified full-support null direction supplies a
sufficient unbounded-likelihood diagnosis, even when EM stopped at a
stable interior point. The check considers the full variable set,
observed-pattern sets, single variables, and reductions of their null
supports; it is not an exhaustive existence test.

Affine differences are column scaled. A singular-value tolerance of
\\100\epsilon\max(n,p)\max(1,d\_{\max})\\, with machine precision
\\\epsilon\\, locates candidate nullspaces. A deterministic normalized
direction must have every coefficient larger in absolute value than
\\p\sqrt{\epsilon}\\ and a residual below that tolerance. A relation
inferred only numerically in an overdetermined subset is reported as
numerical instability, without claiming mathematical nonexistence. A
near-singular iterate at EM's iteration limit is a fitting failure. No
diagnosed obstacle means that no implemented diagnostic failed; it does
not certify finite global existence or uniqueness. These checks concern
the unrestricted likelihood needed for FMI. The constrained substantive
model's own estimates and standard errors do not use the unrestricted
estimate (under lavaan's default information); its chi-square test and
the fit indices built on the unrestricted model's likelihood do. These
stops are errors of classes `dmar_fmi_never_together` and
`dmar_fmi_saturated_singular`, so that a caller can catch them.

## References

Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
missing information from full information maximum likelihood. *Frontiers
in Psychology, 12*, 667802.
[doi:10.3389/fpsyg.2021.667802](https://doi.org/10.3389/fpsyg.2021.667802)

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Orchard, T., & Woodbury, M. A. (1972). A missing information principle:
Theory and applications. In *Proceedings of the Sixth Berkeley Symposium
on Mathematical Statistics and Probability* (Vol. 1, pp. 697–715).
Berkeley, CA: University of California Press.

Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
fraction of missing information from full information maximum
likelihood. *Structural Equation Modeling, 19*(3), 477–494.
[doi:10.1080/10705511.2012.687669](https://doi.org/10.1080/10705511.2012.687669)

## See also

[`missing_patterns`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)
and
[`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md)
for the description of the missing values that comes before any model;
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
and [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) for
FIML fits;
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
for planning the sample size of a reliability interval.

Other missing data:
[`cov_fiml()`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md),
[`df_barnard_rubin()`](https://yelleknek.github.io/DMAR/reference/df_barnard_rubin.md),
[`df_missing()`](https://yelleknek.github.io/DMAR/reference/df_missing.md),
[`missing_coverage()`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md),
[`missing_patterns()`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The six spatial tests of the Holzinger and Swineford data. The revised
# paper form board (t25) and the flags test (t26) were given only at
# Grant-White, so both are missing for the 156 Pasteur pupils.
spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
missing_patterns(holzinger_swineford, items = spatial)
#>  pattern missing_items                     n   proportion n_observed n_missing
#>  111100  t25_paper_form_board_r, t26_flags 156 0.518      4          2        
#>  111111                                    145 0.482      6          0        

fit <- fa_confirmatory(data = holzinger_swineford,
             factors = list(spatial = spatial),
             missing = "fiml", output = "fit")
fmi <- fraction_missing_information(fit)
fmi
#>  term                                             fmi   width_inflation N  
#>  spatial =~ t1_visual_perception                  0.113 1.061           301
#>  spatial =~ t2_cubes                              0.063 1.033           301
#>  spatial =~ t3_paper_form_board                   0.075 1.040           301
#>  spatial =~ t4_lozenges                           0.111 1.061           301
#>  spatial =~ t25_paper_form_board_r                0.512 1.431           301
#>  spatial =~ t26_flags                             0.506 1.423           301
#>  t1_visual_perception ~~ t1_visual_perception     0.214 1.128           301
#>  t2_cubes ~~ t2_cubes                             0.058 1.030           301
#>  t3_paper_form_board ~~ t3_paper_form_board       0.120 1.066           301
#>  t4_lozenges ~~ t4_lozenges                       0.203 1.120           301
#>  t25_paper_form_board_r ~~ t25_paper_form_board_r 0.517 1.438           301
#>  t26_flags ~~ t26_flags                           0.511 1.430           301
#>  t1_visual_perception ~1                          0.000 1.000           301
#>  t2_cubes ~1                                      0.000 1.000           301
#>  t3_paper_form_board ~1                           0.000 1.000           301
#>  t4_lozenges ~1                                   0.000 1.000           301
#>  t25_paper_form_board_r ~1                        0.444 1.341           301
#>  t26_flags ~1                                     0.430 1.325           301
#>  loading_sum_spatial                              0.184 1.107           301
#>  error_sum_spatial                                0.348 1.238           301
#>  omega_spatial                                    0.244 1.150           301
#>  ave_spatial                                      0.221 1.133           301
#>  H_spatial                                        0.216 1.130           301
#>  N_effective N_target estimate se_observed se_complete
#>  267.136     339.157  4.827    0.430       0.405      
#>  282.009     321.270  2.174    0.302       0.292      
#>  278.550     325.259  1.475    0.179       0.172      
#>  267.508     338.685  5.859    0.563       0.530      
#>  147.001     616.331  1.764    0.267       0.186      
#>  148.635     609.555  5.313    0.707       0.497      
#>  236.584     382.955  25.601   3.162       2.803      
#>  283.633     319.431  17.383   1.556       1.510      
#>  264.957     341.947  5.809    0.542       0.509      
#>  239.901     377.660  47.259   5.296       4.728      
#>  145.518     622.612  6.435    0.848       0.589      
#>  147.213     615.442  41.656   5.913       4.136      
#>  301         301      29.615   0.403       0.403      
#>  301         301      24.352   0.271       0.271      
#>  301         301      14.229   0.163       0.163      
#>  301         301      18.003   0.521       0.521      
#>  167.451     541.060  15.747   0.239       0.178      
#>  171.448     528.447  36.602   0.643       0.485      
#>  245.688     368.764  21.412   1.291       1.166      
#>  196.282     461.586  144.144  7.733       6.244      
#>  227.492     398.261  0.761    0.026       0.022      
#>  234.386     386.546  0.352    0.029       0.026      
#>  235.888     384.085  0.775    0.023       0.020      
#> 
#> FMI under the unstructured information (the third estimator of Chen & Savalei, 2021), from 301 cases.
#> The sample FMI is imprecise unless N is several hundred (Chen & Savalei, 2021); ci_method = "percentile" gives an interval.

# The loadings of t25 and t26 have an FMI near .51, close to their 52
# percent of values missing, so each is estimated about as precisely as
# 147 complete cases would estimate it. Omega for the six tests loses
# less, about a quarter of its information: its interval is 1.15 times
# as wide as complete data would give, and the same precision with the
# same missing values would take about 398 cases.
fmi[fmi$term == "omega_spatial",
    c("term", "fmi", "N_effective", "N_target", "width_inflation")]
#>  term          fmi   N_effective N_target width_inflation
#>  omega_spatial 0.244 227.492     398.261  1.150          

# The Hessian estimator, lavaan's default, for comparison.
fraction_missing_information(fit, information = "hessian")$fmi
#>  [1]  0.127929413  0.073126491  0.077115264  0.121462515  0.512742416
#>  [6]  0.512890031  0.177514117  0.068459692  0.031844452  0.159697206
#> [11]  0.523894818  0.519255529 -0.003640621 -0.001612717 -0.002052375
#> [16] -0.003261740  0.457923255  0.442159315  0.201276715  0.300807648
#> [21]  0.246288572  0.231103703  0.225617363
```
