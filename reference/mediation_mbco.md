# Mediation Analysis via Model-Based Constrained Optimization

Tests hypotheses about mediation effects with the model-based
constrained optimization (MBCO) procedure of Tofighi and Kelley (2020):
a likelihood ratio test that compares the full mediation model to a null
model fit subject to the nonlinear constraint that the effect of
interest (an indirect effect, a total effect, or any smooth function of
path coefficients) equals zero. The model is specified in lavaan syntax
but is fit with OpenMx, whose optimizers support the nonlinear equality
constraints the procedure requires (see Details). The function
enumerates the mediation effects implied by the model (the total effect,
the direct effect, the total indirect effect, and every specific
indirect pathway from `x` to `y`), reports each with a confidence
interval, and tests each with the MBCO likelihood ratio statistic and
its *p*-value. Models may include observed or latent variables, one or
many mediators, in parallel or in sequence, and raw data or summary
statistics may be supplied. With a grouping variable the model becomes a
multiple-group SEM: every effect is estimated within each group and the
between-group difference of each effect is tested, which is moderated
mediation with a categorical moderator. With a numeric `moderator`,
interactions stated in the model are probed: conditional effects at
chosen moderator values, the index of moderated mediation, and a joint
constancy test, all with MBCO likelihood ratio inference (see Details).

## Usage

``` r
mediation_mbco(
  model,
  data = NULL,
  S = NULL,
  M = NULL,
  N = NULL,
  group = NULL,
  x = NULL,
  y = NULL,
  moderator = NULL,
  probe_values = NULL,
  hypotheses = NULL,
  ci_method = c("profile_likelihood", "monte_carlo", "wald"),
  conf_level = 0.95,
  B = 10000,
  optimizer = c("SLSQP", "CSOLNP", "NPSOL"),
  seed = NULL
)
```

## Arguments

- model:

  A character string with the model in lavaan syntax (regressions `~`,
  latent variable definitions `=~`, residual covariances `~~`,
  optionally labeled coefficients and defined quantities via `:=`), or a
  fitted lavaan object whose parameter table and data are then reused.
  Parameter labels must be valid R names (e.g., `b1`, not `1b`).

- data:

  A `data.frame` of raw data containing the observed variables. Rows
  with missing values are retained and handled by full information
  maximum likelihood. Supply either `data` or `S` (with `N`), not both.

- S:

  A covariance matrix of the observed variables (with dimnames naming
  the variables), used with `N` and optionally `M` when raw data are not
  available. Complete-data maximum likelihood depends on the data only
  through the sample means and covariance matrix, so an analysis from
  summary statistics reproduces the raw-data analysis exactly (see
  Details). For a multiple-group analysis, a named list of such
  matrices, one per group.

- M:

  Optional named numeric vector of variable means accompanying `S` (for
  multiple groups, a named list of such vectors). When omitted, means of
  zero are used; intercepts are then uninteresting but every effect,
  test, and confidence interval is unaffected.

- N:

  Total sample size. Required with `S`; ignored for raw data. For a
  multiple-group analysis, a vector with one sample size per group in
  `S` (matched by name when named).

- group:

  Optional name (single character string) of the grouping variable in
  `data` for a multiple-group analysis. With summary input, supply `S`,
  `M`, and `N` as lists instead and leave `group` alone (it then merely
  names the constructed grouping column). See the Details section on
  multiple groups, including how lavaan's labeling rules decide which
  parameters differ by group.

- x, y:

  Names (single character strings) of the focal predictor and the
  outcome between which effects are traced. When both are omitted the
  function looks for a unique variable with no incoming regression path
  (the source) and a unique variable with no outgoing regression path
  (the outcome); if either is ambiguous an error asks for `x` and `y`
  explicitly. When the model defines quantities via `:=` or `hypotheses`
  are supplied, `x` and `y` may be omitted entirely and only those
  quantities are estimated and tested.

- moderator:

  Optional name (single character string) of a numeric moderator
  variable, turning on probing: every pathway effect that involves an
  interaction with the moderator is estimated at each probe value, the
  change in the effect per unit moderator is reported (for a pathway
  moderated in one place, the index of moderated mediation), and a
  constancy test asks whether the pathway effect depends on the
  moderator at all. State the interaction in the model syntax with a `:`
  term (e.g., `m ~ x + w + x:w`) or include a precomputed product column
  among the predictors; both are recognized. A categorical moderator
  goes through `group` instead.

- probe_values:

  Numeric vector of moderator values at which moderated effects are
  estimated; names, when given, label the rows. Defaults to the
  moderator's mean and one standard deviation either side (labeled
  `low`, `mean`, `high`), or to the two observed values of a two-valued
  moderator.

- hypotheses:

  Optional named character vector of additional quantities to estimate
  and test against zero, written as R expressions in the parameter
  labels and effect names (e.g.,
  `c(difference = "indirect_via_imagery - indirect_via_repetition")`
  tests the equality of two specific indirect effects). A named *list*
  may mix such scalar expressions with character vectors of several
  expressions; a multi-expression element is tested jointly (all
  constrained to zero at once, with as many degrees of freedom as
  expressions), and its row reports the test alone, with no scalar
  estimate. Each expression of a joint test must restrict the parameters
  beyond what the others do. An expression given twice, a multiple of
  another (`"a"` and `"2*a"`), a combination of others (`"a"`, `"b"`,
  and `"a + b"`; `"a - b"` and `"b - a"`), or an expression that is zero
  wherever the others are (`"a*b"` with `"b"`; `"a*cp"` with `"a"`;
  `"a*b*cp"` with `"a*b"`) would be counted as a degree of freedom while
  restricting nothing, so such a set is refused before any model is fit,
  with a condition of class `dmar_argument_invalid` that names the
  redundant expressions (as `redundant`). Equivalent sets of expressions
  that each restrict, such as `c("a", "b")` and `c("a + b", "a - b")`,
  give the same test. See *Constraints that add no restriction* under
  Details.

- ci_method:

  Confidence interval method for every reported effect:
  `"profile_likelihood"` (default; Neale & Miller, 1997, computed by
  OpenMx), `"monte_carlo"` (simulate the coefficients from their joint
  normal approximation and form percentile intervals of the effect;
  MacKinnon, Lockwood, & Williams, 2004), or `"wald"` (delta method
  symmetric interval; reported for comparison, not recommended for
  indirect effects, whose sampling distributions are skewed).

- conf_level:

  Confidence level. Defaults to 0.95.

- B:

  Number of Monte Carlo draws when `ci_method = "monte_carlo"`, a single
  whole number (not read, and not checked, for the other intervals).
  Defaults to 10000, because a draw costs one evaluation of each effect
  and no refit. The fewest accepted is 400 at the default `conf_level`
  (200 at 0.90, 2000 at 0.99); a smaller `B` is refused before the model
  is fit (class `dmar_count_invalid`). See *Monte Carlo draws* under
  *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- optimizer:

  Optimizer used by OpenMx for all fits: `"SLSQP"` (default),
  `"CSOLNP"`, or `"NPSOL"`. The CRAN build of OpenMx ships SLSQP and
  CSOLNP; NPSOL, the optimizer used in Tofighi and Kelley (2020), is
  available only in the build distributed from the OpenMx website.

- seed:

  Optional integer seed for the Monte Carlo interval, used locally (the
  caller's random number generator state is restored on exit). Default
  `NULL` draws the Monte Carlo replications from the current random
  number stream, which advances that stream.

## Value

A `data.frame` (classes `dmar_mediation_mbco`, `dmar_tbl`) with one row
per effect and columns `pathway` (the traced pathway or defining
expression), `term`, `estimate`, `se` (delta method), `ci_lower`,
`ci_upper`, `lrt` (the MBCO likelihood ratio statistic), `df`,
`p_value`, `delta_aic`, and `delta_bic`. A joint test row (a
several-constraint moderation or hypothesis test) reports the test
columns with `df` equal to the number of constraints, none of which is
redundant (see `hypotheses`); its `estimate`, `se`, and interval are
`NA` because no single number summarizes several constraints.
Attributes: `"conf_level"`, `"ci_method"`, with the Monte Carlo interval
`"B"`, `"B_used"`, `"n_dropped"`, `"B_required"`, and
`"bootstrap_status"` (see *Monte Carlo draws* under Details),
`"optimizer"`, `"x"`, `"y"`, `"N"`, `"deviance"`, `"aic"`, `"bic"`,
`"n_par"` (full model fit information), `"groups"` (the group labels,
reference first; multiple-group fits only), `"R2"` (named vector,
endogenous variables under the full model, per group when grouped),
`"delta_R2"` (matrix of full minus null \\R^2\\, variables by tested
effects), `"moderation"` (for a moderated analysis: the moderator name,
the probe values, and each moderated effect's moderator polynomial,
which is what
[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md)
draws), `"rescale"` (only when badly scaled variables were rescaled for
the fits, see *Badly scaled variables* under Details: `divisors`, the
power of 2 each observed variable was divided by, and `multipliers`, the
factor that turns each parameter of the rescaled fit into the parameter
of the original metric), and `"mx_model"` (the fitted OpenMx full model,
in the rescaled metric when there is a `"rescale"` attribute, an escape
hatch for further OpenMx work). Use `tidy()` and `glance()` for
broom-style views.

## Details

**The MBCO procedure.** A hypothesis about mediation is recast as a
model comparison (Maxwell, Delaney, & Kelley, 2027). The full model
estimates all path coefficients freely. The null model is the same model
estimated subject to the constraint that the effect of interest is zero,
for example \\H_0\\: \beta_1 \beta_2 = 0\\ for the indirect effect of
\\X\\ on \\Y\\ through one mediator. Writing \\D =
-2\\\mathrm{log\\likelihood}\\ for the deviance of each fitted model,
the test statistic is \$\$\mathrm{LRT}\_{\mathrm{MBCO}} =
D\_{\mathrm{null}} - D\_{\mathrm{full}},\$\$ and the *p*-value refers it
to the chi square distribution with degrees of freedom equal to the
number of constraints (here 1). That distribution is the large sample
distribution of the statistic when the null hypothesis is true (Wilks,
1938) and regular: the constraints are smooth and functionally
independent, and the population parameter is identified, lies in the
interior of the parameter space, and is a point near which the null set
is a smooth surface with one dimension fewer per constraint than the
full model. The null hypothesis of no indirect effect is regular at
every point of its null set except those where both of its coefficients
are zero (whatever the values of the other parameters), where the chi
square reference does not hold (see *The reference distribution where
the branches cross* below). The *p*-value is a continuous measure of the
compatibility of the data with the null model, not merely a reject or
fail-to-reject decision. In the simulations of Tofighi and Kelley (2020,
p. 508), with one coefficient of the indirect effect zero, the Type I
error rate of the test stayed close to the nominal .05 for every sample
size (50 to 200) and size of the other coefficient, whereas the
percentile bootstrap, Monte Carlo, and profile likelihood intervals and
the joint significance test rejected too rarely when the other
coefficient was small (0.14) and the sample size was 100 or less.

**Why the constraint is nonlinear.** Null hypotheses about single
parameters are linear constraints: fixing \\\beta_3 = 0\\ restricts the
parameter space to a flat hyperplane, something every structural
equation modeling program does by fixing the parameter and refitting.
The null hypothesis of no mediation is different in kind. The constraint
\\\beta_1 \beta_2 = 0\\ involves a product of parameters, so it is a
nonlinear function of the parameter vector, and its solution set is not
a hyperplane but the union of two hyperplanes: the set where \\\beta_1 =
0\\ (with \\\beta_2\\ free) and the set where \\\beta_2 = 0\\ (with
\\\beta_1\\ free). Mediation is absent in infinitely many ways, and no
single fixed parameter expresses them all: fixing \\\beta_1 = 0\\ alone
tests a more restrictive hypothesis than \\H_0\\: \beta_1 \beta_2 = 0\\.
Maximizing the likelihood over such a constraint set requires an
optimizer that handles general nonlinear equality constraints
(sequential quadratic programming methods such as SLSQP and NPSOL, or
CSOLNP); among R packages for structural equation modeling, OpenMx
provides them, which is why the model is fit there even though it is
specified in lavaan syntax. The same geometry explains why Wald-type
(delta method) tests of a product behave badly near the null: at
\\\beta_1 = \beta_2 = 0\\ the two hyperplanes intersect, the gradient of
\\\beta_1 \beta_2\\ vanishes, and the usual normal approximation for the
product breaks down. The MBCO procedure does not use that approximation,
since it compares maximized likelihoods directly; the chi square
reference of its statistic is itself a large sample approximation,
though, and it too fails at that point, as the next paragraph explains.

**The reference distribution where the branches cross.** Where exactly
one of \\\beta_1\\ and \\\beta_2\\ is zero, the null set near the
population values is the single hyperplane of that coefficient, the null
hypothesis is regular there, and the chi square distribution on 1 degree
of freedom is the large sample distribution of the statistic. At
\\\beta_1 = \beta_2 = 0\\ the two hyperplanes cross and the null
hypothesis is not regular. The best-fitting null model lies on one
branch or the other, so the statistic is the smaller of the two branch
statistics, each the likelihood ratio statistic for fixing one
coefficient at zero. When the estimates of the two coefficients are
independent, as in the single mediator model with observed variables,
the branch statistics are asymptotically the squares \\Z_1^2\\ and
\\Z_2^2\\ of two independent standard normal variables, and the
statistic is asymptotically \\\min(Z_1^2, Z_2^2)\\ rather than a chi
square variable on 1 degree of freedom. It exceeds the .95 quantile of
that chi square distribution with probability \\.05^2 = .0025\\. At
\\\beta_1 = \beta_2 = 0\\, then, the test at the .05 level rejects with
probability near .0025 in large samples, not .05. When the two estimates
are correlated, as they can be in models with latent variables and in
some multiple-group or moderated designs, \\Z_1\\ and \\Z_2\\ are
correlated standard normal variables, and the probability lies between
.0025 and .05, nearer .05 the stronger the correlation (.0049, .0125,
and .0296 for correlations of .3, .6, and .9). The test is conservative
at the crossing in every case: its *p*-value is larger than one computed
from the distribution that holds there, and its power against an
indirect effect whose coefficients are both near zero is lower
accordingly. For a chain of three or more coefficients the same holds
wherever two or more of them are zero, and the test is more conservative
still where more of them are zero (with independent estimates, \\.05^3\\
when three are zero). The function reports the chi square *p*-value
throughout and does not adjust it, because which point of the null set
holds is not known, and a critical value taken from the distribution at
the crossing would make the test reject too often where only one
coefficient is zero. The same qualification applies to every test of a
product of coefficients the function reports (specific indirect effects,
conditional effects at a probe value, and moderation coefficients that
are products), and to a joint test whose constraints are satisfied
together at a point where their gradients are linearly dependent, such
as the crossing of two branches.

**The rate at the crossing in the simulations of Tofighi and Kelley
(2020).** Their simulations report a different rate at \\\beta_1 =
\beta_2 = 0\\. For the single mediator model with \\X\\ binary
(probability .5) and standard normal residuals, 5,000 replications at
each of *N* = 50, 75, 100, and 200, the .05 level, and their test fitted
by the NPSOL optimizer of OpenMx 2.9.6 (p. 507 and its footnote 4), the
empirical Type I error rate of \\\mathrm{LRT}\_{\mathrm{MBCO}}\\ with
both coefficients zero fell within the interval from .025 to .075 that
they took as close enough to .05 for practical purposes, at every sample
size (p. 508; about .03 to .04 in their Figure 3, p. 509), and the other
five methods they studied rejected at rates below .01. The two numbers
describe different things. The .0025 above is a limit: the rate that the
statistic defined by the best-fitting null model approaches as *N* grows
at that point. The published rates are estimates from a simulation of
one design at four sample sizes. Neither is the rate at every sample
size. For the statistic this function computes, the rate in that design
is already close to the limit at these sample sizes: in 400,000
simulated data sets at each *N*, the smaller of the two branch
statistics (the statistic `mediation_mbco()` returns for this model, to
within \\10^{-4}\\ on 630 further data sets of the design at *N* = 50 to
200) exceeded the .95 chi square quantile at rates of .0032, .0029,
.0030, and .0028 for *N* = 50, 75, 100, and 200 (Monte Carlo standard
error .00009), so sample size in this range does not account for the
difference. The statistic of the worse-fitting branch, the one a search
reports when it stops on that branch (see *Local solutions and the
search strategy* below), rejected at .100 to .111 in the same
simulations, and the published rates lie between the two.

**Local solutions and the search strategy.** Because the null set is a
union of surfaces, the constrained deviance surface generally has one
local minimum per branch (one where the mediator does not respond to
\\X\\, one where the outcome does not respond to the mediator, and so on
along longer chains). A constrained optimizer started from the
full-model estimates can converge to whichever branch is nearest rather
than to the branch that fits best. The likelihood ratio test is defined
by the globally best-fitting null model, so `mediation_mbco()` refits
each null model from several starting configurations (the full-model
estimates, and the estimates with each coefficient entering the
constrained effect set to zero in turn), verifies that each candidate
solution actually satisfies the constraint, and keeps the feasible
solution with the smallest deviance. In the memory example below this
matters: for the single mediator model the best-fitting null model sets
the imagery-to-recall path to zero (\\\mathrm{LRT}\_{\mathrm{MBCO}} =
71.31\\, the statistic the example reports), while the branch that sets
the instruction-to-imagery path to zero fits worse
(\\\mathrm{LRT}\_{\mathrm{MBCO}} = 179.02\\). Tofighi and Kelley (2020)
report 175.77 for this example, a value from that worse-fitting branch,
on which their optimizer stopped; their statistic also differs from the
179.02 here because the example runs from the published (rounded)
summary statistics of their Table 1 rather than the full-precision
moments (which give 72.54 and 175.77 for the two branches). Either way
the null model is overwhelmingly incompatible with the data, so the
substantive conclusion is the same.

The null model of a product is fit branch by branch: a product is zero
exactly where one of its factors is, so for each way of choosing one
factor of each constraint, a factor that is a coefficient is fixed at
zero and a factor that is a sum (the conditional coefficient \\a_0 + a_1
w\\ of a moderated path at a probe value \\w\\, for example) is kept as
the constraint of that fit. A label common to every term of a sum is
taken out as a factor (\\a_0 b + a_1 b w = b (a_0 + a_1 w)\\). The null
set is the union of the branches, so when every branch fit is a valid
solution the smallest of them is the constrained minimum, and the
constrained searches from the starting configurations are run only when
a branch fit is not, or when a constraint has no branches (a total
effect, a difference between groups). Such a constraint is also searched
from the full-model estimates moved onto its null set by Gauss-Newton
steps.

**Badly scaled variables.** When the observed variables differ in scale
by orders of magnitude (the ratio of the largest standard deviation to
the smallest is at least 1000, or the largest is at least 1000, or the
smallest at most 0.001, the rule lavaan uses for its `rescale_data`
option), the free parameters do too, and the OpenMx optimizers can stop
far from the maximum. The OpenMx models are then fit to the variables
each divided by a power of 2 near its standard deviation
(\\2^{\mathrm{round}(\log_2 s)}\\ for standard deviation \\s\\, over all
groups), with every parameter in the matching metric (a latent variable
takes the metric of the indicator whose loading is fixed), and every
effect, constraint, and Monte Carlo draw is written in the parameters of
the original metric, so the estimates, standard errors, intervals, and
likelihood ratio statistics are those of the variables as given, and the
deviance, AIC, and BIC are returned to the original metric as well. A
likelihood ratio statistic does not depend on the units of the
variables, and the rescaling changes it only within the precision of the
optimizer. A profile likelihood interval is profiled for the effect in
the metric of the fit and its limits returned to the original metric.
The constraints of a joint test are then handed to the optimizer as
linear combinations whose gradients are orthonormal in the metric of the
fit at the full-model estimates, which leaves their null set unchanged.
Nothing is rescaled when a label shared by parameters in different
metrics (an equality constraint across them) would need two multipliers.
The divisors and the multipliers are returned in the `"rescale"`
attribute. Scales that differ by less than that rule can still slow the
optimizers; the branch fits and the start on the null set described
above are the remedy there, and a variable can be put on a scale near
that of the others before the call (the statistics do not depend on its
units, and the estimates change by the factor applied).

**Effects estimated and tested.** With `x` and `y` resolved, the
function enumerates every directed pathway from `x` to `y` along
regression (`~`) paths. Each pathway through at least one intermediate
variable contributes a specific indirect effect, the product of its path
coefficients, reported as `indirect_via_` followed by the intermediate
variable names. The direct effect is the `x` to `y` coefficient when
that path is in the model. The total indirect effect is the sum of the
specific indirect effects (reported when there are two or more), and the
total effect is the direct effect plus the total indirect effect
(reported when the direct path is in the model). Quantities defined in
the model syntax via `:=` and any `hypotheses` are estimated and tested
the same way, so contrasts of indirect effects, proportions mediated, or
any other smooth function of parameters can be examined; each is tested
against zero, so write an equality of two effects as their difference.
Every reported row carries its estimate, a delta method standard error
(via [`mxSE`](https://rdrr.io/pkg/OpenMx/man/mxSE.html)), the
`ci_method` confidence interval, and the MBCO likelihood ratio test with
its degrees of freedom and *p*-value.

**The parallel two-mediator model.** Tofighi and Kelley (2020) continue
the memory example with imagery and repetition as parallel mediators and
a residual covariance between them, asking whether the indirect effect
through repetition is zero and whether the two specific indirect effects
differ (their Research Questions 2 and 3). That analysis is the
single-mediator call of the examples with the parallel model in place of
the single one and the contrast supplied through `hypotheses`:

    parallel <- "
      imagery    ~ b1*instruction
      repetition ~ b3*instruction
      recall     ~ b2*imagery + b4*repetition + b5*instruction
      imagery ~~ repetition
    "
    mediation_mbco(parallel, S = S_tk, M = M_tk, N = 369,
                   x = "instruction", y = "recall",
                   hypotheses = c(imagery_minus_repetition =
                     "indirect_via_imagery - indirect_via_repetition"),
                   ci_method = "monte_carlo", seed = 113)

Six effects are reported, so six constrained null models are fit, about
twice the cost of the single-mediator analysis. The indirect effect
through repetition is near zero (the paper reports
\\\mathrm{LRT}\_{\mathrm{MBCO}} = 0.083\\, *p* = .773), while the
contrast shows the imagery pathway is larger (the paper reports
\\\mathrm{LRT}\_{\mathrm{MBCO}} = 25.828\\, difference = 2.222, SE =
0.445).

**Choosing the confidence interval.** The default profile likelihood
interval inverts the likelihood ratio test for the effect itself (Neale
& Miller, 1997), so its limits are free to sit asymmetrically about the
estimate, as the skewed sampling distribution of a product of
coefficients calls for, and the interval and the MBCO test tell the same
story. Each bound is a constrained search of its own, which is what
makes it the expensive choice. The Monte Carlo interval is the
inexpensive alternative that also accommodates the skewness: `B`
coefficient vectors are drawn from the joint normal approximation of the
estimates, the effect is evaluated in each draw, and the limits are the
empirical \\(\alpha/2, 1 - \alpha/2)\\ quantiles of those `B` values
(MacKinnon, Lockwood, & Williams, 2004). Ask for it when the profile
searches are slow, when a profile bound fails to converge, or when
comparing with a published analysis that reports one, as Tofighi and
Kelley (2020) do for the memory data in the examples.

**Monte Carlo draws.** The limits of each effect are quantiles of its
`B` drawn values, so the Monte Carlo interval is run only when `B` puts
at least 10 draws beyond each limit (400 at the 95 percent level). A
draw in which an effect is undefined (a division by a draw near zero, or
a function outside its domain) is dropped and counted, never replaced,
and the interval is computed from the rest. The result records `B`
(requested), `B_used` (the number the interval rests on), `n_dropped`,
and `B_required` (the fewest replications a request may ask for; that
count is held to `B`, not to `B_used`) as attributes, and the line under
the table states them, with `bootstrap_status` (`"complete"`,
`"dropped"`, `"exploratory"`, `"insufficient_count"`, or
`"insufficient_share"`). A call gives at most one such warning, whatever
the number of effects. A warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the 2.5 percent of `B` that lies beyond one limit at the 95
percent level, because they could then decide by themselves where that
limit falls. When more than half were dropped, the limits of that effect
are `NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`. The draws are counted effect by effect:
an effect whose interval is withheld does not withhold the others,
`B_used` is the fewest draws any effect's limits rest on, and
`bootstrap_status` is the least favorable status among the effects. When
several effects are affected the one warning lists them (draws dropped,
draws used, and whether the interval is reported), gives in full the
account of the effect with the fewest usable draws, has the least
favorable class among them, and carries the list as the data frame
`effects`. `options(dmar.bootstrap_exploratory = TRUE)` lets a smaller
`B` (2 or more) run for a first look or a demonstration; every such call
gives a warning of class `dmar_bootstrap_exploratory` and the result is
not one to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

The Wald interval is the symmetric delta method interval and is reported
only for comparison.

**Multiple groups (moderated mediation across groups).** With `group`
(or list-form `S`, `M`, `N`), the model is fit as a multiple-group SEM:
the same structure in every group, with group-specific parameters. Every
enumerated effect is then estimated within each group (its term carries
the group label as a suffix), and for each effect the difference from
the reference group (the first group label) is estimated and tested with
the same constrained-optimization machinery, since a difference of two
products is itself a nonlinear function of the parameters. That
difference test is moderated mediation with a categorical moderator:
"does the indirect effect differ across groups?" is exactly the
between-group contrast of the conditional indirect effects. lavaan's
labeling rules decide what varies: a single label such as `b1` on a path
applies to every group and therefore imposes cross-group equality (the
corresponding difference is identically zero and its row is dropped); to
let a path differ by group, leave it unlabeled or give per-group labels
with the vector form `c("b1_f", "b1_m")*x`. With more than two groups,
each non-reference group is compared with the reference group. The call
below fits the simple model in both groups of a data frame `d` whose
`condition` column takes two values, and reports every effect per group
beside its between-group difference; leaving the paths unlabeled lets
them differ by group.

    mediation_mbco("m ~ x
                    y ~ m + x",
                   data = d, group = "condition", x = "x", y = "y",
                   ci_method = "wald")

Each reported row costs its own constrained null model fit, so this
two-group analysis fits nine null models where the single-group analysis
fits three.

**Moderated mediation, probed.** With `moderator`, every regression
coefficient along a pathway becomes a linear function of the moderator
wherever the model contains the matching interaction (stated as `x:w` in
the syntax, or as a product column among the predictors). A pathway
effect, the product of its edge coefficients, is then a polynomial in
the moderator, and the function derives that polynomial symbolically.
Three kinds of rows follow for each moderated effect. First, the
conditional effect at each probe value, each with its own confidence
interval and MBCO test. Second, the polynomial's moderator coefficients:
for a pathway moderated in one place the single such coefficient is
exactly the index of moderated mediation (Hayes, 2015), here tested by
likelihood ratio rather than bootstrap; a pathway moderated in several
places gets one row per power of the moderator. Third, for a pathway
moderated in several places (so the conditional effect is curved in the
moderator), a joint constancy test of all moderator coefficients at
once, with as many degrees of freedom as constraints, asking whether the
pathway effect depends on the moderator at all; as a joint test its row
reports no scalar estimate. Unmoderated effects in the same model keep
their single rows. The null set of a conditional-effect constraint is
again a union of branches (one edge's conditional coefficient or
another's must vanish at the probed value), and the null fits are
started on each branch, and at the model with all interactions removed
for the moderation tests, so the reported statistics reflect the
best-fitting null models. For example:

    mediation_mbco("m ~ x + w + x:w
                    y ~ m + x + w",
                   data = d, x = "x", y = "y", moderator = "w")

[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md)
draws the same conditional effects as curves over the moderator's whole
range with a confidence band, the visual companion to the probe rows.
Bespoke conditional quantities can still be written directly with `:=`
definitions or `hypotheses` when the built-in probing does not cover
them.

**Moderated mediation and mediated moderation.** Moderated mediation
asks whether an indirect effect depends on a moderator \\W\\: the
conditional indirect effect varies with \\w\\ (Muller, Judd, & Yzerbyt,
2005; Preacher, Rucker, & Hayes, 2007). Mediated moderation asks whether
an observed \\X \times W\\ interaction on \\Y\\ is transmitted through
the mediator. The two share their algebra: with first-stage moderation
\\M = a_1 X + a_2 W + a_3 XW + e_M\\ and \\Y = b M + \ldots\\, the
quantity \\a_3 b\\ is at once the index of moderated mediation (Hayes,
2015) and the indirect effect of the product term through \\M\\; what
differs is the question and the reporting emphasis. This function
reports the conditional-indirect-effect framing: the probe rows and the
moderation rows above. For a categorical moderator, use `group`.

**Refusals, warnings, and judgment calls.** The function stops where a
computed number would be mislabeled: a nonrecursive (feedback)
regression structure, categorical-endogenous syntax (thresholds),
explicit `==` constraints, a model with no indirect pathway between `x`
and `y`. It stops as well, with an error of class `dmar_no_maximum`,
when the lavaan fit that the OpenMx fits start from has no maximum, an
error variance running off below zero while the likelihood is flat or
still rising (see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)),
printing the lavaan code of that fit, or none for a fit given as
`model`. It warns where the data make trouble detectable: an endogenous
variable with only two distinct values (a binary mediator or outcome,
where the product of coefficients is not the causal indirect effect), an
interaction term among the predictors that no declared `moderator`
accounts for (declaring the moderator resolves the warning by probing
the moderation), an interaction included without its matching main
effect (the principle of marginality), null models that converge
imperfectly, profile bounds that fail. What it cannot check is left to
the analyst, and stated rather than assumed: the no omitted confounder
assumption, linearity, and the causal direction of the arrows come from
the design, not from the fit.

**Information criteria.** The columns `delta_aic` and `delta_bic` report
AIC and BIC for the null model minus the same criterion for the full
model; positive values favor the full model (the effect improves fit by
more than the parsimony penalty). A scalar equality constraint reduces
the effective number of free parameters by one, so the differences equal
\\\mathrm{LRT}\_{\mathrm{MBCO}} - 2\\ for the AIC and
\\\mathrm{LRT}\_{\mathrm{MBCO}} - \log N\\ for the BIC; for a joint test
the 2 and the \\\log N\\ are multiplied by its degrees of freedom. (The
OpenMx summary counts the same number of estimated parameters in both
models, so its printed AIC difference equals the likelihood ratio
statistic; DMAR counts the constraint against the null model, consistent
with the degrees of freedom of the test.)

**Constraints that add no restriction.** The degrees of freedom of a
test are the number of its constraints, which is the number of
restrictions only when each constraint restricts the parameters beyond
what the others do. Before any OpenMx model is fit, the constraints of
every test, written in the parameter labels, are checked twice. First,
the function forms the Jacobian of the constraints (the matrix of their
derivatives with respect to the free parameters) at the unrestricted
estimates of the lavaan fit, which lie off the null set, and reads its
rank after scaling every column by its largest entry and every row to
unit length, so that the rank does not depend on the units of the
parameters or the scale of the constraints. A singular value below
\\10^{-7}\\ times the largest counts as zero (the same number as the
default `tol` of [`qr`](https://rdrr.io/r/base/qr.html)). The rank is
read again at a fixed point near the estimates, so that an estimate of
exactly zero cannot make independent constraints look dependent, and the
larger rank is used. A test whose constraints have a rank below their
number (a constraint that is a function of the others) is refused (see
`hypotheses`); so is a single constraint whose gradient is zero, which
restricts nothing. The rank is not read on the null set: at the crossing
of the branches of a product constraint its gradient vanishes, and a
rank read there would call every set that contains it dependent.

Second, constraints whose Jacobian has full rank can still restrict the
parameters in fewer directions than their number, when one of them is
zero wherever the others are: `c("a*b", "b")` has the null set \\b =
0\\, `c("a", "a*cp")` the null set \\a = 0\\, and `c("a*b", "a*b*cp")`
the null set \\ab = 0\\, one restriction each. For every constraint of a
joint test the function finds points of the null set of the others, by
Gauss-Newton steps from fixed starting points that cover both signs and
a range of magnitudes of the parameters (in coordinates scaled by the
estimates, and again in coordinates scaled by the derivatives of the
constraints, so that neither the units of the parameters nor the
starting points of one scaling decide the result), and refuses the test
when the constraint is zero at all of them, to within \\10^{-6}\\ of its
typical size, under both scalings. No random numbers are drawn. A set
such as `c("a*b", "a*cp")`, whose null set is the union of \\a = 0\\
(one restriction) and \\b = cp = 0\\ (two), is accepted with 2 degrees
of freedom: the chi square reference on 2 degrees of freedom is right
where the second branch holds, and the test is conservative where only
\\a = 0\\ holds, since the statistic there is asymptotically chi square
on 1 degree of freedom, which exceeds the .95 quantile of the chi square
distribution on 2 degrees of freedom with probability .0144. The checks
do not detect a set none of whose constraints is zero wherever the
others are but whose null set still imposes fewer restrictions than
there are constraints, such as `c("a*b", "a*cp", "b*cp")` (any two of
the three coefficients zero, two restrictions); its test is
conservative. Neither check establishes that the null hypothesis is
regular (see *The reference distribution where the branches cross*
above).

**Change in explained variance.** Following the reporting recommendation
of Tofighi and Kelley (2020), the function computes \\R^2\\ for every
endogenous variable under the full model (the `"R2"` attribute) and the
drop in each \\R^2\\ under every null model (the `"delta_R2"` attribute,
variables by tested effects), so the fit cost of removing an effect can
be read as a change in effect size and not only as a test statistic.

**Summary statistics input.** With complete data the multivariate normal
log likelihood depends on the data only through the sample means and the
sample covariance matrix. When `S` (and optionally `M`) is supplied, the
function therefore constructs an internal data set with exactly those
moments (via [`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) with
`empirical = TRUE`) and proceeds as with raw data; the estimates,
likelihood ratio tests, and confidence intervals are identical to what
the raw data would give, whatever internal data set realizes the
moments. `S` is treated as the unbiased (divisor \\N - 1\\) covariance
matrix, the form reported in articles. This is how a published mediation
analysis can be reproduced, and its hypotheses re-tested with the MBCO
procedure, from a table of descriptive statistics alone.

**Assumptions.** The causal reading of any mediation analysis rests on
the no omitted confounder assumption for the predictor to mediator and
mediator to outcome relations, in addition to the usual distributional
assumptions (residuals multivariate normal, linear relations, no
treatment by mediator interaction); randomizing \\X\\ supports the first
link but not the second (MacKinnon, 2008; Tofighi & Kelley, 2020). With
a binary randomized \\X\\ the variable enters the model as numeric 0/1,
its exogenous variance freely estimated; the normality assumption
applies to the residuals of the endogenous variables.

## References

Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
analysis: Introducing the model-based constrained optimization
procedure. *Psychological Methods, 25*(4), 496–515.
[doi:10.1037/met0000259](https://doi.org/10.1037/met0000259)

Tofighi, D., & Kelley, K. (2020). Indirect effects in sequential
mediation models: Evaluating methods for hypothesis testing and
confidence interval formation. *Multivariate Behavioral Research,
55*(2), 188–210.
[doi:10.1080/00273171.2019.1618545](https://doi.org/10.1080/00273171.2019.1618545)

Hayes, A. F. (2015). An index and test of linear moderated mediation.
*Multivariate Behavioral Research, 50*(1), 1–22.
[doi:10.1080/00273171.2014.962683](https://doi.org/10.1080/00273171.2014.962683)

MacKinnon, D. P. (2008). *Introduction to statistical mediation
analysis*. Erlbaum.

MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
limits for the indirect effect: Distribution of the product and
resampling methods. *Multivariate Behavioral Research, 39*(1), 99–128.
[doi:10.1207/s15327906mbr3901_4](https://doi.org/10.1207/s15327906mbr3901_4)

MacKinnon, D. P., Valente, M. J., & Wurpts, I. C. (2018). Benchmark
validation of statistical models: Application to mediation analysis of
imagery and memory. *Psychological Methods, 23*(4), 654–671.
[doi:10.1037/met0000174](https://doi.org/10.1037/met0000174)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Muller, D., Judd, C. M., & Yzerbyt, V. Y. (2005). When moderation is
mediated and mediation is moderated. *Journal of Personality and Social
Psychology, 89*(6), 852–863.
[doi:10.1037/0022-3514.89.6.852](https://doi.org/10.1037/0022-3514.89.6.852)

Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick, T. R.,
Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes, H. H., & Boker,
S. M. (2016). OpenMx 2.0: Extended structural equation and statistical
modeling. *Psychometrika, 81*(2), 535–549.
[doi:10.1007/s11336-014-9435-8](https://doi.org/10.1007/s11336-014-9435-8)

Neale, M. C., & Miller, M. B. (1997). The use of likelihood-based
confidence intervals in genetic models. *Behavior Genetics, 27*(2),
113–120.
[doi:10.1023/A:1025681223921](https://doi.org/10.1023/A%3A1025681223921)

Preacher, K. J., Rucker, D. D., & Hayes, A. F. (2007). Addressing
moderated mediation hypotheses: Theory, methods, and prescriptions.
*Multivariate Behavioral Research, 42*(1), 185–227.
[doi:10.1080/00273170701341316](https://doi.org/10.1080/00273170701341316)

Wilks, S. S. (1938). The large-sample distribution of the likelihood
ratio for testing composite hypotheses. *Annals of Mathematical
Statistics, 9*(1), 60–62.
[doi:10.1214/aoms/1177732360](https://doi.org/10.1214/aoms/1177732360)

## See also

[`plot_mediation_mbco`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md)
for the conditional-effect display of a moderated analysis;
[`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md) for
the regression-based simple mediation model with bootstrap intervals;
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
and
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)
for planning;
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md)
for the delta method variance.

Other mediation:
[`mediate()`](https://yelleknek.github.io/DMAR/reference/mediate.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Replicate the memory experiment analyses of Tofighi and Kelley
# (2020) from the summary statistics in their Table 1 (data from
# MacKinnon, Valente, & Wurpts, 2018; N = 369). Instruction is 1 for
# imagery rehearsal instructions and 0 for repetition instructions.
# Requires the OpenMx and lavaan packages to be installed.
vars <- c("instruction", "imagery", "repetition", "recall")
sds  <- c(0.50, 2.96, 2.84, 3.40)
R_tk <- matrix(c(1.00,  .62, -.67,  .32,
                  .62, 1.00, -.56,  .51,
                 -.67, -.56, 1.00, -.28,
                  .32,  .51, -.28, 1.00), 4, 4,
               dimnames = list(vars, vars))
S_tk <- outer(sds, sds) * R_tk
M_tk <- c(instruction = 0.51, imagery = 5.66, repetition = 6.08,
          recall = 12.07)

# Single-mediator model: instruction -> imagery -> recall. Every
# reported effect costs its own constrained null model fit in OpenMx,
# which is where the run time goes; the Monte Carlo interval itself
# is inexpensive at any B. The indirect effect is about 2.1 words
# (the paper reports 2.121, SE = 0.276, 95% Monte Carlo CI
# [1.600, 2.682]), and its likelihood ratio statistic of 71.31 is the
# value discussed under Details.
single <- "
  imagery ~ b1*instruction
  recall  ~ b2*imagery + b3*instruction
"
mediation_mbco(single, S = S_tk, M = M_tk, N = 369,
               x = "instruction", y = "recall",
               ci_method = "monte_carlo", seed = 113)
#> Mediation tests via model-based constrained optimization (MBCO)
#>   Effects of instruction on recall; N = 369; optimizer: SLSQP
#>   Full model: deviance = 4040.806567, AIC = 4058.806567, BIC = 4094.003737
#>   R2 (full model): imagery = 0.384, recall = 0.260
#>   Each p_value is from the MBCO likelihood ratio test of the
#>   null model constraining that effect to zero. Causal readings
#>   rest on the no omitted confounder assumption.
#> 
#>  pathway                              term                 estimate se   
#>  instruction -> recall (all pathways) total_effect         2.176    0.335
#>  instruction -> recall                direct_effect        0.042    0.388
#>  instruction -> imagery -> recall     indirect_via_imagery 2.134    0.279
#>  ci_lower ci_upper lrt    df p_value  delta_aic delta_bic
#>  1.522    2.829    39.863 1  < 0.0001 37.863    33.952   
#>  -0.703   0.812    0.012  1  0.9139   -1.988    -5.899   
#>  1.594    2.701    71.306 1  < 0.0001 69.306    65.395   
#> 
#> Confidence level: 95%
#> 
#> Monte Carlo interval: 10000 replications, 10000 used (0 dropped; at least 400 must be requested)

# The parallel two-mediator model of the same paper, with the contrast
# of its two specific indirect effects, is described under Details;
# it is fit the same way from the same summary statistics, adding the
# 'hypotheses' argument. Raw data go in through 'data' rather than
# 'S', 'M', and 'N'; adding 'group' fits the model in every group and
# tests the between-group difference of each effect, which is
# moderated mediation with a categorical moderator (see Details). A
# continuous moderator goes in through 'moderator', and the help page
# for plot_mediation_mbco fits that analysis from a data frame and
# draws the conditional effects it estimates.
```
