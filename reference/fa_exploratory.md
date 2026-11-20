# Exploratory Common-Factor Analysis With Confidence Intervals

Fits a common-factor model, with rotation and uncertainty for its
rotated loadings. This is an exploratory measurement model: every item
can load on every factor. The estimates and rotation-aware standard
errors are computed by lavaan, which is required; DMAR adds the table
and its intervals, the display, the refusals, and the fit measure
policy. The function does not compute principal components. It carries a
name of its own so that attaching lavaan, whose `efa()` it calls, masks
nothing. `efa_k()` is the same function under a shorter name, as
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
is for
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).

## Usage

``` r
fa_exploratory(
  data = NULL,
  nfactors = 1L,
  S = NULL,
  N = NULL,
  items = NULL,
  estimator = "ML",
  missing = "fiml",
  ordered = NULL,
  rotation = "oblimin",
  target = NULL,
  conf_level = 0.95,
  rmsea_conf_level = 0.9,
  seed = NULL,
  point_estimates_only = FALSE
)

efa_k(
  data = NULL,
  nfactors = 1L,
  S = NULL,
  N = NULL,
  items = NULL,
  estimator = "ML",
  missing = "fiml",
  ordered = NULL,
  rotation = "oblimin",
  target = NULL,
  conf_level = 0.95,
  rmsea_conf_level = 0.9,
  seed = NULL,
  point_estimates_only = FALSE
)

# S3 method for class 'dmar_fa_exploratory'
print(x, digits = getOption("dmar.digits", 3L), ...)
```

## Arguments

- data:

  Numeric raw observations. Supply exactly one of `data` or `S`.

- nfactors:

  Number of common factors, a positive integer giving an identified
  exploratory model. With a `target`, its number of factors when
  `nfactors` is not supplied; a supplied value must agree with it.

- S:

  Named positive definite covariance matrix, with `N`. A correlation
  matrix is a covariance matrix with unit variances, and the original
  metric then coincides with the standardized one.

- N:

  Sample size for covariance input.

- items:

  The items to analyze, by name: columns of `data`, or rows and columns
  of `S`, which is then reduced to them, in the order given; only the
  reduced matrix need be finite and positive definite. `NULL`, the
  default, analyzes every column of `data` or of `S`, or, with a
  `target`, the items the target names. Names that cannot be read as at
  least three distinct items are refused with an error of class
  `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

- estimator:

  `"ML"` (default) or `"MLR"` for continuous items; `"ULS"` for the
  minimum residual (minres) solution of Harman and Jones (1966),
  unweighted least squares fitted to the correlation matrix (see
  Details); explicitly request `"WLSMV"` for ordered items. MLR uses the
  same likelihood estimates with Huber-White sandwich errors.

- missing:

  `"fiml"` (default), `"ml"`, or `"listwise"` for continuous items under
  ML or MLR. With ordered items the default fits the items as they stand
  when no item value is missing, as
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  does; with missing values ordered items require an explicit
  `"pairwise"` or `"listwise"` choice, and FIML is not substituted with
  deletion. ULS has no likelihood, so it takes complete raw data, an
  explicit `"listwise"` choice, or a correlation matrix estimated
  elsewhere through `S` (for example the maximum likelihood estimate of
  [`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md));
  with incomplete data the default `"fiml"` is refused rather than
  replaced. Covariance input has no missing-data route.

- ordered:

  `NULL`, `TRUE` for all items, or their names.

- rotation:

  The rotation criterion, by name. Oblique, the factors free to
  correlate: `"oblimin"` (default; direct oblimin with gamma 0, Jennrich
  & Sampson, 1966), `"quartimin"` (Carroll, 1953), `"geomin"` (Yates,
  1987; epsilon 0.001), `"simplimax"` (Kiers, 1994), `"oblimax"`
  (Saunders, 1961), `"infomax"` (McKeon, 1968, as Browne, 2001, presents
  it), and `"bentler"` (Bentler, 1977). Orthogonal, the factors kept
  uncorrelated: `"varimax"` (Kaiser, 1958; with Kaiser row weights, as
  lavaan applies it) and `"quartimax"` (Neuhaus & Wrigley, 1954).
  `"target"` and `"target_orthogonal"` rotate to the partially specified
  target given in `target`, the factors free to correlate (Browne,
  1972b) or kept uncorrelated (Browne, 1972a); see Details. `"none"`
  returns the unrotated solution lavaan identifies, whose factors are
  uncorrelated. Each criterion is applied in the form named here by
  lavaan's gradient projection algorithm from 30 random starts;
  Browne (2001) reviews the family.

- target:

  For `rotation = "target"` or `"target_orthogonal"`, the partially
  specified target, in either of two forms (`NULL`, the default, for
  every other rotation). The usual one is a named list of item keys, the
  form
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  takes for `factors`: each element names a factor and the items
  expected to load on it. An item keyed to a factor is targeted at 0 on
  every factor it is not keyed to and left unspecified on the factors it
  is keyed to (an item keyed to two factors is unspecified on both), and
  an analyzed item keyed to no factor is unspecified on every factor.
  The other is a numeric matrix of items by factors, with the items as
  row names in any order, `NA` for each unspecified loading, and
  elsewhere a target value for the standardized pattern loading,
  strictly between -1 and 1 (see Details); column names, if given, name
  the factors.

- conf_level:

  Confidence level of the intervals for the parameters, the
  communalities, and the determinacy and *H*, default 0.95.

- rmsea_conf_level:

  Confidence level of the RMSEA interval on the fit line and in the
  `fit_measures` attribute (with the gamma hat and adjusted gamma hat
  intervals there, which come from the same noncentrality limits),
  default 0.90, the conventional level, whose limits decide the tests of
  close and not close fit at .05 (see
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md));
  `rmsea_conf_level = 0.95` gives a 95 percent interval. It is separate
  from `conf_level`, which does not reach the RMSEA.

- seed:

  Optional seed for rotation starts; default `NULL` uses the current
  random stream. A supplied seed restores the caller's stream.

- point_estimates_only:

  `FALSE` (default) for the full analysis; `TRUE` for the same estimates
  without their standard errors, tests, and intervals, which are then
  not computed (see Details), for a program that needs only the
  estimates and wants them sooner.

- x:

  For [`print()`](https://rdrr.io/r/base/print.html), a table returned
  by `fa_exploratory()`.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html), the decimal places
  of the display (default `getOption("dmar.digits", 3)`; see
  [`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)).

- ...:

  For [`print()`](https://rdrr.io/r/base/print.html), passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).

## Value

A `dmar_tbl` data frame (classes `dmar_fa_exploratory` and
`dmar_factor_analysis`) with `term`, the parameter labels `lhs`, `op`,
and `rhs`, `estimate`, `se`, `statistic`, `p_value`, `ci_lower`,
`ci_upper`, `conf_level`, and `metric`: the parameters in the original
metric, the standardized solution, the communalities (op
`"communality"`), Hofmann's complexity (op `"complexity"`, with no
interval), and each factor's determinacy and coefficient *H* (ops
`"determinacy"` and `"H"`, the factor as both `lhs` and `rhs`).
Attributes retain `fit` (the lavaan object), `method` (sample size,
estimator, standard error and missing data methods, the rotation and
whether it is orthogonal, the fitted structure, and for a target
rotation the target and the number of loadings it specifies),
`variance_accounted`, `factor_quality`, eligible `fit_measures` (a
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
table, whose RMSEA interval is at `rmsea_conf_level`),
`rmsea_conf_level`, and `point_estimates_only` (`TRUE` or `FALSE`, as
called). `tidy()` preserves the parameter rows; `glance()` reports the
analysis method and sample size, and both carry `point_estimates_only`
as a column; [`print()`](https://rdrr.io/r/base/print.html) shows the
display described in Details.

## Details

The model is fitted to the covariance matrix of the items, with their
means under FIML, and the standardized solution is derived from it
(Cudeck, 1989, on why a correlation matrix is not analyzed as if it were
a covariance matrix); the table carries both metrics. Two routes are
correlation structures by construction: the minimum residual solution,
which is defined on the correlation matrix and is not invariant to the
scale of the items, and the WLSMV analysis of ordered items, which fits
their polychoric correlations and thresholds.

FIML fits the raw incomplete observations directly, under a normal
working likelihood and missingness at random given the included items.
It does not treat an estimated unrestricted covariance matrix as if it
were an ordinary complete-data sample covariance. Factor signs and order
are arbitrary. Wald intervals describe the selected rotated solution,
conditional on the requested factor count; they do not include
uncertainty from choosing that count or rotation.

**Display.** The printed result, and `knit_print()` and
[`as_kable`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
of it, show the solution rather than the parameter list: one row per
item with its standardized pattern loading on each factor and that
loading's interval, its communality with its interval, and its
complexity; above the table the sample, estimator, missing data method,
rotation and its form, and standard error method; beneath it the
variance each factor accounts for, the factor correlations with their
intervals (or the statement that the rotation fixed them at zero), each
factor's determinacy and coefficient *H* with their intervals, and the
fit. The table itself keeps every parameter in both metrics, which
`tidy()` returns, and a subset of it prints as an ordinary DMAR table.

Hofmann's (1978) complexity of an item is the square of the sum of its
squared pattern loadings divided by the sum of its loadings to the
fourth power: 1 for an item that loads on one factor alone, and the
number of factors for an item that loads equally on all of them. The
variance a factor accounts for is its diagonal element of Phi Lambda'
Lambda, with Lambda the standardized pattern loadings and Phi the factor
correlation matrix (the column sum of squared loadings when the factors
are uncorrelated); these sum to the sum of the communalities, and each
is also given as a proportion of the total variance of the items (their
number, in the standardized metric) and of the common variance. No
cumulative column is shown, since the order of the factors is arbitrary.
The `variance_accounted` attribute holds this table. The fit line reads
the `fit_measures` attribute: the chi square statistic with its degrees
of freedom and p-value, the RMSEA with its interval at
`rmsea_conf_level` (90 percent by default), which the line names, the
CFI, TLI, and SRMR, and the version lines that say which indices they
are (see
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md));
under ULS it gives the discrepancy statistic, the NFI, and the SRMR
(below). The measures are computed from the unrotated solution of the
same model: the chi square statistic and the standard indices do not
depend on the rotation, but lavaan's FIML-corrected indices of a rotated
object do, since its rotated parameterization counts more free
parameters than the model has, so the fit of a rotated solution is taken
from the unrotated fit, and every rotation reports the same fit.

**Rotation to a partially specified target.** A partially specified
target states what is expected of some loadings and nothing of the rest:
each loading is either given a target value or left unspecified. The
rotation minimizes the sum of squared differences between the specified
loadings and their target values, in the standardized metric, and asks
nothing of the unspecified ones (Browne, 2001, p. 124). Zeros are the
usual specification, since which loadings should be near zero is easier
to say than how large the others should be (Browne, 2001, p. 124); the
list form of `target` specifies exactly those, a zero for each item on
each factor it is not keyed to. Unlike a confirmatory model, which
forces the specified loadings to their values, the rotation brings them
only as close to their targets as the data allow, so a specified zero
whose rotated loading is far from zero signals a misspecified zero
(Browne, 2001, pp. 124–125), and its interval in the table shows by how
much. A specified value must lie strictly between -1 and 1, and a target
with a value outside that range is refused: the target is stated in the
standardized metric, and such a value is more likely a loading taken
from an item's original metric than an expected standardized loading. No
row weights are applied (Browne, 2001, p. 125, finds the Kaiser and
Cureton-Mulaik row weights unsuited to a partially specified target),
and the rotation is iterated to a gradient tolerance of 1e-8 rather than
lavaan's 1e-5, so that zeros the rotation attains exactly (below) are
attained to about that precision. The factors keep lavaan's names (f1,
f2, ...) in the order of the target's columns, and the display prints
the target beneath the table with the names it gave. Rotation does not
change the model, so the fit is that of every other rotation.

When every specified value is zero, the sign of each factor is
arbitrary, and each factor is reflected, as lavaan does by default, so
that its loadings have a positive sum. Since lavaan reflects all the
factors or none, no factor is reflected when any specified value is
nonzero. A factor whose column holds a nonzero value then takes its sign
from the target (a negative target keeps its sign), but a factor whose
specified values are all zero takes the sign of the random start from
which the rotation reached its minimum. A supplied `seed` fixes that
sign; with `seed = NULL` the loadings of such a factor, and its
correlations with the other factors, can change sign from one run to the
next. Either sign gives the same solution, so such a factor may be read
with its loadings and its correlations reversed together.

The target must identify the rotation, and lavaan does not check that it
does: an unidentified rotation comes back with finite standard errors
and no warning. DMAR checks the conditions of Algina (1980) before
reporting. In oblique rotation each factor's column needs at least m - 1
specified zeros, for m factors, and the standardized unrotated loadings
of the items specified in it must have rank m - 1 (Browne, 2001, p. 115;
Algina, 1980, Theorem 2). A column with a nonzero target value needs one
entry more, m entries with rank m (Theorem 1): with the factor variances
at one, m - 1 restrictions that include a nonzero value need not
identify the rotation, as Algina (1980, pp. 393–394) shows by example,
so DMAR asks for the larger count. In orthogonal rotation fewer suffice:
taken in increasing order of their specified entries, the columns need
at least 0, 1, ..., m - 1 zeros, and again a column with a nonzero value
one entry more (Theorem 3), and the rank of a column's rows is checked
against the same numbers. Two columns with the same specified entries
and values would leave their factors indistinguishable. A target that
fails is refused with an error of class `dmar_target_unidentified`, and
one that cannot be read as a target with class `dmar_target_invalid`;
both inherit `dmar_refusal` and carry a `reason` element (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).
Only a set of rows that is rank deficient to within the precision of the
fitted loadings is refused (a singular value below 1e-6 times the
largest); one that is nearly deficient shows itself in wide intervals.

When a column carries only the minimum number of zeros, the rotation can
attain them exactly. Those loadings are then fixed by the target rather
than estimated (lavaan reports a standard error of zero and a statistic
that means nothing). The table shows each at its target value itself, in
place of the rotated estimate, which differs from it only by the
rotation's convergence error, with `se`, `statistic`, and `p_value` `NA`
and both limits at the target, in both metrics. When every specified
loading is attained this way, as in an oblique rotation with m - 1 zeros
in every column, the rotated solution is the confirmatory factor model
with those loadings fixed at their targets and the factor variances at
one (the factors uncorrelated in an orthogonal rotation).

Myers, Jin, Ahn, Celimli, and Zopluoglu (2015) studied target rotation
by simulation. Three factors were measured by 10 or 20 items (p. 497).
The communalities were low (.20 to .40), high (.60 to .80), or wide (.20
to .80), and minor factors added error to the model (both p. 498). The
target values were in error by up to two standard errors (p. 500). The
errors in the target values had a negligible effect on the accuracy and
the stability of the rotated loadings (pp. 502–503). Specifying more
targets per factor improved both somewhat, the accuracy most when the
communalities were low and the stability most when they were high (p.
502).

In their comparison with geomin, every target value was in error by two
standard errors (pp. 500–501). Target rotation recovered the population
loadings more accurately in all 27 comparisons, one for each number of
targets at each level of communality, though by small amounts (p. 502).
Their measure of accuracy was the absolute value of the bias of a
replication, the signed error averaged over all the loadings (p. 501);
its mean was at or below 0.10 for both rotations, and the largest
difference between them was 0.08 (p. 502). Geomin gave the more stable
solution in 20 of the 27, among them every comparison with low
communalities; target rotation was the more stable with high
communalities and six or more targets per factor, and with wide
communalities and nine or more (pp. 499, 502–503). Their advice on the
number of targets is to specify the target as fully as the theory
allows, particularly when the items are not uniformly highly reliable
indicators of the factors (p. 503). Comparing the rotations, they
conclude that the extra work of specifying a target may be worthwhile
for accuracy; for stability, they expect geomin to do better when the
items are less than uniformly reliable indicators, and target rotation
to do better when the items are uniformly reliable and several targets
are specified (p. 504). Less stability is more sampling variability,
which the intervals of a given analysis show.

The intervals condition on a target fixed before the data were seen.
Browne (2001, pp. 125, 148) describes revising a target in sequence,
leaving unspecified the entries that proved misspecified and specifying
as zeros the loadings that proved near zero; such a sequence is
exploratory, capitalizes on chance, and calls for a confirmatory
analysis of new data. The `method` attribute holds the target (`target`,
items by factors, `NA` where unspecified, with the target's factor
names) and the number of loadings it specifies (`n_specified`).

**Determinacy and H.** The determinacy of a factor is the correlation of
its regression factor score with the factor, the multiple correlation of
the factor with the items (Grice, 2001; the factor determinacy index of
Ferrando & Lorenzo-Seva, 2018): \\\sqrt{\[\Phi \Lambda' R^{-1} \Lambda
\Phi\]\_{kk}}\\, with \\\Lambda\\ the standardized pattern loadings,
\\\Phi\\ the factor correlations, and \\R\\ the item correlation matrix
the model implies. Its square is the largest reliability any weighted
sum of the items reaches as a measure of the factor, which Ferrando and
Lorenzo-Seva call the generalized *H* and
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
the maximal reliability. Coefficient *H* is Hancock and Mueller's (2001)
closed form on the factor's standardized loadings \\l_j\\ alone, \\H =
(1 + (\sum_j l_j^2 / (1 - l_j^2))^{-1})^{-1}\\, applied to one factor of
a multidimensional solution as Rodriguez, Reise, and Haviland (2016)
apply it, so that each item's variance from the other factors counts as
error. With one factor *H* is the squared determinacy; otherwise the two
can differ in either direction, and the squared determinacy is the
quantity to read beside *H*. Both are reported for every factor with a
delta method standard error over lavaan's covariance matrix of the
standardized estimates (the matrix behind the standardized loading and
communality intervals) and Wald limits formed on the logit scale and
transformed back (Browne, 1982; the confidence interval transformation
principle of Kelley, 2007: the limits of a monotone transformation of a
parameter are the transformation of its limits), as
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
forms them by default, so they lie within (0, 1); the `se` column is the
standard error of the index itself, and with one factor it is the
standard error
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
gives its *H* row for the same model fitted the same way (by FIML). For
ordered items both describe the latent responses behind the categories,
not composites of the observed scores (the G-H-latent of Ferrando and
Lorenzo-Seva), and the display says so. Beside each determinacy \\\rho\\
the display and the `factor_quality` attribute give \\2\rho^2 - 1\\
(`min_score_correlation`, with its delta method standard error and the
determinacy's limits mapped through it), the lowest correlation two sets
of factor scores equally consistent with the model can have (Guttman,
1955, as Grice, 2001, p. 435, presents it): zero when \\\rho^2 = .5\\,
negative below. No line is drawn.
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md)
gives the determinacy of the Bartlett score as well, and for each
missing data pattern. The `factor_quality` attribute holds this table.

`estimator = "ULS"` is the minimum residual method: the loadings and
unique variances minimize the sum of squared residuals of the
correlation matrix, the unique variances absorbing the diagonal, so that
the off-diagonal correlations are reproduced as closely as the common
factors allow (Harman & Jones, 1966; `psych::fa(fm = "minres")` gives
the same loadings and communalities, up to the order and signs of the
factors). The estimates need no distributional assumption, and the
correlation matrix is the structure fitted (lavaan's
`correlation = TRUE`), so rescaling an item does not change them. With
raw data the standard errors are lavaan's sandwich for a least squares
estimator with the normal theory covariance matrix of the correlations
(its `se = "robust.sem.nt"`; Browne, 1984), and the rotation enters them
as for ML. A supplied correlation or covariance matrix carries no such
covariance matrix, so with `S` the ULS estimates and communalities are
returned with their standard errors and limits `NA`, and the note under
the table says so. The ULS statistic is in the units of the squared
correlations and has no chi square reference distribution, so in the
`fit_measures` attribute the `chi_square` row holds N - 1 times the ULS
discrepancy with no p-value, the RMSEA rows and the CFI, TLI, RNI, and
IFI are `NA`, and the NFI and the SRMR, which do not rest on that
distribution, are reported (see
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md));
a maximum likelihood fit of the same items gives the fit indices and
their inference.

The table contains original-unit parameters, standardized loadings and
factor correlations, and communalities with intervals. A communality is
one minus the standardized residual variance, so its interval reverses
the limits of that variance. Wald limits are not clipped to \[0, 1\].
Improper or unidentified solutions are refused, with an error of class
`dmar_efa_no_solution` whose `reason` element says whether the model did
not converge, converged to an improper solution (lavaan's
post-estimation check failed: a residual variance below zero, or a
covariance matrix that is not positive definite), or left the
rotation-aware standard errors unavailable;
[`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md)
reads it to report such a count as a row. When a residual variance is
below zero and the likelihood has no maximum at all (the likelihood flat
or still rising as the variance falls, as when too many factors are
extracted), the refusal says so instead (`reason` `"no maximum"`, and
the class `dmar_no_maximum` as well; see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
for the check) and prints the
[`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html) code of the
unrotated fit, which selects `items` from the data or covariance matrix
given, where `bounds = "pos.var"`, lavaan's own default, or fewer
factors are the choices. When the unrestricted FIML likelihood
diagnostic fails but the constrained model and its H0-based errors
remain eligible, parameter results are qualified and H1-dependent fit
indices are withheld. Absence of a diagnosed obstacle does not certify a
unique global maximum.

**Point estimates only.** With `point_estimates_only = TRUE` the table
has the rows and the estimates of the full analysis (the loadings,
factor correlations, and residual variances in both metrics, the
communalities, Hofmann's complexity, and each factor's determinacy and
*H*), and so do the `variance_accounted`, `factor_quality`, and
`fit_measures` attributes, but no standard error, test of a parameter,
or interval is computed: the `se`, `statistic`, `p_value`, `ci_lower`,
and `ci_upper` columns are `NA`, as are the standard errors and limits
of `factor_quality`. lavaan is asked for no standard errors (its
`se = "none"`), so neither the information matrix behind them (with its
sandwich under MLR, WLSMV, and ULS) nor the rotation-aware errors of the
rotated loadings are formed, nor the delta method errors of the
standardized solution, the communalities, the determinacy, and *H*; and
the rotated fit is made without a test statistic, which DMAR never
reads, since the fit measures come from the unrotated fit. That fit
keeps its test, so every fit measure is the one the full analysis
reports, the robust statistic and the FIML-corrected indices under MLR
and the scaled and shifted statistic under WLSMV included. `conf_level`
then sets nothing, since the only interval left is the RMSEA interval of
the fit line, at `rmsea_conf_level`, and the display states no other
level. A loading of a target rotation that is attained exactly is
recognized by lying within 1e-6 of its target, in place of its standard
error of zero, and shown at the target as in the full analysis; in rare
cases (a chance of the order of 1e-5 for each specified loading) an
estimated loading falls that close to its target by chance and is shown
at it, which the full analysis would estimate, within 1e-6 of the same
value. The input and target checks, and the refusals of a solution that
did not converge or is improper and of a likelihood without a maximum
(whose code reproduces the full fit), are made as before. The refusal of
a solution whose rotation-aware standard errors are unavailable is the
full analysis's guard against a solution that is not identified at its
estimates, whose estimates are then not unique; it is kept by inverting
the expected information of the unrotated fit, which is cheap, and, only
when that fails, the information its standard errors would use (the
observed information under FIML); a solution whose information cannot be
inverted is refused with the same class and reason. A failure confined
to the rotation's part of the rotation-aware errors is not detected, so
in that rare case the estimates are returned for a solution the full
analysis refuses. A FIML fit whose standard errors would need the
unavailable unrestricted (H1) information is not refused, since none are
computed; its fit measures are withheld as in the full analysis, with a
warning of the same class whose note speaks of the estimates alone. The
display says that standard errors and intervals were not computed, the
result carries the setting as its `point_estimates_only` attribute, and
the lavaan object in `fit` has no standard errors (and, when rotated, no
test). The time saved depends on the estimator and the data. In one set
of timings with three factors under oblimin, the nine classic Holzinger
and Swineford tests took about a third less time under ML and MLR (0.79
s against 0.51 s under ML) and half or more less under ULS and WLSMV
(0.70 s against 0.27 s for the tests cut into five ordered categories),
and the rotated solution of 18 items with 1,000 cases and 410 patterns
of missing values took 29 s against 17 s under MLR with FIML, and 18 s
against 9 s under ML. What remains is the estimation itself and the fit
measures, chiefly the robust statistic and the FIML-corrected indices of
the unrotated fit, which are kept; an unrotated analysis therefore saves
little or nothing (16.5 s against 16.7 s on those 18 items), and a
program that also fits the unrotated model (for its own fit measures,
say) saves on the rotated fit alone: the two together took 45 s against
34 s.

## References

Algina, J. (1980). A note on identification in the oblique and
orthogonal factor analysis models. *Psychometrika, 45*(3), 393–396.
[doi:10.1007/BF02293911](https://doi.org/10.1007/BF02293911)

Bentler, P. M. (1977). Factor simplicity index and transformations.
*Psychometrika, 42*(2), 277–295.
[doi:10.1007/BF02294054](https://doi.org/10.1007/BF02294054)

Browne, M. W. (1972a). Orthogonal rotation to a partially specified
target. *British Journal of Mathematical and Statistical Psychology,
25*(1), 115–120.
[doi:10.1111/j.2044-8317.1972.tb00482.x](https://doi.org/10.1111/j.2044-8317.1972.tb00482.x)

Browne, M. W. (1972b). Oblique rotation to a partially specified target.
*British Journal of Mathematical and Statistical Psychology, 25*(2),
207–212.
[doi:10.1111/j.2044-8317.1972.tb00492.x](https://doi.org/10.1111/j.2044-8317.1972.tb00492.x)

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
*Topics in applied multivariate analysis* (pp. 72–141). Cambridge, UK:
Cambridge University Press.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*(1), 62–83.
[doi:10.1111/j.2044-8317.1984.tb00789.x](https://doi.org/10.1111/j.2044-8317.1984.tb00789.x)

Browne, M. W. (2001). An overview of analytic rotation in exploratory
factor analysis. *Multivariate Behavioral Research, 36*(1), 111–150.
[doi:10.1207/S15327906MBR3601_05](https://doi.org/10.1207/S15327906MBR3601_05)

Carroll, J. B. (1953). An analytical solution for approximating simple
structure in factor analysis. *Psychometrika, 18*(1), 23–38.
[doi:10.1007/BF02289025](https://doi.org/10.1007/BF02289025)

Cudeck, R. (1989). Analysis of correlation matrices using covariance
structure models. *Psychological Bulletin, 105*(2), 317–327.
[doi:10.1037/0033-2909.105.2.317](https://doi.org/10.1037/0033-2909.105.2.317)

Ferrando, P. J., & Lorenzo-Seva, U. (2018). Assessing the quality and
appropriateness of factor solutions and factor score estimates in
exploratory item factor analysis. *Educational and Psychological
Measurement, 78*(5), 762–780.
[doi:10.1177/0013164417719308](https://doi.org/10.1177/0013164417719308)

Grice, J. W. (2001). Computing and evaluating factor scores.
*Psychological Methods, 6*(4), 430–450.
[doi:10.1037/1082-989X.6.4.430](https://doi.org/10.1037/1082-989X.6.4.430)

Guttman, L. (1955). The determinacy of factor score matrices with
implications for five other basic problems of common-factor theory.
*British Journal of Statistical Psychology, 8*(2), 65–81.
[doi:10.1111/j.2044-8317.1955.tb00321.x](https://doi.org/10.1111/j.2044-8317.1955.tb00321.x)

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Harman, H. H., & Jones, W. H. (1966). Factor analysis by minimizing
residuals (minres). *Psychometrika, 31*(3), 351–368.
[doi:10.1007/BF02289468](https://doi.org/10.1007/BF02289468)

Hofmann, R. J. (1978). Complexity and simplicity as objective indices
descriptive of factor solutions. *Multivariate Behavioral Research,
13*(2), 247–250.
[doi:10.1207/s15327906mbr1302_9](https://doi.org/10.1207/s15327906mbr1302_9)

Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
loadings. *Psychometrika, 38*(4), 593–604.
[doi:10.1007/BF02291497](https://doi.org/10.1007/BF02291497)

Jennrich, R. I., & Sampson, P. F. (1966). Rotation for simple loadings.
*Psychometrika, 31*(3), 313–323.
[doi:10.1007/BF02289465](https://doi.org/10.1007/BF02289465)

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
factor analysis. *Psychometrika, 23*(3), 187–200.
[doi:10.1007/BF02289233](https://doi.org/10.1007/BF02289233)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kiers, H. A. L. (1994). Simplimax: Oblique rotation to an optimal target
with simple structure. *Psychometrika, 59*(4), 567–579.
[doi:10.1007/BF02294392](https://doi.org/10.1007/BF02294392)

McKeon, J. J. (1968). *Rotation for maximum association between factors
and tests*. Unpublished manuscript, Biometric Laboratory, George
Washington University.

Myers, N. D., Jin, Y., Ahn, S., Celimli, S., & Zopluoglu, C. (2015).
Rotation to a partially specified target matrix in exploratory factor
analysis in practice. *Behavior Research Methods, 47*(2), 494–505.
[doi:10.3758/s13428-014-0486-7](https://doi.org/10.3758/s13428-014-0486-7)

Neuhaus, J. O., & Wrigley, C. (1954). The quartimax method: An analytic
approach to orthogonal simple structure. *British Journal of Statistical
Psychology, 7*(2), 81–91.
[doi:10.1111/j.2044-8317.1954.tb00147.x](https://doi.org/10.1111/j.2044-8317.1954.tb00147.x)

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
bifactor models: Calculating and interpreting statistical indices.
*Psychological Methods, 21*(2), 137–150.
[doi:10.1037/met0000045](https://doi.org/10.1037/met0000045)

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

Saunders, D. R. (1961). The rationale for an "oblimax" method of
transformation in factor analysis. *Psychometrika, 26*(3), 317–324.
[doi:10.1007/BF02289800](https://doi.org/10.1007/BF02289800)

Yates, A. (1987). *Multivariate exploratory data analysis: A perspective
on exploratory factor analysis*. State University of New York Press.

## See also

[`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
for the factor count,
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
for the confirmatory model,
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
which takes the lavaan object in the `fit` attribute, and
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
for coefficient *H* of one scale with bootstrap intervals.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
           "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
two <- fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
                      seed = 113)
two                       # the solution: loadings with intervals,
#> Exploratory factor analysis: 2 factors, 6 items, N = 301
#> Estimator: ML; missing data: FIML; rotation: oblimin (oblique); SE: standard
#> 
#>                        item                      f1                      f2
#>        t1_visual_perception    0.546 [0.391, 0.701]    0.229 [0.088, 0.370]
#>                    t2_cubes    0.443 [0.304, 0.582]   0.052 [-0.081, 0.185]
#>                 t4_lozenges    0.766 [0.629, 0.903] -0.072 [-0.125, -0.019]
#>  t6_paragraph_comprehension   0.040 [-0.032, 0.112]    0.833 [0.779, 0.888]
#>                 t7_sentence -0.088 [-0.139, -0.038]    0.899 [0.850, 0.949]
#>             t9_word_meaning    0.085 [0.009, 0.161]    0.805 [0.747, 0.863]
#>           communality complexity
#>  0.426 [0.286, 0.566]      1.341
#>  0.213 [0.103, 0.323]      1.028
#>  0.559 [0.333, 0.784]      1.018
#>  0.716 [0.641, 0.791]      1.005
#>  0.768 [0.690, 0.846]      1.019
#>  0.696 [0.621, 0.771]      1.022
#> (standardized pattern loadings and communalities with 95% intervals, and Hofmann's complexity; tidy(x) has every parameter in both metrics)
#> 
#> Variance accounted for (sums of squared standardized loadings, the diagonal of Phi Lambda' Lambda):
#>  factor sum_sq_loadings proportion_variance proportion_common
#>      f1           1.133               0.189             0.335
#>      f2           2.246               0.374             0.665
#>   total           3.379               0.563             1.000
#> 
#> Factor correlations with 95% intervals:
#>                      f1                   f2
#> f1                    1 0.302 [0.145, 0.459]
#> f2 0.302 [0.145, 0.459]                    1
#> 
#> Factor determinacy and coefficient H with 95% intervals:
#>  factor          determinacy                min_r                    H
#>      f1 0.833 [0.775, 0.878] 0.386 [0.201, 0.541] 0.678 [0.563, 0.775]
#>      f2 0.944 [0.932, 0.955] 0.784 [0.736, 0.824] 0.894 [0.865, 0.917]
#> (determinacy: the correlation of the regression factor score with its factor, whose square is the factor's maximal reliability; min_r: the lowest correlation two sets of scores equally consistent with the model can have, 2 x determinacy^2 - 1 (Guttman, 1955); H: Hancock and Mueller's reliability of the optimally weighted composite from the factor's loadings alone; limits from the logit scale)
#> 
#> Fit: chi square 2.109, df 4, p = 0.7157; RMSEA with 90% interval 0.000 [0.000, 0.064]; CFI 1.000; TLI 1.011; SRMR 0.008
#> Version: robust (auto: the fit uses FIML (missing = "ml"))
#> Chi square: standard; robust measures: FIML corrected indices of Zhang and Savalei (2023), version 6 (normal theory)
#> The fit_measures attribute holds every measure with its version lines.
#> 
#> Confidence level: 95%
#> 
#> Intervals condition on the selected factor count and rotation.
                          # communalities, complexity, variance, fit
attr(two, "variance_accounted")
#>   factor sum_sq_loadings proportion_variance proportion_common
#> 1     f1        1.132905           0.1888174          0.335269
#> 2     f2        2.246187           0.3743645          0.664731
#> 3  total        3.379092           0.5631820          1.000000
attr(two, "factor_quality")  # determinacy and H with intervals
#>   factor determinacy determinacy_se determinacy_ci_lower determinacy_ci_upper
#> 1     f1   0.8325100    0.026150086            0.7748582            0.8777286
#> 2     f2   0.9444978    0.005930719            0.9316611            0.9550396
#>   min_score_correlation min_score_correlation_se min_score_correlation_ci_lower
#> 1             0.3861459               0.08708084                      0.2008103
#> 2             0.7841521               0.02240620                      0.7359848
#>   min_score_correlation_ci_upper         H       H_se H_ci_lower H_ci_upper
#> 1                      0.5408150 0.6781128 0.05463200  0.5632974  0.7748093
#> 2                      0.8242012 0.8936313 0.01334572  0.8645025  0.9170986
generics::tidy(two)[1:6, ]   # every parameter, both metrics
#>                               term lhs op                        rhs   estimate
#> 1       f1 =~ t1_visual_perception  f1 =~       t1_visual_perception  3.8183915
#> 2                   f1 =~ t2_cubes  f1 =~                   t2_cubes  2.0842072
#> 3                f1 =~ t4_lozenges  f1 =~                t4_lozenges  6.9199113
#> 4 f1 =~ t6_paragraph_comprehension  f1 =~ t6_paragraph_comprehension  0.1396080
#> 5                f1 =~ t7_sentence  f1 =~                t7_sentence -0.4550080
#> 6            f1 =~ t9_word_meaning  f1 =~            t9_word_meaning  0.6478122
#>          se statistic      p_value    ci_lower   ci_upper conf_level   metric
#> 1 0.5977834  6.387583 1.685276e-10  2.64675750  4.9900255       0.95 original
#> 2 0.3607920  5.776755 7.615491e-09  1.37706785  2.7913466       0.95 original
#> 3 0.7401286  9.349606 0.000000e+00  5.46928586  8.3705367       0.95 original
#> 4 0.1286618  1.085078 2.778872e-01 -0.11256440  0.3917805       0.95 original
#> 5 0.1322710 -3.439967 5.817852e-04 -0.71425444 -0.1957615       0.95 original
#> 6 0.2968840  2.182038 2.910673e-02  0.06593022  1.2296942       0.95 original
#>   estimator se_method missing h1_diagnostic point_estimates_only
#> 1        ML  standard    fiml          <NA>                FALSE
#> 2        ML  standard    fiml          <NA>                FALSE
#> 3        ML  standard    fiml          <NA>                FALSE
#> 4        ML  standard    fiml          <NA>                FALSE
#> 5        ML  standard    fiml          <NA>                FALSE
#> 6        ML  standard    fiml          <NA>                FALSE

# An orthogonal rotation, its estimates alone: the standard errors,
# tests, and intervals are not computed, which saves time.
fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
               rotation = "varimax", seed = 113, point_estimates_only = TRUE)
#> Exploratory factor analysis: 2 factors, 6 items, N = 301
#> Estimator: ML; missing data: FIML; rotation: varimax (orthogonal), Kaiser row weights; SE: not computed (point_estimates_only = TRUE)
#> 
#>                        item    f1    f2 communality complexity
#>        t1_visual_perception 0.584 0.291       0.426      1.467
#>                    t2_cubes 0.450 0.105       0.213      1.108
#>                 t4_lozenges 0.747 0.022       0.559      1.002
#>  t6_paragraph_comprehension 0.194 0.824       0.716      1.111
#>                 t7_sentence 0.079 0.873       0.768      1.016
#>             t9_word_meaning 0.233 0.801       0.696      1.168
#> (standardized pattern loadings, communalities, and Hofmann's complexity; no intervals, see the note; tidy(x) has every parameter in both metrics)
#> 
#> Variance accounted for (sums of squared standardized loadings, the diagonal of Phi Lambda' Lambda):
#>  factor sum_sq_loadings proportion_variance proportion_common
#>      f1           1.200               0.200             0.355
#>      f2           2.179               0.363             0.645
#>   total           3.379               0.563             1.000
#> 
#> Factor correlations: fixed at zero (orthogonal rotation).
#> 
#> Factor determinacy and coefficient H:
#>  factor determinacy min_r     H
#>      f1       0.821 0.348 0.681
#>      f2       0.935 0.748 0.878
#> (determinacy: the correlation of the regression factor score with its factor, whose square is the factor's maximal reliability; min_r: the lowest correlation two sets of scores equally consistent with the model can have, 2 x determinacy^2 - 1 (Guttman, 1955); H: Hancock and Mueller's reliability of the optimally weighted composite from the factor's loadings alone)
#> 
#> Fit: chi square 2.109, df 4, p = 0.7157; RMSEA with 90% interval 0.000 [0.000, 0.064]; CFI 1.000; TLI 1.011; SRMR 0.008
#> Version: robust (auto: the fit uses FIML (missing = "ml"))
#> Chi square: standard; robust measures: FIML corrected indices of Zhang and Savalei (2023), version 6 (normal theory)
#> The fit_measures attribute holds every measure with its version lines.
#> 
#> Standard errors, tests, and intervals were not computed (point_estimates_only = TRUE); the estimates, the determinacy and H, and the fit measures are those of the full analysis.

# Rotation to a partially specified target, through efa_k(), the same
# function under a shorter name, here of the minimum residual (minres)
# solution: unweighted least squares on the correlation matrix, with
# sandwich standard errors from the raw data. Each item is targeted at
# zero on the factor it is not keyed to, and the items and the number of
# factors are read from the target.
efa_k(holzinger_swineford, rotation = "target", estimator = "ULS",
      seed = 113, target = list(spatial = items[1:3], verbal = items[4:6]))
#> Exploratory factor analysis: 2 factors, 6 items, N = 301
#> Estimator: ULS (minimum residual: unweighted least squares on the correlation matrix); missing data: listwise deletion; rotation: target (oblique), 6 of 12 loadings specified (all 0); SE: sandwich with the normal theory covariance matrix of the correlations (robust.sem.nt)
#> 
#>                        item                      f1                      f2
#>        t1_visual_perception    0.562 [0.393, 0.732]    0.185 [0.066, 0.304]
#>                    t2_cubes    0.459 [0.306, 0.613]   0.014 [-0.118, 0.145]
#>                 t4_lozenges    0.783 [0.624, 0.942] -0.132 [-0.206, -0.059]
#>  t6_paragraph_comprehension   0.038 [-0.036, 0.112]    0.834 [0.774, 0.894]
#>                 t7_sentence -0.089 [-0.153, -0.025]    0.905 [0.853, 0.957]
#>             t9_word_meaning    0.086 [0.012, 0.160]    0.796 [0.739, 0.854]
#>           communality complexity
#>  0.428 [0.260, 0.596]      1.213
#>  0.216 [0.100, 0.332]      1.002
#>  0.554 [0.325, 0.782]      1.057
#>  0.721 [0.637, 0.805]      1.004
#>  0.766 [0.691, 0.842]      1.019
#>  0.693 [0.617, 0.769]      1.023
#> (standardized pattern loadings and communalities with 95% intervals, and Hofmann's complexity; tidy(x) has every parameter in both metrics)
#> 
#> Target (the value where a loading is specified, ? where it is not):
#>                            f1 f2
#> t1_visual_perception        ?  0
#> t2_cubes                    ?  0
#> t4_lozenges                 ?  0
#> t6_paragraph_comprehension  0  ?
#> t7_sentence                 0  ?
#> t9_word_meaning             0  ?
#> Factors: f1 = spatial, f2 = verbal (the table keeps lavaan's names).
#> 
#> Variance accounted for (sums of squared standardized loadings, the diagonal of Phi Lambda' Lambda):
#>  factor sum_sq_loadings proportion_variance proportion_common
#>      f1           1.167               0.195             0.346
#>      f2           2.210               0.368             0.654
#>   total           3.378               0.563             1.000
#> 
#> Factor correlations with 95% intervals:
#>                      f1                   f2
#> f1                    1 0.374 [0.239, 0.509]
#> f2 0.374 [0.239, 0.509]                    1
#> 
#> Factor determinacy and coefficient H with 95% intervals:
#>  factor          determinacy                min_r                    H
#>      f1 0.840 [0.786, 0.883] 0.411 [0.234, 0.558] 0.700 [0.567, 0.806]
#>      f2 0.944 [0.933, 0.954] 0.784 [0.740, 0.821] 0.896 [0.865, 0.920]
#> (determinacy: the correlation of the regression factor score with its factor, whose square is the factor's maximal reliability; min_r: the lowest correlation two sets of scores equally consistent with the model can have, 2 x determinacy^2 - 1 (Guttman, 1955); H: Hancock and Mueller's reliability of the optimally weighted composite from the factor's loadings alone; limits from the logit scale)
#> 
#> Fit: ULS discrepancy statistic 0.373 on 4 df (no reference distribution, so no p-value, RMSEA, CFI, or TLI); NFI 0.999; SRMR 0.008
#> Version: standard (auto: a ULS fit of complete data with the standard test alone, whose statistic lavaan gives no reference distribution and which is the ULS discrepancy of a correlation structure)
#> Chi square: standard
#> The fit_measures attribute holds every measure with its version lines.
#> 
#> Confidence level: 95%
#> 
#> The minimum residual solution is unweighted least squares fitted to the correlation matrix. The ULS statistic has no chi square reference distribution: in the fit_measures attribute the chi_square row is N - 1 times the ULS discrepancy without a p-value, the RMSEA, CFI, TLI, RNI, and IFI are NA, and the NFI and SRMR are reported. Intervals condition on the selected factor count and on a target fixed before the data were seen; a target revised after seeing a rotation makes the analysis exploratory.
```
