# Measurement Invariance Across Groups

Fits the standard ladder of multiple group invariance models for a
measurement model of any number of factors and reports the comparison of
its rungs in one table: configural invariance (same pattern, all
parameters free by group), metric (equal loadings; required before
comparing relations involving the factor), scalar (equal loadings and
intercepts; required before comparing factor means), and strict (equal
residual variances as well; required before comparing observed-score
variances). With ordered indicators a thresholds rung comes first,
because thresholds rather than intercepts carry the location information
(Wu & Estabrook, 2016); with dichotomous indicators the two are not
separately identified, so that rung is folded into metric. Each rung is
tested against the previous with a likelihood ratio test, scaled when
the estimator in force carries a robust test, which is what declaring
ordered items arranges, and the practical-fit changes (delta CFI, delta
RMSEA) are reported alongside, since with large samples the chi square
will flag trivial differences (Cheung & Rensvold, 2002, suggest delta
CFI of about -.01 as a red flag). The fitted lavaan objects come back as
an attribute, so score tests and partial invariance refits do not
require refitting the ladder. Requires lavaan.

## Usage

``` r
measurement_invariance(
  data,
  model = NULL,
  group,
  items = NULL,
  levels = NULL,
  ordered = NULL,
  estimator = "ML",
  missing = "fiml",
  group_partial = NULL,
  parameterization = c("delta", "theta"),
  ...
)
```

## Arguments

- data:

  A `data.frame` with the items and the grouping variable.

- model:

  The measurement model, given either as lavaan model syntax (a single
  string, or a character vector of lines that is collapsed with
  newlines, such as
  `c("visual =~ x1 + x2 + x3", "verbal =~ x4 + x5 + x6")`) or as a named
  list mapping each factor name to a character vector of its item names,
  from which the syntax is built. Exactly one of `model` and `items`
  must be supplied. An indicator of either form that `data` lacks is
  refused with an error of class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

- group:

  Single character string naming the grouping column (two or more
  groups).

- items:

  Character vector (three or more) naming the indicator columns of a
  single factor. A convenience for the one-factor (congeneric) case of
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
  equivalent to `model = "f =~ item_1 + item_2 + ..."`. Exactly one of
  `model` and `items` must be supplied. Items that cannot be read from
  `data` (a name that is not a column, a name given twice, fewer than
  three) are refused with an error of class `dmar_items_invalid`.

- levels:

  Which rungs of the ladder to fit, in order: any leading subset of the
  ladder in force. With continuous indicators the ladder is
  `c("configural", "metric", "scalar", "strict")`; with ordered
  indicators it is
  `c("configural", "thresholds", "metric", "scalar", "strict")`; when
  every indicator is ordered and dichotomous it is
  `c("configural", "metric", "scalar")`, because neither the thresholds
  nor the strict rung is identified for a two-category item (see
  `ordered` below). `NULL` (the default) fits the whole ladder in force.

- ordered:

  Ordered categorical (including binary) items: `TRUE` for every
  indicator in the model, or a character vector naming the ordered ones.
  `NULL` (default) treats all indicators as continuous. Declaring
  ordered items switches the ladder and, unless the estimator already
  carries the mean and variance adjusted test, switches the estimator
  (see `estimator` below). The number of observed categories per
  declared item is counted before the ladder is built, because a
  dichotomous item has a single threshold that is not separately
  identified from the intercept of its underlying response (Millsap &
  Yun-Tein, 2004; Wu & Estabrook, 2016), and its residual variance is
  not a free parameter under either parameterization. When every
  declared item is dichotomous the thresholds rung is folded into
  metric, and the strict rung is dropped as well unless a continuous
  indicator is left to carry a residual variance; when only some
  declared items are dichotomous, both rungs are kept and a message
  names those items, whose thresholds and residual variances the rungs
  leave untested.

- estimator:

  Estimator passed to lavaan. Defaults to `"ML"`; use `"MLR"` for robust
  corrections, in which case the likelihood ratio tests use the scaled
  difference. When `ordered` is supplied, a maximum likelihood or
  generalized least squares (the `"GLS"` label; Browne, 1974) estimator
  is replaced by `"WLSMV"`, and `"DWLS"` and `"ULS"` are replaced by
  `"WLSMV"` and `"ULSMV"`. That second substitution changes only the
  test: in lavaan `"WLSMV"` is `"DWLS"` and `"ULSMV"` is `"ULS"` with
  the mean and variance adjusted statistic, which is what makes the
  rung-to-rung difference test asymptotically valid for ordered data.
  The one categorical estimator left alone is `"WLS"`, the full weight
  matrix estimator of Browne (1984) and Muthén (1984), whose statistic
  is asymptotically chi square already and whose weight matrix a
  substitution would silently change.

- missing:

  Missing data handling passed to lavaan: one of `"fiml"` (the default)
  or `"ml"` (full information maximum likelihood, which fits every case
  with an observed item), `"listwise"` (listwise deletion),
  `"pairwise"`, `"two_stage"`, or `"robust_two_stage"` (two-stage
  maximum likelihood, with normal theory or nonnormality robust standard
  errors and tests; lavaan's `"two.stage"` and `"robust.two.stage"` are
  accepted for these two). Full information maximum likelihood and
  two-stage estimation are not available with the categorical estimator,
  and two-stage estimation needs `estimator = "ML"`; asking for another
  combination is an error. The default analysis is full information
  maximum likelihood whenever a value of a model variable is missing and
  the complete data analysis otherwise, unlike
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  which fits complete data by FIML too: data in which no such value is
  missing, or none but in cases with every model variable missing (which
  both treatments set aside), are fitted as they stand (under FIML the
  robust indices of an MLR ladder would be the corrected ones of Zhang
  and Savalei, 2023, not the complete data ones), and incomplete data by
  FIML. FIML needs continuous items, `estimator = "ML"` or `"MLR"`,
  standard errors other than `se = "robust.sem"`, and a test other than
  the Satorra-Bentler kind passed through `...`; when the default meets
  incomplete items and one of these, the call stops and names the option
  to choose (`missing = "listwise"`, `"pairwise"` for ordered items,
  `estimator = "MLR"`, or, with `estimator = "MLR"`, leaving out `test`)
  rather than deleting cases unasked, and an explicit `missing = "fiml"`
  with `se = "robust.sem"` or such a test stops too, whether or not a
  value is missing, since lavaan would fit the complete cases instead.
  Under FIML a group that never observes a variable of the model, or
  never observes two of them together, is refused before any fit (an
  error of class `dmar_group_coverage`). See *Details*.

- group_partial:

  Character vector of parameters to leave free across groups at every
  rung, passed to lavaan's `group.partial`. Either user-supplied
  parameter labels or lavaan parameter specifications such as
  `"visual =~ x2"` or `"x2 ~1"`. This is the partial invariance case of
  Byrne, Shavelson, and Muthén (1989). `NULL` (default) constrains
  everything the rung calls for, as does an empty vector (`character(0)`
  or lavaan's own default `""`), which frees nothing in lavaan either.

- parameterization:

  Identification of ordered indicators, passed to lavaan: `"delta"` (the
  default, lavaan's) or `"theta"`. Residual variances are free
  parameters only under `"theta"`, so a strict rung requested with
  ordered items switches to `"theta"` with a message. Ignored when no
  item is ordered.

- ...:

  Further arguments passed to
  [`cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html) at every rung (for
  example `std.lv`, `orthogonal`, `cluster`). lavaan's own
  `group.partial` is merged into `group_partial`, and an empty one
  (`NULL`, `character(0)`, or `""`) frees nothing, as in lavaan, so a
  call such as `group.partial = if (free) terms else NULL` runs either
  way; `group.equal`, which the ladder sets rung by rung, is an error
  here (choose the rungs with `levels`).

## Value

A tidy wide `data.frame` (class `dmar_tbl`) with one row per fitted
level and columns `level` (label), `chi_square`, `df`, `p_chi_square`
(exact-fit test), `cfi`, `rmsea`, and, from the second row on, the
step-comparison columns `delta_chi_square`, `delta_df`, `p_value` (the
likelihood ratio test against the previous rung), `delta_cfi`, and
`delta_rmsea`. `delta_chi_square` and `p_value` are `NA` wherever
`delta_df` is zero, since a difference test on zero degrees of freedom
tests nothing. Attributes carry the non-numeric information: `"fits"`,
the named list of fitted lavaan objects, one per level; `"estimator"`,
the estimator actually used; `"ordered"`, `TRUE` when any indicator was
declared ordered; `"test"`, naming the chi square difference test used
(`NA` when a single rung was fit and nothing was compared, and for a ULS
ladder of continuous items without a corrected statistic, which has no
difference test); `"fit_indices"`, the version of the `cfi` and `rmsea`
columns at every rung, `"standard"`, `"robust"`, or `"scaled"`
(`"robust"` is the corrected indices of Zhang and Savalei, 2023, under
FIML and two-stage estimation, and `"scaled"` the fallback of a ladder
of ordered items, or fit by FIML or in two stages, with a scaled
statistic whose robust indices lavaan did not compute at every rung);
`"fit_columns"`, the lines printed under the table that name that
version for a ladder fit by FIML, one that fell back from the robust
indices, or one whose robust indices DMAR computed, and that say why
columns of a ULS ladder of continuous items are `NA`, or why the `rmsea`
of a rung with zero degrees of freedom that does not fit exactly is
(`NULL` otherwise; a two-stage ladder names the version in its own
lines); `"model_test"`, naming the statistic in the `chi_square` column
at every rung (`"standard"`, `"scaled"`, or under the two-stage options
the two-stage residual based or scaled statistic); for a two-stage
ladder, `"two_stage"`, a list with the method, the lines printed under
the table, and, in `rows`, the cases the first stage read in each group
and lavaan's count of them when an auxiliary variable was named, how
many were observed on an auxiliary variable alone, whether the rungs
were computed again with all of them (`corrected`), and whether the
estimates of a rung that constrains across groups were too
(`reestimated`); and `"model"`, the lavaan syntax that was fitted.

## Details

The models are nested by construction, each adding equality constraints
across groups to the previous. With continuous indicators the constraint
sets are none (beyond the configuration), then
`group.equal = "loadings"`, then `c("loadings", "intercepts")`, then
`c("loadings", "intercepts", "residuals")`.

The ordered ladder differs, and the difference is substantive rather
than cosmetic. For an ordered indicator the observed response is a
coarsening of an underlying continuous response at a set of thresholds,
and it is the thresholds, not an intercept, that locate the item on the
latent scale. Constraining loadings while leaving thresholds free across
groups therefore does not deliver what metric invariance is supposed to
deliver, and the accepted sequence (Millsap & Yun-Tein, 2004; Wu &
Estabrook, 2016) constrains thresholds first: `"thresholds"`, then
`c("thresholds", "loadings")` for metric, then
`c("thresholds", "loadings", "intercepts")` for scalar, then adding
`"residuals"` for strict. The intercept rung is not vacuous even when
every indicator is ordered: once thresholds and loadings are
constrained, lavaan frees the underlying-response intercepts in the
non-reference groups, and the scalar rung is what returns them to zero
and lets the latent means be estimated instead. Residual variances,
however, are free parameters only under the theta parameterization, so
`parameterization` switches to `"theta"` when a strict rung is requested
with ordered items.

Dichotomous items are the exception the ordered ladder has to make room
for. A two-category item contributes one threshold, and that threshold,
the intercept of the underlying response, and its residual variance are
not separately identified: the data give one proportion per group per
item, which pins down a single standardized location and nothing else
(Millsap & Yun-Tein, 2004; Wu & Estabrook, 2016). Constraining
thresholds alone across groups is then a reparameterization rather than
a restriction, since lavaan frees the underlying-response intercepts by
exactly as many parameters as the constraint removes, and the residual
variances stay fixed at one in every group whichever parameterization is
in force. So the function counts the observed categories of each
declared ordered indicator before building the ladder, and a message
explains what the count implies. When all of them are dichotomous the
thresholds rung is folded into metric, which constrains thresholds and
loadings together; the strict rung goes too when every indicator in the
model was declared ordered, leaving configural, metric, scalar. A
continuous indicator alongside dichotomous ones keeps the strict rung,
since its residual variance is a free parameter. When only some declared
items are dichotomous the full ordered ladder is kept, since the
polytomous items still carry testable thresholds and residual variances,
and the message names the dichotomous items so their contribution to
those two rungs is not overread. Every rung that is dropped is one that
would have cost zero degrees of freedom.

When the estimator carries a robust test, the difference between two chi
square statistics is not itself chi square distributed.
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html) then
returns the scaled difference test (the Satorra-Bentler or Satorra
correction, chosen by lavaan to match the test in force), and that is
what `delta_chi_square` and `p_value` report; the `"test"` attribute
names the test used. In that case the model-level `chi_square`,
`p_chi_square`, `cfi`, and `rmsea` columns are lavaan's scaled and
robust versions, so `delta_chi_square` will not equal the difference of
consecutive `chi_square` entries. That is a property of scaled
difference testing, not an inconsistency.

One version of the test and one of the indices hold every rung, so that
`delta_cfi` and `delta_rmsea` compare like with like, and the
`"model_test"` and `"fit_indices"` attributes name the versions the
columns hold. Where lavaan does not compute the robust indices at every
rung, every rung falls back to one other version, and a line under the
table says which. For a maximum likelihood fit of continuous items,
whose robust indices estimate the population values of the ML indices
with their small sample bias corrected for nonnormality, it is the
standard indices, which estimate the same population values with the
bias corrected as if the data were normal; the scaled indices, which put
the scaled statistic into the standard formulas, estimate population
values that change with the kurtosis of the data (Savalei, 2018). Beside
the older mean and variance adjusted statistic
(`test = "mean.var.adjusted"`, passed through `...`) lavaan computes no
robust indices, and they are computed here, as
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
computes them: they do not depend on which corrected statistic
accompanies them (Savalei, 2018), so they equal those of the
Satorra-Bentler statistic of the same fit, and a line under the table
says that DMAR computed them. For ordered items, whose standard
statistic is not the basis of an index, it is the scaled indices: lavaan
leaves the robust indices `NA` at the rungs that constrain the
thresholds under the default delta parameterization, as in the ordered
example below. A ladder of continuous items fit by ULS
(`estimator = "ULS"`, `"ULSM"`, or `"ULSMV"`) withholds the ML-derived
robust RMSEA and incremental indices (Savalei, 2018). Invariance to a
common change of units alone does not establish validity outside their
derivation. Eligible scaled indices are reported with a qualification:
they do not estimate the same population discrepancy as ML-based
indices, and ML cutoffs do not automatically apply. The standard ULS
statistic is in the units of squared covariances, so its CFI and RMSEA
are withheld too. At a rung that restricts the means, as the scalar and
strict rungs do through the intercepts, the ULS fit function adds the
squared mean residuals, in the squared units of the items, to the
squared covariance residuals, in the squared units of the covariances,
so the statistic, the loadings, and every index change with the units:
on the first example under `"ULSM"`, the scalar CFI is .9885 for the raw
scores and .9840 with every score divided by 10, and the scalar delta
CFI -.00005 and -.0046, a judgment against the cutoff of Cheung and
Rensvold (2002) that turns on the units. The `cfi` and `rmsea` columns
of such a rung are `NA`, and the `chi_square` column and the difference
test are kept, since their reference distribution does not depend on the
units, though their values do. A ladder fit by `estimator = "ULS"` alone
is tested by Browne's residual based statistic, which has no units, and
its `cfi` and `rmsea` columns, formed from that statistic, are reported
at the rungs whose ULS estimates do not change with the units, those
that leave the means free. At a rung that restricts the means through
parameters the covariances share, as the scalar and strict rungs do
through the loadings, which enter the means through the latent means,
the ULS estimates change with the units, and so do Browne's statistic
evaluated at them and the indices formed from it: on the first example,
the scalar CFI is .8874 for the raw scores and .8861 with every score
divided by 10, and the scalar delta CFI -.0212 and -.0225. Those columns
are `NA` there too, and the `chi_square` column is kept. lavaan's
difference test of a ULS ladder without a corrected statistic refers the
difference of the ULS statistics, in the units of the squared
covariances and not chi square distributed, to the chi square
distribution, so `delta_chi_square` and `p_value` are `NA` at every
rung, and `"test"` is `NA`; `"ULSM"` and `"ULSMV"` give a scaled
difference test. `estimator = "MLR"`, which does not depend on the
units, is the estimator for nonnormal continuous items. When lavaan
reports no scaled statistic at every rung (beside a bootstrap test,
`test = "bollen.stine"`), the `chi_square` column and the difference
test are the standard ones, and the attributes say so. A saturated rung
(zero degrees of freedom and a statistic of zero) fits exactly in every
version and does not decide, since lavaan can leave its robust CFI `NA`
there from a ratio of zeros, or not compute it. Its statistic is zero in
every version, and with it the model term of the CFI, \\T - c \cdot
df\\, whatever the scaling constant \\c\\ (which is 0/0 there), so its
`cfi` is 1 in every version; its RMSEA, \\\sqrt{(T - c \cdot df)/(N
\cdot df)}\\, is 0/0, which no version defines, and it keeps lavaan's 0,
the convention lavaan applies in every version to a model that fits
exactly, so that `delta_cfi` and `delta_rmsea` at the next rung compare
with exact fit. A rung with zero degrees of freedom whose statistic is
not zero, as when an inequality bound binds (lavaan's count of the
degrees of freedom leaves the bound out), does not fit exactly: its
`cfi` is formed from its statistic and decides the version of the
indices (a ladder whose robust CFI lavaan leaves `NA` at such a rung
falls back to another version at every rung, the standard indices for a
maximum likelihood ladder of continuous items), while its `rmsea`, which
divides by the degrees of freedom, is undefined and `NA`, as is
`delta_rmsea` at the next rung, with a line under the table that says
so; it does not decide the version of the test, which no version defines
on zero degrees of freedom.

Which estimators carry such a test is worth being precise about, because
a naive chi square difference on ordered data is not asymptotically
valid. `"MLM"`, `"MLMV"`, `"MLR"`, `"WLSM"`, `"WLSMV"`, `"ULSM"`, and
`"ULSMV"` carry one. `"DWLS"` and `"ULS"` do not, which is why declaring
ordered items promotes them to `"WLSMV"` and `"ULSMV"`: the discrepancy
function and hence the parameter estimates are untouched, and only the
statistic changes. The remaining case is `"WLS"`, the full weight matrix
estimator, whose statistic is asymptotically chi square under the theory
of Browne (1984) and Muthén (1984) and so needs no correction; the
`"test"` attribute reports the standard difference test there. The
sample size that theory asks for is large, which limits full weighted
least squares to large samples, but that is a separate matter from
whether the difference test is the right one.

**Full information maximum likelihood.** With `missing = "fiml"` (the
default when an item value is missing) or `"ml"`, so that the ladder
uses every case with an observed item, as
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
fits the same model by default, the `chi_square` column is the FIML
likelihood ratio statistic (its scaled version with
`estimator = "MLR"`), and the `cfi` and `rmsea` columns are lavaan's
robust versions, the corrected indices of Zhang and Savalei (2023) that
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
reports for the same model: the complete data fit function evaluated at
the FIML estimates, with the saturated FIML estimates in place of the
sample moments, so that they estimate the values the indices would have
had with complete data (lavaan uses their version 6, which assumes
normality, and their version 3, which does not, with
`estimator = "MLR"`). The usual indices computed from the FIML statistic
move with the amount and kind of missing data (Zhang & Savalei, 2023),
so their change from rung to rung is not the change complete data would
show, and a cutoff for a change in fit, such as the delta CFI of about
-.01 of Cheung and Rensvold (2002), is a statement about the latter.
`delta_cfi` and `delta_rmsea` therefore compare the corrected indices,
and one version holds every rung: when lavaan cannot compute them at
some rung with degrees of freedom, every rung reports its scaled indices
(for a fit with a scaled statistic) or its standard ones, and the line
under the table says which. lavaan cannot combine FIML with the
Satorra-Bentler or the mean and variance adjusted statistic
(`test = "satorra.bentler"` or `"mean.var.adjusted"` passed through
`...`), nor with `se = "robust.sem"`: it fits the complete cases
instead, with a warning. Under the default of `missing` such a ladder is
fitted as it stands when no item value is missing and stops otherwise;
asked for by name, FIML with one of them stops before any fit whether or
not a value is missing, as in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
and so does a ladder whose fits lavaan made by another treatment than
the FIML asked for (it fits `estimator = "GLS"`, `"WLS"`, and `"ULS"` in
two stages).

FIML estimates the means, variances, and covariances of the model's
variables within each group, so every group must observe every variable,
and every pair of them together, at least once. A group that does not is
refused by name before any fit, in an error of class
`dmar_group_coverage` that carries the variables and pairs at issue (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md));
lavaan would stop with an internal error. `holzinger_swineford` is such
a data set, by school: `t25_paper_form_board_r` and `t26_flags` were
given in Grant-White only. The message names what to do: leave the
variable (or one of the pair) out of the model, or an auxiliary variable
out of `aux`, for this comparison, so that the ladder compares the
groups on what every group observes, or, when at least two groups would
remain, leave the groups at issue out of `data`. Listwise deletion is no
remedy, since a group that never observes a variable, or never observes
two together, has no case observed on every variable, and the listwise
ladder would stop on that group instead.

Auxiliary variables passed through `...` as lavaan's `aux` enter a FIML
fit as saturated correlates (Graham, 2003), and lavaan makes the mean
and variance of each one ordinary parameters of each group, which the
intercepts and residuals rungs would equate across the groups; the rung
tests would then test whether the auxiliary is distributed alike in the
groups, which says nothing about the measurement model. They are left
free at every rung, beside any `group_partial` terms, so the ladder has
the degrees of freedom of the ladder without the auxiliary. On six
Holzinger and Swineford tests by school with 60 sentence completion
scores deleted at random and general information (whose means are 44.8
and 36.6 in the two schools) as the auxiliary, the scalar rung's
difference test is 23.88 on 4 df and its delta CFI -.031 (24.03 and
-.031 without the auxiliary).

**Two-stage estimation.** With `missing = "two_stage"` every rung is fit
in two stages (Yuan & Bentler, 2000; Savalei & Bentler, 2009): the means
and covariances are first estimated in each group by maximum likelihood
under the saturated normal model (the EM estimates, which use every
partially observed case and any auxiliary variables passed through `...`
as lavaan's `aux`), and the rung is then fit to them by the complete
data maximum likelihood discrepancy function. A case observed on an
auxiliary variable but on none of the model's variables is read by the
first stage and kept, as in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
which says why: lavaan's two-stage estimator leaves such a case out of
its sample size, so every rung is computed again with its sample size in
each group set to the cases the first stage read there. The statistics,
the indices, the difference tests, and the fits of the `"fits"`
attribute then count every case. The second stage weights each group's
discrepancy by its share of the sample size, so a rung that constrains
parameters across groups has estimates that depend on the groups'
counts, and when those grow in different ratios the rung is estimated
again; the configural rung, which constrains nothing across groups,
keeps lavaan's estimates. On six Holzinger and Swineford tests by school
with the whole block missing for 77 pupils at random given general
information (the auxiliary), 130 and 94 cases with an observed test grow
to 156 and 145 read; the common intercept of lozenges at the scalar rung
is 18.00 (SE .697) rather than lavaan's 18.21 (SE .782), and the scalar
rung's difference test is 14.97 on 4 df, against 10.86 from lavaan's own
ladder, 16.46 under FIML with the same auxiliary, and 17.44 after
listwise deletion. A message says how many such cases there are, and a
line under the table names the sample size. lavaan's
[`summary()`](https://rdrr.io/r/base/summary.html) of a fit in the
`"fits"` attribute counts them as used, and
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html) warns
that a difference test between such a fit and one lavaan made of the
same data on its own sample size is invalid. The standard errors are the
two-stage sandwich, and the `chi_square` column of each rung is the
residual based statistic of Savalei and Bentler (2009, Eq. 8), whose
rejection rate their simulations with incomplete normal data found close
to the nominal level, and the statistic
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
reports after two-stage estimation. With `missing = "robust_two_stage"`
the standard errors are robust to nonnormality (Yuan & Bentler, 2000)
and `chi_square` is lavaan's scaled statistic (Savalei & Falk, 2014).
Under both, `cfi` and `rmsea` are lavaan's robust indices, the corrected
indices of Zhang and Savalei (2023): their TS V2 after robust two-stage
estimation, and its normal theory analog after normal theory two-stage
estimation
([`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
explains the difference), so `delta_cfi` and `delta_rmsea` compare
those, and the rung-to-rung test is the scaled difference test of
Satorra and Bentler (2001), which
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html) applies
to two-stage fits. The difference is formed from the two rungs' naive
statistics and their scaling, so `delta_chi_square` is not the
difference of the `chi_square` entries. Two-stage estimation needs
continuous items and `estimator = "ML"`; robustness to nonnormality
comes from `"robust_two_stage"`, not from `estimator = "MLR"`. The
residual based statistic treats every observed variable as random, so a
model given as syntax that regresses on an observed variable needs
`fixed.x = FALSE` passed through `...`: under `"two_stage"` the
statistic is refused for a fit that conditions on such a variable, as it
is for sampling weights, clustered cases, inequality constraints, and a
parameter on a bound, each with a message that says what to do instead.

A rung whose constraints cost no degrees of freedom tests nothing, so
`delta_chi_square` and `p_value` are `NA` when `delta_df` is zero rather
than reporting a difference in chi square that is numerical noise and
can come out negative. The situation is reported in a message. It arises
with dichotomous items (the ladder above avoids the two rungs where it
is structural) and, for instance, with `group_partial` specifications
that free every parameter a rung would constrain.

Failure at a rung does not end the conversation: partial invariance
(freeing the offending parameter) is the usual next step, for which
`group_partial` refits the whole ladder with named parameters free, and
the `"fits"` attribute gives the fitted objects to
[`lavTestScore`](https://rdrr.io/pkg/lavaan/man/lavTestScore.html) for a
score test of which constraint is doing the damage. This function
deliberately reports the standard ladder rather than automating
modification searches.

A rung whose likelihood has no maximum, flat or still rising as an error
variance in some group falls below zero, has no estimates to compare,
only the point where lavaan's optimizer stopped; the function checks
each rung with a negative error variance as
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
describes and stops with an error of class `dmar_no_maximum` that names
the rung, the group, and the item and prints the lavaan code of that
rung's fit, where `bounds = "pos.var"` (also accepted through `...`
here) or a different model are the choices. Three items whose
correlations in one group have a negative product give such a rung when
`std.lv = TRUE` is passed through `...`; under the default scaling by a
marker item, lavaan's optimizer does not converge on such data. A
negative error variance at a maximum of the likelihood is a finite
estimate, and its rung is reported.

A rung that lavaan's optimizer does not converge on has no estimates to
report, and the function stops with an error of class
`dmar_not_converged`, the class
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
gives such a fit, whose message names the rung. The condition carries
the rung as `level` and the lavaan code that reproduces its fit as
`syntax` (the model) and `lavaan_call` (the call, with `my_data` in
place of the data).

## References

Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit
indices for nonnormality. *Multivariate Behavioral Research, 49*(5),
460–470.
[doi:10.1080/00273171.2014.933697](https://doi.org/10.1080/00273171.2014.933697)

Brosseau-Liard, P. E., Savalei, V., & Li, L. (2012). An investigation of
the sample performance of two nonnormality corrections for RMSEA.
*Multivariate Behavioral Research, 47*(6), 904–930.
[doi:10.1080/00273171.2012.715252](https://doi.org/10.1080/00273171.2012.715252)

Browne, M. W. (1974). Generalized least squares estimators in the
analysis of covariance structures. *South African Statistical Journal,
8*, 1–24.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*(1), 62–83.

Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
equivalence of factor covariance and mean structures: The issue of
partial measurement invariance. *Psychological Bulletin, 105*(3),
456–466.

Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit
indexes for testing measurement invariance. *Structural Equation
Modeling, 9*(2), 233–255.
[doi:10.1207/S15328007SEM0902_5](https://doi.org/10.1207/S15328007SEM0902_5)

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525–543.
[doi:10.1007/BF02294825](https://doi.org/10.1007/BF02294825)

Millsap, R. E. (2011). *Statistical approaches to measurement
invariance*. Routledge.

Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance in
ordered-categorical measures. *Multivariate Behavioral Research, 39*(3),
479–515.
[doi:10.1207/s15327906mbr3903_4](https://doi.org/10.1207/s15327906mbr3903_4)

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. [doi:10.1007/BF02296192](https://doi.org/10.1007/BF02296192)

Savalei, V. (2018). On the computation of the RMSEA and CFI from the
mean-and-variance corrected test statistic with nonnormal data in SEM.
*Multivariate Behavioral Research, 53*(3), 419–429.
[doi:10.1080/00273171.2018.1455142](https://doi.org/10.1080/00273171.2018.1455142)

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables. *Structural
Equation Modeling, 16*, 477–497.
[doi:10.1080/10705510903008238](https://doi.org/10.1080/10705510903008238)

Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
robust full information maximum likelihood with incomplete nonnormal
data. *Structural Equation Modeling, 21*(2), 280–302.
[doi:10.1080/10705511.2014.882692](https://doi.org/10.1080/10705511.2014.882692)

Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
analysis models of different levels of invariance for ordered
categorical outcomes. *Psychometrika, 81*(4), 1014–1045.
[doi:10.1007/s11336-016-9506-0](https://doi.org/10.1007/s11336-016-9506-0)

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
*Sociological Methodology, 30*, 165–200.

Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI
following FIML and TS estimation with missing data. *Psychological
Methods, 28*(2), 263–283.
[doi:10.1037/met0000445](https://doi.org/10.1037/met0000445)

## See also

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) and
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
for the single-group measurement models;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for the reliability of the composite the model justifies;
[`compare_cov_structures`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md)
for covariance-structure comparisons outside the factor model;
[`lavTestScore`](https://rdrr.io/pkg/lavaan/man/lavTestScore.html) for
the score test that localizes a failed rung.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
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
# Do the two schools measure the verbal and the reasoning construct
# in the same way? (Holzinger & Swineford, bundled.) The measurement
# model is a named list of factors; for a single factor, name its
# indicators with items = instead, or give lavaan syntax, as the last
# example does.
data(holzinger_swineford)
hs_factors <- list(
  verbal    = c("t6_paragraph_comprehension", "t7_sentence",
                "t9_word_meaning"),
  deduction = c("t20_deduction", "t22_problem_reasoning",
                "t23_series_completion"))

# The whole ladder, configural through metric, scalar, and strict.
# Configural invariance asks whether the same pattern of loadings
# holds at both schools; metric adds the constraint that the loadings
# are equal across schools, and it is the rung that has to hold before
# a relation involving one of these factors is compared across them;
# scalar equates the intercepts as well, which comparing factor means
# requires; strict equates the residual variances. Each row from the
# second on tests its rung against the one below it.
mi <- measurement_invariance(holzinger_swineford, hs_factors,
                             group = "school")
mi
#>  level      chi_square df p_chi_square cfi   rmsea delta_chi_square delta_df
#>  configural 22.505     16 0.1276       0.992 0.052 <NA>             <NA>    
#>  metric     36.346     20 0.0140       0.979 0.074 13.841           4       
#>  scalar     44.217     24 0.0072       0.974 0.075 7.871            4       
#>  strict     54.343     30 0.0042       0.969 0.073 10.126           6       
#>  p_value delta_cfi delta_rmsea
#>  <NA>    <NA>      <NA>       
#>  0.0078  -0.013    0.022      
#>  0.0964  -0.005    0.001      
#>  0.1195  -0.005    -0.001     

# The fitted models travel with the table, so localizing a failed rung
# costs no refitting.
names(attr(mi, "fits"))
#> [1] "configural" "metric"     "scalar"     "strict"    

# The broom verbs: one row per rung of the ladder, and the model-level
# summary (estimator, test flavor, fit index flavor).
generics::tidy(mi)
#>         term chi_square df p_chi_square       cfi      rmsea delta_chi_square
#> 1 configural   22.50471 16  0.127629786 0.9915883 0.05197395               NA
#> 2     metric   36.34593 20  0.014001841 0.9788619 0.07369221        13.841221
#> 3     scalar   44.21741 24  0.007197934 0.9738554 0.07481503         7.871474
#> 4     strict   54.34304 30  0.004209093 0.9685202 0.07342746        10.125635
#>   delta_df    p_value    delta_cfi  delta_rmsea
#> 1       NA         NA           NA           NA
#> 2        4 0.00781944 -0.012726402  0.021718262
#> 3        4 0.09640134 -0.005006486  0.001122817
#> 4        6 0.11945997 -0.005335160 -0.001387562
generics::glance(mi)
#>   n_levels estimator ordered                                test fit_indices
#> 1        4        ML   FALSE standard chi square difference test    standard

# Partial invariance frees one loading across the schools at every
# rung (Byrne, Shavelson, & Muthén, 1989), so the metric rung costs
# one degree of freedom fewer than it did above. Naming a leading
# subset of the rungs in levels = stops the ladder there, which keeps
# this and the next two examples quick; a reported analysis climbs
# the whole ladder.
measurement_invariance(holzinger_swineford, hs_factors,
                       group = "school",
                       group_partial = "verbal =~ t7_sentence",
                       levels = c("configural", "metric"))
#>  level      chi_square df p_chi_square cfi   rmsea delta_chi_square delta_df
#>  configural 22.505     16 0.1276       0.992 0.052 <NA>             <NA>    
#>  metric     30.789     19 0.0426       0.985 0.064 8.284            3       
#>  p_value delta_cfi delta_rmsea
#>  <NA>    <NA>      <NA>       
#>  0.0405  -0.007    0.012      

# With the indicators coded as ordered categories, ordered = TRUE puts
# the thresholds rung first, because thresholds rather than intercepts
# carry the location information there, and makes the rung-to-rung
# tests the scaled difference tests; the "test" attribute names the
# test that was used. The ladder stops at the thresholds rung here;
# metric, scalar, and strict follow it. lavaan computes no robust fit
# indices for the thresholds rung, so both rows report the scaled ones,
# as the line under the table says.
hs_ordered <- holzinger_swineford
for (item in unlist(hs_factors, use.names = FALSE)) {
  hs_ordered[[item]] <- as.integer(cut(
    holzinger_swineford[[item]],
    breaks = quantile(holzinger_swineford[[item]],
                      c(0, .25, .5, .75, 1)),
    include.lowest = TRUE))
}
mi_ordered <- measurement_invariance(hs_ordered, hs_factors,
                                     group = "school", ordered = TRUE,
                                     levels = c("configural",
                                                "thresholds"))
#> Ordered items declared: switching to estimator = "WLSMV".
mi_ordered
#>  level      chi_square df p_chi_square cfi   rmsea delta_chi_square delta_df
#>  configural 18.012     16 0.3232       0.999 0.029 <NA>             <NA>    
#>  thresholds 27.181     22 0.2044       0.996 0.040 10.330           6       
#>  p_value delta_cfi delta_rmsea
#>  <NA>    <NA>      <NA>       
#>  0.1114  -0.002    0.011      
#> 
#> cfi, rmsea: lavaan did not compute its robust versions at every rung, so these columns hold its scaled versions, which put the scaled statistic into the standard formulas; delta_cfi and delta_rmsea compare them
attr(mi_ordered, "test")
#> [1] "scaled chi square difference test (method = \"satorra.2000\")"

# The measurement model can also be given as lavaan syntax, and
# missing = "fiml", the default, fits the ladder by full information
# maximum likelihood when continuous items have incomplete cases; here
# twenty scores on one item are set to missing first. Naming it states
# the choice in the call. With three indicators
# the single factor is just identified in each group, so the
# configural row is saturated (zero degrees of freedom, exact fit)
# and the metric row carries the test. Under FIML the cfi and rmsea
# columns are the corrected indices of Zhang and Savalei (2023), as the
# line under the table says.
hs_missing <- holzinger_swineford
set.seed(113)
hs_missing$t7_sentence[sample(nrow(hs_missing), 20)] <- NA
measurement_invariance(
  hs_missing,
  model = "verbal =~ t6_paragraph_comprehension + t7_sentence +
           t9_word_meaning",
  group = "school", missing = "fiml",
  levels = c("configural", "metric"))
#>  level      chi_square df p_chi_square cfi   rmsea delta_chi_square delta_df
#>  configural 0.000      0  <NA>         1.000 0.000 <NA>             <NA>    
#>  metric     7.671      2  0.0216       0.987 0.140 7.671            2       
#>  p_value delta_cfi delta_rmsea
#>  <NA>    <NA>      <NA>       
#>  0.0216  -0.013    0.140      
#> 
#> cfi, rmsea: the FIML corrected indices of Zhang and Savalei (2023), version 6 (normal theory), which approximate the values complete data would give and are not computed from the chi_square column; delta_cfi and delta_rmsea compare them
```
