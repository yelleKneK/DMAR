# Multiple-Factor Confirmatory Factor Analysis Model

Fits a confirmatory factor analysis model with one or more factors,
where each factor is specified by naming its indicator variables and the
measurement structure is specified by describing what is constrained
(`equal_loading`, `equal_intercept`, `equal_error`) rather than by more
technical terms. The function then reports which classical measurement
structure the description implies (congeneric, essentially
tau-equivalent, tau-equivalent, essentially parallel, or parallel), the
parameter estimates with confidence intervals, fit information, and, per
factor, coefficient omega, the average variance extracted (AVE), and
coefficient *H*, each with a delta method standard error and confidence
interval computed by lavaan from defined parameters (no additional
packages are involved). `cfa_k()` is an alias of this function.

## Usage

``` r
fa_confirmatory(
  data = NULL,
  factors,
  S = NULL,
  N = NULL,
  M = NULL,
  equal_loading = FALSE,
  equal_intercept = FALSE,
  equal_error = FALSE,
  correlated_factors = TRUE,
  meanstructure = NULL,
  estimator = "ML",
  missing = "fiml",
  ordered = NULL,
  se = NULL,
  conf_level = 0.95,
  output = c("verbose", "measurement", "summary", "standardized", "fit"),
  ...,
  wald_df = Inf,
  rmsea_conf_level = 0.9
)

cfa_k(
  data = NULL,
  factors,
  S = NULL,
  N = NULL,
  M = NULL,
  equal_loading = FALSE,
  equal_intercept = FALSE,
  equal_error = FALSE,
  correlated_factors = TRUE,
  meanstructure = NULL,
  estimator = "ML",
  missing = "fiml",
  ordered = NULL,
  se = NULL,
  conf_level = 0.95,
  output = c("verbose", "measurement", "summary", "standardized", "fit"),
  ...,
  wald_df = Inf,
  rmsea_conf_level = 0.9
)
```

## Arguments

- data:

  A raw data matrix or data frame, rows are respondents and columns
  include the items named in `factors`. Supply exactly one of `data` or
  `S`.

- factors:

  Named list. Each element names a factor and gives the character vector
  of its indicator columns (two or more per factor; three or more when
  only one factor is specified). Each item loads on exactly one factor
  (simple structure). An item that is not a column of `data` or an item
  of `S` is refused with an error of class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

- S:

  A symmetric covariance matrix of the items, with dimnames naming the
  items; `N` is then required. Supply exactly one of `data` or `S`.

- N:

  Total sample size. Required with `S`; ignored (inferred from the rows)
  with `data`.

- M:

  Optional named numeric vector of item means, used with a covariance
  matrix to model the mean structure (required there when
  `equal_intercept` is used). Ignored for raw data.

- equal_loading:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the loadings within a factor to a single value (across
  factors nothing is equated). Defaults to `FALSE`.

- equal_intercept:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the item intercepts within a factor to a single value; this
  requires the mean structure (raw data, or `M` with a covariance
  matrix). Defaults to `FALSE`.

- equal_error:

  Logical, or a named logical vector with one element per factor. `TRUE`
  constrains the error variances within a factor to a single value.
  Defaults to `FALSE`.

- correlated_factors:

  Logical. If `TRUE` (default) the factors covary freely; because each
  factor variance is fixed to 1, the `phi` terms for factor pairs are
  the latent correlations. If `FALSE` the factor covariances are fixed
  to zero.

- meanstructure:

  Logical or `NULL`. `NULL` (default) models the mean structure exactly
  when it is needed: under FIML or two-stage estimation of raw data (the
  default, so the item intercepts, the `nu` rows, join the table), when
  any `equal_intercept` is `TRUE`, or when `M` is supplied. Set `TRUE`
  to model intercepts regardless (raw data or `M` required), or `FALSE`
  to suppress them (an error if `equal_intercept` is used).

- estimator:

  Character; estimator passed to lavaan. Must be one of `"ML"` (default;
  maximum likelihood, fully efficient under multivariate normality),
  `"MLR"` (robust maximum likelihood, raw data required: maximum
  likelihood estimates with the test statistic corrected for
  nonnormality, the scaled statistic of Satorra and Bentler, 1994, in
  the form Yuan and Bentler, 2000, gave it for incomplete data; its
  matching Huber-White sandwich standard errors are used when `se` is
  omitted), `"WLS"` (the asymptotic distribution free estimator of
  Browne, 1984; raw data and a large sample required), `"WLSMV"`
  (diagonally weighted least squares with mean- and variance-adjusted
  test statistic; Muthén, 1984; Muthén, du Toit, & Spisic, 1997; the
  standard choice for ordered categorical items, and what `ordered`
  switches to), or `"GLS"` (generalized least squares; Browne, 1974).
  With a robust estimator the reported fit indices are the robust
  versions.

- missing:

  Character; missing-data handling when raw data are supplied: `"fiml"`
  (default; full information maximum likelihood, which uses every case
  with an observed item; lavaan's name `"ml"` is accepted), `"listwise"`
  (listwise deletion), `"pairwise"` (passed to lavaan, for ordered
  items), `"two_stage"` (two-stage maximum likelihood with normal theory
  standard errors and test), or `"robust_two_stage"` (two-stage maximum
  likelihood with standard errors and test robust to nonnormality);
  lavaan's spellings `"two.stage"` and `"robust.two.stage"` are accepted
  for the last two. See *Details* for two-stage estimation. Ignored with
  `S`, except that a two-stage request with `S` is an error, since a
  covariance matrix has no incomplete cases. With `ordered` items,
  `"ml"`/`"fiml"` and the two-stage options are not available; use
  `"pairwise"` or `"listwise"`. FIML needs a maximum likelihood
  estimator of continuous items with normal theory or Huber-White
  standard errors. When the default meets ordered items, a least squares
  estimator (`"WLS"`, `"WLSMV"`, `"GLS"`), `se = "robust.sem"`, or a
  Satorra-Bentler kind of test passed through `...`, the items are
  fitted as they stand if no item value is missing (listwise deletion
  then deletes nothing), and otherwise the call stops and names the
  option to choose, rather than deleting cases unasked. With complete
  data FIML and listwise deletion give the same estimates and the same
  chi square, and for `estimator = "ML"` the same fit indices; the
  standard errors differ slightly, since FIML takes them from the
  observed information, as lavaan does and as Savalei (2010) recommends
  with incomplete data, and listwise deletion from the expected
  information.

- ordered:

  Ordered categorical items: `NULL` (none, the default), `TRUE` (every
  item), or a character vector of item names. Requires raw data. Each
  factor must be all ordered or all continuous. Declaring ordered items
  switches the estimator to `"WLSMV"` (with a message, unless a
  categorical estimator was requested), fits thresholds in place of
  intercepts, and reports each ordered factor's omega on the categorical
  sum score metric (categorical omega with the unrestricted polychoric
  denominator); see *Details*.

- se:

  Standard error method passed to lavaan. The default, `NULL`, follows
  the estimator: model-based `"standard"` errors for ML, GLS and WLS;
  Huber-White sandwich `"robust.huber.white"` errors for MLR; and
  `"robust.sem"` errors for WLSMV. Explicit compatible choices are
  retained. The ML and MLR point estimates use the same likelihood;
  their standard errors and tests use different reference assumptions.
  This difference does not depend on whether optimization is iterative
  or an equivalent solution is computed analytically. `"robust.sem"`
  requires complete data and is refused with FIML; use
  `"robust.huber.white"` to retain incomplete respondents. No
  standard-error request silently changes the missing-data treatment.
  `"none"` requests point estimates only. Under two-stage estimation the
  matching corrected two-stage errors are used; `se` must be `NULL`,
  `"standard"`, or `"none"`. MLR and Huber-White errors require raw
  data. Robust SEM errors from covariance input require ML and a
  matching `NACOV` in `...`.

- conf_level:

  Confidence level for the parameter confidence intervals, including the
  delta method intervals for omega, AVE, and *H*. Defaults to 0.95. It
  does not reach the RMSEA interval, whose level is `rmsea_conf_level`.

- output:

  Format of the returned object:

  `"verbose"`

  :   (default) Parameter estimates with confidence intervals, the
      per-factor defined parameters (`loading_sum`, `error_sum`,
      `omega`, `ave`, `H`), and fit information.

  `"measurement"`

  :   The measurement-property rows only: per factor `omega`, `ave`, and
      `H` (with delta method standard errors and confidence intervals),
      the latent correlation `phi` for every factor pair (with its
      confidence interval), and, for raw data, the heterotrait-monotrait
      ratio `htmt` for every factor pair via
      [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md).

  `"summary"`

  :   The raw [`summary()`](https://rdrr.io/r/base/summary.html) output
      from lavaan (not a data frame).

  `"standardized"`

  :   The standardized parameter estimates from
      [`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html),
      whose defined rows describe the standardized items: `omega` there
      is the reliability of the unit-weighted sum of the standardized
      items, while `ave` and `H` do not depend on the metric of the
      items and equal the rows of the other tables. For a factor with
      `equal_loading` or `equal_error` these rows are computed from each
      item's own standardized loading and error variance (see *The
      standardized solution* under *Details* of `fa_confirmatory`).
      Under `se = "none"` the `se`, `z`, `pvalue`, `ci.lower`, and
      `ci.upper` columns are `NA`.

  `"fit"`

  :   The raw lavaan fit object. Escape hatch for direct lavaan access,
      including likelihood ratio tests between two `fa_confirmatory()`
      fits via
      [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html).

- ...:

  Additional arguments forwarded to
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html) (e.g.,
  `cluster` or `bootstrap`). A `group` argument is an error: see
  *Details* for why, and
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  for the model in several groups.

- wald_df:

  The reference distribution of the Wald intervals and *p*-values of the
  parameter rows: `Inf` (the default), the normal distribution, which
  leaves lavaan's intervals and *p*-values as they are; a single
  positive number, the degrees of freedom of a *t* distribution for
  every parameter row, the intervals and *p*-values recomputed from the
  same estimates and standard errors; or `"missing_information"`, von
  Hippel's (2016) degrees of freedom \\\max\\3, N(1 - \gamma_j) - 1\\\\
  for each loading, from its fraction of missing information
  \\\gamma_j\\ (\\N - 1\\ with complete data), with every other row left
  on the normal reference. A number applies to the tables of
  `output = "verbose"`, `"measurement"`, and `"standardized"`, the
  keyword to that of `"verbose"` only (the `"measurement"` table has no
  loading), and either needs standard errors other than `"none"` and
  `"bootstrap"`; the keyword also needs the fits its calibration covers,
  maximum likelihood (`estimator = "ML"` or `"MLR"`) of continuous items
  by FIML or with complete data, with normal theory or Huber-White
  standard errors. Anything else is an error that says why. It follows
  `...`, so it must be named in full. See *Details*.

- rmsea_conf_level:

  Confidence level of the RMSEA interval of the fit rows, which the
  `rmsea_ci_level` row records. Defaults to 0.90, the conventional
  level, whose limits decide the tests of close and not close fit at .05
  (see
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md));
  `rmsea_conf_level = 0.95` gives a 95 percent interval. It sets the
  level in the table DMAR returns and in its `glance()`;
  `output = "summary"` returns the lavaan summary, which shows lavaan's
  90 percent interval whatever this argument is.

## Value

For `output = "verbose"` (default) and `output = "measurement"`, a
`data.frame` (classes `dmar_fa_confirmatory`, `dmar_cfa_k`, and
`dmar_tbl`) with columns `syntax`, `term`, `estimate`, `se`,
`test_statistic`, `p_value`, `ci_lower`, `ci_upper`. The `"model"`
attribute is a named character vector giving, per factor, the implied
classical structure; the printed header displays it. For
`output = "verbose"` the `"test"` attribute names the statistic in the
`chi_square` row: `"standard"`, `"scaled"` (a robust or categorical
estimator), or, under the two-stage options, the two-stage residual
based or scaled statistic; a two-stage fit also carries the
`"two_stage"` attribute, a list with the method, the lines printed under
the table, and, in `rows`, the cases the first stage read and lavaan's
count of them when an auxiliary variable was named, and the `"fit_rows"`
attribute holds the lines that name the version of the indices after
FIML estimation and the SRMR of a model with a mean structure (`NULL`
when there are none). The tables of `"verbose"` and `"measurement"` also
carry the `"wald_df"` attribute, the degrees of freedom of the *t*
reference of their intervals and *p*-values (`Inf` for the normal one);
under `wald_df = "missing_information"` it holds the degrees of freedom
of each loading in the table, named by term (the other rows keep the
normal reference), the `"wald_df_table"` attribute is the
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
table they came from, and `"wald_df_rows"` names what they are. For
`output = "summary"`, the lavaan summary object; for `"standardized"`,
the standardized solution, with the same `"wald_df"` attribute; for
`"fit"`, the lavaan fit object.

## Details

With `ordered` items, the model is fit by WLSMV to polychoric
correlations with thresholds. A sum score of ordered items lives on the
metric of the observed categories, not on the latent response metric of
the polychoric loadings, so for a factor whose items are ordered the
reported omega is categorical omega, the reliability of the sum of the
observed categories (Green & Yang, 2009), computed from the same fit
with the total variance built from the unrestricted polychoric
correlations (Kelley & Pornprasertmanit, 2016, Eq. 18), the default
denominator of
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md);
the substitution is announced in a message and recorded in the
`omega_metric` attribute, and the delta method interval columns are `NA`
for those rows because the delta method interval describes the latent
response metric. When the data do not determine a polychoric correlation
among a factor's items (a bivariate table with empty cells, as two items
skewed in opposite directions produce), that factor's omega is `NA` with
a warning naming the pair. AVE and coefficient *H* concern the latent
response correlations themselves and are reported unchanged on that
metric. An improper solution (a negative error variance) leaves some of
these properties undefined, and their rows are then `NA` with a warning;
see *Details*. For a bootstrap confidence interval on a categorical
omega, use
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
on the factor's items.

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) and
[`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) are
convenience wrappers around this function for the one and two factor
cases: [`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
takes a vector of items and fits one factor over them, and
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) takes
the items of each of two factors. Both forward every argument here, so
their results are this function's results, with the factors named `f1`
(and `f2`).

**Describing the model instead of naming it.** The classical measurement
structures are nested patterns of within-factor equality constraints
(Lord & Novick, 1968; Graham, 2006):

- congeneric:

  loadings, intercepts, and error variances all free.

- essentially tau-equivalent:

  equal loadings; intercepts and error variances free.

- tau-equivalent:

  equal loadings and equal intercepts; error variances free.

- essentially parallel:

  equal loadings and equal error variances; intercepts free.

- parallel:

  equal loadings, equal intercepts, and equal error variances.

The caller states the constraints; the function reports the implied
name, per factor, in the printed header and in the `"model"` attribute
of the returned table. The distinction between tau-equivalent and
essentially tau-equivalent (and between parallel and essentially
parallel) lives entirely in the mean structure: the covariance structure
of the two members of each pair is identical, so without intercepts in
the model only the "essentially" form can be claimed. That is why
`equal_intercept` requires raw data or `M`: covariances alone cannot
speak to it. A constraint pattern outside the classical list (for
example equal error variances with free loadings) is fit as requested
and labeled descriptively, since it has no conventional name.

Identification fixes each factor variance to 1 (and each factor mean to
0 when the mean structure is modeled), so all loadings are estimated and
within-factor equality constraints are meaningful. Two
`fa_confirmatory()` fits that differ only in descriptor settings are
nested, so `output = "fit"` feeds
[`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html)
directly (with `estimator = "MLR"`, lavaan applies the scaled difference
test).

**Measurement properties.** For factor \\f\\ with unstandardized
loadings \\\lambda_j\\ and error variances \\\psi_j\\ (factor variance
1): coefficient omega \\= (\sum_j \lambda_j)^2 / ((\sum_j \lambda_j)^2 +
\sum_j \psi_j)\\ (McDonald, 1999), the reliability of the unit-weighted
composite; the average variance extracted \\= J^{-1} \sum_j \lambda_j^2
/ (\lambda_j^2 + \psi_j)\\ (Fornell & Larcker, 1981), the mean
proportion of item variance the factor accounts for; and coefficient \\H
= (1 + (\sum_j \lambda_j^2/\psi_j)^{-1})^{-1}\\ (Hancock & Mueller,
2001), the reliability of the optimally weighted composite, which no
single item can drag below its value for any subset. All three are
computed as lavaan defined parameters, so each carries a delta method
standard error and a `conf_level` confidence interval in the same table
as the model parameters. The limits of all three are formed on the logit
scale from that standard error and transformed back to the scale of the
coefficient (Browne, 1982; the confidence interval transformation
principle of Kelley, 2007: the limits of a monotone transformation of a
parameter are the transformation of its limits), as
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
and
[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
form those of *H* by default, so they lie within (0, 1) and are
asymmetric near one, as the sampling distribution of a high reliability
is; for *H* with one factor they are the limits those two functions give
the same fit. The `se` column is the standard error of the coefficient
itself. lavaan's own limits for these defined parameters, the estimate
plus or minus the normal quantile times the standard error, can exceed
one.

**The standardized solution.** With `output = "standardized"` the
defined rows are evaluated at the standardized values: the standardized
loading \\s_j = \lambda_j / (\lambda_j^2 + \psi_j)^{1/2}\\ and error
variance \\e_j = \psi_j / (\lambda_j^2 + \psi_j)\\ of each item take the
places of \\\lambda_j\\ and \\\psi_j\\ in the formulas above (for an
ordered item, those of its latent response). Omega is then the
reliability of the unit-weighted sum of the standardized items, which
differs from the omega of the sum of the items as they are scored
whenever the items' variances differ; the average variance extracted and
*H* are the same on either scale. lavaan evaluates a defined parameter
at the standardized values of the labels it names, and a label that
`equal_loading` or `equal_error` gives several items has one
unstandardized value but a standardized value for each item, which
differ when the items' variances differ; lavaan uses the first item's.
For such a factor the rows are computed here from every item's own
\\s_j\\ and \\e_j\\, with delta method standard errors from their
derivatives with respect to the free parameters and lavaan's covariance
matrix of the estimates (under `se = "bootstrap"`, from lavaan's
bootstrap draws), and with the test statistic and *p*-value, as lavaan
forms them for every other row.
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
does the same for its indices.

**Improper solutions.** When the error variance \\\psi_j\\ of an item is
estimated below zero (a Heywood case; for an ordered item, the residual
variance of its latent response, which the theta parameterization
estimates), the item's communality \\\lambda_j^2 / (\lambda_j^2 +
\psi_j)\\ exceeds one. The average variance extracted averages the
communalities, and coefficient *H*, the largest reliability that a
weighted sum of the items attains, is at least that communality, since
the item weighted alone attains it; for such a factor neither is a
proportion of variance, so its `ave` and `H` rows are `NA`, with a
warning naming the items and their communalities. Omega, the reliability
of the unit-weighted sum, depends on the error variances only through
their sum, and while that sum is positive it remains a proportion and is
reported, as
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
reports it, beside the warning; when the sum is negative the model
implies a negative error variance for the sum score and an omega above
one, and the `omega` row is `NA` as well. Categorical omega is not
defined for an improper ordered solution, so an ordered factor's `omega`
row is `NA` whenever one of its items is improper. On
`holzinger_swineford`, a factor measured by visual perception and
lozenges alone, beside one measured by general information, paragraph
comprehension, and sentence completion, gives visual perception an error
variance below zero and the first factor an *H* of 1.32; that row is
`NA`, while its omega, .71, is reported. A negative error variance can
also mean that the likelihood has no maximum, the variance falling
without bound while the likelihood still rises and lavaan's optimizer
stopping where its steps stop gaining; the function checks for this by
fitting the model again with the variance held further out (see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md))
and, when the log likelihood there is higher or lower by no more than
\\10^{-6}\\ per case, stops for every `output` with an error of class
`dmar_no_maximum` that prints the lavaan code of the fit, where
`bounds = "pos.var"` (also accepted through `...` here) or a different
model are the choices. The lavaan fit (`output = "fit"`) keeps the
defined parameters as lavaan computes them.

**Discriminant validity.** Three complementary readings come from
`output = "measurement"`: (a) the latent correlation `phi` for a factor
pair, with a confidence interval whose upper limit near 1 means the data
cannot distinguish the two factors; (b) the Fornell and Larcker (1981)
comparison, which asks whether each factor's `ave` exceeds the squared
`phi` of its pairs (the factor should share more variance with its own
items than with the other factor); and (c) for raw data, the model-free
`htmt` ratio (Henseler, Ringle, & Sarstedt, 2015). A pair whose ratio
the data leave undefined, as when some items of two factors are never
observed together in a planned missing data design whose forms share
other items, has an `htmt` row of `NA`, with a warning that names it,
while the model based rows are reported: FIML (or two-stage estimation)
estimates the latent correlation from the pairs of their items that are
observed together, through the factor structure. When no item of one
factor is ever observed with an item of the other, as when the two
factors were given to different respondents, the latent correlation
enters the likelihood only through covariances the data never observe,
so it is not identified (the likelihood is flat in it, and lavaan leaves
it at its starting value with no standard error, or under two-stage
estimation fits it to an arbitrary covariance), and `fa_confirmatory()`
refuses the fit, naming the pair: fit each factor alone, whose `omega`,
`ave`, and `H` are identified, or link the forms with items given on
both
([`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md)
reports the cases behind every pair of items). Orthogonal factors
(`correlated_factors = FALSE`) are not refused: their correlation is
fixed at 0 rather than estimated, so nothing is left unidentified, and
since the likelihood then separates into the parts of the factors'
items, each factor's `omega`, `ave`, and `H`, with their standard
errors, are those of its own one factor fit (to the accuracy of the
optimizer). The fit rows of such a fit are lavaan's, whose degrees of
freedom count the covariances no case observes among the moments tested
(lavaan warns that their coverage is zero), though the data cannot test
them; each factor's fit is better judged from its own fit. More
generally, when lavaan cannot compute the standard errors of a fit (its
information matrix cannot be inverted, which suggests a parameter that
is not identified at the estimates), every standard error, test, and
interval is `NA`, with a warning of class `dmar_result_qualified`,
rather than the 0 lavaan reports for the defined quantities. So they are
when lavaan computes them from an information matrix that is singular in
the directions the equality constraints leave free (scaled to a unit
diagonal, its smallest eigenvalue is below \\10^{-6}\\), which it
inverts, sometimes with a warning that the covariance matrix of the
estimates does not appear to be positive definite and sometimes with
none, reporting finite values that do not estimate the sampling
variability of the estimates: the observed information of a FIML fit and
the expected information of a listwise fit or of a covariance matrix
alike, as in
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
and
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md).
The rows report the numbers and their uncertainty; the judgment is the
researcher's.

**Two-stage estimation with incomplete data.** With
`missing = "two_stage"` the model is estimated in two stages (Yuan &
Bentler, 2000; Savalei & Bentler, 2009). The first stage estimates the
means and covariances of the items by maximum likelihood under the
saturated normal model (the EM estimates), from every case with an
observed item and from any auxiliary variables named through lavaan's
`aux` argument, which enter this stage only, so the model itself is
untouched by them. A case observed on an auxiliary variable but on no
item is used as well: the first stage estimates the items and the
auxiliaries jointly, so its auxiliary values inform the item moments,
and its likelihood is summed over every case (Savalei & Bentler, 2009,
Eq. 10). That is what makes the estimates consistent when the items go
missing together at random given the auxiliary variables, the case an
auxiliary is named for. lavaan's two-stage estimator reads such a case
in its first stage but leaves it out of its sample size, and that sample
size turns the per case quantities of the two stages into every standard
error (each too large by the square root of the ratio of the two
counts), the test statistics, the fit indices, and the information
criteria. The fit is therefore computed again with its sample size set
to the cases the first stage read, from the same estimates (in one group
they do not depend on the sample size), so that every number in the
table, and in the lavaan fit of `output = "fit"`, counts them: lavaan's
[`summary()`](https://rdrr.io/r/base/summary.html) of that fit counts
them as used, and
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html) warns
that a difference test between it and a fit lavaan made of the same data
on its own sample size is invalid. A message says how many such cases
there are, and a line under the table names the sample size. When the
fit cannot be computed again, a warning says what rests on lavaan's
sample size instead. On six Holzinger and Swineford tests with the whole
block missing for 77 of 301 pupils, at random given general information
(the auxiliary), the standard error of a verbal loading of 2.840 is .190
rather than lavaan's .220. Cases observed only on auxiliary variables
are kept, as under FIML. Such cases still contribute to the
observed-data likelihood through their auxiliary measurements. The
second stage fits the model to those moments with the complete data
maximum likelihood discrepancy function. The estimates are consistent
when the data are missing at random, as FIML's are, and Savalei and
Bentler (2009) found them only marginally less efficient. What the
second stage cannot supply alone is a standard error or a test, because
it treats the EM moments as the moments of complete data, which makes
its naive chi square and standard errors too optimistic. The standard
errors are therefore the two-stage sandwich of Savalei and Bentler
(2009, Eq. 6), built from the observed information of the saturated
model under missing data, and the `chi_square` row is their residual
based statistic (Eq. 8; Eq. 14 with auxiliary variables), which is
referred to chi square on the model's degrees of freedom. In their
simulations with normal data and 15 or 30 percent of the values missing
it rejected a correct model in no more than 9 percent of the tests at
the .05 level with samples of 200, where the FIML chi square reached
21.9 percent. With `missing = "robust_two_stage"` the standard errors
are the sandwich that is robust to nonnormality (Yuan & Bentler, 2000)
and the `chi_square` row is lavaan's scaled statistic, the pairing
Savalei and Falk (2014) found to outperform robust FIML with incomplete
nonnormal data. Under both, `cfi` and `rmsea` are lavaan's robust
indices, which correct the complete data indices of the second stage for
the missing data as Zhang and Savalei (2023) propose: each subtracts
from the naive statistic *T* an estimate *c* of its expectation under
the model, built from the asymptotic covariance matrix of the stage one
moments. With `missing = "robust_two_stage"` that matrix is the sandwich
of the observed and first order information of the saturated model, and
the indices are their TS V2 (Table 2), the version they recommend after
two-stage estimation, for incomplete nonnormal data as well. With
`missing = "two_stage"` it is the inverse of the observed information,
the form the sandwich takes under normality, so the indices are the
normal theory analog of TS V2, which Zhang and Savalei did not study: it
estimates the same population values and agrees closely with TS V2 when
the data are normal, but it departs from TS V2 as the data depart from
normality, where `"robust_two_stage"` gives TS V2 itself. On nine
Holzinger and Swineford tests with 90 scores on each of two deleted at
random, the two RMSEAs are .0869 and .0859. Zhang and Savalei define
only the RMSEA and the CFI; the `tli` and `nnfi` rows apply the same
correction to the Tucker-Lewis index, and the RMSEA interval inverts the
scaled statistic *T*/*c* and rescales its limits by *c*, both of them
lavaan's extensions. Indices from a FIML fit are never paired with
two-stage estimates. The `AIC` and `BIC` rows are the two-stage
information criteria of Lai (2021), built, like `H0` and `H1`, from the
complete data normal log likelihood evaluated at the EM estimates; Lai
shows that, unlike the FIML criteria, they select in large samples the
model that complete data would select, and they are comparable only with
other two-stage fits of the same variables
([`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
reports the same values). The difference \\-2(H0 - H1)\\ is the naive
statistic, which is not chi square distributed. Two nested two-stage
fits from `output = "fit"` are compared with
[`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html),
which applies the scaled difference test (Satorra & Bentler, 2001). The
`test` attribute names the statistic in the `chi_square` row, and the
lines under the printed table say which statistic and indices it
reports. Two-stage estimation is a normal theory method for raw
continuous data: ordered items, covariance input, and an estimator other
than `"ML"` are refused, robustness coming from `"robust_two_stage"`
rather than from `estimator = "MLR"`. Under `"two_stage"` the residual
based statistic is refused, with a message that points to robust
two-stage estimation, for a fit with sampling weights or clustered
cases, inequality constraints, or a parameter on a bound passed through
`...`, since the statistic assumes independent, equally weighted cases
and a projection on the parameters left free. The residual based
statistic already sums the information of the saturated model over every
case the first stage read, so the correction for cases observed on an
auxiliary variable alone leaves it as it was.

**One group.** The model is fit in one group. Every loading, error
variance, intercept, and factor correlation carries a label, which is
how the equality descriptors and the defined parameters refer to it, and
lavaan holds a label shared by several groups equal across them. With
`group` passed through `...` the groups would therefore share every
measurement parameter, with each factor's variance fixed at 1 (and its
mean at 0) in every group: a model of complete invariance rather than a
model fit by group, whose omega, `ave`, and *H* belong to the pooled
model rather than to any one group. A `group` argument is therefore an
error.
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
fits the model in several groups, beginning with the configural model,
which leaves every parameter free in each group;
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
reads *H* for one factor in one group from its fit, and
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
on one group's rows gives that group's omega.

**Confidence interval conventions.** Parameter rows (including omega,
`ave`, *H*, and `phi`) use `conf_level`, 0.95 by default. The RMSEA
interval uses `rmsea_conf_level`, 0.90 by default, the conventional
level for the RMSEA, because a \\100(1 - 2\alpha)\\ percent interval
decides the one-sided tests of close and not close fit at level
\\\alpha\\ (MacCallum, Browne, & Sugawara, 1996, p. 136); the
`rmsea_ci_level` row records that level, and the line beneath the
printed table names it when it differs from `conf_level`.

The parameter intervals are Wald intervals, the estimate plus or minus a
normal quantile times its standard error, and the *p*-values refer the
ratio of the two (`test_statistic`) to the normal distribution. The
normal reference treats the standard error as known, which it is only in
the limit, and in small samples, more so with missing values, the
intervals are too short (von Hippel, 2016). `wald_df` refers both to the
*t* distribution on that many degrees of freedom instead (von Hippel's
Eq. 29): the limits become the estimate plus or minus the *t* quantile
times the same standard error (for the omega, AVE, and *H* rows, on the
logit scale) and the *p*-value is \\2 P(T\_\nu \> \|z\|)\\; the
`test_statistic` column holds the same ratio under either reference, and
a line under the printed table names the reference distribution. One
number applies to every parameter row, a sensitivity check of the
analyst's own choosing. The fit rows, the RMSEA interval, categorical
omega, and the HTMT rows are untouched.

**Degrees of freedom from the missing information.** With
`wald_df = "missing_information"` each loading gets its own degrees of
freedom by von Hippel's (2016, Eqs. 31 and 34) effective sample size
rule with \\k = 1\\, \$\$\nu_j = \max\\3, N(1 - \gamma_j) - 1\\,\$\$
where \\N\\ is the number of cases the fit used and \\\gamma_j\\ the
parameter's fraction of missing information, estimated from the fit as
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
estimates it (the unstructured estimator of Chen and Savalei, 2021,
which cannot be negative);
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
describes the rule and its floor. A loading of an item that is rarely
missing keeps nearly \\N - 1\\ degrees of freedom, and with complete
data (listwise deletion, or a covariance matrix) every \\\gamma_j\\ is 0
and every \\\nu_j\\ is \\N - 1\\. The `"wald_df"` attribute holds each
parameter's degrees of freedom, the `"wald_df_table"` attribute the
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
table they came from, and the line under the printed table names the
rule and the rows it applied to.

The rule, and \\k = 1\\ rather than \\k\\ equal to the number of items
or von Hippel's complete data rule (his Eq. 32) on \\N - 1\\ or \\N\\
less the number of items, come from a calibration study of these
intervals under full information maximum likelihood (in the package's
source repository, `methodology/missing_data`). It simulated one factor
models as this function writes them (loadings from .8 to .5, the factor
variance fixed at 1), with six items at *N* = 30, 50, 100, and 200 and
with 12 items at *N* = 50 and 100, normal data fit by ML with normal
theory standard errors and nonnormal data (the factor and the errors
drawn from a standardized chi square on 4 degrees of freedom, skewness
1.41 and excess kurtosis 3) fit by MLR with Huber-White standard errors,
with no missing values or with 15 or 30 percent of the values of half
the items missing completely at random or at random, in 2000
replications of each of 60 cells; a repetition of von Hippel's own
design reproduced the coverage he published to within .012. Under a
decision rule fixed before the run (mean absolute distance of the
coverage from .95, over the cells with \\N \le 100\\ and missing values,
below the normal interval's; coverage above .96 in no cell; a mean
distance within .005 of the best rule's; the simplest passing rule
chosen, a cell's coverage being the mean over the target's parameters),
\\k = 1\\ passed for the loadings, with a mean distance of .0171 against
.0247 for the normal reference and no cell's coverage, averaged over its
loadings, above .9502 (single loadings, the largest of 384 with missing
values, reached .9595, with a Monte Carlo standard error of .0044, and
none exceeded .96), and for omega, .0180 against .0242 with no coverage
above .9550. It passed again within the normal theory cells and within
the robust ones taken alone. With nonnormal data the intervals still
covered below .95 (between .908 and .937 for the loadings, averaged over
a cell's loadings), since the robust standard errors themselves are too
small at these sample sizes: the *t* reference narrows the shortfall and
does not remove it.

These coverage figures rest on the data sets that the study's
eligibility rule for the unrestricted (H1) estimates admitted. The
affine-collapse diagnostic sets aside some of those data sets, and the
coverage over the data sets that remain has not been established.
Replaying the eligibility rule on the raw data cannot by itself validate
the refits of the fraction of missing information and the information
diagnostics. The keyword is an approximate, opt-in rule, not a coverage
guarantee.

The keyword is restricted to the study's estimator and parameter types.
The other parameter rows keep the normal reference: the study evaluated
the unstandardized loadings and omega, not factor variances or
correlations, intercepts, error variances, or the sums, AVE, and *H*
built from them, and nothing shows that the rule suits them; omega
passed with its limits formed on its own scale, but its limits here are
formed on the logit scale, where the study found the *t* reference
overcovered (.9645 with normal data), so it keeps the normal reference
as
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)'s
logistic intervals do. It is an error with the `"measurement"` table,
which has no loading, with the standardized solution, with an estimator
other than ML and MLR or ordered items, with pairwise deletion or
two-stage estimation (for which no fraction of missing information is
computed), and with standard errors other than the normal theory and
Huber-White ones (`se = "standard"` and `"robust.huber.white"`;
`"robust.sem"` only under listwise deletion, or from a covariance matrix
with the asymptotic covariance matrix of its moments passed to lavaan as
`NACOV` through `...`, without which lavaan refuses it, since lavaan
computes that sandwich from raw data or a supplied `NACOV`, and under
FIML falls back to listwise deletion, so that `missing = "fiml"` would
not be honored). More generally, it is an error wherever lavaan replaces
the FIML fit that `missing = "fiml"` asks for with a fit to the complete
cases, which it also does for a test statistic of the Satorra-Bentler
kind passed through `...` (`test = "satorra.bentler"`,
`"mean.var.adjusted"`, `"scaled.shifted"`, and their corrected forms),
whether or not a value is missing: `estimator = "MLR"`, whose scaled
statistic is computed under FIML, or the test with
`missing = "listwise"` is the alternative. It is an error, too, where
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
cannot compute the fraction for the data: a pair of items that no case
observes together, or an unrestricted estimate of the means and
covariances that does not exist (see its *Fits that are refused*). Under
FIML the standard errors must come from the observed information
computed from the Hessian of the fitted model, lavaan's default there
and the one the study used; the first order information and the observed
information of the unrestricted model (`observed.information = "h1"`),
which reach lavaan only through `...`, are an error. Under listwise
deletion (or from a covariance matrix) lavaan's normal theory standard
errors use the expected information where the study's cells without
missing values used the observed; at *N* = 30 and 50 they are 1 to 2
percent shorter for the loadings and less than 1 percent shorter for
omega, so the complete data intervals here are, if anything, slightly
shorter than those the study evaluated, whose coverage with normal data,
averaged over a cell's loadings, rose from between .936 and .946 to
between .945 and .949 on \\N - 1\\ degrees of freedom (single loadings
between .934 and .963, where the normal reference reached .962). Several
factors, equal loadings, and auxiliary variables are accepted, and the
rule is applied to each loading as it stands, which the one factor
models of the study support by analogy rather than directly.

**Which version of each fit row.** With complete continuous data and
`estimator = "ML"` the fit rows are the standard ones, under FIML (the
default) as under listwise deletion: with no value missing the FIML
corrected indices described below equal them, and no line under the
table names a version. With a robust or categorical estimator,
`chi_square` is the scaled statistic and `cfi`, `tli`, `nnfi`, `rmsea`,
and the RMSEA interval are lavaan's robust versions, the corrections
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
describes. One version holds all of those rows. For the older mean and
variance adjusted statistic (`test = "mean.var.adjusted"` passed through
`...`), lavaan computes no robust versions, and they are computed here
as
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
computes them: Savalei (2018) shows that the robust RMSEA and CFI do not
depend on which corrected statistic accompanies them, so they equal
those of the Satorra-Bentler statistic of the same fit. Where lavaan did
not compute the robust indices of a maximum likelihood fit of continuous
items with complete data, the rows hold the standard indices, which
estimate the population values of the ML indices with the small sample
bias corrected as if the data were normal, and not the scaled ones,
which put the scaled statistic into the standard formulas and estimate
population values that change with the kurtosis of the data (Savalei,
2018); after FIML or two-stage estimation, and for ordered items or a
least squares estimator, they hold its scaled indices, or else its
standard ones. A model with no degrees of freedom (one factor of three
items, say) fits exactly in every version: its `cfi` is 1, the model
term \\T - c \cdot df\\ of every version being zero, and its `tli`,
`nnfi`, and `rmsea`, which divide by the degrees of freedom and are 0/0,
take lavaan's conventions for a model that fits exactly, 1 and 0,
whatever the version, and no line under the table names a version for
them, under FIML and after two-stage estimation as well (the `srmr` line
and the two-stage lines on the test, the criteria, and the cases are
kept). A model with no degrees of freedom whose statistic is not zero
does not fit exactly, as when an inequality bound binds (lavaan's count
of the degrees of freedom leaves it out): its `cfi` is formed from its
statistic in the version its CFI alone decides, and its `tli`, `nnfi`,
`rmsea`, and RMSEA limits, which divide by the degrees of freedom, are
undefined and `NA`, as
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
reports them, with a line under the table that says so. One factor of
three variables correlated .8, .8, and .5, with `bounds = "pos.var"`
passed through `...`, has a statistic of 32.8 on 0 degrees of freedom
(\\N\\ = 200) and a `cfi` of .925; its `tli`, `nnfi`, `rmsea`, and RMSEA
limits are `NA`, where lavaan reports the values of exact fit, 1 and 0,
whenever the degrees of freedom are zero. Under FIML (`missing = "ml"`
or `"fiml"`) the `chi_square` row is the FIML likelihood ratio statistic
(its scaled version with `estimator = "MLR"`), and those rows are
lavaan's robust versions, which for a FIML fit are the corrected indices
of Zhang and Savalei (2023): they evaluate the complete data fit
function at the FIML estimates, with the saturated FIML estimates in
place of the sample moments, and so approximate the indices complete
data would give, where the usual indices computed from the `chi_square`
row change with the amount of missing data. lavaan uses their version 6,
which assumes normality, and their version 3, which does not, with
`estimator = "MLR"`. Zhang and Savalei define only the RMSEA and the
CFI; the TLI, the NNFI, and the RMSEA interval are lavaan's extension of
the same correction. lavaan cannot combine FIML with the Satorra-Bentler
or the mean and variance adjusted statistic (`test = "satorra.bentler"`
or `"mean.var.adjusted"` passed through `...`); an explicit FIML request
with one of them stops rather than returning a fit to the complete
cases, and the default fits the items as they stand when none is missing
(see `missing`). The `srmr` row is the SRMR of the covariance residuals
alone (lavaan's `srmr_bentler_nomean` when the model has a mean
structure), as
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
reports it, so that the intercepts a FIML fit adds, whose residuals are
zero, do not lower it (lavaan's `srmr` of a model with a mean structure
averages the mean residuals in). With auxiliary variables (lavaan's
`aux` under FIML) it is over the items alone, for the reason
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
gives. A line under the printed table names the version of the indices
after FIML with incomplete data, with `estimator = "MLR"`, and after
two-stage estimation, when they fall back from the robust ones, and when
DMAR computed them.

Common row names under `term`: `lambda_<factor>_<j>` (loadings;
`lambda_<factor>` when equated), `psi_<factor>_<j>` (error variances;
`psi_<factor>` when equated), `nu_<factor>_<j>` (intercepts, when the
mean structure is modeled; `nu_<factor>` when equated), `phi_<factor>`
(factor variance, fixed to 1), `phi_<factor1>_<factor2>` (latent
correlation), the per-factor defined parameters (`loading_sum_<factor>`,
`error_sum_<factor>`, `omega_<factor>`, `ave_<factor>`, `H_<factor>`),
and the fit rows `chi_square`, `df`, `p_chi_square`, `cfi`, `tli`,
`nnfi`, `rmsea`, `rmsea_ci_lower`, `rmsea_ci_upper`, `rmsea_ci_level`,
`srmr`, `AIC`, `BIC`, `H0`, `H1`.

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

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
*Topics in applied multivariate analysis* (pp. 72–141). Cambridge, UK:
Cambridge University Press.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*, 62–83.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Muthén, B., du Toit, S. H. C., & Spisic, D. (1997). *Robust inference
using weighted least squares and quadratic estimating equations in
latent variable modeling with categorical and continuous outcomes*.
Unpublished technical report.

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
standard errors in covariance structure analysis. In A. von Eye & C. C.
Clogg (Eds.), *Latent variables analysis: Applications for developmental
research* (pp. 399–419). Thousand Oaks, CA: Sage.

Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
missing information from full information maximum likelihood. *Frontiers
in Psychology, 12*, 667802.
[doi:10.3389/fpsyg.2021.667802](https://doi.org/10.3389/fpsyg.2021.667802)

Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation
models with unobservable variables and measurement error. *Journal of
Marketing Research, 18*(1), 39–50.

Graham, J. M. (2006). Congeneric and (essentially) tau-equivalent
estimates of score reliability: What they are and how to use them.
*Educational and Psychological Measurement, 66*(6), 930–944.
[doi:10.1177/0013164406288165](https://doi.org/10.1177/0013164406288165)

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*(1), 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for
assessing discriminant validity in variance-based structural equation
modeling. *Journal of the Academy of Marketing Science, 43*(1), 115–135.
[doi:10.1007/s11747-014-0403-8](https://doi.org/10.1007/s11747-014-0403-8)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Lai, K. (2021). Using information criteria under missing data: Full
information maximum likelihood versus two-stage estimation. *Structural
Equation Modeling, 28*(2), 278–291.
[doi:10.1080/10705511.2020.1780925](https://doi.org/10.1080/10705511.2020.1780925)

Lord, F. M., & Novick, M. R. (1968). *Statistical theories of mental
test scores*. Addison-Wesley.

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

McDonald, R. P. (1999). *Test theory: A unified treatment*. Erlbaum.

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. [doi:10.1007/BF02296192](https://doi.org/10.1007/BF02296192)

Savalei, V. (2010). Expected versus observed information in SEM with
incomplete normal and nonnormal data. *Psychological Methods, 15*(4),
352–367. [doi:10.1037/a0020143](https://doi.org/10.1037/a0020143)

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

von Hippel, P. T. (2016). New confidence intervals and bias comparisons
show that maximum likelihood can beat multiple imputation in small
samples. *Structural Equation Modeling, 23*(3), 422–437.
[doi:10.1080/10705511.2015.1047931](https://doi.org/10.1080/10705511.2015.1047931)

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
*Sociological Methodology, 30*, 165–200.

Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI
following FIML and TS estimation with missing data. *Psychological
Methods, 28*(2), 263–283.
[doi:10.1037/met0000445](https://doi.org/10.1037/met0000445)

## See also

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) and
[`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md) for the
one and two factor convenience wrappers;
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
for the bifactor model, a general factor beside group factors, and
[`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md)
for second-order factors over the factors, both fit by the same engine;
[`plot_fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/plot_fa_confirmatory.md)
to display the estimates and the equality question visually;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
and [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) for the
measurement properties as standalone functions;
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
for the across-group analog of these within-factor constraints;
[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html),
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html).

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
data(holzinger_swineford)
hs_factors <- list(
  verbal    = c("t6_paragraph_comprehension", "t7_sentence",
                "t9_word_meaning"),
  deduction = c("t20_deduction", "t22_problem_reasoning",
                "t23_series_completion"))

# Congeneric measurement model for both factors (the default:
# nothing is constrained, and the header names the structure).
fa_confirmatory(holzinger_swineford, hs_factors)
#> Measurement structure, per factor:
#>   verbal: congeneric (no equality constraints)
#>   deduction: congeneric (no equality constraints)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal_1      
#>  verbal =~ t7_sentence                                    lambda_verbal_2      
#>  verbal =~ t9_word_meaning                                lambda_verbal_3      
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal_1         
#>  t7_sentence ~~ t7_sentence                               psi_verbal_2         
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal_3         
#>  t6_paragraph_comprehension ~1                            nu_verbal_1          
#>  t7_sentence ~1                                           nu_verbal_2          
#>  t9_word_meaning ~1                                       nu_verbal_3          
#>  deduction =~ t20_deduction                               lambda_deduction_1   
#>  deduction =~ t22_problem_reasoning                       lambda_deduction_2   
#>  deduction =~ t23_series_completion                       lambda_deduction_3   
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction_2      
#>  t23_series_completion ~~ t23_series_completion           psi_deduction_3      
#>  t20_deduction ~1                                         nu_deduction_1       
#>  t22_problem_reasoning ~1                                 nu_deduction_2       
#>  t23_series_completion ~1                                 nu_deduction_3       
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
#>                                                           chi_square           
#>                                                           df                   
#>                                                           p_chi_square         
#>                                                           cfi                  
#>                                                           tli                  
#>                                                           nnfi                 
#>                                                           rmsea                
#>                                                           rmsea_ci_lower       
#>                                                           rmsea_ci_upper       
#>                                                           rmsea_ci_level       
#>                                                           srmr                 
#>                                                           AIC                  
#>                                                           BIC                  
#>                                                           H0                   
#>                                                           H1                   
#>  estimate  se     test_statistic p_value  ci_lower ci_upper
#>  2.947     0.169  17.387         < 0.0001 2.614    3.279   
#>  4.389     0.249  17.608         < 0.0001 3.900    4.877   
#>  6.489     0.373  17.417         < 0.0001 5.759    7.220   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  3.473     0.419  8.280          < 0.0001 2.651    4.296   
#>  7.295     0.899  8.117          < 0.0001 5.533    9.056   
#>  16.509    2.043  8.082          < 0.0001 12.506   20.512  
#>  9.183     0.201  45.694         < 0.0001 8.789    9.577   
#>  17.362    0.297  58.452         < 0.0001 16.780   17.944  
#>  15.299    0.441  34.667         < 0.0001 14.434   16.164  
#>  11.149    1.152  9.678          < 0.0001 8.891    13.407  
#>  6.738     0.525  12.844         < 0.0001 5.710    7.767   
#>  6.681     0.524  12.751         < 0.0001 5.654    7.708   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  248.254   23.592 10.523         < 0.0001 202.015  294.493 
#>  38.903    4.782  8.135          < 0.0001 29.529   48.276  
#>  38.626    4.773  8.093          < 0.0001 29.272   47.981  
#>  26.890    1.113  24.170         < 0.0001 24.710   29.071  
#>  26.239    0.529  49.579         < 0.0001 25.202   27.276  
#>  18.136    0.526  34.483         < 0.0001 17.105   19.167  
#>  0.730     0.043  16.885         < 0.0001 0.645    0.815   
#>  13.825    0.647  21.376         < 0.0001 12.557   15.092  
#>  27.277    1.998  13.650         < 0.0001 23.360   31.194  
#>  0.875     0.014  64.416         < 0.0001 0.846    0.899   
#>  0.719     0.023  31.557         < 0.0001 0.673    0.762   
#>  0.885     0.012  76.912         < 0.0001 0.860    0.906   
#>  24.568    1.548  15.867         < 0.0001 21.534   27.603  
#>  325.783   23.890 13.637         < 0.0001 278.959  372.606 
#>  0.649     0.036  18.221         < 0.0001 0.577    0.716   
#>  0.469     0.033  14.095         < 0.0001 0.405    0.535   
#>  0.738     0.027  27.128         < 0.0001 0.682    0.788   
#>  13.769    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  8         <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.0880    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.993     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.987     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.987     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.049     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.000     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.092     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.900     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.023     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11767.495 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11837.931 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5864.748 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>           <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)

# Equal loadings: equal_loading = TRUE equates the loadings within
# every factor, and a descriptor can also differ by factor, as here,
# where the verbal loadings are equated and the deduction loadings are
# left free. Because only the covariance structure identifies this
# constraint, the implied structure of the verbal factor is
# essentially tau-equivalent, and the header names each factor's
# structure separately.
fa_confirmatory(holzinger_swineford, hs_factors,
      equal_loading = c(verbal = TRUE, deduction = FALSE))
#> Measurement structure, per factor:
#>   verbal: essentially tau-equivalent (equal loadings)
#>   deduction: congeneric (no equality constraints)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal        
#>  verbal =~ t7_sentence                                    lambda_verbal        
#>  verbal =~ t9_word_meaning                                lambda_verbal        
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal_1         
#>  t7_sentence ~~ t7_sentence                               psi_verbal_2         
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal_3         
#>  t6_paragraph_comprehension ~1                            nu_verbal_1          
#>  t7_sentence ~1                                           nu_verbal_2          
#>  t9_word_meaning ~1                                       nu_verbal_3          
#>  deduction =~ t20_deduction                               lambda_deduction_1   
#>  deduction =~ t22_problem_reasoning                       lambda_deduction_2   
#>  deduction =~ t23_series_completion                       lambda_deduction_3   
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction_2      
#>  t23_series_completion ~~ t23_series_completion           psi_deduction_3      
#>  t20_deduction ~1                                         nu_deduction_1       
#>  t22_problem_reasoning ~1                                 nu_deduction_2       
#>  t23_series_completion ~1                                 nu_deduction_3       
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
#>                                                           chi_square           
#>                                                           df                   
#>                                                           p_chi_square         
#>                                                           cfi                  
#>                                                           tli                  
#>                                                           nnfi                 
#>                                                           rmsea                
#>                                                           rmsea_ci_lower       
#>                                                           rmsea_ci_upper       
#>                                                           rmsea_ci_level       
#>                                                           srmr                 
#>                                                           AIC                  
#>                                                           BIC                  
#>                                                           H0                   
#>                                                           H1                   
#>  estimate  se     test_statistic p_value  ci_lower ci_upper
#>  3.458     0.154  22.478         < 0.0001 3.157    3.760   
#>  3.458     0.154  22.478         < 0.0001 3.157    3.760   
#>  3.458     0.154  22.478         < 0.0001 3.157    3.760   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  1.220     0.533  2.290          0.0220   0.176    2.264   
#>  10.834    1.096  9.882          < 0.0001 8.685    12.983  
#>  31.301    2.711  11.545         < 0.0001 25.987   36.614  
#>  9.183     0.209  43.886         < 0.0001 8.773    9.593   
#>  17.362    0.275  63.094         < 0.0001 16.823   17.901  
#>  15.299    0.379  40.356         < 0.0001 14.556   16.042  
#>  11.181    1.157  9.667          < 0.0001 8.914    13.448  
#>  6.570     0.535  12.277         < 0.0001 5.521    7.619   
#>  6.857     0.530  12.942         < 0.0001 5.818    7.895   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  247.540   23.686 10.451         < 0.0001 201.117  293.963 
#>  41.145    4.971  8.277          < 0.0001 31.402   50.888  
#>  36.245    4.918  7.371          < 0.0001 26.607   45.883  
#>  26.890    1.113  24.170         < 0.0001 24.710   29.071  
#>  26.239    0.529  49.579         < 0.0001 25.202   27.276  
#>  18.136    0.526  34.483         < 0.0001 17.105   19.167  
#>  0.667     0.051  13.176         < 0.0001 0.568    0.766   
#>  10.374    0.462  22.478         < 0.0001 9.470    11.279  
#>  43.354    2.912  14.888         < 0.0001 37.647   49.062  
#>  0.713     0.024  30.295         < 0.0001 0.665    0.757   
#>  0.570     0.020  28.521         < 0.0001 0.530    0.608   
#>  0.919     0.028  32.237         < 0.0001 0.842    0.960   
#>  24.608    1.548  15.895         < 0.0001 21.573   27.642  
#>  324.930   23.906 13.592         < 0.0001 278.076  371.785 
#>  0.651     0.036  18.288         < 0.0001 0.578    0.717   
#>  0.471     0.033  14.167         < 0.0001 0.406    0.536   
#>  0.740     0.027  26.955         < 0.0001 0.683    0.790   
#>  134.162   <NA>   <NA>           <NA>     <NA>     <NA>    
#>  10        <NA>   <NA>           <NA>     <NA>     <NA>    
#>  < 0.0001  <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.846     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.769     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.769     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.203     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.173     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.234     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.900     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.164     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11883.888 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11946.909 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5924.944 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>           <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)

# Measurement properties: omega, ave, and H per factor (each with a
# delta method standard error and confidence interval), the latent
# correlations, and the htmt ratios.
fa_confirmatory(holzinger_swineford, hs_factors, output = "measurement")
#> Measurement structure, per factor:
#>   verbal: congeneric (no equality constraints)
#>   deduction: congeneric (no equality constraints)
#> 
#>  syntax              term                  estimate se    test_statistic
#>  verbal ~~ deduction phi_verbal_deduction  0.730    0.043 16.885        
#>                      omega_verbal          0.875    0.014 64.416        
#>                      ave_verbal            0.719    0.023 31.557        
#>                      H_verbal              0.885    0.012 76.912        
#>                      omega_deduction       0.649    0.036 18.221        
#>                      ave_deduction         0.469    0.033 14.095        
#>                      H_deduction           0.738    0.027 27.128        
#>  verbal ~~ deduction htmt_verbal_deduction 0.731    <NA>  <NA>          
#>  p_value  ci_lower ci_upper
#>  < 0.0001 0.645    0.815   
#>  < 0.0001 0.846    0.899   
#>  < 0.0001 0.673    0.762   
#>  < 0.0001 0.860    0.906   
#>  < 0.0001 0.577    0.716   
#>  < 0.0001 0.405    0.535   
#>  < 0.0001 0.682    0.788   
#>  <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95%

# Equal loadings and intercepts make the items tau-equivalent, and
# equal error variances as well make them parallel, the structure
# the header names here. The mean structure is added because
# equal_intercept asks about it, so the nu terms join the table.
fa_confirmatory(holzinger_swineford, hs_factors, equal_loading = TRUE,
      equal_intercept = TRUE, equal_error = TRUE)
#> Measurement structure, per factor:
#>   verbal: parallel (equal loadings, intercepts, and error variances)
#>   deduction: parallel (equal loadings, intercepts, and error variances)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal        
#>  verbal =~ t7_sentence                                    lambda_verbal        
#>  verbal =~ t9_word_meaning                                lambda_verbal        
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal           
#>  t7_sentence ~~ t7_sentence                               psi_verbal           
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal           
#>  t6_paragraph_comprehension ~1                            nu_verbal            
#>  t7_sentence ~1                                           nu_verbal            
#>  t9_word_meaning ~1                                       nu_verbal            
#>  deduction =~ t20_deduction                               lambda_deduction     
#>  deduction =~ t22_problem_reasoning                       lambda_deduction     
#>  deduction =~ t23_series_completion                       lambda_deduction     
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction        
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction        
#>  t23_series_completion ~~ t23_series_completion           psi_deduction        
#>  t20_deduction ~1                                         nu_deduction         
#>  t22_problem_reasoning ~1                                 nu_deduction         
#>  t23_series_completion ~1                                 nu_deduction         
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
#>                                                           chi_square           
#>                                                           df                   
#>                                                           p_chi_square         
#>                                                           cfi                  
#>                                                           tli                  
#>                                                           nnfi                 
#>                                                           rmsea                
#>                                                           rmsea_ci_lower       
#>                                                           rmsea_ci_upper       
#>                                                           rmsea_ci_level       
#>                                                           srmr                 
#>                                                           AIC                  
#>                                                           BIC                  
#>                                                           H0                   
#>                                                           H1                   
#>  estimate  se     test_statistic p_value  ci_lower ci_upper
#>  3.755     0.274  13.686         < 0.0001 3.217    4.292   
#>  3.755     0.274  13.686         < 0.0001 3.217    4.292   
#>  3.755     0.274  13.686         < 0.0001 3.217    4.292   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  30.411    1.753  17.349         < 0.0001 26.975   33.846  
#>  30.411    1.753  17.349         < 0.0001 26.975   33.846  
#>  30.411    1.753  17.349         < 0.0001 26.975   33.846  
#>  13.948    0.284  49.157         < 0.0001 13.392   14.504  
#>  13.948    0.284  49.157         < 0.0001 13.392   14.504  
#>  13.948    0.284  49.157         < 0.0001 13.392   14.504  
#>  7.559     0.585  12.928         < 0.0001 6.413    8.705   
#>  7.559     0.585  12.928         < 0.0001 6.413    8.705   
#>  7.559     0.585  12.928         < 0.0001 6.413    8.705   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  138.761   7.998  17.349         < 0.0001 123.085  154.437 
#>  138.761   7.998  17.349         < 0.0001 123.085  154.437 
#>  138.761   7.998  17.349         < 0.0001 123.085  154.437 
#>  23.755    0.586  40.532         < 0.0001 22.607   24.904  
#>  23.755    0.586  40.532         < 0.0001 22.607   24.904  
#>  23.755    0.586  40.532         < 0.0001 22.607   24.904  
#>  0.979     0.070  13.915         < 0.0001 0.841    1.117   
#>  11.264    0.823  13.686         < 0.0001 9.651    12.877  
#>  91.233    5.259  17.349         < 0.0001 80.926   101.539 
#>  0.582     0.042  13.930         < 0.0001 0.498    0.661   
#>  0.317     0.037  8.528          < 0.0001 0.249    0.394   
#>  0.582     0.042  13.930         < 0.0001 0.498    0.661   
#>  22.677    1.754  12.928         < 0.0001 19.239   26.115  
#>  416.282   23.994 17.349         < 0.0001 369.255  463.310 
#>  0.553     0.045  12.374         < 0.0001 0.464    0.638   
#>  0.292     0.037  7.815          < 0.0001 0.224    0.370   
#>  0.553     0.045  12.374         < 0.0001 0.464    0.638   
#>  1340.898  <NA>   <NA>           <NA>     <NA>     <NA>    
#>  20        <NA>   <NA>           <NA>     <NA>     <NA>    
#>  < 0.0001  <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.000     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -0.228    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -0.228    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.468     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.447     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.490     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.900     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.759     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  13070.624 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  13096.574 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -6528.312 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5857.863 <NA>   <NA>           <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)

# Ordered categorical items: two correlated factors of three
# four-category items each, generated here. The model is fit by WLSMV
# to polychoric correlations, and each ordered factor's omega is
# reported on the categorical sum score metric (categorical omega),
# which is why the interval columns of those two rows are NA.
set.seed(113)
n <- 200
eta_a <- rnorm(n)
eta_b <- 0.5 * eta_a + sqrt(1 - 0.5^2) * rnorm(n)
lambda <- c(0.6, 0.7, 0.8)
lat <- cbind(outer(eta_a, lambda), outer(eta_b, lambda)) +
  matrix(rnorm(n * 6), n, 6) %*% diag(sqrt(1 - c(lambda, lambda)^2))
likert <- as.data.frame(apply(lat, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1, 0, 1, Inf)))))
names(likert) <- paste0("item_", 1:6)
fa_confirmatory(likert,
      list(scale_a = paste0("item_", 1:3),
           scale_b = paste0("item_", 4:6)),
      ordered = TRUE, output = "measurement")
#> Ordered items declared: switching to estimator = "WLSMV" with robust standard errors.
#> omega_scale_a, omega_scale_b reported on the categorical sum score metric (categorical omega; Green & Yang, 2009; Kelley & Pornprasertmanit, 2016).
#> Measurement structure, per factor:
#>   scale_a: congeneric (no equality constraints)
#>   scale_b: congeneric (no equality constraints)
#> 
#>  syntax             term                 estimate se    test_statistic p_value 
#>  scale_a ~~ scale_b phi_scale_a_scale_b  0.427    0.084 5.089          < 0.0001
#>                     omega_scale_a        0.688    <NA>  <NA>           <NA>    
#>                     ave_scale_a          0.491    0.047 10.511         < 0.0001
#>                     H_scale_a            0.793    0.058 13.695         < 0.0001
#>                     omega_scale_b        0.762    <NA>  <NA>           <NA>    
#>                     ave_scale_b          0.576    0.040 14.362         < 0.0001
#>                     H_scale_b            0.843    0.045 18.788         < 0.0001
#>  scale_a ~~ scale_b htmt_scale_a_scale_b 0.412    <NA>  <NA>           <NA>    
#>  ci_lower ci_upper
#>  0.262    0.591   
#>  <NA>     <NA>    
#>  0.401    0.582   
#>  0.657    0.884   
#>  <NA>     <NA>    
#>  0.496    0.652   
#>  0.735    0.913   
#>  <NA>     <NA>    
#> 
#> Confidence level: 95%

# Two-stage maximum likelihood with incomplete data, here 60 scores on
# one verbal test deleted at random. The chi_square row is the residual
# based statistic of Savalei and Bentler (2009), the cfi and rmsea rows
# are the normal theory analog of the two-stage indices of Zhang and
# Savalei (2023), and the lines under the table say which statistic,
# indices, and criteria are reported.
hs_incomplete <- holzinger_swineford
set.seed(113)
hs_incomplete$t7_sentence[sample(301, 60)] <- NA
fa_confirmatory(hs_incomplete, hs_factors, missing = "two_stage")
#> Measurement structure, per factor:
#>   verbal: congeneric (no equality constraints)
#>   deduction: congeneric (no equality constraints)
#> 
#>  syntax                                                   term                 
#>  verbal =~ t6_paragraph_comprehension                     lambda_verbal_1      
#>  verbal =~ t7_sentence                                    lambda_verbal_2      
#>  verbal =~ t9_word_meaning                                lambda_verbal_3      
#>  verbal ~~ verbal                                         phi_verbal           
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension psi_verbal_1         
#>  t7_sentence ~~ t7_sentence                               psi_verbal_2         
#>  t9_word_meaning ~~ t9_word_meaning                       psi_verbal_3         
#>  t6_paragraph_comprehension ~1                            nu_verbal_1          
#>  t7_sentence ~1                                           nu_verbal_2          
#>  t9_word_meaning ~1                                       nu_verbal_3          
#>  deduction =~ t20_deduction                               lambda_deduction_1   
#>  deduction =~ t22_problem_reasoning                       lambda_deduction_2   
#>  deduction =~ t23_series_completion                       lambda_deduction_3   
#>  deduction ~~ deduction                                   phi_deduction        
#>  t20_deduction ~~ t20_deduction                           psi_deduction_1      
#>  t22_problem_reasoning ~~ t22_problem_reasoning           psi_deduction_2      
#>  t23_series_completion ~~ t23_series_completion           psi_deduction_3      
#>  t20_deduction ~1                                         nu_deduction_1       
#>  t22_problem_reasoning ~1                                 nu_deduction_2       
#>  t23_series_completion ~1                                 nu_deduction_3       
#>  verbal ~~ deduction                                      phi_verbal_deduction 
#>                                                           loading_sum_verbal   
#>                                                           error_sum_verbal     
#>                                                           omega_verbal         
#>                                                           ave_verbal           
#>                                                           H_verbal             
#>                                                           loading_sum_deduction
#>                                                           error_sum_deduction  
#>                                                           omega_deduction      
#>                                                           ave_deduction        
#>                                                           H_deduction          
#>                                                           chi_square           
#>                                                           df                   
#>                                                           p_chi_square         
#>                                                           cfi                  
#>                                                           tli                  
#>                                                           nnfi                 
#>                                                           rmsea                
#>                                                           rmsea_ci_lower       
#>                                                           rmsea_ci_upper       
#>                                                           rmsea_ci_level       
#>                                                           srmr                 
#>                                                           AIC                  
#>                                                           BIC                  
#>                                                           H0                   
#>                                                           H1                   
#>  estimate  se     test_statistic p_value  ci_lower ci_upper
#>  2.959     0.170  17.375         < 0.0001 2.625    3.293   
#>  4.360     0.266  16.403         < 0.0001 3.839    4.881   
#>  6.488     0.373  17.393         < 0.0001 5.757    7.219   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  3.399     0.434  7.829          < 0.0001 2.548    4.250   
#>  7.037     0.954  7.379          < 0.0001 5.168    8.906   
#>  16.527    2.025  8.163          < 0.0001 12.559   20.496  
#>  9.183     0.201  45.694         < 0.0001 8.789    9.577   
#>  17.349    0.307  56.450         < 0.0001 16.747   17.951  
#>  15.299    0.441  34.667         < 0.0001 14.434   16.164  
#>  11.169    1.147  9.737          < 0.0001 8.921    13.418  
#>  6.732     0.525  12.816         < 0.0001 5.702    7.761   
#>  6.681     0.519  12.866         < 0.0001 5.663    7.698   
#>  1.000     0.000  <NA>           <NA>     1.000    1.000   
#>  247.801   23.418 10.581         < 0.0001 201.901  293.700 
#>  38.992    4.787  8.145          < 0.0001 29.609   48.375  
#>  38.628    4.602  8.394          < 0.0001 29.609   47.648  
#>  26.890    1.113  24.170         < 0.0001 24.710   29.071  
#>  26.239    0.529  49.579         < 0.0001 25.202   27.276  
#>  18.136    0.526  34.483         < 0.0001 17.105   19.167  
#>  0.721     0.043  16.717         < 0.0001 0.637    0.806   
#>  13.807    0.652  21.178         < 0.0001 12.530   15.085  
#>  26.963    2.037  13.234         < 0.0001 22.970   30.957  
#>  0.876     0.014  63.527         < 0.0001 0.846    0.901   
#>  0.723     0.023  31.275         < 0.0001 0.675    0.766   
#>  0.887     0.012  76.332         < 0.0001 0.862    0.908   
#>  24.582    1.538  15.982         < 0.0001 21.567   27.596  
#>  325.421   23.625 13.774         < 0.0001 279.117  371.725 
#>  0.650     0.035  18.431         < 0.0001 0.578    0.716   
#>  0.469     0.033  14.243         < 0.0001 0.406    0.534   
#>  0.738     0.027  27.350         < 0.0001 0.682    0.788   
#>  14.301    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  8         <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.0742    <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.990     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.981     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.981     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.059     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.005     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.102     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.900     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  0.026     <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11758.614 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  11829.049 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5860.307 <NA>   <NA>           <NA>     <NA>     <NA>    
#>  -5851.815 <NA>   <NA>           <NA>     <NA>     <NA>    
#> 
#> Confidence level: 95% (90% for the RMSEA interval)
#> 
#> Two-stage maximum likelihood: the model fit to the EM estimates of the means and covariances, with normal theory standard errors (Savalei & Bentler, 2009)
#> 
#> chi_square: the residual based statistic (Savalei & Bentler, 2009)
#> 
#> cfi, rmsea: the normal theory analog of the TS V2 indices of Zhang and Savalei (2023); tli, nnfi, and the RMSEA interval: lavaan's extension of the same correction
#> 
#> AIC, BIC, H0, H1: the complete data normal log likelihood at the EM estimates; AIC and BIC are the two-stage criteria of Lai (2021), comparable only with two-stage fits of the same variables, and -2(H0 - H1) is not a chi square statistic

# Does the equal-loadings description hold? Two fits that differ only
# in a descriptor are nested, so output = "fit" hands them straight to
# lavaan's likelihood ratio test.
fit_free  <- fa_confirmatory(holzinger_swineford, hs_factors, output = "fit")
fit_equal <- fa_confirmatory(holzinger_swineford, hs_factors,
                   equal_loading = TRUE, output = "fit")
lavaan::lavTestLRT(fit_free, fit_equal)
#> 
#> Chi-Squared Difference Test
#> 
#>           Df   AIC   BIC   Chisq Chisq diff Df diff Pr(>Chisq)    
#> fit_free   8 11768 11838  13.769                                  
#> fit_equal 12 11895 11950 149.150     135.38       4  < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
