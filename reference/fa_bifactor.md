# Bifactor Confirmatory Factor Analysis Model

Fits the bifactor model, in which a general factor accounts for what
every item of a scale has in common and two or more group factors
account for what the items of each subscale share beyond it, with all of
the factors uncorrelated (Holzinger & Swineford, 1937; Reise, 2012). The
model is described as in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
by naming the items of each group factor, and is fit by maximum
likelihood to raw data (by full information maximum likelihood, the
default, or after listwise deletion) or to a covariance matrix. The
table reports the loadings, error variances, and intercepts with their
confidence intervals; the indices that describe how much of the total
score and of each subscale score the general factor accounts for, each
with a delta method standard error and a confidence interval (Rodriguez,
Reise, & Haviland, 2016): omega and omega hierarchical (omega_H) for the
total score, omega and omega hierarchical subscale (omega_HS) for each
subscale score, the explained common variance (ECV), coefficient *H*,
and the maximal reliability of each factor, and the percentage of
uncontaminated correlations (PUC); and the fit of the model.

## Usage

``` r
fa_bifactor(
  data = NULL,
  factors,
  S = NULL,
  N = NULL,
  M = NULL,
  general = "general",
  general_only = NULL,
  equal_loading = FALSE,
  equal_intercept = FALSE,
  equal_error = FALSE,
  meanstructure = NULL,
  estimator = "ML",
  missing = "fiml",
  se = NULL,
  conf_level = 0.95,
  rmsea_conf_level = 0.9,
  output = c("verbose", "measurement", "summary", "standardized", "fit"),
  ...,
  wald_df = Inf
)

# S3 method for class 'dmar_fa_bifactor'
tidy(x, ...)

# S3 method for class 'dmar_fa_bifactor'
glance(x, ...)
```

## Arguments

- data:

  A raw data matrix or data frame, rows are respondents and columns
  include the items named in `factors`. Supply exactly one of `data` or
  `S`.

- factors:

  Named list of the group factors. Each element names a group factor and
  gives the character vector of its items; each item is listed under one
  group factor, and the general factor loads on all of them. Two or more
  group factors are needed, each with three or more items (two when
  their loadings are equated with `equal_loading`); see *Details* for
  the structures that are refused. An item, here or in `general_only`,
  that is not a column of `data` or an item of `S` is refused with an
  error of class `dmar_items_invalid` (see
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

- general:

  The name of the general factor, `"general"` by default. It names the
  general factor's rows (`lambda_general_1`, `omega_H_general`).

- general_only:

  `NULL` (the default), or a character vector of items on the general
  factor alone, with no group factor: an item written for the construct
  as a whole, say, or the items of a subscale too small to carry a group
  factor of its own. They count in the general factor's rows (its
  loadings, omega, omega_H, ECV, *H*, maximal reliability, and PUC) and
  in no subscale.

- equal_loading:

  Logical, or a named logical vector with one element for the general
  factor and one for each group factor. `TRUE` constrains the loadings
  within a factor to a single value. Defaults to `FALSE`.

- equal_intercept:

  Logical, or a named logical vector with one element per group factor.
  `TRUE` constrains the intercepts of a group factor's items to a single
  value; this requires the mean structure. The items on the general
  factor alone keep their own. Defaults to `FALSE`.

- equal_error:

  Logical, or a named logical vector with one element per group factor.
  `TRUE` constrains the error variances of a group factor's items to a
  single value; the items on the general factor alone keep their own.
  Defaults to `FALSE`.

- meanstructure:

  Logical or `NULL`. `NULL` (default) models the mean structure exactly
  when it is needed: under FIML (the default with raw data, so the item
  intercepts, the `nu` rows, join the table), when any `equal_intercept`
  is `TRUE`, or when `M` is supplied. `TRUE` models the intercepts
  regardless (raw data or `M` required); `FALSE` suppresses them, which
  FIML does not allow.

- estimator:

  `"ML"` (default), maximum likelihood, or `"MLR"`, maximum likelihood
  with standard errors and a test statistic robust to nonnormality (raw
  data required). Ordered items are not available; see *Details*.

- missing:

  The missing data treatment for raw data: `"fiml"` (default; full
  information maximum likelihood, which uses every case with an observed
  item; lavaan's name `"ml"` is accepted) or `"listwise"`. The default
  is resolved as in
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md):
  with `se = "robust.sem"` or a test statistic FIML cannot compute, the
  items are fit as they stand when no item value is missing, and
  otherwise the call stops and names the option to choose. Ignored with
  `S`.

- se:

  Standard error method passed to lavaan. The default, `NULL`, follows
  the estimator: `"standard"` for ML and `"robust.huber.white"` for MLR.
  `"robust.sem"` requires complete data and is refused with FIML, and
  `"none"` requests point estimates only.

- conf_level:

  Confidence level of the intervals for the parameters and the indices.
  Defaults to 0.95. It does not reach the RMSEA interval, whose level is
  `rmsea_conf_level`.

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

- output:

  Format of the returned object:

  `"verbose"`

  :   (default) The loadings, factor variances (fixed at 1), error
      variances, and intercepts, then the indices, then the fit rows.

  `"measurement"`

  :   The indices alone.

  `"summary"`

  :   The [`summary()`](https://rdrr.io/r/base/summary.html) of the
      lavaan fit (not a data frame).

  `"standardized"`

  :   [`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
      of the fit; see *Details*. Under `se = "none"` its `se`, `z`,
      `pvalue`, `ci.lower`, and `ci.upper` columns are `NA`.

  `"fit"`

  :   The lavaan fit, which
      [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
      [`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
      [`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
      and
      [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html)
      take.

- ...:

  For `fa_bifactor()`, additional arguments forwarded to
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html) (for example
  `cluster`); a `group` argument is an error, as in
  [`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  and so is `ordered`. Unused by `tidy()` and `glance()`.

- wald_df:

  The reference distribution of the Wald intervals and *p*-values: `Inf`
  (the default), the normal distribution, or a single positive number,
  the degrees of freedom of a *t* distribution for every parameter and
  index row (for the indices, on the logit scale), as in
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).
  The keyword `"missing_information"` is refused: the calibration study
  behind it simulated one factor models, and nothing calibrates it for
  the loadings of a bifactor model or for its indices. It follows `...`,
  so it must be named in full.

- x:

  A `dmar_fa_bifactor` table returned by `fa_bifactor()`.

## Value

For `output = "verbose"` (default) and `"measurement"`, a `data.frame`
(classes `dmar_fa_bifactor`, `dmar_tbl`) with columns `syntax`, `term`,
`estimate`, `se`, `test_statistic`, `p_value`, `ci_lower`, and
`ci_upper`. Its attributes record the model: `general`, `factors`, and
`general_only`; `equal`, a list naming the factors whose loadings, error
variances, and intercepts are equated; `has_means`; `omega_metric`
(`"raw_item_sum"`); `improper`, the negative error variances named by
item (empty for a proper solution); and `wald_df`. The verbose table
also carries
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)'s
`test`, `estimator`, `missing`, `fit_indices`, and `fit_rows`. The
printed header describes the model, and the lines under the table say
what the indices describe. For `"summary"`, the lavaan summary; for
`"standardized"`, the standardized solution with the `"wald_df"`
attribute; for `"fit"`, the lavaan fit. `tidy()` gives the parameter and
index rows and `glance()` one row of the fit rows and what produced
them, as for
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).

## Details

**The model.** An item \\i\\ of group factor \\s\\ is \$\$x_i = \nu_i +
\lambda\_{gi} g + \lambda\_{si} s + e_i,\$\$ and an item on the general
factor alone lacks the group term. Every factor's variance is fixed at 1
(and its mean at 0 when the mean structure is modeled), every pair of
factors is uncorrelated, and the first loading of each factor is freed,
so each loading is the covariance of its item with the factor. The
general factor carries what all of the items share; a group factor
carries only what its subscale's items share beyond that, so its
loadings are smaller than those of a factor of the subscale alone, and
they can be near zero or change sign when the general factor takes up
what the subscale measures. Every loading, error variance, and intercept
carries a label, as in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
and the indices are lavaan defined parameters of those labels, so the
fit itself gives each one a delta method standard error. The
standardized solution is derived from this fit of the covariances (or of
the raw data), with the factor variances fixed at 1, not from a fit to
correlations.

**Scope and identification.** A structure outside the model, or one the
data cannot identify, is refused before anything is fit, with an error
of class `dmar_bifactor_unidentified` (which inherits `dmar_refusal`;
see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md))
that names the factors or items at issue and gives the kind of refusal
as its `reason` element. Fewer than two group factors, and an item
listed under two group factors, are outside the model the function fits
(`reason` `"scope"`): a general factor beside a single group factor,
with items on the general factor alone, can be identified, but omega_HS,
ECV, and PUC describe how the common variance divides among subscales,
which takes two or more. Names placed twice or that coincide are refused
with `reason` `"names"`: an item listed twice under one group factor, or
both under a group factor and in `general_only`; the general factor and
a group factor of one name, a factor and an item of one name, or labels
made from the names that coincide. A group factor of one item, or of two
whose loadings are free (they enter the model only through their
product), and two group factors whose general loadings are fixed only up
to a scale are not identified (`reason` `"unidentified"`). With two
group factors and no item on the general factor alone, the covariances
between the items of the two groups fix the general loadings of one
group only up to a factor \\c\\ and those of the other up to \\1/c\\,
and \\c\\ must come from inside a group, whose item variances and
covariances must outnumber the group loadings and error variances it
estimates. Two groups each of three items (or of two items with equal
loadings) leave \\c\\ undetermined, and so the model is not identified;
an item on the general factor alone, a group of four or more items,
equal loadings or equal error variances within a group of three, or
equal general loadings each identify it. Every structure with exactly
two group factors that is fit gives a warning of class
`dmar_bifactor_empirical`: such a model can be empirically
underidentified, since when the item covariances are nearly equal within
each group and between the groups, different general and group loadings
reproduce them equally well (Green & Yang, 2018). Items on the general
factor alone, or items whose loadings differ in strength, protect
against it. The identification checks run once the other arguments and
the data have passed theirs, so the warning comes only with a fit. A
count of parameters cannot see empirical underidentification, and when
lavaan's information matrix cannot be inverted, or is singular (scaled
to a unit diagonal, its smallest eigenvalue is below \\10^{-6}\\), every
standard error, test, and interval is `NA`, the maximal reliability's
included, with a warning. The singularity is checked, as in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
for the observed information of a FIML fit and for the expected
information of a listwise fit or of a covariance matrix alike, since
lavaan can invert either when it is singular.

**The indices.** Write \\\lambda\_{gi}\\ for the general loading of item
\\i\\, \\\lambda\_{si}\\ for its loading on its group factor \\s\\ (zero
for an item on the general factor alone), \\\psi_i\\ for its error
variance, and \\v_i = \lambda\_{gi}^2 + \lambda\_{si}^2 + \psi_i\\ for
its model implied variance. The variance of the total score, the
unweighted sum of the items, is \$\$\sigma^2_X = (\textstyle\sum_i
\lambda\_{gi})^2 + \sum_s (\sum\_{i \in s} \lambda\_{si})^2 + \sum_i
\psi_i,\$\$ and the rows are

- `omega_<general>`:

  Omega for the total score, \\\[(\sum_i \lambda\_{gi})^2 + \sum_s
  (\sum\_{i \in s} \lambda\_{si})^2\] / \sigma^2_X\\: its reliability
  with every factor counted as true score (McDonald, 1999).

- `omega_H_<general>`:

  Omega hierarchical, \\(\sum_i \lambda\_{gi})^2 / \sigma^2_X\\: the
  share of the total score's variance that is due to the general factor
  alone, the group factors' variance counted as error. Its square root
  is the correlation of the total score with the general factor.

- `omega_<s>`:

  Omega for the subscale score of group factor \\s\\, the sum of its
  items, \\(\sum\_{i \in s} \lambda\_{gi})^2 + (\sum\_{i \in s}
  \lambda\_{si})^2\\ over that score's variance, which adds \\\sum\_{i
  \in s} \psi_i\\: its reliability with both factors counted as true
  score.

- `omega_HS_<s>`:

  Omega hierarchical subscale, \\(\sum\_{i \in s} \lambda\_{si})^2\\
  over the same variance: the share of the subscale score's variance
  that is due to its group factor alone, what the subscale score says
  beyond the total.

- `ECV_<general>`, `ECV_<s>`:

  The explained common variance, formed from the standardized loadings
  \\\lambda\_{gi} / \sqrt{v_i}\\ and \\\lambda\_{si} / \sqrt{v_i}\\, as
  [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
  forms it: for the general factor, its share of the common variance of
  all of the items, \\\sum_i \lambda\_{gi}^2 / v_i\\ over \\\sum_i
  (\lambda\_{gi}^2 + \lambda\_{si}^2) / v_i\\; for a group factor's row,
  the general factor's share of the common variance of that factor's
  items, the same sums over them, so one minus it is the group factor's
  own share.

- `H_<general>`, `H_<s>`:

  Coefficient *H* in the form of Rodriguez, Reise, and Haviland (2016):
  the closed form of Hancock and Mueller (2001) over the factor's
  loadings \\\lambda\_{fi}\\, with every other source of an item's
  variance (the other factors' included) counted as its error, \\1 /
  (1 + 1 / \sum_i \lambda\_{fi}^2 / (v_i - \lambda\_{fi}^2))\\.

- `maximal_reliability_<general>`, `maximal_reliability_<s>`:

  The maximal reliability of the factor under the fitted model,
  \\\lambda_f' \Sigma^{-1} \lambda_f\\, with \\\lambda_f\\ the factor's
  loadings (zero for an item off it) and \\\Sigma\\ the model implied
  covariance matrix of the items: the squared multiple correlation of
  the factor with all of the items, the largest reliability any weighted
  sum of them reaches as a measure of that factor, and the square of the
  determinacy of the factor's regression score. *H* treats the variance
  an item shares with the other factors as uncorrelated error, which in
  a bifactor model it is not, so the two differ, in either direction;
  [`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
  says when.

- `PUC_<general>`:

  The percentage of uncontaminated correlations, as a proportion: the
  share of the pairs of items that do not share a group factor, whose
  correlations reflect the general factor alone. It depends on the
  structure, not on the data, so it has no standard error or interval.

Rodriguez, Reise, and Haviland (2016) describe how each is read. In
brief, omega_H against omega for the total score says how much of the
total score's reliable variance is general; omega_HS against omega for a
subscale says how much of the subscale score's reliable variance is its
own; ECV and PUC describe how close the items are to measuring one
thing; and *H* and the maximal reliability describe how well a factor
could be measured by an optimally weighted composite. The table reports
the numbers with their uncertainty; whether a total or subscale score is
worth reporting is the researcher's judgment, and the package draws no
line.

**The metric of omega.** Omega, omega_H, and omega_HS here are the
reliabilities of unit-weighted sums of the items as they are scored, the
coefficients that
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
and `semTools::compRelSEM()` compute.
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
computes them from the standardized solution, for sums of the
standardized items, as Rodriguez, Reise, and Haviland (2016) write them.
The two agree when the items have equal variances and differ otherwise:
on the example below, whose tests' variances range from 8.0 to 153.3,
omega_H is 0.552 for the sum of the scores and 0.577 for the sum of the
standardized scores, and on the simulated example of
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
whose items have nearly equal variances, 0.67305 against 0.67348. ECV is
formed from the standardized loadings, and *H*, the maximal reliability,
and PUC do not depend on the metric of the items, so those rows equal
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)'s
with its model implied denominator. With `output = "standardized"`, the
defined indices are evaluated at the standardized parameters, so the
omega rows there are those of the standardized items. lavaan evaluates a
defined parameter at the standardized values of the labels it names, and
a label that `equal_loading` or `equal_error` gives several items has
one unstandardized value but a standardized value for each item, which
differ when the items' variances differ; lavaan uses the first item's.
Each index that uses such a label is computed here from every item's own
standardized loadings and error variance, as
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
computes its rows, with a delta method standard error from their
derivatives with respect to the free parameters and lavaan's covariance
matrix of the estimates (under `se = "bootstrap"`, from lavaan's
bootstrap draws); the other indices are lavaan's. ECV and *H* then equal
the rows of the other tables.

**Intervals.** The standard errors of omega, omega_H, omega_HS, ECV, and
*H* are lavaan's delta method standard errors of the defined parameters.
The maximal reliability needs a matrix inverse, which lavaan's defined
parameters cannot express, so its standard error is the delta method's
with the derivatives taken numerically over the free parameters and
lavaan's covariance matrix of their estimates, as
[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
forms those of its determinacies. For every index the limits are formed
on the logit scale from that standard error and transformed back
(Browne, 1982; the confidence interval transformation principle of
Kelley, 2007), as
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
forms the limits of omega and *H*, so they lie within (0, 1); the `se`
column is the standard error of the index itself. The sample values of
*H* and of the maximal reliability overstate their population values on
average, the more so for a weakly determined factor (Aguirre-Urreta,
Rönkkö, & McIntosh, 2019);
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
on the fit of `output = "fit"` gives bootstrap intervals for the indices
of
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
and the bootstrap estimate of that bias.

**Improper solutions.** A bifactor model often gives an improper
solution, an error variance estimated below zero, and the estimate is
one of two kinds. It can be a finite maximum of the likelihood: on
`holzinger_swineford`, nine tests in three groups of three (visual
perception, cubes, and lozenges; paragraph comprehension, sentence
completion, and word meaning; addition, counting groups of dots, and
straight and curved capitals) give visual perception an error variance
of -41, and the likelihood falls on either side of that value. Such a
solution is reported. The maximal reliability is not defined for it, and
neither is *H* for a factor on which an item's standardized loading
reaches one, nor an omega whose items' error variances sum below zero;
those rows are `NA`, with a warning that names the items (beside the
warning of class `dmar_heywood_warning` that
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
gives), and the `improper` attribute holds the negative error variances.
The other rows are reported as estimated, and an improper solution's
indices are better read as a sign that the model asks more of the data
than they hold. Or the likelihood can have no maximum at all. With a
fourth group factor on the Woody-McCall mixed fundamentals, revised
paper form board, and flags tests beside those three groups, the
residual covariances of the fourth group's tests, beyond the general
factor, are all negative, which no group factor can reproduce (a single
factor implies three covariances whose product is not negative), and the
likelihood keeps rising, toward a limit, as the error variance of the
paper form board falls without bound. lavaan's optimizer stops at -1163
under FIML, where its steps stop gaining, a point set by its starting
values and tolerance rather than by the data, and every index read there
would describe the optimizer. So when an error variance is negative, the
function fits the model again with it held further out (at ten times its
estimate, or one model implied variance below it, whichever is further),
and when the log likelihood there is higher, or lower by no more than
\\10^{-6}\\ per case (a likelihood that flat leaves the estimate as much
to where the optimizer stopped as a ridge does), it stops, for every
`output` (`"fit"` included), with an error of class `dmar_no_maximum`
(see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md))
that names the test and prints the lavaan code of the fit; there,
`bounds = "pos.var"` (which holds every error variance at zero or above,
and which `...` passes on from here as well) or a different model, such
as one with the covariances of the group's items freed in place of the
group factor, are the choices.

**Ordered items.** `fa_bifactor()` fits continuous items, and `ordered`
is refused. Fit the model to ordered items with
[`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html)
(`ordered = TRUE`, `orthogonal = TRUE`, `std.lv = TRUE`) and pass the
fit to
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
which reports the indices on the metric of the latent responses behind
the categories; the bifactor coefficients of a sum of the observed
categories are not implemented.

**The other bifactor functions.** `output = "fit"` returns the lavaan
fit, so the same model feeds
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
which reports the indices on the standardized metric, without intervals,
and can take the observed composite variance as the denominator, and
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
which gives them bootstrap intervals, an alternative to the delta method
intervals here that does not rest on the normal approximation and that
sets improper replications aside.

Row names under `term`: `lambda_<factor>_<j>` (loadings, the general
factor's numbered over every item, the items of the group factors in
their order and then those of `general_only`; `lambda_<factor>` when
equated), `phi_<factor>` (the factor variances, fixed at 1),
`psi_<factor>_<j>` and `nu_<factor>_<j>` (the error variances and
intercepts of a group factor's items; `psi_<general>_<j>` and
`nu_<general>_<j>` for an item on the general factor alone, numbered as
its general loading), the indices above, and
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)'s
fit rows `chi_square` through `H1`, whose versions and the lines under
the table are as there.

## References

Aguirre-Urreta, M. I., Rönkkö, M., & McIntosh, C. N. (2019). A
cautionary note on the finite sample behavior of maximal reliability.
*Psychological Methods, 24*(2), 236–252.
[doi:10.1037/met0000176](https://doi.org/10.1037/met0000176)

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
*Topics in applied multivariate analysis* (pp. 72–141). Cambridge, UK:
Cambridge University Press.

Green, S., & Yang, Y. (2018). Empirical underidentification with the
bifactor model: A case study. *Educational and Psychological
Measurement, 78*(5), 717–736.
[doi:10.1177/0013164417719947](https://doi.org/10.1177/0013164417719947)

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Holzinger, K. J., & Swineford, F. (1937). The bi-factor method.
*Psychometrika, 2*(1), 41–54.
[doi:10.1007/BF02287965](https://doi.org/10.1007/BF02287965)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

McDonald, R. P. (1999). *Test theory: A unified treatment*. Erlbaum.

Reise, S. P. (2012). The rediscovery of bifactor measurement models.
*Multivariate Behavioral Research, 47*(5), 667–696.
[doi:10.1080/00273171.2012.715555](https://doi.org/10.1080/00273171.2012.715555)

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
bifactor models: Calculating and interpreting statistical indices.
*Psychological Methods, 21*(2), 137–150.
[doi:10.1037/met0000045](https://doi.org/10.1037/met0000045)

## See also

[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
for the correlated factors model of the same items;
[`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md)
for the second-order model whose first-order factors are the group
factors, which is nested in the bifactor model;
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md)
and
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
for the indices of a fitted bifactor model on the standardized metric
and their bootstrap intervals;
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md)
for the factor scores and their determinacy;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
and
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md).

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
data(holzinger_swineford)
hs_groups <- list(
  spatial = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
              "t4_lozenges"),
  verbal  = c("t5_general_information", "t6_paragraph_comprehension",
              "t7_sentence", "t8_word_classification", "t9_word_meaning"),
  memory  = c("t14_word_recognition", "t15_number_recognition",
              "t16_figure_recognition", "t17_object_number",
              "t18_number_figure", "t19_figure_word"))

# A general factor on the 15 tests and three group factors: the
# indices, each with a delta method interval.
fa_bifactor(holzinger_swineford, hs_groups, output = "measurement")
#> Bifactor model: the general factor (general) loads on all 15 items; group factors spatial (4 items), verbal (5 items), memory (6 items)
#> Every factor's variance is fixed at 1 and every pair of factors is uncorrelated
#> 
#>  syntax term                        estimate se    test_statistic p_value 
#>         omega_general               0.876    0.011 79.393         < 0.0001
#>         omega_H_general             0.552    0.044 12.568         < 0.0001
#>         ECV_general                 0.407    0.041 9.844          < 0.0001
#>         H_general                   0.789    0.024 32.200         < 0.0001
#>         maximal_reliability_general 0.681    0.035 19.486         < 0.0001
#>         PUC_general                 0.705    <NA>  <NA>           <NA>    
#>         omega_spatial               0.695    0.037 18.566         < 0.0001
#>         omega_HS_spatial            0.266    0.088 3.022          0.0025  
#>         ECV_spatial                 0.626    0.116 5.411          < 0.0001
#>         H_spatial                   0.393    0.092 4.280          < 0.0001
#>         maximal_reliability_spatial 0.362    0.097 3.712          0.0002  
#>         omega_verbal                0.902    0.010 90.375         < 0.0001
#>         omega_HS_verbal             0.612    0.069 8.877          < 0.0001
#>         ECV_verbal                  0.330    0.075 4.418          < 0.0001
#>         H_verbal                    0.823    0.031 26.450         < 0.0001
#>         maximal_reliability_verbal  0.780    0.046 17.077         < 0.0001
#>         omega_memory                0.728    0.027 27.006         < 0.0001
#>         omega_HS_memory             0.482    0.067 7.204          < 0.0001
#>         ECV_memory                  0.382    0.080 4.776          < 0.0001
#>         H_memory                    0.624    0.043 14.405         < 0.0001
#>         maximal_reliability_memory  0.595    0.050 11.922         < 0.0001
#>  ci_lower ci_upper
#>  0.852    0.896   
#>  0.465    0.636   
#>  0.329    0.490   
#>  0.737    0.833   
#>  0.609    0.745   
#>  <NA>     <NA>    
#>  0.617    0.763   
#>  0.130    0.467   
#>  0.389    0.815   
#>  0.233    0.579   
#>  0.199    0.565   
#>  0.881    0.920   
#>  0.472    0.736   
#>  0.203    0.489   
#>  0.754    0.876   
#>  0.678    0.857   
#>  0.672    0.778   
#>  0.355    0.611   
#>  0.241    0.546   
#>  0.536    0.704   
#>  0.495    0.688   
#> 
#> Confidence level: 95%
#> 
#> omega, omega_H, and omega_HS describe unit-weighted sums of the items as scored (bifactor_indices() gives them for sums of the standardized items); ECV is formed from the standardized loadings, and H and maximal_reliability do not depend on the metric of the items; PUC depends on the structure alone, so it has no interval
#> 
#> For a group factor, omega is the reliability of the total score of its items, omega_HS the share of that score's variance due to the group factor alone, and ECV the general factor's share of the items' common variance; H is the index of Rodriguez, Reise, and Haviland (2016), which counts the other factors' common variance as uncorrelated error, and maximal_reliability is the reliability the optimally weighted items reach as a measure of the factor

# The fit statistics of the same model, one row.
generics::glance(fa_bifactor(holzinger_swineford, hs_groups))
#>   chi_square df     p_value      cfi       tli      rmsea  rmsea_low rmsea_high
#> 1   122.7287 75 0.000423185 0.969631 0.9574834 0.04598074 0.03067577 0.06031208
#>         srmr      AIC      BIC    logLik estimator missing     test fit_indices
#> 1 0.04394018 27584.12 27806.55 -13732.06        ML    fiml standard      robust
#>   rmsea_ci_level h1_diagnostic wald_calibration_note
#> 1            0.9          <NA>                  <NA>

# The lavaan fit feeds bifactor_indices(), whose omega, omega_H, and
# omega_HS describe sums of the standardized tests (0.577 for omega_H,
# against 0.552 for the sum of the scores above).
fit <- fa_bifactor(holzinger_swineford, hs_groups, output = "fit")
bifactor_indices(fit)
#>  factor  ECV   omega omega_H omega_HS PUC   H    
#>  general 0.407 0.878 0.577   <NA>     0.705 0.789
#>  spatial 0.626 0.678 <NA>    0.256    <NA>  0.393
#>  verbal  0.330 0.915 <NA>    0.612    <NA>  0.823
#>  memory  0.382 0.728 <NA>    0.455    <NA>  0.624
#> 
#> Denominator: composite variance implied by the fitted model
#> 
#> For a group factor, ECV is the general factor's share of the common variance of that factor's items (one minus it is the group factor's share), and omega is the reliability of their total score
#> H is the index of Rodriguez, Reise, and Haviland (2016), which counts the other factors' common variance as uncorrelated error; the maximal reliability the fitted model implies is general 0.681, spatial 0.362, verbal 0.780, memory 0.595 (attribute "maximal_reliability"; ci_bifactor_indices() gives each an interval)
```
