# Second-Order Confirmatory Factor Analysis Model

Fits the second-order factor model, in which first-order factors account
for what the items of each subscale share and one or more second-order
factors account for what the first-order factors share. The model is
described as in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
by naming the items of each first-order factor, together with the
first-order factors of each second-order factor, and is fit by maximum
likelihood to raw data (by full information maximum likelihood, the
default, or after listwise deletion) or to a covariance matrix. The
table reports the first-order and second-order loadings, the error
variances, and the intercepts with their confidence intervals; the
indices that describe how much of the total score and of each subscale
score the second-order factor accounts for: omega and omega hierarchical
(omega_H) for the total score, omega and omega hierarchical subscale
(omega_HS) for each subscale score, the explained common variance (ECV),
coefficient *H*, and the maximal reliability, each with a delta method
standard error and a confidence interval, and the percentage of
uncontaminated correlations (PUC), which depends on the structure alone
and so has neither; and the fit of the model. The Schmid-Leiman solution
the estimates imply (Schmid & Leiman, 1957), the loadings of every item
on the second-order factor and on what its first-order factor measures
beyond it, is returned as an attribute.

## Usage

``` r
fa_second_order(
  data = NULL,
  factors,
  second_order,
  S = NULL,
  N = NULL,
  M = NULL,
  equal_loading = FALSE,
  equal_intercept = FALSE,
  equal_error = FALSE,
  equal_second_order_loading = FALSE,
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

# S3 method for class 'dmar_fa_second_order'
tidy(x, ...)

# S3 method for class 'dmar_fa_second_order'
glance(x, ...)
```

## Arguments

- data:

  A raw data matrix or data frame, rows are respondents and columns
  include the items named in `factors`. Supply exactly one of `data` or
  `S`.

- factors:

  Named list of the first-order factors. Each element names a
  first-order factor and gives the character vector of its items, two or
  more; each item is listed under one first-order factor. An item that
  is not a column of `data` or an item of `S` is refused with an error
  of class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

- second_order:

  Named list of the second-order factors. Each element names a
  second-order factor and gives the character vector of its first-order
  factors, named as `factors` names them; each first-order factor is
  listed under one second-order factor. One second-order factor needs
  three or more first-order factors, or two with
  `equal_second_order_loading`; with two or more second-order factors,
  which are always correlated (uncorrelated second-order factors are not
  offered), two first-order factors each are enough. See *Details* for
  the structures that are refused.

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

  Logical, or a named logical vector with one element per first-order
  factor. `TRUE` constrains the loadings of a first-order factor's items
  to a single value. Defaults to `FALSE`.

- equal_intercept:

  Logical, or a named logical vector with one element per first-order
  factor. `TRUE` constrains the intercepts of a first-order factor's
  items to a single value; this requires the mean structure. Defaults to
  `FALSE`.

- equal_error:

  Logical, or a named logical vector with one element per first-order
  factor. `TRUE` constrains the error variances of a first-order
  factor's items to a single value. Defaults to `FALSE`.

- equal_second_order_loading:

  Logical, or a named logical vector with one element per second-order
  factor. `TRUE` constrains the loadings of a second-order factor's
  first-order factors to a single value, so that they correlate equally
  with it. Defaults to `FALSE`.

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
  data required). Ordered items are not available.

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

  :   (default) The first-order loadings, disturbance variances (fixed
      at 1), error variances, and intercepts, the second-order loadings,
      variances (fixed at 1), and correlations, then the indices, then
      the fit rows.

  `"measurement"`

  :   The indices alone.

  `"summary"`

  :   The [`summary()`](https://rdrr.io/r/base/summary.html) of the
      lavaan fit (not a data frame).

  `"standardized"`

  :   [`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
      of the fit, whose rows for the second-order loadings are the
      correlations of the first-order factors with their second-order
      factor; the indices are not among its rows (see *Details*). Under
      `se = "none"` its `se`, `z`, `pvalue`, `ci.lower`, and `ci.upper`
      columns are `NA`.

  `"fit"`

  :   The lavaan fit, which
      [`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html)
      takes to test the model against the bifactor model (see
      *Details*).

- ...:

  For `fa_second_order()`, additional arguments forwarded to
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
  the loadings of a second-order model or for its indices. It follows
  `...`, so it must be named in full.

- x:

  A `dmar_fa_second_order` table returned by `fa_second_order()`.

## Value

For `output = "verbose"` (default) and `"measurement"`, a `data.frame`
(classes `dmar_fa_second_order`, `dmar_tbl`) with columns `syntax`,
`term`, `estimate`, `se`, `test_statistic`, `p_value`, `ci_lower`, and
`ci_upper`. Its attributes record the model: `factors` and
`second_order`; `equal`, a list naming the factors whose loadings, error
variances, intercepts, and second-order loadings are equated;
`has_means`; `structural_df`, the degrees of freedom of the second-order
structure (0 when it is just identified); `omega_metric`
(`"raw_item_sum"`); `improper`, the negative error variances named by
item (empty for a proper solution); `improper_second_order`, the
standardized second-order loadings at one and the second-order
correlations beyond one, named by their terms (empty for a proper
solution); `schmid_leiman`, the Schmid-Leiman solution described above;
and `wald_df`. The verbose table also carries
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

The model is also called the higher-order factor model, and sometimes
the hierarchical factor model. The function is named for the structure
it fits because Schmid and Leiman (1957) and Yung, Thissen, and McLeod
(1999) use "hierarchical factor model" for the other model, the one with
uncorrelated factors in layers, of which the bifactor model fit by
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
is the two-layer case.

**The model.** An item \\i\\ of first-order factor \\f\\ is \$\$x_i =
\nu_i + \lambda_i f + e_i,\$\$ and the first-order factor \\f\\ of
second-order factor \\g\\ is \$\$f = \gamma_f g + u_f,\$\$ where the
disturbance \\u_f\\ is the part of \\f\\ that \\g\\ does not account
for. The model is identified as lavaan's `std.lv = TRUE` identifies it:
the variance of every second-order factor is fixed at 1 (and its mean at
0 when the mean structure is modeled), and so is the variance of every
disturbance, so every first-order loading \\\lambda_i\\ and every
second-order loading \\\gamma_f\\ is free. Both constraints are linear.
The variance of \\f\\ is then \\1 + \gamma_f^2\\, and its correlation
with \\g\\, the standardized second-order loading, is \\\gamma_f /
\sqrt{1 + \gamma_f^2}\\. Two or more second-order factors are correlated
with each other (the `phi` rows), and the function offers correlated
second-order factors only: uncorrelated ones would hold every
first-order factor under one second-order factor uncorrelated with every
first-order factor under another, and would leave a second-order factor
of two first-order factors unidentified
([`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html) with
`std.lv = TRUE` and `orthogonal = TRUE` fits that model). A factor fits
as well reflected: a first-order factor with every one of its loadings
and its second-order loading reversed in sign, a second-order factor
with its second-order loadings and its correlations reversed. Each
first-order factor is reported in the direction in which its loadings
sum to a positive value, and then each second-order factor in the
direction in which its second-order loadings do, so a subscale whose
items are scored in reverse shows positive loadings and a negative
second-order loading. A first-order factor under a second-order factor
with equal loadings is not reflected alone, which would break the
equality, and keeps the direction the common loading gives it. Should a
reflection fail, a warning says so and the factor is reported as the
optimizer left it, which changes the signs of its loadings and nothing
else. Every loading, error variance, and intercept carries a label, as
in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
and the indices are lavaan defined parameters of those labels, so the
fit itself gives each one a delta method standard error. The
standardized solution is derived from this fit of the covariances (or of
the raw data), not from a fit to correlations.

**Scope and identification.** A structure outside the model, or one the
data cannot identify, is refused before anything is fit, with an error
of class `dmar_second_order_unidentified` (which inherits
`dmar_refusal`; see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md))
that names the factors or items at issue and gives the kind of refusal
as its `reason` element, as in
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md).
Not identified (`"unidentified"`): a first-order factor of one item,
whose loading cannot be told apart from its error variance; a
second-order factor of one first-order factor, which cannot be told
apart from it; and one second-order factor over two first-order factors
with free loadings, whose one correlation is the product of the two
standardized loadings (equal loadings identify them). Outside the model
(`"scope"`): an item under two first-order factors, or a first-order
factor under two second-order factors. Names (`"names"`): an item or a
first-order factor listed twice under one factor, a first-order factor
under no second-order factor, and names that coincide. The
identification checks run once the other arguments and the data have
passed theirs. With two or more second-order factors, a second-order
factor of two first-order factors is identified through the correlations
of its first-order factors with those of the others, and is empirically
underidentified when the second-order factors barely correlate. With
incomplete raw data under FIML, first-order factors that no case
observes together (no case has an item of one and an item of the other)
are refused too, as `"unidentified"`, when the pairs of first-order
factors the data do observe together do not identify the second-order
loadings: three first-order factors under one second-order factor with
one pair never observed together leave two correlations for three
loadings, while four with one such pair leave five for four and are fit.
When lavaan's information matrix cannot be inverted, or is singular,
whether it is the observed or the expected information, every standard
error, test, and interval is `NA`, the maximal reliability's included,
with a warning, as in
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md).

**A just identified structural part.** The second-order loadings and
correlations account for the correlations of the first-order factors,
\\k(k - 1)/2\\ of them for \\k\\ first-order factors. Exactly three
first-order factors under one second-order factor have three
correlations and three second-order loadings, and two with equal
loadings have one of each, so the second-order structure is just
identified: it is another parameterization of the correlated factors
model of the same first-order factors, which
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
fits. Whenever that model's three correlations have a positive product
and imply standardized second-order loadings below one in absolute value
(its one correlation is positive, for two factors), the two models fit
exactly alike, with the same chi square on the same degrees of freedom,
and the test of fit says nothing about the second-order structure; a
line under the table says so, and the `structural_df` attribute, the
correlations less the second-order parameters, is 0. On
`holzinger_swineford`, nine tests in three groups of three (visual
perception, cubes, and lozenges; paragraph comprehension, sentence
completion, and word meaning; addition, counting groups of dots, and
straight and curved capitals) give a chi square of 85.31 on 24 degrees
of freedom under both models. Four or more first-order factors leave the
second-order structure something to test.

**The Schmid-Leiman solution.** Substituting the second line of the
model into the first gives \$\$x_i = \nu_i + \lambda_i \gamma_f g +
\lambda_i u_f + e_i,\$\$ a model in which every item loads on the
second-order factor and on its first-order factor's disturbance, all of
them uncorrelated with variance 1: the transformation of Schmid and
Leiman (1957), which turns a second-order solution into one with
uncorrelated factors in layers. In the standardized metric an item's
loading on \\g\\ is \\\lambda^\*\_i \gamma^\*\_f\\ and its loading on
\\u_f\\ is \\\lambda^\*\_i \sqrt{1 - \gamma^{\*2}\_f}\\, with
\\\lambda^\*\_i\\ its standardized first-order loading and
\\\gamma^\*\_f\\ the standardized second-order loading of its factor
(Yung, Thissen, & McLeod, 1999, Eqs. 4 and 10). These are the
`general_loading` and `group_loading` columns of the `schmid_leiman`
attribute, beside `item`, `factor` (the first-order factor),
`second_order`, and `uniqueness`, one minus the sum of their squares.
Within a first-order factor the ratio of the two loadings is
\\\gamma^\*\_f / \sqrt{1 - \gamma^{\*2}\_f}\\ for every item: the
proportionality constraint by which the second-order model differs from
a bifactor model with the same groups, whose general and group loadings
are free (Yung, Thissen, & McLeod, 1999; Chen, West, & Sousa, 2006).

**The indices.** The indices are those of
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
computed from the Schmid-Leiman loadings of the items of each
second-order factor: in the unstandardized metric, the general loading
\\\lambda_i \gamma_f\\ and the group loading \\\lambda_i\\, beside the
error variance \\\psi_i\\, with model implied item variance \\v_i =
\lambda_i^2 (1 + \gamma_f^2) + \psi_i\\. The variance of the total
score, the unweighted sum of the items of a second-order factor \\g\\,
is \$\$\sigma^2_X = (\textstyle\sum_i \lambda_i \gamma\_{f(i)})^2 +
\sum_f (\sum\_{i \in f} \lambda_i)^2 + \sum_i \psi_i,\$\$ and the rows
are

- `omega_<g>`:

  Omega for the total score, the share of \\\sigma^2_X\\ that is not
  error variance: its reliability with every factor counted as true
  score.

- `omega_H_<g>`:

  Omega hierarchical, \\(\sum_i \lambda_i \gamma\_{f(i)})^2 /
  \sigma^2_X\\: the share of the total score's variance that is due to
  the second-order factor alone.

- `omega_<f>`:

  Omega for the subscale score of first-order factor \\f\\, the sum of
  its items, \\(\sum\_{i \in f} \lambda_i)^2 (1 + \gamma_f^2)\\ over
  that score's variance, which adds \\\sum\_{i \in f} \psi_i\\: the
  omega that
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  reports for a factor, here at the estimates of the second-order model.

- `omega_HS_<f>`:

  Omega hierarchical subscale, \\(\sum\_{i \in f} \lambda_i)^2\\ over
  the same variance: the share of the subscale score's variance that is
  due to the disturbance of \\f\\, what the subscale score says beyond
  the second-order factor.

- `ECV_<g>`, `ECV_<f>`:

  The explained common variance, formed from the standardized
  Schmid-Leiman loadings as
  [`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
  forms it: for \\g\\, its share of the common variance of its items;
  for \\f\\, the share of the common variance of the items of \\f\\ that
  is due to \\g\\, which in this model is \\\gamma^{\*2}\_f\\, the
  squared standardized second-order loading, since every item of \\f\\
  divides its common variance in that proportion.

- `H_<g>`, `H_<f>`:

  Coefficient *H* in the form of Rodriguez, Reise, and Haviland (2016),
  the closed form of Hancock and Mueller (2001) over the factor's
  Schmid-Leiman loadings with every other source of an item's variance
  counted as its error.

- `maximal_reliability_<g>`, `maximal_reliability_<f>`:

  The maximal reliability under the fitted model, \\c' \Sigma^{-1} c\\,
  with \\c\\ the covariances of all of the items with the factor and
  \\\Sigma\\ their model implied covariance matrix: the largest
  reliability any weighted sum of the items reaches as a measure of
  \\g\\, or of the disturbance of \\f\\. With several second-order
  factors, the items of the others count too, through the correlations
  of the second-order factors.

- `PUC_<g>`:

  The percentage of uncontaminated correlations, as a proportion: the
  share of the pairs of the items of \\g\\ that do not share a
  first-order factor. It depends on the structure, not on the data, so
  it has no standard error or interval.

With several second-order factors each one's rows describe the sum of
its own items. Rodriguez, Reise, and Haviland (2016) describe how each
index is read; the table reports the numbers with their uncertainty, and
whether a total or subscale score is worth reporting is the researcher's
judgment, on which the package draws no line. Omega, omega_H, and
omega_HS describe unit-weighted sums of the items as they are scored, as
in
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md);
ECV, *H*, the maximal reliability, and PUC do not depend on the metric
of the items. With `output = "standardized"`, lavaan would evaluate the
defined indices at the standardized values of their labels, and a
standardized first-order loading is scaled by the whole variance of its
factor, not by that of the disturbance the indices need, so the index
rows are left out of that output, and lavaan's warning that under
`equal_loading` or `equal_error` their standardized values are ambiguous
is not passed on.

**Intervals.** As in
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md):
the standard errors of omega, omega_H, omega_HS, ECV, and *H* are
lavaan's delta method standard errors of the defined parameters, that of
the maximal reliability is the delta method's with the derivatives taken
numerically, and the limits of every index are formed on the logit scale
from that standard error and transformed back (Browne, 1982; the
confidence interval transformation principle of Kelley, 2007), so they
lie within (0, 1); the `se` column is the standard error of the index
itself.

**The test against the bifactor model.** The second-order model is
nested in the bifactor model whose group factors are its first-order
factors: the bifactor model frees the proportionality constraint above
(Yung, Thissen, & McLeod, 1999; Chen, West, & Sousa, 2006). With
`output = "fit"` from both functions,
`lavaan::lavTestLRT(fa_second_order(..., output = "fit"), fa_bifactor(..., output = "fit"))`
gives the likelihood ratio test of the constraint, on as many degrees of
freedom as there are items less first-order factors (Yung, Thissen, &
McLeod, 1999), as the example shows. When the second-order structure is
just identified the second-order fit is that of the correlated factors
model, and the test compares the bifactor model with it. The comparison
is for one second-order factor; the bifactor model of
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
has one general factor.

**Finding the maximum.** In the parameterization above a standardized
second-order loading reaches one only as \\\gamma_f\\ runs off to
infinity, the loadings \\\lambda_i\\ of its factor shrinking toward zero
beside it, and lavaan's optimizer can stop on the way at a point that is
not the maximum of the likelihood, with other second-order loadings
drifted far from their values there, or fail to converge. The function
therefore fits the model a second time, in a parameterization that
writes each first-order factor on the scale of its second-order factor,
\\f = g + u_f\\, with the variance of \\u_f\\, which is \\1 /
\gamma_f^2\\, free and bounded below at zero, so that a standardized
loading of one is a point the optimizer reaches; the equality
constraints of the call hold exactly in it, as in the first. That fit
starts from the first and, when the first did not converge, came near
the boundary, or fell short, also from lavaan's own starting values; the
largest log likelihood is carried back to the parameterization above, so
the estimates and the fit rows are those of the maximum, and lavaan
accepts each fit only where its gradient is near zero. Log likelihoods
are compared per case (two within 1e-6 per case are taken to be the same
maximum), since the log likelihood is a sum over the cases and the
comparison should not depend on their number. On `holzinger_swineford`,
one second-order factor over the spatial, verbal, and speed tests of the
example and the three tests with missing values (`t24_woody_mccall`,
`t25_paper_form_board_r`, and `t26_flags`) stopped at a chi square of
251.77 in the first parameterization, with three second-order loadings
between 43 and 233; the maximum is 169.77. When the second fit cannot be
made, a warning says that the maximum could not be confirmed and names
any standardized second-order loading of .995 or more as possibly a
point where the optimizer stopped; when neither fit converges, the error
says so and names the standardized second-order loadings the first fit
was taking toward one.

**Improper solutions.** As in
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
a negative error variance makes the maximal reliabilities undefined, and
with them the *H* of a factor on which an item's standardized
Schmid-Leiman loading reaches one and an omega whose items' error
variances sum below zero; those rows are `NA`, with a warning that names
the items (beside the warning of class `dmar_heywood_warning` that
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
gives), and the `improper` attribute holds the negative error variances.
When the likelihood has no maximum at all, flat or still rising as an
error variance falls, the function stops for every `output`, as
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
does, with an error of class `dmar_no_maximum` that prints the lavaan
code of the fit (with a standardized second-order loading at one, that
of the second fit described above, in which the check is made); there,
`bounds = "pos.var"` or a different model are the choices. The second
order has improper solutions of its own, handled in the same way. A
standardized second-order loading of one, a Heywood case at the second
order, leaves the first-order factor's disturbance without variance (a
parameterization that let that variance be negative would estimate it at
or below zero). The loading is at one when the second fit above puts the
variance of \\u_f\\ at its bound of zero (within 1e-6 of it, a
standardized loading of .9999995 or more). The parameterization above
has no point there, so the fit holds the loading within 5e-9 of one,
where the log likelihood is within 1e-6 per case of its maximum, and
that factor's second-order loading and its loadings on its disturbance
have no finite estimate and are `NA`, as are the rows of its
disturbance, omega_HS, ECV, *H*, and the maximal reliability of that
first-order factor. The other rows are those of the maximum, their
standard errors computed with the loading held at one (with
`se = "bootstrap"` they are `NA`, with a warning, since the fit at that
point is computed without the optimizer and a bootstrap of it would not
refit the samples); the fit rows are those of the model at its maximum,
on its own degrees of freedom; and the factor's items have group
loadings of zero in the Schmid-Leiman attribute. With `output = "fit"`,
`"summary"`, or `"standardized"` the object holds the loading at that
point, and a warning says so; with `se = "bootstrap"` that object
carries no standard errors, tests, or intervals (lavaan's
`se = "none"`). Under FIML, lavaan cannot form its robust CFI, TLI,
RMSEA, and GFI for a fit computed at that point, so the `"summary"`
output shows them as `NA`, with one warning in place of lavaan's. A
correlation of two second-order factors beyond one in absolute value,
which leaves the covariance matrix of the latent variables short of
positive definite, leaves every maximal reliability `NA`, since each
reads that correlation. A warning names the factors, and the
`improper_second_order` attribute holds the standardized loadings at one
and the correlations beyond one, named by their terms. Two second-order
factors over two or three first-order factors each are prone to both: of
the ten ways of placing five first-order factors of
`holzinger_swineford` (spatial, verbal, speed, memory, and reasoning
tests), two under one second-order factor and three under the other,
every one puts the standardized second-order loading of the reasoning
factor at one, and seven of them also give a correlation of the
second-order factors beyond one.

**Ordered items.** `fa_second_order()` fits continuous items, and
`ordered` is refused.

Row names under `term`: `lambda_<factor>_<j>` (first-order loadings,
numbered within the factor; `lambda_<factor>` when equated),
`zeta_<factor>` (the disturbance variances, fixed at 1),
`psi_<factor>_<j>` and `nu_<factor>_<j>` (the error variances and
intercepts), `gamma_<second>_<factor>` (second-order loadings;
`gamma_<second>` when equated), `phi_<second>` (the second-order
variances, fixed at 1), `phi_<second>_<second>` (the correlations of the
second-order factors), the indices above, and
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)'s
fit rows `chi_square` through `H1`, whose versions and the lines under
the table are as there.

## References

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
*Topics in applied multivariate analysis* (pp. 72–141). Cambridge, UK:
Cambridge University Press.

Chen, F. F., West, S. G., & Sousa, K. H. (2006). A comparison of
bifactor and second-order models of quality of life. *Multivariate
Behavioral Research, 41*(2), 189–225.
[doi:10.1207/s15327906mbr4102_5](https://doi.org/10.1207/s15327906mbr4102_5)

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
bifactor models: Calculating and interpreting statistical indices.
*Psychological Methods, 21*(2), 137–150.
[doi:10.1037/met0000045](https://doi.org/10.1037/met0000045)

Schmid, J., & Leiman, J. M. (1957). The development of hierarchical
factor solutions. *Psychometrika, 22*(1), 53–61.
[doi:10.1007/BF02289209](https://doi.org/10.1007/BF02289209)

Yung, Y.-F., Thissen, D., & McLeod, L. D. (1999). On the relationship
between the higher-order factor model and the hierarchical factor model.
*Psychometrika, 64*(2), 113–128.
[doi:10.1007/BF02294531](https://doi.org/10.1007/BF02294531)

## See also

[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
for the bifactor model of the same groups, in which the second-order
model is nested;
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
for the correlated factors model of the first-order factors;
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
[`fa_bifactor()`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
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
  spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
  verbal  = c("t6_paragraph_comprehension", "t7_sentence",
              "t9_word_meaning"),
  speed   = c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals"),
  memory  = c("t14_word_recognition", "t15_number_recognition",
              "t17_object_number"))
hs_second <- list(general = c("spatial", "verbal", "speed", "memory"))

# One second-order factor over four first-order factors: the indices,
# each with a delta method interval.
so <- fa_second_order(holzinger_swineford, hs_factors, hs_second,
                      output = "measurement")
so
#> Second-order model: general over spatial (3 items), verbal (3 items), speed (3 items), memory (3 items)
#> Every second-order factor's variance and every first-order factor's disturbance variance is fixed at 1
#> 
#>  syntax term                        estimate se    test_statistic p_value 
#>         omega_general               0.768    0.024 31.802         < 0.0001
#>         omega_H_general             0.445    0.066 6.789          < 0.0001
#>         ECV_general                 0.358    0.041 8.751          < 0.0001
#>         H_general                   0.722    0.033 22.004         < 0.0001
#>         maximal_reliability_general 0.618    0.048 13.006         < 0.0001
#>         PUC_general                 0.818    <NA>  <NA>           <NA>    
#>         omega_spatial               0.639    0.037 17.121         < 0.0001
#>         omega_HS_spatial            0.212    0.094 2.248          0.0246  
#>         ECV_spatial                 0.667    0.145 4.596          < 0.0001
#>         H_spatial                   0.302    0.112 2.694          0.0071  
#>         maximal_reliability_spatial 0.281    0.115 2.451          0.0142  
#>         omega_verbal                0.872    0.014 61.655         < 0.0001
#>         omega_HS_verbal             0.621    0.067 9.277          < 0.0001
#>         ECV_verbal                  0.287    0.077 3.711          0.0002  
#>         H_verbal                    0.760    0.041 18.567         < 0.0001
#>         maximal_reliability_verbal  0.736    0.052 14.252         < 0.0001
#>         omega_speed                 0.676    0.033 20.308         < 0.0001
#>         omega_HS_speed              0.446    0.071 6.258          < 0.0001
#>         ECV_speed                   0.341    0.106 3.227          0.0012  
#>         H_speed                     0.546    0.067 8.207          < 0.0001
#>         maximal_reliability_speed   0.534    0.072 7.424          < 0.0001
#>         omega_memory                0.623    0.044 14.238         < 0.0001
#>         omega_HS_memory             0.499    0.063 7.906          < 0.0001
#>         ECV_memory                  0.198    0.080 2.483          0.0130  
#>         H_memory                    0.549    0.052 10.623         < 0.0001
#>         maximal_reliability_memory  0.545    0.054 10.133         < 0.0001
#>  ci_lower ci_upper
#>  0.717    0.811   
#>  0.323    0.574   
#>  0.282    0.441   
#>  0.653    0.781   
#>  0.522    0.706   
#>  <NA>     <NA>    
#>  0.563    0.708   
#>  0.082    0.449   
#>  0.358    0.879   
#>  0.132    0.551   
#>  0.114    0.543   
#>  0.841    0.897   
#>  0.484    0.741   
#>  0.161    0.458   
#>  0.671    0.831   
#>  0.623    0.824   
#>  0.608    0.738   
#>  0.314    0.586   
#>  0.171    0.565   
#>  0.416    0.671   
#>  0.394    0.669   
#>  0.534    0.704   
#>  0.378    0.621   
#>  0.085    0.398   
#>  0.447    0.647   
#>  0.440    0.647   
#> 
#> Confidence level: 95%
#> 
#> omega, omega_H, and omega_HS describe unit-weighted sums of the items as scored; ECV, H, and maximal_reliability are formed from the Schmid-Leiman loadings (attribute schmid_leiman) and do not depend on the metric of the items; PUC depends on the structure alone, so it has no interval
#> 
#> For a first-order factor, omega is the reliability of the total score of its items, omega_HS the share of that score's variance due to the factor's disturbance (what the second-order factor does not account for), ECV the second-order factor's share of the items' common variance, and H and maximal_reliability describe the disturbance

# The Schmid-Leiman loadings the estimates imply.
attr(so, "schmid_leiman")
#>                                item  factor second_order general_loading
#> 1              t1_visual_perception spatial      general       0.6132207
#> 2                          t2_cubes spatial      general       0.3490749
#> 3                       t4_lozenges spatial      general       0.4928028
#> 4        t6_paragraph_comprehension  verbal      general       0.4570818
#> 5                       t7_sentence  verbal      general       0.4578219
#> 6                   t9_word_meaning  verbal      general       0.4490990
#> 7                      t10_addition   speed      general       0.3362116
#> 8       t12_counting_groups_of_dots   speed      general       0.4176852
#> 9  t13_straight_and_curved_capitals   speed      general       0.3895176
#> 10             t14_word_recognition  memory      general       0.2919927
#> 11           t15_number_recognition  memory      general       0.2602009
#> 12                t17_object_number  memory      general       0.2436182
#>    group_loading uniqueness
#> 1      0.4328373  0.4366123
#> 2      0.2463920  0.8174377
#> 3      0.3478412  0.6361518
#> 4      0.7198484  0.2728945
#> 5      0.7210141  0.2705377
#> 6      0.7072765  0.2980701
#> 7      0.4674189  0.6684813
#> 8      0.5806877  0.4883408
#> 9      0.5415277  0.5550238
#> 10     0.5872058  0.5699296
#> 11     0.5232716  0.6584823
#> 12     0.4899234  0.7006253

# The second-order model is the bifactor model of the same groups with
# the general loadings proportional to the group loadings; the
# likelihood ratio test of that constraint, on 12 - 4 = 8 degrees of
# freedom.
second_order_fit <- fa_second_order(holzinger_swineford, hs_factors,
                                    hs_second, output = "fit")
bifactor_fit <- fa_bifactor(holzinger_swineford, hs_factors,
                            output = "fit")
print_anova(lavaan::lavTestLRT(second_order_fit, bifactor_fit))
#> 
#> Chi-Squared Difference Test
#> 
#> 
#>                  Df      AIC      BIC     Chisq Chisq diff Df diff Pr(>Chisq)
#> bifactor_fit     42 25301.57 25479.52  95.64521         NA      NA       <NA>
#> second_order_fit 50 25330.44 25478.72 140.50574   44.86052       8   < 0.0001
```
