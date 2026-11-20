# Sample Size Planning for Accuracy in Parameter Estimation for Reliability Coefficients

Computes the necessary sample size for the confidence interval on a
population reliability coefficient (coefficient alpha or coefficient
omega, depending on the assumed measurement model) to have expected
width no larger than `width`, or (when `assurance` is supplied) to be no
wider than `width` with the specified probability. This is the accuracy
in parameter estimation (AIPE) counterpart to power based planning for
reliability and is the companion of
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md).
Without an assurance the sample size comes from a closed form for the
expected width of the interval `type` names: the normal theory variance
of coefficient alpha of van Zyl, Neudecker, and Nel (2000) for the
parallel model and for the normal theory interval of the tau equivalent
model (Terry & Kelley, 2012, Eq. 19), and the delta method variance of
the model implied coefficient for the factor analytic interval of the
tau equivalent model (equal loadings) and for the congeneric model
(coefficient omega). When `assurance` is supplied, a Monte Carlo search
starting from the closed form finds the smallest *N* delivering the
requested assurance (Terry & Kelley, 2012).

Coefficient alpha (Guttman, 1945; subsequently popularized by Cronbach,
1951) is planned for under the parallel and tau equivalent (*i.e.*, the
so called True Score) measurement models; coefficient omega (McDonald,
1999) is planned for under the congeneric model. The coefficient is the
one
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
and
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
estimate from the raw item scores: the population, however it is given,
is used as a covariance matrix, so items of unequal variance are planned
for as they are, and a correlation matrix plans for standardized items.

When the coefficient will be estimated from incomplete data by full
information maximum likelihood, the interval at a given *N* is wider
than it would be with complete data. The `fmi` argument plans for that:
it is the fraction of missing information anticipated for the
coefficient, and the plan grows to carry the information that complete
data would have carried (see the section on planning under missing
data).

## Usage

``` r
ss_aipe_reliability(
  model = NULL,
  type = NULL,
  width = NULL,
  true_cov_matrix = NULL,
  conf_level = 0.95,
  assurance = NULL,
  data = NULL,
  n_items = NULL,
  true_item_covariance = NULL,
  true_factor_loading = NULL,
  true_residual_var = NULL,
  G_initial = 500,
  G = 5000,
  start_N = NULL,
  verbose = FALSE,
  fmi = 0,
  seed = NULL
)
```

## Arguments

- model:

  The measurement model assumed for the population. Accepts
  (case-sensitive aliases shown in parentheses): `"Parallel"`
  (`"parallel"`, `"SB"`, `"Spearman Brown"`, `"Spearman-Brown"`, `"sb"`)
  for the strictly parallel items model; `"True Score"`
  (`"True Score Equivalent"`, `"True-Score Equivalent"`, `"Equivalent"`,
  `"Tau Equivalent"`, `"tau-equivalent"`, `"Tau-Equivalent"`,
  `"True-Score"`, `"true-score"`, `"true score"`, `"Cronbach"`,
  `"cronbach"`, `"Chronbach"`, `"alpha"`) for the tau equivalent model
  (in which case the function plans for coefficient alpha); or
  `"Congeneric"` (`"congeneric"`, `"omega"`, `"Omega"`) for the
  congeneric model (in which case the function plans for coefficient
  omega).

- type:

  The method used to construct the confidence interval on the
  reliability coefficient: either `"Factor Analytic"` (McDonald, 1999),
  available for all three measurement models, or `"Normal Theory"` (van
  Zyl, Neudecker, & Nel, 2000), available for the parallel and tau
  equivalent models only (with the congeneric model it is an error).

- width:

  The desired full width of the two-sided confidence interval.

- true_cov_matrix:

  A symmetric population covariance matrix among the items (a planning
  value, not a sample covariance matrix), used as it is to imply the
  population reliability and its sampling distribution; a matrix with a
  negative eigenvalue is refused. A correlation matrix is a covariance
  matrix of standardized items and plans for them.

- conf_level:

  Confidence level (i.e., \\1 - \alpha\\, where \\\alpha\\ is the Type I
  error rate). Default `0.95`.

- assurance:

  Optional probability with which the realized interval is to be no
  wider than `width` (a sample that returns no interval counts as one
  that is wider). When `NULL` (the default), the planner targets the
  *expected* width; when supplied (e.g., 0.80, 0.85, 0.95), the function
  follows the closed form with a Monte Carlo search.

- data:

  A data set from which the population covariance matrix should be
  inferred, from its complete cases; the number of items is the number
  of its columns.

- n_items:

  Number of items (*J* in the formulas of the Details). Needed when the
  population is given by `true_item_covariance` or `true_factor_loading`
  with `true_residual_var`; with `true_cov_matrix` or `data` it is the
  number of columns and need not be given.

- true_item_covariance:

  The presumed common covariance *c* of every pair of items, which with
  `true_residual_var` gives the population: the common loading is the
  square root of *c*, so two items covary *c* and an item's variance is
  *c* plus its `true_residual_var`. It is the inter-item correlation of
  Terry and Kelley (2012) when the items are standardized, that is, when
  `true_residual_var` is \\1 - c\\; otherwise items *j* and *k*
  correlate \\c / \sqrt{(c + \psi^2_j)(c + \psi^2_k)}\\, with
  \\\psi^2_j\\ the residual variance of item *j*. With
  `true_item_covariance = .3` and `true_residual_var = .2`, for example,
  the items correlate .6, not .3. One value.

- true_factor_loading:

  Population factor loading or loadings: one number for the Parallel and
  True Score models, whose items share one loading, and a vector of
  `n_items` loadings for the Congeneric model. A vector for the Parallel
  or True Score model is an error that names the Congeneric model, since
  unequal loadings are what distinguish it.

- true_residual_var:

  Population residual (unique, or error) variance or variances: one
  number for the Parallel model, whose items share one, and a vector of
  `n_items` variances for the True Score and Congeneric models, each
  greater than zero; with the normal theory interval of the True Score
  model, any values that leave the covariance matrix they imply a
  covariance matrix (see Details).

- G_initial:

  Number of Monte Carlo replications in each pass of the initial
  assurance search (default 500): a single whole number of at least 1.

- G:

  Number of Monte Carlo replications in each pass of the final assurance
  verification (default 5000), the count the reported
  `assurance_attained` rests on: a single whole number of at least 1.
  Anything else given for either count is refused with an error of class
  `dmar_count_invalid`. Neither is used without `assurance`.

- start_N:

  Optional starting total sample size for the iterative assurance
  search.

- verbose:

  If `TRUE`, prints the current sample size and empirical assurance at
  each step of the Monte Carlo search.

- fmi:

  The fraction of missing information anticipated for the reliability
  coefficient when it is estimated by full information maximum
  likelihood from incomplete data: a single number in \\\[0, 1)\\, with
  the default `0` for complete data, which returns exactly the complete
  data plan. It is not the proportion of missing values. Writing
  \\\gamma\\ for `fmi`, the plan divides the complete data \\N - 1\\ by
  \\1 - \gamma\\; see the section on planning under missing data for the
  derivation and for where a planning value comes from.

- seed:

  Optional seed for the Monte Carlo assurance search. The default `NULL`
  draws from the current state of the random number generator; a number
  makes the search reproducible, and the caller's generator state is
  restored when the function returns. The closed form (no `assurance`)
  draws nothing.

## Value

A `data.frame` with columns `term` and `value`. Without `assurance` the
data frame has a single row, `"necessary_N"`, giving the necessary *N*.
With `assurance` supplied, the data frame has five rows: `"necessary_N"`
(necessary *N*), `"width"` (echo of the target width), `"assurance"`
(echo of the requested probability), `"assurance_attained"` (the
assurance achieved at the returned *N* in the Monte Carlo verification),
and `"G"` (the number of Monte Carlo replications of the final pass).
When `fmi` is above 0, `necessary_N` is the size under missing data, the
rows are otherwise unchanged, and two attributes record the plan: `fmi`,
the value planned for, and `necessary_N_complete`, the size complete
data would need (for the closed form, its solution rounded up, while the
plan maps the unrounded solution; for an assurance plan, the size at
which the assurance was simulated, so `assurance_attained` refers to
it). With `fmi = 0` neither attribute is added. An assurance plan
carries the number of samples of its final Monte Carlo pass that
returned no interval as the `no_interval` attribute (zero when every
sample returned one), printed under the table when above zero; they are
counted as wide in `assurance_attained`. Both forms carry the class
`dmar_ss_aipe` of the AIPE planners:
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
gives the planned *N*, beside the target width for an assurance plan
(the closed form's table holds only `necessary_N`, so its width is
`NA`), and
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
adds the remaining rows as columns, and the `fmi`,
`necessary_N_complete`, and `no_interval` attributes when the plan has
them (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

**The closed forms.** Each closed form takes the sampling variance of
the coefficient to be \\V / (N - 1)\\ and solves \\2 z \sqrt{V / (N -
1)} = w\\ for *N*, so that \\N = 4 z^2 V / w^2 + 1\\, where *w* is
`width` and *z* the standard normal quantile for `conf_level`. For the
parallel model *V* is the compound symmetry form \\2 (1 - \alpha)^2 J /
(J - 1)\\ of van Zyl et al. (2000), with *J* items, which at a
population of parallel items is also the variance of the factor analytic
interval, so the parallel model has one closed form for both types. For
the normal theory interval of the tau equivalent model *V* is the
general variance of alpha of van Zyl et al. (Terry & Kelley, 2012, Eqs.
7, 8, and 19), evaluated on the population covariance matrix as the
formula defines it. For the congeneric model *V* is the delta method
variance of the model implied omega, \\\omega = (\sum \lambda_j)^2 /
((\sum \lambda_j)^2 + \sum \psi_j^2)\\, with \\\lambda_j\\ the loadings
(`true_factor_loading`) and \\\psi_j^2\\ the error variances
(`true_residual_var`): the gradient of omega taken through the inverse
of the expected information of the single-factor model at the population
(Raykov, 2002; Terry & Kelley, 2012, Eqs. 13 to 16), the variance behind
`reliability_omega(denominator = "model_implied", ci_method = "ml")`,
whose squared standard error is \\V / N\\; taking it as \\V / (N - 1)\\,
as for alpha, adds one case and keeps the planning under missing data
below the same for every model. When the population is given by
`true_cov_matrix` or `data`, the single-factor model is first fit to it
by maximum likelihood. Terry and Kelley (2012, p. 390) had no closed
form for omega: they used the alpha formula as if true score equivalence
held only as a starting value, and found the size whose mean interval
width met the target by an a priori Monte Carlo simulation. That
starting value, taken as the plan for omega, is too large when the
loadings differ among standardized items and can be far too small when
the item variances differ. The delta method form gives the expected
width to first order. The mean width of the factor analytic interval at
the planned *N* meets the target where the samples at that *N* stay away
from an improper solution and from empirical underidentification, and
exceeds it where some do not, whose intervals are far wider, while the
median width stays near the target (within 4 percent of it in the
simulations below, above it in some). Which samples come near depends on
the loadings, the number of items, and *N* together, so no bound on the
loadings alone marks the populations the closed form serves. With three
items the single-factor model is just identified, each squared loading
is a ratio of the covariances (\\\lambda_1^2 = \sigma\_{12} \sigma\_{13}
/ \sigma\_{23}\\), and a sample covariance near zero makes the solution
run away: three loadings of .3 planned for a width of .20 take 687
cases, at which the median width is .202 and the mean .264, and 3
percent of samples give no interval (1,000 samples), while at a width of
.10 (2743 cases) the mean is on target; loadings of .5, .4, and .3 at
.20 (462 cases) average .242, with a median of .206, and three loadings
of .4 (468 cases) .211, with a median of .199 (500 samples each). With
four loadings of .3 the mean width was on target at .20 (527 cases) but
averaged .323 at .30 (235 cases), where a tenth of the samples give no
interval. Where the factor rests on weak loadings the same happens with
more items: loadings of .7, .7, .1, and .1 (omega .46) planned for a
width of .10 take 974 cases, at which the median width is .104 and the
mean .181, and four loadings of .15 take 3437 cases, with a median of
.098 and a mean of .148. In the other populations simulated for this
planner, with four or more items whose loadings share a sign, the mean
was within a few percent of the target (for loadings of opposite sign,
see below). For a population like those above, plan with an assurance
above .5, whose search counts the widths themselves.

For the factor analytic interval of the tau equivalent model *V* is the
same delta method variance with the loadings held equal, the variance
behind
`reliability_alpha(estimator = "model_implied", ci_method = "ml")`: the
single loading's information is the sum of the loading block of the
information above, and its gradient the sum of the loadings' gradients.
That coefficient is the maximum likelihood estimator under the model, so
its variance is no larger than the normal theory variance of alpha, and
it is smaller whenever the error variances differ: with loadings .7 and
error variances .1, .2, 1, and 2, a width of .10 takes 242 cases with
the factor analytic interval and 367 with the normal theory one. With
equal error variances the two coincide.

**Populations that are refused.** The planning population has to be one
the model describes. A population given as `true_cov_matrix` must be a
covariance matrix (no negative eigenvalue), and when the single-factor
model is fit to it, the fit must be proper. Where that model is fit, for
the congeneric model and for the tau equivalent model with the factor
analytic interval, every error variance in `true_residual_var` must be
above zero, as the parallel model requires: a negative one is not a
variance, and a zero one makes an item free of error, the boundary at
which the single-factor solution is improper. The same population given
as `true_cov_matrix` is refused through its fit. With the normal theory
interval of the tau equivalent model nothing is fit: the loading and the
error variances only build the covariance matrix the plan is made from,
which must be a covariance matrix, as `true_cov_matrix` must, and the
plan is the one the same matrix given as `true_cov_matrix` gives. The
congeneric model needs at least three items with a nonzero loading,
without which it is not identified at the population and omega has no
sampling variance. With the factor analytic interval of the parallel and
tau equivalent models, items to which the single-factor model with equal
loadings gives a loading of zero share no common factor in that model,
and the coefficient that interval estimates is zero, on the boundary of
its range, where the closed form does not describe the interval. For a
population given by its loading (`true_item_covariance` or
`true_factor_loading`) that is a loading of zero, and for one given as
`true_cov_matrix` or `data` the loading of the maximum likelihood fit.
For items of equal variance it happens exactly when their covariances
average zero or below (alpha of the population matrix at or below zero);
with unequal variances the fit weighs each covariance by the inverse of
its two items' variances, and the two can disagree either way (three
items of variance 1 that correlate .3 and a fourth of variance 4 that
correlates -.25 with them have alpha -0.276 and a fitted coefficient of
.315, which the tau equivalent model plans for). The parallel model's
closed form is a function of alpha of the population matrix, and it is
refused when alpha is at or below zero; a population whose fitted
coefficient is nonetheless positive is far from parallel, and the
message names the tau equivalent model. A supplied population under
`model = "Parallel"` must have equal diagonal entries and equal
off-diagonal entries, each within \\\sqrt{\epsilon}\max \|S\|\\, where
\\\epsilon\\ is machine precision. This check uses the original item
units. A sample covariance will generally depart from an exact parallel
population; specify the intended population with `true_factor_loading`
and `true_residual_var` when parallelism is a planning assumption rather
than a property of `true_cov_matrix`.

**Where the congeneric closed form is refused.** The delta method takes
omega, a square of the sum of the loadings \\L\\ over the variance of
the composite, to be normal with its first-order variance, which needs
\\L\\ to be large relative to its standard error. When it is not, as
when loadings of opposite sign nearly cancel in the unweighted sum,
omega behaves near zero as the square of a normal variable, and its
first-order variance understates the width of its interval; where the
loadings cancel exactly, the gradient of omega vanishes and the closed
form gives one case. The closed form is therefore checked at the *N* it
gives, or at `n_items + 1` cases when it gives fewer, the fewest at
which the factor analytic interval can be computed (with fewer, the
sample covariance matrix is singular). Write \\\sigma_L\\ for the
standard error of \\\hat L\\ there, \\r = \sigma_L / \|L\|\\, and \\P\\
for the sum of the error variances. Near zero \\\hat\omega\\ is about
\\\hat L^2 / P\\ and its delta method standard error about \\2 \|\hat
L\| \sigma_L / P\\, so with \\\hat L\\ normal, the interval the
estimation functions report, clamped at zero, averages
\$\$\frac{E\[\min(4 z \|X\|, X^2 + 2 z \|X\|)\]}{4 z / r}, \qquad X \sim
N(1 / r, 1),\$\$ times the width planned for. That ratio stays below
one, the clamped interval making the plan conservative, until \\r\\ is
about 1.5 at 95 percent confidence (1.36 at 90 percent and 1.66 at 99),
and it exceeds one beyond; the closed form is refused once it does. In
simulations of the factor analytic interval at the planned *N* for
loadings .7, -.7, .5, and a fourth loading of -.2, -.32, -.35, or -.37
(\\r\\ of .41, 1.07, 1.52, and 2.00 there, the last two refused), the
mean width was .082, .082, .089, and .121 for a target of .10. With a
fourth loading of -.4 (omega .0038, \\r = 3.35\\) the closed form's 49
cases give a mean width of .22, where the interval needs about 95, and
with -.5 the loadings cancel and the delta method gives one case. An
assurance plan is not refused: its search counts the widths of the
intervals themselves, starting at no fewer than `n_items + 1` cases. A
congeneric closed form below `n_items + 1` is floored at that size
before any missing-information inflation, with a warning. This is only a
feasibility minimum. It does not establish that the expected width
target is met at such a small sample size; an assurance plan is
recommended.

Whenever the loadings of a congeneric population differ in sign, the
planner warns. An item that loads opposite to the others lowers omega of
the unweighted sum, as an item keyed in the opposite direction does when
it is not reverse scored. The plan is for omega of the sum of the items
as they are given; reverse score such items before planning, as before
the coefficient is estimated.

The parallel and tau equivalent models need no such check. Their one
loading enters the covariance matrix only through its square, the common
covariance, of which their coefficient is a smooth function, so its
variance does not vanish near zero, and there the clamped interval only
makes the plan conservative: four parallel items that covary .03 (alpha
.11) planned for a width of .30 take 362 cases, at which the factor
analytic interval averages .236 wide.

**The assurance search.** The Monte Carlo search draws samples from the
population covariance matrix and, at each candidate sample size,
computes the realized confidence interval width with the same machinery
the estimation side of the package uses. For `type = "Factor Analytic"`
the single-factor model is fit by maximum likelihood (equal loadings for
the parallel and tau equivalent models, free loadings for the congeneric
model) and the interval is the delta method Wald interval on the model
implied reliability, the interval of
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)`(denominator = "model_implied", ci_method = "ml")`
and of
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)`(estimator = "model_implied", ci_method = "ml")`.
For `type = "Normal Theory"` the interval uses the van Zyl, Neudecker,
and Nel (2000) closed form standard error for the tau equivalent model
and its compound symmetry simplification for the parallel model, the
same closed form behind `reliability_alpha(ci_method = "ml")`. The
congeneric model has no normal theory form, so `type = "Normal Theory"`
is an error there.

The search starts from the closed form size, and with the factor
analytic interval at no fewer than `n_items + 1` cases, the fewest whose
sample covariance matrix the single-factor model can be fit to (a highly
reliable population planned for a wide interval can have a closed form
below that). When that size falls short of the assurance, the size is
raised one case at a time until a pass of `G_initial` samples meets it,
and a final pass of `G` samples confirms the size, again raising it if
needed (Terry & Kelley, 2012). When the closed form size already meets
the assurance, which happens when the assurance is near one half or
below, or when the interval of `type` is narrower than the one the
closed form describes (as the factor analytic interval of the parallel
model can be when the population given as `true_cov_matrix` or `data` is
not parallel), the search moves down: to the size at which the widths it
has drawn would just meet the target (the width of a Wald interval
shrinks as \\1 / \sqrt{N - 1}\\), and from there up one case at a time
to the first size that meets it.

Each sample is treated as the estimation functions treat it, so the
search counts the intervals a user would be given. A sample for which
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
or the model implied
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
would return no interval (a single-factor fit that does not converge,
three items with no single-factor solution, or a solution that has run
away, whose error variances sum below zero so that the model implied
coefficient is above one; McDonald, 2004) gave no interval, and so no
narrow one: it is counted as wider than `width`. So is a sample whose
likelihood has no maximum, whatever the sum, one error variance running
off below zero while the likelihood is flat or still rising, which the
estimation functions refuse with an error of class `dmar_no_maximum`
(see
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)).
The assurance is the probability that the study yields an interval no
wider than `width` (Terry & Kelley, 2012, p. 373), not that probability
among the studies that yield an interval. When some samples of the final
pass returned no interval, the result says how many (the `no_interval`
attribute, printed under the table), since at the planned *N* that share
of studies would end without a coefficient. Each confidence limit is
clamped into \[0, 1\], as the estimation functions clamp it, so a
negative alpha whose interval lies wholly below zero has the interval
\[0, 0\]. A size at which no sample returns an interval stops the search
with the first error.

The cost of the assurance search depends on the interval. The normal
theory intervals are closed forms, so a search at the default
`G_initial = 500` and `G = 5000` finishes in seconds. With the factor
analytic interval a one factor model is fit at every Monte Carlo
iteration and at every candidate sample size the search visits, so a
congeneric plan with an assurance runs for tens of seconds even at
`G_initial = 50` and `G = 200` and for many minutes at the defaults. The
examples on this page therefore stop at the closed form for the
congeneric plan; the assurance version is the same call with `assurance`
supplied (for example, `assurance = .80`). A `seed` makes the search
reproducible and leaves the caller's random number stream as it was;
without one the search draws from, and advances, the current stream.

## Note

Not all of the items can be entered into the function to represent the
population values. For example, either `data` can be used, or
`true_cov_matrix`, or `n_items`, `true_item_covariance`, and
`true_residual_var`, or `n_items`, `true_factor_loading`, and
`true_residual_var`. With a large number of replications (`G`) this
function may take considerable time.

## Planning Under Missing Data

The fraction of missing information (FMI) of a parameter is \$\$\gamma =
1 - \frac{SE_C^2}{SE_O^2},\$\$ where \\SE_O\\ is the standard error of
the estimate from the data as observed and \\SE_C\\ the standard error
the same sample would have given had it been complete (Savalei &
Rhemtulla, 2012, Eq. 4). It is the proportion of the information about
the parameter that the missing values take away, one minus the
asymptotic efficiency of the incomplete data estimator relative to the
complete data one. At a fixed *N*, missing values therefore multiply the
sampling variance of the estimate by \\1 / (1 - \gamma)\\ and the width
of a Wald interval by \\1 / \sqrt{1 - \gamma}\\, the width inflation
factor of Savalei and Rhemtulla (2012); an FMI of .56 makes the interval
1.5 times as wide.

Every closed form on this page has a sampling variance of the form \\V /
(N - 1)\\, so its information about the coefficient is \\(N - 1) / V\\.
With missing data that information becomes \\(N - 1)(1 - \gamma) / V\\,
and the width target is met when \\(N - 1)(1 - \gamma)\\ equals \\N_c -
1\\, where \\N_c\\ is the complete data solution of the closed form
before it is rounded up to a whole number: \$\$N = \frac{N_c - 1}{1 -
\gamma} + 1.\$\$ The factor scales \\N - 1\\ rather than *N* because
\\N - 1\\ is what the planner's information is proportional to; the
target sample size of Savalei and Rhemtulla (2012), \\N / (1 -
\gamma)\\, is the same idea for an information proportional to *N*, and
the two differ by \\\gamma / (1 - \gamma)\\ cases (one case at an FMI of
.5). An FMI of .5 doubles \\N - 1\\, not *N*. The plan is the smallest
whole *N* at or above that value. The complete data sample size, the
whole number \\N_c\\ rounds up to, travels with the result as the
`necessary_N_complete` attribute and is printed under the table. Because
the plan maps the unrounded \\N_c\\, the formula applied to the printed
whole number can give a larger *N* than the plan, by up to about \\1 /
(1 - \gamma)\\ cases (the population of the first example below needs 38
cases with complete data and 42 at an FMI of .1, while the formula
applied to 38 gives 43).

With an `assurance`, the Monte Carlo search runs exactly as it does for
complete data, starting from the complete data closed form, and the size
\\N_e\\ it finds is mapped back with the same formula. The search
therefore takes the distribution of the interval width at *N* under
missing data to be the complete data width distribution at the effective
size \\(N - 1)(1 - \gamma) + 1\\, which keeps the larger relative
variability of the widths of a smaller sample (scaling the complete data
widths at *N* by the width inflation factor would not). This is a
first-order approximation, and
[`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md)
checks it by simulating incomplete data and analyzing them by full
information maximum likelihood (its `fmi` and `missing_prop` arguments).
The `assurance_attained` row is the assurance found at \\N_e\\, the
`necessary_N_complete` attribute is \\N_e\\ itself (a whole number, so
the formula applied to it reproduces the plan), and a `start_N` is read
as a size under missing data and converted to its effective size before
the search.

**Where the FMI comes from.** The FMI belongs to the coefficient, not to
the items, and it is not the proportion of missing values: it depends on
the missing data mechanism, the model, and how strongly the items are
related. Information lost on one item is partly recovered from the items
correlated with it, so the FMI can be lower than the rate of missing
data; under data missing completely at random it is bounded by the
largest rates of univariate and pairwise missingness, and under data
missing at random it can be well above the rate of missing data (Savalei
& Rhemtulla, 2012). At the design stage it comes from a pilot study or a
comparable data set analyzed the way the planned study will be: fit the
single-factor model by full information maximum likelihood (for example
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)`(..., missing = "fiml", output = "fit")`,
with `equal_loading = TRUE` for the tau equivalent coefficient that the
Parallel and True Score plans target) and take the FMI of its `omega_`
defined parameter from
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md).
That FMI refers to the cases the fit used: a case with no observed item
carries no information, and the fit sets it aside (the `N_excluded`
attribute counts such cases). When the planned *N* counts every case
recruited, those who will answer no item included, fold them in as \\1 -
(N\_{used} / N)(1 - \gamma)\\, with \\N\_{used}\\ the cases the fit used
and *N* those recruited;
[`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md)
reports the FMI its simulated data realized in the same way. A sample
FMI is imprecise. Chen and Savalei (2021) found it close to unbiased
even at \\N = 50\\ but variable from sample to sample: precise estimates
took several hundred cases at overall missing rates of 10 to 20 percent
and more than 1,000 at 30 percent. Planning with the upper limit of an
interval for the pilot's FMI, rather than with its point estimate, is
the conservative choice. Von Hippel (2020) recommends the same use of a
pilot's upper limit when deciding how many imputations to run; carrying
it over to sample size planning is our extension.

Two things `fmi` leaves alone. The `data` argument still estimates the
planning covariance matrix from the complete cases of the pilot data
([`na.omit()`](https://rdrr.io/r/stats/na.fail.html)); `fmi` is what
carries the missing data of the planned study. And the planned interval
is the same Wald interval with a normal critical value; the FMI changes
how much information the sample carries, not the reference distribution.

## Warning

In some conditions the factor analytic model fit by
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) (via
lavaan) may fail to converge, and you may see a non-convergence message
from lavaan. The Monte Carlo assurance search counts a non-converged
iteration as an interval wider than `width` and continues, and the
number of such iterations in the final pass is printed under the result;
a sample size at which no iteration converges stops the search with the
first error. Frequent non-convergence usually means the model is poorly
determined by the data, for example because of a small sample size or a
poorly behaved covariance matrix, and it raises the planned *N*, since a
study that yields no interval has not yielded a narrow one.

## References

Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
missing information from full information maximum likelihood. *Frontiers
in Psychology, 12*, 667802.
[doi:10.3389/fpsyg.2021.667802](https://doi.org/10.3389/fpsyg.2021.667802)

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika, 16*(3), 297–334.

Guttman, L. (1945). A basis for analyzing test-retest reliability.
*Psychometrika, 10*(4), 255–282.
[doi:10.1007/BF02288892](https://doi.org/10.1007/BF02288892)

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

McDonald, R. P. (1999). *Test theory: A unified treatment*. Mahwah, NJ:
Lawrence Erlbaum Associates.

McDonald, R. P. (2004). Respecifying improper structures. *Structural
Equation Modeling, 11*(2), 194–209.

Raykov, T. (2002). Analytic estimation of standard error and confidence
interval for scale reliability. *Multivariate Behavioral Research, 37*,
89–103.
[doi:10.1207/S15327906MBR3701_04](https://doi.org/10.1207/S15327906MBR3701_04)

Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
fraction of missing information from full information maximum
likelihood. *Structural Equation Modeling, 19*(3), 477–494.
[doi:10.1080/10705511.2012.687669](https://doi.org/10.1080/10705511.2012.687669)

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution
of the maximum likelihood estimator of Cronbach's alpha. *Psychometrika,
65*(3), 271–280.
[doi:10.1007/BF02296146](https://doi.org/10.1007/BF02296146)

von Hippel, P. T. (2020). How many imputations do you need? A two-stage
calculation using a quadratic rule. *Sociological Methods & Research,
49*(3), 699–718.
[doi:10.1177/0049124117747303](https://doi.org/10.1177/0049124117747303)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md)
to check a plan, under complete or missing data, by simulation.

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Expected confidence interval width (closed form, no Monte Carlo search)
# for six parallel items that correlate .6: standardized items, so
# true_residual_var is 1 - true_item_covariance.
ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
  n_items = 6, true_item_covariance = .6, true_residual_var = .4,
  conf_level = .95, assurance = NULL)
#>  term        value
#>  necessary_N 38   
#> 
#> Confidence level: 95%

# The same plan when the coefficient will be estimated by full information
# maximum likelihood and is expected to lose 30 percent of its information
# to missing values (an FMI of .30 from a pilot study). The complete data
# N minus 1 is divided by 1 minus .30, and the complete data answer prints
# under the table.
ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
  n_items = 6, true_item_covariance = .6, true_residual_var = .4,
  conf_level = .95, fmi = .30)
#>  term        value
#>  necessary_N 54   
#> 
#> Confidence level: 95%
#> 
#> Planned for a fraction of missing information of 0.3 for the coefficient; with complete data, N = 38.

# The assurance cases run a Monte Carlo search. G_initial = 20 and
# G = 100 keep the examples quick; a reported plan deserves the
# defaults of 500 and 5000. The seed makes each search reproducible and
# leaves the random number stream of the session where it was.

# Same population, now targeting an assurance.
ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
  n_items = 6, true_item_covariance = .6, true_residual_var = .4,
  conf_level = .95, assurance = .85, G_initial = 20, G = 100, seed = 113)
#>  term               value
#>  necessary_N        54   
#>  width              0.100
#>  assurance          0.850
#>  assurance_attained 0.860
#>  G                  100  
#> 
#> Confidence level: 95%

# The assurance plan under the same FMI. The search runs at the complete
# data size, the necessary N of the call above (kept as the attribute
# necessary_N_complete), and maps it back.
ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
  n_items = 6, true_item_covariance = .6, true_residual_var = .4,
  conf_level = .95, assurance = .85, G_initial = 20, G = 100, seed = 113,
  fmi = .30)
#>  term               value
#>  necessary_N        77   
#>  width              0.100
#>  assurance          0.850
#>  assurance_attained 0.860
#>  G                  100  
#> 
#> Confidence level: 95%
#> 
#> Planned for a fraction of missing information of 0.3 for the coefficient; with complete data, N = 54, the size at which the assurance was simulated.

# The true score (tau equivalent) model takes true_residual_var as a
# vector of length n_items (number of items) while true_item_covariance
# stays a single value. Error
# variances around 1 - .5 give item variances from .9 to 1.1 and
# correlations near .5, and the plan is for alpha of these items as they
# are, unequal variances and all.
ss_aipe_reliability(model = "True Score", type = "Normal Theory",
  width = .1, n_items = 5, true_item_covariance = .5,
  true_residual_var = c(.4, .5, .6, .4, .6),
  conf_level = .95, assurance = .85, G_initial = 20, G = 100,
  seed = 113)
#>  term               value
#>  necessary_N        137  
#>  width              0.100
#>  assurance          0.850
#>  assurance_attained 0.860
#>  G                  100  
#> 
#> Confidence level: 95%

# Congeneric model, planned from the item loadings and error variances rather
# than from a single covariance. With assurance = NULL the necessary N comes
# from the delta method variance of the model implied omega at the
# population, the variance behind reliability_omega(denominator =
# "model_implied", ci_method = "ml"). Adding an assurance to this plan fits
# a one factor model at every Monte Carlo iteration and is the slow case,
# so the page stops at the closed form (see Details).
ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic", width = .15,
  n_items = 4, true_factor_loading = c(.8, .7, .7, .8),
  true_residual_var = c(.4, .5, .5, .4),
  conf_level = .95, assurance = NULL)
#>  term        value
#>  necessary_N 52   
#> 
#> Confidence level: 95%

# Planning from a presumed population correlation matrix among the items,
# which plans for standardized items.
pop_mat <- rbind(
  c(1.0000000, 0.3813850, 0.4216370, 0.3651484, 0.4472136),
  c(0.3813850, 1.0000000, 0.4020151, 0.3481553, 0.4264014),
  c(0.4216370, 0.4020151, 1.0000000, 0.3849002, 0.4714045),
  c(0.3651484, 0.3481553, 0.3849002, 1.0000000, 0.4082483),
  c(0.4472136, 0.4264014, 0.4714045, 0.4082483, 1.0000000))
ss_aipe_reliability(model = "True Score", type = "Normal Theory", width = .15,
  true_cov_matrix = pop_mat, conf_level = .95, assurance = .85,
  G_initial = 20, G = 100, seed = 113)
#>  term               value
#>  necessary_N        113  
#>  width              0.150
#>  assurance          0.850
#>  assurance_attained 0.890
#>  G                  100  
#> 
#> Confidence level: 95%
```
