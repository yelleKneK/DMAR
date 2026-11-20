# Factor Scores With Standard Errors, Intervals, and Determinacy

Computes, for every row given to a fitted lavaan model and every latent
variable, both standard factor scores: the regression score (Thomson,
1934; Thurstone, 1935) and the Bartlett score (Bartlett, 1937, 1938).
Each comes with a per-case standard error and confidence interval,
beside the number of the factor's items the case answered, the case's
missing data pattern, and the determinacy of each score for that pattern
(the correlation of the score with the factor it stands for). For
ordered items the two scores are their categorical counterparts in
lavaan, the empirical Bayes modal (EBM) score and the maximum likelihood
(ML) score.

## Usage

``` r
factor_scores(fit, data = NULL, conf_level = 0.95)

# S3 method for class 'dmar_factor_scores'
print(x, n = 3L, digits = getOption("dmar.digits", 3L), ...)

# S3 method for class 'dmar_factor_scores'
tidy(x, ...)

# S3 method for class 'dmar_factor_scores'
glance(x, ...)
```

## Arguments

- fit:

  A fitted lavaan model (from
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html),
  [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html), or
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)`(output = "fit")`)
  estimated from raw data. A fit to summary statistics (it holds no
  rows) and a fit that did not converge are refused, as are a two-level
  fit, a fit of continuous items that keeps incomplete cases without
  FIML (two-stage maximum likelihood, pairwise deletion or the available
  cases) and is conditional on exogenous covariates, a fit with
  regressions among the variables lavaan holds as latent in which the
  disturbance of an observed variable regressed on others covaries with
  a factor (an observed mediator, for one), a fit conditional on
  exogenous covariates with regressions and an observed variable that
  lavaan holds among the latent variables but regresses on nothing, and,
  with ordered items, a fit with a direct effect of an exogenous
  covariate on an observed variable that also depends on a factor
  directly, one with a regression on an observed variable that is not an
  exogenous covariate when the variable predicted also depends on a
  factor directly, the predictor is ordered, or, in a row the fit keeps,
  the predictor is missing where the variable predicted is observed and
  depends on a factor through no variable the row has observed, one with
  an observed variable whose disturbance covaries with a factor and
  whose row lavaan keeps among the loadings, one with exogenous
  covariates that is not conditional on them, and any fit under a lavaan
  older than 0.7-3 (see *Details*). A fit with ordered items and
  residual covariances is scored with a warning. For continuous
  unconditional fits, a disturbance related only to an observed
  conditioning variable leaves the factor's regression conditional mean
  valid. Its scores and intervals are retained; Bartlett scores, SEs,
  limits, and determinacy are `NA` with the reason when their stronger
  conditional-unbiasedness requirement fails.

- data:

  Optional but recommended: the data frame given to the model, with the
  same rows in the same order (extra columns are ignored). It is checked
  against the data the model holds (every observed variable present, the
  same number of rows, and the same values in the rows the model used),
  and a frame that differs is refused. With it, the table has one row
  for every row of `data`, including rows the fit did not use (deleted
  listwise, or with a missing group), and their reason and items
  answered are filled in; for a multiple group fit it must hold the
  group variable. Without it, the table runs from row 1 to the last row
  the fit holds, so rows at the end of the data that the fit did not use
  are not represented, and the reason given a row the fit did not use
  cannot say which of its values was missing. The scores are always
  those of the rows the fit was estimated from; scoring new observations
  is not what this function does.

- conf_level:

  Confidence level for the per-case intervals; defaults to 0.95.

- x:

  A `dmar_factor_scores` table returned by `factor_scores()`.

- n:

  For [`print()`](https://rdrr.io/r/base/print.html), the number of rows
  shown for each factor (default 3); the pattern table shows up to
  `4 * n` rows.
  [`knit_print.dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  and
  [`as_kable`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  show the same summary with the default.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html), the decimal places
  of the scores, their standard errors and limits, and the determinacies
  (default `getOption("dmar.digits", 3)`; see
  [`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)).

- ...:

  Unused.

## Value

A `data.frame` (class `dmar_factor_scores`, `dmar_tbl`) with one row for
each row of the data and each factor, ordered by factor and then row, so
that `x$.derived_bartlett[x$factor == "f"]` lines up with the rows of
`data`. Its columns are always the same: `row` (the row of the data),
`group` (the group label; `NA` for a fit without groups), `factor`,
`.derived_regression`, `.derived_regression_se`,
`.derived_regression_ci_lower`, `.derived_regression_ci_upper`,
`.derived_bartlett`, `.derived_bartlett_se`,
`.derived_bartlett_ci_lower`, `.derived_bartlett_ci_upper`,
`items_observed` and `items` (the factor's items the case answered, and
its number of items; a higher order factor's items are those of the
factors it measures, and a loading the model holds at zero, written
`0*x` or `a*x` with `a == 0` directly or through a chain of equalities,
makes no item: the scores and reasons are those of the model without
it), `pattern` (the missing data pattern of the variables scored from,
within the group, numbered from the most complete),
`determinacy_regression`, `determinacy_bartlett` (for the case's
pattern), and `reason` (empty, or why a value is missing or borrowed;
reasons are separated by semicolons). The `.derived_regression` and
`.derived_bartlett` columns hold the EBM and ML scores for ordered
items.

The eight columns that hold the scores, their standard errors, and their
confidence limits carry the prefix `.derived_`, so that a score joined
to the observed data cannot be taken for an observed variable or collide
with one (a data set may well have a column named `regression`). The
`derived_columns` attribute names them by what each holds (`regression`,
`regression_se`, `regression_ci_lower`, `regression_ci_upper`, and the
same four for `bartlett`), so
`x[[attr(x, "derived_columns")[["bartlett"]]]]` is the column of
Bartlett scores. All are present under every supported model.

Attributes: `methods` (the lavaan methods behind the two columns),
`conf_level`, `estimator`, `missing`, `patterns` (one row per group,
pattern, and factor: the number of cases, the items observed, the
variables missing, both determinacies, `min_score_correlation`, and both
standard errors), `determinacy` (per group and factor, for complete
data: both determinacies and `min_score_correlation`), `scored_from`
(`"fit"` or `"indicators"`), `auxiliary` (the auxiliary variables set
aside), `ngroups`, and `lavaan_version`. `tidy()` returns the table as a
plain data frame; `glance()` returns one row with the numbers of rows,
rows scored, factors, groups, and patterns, the two methods, the
estimator, the missing data treatment, `scored_from`, and `conf_level`.

## Details

Both scores are returned because neither is better in general; which one
belongs in a later analysis depends on the role the factor plays there
(see *Which Score Belongs Where*). For relations among factors, the
scores are rarely the right tool at all: a latent variable model,
estimated in one step or by structural after measurement
([`lavaan::sam()`](https://rdrr.io/pkg/lavaan/man/sam.html); Rosseel &
Loh, 2024), estimates them without the bias a regression among scores
carries. The scores are for what only scores can do: describe and rank
individual cases, carry a factor into an analysis outside the model, or
check the model's assumptions case by case. Requires lavaan.

**The Two Scores.** For the observed items \\o\\ of a case, with
\\\Lambda_o\\ their rows of the loading matrix, \\\Theta\_{oo}\\ their
residual covariance matrix, \\\Phi\\ the model implied covariance matrix
of the factors, \\\Sigma\_{oo} = \Lambda_o \Phi \Lambda_o' +
\Theta\_{oo}\\, and \\\mu_o\\ the model implied item means,
\$\$\hat{\eta}\_R = \mathrm{E}(\eta) + \Phi \Lambda_o' \Sigma\_{oo}^{-1}
(y_o - \mu_o)\$\$ is the regression score, the mean of the factor's
distribution given the observed items under the normal factor model, and
\$\$\hat{\eta}\_B = \mathrm{E}(\eta) + (\Lambda_o' \Theta\_{oo}^{-1}
\Lambda_o)^{-1} \Lambda_o' \Theta\_{oo}^{-1} (y_o - \mu_o)\$\$ is the
Bartlett score, the generalized least squares (and, under normality,
maximum likelihood) estimate of the case's factor value treated as a
fixed quantity. The regression score is shrunk toward the factor mean,
by more when the items carry less information, and so has the smaller
mean squared error; the Bartlett score is not shrunk, so its expected
value given the case's factor value is that value (conditionally
unbiased; under missing data this needs a condition on why values are
missing, see *Missing Data*), at the price of a larger variance. Their
standard errors differ in kind: that of the regression score is the
posterior standard deviation, the square root of the diagonal of
\\\Phi - \Phi \Lambda_o' \Sigma\_{oo}^{-1} \Lambda_o \Phi\\, and that of
the Bartlett score is its sampling standard error given the case's
factor value, the square root of the diagonal of \\(\Lambda_o'
\Theta\_{oo}^{-1} \Lambda_o)^{-1}\\. The scores and their standard
errors are those of
[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)
with `se = "standard"`, aligned to the rows of the data (with auxiliary
variables, the same formulas applied to the indicators alone; see
below).

The Bartlett score of a factor exists only when the observed items
identify the factor's value, that is, when the factor's unit vector lies
in the row space of \\\Lambda_o\\, so that its value is a function of
the items' expected values. With simple structure that holds for every
factor with an observed item. With items shared among factors
(cross-loadings, a bifactor model) a pattern can fail it: a case that
answered one item loading on two factors, and nothing else of either,
has items that do not separate the two. The matrix \\\Lambda_o'
\Theta\_{oo}^{-1} \Lambda_o\\ is then singular, and lavaan falls back on
a generalized inverse, whose minimum norm solution is not an estimate of
either factor's value and whose standard error is not the error of
anything; such a factor's Bartlett score, standard error, and
determinacy are `NA` here with the reason. A factor the pattern does
identify keeps its score, which the generalized inverse gives uniquely.

**Which Score Belongs Where.** Using scores in place of factors in a
later regression (factor score regression) with the same type of score
on both sides gives inconsistent estimates of the structural slopes.
Skrondal and Laake (2001) showed that the revised scheme, the regression
score for an explanatory factor and the Bartlett score for a response
factor, gives consistent estimates of every parameter of a latent
regression model, provided the explanatory and the response factors are
fitted as separate measurement models; scores from one joint model of
both are not covered by that result and should not feed a regression
among the same factors, and fitting each factor alone (factorwise) is
consistent only in special cases such as the independent clusters model
with one explanatory factor. Even then the standard errors the second
step reports do not account for the scores having been estimated
(corrected ones may require resampling), no conventional test of fit is
available, and factor score regression works badly for nonlinear models
such as polynomial regression. Devlieger, Mayer, and Rosseel (2016)
compared regression and Bartlett factor score regression, the bias
avoiding method of Skrondal and Laake, and the bias correcting method of
Croon (2002). The bias avoiding method was unbiased for unstandardized
but not for standardized coefficients and had the largest standard error
bias of the four, and Croon's method, with the standard error they
developed for it, was the only suitable alternative to structural
equation modeling, comparable in bias, efficiency, mean squared error,
power, and Type I error rate, though with more standard error bias. The
differences among the methods disappeared when factor score
indeterminacy was low, and they advise assessing it first, which the
determinacy columns here allow. Croon's correction is local structural
after measurement with the maximum likelihood mapping, whose matrix is
the Bartlett score matrix (Rosseel & Loh, 2024), and Rosseel and Loh
recommend local structural after measurement over one-step estimation
when the sample is small to medium or the items have low to moderate
reliability, with as many separate measurement blocks as the model
permits. It is available as
[`lavaan::sam()`](https://rdrr.io/pkg/lavaan/man/sam.html). Rosseel and
Loh developed the framework for complete continuous data, so with
missing data or ordered items its two-step standard errors deserve a
check against a bootstrap.

**Missing Data.** Under full information maximum likelihood
(`missing = "ml"`), lavaan scores each case from its own observed items,
pattern by pattern, which is what \\o\\ denotes above. For continuous
items, the regression score so computed is the posterior weighted
maximum likelihood estimate that Estabrook and Neale (2013) call full
ML, and the Bartlett score so computed is their unweighted ML estimate.
In their simulation (one factor), the two were equally reliable with
complete data, and with missing values the posterior weighted score was
the more reliable (with five items and values missing completely at
random, a mean reliability of .780 against .767, and .752 for the scaled
person mean), which they attribute to its empirical Bayes component. A
case that answered none of a factor's items has no Bartlett score for
that factor; its regression score is kept, but it comes only from the
case's other observed variables, through the model (the items of
correlated factors, a covariate that predicts the factor, an outcome the
factor predicts), it is the factor's mean when none of them bears on the
factor, and its determinacy says how little that is. A case with no
observed variable has no score of either kind. Under listwise deletion
the rows with a missing value were not used by the fit and are not
scored; for continuous items, a fit with `missing = "ml"` scores them,
unless the missing value is an exogenous covariate's. With
`fixed.x = TRUE`, lavaan's default, every fit deletes a row whose
covariate is missing, FIML with `missing = "ml"` included, except one
with `missing = "ml.x"` (which a fit conditional on the covariates does
not allow), and the reason column names that remedy; with ordered items
no fit that lavaan can score keeps such a row. lavaan has no FIML for
ordered items: fitted by WLSMV, they take their parameters from pairwise
or listwise deletion, which assumes the values are missing completely at
random, and under pairwise deletion (`missing = "pairwise"`) the EBM and
ML scores use each case's observed items, so that fit scores the rows
listwise deletion leaves out.

A fit of continuous items that keeps the incomplete cases without FIML
is scored the same way, each case from its own observed items, with the
fit's estimates in the formulas above: a fit by two-stage maximum
likelihood (lavaan's `missing = "two.stage"` or `"robust.two.stage"`, as
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
fits with `missing = "two_stage"`) or by pairwise deletion
(`missing = "pairwise"`, or `"available.cases"`, which for continuous
items gives lavaan the same pairwise complete moments and estimates).
The two-stage estimates are consistent when the values are missing at
random, as the FIML estimates are (Savalei & Bentler, 2009); the
estimates from pairwise deletion require the values to be missing
completely at random. The scores treat the estimates as known, as they
treat the FIML estimates, so the scores, standard errors, and
determinacies of such a fit mean what they mean under FIML and differ
from FIML's only as the estimates differ. A fit without a mean structure
(the default for one group under pairwise deletion) centers its observed
variables at their sample means and its factors at the means lavaan
gives them (`lavInspect(fit, "mean.lv")`): zero, or for a factor a
covariate predicts, the prediction at the covariate's mean, as a FIML
fit has it.
[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)
scores such a fit as complete data, leaving every case with a missing
value without a score, so these scores are computed here, by the closed
forms `lavPredict()` applies pattern by pattern under FIML. Auxiliary
variables of the first stage of a two-stage fit (lavaan's `aux`) are not
in the model and do not enter the scores. Such a fit conditional on
exogenous covariates (`conditional.x = TRUE`) is refused.

Missing values affect the two scores differently. Conditional
unbiasedness of the Bartlett score and the coverage of its Wald interval
are properties of its sampling distribution given the case's factor
value, and with missing values they hold when, given the factors,
whether an item is observed does not depend on the residuals of the
items scored: values missing completely at random, or missing because of
the factors themselves or of variables related to the items only through
the factors. Rubin (1976) showed that sampling distribution inferences
may ignore the process that causes missing data when the missing values
are missing at random and the observed values are observed at random,
whereas likelihood and Bayesian inferences need only the first (with the
parameters of the process distinct from those of the data). Missingness
that depends on the observed values of the factor's own items, the kind
of missing at random that FIML is commonly used for, fails the second:
it selects each pattern on those items' residuals, so the Bartlett score
of a pattern is biased given the factor value, and its interval covers
the factor value less often than stated in that pattern. In a one-factor
model of four items with loadings of .7 and 5,000 cases, the last three
items missing whenever the first fell below -0.3, the Bartlett scores of
the cases left with one item averaged 0.82 below their factor values and
their 95% intervals covered 92.7% of them; with as many values missing
completely at random, -0.01 and 95.3%. The regression score is a
posterior quantity, so missingness that depends only on what the score
conditions on leaves its posterior interval with its average coverage
(94.6% in the first case). Missingness that runs through an auxiliary
variable the scores set aside (see *Auxiliary Variables and Covariates*)
is not of that kind for the regression score, since the pattern then
carries information about the factor beyond the observed items, and the
regression score of a pattern is biased. The Bartlett score is not, but
only when the auxiliary is related to the items through the factors
alone (independent of them given the factors); one that shares variance
with an item beyond the factor, as an earlier wave of the same item
does, selects the pattern on that item's residual, which fails the
condition above. In the same design with the last three items missing
whenever an auxiliary fell below -0.3, an auxiliary that correlated .4
with the factor and with nothing else left the cases with one item
Bartlett scores 0.05 below their factor values on average (within two
Monte Carlo standard errors of zero), with intervals covering 95.6%, and
regression scores 0.17 above them; one that correlated .4 with the
factor and .5 with the first item's residual left them Bartlett scores
0.56 below their factor values, with intervals covering 93.1%. A score
whose error has mean zero given the factor leaves the slopes of a
regression in which it is the response unbiased, and that is the
Bartlett score's place in the revised scheme of factor score regression
above, so missingness that depends on the factor's own observed items
weakens it there. The reliabilities of Estabrook and Neale (2013) quoted
above are for values missing completely at random.

**Determinacy.** The determinacy of a score is its correlation with the
factor, \\\[A \Lambda_o \Phi\]\_{kk} / \sqrt{\[A \Sigma\_{oo} A'\]\_{kk}
\Phi\_{kk}}\\ for the weight matrix \\A\\ of the score, the formula
lavaan uses for `fs.determinacy` and `fs.determinacy.bartlett`,
evaluated here for each missing data pattern. Its square is the
reliability of the score as a measure of the factor. The regression
score has the larger determinacy, since it is the linear function of the
observed items that correlates most with the factor; for a one-factor
model the two are equal, and with uncorrelated residuals their square is
coefficient *H* (Hancock & Mueller, 2001), the reliability of the
optimally weighted composite of the items. Per-pattern determinacy shows
the cost of each missing item. It is computed for continuous items; with
ordered items each case's standard error is given instead. Beside the
regression score's determinacy \\\rho\\, the display and the
`determinacy` and `patterns` attributes report \\2\rho^2 - 1\\
(`min_score_correlation`), the lowest correlation that two sets of
scores, each as consistent with the model as the other, can have
(Guttman, 1955, as Grice, 2001, p. 435, presents it). It is zero when
\\\rho^2 = .5\\, where two such sets can be uncorrelated, and negative
below that, so it says how far the scores of different cases can be
trusted to order the cases on the factor. The package reports the
quantity and draws no line.

**Intervals.** The interval around a regression score is a posterior
(empirical Bayes) interval: its coverage holds on average over cases,
not for each case's factor value. The interval around a Bartlett score
is a Wald interval for the case's factor value. Both treat the model's
estimated parameters as known, so they do not carry the parameters'
sampling variability. For ordered items the EBM interval approximates
the posterior around its mode, and the ML interval is a Wald interval.
The ML score of a factor exists only when the likelihood of the case's
responses has a finite maximum in it. A case whose observed items of a
factor all sit at the same end of their scales (all in the lowest
category, or all in the highest, for positive loadings) has none, since
the likelihood keeps rising as the factor's value moves out. With items
shared among factors (cross-loadings, a bifactor model) the likelihood
can also rise without bound along a direction that moves several factors
together, or stay flat along one that the observed items cannot
separate. The check is exact: a factor has a finite and unique ML score
when no direction along which no observed item's probability falls moves
it, which a nonnegative least squares problem decides (by Farkas'
lemma), and with simple structure it reduces to the rule for the items
at an end. It reads each item's loadings as lavaan's likelihood does,
including those of an observed variable regressed on the factors (an
outcome the factor predicts, or an item with a direct effect of another
factor), which lavaan holds among the regressions. A factor that fails
it has its ML score and interval `NA` with the reason. A factor that
passes it keeps its ML score, also when it shares observed items with a
factor whose likelihood rises without bound: as that factor moves out,
the shared items' probabilities tend to 1 whatever this factor's value,
so its maximum is that of the case's other items, which lavaan's
optimizer reaches within its tolerance (within 1.2e-4 of a direct
maximization of the likelihood for the factor of a cross-loading item in
the tests, and within 1.1% of the standard error in a bifactor model of
12 items). When lavaan's default optimizer stops on its iteration limit
and returns no score for the case, the score comes from lavaan's
`nlminb` optimizer. A higher order factor, which has no items of its
own, has no ML score at all: the likelihood of the responses depends on
the factors it measures, not on it, and lavaan's value is the center of
the small ridge its optimizer adds. Its EBM score, which the factor's
place in the model informs, is kept.

**Ordered Items Need lavaan 0.7-3.** The probability of an ordered
response compares each threshold with the mean of the item's latent
response, its intercept plus its loadings times the factor values. The
EBM and ML scores of lavaan 0.7-2 counted the intercept twice in that
comparison, so they were wrong in any model whose ordered indicators
have intercepts other than zero (a multiple group model with equal
thresholds, for one); lavaan 0.7-3 corrected it, and also changed how
the scores treat the scale factors of the delta parameterization. A fit
with ordered items is therefore refused under an older lavaan, with the
remedy, rather than scored with numbers that may be wrong.

**Multiple Groups.** The scores use each group's own parameters, and the
`group` column names the group. Scores from different groups are on a
common metric only to the extent the measurement model is invariant
across the groups (equal loadings for comparing relations that involve
the scores, and equal intercepts as well for comparing their levels);
see
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md).
Rows whose group is missing are not used by the fit; with `data` they
are kept, without a score.

**Auxiliary Variables and Covariates.** An observed variable that
indicates no factor and enters the model only through covariances with
the indicators is an auxiliary variable carried as a saturated correlate
(Graham, 2003), as the reliability functions of this package carry them.
It helps the estimates under missing data, but under that
parameterization it is uncorrelated with every factor and correlated
with the item residuals, so whatever it shares with the factor is
carried as residual covariance, and a score that conditions on it
subtracts true score variance. When such variables are present the
scores are computed from the indicators alone, with the fit's estimates,
the `scored_from` attribute is `"indicators"`, and a message says so. A
higher order factor, which has no items of its own, then has no Bartlett
score, nor, when the data have missing values, in the scores computed
here for a two-stage fit or one by pairwise deletion
(`missing = "pairwise"` or `"available.cases"`), nor in a fit in which
lavaan gives an observed variable a row of zeros among the loadings (see
*Observed Variables in Regressions*): lavaan forms that score through
the loadings of the factors it measures, a construction not repeated
here, and the reason column says so; its regression score is computed as
for any other factor. A covariate that predicts a factor (`f ~ x`) is
different: the scores condition on it as the model specifies (as on any
observed variable that enters a regression of the model), and Curran,
Cole, Bauer, Hussong, and Gottfredson (2016) found that including such
covariates in the scoring model improved score quality in nearly every
condition of their simulation.

**Observed Variables in Regressions.** An observed variable regressed on
others, such as an item with a direct effect of a covariate (`y1 ~ x`,
the multiple indicators multiple causes model of differential item
functioning) or an outcome the factor predicts (`z ~ f`), is held by
lavaan as a latent variable measured by itself without error, with its
loadings among the regressions and its residual variance among those of
the latent variables. The scores, their standard errors and
determinacies, the check of the ML score, and the check of the solution
take its loadings and residual variance from there, as
[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)
does, so an item with a direct effect enters its factor's Bartlett
score, which is then conditionally unbiased given the factor and the
covariate. An observed variable that lavaan holds the same way but
regresses on nothing, such as a criterion that covaries with a factor
(`f ~~ z`) or a covariate of a factor that also covaries with another
variable (so that lavaan does not count it among the exogenous
covariates), is conditioned on, as an exogenous covariate is: the
regression score is the factor's mean given it and the items, and the
Bartlett score is conditionally unbiased given the factors and it. In a
fit conditional on the covariates (`conditional.x = TRUE`) without
regressions among the variables lavaan holds as latent (no \\B\\
matrix), so is an observed variable regressed only on exogenous
covariates (`z ~ x`): lavaan then keeps the variable's own row among the
loadings, `lavPredict()`'s regression score is the factor's mean given
the items, the variable, and the covariates, and its Bartlett score is
conditionally unbiased given the factors, the variable, and the
covariates, also when the variable covaries with a factor (`f ~ x`,
`z ~ x`, `f ~~ z`;
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html) frees that
covariance by default). In a model with regressions among them,
`lavPredict()` gives a variable regressed on nothing a row of zeros
among the loadings, which leaves its covariances with the factors out of
the scores (with `f ~ x` and `f ~~ z` for three Holzinger and Swineford
tests, its regression score missed the conditional mean by up to 1.9),
so the scores of such a fit are computed here, with the formulas
`lavPredict()` uses; one conditional on the covariates
(`conditional.x = TRUE`) is refused, with the remedy. In a fit without a
mean structure and without exogenous covariates, `lavPredict()` counts
the mean of a variable held this way twice, as its intercept and as the
mean of the latent variable that holds it, which shifts by a constant
every score whose weights involve the variable: the regression score of
every factor related to it and the Bartlett score of a factor that
predicts it (for a factor of three Holzinger and Swineford spatial tests
that predicts a fourth test and correlates .50 with a factor of three
speed tests, by -0.39 and -0.61 for the first factor's regression and
Bartlett scores and by -0.07 for the second factor's regression score);
the scores of such a fit are computed here, with the formulas
`lavPredict()` uses and the means the fit implies. A higher order
factor's Bartlett score, which lavaan forms through the loadings of the
factors it measures, is lavaan's less the constant the double count adds
(-1.09 for a second order factor of three factors of three Holzinger and
Swineford tests each, when the spatial factor predicts the paper form
board test), so it equals that of the fit with a mean structure. In a
model with regressions among the variables lavaan holds as latent, a
variable regressed on others is represented by its measurement equation,
with its disturbance as the error, which is right only when that
disturbance is uncorrelated with the factors and with the variables the
scores condition on, so a fit in which such a variable predicts a
factor, directly or through others (an observed mediator), or has a
residual covariance with a factor is refused: `lavPredict()` then gives
a Bartlett score that is not conditionally unbiased and, when a factor
is involved, a regression score that is not the factor's conditional
mean (with `f ~ a` and `a ~ x` for three Holzinger and Swineford tests,
it missed the conditional mean by up to 1.08, with a correlation of
.92). Carried as a single indicator latent variable (`A =~ a` with
`a ~~ 0*a`), such a variable is scored correctly: written so, lavaan
fixes one scale constraint, the loading or, with `std.lv = TRUE`, the
latent variable's variance, so the refitted model is the same model,
with the refused fit's log likelihood and degrees of freedom (Bollen,
1989, Chapter 7). Written `A =~ 1*a` in a fit with `std.lv = TRUE`, it
fixes both, which is a different model whose scores are not those of the
refused one. In a multiple group fit, `group.equal` holds parameters
equal by kind, and the variable's intercept, residual variance, and
covariances become the latent variable's, so the refit is the same model
only when it carries the refused fit's constraints across the groups:
with equal intercepts, lavaan frees the latent variable's mean after the
first group unless it is fixed at zero in every group (`A ~ 0*1`; for
two schools with equal loadings and intercepts the refit without it had
a log likelihood of -4246.42 on 15 degrees of freedom, against -4246.51
on 16), and the equalities `group.equal` places on residual variances
and covariances, or on those of latent variables, must be given to the
variable's parameters as the refused fit had them and to no others. The
message says which apply.

With ordered items, lavaan's EBM and ML scores leave a direct effect of
an exogenous covariate on an observed variable out of each case's
likelihood, and they take an observed variable that is not an exogenous
covariate as zero where it predicts another. That leaves the scores
unchanged when the variable predicted depends on no factor directly: its
density given its predictors then does not involve the factors, so the
term left out changes each case's likelihood by a factor constant in
them. Such fits are scored: a covariate's effect on an observed outcome
unrelated to the factors (`w ~ x`), and a fully mediated chain of
continuous observed variables (`m ~ f` with `w ~ m`, whose EBM and ML
scores agree with a direct maximization of the model's likelihood for
four Holzinger and Swineford tests cut at their quartiles). A fit in
which the term left out involves the factors is refused: a direct effect
on a variable that also depends on a factor directly (`y1 ~ x` for an
item, whose thresholds are estimated given the covariate, so the scores
are those of a case with the covariate at zero; ML scores off by 0.11 to
0.44 for an item regressed on age), a regression of such a variable on
an observed one (`y4 ~ y3` for local dependence between two items, ML
scores off by up to 0.18 for the four tests; or `w ~ m + f`), a
regression on an ordered item, whose latent response is not observed
(`w ~ y3`), and a regression on an observed variable that a row the fit
keeps lacks while it has the one predicted, where the missing variable
depends on a factor through no variable the row has observed (pairwise
deletion, where the model's likelihood of the variable predicted then
runs through the missing one to the factors). A variable the row has
observed on the way cuts that dependence, and the row is scored: with
`m ~ f`, `w ~ m`, and `s ~ w`, a row lacking `w` but having `m` and `s`
is scored, while a row lacking `m` but having `w` refuses the fit. A fit
is refused, too, for an observed variable lavaan holds among the latent
variables whose disturbance covaries with a factor, regressed on nothing
or, in a fit without regressions among the latent variables, only on
exogenous covariates (`f ~ x` with `w ~ x`, where
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html) frees the
covariance of `f` and `w` by default): the EBM score takes the factor's
prior from its model implied mean given the covariates, which leaves the
variable's value out, so it is the posterior mode given the items and
covariates alone; written as a predictor of the factor (`f ~ x`), a
variable regressed on nothing is an exogenous covariate, which the EBM
score conditions on. A fit with ordered items and exogenous covariates
that is not conditional on them (`conditional.x = FALSE`), which
`lavPredict()` cannot score, is refused as well. The EBM and ML scores
of lavaan also take each case's likelihood as the product of its
observed variables' probabilities given the factors, so they leave
residual covariances among those variables out; a fit with one is scored
with a warning, since its scores are those of the model without it (with
a residual correlation of -0.12 between two of the four tests, the ML
scores of the 160 response patterns differed from the scores that
maximize the model's own likelihood by up to 0.11, 15% of their standard
error, and by 0.03 at the median).

An improper solution is scored with a warning: one with a negative
residual variance (a Heywood case, also of an observed variable
regressed on others), a negative factor variance, a covariance matrix of
the factors that is not positive definite (a factor correlation beyond
1, for one), a residual covariance matrix of the continuous variables
scored from that is not positive definite, or, when none of these holds,
a residual covariance matrix of the structural part of the model that is
not, the conditions lavaan's post-estimation check flags. The regression
weights use the model implied covariance matrix and the Bartlett weights
the loadings and the inverse of the residual covariance matrix, so both
scores inherit the solution. For ordered items lavaan computes no EBM or
ML score for a solution with a negative residual variance, and the
warning and the reason column say so.

## References

Bartlett, M. S. (1937). The statistical conception of mental factors.
*British Journal of Psychology, 28*(1), 97–104.

Bartlett, M. S. (1938). Methods of estimating mental factors. *Nature,
141*(3570), 609–610.

Bollen, K. A. (1989). *Structural equations with latent variables*.
Wiley.

Croon, M. (2002). Using predicted latent scores in general latent
structure models. In G. Marcoulides & I. Moustaki (Eds.), *Latent
variable and latent structure models* (pp. 195–223). Erlbaum.

Curran, P. J., Cole, V., Bauer, D. J., Hussong, A. M., & Gottfredson, N.
(2016). Improving factor score estimation through the use of observed
background characteristics. *Structural Equation Modeling, 23*(6),
827–844.
[doi:10.1080/10705511.2016.1220839](https://doi.org/10.1080/10705511.2016.1220839)

Devlieger, I., Mayer, A., & Rosseel, Y. (2016). Hypothesis testing using
factor score regression: A comparison of four methods. *Educational and
Psychological Measurement, 76*(5), 741–770.
[doi:10.1177/0013164415607618](https://doi.org/10.1177/0013164415607618)

Estabrook, R., & Neale, M. (2013). A comparison of factor score
estimation methods in the presence of missing data: Reliability and an
application to nicotine dependence. *Multivariate Behavioral Research,
48*(1), 1–27.
[doi:10.1080/00273171.2012.730072](https://doi.org/10.1080/00273171.2012.730072)

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

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

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

Rosseel, Y., & Loh, W. W. (2024). A structural after measurement
approach to structural equation modeling. *Psychological Methods,
29*(3), 561–588.
[doi:10.1037/met0000503](https://doi.org/10.1037/met0000503)

Rubin, D. B. (1976). Inference and missing data. *Biometrika, 63*(3),
581–592.
[doi:10.1093/biomet/63.3.581](https://doi.org/10.1093/biomet/63.3.581)

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables. *Structural
Equation Modeling, 16*, 477–497.
[doi:10.1080/10705510903008238](https://doi.org/10.1080/10705510903008238)

Skrondal, A., & Laake, P. (2001). Regression among factor scores.
*Psychometrika, 66*(4), 563–575.
[doi:10.1007/BF02296196](https://doi.org/10.1007/BF02296196)

Thomson, G. H. (1934). The meaning of *i* in the estimate of *g*.
*British Journal of Psychology, 25*, 92–99.

Thurstone, L. L. (1935). *The vectors of mind*. University of Chicago
Press.

## See also

[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
to fit the measurement model;
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
for coefficient *H*;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for the reliability of the unit-weighted composite;
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
before comparing scores across groups;
[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html),
[`lavaan::sam()`](https://rdrr.io/pkg/lavaan/man/sam.html).

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

# Two later spatial tests (t25, t26) were given to the Grant-White
# pupils only, so under full information maximum likelihood each
# Pasteur pupil (rows 1 to 156) is scored from four spatial tests and
# each Grant-White pupil from six.
model <- "
  spatial =~ t1_visual_perception + t2_cubes + t3_paper_form_board +
             t4_lozenges + t25_paper_form_board_r + t26_flags
  verbal  =~ t5_general_information + t6_paragraph_comprehension +
             t7_sentence + t8_word_classification + t9_word_meaning"
fit <- lavaan::cfa(model, data = holzinger_swineford, std.lv = TRUE,
                   missing = "ml")
fs <- factor_scores(fit, data = holzinger_swineford)
fs
#> Factor scores: regression and Bartlett scores, 2 factors (spatial, verbal), 301 rows
#> Fit: ML, missing = "ml"; scored from the fit
#> 
#>  row  factor regression    se           95% CI Bartlett   se           95% CI 
#>    1 spatial     -0.992 0.527  [-2.026, 0.041]   -1.488 0.652 [-2.766, -0.210]
#>    2 spatial     -0.302 0.527  [-1.336, 0.731]   -0.225 0.652  [-1.503, 1.053]
#>    3 spatial     -0.850 0.527  [-1.883, 0.184]   -0.797 0.652  [-2.075, 0.481]
#>    1  verbal     -0.166 0.284  [-0.723, 0.392]   -0.132 0.299  [-0.718, 0.453]
#>    2  verbal     -0.963 0.284 [-1.520, -0.406]   -1.057 0.299 [-1.642, -0.471]
#>    3  verbal     -2.051 0.284 [-2.608, -1.494]   -2.240 0.299 [-2.825, -1.654]
#> (the first 3 rows of each factor)
#> 
#> Determinacy, the correlation of each score with its factor, by missing data pattern; min_r, the lowest correlation two sets of scores equally consistent with the model can have (2 x the regression score's determinacy squared, minus 1; Guttman, 1955):
#>  pattern   n  factor items regression Bartlett min_r
#>        1 145 spatial   6/6      0.888    0.882 0.577
#>        1 145  verbal   5/5      0.959    0.958 0.839
#>        2 156 spatial   4/6      0.850    0.838 0.444
#>        2 156  verbal   5/5      0.959    0.958 0.838
#> 
#> Intervals treat the fitted parameters as known: regression intervals
#> are posterior (empirical Bayes) intervals, Bartlett intervals Wald
#> intervals for each case's factor value. In a later regression, use
#> regression scores for a predictor and Bartlett scores for an outcome,
#> each from its own measurement model (Skrondal & Laake, 2001); for
#> relations among factors, one-step SEM or structural after measurement
#> (lavaan::sam()) is preferred to any regression among scores.
#> 
#> Confidence level: 95%

# One pupil from each school: the spatial scores of the Pasteur pupil
# rest on four tests and have the wider intervals.
fs[fs$row %in% c(1, 157),
   c("row", "factor", ".derived_regression",
     ".derived_regression_ci_lower", ".derived_regression_ci_upper",
     ".derived_bartlett", ".derived_bartlett_ci_lower",
     ".derived_bartlett_ci_upper", "items_observed", "items")]
#>  row factor  .derived_regression .derived_regression_ci_lower
#>  1   spatial -0.992              -2.026                      
#>  157 spatial -1.292              -2.193                      
#>  1   verbal  -0.166              -0.723                      
#>  157 verbal  -0.198              -0.755                      
#>  .derived_regression_ci_upper .derived_bartlett .derived_bartlett_ci_lower
#>  0.041                        -1.488            -2.766                    
#>  -0.391                       -1.729            -2.778                    
#>  0.392                        -0.132            -0.718                    
#>  0.358                        -0.153            -0.739                    
#>  .derived_bartlett_ci_upper items_observed items
#>  -0.210                     4              6    
#>  -0.680                     6              6    
#>  0.453                      5              5    
#>  0.433                      5              5    

# Determinacy and standard errors for each missing data pattern.
attr(fs, "patterns")
#>   group pattern  factor   n items_observed items
#> 1  <NA>       1 spatial 145              6     6
#> 2  <NA>       1  verbal 145              5     5
#> 3  <NA>       2 spatial 156              4     6
#> 4  <NA>       2  verbal 156              5     5
#>                   missing_variables determinacy_regression determinacy_bartlett
#> 1                                                0.8881010            0.8816139
#> 2                                                0.9588085            0.9581213
#> 3 t25_paper_form_board_r, t26_flags              0.8496532            0.8376594
#> 4 t25_paper_form_board_r, t26_flags              0.9587310            0.9581213
#>   min_score_correlation regression_se bartlett_se
#> 1             0.5774469     0.4596483   0.5353493
#> 2             0.8386274     0.2840533   0.2988793
#> 3             0.4438211     0.5273419   0.6520465
#> 4             0.8383304     0.2843146   0.2988793

# The table is ordered by factor and then row, so one factor's column
# lines up with the data.
spatial_bartlett <- fs$.derived_bartlett[fs$factor == "spatial"]
length(spatial_bartlett) == nrow(holzinger_swineford)
#> [1] TRUE

# When scores must stand in for factors in a regression, Skrondal and
# Laake (2001) pair the regression score of the predictor with the
# Bartlett score of the outcome, each from its own measurement model.
# With both factor variances fixed at 1 the slope estimates the latent
# correlation; pairing the regression scores attenuates it. The lm()
# standard errors treat the scores as error free and are not valid.
fit_v <- lavaan::cfa("verbal =~ t5_general_information +
  t6_paragraph_comprehension + t7_sentence + t8_word_classification +
  t9_word_meaning", data = holzinger_swineford, std.lv = TRUE,
  missing = "ml")
fit_s <- lavaan::cfa("spatial =~ t1_visual_perception + t2_cubes +
  t3_paper_form_board + t4_lozenges + t25_paper_form_board_r +
  t26_flags", data = holzinger_swineford, std.lv = TRUE,
  missing = "ml")
fs_v <- factor_scores(fit_v)
fs_s <- factor_scores(fit_s)
c(revised = unname(coef(lm(fs_s$.derived_bartlett ~
                              fs_v$.derived_regression))[2]),
  regression_both = unname(coef(lm(fs_s$.derived_regression ~
                                      fs_v$.derived_regression))[2]),
  latent = lavaan::lavInspect(fit, "cor.lv")["spatial", "verbal"])
#>         revised regression_both          latent 
#>       0.4336799       0.3212779       0.4621871 

# The broom verbs: the table as a plain data frame, and one summary row.
head(generics::tidy(fs))
#>   row group  factor .derived_regression .derived_regression_se
#> 1   1  <NA> spatial          -0.9924726              0.5273419
#> 2   2  <NA> spatial          -0.3023210              0.5273419
#> 3   3  <NA> spatial          -0.8497607              0.5273419
#> 4   4  <NA> spatial           0.6570956              0.5273419
#> 5   5  <NA> spatial          -0.7131244              0.5273419
#> 6   6  <NA> spatial          -0.4098134              0.5273419
#>   .derived_regression_ci_lower .derived_regression_ci_upper .derived_bartlett
#> 1                   -2.0260436                    0.0410985        -1.4876546
#> 2                   -1.3358921                    0.7312500        -0.2251097
#> 3                   -1.8833318                    0.1838104        -0.7966135
#> 4                   -0.3764755                    1.6906667         1.0229176
#> 5                   -1.7466955                    0.3204466        -1.0529895
#> 6                   -1.4433845                    0.6237577        -0.3409758
#>   .derived_bartlett_se .derived_bartlett_ci_lower .derived_bartlett_ci_upper
#> 1            0.6520465                  -2.765642                 -0.2096670
#> 2            0.6520465                  -1.503097                  1.0528779
#> 3            0.6520465                  -2.074601                  0.4813741
#> 4            0.6520465                  -0.255070                  2.3009052
#> 5            0.6520465                  -2.330977                  0.2249981
#> 6            0.6520465                  -1.618963                  0.9370118
#>   items_observed items pattern determinacy_regression determinacy_bartlett
#> 1              4     6       2              0.8496532            0.8376594
#> 2              4     6       2              0.8496532            0.8376594
#> 3              4     6       2              0.8496532            0.8376594
#> 4              4     6       2              0.8496532            0.8376594
#> 5              4     6       2              0.8496532            0.8376594
#> 6              4     6       2              0.8496532            0.8376594
#>   reason
#> 1       
#> 2       
#> 3       
#> 4       
#> 5       
#> 6       
generics::glance(fs)
#>   n_rows n_scored n_factors n_groups n_patterns method_regression
#> 1    301      301         2        1          2        regression
#>   method_bartlett estimator missing scored_from conf_level
#> 1        Bartlett        ML      ml         fit       0.95
```
