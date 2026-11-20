# Categorical Omega for Ordered-Categorical Items, With a Confidence Interval

Estimates categorical omega (\\\omega_C\\; Green & Yang, 2009; Kelley &
Pornprasertmanit, 2016) for a homogeneous composite of
ordered-categorical items and returns a bootstrap confidence interval.
The `denominator` argument selects whether the variance of the summed
categories is built from the unrestricted polychoric correlations
(`"observed"`, the default) or from the polychoric correlations the
fitted single-factor model reproduces (`"model_implied"`, the
coefficient as Green and Yang, 2009, defined it). Incomplete rows are
set aside by listwise deletion (the default) or used through pairwise
deletion or multiple imputation (`missing`); under multiple imputation
the interval is formed by bootstrap then impute (Schomaker & Heumann,
2018).

## Usage

``` r
reliability_omega_categorical(
  data,
  ci_method = c("bca", "percentile", "bootstrap_se", "bootstrap_se_logistic", "none"),
  conf_level = 0.95,
  B = NULL,
  seed = NULL,
  denominator = c("observed", "model_implied"),
  ncpus = 1L,
  missing = c("listwise", "pairwise", "mi"),
  aux = NULL,
  impute = c("pmm", "polr"),
  imputations = 20L,
  progress = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of three or more ordered-categorical
  items (integer codes for the categories), plus any auxiliary variables
  named in `aux`. Only the order of an item's codes matters: its
  categories are scored 0, 1, ..., \\C - 1\\ (see *How the categories
  are scored* in *Details*). A factor whose labels are whole numbers is
  read as those numbers, so a factor labeled 0 to 4 is the item coded 0
  to 4, under every treatment of `missing`; a factor whose labels are
  words is refused, as is an ordered factor whose levels are not in the
  order of the numbers they name. A column read with haven is read as
  its values, with SPSS user-defined missing values made missing.
  Incomplete rows are treated as `missing` says; under the default they
  are set aside. With two items the single-factor model is not
  identified (their one polychoric correlation determines only the
  product of the two loadings), and the function stops.

- ci_method:

  Method for constructing the confidence interval. See *Details*. When
  not supplied, no interval is computed: every interval for categorical
  omega is bootstrap based, and a bootstrap is never run unless
  requested. Ask for `"bca"` (the recommended method) or `"percentile"`
  to obtain one.

- conf_level:

  Confidence level for the interval. Defaults to `0.95`.

- B:

  Number of bootstrap replications. The default, `NULL`, is 10000 under
  `missing = "listwise"` and `"pairwise"`, where BCa, the recommended
  interval, reads its limits from extreme quantiles, and 1000 under
  `missing = "mi"`, where each replication is imputed `imputations`
  times and refitted on each, and the percentile interval is the
  recommended one (see *Run time under multiple imputation*). When no
  item value is missing, `missing = "mi"` has nothing to impute and
  returns the complete data result, whose bootstrap refits once per
  replication, so the default there is 10000, as the message says. When
  a bootstrap is requested, `B` (supplied or default) must be a single
  whole number of at least the count the interval method and the level
  require, as in the rest of the family: at the default `conf_level`,
  400 for `"percentile"` (10 replications beyond each limit), 2000 for
  `"bca"` (50 beyond each nominal limit), and 100 for the two bootstrap
  standard error methods; at 0.99 the first two are 2000 and 10000, so
  the default of 1000 under `missing = "mi"` does not serve a 99 percent
  percentile interval and the refusal names the number. A smaller `B` is
  refused before anything is computed (class `dmar_count_invalid`). When
  fewer replications than required return a value the interval is
  withheld under a warning of class `dmar_bootstrap_insufficient`. See
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md),
  which also describes `options(dmar.bootstrap_exploratory = TRUE)` for
  a first look with a smaller `B`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, under which the resamples are drawn from, and advance, the
  current random number stream; supply an integer for a reproducible
  interval, which is set for the call only, the caller's random number
  state being restored when it returns.

- denominator:

  How the polychoric correlations behind the total variance in the
  denominator are estimated: `"observed"` (the default; the unrestricted
  polychoric correlations) or `"model_implied"` (the correlations
  reproduced by the fitted single-factor model, as in Green and Yang,
  2009). See *Details*. Every interval method applies to both settings,
  since every interval here is bootstrap based and each replication
  recomputes the coefficient with the same denominator.

- ncpus:

  The number of processes over which the bootstrap replications are
  spread. The default `1` runs them in this R session; a larger number
  uses forked workers (boot's `parallel = "multicore"`) or, on Windows,
  which cannot fork, a socket cluster (`"snow"`), whose workers load the
  installed DMAR. boot draws every resample before any replication runs,
  so a given `seed` returns the same interval, counts, and replicates
  whatever `ncpus` is. A session that forks workers and also starts
  background processes through callr or processx needs the environment
  variable `PROCESSX_NOTIFY_OLD_SIGCHLD` set before processx is loaded
  (see
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)).
  Ignored unless a bootstrap method is selected.

- missing:

  How incomplete rows of `data` are treated: `"listwise"` (the default;
  only complete rows are analyzed), `"pairwise"` (each threshold and
  polychoric correlation from the rows that observe its items), or
  `"mi"` (multiple imputation, with a bootstrap then impute interval).
  Listwise and pairwise deletion are consistent only when values are
  missing completely at random; multiple imputation when they are
  missing at random given the items and `aux`. See *Missing data* in
  *Details*.

- aux:

  Optional character vector naming columns of `data` that are auxiliary
  variables: numeric predictors in the imputation model, never items, so
  they inform the imputations without entering the composite (the
  inclusive strategy of Collins, Schafer, & Kam, 2001). Supplying `aux`
  implies `missing = "mi"`, as it implies full information maximum
  likelihood in
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md);
  with `"listwise"` or `"pairwise"` named explicitly it is an error.

- impute:

  The imputation model under `missing = "mi"`: `"pmm"` (the default;
  predictive mean matching) or `"polr"` (the proportional odds model).
  See *The imputation model* in *Details*. Supplying it with another
  `missing` is an error.

- imputations:

  Under `missing = "mi"`, the number of imputations of each bootstrap
  replication (at least 2, default 20), which is also the least number
  of imputations behind the point estimate (whose number is at least 20
  and may be raised by von Hippel's rule; see *The number of
  imputations*). Supplying it with another `missing` is an error.

- progress:

  `NULL` (the default), or a function of two arguments that the
  bootstrap calls as it runs, so that an application can show how far a
  long bootstrap has gone and let its user stop it (a Shiny app, for
  example). It is called as `progress(done, total)`, with the number of
  replications attempted so far and `B`: after each replication when
  `ncpus = 1`, and with more processes after each block of replications
  they complete together, about a fiftieth of `B`
  (`ncpus * ceiling(B / (50 * ncpus))` replications, the last block
  possibly fewer), so that `done` advances by the block. A BCa interval
  takes its acceleration from a regression on the resampling frequencies
  when at least as many replications return a value as there are cases,
  and from the jackknife, one leave-one-out refit per case, when fewer
  do or when that regression fails or gives influence values that are
  not all finite; the jackknife runs at most once in a call, and not at
  all when no interval is formed or every returned replicate has the
  same value. While it runs, the callback is called again after each
  refit (with more processes, after each block of
  `ncpus * ceiling(n / (50 * ncpus))` refits, `n` the number of cases
  resampled), with `done` and `total` both `B`, so it can be called more
  than `B` times, and `done` equal to `total` does not mean that the
  bootstrap has finished. If the callback returns `FALSE`, the bootstrap
  stops before the next replication or block with an error of class
  `dmar_bootstrap_cancelled`, whose elements `done` and `total` say how
  far it got (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md));
  any other value is ignored. An error inside the callback stops the
  call with an error of class `dmar_progress_error`, whose message names
  the callback and which carries `done`, `total`, and, as `parent`, the
  callback's own condition with its class intact, so that an application
  can tell a failing callback from a failing analysis. The callback runs
  outside the random number stream, so a given `seed` returns the same
  result with or without it, and a callback that draws random numbers
  leaves the caller's stream as it found it. With `NULL` nothing is
  called. Ignored unless a bootstrap method is selected.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample \\\omega_C\\), `"se"` (the bootstrap standard deviation across
replications, already on the coefficient scale; `NA` for
`ci_method = "none"`), `"ci_lower"` and `"ci_upper"` (clamped to \[0,
1\]), `"conf_level"`, `"N"` (the rows analyzed: the complete rows under
listwise deletion, otherwise the rows with at least one observed item,
and under multiple imputation with `aux` also the rows with an observed
auxiliary variable), `"N_complete"` (the rows complete on the items;
equal to `"N"` under listwise deletion, and carried so the whole
reliability family returns one shape), and `"n_items"` (number of
items); the rows are the same under every `missing`. Under
`missing = "mi"` the estimate is the mean over the imputations of the
observed data and `se` is the standard error of that mean from the
bootstrap (see *The number of imputations*). Attributes `coefficient`
(`"omega_categorical"`), `denominator` (`"observed"` or
`"model_implied"`), `missing` (`"listwise"`, `"pairwise"`, or `"mi"`;
`"listwise"` also when the items are complete and another treatment had
nothing to do), `aux` (when supplied), `ci_method` (under
`missing = "mi"` the method followed by `"(bootstrap then impute)"`),
and `B` record the computation. The attribute `residual_cor` carries the
local dependence summary of the single-factor fit: a list with `matrix`
(the unrestricted polychoric correlations minus those the model
reproduces, named by item), `largest` (the residual correlation largest
in absolute value, with its sign), `items` (the pair it belongs to), and
`saturated`; a pair whose polychoric correlation the data do not
determine has an `NA` residual and is passed over. With three items the
single-factor model is saturated (zero degrees of freedom) and
reproduces every polychoric correlation, so its residuals are zero and
carry no information about local dependence: `saturated` is then `TRUE`,
`largest` is `NA`, and `items` is `NULL`, rather than naming whichever
pair carries the largest rounding error. The printed table names the
denominator and the largest residual correlation (or says that the model
is saturated) in two lines under the table. A bootstrap call also
carries the attribute `bootstrap`, a list recording the replications the
interval conditions on: `attempted` (`B`), `returned` (the replications
whose value the interval uses; the standard error is their standard
deviation and the percentile and BCa limits are their quantiles),
`failed` (the difference), `reasons` (a named count of why replications
were set aside: `error` when the refit raised an error, `non_finite`
when it was set aside without one, as a fit that lavaan stopped or that
did not converge is, as an improper one is, and, under the observed
denominator, as one that left a polychoric correlation undetermined is,
and, under either, as one whose determined correlations leave the
loadings unidentified is), `error_messages` (the distinct error messages
with their counts, empty when none), `undetermined` (the replications
set aside because the resample left a polychoric correlation
undetermined, which the observed setting needs, or left correlations
undetermined whose absence leaves the loadings unidentified, under
either setting; see *Details*), `improper` (the replications set aside
because their single-factor solution was improper, with a standardized
loading above one, which is an error variance below zero),
`improper_kept` (always zero here: categorical omega is not defined for
an improper solution, so such a replication has no estimate to keep),
`improper_rule` (`"drop"`, for the same reason), `no_solution` (always
zero here: it counts the replications that
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
sets aside because the correlations of three items admit no
single-factor solution, a check the categorical fit does not make),
`not_converged` (zero here: a fit that did not converge is counted among
the `non_finite` ones), `too_few_cases` (always zero here: it counts the
replications that
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
sets aside for too few distinct cases, a rule the categorical fit does
not apply), `caution` (`TRUE` when the replications set aside number at
least \\(1 - C)/2\\ of `B`, with \\C\\ the confidence level, so that
where they would have fallen can decide a limit), `require_half`
(`TRUE`: no interval is formed when fewer than half of the replications
return a value), `required`, `bootstrap_status`, and `B_used` (the
fewest replications a request may ask for, the verdict of the rule, and
the count the limits rest on, as in
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)),
`acceleration` and `acceleration_method` (the BCa acceleration used and
whether it came from the regression or the jackknife; `NA` unless a BCa
interval was formed), `jackknife` (the same `attempted`, `returned`,
`failed`, `reasons`, `error_messages`, `undetermined`, `improper`,
`improper_kept`, `no_solution`, `not_converged`, and `too_few_cases` for
the leave-one-out refits behind the acceleration, all zero when no
jackknife ran), `ncpus` (the number of processes used), `t` (the
estimate of every replication, in order, `NA` where there is none), and
`status` (a factor giving the fate of every replication, with levels
`"proper"`, `"improper"`, `"error"`, `"non_finite"`, `"undetermined"`,
and `"not_converged"`; the interval uses those whose status is
`"proper"`). The attribute `B_used` is the number of replications the
interval uses (`returned`; when fewer than half returned a value, or
fewer than `required`, no interval is formed). The attributes
`n_dropped`, `B_required`, and `bootstrap_status` repeat the counts in
the names every bootstrap result of the package carries (see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)).
The printed table shows these counts in one line, and a warning states
them whenever a replication was set aside or the caution applies.

Under `missing = "mi"` the `bootstrap` attribute has the same elements,
with three differences: `reasons` has a third count, `imputation`, for
the replications whose imputation itself failed (not part of
`non_finite`), whose messages are in `imputation_messages`; `status` has
a seventh level, `"imputation"`; and a replication counted under another
reason is one in which fewer than half of its completed data sets
returned a value, the first that returned none deciding the reason. The
jackknife counts are zero, since BCa is refused. The attribute
`imputation` is a list: `impute`; `imputations_point` (\\M_0\\, the
imputations of the observed data behind the point estimate);
`imputations_returned` (those of them that returned a value, which the
estimate averages); `imputations_set_aside` (the others, counted by
reason: `improper`, `undetermined`, `error`, `non_converged`, and
`non_finite`); `imputations_replicate` (`imputations`);
`imputations_floor` (`max(imputations, 20)`); `imputations_rule` (the
number von Hippel's rule asks for); `maxit` (10, the iterations of the
chained equations); `estimates` (the \\M_0\\ completed data estimates,
`NA` where one was set aside), `status` (the fate of each: `"proper"` or
the reason it was set aside), and `between_sd` (the standard deviation
of those that returned a value); `sd_replicates` (the standard deviation
of the replicate means); `s2_boot`, `s2_imp`, `fmi`, and `fmi_ci_upper`
(the sampling and imputation variances, the fraction of missing
information, and its 95% upper limit); `widening` (\\\sqrt{1 + \gamma /
M_b}\\, the factor by which imputation noise widens the percentile
interval); `replicate_estimates` (the `B` by `imputations` matrix of
completed data estimates, `NA` where one returned no value; `NULL`
without a bootstrap); `replicate_sets_set_aside` (the completed data
sets set aside within the replications that returned a value, counted by
the reasons of `imputations_set_aside`; `NULL` without a bootstrap);
`polr_fallbacks` (how often the proportional odds model fell back to a
multinomial model, in the point estimate's imputations and in the
replications); `ridge` (how often mice applied its ridge penalty,
counted the same way); and `aux_set_aside` (the auxiliary variables left
out of the imputation model as linear functions of other predictors,
`character(0)` when none). The components that need the bootstrap are
`NA` without one. Under `missing = "pairwise"` or `"mi"` the printed
table adds a line naming the treatment and its assumption, and under
`"mi"` a line giving the number of imputations, the fraction of missing
information, and the widening (none when no interval is reported).

## Details

Categorical omega is designed for items measured on an ordered
categorical scale (e.g., Likert items). It uses a probit-link
single-factor model in which each observed item \\X_j\\ is modeled as a
categorization of an underlying continuous response variable
\\X_j^{\*}\\ via thresholds \\t\_{j,c}\\ (Muthén, 1984; Millsap &
Yun-Tein, 2004), fit by diagonally weighted least squares with mean- and
variance-adjusted test statistic (WLSMV), the standard estimator for
ordered categorical items. With the delta parameterization
(\\Var(X_j^{\*}) = 1\\), the population categorical omega is
\$\$\omega_C = \frac{\sum\_{j=1}^{J} \sum\_{j'=1}^{J}
\sigma\_{jj'}\\\left(\lambda\_{j}\lambda\_{j'}\right)}{\sum\_{j=1}^{J}
\sum\_{j'=1}^{J} \sigma\_{jj'}\\\left(\rho\_{X\_{j}^{\*}
X\_{j'}^{\*}}\right)},\$\$ where \\\sigma\_{jj'}(r)\\ is the covariance
of \\(X\_{j}, X\_{j'})\\ computed from a bivariate normal CDF over pairs
of category thresholds and a latent-response correlation \\r\\ (Green &
Yang, 2009, Eq. 13–14, with the full coefficient their Eq. 21; Kelley &
Pornprasertmanit, 2016, Eq. 17–18). The numerator uses the correlations
of the common part, \\\lambda_j \lambda\_{j'}\\, under either setting of
`denominator`; the two settings differ in the correlations \\\rho\\ used
for the total variance. The bivariate normal probabilities, one for each
pair of thresholds of each pair of items, are computed by DMAR itself,
by the method of Drezner and Wesolowsky (1990) as refined by Genz
(2004), which is the algorithm mvtnorm applies in two dimensions. The
two agree to within \\4 \times 10^{-16}\\, except within about
\\10^{-7}\\ of a correlation of one, where mvtnorm loses accuracy and
DMAR's computation does not. Computed together, the probabilities take a
small fraction of the time of the factor model's fit.

- `"observed"` (default):

  The unrestricted polychoric correlations, estimated from the data by a
  saturated bivariate model (Kelley & Pornprasertmanit, 2016, Eq. 18).
  Structure that the single-factor model omits, such as a pair of items
  with shared wording, stays in the total variance, which makes this the
  categorical counterpart of robust omega (see
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)).
  It is the coefficient `MBESS::ci.reliability(type = "categorical")`
  computes, apart from the samples described below in which the data do
  not determine a polychoric correlation.

- `"model_implied"`:

  The latent-response correlations reproduced by the fitted
  single-factor model. This is categorical omega as Green and
  Yang (2009) defined it: their appendix program takes “the reproduced
  (model-based) correlation matrix” as the input to the denominator, so
  correlated errors enter the total variance only when they are
  specified in the model.

When the single-factor model reproduces the polychoric correlations the
two settings coincide, and they diverge as the model misfits. Neither
protects the numerator: shared structure among a few items raises their
loadings under either setting. The result therefore reports the residual
polychoric correlations (the unrestricted minus the reproduced) and
prints the largest under the table, so that a high coefficient resting
on redundant items is visible beside it. When the shared structure sits
on the items that dominate the factor, the factor absorbs it and the
residuals appear among the other items, so the largest residual signals
local dependence more reliably than it locates it.

**What each denominator describes.** Both settings share the numerator,
the variance of the summed categories attributable to the single factor,
and differ in what they take the total variance of the sum to be. The
model implied setting takes the variance the fitted single-factor model
reproduces, so the coefficient describes the sum score as the model
describes it: covariance the model omits is absent from the numerator
and from the denominator alike. The observed setting takes the variance
implied by the polychoric correlations estimated without the
single-factor restriction (still under the latent-response model that
defines polychoric correlations), so covariance the model omits is
counted as variance of the sum that the factor does not explain. The two
estimate the same quantity when the model holds. Kelley and
Pornprasertmanit (2016) used the observed setting because the total
variance of the sum score is a property of the data and need not be
taken from a model that may be wrong; Green and Yang (2009) used the
model implied setting, which is the reliability of the sum score under
the fitted model.

**How the categories are scored.** Categorical omega is the reliability
of a sum in which each item's categories are scored as consecutive
ranks, 0, 1, ..., \\C - 1\\, in the order of the item's codes: Green and
Yang (2009, Eq. 11) and Kelley and Pornprasertmanit (2016, Eq. 16)
define the observed item as \\X\_{ij} = c\\ for \\c = 0, 1, \ldots, C -
1\\, and MBESS and semTools compute the same coefficient. The function
therefore uses the order of the codes in `data`, not their values, and
the coefficient is that of the sum of the codes as given whenever every
item's codes are equally spaced with one step shared by all the items
(codes 0, 2, 4, 6 on every item give the same result as 1, 2, 3, 4).
When they are not, as with codes -2, -1, 1, 2, with items stepping by 2
beside items stepping by 1, or with a category between two used ones
that no respondent chose, the sum of the codes is a different composite,
and a message says so and names the items; the coefficient is still that
of the ranks. In the examples behind this paragraph the coded sum's
reliability was lower by 0.009 to 0.012 for codes -2, -1, 1, 2 on every
item and by 0.034 for two three-category items coded 1, 3, 5 beside two
coded 1, 2, 3, and an unused interior category changed it by about
0.001.

**When the data do not determine a polychoric correlation.** The
observed setting needs every unrestricted polychoric correlation. When
two items are skewed in opposite directions, as in Green and Yang's
(2009) mixed-skew conditions, their bivariate table can have empty cells
in a pattern (every respondent below the top category of one item in the
bottom category of the other) that leaves the polychoric correlation
without a finite standard error. The unrestricted value is then
arbitrary, and the single-factor fit gives the correlation no weight, so
the loadings rest on the correlations the data do determine. Those
identify the loadings only when they connect every item to the others,
directly or through other items, and include at least one cycle of an
odd number of items, such as three items whose three correlations are
all determined (the graph whose edges are the determined pairs must be
connected and not bipartite); an item joined to the rest by a single
determined correlation is then identified as well. When every
correlation between `y1` or `y2` and `y3` or `y4` of the example below
is undetermined, the determined ones fix only the products of the
loadings of `y1` and `y2` and of `y3` and `y4`: the relative sign of the
two pairs is free, and lavaan's default start values can return loadings
of opposite signs and a categorical omega near 0.05, where the
population value is 0.716. In 60 samples of 300 from that design, 13 had
determined correlations that did not identify the loadings: 2 with that
sign reversal (0.05 and 0.09), 8 whose value moved with the start values
by up to 0.06, and 3 whose default start values gave an improper
solution. The function therefore checks the graph. When the determined
correlations identify the loadings, the observed setting stops with an
error naming the pair and pointing to the model implied setting, and the
model implied setting proceeds, says which pair it set aside, and
reports that pair's residual correlation as `NA`. When they do not, both
settings stop with an error naming the pairs and the blocks of items no
determined correlation joins, since the numerator takes the loadings
under either setting. The unrestricted polychoric correlations come from
the single-factor fit itself, which agrees with a second, saturated fit
whenever every correlation is determined. MBESS computes them from such
a saturated fit, which returns 0 for an undetermined correlation and so
removes the pair from the total variance without a message. The example
below uses Green and Yang's design at a sample size where samples rarely
have this problem; at \\N = 500\\ most of them do.

A bootstrap meets the problem far more often than the sample does. The
few respondents in the sparse cells of an oppositely skewed pair are
what determine its polychoric correlation, and a resample omits some of
them as a matter of course, so a sample whose correlations are all
determined still yields many resamples in which one is not. On the data
of the example below, which have no such pair, about half of the
replications of a bootstrap of the observed setting have one (80, 95,
and 97 of 200 in percentile bootstraps with `seed` 1, 2, and 3). Those
replications are set aside, counted in the `bootstrap` attribute
(`undetermined`; see *Value*), named in the printed line and in a
warning, and the interval then conditions on the replications that
redrew the respondents who determine the correlation, which is not the
bootstrap distribution of the coefficient; with more than half set
aside, no interval is reported (see below), and on these data the
observed setting is always near that point. The model implied setting
needs no unrestricted polychoric correlation and sets aside only the
replications whose determined correlations leave the loadings
unidentified (counted as `undetermined` as well), so with items skewed
in opposite directions it keeps more of them, but how many it keeps
depends on the sample size. In a percentile bootstrap of the example's
complete data (\\N = 1000\\, `B = 200`, `seed = 1`) it set aside 1
replication of 200, and over 20 samples of that size from the design
(`B = 50`) none in 15 and at most 18%. Under pairwise deletion of the
example's incomplete copy it set aside 47 of 200: 27 whose loadings were
unidentified, 3 improper, and 17 whose fit failed, so the caution
applies. At \\N = 300\\, in the 45 samples of 60 whose point estimate
was identified and proper, it set aside a median of 42% of the
replications (`B = 50`): 18% where the sample determined every
correlation, 40% where it left one undetermined, and 64% where it left
two, and more than half in 15 of the 45. When fewer than half of the
replications return a value, no interval is reported (see below). Under
either denominator, a replication whose determined correlations do not
identify the loadings is set aside rather than evaluated at arbitrary
starting values.

**When the solution is improper.** With few respondents in the extreme
categories, or an item that nearly duplicates the factor, the
single-factor fit can return a standardized loading above one, a
negative residual variance (a Heywood case). Categorical omega is not
defined for such a solution, since its assembly evaluates the bivariate
normal distribution at the item's communality, and the function stops
with an error naming the item and its loading under either setting of
`denominator`. A bootstrap replication with an improper solution is set
aside and counted (`improper` in the `bootstrap` attribute); these are
the replications in which a standardized loading ran past one (an error
variance fell below zero), so setting them aside can trim the upper tail
of the bootstrap distribution, which is why their number is reported.

Kelley and Pornprasertmanit (2016) found in extensive Monte Carlo
simulation that the bias-corrected and accelerated (BCa) bootstrap
confidence interval for categorical omega achieved acceptable coverage
across a wide variety of threshold patterns, sample sizes, item counts,
and population reliability values, with one documented exception:
coverage dipped somewhat below the acceptable range when the number of
items and the population reliability were both high. They specifically
recommend BCa for categorical omega. Because no bootstrap runs in DMAR
unless the user requests one, the default output is the point estimate
with a message naming the call that produces the recommended interval;
request `ci_method = "bca"` to obtain it. The interval is not a quick
one. Every replication refits the ordered-categorical factor model,
whose sample statistics also supply the unrestricted polychoric
correlations behind the denominator. That is a fraction of a second for
a scale of moderate length, so the recommended interval at the default
`B = 10000` can run for several minutes; when `B` is below the number of
cases, the jackknife behind the BCa acceleration adds one such refit per
case (at the default `B` the acceleration is read from the replicates
themselves at no extra cost; see the `"bca"` entry below). The example
on this page therefore stops at the point estimate; the reported
analysis adds `ci_method = "bca"` to the same call, with `seed` for
reproducibility, and for these oppositely skewed items it bootstraps the
model implied setting, which at this sample size keeps nearly all of the
replications of the complete data (see above), or reports how many
replications of the observed setting were set aside.

**When to use.** Use `reliability_omega_categorical` when items are
ordered-categorical, especially when (a) the number of categories is
small (e.g., two to five), (b) item distributions are skewed, or (c)
threshold patterns differ markedly across items. In Kelley and
Pornprasertmanit's (2016) Study 3, treating ordered items as continuous
and using
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
with the observed total variance in the denominator achieved acceptable
coverage only when threshold patterns were similar across items; in
their experience that condition is rare in practice.

Available confidence interval methods (set via `ci_method`). Every
interval here is bootstrap based: the rows of `data` are resampled with
replacement `B` times (10000 by default, 1000 under `missing = "mi"`;
see *Missing data* below) and categorical omega is recomputed, with the
full WLSMV model refit, on each replication (Efron & Tibshirani, 1993).
Replications whose refit fails or does not converge, most common with
small samples and sparse response categories, are set aside, as are
replications with an improper solution, replications that leave a
polychoric correlation undetermined under the observed setting, and
replications whose determined correlations leave the loadings
unidentified under either setting (all described above), and the
interval is computed from the replications that return a value: the
standard error is their standard deviation, and the percentile and BCa
limits are quantiles of the returned replicates, not of all `B`. How
many were attempted, how many returned a value, and why the rest were
set aside is recorded in the `bootstrap` attribute of the result (see
*Value*), printed in one line under the table, and stated in a warning
whenever a replication was set aside, so an interval that conditions on
a subset of the replications says so and by how much. When fewer than
half of the replications return a value, no interval is reported: the
standard error and the limits are `NA`, and the warning and the printed
line say why. The rest would be a minority of the resamples, chosen by
whether they returned a value, and their spread would describe those
resamples rather than the sampling distribution of the coefficient; this
is the criterion by which the point estimate under `missing = "mi"`
stops when fewer than half of its completed data sets return a value
(Lim & Cheung, 2022, p. 1069). The remedy is a larger sample, which
fills the sparse cells of the bivariate tables. Bootstrap results vary
from run to run; supply `seed` for reproducibility.

- `"bca"`:

  The bias-corrected and accelerated bootstrap, the default and the
  specific recommendation of Kelley and Pornprasertmanit (2016). Where
  the percentile interval reads its limits directly off the empirical
  quantiles of the bootstrap estimates, BCa adjusts the two quantile
  positions for median bias (estimated from the bootstrap distribution)
  and for the rate at which the estimator's variance changes with the
  parameter (the acceleration, estimated by the jackknife), making it
  second-order accurate where the percentile interval is first-order
  accurate (DiCiccio & Efron, 1996). The adjusted quantile positions sit
  farther into the tails than the percentile interval uses, which is why
  the default `B = 10000` is larger than the customary 2000; reduce `B`
  for exploration, not for a reported analysis. The acceleration is
  estimated as boot estimates it: by regressing the bootstrap replicates
  on the resampling frequencies when at least as many replications as
  cases returned a value, at no extra cost, and otherwise, or when that
  regression fails, by the jackknife, one leave-one-out refit per case.
  A leave-one-out refit that fails (when a case is the only respondent
  in one of an item's categories, for example) is counted and set aside,
  and the acceleration is estimated from the refits that returned a
  value. When the acceleration cannot be estimated at all (fewer than
  two leave-one-out refits returned a value), or the bias correction
  cannot be (every returned replicate on one side of the estimate), the
  percentile interval is returned in its place with a warning, and the
  `ci_method` attribute records the substitution as
  `"percentile (BCa acceleration could not be estimated: ...)"`. The
  counts of leave-one-out refits and the acceleration used are in the
  `bootstrap` attribute.

- `"percentile"`:

  Percentile bootstrap: the interval limits are the empirical quantiles
  of the bootstrap estimates.

- `"bootstrap_se"`, `"bootstrap_se_logistic"`:

  Wald intervals using the bootstrap standard deviation as a standard
  error, built on the logit scale for the `_logistic` variant so the
  endpoints respect \[0, 1\]. Under `missing = "mi"` the standard error
  is that of the point estimate (see *The number of imputations* below).

- `"none"`:

  Return only the point estimate.

Under `missing = "mi"` the percentile interval is the recommended one
and BCa is refused (see *Missing data*).

**Missing data.** The `missing` argument sets how incomplete rows of
`data` are treated. Under `"pairwise"` and `"mi"`, rows with no observed
item, which carry no information about the items, are set aside, as
under full information maximum likelihood elsewhere in the family, and
when no item value is missing in the rows that remain both return the
complete data result, with a message. Under multiple imputation a row
with no observed item is kept when an auxiliary variable is observed on
it, and its items are imputed from the auxiliaries, as full information
maximum likelihood elsewhere in the family keeps such a row: setting it
aside would make the analysis conditional on having an observed item,
which biases the estimate when the items go missing together at random
given the auxiliaries.

- `"listwise"` (default):

  Only the rows observed on every column of `data` are analyzed.
  Listwise deletion is consistent when values are missing completely at
  random (MCAR), and even then it discards what the incomplete rows say
  about the items they observe. When fewer than two rows are complete,
  or an item takes a single value in the complete rows, the function
  stops, says so, and names the treatment that uses the incomplete rows.

- `"pairwise"`:

  Each item's thresholds are estimated from the rows that observe it and
  each polychoric correlation from the rows that observe both items
  (lavaan's `missing = "pairwise"`); the single-factor fit, both
  denominators, and the assembly are otherwise unchanged. Pairwise
  deletion is the default treatment of missing data under WLSMV
  estimation in Mplus (Shi, Lee, Fairchild, & Maydeu-Olivares, 2020).
  Asparouhov and Muthén (2010) show it consistent when the missingness
  depends only on covariates in the model; a single-factor model of the
  items has no covariates, so here it is consistent only under MCAR, the
  assumption of listwise deletion. It uses more of the data than
  listwise deletion, not a weaker assumption: when values are missing at
  random (MAR), Shi et al. (2020) found the estimates of pairwise
  deletion severely biased, and they recommend multiple imputation.
  Every interval method applies, each replication resampling the rows
  and refitting by pairwise deletion.

- `"mi"`:

  Multiple imputation, consistent when values are MAR given the items
  and any auxiliary variables named in `aux`. The items are imputed by
  chained equations with mice (van Buuren & Groothuis-Oudshoorn, 2011),
  with every other item and every auxiliary variable as a predictor and
  10 iterations of the chains, and categorical omega is computed on each
  completed data set with the chosen `denominator`; the point estimate
  is their mean (Rubin's rule for a point estimate), and the residual
  polychoric correlations are averaged the same way. An auxiliary
  variable is a predictor in the imputation model only, never an item,
  so it carries information about the missing values without entering
  the composite. The interval is the bootstrap then impute interval of
  Schomaker and Heumann (2018, their Boot MI): `B` bootstrap samples are
  drawn from the incomplete rows, each is imputed `imputations` times,
  the estimate of a replication is the mean of its imputations, and the
  percentile interval of those means is the confidence interval. Because
  the imputation model is refitted on every bootstrap sample, the
  interval carries the sampling variability of the imputation model as
  well as that of the coefficient. Bartlett and Hughes (2020) found that
  imputing and then bootstrapping does not, in general, give valid
  variances when the imputation and analysis models are uncongenial, as
  a regression imputation model and a probit factor model are, while
  bootstrapping and then imputing can. Schomaker and Heumann found the
  coverage of Boot MI near the nominal level in their simulations, and
  found that pooling the estimates of every imputation of every
  replication into one distribution gives intervals that are too wide,
  which DMAR therefore never does. A completed data set of the observed
  data that returns no value is set aside and counted, and the point
  estimate averages the rest; the call stops only when fewer than half
  return a value. A replication computes the same statistic: the mean of
  those of its completed data sets that return a value, set aside when
  fewer than half of them return one or when its imputation fails (see
  *When a completed data set is improper* below).

**When a completed data set is improper.** Imputation can produce an
improper single-factor solution (a standardized loading above one) in a
completed data set when the observed data have none, and the example
below shows how. The polychoric correlation of an item piled into its
top category with one piled into its bottom category is determined by
the few respondents in the cells where the two disagree, and the values
missing at random given `y1` and `y2` are missing most often among the
respondents low on them, who supply those cells. Predictive mean
matching fills each missing value with the observed value of a donor
whose predicted mean is close, and among the observed respondents low on
`y1` and `y2` nearly all chose the bottom category of `y3` (all but one,
in the example's sample). How many imputed respondents land in the rare
cells therefore varies from one completed data set to the next; where
almost none do, the polychoric correlation runs toward one and the
loading of `y3` past it. In 100 completed data sets of the example's
data, 13 were improper: 7 of the 13 in which the imputations put no
respondent in the cells where `y3` is above its bottom category and `y1`
or `y2` below its top, 6 of the 58 in which they put one or two there,
and none of the 29 in which they put three or more. The standard
deviation of the loading of `y3` over those imputations was 0.106,
against 0.046 over 40 imputations under the proportional odds model.
Neither that model nor more donors is a remedy: over 61 samples of the
example's design the proportional odds model gave improper sets in fewer
samples (5 of the 61, against 13) but estimates lower than the complete
data value by 0.066 on average, and predictive mean matching with ten
donors in place of five lowered the average share of improper sets only
from 3% to 2%.

Categorical omega is not defined for an improper solution, so such a
completed data set has no estimate, and Rubin's rules, which average an
estimate from every completed data set, need one. Stopping whenever one
completed data set fails, as Jia and Wu (2019) counted a replication of
their simulation as failed unless all its imputations converged, makes
the method fail by design: with a share \\p\\ of such sets, \\M\\
imputations include one with probability \\1 - (1 - p)^M\\, which
approaches one as von Hippel's rule asks for more imputations. The
function instead sets such a completed data set aside, as it sets aside
a bootstrap replication, and the point estimate averages those that
returned a value, as Lim and Cheung (2022) pooled the imputations that
converged and as lavaan.mi does by default for imputations that did not
converge. Setting aside conditions the estimate on a property of the
estimate, so the number set aside, and why, is reported in a warning and
printed under the table, and the call stops when fewer than half of the
completed data sets return a value, the criterion by which Lim and
Cheung counted an analysis as failed. The improper sets are those with
the largest loadings, so setting them aside lowers the estimate. Over
the 61 samples, 13 had at least one improper completed data set among
20, with a share of 3% on average and 65% at most, and no sample's
complete data or pairwise deletion fit was improper. Against refitting
those sets with every loading bounded at one, setting them aside lowered
the estimate by about 0.034 times the share set aside (0.0017 at 5%, a
tenth of the complete data estimate's sampling standard deviation of
0.0165), and over the 13 samples the two were equally far from the
complete data estimate (root mean square differences 0.0126 and 0.0122).
No completed data set of the observed data is fitted after the
bootstrap, so none can discard a bootstrap already run (see *The number
of imputations*).

A bootstrap replication follows the same rule, set by set: it averages
those of its \\M_b\\ completed data sets that return a value and is set
aside when fewer than half do, because the bootstrap distribution is
that of the statistic computed as the point estimate is computed. The
interval follows it too: when fewer than half of the replications return
a value, none is reported, as for the other treatments (see the
confidence interval methods above). Setting a replication aside whenever
any one of its completed data sets returned no value would be a
different statistic, and the argument above applies to it unchanged: a
replication would be set aside with probability \\1 - (1 - p)^{M_b}\\,
so raising `imputations` would set more replications aside, and the
completed data sets of those replications that did return a value would
be discarded. In a bootstrap of the example's incomplete data under the
model implied setting (`B = 200`, `imputations = 5`, `seed = 1`), in
which a quarter of the completed data sets returned no value, that rule,
applied to the first two, three, four, and five imputations of each
replication, would have set aside 35.5, 45.5, 48.5, and 52.5 percent of
the replications, where the rule used sets aside 12.5, 21, 14.5, and 20
percent. Under the rule used the share set aside goes to zero as \\M_b\\
grows when fewer than half of the completed data sets fail, but not
steadily: an odd \\M_b\\ needs \\(M_b + 1)/2\\ of its sets to return a
value, one more than the even number below it needs of one set fewer, so
it sets aside more replications than either even neighbor (with a fifth
of the sets failing independently, 4.0, 10.4, 2.7, and 5.8 percent at
two to five imputations). Where completed data sets are set aside, an
even number of imputations is the better choice, and the number the
warning names is even. The warning and the `imputation` attribute count
the completed data sets set aside within the replications kept
(`replicate_sets_set_aside`).

**The imputation model.** `impute = "pmm"`, the default, is predictive
mean matching: each missing value is replaced by the observed value of
the same item of a donor whose predicted mean, from a linear regression
on the other variables, is close, with the regression parameters drawn
anew for each imputation. It never imputes a category the item does not
show, though a category rare among the observed values is rarely drawn.
`impute = "polr"` is the proportional odds model, with every other item
entering as a factor (a binary item is imputed by logistic regression,
the same model for two categories). mice draws its categories at the
model's point estimates, without drawing the model's parameters, and Jia
and Wu (2019), pooling its imputations in the usual way, found
acceptable loading estimates but confidence intervals that were too
narrow under missingness at random; under bootstrap then impute the
refit on every bootstrap sample supplies the variability of the
imputation model instead. The default rests on a Monte Carlo simulation
run for this function: 10 items, one factor with loadings from .50 to
.80, \\n = 300\\, and items 6 to 10 missing at random given the mean of
items 1 to 5, about 21% of the values in all. Against the complete data
estimate of the same sample, the mean categorical omega was off by
-0.0035 with five categories (100 replications) and by +0.0015 with
three skewed categories (60 replications) under predictive mean
matching; by -0.019 and +0.035 under the proportional odds model; by
-0.028 and -0.047 under pairwise deletion; and by -0.105 with five
categories under listwise deletion, which with three categories failed
in 55 of the 60 samples. A replication under predictive mean matching
was also about ten times faster (see *Run time* below). The evidence is
one factor structure and one missing data mechanism, so with very rare
categories or another mechanism the comparison could differ.

Every incomplete item is imputed. By default mice stops imputing a
variable whose correlation with another, computed from the rows that
observe both, reaches .999 in absolute value; two items observed
together on only two or three rows, as in a planned missing design with
a small overlap, reach it by arithmetic, and the item left with its
missing values would then be fitted by listwise deletion inside each
completed data set. That check is therefore turned off, an item with a
single observed value is refused before anything is imputed, and an item
mice still does not impute stops the call with its name.

The same check also set aside an auxiliary variable that duplicates
another (a score and its rescaling, age in years and in months). Without
it both reach every item's regression, and when the predictors carry a
second exact dependency (a second such pair, or an auxiliary variable
computed from complete items) mice 3.19.0 removes only one of them and
falls back to a ridge penalty, after which its imputations are unstable:
on the five-item example of the tests, an auxiliary variable computed
from the first three items and its rescaling gave 0.768 with a
between-imputation standard deviation of 0.062, against 0.800 and 0.008
without them. So the auxiliary variables are screened before anything is
imputed: one whose multiple correlation with the columns observed on
every row that observes it (the items and the auxiliary variables kept
before it) reaches .999, mice's own threshold for a pair, is a linear
function of predictors already in the model, carries no information
about the missing values, and is left out, with a message naming it; the
items are not screened, so a planned missing design is still imputed in
full. When mice still applies its ridge penalty, which it does when the
predictors of an item are (nearly) linearly dependent, the warning says
how often, and the `imputation` attribute counts it (`ridge`).

The items are imputed as plain integer codes, whatever their class: a
column read with haven, which carries value labels, and a factor whose
labels are whole numbers are both accepted (see the `data` argument),
and the imputation model is the one for the codes. SPSS user-defined
missing values, which `haven::read_sav(user_na = TRUE)` keeps as values
of a `haven_labelled_spss` column, are missing under every treatment, as
SPSS and haven's [`is.na()`](https://rdrr.io/r/base/NA.html) treat them.

**The number of imputations.** Each bootstrap replication is imputed
`imputations` times (\\M_b\\, 20 by default). The variance of a
replicate mean is the sampling variance of the coefficient plus the
imputation variance divided by \\M_b\\, so imputation noise widens the
percentile interval by about \\\sqrt{1 + \gamma / M_b}\\, where
\\\gamma\\ is the fraction of missing information: by 1.5% at \\\gamma =
.15\\ and by 4.9% at \\\gamma = .50\\ with five imputations. The
function estimates \\\gamma\\ from the replications themselves, by the
one-way analysis of variance of the replicate by imputation estimates
(the variance between replications net of the imputation variance within
them, which is the sampling variance with infinitely many imputations,
against the imputation variance), reports the widening, and warns when
it exceeds 5%, naming the value of `imputations` that would bring it
under; when no interval is reported (fewer than half of the replications
returned a value) neither the warning nor the printed line gives a
widening, and the printed fraction of missing information is marked as
that of the replications that returned one. A replication that sets some
of its completed data sets aside averages \\n_b \< M_b\\ of them, and
the analysis of variance is then the unbalanced one (Searle, Casella, &
McCulloch, 1992, Sec. 3.6), whose interval for the ratio of the
variances is approximate; the widening is \\\sqrt{1 + \gamma \\
\overline{1/n_b}}\\, the mean taken over the replications, and the value
of `imputations` named allows for the share of completed data sets set
aside and is even when any were, whether in the replications kept or in
those set aside for them, but not when only the imputation itself failed
(see *When a completed data set is improper*). Schomaker and Heumann
(2018) found that Boot MI can do well with fewer than five imputations
per replication and covers too often with one. The point estimate
averages \\M_0\\ imputations of the observed data: at least 20 and at
least `imputations`, and more when von Hippel's (2020, Eq. 1) quadratic
rule, \\M = 1 + \frac{1}{2}(\gamma / \mathrm{CV})^2\\ with \\\mathrm{CV}
= .05\\, asks for more at the upper limit of a 95% confidence interval
for \\\gamma\\, as his two stage procedure uses it (here the limit for
the ratio of variance components in the same analysis of variance). His
target is the replicability of a standard error from Rubin's rules,
which bootstrap then impute does not use; the rule is borrowed to keep
the point estimate from changing much if the data were imputed again,
with the bootstrap as the pilot. Since the rule asks for at most 201
imputations (at \\\gamma = 1\\), a call with a bootstrap makes that many
imputations of the observed data, and fits and checks each, before the
bootstrap runs, and the point estimate then uses the first \\M_0\\ of
them, the number fixed by the rule and the floor and never by their
estimates. Nothing found in a completed data set of the observed data
can then discard a finished bootstrap; the imputations the rule does not
use cost at most 181 imputations and fits, less than 1% of the 20000 of
the default bootstrap. The estimate of \\\gamma\\ is imprecise when `B`
is small, and under `impute = "polr"`, whose imputations ignore the
uncertainty of the imputation model's parameters, it is too low. Without
a bootstrap the point estimate averages `max(imputations, 20)`
imputations, which keeps Rubin's relative efficiency, \\(1 + \gamma /
M)^{-1}\\, at .976 or above for \\\gamma\\ up to .5. The `se` row is the
standard error of the point estimate, \\\sqrt{s^2\_{boot} + s^2\_{imp} /
M_0}\\, with \\s^2\_{boot}\\ the sampling variance and \\s^2\_{imp}\\
the imputation variance from the same analysis of variance, so it counts
the noise of the point estimate's own imputations. The analysis of
variance is a construction for this function, not a published procedure
for omega.

**Run time under multiple imputation.** A replication imputes its
resample `imputations` times and refits the factor model on each
completed data set, so the run time grows with `B` times `imputations`.
With 10 items, five categories, and 300 cases, a replication with five
imputations took about half a second under predictive mean matching and
about five seconds under the proportional odds model, so the default
`B = 1000` runs for about eight minutes and well over an hour,
respectively, on one core with five imputations. The default of 20
imputations per replication does about four times that imputation and
fitting work. This larger default reduces imputation noise; it has not
been separately calibrated for interval coverage. The default `B` is
1000 under `missing = "mi"` rather than the 10000 of the other two
treatments (a percentile interval reads its limits from less extreme
quantiles than BCa), and why BCa, whose jackknife would add
`imputations` imputations per case, is refused. `ncpus` spreads the
replications over processes; every resample, and a seed for each
replication's imputations, is drawn before any replication runs, so a
given `seed` returns the same result whatever `ncpus` is.

**Comparison with other packages.** The psych package's `omega` fits a
Schmid-Leiman hierarchical factor model to Pearson correlations or, with
`poly = TRUE`, to polychoric correlations, and does not implement
categorical omega in the sense of Green and Yang (2009), which
`reliability_omega_categorical` computes for ordered-categorical items;
`polychoric` provides polychoric correlation estimation as a separate
tool.

## References

Asparouhov, T., & Muthén, B. (2010). *Weighted least squares estimation
with missing data* (Mplus technical appendix). Muthén & Muthén.

Bartlett, J. W., & Hughes, R. A. (2020). Bootstrap inference for
multiple imputation under uncongeniality and misspecification.
*Statistical Methods in Medical Research, 29*(12), 3533–3546.
[doi:10.1177/0962280220932189](https://doi.org/10.1177/0962280220932189)

Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
inclusive and restrictive strategies in modern missing data procedures.
*Psychological Methods, 6*(4), 330–351.
[doi:10.1037/1082-989X.6.4.330](https://doi.org/10.1037/1082-989X.6.4.330)

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science, 11*(3), 189–228.

Drezner, Z., & Wesolowsky, G. O. (1990). On the computation of the
bivariate normal integral. *Journal of Statistical Computation and
Simulation, 35*(1–2), 101–107.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Genz, A. (2004). Numerical computation of rectangular bivariate and
trivariate normal and t probabilities. *Statistics and Computing,
14*(3), 251–260.
[doi:10.1023/B:STCO.0000035304.20635.31](https://doi.org/10.1023/B%3ASTCO.0000035304.20635.31)

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*, 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

Jia, F., & Wu, W. (2019). Evaluating methods for handling missing
ordinal data in structural equation modeling. *Behavior Research
Methods, 51*, 2337–2355.
[doi:10.3758/s13428-018-1187-4](https://doi.org/10.3758/s13428-018-1187-4)

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Lim, A. J.-M., & Cheung, M. W.-L. (2022). Evaluating FIML and multiple
imputation in joint ordinal-continuous measurements models with missing
data. *Behavior Research Methods, 54*, 1063–1077.
[doi:10.3758/s13428-021-01582-w](https://doi.org/10.3758/s13428-021-01582-w)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance in
ordered-categorical measures. *Multivariate Behavioral Research, 39*(3),
479–515.
[doi:10.1207/s15327906mbr3903_4](https://doi.org/10.1207/s15327906mbr3903_4)

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Schomaker, M., & Heumann, C. (2018). Bootstrap inference when using
multiple imputation. *Statistics in Medicine, 37*(14), 2252–2266.
[doi:10.1002/sim.7654](https://doi.org/10.1002/sim.7654)

Searle, S. R., Casella, G., & McCulloch, C. E. (1992). *Variance
components*. Wiley.

Shi, D., Lee, T., Fairchild, A. J., & Maydeu-Olivares, A. (2020).
Fitting ordinal factor analysis models with missing data: A comparison
between pairwise deletion and multiple imputation. *Educational and
Psychological Measurement, 80*(1), 41–66.
[doi:10.1177/0013164419845039](https://doi.org/10.1177/0013164419845039)

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
imputation by chained equations in R. *Journal of Statistical Software,
45*(3), 1–67.
[doi:10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03)

von Hippel, P. T. (2020). How many imputations do you need? A two-stage
calculation using a quadratic rule. *Sociological Methods & Research,
49*(3), 699–718.
[doi:10.1177/0049124117747303](https://doi.org/10.1177/0049124117747303)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
(use for continuous items),
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md)
(for dichotomous items), `omega`, `polychoric`.

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# The mixed-skew design of Green and Yang (2009, Table 1): four items
# loading 0.8 on one factor, each cut into five categories, two piled
# into the top category and two into the bottom. Their population
# categorical omega for this design is 0.716.
N <- 1000
p_top <- c(.011, .027, .047, .118, .797)
cut_top <- c(-Inf, qnorm(cumsum(p_top)[-5]), Inf)
cut_bottom <- c(-Inf, qnorm(cumsum(rev(p_top))[-5]), Inf)
latent <- outer(rnorm(N), rep(0.8, 4)) +
          matrix(rnorm(N * 4, sd = 0.6), N, 4)
items <- cbind(y1 = cut(latent[, 1], cut_top, labels = FALSE),
               y2 = cut(latent[, 2], cut_top, labels = FALSE),
               y3 = cut(latent[, 3], cut_bottom, labels = FALSE),
               y4 = cut(latent[, 4], cut_bottom, labels = FALSE))

# Categorical omega as Green and Yang defined it: the total variance
# of the summed categories comes from the correlations the fitted
# single-factor model reproduces.
reliability_omega_categorical(data = items, denominator = "model_implied")
#> Categorical omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "bca" (the recommended method) or "percentile"; B = 10000 replications is the default when you do.
#>  term       value
#>  estimate   0.714
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          1000 
#>  N_complete 1000 
#>  n_items    4    
#> 
#> Denominator: variance of the summed categories from the model implied polychoric correlations (Green & Yang, 2009)
#> Largest residual correlation under the single-factor model: 0.028 (y1 with y4)

# The default, categorical omega as Kelley and Pornprasertmanit (2016)
# computed and recommended it: the total variance comes from the
# unrestricted polychoric correlations. The single factor holds in
# this design, so the two agree to within about 0.01; they part ways
# as the model misfits.
reliability_omega_categorical(data = items)
#> Categorical omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "bca" (the recommended method) or "percentile"; B = 10000 replications is the default when you do.
#>  term       value
#>  estimate   0.720
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          1000 
#>  N_complete 1000 
#>  n_items    4    
#> 
#> Denominator: variance of the summed categories from the unrestricted polychoric correlations (Kelley & Pornprasertmanit, 2016)
#> Largest residual correlation under the single-factor model: 0.028 (y1 with y4)

# The same sum score with the items treated as continuous. Robust
# omega falls far below 0.716 here, because items whose categories
# pile up at opposite ends relate to the factor far from linearly,
# which is the case Green and Yang designed their coefficient for.
reliability_omega(data = items)
#> Robust omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "percentile" (or "bca"); B = 10000 replications is the default when you do.
#>  term       value
#>  estimate   0.549
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          1000 
#>  N_complete 1000 
#>  n_items    4    
#> 
#> Denominator: observed composite variance (robust omega)
#> Largest residual correlation under the single-factor model: 0.283 (y3 with y4)

# The same items with values missing at random: respondents low on y1
# and y2 leave y3 and y4 unanswered more often. Listwise deletion, the
# default, keeps only the complete rows, and here that leaves too
# little: the few respondents who determine the correlation of a
# top-piled item with a bottom-piled one (one to five of the 1000 for
# each pair) stand low on y1 and y2, the complete rows keep none of
# them for y2 with y3 or for y1 with y4, and the correlations left do
# not identify the loadings, so a call with the default stops and says
# so (see Details).
z <- as.numeric(scale(items[, "y1"] + items[, "y2"]))
items_mar <- items
items_mar[runif(N) < plogis(-1.5 - 1.2 * z), "y3"] <- NA
items_mar[runif(N) < plogis(-1.5 - 1.2 * z), "y4"] <- NA

# Pairwise deletion estimates each polychoric correlation from the rows
# that observe both items, so every row contributes. It still assumes
# values missing completely at random, as listwise deletion does.
# Multiple imputation, which assumes only that they are missing at
# random, takes longer than an example here can (see Details).
reliability_omega_categorical(data = items_mar,
                              denominator = "model_implied",
                              missing = "pairwise")
#> The data do not determine the polychoric correlation of y2 and y3 (the bivariate table has empty cells). The model implied denominator does not need it, and the correlations the data do determine identify the single-factor loadings; its residual correlation is reported as NA.
#> Categorical omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "bca" (the recommended method) or "percentile"; B = 10000 replications is the default when you do.
#>  term       value
#>  estimate   0.691
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          1000 
#>  N_complete 686  
#>  n_items    4    
#> 
#> Denominator: variance of the summed categories from the model implied polychoric correlations (Green & Yang, 2009)
#> Largest residual correlation under the single-factor model: 0.029 (y2 with y4)
#> 
#> Missing data: pairwise deletion, consistent only when values are missing completely at random, as listwise deletion is (Asparouhov & Muthén, 2010); it uses more of the data, not a weaker assumption, and missing = "mi" is the remedy when values are missing at random
```
