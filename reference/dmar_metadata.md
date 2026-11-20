# Machine-Readable Metadata of a DMAR Result, and the Classes of Its Refusals

Returns, in one list with fixed names, what a program built on DMAR
needs to read beside the numbers of a result table: which estimator,
standard errors, test, and missing data method produced it; which
version of a fit measure was reported and what corrected it; how many
cases and groups it rests on and, after two-stage estimation, how many
rows each stage read; which auxiliary variables entered; whether the
result is fully eligible or qualified by an unavailable unrestricted
(H1) model; and the counts of a bootstrap. Every DMAR table already
carries these as attributes, but under names that differ by family; this
accessor reads them into one stable vocabulary, so that an application
reads `dmar_metadata(x)$eligibility$status` rather than the attribute of
the day. A field a result does not have is `NULL`.

## Usage

``` r
dmar_metadata(x)
```

## Arguments

- x:

  A result table of a DMAR function (a `dmar_tbl`).

## Value

A named list with the fields above; a field the result does not carry is
`NULL`.

## Details

### Fields

- `result`:

  The table's own class, the first element of `class(x)`
  (`"dmar_sem_fit"`, `"dmar_reliability"`, `"dmar_fa_confirmatory"`, and
  so on), and `coefficient`, the reliability coefficient a reliability
  table reports.

- `estimator`, `se_method`, `test`:

  The point estimator, the method of the standard errors, and the test
  statistic actually used, as the function recorded them: lavaan's names
  for a fitted model (`"ML"`, `"robust.huber.white"`,
  `"satorra.bentler (scaled)"`), the reliability family's names for its
  own (`"analytic"`). The fit measure families record lavaan's point
  estimator (`"ML"` for a fit requested with `estimator = "MLR"`, whose
  robustness is in `se_method` and `test`); the information criteria
  record the estimator as requested (`"MLR"`).

- `missing`:

  The missing data method as the function recorded it: lavaan's names
  for a fitted model (`"listwise"`, `"ml"`, `"two.stage"`,
  `"robust.two.stage"`) and DMAR's for its own arguments (`"fiml"`,
  `"two_stage"`, `"robust_two_stage"`, `"pairwise"`, `"none"` or
  `"covariance input"` for a supplied matrix).

- `ci_method`, `wald_df`, `conf_level`, `rmsea_conf_level`,
  `tail_areas`, `denominator`:

  How an interval was formed: the method, the degrees of freedom of a
  Wald interval (`Inf` for the normal), the level, the level of the
  RMSEA interval, unequal tail areas when recorded, and the denominator
  of an omega. `rmsea_conf_level` is the level the table records for its
  RMSEA interval apart from its other intervals, 0.90 by default (the
  fit rows of
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
  [`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
  [`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
  and
  [`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
  whose parameter intervals are at `conf_level`, the fit line of
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md),
  and the fit measures of
  [`irt_grm`](https://yelleknek.github.io/DMAR/reference/irt_grm.md)),
  or else, for a table whose RMSEA interval is its own
  ([`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
  [`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md),
  [`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
  [`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md)),
  the table's level; it is `NULL` for a table with no RMSEA interval.

- `interval_degenerate`:

  For a table of an order-based effect size
  ([`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
  [`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
  [`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md))
  whose estimated variance is 0, so that its two limits equal the
  estimate and are not a confidence interval: `"complete_separation"` or
  `"no_variation"`. `NULL` for a table with an interval.

- `version`, `version_reason`, `correction`:

  For a fit measure table, the version reported (`"standard"`,
  `"scaled"`, `"robust"`), why, and the correction behind a robust
  version, from the `sem_fit` attribute (the information criteria of
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
  record under `sem_information`, read the same way).

- `N`, `groups`, `fit_rows`:

  The cases the result rests on, where the function records them (an
  attribute, or the `N` row of a reliability table), and the number of
  groups; after two-stage estimation, `fit_rows` is the record of the
  rows the first stage read and lavaan counted in each group, how many
  were observed on an auxiliary variable alone, and whether the fit was
  computed again with all of them (`read`, `nobs`, `extra`, `corrected`,
  which is `NA` when no case was observed on an auxiliary variable
  alone, and `reestimated`, whether that computation estimated a model
  that links the groups again).

- `auxiliary`:

  The auxiliary variables of a FIML or two-stage fit, where the function
  recorded them.

- `eligibility`:

  A list with `status`, `code`, and `note`: `"eligible"` when nothing
  qualifies the result (`code` and `note` then `NULL`), and
  `"qualified"` with a code and the text of the warning the function
  gave: `"h1_unavailable"` when the unrestricted (H1) likelihood was
  unavailable and the constrained model's estimates were retained with
  their H0-based uncertainty, and `"two_stage_rows_unavailable"` when a
  two-stage fit whose first stage read cases observed on an auxiliary
  variable alone could not be computed again with every case, so that
  its standard errors and the values that rest on lavaan's N are
  qualified (see *Refusals* for the cases where nothing is returned).

- `bootstrap`:

  For a bootstrap result, a list with `B` (attempted), `used`, `failed`,
  `improper`, and, where the family records them, the finer counts of
  the reliability engine (`returned`, `undetermined`, `not_converged`,
  `too_few_cases`), with `required`, the fewest replications a request
  may ask for, and `status`, one of `"complete"`, `"dropped"`,
  `"exploratory"`, `"insufficient_count"`, or `"insufficient_share"`
  (see
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)).

- `simulation`:

  For the table of a Monte Carlo sensitivity function (the
  `ss_aipe_*_sensitivity()`, `ss_power_*_sensitivity()`, and
  [`ss_seq_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_seq_c_sensitivity.md)
  functions), the record of its run: a list with `attempted` (the
  replications asked for), `returned` (those that produced a confidence
  interval, or a test), `failed` (those that did not), `denominator`
  (the replications the proportion, width, and coverage rows are taken
  over, which is `returned`), `width_target_probability`, and
  `width_target_mc_ci`, with any element the function adds (see the
  section on the record of the simulation on
  [`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md)).
  For the table of a standardized cluster randomized planner
  ([`ss_aipe_crd_es`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd_es.md))
  that simulated a width, the record of the simulation of the design
  reported, with the same three counts, a `denominator` that is
  `attempted` for an assured width, and the elements that page lists
  (`failure_reasons`, `width_mc_ci`, and the totals of the search).
  `NULL` for any other table.

- `transformation`, `methods`, `rotation`, `structure`:

  Family specifics where they exist: the transformation record of
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
  the score methods of
  [`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md)
  as a plain character vector, the rotation and the fitted structure
  (covariance or correlation) of
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md).

- `footer`:

  The lines a table of the fit families prints beneath itself, which
  state the version, the sample, and the qualifications in words. It is
  `NULL` for other tables, including one that prints a line of its own
  (the confidence set of
  [`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md)
  is in that table's `nnt_confidence_set` attribute).

The names are part of DMAR's interface: a field is added when a new
family needs one and is not renamed or removed within a release cycle
without a deprecation.

### Refusals

When DMAR will not run an analysis as asked, it stops with the reason
and what to do instead. Every error DMAR signals on purpose, from the
check of an argument to a fit that did not converge, is a classed
condition: it carries, beside its own class, the class `dmar_refusal`
and a `code` element with the name of its own class without the `dmar_`
prefix, so the refusals can be caught together with
`tryCatch(expr, dmar_refusal = function(e) ...)` or each by its own
class. An error that is not a `dmar_refusal` was not signaled by DMAR's
own checks (it came from R, from another package, or from
[`match.arg()`](https://rdrr.io/r/base/match.arg.html) for a choice that
is not among an argument's values). The classes:

- `dmar_h1_unavailable`:

  The unrestricted (H1) likelihood of a FIML or two-stage fit is
  unavailable (its affine-collapse diagnostic fired), so a fit measure,
  a two-stage criterion, or an EFA uncertainty that needs it is not
  computed; the condition carries the diagnostic as `diagnostic`. The
  same name is the class of the warning a function gives when it retains
  the constrained model's estimates and qualifies them instead;
  `inherits(cnd, "error")` tells the refusal from the qualification.

- `dmar_two_stage_rows_unavailable`:

  A two-stage fit whose first stage read cases observed on an auxiliary
  variable alone could not be computed again with every case, so its fit
  measures are not reported; the condition carries the row counts as
  `rows`.

- `dmar_inadmissible`:

  A fitted model given to
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
  [`sem_residual_fit`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
  [`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
  [`ecvi`](https://yelleknek.github.io/DMAR/reference/ecvi.md), or
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)
  has an inadmissible solution: lavaan's post-estimation check, made on
  the model's own parameters when the fit has auxiliary variables, finds
  a negative variance estimate or a covariance matrix that is not
  positive definite. The condition's `reason` element names each failure
  once, in the order the check makes them: `"negative_variance"`,
  `"cov.lv"` (the covariance matrix of the latent variables), `"theta"`
  (that of the residuals of the observed variables), or `"psi"` (that of
  the residuals of the structural part of the model), the last three the
  names
  [`lavaan::lavInspect()`](https://rdrr.io/pkg/lavaan/man/lavInspect.html)
  gives the matrices; `variables` holds the variables with a negative
  variance estimate; `problems` is a data frame with a row per failure
  and the columns `reason`, `group` (the label of the group, as lavaan
  gives it; `NA` for a fit of one group), `level` (1 or 2 in a
  multilevel fit; `NA` for a fit of one level), `variable` (`NA` for a
  matrix), and `value` (the negative estimate, or the smallest
  eigenvalue of the matrix); and `auxiliary` names the auxiliary
  variables set aside (empty without them). The group here is a label, a
  character string, where the refusals of
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)
  below carry the group's index, a whole number counting the groups in
  lavaan's order.

- `dmar_two_stage_fixed_x`:

  The residual based statistic of a normal theory two-stage fit, which
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  and
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  report, is refused for a fit that conditions on observed exogenous
  variables (`fixed.x = TRUE` or `conditional.x = TRUE`); the condition
  carries those variables as `variables`, and the message says how to
  refit (`fixed.x = FALSE`, with `conditional.x = FALSE` too for a fit
  made with `conditional.x = TRUE`) and on which rows to keep the fit's
  estimates.

- `dmar_exact_fit_bootstrap_scope`:

  A fit that
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)
  does not take. The condition's `reason` element says why:
  `"summary_data"` (a fit to summary statistics), `"ordered"` (ordered
  indicators, named in `variables`), `"estimator"` (an estimator other
  than maximum likelihood, given as `estimator`), `"two_level"`,
  `"cluster"` (a cluster variable, given as `cluster`), `"weights"`
  (sampling weights), `"two_stage"` (two-stage estimation), `"fixed_x"`
  (a fit conditional on observed exogenous variables, named in
  `variables`), `"missing_without_fiml"` (missing values left in the
  analysis data of a fit not made by FIML, as under pairwise deletion,
  with the group's index as `group`), or `"no_test"` (no standard test
  statistic).

- `dmar_pattern_too_small`:

  A missing data pattern with no more rows than observed variables,
  which the pattern-wise transformation of
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)
  cannot use. The condition carries the pattern's index among its
  group's patterns, numbered in the order of their first rows
  (`pattern`), its rows (`n`), the number and the names of the variables
  it observes (`p` and `observed`), and the group's index (`group`).

- `dmar_transformation_undefined`:

  The model based transformation of
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)
  is undefined for a group's data or moments. The condition's `reason`
  element says why: `"not_positive_definite"` (a covariance matrix,
  named in `matrix`, is not positive definite), `"no_mean_structure"`,
  `"nonpositive_variance"` (a saturated variance is not positive),
  `"singular_pattern"` (the model implied covariance matrix of the
  variables a pattern observes, the pattern given as `pattern`, is
  singular), or, for the Newton iteration of the common transformation,
  `"not_converged"`, `"singular_jacobian"`, or `"diverged"` (with the
  step as `step`); `group` is the group's index.

- `dmar_log_ratio_undefined`:

  [`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md)
  or
  [`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md)
  was asked for the log ratio interval of a ratio whose numerator or
  denominator is not positive, where the logarithm is undefined.

- `dmar_bootstrap_needs_data`:

  One of the same two functions was asked for a bootstrap interval from
  summary statistics, which cannot be resampled.

- `dmar_items_invalid`:

  The items to analyze cannot be read from the input, whether they were
  named by `items` (in
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md),
  [`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md),
  [`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md),
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
  and
  [`measurement_alignment`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md)),
  by the factors of
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  and of the functions that fit through it
  ([`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
  [`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
  [`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
  and
  [`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md)),
  or by the model of
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  (a named list of factors or lavaan syntax) or of
  [`sem_sam`](https://yelleknek.github.io/DMAR/reference/sem_sam.md).
  The condition's `reason` element says which check failed: `"form"`
  (not a character vector, or a name is `NA`), `"unknown"` (a name that
  is not a column of `data` or an item of `S`, the empty name included),
  `"duplicate"` (a name given twice), `"count"` (fewer than three
  items), or, with raw data, `"data_names"` (two columns of `data` share
  a name, or, with `items = NULL`, which makes every column an item, a
  column's name is empty or `NA`); and `items` holds the names at issue,
  empty when the failure is not about particular names. A model's
  structure is checked first and without this class (too few items on a
  factor, an item named twice among the factors), so from
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  the functions that fit through it, and the `model` of
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  the reason is `"unknown"`, or `"form"` for a name that is `NA`, and
  from
  [`sem_sam`](https://yelleknek.github.io/DMAR/reference/sem_sam.md)
  `"unknown"` or `"data_names"`; the `items` of
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  and
  [`measurement_alignment`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md)
  can give any reason but `"data_names"`. A target naming items the
  input lacks is `dmar_target_invalid`.

- `dmar_htmt_blocks_invalid`:

  The `blocks` of
  [`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) are not
  disjoint sets of distinct items under distinct names. The condition's
  `reason` element says which check failed: `"construct_names"` (a
  construct name given to two or more blocks), `"repeated"` (an item
  listed more than once in one block), `"shared"` (an item in two or
  more blocks), or `"data_names"` (an item whose name belongs to two or
  more columns of `data`); `blocks` names the blocks and `items` the
  items at issue (empty for `"construct_names"`). The refusal is made
  before any correlation is computed or any resample drawn. A block that
  is not two or more columns of `data` is `dmar_argument_invalid`.

- `dmar_target_invalid`:

  The `target` of a target rotation in
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
  cannot be read as one (or a target was given with another rotation, or
  none with a target rotation); the condition's `reason` element says
  which check failed: `"rotation"`, `"form"`, `"names"`, `"factors"`,
  `"values"`, `"none_specified"`, or `"nfactors"`.

- `dmar_target_unidentified`:

  The target does not identify the rotation: `reason` is `"counts"` (too
  few specified loadings in a column), `"duplicate"` (two columns
  specified alike), or `"rank"` (the unrotated loadings of a column's
  specified items have too low a rank), and `factors` names the columns
  at issue.

- `dmar_group_coverage`:

  Under full information maximum likelihood, a group of
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
  never observes a variable of the model (or an auxiliary variable), or
  never observes two of them together, so the ladder cannot be fitted;
  the condition carries the variables a group never observes as
  `unobserved` (a data frame with columns `group`, `variable`, and
  `n_group`, the cases of the group) and the pairs as `pairs` (columns
  `group`, `variable_1`, and `variable_2`).

- `dmar_bifactor_unidentified`:

  A bifactor structure that
  [`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
  will not fit; the condition carries the factors and items at issue as
  `factors` and `items`, and the kind of refusal as `reason`:
  `"unidentified"` (the data cannot identify the structure), `"names"`
  (an item or factor placed twice, or names that clash), or `"scope"` (a
  structure outside the model the function fits, such as fewer than two
  group factors, which need not be unidentified).

- `dmar_second_order_unidentified`:

  A second-order structure that
  [`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md)
  will not fit; the condition carries the factors and items at issue as
  `factors` and `items`, and the kind of refusal as `reason`, with the
  values of `dmar_bifactor_unidentified` (`"names"` also for a
  first-order factor under no second-order factor, `"scope"` for an item
  on two first-order factors or a first-order factor under two
  second-order factors).

- `dmar_count_invalid`:

  A number of resamples, permutations, or replications that is not a
  single whole number in its admissible range (2.5, a missing or
  infinite value, a vector, or a value below the stated minimum). The
  condition carries the argument's name as `arg` and its limits as `min`
  and `max` (`NULL` when no maximum is stated).

- `dmar_semipartial_r_incoherent`:

  A full-model \\R^2\\ smaller than the square of the semipartial
  correlation it is given with, which no data can produce, in
  [`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
  [`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
  and
  [`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md).
  The condition carries the two values as `R2_full` and `r_sp`.

- `dmar_bca_undefined`:

  [`mediate`](https://yelleknek.github.io/DMAR/reference/mediate.md) was
  asked for the BCa interval of an indirect effect in data where the
  interval is not defined: the indirect effect is not estimable with one
  of the cases left out, so the jackknife acceleration has no value, or
  the bias correction or the acceleration is not finite. The percentile
  interval does not need these quantities.

- `dmar_scores_not_ordered`:

  The observations given to an order-based effect size
  ([`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
  [`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
  [`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md))
  have no common order: an unordered factor, a character vector, an
  ordered factor beside anything but an ordered factor with the same
  levels, time differences in different units, a `Date` beside a
  `POSIXct`, or a class or shape that is not a numeric vector. The
  message says which case applies, what would have been computed, and
  what to supply.

- `dmar_studentized_statistic_undefined`:

  `statistic = "t"` in
  [`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md)
  or
  [`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md)
  for data whose observed standard error is 0 with a nonzero mean
  difference (two constant groups at different values; paired
  differences that are all the same nonzero value), where the
  studentized statistic has no value. The condition carries `reason`
  (`"constant_groups"` or `"constant_differences"`), `alternative`, the
  group sizes as `n_1` and `n_2` or the number of pairs as `n`, and
  `p_value_mean`, the exact *p*-value that `statistic = "mean"` gives.

- `dmar_hedges_correction_undefined`:

  `unbiased = TRUE` in
  [`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md) with a
  control group of 1 or 2, or in
  [`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) with
  \\n_1 + n_2\\ of 2 or 3. On 1 degree of freedom the correction factor
  of Hedges (1981) is 0, so the corrected estimate would be 0 whatever
  the data; on 0 degrees of freedom the factor is not defined, because
  there is no standard deviation. The condition carries `df`.

- `dmar_assurance_invalid`:

  An `assurance` given to
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
  that is neither `NULL` nor a single number strictly between 0 and 1.
  The condition carries the argument's name as `arg` and the value given
  as `value`.

- `dmar_argument_invalid`:

  A value, or a combination of values, that a function cannot use and
  that has no class of its own: the wrong type or length, a number
  outside its range, arguments that cannot be given together or that are
  missing together, a column that is not in the data, data or a fitted
  model the analysis does not cover. It is the class of every argument
  check in the package. The condition carries the name or names of the
  arguments at issue as `arg`. `arg` is `NULL` when the refusal is not
  about a named argument (a property of the data, a search that ran out
  of room), and in some checks whose message names the argument only in
  running text and not at its start. Among others it is given by
  [`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md)
  (a `p_restricted` that is not below `p_full`, squared multiple
  correlations that cannot belong to nested models, an effect stated in
  more than one way); by the partial and semipartial correlation
  functions for a `p` below 2; by
  [`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) and
  [`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md) for
  `unbiased = TRUE` without the sample sizes, and by
  [`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) for a group
  of one with its own standard deviation; by
  [`ss_aipe_sc_ancova_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md)
  when a sample size is to be planned and no `width` is given; and by
  [`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md),
  [`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md),
  [`ss_aipe_cliff_delta_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
  and
  [`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md)
  when both the estimate to plan from and `specified_N` are given.

- `dmar_package_required`:

  The analysis asked for needs a package that is not installed (or, for
  lavaan, a newer version than the one installed). The message names the
  package and how to install it, and the condition carries its name as
  `package`.

- `dmar_computation_failed`:

  The arguments were accepted, and a computation made with them did not
  give a result: a search for a sample size or for a confidence limit
  did not reach its target, a model fitted for the analysis failed or
  gave no proper solution and no more specific class applies
  (`dmar_not_converged`, `dmar_no_maximum`, `dmar_lavaan_error`), a
  refit or a bootstrap failed, or a distribution could not be evaluated
  accurately at the values reached. The message says which and what to
  change.

- `dmar_internal_error`:

  DMAR reached a state its own code should not reach. The arguments are
  not at fault; please report the call (see `BugReports` in the
  package's description).

- `dmar_desired_power_invalid`:

  A `desired_power` that is not a single number strictly between 0 and
  1, in any power planner (the `ss_power_` functions). The condition
  carries the argument's name as `arg` and the value given as `value`.

- `dmar_effect_zero`:

  A sample size was asked of
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)
  or
  [`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md)
  for a planning value under which the null hypothesis is true (a
  squared multiple correlation, or a change in it, of 0), so that power
  equals `alpha_level` at every sample size; or of
  [`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md)
  or
  [`ss_power_rc`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md)
  for a regression coefficient whose effect is 0 at the planning values,
  with a `desired_power` above `alpha_level`.

- `dmar_estimand_undetermined`:

  An estimand asked for by name (`estimand = "total"` or `"partial"`) of
  [`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
  [`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
  or
  [`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
  that the inputs do not determine, such as the total value from an *F*
  statistic, or a total for which no interval is formed. With
  `estimand = "all"` the same estimand is reported as `NA`. The
  condition carries `estimand` and the `reason`.

- `dmar_cost_ratio_invalid`:

  A `cost_ratio` given to
  [`ss_aipe_two_phase_cost`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md)
  that is not a single positive number below 1. The condition carries
  the argument's name as `arg` and the value given as `value`.

- `dmar_two_phase_no_design`:

  Planning values for which
  [`ss_aipe_two_phase_cost`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md)
  cannot represent the sizes of a design as whole numbers (a width too
  narrow for the standard deviation, or a proxy too inexpensive).

- `dmar_no_maximum`:

  A factor analysis or structural equation model fit whose likelihood
  has no maximum: an error variance ran off below zero, and with it held
  further out (at ten times its estimate, or one model implied variance
  below it, whichever is further) the log likelihood is higher, or lower
  by no more than \\10^{-6}\\ per case, so the likelihood is still
  rising or flat along that ridge and lavaan's optimizer stopped at a
  point its starting values and tolerance set; the estimates describe
  the optimizer rather than the data. It comes from
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  (and
  [`cfa_k()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md), and
  [`cfa_2`](https://yelleknek.github.io/DMAR/reference/cfa_2.md)),
  [`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
  [`fa_second_order`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
  (where it is also a `dmar_efa_no_solution`, with `reason`
  `"no maximum"`),
  [`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
  [`measurement_alignment`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
  the single-factor fits of
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  and
  [`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
  and the refits of
  [`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
  and of `reliability_H(fit = )`, for every `output`. It comes as well
  from the fit of the analysis model to the population moments in
  [`ss_aipe_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
  [`ss_power_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
  [`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md),
  and
  [`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md),
  and from the lavaan fit that the OpenMx fits of
  [`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
  start from (DMAR's own, or one given as `model`). The condition
  carries the variables as `variables` (with their groups as `group`),
  their error variances where the optimizer stopped and where the check
  held them as `estimate` and `fixed_at`, the two log likelihoods as
  `loglik`, the check's tolerance and the number of cases as `tolerance`
  and `N`, and the lavaan code of the fit as `syntax` and `lavaan_call`,
  which the message prints (leaving the syntax out when it is too long
  to print, for R cuts an error message at 8170 bytes). Inside a
  bootstrap, such a replication is counted among those that raised an
  error, and the run goes on; in the Monte Carlo of a planner it is set
  aside as a replication that did not converge. A function that reads
  only the chi square of a fit, which does not depend on where on such a
  ridge the optimizer stopped, does not check it:
  [`ss_aipe_rmsea_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea_sensitivity.md)
  and
  [`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md).

- `dmar_null_R2_unsupported`:

  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)
  was given a nonzero `null_R2`. The function plans the omnibus test
  only; the message gives the call to
  [`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md)
  for the test of an increment. The condition carries `null_R2`,
  `true_R2`, `true_f2`, and `p`.

The other refusals, each a `dmar_refusal` as well, are
`dmar_fiml_unbounded` and `dmar_fiml_unstable` (the unrestricted FIML
likelihood has no finite maximum for these data, or is numerically
unstable; the condition carries the certificate as `diagnostic`),
`dmar_efa_no_solution` (an exploratory fit that did not converge, is
improper, or left its standard errors unavailable, as its `reason`
element says; with the reason `"no maximum"` it is also a
`dmar_no_maximum`), `dmar_not_converged` (a fit that did not converge,
from
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
from
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
where it carries the rung as `level` and that rung's lavaan code as
`syntax` and `lavaan_call`, and from the functions that share the check
of
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
those under `dmar_inadmissible`, and from any other function that fits a
model, or is given a fitted one, and reads its estimates; the refusal of
a fit given to
[`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
or [`dmacs`](https://yelleknek.github.io/DMAR/reference/dmacs.md)
carries the function as `fn` and the number of iterations lavaan ran as
`iterations`), `dmar_lavaan_error`, `dmar_too_few_cases`,
`dmar_htmt_undefined`, `dmar_H_improper` (an improper solution given to
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
as a fit), `dmar_fmi_saturated_singular`, and `dmar_fmi_not_inverted`.

A warning is a qualification of a result that is returned. Every warning
DMAR signals on purpose is classed and carries its `code`, so it can be
muffled or read by class from
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html). Two
classes are general: `dmar_result_qualified` (the result is returned and
the message states what to know before using it: an approximation that
was used, an assumption the data do not meet, a quantity reported as
`NA`, cases or replications that were dropped) and
`dmar_argument_ignored` (an argument that was given has no effect on
this call, and the message says why). The others name one kind of
qualification: `dmar_fmi_not_at_maximum`, `dmar_heywood_warning` (an
improper solution: a negative error variance or a standardized loading
beyond one), `dmar_H_marginal_se` (an interval for coefficient H from
the marginal standard errors of the loadings alone), `dmar_nc_F_clamp`
(a lower confidence limit of a noncentrality parameter set to 0; see
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md)),
`dmar_two_stage_aux_only` (cases observed on an auxiliary variable
alone, which lavaan leaves out of its count),
`dmar_negative_reliability`, `dmar_two_stage_n`,
`dmar_semipartial_r_no_R2_full`, `dmar_semipartial_r_rounded` (a
full-model \\R^2\\ below the squared semipartial correlation by no more
than rounding to two decimals allows; the condition carries `R2_full`
and `r_sp`), `dmar_nc_F_clamp_search` (the sample size search of
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)
met the lower-limit clamp of
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md) at
intermediate sample sizes; given once per call, with the count as
`n_evaluations`), `dmar_replications_without_interval`,
`dmar_replications_with_infinite_limit` (replications of
[`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md)
whose interval has an infinite upper limit; given once per call, with
the count as `n_replications`), `dmar_uls_standard_statistic`,
`dmar_bifactor_empirical`, `dmar_smd_undefined`, `dmar_meta_hk_zero_se`
(identical effect sizes under the Hartung-Knapp adjustment; the
condition carries `n_studies`, `estimate`, `se_conventional`,
`lower_conventional`, and `upper_conventional`),
`dmar_interval_degenerate` (an order-based effect size whose estimated
variance is 0, with `reason` `"complete_separation"` or
`"no_variation"`, the `estimate`, the sample sizes `n`, and
`limit_feng_cliff`, the limit Feng and Cliff (2004) give for complete
separation, which the function does not report),
`dmar_bootstrap_replications_dropped` (a bootstrap or Monte Carlo
interval dropped replications in which the statistic was undefined; the
condition carries `B`, `B_used`, `n_dropped`, `B_required`,
`share_dropped`, `caution`, and `status`), `dmar_bootstrap_insufficient`
(more than half of the replications were dropped, or the adjusted
positions of a BCa interval take more replications than were requested,
and the limits are `NA`; the same fields, and for a BCa interval
withheld at its adjusted positions also `positions` and `B_suggest`),
`dmar_bootstrap_exploratory` (the result was computed from fewer
replications than a reported result needs because
`options(dmar.bootstrap_exploratory = TRUE)` is set; the same fields),
and `dmar_bootstrap_note` (a BCa interval replaced by the percentile
interval, or replicates that all have the same value, with every
replication used), all four described in
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).
Each function's help says when.

A refusal is never a table: when a function returns a table, its result
is eligible or qualified, never refused.

The errors `dmar_bootstrap_cancelled` and `dmar_progress_error` are not
refusals. A bootstrap of the reliability family or of
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
stops with the first when its `progress` callback returns `FALSE`, at
the request of the application that called it, and with the second when
the callback itself raises an error. Both carry `done` and `total`, the
replications attempted and the number requested, and the second keeps
the callback's own condition, with its class, as `parent`, so that an
application can tell a failing callback from a failing analysis and
still recognize a condition of its own behind it.

## References

Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
improved confidence interval. *Journal of Modern Applied Statistical
Methods, 3*(2), 322–332.
[doi:10.22237/jmasm/1099267560](https://doi.org/10.22237/jmasm/1099267560)

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.
[doi:10.3102/10769986006002107](https://doi.org/10.3102/10769986006002107)

## See also

[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md) for
the table class and its display,
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)
for `tidy()` and `glance()`, which put the same metadata into columns of
a data frame.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
x <- reliability_alpha(holzinger_swineford[c("t1_visual_perception",
                                             "t2_cubes", "t4_lozenges")])
m <- dmar_metadata(x)
m$coefficient
#> [1] "alpha"
m$ci_method
#> [1] "bonett"
m$eligibility$status
#> [1] "eligible"

# Four items, so that the model has degrees of freedom to test (three
# would make it saturated, which sem_fit_measures() qualifies).
fit <- lavaan::cfa(
  "visual =~ t1_visual_perception + t2_cubes + t3_paper_form_board +
             t4_lozenges",
  data = holzinger_swineford, std.lv = TRUE, estimator = "MLR")
f <- dmar_metadata(sem_fit_measures(fit))
f$version
#> [1] "robust"
f$correction
#> [1] "nonnormality corrections of Brosseau-Liard, Savalei, and Li (2012) and Brosseau-Liard and Savalei (2014)"
```
