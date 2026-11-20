#' Tidy and Glance Methods for DMAR Result Tables
#'
#' @description
#' A DMAR function returns a tidy \code{data.frame} built to be read: one
#' row per quantity, a numeric \code{value} column, and a display layer
#' that rounds the printed values, not the stored ones (see
#' \code{\link{dmar_tbl}}). The tidy verbs \code{\link[generics]{tidy}}
#' and \code{\link[generics]{glance}} give the same numbers in the two
#' shapes a programmer usually wants instead: one row per term with a
#' typed column for each quantity, and a one-row summary of the result as
#' a whole. This page states that contract once, for the confidence
#' interval family, the post hoc family, the contrast tests, and the
#' power-based sample size planners.
#'
#' @details
#' \strong{What the verbs return.}
#' \code{tidy(x)} returns a \code{data.frame} with one row per term, where
#' a term is whatever the family produces one of: a parameter estimate, a
#' contrast, a planned design. Its columns follow the naming convention
#' the broom ecosystem uses, which separates words with dots rather than
#' the underscores DMAR uses everywhere else: \code{term},
#' \code{estimate}, \code{se}, \code{statistic}, \code{p_value},
#' \code{ci_lower}, \code{ci_upper}, and \code{conf_level}. A method
#' reports the subset of those columns its family can fill, plus any
#' column the family genuinely adds, such as \code{p_adjusted} for a
#' multiplicity-adjusted set of comparisons or \code{power} for a sample
#' size planner.
#'
#' \code{glance(x)} returns a one-row \code{data.frame} summarizing the
#' result as a whole, in the same dotted convention: how many comparisons
#' were made, at what confidence level, with which planning inputs. When
#' a result has a single estimand and nothing further to say at the model
#' level, as for a lone effect size and its confidence interval,
#' \code{glance()} coincides with \code{tidy()}. That is expected rather
#' than a defect, since there is no model-level quantity that the single
#' row does not already carry.
#'
#' Neither verb rounds. The \code{dmar_tbl} layer formats what is
#' printed, while \code{tidy()} and \code{glance()} return full
#' precision, which is what makes them the right input to a downstream
#' calculation or plot.
#'
#' \strong{Why \pkg{broom} is not a dependency.}
#' The \code{tidy()} and \code{glance()} generics live in \pkg{generics},
#' a small package that holds the generics and little else. \pkg{broom}
#' imports them from there, and so does DMAR, which registers its methods
#' against \code{generics::tidy} and \code{generics::glance} rather than
#' against \pkg{broom} itself. A user with \pkg{broom} or the tidymodels
#' stack loaded gets DMAR methods on the generic they already call; a user
#' with neither installed can still call \code{generics::tidy()} directly.
#' DMAR never loads \pkg{broom}, and does not need it installed.
#'
#' \strong{The families and the classes they carry.}
#' Each family tags its return with a leading S3 class, ahead of
#' \code{dmar_tbl} and \code{data.frame}, so the verbs dispatch while
#' printing and data-frame behavior are untouched.
#'
#' \tabular{lll}{
#'   \strong{S3 class} \tab \strong{Family} \tab \strong{One \code{tidy()} row is} \cr
#'   \code{dmar_ci_long} \tab confidence intervals, long form \tab an estimate and its limits \cr
#'   \code{dmar_ci_anova} \tab ANOVA effect size intervals \tab an effect and estimand, its estimate and limits \cr
#'   \code{dmar_anova_es} \tab ANOVA effect size estimates \tab an effect and estimand, its estimate \cr
#'   \code{dmar_post_hoc_ci} \tab simultaneous intervals \tab one pairwise or one contrast comparison \cr
#'   \code{dmar_contrast_test} \tab contrast tests \tab one contrast, with its test and its interval \cr
#'   \code{dmar_ss_power} \tab sample size planners \tab a planned size and the power it buys \cr
#'   \code{dmar_ss_power_sensitivity} \tab planner sensitivity studies \tab a planned size and two powers
#' }
#'
#' \strong{The confidence interval family.}
#' Two classes cover the two output shapes.
#'
#' \describe{
#'   \item{\code{dmar_ci_long}}{Long-format interval tables, with rows for
#'     \code{ci_lower} and \code{ci_upper} and, when the function
#'     reports one, an estimate row whose \code{term} is the name of the
#'     parameter. Carried by \code{\link{ci_r}},
#'     \code{\link{ci_smd_c}}, \code{\link{ci_pvaf}},
#'     \code{\link{ci_reg_coef}}, \code{\link{ci_ratio_of_means}}, and
#'     \code{\link{ci_ratio_of_means_paired}}.}
#'   \item{\code{dmar_ci_anova}}{Wide-format ANOVA effect size interval
#'     tables, with one row per effect and columns for the effect name,
#'     the point estimate of each estimand with its limits
#'     (\code{eta_squared_partial}, \code{eta_squared_partial_ci_lower},
#'     \code{eta_squared_partial_ci_upper}, and the same for the total),
#'     and the design metadata. Carried by \code{\link{ci_eta_squared}}
#'     and \code{\link{ci_omega_squared}}, and by
#'     \code{\link{ci_eta_squared_generalized}}, whose one estimate sits
#'     beside \code{ci_lower} and \code{ci_upper}.}
#' }
#'
#' Both produce a \code{data.frame} with \code{term}, \code{estimate},
#' \code{ci_lower}, \code{ci_upper}, and \code{conf_level} (\code{NA}
#' when the object records no level). For \code{dmar_ci_long} it is one
#' row, and \code{glance()} calls \code{tidy()}, since the row is already
#' the whole result. For \code{dmar_ci_anova} it is one row per effect
#' and estimand, the rows of one estimand together: \code{term} names the
#' estimand (\code{"eta_squared_total"}, \code{"eta_squared_partial"}),
#' an \code{effect} column names the effect, a \code{stratum} column
#' follows it for a fit with an \code{Error()} term, and the limits of the
#' estimand are \code{ci_lower} and \code{ci_upper} (\code{NA} where no
#' interval is formed, as for the total of a factorial effect). Its
#' \code{glance()} coincides with \code{tidy()} when that is one row,
#' and otherwise is one row giving the estimands (\code{term}, separated
#' by commas when there are two), the number of effects
#' (\code{n_effects}), \code{N}, and \code{conf_level}.
#' The point estimates of \code{\link{eta_squared}} and
#' \code{\link{omega_squared}} carry \code{dmar_anova_es}, whose
#' \code{tidy()} has the same rows with \code{term}, \code{effect},
#' \code{stratum} (for an \code{aovlist} fit), and \code{estimate}, and
#' whose \code{glance()} gives \code{term}, \code{n_effects}, and
#' \code{N} (\code{NA} for \code{eta_squared()}, which reports none).
#' The two ratio of means functions carry \code{dmar_ci_ratio_of_means}
#' ahead of \code{dmar_ci_long}, because a Fieller confidence set need not
#' be an interval: their one row adds, after the limits,
#' \code{excluded_lower} and \code{excluded_upper} (the ends of the gap of
#' a set that is the complement of an interval, \code{NA} otherwise) and
#' \code{confidence_set}, the form of the set, so that a complement is not
#' read from its limits of \code{-Inf} and \code{Inf} as the whole line
#' (\code{NA} when a bootstrap interval was withheld and the limits are
#' \code{NA}).
#' \code{\link{vargha_delaney_A}} carries \code{dmar_vargha_delaney_A},
#' whose one \code{tidy()} row, named \code{"A"}, adds the test of
#' stochastic equality (\code{se}, \code{statistic}, \code{p_value}) and
#' the two sample sizes; its \code{glance()} coincides with it. Its
#' verbs have \code{df} and \code{test} after \code{p_value} and a
#' logical column \code{interval_degenerate}.
#' \code{\link{cliff_delta}} and
#' \code{\link{probability_of_superiority_paired}} carry
#' \code{dmar_cliff_delta} and
#' \code{dmar_probability_of_superiority_paired}, whose verbs are the
#' default ones with the same logical column, \code{TRUE} when the
#' estimated variance is 0 and the limits are not a confidence interval.
#'
#' A table with no method of its own falls to the default methods for
#' \code{dmar_tbl}. Its \code{tidy()} keeps the confidence level the
#' table records in a \code{conf_level} column: beside the interval of a
#' wide table (\code{NA} when no interval was formed), and, for a long
#' table of \code{term} and \code{value} rows with \code{ci_lower} and
#' \code{ci_upper} rows, on the rows of the interval, its estimate and
#' its two limits, with \code{NA} on every other row (an interval held in
#' rows named for its estimate, as \code{\link{sem_fit_measures}} holds
#' the RMSEA's in \code{rmsea_ci_lower} and \code{rmsea_ci_upper}, is
#' marked the same way, at its own level, except the interval for the
#' standardized mean difference that \code{\link{cles}},
#' \code{\link{nnt_from_smd}}, and \code{\link{proportion_of_superiority}}
#' echo in \code{smd_ci_lower} and \code{smd_ci_upper}, which their own
#' intervals are computed from and whose rows keep \code{NA}). The rows carry the
#' level and say which estimate the limits belong to, so that
#' \code{\link{results_sentence}} on a \code{tidy()} result, or on rows
#' of several \code{tidy()} results stacked with \code{rbind()}, states
#' the level and the estimate of its own interval, as it does on the table
#' itself. In a wide table that carries an interval, the
#' \code{estimate} is the column the interval belongs to (the table's own
#' \code{estimate} column, or its first numeric column that is not a
#' limit, standard error, test statistic, degrees of freedom,
#' \emph{p}-value, confidence level, or count), and the other numeric
#' columns follow under their own names. When that column has a name of
#' its own (\code{kappa} of \code{\link{fleiss_kappa}}, \code{D2} of
#' \code{\link{ci_mahalanobis}}, \code{eta_squared_partial} of
#' \code{\link{simple_effects_AB}}), an \code{estimand} column after
#' \code{estimate} keeps the name, since \code{term} names the row (from
#' the table's label columns, or the row number when it has none) and not
#' the quantity; \code{\link{results_sentence}} labels the sentence of
#' the row with it, as it labels the sentence of the table's own row. The
#' \code{tidy()} of \code{\link{content_validity_index}}, whose estimate is
#' the item's \code{i_cvi}, carries the same column, and so does that of
#' \code{\link{fraction_missing_information}}, whose estimate is the
#' parameter's FMI (\code{"fmi"}) while \code{term} names the parameter.
#' The \code{tidy()} of \code{\link{factor_count}}, whose RMSEA intervals
#' sit in its \code{rmsea_ci_lower} and \code{rmsea_ci_upper} columns beside
#' the other fit measures of each count, keeps the count as
#' \code{estimate} and adds a \code{conf_level} column holding the level of
#' those intervals, 0.90 unless the call set another.
#' The default \code{glance()} of a long table is one row with a column
#' for each term; for a table with an interval it adds an \code{estimand}
#' column, the name of the column that holds the estimate the interval
#' belongs to (the corrected correlation of
#' \code{\link{correction_for_attenuation}}, listed after the observed
#' one, and \code{NA} for a table that reports its limits alone, as
#' \code{\link{ci_rc}} and \code{\link{ci_snr}} do), and a
#' \code{conf_level} column, and names the column of the
#' \emph{p}-value the table's sentence reports, so that
#' \code{\link{results_sentence}} on the row gives the table's own
#' sentence.
#'
#' \strong{Unequal tail areas.} An interval formed from
#' \code{alpha_lower} and \code{alpha_upper} that differ records its
#' level as the coverage \eqn{1 - \alpha_{lower} - \alpha_{upper}} and
#' the two tail areas beside it. Every verb above that returns a
#' \code{conf_level} for such an interval then adds \code{alpha_lower}
#' and \code{alpha_upper} columns after it, holding the two areas on the
#' rows of the interval (\code{NA} on the other rows of a long
#' \code{tidy()}), so that \code{\link{results_sentence}} on the rows
#' states the split, or the one-sided bound, as it does on the table. An
#' equal-tailed interval, the usual one, adds neither column, so its
#' columns are those it has without the record of the split; a stack of
#' an equal-tailed and an unequal-tailed result therefore has different
#' columns, which \code{rbind()} refuses and \code{dplyr::bind_rows()}
#' fills with \code{NA}, read as equal tails.
#'
#' \strong{The post hoc family.}
#' \code{\link{ci_tukey_kramer}}, \code{\link{ci_games_howell}},
#' \code{\link{ci_dunnett_t3}}, \code{\link{ci_scheffe}}, and
#' \code{\link{ci_dunnett}} all carry \code{dmar_post_hoc_ci}. Their
#' source table is wide, with one row per comparison: a \code{contrast}
#' label, a point estimate (\code{mean_difference} for the pairwise and
#' many-to-one procedures, \code{estimate} for Scheffe), a standard
#' error, a test statistic, the \code{ci_lower} and \code{ci_upper}
#' of the simultaneous interval, and the multiplicity-adjusted
#' \code{p_adjusted}.
#' \code{tidy()} maps that to \code{term}, \code{estimate},
#' \code{ci_lower}, \code{ci_upper}, \code{p_adjusted}, and
#' \code{conf_level}, one row per comparison. \code{glance()} describes
#' the family of comparisons as a whole: how many there were, and the
#' simultaneous confidence level they hold jointly.
#'
#' \strong{The contrast tests.}
#' \code{\link{contrast_test}} carries \code{dmar_contrast_test}. Its
#' source table is wide, with one row per contrast: a \code{contrast}
#' label, the estimate \eqn{\hat{\psi} = \sum_i c_i \bar{Y}_i}, its
#' standard error, the \emph{t}-statistic and the degrees of freedom it
#' is referred to, the unadjusted \emph{p}-value, the
#' multiplicity-adjusted \code{p_adjusted}, and the \code{ci_lower} and
#' \code{ci_upper} limits. \code{tidy()} renames those to \code{term},
#' \code{estimate}, \code{ci_lower}, \code{ci_upper}, \code{statistic},
#' \code{df}, \code{p_value}, \code{p_adjusted}, and \code{conf_level},
#' one row per contrast. Both \emph{p}-values are kept, because the pair
#' is what a contrast table is read for: what the contrast would show on
#' its own, and what it shows once the family it belongs to is accounted
#' for.
#'
#' Where a post hoc procedure fixes its adjustment as part of the method,
#' a contrast test chooses one, and the same weights tested under
#' \code{adjust = "none"} and under \code{adjust = "tukey"} are two
#' different inferences. \code{glance()} therefore records the choice
#' alongside the family-level numbers: \code{n_contrasts},
#' \code{adjust}, \code{var_equal}, the smallest adjusted \emph{p}-value
#' \code{p_adjusted_min}, and \code{conf_level}. \code{adjust} and
#' \code{var_equal} name a procedure rather than measure a quantity, so
#' this one-row summary, unlike a DMAR result table, is not numeric
#' throughout.
#' The post hoc and contrast methods, and wide result tables that carry
#' adjustment metadata, retain \code{ci_adjustment} and
#' \code{p_adjustment} in both \code{tidy()} and \code{glance()}. These
#' are distinct because an interval and a test can use different
#' multiplicity procedures. \code{results_sentence()} reports an
#' adjusted \emph{p}-value beside a simultaneous interval only when the
#' recorded procedures match.
#'
#' \strong{The power-based sample size planners.}
#' A planner in the \code{ss_power_*} family returns a long table with a
#' row for the recommended sample size, a row for the realized power, and
#' rows echoing the planning inputs. A planner that reports one size and
#' one power for one design tags its return \code{dmar_ss_power}. This
#' covers the closed-form effect size planners (\code{\link{ss_power_R2}},
#' \code{\link{ss_power_r}}, \code{\link{ss_power_reg_coef}},
#' \code{\link{ss_power_smd}}, \code{\link{ss_power_sem}}), the contrast
#' and ANCOVA planners (\code{\link{ss_power_c}},
#' \code{\link{ss_power_c_ancova}}, \code{\link{ss_power_sc}},
#' \code{\link{ss_power_contrast}},
#' \code{\link{ss_power_equivalence_c}}), the ANOVA and cluster designs
#' (\code{\link{ss_power_one_way_anova}},
#' \code{\link{ss_power_factorial_anova}},
#' \code{\link{ss_power_factorial_ancova}},
#' \code{\link{ss_power_split_plot_anova}},
#' \code{\link{ss_power_rm_anova}},
#' \code{\link{ss_power_mixed_effects}}), and the mediation planner
#' \code{\link{ss_power_indirect_effect}}, whose reported power is the
#' joint power to detect the indirect effect (its row \code{power_ab}) and
#' whose component path powers \code{glance()} carries as extra columns.
#'
#' The size \code{tidy()} reports is the design's planning unit: per
#' group, per cell, per subject, or per cluster. The one-way ANOVA
#' planner, whose natural unit is the total, is summarized by its total
#' \eqn{N}. A design that reports two group sizes reports one of them
#' beside the realized power, falling through to the total \eqn{N} when
#' the per-group sizes are unequal, and \code{glance()} keeps every group
#' size as a column so none is lost. A planner whose result spans several
#' effects, with no single size-and-power summary to give, returns a plain
#' \code{dmar_tbl} and does not gain these verbs at all.
#'
#' The Monte Carlo sensitivity siblings
#' \code{\link{ss_power_R2_sensitivity}} and
#' \code{\link{ss_power_reg_coef_sensitivity}} report two powers at one
#' planned sample size, the empirical (simulated) power and the analytic
#' power, and comparing the two is the object of the study. They carry
#' \code{dmar_ss_power_sensitivity} instead: \code{tidy()} places both
#' powers beside the planned \code{sample_size}, and \code{glance()} adds
#' the simulated distribution of the estimator.
#'
#' \strong{Adding a planner to the family.}
#' A planner opts in by setting \code{dmar_ss_power} as a leading class
#' before routing its return through \code{.as_dmar_tbl()}. The rows the
#' verbs read are named in the internal vectors
#' \code{.SS_POWER_SIZE_TERMS} and \code{.SS_POWER_POWER_TERMS}. A
#' planner whose size or power row is not named there reports \code{NA}
#' rather than failing, so a new row name has to be added to those
#' vectors when a planner introduces one. \code{glance()} leaves out only
#' the two rows \code{tidy()} reads, the size and the power, and keeps
#' every other row as a column, other powers included: the marginal powers
#' of the composite planners (\code{power_a}, \code{power_b}, and
#' \code{power_ab} for parameters labeled \code{a}, \code{b}, and
#' \code{ab} in \code{\link{ss_power_composite_sem}}) and the path powers
#' of \code{\link{ss_power_indirect_effect}}.
#'
#' @param x A DMAR result object carrying one of the classes listed
#'   above.
#' @param \dots Unused, present for consistency with the generics.
#'
#' @return \code{tidy()} returns a \code{data.frame} with one row per
#'   term and broom-convention column names. \code{glance()} returns a
#'   one-row \code{data.frame} summarizing the result as a whole. Both
#'   return values at full precision.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{dmar_tbl}} for the printing layer these tables
#'   share, and the "Reading DMAR result tables" vignette for the wider
#'   output convention.
#'
#' @examples
#' # A single interval: tidy() and glance() coincide, because there is
#' # nothing at the model level the one row does not already carry.
#' res <- ci_r(r = 0.5, N = 50)
#' generics::tidy(res)
#' generics::glance(res)
#'
#' # A family of simultaneous intervals: one tidy() row per comparison,
#' # one glance() row describing the family.
#' set.seed(113)
#' y <- c(rnorm(10, 0), rnorm(10, 1), rnorm(10, 2))
#' g <- factor(rep(c("a", "b", "c"), each = 10))
#' gh <- ci_games_howell(y, group = g)
#' generics::tidy(gh)
#' generics::glance(gh)
#'
#' # A set of contrasts: tidy() keeps both the unadjusted and the
#' # adjusted p-value, and glance() names the adjustment that produced
#' # the second of them.
#' fit <- aov(bdi_post ~ condition, data = depression_bdi)
#' ct <- contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
#' generics::tidy(ct)
#' generics::glance(ct)
#'
#' # A sample size planner: tidy() gives the size and the power it buys,
#' # glance() adds the planning inputs that produced them.
#' plan <- ss_power_smd(true_smd = 0.5, desired_power = 0.80)
#' generics::tidy(plan)
#' generics::glance(plan)
#'
#' @name dmar_tidiers
NULL


# Local %||% helper (also defined in reliability.R / ci_smd.R; the
# repetition keeps each file self-contained and is cheap).
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------------
# The confidence interval family: dmar_ci_long and dmar_ci_anova
# ---------------------------------------------------------------------------

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ci_long <- function(x, ...) {
  # Standard rows that are not the parameter estimate row.
  bound_terms <- c("ci_lower", "ci_upper",
                   "actual_coverage", "achieved_coverage")
  estimate_row <- x[!x$term %in% bound_terms, , drop = FALSE]
  estimate_name <- if (nrow(estimate_row) >= 1L) {
    estimate_row$term[1L]
  } else {
    attr(x, "parameter") %||% "estimate"
  }
  estimate_val <- if (nrow(estimate_row) >= 1L) {
    estimate_row$value[1L]
  } else {
    NA_real_
  }
  .add_tail_columns(data.frame(
    term       = estimate_name,
    estimate   = estimate_val,
    ci_lower   = x$value[x$term == "ci_lower"],
    ci_upper  = x$value[x$term == "ci_upper"],
    conf_level = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  ), x)
}


#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ci_long <- function(x, ...) {
  tidy.dmar_ci_long(x, ...)
}

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ci_anova <- function(x, ...) {
  # The wide ANOVA effect size tables hold, for each effect, one estimate
  # and its limits per estimand: <estimate>, <estimate>_ci_lower, and
  # <estimate>_ci_upper (eta_squared_partial and its limits, from
  # ci_eta_squared() and ci_omega_squared() since 1.1.0), or, in
  # ci_eta_squared_generalized(), one estimate beside ci_lower and
  # ci_upper. tidy() gives one row per effect and estimand: term names the
  # estimand, the effect (and, for a fit with an Error() term, its stratum)
  # has a column of its own, and the limits are ci_lower and ci_upper. The
  # rows of one estimand come together, in the table's order. Through
  # 1.0.1 only the first row was returned, so the tidy() of a factorial
  # fit dropped every effect after the first.
  n <- nrow(x)
  level <- if (!is.null(x$conf_level)) x$conf_level
           else rep(attr(x, "conf_level") %||% NA_real_, n)
  est_cols <- .ci_anova_named_estimates(x)
  pieces <- if (length(est_cols)) {
    lapply(est_cols, function(e)
      list(term = e, estimate = x[[e]],
           lower = x[[paste0(e, "_ci_lower")]],
           upper = x[[paste0(e, "_ci_upper")]]))
  } else {
    e <- .ci_anova_estimate_col(x)
    list(list(term = e,
              estimate = if (e %in% names(x)) x[[e]] else rep(NA_real_, n),
              lower = x$ci_lower, upper = x$ci_upper))
  }
  out <- do.call(rbind, lapply(pieces, function(p) {
    o <- data.frame(term = rep(p$term, n), stringsAsFactors = FALSE)
    o$effect <- if (!is.null(x$effect)) as.character(x$effect)
                else rep(NA_character_, n)
    if (!is.null(x$stratum)) o$stratum <- as.character(x$stratum)
    o$estimate   <- as.numeric(p$estimate)
    o$ci_lower   <- as.numeric(p$lower)
    o$ci_upper   <- as.numeric(p$upper)
    o$conf_level <- as.numeric(level)
    o
  }))
  rownames(out) <- NULL
  .add_tail_columns(out, x)
}


#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ci_anova <- function(x, ...) {
  # One row of tidy() is the whole result, and glance() coincides with
  # tidy(), as for any single estimand. Otherwise the table is summarized
  # by the estimands it reports (term, separated by commas when there are
  # several), how many effects there are, the sample size, and the level;
  # through 1.0.1 its glance() was the first effect's row.
  td <- tidy.dmar_ci_anova(x, ...)
  if (nrow(td) == 1L) return(td)
  N <- if (!is.null(x$N)) unique(x$N) else NA_real_
  .add_tail_columns(data.frame(
    term       = paste(unique(td$term), collapse = ", "),
    n_effects  = nrow(x),
    N          = if (length(N) == 1L) N else NA_real_,
    conf_level = if (length(unique(td$conf_level)) == 1L)
                   td$conf_level[1L] else NA_real_,
    stringsAsFactors = FALSE
  ), x)
}

# The column of an ANOVA effect size table that holds the estimate beside
# ci_lower and ci_upper (ci_eta_squared_generalized()): the first that is
# not the effect, its stratum, a limit, or design metadata.
.ci_anova_estimate_col <- function(x) {
  meta_cols <- c("effect", "stratum", "ci_lower", "ci_upper",
                 "F_value", "df_effect", "df_error", "df_1", "df_2",
                 "N", "conf_level", "method")
  est_col <- setdiff(colnames(x), meta_cols)[1L]
  if (is.null(est_col) || is.na(est_col)) "estimate" else est_col
}

# The columns of an ANOVA effect size table that hold an estimate with its
# limits in <estimate>_ci_lower and <estimate>_ci_upper, in the table's
# order (eta_squared_total, eta_squared_partial); character(0) for a table
# whose limits are ci_lower and ci_upper.
.ci_anova_named_estimates <- function(x) {
  nms <- names(x)
  stems <- sub("_ci_lower$", "", nms[grepl("_ci_lower$", nms)])
  stems[stems %in% nms & paste0(stems, "_ci_upper") %in% nms]
}


# ---------------------------------------------------------------------------
# The ANOVA effect size estimates: dmar_anova_es
# ---------------------------------------------------------------------------

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_anova_es <- function(x, ...) {
  # One row per effect and estimand, as for the intervals: term names the
  # estimand (eta_squared_total, omega_squared_partial), and the effect
  # (and its stratum) a column of its own.
  est_cols <- .anova_es_estimate_cols(x)
  n <- nrow(x)
  out <- do.call(rbind, lapply(est_cols, function(e) {
    o <- data.frame(term = rep(e, n), stringsAsFactors = FALSE)
    o$effect <- if (!is.null(x$effect)) as.character(x$effect)
                else rep(NA_character_, n)
    if (!is.null(x$stratum)) o$stratum <- as.character(x$stratum)
    o$estimate <- as.numeric(x[[e]])
    o
  }))
  rownames(out) <- NULL
  out
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_anova_es <- function(x, ...) {
  td <- tidy.dmar_anova_es(x, ...)
  if (nrow(td) == 1L) return(td)
  N <- if (!is.null(x$N)) unique(x$N) else NA_real_
  data.frame(term = paste(unique(td$term), collapse = ", "),
             n_effects = nrow(x),
             N = if (length(N) == 1L) N else NA_real_,
             stringsAsFactors = FALSE)
}

# The estimate columns of an eta_squared() or omega_squared() table.
.anova_es_estimate_cols <- function(x) {
  grep("^(eta|omega)_squared_(total|partial)$", names(x), value = TRUE)
}


# ---------------------------------------------------------------------------
# The post hoc family: dmar_post_hoc_ci
# ---------------------------------------------------------------------------

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_post_hoc_ci <- function(x, ...) {
  # The point estimate lives in `mean_difference` (Tukey-Kramer, Games-Howell,
  # Dunnett's T3, Dunnett) or `estimate` (Scheffe); pick whichever this
  # object carries.
  est_col <- intersect(c("mean_difference", "estimate"), names(x))[1L]
  out <- data.frame(
    term        = x$contrast,
    estimate    = x[[est_col]],
    ci_lower    = x$ci_lower,
    ci_upper   = x$ci_upper,
    p_adjusted  = x$p_adjusted,
    conf_level  = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
  .dmar_tidy_adjustments(out, x)
}


#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_post_hoc_ci <- function(x, ...) {
  out <- data.frame(
    n_contrasts = nrow(x),
    conf_level  = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
  .dmar_tidy_adjustments(out, x)
}


# ---------------------------------------------------------------------------
# The power-based sample size planners: dmar_ss_power and
# dmar_ss_power_sensitivity
# ---------------------------------------------------------------------------

# The row names the tidy verbs recognize, defined once.
#
# These are the single source of truth for all three uses: tidy() reads the
# sample size and the power from them, and glance() excludes exactly these rows
# from the planning inputs it appends. Keeping one definition is what stops the
# lists from drifting apart, which is how the sample-size list came to hold
# necessary_n_per_group but not specified_n_per_group.
#
# Order is priority order, because intersect() preserves the order of its first
# argument. Per-group rows come before their total-N counterparts on both the
# planning and the specified side: a table carrying both specified_n_per_group
# and specified_N must resolve to the per-group size, the same way one carrying
# necessary_n_per_group and necessary_N does. Reporting per-group when a size
# was planned and total when one was specified would make tidy() mean different
# things for the same design.
.SS_POWER_SIZE_TERMS <- c(# design-specific per-unit sizes: a per-cell or
                          # per-cluster count is the planning unit for a
                          # factorial or cluster design, and is reported
                          # ahead of the design's total. (The repeated
                          # measures planner reports its subjects as the
                          # total, necessary_N or specified_N, below.)
                          "necessary_n_per_cell", "approximate_n_per_cell",
                          "specified_n_per_cell",
                          "n_per_cell",
                          "necessary_n_clusters_per_arm",
                          "specified_n_clusters_per_arm",
                          "necessary_n_clusters",
                          # the cluster randomized planners that solve
                          # for the size of a cluster; in the planners
                          # that solve for both sizes the number of
                          # clusters, above, is the size reported.
                          "necessary_cluster_size",
                          # a crossed design of units and raters
                          # (ss_aipe_crossed) is summarized by its units,
                          # the objects of measurement; the raters, planned
                          # or fixed by the user, follow so that a table
                          # carrying only them still tidies to a size.
                          "necessary_n_units", "necessary_n_raters",
                          "specified_n_raters",
                          # the two-phase planners plan the labeled subsample,
                          # the units that carry the costly target measurement;
                          # the pool they are drawn from is a design input (or,
                          # in the cost planner, a second answer kept beside it)
                          "necessary_n_labeled",
                          # per-group and total sizes shared across the family
                          "necessary_n_per_group", "necessary_N",
                          # the composite ANCOVA planners report a size resolved
                          # under unequal residual variances as approximate_*
                          # rather than necessary_*, because it is the smallest
                          # size at which an approximation reaches the target
                          # rather than one known to attain it (see the section
                          # on unequal residual variances in
                          # ?ss_power_composite_ancova_2group)
                          "approximate_n_per_group", "approximate_total_N",
                          "specified_n_per_group", "specified_N",
                          "specified_n_1", "specified_n_2",
                          # The legacy names sample_size, sample_size_per_group,
                          # N, and n were removed from this registry when the
                          # producers were unified: every planner now emits the
                          # unit-named vocabulary above, so recognizing the old
                          # names would only let a regression slip through.
                          # total sizes: the fallback when a single per-unit
                          # size is undefined, as in an unbalanced design whose
                          # per-group size is NA.
                          "total_N")

# composite_power is the power of a design planned against several effects that
# must all be significant (see ss_power_composite_ancova), which is that
# design's power in the sense tidy() reports. Its approximate_ sibling is the
# same quantity in a design whose groups or cells have unequal residual
# variances, where the composite integral is an approximation and the row name
# says so; both belong here, or a relabeled table would tidy to a missing power
# without complaining.
.SS_POWER_POWER_TERMS <- c("actual_power", "achieved_power",
                           "realized_power", "composite_power",
                           "approximate_composite_power", "power",
                           # ss_power_indirect_effect(): the power to detect
                           # the indirect effect ab; its power_a and power_b
                           # are component powers, not the design's, and stay
                           # glance() columns. It is last, so that a composite
                           # table with a parameter labeled "ab" (whose row
                           # power_ab is the marginal power of that parameter)
                           # is read by its composite_power, and glance()
                           # keeps that table's power_ab as a column.
                           "power_ab")


# The first candidate term present in x, preferring one whose value is not
# missing: an unbalanced design leaves the single per-group size NA, so the
# recognizer must fall through to the total N. Returns the term name, or
# NA_character_ if none is present.
.ss_power_first_hit_term <- function(x, candidates) {
  present <- intersect(candidates, x$term)   # in candidate (priority) order
  if (!length(present)) return(NA_character_)
  vals <- x$value[match(present, x$term)]
  non_na <- present[!is.na(vals)]
  if (length(non_na)) non_na[1L] else present[1L]
}

.ss_power_first_hit <- function(x, candidates) {
  term <- .ss_power_first_hit_term(x, candidates)
  if (is.na(term)) NA_real_ else x$value[x$term == term][1L]
}


# Identify the sample-size row in a long ss_power_* table.
.ss_power_sample_size <- function(x) {
  .ss_power_first_hit(x, .SS_POWER_SIZE_TERMS)
}


# Identify the achieved-power row, if present.
.ss_power_actual_power <- function(x) {
  .ss_power_first_hit(x, .SS_POWER_POWER_TERMS)
}

# The name of the achieved-power row tidy() reads (NA_character_ if none),
# the one row glance() leaves out as the design's power.
.ss_power_power_term <- function(x) {
  .ss_power_first_hit_term(x, .SS_POWER_POWER_TERMS)
}


#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ss_power <- function(x, ...) {
  data.frame(
    term     = "sample_size",
    estimate = .ss_power_sample_size(x),
    power    = .ss_power_actual_power(x),
    stringsAsFactors = FALSE
  )
}


#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ss_power <- function(x, ...) {
  # One row with sample size + power + each input parameter that was
  # echoed in the long output (terms other than the size and power
  # rows themselves).
  size  <- .ss_power_sample_size(x)
  power <- .ss_power_actual_power(x)
  # Exclude only the one size row reported as `estimate` and the one power row
  # reported as `power`, the rows tidy() reads. Other design sizes, for
  # example the second group's n or the total N, are retained as columns so no
  # size is lost for an unbalanced design, and so are other powers: the
  # marginal powers of a composite design (power_a, power_b, and power_ab for
  # parameters labeled a, b, and ab) and the path powers of
  # ss_power_indirect_effect(). Through the first 1.1.0 builds every row named
  # in .SS_POWER_POWER_TERMS was left out, and when power_ab joined it for
  # ss_power_indirect_effect() the glance() of a composite table with a
  # parameter labeled "ab" lost its power_ab column.
  size_term  <- .ss_power_first_hit_term(x, .SS_POWER_SIZE_TERMS)
  power_term <- .ss_power_power_term(x)
  excluded   <- c(size_term, power_term)
  extras_rows <- x[!x$term %in% excluded, , drop = FALSE]
  out <- data.frame(
    term     = "sample_size",
    estimate = size,
    power    = power,
    stringsAsFactors = FALSE
  )
  if (nrow(extras_rows) > 0L) {
    extras <- as.list(extras_rows$value)
    names(extras) <- extras_rows$term
    out <- cbind(out, as.data.frame(extras,
                                    stringsAsFactors = FALSE))
  }
  out
}


# The AIPE planners answer a different question than the power planners: the
# smallest size at which a confidence interval is narrow enough, rather than
# at which a test is powerful enough. Their tidy() therefore reports the
# planned size beside the width the plan was made against. `width` is the
# desired full width the user asked for and nothing else: a planner that
# reports the expected width of its plan, or the width that an assurance
# secures (expected_width, assured_width), keeps that row under its own
# name, and glance() shows it as its own column. A plan made from a budget
# and not from a width has width NA.
.SS_AIPE_WIDTH_TERMS <- "width"

# The composite SEM planner reports one desired width per targeted
# parameter (width_<label>), so a single unsuffixed name may not exist;
# fall back to the first row whose name is "width_" and a label. The
# pct_ci_less_width_<label> and mean_width_<label> rows do not start with
# "width_".
.ss_aipe_width <- function(x) {
  w <- .ss_power_first_hit(x, .SS_AIPE_WIDTH_TERMS)
  if (!is.na(w)) return(w)
  hit <- grep("^width_", x$term)
  if (length(hit)) x$value[hit[1L]] else NA_real_
}

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ss_aipe <- function(x, ...) {
  data.frame(
    term     = "sample_size",
    estimate = .ss_power_sample_size(x),
    width    = .ss_aipe_width(x),
    stringsAsFactors = FALSE
  )
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ss_aipe <- function(x, ...) {
  size  <- .ss_power_sample_size(x)
  width <- .ss_aipe_width(x)
  size_term  <- .ss_power_first_hit_term(x, .SS_POWER_SIZE_TERMS)
  width_term <- .ss_power_first_hit_term(x, .SS_AIPE_WIDTH_TERMS)
  extras_rows <- x[!x$term %in% c(size_term, width_term), , drop = FALSE]
  out <- data.frame(
    term     = "sample_size",
    estimate = size,
    width    = width,
    stringsAsFactors = FALSE
  )
  if (nrow(extras_rows) > 0L) {
    extras <- as.list(extras_rows$value)
    names(extras) <- extras_rows$term
    out <- cbind(out, as.data.frame(extras, stringsAsFactors = FALSE))
  }
  # A plan's scalar attributes that qualify its size travel as columns
  # too, so the one-row summary says what the table's footer says:
  # ss_aipe_reliability() records the fraction of missing information
  # planned for and the size complete data would need (fmi,
  # necessary_N_complete), and, for an assurance plan, the samples of the
  # final pass that returned no interval and were counted as wide
  # (no_interval). Through the first versions of 1.1.0 the summary left
  # them out. ss_aipe_c() records the excess kurtosis of the errors an
  # assurance plan allowed for (true_kurtosis), which travels the same way. A
  # planner that sets none of them is unchanged.
  for (a in c("fmi", "necessary_N_complete", "no_interval", "true_kurtosis")) {
    v <- attr(x, a, exact = TRUE)
    if (is.numeric(v) && length(v) == 1L) out[[a]] <- v
  }
  out
}


# The alpha_lower and alpha_upper columns of an interval formed with
# unequal tail areas (the table's tail_areas attribute; see .tail_areas()
# in R/dmar_output_helpers.R), added to a tidy() or glance() frame `out`
# after its conf_level column, holding the two areas on the rows `rows`
# marks (the interval's rows) and NA on the others; `out` is returned
# unchanged when the table records none, so an equal-tailed interval's
# columns are those it has without the record of the split. Columns, not
# attributes, so that the split stays with its rows when results are
# stacked, as the level does: a
# stack built with rbind() from an equal-tailed result and an unequal one
# stops on their different columns rather than state either split for
# the other's rows, and one built with dplyr::bind_rows() holds NA, read
# as equal tails, on the equal-tailed rows. Not exported.
.add_tail_columns <- function(out, x, rows = rep(TRUE, nrow(out))) {
  tails <- .recorded_tail_areas(x)
  if (is.null(tails) || !nrow(out)) return(out)
  out$alpha_lower <- ifelse(rows, unname(tails[1L]), NA_real_)
  out$alpha_upper <- ifelse(rows, unname(tails[2L]), NA_real_)
  out
}

# Default verbs for any DMAR result table. Every tidy-returning DMAR function
# routes through .as_dmar_tbl(), so these two methods are the floor under the
# whole package. A long table (term beside a single value column) answers
# tidy() with the term/estimate view, with a conf_level column on the rows
# of its interval when it has one and the attributes that name its
# p-value rows (p_terms, and p_sentence, the p-value its sentence
# reports), and glance() with the one-row wide view, which adds, when the
# table has an interval, the name of its estimate's column (estimand) and
# its level (conf_level). A wide
# table (a leading label column beside several typed numeric columns)
# answers through the wide branch: tidy() builds term from the label
# column(s) and takes as estimate the table's estimate (or value) column
# when it has one, else, in a table that carries an interval, the first
# numeric column that is not a limit, standard error, test statistic,
# degrees of freedom, p-value, confidence level, or count (the column the
# interval belongs to, the rule results_sentence() follows, whose name it
# keeps in an estimand column), and else its first numeric column,
# passing the remaining columns through under their own names, with the
# conf_level the table records beside its interval;
# glance() returns the row count with the table's scalar attributes. Through
# 1.0.1 the estimate was always the first numeric column, so a row of
# simple_effects_AB() gave its F statistic as the estimate beside the
# interval for partial eta squared, and a table with an estimate column
# after a count (responder_analysis()) lost the count. Family classes
# (dmar_ss_power, dmar_ss_aipe, the CI and reliability families, and the
# bespoke wide tidiers below) sit ahead of dmar_tbl in the class vector,
# so their richer methods win wherever one exists and these are never
# shadowed by accident. This is what makes the promise in
# vignettes/dmar_output.Rmd, that the verbs answer everywhere, literally
# true.
#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_tbl <- function(x, ...) {
  if (!is.null(x$term) && !is.null(x$value)) {
    out <- data.frame(
      term     = as.character(x$term),
      estimate = as.numeric(x$value),
      stringsAsFactors = FALSE
    )
    # A table with an interval says in its rows which rows form it: a
    # conf_level column holds the level the table records on the
    # interval's estimate (the row its ci_estimate attribute names, or
    # else its first row that is not a limit) and on its two limits, and
    # NA on every other row. The rows carry it, so results_sentence()
    # states the level and the estimate of the rows it is given when
    # tidy() results are stacked with rbind() or subset; an attribute, as
    # the first versions of 1.1.0 used, stays with the first frame of a
    # stack and would be read for the rows of the others.
    # An interval held in rows named for its estimate (<estimate>_ci_lower
    # and <estimate>_ci_upper, as sem_fit_measures() holds the RMSEA's) is
    # marked the same way, at its own level when a <estimate>_ci_level row
    # records one, except the interval the table names in its ci_source
    # attribute: the interval for the standardized mean difference that
    # the intervals of cles(), nnt_from_smd(), and
    # proportion_of_superiority() are computed from, which those tables
    # echo and which stays NA, as it was when its limits were named
    # smd_lower and smd_upper.
    lim <- c("ci_lower", "ci_upper")
    level <- attr(x, "conf_level", exact = TRUE)
    level <- if (is.numeric(level) && length(level) == 1L &&
                 is.finite(level)) level else NA_real_
    row_level <- rep(NA_real_, nrow(out))
    has_ci <- FALSE
    if (all(lim %in% out$term)) {
      ci_rows <- out$term %in% c(.results_sentence_long_estimate(x), lim)
      row_level[ci_rows] <- level
      has_ci <- TRUE
    }
    for (e in setdiff(.results_sentence_row_estimates(x),
                      attr(x, "ci_source", exact = TRUE))) {
      own <- out$estimate[out$term == paste0(e, "_ci_level")]
      ci_rows <- out$term %in% paste0(e, c("", "_ci_lower", "_ci_upper"))
      row_level[ci_rows] <- if (length(own) && is.finite(own[1L])) own[1L]
                            else level
      has_ci <- TRUE
    }
    if (has_ci) {
      out$conf_level <- row_level
      out <- .add_tail_columns(out, x, rows = !is.na(row_level))
    }
    # The p-value rows are named by attributes, which stay with the rows
    # of a table of the same kind; a stack of tables whose p-value rows
    # are named otherwise gives no p-value rather than another's. The
    # ci_source attribute goes with them, for rows whose conf_level is NA
    # throughout (nnt_from_smd() given limits and no level).
    for (a in c("p_terms", "p_sentence", "ci_source")) {
      attr(out, a) <- attr(x, a, exact = TRUE)
    }
    return(out)
  }
  # Wide branch: label columns (character or factor) name the rows; the
  # estimate column is the primary quantity (see the comment above).
  is_num <- vapply(x, is.numeric, logical(1L))
  label_cols <- names(x)[!is_num]
  term <- if (length(label_cols) == 0L) {
    as.character(seq_len(nrow(x)))
  } else {
    do.call(paste, c(lapply(x[label_cols], as.character), sep = ":"))
  }
  num_cols <- names(x)[is_num]
  if (length(num_cols) == 0L) {
    return(data.frame(term = term, stringsAsFactors = FALSE))
  }
  has_ci <- !is.null(.results_sentence_limit_cols(names(x)))
  est_col <- intersect(c("estimate", "value"), num_cols)
  if (!length(est_col) && has_ci) {
    est_col <- num_cols[!grepl(.results_sentence_not_estimate, num_cols)]
  }
  est_col <- if (length(est_col)) est_col[1L] else num_cols[1L]
  out <- data.frame(term = term, estimate = as.numeric(x[[est_col]]),
                    stringsAsFactors = FALSE)
  # The estimate of an interval taken from a column with a name of its
  # own (kappa, D2, eta_squared_partial) keeps that name in an estimand
  # column, since term names the row (from the label columns, or the row
  # number when there are none) and not the quantity: results_sentence()
  # labels the row's sentence with it, as it labels the table's own. A
  # column, not an attribute, so that it stays with its rows in a stack.
  if (has_ci && !est_col %in% c("estimate", "value")) {
    out$estimand <- rep(est_col, nrow(out))
  }
  for (nm in setdiff(num_cols, est_col)) out[[nm]] <- x[[nm]]
  # The level of the interval beside it, NA when the table records none
  # (no interval was formed), so the columns depend on the table's columns
  # alone.
  level <- attr(x, "conf_level", exact = TRUE)
  if (has_ci && is.null(out$conf_level)) {
    out$conf_level <- rep(if (is.numeric(level) && length(level) == 1L)
      level else NA_real_, nrow(out))
  }
  if (has_ci) out <- .add_tail_columns(out, x)
  .dmar_tidy_adjustments(out, x)
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_tbl <- function(x, ...) {
  if (!is.null(x$term) && !is.null(x$value)) {
    # One row, one column per term. Repeated terms (a regions-of-significance
    # table reports one block per group pair) are disambiguated the way base R
    # does, with make.unique suffixes, rather than dropped.
    terms <- make.unique(as.character(x$term), sep = "_")
    vals <- as.list(as.numeric(x$value))
    names(vals) <- terms
    out <- as.data.frame(vals, stringsAsFactors = FALSE)
    # A table with an interval adds the name of the column that holds the
    # interval's estimate (the row its ci_estimate attribute names, or its
    # first row that is not a limit) in an estimand column, and its level
    # in a conf_level column, so that results_sentence() on the row states
    # the level and the estimate of the table's own sentence. Through the
    # first versions of 1.1.0 the row carried neither: its sentence was
    # refused for want of a level, and with the level given it named the
    # first column, the observed correlation of
    # correction_for_attenuation() beside the interval of the corrected
    # one. Columns, so that they stay with the row when glance() results
    # are stacked; the p-value the table's sentence reports (the row its
    # p_sentence attribute names, or else its first p_terms row) is named
    # by the p_sentence attribute of the result, read only where the row
    # has that column.
    lim <- c("ci_lower", "ci_upper")
    if (all(lim %in% x$term)) {
      # A table that reports the limits alone (ci_rc(), ci_snr(),
      # ci_srsnr()) has no estimate, so its estimand is NA; through the
      # first versions of 1.1.0 its glance() stopped with "replacement has
      # 0 rows", where 1.0.1 returned the row of its limits.
      est <- .results_sentence_long_estimate(x)
      if (!"estimand" %in% names(out)) {
        out$estimand <- if (is.null(est)) NA_character_
                        else names(out)[match(est, terms)]
      }
      if (!"conf_level" %in% names(out)) {
        level <- attr(x, "conf_level", exact = TRUE)
        out$conf_level <- if (is.numeric(level) && length(level) == 1L)
          level else NA_real_
      }
      out <- .add_tail_columns(out, x)
      p_own <- attr(x, "p_sentence", exact = TRUE)
      if (!(is.character(p_own) && length(p_own) == 1L &&
            p_own %in% terms)) {
        p_own <- terms[terms %in% attr(x, "p_terms", exact = TRUE)][1L]
      }
      if (!is.na(p_own)) {
        attr(out, "p_sentence") <- names(out)[match(p_own, terms)]
      }
    } else {
      # cles() and nnt_from_smd() hold their intervals in columns named for
      # their estimates, beside the interval for the standardized mean
      # difference they are computed from; the table's ci_source attribute
      # names that one, so that results_sentence() on the row sets it
      # aside, as it did when its limits were named smd_lower and
      # smd_upper.
      attr(out, "ci_source") <- attr(x, "ci_source", exact = TRUE)
    }
    return(out)
  }
  # Wide branch: the row count plus the scalar attributes every wide table
  # may carry (NA when absent, so the shape is predictable).
  .dmar_tidy_adjustments(.add_tail_columns(data.frame(
    n_terms    = nrow(x),
    conf_level = attr(x, "conf_level") %||% NA_real_,
    B_used     = attr(x, "B_used") %||% NA_integer_,
    stringsAsFactors = FALSE
  ), x), x)
}


# ---------------------------------------------------------------------------
# Bespoke wide-table tidiers: the validity and invariance surfaces whose
# rows have real structure (items, ladder rungs, groups). Each producing
# function tags its return with the leading class named here.
# ---------------------------------------------------------------------------

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_content_validity <- function(x, ...) {
  # The estimate is the item's I-CVI, which the term (the item) does not
  # name; the estimand column does, so that results_sentence() gives a row
  # the table's own sentence, "i_cvi = ...", as the default tidy() does for
  # a column with a name of its own.
  out <- data.frame(term = as.character(x$item),
                    estimate = as.numeric(x$i_cvi),
                    estimand = rep("i_cvi", length(x$item)),
                    stringsAsFactors = FALSE)
  for (nm in c("ci_lower", "ci_upper", "kappa", "cvr", "n_experts",
               "n_relevant")) {
    if (!is.null(x[[nm]])) out[[nm]] <- x[[nm]]
  }
  # The level of the intervals beside them, so that a tidy row states the
  # level of its limits.
  if (!is.null(out$ci_lower)) {
    out$conf_level <- rep(attr(x, "conf_level", exact = TRUE) %||% NA_real_,
                          nrow(out))
  }
  out
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_content_validity <- function(x, ...) {
  data.frame(
    s_cvi_ave  = attr(x, "s_cvi_ave") %||% NA_real_,
    s_cvi_ua   = attr(x, "s_cvi_ua") %||% NA_real_,
    n_items    = nrow(x),
    n_experts  = max(x$n_experts, na.rm = TRUE),
    conf_level = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_dmacs <- function(x, ...) {
  out <- data.frame(term = as.character(x$item),
                    estimate = as.numeric(x$dmacs),
                    stringsAsFactors = FALSE)
  for (nm in c("factor_loading_reference", "factor_loading_focal",
               "intercept_reference", "intercept_focal", "sd_pooled")) {
    if (!is.null(x[[nm]])) out[[nm]] <- x[[nm]]
  }
  out
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_dmacs <- function(x, ...) {
  data.frame(
    n_items    = nrow(x),
    reference  = attr(x, "reference") %||% NA_character_,
    focal      = attr(x, "focal") %||% NA_character_,
    mean_focal = attr(x, "mean_focal") %||% NA_real_,
    sd_focal   = attr(x, "sd_focal") %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_measurement_invariance <- function(x, ...) {
  out <- x
  names(out)[names(out) == "level"] <- "term"
  out$term <- as.character(out$term)
  class(out) <- "data.frame"
  attr(out, "fits") <- NULL
  rownames(out) <- NULL
  out
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_measurement_invariance <- function(x, ...) {
  data.frame(
    n_levels    = nrow(x),
    estimator   = attr(x, "estimator") %||% NA_character_,
    ordered     = isTRUE(attr(x, "ordered")),
    test        = attr(x, "test") %||% NA_character_,
    fit_indices = attr(x, "fit_indices") %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_measurement_alignment <- function(x, ...) {
  data.frame(
    term            = as.character(x$group),
    n               = x$n,
    factor_mean     = x$factor_mean,
    factor_variance = x$factor_variance,
    stringsAsFactors = FALSE
  )
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_measurement_alignment <- function(x, ...) {
  data.frame(
    n_groups            = nrow(x),
    alignment           = attr(x, "alignment") %||% NA_character_,
    epsilon             = attr(x, "epsilon") %||% NA_real_,
    simplicity_function = attr(x, "simplicity_function") %||% NA_real_,
    R2_loadings_mean    = mean(unlist(attr(x, "R2_loadings")),
                               na.rm = TRUE),
    R2_intercepts_mean  = mean(unlist(attr(x, "R2_intercepts")),
                               na.rm = TRUE),
    R2_total            = mean(unlist(attr(x, "R2_total")), na.rm = TRUE),
    converged           = isTRUE(attr(x, "converged") == 1) ||
                          isTRUE(attr(x, "converged")),
    n_starts            = attr(x, "n_starts") %||% NA_integer_,
    n_optima            = attr(x, "n_optima") %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

# factor_congruence(): one tidy() row per pair of factors (a factor without
# a partner keeps its row, its term naming the missing partner "(none)",
# as in "X3 | (none)", rather than pasting NA in as text, and its
# coefficient NA), and a glance() row of counts, the two arguments, and
# the mean coefficient over the matched pairs (signed under
# match = "order").

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_factor_congruence <- function(x, ...) {
  none <- function(f) ifelse(is.na(f), "(none)", f)
  data.frame(
    term         = paste(none(x$factor_x), none(x$factor_y), sep = " | "),
    estimate     = as.numeric(x$congruence),
    reflected    = as.logical(x$reflected),
    next_largest = as.numeric(x$next_largest),
    stringsAsFactors = FALSE
  )
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_factor_congruence <- function(x, ...) {
  C <- attr(x, "congruence_matrix", exact = TRUE)
  items <- attr(x, "items", exact = TRUE)
  matched <- !is.na(x$congruence)
  data.frame(
    n_items         = if (is.null(items)) NA_integer_ else length(items),
    n_factors_x     = if (is.null(C)) NA_integer_ else nrow(C),
    n_factors_y     = if (is.null(C)) NA_integer_ else ncol(C),
    n_matched       = sum(matched),
    match           = attr(x, "match", exact = TRUE) %||% NA_character_,
    loadings        = attr(x, "loadings", exact = TRUE) %||% NA_character_,
    mean_congruence = if (any(matched)) mean(x$congruence[matched]) else NA_real_,
    stringsAsFactors = FALSE
  )
}


# The AIPE sensitivity return schema, stated once. Every
# ss_aipe_*_sensitivity() member reports at least these rows: the mean /
# median / SD of the realized interval widths, the proportion of intervals
# at or below the planning width, the tail-specific empirical non-coverage
# of the population value, and the overall empirical Type I error rate.
# The proportion rows are on the 0 to 1 scale, so total_type_I_error is
# the sum of pct_ci_miss_low and pct_ci_miss_high. Beside these, each
# member reports mean_X / median_X / sd_X for its own estimand (mean_icc,
# mean_smd, ...) and echoes its inputs under their unit names (total_N or
# n_per_group for the evaluated size, true_X, supposed_X, width,
# conf_level, and, when an assurance was supplied, assurance). This vector is
# the single source of truth for the family contract; the family-wide
# test in tests/testthat/test-ss_aipe_sensitivity_family.R asserts every
# member against it.
.SS_AIPE_SENS_CORE_TERMS <- c("mean_ci_width", "median_ci_width",
                              "sd_ci_width", "pct_ci_less_width",
                              "pct_ci_miss_low", "pct_ci_miss_high",
                              "total_type_I_error")


# .ss_sens_simulation(G, has_interval, met_width, conf_level)
#
# The record of a sensitivity function's Monte Carlo run, stated once: the
# list every *_sensitivity() function attaches to its table as
# attr(out, "simulation"), so that the number of replications and the
# denominator of the summaries travel with the result.
#
# The rule it encodes (one denominator): the proportion of intervals that
# met the width target (the pct_ci_less_width row), the coverage rows, and the
# width summaries are all taken over the replications that produced an
# interval. A replication without an interval is in no summary row; it is
# counted, and the record gives the number attempted, the number with an
# interval, and the number without.
#   G            : the number of replications attempted (the function's G).
#   has_interval : a logical vector of length G, TRUE for a replication
#                  that returned a confidence interval (for a power
#                  sensitivity function, a test). No NA.
#   met_width    : a logical vector of length G, TRUE for a replication
#                  whose interval met the width target; NA counts as FALSE,
#                  and the entry of a replication without an interval is
#                  not read. NULL for a function with no width target (the
#                  power sensitivity functions, the sequential procedure),
#                  which leaves the two width entries NA.
#   conf_level   : the level of the Monte Carlo interval for the width
#                  proportion, the table's own conf_level (not read when
#                  met_width is NULL).
# Returns a list:
#   attempted                 G, as an integer.
#   returned                  replications with an interval.
#   failed                    replications without one (G - returned).
#   denominator               returned: the number of replications the
#                             proportion rows and the width and coverage
#                             summaries are taken over.
#   width_target_probability  the replications that met the width target
#                             over those with an interval (the
#                             pct_ci_less_width row); NA when none returned an
#                             interval or there is no width target.
#   width_target_mc_ci        the exact binomial (Clopper-Pearson) limits,
#                             from stats::binom.test(), for that proportion
#                             at conf_level. They describe simulation error
#                             only. NA where the proportion is NA.
# A function that records more about its run (the first failure message of
# ss_aipe_reliability_sensitivity(), the infinite widths of
# ss_aipe_cv_sensitivity()) adds its own elements to the list.
# The arguments come from the calling function's own vectors, so a wrong
# shape is a programming error and stops with a plain message.
.ss_sens_simulation <- function(G, has_interval, met_width, conf_level) {
  G <- .dmar_check_count(G, "G")
  if (!is.logical(has_interval) || length(has_interval) != G ||
      anyNA(has_interval)) {
    .dmar_refuse_argument(c(
      "'has_interval' must be a logical vector of length G with no ",
      "missing value."), arg = "has_interval")
  }
  n_met <- NULL
  if (!is.null(met_width)) {
    if (!is.logical(met_width) || length(met_width) != G) {
      .dmar_refuse_argument(
        "'met_width' must be NULL or a logical vector of length G.",
        arg = "met_width")
    }
    n_met <- sum(has_interval & met_width %in% TRUE)
  }
  .ss_sens_counts(G, sum(has_interval), n_met, conf_level)
}

# The list of .ss_sens_simulation() from counts: G attempted, n_returned
# with an interval, and n_met of those that met the width target (NULL
# without a width target).
.ss_sens_counts <- function(G, n_returned, n_met, conf_level) {
  out <- list(attempted = G, returned = n_returned,
              failed = G - n_returned,
              denominator = n_returned,
              width_target_probability = NA_real_,
              width_target_mc_ci = c(NA_real_, NA_real_))
  if (is.null(n_met) || n_returned == 0L) return(out)
  out$width_target_probability <- n_met / n_returned
  out$width_target_mc_ci <- unname(stats::binom.test(
    n_met, n_returned, conf.level = conf_level)$conf.int)
  out
}

# .ss_sens_record(out, G, has_interval, conf_level)
#
# Attaches the record to a finished summary table and returns the table.
# The count that met the width target is read back from the table's own
# pct_ci_less_width row (that proportion times the number of replications with
# an interval), so the record and the row cannot disagree, whichever
# comparison (<= or <) the function's row uses. A row that is not a whole
# number of replications over those with an interval breaks the
# one-denominator rule; that is a programming error and stops. A table
# with no pct_ci_less_width row, or with NA there (no width target was given),
# gets the counts alone.
.ss_sens_record <- function(out, G, has_interval, conf_level = NULL) {
  G <- .dmar_check_count(G, "G")
  if (!is.logical(has_interval) || length(has_interval) != G ||
      anyNA(has_interval)) {
    .dmar_refuse_argument(c(
      "'has_interval' must be a logical vector of length G with no ",
      "missing value."), arg = "has_interval")
  }
  n_returned <- sum(has_interval)
  pct <- out$value[out$term == "pct_ci_less_width"]
  n_met <- NULL
  if (length(pct) == 1L && is.finite(pct) && n_returned > 0L) {
    n_met <- round(pct * n_returned)
    if (abs(n_met - pct * n_returned) > 1e-6) {
      .dmar_refuse_argument(c(
        "The 'pct_ci_less_width' row is not a proportion of the ",
        "replications with an interval."), arg = NULL)
    }
    n_met <- as.integer(n_met)
  }
  attr(out, "simulation") <- .ss_sens_counts(G, n_returned, n_met,
                                             conf_level)
  out
}


# The Monte Carlo sensitivity siblings (ss_power_R2_sensitivity,
# ss_power_reg_coef_sensitivity) report two powers, the empirical (simulated)
# power and the analytic power at the same planned N, and the comparison of the
# two is the point of the study. They therefore carry their own
# dmar_ss_power_sensitivity class: tidy() puts both powers beside the planned
# sample size, and glance() adds the simulated estimator distribution.
.SS_POWER_SENS_HEAD <- c("total_N", "empirical_power", "analytic_power")

.ss_power_sens_head <- function(x) {
  v <- stats::setNames(x$value, x$term)
  data.frame(
    term            = "sensitivity",
    sample_size     = unname(v["total_N"]),
    empirical_power = unname(v["empirical_power"]),
    analytic_power  = unname(v["analytic_power"]),
    stringsAsFactors = FALSE
  )
}


#' @rdname dmar_tidiers
#' @exportS3Method generics::tidy
tidy.dmar_ss_power_sensitivity <- function(x, ...) {
  .ss_power_sens_head(x)
}


#' @rdname dmar_tidiers
#' @exportS3Method generics::glance
glance.dmar_ss_power_sensitivity <- function(x, ...) {
  out    <- .ss_power_sens_head(x)
  v      <- stats::setNames(x$value, x$term)
  extras <- v[setdiff(names(v), .SS_POWER_SENS_HEAD)]
  if (length(extras) > 0L)
    out <- cbind(out, as.data.frame(as.list(extras), stringsAsFactors = FALSE))
  out
}
