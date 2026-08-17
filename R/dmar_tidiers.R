#' Tidy and Glance Methods for DMAR Result Tables
#'
#' @description
#' A DMAR function returns a tidy \code{data.frame} built to be read: one
#' row per quantity, a numeric \code{value} column, and a display layer
#' that rounds sensibly on the way to the console (see
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
#'   \code{dmar_ci_anova} \tab ANOVA effect size intervals \tab an effect size and its limits \cr
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
#'     \code{lower_limit} and \code{upper_limit} and, when the function
#'     reports one, an estimate row whose \code{term} is the name of the
#'     parameter. Carried by \code{\link{ci_r}},
#'     \code{\link{ci_smd_c}}, \code{\link{ci_pvaf}}, and
#'     \code{\link{ci_reg_coef}}.}
#'   \item{\code{dmar_ci_anova}}{Wide-format ANOVA effect size interval
#'     tables, with one row and columns for the effect name, the point
#'     estimate, the limits, and the design metadata. Carried by
#'     \code{\link{ci_eta_squared}},
#'     \code{\link{ci_eta_squared_partial}},
#'     \code{\link{ci_eta_squared_generalized}}, and
#'     \code{\link{ci_omega_squared}}.}
#' }
#'
#' Both produce a one-row \code{data.frame} with \code{term},
#' \code{estimate}, \code{ci_lower}, \code{ci_upper}, and, when the
#' object records it, \code{conf_level}. \code{glance()} on either class
#' calls \code{tidy()}, since the row is already the whole result.
#'
#' \strong{The post hoc family.}
#' \code{\link{ci_tukey_kramer}}, \code{\link{ci_games_howell}},
#' \code{\link{ci_scheffe}}, and \code{\link{ci_dunnett}} all carry
#' \code{dmar_post_hoc_ci}. Their source table is wide, with one row per
#' comparison: a \code{contrast} label, a point estimate
#' (\code{mean_difference} for the pairwise and many-to-one procedures,
#' \code{contrast_value} for Scheffe), a standard error, a test statistic,
#' the \code{lower_limit} and \code{upper_limit} of the simultaneous
#' interval, and the multiplicity-adjusted \code{p_adjusted}.
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
#' joint power to detect the indirect effect and whose component path
#' powers \code{glance()} carries as extra columns.
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
#' vectors when a planner introduces one.
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
#' res <- ci_r(r = 0.5, n = 50)
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
#' plan <- ss_power_smd(smd = 0.5, desired_power = 0.80)
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
  bound_terms <- c("lower_limit", "upper_limit",
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
  data.frame(
    term       = estimate_name,
    estimate   = estimate_val,
    ci_lower   = x$value[x$term == "lower_limit"],
    ci_upper  = x$value[x$term == "upper_limit"],
    conf_level = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
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
  # The wide ANOVA-style tables have columns:
  #   effect, <estimate_name>, lower_limit, upper_limit,
  #   [F_value, df_effect, df_error, N, ...].
  # The estimate column is whichever column is neither "effect" nor
  # a bound nor a metadata column.
  meta_cols <- c("effect", "lower_limit", "upper_limit",
                 "F_value", "df_effect", "df_error", "df_1", "df_2",
                 "N", "conf_level")
  est_col <- setdiff(colnames(x), meta_cols)[1L]
  if (is.na(est_col) || is.null(est_col)) est_col <- "estimate"
  data.frame(
    term       = est_col,
    estimate   = x[[est_col]][1L],
    ci_lower   = x$lower_limit[1L],
    ci_upper  = x$upper_limit[1L],
    conf_level = if (!is.null(x$conf_level)) x$conf_level[1L]
                 else attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
}


#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ci_anova <- function(x, ...) {
  tidy.dmar_ci_anova(x, ...)
}


# ---------------------------------------------------------------------------
# The post hoc family: dmar_post_hoc_ci
# ---------------------------------------------------------------------------

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_post_hoc_ci <- function(x, ...) {
  # The point estimate lives in `mean_difference` (Tukey-Kramer, Games-Howell,
  # Dunnett) or `contrast_value` (Scheffe); pick whichever this object carries.
  est_col <- intersect(c("mean_difference", "contrast_value"), names(x))[1L]
  data.frame(
    term        = x$contrast,
    estimate    = x[[est_col]],
    ci_lower    = x$lower_limit,
    ci_upper   = x$upper_limit,
    p_adjusted  = x$p_adjusted,
    conf_level  = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
}


#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_post_hoc_ci <- function(x, ...) {
  data.frame(
    n_contrasts = nrow(x),
    conf_level  = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
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
.SS_POWER_SIZE_TERMS <- c(# design-specific per-unit sizes: a per-cell,
                          # per-subject, or per-cluster count is the planning
                          # unit for a factorial, repeated measures, or cluster
                          # design, and is reported ahead of the design's total.
                          "necessary_n_per_cell", "approximate_n_per_cell",
                          "specified_n_per_cell",
                          "n_per_cell",
                          "necessary_n_subjects", "specified_n_subjects",
                          "necessary_J_per_arm", "specified_J_per_arm",
                          "necessary_n_clusters",
                          # per-group and total sizes shared across the family
                          "necessary_n_per_group", "necessary_N",
                          # the composite ANCOVA planners report a size resolved
                          # under unequal residual variances as approximate_*
                          # rather than necessary_*, because it is the smallest
                          # size at which an approximation reaches the target
                          # rather than one known to attain it (see the section
                          # on unequal residual variances in
                          # ?ss_power_composite_ancova_2group)
                          "approximate_n_per_group", "approximate_N",
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
                           "approximate_composite_power", "power")


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
  # Exclude only the one size row reported as `estimate` and every power row
  # (a design reports a single power here). Other design sizes, for example the
  # second group's n or the total N, are retained as columns so no size is lost
  # for an unbalanced design.
  size_term <- .ss_power_first_hit_term(x, .SS_POWER_SIZE_TERMS)
  excluded  <- c(size_term, .SS_POWER_POWER_TERMS)
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
# planned size beside the width the plan was made against, read through the
# same size registry and a width registry in planner-priority order (the
# desired width the user asked for, under the names the family uses for it).
.SS_AIPE_WIDTH_TERMS <- c("width", "width_target", "desired_width",
                          "expected_width")

# The composite SEM planner reports one desired width per targeted
# parameter (desired_width_<label>), so a single unsuffixed name may not
# exist; fall back to the first suffixed desired_width_ row.
.ss_aipe_width <- function(x) {
  w <- .ss_power_first_hit(x, .SS_AIPE_WIDTH_TERMS)
  if (!is.na(w)) return(w)
  hit <- grep("^desired_width_", x$term)
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
  out
}


# Default verbs for any DMAR result table. Every tidy-returning DMAR function
# routes through .as_dmar_tbl(), so these two methods are the floor under the
# whole package. A long table (term beside a single value column) answers
# tidy() with the term/estimate view and glance() with the one-row wide
# view. A wide table (a leading label column beside several typed numeric
# columns) answers through the wide branch: tidy() builds term from the
# label column(s) and takes the first numeric column as estimate, passing
# the remaining columns through under their own names, and glance()
# returns the row count with the table's scalar attributes. Family classes
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
    return(data.frame(
      term     = as.character(x$term),
      estimate = as.numeric(x$value),
      stringsAsFactors = FALSE
    ))
  }
  # Wide branch: label columns (character or factor) name the rows; the
  # first numeric column is the primary quantity.
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
  out <- data.frame(term = term, estimate = as.numeric(x[[num_cols[1L]]]),
                    stringsAsFactors = FALSE)
  for (nm in num_cols[-1L]) out[[nm]] <- x[[nm]]
  out
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_tbl <- function(x, ...) {
  if (!is.null(x$term) && !is.null(x$value)) {
    # One row, one column per term. Repeated terms (a regions-of-significance
    # table reports one block per group pair) are disambiguated the way base R
    # does, with make.unique suffixes, rather than dropped.
    vals <- as.list(as.numeric(x$value))
    names(vals) <- make.unique(as.character(x$term), sep = "_")
    return(as.data.frame(vals, stringsAsFactors = FALSE))
  }
  # Wide branch: the row count plus the scalar attributes every wide table
  # may carry (NA when absent, so the shape is predictable).
  data.frame(
    n_terms    = nrow(x),
    conf_level = attr(x, "conf_level") %||% NA_real_,
    B_used     = attr(x, "B_used") %||% NA_integer_,
    stringsAsFactors = FALSE
  )
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
  out <- data.frame(term = as.character(x$item),
                    estimate = as.numeric(x$i_cvi),
                    stringsAsFactors = FALSE)
  for (nm in c("ci_lower", "ci_upper", "kappa", "cvr", "n_experts",
               "n_relevant")) {
    if (!is.null(x[[nm]])) out[[nm]] <- x[[nm]]
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
  for (nm in c("lambda_reference", "lambda_focal", "nu_reference",
               "nu_focal", "sd_pooled")) {
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


# The AIPE sensitivity return schema, stated once. Every
# ss_aipe_*_sensitivity() member reports at least these rows: the mean /
# median / SD of the realized interval widths, the proportion of intervals
# at or below the planning width, the tail-specific empirical non-coverage
# of the population value, and the overall empirical Type I error rate.
# The proportion rows are on the 0 to 1 scale, so total_type_I_error is
# the sum of pct_ci_miss_low and pct_ci_miss_high. Beside these, each
# member reports mean_X / median_X / sd_X for its own estimand (mean_icc,
# mean_smd, ...) and echoes its inputs under their unit names (total_N or
# n_per_group for the evaluated size, true_X, estimated_X, width,
# conf_level, and, when an assurance was supplied, assurance). This vector is
# the single source of truth for the family contract; the family-wide
# test in tests/testthat/test-ss_aipe_sensitivity_family.R asserts every
# member against it.
.SS_AIPE_SENS_CORE_TERMS <- c("mean_ci_width", "median_ci_width",
                              "sd_ci_width", "pct_ci_less_w",
                              "pct_ci_miss_low", "pct_ci_miss_high",
                              "total_type_I_error")


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
