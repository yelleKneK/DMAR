# Internal helpers for the bootstrap intervals of the bifactor indices
# (ci_bifactor_indices()). The bootstrap itself is the reliability
# family's engine, .bootstrap_ci() in R/reliability_internals.R, which
# takes the bifactor point function below as a point function of several
# statistics; the helpers here say which cells of a bifactor_indices()
# table the bootstrap estimates, assemble the rows the fit used, refit the
# model to a resample, and write the lines the printed table carries.
# None is exported.


# The columns of a bifactor_indices() table that hold indices, in the
# table's order.
.BIFACTOR_INDEX_COLUMNS <- c("ECV", "omega", "omega_H", "omega_HS", "PUC",
                             "H")


# .bifactor_cells(tab)
#
# The cells of a bifactor_indices() table that apply to its rows, in the
# order of its rows and columns, each row followed by the maximal
# reliability of its factor (the table's "maximal_reliability"
# attribute): for the general factor (the first row) ECV, omega,
# omega_H, PUC, H, and maximal_reliability; for each group factor ECV,
# omega, omega_HS, H, and maximal_reliability. The set follows from the
# structure alone (which row is the general factor), never from which
# values happen to be finite, so the rows of ci_bifactor_indices() depend
# on the model and on nothing else.
#
# Inputs:
#   tab : a bifactor_indices() table.
# Returns:
#   A data frame with columns `factor` and `index`.
.bifactor_cells <- function(tab) {
  f <- as.character(tab$factor)
  general <- c("ECV", "omega", "omega_H", "PUC", "H", "maximal_reliability")
  group <- c("ECV", "omega", "omega_HS", "H", "maximal_reliability")
  data.frame(factor = c(rep(f[1L], length(general)),
                        rep(f[-1L], each = length(group))),
             index = c(general, rep(group, length(f) - 1L)),
             stringsAsFactors = FALSE)
}


# .bifactor_cell_values(tab, cells)
#
# The values of a bifactor_indices() table at the given cells, as a plain
# numeric vector in the order of `cells` (NA for a cell the table does not
# hold). A maximal_reliability cell is read from the table's attribute of
# that name.
.bifactor_cell_values <- function(tab, cells) {
  rows <- match(cells$factor, as.character(tab$factor))
  maximal <- attr(tab, "maximal_reliability", exact = TRUE)
  vapply(seq_len(nrow(cells)), function(k) {
    if (cells$index[k] == "maximal_reliability") {
      v <- if (is.null(maximal)) NA_real_ else maximal[cells$factor[k]]
      return(if (length(v) == 1L) as.numeric(v) else NA_real_)
    }
    col <- tab[[cells$index[k]]]
    if (is.null(col) || is.na(rows[k])) NA_real_ else as.numeric(col[rows[k]])
  }, numeric(1))
}


# .bifactor_maximal_cells(cells)
#
# Which of `cells` are maximal reliabilities, a logical vector.
.bifactor_maximal_cells <- function(cells) {
  cells$index == "maximal_reliability"
}


# .bifactor_boot_pool(fit)
#
# The rows the fit used, which the bootstrap resamples: the data lavaan
# stored with the fit, less any row with no observed value on the model's
# variables. Under listwise deletion (lavaan's default, and the default
# for ordered items) lavaan stores only the complete rows; under
# missing = "ml" it stores every row of the data, including a row with
# nothing observed, which it reads but leaves out of the sample, so that
# row is dropped here and every other row keeps its pattern of observed
# and missing values. The count is checked against the cases lavaan says
# it used, so the bootstrap and the point estimate describe the same
# sample.
#
# Inputs:
#   fit : a single-group lavaan fit to raw data.
# Returns:
#   A list with `data` (a data frame of the model's observed variables,
#   one row per case the fit used) and `incomplete` (the number of those
#   rows with a missing value).
.bifactor_boot_pool <- function(fit) {
  X <- tryCatch(lavaan::lavInspect(fit, "data"), error = function(e) NULL)
  if (!is.matrix(X) || !nrow(X) || !ncol(X)) {
    .dmar_refuse_argument(c(
      "'fit' carries no raw data, as a fit to a covariance matrix ",
      "does not; a bootstrap resamples the cases, so refit the model ",
      "to the raw data (lavaan::cfa(..., data = )) and pass that fit."),
      arg = "fit")
  }
  seen <- rowSums(!is.na(X))
  pool <- as.data.frame(X[seen > 0L, , drop = FALSE])
  rownames(pool) <- NULL
  used <- as.integer(sum(lavaan::lavInspect(fit, "nobs")))
  if (nrow(pool) != used) {
    .dmar_refuse_argument(
      sprintf(paste0(
   "The rows stored with 'fit' (%d with an observed value) are not ",
   "the %d cases the fit used, so a bootstrap of them would not ",
   "describe this fit."), nrow(pool), used), arg = "fit")
  }
  list(data = pool, incomplete = sum(rowSums(is.na(pool)) > 0L))
}


# .bifactor_refit_options(fit, denominator)
#
# The fit's own options (estimator, missing data treatment, scaling,
# ordered thresholds, optimizer, ...) for refitting its model to a
# resample, with the parts no replication reads switched off: standard
# errors and test statistics (robust maximum likelihood and WLSMV have
# the estimates of maximum likelihood and DWLS and differ only in those,
# so their replications are the cheaper fits, with the same estimates),
# the baseline model, and, unless the observed denominator needs the
# unrestricted (h1) covariance matrix (under missing = "ml" the full
# information estimate behind lavaan's sample statistics), the
# unrestricted model.
.bifactor_refit_options <- function(fit, denominator) {
  opts <- lavaan::lavInspect(fit, "options")
  opts$se <- "none"
  opts$test <- "none"
  opts$baseline <- FALSE
  if (denominator == "model_implied") opts$h1 <- FALSE
  opts
}


# .bifactor_point_fn(model, options, ordered, categories, general,
#                    denominator, cells)
#
# The point function the bootstrap engine applies to each resample, in
# the form .bootstrap_ci() documents for a point function of several
# statistics: the model is refit to the resample with the fit's parameter
# table and options (default starting values, as for the fit itself), and
# the indices at `cells` are read from bifactor_indices() on the refit.
# Built in its own small frame, so a socket worker receives these values
# and not the frame of ci_bifactor_indices(), which holds the fit.
#
# A refit that did not converge returns NA with the reason
# "not_converged". A refit whose likelihood has no maximum (an error
# variance that runs off below zero, the likelihood flat or still rising,
# where lavaan's optimizer stopped somewhere along the way;
# .lavaan_no_maximum()) stops with a short
# refusal of class dmar_no_maximum, which the engine counts among the
# replications that raised an error and quotes, naming the variable
# (.no_max_replication()). A refit whose solution is improper
# (bifactor_indices() flags a standardized loading beyond one or a
# negative standardized residual variance, or lavaan's post-estimation
# check fails, made on the model's own parameters when the fit has
# auxiliary variables, which the refit inherits through the fit's
# options; .sem_admissible(), R/sem_fit_internals.R) returns its
# indices with the reason "improper", which the engine uses under
# improper = "keep" and sets aside under "drop"; when such a solution
# leaves an index that is not finite or lies outside [0, 1] (an omega
# above one is a composite whose implied error variance is negative, and
# H is NA for a factor with a standardized loading at or past one), it
# yields no coefficient, and NA with the same reason sets it aside under
# either rule. The maximal reliabilities are not among the indices that
# test: bifactor_indices() gives none for an improper solution, so they
# are NA there and are marked not defined for the replication (the
# engine's dmar_undefined attribute), and the replication is used or set
# aside on its other indices alone. Their intervals therefore rest on
# the proper replications, the same ones under either rule. A resample
# in which an ordered item lost a category raises an error that says so,
# since the item's thresholds cannot all be refit; the engine counts it
# with that message. lavaan's warnings are muffled inside a replication,
# since its fate is what the engine reports.
#
# Inputs:
#   model       : the fit's parameter table without its estimates,
#                 standard errors, and starting values.
#   options     : .bifactor_refit_options().
#   ordered     : the fit's ordered variables (character(0) for none).
#   categories  : for each ordered variable, the categories observed in
#                 the rows the fit used.
#   general     : the name of the general factor.
#   denominator : "model_implied" or "observed".
#   cells       : the cells to return (.bifactor_cells(), PUC left out,
#                 and the maximal reliabilities too when the point
#                 estimate is improper and has none).
# Returns:
#   A function of a resample returning nrow(cells) numbers.
.bifactor_point_fn <- function(model, options, ordered, categories, general,
                               denominator, cells) {
  force(model)
  force(options)
  force(ordered)
  force(categories)
  force(general)
  force(denominator)
  force(cells)
  function(d) {
    d <- as.data.frame(d)
    for (v in ordered) {
      lost <- setdiff(categories[[v]], d[[v]])
      if (length(lost)) {
        .dmar_refuse_argument(
          sprintf(paste0("the resample drew no response in category %s ",
                         "of the ordered item %s, so its thresholds ",
                         "cannot all be refit"),
                  format(lost[1L]), v), arg = NULL)
      }
    }
    muffle <- function(w) invokeRestart("muffleWarning")
    none <- rep(NA_real_, nrow(cells))
    refit <- withCallingHandlers(
      lavaan::lavaan(model = model, data = d,
                     ordered = if (length(ordered)) ordered,
                     slotOptions = options),
      warning = muffle)
    if (!isTRUE(lavaan::lavInspect(refit, "converged"))) {
      attr(none, "dmar_reason") <- "not_converged"
      return(none)
    }
    # A refit whose likelihood has no maximum stops, and the engine counts
    # it among the replications that raised an error
    # (.no_max_replication()).
    no_max <- .lavaan_no_maximum(refit)
    if (!is.null(no_max)) .no_max_replication(no_max)
    tab <- withCallingHandlers(
      bifactor_indices(refit, general = general, denominator = denominator),
      warning = muffle)
    post_ok <- withCallingHandlers(.sem_admissible(refit, explain = FALSE),
                                   warning = muffle)
    vals <- .bifactor_cell_values(tab, cells)
    maximal <- .bifactor_maximal_cells(cells)
    if (isTRUE(attr(tab, "improper")) || !isTRUE(post_ok)) {
      other <- vals[!maximal]
      if (!all(is.finite(other)) || any(other < 0 | other > 1)) {
        vals <- none
      }
      vals[maximal] <- NA_real_
      attr(vals, "dmar_reason") <- "improper"
    }
    if (any(maximal & !is.finite(vals))) {
      attr(vals, "dmar_undefined") <- maximal & !is.finite(vals)
    }
    vals
  }
}


# .bifactor_improper_note(improper)
#
# The sentence the engine's warning gives about the replications set
# aside as improper, in the words of a bifactor model (the engine's own
# sentence speaks of a single-factor loading).
.bifactor_improper_note <- function(improper) {
  paste0(
    "An improper replication is a bifactor solution with a standardized ",
    "loading beyond one or a negative residual variance; setting these ",
    "aside conditions the intervals on admissible solutions",
    if (improper == "keep") {
      paste0("; under improper = \"keep\" an improper replication is set ",
             "aside only when an index other than the maximal ",
             "reliability is not finite or falls outside [0, 1], and since ",
             "the maximal reliability is not defined for an improper ",
             "solution, its intervals rest on the proper replications ",
             "alone")
    } else {
      ""
    })
}


# .bifactor_definition_lines(maximal = NULL)
#
# What the cells of a bifactor table that are most often misread mean:
# ECV and omega on a group factor's row describe that factor's items (the
# general factor's share of their common variance, and the reliability of
# their total score), and H is the index of Rodriguez, Reise, and
# Haviland (2016), not the maximal reliability. With `maximal` (the
# "maximal_reliability" attribute of bifactor_indices()) the second line
# gives the maximal reliability of each factor, as a point estimate, and
# names ci_bifactor_indices() for its interval; without it (the table of
# ci_bifactor_indices(), whose maximal_reliability rows carry the
# intervals), it says what those rows are.
#
# Inputs:
#   maximal : NULL, or a numeric vector named by factor (NA throughout
#             for an improper solution).
# Returns:
#   Two strings, named "bifactor_definitions".
.bifactor_definition_lines <- function(maximal = NULL) {
  ecv <- paste0(
    "For a group factor, ECV is the general factor's share of the ",
    "common variance of that factor's items (one minus it is the group ",
    "factor's share), and omega is the reliability of their total score")
  h <- paste0(
    "H is the index of Rodriguez, Reise, and Haviland (2016), which ",
    "counts the other factors' common variance as uncorrelated error; ")
  h <- if (is.null(maximal)) {
    paste0(h, "maximal_reliability is the reliability the optimally ",
           "weighted items reach as a measure of the factor under the ",
           "fitted model")
  } else if (!any(is.finite(maximal))) {
    paste0(h, "the maximal reliability the fitted model implies is not ",
           "defined for an improper solution")
  } else {
    paste0(h, "the maximal reliability the fitted model implies is ",
           paste(names(maximal), formatC(maximal, format = "f", digits = 3),
                 collapse = ", "),
           " (attribute \"maximal_reliability\"; ci_bifactor_indices() ",
           "gives each an interval)")
  }
  stats::setNames(c(ecv, h), rep("bifactor_definitions", 2L))
}


# .bifactor_indices_footer(x)
#
# The lines a bifactor_indices() table prints under itself, beside the
# denominator line .dmar_tbl_footer() (R/dmar_output_helpers.R) writes:
# .bifactor_definition_lines() with the table's maximal reliabilities.
# A table is recognized by its "maximal_reliability" attribute, which
# only bifactor_indices() sets; a subset of the table, which loses its
# attributes, prints without the lines.
#
# Inputs:
#   x : a table.
# Returns:
#   A character vector named "bifactor_definitions", or character(0).
.bifactor_indices_footer <- function(x) {
  maximal <- attr(x, "maximal_reliability", exact = TRUE)
  if (inherits(x, "dmar_ci_bifactor") || !is.numeric(maximal) ||
      is.null(names(maximal))) {
    return(character(0))
  }
  .bifactor_definition_lines(maximal)
}


# .ci_bifactor_footer(x)
#
# The lines a ci_bifactor_indices() table prints under itself, beside the
# confidence level, the denominator, and the bootstrap counts that
# .dmar_tbl_footer() (R/dmar_output_helpers.R) writes for every table: how
# the intervals were formed, over which rows, what was done with improper
# replications, which intervals fell back from BCa to the percentile
# interval, why PUC has none, and, when the maximal reliability rests on
# fewer replications than the other indices (the improper ones kept
# under improper = "keep" have none) or has no estimate (an improper
# point estimate), on how many or why; then the bootstrap estimates of
# the bias of H and of the maximal reliability (.bifactor_bias_line());
# then what ECV and omega mean for a group factor and which H the table
# carries. .dmar_tbl_footer() calls this for every table, so the console
# and a knitted report show the same text.
#
# Inputs:
#   x : a table (anything not of class dmar_ci_bifactor gives nothing).
# Returns:
#   A character vector, the line named "bifactor_interval", the line
#   named "bifactor_bias" when there is a bias to give, and then the two
#   .bifactor_definition_lines() (in the form that describes the
#   maximal_reliability rows), or character(0).
.ci_bifactor_footer <- function(x) {
  if (!inherits(x, "dmar_ci_bifactor")) return(character(0))
  requested <- attr(x, "ci_method_requested", exact = TRUE)
  if (is.null(requested)) return(character(0))
  method <- attr(x, "ci_method", exact = TRUE)
  bs <- attr(x, "bootstrap", exact = TRUE)
  rule <- if (is.list(bs)) bs$improper_rule else NULL
  N <- attr(x, "N", exact = TRUE)
  incomplete <- attr(x, "N_incomplete", exact = TRUE)
  line <- sprintf(
    "Intervals: %s bootstrap over the %s rows the fit used%s",
    if (requested == "bca") {
      "bias-corrected and accelerated (BCa)"
    } else {
      "percentile"
    },
    format(N),
    if (isTRUE(incomplete > 0L)) {
      sprintf(paste0(" (%d of them incomplete, resampled with their ",
                     "missing data patterns)"), incomplete)
    } else {
      ""
    })
  if (identical(rule, "drop")) {
    line <- paste0(line, ", with improper replications set aside")
  } else if (identical(rule, "keep")) {
    line <- paste0(line, ", with improper replications kept where every ",
                   "index but the maximal reliability is a proportion")
  }
  fell_back <- names(method)[startsWith(method, "percentile (")]
  if (requested == "bca" && length(fell_back)) {
    line <- paste0(line, sprintf(
      "; percentile for %s, where the BCa interval could not be formed",
      paste(fell_back, collapse = ", ")))
  }
  if ("PUC" %in% x$index) {
    line <- paste0(line, "; PUC depends only on which items share a ",
                   "group factor, so it has no interval")
  }
  maximal_rows <- x$index == "maximal_reliability"
  if (any(maximal_rows) && isTRUE(attr(x, "improper", exact = TRUE))) {
    line <- paste0(line, "; the maximal reliability is not defined for an ",
                   "improper solution, that of the fit included, so its ",
                   "rows have no estimate and no interval")
  } else if (any(maximal_rows) && is.list(bs) &&
             isTRUE(bs$improper_kept > 0L)) {
    t <- bs$t
    col <- if (is.matrix(t)) {
      which(endsWith(colnames(t), ":maximal_reliability"))[1L]
    } else {
      NA_integer_
    }
    if (!is.na(col)) {
      line <- paste0(line, sprintf(paste0(
        "; the maximal reliability is not defined for an improper ",
        "solution, so its intervals rest on the %d proper replications ",
        "of the %d used"), sum(is.finite(t[, col])), bs$returned))
    }
  }
  c(bifactor_interval = line,
    .bifactor_bias_line(if (is.list(bs)) bs$bias),
    .bifactor_definition_lines())
}


# .bifactor_bias_line(bias)
#
# The line a ci_bifactor_indices() table prints with the bootstrap
# estimates of the bias of H and of the maximal reliability of each
# factor, the two indices whose sample values overstate their population
# values on average (the help of ci_bifactor_indices() says why). The
# other indices' estimates are in the table's "bootstrap" attribute
# only.
#
# Inputs:
#   bias : the `bias` element of the "bootstrap" attribute, a numeric
#          vector named "factor:index", or NULL.
# Returns:
#   One string named "bifactor_bias", or character(0) when there is no
#   finite estimate to give.
.bifactor_bias_line <- function(bias) {
  if (!is.numeric(bias) || is.null(names(bias))) return(character(0))
  part <- function(index) {
    b <- bias[endsWith(names(bias), paste0(":", index)) & is.finite(bias)]
    if (!length(b)) return(NULL)
    paste0(index, " ",
           paste(sub(paste0(":", index, "$"), "", names(b)),
                 formatC(b, digits = 3, format = "g"), collapse = ", "))
  }
  parts <- c(part("H"), part("maximal_reliability"))
  if (!length(parts)) return(character(0))
  c(bifactor_bias = paste0(
    "Bootstrap estimates of bias (mean of the replicates used minus the ",
    "estimate): ", paste(parts, collapse = "; ")))
}
