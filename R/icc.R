# Internal: type label aliases.
.canonicalize_icc_type <- function(type) {
  alias <- c(
    "1" = "ICC(1,1)", "1k" = "ICC(1,k)", "1,1" = "ICC(1,1)", "1,k" = "ICC(1,k)",
    "2" = "ICC(2,1)", "2k" = "ICC(2,k)", "2,1" = "ICC(2,1)", "2,k" = "ICC(2,k)",
    "3" = "ICC(3,1)", "3k" = "ICC(3,k)", "3,1" = "ICC(3,1)", "3,k" = "ICC(3,k)",
    "ICC(1,1)" = "ICC(1,1)", "ICC(1,k)" = "ICC(1,k)",
    "ICC(2,1)" = "ICC(2,1)", "ICC(2,k)" = "ICC(2,k)",
    "ICC(3,1)" = "ICC(3,1)", "ICC(3,k)" = "ICC(3,k)",
    "all" = "all"
  )
  out <- unname(alias[as.character(type)])
  if (any(is.na(out))) {
    bad <- unique(type[is.na(out)])
    stop("Unrecognized ICC type(s): ", paste(bad, collapse = ", "),
         ". Valid: 1, 2, 3, 1k, 2k, 3k, ICC(1,1), ..., ICC(3,k), all.",
         call. = FALSE)
  }
  if (any(out == "all")) out <- c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)",
                                  "ICC(1,k)", "ICC(2,k)", "ICC(3,k)")
  unique(out)
}


# Internal: ANOVA decomposition for an n x k matrix of ratings.
# Returns the mean squares and degrees of freedom needed for ICC formulas.
.icc_decomposition <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'x' must be a numeric matrix or data.frame (rows = subjects, ",
         "columns = raters/measurements).", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("Missing values are not supported by icc(); remove or impute first.",
         call. = FALSE)
  }
  n <- nrow(x); k <- ncol(x)
  if (n < 2L || k < 2L) {
    stop("'x' must have at least 2 rows (subjects) and 2 columns (raters).",
         call. = FALSE)
  }

  grand_mean <- mean(x)
  row_means  <- rowMeans(x)
  col_means  <- colMeans(x)

  ss_between_subjects <- k * sum((row_means - grand_mean)^2)
  ss_between_raters   <- n * sum((col_means - grand_mean)^2)
  ss_total            <- sum((x - grand_mean)^2)
  ss_within_subjects  <- ss_total - ss_between_subjects   # one-way decomposition
  ss_residual         <- ss_total - ss_between_subjects - ss_between_raters

  ms_R  <- ss_between_subjects / (n - 1L)                  # MS_subjects
  ms_C  <- ss_between_raters   / (k - 1L)                  # MS_raters
  ms_E  <- ss_residual         / ((n - 1L) * (k - 1L))     # MS_residual (two-way)
  ms_W  <- ss_within_subjects  / (n * (k - 1L))            # MS_within (one-way)

  list(n = n, k = k,
       ms_R = ms_R, ms_C = ms_C, ms_E = ms_E, ms_W = ms_W)
}


# Internal: compute one ICC type and its F-distribution-based CI. Formulas from
# Shrout & Fleiss (1979), Table 4 (for the F tests) and pages 425--426 (for the CIs).
.icc_one_type <- function(d, type, conf_level) {
  alpha <- 1 - conf_level
  n  <- d$n; k <- d$k
  ms_R <- d$ms_R; ms_C <- d$ms_C; ms_E <- d$ms_E; ms_W <- d$ms_W

  switch(type,
    "ICC(1,1)" = {
      value   <- (ms_R - ms_W) / (ms_R + (k - 1) * ms_W)
      F_obs   <- ms_R / ms_W
      df_1    <- n - 1
      df_2    <- n * (k - 1)
      F_L     <- F_obs / stats::qf(1 - alpha / 2, df_1, df_2)
      F_U     <- F_obs * stats::qf(1 - alpha / 2, df_2, df_1)
      lo      <- (F_L - 1) / (F_L + k - 1)
      hi      <- (F_U - 1) / (F_U + k - 1)
      p_val   <- stats::pf(F_obs, df_1, df_2, lower.tail = FALSE)
    },
    "ICC(1,k)" = {
      value   <- (ms_R - ms_W) / ms_R
      F_obs   <- ms_R / ms_W
      df_1    <- n - 1
      df_2    <- n * (k - 1)
      F_L     <- F_obs / stats::qf(1 - alpha / 2, df_1, df_2)
      F_U     <- F_obs * stats::qf(1 - alpha / 2, df_2, df_1)
      lo      <- 1 - 1 / F_L
      hi      <- 1 - 1 / F_U
      p_val   <- stats::pf(F_obs, df_1, df_2, lower.tail = FALSE)
    },
    "ICC(2,1)" = {
      value   <- (ms_R - ms_E) /
                 (ms_R + (k - 1) * ms_E + k * (ms_C - ms_E) / n)
      F_obs   <- ms_R / ms_E
      df_1    <- n - 1
      df_2    <- (n - 1) * (k - 1)
      # Shrout & Fleiss (1979) p. 426 formulas for ICC(2,1) CI (asymmetric).
      Fj      <- ms_C / ms_E
      vn      <- (k - 1) * (n - 1) *
                 ((k * value * Fj + n * (1 + (k - 1) * value) - k * value)^2)
      vd      <- (n - 1) * k^2 * value^2 * Fj^2 +
                 (n * (1 + (k - 1) * value) - k * value)^2
      v       <- vn / vd
      F_L     <- stats::qf(1 - alpha / 2, df_1, v)
      F_U     <- stats::qf(1 - alpha / 2, v, df_1)
      lo      <- n * (ms_R - F_L * ms_E) /
                 (F_L * (k * ms_C + (k * n - k - n) * ms_E) + n * ms_R)
      hi      <- n * (F_U * ms_R - ms_E) /
                 (k * ms_C + (k * n - k - n) * ms_E + n * F_U * ms_R)
      p_val   <- stats::pf(F_obs, df_1, df_2, lower.tail = FALSE)
    },
    "ICC(2,k)" = {
      icc1    <- .icc_one_type(d, "ICC(2,1)", conf_level)
      single  <- icc1$value
      lo_s    <- icc1$lower_limit
      hi_s    <- icc1$upper_limit
      value   <- (k * single)              / (1 + (k - 1) * single)
      lo      <- (k * lo_s)                / (1 + (k - 1) * lo_s)
      hi      <- (k * hi_s)                / (1 + (k - 1) * hi_s)
      F_obs   <- icc1$F_value
      df_1    <- icc1$df_1
      df_2    <- icc1$df_2
      p_val   <- icc1$p_value
    },
    "ICC(3,1)" = {
      value   <- (ms_R - ms_E) / (ms_R + (k - 1) * ms_E)
      F_obs   <- ms_R / ms_E
      df_1    <- n - 1
      df_2    <- (n - 1) * (k - 1)
      F_L     <- F_obs / stats::qf(1 - alpha / 2, df_1, df_2)
      F_U     <- F_obs * stats::qf(1 - alpha / 2, df_2, df_1)
      lo      <- (F_L - 1) / (F_L + k - 1)
      hi      <- (F_U - 1) / (F_U + k - 1)
      p_val   <- stats::pf(F_obs, df_1, df_2, lower.tail = FALSE)
    },
    "ICC(3,k)" = {
      value   <- (ms_R - ms_E) / ms_R
      F_obs   <- ms_R / ms_E
      df_1    <- n - 1
      df_2    <- (n - 1) * (k - 1)
      F_L     <- F_obs / stats::qf(1 - alpha / 2, df_1, df_2)
      F_U     <- F_obs * stats::qf(1 - alpha / 2, df_2, df_1)
      lo      <- 1 - 1 / F_L
      hi      <- 1 - 1 / F_U
      p_val   <- stats::pf(F_obs, df_1, df_2, lower.tail = FALSE)
    },
    stop("Unknown ICC type: ", type, call. = FALSE)
  )

  data.frame(
    type        = type,
    value       = value,
    lower_limit = lo,
    upper_limit = hi,
    F_value     = F_obs,
    df_1        = df_1,
    df_2        = df_2,
    p_value     = p_val,
    stringsAsFactors = FALSE,
    row.names   = NULL
  )
}


#' Intraclass Correlation Coefficients With Confidence Intervals
#'
#' Computes one or more of the six standard intraclass correlation coefficients
#' (ICC) of Shrout and Fleiss (1979) for an \eqn{n \times k} matrix of ratings
#' (rows = subjects, columns = raters/measurements), along with the
#' \emph{F}-distribution-based confidence interval for each.
#'
#' @param x An \eqn{n \times k} numeric matrix or \code{data.frame}: each row
#'   is a subject (target) and each column is a rater (or repeated
#'   measurement). Missing values are not supported; remove or impute first.
#' @param type Which ICC variant(s) to return. Aliases include \code{"1"},
#'   \code{"2"}, \code{"3"} (single rater versions of types 1, 2, 3) and
#'   \code{"1k"}, \code{"2k"}, \code{"3k"} (average-of-\eqn{k}-raters
#'   versions). The full names \code{"ICC(1,1)"}, \code{"ICC(2,1)"},
#'   \code{"ICC(3,1)"}, \code{"ICC(1,k)"}, \code{"ICC(2,k)"},
#'   \code{"ICC(3,k)"} are also accepted, as is \code{"all"}. Vectors are
#'   accepted; a row is returned for each type.
#' @param conf_level Confidence level for the interval (default \code{0.95}).
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   requested ICC type and columns \code{type}, \code{value} (point
#'   estimate), \code{lower_limit},
#'   \code{upper_limit}, \code{F_value}, \code{df_1}, \code{df_2}, and
#'   \code{p_value} (for the implicit null \eqn{H_0\!: \mathrm{ICC} = 0}).
#'
#' @details The six variants follow Shrout and Fleiss's (1979) classification:
#' \itemize{
#'   \item \code{ICC(1,1)}: one-way random model, single rater. Each subject
#'     is rated by a different (random) rater; appropriate when raters are
#'     not crossed with subjects.
#'   \item \code{ICC(1,k)}: one-way random model, average of \eqn{k} raters.
#'   \item \code{ICC(2,1)}: two-way random model, single rater, absolute
#'     agreement. Subjects \eqn{\times} raters fully crossed; both effects
#'     random.
#'   \item \code{ICC(2,k)}: two-way random model, average of \eqn{k} raters,
#'     absolute agreement.
#'   \item \code{ICC(3,1)}: two-way mixed model, single rater, consistency.
#'     Raters are fixed (the only ones of interest); rates differences in
#'     mean across raters are not penalized.
#'   \item \code{ICC(3,k)}: two-way mixed model, average of \eqn{k} raters,
#'     consistency.
#' }
#' Confidence intervals follow the \emph{F}-distribution-based formulas of
#' Shrout and Fleiss (1979, pp.\ 425--426); the \code{ICC(2,1)} interval uses
#' their asymmetric approximate-\emph{df} formulation.
#'
#' @references
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
#'   assessing rater reliability. \emph{Psychological Bulletin, 86}(2),
#'   420--428.
#'
#' McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some
#'   intraclass correlation coefficients. \emph{Psychological Methods, 1}(1),
#'   30--46. \doi{10.1037/1082-989X.1.1.30}
#'
#' @examples
#' # Shrout & Fleiss (1979), Table 1: 4 raters, 6 targets.
#' shrout_fleiss <- matrix(
#'   c(9, 2, 5, 8,
#'     6, 1, 3, 2,
#'     8, 4, 6, 8,
#'     7, 1, 2, 6,
#'     10, 5, 6, 9,
#'     6, 2, 4, 7),
#'   nrow = 6, byrow = TRUE,
#'   dimnames = list(paste0("Target_", 1:6),
#'                   paste0("Judge_",  1:4))
#' )
#'
#' icc(shrout_fleiss, type = "all")
#'
#' # Just the two-way mixed-model single-rater consistency ICC:
#' icc(shrout_fleiss, type = "ICC(3,1)")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_r}}, \code{\link{descriptives}}
#'
#' @keywords htest multivariate
#'
#' @family reliability
#'
#' @export
#' @import stats
icc <- function(x, type = "ICC(2,1)", conf_level = 0.95) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  types <- .canonicalize_icc_type(type)
  d     <- .icc_decomposition(x)

  rows <- lapply(types, function(t) .icc_one_type(d, t, conf_level))
  out  <- do.call(rbind, rows)
  rownames(out) <- NULL
  .as_dmar_tbl(out, conf_level = conf_level)
}
