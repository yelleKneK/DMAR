#' Responder Analysis: Who Cleared the Threshold, by Group
#'
#' The clinical and behavioral endpoint that mean differences hide: the
#' proportion of each group whose outcome reaches a meaningful threshold
#' (a minimal clinically important difference, a remission cut, a mastery
#' criterion). For each group the function reports the responder count and
#' proportion with a Wilson confidence interval; with exactly two groups it
#' adds the risk difference with the Newcombe (1998) score-based hybrid
#' interval and the number needed to treat; and across any number of
#' groups it reports the omnibus chi square test of equal responder
#' proportions. An optional \code{sweep} repeats the analysis over a grid
#' of thresholds, making explicit how conclusions depend on where the line
#' is drawn, disclosing the threshold dependence that any
#' single-threshold claim leaves implicit.
#'
#' @param x Numeric vector of outcomes (for example, change scores).
#' @param group Group labels, one per observation (coerced to factor; the
#'   first level is the reference for the two-group difference).
#' @param threshold The cut defining response.
#' @param direction \code{"ge"} (default): a responder has
#'   \code{x >= threshold}; \code{"le"}: \code{x <= threshold} (for
#'   outcomes where lower is better).
#' @param conf_level Confidence level for all intervals. Defaults to 0.95.
#' @param sweep Optional numeric vector of additional thresholds; the
#'   analysis is repeated at each and stacked with a leading
#'   \code{threshold} column.
#'
#' @details
#' Per-group intervals are Wilson score intervals
#' (\code{\link{ci_proportion}}). The two-group risk difference uses
#' Newcombe's method 10: the difference interval is assembled from the two
#' Wilson limits, which keeps it inside [-1, 1] and well behaved at
#' boundary counts. The number needed to treat is \eqn{1/|\Delta|}, with
#' its interval from the inverted difference limits when the difference
#' interval excludes zero; when it includes zero the NNT interval is
#' reported as \code{NA} (the interval is disjoint and an interval on the
#' NNT scale would mislead; Altman, 1998). Dichotomizing throws away
#' information, so a responder analysis complements, never replaces, the
#' analysis of the continuous outcome (Maxwell, Delaney, & Kelley, 2027).
#'
#' @return A tidy wide \code{data.frame} (class \code{dmar_tbl}). One row
#'   per group with \code{group}, \code{n}, \code{responders},
#'   \code{estimate} (the proportion), \code{lower_limit},
#'   \code{upper_limit}; with two groups, a \code{difference} row
#'   (second level minus first) and an \code{nnt} row; and a final
#'   \code{omnibus} row carrying \code{chi_square}, \code{df}, and
#'   \code{p_value} (columns that are \code{NA} on the other rows). When
#'   \code{sweep} is supplied, the same table is stacked per threshold
#'   with a leading \code{threshold} column.
#'
#' @references
#' Altman, D. G. (1998). Confidence intervals for the number needed to
#'   treat. \emph{BMJ, 317}(7168), 1309--1312. \doi{10.1136/bmj.317.7168.1309}
#'
#' Newcombe, R. G. (1998). Interval estimation for the difference between
#'   independent proportions: Comparison of eleven methods.
#'   \emph{Statistics in Medicine, 17}(8), 873--890.
#'   \doi{10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_proportion}} for the per-group interval;
#'   \code{\link{nnt_from_smd}} for the model-based route to the number
#'   needed to treat from a standardized mean difference;
#'   \code{\link{cliff_delta}} and \code{\link{proportion_of_superiority}}
#'   for dominance-style effect sizes on the continuous outcome.
#'
#' @family effect size estimates
#'
#' @keywords htest design
#'
#' @examples
#' # A two-arm trial: change scores, response defined as a gain of 10+.
#' set.seed(113)
#' change <- c(rnorm(60, 8, 9), rnorm(60, 13, 9))
#' arm    <- rep(c("control", "treatment"), each = 60)
#' responder_analysis(change, arm, threshold = 10)
#'
#' # How threshold-dependent is that conclusion?
#' responder_analysis(change, arm, threshold = 10, sweep = c(5, 15))
#'
#' @export
#' @importFrom stats chisq.test qnorm
responder_analysis <- function(x, group, threshold,
                               direction = c("ge", "le"),
                               conf_level = 0.95, sweep = NULL) {
  direction <- match.arg(direction)
  if (!is.numeric(x) || length(x) < 4L) {
    stop("'x' must be a numeric vector with at least 4 observations.",
         call. = FALSE)
  }
  if (length(group) != length(x)) {
    stop("'group' must have one label per observation.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("'threshold' must be a single number.", call. = FALSE)
  }
  if (!is.null(sweep) && (!is.numeric(sweep) || anyNA(sweep))) {
    stop("'sweep' must be a numeric vector of thresholds.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  ok <- !is.na(x) & !is.na(group)
  x <- x[ok]; g <- droplevels(factor(group[ok]))
  if (nlevels(g) < 2L) {
    stop("'group' must have at least two groups after removing missing ",
         "values.", call. = FALSE)
  }

  one_threshold <- function(thr) {
    resp <- if (direction == "ge") x >= thr else x <= thr
    n_g <- tapply(resp, g, length)
    r_g <- tapply(resp, g, sum)
    p_g <- r_g / n_g
    z   <- qnorm(1 - (1 - conf_level) / 2)
    wilson <- function(r, n) {
      p <- r / n; den <- 1 + z^2 / n
      ctr <- (p + z^2 / (2 * n)) / den
      hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
      c(ctr - hw, ctr + hw)
    }
    lims <- t(vapply(seq_along(n_g), function(i) wilson(r_g[i], n_g[i]),
                     numeric(2)))
    out <- data.frame(
      group       = names(n_g),
      n           = as.numeric(n_g),
      responders  = as.numeric(r_g),
      estimate    = as.numeric(p_g),
      lower_limit = lims[, 1],
      upper_limit = lims[, 2],
      chi_square  = NA_real_,
      df          = NA_real_,
      p_value     = NA_real_,
      stringsAsFactors = FALSE
    )

    if (nlevels(g) == 2L) {
      d  <- p_g[2] - p_g[1]
      l1 <- wilson(r_g[1], n_g[1]); l2 <- wilson(r_g[2], n_g[2])
      # Newcombe method 10: combine the per-group score limits.
      lo <- d - sqrt((p_g[2] - l2[1])^2 + (l1[2] - p_g[1])^2)
      hi <- d + sqrt((l2[2] - p_g[2])^2 + (p_g[1] - l1[1])^2)
      out <- rbind(out, data.frame(
        group = "difference", n = NA_real_, responders = NA_real_,
        estimate = as.numeric(d), lower_limit = lo, upper_limit = hi,
        chi_square = NA_real_, df = NA_real_, p_value = NA_real_,
        stringsAsFactors = FALSE))
      excl0 <- lo > 0 || hi < 0
      out <- rbind(out, data.frame(
        group = "nnt", n = NA_real_, responders = NA_real_,
        estimate = 1 / abs(d),
        lower_limit = if (excl0) 1 / max(abs(lo), abs(hi)) else NA_real_,
        upper_limit = if (excl0) 1 / min(abs(lo), abs(hi)) else NA_real_,
        chi_square = NA_real_, df = NA_real_, p_value = NA_real_,
        stringsAsFactors = FALSE))
    }

    tab <- table(g, factor(resp, levels = c(FALSE, TRUE)))
    om  <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    rbind(out, data.frame(
      group = "omnibus", n = NA_real_, responders = NA_real_,
      estimate = NA_real_, lower_limit = NA_real_, upper_limit = NA_real_,
      chi_square = unname(om$statistic), df = unname(om$parameter),
      p_value = om$p.value, stringsAsFactors = FALSE))
  }

  if (is.null(sweep)) {
    res <- one_threshold(threshold)
  } else {
    thr_all <- unique(c(threshold, sweep))
    res <- do.call(rbind, lapply(thr_all, function(t) {
      cbind(threshold = t, one_threshold(t))
    }))
    rownames(res) <- NULL
  }
  .as_dmar_tbl(res, conf_level = conf_level)
}
