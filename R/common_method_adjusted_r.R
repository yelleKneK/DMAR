#' Correlations Adjusted for Common Method Variance, With Intervals and Tests
#'
#' @description
#' When every variable in a study is measured the same way, often by one
#' respondent on one questionnaire, part of each correlation can come from
#' the shared method rather than from the constructs. The marker variable
#' technique of Lindell and Whitney (2001) takes a correlation that should
#' be zero on substantive grounds, the correlation \eqn{r_M} of a marker
#' variable with the substantive variables, as an estimate of the method
#' effect common to every correlation, and removes it from each
#' substantive correlation. \code{common_method_adjusted_r()} returns, for
#' one or more values of \eqn{r_M}, every adjusted correlation with its
#' confidence interval and Lindell and Whitney's \emph{t} test, so that a
#' reader can see whether a relationship survives the adjustment and how
#' that conclusion depends on the marker correlation assumed.
#'
#' The marker correlation is itself an estimate, and the proxy used when
#' no marker was measured (the smallest positive correlation among the
#' substantive variables) is a choice. Reading the adjusted correlations
#' over a grid of values of \eqn{r_M}, from no adjustment to a generous
#' one, is the treatment of that uncertainty: a correlation whose interval
#' excludes zero across the plausible range does not depend on the value
#' chosen. \code{\link{common_method_marker}} reports the adjusted matrix
#' for a single marker correlation without intervals or tests, since it
#' does not take \eqn{N}; this function takes \eqn{N} and gives both.
#'
#' @param R A correlation matrix among the substantive variables (not the
#'   marker): symmetric, with a unit diagonal. Its column (or row) names
#'   name the variables; an unnamed matrix gets \code{V1}, \code{V2}, and
#'   so on.
#' @param N The sample size on which \code{R} rests, at least 4.
#' @param marker_r One or more marker correlations \eqn{r_M}, each in
#'   \eqn{(-1, 1)}, in the order the rows should follow. When \code{NULL}
#'   (the default) the smallest positive off-diagonal correlation in
#'   \code{R} is used as the proxy marker, as in
#'   \code{\link{common_method_marker}}.
#' @param conf_level Confidence level of the intervals (default 0.95).
#'
#' @details
#' \strong{The adjustment.} Each substantive correlation \eqn{r_{ij}} is
#' adjusted as
#' \deqn{r^{A}_{ij} = \frac{r_{ij} - r_M}{1 - r_M}}
#' (Lindell & Whitney, 2001, Equation 4), the partial correlation of
#' \eqn{i} and \eqn{j} with the method effect removed when that effect
#' contributes the same \eqn{r_M} to every correlation. The adjusted
#' matrix for each value of \code{marker_r} is the \code{"adjusted"}
#' attribute; for a single value it equals the \code{"adjusted"}
#' attribute of \code{\link{common_method_marker}}.
#'
#' \strong{The test.} Lindell and Whitney (2001, Equation 5) test each
#' adjusted correlation as a partial correlation with one variable
#' partialed out,
#' \deqn{t = \frac{r^{A}}{\sqrt{(1 - (r^{A})^2) / (N - 3)}},}
#' referred to the \emph{t} distribution with \eqn{N - 3} degrees of
#' freedom, with a two-sided \emph{p}-value.
#'
#' \strong{The interval.} With \eqn{r_M} treated as known, the adjustment
#' is an increasing linear function of the population correlation
#' \eqn{\rho_{ij}}, so a confidence interval for \eqn{\rho_{ij}} maps
#' into one for the adjusted correlation with the same coverage: each
#' limit of the Fisher's \emph{Z} interval for \eqn{\rho_{ij}} (the
#' interval of \code{\link{ci_r}}, based on the variance \eqn{1/(N - 3)})
#' is adjusted by the same equation as the correlation. A lower limit
#' that the map carries below \eqn{-1} is set to \eqn{-1}, the bound of a
#' correlation. Because the map is increasing, that happens exactly when
#' the lower limit for \eqn{\rho_{ij}} falls below \eqn{2 r_M - 1},
#' whatever the sign of \eqn{r_{ij}} or of the adjusted correlation:
#' never when \eqn{r_M \le 0}; for \eqn{0 < r_M < .5} only when the
#' lower limit for \eqn{\rho_{ij}} is negative, as it is for a negative
#' correlation, one near zero, or one from a small sample (with
#' \eqn{r_{ij} = .35}, \eqn{N = 5}, and \eqn{r_M = .3} the adjusted
#' correlation is .071 and its lower limit is set to \eqn{-1}); and for
#' \eqn{r_M \ge .5} even when that limit is positive.
#' The interval is conditional on \eqn{r_M}: it carries the sampling
#' variability of \eqn{r_{ij}} and none of that of the marker
#' correlation, which is what the grid over \code{marker_r} is for.
#'
#' An adjusted correlation at or beyond \eqn{\pm 1} lies outside the
#' range of a correlation. It arises from a correlation below
#' \eqn{2 r_M - 1}, which a method effect of \eqn{r_M} common to every
#' correlation cannot produce (or from a correlation of exactly 1), and
#' it gets no interval and no test: its \code{ci_lower},
#' \code{ci_upper}, \code{t_value}, and \code{p_value} are \code{NA}.
#'
#' The interval and the test are different approximations. The interval
#' excludes zero exactly when the Fisher's \emph{Z} interval for
#' \eqn{\rho_{ij}} excludes \eqn{r_M}, while the test treats the adjusted
#' correlation as a sample partial correlation, so near the boundary of
#' significance the two can disagree. The test is kept for parity with
#' the paper; the interval is what to report beside the estimate. Neither
#' is adjusted for the number of pairs or of marker values.
#'
#' \strong{Reading the results.} The proxy marker is one of the
#' substantive correlations, so the pair that supplied it is adjusted to
#' exactly zero and its test says nothing about that pair. The technique
#' presumes the variables are reflected so that their intercorrelations
#' are positive; a negative correlation is pushed further from zero rather
#' than toward it, so reverse code as needed before computing \code{R}.
#' A correlation that remains different from zero after the adjustment is
#' consistent with a relationship that is not only a method artifact, under
#' the assumption that the method effect is the same for every pair; it is
#' not evidence that no method variance is present.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   value of \code{marker_r} and pair of variables, the pairs in the
#'   order (1, 2), (1, 3), ..., (2, 3), ... within each marker value, and
#'   the columns \code{marker_r}, \code{variable_1}, \code{variable_2},
#'   \code{r} (the unadjusted correlation), \code{estimate} (the adjusted
#'   correlation), \code{ci_lower}, \code{ci_upper}, \code{t_value},
#'   \code{df} (\eqn{N - 3}), and \code{p_value}. The attributes are
#'   \code{adjusted} (a list of the adjusted correlation matrices, one per
#'   value of \code{marker_r}, named by it), \code{N}, \code{conf_level},
#'   and \code{marker_rule} (\code{"smallest_positive"} when the proxy was
#'   used, \code{"supplied"} otherwise). \code{tidy()} returns one row per
#'   pair and marker value, with \code{term} naming the pair
#'   (\code{"variable_1:variable_2"}) and the columns \code{marker_r},
#'   \code{r}, \code{estimate}, \code{ci_lower}, \code{ci_upper},
#'   \code{statistic} (the \emph{t}), \code{df}, \code{p_value}, and
#'   \code{conf_level}; \code{glance()} returns \code{n_pairs},
#'   \code{n_markers}, \code{N}, and \code{conf_level}.
#'
#' @references
#' Lindell, M. K., & Whitney, D. J. (2001). Accounting for common method
#'   variance in cross-sectional research designs. \emph{Journal of Applied
#'   Psychology, 86}(1), 114--121. \doi{10.1037/0021-9010.86.1.114}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{common_method_marker}} for the adjusted matrix and
#'   its summary at one marker correlation;
#'   \code{\link{common_method_single_factor}} for the single-factor
#'   screen; \code{\link{ci_r}} for the interval of an unadjusted
#'   correlation.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Visual, textual, and speed composites of the nine Holzinger and
#' # Swineford tests in the primary set, each the mean of three tests.
#' hs <- holzinger_swineford
#' scores <- data.frame(
#'   visual  = rowMeans(hs[c("t1_visual_perception", "t2_cubes",
#'                           "t4_lozenges")]),
#'   textual = rowMeans(hs[c("t6_paragraph_comprehension", "t7_sentence",
#'                           "t9_word_meaning")]),
#'   speed   = rowMeans(hs[c("t10_addition", "t12_counting_groups_of_dots",
#'                           "t13_straight_and_curved_capitals")]))
#' R <- stats::cor(scores)
#'
#' # The smallest positive correlation as the proxy marker; the pair that
#' # supplied it (textual with speed) is adjusted to zero. Visual with
#' # speed shows the interval and the test parting near the boundary: the
#' # interval includes zero while the p-value is below .05.
#' common_method_adjusted_r(R, N = nrow(scores))
#'
#' # A grid of marker correlations, from no adjustment to a large one.
#' grid <- common_method_adjusted_r(R, N = nrow(scores),
#'                                  marker_r = c(0, 0.1, 0.2, 0.3))
#' grid[grid$variable_1 == "visual" & grid$variable_2 == "textual", ]
#' attr(grid, "adjusted")[["0.2"]]
#'
#' @export
common_method_adjusted_r <- function(R, N, marker_r = NULL,
                                     conf_level = 0.95) {
  if (!is.matrix(R) || !is.numeric(R) || nrow(R) < 2L ||
      nrow(R) != ncol(R) || anyNA(R) || !isSymmetric(unname(R)) ||
      any(abs(R) > 1 + 1e-8) || any(abs(diag(R) - 1) > 1e-8)) {
    .dmar_refuse_argument(c(
      "'R' must be a correlation matrix of two or more variables: ",
      "numeric, symmetric, with a unit diagonal and no missing values."),
      arg = "R")
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 4) {
    .dmar_refuse_argument(c(
      "'N' must be a single sample size of at least 4; the test has ",
      "N - 3 degrees of freedom and the interval rests on the ",
      "variance 1/(N - 3)."), arg = "N")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number between 0 and 1.",
      arg = "conf_level")
  }
  p <- nrow(R)
  vars <- colnames(R) %||% rownames(R) %||% paste0("V", seq_len(p))
  dimnames(R) <- list(vars, vars)

  # The pairs (i, j), i < j, in the order (1, 2), (1, 3), ..., (2, 3).
  pairs <- which(upper.tri(R), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  r <- R[pairs]

  marker_rule <- "supplied"
  if (is.null(marker_r)) {
    # The proxy marker of common_method_marker(): the smallest positive
    # off-diagonal correlation.
    positive <- r[r > 0]
    if (!length(positive)) {
      .dmar_refuse_argument(c(
        "No positive correlation is available to serve as a proxy ",
        "marker; supply 'marker_r'."), arg = "marker_r")
    }
    marker_r <- min(positive)
    marker_rule <- "smallest_positive"
  }
  if (!is.numeric(marker_r) || length(marker_r) < 1L ||
      any(!is.finite(marker_r)) || any(marker_r <= -1 | marker_r >= 1)) {
    .dmar_refuse_argument(
      "'marker_r' must hold one or more numbers in (-1, 1).", arg = "marker_r")
  }

  # Fisher's Z interval for each population correlation (as in ci_r()),
  # mapped below through the same adjustment as the correlation.
  df <- N - 3
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  z <- atanh(r)
  rho_lower <- tanh(z - crit / sqrt(df))
  rho_upper <- tanh(z + crit / sqrt(df))

  rows <- lapply(marker_r, function(m) {
    estimate <- (r - m) / (1 - m)
    # The map carries a lower limit below -1 exactly when the lower limit
    # for rho is below 2 m - 1, whatever the sign of r (for m > 0 a small
    # sample is enough); it is cut at the bound of a correlation. An
    # adjusted correlation at or beyond -1 or 1 lies outside the range of
    # a correlation and gets no limits and no test.
    inside <- abs(estimate) < 1
    lower <- rep(NA_real_, length(estimate))
    upper <- rep(NA_real_, length(estimate))
    t_value <- rep(NA_real_, length(estimate))
    lower[inside] <- pmax((rho_lower[inside] - m) / (1 - m), -1)
    upper[inside] <- (rho_upper[inside] - m) / (1 - m)
    t_value[inside] <- estimate[inside] /
      sqrt((1 - estimate[inside]^2) / df)
    data.frame(marker_r   = m,
               variable_1 = vars[pairs[, 1L]],
               variable_2 = vars[pairs[, 2L]],
               r          = r,
               estimate   = estimate,
               ci_lower   = lower,
               ci_upper   = upper,
               t_value    = t_value,
               df         = df,
               p_value    = 2 * stats::pt(-abs(t_value), df),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  adjusted <- lapply(marker_r, function(m) {
    a <- (R - m) / (1 - m)
    diag(a) <- 1
    a
  })
  names(adjusted) <- as.character(marker_r)

  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      subclass = "dmar_common_method_adjusted_r")
  attr(out, "adjusted") <- adjusted
  attr(out, "N") <- N
  attr(out, "marker_rule") <- marker_rule
  out
}


# The tidy() and glance() methods for common_method_adjusted_r(), documented
# on its own page. The bespoke pair exists because the generic wide tidier
# would take the first numeric column, marker_r, as the estimate (and then
# overwrite it with the estimate column), dropping the marker value that
# tells apart the rows of one pair.

#' @rdname common_method_adjusted_r
#' @param x A \code{dmar_common_method_adjusted_r} table returned by
#'   \code{common_method_adjusted_r()}.
#' @param ... Unused.
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_common_method_adjusted_r <- function(x, ...) {
  data.frame(term       = paste(x$variable_1, x$variable_2, sep = ":"),
             marker_r   = x$marker_r,
             r          = x$r,
             estimate   = x$estimate,
             ci_lower   = x$ci_lower,
             ci_upper   = x$ci_upper,
             statistic  = x$t_value,
             df         = x$df,
             p_value    = x$p_value,
             conf_level = attr(x, "conf_level") %||% NA_real_,
             stringsAsFactors = FALSE,
             row.names  = NULL)
}

#' @rdname common_method_adjusted_r
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_common_method_adjusted_r <- function(x, ...) {
  data.frame(n_pairs    = length(unique(paste(x$variable_1, x$variable_2,
                                              sep = ":"))),
             n_markers  = length(attr(x, "adjusted")),
             N          = attr(x, "N") %||% NA_real_,
             conf_level = attr(x, "conf_level") %||% NA_real_,
             stringsAsFactors = FALSE,
             row.names  = NULL)
}
