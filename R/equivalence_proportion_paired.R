#' Equivalence and Noninferiority Tests for Two Paired Proportions
#'
#' Tests whether the difference between two proportions measured on the
#' same units, \eqn{\psi = p_x - p_y}, lies inside equivalence bounds or
#' above a noninferiority margin: whether a cheaper diagnostic procedure
#' finds a condition about as often as the standard one on the same
#' patients, or whether a short form of a scale classifies respondents
#' about as often as the full form. The tests are the score tests of Nam
#' (1997), which take the place of the \emph{t} tests of
#' \code{\link{equivalence_c}}, and the decisions are read from the
#' score interval of \code{\link{ci_proportion_paired}} (Tango, 1998),
#' with which they agree.
#'
#' @details
#' \strong{The tests.} With the score statistic \eqn{Z(\psi_0)} described
#' on \code{\link{ci_proportion_paired}}, the noninferiority null
#' hypothesis \eqn{H_0: \psi \le -\delta_L} is rejected in favor of
#' \eqn{\psi > -\delta_L} when \eqn{Z(-\delta_L) > z_{1-\alpha}}, and the
#' nonsuperiority null hypothesis \eqn{H_0: \psi \ge \delta_U} is
#' rejected in favor of \eqn{\psi < \delta_U} when
#' \eqn{Z(\delta_U) < -z_{1-\alpha}}. Equivalence, \eqn{-\delta_L < \psi
#' < \delta_U}, is claimed when both are rejected. Because the claim
#' requires both rejections, the two one-sided tests form an
#' intersection-union test (Berger & Hsu, 1996): each is carried out at
#' level \eqn{\alpha} and the equivalence claim has level \eqn{\alpha},
#' with no adjustment for there being two tests. The \emph{p}-value of the
#' equivalence test, \code{p_tost}, is the larger of the two one-sided
#' \emph{p}-values.
#'
#' \strong{Why the score statistic and not McNemar's.} McNemar's test
#' estimates the variance of the per-unit difference as it would be if
#' \eqn{\psi} were 0, so it tests no difference and nothing else. A test
#' against a margin needs the variance at \eqn{\psi = -\delta_L} (or
#' \eqn{\delta_U}), which depends on a discordant probability that the
#' hypothesis does not fix; the score statistic replaces that probability
#' with its maximum likelihood estimate under the hypothesis (Nam, 1997).
#' At \eqn{\psi_0 = 0} the score statistic is McNemar's, so the tests here
#' generalize it. Nam (1997) found that a Wald-type statistic, with the
#' variance estimated from the observed proportions instead, rejects more
#' often than its nominal level, while the score test holds close to it.
#'
#' \strong{Reading the interval.} \eqn{Z} decreases in \eqn{\psi_0}, and
#' the score interval is the set of differences the score test does not
#' reject. So the noninferiority test rejects at level \eqn{\alpha}
#' exactly when the lower limit of the 100(1 - 2\eqn{\alpha})\% score
#' interval exceeds \eqn{-\delta_L}, and the nonsuperiority test exactly
#' when its upper limit is below \eqn{\delta_U}. The flags are computed
#' from that interval, which the table reports.
#'
#' \strong{The flags and the verdict.} The flags read the interval
#' against the bounds:
#' \itemize{
#'   \item \code{equivalent}: the interval lies inside
#'     \eqn{(-\delta_L, \delta_U)}.
#'   \item \code{noninferior}: the lower limit exceeds \eqn{-\delta_L}.
#'   \item \code{superior}: the lower limit exceeds \eqn{\delta_U}.
#'   \item \code{inferior}: the upper limit is below \eqn{-\delta_L}.
#'   \item \code{nonsuperior}: the upper limit is below \eqn{\delta_U}, so
#'     a difference in favor of \code{x} as large as \eqn{\delta_U} is
#'     ruled out. It is the mirror image of noninferiority.
#' }
#' With both bounds, the \code{"verdict"} attribute names the most
#' informative statement the interval supports: \code{"Equivalent"},
#' \code{"Superior"}, \code{"Inferior"}, \code{"Noninferior only"} (an
#' interval above \eqn{-\delta_L} that contains \eqn{\delta_U}),
#' \code{"Not relevantly superior only"} (an interval below
#' \eqn{\delta_U} that contains \eqn{-\delta_L}), or
#' \code{"Inconclusive"} (an interval that contains both bounds). An
#' inconclusive result is a statement about precision, not evidence of a
#' difference.
#'
#' \strong{Noninferiority alone.} A call with \code{prop_diff_lower} and no
#' \code{prop_diff_upper} asks only whether \eqn{\psi > -\delta_L}, and no
#' upper bound is assumed for it. The rows that would need one,
#' \code{z_upper}, \code{p_upper}, \code{p_tost}, \code{prop_diff_upper},
#' \code{equivalent}, \code{superior}, and \code{nonsuperior}, are
#' \code{NA}, so the table keeps its shape; the verdict is
#' \code{"Noninferior"} (the lower limit exceeds \eqn{-\delta_L}),
#' \code{"Inferior"} (the upper limit is below \eqn{-\delta_L}), or
#' \code{"Inconclusive"}; and the \emph{p}-value a results sentence
#' reports is \code{p_noninferiority}. This is the noninferiority test of
#' \code{\link{equivalence_c}}, which reads \eqn{-\delta_L} alone. A
#' call with \code{prop_diff_upper} alone gives the symmetric bounds
#' \eqn{(-\delta_U, \delta_U)}, as \code{\link{equivalence_c}} does.
#'
#' \strong{Choosing the bounds.} The bounds are fixed before the data are
#' examined, on substantive grounds, as the smallest difference in the
#' two proportions that would matter, and on the scale of the difference
#' itself: \code{prop_diff_lower = 0.05} is five percentage points.
#'
#' @inheritParams ci_proportion_paired
#' @param prop_diff_lower,prop_diff_upper Bounds on the scale of the
#'   difference in proportions, each a single number in (0, 1); the
#'   equivalence region is \eqn{(-\delta_L, \delta_U)}, with
#'   \code{prop_diff_lower} as \eqn{\delta_L} and \code{prop_diff_upper}
#'   as \eqn{\delta_U}. Supply at least one. With
#'   \code{prop_diff_upper} alone the bounds are symmetric. With
#'   \code{prop_diff_lower} alone the call is a noninferiority test with
#'   margin \eqn{\delta_L}, and the rows that would need an upper bound
#'   are \code{NA} (see Details).
#' @param alpha_level One-sided significance level of each test. Defaults
#'   to 0.05, so the decisions are read from the 90\% score interval.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{difference} (\eqn{\hat\psi}{psi-hat}); the counts \code{n_10},
#'   \code{n_01}, and \code{N}; the score statistics \code{z_lower}
#'   (\eqn{Z(-\delta_L)}) and \code{z_upper} (\eqn{Z(\delta_U)}) and their
#'   one-sided \emph{p}-values \code{p_lower} and \code{p_upper}; the TOST
#'   \emph{p}-value \code{p_tost}; the noninferiority \emph{p}-value
#'   \code{p_noninferiority} (equal to \code{p_lower} by construction);
#'   the 100(1 - 2\eqn{\alpha})\% score limits \code{ci_lower} and
#'   \code{ci_upper}; the bounds \code{prop_diff_lower} (stored signed, as
#'   \eqn{-\delta_L}) and \code{prop_diff_upper}; and the flags
#'   \code{equivalent}, \code{noninferior}, \code{superior},
#'   \code{inferior}, and \code{nonsuperior} (1 = declared, 0 = not). A
#'   noninferiority test alone, with no \code{prop_diff_upper}, has the same
#'   rows, with \code{NA} in those that would need an upper bound. The
#'   attributes are \code{conf_level} (\eqn{1 - 2\alpha}), the
#'   \code{verdict}, and \code{p_sentence}, which names the
#'   \emph{p}-value \code{\link{results_sentence}} reports beside the
#'   interval: \code{p_tost}, or \code{p_noninferiority} for a
#'   noninferiority test alone.
#'
#' @references
#' Berger, R. L., & Hsu, J. C. (1996). Bioequivalence trials,
#'   intersection-union tests and equivalence confidence sets.
#'   \emph{Statistical Science, 11}(4), 283--319.
#'   \doi{10.1214/ss/1032280304}
#'
#' McNeil, B. J., Keeler, E., & Adelstein, S. J. (1975). Primer on
#'   certain elements of medical decision making. \emph{New England
#'   Journal of Medicine, 293}(5), 211--215.
#'   \doi{10.1056/NEJM197507312930501}
#'
#' Nam, J.-m. (1997). Establishing equivalence of two treatments and
#'   sample size requirements in matched-pairs design. \emph{Biometrics,
#'   53}(4), 1422--1430. \doi{10.2307/2533508}
#'
#' Tango, T. (1998). Equivalence test and confidence interval for the
#'   difference in proportions for the paired-sample design.
#'   \emph{Statistics in Medicine, 17}(8), 891--908.
#'   \doi{10.1002/(SICI)1097-0258(19980430)17:8<891::AID-SIM780>3.0.CO;2-B}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_proportion_paired}} for the interval alone;
#'   \code{\link{equivalence_c}} for a linear contrast of means.
#'
#' @family equivalence testing
#'
#' @family hypothesis tests
#'
#' @keywords htest
#'
#' @examples
#' # Liver scan data from McNeil, Keeler, and Adelstein (1975), the
#' # example of Nam (1997, p. 1427): 344 patients classified by scan (x)
#' # and by pathology (y), 32 positive by scan only and 27 by pathology
#' # only. Does the scan understate the rate of disease by less than five
#' # percentage points? That is a noninferiority question, so only
#' # prop_diff_lower is given, and the rows that would need an upper bound
#' # are NA. z_lower is the score statistic Nam reports, 2.795; the table
#' # rounds it for display and stores the full value.
#' res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
#'                                      prop_diff_lower = 0.05)
#' res
#' res$value[res$term == "z_lower"]
#' attr(res, "verdict")
#'
#' # Agreement within five points in both directions is a stronger
#' # claim. The interval rules out a shortfall of five points but not an
#' # excess of five, so the verdict is noninferiority only.
#' res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
#'                                      prop_diff_lower = 0.05,
#'                                      prop_diff_upper = 0.05)
#' attr(res, "verdict")
#'
#' # A short form (x) and the full form (y) of a screening scale given to
#' # the same 400 respondents: 38 screen positive on the short form only
#' # and 52 on the full form only. A shortfall of five points is not ruled
#' # out, but an excess of five is.
#' res <- equivalence_proportion_paired(n_10 = 38, n_01 = 52, N = 400,
#'                                      prop_diff_lower = 0.05,
#'                                      prop_diff_upper = 0.05)
#' res
#' attr(res, "verdict")
#'
#' @export
equivalence_proportion_paired <- function(x = NULL, y = NULL, n_10 = NULL,
                                          n_01 = NULL, N = NULL,
                                          prop_diff_lower = NULL,
                                          prop_diff_upper = NULL,
                                          alpha_level = 0.05) {
  if (is.null(prop_diff_lower) && is.null(prop_diff_upper))
    .dmar_refuse_argument(c(
      "Supply 'prop_diff_lower', 'prop_diff_upper', or both (the bounds on ",
      "the difference in proportions)."),
      arg = c("prop_diff_lower", "prop_diff_upper"))
  # prop_diff_lower alone asks the noninferiority question and nothing more,
  # so no upper bound is invented for it (the rows that would need one
  # are NA); prop_diff_upper alone gives symmetric bounds, as in
  # equivalence_c().
  noninferiority_only <- is.null(prop_diff_upper)
  for (nm in c("prop_diff_lower", "prop_diff_upper")) {
    v <- get(nm)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L ||
                        !is.finite(v) || v <= 0 || v >= 1))
      .dmar_refuse_argument(c(
        "'", nm, "' must be a single number in (0, 1), a bound on the ",
        "difference in proportions."), arg = NULL)
  }
  if (is.null(prop_diff_lower)) prop_diff_lower <- prop_diff_upper
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 0.5)
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 0.5).", arg = "alpha_level")
  k <- .paired_rd_counts(x, y, n_10, n_01, N)
  n_10 <- k[["n_10"]]
  n_01 <- k[["n_01"]]
  N    <- k[["N"]]
  est <- (n_10 - n_01) / N

  z_lower <- .paired_rd_score(n_10, n_01, N, -prop_diff_lower)
  p_lower <- stats::pnorm(z_lower, lower.tail = FALSE)
  # The noninferiority test is exactly the lower one-sided test, so it is
  # reported rather than recomputed, as equivalence_c() does.
  p_noninferiority <- p_lower

  lims <- .paired_rd_score_limits(n_10, n_01, N,
                                  stats::qnorm(1 - alpha_level))
  noninferior <- as.integer(lims[1L] > -prop_diff_lower)
  inferior    <- as.integer(lims[2L] < -prop_diff_lower)

  if (noninferiority_only) {
    z_upper <- p_upper <- p_tost <- NA_real_
    equivalent <- superior <- nonsuperior <- NA_integer_
    verdict <-
      if (noninferior == 1L) "Noninferior"
      else if (inferior == 1L) "Inferior"
      else "Inconclusive"
  } else {
    z_upper <- .paired_rd_score(n_10, n_01, N, prop_diff_upper)
    p_upper <- stats::pnorm(z_upper, lower.tail = TRUE)
    p_tost  <- max(p_lower, p_upper)
    equivalent  <- as.integer(noninferior == 1L && lims[2L] < prop_diff_upper)
    superior    <- as.integer(lims[1L] > prop_diff_upper)
    nonsuperior <- as.integer(lims[2L] < prop_diff_upper)
    # The order of the checks matters: superiority implies
    # noninferiority, and inferiority nonsuperiority, so the more
    # informative label wins.
    verdict <-
      if (equivalent == 1L) "Equivalent"
      else if (superior == 1L) "Superior"
      else if (inferior == 1L) "Inferior"
      else if (noninferior == 1L) "Noninferior only"
      else if (nonsuperior == 1L) "Not relevantly superior only"
      else "Inconclusive"
  }

  out <- data.frame(
    term  = c("difference", "n_10", "n_01", "N",
              "z_lower", "z_upper",
              "p_lower", "p_upper", "p_tost", "p_noninferiority",
              "ci_lower", "ci_upper",
              "prop_diff_lower", "prop_diff_upper",
              "equivalent", "noninferior", "superior", "inferior",
              "nonsuperior"),
    value = c(est, n_10, n_01, N,
              z_lower, z_upper,
              p_lower, p_upper, p_tost, p_noninferiority,
              lims[1L], lims[2L],
              -prop_diff_lower,
              if (noninferiority_only) NA_real_ else prop_diff_upper,
              equivalent, noninferior, superior, inferior,
              nonsuperior),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = 1 - 2 * alpha_level,
                      p_terms = c("p_lower", "p_upper", "p_tost",
                                  "p_noninferiority"))
  attr(out, "verdict") <- verdict
  # The p-value a results sentence reports beside the 1 - 2 alpha
  # interval is that of the question asked, equivalence or noninferiority
  # alone, not the first p_terms row.
  attr(out, "p_sentence") <- if (noninferiority_only) {
    "p_noninferiority"
  } else {
    "p_tost"
  }
  out
}
