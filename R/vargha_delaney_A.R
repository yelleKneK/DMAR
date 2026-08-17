#' Vargha and Delaney's \emph{A} (Stochastic-Superiority Effect Size)
#'
#' Computes Vargha and Delaney's (2000) \emph{A}, the probability that a
#' randomly drawn observation from the first sample exceeds a randomly drawn
#' observation from the second (with tied pairs counted as half), together
#' with its asymptotic standard error from DeLong, DeLong, and Clarke-Pearson
#' (1988) and a confidence interval on the population \emph{A}. \emph{A} is
#' equivalent to the receiver-operating-characteristic area-under-the-curve
#' (AUC) and to the common-language effect size (McGraw & Wong, 1992) under a
#' continuous-response assumption, and is a robust, scale-free ordinal effect
#' size that does not require equal variances or normality.
#'
#' @param x Either a numeric vector of observations from group 1, or a
#'   two-sided formula of the form \code{outcome ~ group} (in which case
#'   \code{data} must be supplied and the grouping variable must have exactly
#'   two levels).
#' @param y Numeric vector of observations from group 2. Ignored when \code{x}
#'   is a formula.
#' @param data Optional data frame containing the variables named in the
#'   formula \code{x}.
#' @param conf_level Confidence coverage for a symmetric interval (default
#'   \code{0.95}).
#' @param ci_method Either \code{"logit"} (default; Wald interval on the
#'   logit of \emph{A} with back-transformation, recommended for finite
#'   samples; Newcombe, 2006b) or \code{"wald"} (untransformed Wald on the
#'   original scale, clipped to [0, 1]).
#'
#' @return A one-row \code{data.frame} with columns \code{A} (point
#'   estimate), \code{se} (DeLong-DeLong-Clarke-Pearson standard error),
#'   \code{lower_limit} and \code{upper_limit} (confidence limits at
#'   \code{conf_level}), \code{z_value} and \code{p_value} (Wald test of
#'   \eqn{H_0\!: A = 0.5}, i.e., stochastic equality), \code{n_1} and
#'   \code{n_2} (group sample sizes), and \code{ci_method}.
#'
#' @details
#' \strong{Definition.} For independent samples \eqn{X_1, \ldots, X_{n_1}} and
#' \eqn{Y_1, \ldots, Y_{n_2}},
#' \deqn{A = \Pr(X > Y) + \tfrac{1}{2}\,\Pr(X = Y).}
#' Values of \eqn{A = 0.5} indicate stochastic equality; \eqn{A > 0.5}
#' indicates that group 1 tends to score higher. Qualitative magnitude
#' labels for \emph{A} are not reported here, in keeping with the DMAR
#' convention of reporting effect sizes as numbers with confidence
#' intervals.
#'
#' \strong{Sample estimate.} Equivalent rank-based computation:
#' \deqn{\hat A = \frac{\bar R_X - (n_1 + 1)/2}{n_2},}
#' where \eqn{\bar R_X} is the mean rank of the first sample in the pooled
#' ranking with mid-ranks for ties (Vargha & Delaney, 2000, p. 109).
#' Equivalently, \eqn{\hat A = U / (n_1 n_2)}, with \eqn{U} the Mann-Whitney
#' \emph{U}-statistic counting \eqn{X_i > Y_j} (tied pairs at \eqn{1/2}).
#'
#' \strong{Standard error.} The function uses the DeLong-DeLong-Clarke-Pearson
#' (1988) U-statistic variance estimator, which is unbiased under sampling
#' from any joint distribution (no parametric or homoscedasticity assumption).
#' Defining the placement components
#' \deqn{V_{10}(X_i) = \frac{1}{n_2}\sum_{j} \psi(X_i, Y_j), \qquad
#'         V_{01}(Y_j) = \frac{1}{n_1}\sum_{i} \psi(X_i, Y_j),}
#' with \eqn{\psi(x, y) = 1, \tfrac{1}{2}, 0} as \eqn{x > y, =, <}, the
#' variance estimate is
#' \deqn{\widehat{\mathrm{Var}}(\hat A) = \frac{S^2_{10}}{n_1} +
#'                                        \frac{S^2_{01}}{n_2},}
#' where \eqn{S^2_{10}} and \eqn{S^2_{01}} are the sample variances of the
#' \eqn{V_{10}} and \eqn{V_{01}} placement components. This is identical to
#' the (single-curve) DeLong AUC variance and is the standard nonparametric
#' variance for the Mann-Whitney functional (Brunner & Munzel, 2000).
#'
#' \strong{Confidence interval.} \code{ci_method = "logit"} (the default)
#' constructs a Wald interval on \eqn{\mathrm{logit}(A) = \log\{A/(1-A)\}}
#' using the delta method standard error
#' \eqn{\widehat{\mathrm{SE}}(\hat A)/\{\hat A(1 - \hat A)\}} and back-
#' transforms with the inverse logit. Newcombe (2006a, 2006b) showed in
#' extensive coverage simulations that logit-Wald has notably better small-
#' sample coverage than untransformed Wald, while remaining simple and free
#' of iteration. \code{ci_method = "wald"} returns the untransformed Wald
#' interval, clipped to \eqn{[0, 1]}.
#'
#' @references
#' Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher problem:
#'   Asymptotic theory and a small-sample approximation. \emph{Biometrical
#'   Journal, 42}(1), 17--25.
#'   \doi{10.1002/(SICI)1521-4036(200001)42:1<17::AID-BIMJ17>3.0.CO;2-U}
#'
#' DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing the
#'   areas under two or more correlated receiver operating characteristic
#'   curves: A nonparametric approach. \emph{Biometrics, 44}(3), 837--845.
#'
#' Hanley, J. A., & McNeil, B. J. (1982). The meaning and use of the area
#'   under a receiver operating characteristic (ROC) curve. \emph{Radiology,
#'   143}(1), 29--36.
#'
#' McGraw, K. O., & Wong, S. P. (1992). A common language effect size
#'   statistic. \emph{Psychological Bulletin, 111}(2), 361--365.
#'   \doi{10.1037/0033-2909.111.2.361}
#'
#' Newcombe, R. G. (2006a). Confidence intervals for an effect size measure
#'   based on the Mann-Whitney statistic. Part 1: General issues and
#'   tail-area-based methods. \emph{Statistics in Medicine, 25}(4), 543--557.
#'   \doi{10.1002/sim.2323}
#'
#' Newcombe, R. G. (2006b). Confidence intervals for an effect size measure
#'   based on the Mann-Whitney statistic. Part 2: Asymptotic methods and
#'   evaluation. \emph{Statistics in Medicine, 25}(4), 559--573.
#'   \doi{10.1002/sim.2324}
#'
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL
#'   common language effect size statistics of McGraw and Wong.
#'   \emph{Journal of Educational and Behavioral Statistics, 25}(2),
#'   101--132. \doi{10.3102/10769986025002101}
#'
#' @examples
#' # Two numeric vectors.
#' set.seed(113)
#' x <- rnorm(40, mean = 0.6)
#' y <- rnorm(40, mean = 0)
#' vargha_delaney_A(x, y)
#'
#' # Formula interface on the pygmalion field experiment. The first factor
#' # level (Control) forms group 1, so A below 0.5 says a randomly drawn
#' # control child tends to score below a child from the bloomer group.
#' vargha_delaney_A(iq_8 ~ treatment, data = pygmalion)
#'
#' # Wald (untransformed) interval rather than logit.
#' vargha_delaney_A(x, y, ci_method = "wald")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{smd}}, \code{\link{ci_smd}}, \code{\link{cohen_kappa}}
#'
#' @concept probability of superiority
#' @concept Vargha-Delaney A
#'
#' @keywords htest nonparametric
#'
#' @export
vargha_delaney_A <- function(x, y = NULL, data = NULL,
                             conf_level = 0.95,
                             ci_method  = c("logit", "wald")) {

  ci_method <- match.arg(ci_method)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  if (inherits(x, "formula")) {
    if (is.null(data)) {
      stop("'data' must be supplied when 'x' is a formula.", call. = FALSE)
    }
    mf <- stats::model.frame(x, data = data, na.action = stats::na.omit)
    if (ncol(mf) != 2L) {
      stop("Formula must have one outcome and one grouping variable.",
           call. = FALSE)
    }
    y_vals <- mf[[1L]]
    g_vals <- as.factor(mf[[2L]])
    lvl    <- levels(g_vals)
    if (length(lvl) != 2L) {
      stop("Grouping variable must have exactly two levels.", call. = FALSE)
    }
    x_vec <- as.numeric(y_vals[g_vals == lvl[1L]])
    y_vec <- as.numeric(y_vals[g_vals == lvl[2L]])
  } else {
    if (is.null(y)) stop("Supply 'y' (or a formula in 'x' with 'data').",
                         call. = FALSE)
    x_vec <- as.numeric(x); x_vec <- x_vec[!is.na(x_vec)]
    y_vec <- as.numeric(y); y_vec <- y_vec[!is.na(y_vec)]
  }

  n_1 <- length(x_vec)
  n_2 <- length(y_vec)
  if (n_1 < 2L || n_2 < 2L) {
    stop("Each sample must have at least 2 non-missing observations.",
         call. = FALSE)
  }

  # --- Point estimate via mean ranks (handles ties via mid-ranks).
  combined       <- c(x_vec, y_vec)
  ranks_combined <- rank(combined, ties.method = "average")
  R_bar_x        <- mean(ranks_combined[seq_len(n_1)])
  A              <- (R_bar_x - (n_1 + 1) / 2) / n_2

  # --- DeLong-DeLong-Clarke-Pearson (1988) variance via placement components.
  psi <- outer(x_vec, y_vec,
               function(a, b) (a > b) + 0.5 * (a == b))
  V_10  <- rowMeans(psi)   # length n_1: P(X_i > Y) + 0.5 P(X_i = Y) | X_i
  V_01  <- colMeans(psi)   # length n_2: P(X > Y_j) + 0.5 P(X = Y_j) | Y_j
  S2_10 <- stats::var(V_10)
  S2_01 <- stats::var(V_01)
  var_A <- S2_10 / n_1 + S2_01 / n_2
  se_A  <- sqrt(max(0, var_A))

  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (ci_method == "logit" && A > 0 && A < 1 && se_A > 0) {
    logit_A    <- log(A / (1 - A))
    se_logit_A <- se_A / (A * (1 - A))
    lo_logit   <- logit_A - z_crit * se_logit_A
    hi_logit   <- logit_A + z_crit * se_logit_A
    lower_limit <- 1 / (1 + exp(-lo_logit))
    upper_limit <- 1 / (1 + exp(-hi_logit))
  } else {
    lower_limit <- max(0, A - z_crit * se_A)
    upper_limit <- min(1, A + z_crit * se_A)
  }

  if (se_A > 0) {
    z_value <- (A - 0.5) / se_A
    p_value <- 2 * stats::pnorm(-abs(z_value))
  } else {
    z_value <- NA_real_
    p_value <- NA_real_
  }

  data.frame(
    A           = A,
    se          = se_A,
    lower_limit = lower_limit,
    upper_limit = upper_limit,
    z_value     = z_value,
    p_value     = p_value,
    n_1         = n_1,
    n_2         = n_2,
    ci_method   = ci_method,
    stringsAsFactors = FALSE,
    row.names   = NULL
  )
}
