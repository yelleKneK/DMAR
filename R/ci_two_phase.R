#' Confidence Interval for a Mean Under Two-Phase Sampling With a Proxy Measurement
#'
#' @description
#' The confidence interval for the population mean of a target measurement
#' that is costly to obtain and so is observed on only a random subsample of
#' units, when an inexpensive proxy measurement is observed on every unit of
#' a larger pool. The design is two-phase, or double, sampling (Neyman,
#' 1938): the first phase draws the pool and measures the proxy on all of
#' it, and the second draws a simple random subsample of the pool, the
#' labeled units, and measures the target on them. The regression estimator
#' of the target mean, the row \code{mean_Y} of the table,
#' \deqn{\hat{\mu}_Y = \bar y_n + \hat b\,(\bar x_M - \bar x_n),}
#' adjusts the labeled units' target mean \eqn{\bar y_n} by how far their
#' proxy mean \eqn{\bar x_n} sits from the proxy mean of the whole pool
#' \eqn{\bar x_M}, through the least squares slope \eqn{\hat b} of the
#' target on the proxy in the labeled units. When the two measurements are
#' strongly correlated, the interval is much narrower than the one the
#' labeled units give alone, which is the point of the design.
#'
#' @param y The target measurement on the \eqn{n} labeled units (numeric,
#'   or logical for a binary measurement, taken as 0 and 1).
#' @param x The proxy measurement on the same \eqn{n} units, aligned with
#'   \code{y}.
#' @param x_unlabeled The proxy measurement on the units of the pool that
#'   have no target value, \eqn{N = M - n} of them, where \eqn{M} is the
#'   size of the pool.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param method \code{"regression"} (default), the regression estimator
#'   above; or \code{"prediction_powered"}, the prediction-powered estimate
#'   of Angelopoulos et al. (2023), which gives the proxy the fixed weight
#'   \eqn{\lambda = 1}, the row \code{proxy_weight} (see Details).
#' @param cluster,cluster_unlabeled Optional cluster identifiers for the
#'   labeled and the unlabeled units (for example the respondent who
#'   contributed several units, or the prompt that generated several
#'   items), supplied together or not at all. When given, the standard error
#'   is computed from the cluster totals of each unit's contribution to the
#'   estimate, and the limits use the \emph{t} distribution on
#'   \eqn{Q_n - 1} degrees of freedom, where \eqn{Q_n} is the number of
#'   clusters among the labeled units (see Details).
#'
#' @details
#' \strong{The estimand.} The interval is for \eqn{\mu_Y}, the mean of the
#' target measurement in the population from which the pool was drawn at
#' random, the labeled units being a simple random subsample of the pool.
#' The proxy enters only as a covariate. It need not be on the target's
#' scale, unbiased for it, or even measured in its units, and whatever its
#' bias, that bias does not enter the estimate; what the proxy contributes is
#' its correlation \eqn{\rho} with the target. In the evaluation of a system
#' whose outputs an automated judge labels, for example, the target is a
#' human assessor's label of an output, the proxy is the judge's label of the
#' same output, and the pool is the set of outputs the judge labels.
#'
#' \strong{Standard error.} To first order,
#' \deqn{\mathrm{Var}(\hat{\mu}_Y) = \frac{\sigma_Y^2 (1 - \rho^2)}{n} +
#'   \frac{\rho^2 \sigma_Y^2}{M},}
#' the first term the residual variation the regression leaves among the
#' \eqn{n} labeled units and the second the sampling variation of the pool
#' itself. The reported standard error estimates
#' \eqn{\sigma_e^2/n + \beta^2\sigma_X^2/M}, with the residual mean square
#' of the target on the proxy in the labeled units (on \eqn{n - 2} degrees
#' of freedom) for \eqn{\sigma_e^2} and \eqn{\hat b^2} times the proxy's
#' variance over the whole pool for \eqn{\beta^2 \sigma_X^2}, and the limits
#' use the \emph{t} distribution on the \eqn{n - 2} degrees of freedom of
#' the first term, which dominates when the pool is large. If the pool is
#' itself the whole population of interest, a fixed collection of units
#' every one of which carries the proxy, the estimand is the pool's own
#' target mean and the reported standard error is conservative for it: the
#' second term is then not needed, and the first can be multiplied by
#' \eqn{1 - n/M}.
#'
#' \strong{Relation to prediction-powered inference.} Writing
#' \eqn{\bar x_N} for the proxy mean of the unlabeled units,
#' \eqn{\hat{\mu}_Y = \lambda \bar x_N + (\bar y_n - \lambda \bar x_n)} with
#' \eqn{\lambda = \hat b N / M}. With \eqn{\lambda = 1} in its place this is
#' the prediction-powered estimate of Angelopoulos et al. (2023): the proxy
#' mean of the unlabeled units, corrected by the mean target-minus-proxy
#' difference of the labeled units. That estimate is unbiased for any
#' proxy, but it uses the proxy efficiently only when the proxy is on the
#' target's scale with a slope near one, while the regression estimator
#' learns the weight from the data; the \code{proxy_weight} row reports
#' \eqn{\lambda}, the weight each method applied. \code{method = "prediction_powered"} reports
#' the interval with standard error
#' \eqn{\{s_{Y - X}^2/n + s_{X,N}^2/N\}^{1/2}} and the \emph{t}
#' distribution on \eqn{n - 1} degrees of freedom; Angelopoulos et al.
#' state it with the normal quantile, which the \emph{t} interval approaches
#' as \eqn{n} grows.
#'
#' \strong{Clusters.} With \code{cluster} and \code{cluster_unlabeled}, each
#' unit's contribution to the linearized estimate is summed within its
#' cluster, and the standard error is
#' \eqn{\{Q_n/(Q_n - 1) \sum_q T_q^2\}^{1/2}}, the sum running over every
#' cluster of the pool, \eqn{T_q} the total of cluster \eqn{q}, and
#' \eqn{Q_n} the number of clusters among the labeled units. For the
#' regression estimate a labeled unit contributes
#' \eqn{e_i/n + \hat b(x_i - \bar x_M)/M} and an unlabeled unit
#' \eqn{\hat b(x_j - \bar x_M)/M}, where \eqn{e_i} is its residual; for
#' the prediction-powered estimate a labeled unit contributes its
#' target-minus-proxy difference less the mean difference, over \eqn{n},
#' and an unlabeled unit its proxy less \eqn{\bar x_N}, over \eqn{N}. The
#' term that dominates the variance, the residual term, comes from the
#' labeled clusters alone, so they, and not all the clusters of the pool,
#' set both the degrees of freedom, \eqn{Q_n - 1}, and the small-sample
#' factor \eqn{Q_n/(Q_n - 1)}. When the labeled units come from a few
#' clusters this matters: the interval then rests on those few, however
#' many clusters the pool has. With one unit per cluster the labeled part
#' of the prediction-powered standard error reduces to
#' \eqn{s_{Y - X}^2/n}, its form without clusters. At least two labeled
#' clusters are required.
#'
#' \strong{A constant proxy or target.} The regression estimate needs both
#' measurements to vary among the labeled units. A proxy that takes one
#' value there leaves the slope undefined, and a target that takes one value
#' (a binary target that is 1 on every labeled unit, for example) makes the
#' residual variation and the slope both zero, so the interval would have
#' zero width. Each is refused with an error rather than reported as an
#' interval; more labeled units are the remedy. The prediction-powered
#' estimate does not rest on the regression and is defined in both cases.
#'
#' \strong{Missing values.} A labeled unit enters only when its target,
#' its proxy, and (if given) its cluster are all present, and an unlabeled
#' unit only when its proxy (and cluster) is present; the others are set
#' aside, which leaves the design intact when values are missing completely
#' at random. The \code{n_labeled} and \code{n_pool} rows count the units
#' used.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}), its \code{method}
#'   attribute naming the method, with rows
#'   \item{mean_Y}{the estimate of \eqn{\mu_Y}, the population mean of the
#'     target.}
#'   \item{se, df}{its standard error and the degrees of freedom of the
#'     \emph{t} reference distribution.}
#'   \item{ci_lower, ci_upper}{the confidence limits.}
#'   \item{slope}{the least squares slope \eqn{\hat b} of the target on the
#'     proxy in the labeled units. The regression estimate applies it; the
#'     prediction-powered estimate does not, and the row is reported for
#'     comparison (\code{NA} when the proxy is constant on the labeled
#'     units).}
#'   \item{proxy_weight}{the weight \eqn{\lambda} the estimate gives the
#'     proxy mean of the unlabeled units: \eqn{\hat b N / M} for the
#'     regression estimate and 1 for the prediction-powered estimate.}
#'   \item{correlation}{the correlation of target and proxy in the labeled
#'     units.}
#'   \item{proxy_bias}{the proxy mean minus the target mean in the labeled
#'     units, the proxy's bias when the two share a scale. It does not enter
#'     the regression estimate.}
#'   \item{se_ratio_target_only}{the standard error divided by that of the
#'     labeled units' target mean alone (\eqn{s_Y/\sqrt n}, or its cluster
#'     counterpart), so values below 1 are the precision the proxy buys.}
#'   \item{n_labeled, n_pool}{the numbers of labeled units and of units in
#'     the pool, \eqn{n} and \eqn{M}.}
#'
#' @references
#' Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., & Zrnic, T.
#'   (2023). Prediction-powered inference. \emph{Science, 382}(6671),
#'   669--674. \doi{10.1126/science.adi6000}
#'
#' Neyman, J. (1938). Contribution to the theory of sampling human
#'   populations. \emph{Journal of the American Statistical Association,
#'   33}(201), 101--116. \doi{10.1080/01621459.1938.10503378}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_two_phase}} and
#'   \code{\link{ss_aipe_two_phase_cost}} to plan the numbers of labeled and
#'   pool units, \code{\link{ss_aipe_two_phase_sensitivity}} to evaluate a
#'   plan, and \code{\link{ss_seq_two_phase}} for the sequential rule that
#'   needs no planning values.
#'
#' @family confidence intervals
#'
#' @keywords htest design
#'
#' @examples
#' # A brief screening scale (the proxy) is scored for all 2,000
#' # respondents of a survey; a structured interview (the target) is given
#' # to a random 150 of them.
#' set.seed(113)
#' screen    <- rnorm(2000, mean = 50, sd = 10)
#' interview <- 20 + 0.6 * screen + rnorm(2000, sd = 5)
#' labeled   <- sample(2000, 150)
#' ci_two_phase(y = interview[labeled], x = screen[labeled],
#'              x_unlabeled = screen[-labeled])
#'
#' # The screener is not on the interview's scale: the interview rises
#' # about 0.6 points per screener point. The prediction-powered estimate
#' # gives the screener a weight of 1 regardless, so it stays unbiased but
#' # its interval is wider than the regression estimator's.
#' ci_two_phase(y = interview[labeled], x = screen[labeled],
#'              x_unlabeled = screen[-labeled],
#'              method = "prediction_powered")
#'
#' @export
ci_two_phase <- function(y, x, x_unlabeled, conf_level = 0.95,
                         method = c("regression", "prediction_powered"),
                         cluster = NULL, cluster_unlabeled = NULL) {
  .dmar_check_conf_level(conf_level)
  method <- match.arg(method)
  .two_phase_check_planning(conf_level = conf_level)
  d <- .two_phase_prepare(y, x, x_unlabeled, cluster, cluster_unlabeled,
                          min_unlabeled = if (method == "regression") 1L else 2L)
  y <- d$y; x <- d$x; xu <- d$x_unlabeled
  n <- length(y); N <- length(xu); M <- n + N
  ybar <- mean(y); xbar_n <- mean(x); xbar_N <- mean(xu)
  xbar_M <- (n * xbar_n + N * xbar_N) / M

  sd_x <- stats::sd(x); sd_y <- stats::sd(y)
  slope <- if (sd_x > 0) stats::cov(y, x) / sd_x^2 else NA_real_
  correlation <- if (sd_x > 0 && sd_y > 0) stats::cor(y, x) else NA_real_
  clustered <- !is.null(d$cluster)

  if (method == "regression") {
    if (is.na(slope))
      .dmar_refuse_argument(c(
        "The proxy takes a single value on the labeled units, so its ",
        "relation to the target cannot be estimated and the regression ",
        "estimator is not defined."), arg = NULL)
    if (all(y == y[1L]))
      .dmar_refuse_argument(c(
        "The target takes a single value on the labeled units (", y[1L],
        "), so the residual variation and the slope are both zero and ",
        "the regression interval would have zero width. Label more ",
        "units; a target that rarely departs from one value needs more ",
        "labeled units before its variation shows."), arg = NULL)
    fast <- .ci_two_phase_fast(y, x, N, xbar_N, sum((xu - xbar_N)^2),
                               conf_level)
    est <- fast[["mean_Y"]]
    se  <- fast[["se"]]
    df  <- fast[["df"]]
    lam <- slope * N / M
    # Each unit's contribution to the linearized estimate: the residual
    # share of a labeled unit plus every unit's share of the pool mean.
    infl_lab <- (y - ybar - slope * (x - xbar_n)) / n + slope * (x - xbar_M) / M
    infl_unl <- slope * (xu - xbar_M) / M
  } else {
    lam <- 1
    est <- xbar_N + (ybar - xbar_n)
    se  <- sqrt(stats::var(y - x) / n + stats::var(xu) / N)
    df  <- n - 1
    infl_lab <- (y - x - (ybar - xbar_n)) / n
    infl_unl <- (xu - xbar_N) / N
  }

  if (clustered) {
    totals <- tapply(c(infl_lab, infl_unl), c(d$cluster, d$cluster_unlabeled),
                     sum)
    totals_target <- tapply((y - ybar) / n, d$cluster, sum)
    Q_lab <- length(totals_target)
    if (Q_lab < 2L)
      .dmar_refuse_argument(c(
        "The labeled units span fewer than two clusters, so a cluster ",
        "standard error cannot be formed."), arg = NULL)
    # The labeled contributions sum to zero over the labeled units (the
    # residuals, or the centered differences), so their cluster totals sum
    # to zero over the Q_lab labeled clusters, and the sum of squared totals
    # understates their variance by about (Q_lab - 1) / Q_lab. That part
    # dominates the variance, so Q_lab sets both the factor and the degrees
    # of freedom; clusters that hold only unlabeled units carry the small
    # pool term alone and add no information about the residual variation.
    se <- sqrt(Q_lab / (Q_lab - 1) * sum(totals^2))
    df <- Q_lab - 1
    se_target <- sqrt(Q_lab / (Q_lab - 1) * sum(totals_target^2))
  } else {
    se_target <- sd_y / sqrt(n)
  }

  q <- stats::qt(1 - (1 - conf_level) / 2, df)
  # The estimate is named for its estimand, the mean of the target Y, and
  # the weight the estimate gives the proxy for what it is (Ken's ruling of
  # 2026-11-11).
  out <- data.frame(
    term  = c("mean_Y", "se", "df", "ci_lower", "ci_upper",
              "slope", "proxy_weight", "correlation", "proxy_bias",
              "se_ratio_target_only", "n_labeled", "n_pool"),
    value = c(est, se, df, est - q * se, est + q * se,
              slope, lam, correlation, xbar_n - ybar,
              if (se_target > 0) se / se_target else NA_real_, n, M),
    stringsAsFactors = FALSE, row.names = NULL
  )
  attr(out, "method") <- method
  .as_dmar_tbl(out, conf_level = conf_level)
}
