#' Sample Size Planning From the AIPE Perspective for the Least Costly Two-Phase Design
#'
#' @description
#' Takes planning values for a two-phase design (the standard deviation of a
#' costly target measurement, its correlation with an inexpensive proxy, and
#' the cost of a proxy relative to a target measurement) and a desired
#' confidence interval width for the target mean. Returns the number of
#' labeled units \eqn{n}, which carry both measurements, and the number of
#' pool units \eqn{M}, which carry the proxy, of the least costly design in
#' whole numbers whose interval from \code{\link{ci_two_phase}} has an
#' expected width no larger than desired. The table also gives the cost
#' relative to measuring the target alone and the correlation the proxy
#' must have before it can save anything.
#'
#' @param true_sd_Y Planning value of the standard deviation of the target
#'   measurement, \eqn{\sigma_Y}.
#' @param true_cor Planning value of the correlation between the target and the
#'   proxy, \eqn{\rho}. Only its absolute value matters.
#' @param width The desired full width of the confidence interval (the
#'   upper limit minus the lower limit), in the units of the target.
#' @param cost_ratio The cost of one proxy measurement relative to one
#'   target measurement, \eqn{r = c_P / c_T}, a number in (0, 1).
#' @param conf_level Confidence level (default \code{0.95}).
#'
#' @details
#' \strong{Estimand and interval.} The quantity to be estimated is
#' \eqn{\mu_Y}, the population mean of the target measurement. The sampling
#' unit is one unit of the pool: \eqn{M} units are drawn at random, the
#' proxy is observed on all of them, and the target on a simple random
#' subsample of \eqn{n} of them. The interval planned for is the one
#' \code{\link{ci_two_phase}} reports for the regression estimator of
#' double sampling, a Wald-type interval whose standard error estimates
#' \eqn{\sigma_Y^2(1 - \rho^2)/n + \rho^2\sigma_Y^2/M} and whose reference
#' distribution is \emph{t} on \eqn{n - 2} degrees of freedom. This is a
#' plan for the expected width, with no assurance: about half of all
#' realized intervals are wider than \code{width}.
#' \code{\link{ss_aipe_two_phase}} adds an assurance for a given pool.
#'
#' \strong{The cost model.} Every unit of the pool, the labeled ones
#' included, costs \eqn{c_P} for its proxy, and each labeled unit costs
#' \eqn{c_T} more for its target, so a design costs
#' \eqn{c_T n + c_P M}, which is \eqn{n + rM} in units of one target
#' measurement. In an evaluation of a system whose outputs an automated
#' judge labels, for example, \eqn{c_T} is the cost of a human assessor's
#' label and \eqn{c_P} that of the judge's. The target alone on \eqn{n}
#' units costs \eqn{n}.
#'
#' \strong{Which designs are admissible.} A two-phase design has
#' \eqn{n \ge 3} and \eqn{M \ge n + 1}, the least
#' \code{\link{ci_two_phase}} accepts, and meets the width when, with
#' \eqn{q = 1 - (1 - \mathrm{conf\_level})/2} and \eqn{\omega} the desired width
#' (\code{width}),
#' \deqn{2\, t_{q,\, n-2}\, \sigma_Y \left\{\frac{1 - \rho^2}{n} +
#'   \frac{\rho^2}{M}\right\}^{1/2} \le \omega,}
#' the criterion \code{\link{ss_aipe_two_phase}} uses for a given pool, so
#' that planner, given \code{n_pool = necessary_n_pool}, returns
#' \code{necessary_n_labeled}. The left side sets the two estimates in the
#' standard error at their planning values, which overstates the expected
#' width slightly (see \code{\link{ss_aipe_two_phase}}). The alternative
#' without a proxy is the \emph{t} interval for one mean on the smallest
#' \eqn{n} with
#' \eqn{2\, t_{q,\, n-1}\, c_4(n)\, \sigma_Y / \sqrt{n} \le \omega}, where
#' \eqn{c_4(n)\sigma_Y} is the expected sample standard deviation under
#' normality; this is the \code{necessary_n_target_only} row here and in
#' \code{\link{ss_aipe_two_phase}}. Because the \emph{t} quantile is part of
#' the requirement, a design with very few labeled units is admissible only
#' when it is wide enough to pay for its quantile (12.7 at a 95\% level
#' with 3 labeled units), and no separate minimum is imposed.
#'
#' \strong{The search.} For each \eqn{n} the smallest admissible pool
#' follows from the inequality above, so the search is over \eqn{n} alone.
#' Since \eqn{t_{q,\, n-2}} is at least the normal quantile \eqn{z_q}, no
#' design with \eqn{n} labeled units can cost less than
#' \eqn{n + r \max\{n + 1,\; b/(1 - a/n)\}} for \eqn{n > a}, and none
#' exists for \eqn{n \le a}, where \eqn{a = (1 - \rho^2) s} and
#' \eqn{b = \rho^2 s} with \eqn{s = (2 z_q \sigma_Y/\omega)^2}. That bound
#' is convex in \eqn{n}, so the function starts at its minimizer with the
#' target alone as the design to beat, examines every \eqn{n} above and
#' then below the minimizer for as long as the bound is less than the cost
#' of the best design found, and returns that design. No design in whole
#' numbers costs less: every \eqn{n} not examined has a bound at least as
#' large as the cost returned. A two-phase design is returned only when it
#' costs less than the target alone, and among designs of equal cost the
#' first found is kept, so a tie with the target alone returns the target
#' alone. The search and its bound were derived for DMAR; in 2,200
#' planning sets drawn at random (correlations of either sign from .02 to
#' .999, cost ratios from .001 to .98, four confidence levels) the result
#' equaled an exhaustive search over pairs \eqn{(n, M)} every time.
#'
#' \strong{The fractional solution, for comparison.} If the sizes could be
#' fractional and the normal quantile were used, minimizing
#' \eqn{n + rM} subject to \eqn{a/n + b/M = 1} (a Lagrange multiplier
#' argument) gives \eqn{n = a + \sqrt{rab}} and
#' \eqn{M = b + \sqrt{ab/r}}, at a cost of
#' \eqn{(\sqrt{1 - \rho^2} + |\rho|\sqrt{r})^2} times the \eqn{s} units the
#' target alone would need. That factor is the
#' \code{relative_cost_continuous} row. It is below 1 exactly when
#' \eqn{|\rho| > 2\sqrt{r}/(1 + r)}, the \code{threshold_cor} row; at or
#' below that correlation the row is 1. The two sizes are in the ratio
#' \eqn{M/n = \sqrt{b/(ra)}}: each is proportional to its share of the
#' variance and inversely proportional to the square root of its cost.
#' Neyman (1938, Equations 49, 52, and 53, pp. 110--111) derived that
#' square root rule for double sampling in which the inexpensive variable
#' forms strata and the variance is minimized for a fixed total cost; the
#' form here, with a regression on the proxy and the cost minimized for a
#' fixed width, is the same allocation problem stated the other way
#' round. For the regression estimator itself, Tamhane (1978, Equation
#' 2.2, p. 421) gives the same ratio of the two sizes for a fixed budget,
#' and the correlation above which double sampling is no worse than the
#' target alone, \eqn{\rho^2 \ge 4 c_1 c_2/(c_1 + c_2)^2} with
#' \eqn{c_1} and \eqn{c_2} the costs of the proxy and of the target,
#' which is \code{threshold_cor}; the result is for a bivariate normal
#' population with a known covariance matrix. The fractional solution is
#' a lower bound: no design in whole numbers that meets the
#' width costs less than \code{relative_cost_continuous} times \eqn{s}
#' target measurements. The design returned costs more than that bound
#' because sizes are whole numbers and because the interval uses the
#' \emph{t} quantile (a median of 7\% more, and far more when only a few
#' units are needed, in the 891 random planning sets above in which a
#' two-phase design was returned).
#'
#' \strong{What can happen at the edges.}
#' \itemize{
#'   \item A correlation at or below \code{threshold_cor}: the target alone
#'     is returned (\code{necessary_n_pool} is 0 and \code{relative_cost}
#'     is 1). In the random planning sets above no two-phase design was
#'     ever the least costly there.
#'   \item A correlation above \code{threshold_cor}: the threshold is exact
#'     for fractional sizes and the normal quantile only. The target alone
#'     can still be the least costly design in whole numbers, most often
#'     when few units are needed or the correlation is a little above the
#'     threshold (376 of the 1,267 random planning sets above the
#'     threshold), and it is then returned.
#'   \item \code{true_cor = 0}: the proxy carries no information and the target
#'     alone is returned.
#'   \item A proxy that costs as much as the target or more
#'     (\code{cost_ratio} of 1 or more) never saves anything and is refused
#'     with an error of class \code{dmar_cost_ratio_invalid}.
#'   \item A very wide \code{width}: the target alone needs at least 2
#'     units, the fewest from which a standard deviation can be estimated.
#'   \item Sizes that cannot be represented as whole numbers (above
#'     \eqn{2^{52}}, about 4.5e15), or a search that would examine more
#'     than \eqn{10^8} values of \eqn{n}: refused with an error of class
#'     \code{dmar_two_phase_no_design}. Neither arises for a width within
#'     a factor of a million of \code{true_sd_Y}.
#' }
#'
#' @return A \code{data.frame} (classes \code{dmar_ss_aipe} and
#'   \code{dmar_tbl}) with rows
#'   \item{necessary_n_labeled}{the number of labeled units, \eqn{n}; when
#'     \code{necessary_n_pool} is 0, the number of units for the target
#'     alone.}
#'   \item{necessary_n_pool}{the number of units in the pool, \eqn{M},
#'     every one of which carries the proxy; 0 when the design measures
#'     the target alone.}
#'   \item{relative_cost}{the cost of the returned design, \eqn{n + rM},
#'     over the cost of the target alone, \code{necessary_n_target_only},
#'     both in whole numbers under the criteria above; 1 when the returned
#'     design is the target alone.}
#'   \item{relative_cost_continuous}{the same ratio for fractional sizes
#'     under the normal quantile,
#'     \eqn{(\sqrt{1 - \rho^2} + |\rho|\sqrt{r})^2}, or 1 at or below
#'     \code{threshold_cor}; reported for comparison.}
#'   \item{expected_width}{the planned full width at the returned sizes
#'     and the planning values: the left side of the two-phase inequality
#'     above, or of the target-alone inequality when the pool is 0. It is
#'     at most \code{width}.}
#'   \item{necessary_n_target_only}{the number of units the target alone
#'     needs for the same width.}
#'   \item{threshold_cor}{the absolute correlation above which the proxy
#'     lowers the cost of the fractional solution,
#'     \eqn{2\sqrt{r}/(1 + r)}.}
#'   \item{width, true_sd_Y, true_cor, cost_ratio, conf_level}{the inputs.}
#'
#' @references
#' Neyman, J. (1938). Contribution to the theory of sampling human
#'   populations. \emph{Journal of the American Statistical Association,
#'   33}(201), 101--116. \doi{10.1080/01621459.1938.10503378}
#'
#' Tamhane, A. C. (1978). Inference based on regression estimator in
#'   double sampling. \emph{Biometrika, 65}(2), 419--427.
#'   \doi{10.1093/biomet/65.2.419}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_two_phase}} for the labeled units given a
#'   pool, with or without an assurance, \code{\link{ci_two_phase}} for the
#'   interval, and \code{\link{ss_aipe_two_phase_sensitivity}} to evaluate
#'   a design by simulation.
#'
#' @family AIPE sample size planning
#'
#' @keywords design
#'
#' @examples
#' # A proxy measurement costs 1% of a target measurement. Planning values
#' # of 0.43 for the target's standard deviation and 0.55 for its
#' # correlation with the proxy; a 90% interval no wider than 0.10. The
#' # least costly design labels 151 units of a pool of 957 and costs 79%
#' # of what the 202 units of the target alone would.
#' ss_aipe_two_phase_cost(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
#'                        cost_ratio = 0.01, conf_level = 0.90)
#'
#' # At 10% of the cost the same proxy no longer pays: a correlation of 0.55
#' # is below the threshold, and the design measures the target alone.
#' ss_aipe_two_phase_cost(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
#'                        cost_ratio = 0.10, conf_level = 0.90)
#'
#' # A worked example in standard deviation units: a proxy that correlates
#' # 0.90 with the target at 5% of its cost, and a 95% interval of width
#' # 0.20. The design is 109 labeled units in a pool of 1011, a cost of
#' # 109 + 0.05 * 1011 = 159.55 target measurements against 387.
#' plan <- ss_aipe_two_phase_cost(true_sd_Y = 1, true_cor = 0.90, width = 0.20,
#'                                cost_ratio = 0.05)
#' plan
#'
#' # The companion planner, given that pool, returns the same number of
#' # labeled units.
#' ss_aipe_two_phase(true_sd_Y = 1, true_cor = 0.90, width = 0.20, n_pool = 1011)
#'
#' @export
ss_aipe_two_phase_cost <- function(true_sd_Y, true_cor, width, cost_ratio,
                                   conf_level = 0.95) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  .two_phase_check_planning(sd_y = true_sd_Y, rho = true_cor, width = width,
                            conf_level = conf_level)
  if (!is.numeric(cost_ratio) || length(cost_ratio) != 1L ||
      is.na(cost_ratio) || cost_ratio <= 0)
    .dmar_refuse(c(
      "'cost_ratio' must be a single positive number below 1: the cost of ",
      "a proxy measurement relative to a target measurement. A ratio of 0 ",
      "or less is not a cost; give the proxy's cost divided by the ",
      "target's, for example cost_ratio = 0.05."),
      "cost_ratio_invalid", arg = "cost_ratio", value = cost_ratio)
  if (cost_ratio >= 1)
    .dmar_refuse(c(
      "'cost_ratio' must be below 1; ", format(cost_ratio), " was given. ",
      "A proxy measurement that costs as much as a target measurement or ",
      "more never lowers the cost of a design; measuring the target alone ",
      "is then the least costly, and ss_aipe_two_phase() gives that ",
      "sample size in its necessary_n_target_only row."),
      "cost_ratio_invalid", arg = "cost_ratio", value = cost_ratio)

  q <- 1 - (1 - conf_level) / 2
  z <- stats::qnorm(q)
  r <- cost_ratio
  # sigma_Y^2 / V, the number of units the target alone needs under the
  # normal quantile before rounding, formed as one ratio so that a very
  # small or very large 'true_sd_Y' with a width on its scale neither underflows
  # nor overflows. a and b are its two parts: a design (n, M) meets the
  # width under the normal quantile exactly when a / n + b / M <= 1.
  s2 <- (2 * z * true_sd_Y / width)^2
  a <- (1 - true_cor^2) * s2
  b <- true_cor^2 * s2
  threshold <- 2 * sqrt(r) / (1 + r)
  # The fractional solution (Neyman's allocation), for comparison only: its
  # cost relative to the target alone, and its pool size, which is what
  # the representable-size guard looks at.
  sd_e <- sqrt(1 - true_cor^2)
  sd_b <- abs(true_cor)
  pays_continuous <- sd_b > threshold
  relative_cost_continuous <-
    if (pays_continuous) (sd_e + sd_b * sqrt(r))^2 else 1
  M_c <- sd_b * (sd_e + sd_b * sqrt(r)) * s2 / sqrt(r)

  # Whole numbers are exact in a double up to 2^52; beyond that a size
  # cannot be rounded up or raised by one.
  no_design <- function()
    .dmar_refuse(c(
      "No design in whole numbers can be represented for these planning ",
      "values: the width is too narrow for 'true_sd_Y', or the proxy too ",
      "inexpensive, for sizes below 2^52 (about 4.5e15)."),
      "two_phase_no_design")
  if (!is.finite(s2) || s2 > 2^52 ||
      (pays_continuous && (!is.finite(M_c) || M_c > 2^52)))
    no_design()

  # The interval ci_two_phase() reports, at the planning values: full width
  # 2 t(q; n - 2) sd_y sqrt((1 - rho^2) / n + rho^2 / M). A design is
  # feasible when that is no wider than 'width'. Vectorized in n and M.
  width_at <- function(n, M)
    2 * stats::qt(q, n - 2) * true_sd_Y *
      sqrt((1 - true_cor^2) / n + true_cor^2 / M)

  # The target alone: the criterion and the search of ss_aipe_two_phase()
  # for its necessary_n_target_only row (the expected width of the t
  # interval for one mean, exact under normality through c4(n)), so the
  # two functions give one number. The search stops with a plain error
  # when it passes 2^52; that is the same condition as the guard above.
  width_target <- function(n)
    2 * stats::qt(q, n - 1) * .c4(n) * true_sd_Y / sqrt(n)
  meets_target <- function(n) width_target(n) <= width
  n_target <- tryCatch(
    .two_phase_search(meets_target, meets_target,
                      function(lo, hi) meets_target(lo), s2, 2),
    error = function(e) no_design())

  # The smallest pool for each n of a vector, or NA where no pool of a
  # representable size gives the width with n labeled units. The division
  # can be off by a unit through rounding, so the pool is raised while the
  # width is not met and lowered while one unit fewer still meets it.
  least_pool <- function(n) {
    rem <- (width / (2 * stats::qt(q, n - 2) * true_sd_Y))^2 -
      (1 - true_cor^2) / n
    M <- rep(NA_real_, length(n))
    ok <- which(rem > 0)
    if (length(ok) == 0L) return(M)
    M[ok] <- pmax(n[ok] + 1, ceiling(true_cor^2 / rem[ok]))
    M[which(M > 2^52)] <- NA_real_
    repeat {
      up <- which(!is.na(M) & width_at(n, M) > width)
      if (length(up) == 0L) break
      M[up] <- M[up] + 1
      M[which(M > 2^52)] <- NA_real_
    }
    repeat {
      down <- which(!is.na(M) & M - 1 > n)
      down <- down[width_at(n[down], M[down] - 1) <= width]
      if (length(down) == 0L) break
      M[down] <- M[down] - 1
    }
    M
  }

  # A lower bound for the cost of any design with n labeled units. The t
  # quantile is at least the normal quantile, so a feasible (n, M) has
  # a / n + b / M <= 1: none exists for n <= a, and for n > a the pool is
  # at least b / (1 - a / n), and at least n + 1. The bound
  #   lower(n) = n + r max(n + 1, b / (1 - a / n))
  # is convex in n > a (the larger of a convex and a linear function), so
  # it rises on either side of its minimizer.
  lower <- function(n) {
    out <- rep(Inf, length(n))
    ok <- which(n > a)
    out[ok] <- n[ok] + r * pmax(n[ok] + 1, b / (1 - a / n[ok]))
    out
  }
  # The minimizer of the bound: that of n + r b / (1 - a / n), which is
  # a + sqrt(r a b), unless the pool there is below n + 1; the minimizer is
  # then the n at which b / (1 - a / n) = n + 1, the larger root of
  # n^2 - (a + b - 1) n - a = 0, since the bound falls before that point
  # and rises with slope 1 + r after it.
  x <- a + sqrt(r * a * b)
  if (!isTRUE(b / (1 - a / x) >= x + 1))
    x <- ((a + b - 1) + sqrt((a + b - 1)^2 + 4 * a)) / 2
  n_c <- max(3, ceiling(x), floor(a) + 1)

  # The search. The best design so far starts as the target alone. From
  # n_c the scan goes up while the bound is below the best cost so far,
  # then down from n_c - 1 on the same rule. Because the bound is convex
  # and n_c is at or just above its minimizer, every n not scanned has a
  # bound, and so a cost, of at least the best cost found: the result is
  # the least costly design in whole numbers. A design replaces the best
  # so far only when it is cheaper by more than a relative 1e-12, so a tie
  # keeps the design found first, and a tie with the target alone keeps
  # the target alone. The n are evaluated in blocks (qt() is vectorized);
  # the block is walked in scan order, so the result is that of a scan one
  # n at a time.
  state <- new.env(parent = emptyenv())
  state$n <- n_target
  state$M <- 0
  state$cost <- n_target
  state$evaluations <- 0
  scan <- function(from, step) {
    block <- 1
    repeat {
      n <- from + step * (seq_len(block) - 1)
      n <- n[n >= 3 & n <= 2^52]
      if (length(n) == 0L) return(invisible(NULL))
      bound <- lower(n)
      stop_at <- which(bound >= state$cost)
      if (length(stop_at) > 0L) n <- n[seq_len(stop_at[1L] - 1L)]
      if (length(n) > 0L) {
        state$evaluations <- state$evaluations + length(n)
        if (state$evaluations > 1e8) no_design()
        M <- least_pool(n)
        cost <- n + r * M
        for (i in seq_along(n)) {
          if (bound[i] >= state$cost) return(invisible(NULL))
          if (!is.na(cost[i]) && cost[i] < state$cost * (1 - 1e-12)) {
            state$n <- n[i]
            state$M <- M[i]
            state$cost <- cost[i]
          }
        }
      }
      if (length(stop_at) > 0L) return(invisible(NULL))
      from <- from + step * block
      block <- min(2 * block, 65536)
    }
  }
  scan(n_c, 1)
  scan(n_c - 1, -1)

  n <- state$n
  M <- state$M
  relative_cost <- if (M > 0) (n + r * M) / n_target else 1
  expected_width <- if (M > 0) width_at(n, M) else width_target(n)

  out <- data.frame(
    term  = c("necessary_n_labeled", "necessary_n_pool", "relative_cost",
              "relative_cost_continuous", "expected_width",
              "necessary_n_target_only", "threshold_cor", "width", "true_sd_Y",
              "true_cor", "cost_ratio", "conf_level"),
    value = c(n, M, relative_cost, relative_cost_continuous, expected_width,
              n_target, threshold, width, true_sd_Y, true_cor, r, conf_level),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
