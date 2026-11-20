#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Mean Under Two-Phase Sampling
#'
#' @description
#' Evaluates a two-phase design by simulation: how wide the confidence
#' intervals from \code{\link{ci_two_phase}} turn out to be, how often they
#' are no wider than planned, and how often they miss the population mean
#' of the target, when the design is carried out in a population that may
#' differ from the one the plan assumed. The design is either planned by
#' \code{\link{ss_aipe_two_phase}} from planning values, so that the
#' distance between planning values and population values is what is
#' studied, or fixed by the user.
#'
#' @param true_sd_Y,true_cor The population standard deviation of the
#'   target measurement and its correlation with the proxy, from which the
#'   data are generated when \code{population} is \code{NULL}.
#' @param supposed_sd_Y,supposed_cor Planning values passed to
#'   \code{\link{ss_aipe_two_phase}} to choose the number of labeled units;
#'   supply both, or \code{specified_n_labeled} instead.
#' @param specified_n_labeled The number of labeled units to evaluate, in
#'   place of a planned one.
#' @param n_pool The number of units in the pool, \eqn{M}, a whole number
#'   larger than the number of labeled units.
#' @param width The desired full width of the confidence interval, the
#'   planning target and the yardstick for \code{pct_ci_less_width}.
#' @param population \code{NULL} (default) for a bivariate normal
#'   population, or a data frame of observed (target, proxy) pairs with
#'   numeric columns \code{y} and \code{x}, to be resampled with replacement
#'   as the population; \code{true_sd_Y} and \code{true_cor} are then
#'   computed from it and must not be supplied.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance,true_kurtosis Passed to \code{\link{ss_aipe_two_phase}} when
#'   the design is planned; see that page.
#' @param G The number of Monte Carlo replications (default 1000).
#' @param print_iter Logical; print the replication number as the loop
#'   runs.
#' @param filename Optional path for a comma separated file recording every
#'   replication (the estimate \code{mean_Y}, the two confidence limits, the interval
#'   width, and two indicators of whether the interval missed the population
#'   mean below or above): nothing is written when \code{filename} is
#'   \code{NULL} (the default), a new file with a header row is created
#'   otherwise, an existing file at that path is appended to, and a
#'   throwaway run should point it at \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (default) to draw from the current random number
#'   stream, or a number that makes the run reproducible without disturbing
#'   the stream of the calling script.
#'
#' @details
#' \strong{The data generating mechanism.} With \code{population = NULL},
#' each replication draws the \eqn{n} labeled units as (target, proxy)
#' pairs from a bivariate normal distribution with standard deviations
#' \code{true_sd_Y} and 1, correlation \code{true_cor}, and means 0, and the
#' \eqn{M - n} unlabeled units' proxies from the same proxy distribution.
#' The interval's width and coverage do not depend on the two means or on
#' the proxy's scale, so these choices lose no generality; the population
#' target mean the intervals are scored against is then 0. Because the
#' regression estimator uses the unlabeled proxies only through their mean
#' and their sum of squares, the simulation draws those two from their
#' exact joint distribution (independent normal and scaled chi square)
#' rather than drawing every unit of the pool, which keeps a large pool
#' cheap. With a \code{population} data frame, each replication resamples
#' the labeled pairs and the unlabeled proxies with replacement from its
#' rows, a population that keeps whatever skew, bounds, or nonlinearity
#' the observed pairs have; the target mean scored against is then the mean
#' of its \code{y} column. In the evaluation of a system whose outputs an
#' automated judge labels, for example, a pilot set of outputs carrying both
#' a human assessor's label and the judge's makes such a population.
#'
#' \strong{What is evaluated.} Each replication computes the regression
#' interval of \code{\link{ci_two_phase}} (the interval the planner plans
#' for) and records its width and whether it missed the population target
#' mean below or above. A replication whose labeled units all share one
#' proxy value has no regression interval, and one whose labeled units all
#' share one target value would have an interval of zero width, which
#' \code{\link{ci_two_phase}} refuses. Either can happen only with a
#' \code{population} data frame, chiefly one whose proxy takes few distinct
#' values or whose target is binary with a proportion near 0 or 1. Such
#' replications are dropped, with a single warning that counts each kind,
#' and the number kept is \code{attr(x, "simulation")$returned} (see the
#' section on the record of the simulation). Many dropped
#' replications say that the design has too few labeled units for that
#' population, and the summaries, which describe the replications that
#' produced an interval, then flatter the design.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with \code{term} and
#'   \code{value} columns: \code{mean_mean_Y}, \code{median_mean_Y}, and
#'   \code{sd_mean_Y} for the estimates of the target mean; the realized
#'   widths (\code{mean_ci_width}, \code{median_ci_width},
#'   \code{sd_ci_width}); the proportion of intervals no wider than
#'   \code{width} (\code{pct_ci_less_width}); the proportions that missed the
#'   population mean from above (\code{pct_ci_miss_low}, the mean below the
#'   lower limit) and from below (\code{pct_ci_miss_high}) and their sum
#'   (\code{total_type_I_error}), all on the 0 to 1 scale; then
#'   \code{n_labeled}, \code{total_N} (the pool size \eqn{M}, the design's
#'   total number of units, which the rest of the two-phase family calls
#'   \code{n_pool}; the sensitivity family's shared schema names the size
#'   of a whole design \code{total_N}),
#'   \code{true_mean_Y} (the population target mean), \code{true_sd_Y},
#'   \code{true_cor}, \code{supposed_sd_Y} and \code{supposed_cor}
#'   (\code{NA} when \code{specified_n_labeled} was given), \code{width},
#'   \code{conf_level}, and, only when one was supplied, \code{assurance}
#'   and \code{true_kurtosis}. The number of replications used is in the
#'   \code{simulation} attribute (its element \code{returned}), described
#'   below.
#'
#' @references
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Neyman, J. (1938). Contribution to the theory of sampling human
#'   populations. \emph{Journal of the American Statistical Association,
#'   33}(201), 101--116. \doi{10.1080/01621459.1938.10503378}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_two_phase}}, \code{\link{ci_two_phase}},
#'   \code{\link{ss_aipe_two_phase_cost}}
#'
#' @family AIPE sample size planning
#'
#' @keywords design htest
#'
#' @examples
#' # A plan made with a correlation of 0.55, carried out where the
#' # correlation is only 0.45. G = 200 keeps the example quick; a reported
#' # evaluation deserves the default G = 1000 or more.
#' ss_aipe_two_phase_sensitivity(true_sd_Y = 0.43, true_cor = 0.45,
#'                               supposed_sd_Y = 0.43, supposed_cor = 0.55,
#'                               n_pool = 20000, width = 0.10,
#'                               conf_level = 0.90, G = 200, seed = 113)
#'
#' # A fixed design, 150 labeled units in a pool of 5,000, evaluated in a
#' # population of binary (target, proxy) pairs resampled from a pilot of
#' # 300 units: how often is the 90% interval narrower than 0.085?
#' set.seed(113)
#' proxy  <- rbinom(300, 1, 0.6)
#' target <- ifelse(runif(300) < 0.8, proxy, rbinom(300, 1, 0.6))
#' ss_aipe_two_phase_sensitivity(specified_n_labeled = 150, n_pool = 5000,
#'                               width = 0.085, conf_level = 0.90,
#'                               population = data.frame(y = target, x = proxy),
#'                               G = 200, seed = 113)
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_two_phase_sensitivity <- function(true_sd_Y = NULL, true_cor = NULL,
                                          supposed_sd_Y = NULL,
                                          supposed_cor = NULL,
                                          specified_n_labeled = NULL,
                                          n_pool, width, population = NULL,
                                          conf_level = 0.95, assurance = NULL,
                                          true_kurtosis = 0, G = 1000,
                                          print_iter = FALSE, filename = NULL,
                                          seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  planned <- !is.null(supposed_sd_Y) || !is.null(supposed_cor)
  if (planned && !is.null(specified_n_labeled))
    .dmar_refuse_argument(c(
      "Supply 'supposed_sd_Y' and 'supposed_cor', or ",
      "'specified_n_labeled', but not both."),
      arg = c("supposed_sd_Y", "supposed_cor", "specified_n_labeled"))
  if (!planned && is.null(specified_n_labeled))
    .dmar_refuse_argument(c(
      "Supply either 'supposed_sd_Y' and 'supposed_cor' (to plan the ",
      "number of labeled units) or 'specified_n_labeled'."),
      arg = c("supposed_sd_Y", "supposed_cor", "specified_n_labeled"))
  if (planned && (is.null(supposed_sd_Y) || is.null(supposed_cor)))
    .dmar_refuse_argument(c(
      "Supply both 'supposed_sd_Y' and 'supposed_cor' to plan the ",
      "number of labeled units."), arg = c("supposed_sd_Y", "supposed_cor"))
  .two_phase_check_planning(sd_y = supposed_sd_Y, rho = supposed_cor,
                            width = width, conf_level = conf_level,
                            assurance = assurance, true_kurtosis = true_kurtosis,
                            sd_name = "supposed_sd_Y",
                            rho_name = "supposed_cor")
  if (!is.numeric(n_pool) || length(n_pool) != 1L || !is.finite(n_pool) ||
      n_pool < 4 || n_pool != round(n_pool))
    .dmar_refuse_argument(
      "'n_pool' must be a single finite whole number of at least 4.",
      arg = "n_pool")
  G <- .dmar_check_count(G, "G")
  .check_filename(filename)

  if (is.null(population)) {
    .two_phase_check_planning(sd_y = true_sd_Y, rho = true_cor,
                              sd_name = "true_sd_Y", rho_name = "true_cor")
    if (is.null(true_sd_Y) || is.null(true_cor))
      .dmar_refuse_argument(c(
        "Supply 'true_sd_Y' and 'true_cor' (the population that generates ",
        "the data), or a 'population' data frame."),
        arg = c("true_sd_Y", "true_cor", "population"))
    true_mean_Y <- 0
  } else {
    if (!is.null(true_sd_Y) || !is.null(true_cor))
      .dmar_refuse_argument(c(
        "With a 'population' data frame, 'true_sd_Y' and 'true_cor' are ",
        "computed from it; leave them NULL."),
        arg = c("population", "true_sd_Y", "true_cor"))
    if (!is.data.frame(population) || !all(c("y", "x") %in% names(population)))
      .dmar_refuse_argument(c(
        "'population' must be a data frame with columns 'y' (the target) ",
        "and 'x' (the proxy)."), arg = "population")
    pop_y <- population$y; pop_x <- population$x
    if (is.logical(pop_y)) pop_y <- as.numeric(pop_y)
    if (is.logical(pop_x)) pop_x <- as.numeric(pop_x)
    if (!is.numeric(pop_y) || !is.numeric(pop_x) ||
        anyNA(pop_y) || anyNA(pop_x) ||
        any(!is.finite(pop_y)) || any(!is.finite(pop_x)))
      .dmar_refuse_argument(c(
        "The 'y' and 'x' columns of 'population' must be numeric with no ",
        "missing or infinite values."), arg = "population")
    P <- length(pop_y)
    if (P < 3L || stats::sd(pop_x) == 0 || stats::sd(pop_y) == 0)
      .dmar_refuse_argument(c(
        "'population' needs at least 3 rows and must vary in both 'y' ",
        "and 'x'."), arg = "population")
    true_mean_Y <- mean(pop_y)
    # The population is the finite set being resampled, so its moments
    # take the divisor P.
    true_sd_Y  <- sqrt(mean((pop_y - true_mean_Y)^2))
    true_cor   <- stats::cor(pop_y, pop_x)
  }

  if (planned) {
    plan <- ss_aipe_two_phase(true_sd_Y = supposed_sd_Y,
                              true_cor = supposed_cor,
                              width = width, n_pool = n_pool,
                              conf_level = conf_level, assurance = assurance,
                              true_kurtosis = true_kurtosis)
    n <- plan$value[plan$term == "necessary_n_labeled"]
  } else {
    n <- specified_n_labeled
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 3 ||
        n != round(n))
      .dmar_refuse_argument(c(
        "'specified_n_labeled' must be a single whole number of at ",
        "least 3."), arg = "specified_n_labeled")
  }
  if (n >= n_pool)
    .dmar_refuse_argument(c(
      "The design has ", n, " labeled units in a pool of ", n_pool,
      ", which leaves no unit with the proxy alone; 'n_pool' must exceed ",
      "the number of labeled units."), arg = "n_pool")
  n_u <- n_pool - n

  .dmar_local_seed(seed)
  mean_Y <- lower <- upper <- rep(NA_real_, G)
  # Replications without a usable interval, by reason: the proxy constant on
  # the labeled units (no regression) or the target constant on them (an
  # interval of zero width). Only a resampled population can produce either.
  constant_x <- constant_y <- 0L
  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    if (is.null(population)) {
      x <- stats::rnorm(n)
      y <- true_sd_Y * (true_cor * x + sqrt(1 - true_cor^2) * stats::rnorm(n))
      mean_u <- stats::rnorm(1L, sd = 1 / sqrt(n_u))
      ss_u <- if (n_u > 1) stats::rchisq(1L, df = n_u - 1) else 0
    } else {
      idx <- sample.int(P, n, replace = TRUE)
      y <- pop_y[idx]; x <- pop_x[idx]
      xu <- pop_x[sample.int(P, n_u, replace = TRUE)]
      mean_u <- mean(xu)
      ss_u <- sum((xu - mean_u)^2)
      if (all(x == x[1L])) {
        constant_x <- constant_x + 1L
        next
      }
      if (all(y == y[1L])) {
        constant_y <- constant_y + 1L
        next
      }
    }
    r <- .ci_two_phase_fast(y, x, n_u, mean_u, ss_u, conf_level)
    mean_Y[g] <- r[["mean_Y"]]
    lower[g] <- r[["ci_lower"]]
    upper[g] <- r[["ci_upper"]]
  }

  ci_width <- upper - lower
  miss_low <- true_mean_Y < lower
  miss_high <- true_mean_Y > upper
  if (!is.null(filename)) {
    per_rep <- data.frame(mean_Y = mean_Y, ci_lower = lower, ci_upper = upper,
                          ci_width = ci_width, type_I_lower = miss_low,
                          type_I_upper = miss_high)
    .write_sensitivity_csv(per_rep, filename)
  }

  ok <- !is.na(ci_width)
  G_used <- sum(ok)
  # Replications whose interval has a missing limit for any reason other
  # than the two above are counted too, so the message accounts for every
  # replication dropped, and only the constant draws are put down to the
  # number of labeled units. (A population whose values are so large that
  # their sums pass the largest double gives such intervals.) G is an
  # integer (.dmar_check_count()), and so is every count here, so each
  # prints in full, not with an exponent.
  other <- G - G_used - constant_x - constant_y
  reasons <- c(
    if (constant_x > 0L)
      paste0(constant_x, " drew labeled units that all share one proxy ",
             "value, which leaves no regression interval"),
    if (constant_y > 0L)
      paste0(constant_y, " drew labeled units that all share one target ",
             "value, which would give an interval of zero width"),
    if (other > 0L)
      paste0(other, " gave an interval with a missing limit"))
  if (G_used == 0L) {
    none <- c("No replication produced an interval: of ", G,
              " replications, ", paste(reasons, collapse = ", and "), ".")
    # Constant draws are a property of the design for this population, so
    # the arguments are at fault; a missing limit for another reason is a
    # computation that gave no result with arguments that were accepted.
    if (other == 0L)
      .dmar_refuse_argument(c(none, " The design has too few labeled ",
                              "units for this population."), arg = NULL)
    .dmar_refuse(none, "computation_failed")
  }
  if (G_used < G)
    .dmar_warn(c(
      "Of ", G, " replications, ", paste(reasons, collapse = ", and "),
      "; they were dropped, and the summaries use the remaining ",
      G_used, "."), "result_qualified")
  mean_Y <- mean_Y[ok]; ci_width <- ci_width[ok]
  miss_low <- miss_low[ok]; miss_high <- miss_high[ok]

  out <- data.frame(
    term  = c("mean_mean_Y", "median_mean_Y", "sd_mean_Y",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_labeled", "total_N",
              "true_mean_Y", "true_sd_Y", "true_cor",
              "supposed_sd_Y", "supposed_cor", "width", "conf_level",
              if (!is.null(assurance)) c("assurance", "true_kurtosis")),
    value = c(mean(mean_Y), stats::median(mean_Y),
              if (G_used > 1L) stats::sd(mean_Y) else NA_real_,
              mean(ci_width), stats::median(ci_width),
              if (G_used > 1L) stats::sd(ci_width) else NA_real_,
              mean(ci_width <= width),
              mean(miss_low), mean(miss_high), mean(miss_low | miss_high),
              n, n_pool,
              true_mean_Y, true_sd_Y, true_cor,
              if (planned) supposed_sd_Y else NA_real_,
              if (planned) supposed_cor else NA_real_,
              width, conf_level,
              if (!is.null(assurance)) c(assurance, true_kurtosis)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G, ok,
                  conf_level)
}
