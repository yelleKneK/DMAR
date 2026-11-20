# Monte Carlo sensitivity analysis for AIPE planning of a single mean.
#' Sensitivity Analysis for Sample Size Planning for a Mean From the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description
#' Evaluates a sample size plan for the confidence interval for a mean by
#' Monte Carlo simulation. Each of \code{G} replications draws a sample of
#' the planned (or a specified) size from the data generating population,
#' forms the \emph{t} interval, and records its width and whether it missed
#' the true mean. The summary shows what a plan from
#' \code{\link{ss_aipe_mean}} delivers when the planning value of the
#' standard deviation is wrong, when the population is not normal, or both.
#' The population is either normal, with mean \code{true_mean} and standard
#' deviation \code{true_sd}, or the empirical distribution of a vector of
#' observed scores, resampled with replacement, which is how a plan is
#' checked against the skewness, tails, and discreteness of real data such
#' as a pilot's difference scores.
#'
#' @param true_sd The data generating standard deviation of a normal
#'   population. Supply it, or \code{population}, but not both.
#' @param supposed_sd The supposed standard deviation, the planning value
#'   passed to \code{\link{ss_aipe_mean}} to plan the sample size (not an
#'   estimate computed from data). Supply this or
#'   \code{specified_N}, but not both.
#' @param specified_N A sample size to evaluate in place of a planned one;
#'   the planner is then bypassed. A single whole number of at least 2.
#' @param width The desired full width of the confidence interval: the
#'   planning target, and the bound against which \code{pct_ci_less_width}
#'   counts the realized widths.
#' @param true_mean The data generating mean of a normal population
#'   (default 0). The width does not depend on it; the misses are counted
#'   against it. Not supplied with \code{population}.
#' @param population \code{NULL} (the default) for a normal population, or
#'   a numeric vector of scores whose empirical distribution is the
#'   population: each replication draws its sample from these scores with
#'   replacement. The true mean is then their mean and the true standard
#'   deviation their standard deviation with divisor equal to their number,
#'   the standard deviation of the distribution resampled. Missing values
#'   are refused rather than dropped.
#' @param conf_level The confidence level of the interval (default
#'   \code{0.95}).
#' @param assurance \code{NULL} (the default) or the assurance passed to
#'   \code{\link{ss_aipe_mean}} when it plans the sample size; when one is
#'   supplied it is also echoed, as the value \code{pct_ci_less_width} is read
#'   against.
#' @param true_kurtosis The planning value of the excess kurtosis passed to
#'   \code{\link{ss_aipe_mean}} (default 0), on the scale
#'   \code{\link{kurtosis}} reports. It changes the plan only together with
#'   \code{supposed_sd} and an \code{assurance}; the function warns when a
#'   nonzero value would go unused. As there, it must be greater than -2;
#'   the \code{true_kurtosis} argument of \code{\link{ss_aipe_mean}} says what
#'   to pass for a pilot whose \code{kurtosis()} is -2 or less.
#' @param G The number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print the index of each replication as it
#'   runs, which is helpful in a long run; default \code{FALSE}.
#' @param filename Optional path of a CSV file to receive the
#'   per-replication results (the sample mean and standard deviation, the
#'   interval limits and width, and the two tail misses), appended when the
#'   file already exists and created otherwise. The default \code{NULL}
#'   writes nothing; a throwaway run that wants the file should point it
#'   at \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the run reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}.
#'   The rows \code{mean_mean}, \code{median_mean}, and \code{sd_mean}
#'   summarize the \code{G} sample means; \code{mean_ci_width},
#'   \code{median_ci_width}, and \code{sd_ci_width} summarize the interval
#'   widths; \code{pct_ci_less_width} is the proportion of intervals no wider
#'   than \code{width}, to be read against the assurance of the plan (or
#'   against a value near one half for an expected width plan);
#'   \code{pct_ci_miss_low} is the proportion of intervals lying entirely
#'   above the true mean (the true mean below the lower limit) and
#'   \code{pct_ci_miss_high} the proportion lying entirely below it, with
#'   \code{total_type_I_error} their sum. These proportions are on the 0 to
#'   1 scale. The input echoes are \code{total_N} (the sample size
#'   evaluated), \code{true_sd} and \code{true_mean} (computed from
#'   \code{population} when one is supplied), \code{supposed_sd} (NA when
#'   \code{specified_N} was supplied instead), \code{width},
#'   \code{conf_level}, and, when an assurance was supplied,
#'   \code{assurance}. The attribute \code{population} records whether the
#'   population was \code{"normal"} or \code{"resampled"}.
#'
#' @details
#' \strong{What each replication does.} A sample of \code{total_N} scores
#' is drawn, and the interval
#' \eqn{\bar{y} \pm t_{1-\alpha/2,\,N-1}\,
#' S/\sqrt{N}}{ybar +/- t(1 - alpha/2, N - 1) S / sqrt(N)} is formed,
#' the one \code{t.test()} reports and the one the plan was made for. Every
#' replication returns an interval, so none is dropped. A resample whose
#' scores are all equal, which a discrete population allows at a small
#' \emph{N}, gives an interval of zero width, and it is counted as it
#' falls.
#'
#' \strong{What to read.} The width rows against \code{width} show the
#' effect of a planning value of the standard deviation that differs from
#' the true one: the width is proportional to \eqn{S}, so a population
#' standard deviation larger than the planning value widens every
#' interval. \code{pct_ci_less_width} against the plan's assurance shows
#' whether the plan delivered the probability it promised; with a resampled
#' population it is the check, free of the two-moment approximation, on
#' the planner's \code{true_kurtosis} adjustment. The two miss rows, each
#' against \eqn{\alpha/2}{alpha/2}, show the coverage of the \emph{t}
#' interval in the population at hand: under normality both tails miss at
#' the nominal rate, and a skewed population misses more often in one tail
#' than in the other, which the plan does not address.
#'
#' \strong{A pilot as the population.} Passing a pilot sample's scores as
#' \code{population} and its standard deviation as \code{supposed_sd}
#' evaluates the plan in a population with exactly the pilot's shape, so
#' the planning value is correct (up to the divisor) and what the
#' simulation isolates is the effect of the shape. The pilot is itself a
#' sample, so its shape, and above all its tails, carries sampling error of
#' its own.
#'
#' Sensitivity analysis by simulation, in which the population value that
#' generates the data and the value the plan assumed are set separately,
#' was introduced for the standardized mean difference by Kelley and
#' Rausch (2006); the \code{ss_aipe_*_sensitivity()} family carries it to
#' the other AIPE planners, with a common set of rows.
#'
#' @references
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_mean}}, \code{\link{kurtosis}},
#'   \code{\link{ss_aipe_c_sensitivity}}
#'
#' @examples
#' # The pilot of the ss_aipe_mean() page: the IQ gains of the 246 children
#' # in the control group of the pygmalion data.
#' gain <- pygmalion$iq_gain[pygmalion$treatment == "Control"]
#'
#' # A plan for 90 percent assurance made under normal theory, evaluated by
#' # resampling the pilot's own gains, whose tails are heavier than the
#' # normal's: the share of intervals no wider than 4 points (the row
#' # pct_ci_less_width) falls short of .90.
#' ss_aipe_mean_sensitivity(supposed_sd = sd(gain), width = 4,
#'                          population = gain, assurance = 0.90, seed = 113)
#'
#' # The plan made with the pilot's excess kurtosis is six children larger,
#' # and the share comes close to .90.
#' ss_aipe_mean_sensitivity(supposed_sd = sd(gain), width = 4,
#'                          population = gain, assurance = 0.90,
#'                          true_kurtosis = kurtosis(gain), seed = 113)
#'
#' # A normal population whose standard deviation is 15 when the plan
#' # assumed 13.3: the realized intervals run wider than the target.
#' ss_aipe_mean_sensitivity(true_sd = 15, supposed_sd = 13.3, width = 4,
#'                          seed = 113)
#'
#' @keywords design htest
#'
#' @family AIPE sample size planning
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_mean_sensitivity <- function(true_sd = NULL, supposed_sd = NULL,
                                     specified_N = NULL, width,
                                     true_mean = 0, population = NULL,
                                     conf_level = 0.95, assurance = NULL,
                                     true_kurtosis = 0, G = 1000,
                                     print_iter = FALSE, filename = NULL,
                                     seed = NULL) {
  if (is.null(supposed_sd) == is.null(specified_N))
    .dmar_refuse_argument(c(
      "Supply exactly one of 'supposed_sd' (to plan N with ",
      "ss_aipe_mean()) and 'specified_N' (to evaluate a given N)."),
      arg = c("supposed_sd", "specified_N"))

  # The data generating population: normal, or the empirical distribution of
  # the scores supplied, whose mean and standard deviation (divisor equal to
  # their number) are then the true values the echoes report.
  if (is.null(population)) {
    if (is.null(true_sd))
      .dmar_refuse_argument(c(
        "Supply 'true_sd' for a normal population or 'population', a ",
        "vector of scores to resample."), arg = c("true_sd", "population"))
    if (!is.numeric(true_sd) || length(true_sd) != 1L ||
        !is.finite(true_sd) || true_sd <= 0)
      .dmar_refuse_argument(
        "'true_sd' must be a single finite positive number.", arg = "true_sd")
    if (!is.numeric(true_mean) || length(true_mean) != 1L ||
        !is.finite(true_mean))
      .dmar_refuse_argument(
        "'true_mean' must be a single finite number.", arg = "true_mean")
    mu    <- true_mean
    sigma <- true_sd
  } else {
    if (!is.null(true_sd))
      .dmar_refuse_argument(c(
        "Supply 'true_sd' for a normal population or 'population' for a ",
        "resampled one, not both."), arg = c("true_sd", "population"))
    if (!missing(true_mean))
      .dmar_refuse_argument(c(
        "'true_mean' describes a normal population; with 'population' ",
        "the true mean is the mean of the scores supplied."),
        arg = c("true_mean", "population"))
    if (!is.numeric(population) || length(population) < 2L)
      .dmar_refuse_argument(
        "'population' must be a numeric vector of at least two scores.",
        arg = "population")
    if (!all(is.finite(population)))
      .dmar_refuse_argument(c(
        "'population' must not contain missing or infinite values; ",
        "remove them (for example with na.omit()) before resampling."),
        arg = "population")
    mu    <- mean(population)
    sigma <- sqrt(mean((population - mu)^2))
    if (sigma == 0)
      .dmar_refuse_argument(
        "'population' must not be constant.", arg = "population")
  }

  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) ||
      width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single finite positive number.", arg = "width")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument(
      "'conf_level' must be a single number strictly between 0 and 1.",
      arg = "conf_level")
  G <- .dmar_check_count(G, "G")
  if (!is.logical(print_iter) || length(print_iter) != 1L || is.na(print_iter))
    .dmar_refuse_argument(
      "'print_iter' must be TRUE or FALSE.", arg = "print_iter")
  .check_filename(filename)

  if (!is.null(supposed_sd)) {
    if (!is.numeric(supposed_sd) || length(supposed_sd) != 1L ||
        !is.finite(supposed_sd) || supposed_sd <= 0)
      .dmar_refuse_argument(
        "'supposed_sd' must be a single finite positive number.",
        arg = "supposed_sd")
    plan <- ss_aipe_mean(true_sd = supposed_sd, width = width,
                         conf_level = conf_level, assurance = assurance,
                         true_kurtosis = true_kurtosis)
    n <- plan$value[plan$term == "necessary_N"]
  } else {
    if (!is.numeric(specified_N) || length(specified_N) != 1L ||
        !is.finite(specified_N) || specified_N < 2 ||
        specified_N != floor(specified_N))
      .dmar_refuse_argument(
        "'specified_N' must be a single whole number of at least 2.",
        arg = "specified_N")
    if (!is.null(assurance) &&
        (!is.numeric(assurance) || length(assurance) != 1L ||
         !is.finite(assurance) || assurance <= 0.5 || assurance >= 1))
      .dmar_refuse_argument(c(
        "'assurance' must be NULL or a single number strictly between ",
        "0.5 and 1."), arg = "assurance")
    if (!is.numeric(true_kurtosis) || length(true_kurtosis) != 1L ||
        !is.finite(true_kurtosis) || true_kurtosis <= -2)
      .dmar_refuse_argument(c(
        "'true_kurtosis' must be a single finite number greater than -2 (the ",
        "excess kurtosis, 0 for a normal population); for a pilot whose ",
        "kurtosis() is -2 or less, see the 'true_kurtosis' argument of ",
        "?ss_aipe_mean."), arg = "true_kurtosis")
    if (true_kurtosis != 0)
      .dmar_warn(c(
        "'true_kurtosis' enters only the planner's assurance adjustment and ",
        "is unused when 'specified_N' is supplied."), "argument_ignored")
    n <- specified_N
  }

  .dmar_local_seed(seed)

  draw <- if (is.null(population)) {
    function() stats::rnorm(n, mean = mu, sd = sigma)
  } else {
    function() population[sample.int(length(population), n, replace = TRUE)]
  }
  mean_hat <- numeric(G)
  sd_hat   <- numeric(G)
  for (g in seq_len(G)) {
    if (print_iter) cat(g, "\n")
    y <- draw()
    mean_hat[g] <- mean(y)
    sd_hat[g]   <- stats::sd(y)
  }

  half     <- stats::qt(1 - (1 - conf_level) / 2, n - 1) * sd_hat / sqrt(n)
  ci_lower <- mean_hat - half
  ci_upper <- mean_hat + half
  ci_width <- 2 * half
  # The family's tails: a miss low is an interval lying entirely above the
  # true mean (the truth below its lower limit), a miss high the reverse.
  miss_low  <- mu < ci_lower
  miss_high <- mu > ci_upper

  if (!is.null(filename)) {
    per_rep <- data.frame(mean_hat = mean_hat, sd_hat = sd_hat,
                          ci_lower = ci_lower, ci_upper = ci_upper,
                          ci_width = ci_width, type_I_lower = miss_low,
                          type_I_upper = miss_high)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_mean", "median_mean", "sd_mean",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "true_sd", "true_mean", "supposed_sd",
              "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(mean_hat), stats::median(mean_hat), stats::sd(mean_hat),
              mean(ci_width), stats::median(ci_width), stats::sd(ci_width),
              mean(ci_width <= width),
              mean(miss_low), mean(miss_high), mean(miss_low | miss_high),
              n, sigma, mu,
              if (is.null(supposed_sd)) NA_real_ else supposed_sd,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "population") <- if (is.null(population)) "normal" else "resampled"
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                  !is.na(ci_width), conf_level)
}
