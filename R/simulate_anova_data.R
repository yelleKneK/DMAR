#' Simulate Data From a One-Way Fixed-Effects ANOVA Model
#'
#' Generates random data appropriate for a one-way fixed-effects analysis of
#' variance. Each group's observations are drawn from a normal distribution
#' with that group's population mean and a common (or per-group) standard
#' deviation. Per-group sample sizes may be equal or unequal.
#'
#' @param true_mean A numeric vector of length \code{a} giving the population mean
#'   in each group.
#' @param true_sd Within-group population standard deviation. Either a single
#'   number (homoscedastic; common across groups) or a numeric vector of
#'   length \code{a} (heteroscedastic; one SD per group).
#' @param a The number of fixed levels of the grouping factor (per the
#'   convention used throughout DMAR for fixed-factor designs).
#' @param n_per_group A single number (equal sample size per group) or a numeric
#'   vector of length \code{a} giving the sample size in each group.
#' @param seed Optional integer random seed for reproducibility (default
#'   \code{NULL}; supply an integer such as \code{113} for reproducible
#'   output).
#'
#' @return A long-format \code{data.frame} with one row per simulated
#'   subject and two columns:
#'   \describe{
#'     \item{\code{group}}{A factor with \code{a} levels.}
#'     \item{\code{y}}{Numeric simulated outcome.}
#'   }
#'
#' @details
#' The fixed-effects ANOVA model assumes group-specific means and a common
#' within-group variance. Setting \code{true_sd} to a vector relaxes the
#' homoscedasticity assumption; in that case the simulated data violate the
#' standard ANOVA assumption, which allows studying the robustness of the
#' ANOVA \emph{F} test or the performance of Welch-style alternatives.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{simulate_ancova_data}},
#'   \code{\link{simulate_regression_data}}, \code{\link{contrast_test}},
#'   \code{\link{ss_power_contrast}}
#'
#' @family data simulators
#'
#' @keywords design datagen
#'
#' @examples
#' # Three-group ANOVA, equal n per group.
#' set.seed(113)
#' d <- simulate_anova_data(true_mean = c(50, 55, 60), true_sd = 8, a = 3, n_per_group = 30)
#' aggregate(y ~ group, data = d, FUN = mean)
#'
#' # Same design with unequal n per group.
#' simulate_anova_data(true_mean = c(50, 55, 60), true_sd = 8, a = 3,
#'                     n_per_group = c(40, 30, 20), seed = 113)
#'
#' # Heteroscedastic case: each group has its own SD.
#' simulate_anova_data(true_mean = c(50, 55, 60), true_sd = c(5, 8, 12),
#'                     a = 3, n_per_group = 30, seed = 113)
#'
#' @export
#' @import stats
simulate_anova_data <- function(true_mean, true_sd, a, n_per_group, seed = NULL) {

  if (!is.numeric(a) || length(a) != 1L || a < 2L || a != round(a)) {
    .dmar_refuse_argument("'a' must be a single integer >= 2.", arg = "a")
  }
  if (!is.numeric(true_mean) || length(true_mean) != a) {
    .dmar_refuse_argument(
      "'true_mean' must be a numeric vector of length 'a' (one population mean per group).",
      arg = c("true_mean", "a"))
  }

  # Sigma: scalar (homoscedastic) or vector of length a (heteroscedastic).
  if (length(true_sd) == 1L) {
    if (!is.numeric(true_sd) || true_sd <= 0) {
      .dmar_refuse_argument(c(
        "'true_sd' must be a single positive number, or a numeric vector ",
        "of length 'a' of positive numbers."), arg = c("true_sd", "a"))
    }
    sigma_per_group <- rep(as.numeric(true_sd), a)
  } else if (length(true_sd) == a) {
    if (!is.numeric(true_sd) || any(true_sd <= 0)) {
      .dmar_refuse_argument(
        "Every element of 'true_sd' must be a positive number.", arg = "true_sd")
    }
    sigma_per_group <- as.numeric(true_sd)
  } else {
    .dmar_refuse_argument(
      "'true_sd' must be a single number or a numeric vector of length 'a'.",
      arg = c("true_sd", "a"))
  }

  # n_per_group: scalar (equal) or vector of length a (unequal)
  if (length(n_per_group) == 1L) {
    if (!is.numeric(n_per_group) || n_per_group < 2L || n_per_group != round(n_per_group)) {
      .dmar_refuse_argument(
        "'n_per_group' must be a positive integer >= 2.", arg = "n_per_group")
    }
    n_per_group <- rep(as.integer(n_per_group), a)
  } else if (length(n_per_group) == a) {
    if (!is.numeric(n_per_group) || any(n_per_group < 2L) || any(n_per_group != round(n_per_group))) {
      .dmar_refuse_argument(
        "Each entry of 'n_per_group' must be a positive integer >= 2.",
        arg = "n_per_group")
    }
    n_per_group <- as.integer(n_per_group)
  } else {
    .dmar_refuse_argument(
      "'n_per_group' must be a single number or a numeric vector of length 'a'.",
      arg = c("n_per_group", "a"))
  }

  .dmar_local_seed(seed)

  N_total <- sum(n_per_group)
  group   <- factor(rep(seq_len(a), times = n_per_group),
                    levels = seq_len(a),
                    labels = as.character(seq_len(a)))
  y       <- numeric(N_total)
  pos     <- 1L
  for (j in seq_len(a)) {
    n_j <- n_per_group[j]
    y[pos:(pos + n_j - 1L)] <- stats::rnorm(n_j, mean = true_mean[j],
                                            sd   = sigma_per_group[j])
    pos <- pos + n_j
  }
  data.frame(group = group, y = y, stringsAsFactors = FALSE)
}
