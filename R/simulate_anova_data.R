#' Simulate Data From a One-Way Fixed-Effects ANOVA Model
#'
#' Generates random data appropriate for a one-way fixed-effects analysis of
#' variance. Each group's observations are drawn from a normal distribution
#' with that group's population mean and a common (or per-group) standard
#' deviation. Per-group sample sizes may be equal or unequal.
#'
#' @param mu A numeric vector of length \code{a} giving the population mean
#'   in each group.
#' @param sigma Within-group population standard deviation. Either a single
#'   number (homoscedastic; common across groups) or a numeric vector of
#'   length \code{a} (heteroscedastic; one SD per group).
#' @param a The number of fixed levels of the grouping factor (per the
#'   convention used throughout DMAR for fixed-factor designs).
#' @param n A single number (equal sample size per group) or a numeric
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
#' within-group variance. Setting \code{sigma} to a vector relaxes the
#' homoscedasticity assumption; in that case the simulated data violate the
#' standard ANOVA assumption (a useful feature for studying robustness or
#' the performance of Welch-style alternatives).
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
#' d <- simulate_anova_data(mu = c(50, 55, 60), sigma = 8, a = 3, n = 30)
#' aggregate(y ~ group, data = d, FUN = mean)
#'
#' # Same design with unequal n per group.
#' simulate_anova_data(mu = c(50, 55, 60), sigma = 8, a = 3,
#'                     n = c(40, 30, 20), seed = 113)
#'
#' # Heteroscedastic case: each group has its own SD.
#' simulate_anova_data(mu = c(50, 55, 60), sigma = c(5, 8, 12),
#'                     a = 3, n = 30, seed = 113)
#'
#' @export
#' @import stats
simulate_anova_data <- function(mu, sigma, a, n, seed = NULL) {

  if (!is.numeric(a) || length(a) != 1L || a < 2L || a != round(a)) {
    stop("'a' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(mu) || length(mu) != a) {
    stop("'mu' must be a numeric vector of length 'a' (one population mean per group).",
         call. = FALSE)
  }

  # Sigma: scalar (homoscedastic) or vector of length a (heteroscedastic).
  if (length(sigma) == 1L) {
    if (!is.numeric(sigma) || sigma <= 0) {
      stop("'sigma' must be a single positive number, or a numeric vector ",
           "of length 'a' of positive numbers.", call. = FALSE)
    }
    sigma_per_group <- rep(as.numeric(sigma), a)
  } else if (length(sigma) == a) {
    if (!is.numeric(sigma) || any(sigma <= 0)) {
      stop("Every element of 'sigma' must be a positive number.",
           call. = FALSE)
    }
    sigma_per_group <- as.numeric(sigma)
  } else {
    stop("'sigma' must be a single number or a numeric vector of length 'a'.",
         call. = FALSE)
  }

  # n: scalar (equal) or vector of length a (unequal)
  if (length(n) == 1L) {
    if (!is.numeric(n) || n < 2L || n != round(n)) {
      stop("'n' must be a positive integer >= 2.", call. = FALSE)
    }
    n_per_group <- rep(as.integer(n), a)
  } else if (length(n) == a) {
    if (!is.numeric(n) || any(n < 2L) || any(n != round(n))) {
      stop("Each entry of 'n' must be a positive integer >= 2.",
           call. = FALSE)
    }
    n_per_group <- as.integer(n)
  } else {
    stop("'n' must be a single number or a numeric vector of length 'a'.",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      .old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  N_total <- sum(n_per_group)
  group   <- factor(rep(seq_len(a), times = n_per_group),
                    levels = seq_len(a),
                    labels = as.character(seq_len(a)))
  y       <- numeric(N_total)
  pos     <- 1L
  for (j in seq_len(a)) {
    n_j <- n_per_group[j]
    y[pos:(pos + n_j - 1L)] <- stats::rnorm(n_j, mean = mu[j],
                                            sd   = sigma_per_group[j])
    pos <- pos + n_j
  }
  data.frame(group = group, y = y, stringsAsFactors = FALSE)
}
