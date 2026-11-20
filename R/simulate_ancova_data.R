#' Simulate Data From a One-Covariate ANCOVA Model
#'
#' Generates random data appropriate for an analysis of covariance with one
#' continuous outcome (\eqn{Y}) and one continuous covariate (\eqn{X}) crossed
#' with \eqn{a} fixed groups. The covariate is treated as a random variable;
#' \eqn{Y} and \eqn{X} are jointly multivariate normal within each group. Both
#' the randomized-design case (population covariate mean common across groups)
#' and the non-randomized / preexisting-groups case (population covariate
#' means differ across groups; \eqn{Y}-on-\eqn{X} correlation may also differ)
#' are supported. Per-group sample sizes may be equal or unequal.
#'
#' @param true_mean_Y A numeric vector of length \code{a} giving the population mean
#'   of \eqn{Y} in each group.
#' @param true_mean_X When \code{randomized = TRUE}, a single number giving the
#'   common population mean of \eqn{X} across all groups. When
#'   \code{randomized = FALSE}, a numeric vector of length \code{a} giving
#'   the population covariate mean in each group.
#' @param true_sd_Y The population standard deviation of \eqn{Y} (assumed
#'   common across groups).
#' @param true_sd_X The population standard deviation of \eqn{X} (assumed
#'   common across groups).
#' @param true_cor The population correlation between \eqn{Y} and \eqn{X}. When
#'   \code{randomized = TRUE}, must be a single number (see Details). When
#'   \code{randomized = FALSE}, may be a single number (recycled to all
#'   \code{a} groups) or a numeric vector of length \code{a} giving a
#'   distinct correlation in each group.
#' @param a The number of fixed levels of the grouping factor (i.e., the
#'   number of conditions in a fixed-effects ANCOVA design). Use this
#'   argument when groups are the design levels of interest. (For
#'   sample-selected groups, e.g., classrooms or schools randomly drawn
#'   from a population, the convention in this package is to use
#'   \code{J} instead.)
#' @param n_per_group A single number (equal sample size per group) or a numeric
#'   vector of length \code{a} giving the sample size in each group.
#' @param randomized Logical. \code{TRUE} (the default) for a randomized
#'   design (random assignment of subjects to groups, so that the
#'   population covariate mean is the same in every group). \code{FALSE}
#'   for a non-randomized / preexisting-groups design.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A long-format \code{data.frame} with one row per simulated subject
#'   and three columns:
#'   \describe{
#'     \item{\code{group}}{A factor with \code{a} levels (\code{"1"}, ...,
#'       \code{as.character(a)}) identifying each subject's group.}
#'     \item{\code{y}}{Numeric simulated outcome.}
#'     \item{\code{x}}{Numeric simulated covariate.}
#'   }
#'   This format is directly usable with \code{aov()}, \code{lm()}, and other
#'   model-fitting functions.
#'
#' @details
#' Each group's \eqn{(Y, X)} pairs are drawn from a bivariate normal
#' distribution with mean \eqn{(\mu_{Y,j}, \mu_{X,j})} and covariance
#' \deqn{\Sigma_j = \begin{pmatrix} \sigma_Y^2 & \rho_j\,\sigma_Y\,\sigma_X \\ \rho_j\,\sigma_Y\,\sigma_X & \sigma_X^2 \end{pmatrix}.}
#'
#' \strong{Why \code{true_cor} must be a single number when \code{randomized = TRUE}.}
#' Random assignment forms each group as an exchangeable random sample from
#' the same population. The bivariate distribution of \eqn{(Y, X)} is
#' therefore the same in every group, including the correlation. Allowing
#' \code{true_cor} to differ across groups would silently break that
#' interpretation and produce data that no randomized design could plausibly
#' have generated. The function therefore stops with an error in that case;
#' use \code{randomized = FALSE} if you genuinely want group-specific
#' correlations.
#'
#' \strong{Convention on group labels.} The argument \code{a} is used here
#' (and throughout DMAR's experimental-design functions) for the number
#' of \emph{fixed} levels of a designed factor, the levels you intend to
#' compare. The letter \code{J} is reserved for the number of \emph{sample-
#' selected} groups, e.g., when classrooms or schools are randomly sampled
#' from a population (a random-effects context).
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link[MASS]{mvrnorm}}, \code{\link{ci_c_ancova}},
#'   \code{\link{ss_aipe_c_ancova}}
#'
#' @family data simulators
#'
#' @keywords design datagen
#'
#' @examples
#' # 1. Randomized design, two groups, equal n.
#' set.seed(113)
#' simple <- simulate_ancova_data(
#'   true_mean_Y = c(3, 5),
#'   true_mean_X = 10,
#'   true_sd_Y    = 1,
#'   true_sd_X    = 2,
#'   true_cor   = 0.8,
#'   a          = 2,
#'   n_per_group = 20
#' )
#' head(simple)
#'
#' # 2. Four preexisting groups with different correlations and unequal n.
#' #    The first two groups share true_cor = 0.30; the second two share a larger
#' #    true_cor = 0.60. The four groups enter the data generating mechanism
#' #    only through their per-group means and correlations, as four
#' #    distinct populations being sampled. In a downstream
#' #    analysis these four groups could be cross-classified as a 2 x 2
#' #    factorial design (e.g., the first factor distinguishing groups 1-2
#' #    from groups 3-4, and the second factor distinguishing groups 1, 3
#' #    from groups 2, 4) and analyzed via factorial ANCOVA.
#' set.seed(113)
#' preexisting <- simulate_ancova_data(
#'   true_mean_Y = c(50, 55, 60, 65),
#'   true_mean_X = c(10, 12, 11, 13),
#'   true_sd_Y    = 8,
#'   true_sd_X    = 3,
#'   true_cor   = c(0.30, 0.30, 0.60, 0.60),
#'   a          = 4,
#'   n_per_group = c(40, 35, 45, 30),
#'   randomized = FALSE
#' )
#' aggregate(cbind(y, x) ~ group, data = preexisting,
#'           FUN = function(z) round(c(mean = mean(z), sd = sd(z)), 2))
#'
#' @export
#' @importFrom MASS mvrnorm
simulate_ancova_data <- function(true_mean_Y, true_mean_X, true_sd_Y, true_sd_X, true_cor,
                                 a, n_per_group, randomized = TRUE,
                                 seed = NULL) {

  # ---------- Validate scalar arguments ----------
  if (!is.numeric(a) || length(a) != 1L || a < 2L || a != round(a)) {
    .dmar_refuse_argument(
      "'a' must be a single integer >= 2 (the number of fixed groups).",
      arg = "a")
  }
  if (!is.numeric(true_sd_Y) || length(true_sd_Y) != 1L || true_sd_Y <= 0) {
    .dmar_refuse_argument(
      "'true_sd_Y' must be a single positive number.", arg = "true_sd_Y")
  }
  if (!is.numeric(true_sd_X) || length(true_sd_X) != 1L || true_sd_X <= 0) {
    .dmar_refuse_argument(
      "'true_sd_X' must be a single positive number.", arg = "true_sd_X")
  }
  if (!is.logical(randomized) || length(randomized) != 1L) {
    .dmar_refuse_argument(
      "'randomized' must be a single logical value.", arg = "randomized")
  }

  # ---------- Validate vector arguments ----------
  if (!is.numeric(true_mean_Y) || length(true_mean_Y) != a) {
    .dmar_refuse_argument(
      "'true_mean_Y' must be a numeric vector of length 'a' (one population mean of Y per group).",
      arg = c("true_mean_Y", "a"))
  }

  if (isTRUE(randomized)) {
    if (length(true_mean_X) != 1L) {
      .dmar_refuse_argument(c(
        "Under 'randomized = TRUE', 'true_mean_X' must be a single number ",
        "(the common population covariate mean across groups)."), arg = "true_mean_X")
    }
    if (length(true_cor) != 1L) {
      .dmar_refuse_argument(c(
        "Under 'randomized = TRUE', 'true_cor' must be a single number. ",
        "Random assignment forms each group as a sample from a common ",
        "population, so the Y-X correlation is by definition the same ",
        "in every group; allowing it to differ would contradict the ",
        "randomized-design assumption. Use 'randomized = FALSE' to ",
        "simulate group-specific correlations."), arg = "true_cor")
    }
    mu_x_per_group  <- rep(as.numeric(true_mean_X),  a)
    rho_per_group   <- rep(as.numeric(true_cor),   a)
  } else {
    if (length(true_mean_X) != a) {
      .dmar_refuse_argument(c(
        "Under 'randomized = FALSE', 'true_mean_X' must be a numeric vector of length 'a' ",
        "(one population covariate mean per group)."), arg = c("true_mean_X", "a"))
    }
    if (length(true_cor) == 1L) {
      rho_per_group <- rep(as.numeric(true_cor), a)
    } else if (length(true_cor) == a) {
      rho_per_group <- as.numeric(true_cor)
    } else {
      .dmar_refuse_argument(c(
        "Under 'randomized = FALSE', 'true_cor' must be either a single number ",
        "(common correlation, recycled across all groups) or a numeric ",
        "vector of length 'a' (one correlation per group)."),
        arg = c("true_cor", "a"))
    }
    mu_x_per_group <- as.numeric(true_mean_X)
  }

  if (any(abs(rho_per_group) >= 1)) {
    .dmar_refuse_argument(
      "All correlations in 'true_cor' must lie strictly between -1 and 1.",
      arg = "true_cor")
  }

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
    .dmar_refuse_argument(c(
      "'n_per_group' must be a single number (equal sizes) or a numeric vector ",
      "of length 'a' (one sample size per group)."),
      arg = c("n_per_group", "a"))
  }

  # ---------- Generate the data ----------
  if (!requireNamespace("MASS", quietly = TRUE)) {
    .dmar_refuse(
      "Package 'MASS' is required. Install with install.packages(\"MASS\").",
      "package_required", package = "MASS")
  }
  .dmar_local_seed(seed)

  N_total  <- sum(n_per_group)
  group    <- factor(rep(seq_len(a), times = n_per_group),
                     levels = seq_len(a),
                     labels = as.character(seq_len(a)))
  y        <- numeric(N_total)
  x        <- numeric(N_total)

  pos <- 1L
  for (j in seq_len(a)) {
    n_j     <- n_per_group[j]
    cov_yx  <- rho_per_group[j] * true_sd_Y * true_sd_X
    Sigma   <- matrix(c(true_sd_Y^2, cov_yx,
                        cov_yx,    true_sd_X^2),
                      nrow = 2)
    sample_j <- MASS::mvrnorm(
      n_j,
      mu    = c(true_mean_Y[j], mu_x_per_group[j]),
      Sigma = Sigma
    )
    idx        <- pos:(pos + n_j - 1L)
    y[idx]     <- sample_j[, 1L]
    x[idx]     <- sample_j[, 2L]
    pos        <- pos + n_j
  }

  data.frame(group = group, y = y, x = x, stringsAsFactors = FALSE)
}
