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
#' @param mu_y A numeric vector of length \code{a} giving the population mean
#'   of \eqn{Y} in each group.
#' @param mu_x When \code{randomized = TRUE}, a single number giving the
#'   common population mean of \eqn{X} across all groups. When
#'   \code{randomized = FALSE}, a numeric vector of length \code{a} giving
#'   the population covariate mean in each group.
#' @param sigma_y The population standard deviation of \eqn{Y} (assumed
#'   common across groups).
#' @param sigma_x The population standard deviation of \eqn{X} (assumed
#'   common across groups).
#' @param rho The population correlation between \eqn{Y} and \eqn{X}. When
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
#' @param n A single number (equal sample size per group) or a numeric
#'   vector of length \code{a} giving the sample size in each group.
#' @param randomized Logical. \code{TRUE} (the default) for a randomized
#'   design (random assignment of subjects to groups, so that the
#'   population covariate mean is the same in every group). \code{FALSE}
#'   for a non-randomized / preexisting-groups design.
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
#' \strong{Why \code{rho} must be a single number when \code{randomized = TRUE}.}
#' Random assignment forms each group as an exchangeable random sample from
#' the same population. The bivariate distribution of \eqn{(Y, X)} is
#' therefore the same in every group, including the correlation. Allowing
#' \code{rho} to differ across groups would silently break that
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
#'   mu_y       = c(3, 5),
#'   mu_x       = 10,
#'   sigma_y    = 1,
#'   sigma_x    = 2,
#'   rho        = 0.8,
#'   a          = 2,
#'   n          = 20
#' )
#' head(simple)
#'
#' # 2. Four preexisting groups with different correlations and unequal n.
#' #    The first two groups share rho = 0.30; the second two share a larger
#' #    rho = 0.60. The four groups are not used in the data-generation
#' #    machinery beyond their per-group means and correlations -- they are
#' #    just four distinct populations being sampled. In a downstream
#' #    analysis these four groups could be cross-classified as a 2 x 2
#' #    factorial design (e.g., the first factor distinguishing groups 1-2
#' #    from groups 3-4, and the second factor distinguishing groups 1, 3
#' #    from groups 2, 4) and analyzed via factorial ANCOVA.
#' set.seed(113)
#' preexisting <- simulate_ancova_data(
#'   mu_y       = c(50, 55, 60, 65),
#'   mu_x       = c(10, 12, 11, 13),
#'   sigma_y    = 8,
#'   sigma_x    = 3,
#'   rho        = c(0.30, 0.30, 0.60, 0.60),
#'   a          = 4,
#'   n          = c(40, 35, 45, 30),
#'   randomized = FALSE
#' )
#' aggregate(cbind(y, x) ~ group, data = preexisting,
#'           FUN = function(z) round(c(mean = mean(z), sd = sd(z)), 2))
#'
#' @export
#' @importFrom MASS mvrnorm
simulate_ancova_data <- function(mu_y, mu_x, sigma_y, sigma_x, rho,
                                 a, n, randomized = TRUE) {

  # ---------- Validate scalar arguments ----------
  if (!is.numeric(a) || length(a) != 1L || a < 2L || a != round(a)) {
    stop("'a' must be a single integer >= 2 (the number of fixed groups).",
         call. = FALSE)
  }
  if (!is.numeric(sigma_y) || length(sigma_y) != 1L || sigma_y <= 0) {
    stop("'sigma_y' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(sigma_x) || length(sigma_x) != 1L || sigma_x <= 0) {
    stop("'sigma_x' must be a single positive number.", call. = FALSE)
  }
  if (!is.logical(randomized) || length(randomized) != 1L) {
    stop("'randomized' must be a single logical value.", call. = FALSE)
  }

  # ---------- Validate vector arguments ----------
  if (!is.numeric(mu_y) || length(mu_y) != a) {
    stop("'mu_y' must be a numeric vector of length 'a' (one population mean of Y per group).",
         call. = FALSE)
  }

  if (isTRUE(randomized)) {
    if (length(mu_x) != 1L) {
      stop("Under 'randomized = TRUE', 'mu_x' must be a single number ",
           "(the common population covariate mean across groups).",
           call. = FALSE)
    }
    if (length(rho) != 1L) {
      stop("Under 'randomized = TRUE', 'rho' must be a single number. ",
           "Random assignment forms each group as a sample from a common ",
           "population, so the Y-X correlation is by definition the same ",
           "in every group; allowing it to differ would contradict the ",
           "randomized-design assumption. Use 'randomized = FALSE' to ",
           "simulate group-specific correlations.", call. = FALSE)
    }
    mu_x_per_group  <- rep(as.numeric(mu_x),  a)
    rho_per_group   <- rep(as.numeric(rho),   a)
  } else {
    if (length(mu_x) != a) {
      stop("Under 'randomized = FALSE', 'mu_x' must be a numeric vector of length 'a' ",
           "(one population covariate mean per group).", call. = FALSE)
    }
    if (length(rho) == 1L) {
      rho_per_group <- rep(as.numeric(rho), a)
    } else if (length(rho) == a) {
      rho_per_group <- as.numeric(rho)
    } else {
      stop("Under 'randomized = FALSE', 'rho' must be either a single number ",
           "(common correlation, recycled across all groups) or a numeric ",
           "vector of length 'a' (one correlation per group).", call. = FALSE)
    }
    mu_x_per_group <- as.numeric(mu_x)
  }

  if (any(abs(rho_per_group) >= 1)) {
    stop("All correlations in 'rho' must lie strictly between -1 and 1.",
         call. = FALSE)
  }

  if (length(n) == 1L) {
    if (!is.numeric(n) || n < 2L || n != round(n)) {
      stop("'n' must be a positive integer >= 2.", call. = FALSE)
    }
    n_per_group <- rep(as.integer(n), a)
  } else if (length(n) == a) {
    if (!is.numeric(n) || any(n < 2L) || any(n != round(n))) {
      stop("Each entry of 'n' must be a positive integer >= 2.", call. = FALSE)
    }
    n_per_group <- as.integer(n)
  } else {
    stop("'n' must be a single number (equal n per group) or a numeric vector ",
         "of length 'a' (one sample size per group).", call. = FALSE)
  }

  # ---------- Generate the data ----------
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install with install.packages(\"MASS\").",
         call. = FALSE)
  }

  N_total  <- sum(n_per_group)
  group    <- factor(rep(seq_len(a), times = n_per_group),
                     levels = seq_len(a),
                     labels = as.character(seq_len(a)))
  y        <- numeric(N_total)
  x        <- numeric(N_total)

  pos <- 1L
  for (j in seq_len(a)) {
    n_j     <- n_per_group[j]
    cov_yx  <- rho_per_group[j] * sigma_y * sigma_x
    Sigma   <- matrix(c(sigma_y^2, cov_yx,
                        cov_yx,    sigma_x^2),
                      nrow = 2)
    sample_j <- MASS::mvrnorm(
      n_j,
      mu    = c(mu_y[j], mu_x_per_group[j]),
      Sigma = Sigma
    )
    idx        <- pos:(pos + n_j - 1L)
    y[idx]     <- sample_j[, 1L]
    x[idx]     <- sample_j[, 2L]
    pos        <- pos + n_j
  }

  data.frame(group = group, y = y, x = x, stringsAsFactors = FALSE)
}
