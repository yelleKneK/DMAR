#' Simulate Data From a Multivariate Normal Multiple-Regression Model
#'
#' Generates random data \eqn{(Y, X_1, \ldots, X_p)} jointly multivariate
#' normal with user-specified marginal means, marginal SDs, and full
#' correlation structure. Useful as a backbone for sensitivity analyses,
#' Monte Carlo studies of regression sample size methods, and pedagogical
#' demonstrations.
#'
#' @param N The total sample size (a positive integer \eqn{\ge p + 2}).
#' @param p The number of predictor variables.
#' @param rho_YX A numeric vector of length \code{p} giving the population
#'   correlations between \eqn{Y} and each predictor \eqn{X_j}.
#' @param rho_XX A \eqn{p \times p} symmetric correlation matrix for the
#'   predictors. Defaults to the identity matrix (orthogonal predictors).
#' @param mu_Y The population mean of \eqn{Y} (default \code{0}).
#' @param mu_X A numeric vector of length \code{p} giving the population
#'   mean of each predictor (default \code{0}, recycled across predictors).
#' @param sigma_Y The population standard deviation of \eqn{Y} (default
#'   \code{1}, in which case the simulated \eqn{Y} is on the standardized
#'   scale).
#' @param sigma_X A single number or a numeric vector of length \code{p}
#'   giving the population standard deviation of each predictor (default
#'   \code{1}, standardized).
#' @param seed Optional integer random seed for reproducibility (default
#'   \code{NULL}).
#' @param column_names Optional character vector of length \code{p + 1}
#'   giving column names for the returned \code{data.frame}; defaults to
#'   \code{c("y", "x1", "x2", ..., "xp")}.
#'
#' @return A \code{data.frame} with \code{N} rows and \code{p + 1} columns:
#'   the outcome \eqn{Y} (first column) followed by predictors
#'   \eqn{X_1, \ldots, X_p}.
#'
#' @details
#' Internally the joint correlation matrix is assembled as
#' \deqn{R = \begin{pmatrix} 1 & \rho_{YX}^\top \\ \rho_{YX} & R_{XX} \end{pmatrix},}
#' converted to a covariance matrix via the supplied SDs, and \code{N}
#' draws are taken using \code{\link[MASS]{mvrnorm}}. The resulting
#' \eqn{Y} and predictors satisfy the requested marginal means and
#' standard deviations and (in expectation) the requested correlation
#' structure.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{simulate_ancova_data}},
#'   \code{\link{simulate_anova_data}}, \code{\link{ss_aipe_R2}},
#'   \code{\link{ss_aipe_reg_coef}}, \code{\link{ci_R2}}
#'
#' @family data simulators
#'
#' @keywords design datagen
#'
#' @examples
#' # Five orthogonal predictors, each correlating .30 with Y.
#' set.seed(113)
#' d <- simulate_regression_data(
#'   N      = 200,
#'   p      = 5,
#'   rho_YX = rep(0.30, 5)
#' )
#' summary(lm(y ~ ., data = d))$r.squared   # ~ 5 * 0.30^2 = 0.45
#'
#' # Predictors with shared structure (exchangeable correlation matrix).
#' rho_XX <- matrix(0.5, nrow = 5, ncol = 5); diag(rho_XX) <- 1
#' simulate_regression_data(
#'   N      = 300,
#'   p      = 5,
#'   rho_YX = c(.50, .40, .30, .20, .10),
#'   rho_XX = rho_XX,
#'   seed   = 113
#' )[1:3, ]
#'
#' @export
#' @importFrom MASS mvrnorm
simulate_regression_data <- function(N, p, rho_YX,
                                     rho_XX       = NULL,
                                     mu_Y         = 0,
                                     mu_X         = 0,
                                     sigma_Y      = 1,
                                     sigma_X      = 1,
                                     seed         = NULL,
                                     column_names = NULL) {

  if (!is.numeric(N) || length(N) != 1L || N != round(N) || N < 2L) {
    stop("'N' must be a positive integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1L || p != round(p) || p < 1L) {
    stop("'p' must be a positive integer.", call. = FALSE)
  }
  if (N <= p + 1L) {
    stop("'N' must be > p + 1 (and ideally much larger).", call. = FALSE)
  }
  if (!is.numeric(rho_YX) || length(rho_YX) != p) {
    stop("'rho_YX' must be a numeric vector of length 'p'.", call. = FALSE)
  }
  if (any(abs(rho_YX) >= 1)) {
    stop("All entries of 'rho_YX' must lie strictly in (-1, 1).",
         call. = FALSE)
  }

  if (is.null(rho_XX)) {
    rho_XX <- diag(p)
  } else {
    rho_XX <- as.matrix(rho_XX)
    if (!isTRUE(all.equal(dim(rho_XX), c(p, p)))) {
      stop("'rho_XX' must be a p x p numeric matrix.", call. = FALSE)
    }
    if (!isTRUE(all.equal(rho_XX, t(rho_XX), tolerance = 1e-8))) {
      stop("'rho_XX' must be symmetric.", call. = FALSE)
    }
    if (any(abs(rho_XX[lower.tri(rho_XX)]) >= 1)) {
      stop("All off-diagonal entries of 'rho_XX' must be in (-1, 1).",
           call. = FALSE)
    }
  }

  # mu_X / sigma_X: scalar -> recycle to length p.
  if (length(mu_X)    == 1L) mu_X    <- rep(as.numeric(mu_X),    p)
  if (length(sigma_X) == 1L) sigma_X <- rep(as.numeric(sigma_X), p)
  if (length(mu_X)    != p) stop("'mu_X' must have length 1 or 'p'.", call. = FALSE)
  if (length(sigma_X) != p) stop("'sigma_X' must have length 1 or 'p'.", call. = FALSE)
  if (any(sigma_X <= 0))    stop("All entries of 'sigma_X' must be positive.", call. = FALSE)
  if (sigma_Y <= 0)         stop("'sigma_Y' must be positive.", call. = FALSE)

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install with install.packages(\"MASS\").",
         call. = FALSE)
  }

  # Assemble joint correlation matrix and convert to covariance matrix.
  R <- rbind(c(1, rho_YX), cbind(rho_YX, rho_XX))
  sd_vec <- c(sigma_Y, sigma_X)
  Sigma  <- diag(sd_vec) %*% R %*% diag(sd_vec)
  mu_vec <- c(mu_Y, mu_X)

  # Verify the joint matrix is positive definite, since combining a
  # user-supplied rho_YX with rho_XX can produce inconsistent structures.
  eigs <- eigen(Sigma, only.values = TRUE)$values
  if (min(eigs) <= 0) {
    stop("The joint covariance matrix implied by 'rho_YX' and 'rho_XX' is ",
         "not positive definite. Check that the requested correlation ",
         "structure is internally consistent.", call. = FALSE)
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
  M <- MASS::mvrnorm(N, mu = mu_vec, Sigma = Sigma)

  if (is.null(column_names)) {
    column_names <- c("y", paste0("x", seq_len(p)))
  } else if (length(column_names) != p + 1L) {
    stop("'column_names' must have length p + 1 (Y followed by the p X variables).",
         call. = FALSE)
  }
  out <- as.data.frame(M)
  names(out) <- column_names
  out
}
