#' Simulate Data From a Multivariate Normal Multiple-Regression Model
#'
#' Generates random data \eqn{(Y, X_1, \ldots, X_p)} jointly multivariate
#' normal with user-specified marginal means, marginal SDs, and full
#' correlation structure. Useful for sensitivity analyses, Monte Carlo
#' studies of regression sample size methods, and teaching demonstrations.
#'
#' @param N The total sample size (a positive integer \eqn{\ge p + 2}).
#' @param p The number of predictor variables.
#' @param true_cor_YX A numeric vector of length \code{p} giving the population
#'   correlations between \eqn{Y} and each predictor \eqn{X_j}.
#' @param true_cor_XX A \eqn{p \times p} symmetric correlation matrix for the
#'   predictors. Defaults to the identity matrix (orthogonal predictors).
#' @param true_mean_Y The population mean of \eqn{Y} (default \code{0}).
#' @param true_mean_X A numeric vector of length \code{p} giving the population
#'   mean of each predictor (default \code{0}, recycled across predictors).
#' @param true_sd_Y The population standard deviation of \eqn{Y} (default
#'   \code{1}, in which case the simulated \eqn{Y} is on the standardized
#'   scale).
#' @param true_sd_X A single number or a numeric vector of length \code{p}
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
#' where \eqn{\rho_{YX}} is \code{true_cor_YX} and \eqn{R_{XX}} is
#' \code{true_cor_XX},
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
#'   true_cor_YX = rep(0.30, 5)
#' )
#' summary(lm(y ~ ., data = d))$r.squared   # about 0.45: five predictors, each 0.30^2
#'
#' # Predictors with shared structure (exchangeable correlation matrix).
#' true_cor_XX <- matrix(0.5, nrow = 5, ncol = 5); diag(true_cor_XX) <- 1
#' simulate_regression_data(
#'   N      = 300,
#'   p      = 5,
#'   true_cor_YX = c(.50, .40, .30, .20, .10),
#'   true_cor_XX = true_cor_XX,
#'   seed   = 113
#' )[1:3, ]
#'
#' @export
#' @importFrom MASS mvrnorm
simulate_regression_data <- function(N, p, true_cor_YX,
                                     true_cor_XX  = NULL,
                                     true_mean_Y  = 0,
                                     true_mean_X  = 0,
                                     true_sd_Y    = 1,
                                     true_sd_X    = 1,
                                     seed         = NULL,
                                     column_names = NULL) {

  if (!is.numeric(N) || length(N) != 1L || N != round(N) || N < 2L) {
    .dmar_refuse_argument("'N' must be a positive integer >= 2.", arg = "N")
  }
  if (!is.numeric(p) || length(p) != 1L || p != round(p) || p < 1L) {
    .dmar_refuse_argument("'p' must be a positive integer.", arg = "p")
  }
  if (N <= p + 1L) {
    .dmar_refuse_argument(
      "'N' must be > p + 1 (and ideally much larger).", arg = "N")
  }
  if (!is.numeric(true_cor_YX) || length(true_cor_YX) != p) {
    .dmar_refuse_argument(
      "'true_cor_YX' must be a numeric vector of length 'p'.",
      arg = c("true_cor_YX", "p"))
  }
  if (any(abs(true_cor_YX) >= 1)) {
    .dmar_refuse_argument(
      "All entries of 'true_cor_YX' must lie strictly in (-1, 1).", arg = "true_cor_YX")
  }

  if (is.null(true_cor_XX)) {
    true_cor_XX <- diag(p)
  } else {
    true_cor_XX <- as.matrix(true_cor_XX)
    if (!isTRUE(all.equal(dim(true_cor_XX), c(p, p)))) {
      .dmar_refuse_argument(
        "'true_cor_XX' must be a p x p numeric matrix.", arg = "true_cor_XX")
    }
    if (!isTRUE(all.equal(true_cor_XX, t(true_cor_XX), tolerance = 1e-8))) {
      .dmar_refuse_argument("'true_cor_XX' must be symmetric.", arg = "true_cor_XX")
    }
    if (any(abs(true_cor_XX[lower.tri(true_cor_XX)]) >= 1)) {
      .dmar_refuse_argument(
        "All off-diagonal entries of 'true_cor_XX' must be in (-1, 1).",
        arg = "true_cor_XX")
    }
  }

  # true_mean_X / true_sd_X: scalar -> recycle to length p.
  if (length(true_mean_X) == 1L) true_mean_X <- rep(as.numeric(true_mean_X), p)
  if (length(true_sd_X)   == 1L) true_sd_X   <- rep(as.numeric(true_sd_X),   p)
  if (length(true_mean_X) != p) .dmar_refuse_argument(
    "'true_mean_X' must have length 1 or 'p'.", arg = c("true_mean_X", "p"))
  if (length(true_sd_X) != p) .dmar_refuse_argument(
    "'true_sd_X' must have length 1 or 'p'.", arg = c("true_sd_X", "p"))
  if (any(true_sd_X <= 0))    .dmar_refuse_argument(
    "All entries of 'true_sd_X' must be positive.", arg = "true_sd_X")
  if (true_sd_Y <= 0)         .dmar_refuse_argument(
    "'true_sd_Y' must be positive.", arg = "true_sd_Y")

  if (!requireNamespace("MASS", quietly = TRUE)) {
    .dmar_refuse(
      "Package 'MASS' is required. Install with install.packages(\"MASS\").",
      "package_required", package = "MASS")
  }

  # Assemble joint correlation matrix and convert to covariance matrix.
  R <- rbind(c(1, true_cor_YX), cbind(true_cor_YX, true_cor_XX))
  sd_vec <- c(true_sd_Y, true_sd_X)
  Sigma  <- diag(sd_vec) %*% R %*% diag(sd_vec)
  mu_vec <- c(true_mean_Y, true_mean_X)

  # Verify the joint matrix is positive definite, since combining a
  # user-supplied true_cor_YX with true_cor_XX can produce inconsistent structures.
  eigs <- eigen(Sigma, only.values = TRUE)$values
  if (min(eigs) <= 0) {
    .dmar_refuse_argument(c(
      "The joint covariance matrix implied by 'true_cor_YX' and 'true_cor_XX' is ",
      "not positive definite. Check that the requested correlation ",
      "structure is internally consistent."), arg = c("true_cor_YX", "true_cor_XX"))
  }

  .dmar_local_seed(seed)
  M <- MASS::mvrnorm(N, mu = mu_vec, Sigma = Sigma)

  if (is.null(column_names)) {
    column_names <- c("y", paste0("x", seq_len(p)))
  } else if (length(column_names) != p + 1L) {
    .dmar_refuse_argument(
      "'column_names' must have length p + 1 (Y followed by the p X variables).",
      arg = "column_names")
  }
  out <- as.data.frame(M)
  names(out) <- column_names
  out
}
