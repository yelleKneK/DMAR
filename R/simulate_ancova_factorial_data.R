# Internal: assemble cell labels in a stable lexicographic order across
# factors A, B, C, D. expand.grid varies the *first* argument fastest, which
# matches what most readers expect for textbook factorial-design tables.
.factorial_cells <- function(a, b, c, d) {
  if (b == 1L && c == 1L && d == 1L) {
    cells <- data.frame(A = factor(seq_len(a)), stringsAsFactors = FALSE)
  } else if (c == 1L && d == 1L) {
    cells <- expand.grid(A = factor(seq_len(a)),
                         B = factor(seq_len(b)),
                         stringsAsFactors = FALSE)
  } else if (d == 1L) {
    cells <- expand.grid(A = factor(seq_len(a)),
                         B = factor(seq_len(b)),
                         C = factor(seq_len(c)),
                         stringsAsFactors = FALSE)
  } else {
    cells <- expand.grid(A = factor(seq_len(a)),
                         B = factor(seq_len(b)),
                         C = factor(seq_len(c)),
                         D = factor(seq_len(d)),
                         stringsAsFactors = FALSE)
  }
  cells
}


#' Simulate Data From a Factorial ANCOVA Design (up to Four Factors, Any
#' Number of Covariates)
#'
#' Generates random data appropriate for an analysis of covariance with up
#' to four crossed fixed factors (\eqn{A, B, C, D}) and one or more
#' continuous covariates. Within every cell the outcome \eqn{Y} and the
#' covariates \eqn{X_1, \ldots, X_q} are jointly multivariate normal with a
#' common (homogeneous) within-cell covariance structure, the standard
#' assumption underlying classical ANCOVA. Per-cell sample sizes may be
#' equal or unequal.
#'
#' @param a Number of levels of the first factor (must be at least \code{2}).
#' @param b Number of levels of the second factor (default \code{1}; use \code{1} if absent).
#' @param c Number of levels of the third factor (default \code{1}; use \code{1} if absent).
#' @param d Number of levels of the fourth factor (default \code{1}; use \code{1} if absent).
#'   So \code{a = 3, b = 2, c = 1, d = 1} specifies a 3 \eqn{\times} 2 design.
#' @param n_covariates Integer \eqn{\ge 1}: the number of continuous
#'   covariates. Default \code{1}.
#' @param mu_y Numeric vector of length \eqn{a \cdot b \cdot c \cdot d}
#'   giving the population cell means of \eqn{Y}, in the same row order as
#'   \code{\link[base]{expand.grid}(A, B, C, D)} (which varies factor
#'   \code{A} fastest, then \code{B}, then \code{C}, then \code{D}).
#' @param mu_x Numeric matrix of dimension
#'   \eqn{(a \cdot b \cdot c \cdot d) \times q} giving the population cell
#'   means of each covariate (rows = cells in the same order as
#'   \code{mu_y}; columns = covariates). When \code{randomized = TRUE}
#'   every column must be constant across cells.
#' @param sigma_y Within-cell standard deviation of \eqn{Y} (a single
#'   positive number, common across cells).
#' @param sigma_x Within-cell standard deviations of the covariates.
#'   Either a single number (recycled to all \code{n_covariates}) or a
#'   numeric vector of length \code{n_covariates}.
#' @param rho_y_x Within-cell correlations between \eqn{Y} and each
#'   covariate. Either a single number (recycled to all
#'   \code{n_covariates}) or a numeric vector of length
#'   \code{n_covariates}.
#' @param rho_x_x Within-cell correlation matrix among the covariates,
#'   \eqn{q \times q}. Default \code{NULL}, interpreted as the
#'   \code{n_covariates}-dimensional identity.
#' @param n A single number (equal sample size per cell) or a numeric
#'   vector of length \eqn{a \cdot b \cdot c \cdot d} giving per-cell
#'   sample sizes (in the same row order as \code{mu_y}).
#' @param randomized Logical. \code{TRUE} (the default) for a randomized
#'   design, each cell is sampled from the same population covariate
#'   distribution, so \code{mu_x} must be constant across cells. \code{FALSE}
#'   for non-randomized / preexisting-groups designs in which the cells'
#'   population covariate means may differ.
#'
#' @return A long-format \code{data.frame} with one row per simulated
#'   subject and the following columns:
#'   \describe{
#'     \item{\code{A}, \code{B}, \code{C}, \code{D}}{Factor columns for
#'       each present design factor (omitted when the factor is absent,
#'       i.e., its corresponding \code{a/b/c/d} argument is \code{1}).}
#'     \item{\code{x1}, \code{x2}, ..., \code{x<q>}}{Numeric simulated
#'       covariates.}
#'     \item{\code{y}}{Numeric simulated outcome.}
#'   }
#'
#' @details
#' This is the factorial generalization of \code{\link{simulate_ancova_data}}
#' (which is the special case \code{b = c = d = 1, n_covariates = 1}). All
#' cells share the same within-cell covariance structure (homogeneity of
#' regression, the classical ANCOVA assumption); the difference between
#' randomized and non-randomized designs lies entirely in whether the cell
#' covariate means are constrained to be equal.
#'
#' \strong{Why \code{mu_x} must be constant across cells when
#' \code{randomized = TRUE}.} Random assignment forms each cell as an
#' exchangeable random sample from the same joint distribution of
#' covariates and outcome. Cell-specific covariate means would silently
#' break that interpretation. The function checks the constraint and
#' stops with an informative error if it is violated. (The same logic
#' that powers \code{\link{simulate_ancova_data}}.)
#'
#' \strong{Cell ordering.} The function uses
#' \code{\link[base]{expand.grid}}'s convention, factor \code{A} varies
#' fastest, then \code{B}, then \code{C}, then \code{D}. So for a
#' \eqn{2 \times 3} design, the six cells of \code{mu_y} are
#' \eqn{(A_1 B_1), (A_2 B_1), (A_1 B_2), (A_2 B_2), (A_1 B_3), (A_2 B_3)}.
#' If you build the cell specification by passing the factor levels to
#' \code{expand.grid} in the same order, the row indexing automatically
#' matches.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{simulate_ancova_data}} (one factor, one covariate
#'   special case), \code{\link{simulate_anova_data}},
#'   \code{\link{simulate_regression_data}}
#'
#' @family data simulators
#'
#' @keywords design datagen
#'
#' @examples
#' # 1. 2 x 2 randomized design, single covariate.
#' set.seed(113)
#' design_2x2 <- expand.grid(A = factor(1:2), B = factor(1:2))
#' design_2x2$mu_y <- c(50, 60, 55, 65)   # cell means in expand.grid order
#'
#' d1 <- simulate_ancova_factorial_data(
#'   a            = 2, b = 2,
#'   mu_y         = design_2x2$mu_y,
#'   mu_x         = matrix(10, nrow = 4, ncol = 1),  # constant covariate mean
#'   sigma_y      = 8,
#'   sigma_x      = 3,
#'   rho_y_x      = 0.40,
#'   n            = 30
#' )
#' aggregate(y ~ A + B, data = d1, FUN = mean)
#'
#' # 2. 3 x 2 nonrandomized design with two covariates and unequal n.
#' set.seed(113)
#' a <- 3; b <- 2; q <- 2
#' n_cells <- a * b
#' d2 <- simulate_ancova_factorial_data(
#'   a            = a, b = b,
#'   n_covariates = q,
#'   mu_y         = c(50, 55, 60,  52, 58, 64),
#'   mu_x         = matrix(c(10, 11, 12,  9, 10, 11,
#'                            5,  6,  7,  4,  5,  6),
#'                          nrow = n_cells, ncol = q),
#'   sigma_y      = 8,
#'   sigma_x      = c(3, 2),
#'   rho_y_x      = c(0.40, 0.25),
#'   rho_x_x      = matrix(c(1,   0.3,
#'                           0.3, 1),
#'                         nrow = 2),
#'   n            = c(30, 25, 20,  35, 30, 25),
#'   randomized   = FALSE
#' )
#' head(d2)
#'
#' # 3. 2 x 2 x 2 randomized design with one covariate.
#' set.seed(113)
#' d3 <- simulate_ancova_factorial_data(
#'   a       = 2, b = 2, c = 2,
#'   mu_y    = c(50, 55, 52, 57,  53, 58, 55, 60),  # 8 cells in A-fastest order
#'   mu_x    = matrix(10, nrow = 8, ncol = 1),
#'   sigma_y = 8,
#'   sigma_x = 3,
#'   rho_y_x = 0.40,
#'   n       = 25
#' )
#' table(d3$A, d3$B, d3$C)  # 25 per cell, 200 total
#'
#' @export
#' @importFrom MASS mvrnorm
simulate_ancova_factorial_data <- function(a, b = 1, c = 1, d = 1,
                                           n_covariates = 1,
                                           mu_y, mu_x,
                                           sigma_y, sigma_x, rho_y_x,
                                           rho_x_x = NULL,
                                           n,
                                           randomized = TRUE) {

  # ---------- Validate factor counts ----------
  for (nm in c("a", "b", "c", "d")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || val < 1L || val != round(val)) {
      stop("'", nm, "' must be a single positive integer.", call. = FALSE)
    }
  }
  if (a < 2L) stop("'a' must be at least 2.", call. = FALSE)

  # ---------- Validate covariate count ----------
  if (!is.numeric(n_covariates) || length(n_covariates) != 1L ||
      n_covariates < 1L || n_covariates != round(n_covariates)) {
    stop("'n_covariates' must be a single positive integer.", call. = FALSE)
  }
  q <- as.integer(n_covariates)

  # ---------- Cell layout ----------
  cells <- .factorial_cells(a, b, c, d)
  K <- nrow(cells)
  if (K != a * b * c * d) {
    stop("Internal error: cell count mismatch.", call. = FALSE)
  }

  # ---------- Validate mu_y ----------
  if (!is.numeric(mu_y) || length(mu_y) != K) {
    stop("'mu_y' must be a numeric vector of length a * b * c * d (= ", K, ").",
         call. = FALSE)
  }

  # ---------- Validate mu_x ----------
  mu_x <- as.matrix(mu_x)
  if (!is.numeric(mu_x) || nrow(mu_x) != K || ncol(mu_x) != q) {
    stop("'mu_x' must be a numeric matrix with ", K, " rows (cells) and ",
         q, " columns (covariates).", call. = FALSE)
  }

  if (isTRUE(randomized)) {
    constant_per_col <- apply(mu_x, 2, function(col) length(unique(col)) == 1L)
    if (!all(constant_per_col)) {
      stop("Under 'randomized = TRUE', every column of 'mu_x' must be constant ",
           "across cells. Random assignment forms each cell as a sample from a ",
           "common joint distribution of covariates and outcome, so the cell ",
           "covariate means cannot differ. Use 'randomized = FALSE' for a ",
           "non-randomized / preexisting-groups design.", call. = FALSE)
    }
  }

  # ---------- Validate sigma_y, sigma_x ----------
  if (!is.numeric(sigma_y) || length(sigma_y) != 1L || sigma_y <= 0) {
    stop("'sigma_y' must be a single positive number.", call. = FALSE)
  }
  if (length(sigma_x) == 1L) sigma_x <- rep(as.numeric(sigma_x), q)
  if (length(sigma_x) != q || any(sigma_x <= 0)) {
    stop("'sigma_x' must be a single positive number or a numeric vector of ",
         "length 'n_covariates' of positive numbers.", call. = FALSE)
  }

  # ---------- Validate rho_y_x ----------
  if (length(rho_y_x) == 1L) rho_y_x <- rep(as.numeric(rho_y_x), q)
  if (length(rho_y_x) != q || any(abs(rho_y_x) >= 1)) {
    stop("'rho_y_x' must be a single number or a numeric vector of length ",
         "'n_covariates', all entries strictly in (-1, 1).", call. = FALSE)
  }

  # ---------- Validate rho_x_x ----------
  if (is.null(rho_x_x)) {
    rho_x_x <- diag(q)
  } else {
    rho_x_x <- as.matrix(rho_x_x)
    if (!isTRUE(all.equal(dim(rho_x_x), c(q, q)))) {
      stop("'rho_x_x' must be a ", q, " x ", q, " correlation matrix.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(rho_x_x, t(rho_x_x), tolerance = 1e-8))) {
      stop("'rho_x_x' must be symmetric.", call. = FALSE)
    }
    if (any(abs(rho_x_x[lower.tri(rho_x_x)]) >= 1)) {
      stop("Off-diagonal entries of 'rho_x_x' must be in (-1, 1).",
           call. = FALSE)
    }
  }

  # ---------- Build joint within-cell covariance ----------
  # First row/column is Y; remaining rows/columns are X1..Xq.
  R <- rbind(c(1, rho_y_x), cbind(rho_y_x, rho_x_x))
  sd_vec <- c(sigma_y, sigma_x)
  Sigma  <- diag(sd_vec) %*% R %*% diag(sd_vec)

  if (min(eigen(Sigma, only.values = TRUE)$values) <= 0) {
    stop("The implied within-cell covariance matrix is not positive definite. ",
         "Check that 'rho_y_x' and 'rho_x_x' are jointly consistent.",
         call. = FALSE)
  }

  # ---------- Validate n ----------
  if (length(n) == 1L) {
    if (!is.numeric(n) || n < 2L || n != round(n)) {
      stop("'n' must be a positive integer >= 2.", call. = FALSE)
    }
    n_per_cell <- rep(as.integer(n), K)
  } else if (length(n) == K) {
    if (!is.numeric(n) || any(n < 2L) || any(n != round(n))) {
      stop("Each entry of 'n' must be a positive integer >= 2.", call. = FALSE)
    }
    n_per_cell <- as.integer(n)
  } else {
    stop("'n' must be a single number or a numeric vector of length ", K, ".",
         call. = FALSE)
  }

  # ---------- Generate ----------
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.", call. = FALSE)
  }

  N_total <- sum(n_per_cell)
  out_y   <- numeric(N_total)
  out_x   <- matrix(NA_real_, nrow = N_total, ncol = q)
  cell_id <- integer(N_total)
  pos     <- 1L

  for (k in seq_len(K)) {
    n_k <- n_per_cell[k]
    mu_k <- c(mu_y[k], mu_x[k, ])
    draw <- MASS::mvrnorm(n_k, mu = mu_k, Sigma = Sigma)
    if (n_k == 1L) draw <- matrix(draw, nrow = 1L)
    idx <- pos:(pos + n_k - 1L)
    out_y[idx]    <- draw[, 1L]
    out_x[idx, ]  <- draw[, -1L, drop = FALSE]
    cell_id[idx]  <- k
    pos           <- pos + n_k
  }

  # ---------- Assemble long-format data.frame ----------
  long <- cells[cell_id, , drop = FALSE]
  rownames(long) <- NULL
  for (j in seq_len(q)) {
    long[[paste0("x", j)]] <- out_x[, j]
  }
  long$y <- out_y
  long
}
