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
#' @param true_mean_Y Numeric vector of length \eqn{a \cdot b \cdot c \cdot d}
#'   giving the population cell means of \eqn{Y}, in the same row order as
#'   \code{\link[base]{expand.grid}(A, B, C, D)} (which varies factor
#'   \code{A} fastest, then \code{B}, then \code{C}, then \code{D}).
#' @param true_mean_X Numeric matrix of dimension
#'   \eqn{(a \cdot b \cdot c \cdot d) \times q} giving the population cell
#'   means of each covariate (rows = cells in the same order as
#'   \code{true_mean_Y}; columns = covariates). When \code{randomized = TRUE}
#'   every column must be constant across cells.
#' @param true_sd_Y Within-cell standard deviation of \eqn{Y} (a single
#'   positive number, common across cells).
#' @param true_sd_X Within-cell standard deviations of the covariates.
#'   Either a single number (recycled to all \code{n_covariates}) or a
#'   numeric vector of length \code{n_covariates}.
#' @param true_cor_YX Within-cell correlations between \eqn{Y} and each
#'   covariate. Either a single number (recycled to all
#'   \code{n_covariates}) or a numeric vector of length
#'   \code{n_covariates}.
#' @param true_cor_XX Within-cell correlation matrix among the covariates,
#'   \eqn{q \times q}. Default \code{NULL}, interpreted as the
#'   \code{n_covariates}-dimensional identity.
#' @param n_per_cell A single number (equal sample size per cell) or a numeric
#'   vector of length \eqn{a \cdot b \cdot c \cdot d} giving per-cell
#'   sample sizes (in the same row order as \code{true_mean_Y}).
#' @param randomized Logical. \code{TRUE} (the default) for a randomized
#'   design, each cell is sampled from the same population covariate
#'   distribution, so \code{true_mean_X} must be constant across cells. \code{FALSE}
#'   for non-randomized / preexisting-groups designs in which the cells'
#'   population covariate means may differ.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
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
#' \strong{Why \code{true_mean_X} must be constant across cells when
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
#' \eqn{2 \times 3} design, the six cells of \code{true_mean_Y} are
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
#' design_2x2$true_mean_Y <- c(50, 60, 55, 65)   # cell means in expand.grid order
#'
#' d1 <- simulate_ancova_factorial_data(
#'   a            = 2, b = 2,
#'   true_mean_Y  = design_2x2$true_mean_Y,
#'   true_mean_X  = matrix(10, nrow = 4, ncol = 1),  # constant covariate mean
#'   true_sd_Y      = 8,
#'   true_sd_X      = 3,
#'   true_cor_YX  = 0.40,
#'   n_per_cell   = 30
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
#'   true_mean_Y  = c(50, 55, 60,  52, 58, 64),
#'   true_mean_X  = matrix(c(10, 11, 12,  9, 10, 11,
#'                            5,  6,  7,  4,  5,  6),
#'                          nrow = n_cells, ncol = q),
#'   true_sd_Y      = 8,
#'   true_sd_X      = c(3, 2),
#'   true_cor_YX  = c(0.40, 0.25),
#'   true_cor_XX  = matrix(c(1,   0.3,
#'                           0.3, 1),
#'                         nrow = 2),
#'   n_per_cell   = c(30, 25, 20,  35, 30, 25),
#'   randomized   = FALSE
#' )
#' head(d2)
#'
#' # 3. 2 x 2 x 2 randomized design with one covariate.
#' set.seed(113)
#' d3 <- simulate_ancova_factorial_data(
#'   a       = 2, b = 2, c = 2,
#'   true_mean_Y = c(50, 55, 52, 57,  53, 58, 55, 60),  # 8 cells in A-fastest order
#'   true_mean_X = matrix(10, nrow = 8, ncol = 1),
#'   true_sd_Y = 8,
#'   true_sd_X = 3,
#'   true_cor_YX = 0.40,
#'   n_per_cell = 25
#' )
#' table(d3$A, d3$B, d3$C)  # 25 per cell, 200 total
#'
#' @export
#' @importFrom MASS mvrnorm
simulate_ancova_factorial_data <- function(a, b = 1, c = 1, d = 1,
                                           n_covariates = 1,
                                           true_mean_Y, true_mean_X,
                                           true_sd_Y, true_sd_X, true_cor_YX,
                                           true_cor_XX = NULL,
                                           n_per_cell,
                                           randomized = TRUE,
                                           seed = NULL) {

  # ---------- Validate factor counts ----------
  for (nm in c("a", "b", "c", "d")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || val < 1L || val != round(val)) {
      .dmar_refuse_argument(c(
        "'", nm, "' must be a single positive integer."), arg = NULL)
    }
  }
  if (a < 2L) .dmar_refuse_argument("'a' must be at least 2.", arg = "a")

  # ---------- Validate covariate count ----------
  if (!is.numeric(n_covariates) || length(n_covariates) != 1L ||
      n_covariates < 1L || n_covariates != round(n_covariates)) {
    .dmar_refuse_argument(
      "'n_covariates' must be a single positive integer.", arg = "n_covariates")
  }
  q <- as.integer(n_covariates)

  # ---------- Cell layout ----------
  cells <- .factorial_cells(a, b, c, d)
  K <- nrow(cells)
  if (K != a * b * c * d) {
    .dmar_refuse("Internal error: cell count mismatch.", "internal_error")
  }

  # ---------- Validate true_mean_Y ----------
  if (!is.numeric(true_mean_Y) || length(true_mean_Y) != K) {
    .dmar_refuse_argument(c(
      "'true_mean_Y' must be a numeric vector of length a * b * c * d (= ", K, ")."),
      arg = "true_mean_Y")
  }

  # ---------- Validate true_mean_X ----------
  true_mean_X <- as.matrix(true_mean_X)
  if (!is.numeric(true_mean_X) || nrow(true_mean_X) != K || ncol(true_mean_X) != q) {
    .dmar_refuse_argument(c(
      "'true_mean_X' must be a numeric matrix with ", K, " rows (cells) and ",
      q, " columns (covariates)."), arg = "true_mean_X")
  }

  if (isTRUE(randomized)) {
    constant_per_col <- apply(true_mean_X, 2, function(col) length(unique(col)) == 1L)
    if (!all(constant_per_col)) {
      .dmar_refuse_argument(c(
        "Under 'randomized = TRUE', every column of 'true_mean_X' must be constant ",
        "across cells. Random assignment forms each cell as a sample from a ",
        "common joint distribution of covariates and outcome, so the cell ",
        "covariate means cannot differ. Use 'randomized = FALSE' for a ",
        "non-randomized / preexisting-groups design."), arg = "true_mean_X")
    }
  }

  # ---------- Validate true_sd_Y, true_sd_X ----------
  if (!is.numeric(true_sd_Y) || length(true_sd_Y) != 1L || true_sd_Y <= 0) {
    .dmar_refuse_argument(
      "'true_sd_Y' must be a single positive number.", arg = "true_sd_Y")
  }
  if (length(true_sd_X) == 1L) true_sd_X <- rep(as.numeric(true_sd_X), q)
  if (length(true_sd_X) != q || any(true_sd_X <= 0)) {
    .dmar_refuse_argument(c(
      "'true_sd_X' must be a single positive number or a numeric vector of ",
      "length 'n_covariates' of positive numbers."),
      arg = c("true_sd_X", "n_covariates"))
  }

  # ---------- Validate true_cor_YX ----------
  if (length(true_cor_YX) == 1L) true_cor_YX <- rep(as.numeric(true_cor_YX), q)
  if (length(true_cor_YX) != q || any(abs(true_cor_YX) >= 1)) {
    .dmar_refuse_argument(c(
      "'true_cor_YX' must be a single number or a numeric vector of length ",
      "'n_covariates', all entries strictly in (-1, 1)."),
      arg = c("true_cor_YX", "n_covariates"))
  }

  # ---------- Validate true_cor_XX ----------
  if (is.null(true_cor_XX)) {
    true_cor_XX <- diag(q)
  } else {
    true_cor_XX <- as.matrix(true_cor_XX)
    if (!isTRUE(all.equal(dim(true_cor_XX), c(q, q)))) {
      .dmar_refuse_argument(c(
        "'true_cor_XX' must be a ", q, " x ", q, " correlation matrix."),
        arg = "true_cor_XX")
    }
    if (!isTRUE(all.equal(true_cor_XX, t(true_cor_XX), tolerance = 1e-8))) {
      .dmar_refuse_argument("'true_cor_XX' must be symmetric.", arg = "true_cor_XX")
    }
    if (any(abs(true_cor_XX[lower.tri(true_cor_XX)]) >= 1)) {
      .dmar_refuse_argument(
        "Off-diagonal entries of 'true_cor_XX' must be in (-1, 1).",
        arg = "true_cor_XX")
    }
  }

  # ---------- Build joint within-cell covariance ----------
  # First row/column is Y; remaining rows/columns are X1..Xq.
  R <- rbind(c(1, true_cor_YX), cbind(true_cor_YX, true_cor_XX))
  sd_vec <- c(true_sd_Y, true_sd_X)
  Sigma  <- diag(sd_vec) %*% R %*% diag(sd_vec)

  if (min(eigen(Sigma, only.values = TRUE)$values) <= 0) {
    .dmar_refuse_argument(c(
      "The implied within-cell covariance matrix is not positive definite. ",
      "Check that 'true_cor_YX' and 'true_cor_XX' are jointly consistent."),
      arg = c("true_cor_YX", "true_cor_XX"))
  }

  # ---------- Validate n_per_cell ----------
  if (length(n_per_cell) == 1L) {
    if (!is.numeric(n_per_cell) || n_per_cell < 2L || n_per_cell != round(n_per_cell)) {
      .dmar_refuse_argument(
        "'n_per_cell' must be a positive integer >= 2.", arg = "n_per_cell")
    }
    n_per_cell <- rep(as.integer(n_per_cell), K)
  } else if (length(n_per_cell) == K) {
    if (!is.numeric(n_per_cell) || any(n_per_cell < 2L) || any(n_per_cell != round(n_per_cell))) {
      .dmar_refuse_argument(
        "Each entry of 'n_per_cell' must be a positive integer >= 2.",
        arg = "n_per_cell")
    }
    n_per_cell <- as.integer(n_per_cell)
  } else {
    .dmar_refuse_argument(c(
      "'n_per_cell' must be a single number or a numeric vector of length ", K, "."),
      arg = "n_per_cell")
  }

  # ---------- Generate ----------
  if (!requireNamespace("MASS", quietly = TRUE)) {
    .dmar_refuse(
      "Package 'MASS' is required.", "package_required", package = "MASS")
  }
  .dmar_local_seed(seed)

  N_total <- sum(n_per_cell)
  out_y   <- numeric(N_total)
  out_x   <- matrix(NA_real_, nrow = N_total, ncol = q)
  cell_id <- integer(N_total)
  pos     <- 1L

  for (k in seq_len(K)) {
    n_k <- n_per_cell[k]
    mu_k <- c(true_mean_Y[k], true_mean_X[k, ])
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
