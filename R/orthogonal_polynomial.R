# Orthogonal-polynomial (trend) contrast coefficients for a quantitative factor.
#' Orthogonal-Polynomial (Trend) Contrast Coefficients
#'
#' Builds the set of orthogonal-polynomial contrasts (linear,
#' quadratic, cubic, \ldots) for a factor whose \eqn{a} levels are
#' quantitative, so that between-group variation can be decomposed into
#' independent trend components. The columns are mutually orthogonal and
#' each sums to zero (so each is a contrast on the group means). By
#' default the coefficients are returned in \emph{orthonormal} form,
#' meaning each column also has unit length, \eqn{\sum_i c_i^2 = 1}; the
#' alternative \code{type = "integer"} rescales every column to the
#' small whole numbers used in the published orthogonal-polynomial
#' tables, which are easier to read by eye and match hand computation.
#' The per-trend sum of squared coefficients \eqn{\sum_i c_i^2} is
#' carried on the returned object and shown when it is printed.
#'
#' @param levels One of: a single integer giving the number of levels
#'   \eqn{a} (labels default to \code{"L1"}, \code{"L2"}, \ldots, and
#'   the levels are treated as equally spaced); a character or factor
#'   vector of level labels (equally spaced unless \code{scores} is
#'   given); or a numeric vector of the quantitative level values
#'   themselves, which are used both as the labels and as the spacing.
#' @param scores Optional numeric vector of length \eqn{a} giving the
#'   quantitative position of each level. Supply this when the labels
#'   are non-numeric but the spacing is unequal (e.g., doses 0, 1, 2,
#'   and 4). Defaults to equally spaced positions \code{1:a}.
#' @param type Either \code{"orthonormal"} (default) for unit-length
#'   columns (\eqn{\sum_i c_i^2 = 1}, identical to
#'   \code{\link[stats]{contr.poly}}) or \code{"integer"} for the
#'   minimal whole-number coefficients of the classic trend-coefficient
#'   table. \code{"integer"} requires equally spaced \code{scores}.
#' @param degree Highest-order trend to return, an integer between
#'   \code{1} and \eqn{a - 1}. Defaults to \eqn{a - 1} (the full set the
#'   design can support).
#'
#' @return A numeric \eqn{a \times \text{degree}} matrix with row names
#'   = the level labels and column names = the trend names
#'   (\code{"linear"}, \code{"quadratic"}, \code{"cubic"},
#'   \code{"quartic"}, \ldots). The matrix can be assigned directly to
#'   \code{contrasts(factor)} or passed to
#'   \code{\link{contrast_test}}, \code{\link{is_orthogonal_set}}, or
#'   \code{\link{ci_c}}. The per-column sum of squared coefficients is
#'   stored as \code{attr(*, "sum_sq")} (a named numeric vector), and
#'   the level spacing and \code{type} are stored as
#'   \code{attr(*, "scores")} and \code{attr(*, "type")}. The object
#'   carries class \code{"orthogonal_polynomial"} so that printing shows
#'   the coefficients alongside \eqn{\sum_i c_i^2}. When printed, the object
#'   is shown in the textbook Table A.10 orientation (trends in rows, levels
#'   in columns) with the \eqn{\sum_i c_i^2} values as a final column; the
#'   stored matrix is the transpose of that display (levels in rows) so it can
#'   be assigned directly to \code{contrasts()}.
#'
#' @details
#' \strong{What a trend contrast is.} When the levels of a factor are
#' quantitative and ordered (minutes of study, dose, day), the omnibus
#' between-group variation can be partitioned into a linear trend (does
#' the mean rise or fall steadily?), a quadratic trend (is there
#' curvature?), a cubic trend (an S-shape?), and so on, up to order
#' \eqn{a - 1}. Each trend is a single 1-\emph{df} contrast on the group
#' means, \eqn{\hat\psi = \sum_i c_i \bar Y_i}, and because the contrasts
#' are mutually orthogonal their sums of squares add up to the omnibus
#' between-group sum of squares exactly.
#'
#' \strong{Orthonormal versus integer scaling.} The two \code{type}
#' values return the \emph{same} trends (same directions, same
#' hypotheses, identical \emph{F}, \emph{t}, and \emph{p} for a given
#' trend); they differ only in how each column is scaled.
#' \itemize{
#'   \item \code{"orthonormal"} normalizes each column to unit length,
#'     \eqn{\sum_i c_i^2 = 1}. The trend portion of the design is then an
#'     orthonormal basis, the sum of squares for a trend is simply
#'     \eqn{\hat\psi^2}, and in a fitted \code{\link[stats]{lm}} the
#'     \code{.linear} / \code{.quadratic} coefficients are the trend
#'     estimates on a common scale. This is the form
#'     \code{\link[stats]{lm}} and \code{\link[stats]{aov}} use
#'     internally and is defined for any spacing, equal or unequal.
#'   \item \code{"integer"} rescales each column to the smallest whole
#'     numbers with the same ratios, reproducing the published
#'     orthogonal-polynomial coefficient table (e.g., for \eqn{a = 4} the
#'     linear contrast is \eqn{-3, -1, 1, 3}). These are easy to read and
#'     to compute with by hand, but the columns no longer share a common
#'     length: \eqn{\sum_i c_i^2} varies by trend (for \eqn{a = 4}: 20,
#'     4, and 20), so the sum of squares for a trend is
#'     \eqn{n\, \hat\psi^2 / \sum_i c_i^2}. Integer coefficients exist
#'     only for equally spaced levels.
#' }
#' The reported \eqn{\sum_i c_i^2} (shown on printing and stored in
#' \code{attr(*, "sum_sq")}) is exactly the divisor in that
#' sum-of-squares formula; for the orthonormal form it is 1 for every
#' trend.
#'
#' \strong{Unequal spacing.} With unequally spaced \code{scores} the
#' orthonormal polynomials are still uniquely defined and are returned by
#' the \code{"orthonormal"} type. Whole-number coefficients generally do
#' not exist in that case, so \code{type = "integer"} is an error.
#'
#' \strong{Relation to base R.} For equally spaced levels the
#' orthonormal output is identical to \code{\link[stats]{contr.poly}(a)};
#' \code{orthogonal_polynomial} adds the integer table form, an explicit
#' \code{scores} argument for unequal spacing, meaningful trend names,
#' and the \eqn{\sum_i c_i^2} report.
#'
#' @references
#' Fisher, R. A., & Yates, F. (1953). \emph{Statistical tables for
#'   biological, agricultural and medical research} (4th ed.). Oliver
#'   and Boyd. (Origin of the tabulated integer coefficients.)
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 6 on trend analysis; Appendix
#'   Table A.10 reports these coefficients.)
#'
#' @seealso \code{\link[stats]{contr.poly}}, \code{\link{effects_coding}},
#'   \code{\link{helmert_coding}}, \code{\link{is_orthogonal_set}},
#'   \code{\link{contrast_test}}, \code{\link{ci_c}}
#'
#' @examples
#' # Trend (orthogonal-polynomial) analysis decomposes the omnibus effect of a
#' # quantitative factor (dose, time, trial block, stimulus intensity) into
#' # independent linear, quadratic, cubic, ... components, as in the trend
#' # analysis of Maxwell, Delaney, and Kelley (2027, Chapter 6).
#'
#' # 1. The a = 4 coefficient table in the integer form of the textbook
#' #    appendix (Table A.10). Printing puts the trends in rows and appends a
#' #    final sum-of-squared-coefficients (sum c^2) column: linear -3 -1 1 3,
#' #    and per-trend sum c^2 of 20, 4, 20.
#' orthogonal_polynomial(4, type = "integer")
#'
#' # 2. The default orthonormal form (identical to stats::contr.poly(4)); every
#' #    trend has sum c^2 = 1, the scale lm() and aov() use internally.
#' orthogonal_polynomial(4)
#'
#' # 3. A worked trend analysis in the style of MDK Chapter 6. An outcome is
#' #    measured at a = 5 equally spaced levels of a quantitative factor (say
#' #    stimulus intensity 1..5), n = 8 per level, from a population with a
#' #    strong linear and a mild quadratic trend. Assigning the trend contrasts
#' #    to the factor makes lm() report each trend test directly.
#' set.seed(113)
#' intensity <- factor(rep(1:5, each = 8))
#' contrasts(intensity) <- orthogonal_polynomial(levels(intensity))
#' y <- c(10, 9, 7, 6, 6)[as.integer(intensity)] +
#'        rnorm(length(intensity), sd = 1.2)
#' round(coef(summary(lm(y ~ intensity))), 3)
#'
#' # 4. The same trends through DMAR's contrast_test(), which reports each
#' #    trend's estimate, standard error, t, p, and confidence interval from a
#' #    fitted one-way model. The integer columns are the contrast weights.
#' op <- orthogonal_polynomial(5, type = "integer")
#' contrast_test(aov(y ~ intensity),
#'               contrasts = list(linear    = op[, "linear"],
#'                                quadratic = op[, "quadratic"],
#'                                cubic     = op[, "cubic"]))
#'
#' # 5. Unequally spaced doses (0, 1, 2, 4 mg): integer coefficients no longer
#' #    exist, but the orthonormal trends remain uniquely defined.
#' orthogonal_polynomial(c("0 mg", "1 mg", "2 mg", "4 mg"),
#'                       scores = c(0, 1, 2, 4))
#'
#' # 6. Confirm a returned set is mutually orthogonal.
#' is_orthogonal_set(orthogonal_polynomial(5, type = "integer"))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family design utilities
#'
#' @export
#' @import stats

orthogonal_polynomial <- function(levels, scores = NULL,
                                  type = c("orthonormal", "integer"),
                                  degree = NULL) {
  type <- match.arg(type)

  # Resolve level labels and their quantitative positions (scores).
  if (is.numeric(levels) && length(levels) > 1L) {
    labs   <- as.character(levels)
    scores <- if (is.null(scores)) as.numeric(levels) else scores
  } else {
    labs <- .normalize_levels(levels)
    if (is.null(scores)) scores <- seq_along(labs)
  }
  a <- length(labs)
  if (a < 3L)
    stop("Need at least 3 levels to define a polynomial (linear) trend.")
  if (!is.numeric(scores) || length(scores) != a)
    stop(sprintf("'scores' must be a numeric vector of length %d.", a))
  if (anyDuplicated(scores))
    stop("'scores' must be distinct.")

  if (is.null(degree)) degree <- a - 1L
  if (length(degree) != 1L || degree < 1L || degree > a - 1L)
    stop(sprintf("'degree' must be a single integer between 1 and %d.", a - 1L))

  d <- diff(sort(scores))
  equally_spaced <- max(abs(d - d[1L])) < 1e-9 * max(abs(d))

  if (type == "integer") {
    if (!equally_spaced)
      stop("type = \"integer\" requires equally spaced 'scores'; ",
           "use type = \"orthonormal\" for unequal spacing.")
    ord <- order(scores)
    P   <- stats::contr.poly(a)                 # orthonormal, ascending order
    M0  <- apply(P, 2L, .poly_min_int)
    if (is.null(dim(M0))) M0 <- matrix(M0, nrow = a)
    M <- matrix(0, a, a - 1L)
    M[ord, ] <- M0                              # restore the input order
  } else {
    P <- stats::poly(scores, degree = a - 1L)
    M <- unclass(P) / rep(sqrt(colSums(P^2)), each = a)
    attributes(M) <- NULL
    dim(M) <- c(a, a - 1L)
  }

  M <- M[, seq_len(degree), drop = FALSE]
  dimnames(M) <- list(labs, .trend_names(degree))

  structure(
    M,
    sum_sq = stats::setNames(colSums(M^2), colnames(M)),
    scores = stats::setNames(scores, labs),
    type   = type,
    class  = c("orthogonal_polynomial", "matrix", "array")
  )
}

#' @export
print.orthogonal_polynomial <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  type   <- attr(x, "type")
  sum_sq <- attr(x, "sum_sq")

  body <- unclass(x)
  attr(body, "sum_sq") <- attr(body, "scores") <- attr(body, "type") <- NULL

  # Display in the textbook (MDK Table A.10) orientation: trends in rows,
  # levels in columns, with the per-trend sum of squared coefficients as the
  # final column rather than an offset block. The stored object is the
  # transpose of this (levels x trends) so it assigns directly to contrasts().
  disp <- t(body)                                   # trends x levels
  cat(sprintf("Orthogonal-polynomial (trend) contrasts (%s), %d level%s\n",
              type, nrow(x), if (nrow(x) == 1L) "" else "s"))
  cat("Layout: trends in rows, levels in columns (MDK Table A.10);",
      "final column is sum c^2.\n")
  cat("Stored object is the transpose (levels x trends) for contrasts().\n\n")

  if (identical(type, "integer")) {
    tab <- cbind(round(disp), "sum c^2" = round(sum_sq))
  } else {
    tab <- cbind(round(disp, digits), "sum c^2" = round(sum_sq, digits))
  }
  print(tab, ...)
  invisible(x)
}

# ---- internal helpers ----

# Greatest common divisor of two integers, and of an integer vector.
.gcd2 <- function(a, b) {
  while (b != 0) { t <- b; b <- a %% b; a <- t }
  abs(a)
}
.gcd_vec <- function(x) {
  x <- x[x != 0]
  if (!length(x)) return(1)
  Reduce(.gcd2, x)
}

# Minimal whole-number representation of one orthogonal-polynomial
# column: divide by the smallest nonzero magnitude, find the smallest
# integer multiplier that clears the remaining fractions, then reduce by
# the gcd. Errors if no integer scaling is found (e.g., unequal spacing
# slipped through), pointing the user at the orthonormal form.
.poly_min_int <- function(col, tol = 1e-6, max_m = 100000L) {
  nz <- abs(col)[abs(col) > tol]
  w  <- col / min(nz)
  m  <- NA_integer_
  for (cand in seq_len(max_m)) {
    if (all(abs(cand * w - round(cand * w)) < tol)) { m <- cand; break }
  }
  if (is.na(m))
    stop("Could not find an integer scaling for a trend column; ",
         "use type = \"orthonormal\".")
  v <- round(m * w)
  as.numeric(v / .gcd_vec(v))
}

# Human-readable trend names for the first 'degree' polynomial terms.
.trend_names <- function(degree) {
  base <- c("linear", "quadratic", "cubic", "quartic", "quintic",
            "sextic", "septic", "octic", "nonic", "decic")
  vapply(seq_len(degree),
         function(k) if (k <= length(base)) base[k] else paste0("degree_", k),
         character(1))
}
