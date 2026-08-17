# Tidy descriptive statistics for one or more variables.
#' Descriptive Statistics for One or More Variables
#'
#' Computes a compact set of descriptive statistics useful for data screening and
#' psychometric work (including skewness and excess kurtosis), with an optional
#' correlation matrix for the numeric variables.
#'
#' @param x A \code{data.frame}, tibble, or \code{matrix} containing the variables
#'   to summarize. A \code{matrix} is coerced to a \code{data.frame}.
#' @param correlations Logical. If \code{TRUE}, also return a correlation matrix for
#'   the numeric variables (see Details).
#' @param listwise Logical. If \code{TRUE}, apply listwise deletion (drop every row
#'   containing any \code{NA}) to \code{x} \emph{before} computing any statistics.
#'   If \code{FALSE} (the default), each variable is summarized using its own
#'   available (non-missing) observations.
#'
#' @return A list with two elements, returned in a stable shape regardless of the
#'   \code{correlations} argument:
#'   \describe{
#'     \item{\code{descriptives}}{A \code{data.frame} with one row per input
#'       variable and columns \code{variable}, \code{type}, \code{n},
#'       \code{n_missing}, \code{prop_missing}, \code{mean}, \code{median},
#'       \code{sd}, \code{min}, \code{max}, \code{q25}, \code{q75},
#'       \code{skewness}, and \code{kurtosis} (excess kurtosis). For non-numeric
#'       variables, the numeric summary columns are \code{NA}.}
#'     \item{\code{correlations}}{A \eqn{p \times p} correlation matrix among the
#'       numeric variables, or \code{NULL} if \code{correlations = FALSE} or fewer
#'       than two numeric variables are available. This is returned as a plain
#'       matrix (not a table), because its natural structure is a symmetric
#'       two-dimensional array; that format reads well and plugs directly into
#'       downstream tools such as \code{\link[stats]{cov2cor}} or factor analysis
#'       and SEM software.}
#'   }
#'
#' @details
#' \strong{Skewness and kurtosis.} The reported values use the bias-corrected
#' formulas commonly referred to as SAS/SPSS Type 2:
#' \deqn{\mathrm{skewness} = \frac{n}{(n-1)(n-2)} \sum_{i=1}^{n} \left(\frac{x_i - \bar{x}}{s}\right)^3,}
#' \deqn{\mathrm{kurtosis} = \frac{n(n+1)}{(n-1)(n-2)(n-3)} \sum_{i=1}^{n} \left(\frac{x_i - \bar{x}}{s}\right)^4 - \frac{3(n-1)^2}{(n-2)(n-3)},}
#' where \eqn{s} is the (divisor-\eqn{n-1}) sample standard deviation. Kurtosis is
#' reported as excess kurtosis, so a normal distribution has an expected value of 0.
#' As a rough guide to deciding whether normal-theory inference (e.g., maximum
#' likelihood in factor analysis or SEM) is defensible, values of
#' \eqn{|\mathrm{skewness}| > 2} or \eqn{|\mathrm{kurtosis}| > 7} are frequently
#' flagged as problematic.
#'
#' \strong{Variable-type handling.} Each column is classified as
#' \code{"numeric"}, \code{"integer"}, \code{"logical"}, \code{"factor"},
#' \code{"character"}, or whatever its first class is otherwise. Columns of type
#' \code{"numeric"}, \code{"integer"}, and \code{"logical"} receive full
#' distributional summaries (logicals are coerced so that \code{mean} is the
#' proportion of \code{TRUE}s). All other columns report only \code{n},
#' \code{n_missing}, and \code{prop_missing}; their remaining columns are
#' \code{NA}.
#'
#' \strong{Correlations.} The correlation matrix, when requested, uses Pearson
#' correlations with \code{use = "pairwise.complete.obs"} so that each pairwise
#' correlation uses the maximum information available. Only numeric-type variables
#' are included.
#'
#' @examples
#' # Four cognitive tests from the Holzinger and Swineford (1939) study
#' # (301 children in two schools).
#' hs_tests <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
#'                                     "t4_lozenges",
#'                                     "t6_paragraph_comprehension")]
#' descriptives(hs_tests)
#'
#' # Include a correlation matrix (useful during scale development).
#' descriptives(hs_tests, correlations = TRUE)
#'
#' # Mixed-type data: numeric summaries for the test scores, and
#' # type = "factor" (with NA numeric columns) for school.
#' descriptives(holzinger_swineford[, c("school", "t1_visual_perception",
#'                                      "t2_cubes")])
#'
#' # Data with missing values: the revised paper form board and flags tests
#' # were administered only in the Grant-White school, so 156 of the 301
#' # children have no score. Per-variable N and missingness are reported.
#' hs_partial <- holzinger_swineford[, c("t1_visual_perception",
#'                                       "t25_paper_form_board_r",
#'                                       "t26_flags")]
#' descriptives(hs_partial)
#'
#' # The same data with listwise deletion applied first.
#' descriptives(hs_partial, listwise = TRUE)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link[stats]{cor}}, \code{\link[stats]{quantile}}
#'
#' @keywords univar multivariate
#'
#' @family descriptive statistics
#'
#' @export
#' @import stats

descriptives <- function(x, correlations = FALSE, listwise = FALSE) {
  # Normalize input to a data.frame (tibble inherits from data.frame, so it passes
  # through; a matrix is coerced column by column).
  if (is.matrix(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) stop("'x' must be a data frame, tibble, or matrix.")
  if (ncol(x) == 0) stop("'x' has no columns to summarize.")

  # Optional listwise deletion. This removes rows with any NA before computing
  # any statistics, so per-variable 'n_missing' values will be zero afterward.
  if (isTRUE(listwise)) {
    x <- x[stats::complete.cases(x), , drop = FALSE]
    if (nrow(x) == 0) stop("No complete cases remain after listwise deletion.")
  }

  # Classify each variable so we can branch numeric vs. non-numeric summaries.
  var_names  <- names(x)
  var_types  <- vapply(x, .variable_type, character(1))
  is_numlike <- var_types %in% c("numeric", "integer", "logical")

  # Build per-variable summary rows and stack them into a single data.frame.
  per_variable <- lapply(seq_along(x), function(i) {
    col <- x[[i]]
    if (is_numlike[i]) {
      .numeric_stats(as.numeric(col))
    } else {
      .nonnumeric_stats(col)
    }
  })
  stats_df <- do.call(rbind, per_variable)

  descriptives_tbl <- data.frame(
    variable = var_names,
    type     = var_types,
    stats_df,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  # Correlation matrix (numeric-like variables only).
  correlation_matrix <- NULL
  if (isTRUE(correlations)) {
    numeric_x <- x[, is_numlike, drop = FALSE]
    if (ncol(numeric_x) < 2) {
      warning("A correlation matrix requires at least two numeric variables; returning NULL.")
    } else {
      numeric_x[] <- lapply(numeric_x, as.numeric)
      correlation_matrix <- stats::cor(numeric_x, use = "pairwise.complete.obs")
    }
  }

  list(
    descriptives = descriptives_tbl,
    correlations = correlation_matrix
  )
}


# Internal helpers. These are file-local and use the leading-dot naming convention
# reserved in this package for functions that are not exported.

.variable_type <- function(col) {
  if (is.factor(col))    return("factor")
  if (is.logical(col))   return("logical")
  if (is.integer(col))   return("integer")
  if (is.numeric(col))   return("numeric")
  if (is.character(col)) return("character")
  class(col)[1]
}

.numeric_stats <- function(x) {
  x_valid <- x[!is.na(x)]
  n_valid <- length(x_valid)
  n_miss  <- sum(is.na(x))
  n_total <- n_valid + n_miss

  # All-missing or empty: return the missingness bookkeeping with NA stats.
  if (n_valid == 0) {
    return(data.frame(
      n = 0L, n_missing = n_miss, prop_missing = if (n_total == 0) NA_real_ else n_miss / n_total,
      mean = NA_real_, median = NA_real_, sd = NA_real_,
      min = NA_real_, max = NA_real_, q25 = NA_real_, q75 = NA_real_,
      skewness = NA_real_, kurtosis = NA_real_
    ))
  }

  qs <- stats::quantile(x_valid, probs = c(.25, .75), names = FALSE)

  data.frame(
    n            = n_valid,
    n_missing    = n_miss,
    prop_missing = n_miss / n_total,
    mean         = mean(x_valid),
    median       = stats::median(x_valid),
    sd           = stats::sd(x_valid),
    min          = min(x_valid),
    max          = max(x_valid),
    q25          = qs[1],
    q75          = qs[2],
    skewness     = skewness(x_valid),
    kurtosis     = kurtosis(x_valid)
  )
}

.nonnumeric_stats <- function(col) {
  n_miss  <- sum(is.na(col))
  n_valid <- length(col) - n_miss
  n_total <- n_valid + n_miss
  data.frame(
    n            = n_valid,
    n_missing    = n_miss,
    prop_missing = if (n_total == 0) NA_real_ else n_miss / n_total,
    mean = NA_real_, median = NA_real_, sd = NA_real_,
    min = NA_real_, max = NA_real_, q25 = NA_real_, q75 = NA_real_,
    skewness = NA_real_, kurtosis = NA_real_
  )
}

# skewness() and kurtosis() are now exported public functions; see
# R/skewness.R and R/kurtosis.R.
