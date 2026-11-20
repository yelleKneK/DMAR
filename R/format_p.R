#' Format \emph{p}-Values for Display at a Fixed Number of Decimal Places
#'
#' Render a vector of \emph{p}-values as character strings at a fixed
#' number of decimal places (default 4) with a floor label
#' \dQuote{< 10^(-digits_p)} for values too small to express. Never uses
#' scientific notation. This is the package-wide convention for
#' displaying \emph{p}-values used by the \code{\link{dmar_tbl}} print
#' layer, by \code{\link{correlations_test}}, and by the helper
#' display functions \code{\link{print_anova}} and
#' \code{\link{print_summary}}. Exposing it as a user-facing function
#' lets analysts apply the same convention to ad hoc \emph{p}-values
#' that they want to report in prose or in a manually constructed
#' table.
#'
#' The function does not modify the input value; the underlying
#' numeric \emph{p}-value retains full precision and can still be
#' indexed out of whatever object holds it. Only the returned display
#' string is rounded.
#'
#' @param p Numeric vector of \emph{p}-values. \code{NA} values pass
#'   through as \code{NA_character_}. The function only formats: it shows
#'   each value as given and does not check that it is a probability. A
#'   value below 0 (and \code{-Inf}) is below the floor and shows as the
#'   floor label, a value above 1 is shown at \code{digits_p} decimal
#'   places (1.1 as \dQuote{1.1000}), \code{Inf} is shown as Inf (with
#'   the leading spaces \code{\link[base]{formatC}} gives it), and
#'   \code{NaN} is returned as \code{NA_character_}.
#' @param digits_p Number of decimal places, a single whole number from 1
#'   to 324, the most decimal places \code{\link[base]{formatC}} writes
#'   (\code{3} and \code{3L} are accepted; \code{3.9} is refused with a
#'   condition of class \code{dmar_argument_invalid}, rather than cut to
#'   3). Default \code{4L}. Values strictly below
#'   \eqn{10^{-\mathrm{digits\_p}}} render as the floor label
#'   \dQuote{< 10^(-digits_p)} (for example \dQuote{< 0.0001} when
#'   \code{digits_p = 4}).
#'
#' @return A character vector of the same length as \code{p}.
#'
#' @author Ken Kelley
#'
#' @seealso \code{\link{print_anova}}, \code{\link{print_summary}},
#'   \code{\link{dmar_tbl}}.
#'
#' @examples
#' # Round to four decimals with a "< 0.0001" floor:
#' format_p(c(0.5, 0.0234, 0.0001234, 1e-10, NA))
#'
#' # Six decimals when more precision is wanted:
#' format_p(0.0001234, digits_p = 6)
#'
#' # Use inline in prose for a publication-style summary:
#' fit <- lm(weight ~ Time + Diet, data = ChickWeight)
#' p_time <- summary(fit)$coefficients["Time", "Pr(>|t|)"]
#' paste0("The Time coefficient was significant (p = ", format_p(p_time), ").")
#'
#' @export
format_p <- function(p, digits_p = 4L) {
  if (!is.numeric(p)) {
    .dmar_refuse_argument("'p' must be numeric.", arg = "p")
  }
  # A whole number is required: through the first builds of 1.1.0 a
  # fractional digits_p was cut by as.integer() (3.9 showed three decimal
  # places). The comparison with round() is exact, so 3 and 3L pass. The
  # largest is 324, the most decimal places formatC() writes (and the
  # place of the smallest positive double); a larger value was cut to 324
  # by formatC() with a warning, or, above the largest integer, became NA
  # in as.integer() with a warning.
  if (!is.numeric(digits_p) || length(digits_p) != 1L ||
      !is.finite(digits_p) || digits_p < 1L || digits_p > 324L ||
      digits_p != round(digits_p)) {
    .dmar_refuse_argument(c(
      "'digits_p' must be a single positive integer, a whole number of ",
      "decimal places from 1 to 324; ", .dmar_given(digits_p), "."),
      arg = "digits_p", value = digits_p)
  }
  dp     <- as.integer(digits_p)
  thresh <- 10^(-dp)
  floor_label <- paste0("< ", formatC(thresh, format = "f", digits = dp))
  ifelse(is.na(p), NA_character_,
  ifelse(p < thresh, floor_label,
         formatC(p, format = "f", digits = dp)))
}
