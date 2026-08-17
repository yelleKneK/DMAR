#' Format \emph{p}-values for Display the DMAR Way
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
#'   through as \code{NA_character_}.
#' @param digits_p Integer number of decimal places. Default
#'   \code{4L}. Values strictly below \eqn{10^{-\mathrm{digits\_p}}}
#'   render as the floor label \dQuote{< 10^(-digits_p)} (for example
#'   \dQuote{< 0.0001} when \code{digits_p = 4}).
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
    stop("'p' must be numeric.", call. = FALSE)
  }
  if (!is.numeric(digits_p) || length(digits_p) != 1L ||
      !is.finite(digits_p) || digits_p < 1L) {
    stop("'digits_p' must be a single positive integer.", call. = FALSE)
  }
  dp     <- as.integer(digits_p)
  thresh <- 10^(-dp)
  floor_label <- paste0("< ", formatC(thresh, format = "f", digits = dp))
  ifelse(is.na(p), NA_character_,
  ifelse(p < thresh, floor_label,
         formatC(p, format = "f", digits = dp)))
}
