#' Print a Model Comparison or ANOVA Table With DMAR \emph{p}-value Formatting
#'
#' Pretty-print an ANOVA-like object (the output of \code{stats::anova},
#' \code{car::Anova}, \code{lmerTest::anova}, etc.) with \emph{p}-values
#' formatted at a fixed number of decimal places (default 4) and with a
#' \dQuote{< 10^(-digits_p)} floor for values too small to express. The
#' default behavior of \code{print.anova} routes \emph{p}-values through
#' \code{stats::format.pval}, which applies its own digit rule
#' (\code{max(1L, getOption("digits") - 2L)}) and switches to scientific
#' notation for tiny values. \code{print_anova()} sidesteps that by
#' converting the \emph{p}-value columns to character strings up front
#' and printing as a data frame.
#'
#' The returned object is the input \code{x} invisibly, unchanged: the
#' underlying numeric \emph{p}-values retain full precision and can
#' still be indexed (for example as \code{x[["Pr(>F)"]]}).
#'
#' Any column whose name starts with \code{Pr(} is formatted as a
#' \emph{p}-value column. Other columns print at whatever
#' \code{getOption("digits")} dictates (so set
#' \code{options(digits = 4)} for a uniformly compact display).
#'
#' @param x An ANOVA-like data frame with one or more \code{Pr(...)}
#'   columns. Accepts \code{anova} objects from \code{stats::anova},
#'   \code{car::Anova}, \code{car::Manova}, and \code{lmerTest::anova}.
#' @param digits_p Integer number of decimal places for the
#'   \emph{p}-value column(s). Default \code{4L}.
#'
#' @return The input \code{x}, invisibly and unchanged.
#'
#' @author Ken Kelley
#'
#' @seealso \code{\link{format_p}}, \code{\link{print_summary}}.
#'
#' @examples
#' fit <- lm(weight ~ Time + Diet, data = ChickWeight)
#' print_anova(anova(fit))
#'
#' print_anova(car::Anova(fit, type = "III"))
#'
#' # Underlying numeric p-values are untouched:
#' a <- anova(fit)
#' print_anova(a)
#' a[["Pr(>F)"]]   # full-precision doubles
#'
#' @export
print_anova <- function(x, digits_p = 4L) {
  hdg <- attr(x, "heading")
  if (!is.null(hdg)) {
    cat(hdg, sep = "\n")
    cat("\n")
  }
  df <- as.data.frame(x)
  pcols <- grep("^Pr\\(", colnames(df), value = TRUE)
  for (col in pcols) {
    df[[col]] <- format_p(df[[col]], digits_p = digits_p)
  }
  print(df, right = TRUE)
  invisible(x)
}
