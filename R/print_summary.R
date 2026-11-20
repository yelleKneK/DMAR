#' Print a Model Summary With DMAR \emph{p}-value Formatting
#'
#' Pretty-print a model summary (the output of \code{summary.lm},
#' \code{summary.glm}, or \code{summary} on an \pkg{lme4} or
#' \pkg{lmerTest} fit) with \emph{p}-values formatted at a fixed
#' number of decimal places (default 4) and with a
#' \dQuote{< 10^(-digits_p)} floor for values too small to express.
#' The default \code{print.summary.lm} / \code{print.summary.merMod}
#' routes \emph{p}-values through \code{stats::format.pval}, which
#' applies its own digit rule and switches to scientific notation for
#' tiny values. \code{print_summary()} sidesteps that by converting
#' the \emph{p}-value columns to character strings up front and
#' printing as a data frame.
#'
#' For a linear model, the function prints the coefficient table, the
#' residual standard error and degrees of freedom, the multiple and
#' adjusted \eqn{R^2}, and the omnibus \emph{F} test and its
#' \emph{p}-value. For a mixed-effects model fit through
#' \pkg{lme4} / \pkg{lmerTest}, the function prints the random-effect
#' variances (from \code{lme4::VarCorr}) and the fixed-effect
#' coefficient table.
#'
#' The returned object is the model summary, invisibly and unchanged:
#' the underlying numeric \emph{p}-values retain full precision and
#' can still be indexed (for example as
#' \code{coef(summary(fit))[, "Pr(>|t|)"]}).
#'
#' @param fit A fitted model object with a \code{summary} method that
#'   returns coefficients via \code{coef(summary(fit))}, including a
#'   \code{Pr(...)} column. Tested with \code{lm}, \code{glm},
#'   \code{lme4::lmer}, and \code{lmerTest::lmer}.
#' @param digits_p Integer number of decimal places for the
#'   \emph{p}-value column(s). Default \code{4L}.
#'
#' @return The model summary, invisibly and unchanged.
#'
#' @author Ken Kelley
#'
#' @seealso \code{\link{format_p}}, \code{\link{print_anova}}.
#'
#' @examples
#' fit_lm <- lm(weight ~ Time + Diet, data = ChickWeight)
#' print_summary(fit_lm)
#'
#' fit_lmer <- lme4::lmer(weight ~ Time + (1 | Chick), data = ChickWeight)
#' print_summary(fit_lmer)
#'
#' # Underlying numeric p-values are untouched:
#' sm <- summary(fit_lm)
#' sm$coefficients[, "Pr(>|t|)"]   # full-precision doubles
#'
#' @export
print_summary <- function(fit, digits_p = 4L) {
  sm <- summary(fit)
  coefs <- tryCatch(coef(sm), error = function(e) NULL)
  if (is.null(coefs)) {
    print(sm)
    return(invisible(sm))
  }
  df <- as.data.frame(coefs)
  pcols <- grep("^Pr\\(", colnames(df), value = TRUE)
  for (col in pcols) {
    df[[col]] <- format_p(df[[col]], digits_p = digits_p)
  }
  is_mer <- inherits(fit, "merMod")
  if (is_mer && requireNamespace("lme4", quietly = TRUE)) {
    cat("Random effects:\n")
    print(lme4::VarCorr(fit))
    cat("\nFixed effects:\n")
  } else {
    cat("Coefficients:\n")
  }
  print(df, right = TRUE)
  if (!is.null(sm$sigma) && !is_mer) {
    cat("\nResidual standard error: ", formatC(sm$sigma, digits = 4),
        " on ", sm$df[2L], " degrees of freedom\n", sep = "")
  }
  if (!is.null(sm$r.squared)) {
    cat("Multiple R-squared: ", formatC(sm$r.squared, digits = 4),
        ",  Adjusted R-squared: ", formatC(sm$adj.r.squared, digits = 4),
        "\n", sep = "")
    if (!is.null(sm$fstatistic)) {
      fst <- sm$fstatistic
      pv  <- pf(fst[1L], fst[2L], fst[3L], lower.tail = FALSE)
      cat("F-statistic: ", formatC(fst[1L], digits = 4),
          " on ", fst[2L], " and ", fst[3L], " DF, p-value: ",
          format_p(pv, digits_p = digits_p), "\n", sep = "")
    }
  }
  invisible(sm)
}
