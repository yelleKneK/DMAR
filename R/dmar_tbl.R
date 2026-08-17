#' Printing for DMAR Result Tables
#'
#' Most \pkg{DMAR} estimation and testing functions return a tidy
#' \code{data.frame} with a \code{term} column and one or more numeric
#' columns. Because a single numeric column often holds quantities on
#' very different scales (for example, whole-number degrees of freedom
#' alongside an \emph{F} statistic, an effect size, and a small
#' \emph{p}-value), the base \code{\link[base]{print.data.frame}}
#' method formats the whole column with one common format and is easily
#' pushed into scientific notation with many trailing digits. The
#' \code{dmar_tbl} class supplies \code{print} and \code{format}
#' methods that format each value on its own terms: whole numbers (such
#' as degrees of freedom and sample sizes) print without a decimal
#' part, other values print to a small number of significant figures,
#' and scientific notation is reserved for magnitudes where it is the
#' clearer choice (for example, a very small \emph{p}-value).
#'
#' The stored numeric values are never rounded; only their display
#' changes, so downstream arithmetic on the returned object (confidence
#' interval widths, further calculations) uses full precision.
#'
#' The same formatting applies to every \code{dmar_tbl}, whether the
#' table is long (a \code{term} column beside a single \code{value} or
#' \code{estimate} column) or wide (a leading label column such as
#' \code{term}, \code{effect}, or \code{sample_type} beside several
#' typed numeric columns). Each numeric column is formatted on its own
#' terms, so the shape of the table does not matter.
#'
#' Display precision is controlled by the \code{digits} argument or,
#' globally, by \code{options(dmar.digits = )}. The default is 3
#' significant figures. (The option is dot-named because that is the R
#' convention for package options, for example \code{dplyr.width} and
#' \code{knitr.table.format}; it is not a function or argument name and
#' so is outside the package's snake_case rule.)
#'
#' \emph{p}-values are shown to a fixed number of decimal places (four
#' by default, set by \code{digits_p}) rather than to significant
#' figures, which is the conventional way to report them. A
#' \emph{p}-value smaller than the smallest magnitude those decimals
#' can represent prints as \dQuote{< 0.0001} instead of rounding to
#' \code{0.0000}. A column is treated as holding \emph{p}-values when
#' it is named \code{p_value}, \code{p.value}, \code{p_adjusted} (the
#' multiplicity-adjusted case, as in \code{\link{ci_dunnett}}), or
#' \code{p_chi_square} (the exact-fit test of a fitted model, as in
#' \code{\link{measurement_invariance}}); in a long-format
#' table whose quantities share a single \code{value} column, the rows
#' to format this way are named by the producing function through a
#' \code{p_terms} attribute.
#'
#' A few quantities read better at a fixed number of decimal places than
#' at significant figures even though they are not \emph{p}-values:
#' information criteria such as AIC and BIC, and log-likelihoods, where a
#' model comparison difference of a few points would be rounded away by
#' three significant figures (an AIC of 2284.830 would otherwise print as
#' 2280). The producing function names these rows through a
#' \code{fixed_terms} attribute, and they print to \code{digits_fixed}
#' decimal places (three by default).
#'
#' To see more precision than the display shows, raise \code{digits}
#' (for example \code{print(x, digits = 8)}) or read the columns
#' directly, since the stored values are never rounded: \code{x$value}
#' or \code{x[["p_value"]]} returns the numbers at full precision.
#'
#' @section Using the result in your own code: You do not need to know
#'   anything about S3 classes to use a \code{dmar_tbl}. It is an
#'   ordinary \code{data.frame} with a print method, so everything you
#'   already do with a data frame works: \code{x$value} pulls the
#'   numeric column, \code{x[x$term == "smd", ]} selects a row, and the
#'   full-precision numbers are right there for any further calculation.
#'   Three common needs:
#'   \itemize{
#'     \item \emph{Read one number.} Index it like any data frame, for
#'       example \code{x$value[x$term == "upper_limit"]}. The display
#'       rounds; the stored value does not, so this returns the number
#'       at full precision.
#'     \item \emph{See more (or fewer) digits.} Use \code{print(x,
#'       digits = 6)} for a single table, or \code{options(dmar.digits =
#'       6)} for the rest of the session.
#'     \item \emph{Hand the result to other tools.} \code{tidy(x)}
#'       returns a one-row-per-term table in DMAR's own column
#'       vocabulary (\code{term}, \code{estimate}, \code{se},
#'       \code{statistic}, \code{p_value}, \code{ci_lower},
#'       \code{ci_upper}, and so on), which is the convenient
#'       \dQuote{wide} view for plotting or joining; \code{glance(x)}
#'       returns a one-row model-level summary. Both come from the
#'       \pkg{generics} package and need no extra setup. For a
#'       single-estimand result such as \code{\link{ci_smd}} the table is
#'       already one row, so \code{glance()} coincides with \code{tidy()}
#'       (there are no extra model-level statistics to report); for a
#'       multi-row result such as \code{\link{mlmr}} they differ.
#'   }
#'
#' @param x A \code{dmar_tbl} object (a \code{data.frame} returned
#'   by a DMAR function).
#' @param digits Number of significant figures for non-integer values.
#'   Defaults to \code{getOption("dmar.digits", 3L)}.
#' @param digits_p Number of decimal places for \emph{p}-values.
#'   Defaults to 4. A \emph{p}-value below \code{10^(-digits_p)} prints
#'   as \dQuote{< 0.0001} (with the threshold tracking \code{digits_p}).
#' @param digits_fixed Number of decimal places for \code{fixed_terms}
#'   rows (information criteria such as AIC and BIC, and log-likelihoods).
#'   Defaults to 3.
#' @param ... Additional arguments passed to
#'   \code{\link[base]{print.data.frame}}.
#'
#' @return \code{print.dmar_tbl} returns \code{x} invisibly.
#'   \code{format.dmar_tbl} returns a \code{data.frame} whose numeric
#'   columns have been formatted to character for display.
#'
#' @seealso \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} for the wide one-row-per-term and the
#'   one-row summary views. For a gentle, non-technical tour of how to
#'   read and use DMAR result tables, see the \dQuote{Reading DMAR result
#'   tables} vignette: \code{vignette("dmar_output", package = "DMAR")}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Every DMAR estimation function returns a table that prints this way.
#' x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
#' x                       # rounded for reading; sample sizes have no decimals
#'
#' # The stored numbers keep full precision; only the display rounds.
#' x$value[x$term == "smd"]
#' print(x, digits = 8)    # ask the display for more digits
#'
#' # Pull a single number out, exactly as you would from a data frame.
#' x$value[x$term == "upper_limit"]
#'
#' # The broom verbs give the programmer-friendly wide and summary views.
#' generics::tidy(x)
#' generics::glance(x)
#'
#' # The same display rules apply to wide tables (several typed columns),
#' # for example an effect size with its confidence interval per effect.
#' ci_eta_squared(aov(iq_8 ~ treatment, data = pygmalion))
#'
#' @keywords internal
#' @name dmar_tbl
NULL

# Element-wise display formatter shared by format.dmar_tbl / print.dmar_tbl.
# Whole numbers (e.g., degrees of freedom, sample sizes) print without a
# decimal part; other finite values print to `digits` significant figures,
# with base R deciding fixed vs. scientific notation per value (so a very
# small p-value reads as 1.01e-05 while a mean of 1234.5 reads as 1230).
# Not exported; consumed only by the print/format methods below.
.format_dmar_value <- function(v, digits = 3L) {
  vapply(v, function(x) {
    if (is.na(x)) return(NA_character_)
    if (is.finite(x) && x == round(x) && abs(x) < 1e15)
      return(format(x, scientific = FALSE, trim = TRUE))
    format(signif(x, digits), digits = digits, trim = TRUE)
  }, character(1L))
}

# Fixed-decimal p-value formatter shared by format.dmar_tbl. p-values read
# best to a fixed number of decimal places (the convention is four), not to
# significant figures, so they never slip into scientific notation. A value
# below the smallest representable magnitude (10^(-digits_p)) prints as
# "< 0.0001" rather than rounding to "0.0000". Not exported.
.format_dmar_pvalue <- function(p, digits_p = 4L) {
  thresh <- 10^(-digits_p)
  floor_label <- paste0("< ", formatC(thresh, format = "f", digits = digits_p))
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < thresh) return(floor_label)
    formatC(x, format = "f", digits = digits_p)
  }, character(1L))
}

# Fixed-decimal formatter for non-p-value quantities that read best at a set
# number of decimal places rather than significant figures. Information
# criteria (AIC, BIC) and log-likelihoods are the motivating case: under 3
# significant figures a model comparison difference of a few points is lost
# (e.g., AIC 2284.830 would print as "2280", and two log-likelihoods that
# differ in the first decimal would collapse to the same display). Unlike the
# p-value formatter there is no floor and no scientific fallback. Not exported;
# consumed via the fixed_terms attribute by format.dmar_tbl.
.format_dmar_fixed <- function(v, digits_fixed = 3L) {
  vapply(v, function(x) {
    if (is.na(x)) return(NA_character_)
    formatC(x, format = "f", digits = digits_fixed)
  }, character(1L))
}

# Tag a result table as a dmar_tbl so it prints via print.dmar_tbl
# (whole numbers without a decimal part, other values to `digits`
# significant figures, p-values to fixed decimals). The rule is uniform:
# every tidy-returning DMAR function routes its result through this helper,
# whether the table is long (a `term` column plus a single `value` /
# `estimate` column) or wide (a leading label column such as `term`,
# `effect`, or `sample_type` plus one or more numeric columns). The format
# method works on every numeric column, so the shape of the table does not
# matter.
#
# `conf_level`, when supplied, is recorded for the print footer. `p_terms`,
# when supplied, names the rows of the long-format primary column that hold
# p-values (a wide table needs no p_terms for a column named p_value,
# p.value, p_adjusted, or p_chi_square, as the format method finds those by
# name). `fixed_terms`, when supplied, names
# rows of the primary column to print to a fixed number of decimal places
# instead of significant figures (information criteria such as AIC / BIC).
# The dmar_tbl class is inserted just before `data.frame`, so any leading
# subclass (e.g., dmar_ci_long, dmar_ss_power) keeps dispatching its own
# tidy() / glance() methods while print falls through to print.dmar_tbl.
# Idempotent. Not exported; this is the one-line opt-in every tidy-returning
# DMAR function calls on its result just before returning it.
.as_dmar_tbl <- function(out, conf_level = NULL, p_terms = NULL,
                         fixed_terms = NULL, ss_type = NULL,
                         subclass = NULL) {
  if (!is.null(conf_level))  attr(out, "conf_level")  <- conf_level
  if (!is.null(p_terms))     attr(out, "p_terms")     <- p_terms
  if (!is.null(fixed_terms)) attr(out, "fixed_terms") <- fixed_terms
  if (!is.null(ss_type))     attr(out, "ss_type")     <- ss_type
  cls <- class(out)
  if (!("dmar_tbl" %in% cls))
    class(out) <- c(setdiff(cls, "data.frame"), "dmar_tbl", "data.frame")
  # A leading dispatch subclass (e.g. dmar_ss_aipe, dmar_ss_power) rides in
  # front of dmar_tbl so family tidiers dispatch while printing falls
  # through to print.dmar_tbl; idempotent like the dmar_tbl insertion.
  if (!is.null(subclass) && !(subclass %in% class(out)))
    class(out) <- c(subclass, class(out))
  out
}

#' @rdname dmar_tbl
#' @export
format.dmar_tbl <- function(x, digits = getOption("dmar.digits", 3L),
                            digits_p = 4L, digits_fixed = 3L, ...) {
  raw <- x
  class(raw) <- "data.frame"
  out <- raw
  num <- vapply(out, is.numeric, logical(1L))
  out[num] <- lapply(out[num], .format_dmar_value, digits = digits)

  # The "primary" numeric column carries the per-term quantities in a long
  # table: `value`, or `estimate` for coefficient-style tables. The p_terms
  # and fixed_terms attributes name rows of this column to reformat.
  primary <- intersect(c("value", "estimate"), names(out)[num])
  primary <- if (length(primary)) primary[1L] else NA_character_

  # p-values get fixed decimals, overriding the significant-figure format
  # above. They are detected two ways: a dedicated column (p_value /
  # p.value / p_adjusted / p_chi_square) in a wide table, or named rows of
  # the primary column in a long table (the producing function lists them in
  # a p_terms attribute). The names are package-wide conventions, so a
  # function that reports an exact-fit p-value as p_chi_square, as
  # measurement_invariance() does, formats it the same way as the p_value
  # column beside it.
  for (nm in intersect(c("p_value", "p.value", "p_adjusted", "p_chi_square"),
                       names(out)[num]))
    out[[nm]] <- .format_dmar_pvalue(raw[[nm]], digits_p = digits_p)
  p_terms <- attr(x, "p_terms")
  if (!is.null(p_terms) && "term" %in% names(out) && !is.na(primary)) {
    idx <- raw$term %in% p_terms
    out[[primary]][idx] <- .format_dmar_pvalue(raw[[primary]][idx],
                                               digits_p = digits_p)
  }

  # fixed_terms get a fixed number of decimal places (digits_fixed) instead of
  # significant figures: information criteria such as AIC / BIC, where a
  # model comparison difference of a few points would otherwise round away.
  fixed_terms <- attr(x, "fixed_terms")
  if (!is.null(fixed_terms) && "term" %in% names(out) && !is.na(primary)) {
    idx <- raw$term %in% fixed_terms
    out[[primary]][idx] <- .format_dmar_fixed(raw[[primary]][idx],
                                              digits_fixed = digits_fixed)
  }
  out
}

#' @rdname dmar_tbl
#' @export
print.dmar_tbl <- function(x, digits = getOption("dmar.digits", 3L),
                           digits_p = 4L, digits_fixed = 3L, ...) {
  disp <- format(x, digits = digits, digits_p = digits_p,
                 digits_fixed = digits_fixed)
  # row.names / right are passed through `...` (the documented contract), but
  # default to FALSE for the look. Setting them via do.call rather than as
  # fixed arguments lets a caller override them, e.g. print(x, row.names = TRUE),
  # without colliding with a hard-coded value.
  dots <- list(...)
  if (is.null(dots[["row.names"]])) dots[["row.names"]] <- FALSE
  if (is.null(dots[["right"]]))     dots[["right"]]     <- FALSE
  do.call(print.data.frame, c(list(disp), dots))
  ss_type <- attr(x, "ss_type")
  if (!is.null(ss_type)) {
    roman <- switch(as.character(ss_type),
                    "1" = "I", "2" = "II", "3" = "III", as.character(ss_type))
    cat(sprintf("\nSum of squares: Type %s\n", roman))
  }
  conf_level <- attr(x, "conf_level")
  if (!is.null(conf_level))
    cat(sprintf("\nConfidence level: %g%%\n", 100 * conf_level))
  invisible(x)
}
