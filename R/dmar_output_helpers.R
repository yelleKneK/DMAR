#' Publication-Ready Display of DMAR Result Tables
#'
#' The \code{\link{dmar_tbl}} print layer formats a result table for the
#' console. These helpers carry the same formatting into the two places a
#' researcher writes up an analysis: a knitted report and a results sentence.
#'
#' \code{knit_print.dmar_tbl} is the \pkg{knitr} print method, so a
#' \code{dmar_tbl} dropped into an R Markdown chunk renders as a formatted
#' \code{\link[knitr]{kable}} (sensible rounding, whole-number sample sizes,
#' \emph{p}-values to fixed decimals) rather than as a raw dump of doubles.
#' \code{as_kable} is the explicit form of the same rendering: it returns the
#' \code{knitr_kable} object so the caller can pipe it into
#' further styling or embed it in a larger document. Both reuse
#' \code{\link{format.dmar_tbl}}, so what a reader sees in a report matches
#' what they saw at the console, and neither rounds the stored numbers.
#'
#' \code{results_sentence} turns a table that carries a confidence interval into the
#' one sentence an author puts in a results section, for example
#' \dQuote{smd = 0.50, 95\% CI [0.10, 0.90]}. It reads the estimate and its
#' limits from the numeric columns at full precision and formats them for the
#' sentence, so the reported numbers are exact to the requested decimals rather
#' than transcribed from the rounded console display.
#'
#' @param x A \code{dmar_tbl} object (a \code{data.frame} returned by a DMAR
#'   function). For \code{results_sentence}, a \code{dmar_tbl} that carries a
#'   confidence interval, either as \code{lower_limit} / \code{upper_limit}
#'   rows of a long table or as \code{lower_limit} / \code{upper_limit}
#'   columns of a wide table.
#' @param format Passed to \code{\link[knitr]{kable}} as its \code{format}
#'   argument (for example \code{"html"}, \code{"latex"}, \code{"pipe"}). The
#'   default \code{NULL} lets \pkg{knitr} choose based on the output context.
#' @param label For \code{results_sentence}, the label that leads the sentence
#'   (for example \code{"Cohen's d"}). The default \code{NULL} uses the name of
#'   the point-estimate term.
#' @param digits For \code{results_sentence}, the number of decimal places for
#'   the estimate and the interval limits. Default 2.
#' @param ... Additional arguments. For \code{knit_print.dmar_tbl} and
#'   \code{as_kable.dmar_tbl}, passed to \code{\link[knitr]{kable}}.
#'
#' @return \code{knit_print.dmar_tbl} returns a \code{knit_asis} object (the
#'   rendered table) for \pkg{knitr} to place in the document.
#'   \code{as_kable} returns a \code{knitr_kable} object.
#'   \code{results_sentence} returns a length-one character string.
#'
#' @seealso \code{\link{dmar_tbl}} for the console print layer and
#'   \code{\link{format_p}} for the \emph{p}-value convention these helpers
#'   reuse.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
#'
#' # A knitr_kable that keeps every column, ready for a report.
#' as_kable(x)
#'
#' # The sentence an author writes in a results section.
#' results_sentence(x, label = "Cohen's d")
#'
#' # Wide tables (one interval per row) work the same way.
#' results_sentence(ci_R2(R2 = 0.25, N = 100, p = 5))
#'
#' @name dmar_output_helpers
NULL

#' @rdname dmar_output_helpers
#' @exportS3Method knitr::knit_print
knit_print.dmar_tbl <- function(x, ...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    # Without knitr there is nothing to render into; fall back to the console
    # print method so the object still shows something useful.
    print(x)
    return(invisible(x))
  }
  disp <- format(x)
  tbl <- knitr::kable(disp, row.names = FALSE, ...)

  # Append the same footer the console print method shows (sum-of-squares type
  # and confidence level) as a trailing note beneath the table.
  cap <- .dmar_tbl_footer(x)
  body <- if (length(cap)) c(tbl, "", cap) else tbl
  knitr::knit_print(knitr::asis_output(paste(body, collapse = "\n")))
}

#' Coerce a DMAR Result Table to a \code{knitr::kable}
#'
#' @rdname dmar_output_helpers
#' @export
as_kable <- function(x, ...) {
  UseMethod("as_kable")
}

#' @rdname dmar_output_helpers
#' @export
as_kable.dmar_tbl <- function(x, format = NULL, ...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for as_kable(). ",
         "Install with install.packages('knitr').", call. = FALSE)
  }
  # Reuse the console formatter so every numeric column (long or wide) is
  # rendered on its own terms, then hand the character table to knitr::kable.
  disp <- format(x)
  cap <- .dmar_tbl_footer(x)
  caption <- if (length(cap)) paste(cap, collapse = "; ") else NULL
  tbl <- knitr::kable(disp, format = format, row.names = FALSE,
                      caption = caption, ...)
  tbl
}

#' Write a Publication-Style Results Sentence From a DMAR Result Table
#'
#' @rdname dmar_output_helpers
#' @export
results_sentence <- function(x, label = NULL, digits = 2) {
  if (!inherits(x, "dmar_tbl") && !is.data.frame(x)) {
    stop("'x' must be a dmar_tbl (a DMAR result table).", call. = FALSE)
  }
  nms <- names(x)

  if (all(c("lower_limit", "upper_limit") %in% nms) &&
      !("term" %in% nms && "value" %in% nms)) {
    # Wide table: lower_limit / upper_limit are columns. The estimate is the
    # first non-limit numeric column (the primary column, e.g. eta_squared).
    if (nrow(x) != 1L) {
      stop("results_sentence() needs a single-row table; this table has ",
           nrow(x), " rows. Select one row first, e.g. x[x$effect == ..., ].",
           call. = FALSE)
    }
    num <- vapply(x, is.numeric, logical(1L))
    est_col <- setdiff(names(x)[num],
                       c("lower_limit", "upper_limit", "N",
                         "p_value", "p.value", "p_adjusted"))
    if (!length(est_col))
      stop("Could not find a point-estimate column in this table.",
           call. = FALSE)
    est_col <- est_col[1L]
    est <- x[[est_col]][1L]
    lo  <- x[["lower_limit"]][1L]
    hi  <- x[["upper_limit"]][1L]
    est_label <- if (is.null(label)) est_col else label
    p_val <- .results_sentence_p_wide(x)
  } else if ("term" %in% nms) {
    # Long table: lower_limit / upper_limit are rows of `term`, beside a point-
    # estimate row, in the primary numeric column (`value` or `estimate`).
    primary <- intersect(c("value", "estimate"),
                         names(x)[vapply(x, is.numeric, logical(1L))])
    if (!length(primary))
      stop("results_sentence() needs a numeric 'value' or 'estimate' column.",
           call. = FALSE)
    primary <- primary[1L]
    if (!all(c("lower_limit", "upper_limit") %in% x$term))
      stop("This table carries no confidence interval ",
           "(no 'lower_limit'/'upper_limit' rows), so no sentence can be ",
           "written. results_sentence() needs an interval-carrying table.",
           call. = FALSE)
    lo <- x[[primary]][x$term == "lower_limit"][1L]
    hi <- x[[primary]][x$term == "upper_limit"][1L]
    est_term <- setdiff(x$term, c("lower_limit", "upper_limit"))
    if (!length(est_term))
      stop("Could not find a point-estimate row in this table.", call. = FALSE)
    est_term <- est_term[1L]
    est <- x[[primary]][x$term == est_term][1L]
    est_label <- if (is.null(label)) est_term else label
    p_val <- .results_sentence_p_long(x, primary)
  } else {
    stop("This table carries no confidence interval, so no sentence can be ",
         "written. results_sentence() needs an interval-carrying table.",
         call. = FALSE)
  }

  conf_level <- attr(x, "conf_level")
  if (is.null(conf_level)) conf_level <- 0.95
  pct <- .format_conf_pct(conf_level)

  fmt <- function(v) formatC(v, format = "f", digits = digits)
  out <- sprintf("%s = %s, %s%% CI [%s, %s]",
                 est_label, fmt(est), pct, fmt(lo), fmt(hi))
  if (!is.null(p_val) && !is.na(p_val))
    out <- paste0(out, ", p = ", format_p(p_val))
  out
}

# Footer lines shared by knit_print / as_kable: the sum-of-squares type and the
# confidence level, matching what print.dmar_tbl writes below the table. Returns
# a character vector (possibly empty). Not exported.
.dmar_tbl_footer <- function(x) {
  lines <- character(0)
  ss_type <- attr(x, "ss_type")
  if (!is.null(ss_type)) {
    roman <- switch(as.character(ss_type),
                    "1" = "I", "2" = "II", "3" = "III", as.character(ss_type))
    lines <- c(lines, sprintf("Sum of squares: Type %s", roman))
  }
  conf_level <- attr(x, "conf_level")
  if (!is.null(conf_level))
    lines <- c(lines, sprintf("Confidence level: %s%%",
                              .format_conf_pct(conf_level)))
  lines
}

# Format a confidence level (a proportion such as 0.95) as a percent for
# display, dropping a trailing ".0" so 0.95 reads "95" not "95.0". Not exported.
.format_conf_pct <- function(conf_level) {
  # Keep 0.95 as "95" and 0.975 as "97.5", with no trailing ".0" and no width
  # padding (format(scientific = FALSE, trim = TRUE) does both).
  format(100 * conf_level, scientific = FALSE, trim = TRUE)
}

# Pull a p-value out of a wide interval-carrying table, if one is present. The
# column is detected by the same names the format method uses. Not exported.
.results_sentence_p_wide <- function(x) {
  pcol <- intersect(c("p_value", "p.value", "p_adjusted"), names(x))
  if (length(pcol)) x[[pcol[1L]]][1L] else NULL
}

# Pull a p-value out of a long interval-carrying table, if one is present:
# either a dedicated p-value column or a p_terms-named row of the primary
# column. Not exported.
.results_sentence_p_long <- function(x, primary) {
  pcol <- intersect(c("p_value", "p.value", "p_adjusted"), names(x))
  if (length(pcol)) return(x[[pcol[1L]]][1L])
  p_terms <- attr(x, "p_terms")
  if (!is.null(p_terms)) {
    hit <- x$term %in% p_terms
    if (any(hit)) return(x[[primary]][hit][1L])
  }
  NULL
}
