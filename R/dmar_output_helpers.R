#' Publication-Ready Display of DMAR Result Tables
#'
#' The \code{\link{dmar_tbl}} print layer formats a result table for the
#' console. These helpers carry the same formatting into the two places a
#' researcher writes up an analysis: a knitted report and a results sentence.
#'
#' \code{knit_print.dmar_tbl} is the \pkg{knitr} print method, so a
#' \code{dmar_tbl} dropped into an R Markdown chunk renders as a formatted
#' \code{\link[knitr]{kable}} (values to a fixed number of decimal places,
#' whole-number sample sizes without decimals, \emph{p}-values to fixed
#' decimals) rather than as unrounded doubles.
#' \code{as_kable} is the explicit form of the same rendering: it returns the
#' \code{knitr_kable} object so the caller can pipe it into
#' further styling or embed it in a larger document. Both show the table
#' the console shows, formatted by \code{\link{format.dmar_tbl}} (for
#' most results the whole table; for a few, such as
#' \code{\link{fraction_missing_information}}, which leaves out the
#' columns that carry nothing for the call, and
#' \code{\link{factor_scores}}, which shows the first rows of each
#' factor, the same view the console prints), together with every line the
#' console writes with it: the sum-of-squares type, the confidence level,
#' how an omega's denominator was estimated and its largest residual
#' correlation, the bootstrap replication counts, the missing data a
#' sample size plan assumed, and the lines a result of its own adds (which
#' version of a fit measure was reported and why, what a model comparison
#' checked, how a fraction of missing information was estimated, which
#' rule set a table of degrees of freedom). What a
#' reader sees in a report therefore matches what they saw at the
#' console, and neither rounds the stored numbers.
#' \code{knit_print.dmar_tbl} sets each of those lines as its own
#' paragraph beneath the table (a line that introduces the table, as the
#' summary of \code{\link{factor_scores}}, the measurement structure of
#' \code{\link{fa_confirmatory}}, and the model summary of
#' \code{\link{mediation_mbco}} do, comes before it, and a small table
#' among the lines is set as a code block so its columns stay aligned),
#' and sets the whole off with a blank line before and after it, so that
#' two tables printed in one chunk stay apart;
#' \code{as_kable} joins them into the table's caption. Both write those
#' lines in the markup of the format \pkg{knitr} renders the table in,
#' escaped for LaTeX (a document knitted from \code{.Rnw}, or
#' \code{format = "latex"}) and for HTML, so that the percent sign of a
#' confidence level or an underscore in an item name prints as itself;
#' \pkg{knitr} escapes a table's cells and column names but writes a
#' caption as it is given.
#'
#' \code{results_sentence} turns a table that carries a confidence interval into the
#' one sentence an author puts in a results section, for example
#' \dQuote{smd = 0.50, 95\% CI [0.10, 0.90]}. It reads the estimate and its
#' limits from the numeric columns at full precision and formats them for the
#' sentence, so the reported numbers are exact to the requested decimals rather
#' than transcribed from the rounded console display. In a wide table the
#' estimate is the column the row names in an \code{estimand} column, as
#' the \code{glance()} row of a long table names the column of the
#' estimate its interval belongs to (see \code{\link{dmar_tidiers}}),
#' else the \code{estimate} (or \code{value}) column when there is
#' one, and otherwise the first numeric column that is not a limit, a
#' standard error, a test statistic with its degrees of freedom and
#' \emph{p}-value, or a count; a \emph{p}-value the table carries is
#' appended, to the four decimals of \code{\link{format_p}}, and one
#' below that floor as the bound it is (\dQuote{p < 0.0001}). In a long
#' table, whose interval is held in
#' \code{ci_lower} and \code{ci_upper} rows, the estimate is the row
#' the table names as the one the interval belongs to, when it names one
#' (the corrected correlation of \code{\link{correction_for_attenuation}},
#' whose table lists the observed correlation first, and the proportion
#' of superiority of \code{\link{proportion_of_superiority}}, listed after
#' the standardized mean difference it is computed from), and otherwise
#' its first row that is not a limit. Its \emph{p}-value is likewise the
#' row the table names as the \emph{p}-value of its test, when it names
#' one, and otherwise its first \emph{p}-value row: the equivalence tables
#' of \code{\link{equivalence_r}}, \code{\link{equivalence_smd}}, and
#' \code{\link{equivalence_c}} list the \emph{p}-values of the two
#' one-sided tests first and name the TOST \emph{p}-value, the test their
#' \eqn{1 - 2\alpha} interval decides. A long table may instead hold an
#' estimate's interval in rows of its own, named for the estimate:
#' \code{\link{sem_fit_measures}}
#' reports \code{rmsea} beside \code{rmsea_ci_lower} and
#' \code{rmsea_ci_upper}, and \code{gamma_hat} and
#' \code{adjusted_gamma_hat} the same way. The sentence then reports the
#' estimate named by \code{term}, or the first such estimate in the table
#' (the RMSEA there), with no \emph{p}-value, since the
#' \emph{p}-values of such a table belong to other tests (of exact, close,
#' and not close fit). \code{term} also reaches the RMSEA rows of a
#' \code{\link{fa_confirmatory}} table. The confidence level is the
#' interval's own level when the table records it (\code{fa_confirmatory()}
#' records the level of its RMSEA interval, its \code{rmsea_conf_level},
#' in the \code{rmsea_ci_level} row and in an attribute of the same name,
#' which is read when a selection of rows leaves that row out; a table
#' stacked with \code{rbind()} keeps the attributes of its first table
#' only, so keep the \code{rmsea_ci_level} row in each selection when
#' stacking tables fit at different levels), and
#' otherwise the table's
#' \code{conf_level} attribute, or its \code{conf_level} row or column (as
#' the reliability coefficients and the \code{tidy()} and \code{glance()}
#' verbs carry it).
#' The \code{tidy()} of a long table holds the level in a
#' \code{conf_level} column on the rows of its interval, the estimate and
#' the two limits, and \code{NA} on the others; the sentence reads the
#' level and the estimate from those rows, so that rows of several
#' \code{tidy()} results stacked with \code{rbind()} give the sentence of
#' their own table (a frame's attributes are those of the first frame of
#' a stack). Every DMAR table that forms an interval records its level, and a
#' sentence states no level the table does not record: a table that
#' records none, such as one assembled by hand, is refused unless
#' \code{conf_level} gives the level its interval was formed at, and a
#' \code{conf_level} that differs from the level the table records is
#' refused.
#'
#' A table whose two limits are not a confidence interval gets a sentence
#' that says so. \code{\link{vargha_delaney_A}}, \code{\link{cliff_delta}},
#' and \code{\link{probability_of_superiority_paired}} report both limits
#' equal to the estimate when the estimated variance is 0 (completely
#' separated groups, or no variation) and record the case in the table's
#' \code{interval_degenerate} attribute. The sentence then gives the
#' estimate and the reason and no interval, \dQuote{A = 1.00 (no
#' confidence interval: the groups are completely separated, so the
#' estimated variance is 0)} or \dQuote{A = 0.50 (no confidence interval:
#' every score is equal, so the estimated variance is 0)}, with
#' \dQuote{every pair favors the same member} and \dQuote{every pair is
#' tied} for the paired function. A \emph{p}-value is appended only when
#' the table has one (the permutation test of
#' \code{\link{vargha_delaney_A}}). The \code{tidy()} and \code{glance()}
#' rows of these tables carry the mark in a logical column
#' \code{interval_degenerate}, which says that the variance is 0 and not
#' which case gave it, so their sentence reads \dQuote{A = 1.00 (no
#' confidence interval: the estimated variance is 0)}.
#'
#' The level a sentence states is the interval's coverage, which for an
#' interval formed from tail areas \code{alpha_lower} and
#' \code{alpha_upper} is \eqn{1 - \alpha_{lower} - \alpha_{upper}}
#' whatever the split. A sentence in the form \dQuote{95\% CI [L, U]}
#' is read as the equal-tailed interval, so an interval whose tail areas
#' differ says so, \dQuote{smd = 0.50, 95\% CI [0.03, 0.85] (unequal
#' tail areas: 0.01 below, 0.04 above)}, and an interval with a tail
#' area of 0 is written as the one-sided confidence bound it is,
#' \dQuote{smd = 0.50, one-sided 95\% upper confidence bound 0.83},
#' since its other limit is infinite or the bound of the parameter (0 for
#' an RMSEA or a squared multiple correlation) and carries no
#' information. The split is read from the table's record of it, never
#' from a limit: an equal-tailed interval whose lower limit is 0 because
#' the statistic fell below its critical value (Steiger, 2004) is written
#' as the two-sided interval it is. The line beneath a printed table
#' states the same split, for example \dQuote{Confidence level: 95\%
#' (unequal tail areas: 0.01 below the lower limit, 0.04 above the upper
#' limit)}, and the \code{tidy()} and \code{glance()} rows of such an
#' interval carry it in \code{alpha_lower} and \code{alpha_upper}
#' columns (see \code{\link{dmar_tidiers}}).
#'
#' A confidence set that is not an interval is written as what it is. The
#' Fieller set of \code{\link{ci_ratio_of_means}} and
#' \code{\link{ci_ratio_of_means_paired}} can be the complement of an
#' interval or the whole line, and with unequal tail areas an interval
#' unbounded on one side or the union of two intervals. The sentence reads
#' the form from the table's \code{confidence_set} attribute (or the
#' \code{confidence_set} column of its \code{tidy()} row) and the ends of
#' the excluded gap from its \code{excluded_lower} and
#' \code{excluded_upper} rows (or columns), and writes, for example,
#' \dQuote{ratio = 0.32, 99.5\% confidence set: every value except those
#' between 0.90 and 147.05}, \dQuote{ratio = 0.32, 99.9\% confidence set:
#' the whole line, which excludes no value}, or \dQuote{ratio = 0.32, 97.8\%
#' confidence set: the union of [-0.61, 0.93] and [9.55, Inf) (unequal
#' tail areas: 0.02 below, 0.002 above)}, since \dQuote{99.5\% CI [-Inf,
#' Inf]} would read as the whole line and hide the gap. An infinite end is
#' written open and a finite one closed. A bounded Fieller interval is
#' written as any other interval.
#'
#' A table whose limits
#' are \code{NA}, because no interval was formed (as from
#' \code{\link{krippendorff_alpha}} without a bootstrap), has no sentence to
#' write and gives an error.
#'
#' @param x A \code{dmar_tbl} object (a \code{data.frame} returned by a DMAR
#'   function). For \code{results_sentence}, a \code{dmar_tbl} that carries a
#'   confidence interval, either as \code{ci_lower} / \code{ci_upper}
#'   rows of a long table or as \code{ci_lower} / \code{ci_upper}
#'   columns of a one-row wide
#'   table, such as one row of \code{\link{common_method_adjusted_r}} or
#'   of a \code{tidy()} result, or as \code{<estimate>_ci_lower} /
#'   \code{<estimate>_ci_upper} rows beside the estimate's row, as in a
#'   \code{\link{sem_fit_measures}} table, or as
#'   \code{<estimate>_ci_lower} / \code{<estimate>_ci_upper} columns
#'   beside the estimate's column of a one-row wide table, as in one row
#'   of \code{\link{ci_eta_squared}} or \code{\link{ci_omega_squared}}.
#' @param format Passed to \code{\link[knitr]{kable}} as its \code{format}
#'   argument (for example \code{"html"}, \code{"latex"}, \code{"pipe"}). The
#'   default \code{NULL} lets \pkg{knitr} choose based on the output context.
#' @param label For \code{results_sentence}, the label that leads the sentence
#'   (for example \code{"Cohen's d"}). The default \code{NULL} uses the name of
#'   the point-estimate term: the estimate's row of a long table, the
#'   quantity a \code{tidy()} row names in its \code{estimand} column
#'   (see \code{\link{dmar_tidiers}}), the row's \code{term} beside a
#'   wide table's \code{estimate} or \code{value} column, and otherwise
#'   the estimate column's name.
#' @param digits For \code{results_sentence}, the number of decimal places for
#'   the estimate and the interval limits. Default 2. Like the \code{digits}
#'   of the table display (see \code{\link{dmar_tbl}}), it counts decimal
#'   places, but it keeps a default of its own and does not follow
#'   \code{options(dmar.digits = )}, so a sentence written for a manuscript
#'   does not change with the console display. As in the table, a value
#'   that rounds to zero is written without a minus sign (0.00, not
#'   -0.00).
#' @param term For \code{results_sentence}, the estimate whose interval a
#'   table holds in \code{<term>_ci_lower} and \code{<term>_ci_upper}
#'   rows, for example \code{"gamma_hat"} in a
#'   \code{\link{sem_fit_measures}} table or \code{"rmsea"} in a
#'   \code{\link{fa_confirmatory}} table, or in columns of those names,
#'   for example \code{"eta_squared_total"} in a row of
#'   \code{\link{ci_eta_squared}}. The default \code{NULL} reads the
#'   table's interval as described above, taking the first estimate with
#'   such rows when the table holds its interval in no other form, and,
#'   of such columns, the one estimate whose limits are not \code{NA} (a
#'   row with several is refused unless \code{term} names one). Either
#'   way the interval for the standardized mean difference that
#'   \code{\link{cles}} and \code{\link{nnt_from_smd}} echo beside their
#'   own, which theirs are computed from, is taken only when \code{term}
#'   is \code{"smd"}.
#' @param conf_level For \code{results_sentence}, the confidence level of
#'   the table's interval, for a table that records none (one assembled
#'   by hand, say). The default \code{NULL} reads the level the table
#'   records and refuses a table that records none, rather than assuming
#'   one; a value that differs from the level the table records is
#'   refused.
#' @param ... Additional arguments. For \code{knit_print.dmar_tbl} and
#'   \code{as_kable.dmar_tbl}, passed to \code{\link[knitr]{kable}}
#'   (other than \code{caption}, which \code{as_kable} writes itself).
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
#' @references
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and
#'   contrast analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
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
#' # So does one row of a tidy() result, here the same interval.
#' results_sentence(generics::tidy(x))
#'
#' # A table that holds an interval in rows named for its estimate: the
#' # RMSEA of a model with chi square 85.3 on 24 df and N = 301, with its
#' # 90 percent interval (the default level of every RMSEA interval), and
#' # then gamma hat.
#' fm <- sem_fit_measures(chisq = 85.3, df = 24, N = 301, p = 9)
#' results_sentence(fm, label = "RMSEA", digits = 3)
#' results_sentence(fm, term = "gamma_hat", digits = 3)
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
  disp <- format(.dmar_tbl_display(x))
  tbl <- knitr::kable(disp, row.names = FALSE, ...)

  # The lines the console writes with the table (see .dmar_tbl_footer()),
  # one paragraph per line so Markdown does not run them together: those
  # that introduce the table before it, the rest beneath it. They are set
  # in the markup of the format knitr chose for the table, so in an Rnw
  # document (LaTeX) the "%" of "Confidence level: 95%" does not start a
  # comment and an underscore in an item name does not stop the
  # compilation.
  fmt <- attr(tbl, "format")
  parts <- .dmar_footer_parts(.dmar_tbl_footer(x))
  body <- c(.dmar_footer_markdown(parts$header, fmt),
            paste(tbl, collapse = "\n"),
            .dmar_footer_markdown(parts$body, fmt))
  # A blank line before and after the whole output, as knitr sets off a
  # kable of its own: knitr writes consecutive outputs of a chunk with
  # nothing between them, so without it two tables printed in one chunk
  # (with echo = FALSE) ran together, the second table's header row
  # appended to the first table's last footer line, which pandoc then read
  # as a header cell. Through 1.0.1 that was so.
  knitr::knit_print(knitr::asis_output(
    paste0("\n\n", paste(body, collapse = "\n\n"), "\n\n")))
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
    .dmar_refuse_argument(c(
      "Package 'knitr' is required for as_kable(). ",
      "Install with install.packages('knitr')."), arg = NULL)
  }
  # The caption is the table's own lines; a caption passed through `...`
  # would collide with it (knitr::kable() would see two), so it is refused
  # in words rather than by R's argument matching.
  if ("caption" %in% names(list(...))) {
    .dmar_refuse_argument(c(
      "as_kable() writes the caption itself, from the lines the console ",
      "prints with the table. For a caption of your own, call ",
      "knitr::kable() on format(x)."), arg = NULL)
  }
  # Reuse the console formatter so every numeric column (long or wide) is
  # rendered on its own terms, then hand the character table to knitr::kable.
  # The table is the one the console shows, and every line the console
  # writes with it goes into the caption.
  disp <- format(.dmar_tbl_display(x))
  tbl <- knitr::kable(disp, format = format, row.names = FALSE, ...)
  caption <- .dmar_footer_caption(.dmar_tbl_footer(x))
  if (is.null(caption)) return(tbl)
  # knitr escapes the cells and the column names but writes a caption as
  # it is given, so the caption is escaped here for the format knitr
  # renders, read off the table above (for format = NULL knitr chooses it
  # from the document being knitted). Through 1.0.1 it was not: in LaTeX
  # the "%" of "Confidence level: 95%" started a comment that swallowed
  # the caption's closing brace, so the document did not compile, and an
  # underscore in an item name stopped it as well.
  fmt <- attr(tbl, "format")
  knitr::kable(disp, format = fmt, row.names = FALSE,
               caption = .dmar_footer_escape(caption, fmt), ...)
}

#' Write a Publication-Style Results Sentence From a DMAR Result Table
#'
#' @rdname dmar_output_helpers
#' @export
results_sentence <- function(x, label = NULL, digits = 2, term = NULL,
                             conf_level = NULL) {
  if (!inherits(x, "dmar_tbl") && !is.data.frame(x)) {
    .dmar_refuse_argument(
      "'x' must be a dmar_tbl (a DMAR result table).", arg = "x")
  }
  if (!is.null(term) &&
      (!is.character(term) || length(term) != 1L || is.na(term) ||
       !nzchar(term))) {
    .dmar_refuse_argument(
      "'term' must be NULL or a single term of the table.", arg = "term")
  }
  if (!is.null(conf_level) &&
      (!is.numeric(conf_level) || length(conf_level) != 1L ||
       !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)) {
    .dmar_refuse_argument(
      "'conf_level' must be NULL or a single number in (0, 1).",
      arg = "conf_level")
  }
  given_level <- conf_level
  ignore_attr <- FALSE
  # A wide table that holds the interval of each of its estimates in
  # columns named for the estimate, <estimate>_ci_lower and
  # <estimate>_ci_upper (the total and the partial value of
  # ci_eta_squared() and ci_omega_squared()), is read as the table of the
  # one estimate `term` names, or of the one whose limits are not NA.
  named <- .results_sentence_named_interval(x, term)
  x <- named$x
  if (named$used_term) term <- NULL
  nms <- names(x)
  # A long table holds its interval in ci_lower / ci_upper rows of
  # `term`; a wide table holds it in the ci_lower / ci_upper pair of
  # columns, as the coefficient tables and the tidy() verbs do.
  long <- "term" %in% nms &&
    all(c("ci_lower", "ci_upper") %in% x$term)
  lim_cols <- .results_sentence_limit_cols(nms)
  # A table can also hold an estimate's interval in rows of its own,
  # <estimate>_ci_lower and <estimate>_ci_upper beside the estimate's row,
  # as sem_fit_measures() holds the RMSEA's (and gamma hat's) and fa_confirmatory()
  # the RMSEA's. Through 1.0.1 such a table had no sentence.
  row_est <- .results_sentence_row_estimates(x)
  if (!is.null(term) && !term %in% row_est) {
    .dmar_refuse_argument(c(
      "term = \"", term, "\" names no estimate whose interval this ",
      "table holds in ", term, "_ci_lower and ", term, "_ci_upper rows; ",
      if (length(row_est)) {
        paste0("the estimates that have such rows are ",
               paste0("\"", row_est, "\"", collapse = ", "), ".")
      } else {
        "the table has no such rows."
      }), arg = NULL)
  }
  rows_form <- !is.null(term) ||
    (!long && is.null(lim_cols) && length(row_est) > 0L)
  level_row <- NULL
  # The unequal tail areas of the interval, when its rows carry them in
  # alpha_lower and alpha_upper columns (a tidy() or glance() row), read
  # from the interval's row; NULL when the rows have no such columns, and
  # then the table's tail_areas attribute is read below.
  tails_row <- NULL
  tail_cols <- all(c("alpha_lower", "alpha_upper") %in% nms) &&
    is.numeric(x$alpha_lower) && is.numeric(x$alpha_upper)
  # A confidence set that need not be an interval (the Fieller set of the
  # ratio of means functions), read below from a table that records one;
  # NULL for every other table, whose sentence it leaves as it was.
  set <- NULL
  gap_terms <- c("excluded_lower", "excluded_upper")

  if (rows_form) {
    # The estimate is `term`, or else the first in the table's order that
    # has the two rows (the RMSEA of sem_fit_measures()), with the order
    # of .results_sentence_row_estimates(). A p-value the
    # table carries belongs to another test (exact, close, or not close
    # fit), so none is appended.
    est_term <- if (is.null(term)) row_est[1L] else term
    primary <- intersect(c("value", "estimate"),
                         names(x)[vapply(x, is.numeric, logical(1L))])
    if (!length(primary))
      .dmar_refuse_argument(
        "results_sentence() needs a numeric 'value' or 'estimate' column.",
        arg = NULL)
    primary <- primary[1L]
    row_of <- function(t) x[[primary]][x$term == t][1L]
    est <- row_of(est_term)
    lo  <- row_of(paste0(est_term, "_ci_lower"))
    hi  <- row_of(paste0(est_term, "_ci_upper"))
    est_label <- if (is.null(label)) est_term else label
    p_val <- NULL
    # The interval's own level, when a row records it (fa_confirmatory()'s
    # rmsea_ci_level, the level of its RMSEA interval), comes before the
    # table's.
    lv <- paste0(est_term, "_ci_level")
    if (lv %in% x$term) level_row <- row_of(lv)
    conf_col <- if ("conf_level" %in% x$term) row_of("conf_level")
    # The tidy() of such a table holds the level in a conf_level column on
    # the interval's rows, which comes before an attribute a stack keeps
    # from its first frame.
    if ("conf_level" %in% names(x) && is.numeric(x$conf_level)) {
      conf_col <- x$conf_level[x$term == paste0(est_term, "_ci_lower")][1L]
      ignore_attr <- TRUE
    }
    # A table whose RMSEA interval is at a level of its own also records
    # that level in its rmsea_conf_level attribute, which stays with the
    # table when rows are selected and the rmsea_ci_level row is left out.
    if (is.null(level_row) && !ignore_attr && identical(est_term, "rmsea")) {
      own <- attr(x, "rmsea_conf_level", exact = TRUE)
      if (is.numeric(own) && length(own) == 1L && is.finite(own))
        level_row <- own
    }
    if (tail_cols) {
      on <- x$term == paste0(est_term, "_ci_lower")
      tails_row <- c(x$alpha_lower[on][1L], x$alpha_upper[on][1L])
    }
  } else if (!long && !is.null(lim_cols)) {
    # Wide table: one row, the limits in `lim_cols`. The estimate is the
    # `estimate` (or `value`) column when there is one, and otherwise the
    # first numeric column that is not a limit, a standard error, a test
    # statistic, its degrees of freedom or p-value, or a count (the
    # primary column, e.g. eta_squared).
    if (nrow(x) != 1L) {
      .dmar_refuse_argument(c(
        "results_sentence() needs a single-row table; this table has ",
        nrow(x), " rows. Select one row first, e.g. x[x$effect == ..., ]",
        if (length(row_est)) {
          paste0(", or name the estimate whose interval the table holds ",
                 "in rows of its own, e.g. term = \"", row_est[1L], "\"")
        }, "."), arg = NULL)
    }
    num <- names(x)[vapply(x, is.numeric, logical(1L))]
    # A glance() row of a long table names the column that holds its
    # interval's estimate in an estimand column (the corrected correlation
    # of correction_for_attenuation(), listed after the observed one).
    estimand <- .results_sentence_estimand(x)
    est_col <- if (!is.null(estimand) && estimand %in% num &&
                   !estimand %in% lim_cols) {
      estimand
    } else {
      intersect(c("estimate", "value"), num)
    }
    if (!length(est_col)) {
      est_col <- num[!grepl(.results_sentence_not_estimate, num)]
    }
    if (!length(est_col))
      .dmar_refuse_argument(
        "Could not find a point-estimate column in this table.", arg = NULL)
    est_col <- est_col[1L]
    est <- x[[est_col]][1L]
    lo  <- x[[lim_cols[1L]]][1L]
    hi  <- x[[lim_cols[2L]]][1L]
    # The default label names the estimate: the quantity a tidy() row
    # names in its estimand column (the name of the column the default
    # tidy() took its estimate from, kappa or D2, where its term holds the
    # row number or the values of the label columns), else the row's term
    # when a generic estimate column sits beside one (a tidy() row, a
    # coefficient), and otherwise the column's own name.
    est_label <- if (!is.null(label)) {
      label
    } else if (!is.null(estimand) && est_col %in% c("estimate", "value")) {
      estimand
    } else if (est_col %in% c("estimate", "value") && "term" %in% nms &&
               !is.numeric(x$term)) {
      as.character(x$term[1L])
    } else {
      est_col
    }
    p_val <- .results_sentence_p_wide(x)
    conf_col <- if ("conf_level" %in% num) x$conf_level[1L]
    if (tail_cols) tails_row <- c(x$alpha_lower[1L], x$alpha_upper[1L])
    # The tidy() row of a ratio of means table names the form of its set
    # in a confidence_set column, beside the ends of its gap.
    if (all(gap_terms %in% nms) && is.numeric(x$excluded_lower) &&
        is.numeric(x$excluded_upper)) {
      set <- .results_sentence_set(x[["confidence_set"]],
                                   c(x$excluded_lower[1L],
                                     x$excluded_upper[1L]), lo, hi)
    }
  } else if (long) {
    # Long table: ci_lower / ci_upper are rows of `term`, beside a point-
    # estimate row, in the primary numeric column (`value` or `estimate`).
    primary <- intersect(c("value", "estimate"),
                         names(x)[vapply(x, is.numeric, logical(1L))])
    if (!length(primary))
      .dmar_refuse_argument(
        "results_sentence() needs a numeric 'value' or 'estimate' column.",
        arg = NULL)
    primary <- primary[1L]
    lo <- x[[primary]][x$term == "ci_lower"][1L]
    hi <- x[[primary]][x$term == "ci_upper"][1L]
    # The estimate is the row the table names as the one its interval
    # belongs to (the ci_estimate attribute), and otherwise its first row
    # that is not a limit. Through 1.0.1 it was always the first, so the
    # interval of correction_for_attenuation()'s corrected correlation was
    # reported beside the observed correlation listed before it. The
    # tidy() of a long table names it in its rows instead: a conf_level
    # column holds the level on the estimate and the two limits, and NA on
    # every other row, so that a sentence from rows of stacked tidy()
    # results reads the level and the estimate of those rows, and not the
    # attributes of the first frame of the stack.
    row_levels <- "conf_level" %in% names(x) && is.numeric(x$conf_level)
    if (row_levels) {
      lim_level <- x$conf_level[x$term == "ci_lower"][1L]
      not_est <- c("ci_lower", "ci_upper",
                   paste0(row_est, rep(c("", "_ci_lower", "_ci_upper"),
                                       each = length(row_est))))
      marked <- setdiff(x$term[!is.na(x$conf_level)], not_est)
      if (!is.na(lim_level) && !length(marked) &&
          length(setdiff(x$term, not_est))) {
        .dmar_refuse_argument(c(
          "The rows of this table carry the level of its interval but ",
          "not its estimate (no row but the limits has a conf_level); ",
          "keep the estimate's row with its limits."), arg = NULL)
      }
    }
    est_term <- if (row_levels && length(marked)) {
      marked[1L]
    } else {
      .results_sentence_long_estimate(x)
    }
    if (is.null(est_term))
      .dmar_refuse_argument(
        "Could not find a point-estimate row in this table.", arg = NULL)
    est <- x[[primary]][x$term == est_term][1L]
    est_label <- if (is.null(label)) est_term else label
    p_val <- .results_sentence_p_long(x, primary)
    conf_col <- if ("conf_level" %in% x$term) {
      x[[primary]][x$term == "conf_level"][1L]
    } else if (row_levels) {
      lim_level
    }
    # The rows' own level comes before an attribute, which a stack of
    # frames keeps from its first.
    if (row_levels) ignore_attr <- TRUE
    if (tail_cols) {
      on <- x$term == "ci_lower"
      tails_row <- c(x$alpha_lower[on][1L], x$alpha_upper[on][1L])
    }
    # A ratio of means table names the form of its set in its
    # confidence_set attribute and holds the ends of its gap in rows.
    if (all(gap_terms %in% x$term)) {
      set <- .results_sentence_set(
        if (!ignore_attr) attr(x, "confidence_set", exact = TRUE),
        c(x[[primary]][x$term == gap_terms[1L]][1L],
          x[[primary]][x$term == gap_terms[2L]][1L]), lo, hi)
    }
  } else {
    .dmar_refuse_argument(c(
      "results_sentence() reads a confidence interval from ",
      "ci_lower and ci_upper columns ",
      "of a one-row table, from ci_lower and ci_upper rows of ",
      "a term column, or from <estimate>_ci_lower and ",
      "<estimate>_ci_upper rows beside the estimate's row, and this ",
      "table has none of these. For an interval held under other ",
      "names, write the sentence from the columns directly."), arg = NULL)
  }
  # A table whose shape does not depend on its arguments keeps its limits
  # NA when no interval was formed (krippendorff_alpha() with
  # boot = FALSE); "95% CI [NA, NA]" is not a sentence anyone reports.
  if (is.na(lo) || is.na(hi)) {
    .dmar_refuse_argument(c(
      "This table's confidence interval limits are NA (no interval was ",
      "formed, for example because no bootstrap was run), so no ",
      "sentence can be written."), arg = NULL)
  }

  # The confidence level is the interval's own level row, when the table
  # has one, or else the table's attribute, or else its conf_level row
  # (the reliability family) or column (a tidy() row). Through 1.0.1 only
  # the attribute was read, so a reliability coefficient's 90% interval
  # was written as a 95% one. A table that records no level is refused
  # unless the caller gives it, since the sentence states no level the
  # table does not record; through 1.0.1 it was taken to be 95%, and the
  # 90% intervals of ci_eta_squared() and vargha_delaney_A(), which
  # recorded none, were written "95% CI".
  recorded <- .results_sentence_level(
    level_row,
    if (!ignore_attr) attr(x, "conf_level", exact = TRUE),
    conf_col)
  if (is.null(recorded)) {
    if (is.null(given_level)) {
      .dmar_refuse_argument(c(
        "This table records no confidence level (no conf_level ",
        "attribute, row, or column), so the sentence cannot state the ",
        "level of its interval. Give the level the interval was formed ",
        "at, e.g. results_sentence(x, conf_level = 0.95)."), arg = NULL)
    }
    conf_level <- given_level
  } else {
    if (!is.null(given_level) && abs(given_level - recorded) > 1e-8) {
      .dmar_refuse_argument(c(
        "conf_level = ", format(given_level), " is not the level this ",
        "table records for its interval, ", .format_conf_pct(recorded),
        "%. Leave conf_level as NULL to report the recorded level."),
        arg = NULL)
    }
    conf_level <- recorded
  }
  pct <- .format_conf_pct(conf_level)

  # An interval with unequal tail areas says so, since "95% CI [L, U]" is
  # read as the equal-tailed interval, and a one-sided interval (a tail
  # area of 0) is written as the confidence bound it is, whose other limit
  # (infinite, or the parameter's own bound, 0 for an RMSEA or an R2)
  # carries no information. The areas are those the rows carry, else the
  # table's attribute; one-sidedness is read from a tail area of 0, never
  # from a limit, since an equal-tailed interval can have a lower limit of
  # 0 because the statistic fell below its critical value (Steiger, 2004).
  # Through 1.0.1, and in the first versions of 1.1.0, a .01/.04 interval
  # and a one-sided one were both written "95% CI [L, U]".
  tails <- if (!is.null(tails_row)) {
    .tail_areas_valid(tails_row)
  } else if (!ignore_attr) {
    .recorded_tail_areas(x)
  }
  # A value that rounds to zero is written without a minus sign, as the
  # table writes it: a lower limit of -0.004 reads 0.00, where through
  # 1.0.x it read -0.00 (an infinite limit keeps its sign).
  fmt <- function(v) .dmar_fixed_text(v, digits)
  # A confidence set that is not an interval is written as what it is: a
  # complement of an interval as every value outside its gap, the whole
  # line as such, and a union of two intervals as its two pieces. So is an
  # interval with an infinite limit on a side whose tail area is not zero,
  # which is not the one-sided bound written below. Read as an interval,
  # the complement in Fieller's own example would be written "99.5% CI
  # [-Inf, Inf]", which reads as the whole line and hides the gap.
  # A table of an order-based effect size whose estimated variance is 0
  # (vargha_delaney_A(), cliff_delta(), probability_of_superiority_paired()
  # with completely separated groups or no variation) reports both limits
  # equal to the estimate, and they are not a confidence interval. The
  # sentence gives the estimate and says why there is no interval; through
  # the first versions of 1.1.0 it read "A = 1.00, 95% CI [1.00, 1.00]".
  # The table records the case in its interval_degenerate attribute, and
  # its tidy() and glance() rows in a logical column of that name.
  degenerate <- .results_sentence_degenerate(x, ignore_attr)
  out <- if (!is.null(degenerate)) {
    sprintf("%s = %s (no confidence interval: %s)",
            est_label, fmt(est), degenerate)
  } else if (.results_sentence_set_shown(set, lo, hi, tails)) {
    .results_sentence_set_text(set, est_label, est, pct, lo, hi, tails, fmt)
  } else if (is.null(tails)) {
    sprintf("%s = %s, %s%% CI [%s, %s]",
            est_label, fmt(est), pct, fmt(lo), fmt(hi))
  } else if (tails[[1L]] == 0) {
    sprintf("%s = %s, one-sided %s%% upper confidence bound %s",
            est_label, fmt(est), pct, fmt(hi))
  } else if (tails[[2L]] == 0) {
    sprintf("%s = %s, one-sided %s%% lower confidence bound %s",
            est_label, fmt(est), pct, fmt(lo))
  } else {
    sprintf("%s = %s, %s%% CI [%s, %s] (unequal tail areas: %s below, %s above)",
            est_label, fmt(est), pct, fmt(lo), fmt(hi),
            .format_tail(tails[[1L]]), .format_tail(tails[[2L]]))
  }
  # A p-value below format_p()'s floor is written as the bound it is,
  # "p < 0.0001"; through 1.0.1 the sentence read "p = < 0.0001".
  if (!is.null(p_val) && !is.na(p_val)) {
    p_txt <- format_p(p_val)
    p_label <- attr(p_val, "inference_label", exact = TRUE) %||% "p"
    out <- paste0(out, ", ", p_label, if (grepl("^[<>]", p_txt)) " " else " = ", p_txt)
  } else if (.results_sentence_simultaneous(x)) {
    out <- paste0(out, "; no matching multiplicity-adjusted p-value is available")
  }
  out <- paste0(out, .simultaneous_phrase(x))
  out
}

# ---- The lines that go with a table ------------------------------------
#
# .dmar_tbl_footer(x, ...) returns the lines every display of a dmar_tbl
# writes with the table, in display order, as a character vector named by
# the block each line belongs to. print.dmar_tbl() (and the print methods
# that share its engine, .dmar_tbl_write() in R/dmar_tbl.R),
# knit_print.dmar_tbl(), and as_kable.dmar_tbl() take their lines from here
# and nowhere else, so the console and a knitted report cannot drift apart.
#
# It is an internal S3 generic. The default method writes the lines any
# table can carry, read from its attributes: the sum-of-squares type, the
# confidence level, the t reference distribution of Wald intervals given
# a finite wald_df, how the total variance in an omega denominator was
# estimated, the largest residual correlation of a single-factor
# reliability fit (or that the fit is saturated), what the cells of a
# table of bifactor indices mean and how their intervals were formed,
# the bootstrap replication counts, the number of replications of a
# sensitivity analysis, a
# missing data treatment whose assumption the reader needs (pairwise
# deletion, multiple imputation), what a two-stage fit reports, the
# lines of a coefficient H table, the missing data a sample size plan
# assumed, the excess kurtosis of the errors an assurance plan allowed
# for, and how an assurance plan found the probability it reports. A
# class with lines of its own (dmar_sem_fit, dmar_sem_information,
# dmar_compare_sem, dmar_fmi, dmar_factor_scores, dmar_fa_confirmatory,
# dmar_fa_bifactor, dmar_fa_second_order, dmar_measurement_invariance,
# dmar_ave, dmar_mediation_mbco, dmar_df_missing, dmar_df_barnard_rubin,
# dmar_nnt, dmar_responder, dmar_vargha_delaney_A, dmar_cliff_delta,
# dmar_probability_of_superiority_paired, dmar_ss_aipe_cv_sensitivity,
# dmar_ss_aipe_partial_r_sensitivity, dmar_mr)
# has a method in its own file that keeps the text in a helper of its own
# and returns it with NextMethod()'s lines, through .dmar_footer_add(),
# and is registered with
# an @exportS3Method tag (S3method() in NAMESPACE registers a method of a
# generic the namespace defines without exporting either), so that the
# generic dispatches from any caller, a test or another package's ::: call
# included, and not only from inside the namespace. A generic was chosen
# over a footer attribute because several classes build their lines at
# display time from attributes (the FMI estimator, the rows a factor
# score summary shows), because a subset of a table that loses its
# attributes then falls back to the default lines without special code,
# and because the returned object carries nothing it did not carry before.
#
# The names. Consecutive lines with the same name form one block, which
# the console prints without blank lines between its lines; each block
# starts after a blank line, except that a "residual_cor" line stays with
# the "denominator" line before it and the "table_note" lines sit
# directly under the table, wherever a method put them. Lines named
# "header" introduce the table: the console prints them above it,
# knit_print() sets them as paragraphs before it, and as_kable() puts
# them first in the caption. A name ending in "_table"
# marks the lines of a small table laid out for the console; they belong
# to the block named without the suffix, knit_print() sets them as a
# code block so their columns stay aligned, and as_kable() closes
# up their spacing. A line may hold "\n" where the console wants a line
# break inside one paragraph, which Markdown and a caption read as a
# space. Neither generic nor any method is exported.
.dmar_tbl_footer <- function(x, ...) UseMethod(".dmar_tbl_footer")

#' @exportS3Method .dmar_tbl_footer default
#' @noRd
.dmar_tbl_footer.default <- function(x, ...) {
  lines <- character(0)
  # The record of a sensitivity run (.ss_sens_simulation() in
  # R/dmar_tidiers.R): the number of replications, and, when some returned
  # no interval (no test, for a power sensitivity table), the three counts,
  # since the summary rows describe only the replications with one.
  simulation <- attr(x, "simulation", exact = TRUE)
  if (is.list(simulation) && !is.null(simulation$attempted)) {
    lines["simulation"] <- if (isTRUE(simulation$failed > 0L)) {
      sprintf("Replications: %d attempted, %d with %s, %d without",
              simulation$attempted, simulation$returned,
              if (inherits(x, "dmar_ss_power_sensitivity")) "a test" else
                "an interval",
              simulation$failed)
    } else {
      sprintf("Replications: %d", simulation$attempted)
    }
  }
  ss_type <- attr(x, "ss_type", exact = TRUE)
  if (!is.null(ss_type)) {
    roman <- switch(as.character(ss_type),
                    "1" = "I", "2" = "II", "3" = "III", as.character(ss_type))
    lines["ss_type"] <- sprintf("Sum of squares: Type %s", roman)
  }
  conf_level <- attr(x, "conf_level", exact = TRUE)
  if (!is.null(conf_level))
    lines["conf_level"] <- paste0(
      sprintf("Confidence level: %s%%", .format_conf_pct(conf_level)),
      .simultaneous_phrase(x),
      .tail_phrase(.recorded_tail_areas(x), conf_level),
      .rmsea_level_phrase(x, conf_level))
  # Wald intervals referred to a t distribution (a finite wald_df
  # attribute: the reliability coefficients and fa_confirmatory()) say so; the line
  # is built by .wald_df_footer_line() in
  # R/missing_information_internals.R.
  wd_line <- .wald_df_footer_line(x)
  if (!is.null(wd_line)) lines["wald_df"] <- wd_line
  calibration <- attr(x, "wald_calibration_note", exact = TRUE)
  if (!is.null(calibration)) lines["wald_calibration"] <- calibration
  # An omega coefficient names how its denominator was estimated, so the
  # table says which coefficient it is without the reader having to know
  # the default, and a reliability fit names its largest residual
  # correlation, the local dependence that no denominator corrects (see
  # .denominator_phrase() and .residual_cor_summary() in
  # R/reliability_internals.R).
  denominator <- attr(x, "denominator", exact = TRUE)
  if (is.character(denominator) && length(denominator) == 1L) {
    lines["denominator"] <- sprintf(
      "Denominator: %s",
      .denominator_phrase(denominator,
                          attr(x, "coefficient", exact = TRUE),
                          attr(x, "metric", exact = TRUE)))
  }
  rc <- attr(x, "residual_cor", exact = TRUE)
  if (is.list(rc) && isTRUE(rc$saturated)) {
    lines["residual_cor"] <- paste0(
      "With three items the single-factor model is saturated (df = 0), ",
      "so its residual correlations are zero and carry no information ",
      "about local dependence")
  } else if (is.list(rc) && length(rc$largest) == 1L &&
             is.finite(rc$largest)) {
    lines["residual_cor"] <- sprintf(paste0(
      "Largest residual correlation under the single-factor model: ",
      "%s (%s with %s)"),
      formatC(rc$largest, format = "f", digits = 3),
      rc$items[1L], rc$items[2L])
  }
  # The bifactor indices say what ECV and omega mean for a group factor
  # and which H they report, with the maximal reliability of each factor,
  # and their intervals say how they were formed, over which rows, and
  # why PUC has none; the lines are built by .bifactor_indices_footer()
  # and .ci_bifactor_footer() in R/bifactor_internals.R and are empty for
  # every other table.
  lines <- c(lines, .bifactor_indices_footer(x), .ci_bifactor_footer(x))
  # A bootstrap interval says in one line how many replications it
  # conditions on, why the rest were set aside, and how many came from
  # an improper solution (the reliability family's `bootstrap`
  # attribute). The line is built by .bootstrap_footer_line() in
  # R/reliability_internals.R, beside the engine that records the
  # counts (.bootstrap_ci()).
  bs_line <- .bootstrap_footer_line(x)
  if (!is.null(bs_line)) lines["bootstrap"] <- bs_line
  # The counts of a bootstrap or Monte Carlo run (.dmar_boot_footer() in
  # R/dmar_bootstrap.R), for a table whose family has no bootstrap line of
  # its own: mediate(), the kappa functions, krippendorff_alpha(), htmt(),
  # R2_mixed_effects(), and ci_eta_squared_generalized(). The reliability
  # family and average_variance_extracted() print theirs from the
  # `bootstrap` attribute, and the classes named below add the line in
  # their own methods.
  if (is.null(attr(x, "bootstrap", exact = TRUE)) &&
      !inherits(x, c("dmar_ci_ratio_of_means", "dmar_mediation_mbco",
                     "dmar_fmi", "dmar_parallel_analysis",
                     "dmar_factor_count", "dmar_sem_fit",
                     "dmar_ci_bifactor"))) {
    # A planner that runs B replications inside each simulated study shows
    # the line only for an exploratory run.
    boot_line <- if (identical(attr(x, "bootstrap_scope", exact = TRUE),
                               "per_study")) {
      .dmar_boot_footer_per_study(x)
    } else {
      .dmar_boot_footer(x)
    }
    if (!is.null(boot_line)) lines["bootstrap_counts"] <- boot_line
  }
  effect_inference <- attr(x, "effect_size_inference", exact = TRUE)
  if (!is.null(effect_inference)) lines["effect_size_inference"] <- effect_inference
  h1_diagnostic <- attr(x, "h1_diagnostic", exact = TRUE)
  if (!is.null(h1_diagnostic)) lines["h1_diagnostic"] <- h1_diagnostic
  # Pairwise deletion and multiple imputation name their assumption, and
  # an imputation the counts behind its estimate and interval; the lines
  # are built by .missing_treatment_footer_lines() in
  # R/multiple_imputation_internals.R and are empty for every other
  # table, listwise deletion and full information maximum likelihood
  # included.
  mt_lines <- .missing_treatment_footer_lines(x)
  if (length(mt_lines)) lines[names(mt_lines)] <- mt_lines
  # A model estimated in two stages says which statistic, indices, and
  # information criteria it reports (fa_confirmatory() and measurement_invariance();
  # the lines are built by .ts_footer() in R/two_stage_internals.R).
  ts_lines <- .two_stage_footer_lines(x)
  if (length(ts_lines)) lines[names(ts_lines)] <- ts_lines
  # A table of coefficient H adds its own lines (which factor of a fitted
  # model, the latent response metric, the standard error's basis, the
  # bootstrap bias); see .H_footer_lines() in R/reliability_internals.R.
  lines <- c(lines, .H_footer_lines(x))
  # A sample size plan made for incomplete data names the fraction of
  # missing information it assumed and the complete data N beside it, and
  # a sensitivity analysis names the missingness it simulated (see
  # .aipe_missing_phrase() in R/ss_aipe_missing_internals.R). Only a table
  # that carries those attributes under their exact names is asked: attr()
  # matches names partially, so a fraction_missing_information() table,
  # whose fmi_max attribute "fmi" would otherwise match, is left out.
  if (!is.null(attr(x, "fmi", exact = TRUE)) ||
      !is.null(attr(x, "missing_prop", exact = TRUE))) {
    missing_line <- .aipe_missing_phrase(x)
    if (!is.null(missing_line)) lines["missing_data"] <- missing_line
  }
  # An assurance plan of ss_aipe_reliability() whose final Monte Carlo pass
  # had samples that returned no interval, counted as wide, says how many
  # (see .ss_aipe_no_interval_line() in R/ss_aipe_reliability.R).
  if (!is.null(attr(x, "no_interval", exact = TRUE))) {
    ni_line <- .ss_aipe_no_interval_line(x)
    if (!is.null(ni_line)) lines["no_interval"] <- ni_line
  }
  # An assurance plan of ss_aipe_c() made for errors that are not normal
  # names the excess kurtosis its assurance allows for (see
  # .ss_aipe_kurtosis_line() in R/ss_aipe_c.R).
  if (!is.null(attr(x, "true_kurtosis", exact = TRUE))) {
    k_line <- .ss_aipe_kurtosis_line(x)
    if (!is.null(k_line)) lines["kurtosis"] <- k_line
  }
  # An assurance plan of ss_aipe_smd(), ss_aipe_sc(), ss_aipe_sc_ancova(),
  # ss_aipe_reg_coef(), ss_aipe_rc(), or ss_aipe_src() names its
  # assurance_method and says how assurance_attained was found, or that
  # the published approximation does not compute it (see
  # .ss_aipe_assurance_line() in R/ss_aipe_assurance_internals.R).
  if (!is.null(attr(x, "assurance_method", exact = TRUE))) {
    a_line <- .ss_aipe_assurance_line(x)
    if (!is.null(a_line)) lines["assurance"] <- a_line
  }
  lines
}

# .dmar_tbl_display(x, ...) returns the table every display shows: by
# default the table itself. A class whose display leaves out columns or
# rows returns that view instead (dmar_fmi drops the columns that carry
# nothing for the call; dmar_factor_scores shows the first rows of each
# factor), so a knitted table has the columns and rows the console shows.
# The value is a dmar_tbl, whose format() method makes the display
# strings. An internal S3 generic like .dmar_tbl_footer(); not exported.
.dmar_tbl_display <- function(x, ...) UseMethod(".dmar_tbl_display")

#' @exportS3Method .dmar_tbl_display default
#' @noRd
.dmar_tbl_display.default <- function(x, ...) x

# The level of an RMSEA interval held in the rows of a table whose other
# intervals are at another level (the fit rows of fa_confirmatory(),
# fa_bifactor(), and fa_second_order(), whose rmsea_conf_level attribute
# records it), for the line that states the table's level: " (90% for the
# RMSEA interval)", or "" when the table records no such level, holds no
# RMSEA limits, or forms them at its own level.
.rmsea_level_phrase <- function(x, conf_level) {
  level <- attr(x, "rmsea_conf_level", exact = TRUE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      isTRUE(all.equal(level, as.numeric(conf_level)))) return("")
  col <- intersect(c("estimate", "value"), names(x))
  if (!length(col) || is.null(x[["term"]])) return("")
  limits <- x[[col[1L]]][x[["term"]] %in% c("rmsea_ci_lower", "rmsea_ci_upper")]
  if (!length(limits) || !any(is.finite(limits))) return("")
  sprintf(" (%s%% for the RMSEA interval)", .format_conf_pct(level))
}

# Appends a class's lines (a character vector, possibly NULL or empty) to
# the lines from NextMethod(), all named `block`. Used by the class methods
# of .dmar_tbl_footer().
.dmar_footer_add <- function(lines, new, block) {
  if (!length(new)) return(lines)
  c(lines, stats::setNames(as.character(new), rep(block, length(new))))
}

# Splits footer lines into those that introduce the table ("header") and
# those beneath it, the "table_note" lines first. Returns a list with
# elements `header` and `body`.
.dmar_footer_parts <- function(lines) {
  kind <- names(lines) %||% rep("", length(lines))
  head <- kind == "header"
  note <- kind == "table_note"
  list(header = lines[head], body = c(lines[note], lines[!head & !note]))
}

# The block each footer line belongs to: its name, without a "_table"
# suffix.
.dmar_footer_block <- function(lines) {
  sub("_table$", "", names(lines) %||% rep("", length(lines)))
}

# Whether each line starts a new block, which the console sets off with a
# blank line: a line starts one unless it continues the block of the line
# before it, is a "residual_cor" line under a "denominator" line, or is a
# "table_note" directly under the table.
.dmar_footer_breaks <- function(lines) {
  block <- .dmar_footer_block(lines)
  n <- length(block)
  if (!n) return(logical(0))
  prev <- c(NA_character_, block[-n])
  brk <- is.na(prev) | block != prev
  brk[block == "residual_cor" & prev %in% "denominator"] <- FALSE
  brk[1L] <- block[1L] != "table_note"
  brk
}

# Writes footer lines at the console, a blank line before each block.
.dmar_footer_cat <- function(lines) {
  brk <- .dmar_footer_breaks(lines)
  for (k in seq_along(lines)) {
    cat(if (brk[k]) "\n", lines[[k]], "\n", sep = "")
  }
  invisible(lines)
}

# The markup knit_print() writes for footer lines: one paragraph per
# line, and each run of "_table" lines as one block of fixed-width text.
# `format` is the format of the kable the lines go with: in Markdown (the
# default, and every format but "latex" and "html") a paragraph is the
# line itself and a block is a fenced code block; in LaTeX (an Rnw
# document) a paragraph is the line escaped for text mode and a block a
# verbatim environment; in HTML a paragraph is a <p> element and a block a
# <pre> element, their text escaped. Returns a character vector of chunks
# for the caller to separate with blank lines.
.dmar_footer_markdown <- function(lines, format = "pipe") {
  if (!length(lines)) return(character(0))
  format <- if (identical(format, "latex") || identical(format, "html")) {
    format
  } else {
    "markdown"
  }
  fence <- switch(format,
                  latex    = c("\\begin{verbatim}", "\\end{verbatim}"),
                  html     = c("<pre>", "</pre>"),
                  markdown = c("```", "```"))
  is_table <- grepl("_table$", names(lines))
  run <- cumsum(c(TRUE, is_table[-1L] != is_table[-length(is_table)]))
  out <- character(0)
  for (r in unique(run)) {
    k <- which(run == r)
    if (is_table[k[1L]]) {
      txt <- unname(lines[k])
      if (format == "html") txt <- .escape_html(txt)
      out <- c(out, paste(c(fence[1L], txt, fence[2L]), collapse = "\n"))
    } else {
      txt <- .dmar_footer_escape(trimws(unname(lines[k]), which = "left"),
                                 format)
      if (format == "html") txt <- paste0("<p>", txt, "</p>")
      out <- c(out, txt)
    }
  }
  out
}

# Escapes footer text (a caption, or the paragraphs knit_print() writes)
# for the format of the kable it goes with: text-mode LaTeX through
# .escape_latex() and HTML through .escape_html() (both in
# R/correlations_test.R), and any other format unchanged, since Markdown
# reads the lines as they are.
.dmar_footer_escape <- function(txt, format) {
  if (identical(format, "latex")) return(.escape_latex(txt))
  if (identical(format, "html")) return(.escape_html(txt))
  txt
}

# The caption as_kable() writes from footer lines (NULL when there are
# none): the lines in order, the lines that introduce the table first,
# each line's spacing closed up (a caption is one line of text), joined
# with "; ", or with a space after a line that ends a sentence or
# introduces what follows.
.dmar_footer_caption <- function(lines) {
  if (!length(lines)) return(NULL)
  parts <- .dmar_footer_parts(lines)
  lines <- c(parts$header, parts$body)
  txt <- trimws(gsub("\\s+", " ", unname(lines)))
  txt <- txt[nzchar(txt)]
  if (!length(txt)) return(NULL)
  sep <- ifelse(grepl("[.:]$", txt), " ", "; ")
  paste0(paste0(txt[-length(txt)], sep[-length(txt)], collapse = ""),
         txt[length(txt)])
}

# Format a confidence level (a proportion such as 0.95) as a percent for
# display, with at least `digits` significant digits and no trailing ".0"
# or width padding (format(scientific = FALSE, trim = TRUE) does both), so
# 0.95 reads "95" and 0.975 reads "97.5". A level below 1 is never shown
# as 100: digits are added while the rounded display would read 100 (the
# level of each of 28 Bonferroni intervals at a family level of .999,
# 1 - .001 / 28 = .99996429, reads "99.996", where four significant digits
# gave "100"), up to the 15 a double holds, so only a level within about
# 1e-15 of 1 can still read 100. The digits are given here, not taken
# from getOption("digits"), under which options(digits = 2) wrote a 97.5%
# interval as "98%". The text shown follows getOption("OutDec"), as the
# numbers of the table and of its sentence do ("97,5" under OutDec = ","),
# so whether the display reads 100 is decided on a copy written with "."
# as the decimal mark, which as.numeric() reads under any OutDec (a copy
# with "," read as NA and stopped the loop). Display only; the stored
# level is never rounded. Not exported.
.format_conf_pct <- function(conf_level, digits = 7L) {
  .format_below_top(100 * as.numeric(conf_level), top = 100, digits = digits)
}

# The same display of a confidence level kept as a proportion, for output
# that has always written it so ("level = 0.975" in summary() of mlmr()):
# at least `digits` significant digits whatever getOption("digits") is,
# and never "1" for a level below 1. Not exported.
.format_conf_prop <- function(conf_level, digits = 7L) {
  .format_below_top(as.numeric(conf_level), top = 1, digits = digits)
}

# The engine of .format_conf_pct() and .format_conf_prop(): `x` written
# with at least `digits` significant digits, more while a value below
# `top` would read `top`, up to 15. Not exported.
.format_below_top <- function(x, top, digits) {
  vapply(x, function(p) {
    d <- as.integer(digits)
    reads_top <- function(d) {
      as.numeric(format(p, digits = d, scientific = FALSE, trim = TRUE,
                        decimal.mark = ".")) >= top
    }
    while (is.finite(p) && p < top && d < 15L && reads_top(d)) {
      d <- d + 1L
    }
    format(p, digits = d, scientific = FALSE, trim = TRUE)
  }, character(1L))
}

# The column labels of a two-sided interval at `conf_level` from
# confint(): the percentiles of its two limits with a space before the
# percent sign, written as stats::confint.default() writes them (three
# significant digits for the smaller of the two tail percentiles, which
# sets the decimals of both, so 0.975 reads "1.25 %" and "98.75 %"), not
# with getOption("digits"), under which options(digits = 2) labeled them
# "1.3 %" and "99 %". Display only. Not exported.
.ci_quantile_labels <- function(conf_level) {
  a <- (1 - conf_level) / 2
  paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE,
               digits = 3L), "%")
}

# What the confidence level line adds for intervals that hold jointly at
# the recorded level by an adjustment of each one (the `simultaneous`
# attribute of contrast_test(adjust = "bonferroni") and of
# pairwise_within(bonferroni_ci = TRUE): the method, the number of
# comparisons, and the level of each interval, shown to at least four
# significant digits, 98.33 rather than 98.33333, and to more when four
# would round a level below 100% to 100 (99.996 for 28 intervals at a
# family level of .999); "" for any other table. Not exported.
.simultaneous_phrase <- function(x) {
  s <- attr(x, "simultaneous", exact = TRUE)
  if (!is.list(s) || is.null(s$method) || is.null(s$n_comparisons)) {
    method <- .results_sentence_adjustment(x, "ci_adjustment")
    labels <- c(bonferroni = "Bonferroni", tukey = "Tukey",
                 scheffe = "Scheffe", dunnett = "Dunnett", games_howell = "Games-Howell",
                 dunnett_t3 = "Dunnett's T3")
    if (is.null(method) || !method %in% names(labels)) return("")
    return(paste0(", simultaneous interval (", labels[[method]], ")"))
  }
  each <- if (is.numeric(s$each) && length(s$each) == 1L) {
    sprintf(": each interval at %s%%", .format_conf_pct(s$each, digits = 4L))
  } else ""
  sprintf(", simultaneous across the %d comparisons (%s%s)",
          as.integer(s$n_comparisons), s$method, each)
}

# The tail areas of an interval formed from alpha_lower and alpha_upper,
# c(lower = alpha_lower, upper = alpha_upper), when the two differ (by
# more than 1e-12), and NULL when they are equal or either is not a single
# finite number. The confidence coefficient of an interval [L, U] with
# tail areas a_L below L and a_U above U is 1 - a_L - a_U whatever the
# split (Neyman, 1937; Casella & Berger, 2002, sec. 9.2), and that is the
# level a table records; the split is recorded beside it only when it is
# unequal, because "95% CI [L, U]" is read as the equal-tailed interval,
# and an interval with a tail area of 0 is a one-sided confidence bound
# whose other limit (an infinite one, or the parameter's own bound, such
# as 0 for an RMSEA or a squared multiple correlation) carries no
# information. Each constructor that takes alpha_lower and alpha_upper
# stores the result in a tail_areas attribute. Not exported.
.tail_areas <- function(alpha_lower, alpha_upper) {
  if (!is.numeric(alpha_lower) || !is.numeric(alpha_upper) ||
      length(alpha_lower) != 1L || length(alpha_upper) != 1L ||
      !is.finite(alpha_lower) || !is.finite(alpha_upper) ||
      abs(alpha_lower - alpha_upper) <= 1e-12) {
    return(NULL)
  }
  c(lower = alpha_lower, upper = alpha_upper)
}

# The unequal tail areas a table records (its tail_areas attribute), or
# NULL. Not exported.
.recorded_tail_areas <- function(x) {
  .tail_areas_valid(attr(x, "tail_areas", exact = TRUE))
}

# `tails` when it is a pair of finite tail areas that differ, else NULL.
# Not exported.
.tail_areas_valid <- function(tails) {
  if (!is.numeric(tails) || length(tails) != 2L) return(NULL)
  .tail_areas(unname(tails[1L]), unname(tails[2L]))
}

# A tail area for display, to four significant digits (0.01, 0.025,
# 0.01667), never rounded to 0 when it is positive, and not dependent on
# getOption("digits"). Not exported.
.format_tail <- function(a) {
  format(a, digits = 4L, scientific = FALSE, trim = TRUE)
}

# What the confidence level line adds for an interval with unequal tail
# areas (`tails`, from .recorded_tail_areas()): the two areas, or, when
# one is 0, that the interval is one-sided and which limit is the
# confidence bound at the level; "" when the tails are equal or not
# recorded. Not exported.
.tail_phrase <- function(tails, conf_level) {
  if (is.null(tails)) return("")
  lo <- unname(tails[1L])
  hi <- unname(tails[2L])
  pct <- .format_conf_pct(conf_level)
  if (lo == 0) {
    sprintf(paste0(" (one-sided: the upper limit is the %s%% upper ",
                   "confidence bound, with tail area %s above it and ",
                   "none below the lower limit)"), pct, .format_tail(hi))
  } else if (hi == 0) {
    sprintf(paste0(" (one-sided: the lower limit is the %s%% lower ",
                   "confidence bound, with tail area %s below it and ",
                   "none above the upper limit)"), pct, .format_tail(lo))
  } else {
    sprintf(paste0(" (unequal tail areas: %s below the lower limit, %s ",
                   "above the upper limit)"),
            .format_tail(lo), .format_tail(hi))
  }
}

# The level a table records for its interval: the interval's own level
# row (`level_row`), else the table's conf_level attribute (`attr_level`),
# else its conf_level row or column (`col_level`), each taken only when it
# is a single finite number; NULL when the table records none. Not
# exported.
.results_sentence_level <- function(level_row, attr_level, col_level) {
  for (v in list(level_row, attr_level, col_level)) {
    if (is.numeric(v) && length(v) >= 1L && is.finite(v[1L])) return(v[1L])
  }
  NULL
}

# The estimate a long table's ci_lower and ci_upper rows belong to:
# the term the table names in its ci_estimate attribute, when that is one
# of its rows, and otherwise its first row that is not a limit; NULL when
# there is none. A function whose first row is not the interval's
# estimate names it (correction_for_attenuation(), whose table lists the
# observed correlation before the corrected one the interval is for;
# proportion_of_superiority(), which lists the standardized mean
# difference first). Not exported.
.results_sentence_long_estimate <- function(x) {
  terms <- as.character(x$term)
  own <- attr(x, "ci_estimate", exact = TRUE)
  if (is.character(own) && length(own) == 1L && own %in% terms) return(own)
  rest <- setdiff(terms, c("ci_lower", "ci_upper"))
  if (length(rest)) rest[1L] else NULL
}

# The confidence set of a table whose set need not be an interval, as the
# Fieller set of ci_ratio_of_means() and ci_ratio_of_means_paired() need
# not: list(shape, gap), with `shape` the form the table names (its
# confidence_set attribute, or the confidence_set column of its tidy()
# row) and `gap` the ends of its excluded gap. A table that keeps the gap
# but has lost the attribute, as a subset of rows does, has its form read
# from the gap and the limits `lo` and `hi` instead. Not exported.
.results_sentence_set <- function(shape, gap, lo, hi) {
  if (is.factor(shape)) shape <- as.character(shape)
  gap <- as.numeric(gap)
  if (!is.character(shape) || length(shape) < 1L || is.na(shape[1L])) {
    unbounded <- identical(c(lo, hi), c(-Inf, Inf))
    shape <- if (!anyNA(gap)) {
      if (unbounded) "complement of an interval" else "union of two intervals"
    } else if (unbounded) {
      "whole line"
    } else {
      "interval"
    }
  }
  list(shape = shape[1L], gap = gap)
}

# Whether results_sentence() writes the set `set` (from
# .results_sentence_set()) in its own form: always when it is not an
# interval, and for an interval when a limit is infinite on a side whose
# tail area is not zero (an interval with a one-sided bound is written as
# that bound). `tails` are the recorded unequal tail areas, NULL when they
# are equal and so both positive. Not exported.
.results_sentence_set_shown <- function(set, lo, hi, tails) {
  if (is.null(set)) return(FALSE)
  if (!identical(set$shape, "interval")) return(TRUE)
  a <- if (is.null(tails)) c(1, 1) else unname(tails)
  (a[1L] > 0 && is.infinite(lo)) || (a[2L] > 0 && is.infinite(hi))
}

# The sentence for such a set. An infinite end is written open, "(-Inf"
# or "Inf)", and a finite one closed, since the set contains its finite
# ends; the gap of a complement is written in words. Unequal tail areas
# are stated as results_sentence() states them for an interval, and a
# tail area of 0 as one-sided. Not exported.
.results_sentence_set_text <- function(set, label, est, pct, lo, hi, tails,
                                       fmt) {
  end <- function(v, left) {
    if (is.infinite(v)) return(if (left) "(-Inf" else "Inf)")
    if (left) paste0("[", fmt(v)) else paste0(fmt(v), "]")
  }
  piece <- function(a, b) paste0(end(a, TRUE), ", ", end(b, FALSE))
  gap <- set$gap
  what <- switch(set$shape,
    "complement of an interval" = sprintf(
      "confidence set: every value except those between %s and %s",
      fmt(gap[1L]), fmt(gap[2L])),
    "whole line" = "confidence set: the whole line, which excludes no value",
    "union of two intervals" = sprintf("confidence set: the union of %s and %s",
                                       piece(lo, gap[1L]), piece(gap[2L], hi)),
    paste("CI", piece(lo, hi)))
  out <- sprintf("%s = %s, %s%% %s", label, fmt(est), pct, what)
  if (is.null(tails)) return(out)
  a <- unname(tails)
  paste0(out, if (a[1L] == 0) {
    sprintf(" (one-sided: tail area %s above, none below)", .format_tail(a[2L]))
  } else if (a[2L] == 0) {
    sprintf(" (one-sided: tail area %s below, none above)", .format_tail(a[1L]))
  } else {
    sprintf(" (unequal tail areas: %s below, %s above)",
            .format_tail(a[1L]), .format_tail(a[2L]))
  })
}

# The quantity a one-row wide table names in an estimand column, a single
# non-missing string: the name of the column the default tidy() took the
# row's estimate from, or the column of a glance() row that holds the
# estimate of its interval; NULL when it names none. Not exported.
.results_sentence_estimand <- function(x) {
  e <- x[["estimand"]]
  if (is.factor(e)) e <- as.character(e)
  if (is.character(e) && length(e) >= 1L && !is.na(e[1L]) &&
      nzchar(e[1L])) e[1L] else NULL
}

# A wide table without ci_lower and ci_upper columns whose estimates hold
# their limits in columns named for them, <estimate>_ci_lower and
# <estimate>_ci_upper, reduced to the columns of one estimate, its limits
# renamed ci_lower and ci_upper and its name in an estimand column, which
# labels the sentence. The estimate is `term` when it names one (used_term
# is then TRUE), and otherwise the one estimate whose limits are not all
# NA; a table with several such estimates is refused unless `term` names
# one. An estimate the table names in its ci_source attribute (the
# standardized mean difference in the glance() row of cles() and
# nnt_from_smd(), whose interval the others are computed from) is set
# aside unless `term` names it or no other estimate has an interval. Any
# other table is returned as it is. Not exported.
.results_sentence_named_interval <- function(x, term) {
  nms <- names(x)
  unchanged <- list(x = x, used_term = FALSE)
  if (all(c("ci_lower", "ci_upper") %in% nms) || "term" %in% nms)
    return(unchanged)
  stems <- sub("_ci_lower$", "", nms[grepl("_ci_lower$", nms)])
  stems <- stems[stems %in% nms & paste0(stems, "_ci_upper") %in% nms]
  stems <- stems[vapply(stems, function(e) is.numeric(x[[e]]), logical(1L))]
  if (!length(stems)) return(unchanged)
  src <- attr(x, "ci_source", exact = TRUE)
  if (is.character(src))
    stems <- c(setdiff(stems, src), intersect(stems, src))
  used_term <- FALSE
  if (!is.null(term)) {
    if (!term %in% stems) return(unchanged)
    pick <- term
    used_term <- TRUE
  } else {
    formed <- stems[vapply(stems, function(e)
      any(!is.na(x[[paste0(e, "_ci_lower")]]) &
            !is.na(x[[paste0(e, "_ci_upper")]])), logical(1L))]
    if (is.character(src) && any(!formed %in% src))
      formed <- setdiff(formed, src)
    if (length(formed) > 1L) {
      .dmar_refuse_argument(c(
        "This table holds the intervals of several estimates (",
        paste0("\"", formed, "\"", collapse = ", "), "); name the one ",
        "for the sentence, e.g. term = \"", formed[length(formed)], "\"."),
        arg = "term")
    }
    pick <- if (length(formed)) formed else stems[1L]
  }
  others <- unlist(lapply(setdiff(stems, pick), function(e)
    paste0(e, c("", "_ci_lower", "_ci_upper"))))
  y <- x[, setdiff(nms, others), drop = FALSE]
  names(y)[names(y) == paste0(pick, "_ci_lower")] <- "ci_lower"
  names(y)[names(y) == paste0(pick, "_ci_upper")] <- "ci_upper"
  y$estimand <- rep(pick, nrow(y))
  for (a in c("conf_level", "tail_areas"))
    attr(y, a) <- attr(x, a, exact = TRUE)
  list(x = y, used_term = used_term)
}

# The pair of columns that holds a wide table's interval, ci_lower and
# ci_upper; NULL when the table does not have both. Not exported.
.results_sentence_limit_cols <- function(nms) {
  pair <- c("ci_lower", "ci_upper")
  if (all(pair %in% nms)) pair else NULL
}

# The terms of a table that hold their interval in rows of their own,
# <term>_ci_lower and <term>_ci_upper beside the term's row (the rmsea,
# gamma_hat, and adjusted_gamma_hat of sem_fit_measures(), the rmsea of
# fa_confirmatory()), in the table's order, with two exceptions. The
# estimate a table names in its ci_source attribute comes last: cles(),
# nnt_from_smd(), and proportion_of_superiority() echo the interval for
# the standardized mean difference their own intervals are computed from
# (smd_ci_lower, smd_ci_upper) ahead of those intervals, and name it
# there. And in rows that carry their level in a conf_level column (a
# tidy() result, or a stack or subset of such results), the estimates
# whose limit rows carry a level come first, since tidy() gives the
# ci_source interval none and a subset of rows loses the attribute.
# character(0) when there are none or the table has no term column. Not
# exported.
.results_sentence_row_estimates <- function(x) {
  if (!"term" %in% names(x) || is.numeric(x$term)) return(character(0))
  terms <- as.character(x$term)
  est <- terms[paste0(terms, "_ci_lower") %in% terms &
                 paste0(terms, "_ci_upper") %in% terms]
  src <- attr(x, "ci_source", exact = TRUE)
  if (is.character(src)) est <- c(setdiff(est, src), intersect(est, src))
  if (length(est) && is.numeric(x[["conf_level"]])) {
    marked <- !is.na(x$conf_level[match(paste0(est, "_ci_lower"), terms)])
    est <- c(est[marked], est[!marked])
  }
  est
}

# The names of the numeric columns of a wide table that are never its
# point estimate: the limits, a standard error, a test statistic with its
# degrees of freedom and p-value, a confidence level and its tail areas,
# and a count or sample size. Not exported.
.results_sentence_not_estimate <- paste0(
  "^(ci_lower|ci_upper|conf_level|alpha_lower|",
  "alpha_upper|se|statistic|",
  "(t|z|F|q)_(value|statistic)|chi_square|p_value|p\\.value|p_adjusted|",
  "df|df_.+|N|N_.+|n|n_.+)$")

# Why a table of an order-based effect size has no confidence interval,
# in the words results_sentence() writes, or NULL when it has one. The
# table's interval_degenerate attribute names the case
# ("complete_separation" or "no_variation"), and the words are those of
# the function that made the table, pairs for
# probability_of_superiority_paired() and groups for the other two. Rows
# that carry a logical interval_degenerate column (a tidy() or glance()
# result, or a stack of them) are read from the column: it says that the
# estimated variance is 0, not which case gave it. `ignore_attr` is TRUE
# when the rows carry their own level, so that the attribute of the first
# frame of a stack is not read for the rows of another. Not exported.
.results_sentence_degenerate <- function(x, ignore_attr = FALSE) {
  mark <- x[["interval_degenerate"]]
  if (is.logical(mark)) {
    on <- if ("term" %in% names(x) && "ci_lower" %in% x$term) {
      mark[x$term == "ci_lower"][1L]
    } else {
      mark[1L]
    }
    return(if (isTRUE(on)) "the estimated variance is 0")
  }
  reason <- if (!ignore_attr) attr(x, "interval_degenerate", exact = TRUE)
  if (!is.character(reason) || length(reason) != 1L) return(NULL)
  paired <- inherits(x, "dmar_probability_of_superiority_paired")
  what <- switch(reason,
    complete_separation = if (paired) {
      "every pair favors the same member"
    } else {
      "the groups are completely separated"
    },
    no_variation = if (paired) "every pair is tied" else
      "every score is equal",
    NULL)
  if (is.null(what)) return(NULL)
  paste0(what, ", so the estimated variance is 0")
}

# Pull a p-value out of a wide interval-carrying table, if one is present:
# the column its p_sentence attribute names, when the row has it (the TOST
# p-value of the glance() row of an equivalence table), and otherwise the
# first column with one of the names the format method uses. Not exported.
.results_sentence_p_wide <- function(x) {
  simultaneous <- .results_sentence_simultaneous(x)
  if (simultaneous) {
    if (!.results_sentence_matching_adjustment(x) || !"p_adjusted" %in% names(x)) return(NULL)
    return(.results_sentence_named_p(x$p_adjusted[1L], "multiplicity-adjusted p"))
  }
  own <- attr(x, "p_sentence", exact = TRUE)
  if (is.character(own) && length(own) == 1L && own %in% names(x) &&
      is.numeric(x[[own]]) && !is.na(x[[own]][1L])) {
    return(.results_sentence_named_p(x[[own]][1L], if (identical(own, "p_tost")) "TOST p" else if (own %in% c("p", "p_value", "p.value")) "unadjusted p" else own))
  }
  pcol <- intersect(c("p_value", "p.value", "p_adjusted"), names(x))
  if (!length(pcol)) return(NULL)
  .results_sentence_named_p(x[[pcol[1L]]][1L],
    if (identical(pcol[1L], "p_adjusted")) "multiplicity-adjusted p" else "unadjusted p")
}

# Pull a p-value out of a long interval-carrying table, if one is present:
# either a dedicated p-value column, the row the table names in its
# p_sentence attribute (the TOST p-value of the equivalence tables, which
# list the p-values of the two one-sided tests first), or the first
# p_terms-named row of the primary column. Not exported.
.results_sentence_p_long <- function(x, primary) {
  if (.results_sentence_simultaneous(x)) {
    if (!.results_sentence_matching_adjustment(x)) return(NULL)
    if ("p_adjusted" %in% names(x)) return(.results_sentence_named_p(x$p_adjusted[1L], "multiplicity-adjusted p"))
    if ("p_adjusted" %in% x$term) return(.results_sentence_named_p(x[[primary]][x$term == "p_adjusted"][1L], "multiplicity-adjusted p"))
    return(NULL)
  }
  own <- attr(x, "p_sentence", exact = TRUE)
  if (is.character(own) && length(own) == 1L && own %in% x$term) {
    return(.results_sentence_named_p(x[[primary]][x$term == own][1L], if (identical(own, "p_tost")) "TOST p" else if (own %in% c("p", "p_value", "p.value")) "unadjusted p" else own))
  }
  pcol <- intersect(c("p_value", "p.value", "p_adjusted"), names(x))
  if (length(pcol)) return(.results_sentence_p_wide(x))
  p_terms <- attr(x, "p_terms", exact = TRUE)
  if (!is.null(p_terms)) {
    hit <- x$term %in% p_terms
    if (any(hit)) return(.results_sentence_named_p(x[[primary]][hit][1L], "unadjusted p"))
  }
  NULL
}

# Match the p-value to the interval's family, not merely the first column.
.results_sentence_simultaneous <- function(x) {
  is.list(attr(x, "simultaneous", exact = TRUE)) ||
    (.results_sentence_adjustment(x, "ci_adjustment") %||% "none") %in% c("bonferroni", "tukey", "scheffe", "dunnett", "games_howell", "dunnett_t3") ||
    (attr(x, "adjust", exact = TRUE) %||% "none") %in%
      c("bonferroni", "tukey", "scheffe", "dunnett") ||
    inherits(x, "dmar_post_hoc_ci")
}

.results_sentence_adjustment <- function(x, name) {
  own <- x[[name]]
  if (is.character(own)) own <- unique(own)
  if (is.null(own)) own <- attr(x, name, exact = TRUE)
  if (!is.character(own) || length(own) != 1L || is.na(own)) return(NULL)
  tolower(own)
}

.results_sentence_matching_adjustment <- function(x) {
  ci <- .results_sentence_adjustment(x, "ci_adjustment") %||%
    tolower(attr(x, "simultaneous", exact = TRUE)$method %||% attr(x, "adjust", exact = TRUE))
  p <- .results_sentence_adjustment(x, "p_adjustment") %||% attr(x, "adjust", exact = TRUE)
  length(ci) == 1L && length(p) == 1L && !is.na(ci) && !is.na(p) &&
    identical(tolower(ci), tolower(p))
}

# Carry adjustment provenance per row through tidy(), including row stacks.
.dmar_tidy_adjustments <- function(out, x) {
  for (name in c("ci_adjustment", "p_adjustment")) {
    value <- attr(x, name, exact = TRUE)
    if (is.character(value) && length(value) == 1L) out[[name]] <- rep(value, nrow(out))
  }
  out
}

.results_sentence_named_p <- function(p, label) {
  attr(p, "inference_label") <- label
  p
}
