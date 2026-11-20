#' Missing Data Patterns
#'
#' Tabulates the distinct patterns of observed and missing values across a
#' set of items: which items each pattern leaves missing, and how many
#' cases show it. The table is the first description of incomplete data,
#' the one a reader needs before any model is fit. It shows whether the
#' missing values are scattered across many patterns or concentrated in a
#' few, whether they arrive together (two tests given to only one school
#' share one pattern), and how many cases are complete, which is all
#' listwise deletion would keep. Full information maximum likelihood
#' (FIML) uses every pattern, and the fraction of missing information
#' (\code{\link{fraction_missing_information}}) then says what the
#' patterns cost each parameter.
#'
#' @param data A \code{data.frame}.
#' @param items Optional character vector naming the columns of
#'   \code{data} to describe; \code{NULL} (the default) uses every column.
#'   Any column type is allowed, since only whether a value is \code{NA}
#'   matters.
#'
#' @details
#' A pattern is written as a string of 1s and 0s, one character per item
#' in the order of \code{items}, with 1 for an observed value and 0 for a
#' missing one (the convention of \code{mice::md.pattern()}). The rows are
#' sorted by the number of cases, largest first, then by the number of
#' items missing, fewest first, so the complete cases usually lead and
#' the rare patterns trail.
#'
#' The pattern in which every item is missing is kept as a row
#' (\code{n_observed = 0}) when it occurs. \pkg{lavaan} sets such cases
#' aside before fitting, without reporting them among its own patterns,
#' so the number of cases a FIML fit uses is \code{N} minus
#' \code{N_all_missing}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   distinct pattern and the same columns whatever the number of items:
#'   \code{pattern} (the string of 1s and 0s), \code{missing_items} (the
#'   names of the items the pattern leaves missing, separated by commas;
#'   \code{""} for the complete pattern), \code{n} (the number of cases),
#'   \code{proportion} (of all rows of \code{data}), \code{n_observed},
#'   and \code{n_missing} (the numbers of items observed and missing in
#'   the pattern). The attributes are \code{matrix} (the logical pattern
#'   matrix, \code{TRUE} for observed, one row per pattern in the order of
#'   the table, the layout of \code{lavaan::lavInspect(fit, "patterns")}),
#'   \code{N} (the rows of \code{data}), \code{N_complete},
#'   \code{N_all_missing}, \code{n_patterns}, and \code{items}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{missing_coverage}} for the proportion of cases
#'   observed on each pair of items; \code{\link{fraction_missing_information}}
#'   for what the missing values cost each parameter of a FIML fit.
#'
#' @family missing data
#'
#' @keywords manip
#'
#' @examples
#' # The six spatial tests of the Holzinger and Swineford data: two
#' # patterns, the 145 Grant-White pupils with every test and the 156
#' # Pasteur pupils, who were not given the revised paper form board
#' # (t25) or the flags test (t26).
#' spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
#'              "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
#' missing_patterns(holzinger_swineford, items = spatial)
#'
#' # Values missing completely at random on two more tests add patterns.
#' set.seed(113)
#' d <- holzinger_swineford[spatial]
#' d$t1_visual_perception[runif(nrow(d)) < .10] <- NA
#' d$t2_cubes[runif(nrow(d)) < .10] <- NA
#' pat <- missing_patterns(d)
#' pat
#' attr(pat, "N_complete")
#'
#' @export
missing_patterns <- function(data, items = NULL) {
  observed <- .missing_item_frame(data, items)
  items <- colnames(observed)
  N <- nrow(observed)
  key <- apply(observed, 1L, function(r) paste(as.integer(r), collapse = ""))
  first <- !duplicated(key)
  pattern <- key[first]
  pattern_matrix <- observed[first, , drop = FALSE]
  n <- as.integer(tabulate(match(key, pattern), nbins = length(pattern)))
  n_observed <- as.integer(rowSums(pattern_matrix))
  n_missing <- length(items) - n_observed
  missing_items <- apply(pattern_matrix, 1L, function(r)
    paste(items[!r], collapse = ", "))

  # Largest first, then fewest items missing, then the pattern string in
  # decreasing order (more items observed earlier in the item order), so
  # the order never depends on the order of the rows of `data`.
  ord <- order(-n, n_missing, -xtfrm(pattern))
  out <- data.frame(
    pattern       = pattern[ord],
    missing_items = unname(missing_items[ord]),
    n             = n[ord],
    proportion    = n[ord] / N,
    n_observed    = n_observed[ord],
    n_missing     = as.integer(n_missing[ord]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  pattern_matrix <- pattern_matrix[ord, , drop = FALSE]
  rownames(pattern_matrix) <- NULL

  out <- .as_dmar_tbl(out)
  attr(out, "matrix") <- pattern_matrix
  attr(out, "N") <- N
  attr(out, "N_complete") <- sum(rowSums(observed) == length(items))
  attr(out, "N_all_missing") <- sum(rowSums(observed) == 0L)
  attr(out, "n_patterns") <- length(pattern)
  attr(out, "items") <- items
  out
}
