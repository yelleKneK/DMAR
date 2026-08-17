# Effects-coding (deviation-coding) and Helmert-coding contrast matrices.
#' Effects-Coding Contrast Matrix for a Factor
#'
#' Builds the effects-coding (also called \emph{deviation coding} or
#' \emph{sum-to-zero}) contrast matrix for a factor with \eqn{a}
#' levels. Each non-reference level contrasts with the grand mean
#' (rather than with a reference category as in dummy coding). The
#' returned matrix has rows = levels and columns named after the
#' levels, replacing the numeric column names produced by
#' \code{stats::contr.sum()}.
#'
#' @param levels Either an integer giving the number of levels or a
#'   character / factor vector giving the level labels. If integer,
#'   the labels default to \code{"L1"}, \code{"L2"}, \ldots
#' @param reference Optional character name of the reference level
#'   (whose coefficients are all \eqn{-1}). Default is the last
#'   level.
#'
#' @return A numeric \eqn{a \times (a - 1)} matrix with row names =
#'   the factor levels and column names = the non-reference levels.
#'   Suitable for assignment to \code{contrasts(factor)} or use in
#'   manual contrast construction.
#'
#' @details
#' \strong{Why effects coding.} Effects coding gives the regression
#' intercept the interpretation of the grand mean (rather than the
#' reference-category mean), and each slope coefficient becomes the
#' deviation of that level's mean from the grand mean (rather than
#' the difference vs the reference category). For balanced designs
#' the effect coefficients are orthogonal to the intercept.
#'
#' \strong{Equivalent to.} \code{stats::contr.sum()} but with
#' meaningful column names (the level labels), which is what is lost
#' in the base-R implementation.
#'
#' @references
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003).
#'   \emph{Applied multiple regression/correlation analysis for the
#'   behavioral sciences} (3rd ed.). Lawrence Erlbaum. (See Chapter
#'   8.)
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapters 4, 7.)
#'
#' @seealso \code{\link[stats]{contr.sum}}, \code{\link{helmert_coding}},
#'   \code{\link{is_orthogonal_set}}
#'
#' @examples
#' # 1. Effects coding for a 4-level factor:
#' effects_coding(c("low", "med", "high", "very_high"))
#'
#' # 2. With "low" as the reference category:
#' effects_coding(c("low", "med", "high", "very_high"),
#'                reference = "low")
#'
#' # 3. Assign to a factor for modeling:
#' f <- factor(c("a", "b", "c", "a", "b", "c"))
#' contrasts(f) <- effects_coding(levels(f))
#' model.matrix(~ f)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family design utilities
#'
#' @export

effects_coding <- function(levels, reference = NULL) {
  labs <- .normalize_levels(levels)
  a <- length(labs)
  if (a < 2L) stop("Need at least 2 levels.")
  if (is.null(reference)) reference <- labs[a]
  if (!reference %in% labs)
    stop(sprintf("'reference' = '%s' is not among the levels.", reference))

  other <- setdiff(labs, reference)
  M <- matrix(0, nrow = a, ncol = a - 1L,
              dimnames = list(labs, other))
  for (j in seq_along(other)) M[other[j], j] <- 1
  M[reference, ] <- -1
  M
}

#' Helmert-Coding Contrast Matrix for a Factor
#'
#' Builds the Helmert-coding contrast matrix for a factor with \eqn{a}
#' levels. The \eqn{k}-th column contrasts the \eqn{(k + 1)}-th level
#' against the average of all preceding levels, giving a fully
#' orthogonal set under equal sample sizes. The returned matrix has columns named after the
#' contrasted level rather than the numeric column names produced by
#' \code{stats::contr.helmert()}.
#'
#' @param levels Either an integer giving the number of levels or a
#'   character / factor vector giving the level labels. If integer,
#'   the labels default to \code{"L1"}, \code{"L2"}, \ldots
#'
#' @return A numeric \eqn{a \times (a - 1)} matrix with row names =
#'   the factor levels and column names of the form
#'   \code{"L2_vs_prior"}, \code{"L3_vs_prior"}, \ldots
#'
#' @details
#' \strong{Why Helmert.} Helmert contrasts are the canonical
#' "sequential" orthogonal contrast set: under equal-\eqn{n}, every
#' column is orthogonal to every other column and to the intercept.
#' They are useful when the factor has a natural ordering and the
#' research questions are "does the \eqn{k}-th level differ from the
#' average of the preceding levels?"
#'
#' \strong{Equivalent to.} \code{stats::contr.helmert()} but with
#' interpretable column names.
#'
#' @references
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003).
#'   \emph{Applied multiple regression/correlation analysis for the
#'   behavioral sciences} (3rd ed.). Lawrence Erlbaum.
#'
#' @seealso \code{\link[stats]{contr.helmert}},
#'   \code{\link{effects_coding}}, \code{\link{is_orthogonal_set}}
#'
#' @examples
#' # 1. Helmert coding for a 4-level factor:
#' helmert_coding(c("baseline", "week1", "week2", "week3"))
#'
#' # 2. Confirm orthogonality:
#' M <- helmert_coding(4)
#' is_orthogonal_set(M)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family design utilities
#'
#' @export

helmert_coding <- function(levels) {
  labs <- .normalize_levels(levels)
  a <- length(labs)
  if (a < 2L) stop("Need at least 2 levels.")
  M <- stats::contr.helmert(a)
  rownames(M) <- labs
  colnames(M) <- paste0(labs[-1L], "_vs_prior")
  M
}

# ---- internal ----

.normalize_levels <- function(levels) {
  if (is.factor(levels)) return(levels(levels))
  if (is.character(levels)) return(levels)
  if (is.numeric(levels) && length(levels) == 1L && levels >= 2L)
    return(paste0("L", seq_len(levels)))
  stop("'levels' must be a character vector, factor, or single positive integer.")
}
