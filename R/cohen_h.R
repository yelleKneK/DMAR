#' Cohen's \emph{h} Effect Size for a Difference Between Two Proportions
#'
#' Computes Cohen's \emph{h}, the effect size for the difference between two
#' proportions on the arcsine (variance-stabilizing) scale,
#' \deqn{h = \varphi_1 - \varphi_2, \qquad \varphi_i = 2\,\arcsin\!\sqrt{p_i}.}
#' The arcsine transform spaces proportions so that a given \emph{h} carries the
#' same detectability wherever the proportions sit, which a raw difference
#' \eqn{p_1 - p_2} does not: a shift from .01 to .05 is easier to detect than one
#' from .41 to .45, and \emph{h} reflects that while the raw difference does not.
#' Cohen's \emph{h} is the proportion analogue of the standardized mean
#' difference (\code{\link{smd}}): the effect size on which power and sample size
#' planning for a difference between two proportions is conventionally based.
#'
#' \emph{Cohen's h is a population quantity}: supplied with population
#' proportions it returns the population value, supplied with sample proportions
#' it returns the corresponding sample value. It is signed, positive when
#' \code{p1} exceeds \code{p2}; its magnitude \code{abs()} is the size of the
#' effect irrespective of direction.
#'
#' @param p1,p2 The two proportions, each in \eqn{[0, 1]}. \code{h} is
#'   \eqn{\varphi(p_1) - \varphi(p_2)}, so it is positive when \code{p1} is the
#'   larger.
#'
#' @return A 1-row \code{data.frame} with columns \code{term} and \code{value};
#'   \code{term} is \code{"cohen_h"} and \code{value} is the signed effect size.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum. (See Chapter 6.)
#'
#' @seealso \code{\link{smd}} for the standardized mean difference,
#'   \code{\link{cohen_f}}, \code{\link{ci_proportion}}
#'
#' @examples
#' # A shift from .40 to .55.
#' cohen_h(p1 = 0.55, p2 = 0.40)
#'
#' # The same raw difference near the floor is a larger h, since a difference is
#' # easier to detect where the proportions are small.
#' cohen_h(p1 = 0.20, p2 = 0.05)
#'
#' # Its magnitude is the size irrespective of direction.
#' abs(cohen_h(p1 = 0.40, p2 = 0.55)$value)
#'
#' @export
cohen_h <- function(p1, p2) {
  if (!is.numeric(p1) || length(p1) != 1L || is.na(p1) || p1 < 0 || p1 > 1) {
    stop("'p1' must be a single proportion in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(p2) || length(p2) != 1L || is.na(p2) || p2 < 0 || p2 > 1) {
    stop("'p2' must be a single proportion in [0, 1].", call. = FALSE)
  }
  phi <- function(p) 2 * asin(sqrt(p))
  h <- phi(p1) - phi(p2)
  .as_dmar_tbl(data.frame(term = "cohen_h", value = h))
}
