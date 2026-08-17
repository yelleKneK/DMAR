#' Cohen's (1968) Psychiatric Diagnosis Agreement Table
#'
#' The illustrative agreement matrix from Cohen's (1968) weighted kappa
#' paper, Table 1: two judges independently assign \emph{N} = 200 cases
#' to three diagnostic categories (personality disorder, neurosis,
#' psychosis). The data set reproduces the printed table cell for cell
#' and in its original layout, Judge B indexing the rows and Judge A the
#' columns, one row per cell of the 3 x 3 matrix. Each cell carries the
#' three quantities Cohen prints: the ratio-scaled disagreement weight,
#' the chance-expected proportion (his parenthetical values), and the
#' observed proportion; the raw frequency is the observed proportion
#' times \emph{N}.
#'
#' @format A data frame with 9 observations (one per cell of the 3 x 3
#' agreement matrix) on 6 variables.
#' \describe{
#'   \item{\code{judge_b}}{Factor: Judge B's diagnostic category (the
#'     table's rows), with levels \code{Personality disorder},
#'     \code{Neurosis}, \code{Psychosis}.}
#'   \item{\code{judge_a}}{Factor: Judge A's diagnostic category (the
#'     table's columns), same levels.}
#'   \item{\code{frequency}}{Number of the 200 cases jointly assigned
#'     to the cell.}
#'   \item{\code{disagreement_weight}}{Cohen's ratio-scaled
#'     disagreement weight \eqn{v_{ij}} for the cell: 0 on the
#'     agreement diagonal, 1 for a personality disorder-neurosis
#'     confusion, 3 for personality disorder-psychosis, and 6 for
#'     neurosis-psychosis, the confusion the illustration treats as
#'     gravest.}
#'   \item{\code{observed_proportion}}{Observed proportion of cases in
#'     the cell, \code{frequency / 200}.}
#'   \item{\code{expected_proportion}}{Chance-expected proportion of
#'     cases in the cell, the product of the cell's row (Judge B) and
#'     column (Judge A) marginal proportions; the parenthetical values
#'     in Cohen's Table 1.}
#' }
#'
#' @details
#' Cohen built this table to make a point that is easy to miss:
#' weighted kappa is fully chance corrected, and it can be
#' \emph{smaller} than unweighted kappa on the same data. Here the
#' judges disagree far less than chance expectation in the mildly
#' weighted personality disorder-neurosis cells but at about the chance
#' level in the heavily weighted neurosis-psychosis cells, so
#' \eqn{\kappa = .492} while \eqn{\kappa_W = .348}: they disagree least
#' where it matters least. Interchanging the 6 and 1 weights reverses
#' the conclusion (\eqn{\kappa_W = .574}).
#'
#' The reconstruction was verified against every quantity Cohen
#' computes from the table: the marginals (.50/.30/.20 for Judge B,
#' .60/.30/.10 for Judge A), the chance-expected cell proportions, the
#' weighted disagreement sums \eqn{q'_o = .90} and \eqn{q'_c = 1.38},
#' \eqn{\kappa = .492}, \eqn{\kappa_W = .348}, and his Formula 10 and
#' 13 standard errors (.0901 and .0916). The \code{\link{cohen_kappa}}
#' help page replicates the full set of analyses, and the weighted
#' kappa vignette works the illustration end to end, including the
#' orientation of the printed weight display in the paper's
#' asymmetric-weight validity reinterpretation.
#'
#' @source
#' Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision
#' for scaled disagreement or partial credit. \emph{Psychological
#' Bulletin, 70}(4), 213--220 (Table 1, p. 214).
#'
#' @references
#' Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision
#'   for scaled disagreement or partial credit. \emph{Psychological
#'   Bulletin, 70}(4), 213--220. \doi{10.1037/h0026256}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cohen_kappa}}, which analyzes this table in its
#' examples.
#'
#' @examples
#' data(diagnosis_agreement)
#'
#' # Rebuild Cohen's Table 1 layout (Judge B in rows, Judge A in columns).
#' xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
#'
#' # Unweighted and weighted kappa, reproducing kappa = .492 and
#' # kappa_w = .348.
#' tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
#' v   <- xtabs(disagreement_weight ~ judge_b + judge_a,
#'              data = diagnosis_agreement)
#' cohen_kappa(table = tab)
#' cohen_kappa(table = tab, weights = unclass(v),
#'             weight_scaling = "disagreement")
#'
#' # Cohen's Formula 8 directly from the per-cell quantities.
#' with(diagnosis_agreement,
#'      1 - sum(disagreement_weight * observed_proportion) /
#'          sum(disagreement_weight * expected_proportion))
#'
#' @keywords datasets
#' @family reliability
"diagnosis_agreement"
