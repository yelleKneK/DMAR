# Asymptotic variance of the robust trimmed SMD for AIPE planning.
#' Asymptotic Variance of the Robust Trimmed SMD
#'
#' Computes the asymptotic variance of the Algina-Keselman-Penfield
#' (2005) robust standardized mean difference under the Yuen (1974)
#' trimmed-mean framework, suitable for AIPE sample size planning
#' for robust effect sizes (Keselman, Algina, Lix, Wilcox, & Deering,
#' 2008).
#'
#' @param population_smd_trimmed Anticipated population value of the
#'   robust trimmed SMD \eqn{\delta_R}.
#' @param n_1,n_2 Per-group sample sizes.
#' @param trim Proportion to trim and Winsorize. Default \code{0.20}.
#'
#' @return A 1-row \code{data.frame} with columns \code{term}
#'   (\code{"var_smd_trimmed"}) and \code{value} (the variance).
#'
#' @details
#' \strong{Variance formula.} Under random sampling with trimming
#' proportion \eqn{\gamma} from each tail, the variance of the
#' trimmed-mean difference scales by \eqn{1 / h_j} (where \eqn{h_j =
#' n_j - 2 \lfloor \gamma n_j \rfloor} is the number of retained
#' observations in group \eqn{j}) rather than \eqn{1 / n_j}. The
#' large-sample variance of the standardized version, written on the
#' \eqn{d_R} scale, is
#' \deqn{\mathrm{Var}(\hat d_R) \;\approx\;
#'   \frac{h_1 + h_2}{h_1 h_2} +
#'   \frac{\delta_R^2}{2 (h_1 + h_2)}.}
#' For \eqn{\gamma = 0} this reduces to the standard Hedges-Olkin
#' (1985) variance of Cohen's \emph{d}.
#'
#' \strong{When to use.} For AIPE planning of a robust effect size
#' study, use \code{var_smd_trimmed()} in place of
#' \code{\link{var_smd}()}. Pair with \code{\link{smd_trimmed}()} for
#' the point estimate and noncentral \emph{t} CI.
#'
#' @references
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An
#'   alternative to Cohen's standardized mean difference effect
#'   size: A robust parameter and confidence interval in the two
#'   independent groups case. \emph{Psychological Methods, 10}(3),
#'   317--328. \doi{10.1037/1082-989X.10.3.317}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Keselman, H. J., Algina, J., Lix, L. M., Wilcox, R. R., & Deering,
#'   K. N. (2008). A generally robust approach for testing hypotheses
#'   and setting confidence intervals for effect sizes.
#'   \emph{Psychological Methods, 13}(2), 110--129.
#'   \doi{10.1037/1082-989X.13.2.110}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Yuen, K. K. (1974). The two-sample trimmed \emph{t} for unequal
#'   population variances. \emph{Biometrika, 61}(1), 165--170.
#'
#' @seealso \code{\link{smd_trimmed}}, \code{\link{var_smd}},
#'   \code{\link{ss_aipe_smd}}
#'
#' @examples
#' # 1. Population delta_R = 0.5, n = 30 per group, 20% trim:
#' var_smd_trimmed(population_smd_trimmed = 0.5, n_1 = 30, n_2 = 30)
#'
#' # 2. Variance scales by h / n via the trimming proportion:
#' var_smd_trimmed(0.5, 30, 30, trim = 0.00)$value
#' var_smd_trimmed(0.5, 30, 30, trim = 0.20)$value
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_smd_trimmed <- function(population_smd_trimmed, n_1, n_2,
                            trim = 0.20) {
  for (nm in c("population_smd_trimmed", "n_1", "n_2", "trim")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L)
      stop(sprintf("'%s' must be a single numeric value.", nm))
  }
  if (n_1 < 4 || n_2 < 4) stop("'n_1' and 'n_2' must each be >= 4.")
  if (trim < 0 || trim >= 0.5) stop("'trim' must be in [0, 0.5).")

  h_1 <- n_1 - 2L * floor(trim * n_1)
  h_2 <- n_2 - 2L * floor(trim * n_2)
  v <- (h_1 + h_2) / (h_1 * h_2) +
       population_smd_trimmed^2 / (2 * (h_1 + h_2))

  out <- data.frame(term  = "var_smd_trimmed",
                    value = v,
                    stringsAsFactors = FALSE,
                    row.names = NULL)
  .as_dmar_tbl(out)
}
