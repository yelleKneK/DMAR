#' Asymptotic Variance of the Intraclass Correlation Coefficient
#'
#' Computes the asymptotic (large-sample) variance of the intraclass
#' correlation coefficient (ICC) for any of the six classical Shrout-Fleiss
#' (1979) forms, given a value of the population ICC, the number of subjects
#' \eqn{n}, and the number of raters \eqn{k}. The single-rater forms use
#' Smith's (1956) one-way-ANOVA asymptotic variance, equivalently the
#' Fisher (1925) information-matrix result, and the average-of-\eqn{k}
#' forms use a delta method transformation through the Spearman-Brown
#' relation.
#'
#' @param rho The intraclass correlation coefficient at which the
#'   asymptotic variance is evaluated, on the scale matching \code{type}:
#'   for the single-rater types (\code{"ICC(1,1)"}, \code{"ICC(2,1)"},
#'   \code{"ICC(3,1)"}) supply the single-rater ICC; for the
#'   average-of-\eqn{k} types (\code{"ICC(1,k)"}, \code{"ICC(2,k)"},
#'   \code{"ICC(3,k)"}) supply the average-of-\eqn{k} ICC. Must lie in
#'   \eqn{[0, 1]}.
#'
#'   \strong{Population value or sample value?} The formula is derived
#'   in terms of the unknown population value \eqn{\rho}. In applied
#'   work \eqn{\rho} is never known, so two conventions are common:
#'   (i) for \emph{prospective} variance calculation (e.g., when
#'   planning a study and asking how precise an estimator will be at
#'   plausible truths), supply your \emph{anticipated} population value
#'   motivated by prior literature, theory, or a pilot; (ii) for
#'   \emph{post hoc} variance calculation (e.g., constructing a Wald
#'   standard error or a meta-analytic weight from an observed sample),
#'   supply the \emph{sample} estimate \eqn{\hat\rho} as a plug-in for
#'   the population value. The two postures look identical at the call
#'   site but have different conceptual content; the plug-in case (ii)
#'   yields a consistent (not exact) estimator of the variance.
#' @param n Number of subjects (targets) rated.
#' @param k Number of raters (or repeated measurements) per subject.
#' @param type Which Shrout-Fleiss ICC variant the population value
#'   represents. One of \code{"ICC(1,1)"}, \code{"ICC(2,1)"},
#'   \code{"ICC(3,1)"}, \code{"ICC(1,k)"}, \code{"ICC(2,k)"},
#'   \code{"ICC(3,k)"}. The shorthand aliases used in \code{\link{icc}}
#'   (\code{"1"}, \code{"2"}, \code{"3"}, \code{"1k"}, \code{"2k"},
#'   \code{"3k"}) are also accepted.
#'
#' @return A one-row \code{data.frame} with columns \code{term} (always
#'   \code{"var_icc"}) and \code{value} (the asymptotic variance).
#'
#' @details
#' \strong{Single-rater forms.} For the single-rater intraclass correlation
#' in a balanced design with \eqn{n} subjects and \eqn{k} measurements per
#' subject, the standard large-sample variance derived from the Fisher
#' information matrix is (Smith, 1956; Donner, 1986; Searle, 1971, ch. 11):
#' \deqn{\mathrm{Var}(\hat\rho) \;\approx\;
#'         \frac{2 (1 - \rho)^2 \bigl(1 + (k-1)\rho\bigr)^2}
#'              {n\,k\,(k - 1)}.}
#' This expression is exact under the one-way random-effects model
#' (\code{ICC(1,1)}) and serves as the standard large-sample approximation
#' for the two-way single-rater forms (\code{ICC(2,1)}, \code{ICC(3,1)})
#' as well; differences with the exact two-way variance vanish at order
#' \eqn{1/n} (Bonett, 2002; Burdick & Graybill, 1992). For small \eqn{n} or
#' moderate \eqn{k} where exact two-way inference matters, prefer the
#' \emph{F}-distribution-based confidence intervals returned by
#' \code{\link{icc}}, which follow Shrout and Fleiss (1979) directly.
#'
#' \strong{Average-of-\emph{k} forms.} Applying the Spearman-Brown
#' transformation \eqn{\rho_k = k\rho / [1 + (k - 1)\rho]} together with the
#' delta method gives the closed-form
#' \deqn{\mathrm{Var}(\hat\rho_k) \;\approx\;
#'         \frac{2\,k\,(1 - \rho_k)^2}{n\,(k - 1)},}
#' expressed directly in the average-level ICC \eqn{\rho_k} (so the user
#' need not invert Spearman-Brown when working with reliability of
#' composites). The reduction to this form follows from the substitutions
#' \eqn{1 - \rho = k(1 - \rho_k)/[k - (k-1)\rho_k]} and
#' \eqn{1 + (k-1)\rho = k / [k - (k-1)\rho_k]}.
#'
#' \strong{Use cases.} The asymptotic variance is the natural ingredient for
#' Wald-style inference, sample size planning for the width of an ICC
#' confidence interval (compare with Bonett, 2002, which uses a Fisher-style
#' transformation), and meta-analytic synthesis of ICCs across studies (the
#' inverse of \code{value} weights each study). For confidence intervals
#' themselves, prefer \code{\link{icc}}, which uses the exact
#' \emph{F}-distribution inversion of Shrout and Fleiss (1979, pp. 425--426).
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for estimating intraclass
#'   correlations with desired precision. \emph{Statistics in Medicine,
#'   21}(9), 1331--1335. \doi{10.1002/sim.1108}
#'
#' Burdick, R. K., & Graybill, F. A. (1992). \emph{Confidence Intervals on
#'   Variance Components}. Marcel Dekker.
#'
#' Donner, A. (1986). A review of inference procedures for the intraclass
#'   correlation coefficient in the one-way random effects model.
#'   \emph{International Statistical Review, 54}(1), 67--82.
#'
#' Fisher, R. A. (1925). \emph{Statistical Methods for Research Workers}.
#'   Oliver & Boyd.
#'
#' McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some
#'   intraclass correlation coefficients. \emph{Psychological Methods,
#'   1}(1), 30--46. \doi{10.1037/1082-989X.1.1.30}
#'
#' Searle, S. R. (1971). \emph{Linear Models}. Wiley.
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
#'   assessing rater reliability. \emph{Psychological Bulletin, 86}(2),
#'   420--428.
#'
#' Smith, C. A. B. (1956). On the estimation of intraclass correlation.
#'   \emph{Annals of Human Genetics, 21}(4), 363--373.
#'
#' @examples
#' # Single-rater one-way ICC at rho = .60 with 30 subjects and 4 raters.
#' var_icc(rho = 0.60, n = 30, k = 4, type = "ICC(1,1)")
#'
#' # Same study, but expressed at the average-of-4-rater level.
#' #     The Spearman-Brown-transformed population value is
#' #     rho_k = 4*.6 / (1 + 3*.6) = 0.857
#' var_icc(rho = 0.857, n = 30, k = 4, type = "ICC(1,k)")
#'
#' # Two-way mixed-model consistency ICC (uses the same one-way
#' #     asymptotic variance as a large-sample approximation; see Details).
#' var_icc(rho = 0.60, n = 30, k = 4, type = "ICC(3,1)")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{icc}}, \code{\link{ss_aipe_reliability}},
#'   \code{\link{var_R2}}
#'
#' @keywords design univar
#'
#' @export
var_icc <- function(rho, n, k, type = "ICC(1,1)") {
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) ||
      rho < 0 || rho > 1) {
    stop("'rho' must be a single number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n < 2L) {
    stop("'n' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k < 2L) {
    stop("'k' must be a single integer >= 2.", call. = FALSE)
  }
  if (length(type) != 1L) {
    stop("'type' must be a single ICC label; see ?icc for valid values.",
         call. = FALSE)
  }

  alias <- c(
    "1"        = "ICC(1,1)", "2"        = "ICC(2,1)", "3"       = "ICC(3,1)",
    "1k"       = "ICC(1,k)", "2k"       = "ICC(2,k)", "3k"      = "ICC(3,k)",
    "1,1"      = "ICC(1,1)", "2,1"      = "ICC(2,1)", "3,1"     = "ICC(3,1)",
    "1,k"      = "ICC(1,k)", "2,k"      = "ICC(2,k)", "3,k"     = "ICC(3,k)",
    "ICC(1,1)" = "ICC(1,1)", "ICC(1,k)" = "ICC(1,k)",
    "ICC(2,1)" = "ICC(2,1)", "ICC(2,k)" = "ICC(2,k)",
    "ICC(3,1)" = "ICC(3,1)", "ICC(3,k)" = "ICC(3,k)"
  )
  type_full <- unname(alias[as.character(type)])
  if (is.na(type_full)) {
    stop("Unrecognized ICC type: ", type,
         ". Valid: 1, 2, 3, 1k, 2k, 3k, ICC(1,1), ..., ICC(3,k).",
         call. = FALSE)
  }

  single_rater <- type_full %in% c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)")
  if (single_rater) {
    value <- 2 * (1 - rho)^2 * (1 + (k - 1) * rho)^2 / (n * k * (k - 1))
  } else {
    value <- 2 * k * (1 - rho)^2 / (n * (k - 1))
  }

  out <- data.frame(term = "var_icc", value = value,
                    stringsAsFactors = FALSE, row.names = NULL)
  .as_dmar_tbl(out)
}
