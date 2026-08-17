# Critical value for the Bryant-Paulson generalized studentized range (ANCOVA)
#' Provides the Critical Value for the Bryant--Paulson ANCOVA Multiple-Comparison Procedure
#'
#' Computes the critical value of the Bryant--Paulson generalized studentized
#' range, the reference distribution for multiple comparisons of adjusted
#' means in an analysis of covariance with random covariates. The single-step
#' value is the multiplier for simultaneous confidence intervals on pairwise
#' differences of adjusted means; the \code{"duncan"} option instead returns
#' the significant range of the stepwise Duncan multiple-range procedure.
#'
#' @param alpha_level Type I error rate (i.e., the false positive rate). As with
#'   \code{\link{cv_tukey_hsd}}, the full \code{alpha_level} applies to the upper tail
#'   of the (non-negative) generalized studentized range distribution; it is
#'   not split between two tails (see Details).
#' @param df The ANCOVA error degrees of freedom (a positive number; in a
#'   one-way ANCOVA, \code{df = N - groups - covariates}).
#' @param groups The number of groups whose adjusted means are being compared
#'   (an integer of at least 2).
#' @param covariates The number of \emph{random} covariates in the ANCOVA, the
#'   parameter \eqn{p} of the Bryant--Paulson distribution (a non-negative
#'   integer). With \code{covariates = 0} the critical value reduces to the
#'   ordinary studentized range and the result equals \eqn{\sqrt2} times
#'   \code{\link{cv_tukey_hsd}} (see Details).
#' @param procedure One of \code{"tukey"} (default) or \code{"duncan"}.
#'   \code{"tukey"} returns the single-step \emph{simultaneous} critical value
#'   \eqn{q_{\alpha;p,k,\nu}} used for familywise (Tukey--Kramer-type)
#'   confidence intervals; \code{"duncan"} returns the stepwise Duncan
#'   multiple-range significant range \eqn{r_{\alpha;p,k,\nu}} tabled by Bryant
#'   and Bruvold (1980, Table 2) (see Details).
#' @param verbose Provides extra information (the tail areas) about the
#'   critical value.
#'
#' @return Returns the critical value in a output style (a
#'   \code{data.frame} with class \code{dmar_tbl} and one row per critical
#'   value, following the format used by \code{\link{cv_tukey_hsd}} and
#'   \code{\link{cv_t}}). The \code{value} is on the \emph{studentized-range
#'   scale} (the scale on which Bryant and Paulson tabulate their critical
#'   values and on which \code{\link{ci_c_ancova_bp}} uses them). When
#'   \code{verbose = TRUE} and \code{procedure = "tukey"}, the upper- and
#'   lower-tail areas of the Bryant--Paulson distribution at the critical value
#'   are also returned; for \code{procedure = "duncan"} the tail areas are
#'   \code{NA} because the significant range is a stepwise quantity rather than
#'   a single quantile.
#'
#' @details
#' The Bryant--Paulson procedure is the analysis-of-covariance generalization
#' of Tukey's method (\code{\link{cv_tukey_hsd}}) for comparing adjusted means
#' when the covariate is \emph{random}. Because the covariate adjustment must
#' be estimated, the studentized range of adjusted means is stochastically
#' larger than the ordinary studentized range, so the Bryant--Paulson critical
#' value exceeds Tukey's; using the latter would give intervals that are too
#' narrow and a familywise error rate above \code{alpha_level}. The single-step
#' (\code{procedure = "tukey"}) value is the multiplier for a family of
#' simultaneous confidence intervals on the pairwise differences of adjusted
#' means that jointly hold at level \eqn{1 - \alpha}. Maxwell, Delaney, and
#' Kelley (2027, Chapter 9) develop multiple comparisons of adjusted means in
#' the analysis of covariance, the setting this critical value serves.
#'
#' The reference distribution is the Bryant--Paulson generalized studentized
#' range, implemented in \code{\link{qbryant_paulson}}. Its quantiles are not a
#' standard base-R distribution and are not the multivariate \emph{t} quantiles
#' that \code{\link{cv_dunnett}} and \code{\link{cv_smm}} obtain from
#' \pkg{mvtnorm}; they are computed directly by \code{\link{qbryant_paulson}},
#' so this function depends on neither base-R nor \pkg{mvtnorm} multiple-mean
#' machinery.
#'
#' \strong{Scale.} The returned \code{value} is on the studentized-range scale,
#' \eqn{q_{\alpha;p,k,\nu}}, the scale of Bryant and Paulson's (1976) and
#' Bryant and Bruvold's (1980) tables and of Eq. (2.4) of the latter. A pair of
#' adjusted means is declared different when
#' \eqn{|\hat\theta_i - \hat\theta_j| > q_{\alpha;p,k,\nu}\,\hat\sigma_{y\mid x}\sqrt{1/n}}.
#' This differs from \code{\link{cv_tukey_hsd}}, which divides its value by
#' \eqn{\sqrt2} to report on the pairwise mean-difference scale; divide the
#' value here by \eqn{\sqrt2} to obtain that scale. With \code{covariates = 0},
#' \code{cv_bryant_paulson} returns exactly \eqn{\sqrt2 \times}
#' \code{cv_tukey_hsd}.
#'
#' \strong{Duncan multiple-range.} For \code{procedure = "duncan"} the value is
#' the \dQuote{significant range} of Duncan's stepwise test as extended to
#' ANCOVA by Bryant and Bruvold (1980, Section 4). With variable protection
#' levels \eqn{\alpha_k = 1 - (1-\alpha)^{k-1}},
#' \deqn{r_{\alpha;p,2,\nu} = q_{\alpha;p,2,\nu}, \qquad
#'       r_{\alpha;p,k,\nu} = \max\{\, r_{\alpha;p,k-1,\nu},\;
#'       q_{\alpha_k;p,k,\nu} \,\}, \quad k > 2.}
#' These are the values in Bryant and Bruvold's Table 2, reproduced by this
#' function to the tabled two-decimal precision (see the package tests).
#'
#' @references
#' Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method of
#'   multiple comparisons to experimental designs with random concomitant
#'   variables. \emph{Biometrika, 63}, 631--638.
#'
#' Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures in
#'   the analysis of covariance. \emph{Journal of the American Statistical
#'   Association, 75}(372), 874--880. \doi{10.2307/2287175}
#'
#' Duncan, D. B. (1955). Multiple range and multiple F tests.
#'   \emph{Biometrics, 11}, 1--42.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9, where multiple comparisons of
#'   adjusted means in the analysis of covariance are developed; Appendix
#'   Table A.8 reports these critical values.)
#'
#' @examples
#' # Multiple comparisons of adjusted means in ANCOVA (the setting of Maxwell,
#' # Delaney, and Kelley, 2027, Chapter 9), using the worked example of Bryant
#' # and Bruvold (1980): 6 panels, 1 random covariate, 14 error df,
#' # alpha_level = .05. The single-step value is the multiplier for simultaneous
#' # confidence intervals on the pairwise differences of adjusted means. With a
#' # covariate present there is no closed form for it: the Bryant-Paulson
#' # distribution function is integrated numerically and then inverted, which
#' # takes about half a second, so the call is shown here rather than run.
#' # cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 1)
#' # It returns 4.83, the entry in Table 1 of Bryant and Paulson (1976) and the
#' # multiplier behind the simultaneous intervals of the 1980 worked example.
#'
#' # With no covariates the reference distribution is the ordinary studentized
#' # range, which base R supplies directly, so these calls are quick and do run
#' # here. The critical value is sqrt(2) times the Tukey HSD critical value.
#' cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 0)$value
#' sqrt(2) * cv_tukey_hsd(alpha_level = .05, df = 14, groups = 6)$value
#'
#' # The stepwise Duncan multiple-range significant range comes from
#' # procedure = "duncan". With covariates = 0 it reduces to Duncan's (1955)
#' # own significant studentized range, 3.37 for a stretch of 6 groups on 14
#' # error degrees of freedom.
#' cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 0,
#'                   procedure = "duncan")
#' # One random covariate raises that range to 3.50, the entry in Table 2 of
#' # Bryant and Bruvold (1980). Being stepwise, it inverts a separate
#' # Bryant-Paulson quantile for every stretch from 2 to 6 groups, so it costs
#' # several seconds and is not run here either. The package tests and the
#' # tests check it against the paper's Section 4 example.
#' # cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 1,
#' #                   procedure = "duncan")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_tukey_hsd}}, \code{\link{cv_scheffe}},
#'   \code{\link{qbryant_paulson}}, \code{\link{ci_c_ancova_bp}}
#'
#' @keywords htest design
#'
#' @family critical values
#'
#' @export

cv_bryant_paulson <- function(alpha_level, df, groups, covariates = 1,
                              procedure = c("tukey", "duncan"),
                              verbose = TRUE) {
  procedure <- match.arg(procedure)

  # Argument checks (mirroring cv_tukey_hsd).
  if (missing(df)) stop("You must specify the degrees of freedom (i.e., 'df'), which must be a positive number.")
  if (!df > 0) stop("You must specify a positive value for the degrees of freedom.")
  if (missing(groups)) stop("You must specify the number of groups (i.e., 'groups').")
  if (groups < 2) stop("This function is appropriate for situations in which there are two or more groups.")
  if (covariates < 0 || covariates != round(covariates)) stop("'covariates' must be a non-negative integer (the number of covariates).")
  if (missing(alpha_level)) stop("You must specify 'alpha_level'.")
  if (alpha_level <= 0 || alpha_level >= 1) stop("Specify 'alpha_level' to be greater than zero and less than 1.")

  if (procedure == "tukey") {
    # Single-step simultaneous critical value q_{alpha_level; p, k, nu}: the upper
    # alpha_level point of the Bryant-Paulson Q_p distribution (Bryant & Paulson,
    # 1976, Eq. 7 / Table 1).  It is the multiplier in their simultaneous
    # confidence interval, Eq. (8).
    value <- qbryant_paulson(1 - alpha_level, num_covariates = covariates,
                             num_groups = groups, df = df)
    area_less    <- pbryant_paulson(value, num_covariates = covariates,
                                    num_groups = groups, df = df)
    area_greater <- pbryant_paulson(value, num_covariates = covariates,
                                    num_groups = groups, df = df,
                                    lower_tail = FALSE)
  } else {
    # Duncan multiple-range significant range r_{alpha_level; p, k, nu}.
    # Bryant & Paulson (1976, Sec. 5) noted the Duncan extension "would
    # require further computation"; it was carried out by Bryant & Bruvold
    # (1980, Sec. 4, and their Table 2), who define, with variable protection
    # levels alpha_k = 1 - (1-alpha_level)^(k-1):
    #     r_{alpha_level;p,2,nu} = q_{alpha_level;p,2,nu}
    #     r_{alpha_level;p,k,nu} = max( r_{alpha_level;p,k-1,nu}, q_{alpha_k;p,k,nu} ),  k>2,
    # i.e. a running maximum over the single-step q values evaluated at the
    # stretch-dependent protection levels.
    value <- qbryant_paulson(1 - alpha_level, num_covariates = covariates,
                             num_groups = 2, df = df)
    if (groups >= 3) for (j in 3:groups) {
      qj <- qbryant_paulson((1 - alpha_level)^(j - 1), num_covariates = covariates,
                            num_groups = j, df = df)
      value <- max(value, qj)
    }
    area_less <- area_greater <- NA_real_
  }

  term <- "upper_cv"
  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term = term, value = value,
                                   area_less = area_less,
                                   area_greater = area_greater)))
  }
  .as_dmar_tbl(data.frame(term = term, value = value))
}
