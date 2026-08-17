#' Unbiased and Adjusted Estimators of the Population Squared Multiple Correlation
#'
#' @description
#' Estimates the population squared multiple correlation coefficient
#' \eqn{\rho^2} from an observed sample \eqn{R^2}, correcting the well-known
#' positive (upward) bias of \eqn{R^2}. Two estimators are available: the
#' (essentially) unbiased Olkin and Pratt (1958) estimator (the default), and
#' the classic Ezekiel (1930) adjusted-\eqn{R^2} shrinkage formula reported by
#' \code{\link[stats]{summary.lm}} as \code{adj.r.squared}. This is the inverse
#' direction of \code{\link{expected_R2}}, which gives the forward expectation
#' \eqn{E[R^2 \mid \rho^2]}.
#'
#' @param R2 Observed sample squared multiple correlation coefficient, in
#'   \eqn{[0, 1]}.
#' @param N Sample size.
#' @param p Number of predictor variables.
#' @param method Which estimator to compute: \code{"olkin_pratt"} (default) for
#'   the (essentially) unbiased Olkin-Pratt (1958) estimator, or
#'   \code{"ezekiel"} for the Ezekiel (1930) adjusted \eqn{R^2} (the value
#'   \code{\link[stats]{summary.lm}} reports as \code{adj.r.squared}).
#'
#' @details
#' The sample \eqn{R^2} overestimates \eqn{\rho^2}; the bias is larger for
#' smaller samples and for more predictors. Two corrections are offered.
#'
#' The Ezekiel (1930) adjusted estimator is
#' \deqn{\hat\rho^2_{\mathrm{Ezekiel}} = 1 - \frac{N - 1}{N - p - 1}\,(1 - R^2).}
#' It reduces the bias but is not unbiased; it is exactly the quantity
#' \code{summary(lm(...))$adj.r.squared} reports.
#'
#' The Olkin and Pratt (1958) estimator is (essentially) unbiased:
#' \deqn{\hat\rho^2_{\mathrm{OP}} = 1 - \frac{N - 3}{N - p - 1}\,(1 - R^2)\;
#'        {}_2F_1\!\left(1, 1; \frac{N - p + 1}{2}; 1 - R^2\right),}
#' where \eqn{{}_2F_1} is the Gaussian hypergeometric function (the same
#' function \code{\link{expected_R2}} uses for the forward direction; see Stuart,
#' Ord, & Arnold, 1999, section 28). Both estimators can fall below 0 for very
#' small \eqn{R^2}; that is expected behavior for a bias-corrected estimator and
#' is not truncated here (matching \code{adj.r.squared}, which is also allowed to
#' be negative).
#'
#' @return A 1-row \code{data.frame} (class \code{dmar_tbl}) with columns
#'   \code{term} and \code{value}. The \code{term} is
#'   \code{"unbiased_population_R2"} when \code{method = "olkin_pratt"} and
#'   \code{"adjusted_population_R2"} when \code{method = "ezekiel"}; \code{value}
#'   is the corresponding estimate of \eqn{\rho^2}.
#'
#' @references
#' Ezekiel, M. (1930). \emph{Methods of correlation analysis}. Wiley.
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple correlation
#'   coefficient: Accuracy in parameter estimation via narrow confidence
#'   intervals. \emph{Multivariate Behavioral Research, 43}(4), 524--555.
#'   \doi{10.1080/00273170802490632}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain
#'   correlation coefficients. \emph{The Annals of Mathematical Statistics,
#'   29}(1), 201--211.
#'
#' Stuart, A., Ord, J. K., & Arnold, S. (1999). \emph{Kendall's advanced
#'   theory of statistics, volume 2A: Classical inference and the linear
#'   model} (6th ed.). Arnold.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{expected_R2}} for the forward expectation
#'   \eqn{E[R^2 \mid \rho^2]}, and \code{\link{ci_R2}}, \code{\link{var_R2}},
#'   \code{\link{ss_aipe_R2}} for interval, variance, and planning tools on the
#'   same effect size.
#'
#' @examples
#' # An observed R^2 = .50 with N = 50 and p = 5 predictors overstates rho^2.
#' # The Olkin-Pratt (essentially unbiased) estimate is the default:
#' unbiased_R2(R2 = .50, N = 50, p = 5)
#'
#' # The Ezekiel adjusted R^2 (what summary(lm) reports) over-shrinks slightly,
#' # so it typically sits a little below the Olkin-Pratt value:
#' unbiased_R2(R2 = .50, N = 50, p = 5, method = "ezekiel")
#'
#' # The Ezekiel option reproduces summary(lm)$adj.r.squared exactly.
#' set.seed(113)
#' d   <- as.data.frame(matrix(rnorm(50 * 6), 50, 6))
#' fit <- lm(V1 ~ ., data = d)
#' s   <- summary(fit)
#' unbiased_R2(R2 = s$r.squared, N = 50, p = 5, method = "ezekiel")$value
#' s$adj.r.squared
#'
#' # The bias (and so the correction) shrinks as N grows for fixed R^2 and p.
#' unbiased_R2(.50, 50, 5)
#' unbiased_R2(.50, 500, 5)
#'
#' @export
unbiased_R2 <- function(R2, N, p, method = c("olkin_pratt", "ezekiel")) {
  method <- match.arg(method)
  if (!is.numeric(R2) || length(R2) != 1L || is.na(R2) || R2 < 0 || R2 > 1)
    stop("'R2' must be a single number in [0, 1].", call. = FALSE)
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 1 || p != round(p))
    stop("'p' must be a single positive integer (the number of predictors).",
         call. = FALSE)
  if (!is.numeric(N) || length(N) != 1L || is.na(N) || N != round(N) || N < p + 2)
    stop("'N' must be a single integer greater than 'p' + 1.", call. = FALSE)

  if (method == "ezekiel") {
    value <- 1 - ((N - 1) / (N - p - 1)) * (1 - R2)
    term  <- "adjusted_population_R2"
  } else {
    # The 2F1 argument is 1 - R2; at R2 == 0 it equals 1, where the series is
    # nudged off the boundary to keep the hypergeometric value finite (Olkin &
    # Pratt, 1958). gsl::hyperg_2F1(1, 1, (N - p + 1) / 2, 1 - R2x) is another
    # way to obtain it; DMAR computes the 2F1 in base R via .hyperg_2F1() (no
    # GSL system dependency; see R/R2_internals.R), which is also more accurate
    # than gsl as the argument approaches 1.
    R2x   <- if (R2 == 0) 9.9e-12 else R2
    value <- 1 - ((N - 3) / (N - p - 1)) * (1 - R2x) *
      .hyperg_2F1(1, 1, (N - p + 1) / 2, 1 - R2x)
    term  <- "unbiased_population_R2"
  }
  .as_dmar_tbl(data.frame(term = term, value = value))
}
