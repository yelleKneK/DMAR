#' Kuder-Richardson Formula 20 (KR-20) With a Confidence Interval
#'
#' @description
#' Estimates Kuder-Richardson formula 20 (Kuder & Richardson, 1937) for a
#' homogeneous composite scored on dichotomous (0/1) items, and returns a
#' confidence interval for the population coefficient.
#'
#' @details
#' Kuder and Richardson's (1937) formula 20 for a \eqn{J}-item composite
#' of binary items is
#' \deqn{KR_{20} = \frac{J}{J-1}\left(1 - \frac{\sum_{j} p_{j} q_{j}}{s_{Y}^{2}}\right),}
#' where \eqn{p_{j}} is the proportion of respondents endorsing item
#' \emph{j} (i.e., scoring 1), \eqn{q_{j} = 1 - p_{j}}, and
#' \eqn{s_{Y}^{2}} is the variance of the composite score. For
#' dichotomous items \eqn{p_{j} q_{j}} is the item variance, so KR-20 is
#' algebraically identical to coefficient \eqn{\alpha} (Guttman, 1945;
#' Cronbach, 1951) computed from the item covariance matrix. KR-20
#' predates coefficient \eqn{\alpha} by 14 years and Feldt's (1965)
#' \emph{F}-distribution interval was derived specifically for the
#' sampling distribution of KR-20.
#'
#' Because KR-20 is a special case of coefficient \eqn{\alpha} (Guttman,
#' 1945; Cronbach, 1951), the same considerations apply: KR-20 equals
#' the population reliability of the composite under essential
#' \eqn{\tau}-equivalence (i.e., equal factor loadings) and serves as a
#' lower bound otherwise. For dichotomous items see also
#' \code{\link{reliability_omega_categorical}} (categorical omega), which treats
#' the items via a probit-link single-factor model and does not assume
#' equal loadings.
#'
#' Available confidence interval methods (set via \code{ci_method}) are
#' the same as for \code{\link{reliability_alpha}}; see that function's
#' \emph{Details} for full descriptions. The default for KR-20 is
#' \code{"feldt"}, the \emph{F}-distribution interval originally
#' developed for KR-20.
#'
#' The menu also includes the nonparametric bootstrap intervals
#' \code{"percentile"}, \code{"bca"}, \code{"bootstrap_se"}, and
#' \code{"bootstrap_se_logistic"}, which resample the rows of
#' \code{data} with replacement \code{B} times and recompute KR-20 on
#' each replication (Efron & Tibshirani, 1993). They are worth the cost
#' when the closed forms are least trustworthy, which for dichotomous
#' items means highly unbalanced item difficulties or a sample size too
#' small for the normal-theory derivations behind \code{"feldt"} and
#' \code{"bonett"}. No bootstrap runs unless \code{ci_method} asks for
#' one; when it does, the default is \code{B = 10000} replications, and
#' supplying \code{seed} makes the interval reproducible.
#'
#' \strong{Comparison with other packages.} The \pkg{psych} package
#' computes the same quantity via \code{\link[psych]{alpha}} (since
#' \eqn{\alpha} on 0/1 data \emph{is} KR-20). \code{reliability_kr20}
#' restricts input to raw 0/1 data so the dichotomous-items assumption
#' cannot be quietly violated, presents the historical formula in the
#' documentation, and accompanies the point estimate with a confidence
#' interval drawn from the methods compared in Kelley and Pornprasertmanit
#' (2016).
#'
#' @param data A numeric matrix or data frame of 0/1 item scores (rows
#'   are respondents, columns are items). Rows with any missing values
#'   are listwise-deleted. Non-binary values trigger an error.
#' @param ci_method Method for constructing the confidence interval; see
#'   \code{\link{reliability_alpha}} for the full list. Defaults to
#'   \code{"feldt"}.
#' @param conf_level Confidence level for the interval (1 - Type I error
#'   rate). Defaults to \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected. Defaults to \code{10000}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, which leaves the user's current RNG
#'   state intact; supply an integer for reproducibility.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#' and rows
#' \code{"estimate"} (sample KR-20),
#' \code{"se"} (standard error on the coefficient scale, \code{NA} for
#' methods that do not produce one; for the transformation-based
#' intervals \code{"fisher"}, \code{"bonett"}, and
#' \code{"hakstian_whalen"} it is the delta method back-transform of the
#' transformation-scale standard error),
#' \code{"se_transformed"} (only for those transformation-based
#' intervals: the standard error on the transformation scale, with the
#' scale named in the attribute \code{se_transform_scale}:
#' \code{"fisher_z"}, \code{"log(1-alpha)"}, or \code{"cube_root"}),
#' \code{"lower_limit"} and \code{"upper_limit"} (clamped to [0, 1]),
#' \code{"conf_level"}, \code{"N"} (effective sample size after
#' listwise deletion), \code{"N_complete"} (the complete cases; equal
#' to \code{"N"} here, and carried so the whole reliability family
#' returns one shape), and \code{"J"} (number of items). Attributes
#' \code{coefficient} (\code{"kr20"}) and \code{ci_method} record the
#' computation; bootstrap calls also record \code{B}.
#'
#' @references
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure
#'   of tests. \emph{Psychometrika, 16}(3), 297--334.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Feldt, L. S. (1965). The approximate sampling distribution of
#'   Kuder-Richardson reliability coefficient twenty.
#'   \emph{Psychometrika, 30}, 357--370.
#'
#' Guttman, L. (1945). A basis for analyzing test-retest reliability.
#'   \emph{Psychometrika, 10}(4), 255--282.
#'
#' Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
#'   formation for reliability coefficients of homogeneous measurement
#'   instruments. \emph{Methodology, 8}, 39--50.
#'   \doi{10.1027/1614-2241/a000036}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Kuder, G. F., & Richardson, M. W. (1937). The theory of the estimation
#'   of test reliability. \emph{Psychometrika, 2}, 151--160.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{reliability}} (general wrapper),
#' \code{\link{reliability_alpha}},
#' \code{\link{reliability_omega_categorical}},
#' \code{\link[psych]{alpha}}.
#'
#' @examples
#' set.seed(113)
#' # Ten dichotomous items with a single underlying ability.
#' N <- 300
#' J <- 10
#' ability <- rnorm(N)
#' loadings <- rep(0.6, J)
#' latent <- outer(ability, loadings) +
#'           matrix(rnorm(N * J, sd = sqrt(1 - 0.6^2)), N, J)
#' items <- (latent > 0) * 1
#' colnames(items) <- paste0("y", seq_len(J))
#'
#' reliability_kr20(data = items)
#' reliability_kr20(data = items, ci_method = "bonett")
#'
#' # A bootstrap interval recomputes KR-20 on each of B resamples of the
#' # rows, so it is shown rather than run; the call is
#' #   reliability_kr20(data = items, ci_method = "percentile",
#' #                    B = 10000, seed = 113)
#' # and the default B = 10000 is what a reported interval deserves.
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_kr20 <- function(data,
                             ci_method = c("feldt", "bonett", "fisher",
                                           "hakstian_whalen",
                                           "ml", "ml_logistic",
                                           "adf", "adf_logistic",
                                           "bootstrap_se",
                                           "bootstrap_se_logistic",
                                           "percentile", "bca", "none"),
                             conf_level = 0.95,
                             B = 10000,
                             seed = NULL) {
  ci_method <- match.arg(ci_method)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (missing(data) || is.null(data)) {
    stop("'data' (a 0/1 matrix or data frame of item scores) is required.",
         call. = FALSE)
  }

  data <- as.data.frame(data)
  data <- data[stats::complete.cases(data), , drop = FALSE]
  data_mat <- as.matrix(data)
  if (!all(data_mat %in% c(0, 1))) {
    stop("KR-20 requires items scored 0/1; non-binary values were found.",
         call. = FALSE)
  }
  if (ncol(data_mat) < 2L) {
    stop("At least two items are required.", call. = FALSE)
  }

  N <- nrow(data_mat)
  J <- ncol(data_mat)
  S <- stats::cov(data_mat)
  estimate <- .kr20_from_data(data_mat)

  ci <- switch(
    ci_method,
    none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
    feldt = .ci_feldt(estimate, N, J, conf_level),
    fisher = .ci_fisher(estimate, N, conf_level),
    bonett = .ci_bonett(estimate, N, J, conf_level),
    hakstian_whalen = .ci_hakstian_whalen(estimate, N, J, conf_level),
    ml = .ci_alpha_ml(estimate, S, N, J, conf_level, logistic = FALSE),
    ml_logistic = .ci_alpha_ml(estimate, S, N, J, conf_level, logistic = TRUE),
    adf = .ci_alpha_adf(estimate, data_mat, conf_level, logistic = FALSE),
    adf_logistic = .ci_alpha_adf(estimate, data_mat, conf_level,
                                 logistic = TRUE),
    bootstrap_se = ,
    bootstrap_se_logistic = ,
    percentile = ,
    bca = .bootstrap_ci(
      data = data_mat,
      point_fn = function(d) .kr20_from_data(d),
      B = B, conf_level = conf_level, kind = ci_method, seed = seed
    )
  )

  .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = N, J = J,
    coefficient = "kr20",
    ci_method = ci_method,
    B = if (ci_method %in% c("bootstrap_se", "bootstrap_se_logistic",
                             "percentile", "bca")) B else NA_integer_,
    se_transformed = ci$se_transformed,
    se_transform_scale = ci$se_transform_scale
  )
}
