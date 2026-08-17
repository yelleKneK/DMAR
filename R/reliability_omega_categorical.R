#' Categorical Omega for Ordered-Categorical Items, With a Confidence Interval
#'
#' @description
#' Estimates categorical omega (\eqn{\omega_C}; Green & Yang, 2009;
#' Kelley & Pornprasertmanit, 2016) for a homogeneous composite of
#' ordered-categorical items and returns a bootstrap confidence interval.
#'
#' @details
#' Categorical omega is designed for items measured on an ordered
#' categorical scale (e.g., Likert items). It uses a probit-link
#' single-factor model in which each observed item \eqn{X_j} is modeled
#' as a categorization of an underlying continuous response variable
#' \eqn{X_j^{*}} via thresholds \eqn{t_{j,c}} (Muthén, 1984; Millsap &
#' Yun-Tein, 2004), fit by diagonally weighted least squares with mean-
#' and variance-adjusted test statistic (WLSMV), the standard estimator
#' for ordered categorical items. With the delta parameterization
#' (\eqn{Var(X_j^{*}) = 1}), the population categorical omega is
#' \deqn{\omega_C = \frac{\sum_{j=1}^{J} \sum_{j'=1}^{J} \sigma_{jj'}\!\left(\lambda_{j}\lambda_{j'}\right)}{\sum_{j=1}^{J} \sum_{j'=1}^{J} \sigma_{jj'}\!\left(\rho_{X_{j}^{*} X_{j'}^{*}}\right)},}
#' where \eqn{\sigma_{jj'}(r)} is the model implied covariance of
#' \eqn{(X_{j}, X_{j'})} computed from a bivariate normal CDF over pairs
#' of category thresholds and a correlation \eqn{r} (Green & Yang, 2009,
#' Eq. 13--14, with the full coefficient their Eq. 21; Kelley &
#' Pornprasertmanit, 2016, Eq. 17--18). The
#' numerator uses model implied polychoric correlations
#' (\eqn{\lambda_j \lambda_{j'}}), while the denominator uses observed
#' polychoric correlations estimated from the data via a saturated
#' bivariate model.
#'
#' Kelley and Pornprasertmanit (2016) found in extensive Monte Carlo
#' simulation that the bias-corrected and accelerated (BCa) bootstrap
#' confidence interval for categorical omega achieved acceptable coverage
#' across a wide variety of threshold patterns, sample sizes, item
#' counts, and population reliability values, with one documented
#' exception: coverage dipped somewhat below the acceptable range when
#' the number of items and the population reliability were both high.
#' They specifically recommend BCa for categorical omega. Because no bootstrap runs in
#' \pkg{DMAR} unless the user requests one, the default output is the
#' point estimate with a message naming the call that produces the
#' recommended interval; request \code{ci_method = "bca"} to obtain it.
#'
#' \strong{When to use.} Use \code{reliability_omega_categorical} when items are
#' ordered-categorical, especially when (a) the number of categories is
#' small (e.g., two to five), (b) item distributions are skewed, or (c)
#' threshold patterns differ markedly across items. In Kelley and
#' Pornprasertmanit's (2016) Study 3, treating ordered items as
#' continuous and using \code{\link{reliability_omega}} with the
#' observed total variance in the denominator achieved
#' acceptable coverage only when threshold patterns were similar across
#' items; in their experience that condition is rare in practice.
#'
#' Available confidence interval methods (set via \code{ci_method}).
#' Every interval here is bootstrap based: the rows of \code{data} are
#' resampled with replacement \code{B} times (10000 by default) and
#' categorical omega is recomputed, with the full WLSMV model refit, on
#' each replication (Efron & Tibshirani, 1993). Replications whose
#' refit fails or does not converge, most common with small samples and
#' sparse response categories, are dropped, and the interval is
#' computed from the replications that return a value. Bootstrap
#' results vary from run to run; supply \code{seed} for
#' reproducibility.
#' \describe{
#'   \item{\code{"bca"}}{The bias-corrected and accelerated bootstrap,
#'   the default and the specific recommendation of Kelley and
#'   Pornprasertmanit (2016). Where the percentile interval reads its
#'   limits directly off the empirical quantiles of the bootstrap
#'   estimates, BCa adjusts the two quantile positions for median bias
#'   (estimated from the bootstrap distribution) and for the rate at
#'   which the estimator's variance changes with the parameter (the
#'   acceleration, estimated by the jackknife), making it second-order
#'   accurate where the percentile interval is first-order accurate
#'   (DiCiccio & Efron, 1996). The adjusted quantile positions sit
#'   farther into the tails than the percentile interval uses, which
#'   is why the default \code{B = 10000} is larger than the customary
#'   2000; reduce \code{B} for exploration, not for a reported
#'   analysis.}
#'   \item{\code{"percentile"}}{Percentile bootstrap: the interval
#'   limits are the empirical quantiles of the bootstrap estimates.}
#'   \item{\code{"bootstrap_se"}, \code{"bootstrap_se_logistic"}}{Wald
#'   intervals using the bootstrap standard deviation as a standard
#'   error, built on the logit scale for the \code{_logistic} variant
#'   so the endpoints respect [0, 1].}
#'   \item{\code{"none"}}{Return only the point estimate.}
#' }
#'
#' \strong{Comparison with other packages.} The \pkg{psych} package's
#' \code{\link[psych]{omega}} fits a Schmid-Leiman hierarchical factor
#' model on continuous (or treated-as-continuous) items and does not
#' implement categorical omega in the sense of Green and Yang (2009).
#' For ordered-categorical items \code{reliability_omega_categorical} is the
#' appropriate choice; \code{\link[psych]{polychoric}} provides
#' polychoric correlation estimation as a separate tool.
#'
#' @param data A numeric matrix or data frame of ordered-categorical item
#'   scores (integer codes for the categories). Rows with any missing
#'   values are listwise-deleted.
#' @param ci_method Method for constructing the confidence interval. See
#'   \emph{Details}. When not supplied, no interval is computed: every
#'   interval for categorical omega is bootstrap based, and a bootstrap
#'   is never run unless requested. Ask for \code{"bca"} (the
#'   recommended method) or \code{"percentile"} to obtain one.
#' @param conf_level Confidence level for the interval. Defaults to
#'   \code{0.95}.
#' @param B Number of bootstrap replications. Defaults to \code{10000}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, which leaves the user's current RNG state intact; supply an integer for reproducibility.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#' and rows
#' \code{"estimate"} (sample \eqn{\omega_C}),
#' \code{"se"} (the bootstrap standard deviation across replications,
#' already on the coefficient scale; \code{NA} for
#' \code{ci_method = "none"}),
#' \code{"lower_limit"} and \code{"upper_limit"} (clamped to [0, 1]),
#' \code{"conf_level"}, \code{"N"}, \code{"N_complete"} (the complete
#' cases; equal to \code{"N"} here, and carried so the whole
#' reliability family returns one shape), and \code{"J"}. Attributes
#' \code{coefficient} (\code{"omega_categorical"}), \code{ci_method}, and \code{B}
#' record the computation.
#'
#' @references
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#'   \emph{Statistical Science, 11}(3), 189--228.
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores
#'   using structural equation modeling: An alternative to coefficient
#'   alpha. \emph{Psychometrika, 74}, 155--167.
#'   \doi{10.1007/s11336-008-9099-3}
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
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance
#'   in ordered-categorical measures.
#'   \emph{Multivariate Behavioral Research, 39}(3), 479--515.
#'   \doi{10.1207/s15327906mbr3903_4}
#'
#' Muthén, B. (1984). A general structural equation model with
#'   dichotomous, ordered categorical, and continuous latent variable
#'   indicators. \emph{Psychometrika, 49}(1), 115--132.
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
#' \code{\link{reliability_omega}} (use for continuous items),
#' \code{\link{reliability_kr20}} (for dichotomous items),
#' \code{\link[psych]{omega}},
#' \code{\link[psych]{polychoric}}.
#'
#' @examples
#' set.seed(113)
#' # Six 5-category items with a single latent factor.
#' N <- 200
#' J <- 6
#' loadings <- rep(0.7, J)
#' eta <- rnorm(N)
#' latent <- outer(eta, loadings) +
#'           matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
#' items <- apply(latent, 2, function(x)
#'   as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
#'                  labels = FALSE)))
#' colnames(items) <- paste0("y", seq_len(J))
#'
#' # Default: point estimate only, with a message naming the call that
#' # produces the recommended interval.
#' reliability_omega_categorical(data = items)
#'
#' # The same items treated as continuous, for contrast. That call fits a
#' # second factor analysis model, so it is shown rather than run:
#' #   reliability_omega(data = items)
#' # With five categories and thresholds spread across the latent scale
#' # the two coefficients nearly agree on these data; the gap widens as
#' # the categories get coarser and as the thresholds move into the
#' # tails, which is where the categorical coefficient is worth its cost.
#'
#' # Every interval for categorical omega is bootstrap based, and each
#' # replication refits the ordered-categorical factor analysis model.
#' # That refitting is why it is not run here; the call is
#' #   reliability_omega_categorical(data = items, ci_method = "bca",
#' #                                 B = 10000, seed = 113)
#' # and a reported interval deserves the BCa method at the default
#' # B = 10000.
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability_omega_categorical <- function(data,
                                ci_method = c("bca", "percentile",
                                              "bootstrap_se",
                                              "bootstrap_se_logistic",
                                              "none"),
                                conf_level = 0.95,
                                B = 10000,
                                seed = NULL) {
  ci_supplied <- !missing(ci_method)
  ci_method <- match.arg(ci_method)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (missing(data) || is.null(data)) {
    stop("'data' (a matrix or data frame of ordered-categorical item ",
         "scores) is required.", call. = FALSE)
  }

  data <- as.data.frame(data)
  data <- data[stats::complete.cases(data), , drop = FALSE]
  if (ncol(data) < 2L) {
    stop("At least two items are required.", call. = FALSE)
  }
  for (j in seq_len(ncol(data))) {
    if (!all(data[[j]] == as.integer(data[[j]]))) {
      stop("Categorical omega requires items coded as integer category ",
           "values; column ", j, " contains non-integer values.",
           call. = FALSE)
    }
  }

  N <- nrow(data)
  J <- ncol(data)

  if (!ci_supplied) {
    # Every interval for categorical omega is bootstrap based, and no
    # bootstrap runs unless the user asks for one. The message waits
    # until after input validation so a call that is about to error does
    # not first receive advice.
    ci_method <- "none"
    message("Categorical omega is reported without a confidence interval ",
            "by default because its interval is bootstrap based. Request ",
            "it with ci_method = \"bca\" (the recommended method) or ",
            "\"percentile\"; B = 10000 replications is the default when ",
            "you do.")
  }

  point_fn <- function(d) .omega_c_from_data(d)
  estimate <- point_fn(data)
  if (is.na(estimate)) {
    stop("The single-factor categorical CFA model did not converge.",
         call. = FALSE)
  }

  ci <- switch(
    ci_method,
    none = list(se = NA_real_, lower = NA_real_, upper = NA_real_),
    bootstrap_se = ,
    bootstrap_se_logistic = ,
    percentile = ,
    bca = .bootstrap_ci(
      data = data,
      point_fn = point_fn,
      B = B, conf_level = conf_level, kind = ci_method, seed = seed
    )
  )

  .relia_result(
    estimate = estimate, se = ci$se,
    lower = ci$lower, upper = ci$upper,
    conf_level = conf_level, N = N, J = J,
    coefficient = "omega_categorical",
    ci_method = ci_method,
    B = if (ci_method != "none") B else NA_integer_
  )
}

