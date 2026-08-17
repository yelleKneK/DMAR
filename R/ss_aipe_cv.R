#' Sample Size Planning for the Coefficient of Variation Given the Goal of Accuracy in Parameter Estimation Approach to Sample Size Planning
#'
#' @description
#' Determines the necessary sample size so that the expected confidence interval width for the coefficient
#' of variation will be sufficiently narrow, optionally with a desired degree of certainty that the interval
#' will not be wider than desired. The population coefficient of variation may be given directly
#' as \code{C_of_V} or through \code{mu} and \code{sigma}, in which case \code{C_of_V} is taken
#' as \code{sigma / mu}. The value of \code{C_of_V} should be positive.
#'
#' @param C_of_V Population coefficient of variation on which the sample size procedure is based
#' @param width Desired (full) width of the confidence interval
#' @param conf_level Confidence interval coverage; 1-Type I error rate
#' @param assurance Value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g., .80, .90, .95)
#' @param mu Population mean (specified with \code{sigma} when \code{C_of_V} is not specified)
#' @param sigma Population standard deviation (specified with \code{mu} when \code{C_of_V} is not specified)
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit
#' @param \dots For modifying parameters of functions this function calls
#'
#' @return
#' Returns the necessary sample size given the input specifications.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
#'   variation with minimum risk: A sequential method for minimizing
#'   sampling error and study cost.
#'   \emph{Multivariate Behavioral Research, 51}(5), 627--648.
#'   \doi{10.1080/00273171.2016.1203279}
#'
#' Kelley, K. (2007). Sample size planning for the coefficient of variation
#'   from the accuracy in parameter estimation approach.
#'   \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_aipe_cv_sensitivity}}, \code{\link{cv}}
#'
#' @examples
#' # Suppose one wishes to have a confidence interval with an expected width of .10
#' # for a 99% confidence interval when the population coefficient of variation is .10.
#' ss_aipe_cv(C_of_V = .1, width = .1, conf_level = .99)
#'
#' # The same planning problem parameterized by the population mean and standard
#' # deviation: mu = 10 and sigma = 1 imply the same coefficient of variation, .10.
#' ss_aipe_cv(mu = 10, sigma = 1, width = .1, conf_level = .99)
#'
#' # Ensuring that the confidence interval will be sufficiently narrow with a 99\%
#' # certainty for the situation above.
#' ss_aipe_cv(C_of_V = .1, width = .1, conf_level = .99, assurance = .99)
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_cv <- function(C_of_V = NULL, width = NULL, conf_level = .95, assurance = NULL, mu = NULL, sigma = NULL, alpha_lower = NULL, alpha_upper = NULL, ...) {
  # The coefficient of variation can be supplied directly as 'C_of_V' or
  # through the population mean and standard deviation, in which case
  # C_of_V = sigma / mu. Exactly one of the two parameterizations is allowed;
  # accepting both would let them disagree silently.
  if (!is.null(C_of_V) && (!is.null(mu) || !is.null(sigma)))
    stop("Specify either 'C_of_V' or both 'mu' and 'sigma'; you cannot mix the two parameterizations.", call. = FALSE)
  if (is.null(C_of_V)) {
    if (is.null(mu) || is.null(sigma))
      stop("Specify either 'C_of_V' or both 'mu' and 'sigma'.", call. = FALSE)
    if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) ||
        !is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma))
      stop("'mu' and 'sigma' must each be a single finite number.", call. = FALSE)
    C_of_V <- sigma / mu
  }

  if (C_of_V <= 0) stop("The coefficient of variation ('C_of_V', or the ratio 'sigma' / 'mu') should be positive (see Chattopadhyay and Kelley, 2016)")

  # A non-positive or non-finite 'width' would make the search target
  # unreachable (the interval width can never fall to or below it), so the
  # increment-by-one search would run without terminating. Validate at entry.
  if (is.null(width) || !is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
    stop("'width' must be a single finite positive number.", call. = FALSE)

  if (!is.null(conf_level) && (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1))
    stop("'conf_level' must be NULL or a single number strictly between 0 and 1.", call. = FALSE)

  if (is.null(conf_level)) {
    if (alpha_lower >= 1 || alpha_lower < 0) stop("\'alpha_lower\' is not correctly specified.")
    if (alpha_upper >= 1 || alpha_upper < 0) stop("\'alpha_upper\' is not correctly specified.")
  }

  if (!is.null(conf_level)) {
    if (!is.null(alpha_lower) || !is.null(alpha_upper)) stop("Since \'conf_level\' is specified, \'alpha_lower\' and \'alpha_upper\' should be \'NULL\'.")
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (!is.null(assurance)) {
    if ((assurance <= 0) || (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than .50 and less than 1.", call. = FALSE)
    if (assurance <= .50) stop("The 'assurance' should be > .5 (but less than 1).", call. = FALSE)
  }

  #####
  minimal_N <- 4

  Lim_0 <- ci_cv(cv = cv(cv = C_of_V, N = minimal_N, unbiased = TRUE)$value, n = minimal_N, alpha_lower = alpha_lower, alpha_upper = alpha_upper, conf_level = NULL)
  Current_Width <- Lim_0$value[Lim_0$term == "upper_limit"] - Lim_0$value[Lim_0$term == "lower_limit"]
  dif <- Current_Width - width
  N_0 <- minimal_N
  # Seed N with the minimal sample size so that if the requested width is
  # already attained at minimal_N (dif <= 0), the loop never runs and N is
  # still defined. Without this, N would be looked up in the enclosing scopes
  # and either error ("object 'N' not found") or silently return a stray
  # workspace variable named N.
  N <- minimal_N

  # The iterative search calls ci_cv() (and through it conf_limits_nct())
  # once per candidate N. Because the noncentrality parameter grows with
  # sqrt(N), it commonly exceeds 37.62, the value beyond which R's noncentral
  # t is inaccurate, and conf_limits_nct() warns on each such call. Surfacing
  # that warning once per iteration produces dozens of identical messages, so
  # we muffle the per-iteration warnings, count them, and emit a single
  # summary warning after the search completes.
  .ncp_count <- 0L
  # An explicit cap on the increment-by-one search so a pathological input
  # cannot loop without terminating.
  .iter <- 0L
  max_iter <- 100000L
  withCallingHandlers(
    while (dif > 0) {
      .iter <- .iter + 1L
      if (.iter > max_iter)
        stop("The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible.", call. = FALSE)
      N <- N_0 + 1
      CI_for_CV <- ci_cv(cv = cv(cv = C_of_V, N = N, unbiased = TRUE)$value, n = N, alpha_lower = alpha_lower, alpha_upper = alpha_upper, conf_level = NULL)
      Current_Width <- CI_for_CV$value[CI_for_CV$term == "upper_limit"] - CI_for_CV$value[CI_for_CV$term == "lower_limit"]
      dif <- Current_Width - width
      N_0 <- N
    },
    warning = function(w) {
      if (grepl("noncentrality parameter exceeds 37.62", conditionMessage(w))) {
        .ncp_count <<- .ncp_count + 1L
        invokeRestart("muffleWarning")
      }
    }
  )
  if (.ncp_count > 0L) {
    warning(sprintf(
      "During the iterative sample size search, the noncentrality parameter exceeded 37.62 in magnitude (the limit of R's noncentral t accuracy) in %d intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_nct.",
      .ncp_count
    ), call. = FALSE)
  }

  #############################################################################################

  # Now, incorporate a desired degree of certainty.
  if (!is.null(assurance)) {
    # Here the noncentral value that will be exceeded only (1-assurance)100% of the time is found (which leads to confidence intervals wider than desired).
    beyond_CV_NCP <- suppressWarnings(qt(p = 1 - assurance, df = N - 1, ncp = sqrt(N) / C_of_V, lower.tail = TRUE, log.p = FALSE))

    # Now the noncentral parameter is transformed into a coefficient of variation.
    Lim_for_Certainty <- sqrt(N) / beyond_CV_NCP

    # Now calculate sample size using the value not to be exceeded more than (1-assurance)100% of the time.
    N_gamma <- ss_aipe_cv(C_of_V = cv(cv = Lim_for_Certainty, N = N, unbiased = TRUE)$value, width = width, alpha_lower = alpha_lower, alpha_upper = alpha_upper, conf_level = NULL, assurance = NULL)[1, 2]
  }

  if (is.null(assurance)) {
    return(.as_dmar_tbl(data.frame(term = "necessary_N", value = N), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }

  if (!is.null(assurance)) {
    return(.as_dmar_tbl(data.frame(term = "necessary_N", value = N_gamma), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }
}
