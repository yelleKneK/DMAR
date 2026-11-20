#' Sample Size Planning for a Contrast in Randomized ANCOVA From the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description
#' Plans the sample size per group so that the confidence interval for an
#' unstandardized contrast in a one-covariate randomized ANCOVA is
#' sufficiently narrow, following the accuracy in parameter estimation (AIPE)
#' approach. To the extent the covariate correlates with the response, the
#' covariate adjustment shrinks the error variance, so the desired precision
#' is reached with a smaller sample size than the corresponding ANOVA design
#' requires.
#'
#' @param true_var_error_ancova The population error variance of the ANCOVA model (i.e., the mean square within of the ANCOVA model)
#' @param true_var_error_anova The population error variance of the ANOVA model (i.e., the mean square within of the ANOVA model)
#' @param true_cor The population correlation coefficient of the response and the covariate
#' @param c_weights The contrast weights
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#'
#' @details
#' Either the error variance of the ANCOVA model or of the ANOVA model can be used to plan the appropriate
#' sample size per group. When using the error variance of the ANOVA model to plan sample size, the correlation
#' coefficient of the response and the covariate is also needed.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}:
#' \item{necessary_n_per_group}{The necessary sample size \emph{per group}}
#'
#' @references
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample-size planning.
#'   \emph{Evaluation & the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#'   ANCOVA and ANOVA contrasts: Sample size planning via narrow
#'   confidence intervals.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 65},
#'   350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_c_ancova}}, \code{\link{ci_sc_ancova}}, \code{\link{ss_aipe_c}}
#'
#' @examples
#' # Suppose the population error variance of some three-group ANOVA model
#' # is believed to be 40, and the population correlation coefficient
#' # of the response and the covariate is 0.22. The researcher is
#' # interested in the difference between the mean of group 1 and
#' # the average of means of group 2 and 3. To plan the sample size so
#' # that, with 90 percent certainty, the obtained 95 percent full
#' # confidence interval width is no wider than 3:
#'
#' ss_aipe_c_ancova(true_var_error_anova = 40, true_cor = .22, c_weights = c(1, -0.5, -0.5),
#'                  width = 3, assurance = .90)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_c_ancova <- function(true_var_error_ancova = NULL, true_var_error_anova = NULL, true_cor = NULL, c_weights, width, conf_level = .95,
                             assurance = NULL) {
  if (is.null(true_var_error_ancova) && is.null(true_var_error_anova)) .dmar_refuse_argument(
    "Please specify either the ANCOVA error variance, or both the ANOVA error variance and the correlation coefficient",
    arg = NULL)

  if (is.null(true_var_error_ancova)) {
    if (is.null(true_var_error_anova) || is.null(true_cor)) .dmar_refuse_argument(
      "Please specify either the ANCOVA error variance, or both the ANOVA error variance and the correlation coefficient",
      arg = NULL)
    if (!is.numeric(true_cor) || length(true_cor) != 1L || !is.finite(true_cor) || abs(true_cor) >= 1)
      .dmar_refuse_argument(
        "'true_cor' must be a single number strictly between -1 and 1.",
        arg = "true_cor")
    if (!is.numeric(true_var_error_anova) || length(true_var_error_anova) != 1L || !is.finite(true_var_error_anova) || true_var_error_anova <= 0)
      .dmar_refuse_argument(
        "'true_var_error_anova' must be a single finite positive number.",
        arg = "true_var_error_anova")
    error_var <- true_var_error_anova * (1 - true_cor^2)
  }
  if (!is.null(true_var_error_ancova)) {
    if (!is.null(true_var_error_anova) || !is.null(true_cor)) .dmar_refuse_argument(
      "Since you input the ANCOVA error variance, do not input the ANOVA error variance and the correlation coefficient",
      arg = NULL)
    if (!is.numeric(true_var_error_ancova) || length(true_var_error_ancova) != 1L || !is.finite(true_var_error_ancova) || true_var_error_ancova <= 0)
      .dmar_refuse_argument(
        "'true_var_error_ancova' must be a single finite positive number.",
        arg = "true_var_error_ancova")
    error_var <- true_var_error_ancova
  }

  if (abs(sum(c_weights)) > 1e-8) .dmar_refuse_argument(
    "'c_weights' must sum to zero.", arg = "c_weights")
  # All-zero weights define no contrast; they gave a plan of 2 per group.
  .dmar_check_contrast_nonzero(c_weights)
  if (sum(c_weights[c_weights > 0]) > 1) .dmar_refuse_argument(
    "The positive values of 'c_weights' must not sum to more than 1; give the weights as fractions, such as c(1/2, 1/2, -1).", arg = "c_weights")

  # Validate the planning targets at entry so a boundary input errors clearly
  # rather than sending qnorm()/qt() out of domain and crashing the fixed-point
  # search with a missing value.
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single finite positive number.", arg = "width")
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument(
      "'conf_level' must be a single number strictly between 0 and 1.",
      arg = "conf_level")
  if (!is.null(assurance) && (!is.numeric(assurance) || length(assurance) != 1L || !is.finite(assurance) || assurance <= 0 || assurance >= 1))
    .dmar_refuse_argument(
      "'assurance' must be NULL or a single number strictly between 0 and 1.",
      arg = "assurance")

  alpha <- 1 - conf_level
  J <- length(c_weights)
  sigma <- sqrt(error_var)

  # The admissible minimum per-group sample size: the ANCOVA error degrees of
  # freedom are J * (n - 1) - 1 (one covariate), so n must be at least 2 for a
  # positive error df. A target so wide that the closed-form n falls below 2 is
  # met at n = 2; return that rather than a degenerate n = 1 (or a crash).
  min_n <- 2L
  max_iter <- 10000L

  # Floor the per-group n at 2 when forming the t / chi square degrees of
  # freedom so that J * (n - 1) - 1 stays positive; for a very wide target the
  # z-based start can be below 1 per group, which would otherwise pass a
  # nonpositive df to qt()/qchisq() and produce NaN.
  df_of <- function(n) max(n, 2) * J - J - 1

  n <- (sigma^2 * 4 * (qnorm(1 - alpha / 2))^2 * sum(c_weights^2)) / width^2
  tol <- 1e-6
  dif <- tol + 1

  if (is.null(assurance)) {
    iter <- 0L
    while (dif > tol) {
      iter <- iter + 1L
      if (iter > max_iter)
        .dmar_refuse(c(
          "The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible."),
          "computation_failed")
      n_p <- n
      n <- (sigma^2 * 4 * (qt(1 - alpha / 2, df_of(n)))^2 * sum(c_weights^2)) / width^2
      dif <- abs(n - n_p)
    }
    n <- max(ceiling(n), min_n)
    return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }

  if (!is.null(assurance)) {
    iter <- 0L
    while (dif > tol) {
      iter <- iter + 1L
      if (iter > max_iter)
        .dmar_refuse(c(
          "The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible."),
          "computation_failed")
      n_p <- n
      df_n <- df_of(n)
      n <- ((sigma^2 * 4 * (qt(1 - alpha / 2, df_n))^2 * sum(c_weights^2)) / width^2) * (qchisq(assurance, df_n) / df_n)
      dif <- abs(n - n_p)
    }
    n <- max(ceiling(n), min_n)
    return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }
}
