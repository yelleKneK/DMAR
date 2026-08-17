#' Sample Size Planning for an ANOVA Contrast From the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description
#' Plans the sample size \emph{per group} so that the confidence interval for
#' an unstandardized contrast of means in a fixed effects analysis of variance
#' is sufficiently narrow, following the accuracy in parameter estimation
#' (AIPE) approach: the design goal is a contrast estimated with the precision
#' the research question requires, not merely one detected as nonzero.
#' AIPE sample size planning for ANOVA and ANCOVA contrasts is developed
#' in Lai and Kelley (2012).
#'
#' @param error_variance The common error variance; i.e., the mean square error
#' @param c_weights The contrast weights
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#' @param MSwithin An alias for \code{error_variance}
#' @param SD The standard deviation of the common error in ANOVA model
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}:
#' \item{necessary_n_per_group}{the necessary sample size \emph{per group}}
#'
#' @references
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample size planning.
#'   \emph{Evaluation and the Health Professions, 26}(3), 258--287.
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
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Be sure to use the error variance and not its square root (i.e., the standard deviation of the errors).
#'
#' @seealso
#' \code{\link{ss_aipe_sc}}, \code{\link{ss_aipe_c_ancova}}, \code{\link{ci_c}}
#'
#' @examples
#' # Suppose the population error variance of some three-group ANOVA model
#' # is believed to be 40. The researcher is interested in the difference
#' # between the mean of group 1 and the average of means of group 2 and 3.
#' # To plan the sample size so that, with 90 percent certainty, the
#' # obtained 95 percent full confidence interval width is no wider than 3:
#'
#' ss_aipe_c(error_variance = 40, c_weights = c(1, -0.5, -0.5),
#'           width = 3, assurance = .90)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_c <- function(error_variance = NULL, c_weights, width, conf_level = .95, assurance = NULL, MSwithin = NULL,
                      SD = NULL, ...) {
    ####################################################################
    if (is.null(error_variance) && is.null(MSwithin) && is.null(SD)) stop("You must specify the variability of the contrast (e.g., 'error_variance', 'MSwithin', or 'SD'")

    if (!is.null(error_variance) && !is.null(MSwithin)) {
        if (error_variance != MSwithin) stop("You provided discrepant information about the estimated standard deviation of the contrast")
    }
    if (!is.null(error_variance) && !is.null(SD)) {
        if (error_variance != SD^2) stop("You provided discrepant information about the estimated standard deviation of the contrast")
    }
    if (!is.null(MSwithin) && !is.null(SD)) {
        if (MSwithin != SD^2) stop("You provided discrepant information about the estimated standard deviation of the contrast")
    }

    if (is.null(error_variance) && !is.null(MSwithin)) error_variance <- MSwithin
    if (is.null(error_variance) && !is.null(SD)) error_variance <- SD^2

    if (!is.numeric(error_variance) || length(error_variance) != 1L || !is.finite(error_variance) || error_variance <= 0)
        stop("The error variance (i.e., 'error_variance', 'MSwithin', or 'SD^2') must be a single finite positive number.", call. = FALSE)

    if (abs(sum(c_weights)) > 1e-8) stop("The sum of the coefficients must be zero")
    if (sum(c_weights[c_weights > 0]) > 1) stop("Please use fractions to specify the contrast weights")

    # Validate the planning targets at entry so that a boundary input produces a
    # clear error rather than a NaN or a nonsensical N. A negative or infinite
    # width, or a conf_level outside (0, 1), would otherwise send qnorm()/qt()
    # out of domain and crash the fixed-point search with a missing value.
    if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
        stop("'width' must be a single finite positive number.", call. = FALSE)
    if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
        stop("'conf_level' must be a single number strictly between 0 and 1.", call. = FALSE)
    if (!is.null(assurance) && (!is.numeric(assurance) || length(assurance) != 1L || !is.finite(assurance) || assurance <= 0 || assurance >= 1))
        stop("'assurance' must be NULL or a single number strictly between 0 and 1.", call. = FALSE)

    #####################################################################
    alpha <- 1 - conf_level
    J <- length(c_weights)
    sigma <- sqrt(error_variance)

    # The admissible minimum per-group sample size: the ANOVA error degrees of
    # freedom are J * (n - 1), so n must be at least 2 for a positive error df.
    # A target so wide that the closed-form n falls below 2 is met at n = 2, the
    # smallest design that can estimate the error variance; return that rather
    # than a degenerate n = 1.
    min_n <- 2L
    # An explicit cap on the fixed-point search so a pathological input cannot
    # loop without terminating.
    max_iter <- 10000L

    n <- (sigma^2 * 4 * (qnorm(1 - alpha / 2))^2 * sum(c_weights^2)) / width^2
    tol <- 1e-6
    dif <- tol + 1

    if (is.null(assurance)) {
        iter <- 0L
        while (dif > tol) {
            iter <- iter + 1L
            if (iter > max_iter)
                stop("The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible.", call. = FALSE)
            n_p <- n
            # Floor the per-group n at 2 when forming the t degrees of freedom so
            # that n * J - J = J * (n - 1) stays positive; for a very wide target
            # the z-based start can be below 1 per group, which would otherwise
            # pass negative df to qt() and produce NaN.
            df <- max(n, 2) * J - J
            n <- (sigma^2 * 4 * (qt(1 - alpha / 2, df))^2 * sum(c_weights^2)) / width^2
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
                stop("The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible.", call. = FALSE)
            n_p <- n
            # See the note above: floor the per-group n at 2 in the degrees of
            # freedom so qt() and qchisq() never receive a nonpositive df.
            df <- max(n, 2) * J - J
            n <- ((sigma^2 * 4 * (qt(1 - alpha / 2, df))^2 * sum(c_weights^2)) / width^2) * (qchisq(assurance, df) / df)
            dif <- abs(n - n_p)
        }
        n <- max(ceiling(n), min_n)
        return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
}
