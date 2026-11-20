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
#' @param true_var_error The population error variance: the variance
#'   of the errors common to the groups, which the mean square error (the
#'   mean square within groups) of an analysis of variance estimates. It is
#'   a variance, not a standard deviation: a planning value stated as an
#'   error standard deviation is squared before it is given. Must be a
#'   single finite positive number. This is the one way to state the error
#'   variability, as in \code{\link{ss_aipe_c_sensitivity}} and
#'   \code{\link{ss_power_contrast}}.
#' @param c_weights The contrast weights
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage, (i.e., 1 - Type I error rate). A level above .999999998, whose tail areas are below 1e-9, is refused (see "Small tail areas" in \code{\link{ci_nc_t}}).
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#' @param true_kurtosis The population excess kurtosis of the errors (their fourth
#'   standardized moment less 3), on the scale \code{\link{kurtosis}}
#'   reports, so the value it gives for the residuals of a pilot study can be
#'   passed directly. It is used only with \code{assurance}, through the
#'   effective degrees of freedom described in Details: a negative value
#'   (lighter tails than the normal, as bounded outcomes such as proportions
#'   and rating scale scores often have) makes the sample variance less
#'   variable than normal theory implies, and a positive value (heavier tails)
#'   makes it more variable. The default \code{0} is the normal theory
#'   value. Must be a single finite number
#'   greater than -2; no distribution has a smaller excess kurtosis, and only
#'   a symmetric two-point distribution attains -2.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}:
#' \item{necessary_n_per_group}{the necessary sample size \emph{per group}}
#'
#' A plan made with \code{assurance} and a nonzero \code{true_kurtosis} keeps
#' that one row, carries the value as the attribute \code{true_kurtosis}, prints
#' a line under the table saying that the assurance allows for that excess
#' kurtosis of the errors, and reports it as a \code{true_kurtosis} column of
#' \code{generics::glance()}. With the default \code{true_kurtosis = 0}, or
#' without \code{assurance}, no attribute is added.
#'
#' @details
#' Without \code{assurance}, the planned \emph{n} per group is the smallest
#' integer (and at least 2) satisfying
#' \deqn{n \ge \frac{4\, t^2_{1-\alpha/2,\, J(n-1)}\, \sigma^2 \sum_j c_j^2}{\omega^2},}
#' where \eqn{\sigma^2} is the error variance (\code{true_var_error}),
#' \eqn{c_j} are the contrast
#' weights of the \eqn{J} groups, and \eqn{\omega} is the desired full
#' width. That is the width the interval would have if the pooled sample
#' variance \eqn{S^2} equaled \eqn{\sigma^2}; because \eqn{S^2} varies from
#' sample to sample, the realized width exceeds \eqn{\omega} in a large
#' proportion of studies, a proportion that approaches one half as \emph{n}
#' grows. With \code{assurance} equal to \eqn{\gamma}, the right side is
#' multiplied by \eqn{\chi^2_{\gamma,\, \nu} / \nu}, the \eqn{\gamma}
#' quantile of a chi square on \eqn{\nu} degrees of freedom divided by
#' \eqn{\nu}, which is the \eqn{\gamma} quantile of \eqn{S^2/\sigma^2} when
#' \eqn{\nu S^2/\sigma^2} is chi square on \eqn{\nu} degrees of freedom; the
#' realized width is then at most \eqn{\omega} with probability
#' \eqn{\gamma}. For normal errors \eqn{\nu = J(n - 1)} and the chi square
#' law is exact.
#'
#' For errors that are not normal, the variability of \eqn{S^2} depends on
#' their excess kurtosis \eqn{\kappa}. With \eqn{J} groups of \eqn{n}, the
#' pooled variance has
#' \eqn{\mathrm{var}(S^2) = \sigma^4 [2/\{J(n - 1)\} + \kappa/(nJ)]}, so,
#' for an assurance above one half, a plan made under normality delivers
#' more assurance than requested when \eqn{\kappa < 0} and less when
#' \eqn{\kappa > 0}. The \code{true_kurtosis} argument replaces \eqn{J(n - 1)}
#' in the chi square quantile with the effective degrees of freedom
#' \deqn{\nu = \frac{2}{2/\{J(n - 1)\} + \kappa/(nJ)},}
#' the degrees of freedom for which \eqn{\chi^2_\nu/\nu} has the variance
#' of \eqn{S^2/\sigma^2}; \eqn{\kappa = 0} gives \eqn{\nu = J(n - 1)} and
#' the normal theory plan, and the degrees of freedom of the \emph{t}
#' quantile are not changed. This is the approximation of
#' Box (1953, p. 319, Eq. 3): a variance estimate on \eqn{\phi}
#' degrees of freedom from a population with excess kurtosis \eqn{\gamma_2}
#' has variance \eqn{2\sigma^4/(d\phi)}, with
#' \eqn{d = [1 + \phi\gamma_2/\{2(\phi + 1)\}]^{-1}}, and is distributed
#' approximately as a normal theory estimate on \eqn{d\phi} degrees of
#' freedom; with \eqn{\phi = n - 1} in each group, \eqn{\nu} is \eqn{d}
#' times \eqn{J(n - 1)}. The approximation matches two moments and leaves
#' the shape to the chi square. It requires errors with a finite fourth
#' moment, and in Box's sampling experiments (p. 325) the large sample
#' result was approached more slowly when the kurtosis was large. The
#' kurtosis of a pilot study's residuals is itself an imprecise estimate in
#' samples of the size pilots have, so the value supplied is a planning
#' value, as the error variance is. The sensitivity companion,
#' \code{\link{ss_aipe_c_sensitivity}}, simulates normal errors and has no
#' \code{true_kurtosis} argument, so it evaluates a plan under the normal theory
#' law, not under the errors a nonzero \code{true_kurtosis} allows for.
#'
#' A Monte Carlo check of the approximation used two groups,
#' \code{c_weights = c(1, -1)}, an error variance of 1, a width of .8, and
#' an assurance of .80, with 10,000 replications per condition. With
#' uniform errors (excess kurtosis -1.2), the plan made with the kurtosis
#' supplied (53 per group) met the width target in a proportion .811 of
#' the replications, where the plan made under normality (55 per group)
#' met it in .915; with \emph{t} errors on 10 degrees of freedom (excess
#' kurtosis 1), the plan with the kurtosis supplied (56 per group) met it
#' in .811, where the normal theory plan (55 per group) met it in .776,
#' short of the .80 requested. The study's script and results are kept with
#' the package's development records and are not distributed with the
#' package.
#'
#' @references
#' Box, G. E. P. (1953). Non-normality and tests on variances.
#'   \emph{Biometrika, 40}(3--4), 318--335.
#'   \doi{10.1093/biomet/40.3-4.318}
#'
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
#' ss_aipe_c(true_var_error = 40, c_weights = c(1, -0.5, -0.5),
#'           width = 3, assurance = .90)
#'
#' # The planning values can come from a pilot study. The posttest scores in
#' # depression_bdi have an error variance near 40, and their residuals are
#' # lighter tailed than the normal, as scores on a bounded inventory often
#' # are. Thirty residuals estimate the kurtosis imprecisely, so it is a
#' # planning value to weigh, as the error variance is.
#' pilot <- lm(bdi_post ~ condition, data = depression_bdi)
#' sigma(pilot)^2
#' kurtosis(residuals(pilot))
#' ss_aipe_c(true_var_error = sigma(pilot)^2, c_weights = c(1, -0.5, -0.5),
#'           width = 3, assurance = .90, true_kurtosis = kurtosis(residuals(pilot)))
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_c <- function(true_var_error, c_weights, width, conf_level = .95, assurance = NULL,
                      true_kurtosis = 0) {
    ####################################################################
    # The error variance is given one way, as true_var_error, the name
    # ss_aipe_c_sensitivity() and ss_power_contrast() use for it.
    if (missing(true_var_error) || is.null(true_var_error))
        .dmar_refuse_argument(c(
            "'true_var_error' must be given: the population error ",
            "variance, which the mean square error (the mean square within ",
            "groups) of an analysis of variance estimates. The width of the ",
            "interval for the contrast is proportional to its square root, ",
            "so no sample size can be planned without it. It is a variance: ",
            "if the planning value is an error standard deviation, give its ",
            "square."),
            arg = "true_var_error")
    if (!is.numeric(true_var_error) || length(true_var_error) != 1L || !is.finite(true_var_error) || true_var_error <= 0)
        .dmar_refuse_argument(c(
            "'true_var_error' must be a single finite positive number: ",
            "the population error variance (the mean square error estimates ",
            "it). A variance is positive, and the plan is for one contrast ",
            "under one common error variance. If the planning value is an ",
            "error standard deviation, give its square."),
            arg = "true_var_error", value = true_var_error)

    if (abs(sum(c_weights)) > 1e-8) .dmar_refuse_argument(
      "'c_weights' must sum to zero.", arg = "c_weights")
    # All-zero weights define no contrast; they gave a plan of 2 per group.
    .dmar_check_contrast_nonzero(c_weights)
    if (sum(c_weights[c_weights > 0]) > 1) .dmar_refuse_argument(
      "The positive values of 'c_weights' must not sum to more than 1; give the weights as fractions, such as c(1/2, 1/2, -1).", arg = "c_weights")

    # Validate the planning targets at entry so that a boundary input produces a
    # clear error rather than a NaN or a nonsensical N. A negative or infinite
    # width, or a conf_level outside (0, 1), would otherwise send qnorm()/qt()
    # out of domain and crash the fixed-point search with a missing value.
    if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
        .dmar_refuse_argument(
          "'width' must be a single finite positive number.", arg = "width")
    if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
        .dmar_refuse_argument(
          "'conf_level' must be a single number strictly between 0 and 1.",
          arg = "conf_level")
    # The floor on a tail area (R/dmar_conditions.R), which ci_c(), the
    # interval planned for, keeps.
    .dmar_check_tail_floor((1 - conf_level) / 2, (1 - conf_level) / 2,
                           conf_level = conf_level)
    if (!is.null(assurance) && (!is.numeric(assurance) || length(assurance) != 1L || !is.finite(assurance) || assurance <= 0 || assurance >= 1))
        .dmar_refuse_argument(
          "'assurance' must be NULL or a single number strictly between 0 and 1.",
          arg = "assurance")
    # Excess kurtosis is at least -2 for every distribution (-2 only for a
    # symmetric two-point law); a finite value is required because the
    # effective degrees of freedom below need a finite fourth moment. It is
    # validated even without an assurance, where it has no effect, so that a
    # mistyped value is caught rather than silently ignored.
    if (!is.numeric(true_kurtosis) || length(true_kurtosis) != 1L || !is.finite(true_kurtosis) || true_kurtosis <= -2)
        .dmar_refuse_argument(
          "'true_kurtosis' (the excess kurtosis of the errors) must be a single finite number greater than -2.",
          arg = "true_kurtosis")

    #####################################################################
    alpha <- 1 - conf_level
    J <- length(c_weights)
    sigma <- sqrt(true_var_error)

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
                .dmar_refuse(c(
                  "The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible."),
                  "computation_failed")
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
                .dmar_refuse(c(
                  "The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible."),
                  "computation_failed")
            n_p <- n
            # See the note above: floor the per-group n at 2 in the degrees of
            # freedom so qt() and qchisq() never receive a nonpositive df.
            n_df <- max(n, 2)
            df <- n_df * J - J
            # The effective degrees of freedom of the pooled variance under
            # errors with excess kurtosis kappa, Box's (1953, p. 319, Eq. 3)
            # d times df: var(S^2) / sigma^4 = 2 / df + kappa / (n J), and the
            # chi square on nu df has the same variance after scaling by nu, so
            # nu = 2 / (2 / df + kappa / (n J)) = df / (1 + kappa df / (2 n J)).
            # The second form is the one computed: at kappa = 0 the divisor is
            # exactly 1, so nu is df to the last bit and the plan is the
            # normal theory plan of earlier releases, not merely close to it.
            # kappa > -2 keeps the divisor positive, since df / (n J) < 1.
            nu <- df / (1 + true_kurtosis * df / (2 * n_df * J))
            n <- ((sigma^2 * 4 * (qt(1 - alpha / 2, df))^2 * sum(c_weights^2)) / width^2) * (qchisq(assurance, nu) / nu)
            dif <- abs(n - n_p)
        }
        n <- max(ceiling(n), min_n)
        out <- .as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe")
        # A plan made for errors that are not normal records the excess
        # kurtosis it allowed for, which .ss_aipe_kurtosis_line() writes
        # under the table and glance() reports as a column. With
        # true_kurtosis = 0 nothing is added, so the object is the one earlier
        # releases returned.
        if (true_kurtosis != 0) attr(out, "true_kurtosis") <- true_kurtosis
        return(out)
    }
}

# .ss_aipe_kurtosis_line(x)
#
# The line written under an assurance plan of ss_aipe_c() made with a
# nonzero true_kurtosis (the `true_kurtosis` attribute), or NULL. Its caller is
# .dmar_tbl_footer.default() in R/dmar_output_helpers.R, through which
# print.dmar_tbl(), knit_print(), and as_kable() all reach it. The number
# is formatted for reading only; the attribute keeps full precision.
.ss_aipe_kurtosis_line <- function(x) {
    k <- attr(x, "true_kurtosis", exact = TRUE)
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k == 0) return(NULL)
    paste0("The assurance allows for errors with an excess kurtosis of ",
           trimws(formatC(k, digits = 3, format = "fg")),
           " (normal errors have 0).")
}
