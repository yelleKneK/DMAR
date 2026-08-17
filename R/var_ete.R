#' Variance of the Estimated Treatment Effect in Two-Group ANCOVA With
#' Heterogeneous Slopes
#'
#' Computes the variance of the estimated treatment effect (ETE) at a
#' chosen covariate value in a two-group analysis of covariance with
#' heterogeneity of regression and a random covariate, following Li,
#' McLouth, and Delaney (2020). When the two groups' slopes differ, the
#' treatment effect is a function of the covariate, and its sampling
#' variance at the sample grand mean, one standard deviation from the
#' mean, or a fixed covariate value must account for the covariate
#' being a random variable rather than a set of fixed constants; the
#' fixed-constant formulas understate or misstate that variability.
#' This is the reimplementation of \code{var.ete()} from \pkg{MBESS},
#' contributed there by Li Li.
#'
#' @param sigma2 Residual error variance: the population value when
#'   \code{type = "population"}, the sample estimate when
#'   \code{type = "sample"}.
#' @param sigma2_Z Variance of the random covariate: population value
#'   or sample estimate, matching \code{type}.
#' @param n_1,n_2 Sample sizes of the two groups (each must exceed 3;
#'   the formulas involve \eqn{n - 3} in denominators).
#' @param beta_1,beta_2 Slopes of the covariate in group 1 and
#'   group 2: population values or sample estimates, matching
#'   \code{type}.
#' @param mu_Z Mean of the covariate (population value or sample mean,
#'   matching \code{type}). Defaults to 0. Used when
#'   \code{covariate_value = "fixed"}.
#' @param fixed_value The fixed covariate value at which the treatment
#'   effect is assessed when \code{covariate_value = "fixed"}.
#'   Defaults to 0.
#' @param type \code{"sample"} (default) for the unbiased estimate of
#'   the variance from sample slopes and variances, or
#'   \code{"population"} for the variance from population values.
#' @param covariate_value Where the treatment effect is assessed:
#'   \code{"sample_mean"} (default) at the sample grand mean of the
#'   covariate, \code{"sd"} at the grand mean plus or minus one sample
#'   standard deviation, or \code{"fixed"} at \code{fixed_value}.
#'
#' @details
#' Randomized experiments with a covariate commonly probe the simple
#' treatment effect at the grand mean and one standard deviation
#' either side of it when the slopes differ across groups. The
#' variance expressions here treat the covariate as normally
#' distributed rather than fixed, which Li, McLouth, and Delaney
#' (2020) show can change the estimated standard error substantially
#' when heterogeneity of regression is strong. The square root of the
#' returned value is the standard error used for a confidence interval
#' or test of the treatment effect at the chosen covariate value.
#'
#' At the sample grand mean of the covariate, writing
#' \eqn{N = n_1 + n_2}, the population variance (their Equation 10) is
#' \deqn{\mathrm{Var} = \sigma^2 C_0
#'   + \frac{(\beta_1 - \beta_2)^2 \sigma^2_Z}{N}, \qquad
#'   C_0 = \frac{1}{n_1} + \frac{1}{n_2}
#'   + \frac{n_2}{N n_1 (n_1 - 3)} + \frac{n_1}{N n_2 (n_2 - 3)},}
#' and with \code{type = "sample"} the returned value is their unbiased
#' estimator (Equation C.7), which subtracts
#' \eqn{\sigma^2 \{(N-3)/(n_1-3) + (N-3)/(n_2-3)\} / \{N (N-1)\}}
#' so that plugging in sample estimates does not overstate the variance.
#' The \code{covariate_value = "sd"} expressions are their Equations 12
#' and C.9, which add the variance contribution of estimating the
#' covariate's standard deviation, and the \code{"fixed"} expressions are
#' their Equations 14 and C.10.
#'
#' The two \code{"fixed"} estimands differ in where the deviation of
#' \code{fixed_value} is measured from. With \code{type = "sample"} the
#' deviation is taken from the sample grand mean (Equation C.10), so
#' \code{mu_Z} should be the sample mean of the covariate. With
#' \code{type = "population"} the deviation is taken from a known
#' population mean (Equation 14); evaluating the treatment effect at that
#' known mean itself, as in the paper's worked example, sets
#' \code{fixed_value = mu_Z}, which zeroes the deviation term.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row,
#'   \code{term = "var_ete"}, whose \code{value} is the variance of
#'   the estimated treatment effect at the chosen covariate value. The
#'   \code{type} and \code{covariate_value} choices are recorded as
#'   attributes of the same names.
#'
#' @references
#' Kelley, K. (2007a). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of
#'   Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2007b). Methods for the behavioral, educational, and
#'   social sciences: An R package. \emph{Behavior Research Methods,
#'   39}(4), 979--984. \doi{10.3758/BF03192993}
#'
#' Li, L., McLouth, C. J., & Delaney, H. D. (2020). Analysis of
#'   covariance in randomized experiments with heterogeneity of
#'   regression and a random covariate: The variance of the estimated
#'   treatment effect at selected covariate values. \emph{Multivariate
#'   Behavioral Research, 55}(6), 926--940.
#'   \doi{10.1080/00273171.2019.1693953}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9 on heterogeneity of
#'   regression.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ancova}} for the model whose treatment effect
#'   this variance describes; \code{\link{regions_of_significance}}
#'   for the companion question of where a moderated effect is
#'   distinguishable from zero.
#'
#' @family variance utilities
#'
#' @keywords htest
#'
#' @examples
#' # Pygmalion data (Maxwell, Delaney, & Kelley, 2027): the treatment
#' # effect of the "Bloomer" expectation at the covariate grand mean,
#' # with heterogeneous pre-IQ slopes.
#' data(pygmalion)
#' fit <- lm(iq_8 ~ iq_pre * treatment, data = pygmalion)
#' s2  <- sum(residuals(fit)^2) / fit$df.residual
#' var_ete(sigma2 = s2, sigma2_Z = var(pygmalion$iq_pre),
#'         n_1 = sum(pygmalion$treatment == "Bloomer"),
#'         n_2 = sum(pygmalion$treatment == "Control"),
#'         beta_1 = coef(fit)["iq_pre"] + coef(fit)["iq_pre:treatmentBloomer"],
#'         beta_2 = coef(fit)["iq_pre"])
#'
#' @export
var_ete <- function(sigma2, sigma2_Z, n_1, n_2, beta_1, beta_2,
                    mu_Z = 0, fixed_value = 0,
                    type = c("sample", "population"),
                    covariate_value = c("sample_mean", "sd", "fixed")) {
  type <- match.arg(type)
  covariate_value <- match.arg(covariate_value)
  for (nm in c("sigma2", "sigma2_Z", "n_1", "n_2", "beta_1", "beta_2",
               "mu_Z", "fixed_value")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
      stop("'", nm, "' must be a single non-missing number.",
           call. = FALSE)
    }
  }
  if (sigma2 <= 0 || sigma2_Z <= 0) {
    stop("'sigma2' and 'sigma2_Z' must be positive.", call. = FALSE)
  }
  if (n_1 <= 3 || n_2 <= 3 || n_1 != round(n_1) || n_2 != round(n_2)) {
    stop("'n_1' and 'n_2' must be whole numbers greater than 3.",
         call. = FALSE)
  }
  n_1 <- as.numeric(n_1)
  n_2 <- as.numeric(n_2)
  N <- n_1 + n_2
  # Shared building blocks of Li, McLouth, and Delaney (2020).
  base <- 1 / n_1 + 1 / n_2 + n_2 / (N * n_1 * (n_1 - 3)) +
    n_1 / (N * n_2 * (n_2 - 3))
  tail_sum <- (N - 3) / (n_1 - 3) + (N - 3) / (n_2 - 3)
  slope2 <- (beta_1 - beta_2)^2

  if (covariate_value == "sample_mean") {
    variance <- sigma2 * base + sigma2_Z * slope2 / N
    if (type == "sample") {
      variance <- variance - sigma2 / N * (tail_sum / (N - 1))
    }
  } else if (covariate_value == "sd") {
    c_1 <- 1 / N + 1 - 2 / (N - 1) *
      exp(2 * (lgamma(N / 2) - lgamma((N - 1) / 2)))
    c_0 <- base + tail_sum / (N - 1)
    variance <- sigma2 * c_0 + sigma2_Z * slope2 * c_1
    if (type == "sample") {
      variance <- variance - c_1 * sigma2 * tail_sum / (N - 1)
    }
  } else {
    c_1 <- (n_1 - 2) / (n_1 * (n_1 - 3)) + (n_2 - 2) / (n_2 * (n_2 - 3))
    c_2 <- 1 / (n_1 - 3) + 1 / (n_2 - 3)
    dev2 <- (fixed_value - mu_Z)^2
    if (type == "population") {
      variance <- sigma2 * (c_1 + dev2 * c_2 / sigma2_Z)
    } else {
      variance <- sigma2 *
        (c_1 + dev2 * c_2 / sigma2_Z * (N - 3) / (N - 1) - c_2 / N)
    }
  }

  out <- data.frame(term = "var_ete", value = unname(variance),
                    stringsAsFactors = FALSE)
  out <- .as_dmar_tbl(out)
  attr(out, "type") <- type
  attr(out, "covariate_value") <- covariate_value
  out
}
