#' Sample Size Planning for Accuracy in Parameter Estimation for Reliability Coefficients
#'
#' @description
#' Computes the necessary sample size for the confidence interval on a
#' population reliability coefficient (coefficient alpha or coefficient
#' omega, depending on the assumed measurement model) to have expected
#' width no larger than \code{width}, or (when \code{assurance} is
#' supplied) to be no wider than \code{width} with the specified
#' probability. This is the accuracy in parameter estimation (AIPE)
#' counterpart to power based planning for reliability and is the
#' companion of \code{\link{reliability}}. The closed form planner uses
#' the Terry and Kelley (2012) formulas under the user-selected measurement
#' model and confidence interval type; when \code{assurance} is supplied,
#' the function follows the closed form with an internal Monte Carlo
#' simulation to find the smallest \emph{N} delivering the requested
#' assurance.
#'
#' Coefficient alpha (Guttman, 1945; subsequently popularized by
#' Cronbach, 1951) is returned for the parallel and tau equivalent
#' (\emph{i.e.}, the so called True Score) measurement models;
#' coefficient omega (McDonald, 1999) is returned for the congeneric
#' model.
#'
#' @details
#' The Monte Carlo assurance search simulates covariance matrices from
#' the population implied by the inputs and, at each candidate sample
#' size, computes the realized confidence interval width with the same
#' machinery the estimation side of the package uses. For
#' \code{type = "Factor Analytic"} the single-factor model is fit by
#' maximum likelihood (equal loadings for the parallel and tau
#' equivalent models, free loadings for the congeneric model) and the
#' interval is the delta method Wald interval on the model implied
#' reliability, the interval of
#' \code{\link{reliability_omega}(denominator = "model_implied",
#' ci_method = "ml")} and of
#' \code{\link{reliability_alpha}(estimator = "model_implied",
#' ci_method = "ml")}. For \code{type = "Normal Theory"} the interval
#' uses the van Zyl, Neudecker, and Nel (2000) closed form standard
#' error for the tau equivalent model and its compound symmetry
#' simplification for the parallel model, the same closed form behind
#' \code{reliability_alpha(ci_method = "ml")}. The congeneric model has
#' no normal theory form, so \code{type = "Normal Theory"} is an error
#' there.
#'
#' @param model The measurement model assumed for the population.
#'   Accepts (case-sensitive aliases shown in parentheses):
#'   \code{"Parallel"} (\code{"parallel"}, \code{"SB"},
#'   \code{"Spearman Brown"}, \code{"Spearman-Brown"}, \code{"sb"}) for
#'   the strictly parallel items model; \code{"True Score"}
#'   (\code{"True Score Equivalent"}, \code{"True-Score Equivalent"},
#'   \code{"Equivalent"}, \code{"Tau Equivalent"}, \code{"tau-equivalent"},
#'   \code{"Tau-Equivalent"}, \code{"True-Score"}, \code{"true-score"},
#'   \code{"true score"}, \code{"Cronbach"}, \code{"cronbach"},
#'   \code{"Chronbach"}, \code{"alpha"}) for the tau equivalent model
#'   (in which case the function plans for coefficient alpha); or
#'   \code{"Congeneric"} (\code{"congeneric"}, \code{"omega"},
#'   \code{"Omega"}) for the congeneric model (in which case the function
#'   plans for coefficient omega).
#' @param type The method used to construct the confidence interval on the
#'   reliability coefficient: either \code{"Factor Analytic"} (McDonald,
#'   1999), available for all three measurement models, or
#'   \code{"Normal Theory"} (van Zyl, Neudecker, & Nel, 2000), available
#'   for the parallel and tau equivalent models only.
#' @param width The desired full width of the two-sided confidence interval.
#' @param S A symmetric population covariance (or correlation) matrix among
#'   the items, used to imply the population reliability and its sampling
#'   distribution.
#' @param conf_level Confidence level (i.e., \eqn{1 - \alpha}, where
#'   \eqn{\alpha} is the Type I error rate). Default \code{0.95}.
#' @param assurance Optional probability with which the realized interval
#'   is to be no wider than \code{width}. When \code{NULL} (the default),
#'   the planner targets the \emph{expected} width; when supplied
#'   (e.g., 0.80, 0.85, 0.95), the function follows the closed form with a
#'   Monte Carlo search.
#' @param data A data set from which the population covariance matrix
#'   should be inferred.
#' @param i Number of items.
#' @param cor_est The presumed inter-item correlation. One value for the
#'   parallel and tau equivalent models.
#' @param lambda Vector of population factor loadings.
#' @param psi_square Vector of population unique (error) variances.
#' @param initial_iter Number of Monte Carlo iterations used in the
#'   initial assurance search.
#' @param final_iter Number of Monte Carlo iterations used in the final
#'   assurance verification.
#' @param start_ss Optional starting sample size for the iterative
#'   assurance search.
#' @param verbose If \code{TRUE}, prints the current sample size and
#'   empirical assurance at each step of the Monte Carlo search.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}.
#' Without \code{assurance} the data frame has a single row,
#' \code{"necessary_N"}, giving the necessary \emph{N}. With \code{assurance}
#' supplied, the data frame has five rows: \code{"necessary_N"} (necessary
#' \emph{N}), \code{"width"} (echo of the target width),
#' \code{"specified_assurance"} (echo of the requested probability),
#' \code{"empirical_assurance"} (the assurance achieved at the returned
#' \emph{N} in the Monte Carlo verification), and \code{"final_iter"}
#' (number of Monte Carlo iterations used).
#'
#' @references
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
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. Mahwah, NJ:
#'   Lawrence Erlbaum Associates.
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{British Journal of Mathematical
#'   and Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution of
#'   the maximum likelihood estimator of Cronbach's alpha.
#'   \emph{Psychometrika, 65}(3), 271--280. \doi{10.1007/BF02296146}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note Not all of the items can be entered into the function to represent the population values.
#' For example, either 'data' can be used, or \code{S}, or  \code{i}, \code{cor_est}, and \code{psi_square},
#' or \code{i}, \code{lambda}, and \code{psi_square}. With a large number of iterations (\code{final_iter})
#' this function may take considerable time.
#'
#' @section Warning: In some conditions the factor analytic model fit by
#' \code{\link{cfa_1}} (via \pkg{lavaan}) may fail to converge, and you may
#' see a non-convergence message from \pkg{lavaan}. The Monte Carlo assurance
#' search treats a non-converged iteration as missing and continues, so a few
#' such messages do not invalidate the result. Frequent non-convergence
#' usually means the model is poorly determined by the data, for example
#' because of a small sample size, a low number of iterations, or a poorly
#' behaved covariance matrix.
#'
#' @seealso \code{\link{reliability}}, \code{\link{cfa_1}}
#'
#' @examples
#' # Expected confidence interval width (closed form, no Monte Carlo search).
#' ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
#'   i = 6, cor_est = .3, psi_square = .2, conf_level = .95, assurance = NULL)
#'
#' # The assurance cases run a Monte Carlo search; the iteration counts below are
#' # reduced so the example runs quickly. Raise initial_iter and final_iter for a
#' # production plan.
#' set.seed(113)
#'
#' # Same population, now targeting an assurance.
#' ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
#'   i = 6, cor_est = .3, psi_square = .2, conf_level = .95, assurance = .85,
#'   initial_iter = 50, final_iter = 200)
#'
#' # The true score (tau equivalent) model takes psi_square as a vector of
#' # length i (number of items) while cor_est stays a single value. Its
#' # assurance search is the slowest of the normal theory calls on this page,
#' # and the S matrix example at the end already plans for the true score
#' # model, so this one is shown rather than run:
#' #   ss_aipe_reliability(model = "True Score", type = "Normal Theory",
#' #     width = .1, i = 5, cor_est = .3, psi_square = c(.2, .3, .3, .2, .3),
#' #     conf_level = .95, assurance = .85, initial_iter = 50,
#' #     final_iter = 200)
#'
#' # Congeneric model, planned from the item loadings and error variances rather
#' # than from a single correlation. With assurance = NULL the necessary N comes
#' # from the closed form expected width evaluated at the implied population
#' # correlation matrix, so type does not enter the answer; type selects the
#' # interval that the Monte Carlo assurance search evaluates.
#' ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic", width = .15,
#'   i = 4, lambda = c(.8, .7, .7, .8), psi_square = c(.4, .5, .5, .4),
#'   conf_level = .95, assurance = NULL)
#'
#' # Adding an assurance to that congeneric plan is the expensive case: with the
#' # factor analytic interval a one factor model is fit at every Monte Carlo
#' # iteration and at every candidate sample size the search visits, so it runs
#' # for tens of seconds at the reduced counts used above and for many minutes at
#' # the defaults. That is why it is not run here; the call is
#' #   ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
#' #     width = .15, i = 4, lambda = c(.8, .7, .7, .8),
#' #     psi_square = c(.4, .5, .5, .4), conf_level = .95, assurance = .80,
#' #     initial_iter = 50, final_iter = 200)
#'
#' # Planning from a presumed population correlation matrix among the items.
#' pop_mat <- rbind(
#'   c(1.0000000, 0.3813850, 0.4216370, 0.3651484, 0.4472136),
#'   c(0.3813850, 1.0000000, 0.4020151, 0.3481553, 0.4264014),
#'   c(0.4216370, 0.4020151, 1.0000000, 0.3849002, 0.4714045),
#'   c(0.3651484, 0.3481553, 0.3849002, 1.0000000, 0.4082483),
#'   c(0.4472136, 0.4264014, 0.4714045, 0.4082483, 1.0000000))
#' ss_aipe_reliability(model = "True Score", type = "Normal Theory", width = .15,
#'   S = pop_mat, conf_level = .95, assurance = .85, initial_iter = 50,
#'   final_iter = 200)
#'
#' @keywords design
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_reliability <- function(model = NULL, type = NULL, width = NULL, S = NULL,
                                conf_level = 0.95, assurance = NULL, data = NULL, i = NULL,
                                cor_est = NULL, lambda = NULL, psi_square = NULL, initial_iter = 500,
                                final_iter = 5000, start_ss = NULL, verbose = FALSE) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")

  model_type1 <- c(
    "Parallel", "SB", "Spearman Brown", "Spearman-Brown",
    "sb", "parallel"
  )
  model_type2 <- c(
    "True Score", "True Score Equivalent", "True-Score Equivalent",
    "Equivalent", "Tau Equivalent", "Chronbach", "Cronbach",
    "cronbach", "alpha", "true score", "tau-equivalent",
    "Tau-Equivalent", "True-Score", "true-score"
  )
  model_type3 <- c("Congeneric", "congeneric", "omega", "Omega")
  if (sum(model == model_type1, model == model_type2, model ==
    model_type3) != 1) {
    stop("Assign one and only one of the three types of models to 'model'.")
  }
  if (sum(model == model_type1) == 1) {
    Model_to_Use <- "Parallel"
  }
  if (sum(model == model_type2) == 1) {
    Model_to_Use <- "True Score"
  }
  if (sum(model == model_type3) == 1) {
    Model_to_Use <- "Congeneric"
  }
  type_1 <- c(
    "Normal Theory", "Normal theory", "normal theory",
    "nt", "NT"
  )
  type_2 <- c(
    "Factor Analytic", "factor analytic", "Factor analytic",
    "factor Analytic"
  )
  if (sum(type == type_1, type == type_2) != 1) {
    stop("Assign either Factor Analytic or Normal Theory to 'type'.")
  }
  if (sum(type == type_1) == 1) {
    Type_to_Use <- "Normal Theory"
  }
  if (sum(type == type_2) == 1) {
    Type_to_Use <- "Factor Analytic"
  }
  if (sum(model == model_type1) == 1) {
    if (!is.null(cor_est)) {
      if (!is.null(lambda)) {
        stop("Please enter either cor_est or lambda, but not both.")
      }
      if (is.null(psi_square)) {
        stop("Problem: please enter all of the necessary information: `i', `cor_est' or `lambda', and `psi.square.'")
      }
      if (psi_square <= 0) {
        stop("Problem: `psi_square' must be greater than zero")
      }
      if (length(psi_square) > 1) {
        stop("Problem: 'psi_square' must be one number for the Parallel model. If you want to enter multiple 'psi_square' values, you must use either `True Score' or `Congeneric.'")
      }
      if (is.null(i)) {
        stop("Problem: please enter all of the necessary information: i, cor_est or lambda, and psi.square.")
      }
      if (i <= 0) {
        stop("Problem: `i' must be greater than zero")
      }
      if (length(cor_est) >= 2) {
        stop("Problem: you can only enter one 'cor_est' value.")
      }
      lambda_1 <- sqrt(cor_est)
      v_lambda <- matrix(data = lambda_1, nrow = 1, ncol = i)
      v_Psi <- matrix(data = psi_square, nrow = 1, ncol = i)
      temp_mat <- covmat_from_cfa(
        lambda = v_lambda, psi_squared = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$population_cov
    }
    if (!is.null(lambda)) {
      if (!is.null(cor_est)) {
        stop("Problem: please enter either cor_est or lambda, but not both.")
      }
      if (is.null(psi_square)) {
        stop("Problem: please enter all of the necessary information: i, cor_est or lambda, and psi.square.")
      }
      if (psi_square <= 0) {
        stop("Problem: `psi_square' must be greater than zero.")
      }
      if (length(psi_square) > 1) {
        stop("Problem: 'psi_square' must be one number for the Parallel model. If you want to enter multiple 'psi_square' values, you must use either the True Score or Congeneric model.")
      }
      if (is.null(i)) {
        stop("Problem: please enter all of the necessary information: i, cor_est or lambda, and psi.square.")
      }
      if (i <= 0) {
        stop("Problem: `i' must be greater than zero.")
      }
      v_lambda <- matrix(data = lambda, nrow = 1, ncol = i)
      v_Psi <- matrix(data = psi_square, nrow = 1, ncol = i)
      temp_mat <- covmat_from_cfa(
        lambda = v_lambda, psi_squared = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$population_cov
    }
    if (!is.null(S)) {
      if (!isSymmetric(S, tol = 1e-05)) {
        stop("Input a symmetric covariance or correlation matrix, 'S'")
      }
      cov_mat <- S
      i <- dim(S)[1]
    }
    if (!is.null(data)) {
      data <- na.omit(data)
      cov_mat <- var(data, y = NULL, na.rm = TRUE)
    }
    sigma_jj <- sum(diag(cov_mat))
    sigma2_Y <- sum(cov_mat)
    p <- ncol(cov_mat)
    alpha <- (p / (p - 1)) * (1 - sigma_jj / sigma2_Y)
    k <- 1 - alpha
    Crit_Value <- qnorm(1 - (1 - conf_level) / 2)
    Nec_N <- as.numeric(((Crit_Value^2) * 8 * (k^2) * p) / ((width^2) *
      (p - 1)) + 1)
  }
  if (sum(model == model_type2) == 1) {
    if (!is.null(cor_est)) {
      if (!is.null(lambda)) {
        stop("Problem: please enter either `cor_est' or `lambda', but not both")
      }
      if (is.null(psi_square)) {
        stop("Problem: please enter a vector of `psi_square' values")
      }
      if (is.null(i)) {
        stop("Problem: please enter a value for i")
      }
      if ((i) != length(psi_square)) {
        stop("The number of values entered in the psi_square vector should be the same as the quantity entered for i")
      }
      if (i <= 0) {
        stop("Problem: `i' must be greater than zero")
      }
      if (length(cor_est) >= 2) {
        stop("You can only enter one 'cor_est' value.")
      }
      lambda_1 <- sqrt(cor_est)
      lambda_vector <- rep(lambda_1, times = i)
      Population_Cov <- covmat_from_cfa(
        lambda = lambda_vector,
        psi_squared = psi_square
      )$population_cov
      cor_mat <- cov2cor(Population_Cov)
      cov_mat <- cor_mat
    }
    if (!is.null(lambda)) {
      if (!is.null(cor_est)) {
        stop("Please enter either cor_est or lambda, but not both")
      }
      if (is.null(psi_square)) {
        stop("Please enter all of the necessary information: i, cor_est or lambda, and psi_square")
      }
      if (length(psi_square) != (i)) {
        stop("The number of values entered in the psi_square vector should be the same as the quantity entered for i")
      }
      if (is.null(i)) {
        stop("You need to enter a quantity for i if you enter lambda and psi_square")
      }
      if (i <= 0) {
        stop("i must be greater than zero")
      }
      if (length(lambda) > 1) {
        stop("'lambda' must be one number for the True Score model. If you want to enter multiple 'lambda' values, you must use the Congeneric model.")
      }
      lambda_vector <- rep(lambda, times = i)
      Population_Cov <- covmat_from_cfa(
        lambda = lambda_vector,
        psi_squared = psi_square
      )$population_cov
      cor_mat <- cov2cor(Population_Cov)
      cov_mat <- cor_mat
    }
    if (!is.null(S)) {
      if (!isSymmetric(S, tol = 1e-05)) {
        stop("Input a symmetric covariance or correlation matrix 'S'")
      }
      cor_mat <- cov2cor(S)
      cov_mat <- cor_mat
      i <- dim(S)[1]
    }
    if (!is.null(data)) {
      data <- na.omit(data)
      cor_mat <- cor(data)
      cov_mat <- cor_mat
    }
    p <- ncol(cor_mat)
    j <- cbind(rep(1, times = p))
    Crit_Value <- qnorm(1 - (1 - conf_level) / 2)
    step_1 <- (p^2 / (p - 1)^2)
    gamma_1 <- 2 / ((t(j) %*% cor_mat %*% j)^3)
    gamma_2_1_1 <- (t(j) %*% cor_mat %*% j)
    gamma_2_1_2 <- ((sum(diag(cor_mat %*% cor_mat))) + (sum(diag(cor_mat)))^2)
    gamma_2_1 <- gamma_2_1_1 * gamma_2_1_2
    gamma_2_2 <- 2 * (sum(diag(cor_mat))) * (t(j) %*% (cor_mat %*%
      cor_mat) %*% j)
    gamma_2 <- gamma_2_1 - gamma_2_2
    gamma_final <- gamma_1 * gamma_2
    Nec_N <- as.numeric(((step_1 * gamma_final) / ((width / (2 *
      Crit_Value))^2)) + 1)
  }
  if (sum(model == model_type3) == 1) {
    if (!is.null(cor_est)) {
      if (!is.null(lambda)) {
        stop("Please enter either cor_est or lambda, but not both")
      }
      if (is.null(psi_square)) {
        stop("Please enter all of the necessary information: i, cor_est or lambda, and psi_square")
      }
      if (is.null(i)) {
        stop("Please enter all of the necessary information: i, cor_est or lambda, and psi_square")
      }
      if (i <= 0) {
        stop("i must be greater than zero")
      }
      if ((i) != length(psi_square)) {
        stop("The number of values entered in the psi_square vector should be the same as the quantity entered for i")
      }
      if (length(cor_est) >= 2) {
        stop("If you have multiple values for 'cor_est', please put them as 'lambda' values. The square root of cor_est equals lambda.")
      }
      # print("You entered one value for 'cor_est' with the Congeneric model. This model allows for multiple 'lambda' values. If you have only one 'cor_est' or 'lambda' value, you might want to use the True Score model.")
      lambda_1 <- sqrt(cor_est)
      v_lambda <- matrix(data = lambda_1, nrow = 1, ncol = i)
      v_Psi <- matrix(data = psi_square, nrow = 1, ncol = i)
      temp_mat <- covmat_from_cfa(
        lambda = v_lambda, psi_squared = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$population_cov
      cor_mat <- cov2cor(cov_mat)
    }
    if (!is.null(lambda)) {
      if (!is.null(cor_est)) {
        stop("Please enter either cor_est or lambda, but not both")
      }
      if (is.null(psi_square)) {
        stop("Please enter all of the necessary information: i, cor_est or lambda, and psi_square")
      }
      if (is.null(i)) {
        stop("Please enter all of the necessary information: i, cor_est or lambda, and psi_square")
      }
      if (i <= 0) {
        stop("i must be greater than zero")
      }
      if ((i) != length(lambda)) {
        stop("The number of values entered in the lambda vector should be the same as the quantity entered for i")
      }
      if ((i) != length(psi_square)) {
        stop("The number of values entered in the psi_square vector should be the same as the quantity entered for i")
      }
      v_lambda <- matrix(data = lambda, nrow = 1, ncol = i)
      v_Psi <- matrix(data = psi_square, nrow = 1, ncol = i)
      temp_mat <- covmat_from_cfa(
        lambda = v_lambda, psi_squared = v_Psi,
        tol_det = 1e-05
      )
      cov_mat <- Sigma <- temp_mat$population_cov
      cor_mat <- cov2cor(cov_mat)
    }
    if (!is.null(data)) {
      data <- na.omit(data)
      cov_mat <- var(data, y = NULL, na.rm = TRUE)
      cor_mat <- cov2cor(cov_mat)
    }
    if (!is.null(S)) {
      if (!isSymmetric(S, tol = 1e-05)) {
        stop("Input a symmetric covariance or correlation matrix 'S'")
      }
      cor_mat <- cov2cor(S)
      cov_mat <- cor_mat
      i <- dim(S)[1]
    }
    p <- ncol(cor_mat)
    j <- cbind(rep(1, times = p))
    Crit_Value <- qnorm(1 - (1 - conf_level) / 2)
    step_1 <- (p^2 / (p - 1)^2)
    gamma_1 <- 2 / ((t(j) %*% cor_mat %*% j)^3)
    gamma_2_1_1 <- (t(j) %*% cor_mat %*% j)
    gamma_2_1_2 <- ((sum(diag(cor_mat %*% cor_mat))) + (sum(diag(cor_mat)))^2)
    gamma_2_1 <- gamma_2_1_1 * gamma_2_1_2
    gamma_2_2 <- 2 * (sum(diag(cor_mat))) * (t(j) %*% (cor_mat %*%
      cor_mat) %*% j)
    gamma_2 <- gamma_2_1 - gamma_2_2
    gamma_final <- gamma_1 * gamma_2
    Nec_N <- as.numeric(((step_1 * gamma_final) / ((width / (2 *
      Crit_Value))^2)) + 1)
  }
  if (!is.null(assurance)) {
    if (assurance > 1) {
      assurance <- assurance / 100
    }
    # print("An a priori Monte Carlo simulation study has been started so that the exact sample size for the requested condition can be determined. Please be patient, as this process may take several minutes (or longer given the computer and condition).")
    initial_assurance_N <- ceiling(Nec_N) + 1
    if (sum(model == model_type1) == 1) {
      Model_to_Use <- "Parallel"
    }
    if (sum(model == model_type2) == 1) {
      Model_to_Use <- "True Score"
    }
    if (sum(model == model_type3) == 1) {
      Model_to_Use <- "Congeneric"
    }
    n_i <- initial_assurance_N
    if (is.null(start_ss)) {
      Difference <- -1
      while (Difference < 0) {
        CI_Result <- rep(0, initial_iter)
        for (iters in 1:initial_iter) {
          sim_data <- var(MASS::mvrnorm(n = n_i, mu = rep(
            0,
            i
          ), Sigma = cov_mat, tol = 1e-06, empirical = FALSE))
          CI_Result_raw <- try(.ss_aipe_reliability_ci(
            S = sim_data,
            N = n_i, model = Model_to_Use, type = Type_to_Use,
            conf_level = conf_level
          ), silent = TRUE)
          if (inherits(CI_Result_raw, "try-error")) {
            CI_Result[iters] <- NA
          }
          if (!inherits(CI_Result_raw, "try-error")) {
            CI_Result[iters] <- CI_Result_raw$upper -
              CI_Result_raw$lower
          }
          Difference <- mean(na.omit(CI_Result) < width) - assurance
          if (verbose == TRUE) print(paste("The current assurance is", mean(na.omit(CI_Result) < width), "at the current sample size of", n_i, "you are in the 'initial_iter' stage."))
          if (Difference < 0) break
        }
        if (Difference < 0) {
          n_i <- n_i + 1
        }
      }
      initial_n_i <- n_i
    }
    if (!is.null(start_ss)) {
      CI_Result <- rep(0, final_iter)
      n_i <- start_ss
      Difference <- 0
      for (iters in 1:final_iter) {
        sim_data <- var(MASS::mvrnorm(n = n_i, mu = rep(
          0,
          i
        ), Sigma = cov_mat, tol = 1e-06, empirical = FALSE))
        CI_Result_raw <- try(.ss_aipe_reliability_ci(
          S = sim_data,
          N = n_i, model = Model_to_Use, type = Type_to_Use,
          conf_level = conf_level
        ), silent = TRUE)
        if (inherits(CI_Result_raw, "try-error")) {
          CI_Result[iters] <- NA
        }
        if (!inherits(CI_Result_raw, "try-error")) {
          CI_Result[iters] <- CI_Result_raw$upper -
            CI_Result_raw$lower
        }
        Difference <- mean(na.omit(CI_Result) < width) - assurance
        if (verbose == TRUE) print(paste("The current assurance is", mean(na.omit(CI_Result) < width), "at the current sample size of", n_i, "you are in the 'final_iter' stage."))
        if (Difference < 0) break
      }

      if (Difference > 0) {
        while (Difference > 0) {
          CI_Result <- rep(1, final_iter)
          for (iters in 1:final_iter) {
            sim_data <- var(MASS::mvrnorm(n = n_i, mu = rep(
              0,
              i
            ), Sigma = cov_mat, tol = 1e-06, empirical = FALSE))
            CI_Result_raw <- try(.ss_aipe_reliability_ci(
              S = sim_data,
              N = n_i, model = Model_to_Use, type = Type_to_Use,
              conf_level = conf_level
            ), silent = TRUE)
            if (inherits(CI_Result_raw, "try-error")) {
              CI_Result[iters] <- NA
            }
            if (!inherits(CI_Result_raw, "try-error")) {
              CI_Result[iters] <- CI_Result_raw$upper -
                CI_Result_raw$lower
            }
            Difference <- mean(na.omit(CI_Result) < width) - assurance
            if (verbose == TRUE) print(paste("The current assurance is", mean(na.omit(CI_Result) < width), "at the current sample size of", n_i, "you are in the 'final_iter' stage."))
            if (Difference > 0) break
          }

          if (Difference > 0) {
            n_i <- n_i - 1
          }
        }
      }
    }
    Difference <- -1
    while (Difference < 0) {
      CI_Result <- rep(0, final_iter)
      for (iters in 1:final_iter) {
        sim_data <- var(MASS::mvrnorm(n = n_i, mu = rep(
          0,
          i
        ), Sigma = cov_mat, tol = 1e-06, empirical = FALSE))
        CI_Result_raw <- try(.ss_aipe_reliability_ci(
          S = sim_data,
          N = n_i, model = Model_to_Use, type = Type_to_Use,
          conf_level = conf_level
        ), silent = TRUE)
        if (inherits(CI_Result_raw, "try-error")) {
          CI_Result[iters] <- NA
        }
        if (!inherits(CI_Result_raw, "try-error")) {
          CI_Result[iters] <- CI_Result_raw$upper -
            CI_Result_raw$lower
        }
        Difference <- mean(na.omit(CI_Result) < width) - assurance
        if (verbose == TRUE) print(paste("The current assurance is", mean(na.omit(CI_Result) < width), "at the current sample size of", n_i, "; you are in the 'final_iter' stage."))
        if (Difference < 0) break
      }

      if (Difference < 0) {
        n_i <- n_i + 1
      }
    }
    Nec_N_assurance <- n_i
    empirical_assurance <- mean(na.omit(CI_Result) < width)
    # print(paste("A sample size of", n_i, "leads to an empirical assurance of", round(mean(na.omit(CI_Result) < width), 3)))
  }
  if (is.null(assurance)) {
    return(.as_dmar_tbl(data.frame(term = "necessary_N", value = ceiling(Nec_N)), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }
  if (!is.null(assurance)) {
    return(.as_dmar_tbl(data.frame(
      term = c("necessary_N", "width", "specified_assurance", "empirical_assurance", "final_iter"),
      value = c(ceiling(Nec_N_assurance), width, assurance, empirical_assurance, final_iter)
    ), conf_level = conf_level))
  }
}


# .ss_aipe_reliability_ci(S, N, model, type, conf_level)
#
# Confidence interval for the reliability coefficient implied by a sample
# covariance (or correlation) matrix, used by the Monte Carlo assurance
# search in ss_aipe_reliability(). The caller has already normalized
# `model` to one of "Parallel", "True Score", or "Congeneric" and `type`
# to "Factor Analytic" or "Normal Theory".
#
# The computation is the package's own reliability machinery, the same
# engine the estimation side of DMAR uses:
#
#   - "Factor Analytic": .omega_fit_cfa() fits the single-factor model by
#     maximum likelihood (equal loadings for the parallel and tau
#     equivalent models, free loadings for the congeneric model) and
#     .ci_omega_delta() builds the delta method Wald interval on the
#     omega defined parameter. This is the interval of
#     reliability_omega(denominator = "model_implied", ci_method = "ml")
#     and of reliability_alpha(estimator = "model_implied",
#     ci_method = "ml"). Note that the factor analytic treatment of the
#     parallel model matches the historical planner engine: equal
#     loadings are imposed but the error variances stay free, so the
#     parallel and tau equivalent models share the same factor analytic
#     interval.
#
#   - "Normal Theory": the van Zyl, Neudecker, and Nel (2000) closed
#     form standard error from .ci_alpha_ml() for the tau equivalent
#     model, and its compound symmetry simplification,
#     Var(alpha_hat) = 2 (1 - alpha)^2 J / ((N - 1) (J - 1)), for the
#     parallel model. The congeneric model has no normal theory form
#     and is an error.
#
# Both interval types are clamped to [0, 1], matching the estimation
# functions (reliability is bounded by definition, and the planner
# compares realized widths against a target width on that scale).
#
# Inputs:
#   S          : symmetric covariance or correlation matrix among items.
#   N          : sample size that generated S.
#   model      : "Parallel", "True Score", or "Congeneric".
#   type       : "Factor Analytic" or "Normal Theory".
#   conf_level : confidence level.
#
# Returns:
#   A list with `estimate`, `se`, `lower`, and `upper`. Stops (rather
#   than returning NA) when the factor analytic fit does not converge;
#   the planner's try() treats that iteration as missing.
.ss_aipe_reliability_ci <- function(S, N, model, type, conf_level) {
  if (type == "Factor Analytic") {
    equal_loadings <- model %in% c("Parallel", "True Score")
    fit <- .omega_fit_cfa(S = S, N = as.numeric(N),
                          equal_loadings = equal_loadings, se = "standard")
    if (!fit$converged) {
      stop("The single-factor model for coefficient omega did not converge.",
           call. = FALSE)
    }
    cl <- .ci_omega_delta(fit, conf_level = conf_level, logistic = FALSE)
    return(list(estimate = fit$omega, se = cl$se,
                lower = max(cl$lower, 0), upper = min(cl$upper, 1)))
  }
  if (type == "Normal Theory") {
    if (model == "Congeneric") {
      stop("The Congeneric model can not be used with `Normal Theory.'")
    }
    p <- ncol(S)
    estimate <- .alpha_from_S(S)
    if (model == "Parallel") {
      crit <- stats::qnorm(1 - (1 - conf_level) / 2)
      se <- sqrt((2 * (1 - estimate)^2 * p) / ((N - 1) * (p - 1)))
      cl <- list(se = se,
                 lower = estimate - crit * se,
                 upper = estimate + crit * se)
    } else {
      cl <- .ci_alpha_ml(estimate, S = S, N = N, J = p,
                         conf_level = conf_level, logistic = FALSE)
    }
    return(list(estimate = estimate, se = cl$se,
                lower = max(cl$lower, 0), upper = min(cl$upper, 1)))
  }
  stop("Unknown interval type: ", type, call. = FALSE)
}
