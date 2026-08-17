#' Sample Size Planning for a Single Regression Coefficient (AIPE)
#'
#' @description
#' Computes the necessary sample size for the confidence interval on a
#' targeted regression coefficient \eqn{\beta_j} (or its standardized
#' counterpart) in a multiple regression with \emph{p} predictors to be no
#' wider than a user-specified value. This is the accuracy in parameter
#' estimation (AIPE) framework of Kelley and Maxwell (2003), targeted at a
#' specific coefficient rather than at the omnibus \eqn{R^2}. The
#' \code{noncentral = TRUE} variant inverts the noncentral \emph{t}
#' distribution of the standardized \eqn{b_j} under joint multivariate
#' normality of the predictors; the default central-\emph{t} variant uses
#' a closed form Wald style approximation. Optionally, supplying
#' \code{assurance} returns the larger \emph{N} that guarantees the
#' realized width with the specified probability rather than just on
#' average.
#'
#' @param rho2_Y_X Population value of \eqn{\rho^2_{Y \cdot X_1, \ldots, X_p}},
#'   the squared multiple correlation of the outcome \emph{Y} with all
#'   \emph{p} predictors.
#' @param rho2_j_X_without_j Population value of
#'   \eqn{\rho^2_{X_j \cdot X_{-j}}}, the squared multiple correlation when
#'   the \emph{j}th predictor is regressed on the remaining \eqn{p - 1}
#'   predictors. Quantifies the multicollinearity faced by the targeted
#'   coefficient.
#' @param p The number of predictor variables.
#' @param b_j The (unstandardized) regression coefficient for the \emph{j}th
#'   predictor, the predictor of interest.
#' @param width Desired (full) width of the two-sided confidence interval
#'   on \eqn{\beta_j}.
#' @param which_width Which portion of the confidence interval
#'   \code{width} refers to. Only \code{"Full"} is currently implemented.
#' @param sigma_Y Population standard deviation of \emph{Y}.
#' @param sigma_X Population standard deviation of the \emph{j}th predictor.
#' @param rho_XX Population correlation matrix for the \emph{p} predictor
#'   variables. If supplied with \code{rho_YX}, \code{rho2_Y_X} and
#'   \code{rho2_j_X_without_j} are derived from the covariance structure.
#' @param rho_YX Length-\emph{p} vector of population correlations between
#'   \emph{Y} and the \emph{p} predictors.
#' @param which_predictor Which of the \emph{p} predictors is the targeted
#'   coefficient.
#' @param noncentral If \code{TRUE}, plans using the exact noncentral
#'   \emph{t} sampling distribution of the standardized \eqn{b_j} (Kelley,
#'   2007). If \code{FALSE} (the default), uses the central \emph{t}
#'   approximation. The noncentral path requires
#'   \code{sigma_Y = sigma_X = 1} (\emph{i.e.}, a standardized solution).
#' @param alpha_lower Type I error rate for the lower confidence limit.
#' @param alpha_upper Type I error rate for the upper confidence limit.
#' @param conf_level Confidence level (i.e., \eqn{1 - \alpha}, where
#'   \eqn{\alpha} is the Type I error rate). Default \code{0.95}. Mutually
#'   exclusive with \code{alpha_lower} and \code{alpha_upper}.
#' @param assurance Optional probability with which the realized confidence
#'   interval is to be no wider than \code{width}. When \code{NULL} (the
#'   default), the planning targets the \emph{expected} width.
#'
#' @details
#' \strong{Calling conventions.} The function offers several mutually
#' exclusive ways to supply the population information needed to plan
#' \eqn{N}; the user picks the one most aligned with their available
#' planning values. Specify exactly one of:
#' \itemize{
#'   \item \emph{Covariance structure path.} Supply \code{rho_XX} and
#'     \code{rho_YX} along with \code{which_predictor}. The function
#'     derives \eqn{\rho^2_{Y \cdot X}} and
#'     \eqn{\rho^2_{X_j \cdot X_{-j}}} from the covariance structure and
#'     also computes the population \eqn{b_j} for consistency checks.
#'   \item \emph{Squared multiple correlations path.} Supply
#'     \code{rho2_Y_X}, \code{rho2_j_X_without_j}, \code{p}, and
#'     \code{b_j} directly when the user already has these planning values
#'     from prior research and does not need to specify the full
#'     covariance structure.
#'   \item \emph{Standardized solution.} For the \code{noncentral = TRUE}
#'     path, set \code{sigma_Y = sigma_X = 1}; the result returns the
#'     standardized regression coefficient sample size.
#' }
#'
#' \strong{Noncentral vs.\ central planning.} The central \emph{t} closed
#' form is fast and adequate at moderate to large \emph{N}; the
#' \code{noncentral = TRUE} path additionally accounts for the noncentral
#' \emph{t} sampling distribution of the standardized \eqn{b_j} and is
#' preferred when planning at small to moderate \emph{N} or when reporting
#' planning that will be matched against the noncentral CI from
#' \code{\link{ci_reg_coef}}.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} value is \code{"necessary_N"} and \code{value} is the
#' necessary total sample size \emph{N} given the input specifications.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 4 on individual
#'   comparisons of means and Chapter 6 on trend analysis.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_aipe_reg_coef_sensitivity}}, \code{\link{conf_limits_nct}}
#'
#' @examples
#' # 1. Covariance structure path: supply the population correlation
#' #    matrix and the population YX cross-correlations. Five predictors
#' #    in an exchangeable structure (all pairwise correlations 0.5,
#' #    all Y-X correlations 0.3).
#' rho_YX <- c(.3, .3, .3, .3, .3)
#' rho_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5),
#'                 c(.5, .5, 1, .5, .5), c(.5, .5, .5, 1, .5),
#'                 c(.5, .5, .5, .5, 1))
#'
#' # Closed-form (central t) planning, targeting the first predictor's
#' # standardized coefficient.
#' ss_aipe_reg_coef(width = .10, which_width = "Full",
#'                  sigma_Y = 1, sigma_X = 1,
#'                  rho_XX = rho_XX, rho_YX = rho_YX,
#'                  which_predictor = 1, noncentral = FALSE,
#'                  conf_level = .95)
#'
#' # Adding assurance (.85): the realized CI width will be no larger than
#' # 0.10 in 85 percent of replications. Required N grows accordingly.
#' ss_aipe_reg_coef(width = .10, which_width = "Full",
#'                  sigma_Y = 1, sigma_X = 1,
#'                  rho_XX = rho_XX, rho_YX = rho_YX,
#'                  which_predictor = 1, noncentral = FALSE,
#'                  conf_level = .95, assurance = .85)
#'
#' # Exact noncentral t planning. Required N differs at small to
#' # moderate samples.
#' ss_aipe_reg_coef(width = .10, which_width = "Full",
#'                  sigma_Y = 1, sigma_X = 1,
#'                  rho_XX = rho_XX, rho_YX = rho_YX,
#'                  which_predictor = 1, noncentral = TRUE,
#'                  conf_level = .95)
#'
#' # 2. Squared multiple correlations path: when the user has planning
#' #    values for rho^2_Y.X and rho^2_j.X_-j directly (e.g., from a
#' #    prior power analysis), without specifying the full covariance
#' #    structure. b_j is required on this path.
#' ss_aipe_reg_coef(rho2_Y_X = 0.30, rho2_j_X_without_j = 0.20,
#'                  p = 5, b_j = 0.25,
#'                  width = .15, which_width = "Full",
#'                  sigma_Y = 1, sigma_X = 1,
#'                  noncentral = FALSE, conf_level = .95)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_reg_coef <- function(rho2_Y_X = NULL, rho2_j_X_without_j = NULL, p = NULL,
                             b_j = NULL, width, which_width = "Full", sigma_Y = 1, sigma_X = 1,
                             rho_XX = NULL, rho_YX = NULL, which_predictor = NULL, noncentral = FALSE,
                             alpha_lower = NULL, alpha_upper = NULL, conf_level = 0.95,
                             assurance = NULL) {

  if (!is.null(p)) {
    if (p == 1) {
      if (!is.null(rho_XX)) rho_XX <- as.matrix(rho_XX)
      if (!is.null(rho_YX)) rho_YX <- as.matrix(rho_YX)
      if (!is.null(rho2_j_X_without_j)) {
        if (rho2_j_X_without_j != 0) stop("If p=1, how could 'rho2_j_X_without_j' be nonzero?")
      }
      if (!is.null(rho_XX)) {
        if (dim(rho_XX) != c(1, 1)) stop("If p=1, 'rho_XX' should be a 1 x 1 matrix?")
        if (dim(rho_YX) != c(1, 1)) stop("If p=1, how can 'rho_XX' not be a 1 x 1 matrix?")
      }
    }
  }
  # Thanks to Jan Herman for modificatoin of the function just above here to allow p=1 to work properly.

  Expected_R2 <- function(population_R2, N, p) {
    # gsl::hyperg_2F1(1, 1, 0.5 * (N + 1), population_R2) is another way to
    # obtain this 2F1 value; computed in base R via .hyperg_2F1() (no GSL
    # system dependency; see R/R2_internals.R).
    Value <- 1 - ((N - p - 1) / (N - 1)) * (1 - population_R2) *
      .hyperg_2F1(1, 1, 0.5 * (N + 1), population_R2)
    Value <- max(0, Value)
    return(Value)
  }
  if (noncentral == TRUE && is.null(sigma_X)) {
    sigma_X <- 1
  }
  if (noncentral == TRUE && is.null(sigma_Y)) {
    sigma_Y <- 1
  }
  if (noncentral == TRUE && sigma_Y != 1) {
    stop("Since you've specified 'noncentral=TRUE', all variances should be one (your 'sigma_Y' (i.e., the Y standard deviation) is not one), for a standardized solution.")
  }
  if (noncentral == TRUE && sigma_X != 1) {
    stop("Since you've specified 'noncentral=TRUE', all variances should be one (your 'sigma_X' (i.e., the X standard deviation) is not one), for a standardized solution.")
  }
  if (is.null(p) && is.null(rho_XX)) {
    stop("Since rho_XX is not specified, you must specify 'p'.")
  }
  if (!is.null(rho_XX)) {
    if (!(sum(round(rho_XX, 5) == round(t(rho_XX), 5)) ==
      dim(as.matrix(rho_XX))[1] * dim(as.matrix(rho_XX))[2])) {
      stop("The correlation matrix, 'rho_XX' should be symmetric.")
    }
  }
  if (is.null(p) && !is.null(rho_XX)) {
    p <- dim(rho_XX)[1]
  }
  char.expand(which_width, c("Full", "Lower", "Upper"), nomatch = stop("Problems with 'which_width' specification. You must choose either 'Full', 'Lower', or 'Upper'.",
    call. = FALSE
  ))
  if (which_width == "Lower" || which_width == "Upper") {
    stop("At the present time, only the 'which_width' of 'Full' is implemented.",
      call. = FALSE
    )
  }
  if (is.null(conf_level)) {
    if (!is.numeric(alpha_lower) || !is.numeric(alpha_upper)) {
      stop("Since 'conf_level' is not specified, you need to correctly specify 'alpha_lower' and 'alpha_upper'.")
    }
    if (alpha_lower < 0 || alpha_lower >= 1 || alpha_upper <
      0 | alpha_upper >= 1) {
      stop("You have not correctly specified 'alpha_lower' and/or 'alpha_upper'.")
    }
  }
  if (!is.null(conf_level)) {
    if (!is.null(alpha_lower) || !is.null(alpha_upper)) {
      stop("Since 'conf_level' is specified, do not specify 'alpha_lower' and 'alpha_upper'.")
    }
    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  }
  if (!is.null(assurance)) {
    if ((assurance <= 0) || (assurance >=
      1)) {
      stop("The 'assurance' must either be NULL or some value greater than .50 and less than 1.",
        call. = FALSE
      )
    }
    if (assurance <= 0.5) {
      stop("The 'assurance' should be > .5 (but less than 1).",
        call. = FALSE
      )
    }
  }
  if ((!is.null(rho2_j_X_without_j) && !is.null(rho2_Y_X)) &
    (!is.null(rho_XX) & !is.null(rho_YX))) {
    stop("Since 'rho2_j_X_without_j' and 'rho2_Y_X' are specified, do not specify 'rho_XX' or 'rho_YX' (or vice versa).",
      call. = FALSE
    )
  }
  if (!is.null(rho_XX) && !is.null(rho_YX)) {
    if (!is.null(rho2_j_X_without_j) || !is.null(rho2_Y_X)) {
      stop("Since 'rho_XX' and 'rho_YX' are specified, do not specify 'rho2_Y_X' or 'rho2_j_X_without_j'.",
        call. = FALSE
      )
    }
    rho2_Y_X <- (rho_YX %*% solve(rho_XX) %*% rho_YX)
    rho2_j_X_without_j <- 1 - 1 / solve(rho_XX)[
      which_predictor,
      which_predictor
    ]
    # if-statement added by Jan Herman to cope with the p=1 case
    if (p > 1) {
      rho2_Y_X_without_j <- (rho_YX[-which_predictor] %*% solve(rho_XX[
        -which_predictor,
        -which_predictor
      ]) %*% rho_YX[-which_predictor])
      b_j_tmp <- (solve(rho_XX) %*% rho_YX)[which_predictor]
    } else {
      rho2_Y_X_without_j <- as.matrix(0) # If p=1 then rho2_Y_X_without_j will be zero
    }
    if (!is.null(b_j)) {
      if (!isTRUE(all.equal(as.numeric(b_j), as.numeric(b_j_tmp)))) {
        stop("The covariance structure implied regression coefficient and 'b_j' are not equal; this is a problem.")
      }
    }
  }
  if (!is.null(rho2_j_X_without_j) && !is.null(rho2_Y_X) && (is.null(rho_XX) |
    is.null(rho_YX))) {
    if (is.null(b_j)) {
      stop("Since 'rho_XX' and 'rho_YX' are not specified, implying 'rho2_j_X_without_j' and 'rho2_Y_X' are specified, 'b_j' must also be specified.",
        call. = FALSE
      )
    }
  }
  n0 <- (qnorm(1 - (alpha_lower + alpha_upper) / 2) / (width *
    0.5))^2 * (1 - rho2_Y_X) / (1 - rho2_j_X_without_j) * (sigma_Y^2 / sigma_X^2) +
    p + 1
  n1 <- max(ceiling(n0) - 10, 2 * p)
  Diff <- 1
  while (Diff > 0) {
    n1 <- n1 + 1
    E_rho2_Y_X <- Expected_R2(
      population_R2 = rho2_Y_X, N = n1,
      p = p
    )
    E_rho2_j_X_without_j <- Expected_R2(
      population_R2 = rho2_j_X_without_j,
      N = n1, p = p
    )
    CV_i <- (qt(1 - (alpha_lower + alpha_upper) / 2, n1 - p -
      1))
    SD_i <- sqrt(((1 - E_rho2_Y_X) / ((1 - E_rho2_j_X_without_j) *
      (n1 - p - 1)))) * (sigma_Y / sigma_X)
    Current_Width <- CV_i * SD_i * 2
    Diff <- Current_Width - width
  }
  N <- n1
  if (!is.null(assurance)) {
    if ((assurance <= 0) || (assurance >=
      1)) {
      stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.",
        call. = FALSE
      )
    }
    if (assurance <= 0.5) {
      stop("The 'assurance' should be > .5 (but less than 1).",
        call. = FALSE
      )
    }
    E_rho2_Y_X <- Expected_R2(
      population_R2 = rho2_Y_X, N = N,
      p = p
    )
    E_rho2_j_X_without_j <- Expected_R2(
      population_R2 = rho2_j_X_without_j,
      N = N, p = p
    )
    N_M <- (qt(1 - (alpha_lower + alpha_upper) / 2, N - p -
      1) / (width * 0.5))^2 * ((1 - E_rho2_Y_X) / (1 - E_rho2_j_X_without_j)) *
      (sigma_Y^2 / sigma_X^2) * (qchisq(
        assurance,
        N - 1
      ) / (N - p - 1)) + p + 1
    N_M <- ceiling(N_M)
  }
  if (noncentral == FALSE && is.null(assurance)) {
    return(.as_dmar_tbl(data.frame(term = "necessary_N", value = N), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }
  if (noncentral == FALSE && !is.null(assurance)) {
    return(.as_dmar_tbl(data.frame(term = "necessary_N", value = N_M), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }
  if (noncentral == TRUE) {
    if (is.null(b_j)) {
      b_j <- (solve(rho_XX) %*% rho_YX)[which_predictor]
      if (is.null(b_j)) stop("b_j must be specified directly or obtained from other combinations of parameters.")
    }
    n2 <- max(N - 6, 2 * p + 1)
    Diff <- 1
    while (Diff > 0) {
      n2 <- n2 + 1
      E_rho2_Y_X <- Expected_R2(
        population_R2 = rho2_Y_X,
        N = n2, p = p
      )
      E_rho2_j_X_without_j <- Expected_R2(
        population_R2 = rho2_j_X_without_j,
        N = n2, p = p
      )
      CI_Result_NC <- ci_reg_coef(
        b_j = b_j, SE_b_j = NULL,
        s_Y = sigma_Y, s_X = sigma_X, N = n2, p = p,
        R2_Y_X = E_rho2_Y_X, R2_j_X_without_j = E_rho2_j_X_without_j,
        conf_level = NULL, R2_Y_X_without_j = NULL, t_value = NULL,
        alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        noncentral = TRUE
      )
      current_width <- CI_Result_NC[which(CI_Result_NC$term == "upper_limit"), 2] - CI_Result_NC[which(CI_Result_NC$term == "lower_limit"), 2]
      Diff <- current_width - width
    }
    N_NC <- n2
    if (noncentral == TRUE && is.null(assurance)) {
      return(.as_dmar_tbl(data.frame(term = "necessary_N", value = N_NC), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
    if (!is.null(assurance)) {
      E_rho2_Y_X <- Expected_R2(
        population_R2 = rho2_Y_X,
        N = N_NC, p = p
      )
      E_rho2_j_X_without_j <- Expected_R2(
        population_R2 = rho2_j_X_without_j,
        N = N_NC, p = p
      )
      SE_b_M <- sqrt((((1 - E_rho2_Y_X) * sigma_Y^2) / (((1 -
        E_rho2_j_X_without_j) * (N_NC - p - 1)) * sigma_X^2)))
      CI_R2 <- ci_R2(
        R2 = E_rho2_Y_X, df_1 = p, df_2 = N_NC -
          p - 1, conf_level = assurance, F_value = NULL,
        N = NULL, p = NULL, alpha_lower = NULL, alpha_upper = NULL
      )[1, 2]
      n3 <- N_NC
      Diff <- 1
      while (Diff > 0) {
        n3 <- n3 + 1
        E_rho2_Y_X <- Expected_R2(
          population_R2 = CI_R2,
          N = n3, p = p
        )
        E_rho2_j_X_without_j <- Expected_R2(
          population_R2 = rho2_j_X_without_j,
          N = n3, p = p
        )
        CI_Result_NC <- ci_reg_coef(
          b_j = b_j, SE_b_j = NULL,
          s_Y = sigma_Y, s_X = sigma_X, N = n3, p = p,
          R2_Y_X = E_rho2_Y_X, R2_j_X_without_j = E_rho2_j_X_without_j,
          conf_level = NULL, R2_Y_X_without_j = NULL,
          t_value = NULL, alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, noncentral = TRUE
        )
        current_width <- CI_Result_NC[which(CI_Result_NC$term == "upper_limit"), 2] - CI_Result_NC[which(CI_Result_NC$term == "lower_limit"), 2]
        Diff <- current_width - width
      }
      N_NC_M <- n3
      if (noncentral == TRUE && !is.null(assurance)) {
        return(.as_dmar_tbl(data.frame(term = "necessary_N", value = N_NC_M), conf_level = conf_level, subclass = "dmar_ss_aipe"))
      }
    }
  }
}
