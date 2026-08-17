#' Sample Size Planning for Accuracy in Parameter Estimation for the Multiple Correlation Coefficient
#'
#' @description
#' Determines necessary sample size for the multiple correlation coefficient so that the confidence
#' interval for the population multiple correlation coefficient is sufficiently narrow. Optionally,
#' there is a certainty parameter that allows one to be a specified percent certain that the observed
#' interval will be no wider than desired.
#'
#' @param population_R2 Value of the population multiple correlation coefficient
#' @param conf_level Confidence interval level (e.g., .95, .99, .90); 1-Type I error rate
#' @param width Width of the confidence interval (see \code{which_width})
#' @param random_predictors Whether or not the predictor variables are random (set to \code{TRUE}) or are fixed (set to \code{FALSE})
#' @param which_width Defines the width that \code{width} refers to
#' @param p The number of predictor variables
#' @param assurance Value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g, .80, .90, .95)
#' @param verify_ss Evaluates numerically via an internal Monte Carlo simulation the exact sample size given the specifications
#' @param tol The tolerance of the iterative function \code{conf_limits_nct} for convergence
#' @param \dots For modifying the parameters of functions this function calls upon
#'
#' @details
#' This function determines a necessary sample size so that the expected confidence interval width for the
#' squared multiple correlation coefficient is sufficiently narrow (when \code{assurance=NULL})
#' so that the obtained confidence interval is no larger than the value specified with some desired degree
#' of certainty (i.e., a probability that the obtained width is less than the specified width). The method
#' depends on whether or not the regressors are regarded as fixed or random. This is the case because the
#' distribution theory for the two cases is different and thus the confidence interval procedure is
#' conditional on the type of regressors. The default methods are approximate but can be made exact with the
#' specification of \code{verify_ss=TRUE}, which performs an a priori Monte Carlo simulation study.
#' Kelley (2008) and Kelley & Maxwell (2008) detail the methods used in the function, with the former focusing
#' on random regressors and the latter on fixed regressors.
#'
#' It is recommended that the option \code{verify_ss} should always be used! Doing so uses the method implied
#' sample size as an estimate and then evaluates with an internal Monte Carlo simulation
#' (i.e., via "brute-force" methods) the exact sample size given the goals specified. When \code{verify_ss=TRUE},
#' the default number of iterations is 10,000 but this can be changed by specifying G=5000 (or some other value;
#' 10000 is the recommended) When \code{verify_ss=TRUE} is specified, an internal function \code{verify_ss_aipe_r2}
#' calls upon the \code{ss_aipe_R2_sensitivity} function for purposes of the internal Monte Carlo simulation
#' study. See the \code{verify_ss_aipe_r2} function for arguments that can be passed from \code{ss_aipe_R2}
#' to \code{verify_ss_aipe_r2}.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} value is \code{"necessary_N"} and \code{value} is the
#' necessary total sample size \emph{N} given the input specifications.
#'
#' @references
#' Algina, J. & Olejnik, S. (2000). Determining sample size for accurate estimation of the squared
#' multiple correlation coefficient. \emph{Multivariate Behavioral Research, 35}, 119--137.
#'   \doi{10.1207/s15327906mbr3501_5}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43},
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
#'   interval estimation, power calculations, sample size estimation, and
#'   hypothesis testing in multiple regression. \emph{Behavior Research
#'   Methods, Instruments, & Computers, 24}(4), 581--582.
#'   \doi{10.3758/BF03203611}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' With \code{verify_ss = TRUE} the function can take some time to converge
#' (e.g., several minutes to a quarter hour) because the closed form
#' approximation is followed by an a priori Monte Carlo simulation. The
#' default \code{verify_ss = FALSE} returns the closed form approximation
#' only and is essentially instantaneous.
#'
#' @seealso
#' \code{\link{ci_R2}}, \code{\link{conf_limits_nct}}, \code{\link{ss_aipe_R2_sensitivity}}
#'
#' @examples
#' # 1. Closed form planner under random predictors (the typical case).
#' #    Sample size sufficient for the expected CI width on rho^2 to be .10.
#' ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
#'            which_width = "Full", p = 5, random_predictors = TRUE)
#'
#' # 2. The same target under fixed predictors (planned dosing levels,
#' #    factorial covariates, and the like) needs a smaller N, and adding an
#' #    assurance of .85, so that the realized width is no larger than the
#' #    target in 85 percent of replications rather than only on average,
#' #    needs a larger one. Each is another pass of the same iterative search
#' #    over N, so the two calls are shown here rather than run:
#' #
#' #    ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
#' #               which_width = "Full", p = 5, random_predictors = FALSE)
#' #
#' #    ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
#' #               which_width = "Full", p = 5, assurance = .85,
#' #               random_predictors = TRUE)
#'
#' # 3. verify_ss = TRUE follows the closed form approximation with an a
#' #    priori Monte Carlo simulation of the realized width, starting from
#' #    the closed form answer and returning the sample size the simulation
#' #    settles on, which is what a plan meant to be defended deserves. G is
#' #    the number of replications in that simulation; 10000 is the default
#' #    and the recommendation. The call runs for several minutes, so it too
#' #    is shown rather than run:
#' #
#' #    set.seed(113)
#' #    ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
#' #               which_width = "Full", p = 5, random_predictors = TRUE,
#' #               verify_ss = TRUE, G = 10000)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_R2 <- function(population_R2 = NULL, conf_level = 0.95, width = NULL, random_predictors = TRUE,
                       which_width = "Full", p = NULL, assurance = NULL,
                       verify_ss = FALSE, tol = 1e-09, ...) {
  # The iterative search calls ci_R2() (and through it conf_limits_ncf())
  # repeatedly with intermediate trial values of N. For early values of N the
  # noncentral F lower limit is often clamped to 0, which conf_limits_ncf
  # signals via a warning. Surfacing that warning once per iteration produces
  # dozens of identical messages. We muffle the per-iteration warnings with
  # withCallingHandlers, count them, and emit a single summary warning at the
  # end (via on.exit so it fires on any return path).
  .clamp_count <- 0L
  on.exit({
    if (.clamp_count > 0L) {
      warning(sprintf(
        "During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in %d intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_ncf for the meaning of the clamp.",
        .clamp_count
      ), call. = FALSE)
    }
  }, add = TRUE)

  withCallingHandlers({
  verify_ss_aipe_r2 <- function(population_R2 = NULL, conf_level = .95, width = NULL, random_predictors = TRUE,
                                which_width = "Full", p = NULL, n = NULL, assurance = NULL, g = 500, G = 10000,
                                print_iter = FALSE, ...) {
    # Based on method called for Sample size.
    Summary_of_SS <- ss_aipe_R2_sensitivity(
      specified_N = n,
      true_R2 = population_R2, estimated_R2 = NULL, w = width, p = p,
      random_predictors = random_predictors, assurance = assurance,
      conf_level = conf_level, generate_random_predictors = random_predictors,
      rho_yx = .3, rho_xx = .3, G = g, print_iter = print_iter
    )

    i <- 1

    if (is.null(assurance)) {
      E_Width_Info <- as.matrix(cbind(E.Width.CI = Summary_of_SS[which(Summary_of_SS$term == "mean_ci_width"), 2], N = n))
      w_Bigger <- E_Width_Info[1, 1] > width
      w_Smaller <- E_Width_Info[1, 1] < width

      while (i != 0) {
        i <- i + 1

        if (w_Bigger == TRUE) n <- n + 1
        if (w_Smaller == TRUE) n <- n - 1

        Summary_of_SS <- ss_aipe_R2_sensitivity(
          specified_N = n, true_R2 = population_R2,
          estimated_R2 = NULL, w = width, p = p, random_predictors = random_predictors,
          assurance = NULL, conf_level = conf_level,
          generate_random_predictors = random_predictors,
          rho_yx = .3, rho_xx = .3, G = g, print_iter = print_iter
        )

        E_Width_Info <- rbind(E_Width_Info, c(E.Width.CI = Summary_of_SS[which(Summary_of_SS$term == "mean_ci_width"), 2], N = n))

        if (w_Bigger == FALSE) w_Bigger <- E_Width_Info[i, 1] > width
        if (w_Smaller == FALSE) w_Smaller <- E_Width_Info[i, 1] < width

        if ((w_Bigger == TRUE) && (w_Smaller == TRUE)) break()
      }

      # Now get exact sample size (using previous results).
      n <- E_Width_Info[nrow(E_Width_Info), 2]

      # Based on method called for Sample size.
      Summary_of_SS <- ss_aipe_R2_sensitivity(
        specified_N = n,
        true_R2 = population_R2, estimated_R2 = NULL, w = width, p = p,
        random_predictors = random_predictors, assurance = assurance,
        conf_level = conf_level, generate_random_predictors = random_predictors,
        rho_yx = .3, rho_xx = .3, G = G, print_iter = print_iter
      )

      E_Width_Info <- as.matrix(cbind(E.Width.CI = Summary_of_SS[which(Summary_of_SS$term == "mean_ci_width"), 2], N = n))

      i <- 1

      w_Bigger <- E_Width_Info[1, 1] > width
      w_Smaller <- E_Width_Info[1, 1] < width

      while (i != 0) {
        i <- i + 1

        if (w_Bigger == TRUE) n <- n + 1
        if (w_Smaller == TRUE) n <- n - 1

        Summary_of_SS <- ss_aipe_R2_sensitivity(
          specified_N = n, true_R2 = population_R2,
          estimated_R2 = NULL, w = width, p = p, random_predictors = random_predictors,
          assurance = NULL, conf_level = conf_level,
          generate_random_predictors = random_predictors,
          rho_yx = .3, rho_xx = .3, G = G, print_iter = print_iter
        )

        E_Width_Info <- rbind(E_Width_Info, c(E.Width.CI = Summary_of_SS[which(Summary_of_SS$term == "mean_ci_width"), 2], N = n))

        if (w_Bigger == FALSE) w_Bigger <- E_Width_Info[i, 1] > width
        if (w_Smaller == FALSE) w_Smaller <- E_Width_Info[i, 1] < width

        if ((w_Bigger == TRUE) && (w_Smaller == TRUE)) break()
      }

      Contending <- (E_Width_Info[, 1] <= width)
      To_Use <- min(E_Width_Info[Contending, 2])
      Here <- E_Width_Info[which(E_Width_Info[, 2] == To_Use), ]

      # Result_Full <- list(Required.Sample.Size = as.numeric(Here[2]), Expected.Width = as.numeric(Here[1]))
      Result_Full <- as.numeric(Here[2])

      return(Result_Full)
    }
    ############################################################################
    ############################################################################
    # The following is used if there is a desired degree of certainty specified.
    ############################################################################
    ############################################################################

    if (!is.null(assurance)) {
      E_Gamma_Info <- as.matrix(cbind(E.Gamma = Summary_of_SS[which(Summary_of_SS$term == "pct_ci_less_w"), 2], N = n))
      g_Bigger <- E_Gamma_Info[1, 1] > assurance
      g_Smaller <- E_Gamma_Info[1, 1] < assurance


      while (i != 0) {
        i <- i + 1

        if (g_Bigger == TRUE) n <- n - 1
        if (g_Smaller == TRUE) n <- n + 1

        Summary_of_SS <- ss_aipe_R2_sensitivity(
          specified_N = n, true_R2 = population_R2,
          estimated_R2 = NULL, w = width, p = p, random_predictors = random_predictors,
          assurance = NULL, conf_level = conf_level,
          generate_random_predictors = random_predictors,
          rho_yx = .3, rho_xx = .3, G = g, print_iter = print_iter
        )

        E_Gamma_Info <- rbind(E_Gamma_Info, c(E.Gamma = Summary_of_SS[which(Summary_of_SS$term == "pct_ci_less_w"), 2], N = n))

        if (g_Bigger == FALSE) g_Bigger <- E_Gamma_Info[i, 1] > assurance
        if (g_Smaller == FALSE) g_Smaller <- E_Gamma_Info[i, 1] < assurance

        if ((g_Bigger == TRUE) && (g_Smaller == TRUE)) break()
      }

      # Now get exact sample size (using previous results).
      n <- E_Gamma_Info[nrow(E_Gamma_Info), 2]

      # Based on method called for Sample size.
      Summary_of_SS <- ss_aipe_R2_sensitivity(
        specified_N = n,
        true_R2 = population_R2, estimated_R2 = NULL, w = width, p = p,
        random_predictors = random_predictors, assurance = assurance,
        conf_level = conf_level, generate_random_predictors = random_predictors,
        rho_yx = .3, rho_xx = .3, G = G, print_iter = print_iter
      )

      E_Gamma_Info <- as.matrix(cbind(E.Gamma = Summary_of_SS[which(Summary_of_SS$term == "pct_ci_less_w"), 2], N = n))

      i <- 1

      g_Bigger <- E_Gamma_Info[1, 1] > assurance
      g_Smaller <- E_Gamma_Info[1, 1] < assurance

      while (i != 0) {
        i <- i + 1

        if (g_Bigger == TRUE) n <- n - 1
        if (g_Smaller == TRUE) n <- n + 1

        Summary_of_SS <- ss_aipe_R2_sensitivity(
          specified_N = n, true_R2 = population_R2,
          estimated_R2 = NULL, w = width, p = p, random_predictors = random_predictors,
          assurance = NULL, conf_level = conf_level,
          generate_random_predictors = random_predictors,
          rho_yx = .3, rho_xx = .3, G = G, print_iter = print_iter
        )

        E_Gamma_Info <- rbind(E_Gamma_Info, c(E.Gamma = Summary_of_SS[which(Summary_of_SS$term == "pct_ci_less_w"), 2], N = n))

        if (g_Bigger == FALSE) g_Bigger <- E_Gamma_Info[i, 1] > assurance
        if (g_Smaller == FALSE) g_Smaller <- E_Gamma_Info[i, 1] < assurance

        if ((g_Bigger == TRUE) && (g_Smaller == TRUE)) break()
      }

      Contending <- (E_Gamma_Info[, 1] >= assurance)
      To_Use <- min(E_Gamma_Info[Contending, 2])
      Here <- E_Gamma_Info[which(E_Gamma_Info[, 2] == To_Use), ]

      Result_Full <- as.numeric(Here[2])
      return(Result_Full)
    }
  }


  ################# ss_aipe_R2 #################

  char.expand(which_width, c("Full", "Lower", "Upper"), nomatch = stop("Problems with 'which_width' specification. You must choose either 'Full', 'Lower', or 'Upper'.",
    call. = FALSE
  ))
  if (is.null(p)) {
    stop("You need to specify 'p', the number of predictors.")
  }
  if (is.null(population_R2)) {
    stop("You need to specify the population squared multiple correlation coefficient, 'population_R2'.")
  }
  if (!is.null(assurance)) {
    if (assurance <= 0.49 || assurance > 1) {
      stop("The 'assurance' must be between .50 and 1.")
    }
  }
  Expected_R2 <- function(population_R2, N, p) {
    # gsl::hyperg_2F1(1, 1, 0.5 * (N + 1), population_R2) is another way to
    # obtain this 2F1 value; computed in base R via .hyperg_2F1() (no GSL
    # system dependency; see R/R2_internals.R).
    Value <- 1 - ((N - p - 1) / (N - 1)) * (1 - population_R2) *
      .hyperg_2F1(1, 1, 0.5 * (N + 1), population_R2)
    Value <- max(0, Value)
    return(Value)
  }
  To_Find_Pop_R2_Given_Expectation <- function(E_R2_low = conf_limit_desired_certainty_lower,
                                               E_R2_up = conf_limit_desired_certainty_upper, N = N,
                                               p = p) {
    True_vals <- seq(0.001, 0.999, 0.001)
    Exp_vals <- rep(NA, length(True_vals))
    for (i in 1:length(True_vals)) {
      Exp_vals[i] <- Expected_R2(True_vals[i], N, p)
    }
    match(round(E_R2_low, 2), round(Exp_vals, 3))
    return(c(mean(True_vals[match(round(E_R2_low, 2), round(
      Exp_vals,
      3
    ))]), mean(True_vals[match(round(E_R2_up, 2), round(
      Exp_vals,
      3
    ))])))
  }
  alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  N_0 <- p + 1 + p
  Continue <- TRUE
  while (Continue == TRUE) {
    CI_0 <- ci_R2(
      R2 = Expected_R2(
        population_R2 = population_R2,
        N = N_0, p = p
      ), conf_level = NULL, alpha_lower = alpha_lower,
      alpha_upper = alpha_upper, N = N_0, p = p, random_predictors = FALSE
    )
    Continue <- sum(is.na(CI_0$value)) > 0
    N_0 <- N_0 + 1
  }
  if (which_width == "Full") {
    N_1 <- N_0 + 1
    # Inner-loop fast path: skip ci_R2's data.frame wrapping; the
    # fixed-predictors helper returns just c(lower, upper).
    R2_hat <- Expected_R2(population_R2 = population_R2,
                          N = N_1, p = p)
    lims <- .ci_R2_fixed_limits_fast(R2_hat, N_1, p, alpha_lower,
                                     alpha_upper, tol = tol)
    w_F <- lims[2L] - lims[1L]
    Diff <- w_F - width
    while (Diff > tol || is.na(Diff)) {
      N_1 <- N_1 + 1
      R2_hat <- Expected_R2(population_R2 = population_R2,
                            N = N_1, p = p)
      lims <- .ci_R2_fixed_limits_fast(R2_hat, N_1, p, alpha_lower,
                                       alpha_upper, tol = tol)
      w_F <- lims[2L] - lims[1L]
      Diff <- w_F - width
    }
    if (random_predictors == TRUE) {
      # Inner-loop fast path: call .ci_R2_random_limits_fast which
      # returns just c(lower, upper) and skips the per-iteration
      # data.frame and class-tagging overhead of ci_R2(). Falls
      # back to ci_R2() if the fast path detects a degenerate
      # one-sided case.
      R2_hat <- Expected_R2(population_R2 = population_R2,
                            N = N_1, p = p)
      lims <- .ci_R2_random_limits_fast(R2_hat, N_1, p,
                                        alpha_lower, alpha_upper)
      if (anyNA(lims)) {
        # Defer to the canonical ci_R2() path for edge cases.
        CI_1 <- ci_R2(
          R2 = R2_hat, alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, N = N_1, p = p,
          random_predictors = TRUE
        )
        lims <- c(CI_1[1, 2], CI_1[3, 2])
      }
      w_F <- lims[2L] - lims[1L]
      Diff <- w_F - width
      while (Diff > tol || is.na(Diff)) {
        N_1 <- N_1 + 1
        R2_hat <- Expected_R2(population_R2 = population_R2,
                              N = N_1, p = p)
        lims <- .ci_R2_random_limits_fast(R2_hat, N_1, p,
                                          alpha_lower, alpha_upper)
        if (anyNA(lims)) {
          CI_1 <- ci_R2(
            R2 = R2_hat, alpha_lower = alpha_lower,
            alpha_upper = alpha_upper, N = N_1, p = p,
            random_predictors = TRUE
          )
          lims <- c(CI_1[1, 2], CI_1[3, 2])
        }
        w_F <- lims[2L] - lims[1L]
        Diff <- w_F - width
      }
    }
    Result_Full <- data.frame(term = "necessary_N", value = N_1)
    if (verify_ss == FALSE) {
      if (is.null(assurance) == TRUE) {
        return(.as_dmar_tbl(Result_Full, conf_level = conf_level, subclass = "dmar_ss_aipe"))
      }
    }
    if (verify_ss == TRUE) {
      Result_Full <- data.frame(term = "necessary_N", value = verify_ss_aipe_r2(
        population_R2 = population_R2,
        conf_level = conf_level, width = width, random_predictors = random_predictors,
        which_width = "Full", p = p, n = Result_Full[1, 2],
        assurance = assurance, ...
      ))
      if (is.null(assurance) == TRUE) {
        return(.as_dmar_tbl(Result_Full, conf_level = conf_level, subclass = "dmar_ss_aipe"))
      }
    }
    if (!is.null(assurance)) {
      conf_limit_desired_certainty_lower <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ), alpha_lower = 1 - assurance,
        alpha_upper = 0, N = N_1, p = p, random_predictors = random_predictors
      )[1, 2]
      conf_limit_desired_certainty_upper <- ci_R2(R2 = Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      ), alpha_lower = 0, alpha_upper = 1 -
        assurance, N = N_1, p = p, random_predictors = random_predictors)[3, 2]
      N_2 <- N_1
      CI_2 <- ci_R2(
        R2 = conf_limit_desired_certainty_upper,
        alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_2, p = p, random_predictors = random_predictors
      )
      w_F <- CI_2[3, 2] - CI_2[1, 2]
      Diff <- w_F - width
      while (Diff > tol || is.na(Diff)) {
        N_2 <- N_2 + 1
        conf_limit_desired_certainty_upper_i <- ci_R2(R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_2, p = p
        ), alpha_lower = 0, alpha_upper = 1 -
          assurance, N = N_2, p = p, random_predictors = random_predictors)[3, 2]
        CI_2 <- ci_R2(
          R2 = conf_limit_desired_certainty_upper_i,
          alpha_lower = alpha_lower, alpha_upper = alpha_upper,
          N = N_2, p = p, random_predictors = random_predictors
        )
        w_F <- CI_2[3, 2] - CI_2[1, 2]
        Diff <- w_F - width
      }
      N_Upper_Conf_Lim <- N_2
      Ex_Width_Upper <- w_F
      N_3 <- N_1
      CI_3 <- ci_R2(
        R2 = conf_limit_desired_certainty_lower,
        alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_3, p = p, random_predictors = random_predictors
      )
      w_F <- CI_3[3, 2] - CI_3[1, 2]
      Diff <- w_F - width
      while (Diff > tol || is.na(Diff)) {
        N_3 <- N_3 + 1
        conf_limit_desired_certainty_lower_i <- ci_R2(
          R2 = Expected_R2(
            population_R2 = population_R2,
            N = N_3, p = p
          ), alpha_lower = 1 - assurance,
          alpha_upper = 0, N = N_3, p = p, random_predictors = random_predictors
        )[1, 2]
        CI_3 <- ci_R2(
          R2 = conf_limit_desired_certainty_lower_i,
          alpha_lower = alpha_lower, alpha_upper = alpha_upper,
          N = N_3, p = p, random_predictors = random_predictors
        )
        w_F <- CI_3[3, 2] - CI_3[1, 2]
        Diff <- w_F - width
      }
      N_Lower_Conf_Lim <- N_3
      Ex_Width_Lower <- w_F
      if (N_Upper_Conf_Lim >= N_Lower_Conf_Lim) {
        Result_Full <- data.frame(term = "necessary_N", value = N_Upper_Conf_Lim)
      }
      if (N_Lower_Conf_Lim > N_Upper_Conf_Lim) {
        Result_Full <- data.frame(term = "necessary_N", value = N_Lower_Conf_Lim)
      }
      CI_WIDTH_R2 <- function(R2, conf_level, N, p, random_predictors) {
        Lims <- ci_R2(
          R2 = Expected_R2(
            population_R2 = R2,
            N, p
          ), conf_level = conf_level, N = N, p = p,
          random_predictors = random_predictors
        )
        Lims[3, 2] - Lims[1, 2]
      }
      Pop_Values <- To_Find_Pop_R2_Given_Expectation(
        E_R2_low = conf_limit_desired_certainty_lower,
        E_R2_up = conf_limit_desired_certainty_upper,
        N = N_1, p = p
      )
      optimize_result <- optimize(
        f = CI_WIDTH_R2, interval = c(max(Pop_Values[1] *
          0.98, 1e-06), min(Pop_Values[2] * 1.02, 0.999999)),
        maximum = TRUE, tol = .Machine$double.eps^0.5,
        N = N_1, p = p, conf_level = conf_level, random_predictors = random_predictors
      )$maximum
      optimize_result <- Expected_R2(
        population_R2 = optimize_result,
        N_1, p
      )
      N_Op <- NULL
      N_Op_Up <- NULL
      N_Op_Low <- NULL
      if ((round(conf_limit_desired_certainty_lower, 4) <
        round(optimize_result, 4)) & (round(
        optimize_result,
        4
      ) < round(
        conf_limit_desired_certainty_upper,
        4
      ))) {
        if (random_predictors == TRUE) {
          conf_limit_desired_certainty_lower <- ci_R2(
            R2 = optimize_result,
            alpha_lower = 1 - assurance, alpha_upper = 0,
            N = N_1, p = p, random_predictors = random_predictors
          )[1, 2]
          conf_limit_desired_certainty_upper <- ci_R2(
            R2 = optimize_result,
            alpha_lower = 0, alpha_upper = 1 - assurance,
            N = N_1, p = p, random_predictors = random_predictors
          )[3, 2]
          if (conf_limit_desired_certainty_lower > optimize_result) {
            conf_limit_desired_certainty_lower <- optimize_result
          }
          N_Op_Low <- N_1
          CI_Op_Low <- ci_R2(
            R2 = conf_limit_desired_certainty_lower,
            alpha_lower = alpha_lower, alpha_upper = alpha_upper,
            N = N_Op_Low, p = p, random_predictors = random_predictors
          )
          w_F <- CI_Op_Low[3, 2] - CI_Op_Low[1, 2]
          Diff <- w_F - width
          while (Diff > tol || is.na(Diff)) {
            N_Op_Low <- N_Op_Low + 1
            CI_Op_Low <- ci_R2(
              R2 = conf_limit_desired_certainty_lower,
              alpha_lower = alpha_lower, alpha_upper = alpha_upper,
              N = N_Op_Low, p = p, random_predictors = random_predictors
            )
            w_F <- CI_Op_Low[3, 2] - CI_Op_Low[1, 2]
            Diff <- w_F - width
          }
          w_Op_Width_Low <- w_F
          if (conf_limit_desired_certainty_upper < optimize_result) {
            conf_limit_desired_certainty_upper <- optimize_result
          }
          N_Op_Up <- N_1
          CI_Op_Up <- ci_R2(
            R2 = conf_limit_desired_certainty_upper,
            alpha_lower = alpha_lower, alpha_upper = alpha_upper,
            N = N_Op_Up, p = p, random_predictors = random_predictors
          )
          w_F <- CI_Op_Low[3, 2] - CI_Op_Low[1, 2]
          Diff <- w_F - width
          while (Diff > tol || is.na(Diff)) {
            N_Op_Up <- N_Op_Up + 1
            CI_Op_Up <- ci_R2(
              R2 = conf_limit_desired_certainty_upper,
              alpha_lower = alpha_lower, alpha_upper = alpha_upper,
              N = N_Op_Up, p = p, random_predictors = random_predictors
            )
            w_F <- CI_Op_Up[3, 2] - CI_Op_Up[1, 2]
            Diff <- w_F - width
          }
          w_Op_Width_Up <- w_F
          if (N_1 == N_Lower_Conf_Lim && N_1 == N_Upper_Conf_Lim &
            N_1 == N_Op_Up & N_1 == N_Op_Low) {
            N_Op <- N_1
            CI_Op <- ci_R2(
              R2 = optimize_result, alpha_lower = alpha_lower,
              alpha_upper = alpha_upper, N = N_Op, p = p,
              random_predictors = random_predictors
            )
            w_F <- CI_Op[3, 2] - CI_Op[1, 2]
            Diff <- w_F - width
            while (Diff > tol || is.na(Diff)) {
              N_Op <- N_Op + 1
              CI_Op <- ci_R2(
                R2 = optimize_result, alpha_lower = alpha_lower,
                alpha_upper = alpha_upper, N = N_Op,
                p = p, random_predictors = random_predictors
              )
              w_F <- CI_Op[3, 2] - CI_Op[1, 2]
              Diff <- w_F - width
            }
            w_Op <- w_F
          }
          Necessary_N <- max(N_1, N_Lower_Conf_Lim, N_Upper_Conf_Lim,
            N_Op_Low, N_Op_Up, N_Op,
            na.rm = TRUE
          )
          For_Ex_Width <- ci_R2(
            R2 = Expected_R2(
              population_R2 = population_R2,
              N = Necessary_N, p = p
            ), alpha_lower = alpha_lower,
            alpha_upper = alpha_upper, N = Necessary_N,
            p = p, random_predictors = random_predictors
          )
          Result_Full <- data.frame(term = "necessary_N", value = Necessary_N)
        }
        if (random_predictors == FALSE) {
          Prob_F_Max_Width <- pf(q = convert_R2_f(
            R2 = optimize_result,
            p = p, N = N_1
          )[1, 2], df1 = p, df2 = N_1 - p -
            1, ncp = convert_R2_lambda(
            R2 = population_R2,
            N = N_1
          )[1, 2])
          low_lim_F <- qf(p = max(0, (Prob_F_Max_Width -
            (1 - assurance) / 2)), df1 = p, df2 = N_1 -
            p - 1, ncp = convert_R2_lambda(
            R2 = population_R2,
            N = N_1
          )[1, 2])
          up_lim_F <- qf(
            p = min((Prob_F_Max_Width +
              (1 - assurance) / 2), 1), df1 = p,
            df2 = N_1 - p - 1, ncp = convert_R2_lambda(
              R2 = population_R2,
              N = N_1
            )[1, 2]
          )
          low_R2 <- convert_f_R2(
            F_value = low_lim_F, df_1 = p,
            df_2 = N_1 - p - 1
          )[1, 2]
          up_R2 <- convert_f_R2(
            F_value = up_lim_F, df_1 = p,
            df_2 = N_1 - p - 1
          )[1, 2]
          n_1 <- N_1 - 1
          if (up_R2 > 0 && up_R2 < 1) {
            N_4 <- n_1 + 1
            CI_4 <- ci_R2(
              R2 = Expected_R2(
                population_R2 = up_R2,
                N = N_4, p = p
              ), alpha_lower = alpha_lower,
              alpha_upper = alpha_upper, N = N_4, p = p,
              random_predictors = random_predictors
            )
            w_F <- CI_4[3, 2] - CI_4[1, 2]
            Diff <- w_F - width
            while (Diff > tol || is.na(Diff)) {
              N_4 <- N_4 + 1
              CI_4 <- ci_R2(
                R2 = Expected_R2(
                  population_R2 = up_R2,
                  N = N_4, p = p
                ), alpha_lower = alpha_lower,
                alpha_upper = alpha_upper, N = N_4, p = p,
                random_predictors = random_predictors
              )
              w_F <- CI_4[3, 2] - CI_4[1, 2]
              Diff <- w_F - width
            }
            n_up_R2 <- N_4
            Ex_Width_Upper_after_optim <- w_F
          }
          if (low_R2 > 0 && low_R2 < 1) {
            N_5 <- n_1 + 1
            CI_5 <- ci_R2(
              R2 = Expected_R2(
                population_R2 = low_R2,
                N = N_5, p = p
              ), alpha_lower = alpha_lower,
              alpha_upper = alpha_upper, N = N_5, p = p,
              random_predictors = random_predictors
            )
            w_F <- CI_5[3, 2] - CI_5[1, 2]
            Diff <- w_F - width
            while (Diff > tol || is.na(Diff)) {
              N_5 <- N_5 + 1
              CI_5 <- ci_R2(
                R2 = Expected_R2(
                  population_R2 = low_R2,
                  N = N_5, p = p
                ), alpha_lower = alpha_lower,
                alpha_upper = alpha_upper, N = N_5, p = p,
                random_predictors = random_predictors
              )
              w_F <- CI_5[3, 2] - CI_5[1, 2]
              Diff <- w_F - width
            }
            n_low_R2 <- N_5
            Ex_Width_Lower_after_optim <- w_F
          }
          if (low_R2 == 0 || low_R2 == 1) {
            n_low_R2 <- N_1
          }
          if (up_R2 == 0 || up_R2 == 1) {
            n_up_R2 <- N_1
          }
          if ((n_low_R2 >= n_up_R2) && (n_low_R2 > Result_Full[1, 2])) {
            Result_Full <- data.frame(term = "necessary_N", value = n_low_R2)
          }
          if ((n_up_R2 > n_low_R2) && (n_up_R2 > Result_Full[1, 2])) {
            Result_Full <- data.frame(term = "necessary_N", value = n_up_R2)
          }
        }
      }
      if (verify_ss == FALSE) {
        return(.as_dmar_tbl(Result_Full, conf_level = conf_level, subclass = "dmar_ss_aipe"))
      }
      if (verify_ss == TRUE) {
        Result_Full <- data.frame(term = "necessary_N", value = verify_ss_aipe_r2(
          population_R2 = population_R2,
          conf_level = conf_level, width = width, random_predictors = random_predictors,
          which_width = "Full", p = p, n = Result_Full[1, 2],
          assurance = assurance,
          ...
        ))
        return(.as_dmar_tbl(Result_Full, conf_level = conf_level, subclass = "dmar_ss_aipe"))
      }
    }
  }









  if (which_width == "Lower") {
    if (verify_ss == TRUE) {
      stop("verify_ss can only be used with 'which_width = \"Full\"'.")
    }
    N_1 <- N_0 + 1
    CI_1 <- ci_R2(
      R2 = Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
      N = N_1, p = p, random_predictors = FALSE
    )
    w_L <- Expected_R2(
      population_R2 = population_R2, N = N_1,
      p = p
    ) - CI_1[1, 2]
    Diff <- w_L - width
    while (Diff > tol || is.na(Diff)) {
      N_1 <- N_1 + 1
      CI_1 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_1, p = p, random_predictors = FALSE
      )
      w_L <- Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      ) - CI_1[1, 2]
      Diff <- w_L - width
    }
    if (random_predictors == TRUE) {
      CI_1 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_1, p = p, random_predictors = TRUE
      )
      w_L <- Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      ) - CI_1[1, 2]
      Diff <- w_L - width
      while (Diff > tol || is.na(Diff)) {
        N_1 <- N_1 + 1
        CI_1 <- ci_R2(
          R2 = Expected_R2(
            population_R2 = population_R2,
            N = N_1, p = p
          ), alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, N = N_1, p = p,
          random_predictors = TRUE
        )
        w_L <- Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ) - CI_1[1, 2]
        Diff <- w_L - width
      }
    }
    if (is.null(assurance)) {
      Result_Low <- data.frame(term = "necessary_N", value = N_1)
      return(.as_dmar_tbl(Result_Low, conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
    if (!is.null(assurance)) {
      conf_limit_desired_certainty <- ci_R2(R2 = Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      ), alpha_lower = 0, alpha_upper = 1 -
        assurance, N = N_1, p = p, random_predictors = random_predictors)[3, 2]
      N_2 <- N_1
      CI_2 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = conf_limit_desired_certainty,
          N = N_2, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_2, p = p, random_predictors = random_predictors
      )
      w_L <- conf_limit_desired_certainty - CI_2[1, 2]
      Diff <- w_L - width
      while (Diff > tol) {
        N_2 <- N_2 + 1
        CI_2 <- ci_R2(
          R2 = Expected_R2(
            population_R2 = conf_limit_desired_certainty,
            N = N_2, p = p
          ), alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, N = N_2, p = p,
          random_predictors = random_predictors
        )
        w_L <- conf_limit_desired_certainty - CI_2[1, 2]
        Diff <- w_L - width
      }
      w_L_M_2 <- Expected_R2(
        population_R2 = population_R2,
        N = N_2, p = p
      ) - ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_2, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_2, p = p, random_predictors = random_predictors
      )[1, 2]
      Result_Low_2 <- data.frame(term = "necessary_N", value = N_2)
      conf_limit_desired_certainty <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ), alpha_lower = 1 - assurance,
        alpha_upper = 0, N = N_1, p = p, random_predictors = random_predictors
      )[1, 2]
      N_3 <- N_1
      CI_2 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = conf_limit_desired_certainty,
          N = N_3, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_3, p = p, random_predictors = random_predictors
      )
      w_L <- conf_limit_desired_certainty - CI_2[1, 2]
      Diff <- w_L - width
      while (Diff > tol) {
        N_3 <- N_3 + 1
        CI_2 <- ci_R2(
          R2 = Expected_R2(
            population_R2 = conf_limit_desired_certainty,
            N = N_3, p = p
          ), alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, N = N_3, p = p,
          random_predictors = random_predictors
        )
        w_L <- conf_limit_desired_certainty - CI_2[1, 2]
        Diff <- w_L - width
      }
      w_L_M_3 <- Expected_R2(
        population_R2 = population_R2,
        N = N_3, p = p
      ) - ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_3, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_3, p = p, random_predictors = random_predictors
      )[1, 2]
      Result_Low_3 <- data.frame(term = "necessary_N", value = N_3)
      if (N_2 >= N_3) {
        Result_Low <- Result_Low_2
      }
      if (N_3 > N_2) {
        Result_Low <- Result_Low_3
      }
      return(.as_dmar_tbl(Result_Low, conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
  }
  if (which_width == "Upper") {
    if (verify_ss == TRUE) {
      stop("verify_ss can only be used with 'which_width = \"Full\"'.")
    }
    N_1 <- N_0 + 1
    CI_1 <- ci_R2(
      R2 = Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
      N = N_1, p = p, random_predictors = FALSE
    )
    w_U <- CI_1[3, 2] - Expected_R2(
      population_R2 = population_R2,
      N = N_1, p = p
    )
    Diff <- w_U - width
    while (Diff > tol || is.na(Diff)) {
      N_1 <- N_1 + 1
      CI_1 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_1, p = p, random_predictors = FALSE
      )
      w_U <- CI_1[3, 2] - Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      )
      Diff <- w_U - width
    }
    if (random_predictors == TRUE) {
      CI_1 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_1, p = p, random_predictors = TRUE
      )
      w_U <- CI_1[3, 2] - Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      )
      Diff <- w_U - width
      while (Diff > tol || is.na(Diff)) {
        N_1 <- N_1 + 1
        CI_1 <- ci_R2(
          R2 = Expected_R2(
            population_R2 = population_R2,
            N = N_1, p = p
          ), alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, N = N_1, p = p,
          random_predictors = TRUE
        )
        w_U <- CI_1[3, 2] - Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        )
        Diff <- w_U - width
      }
    }
    if (is.null(assurance)) {
      Result_Up <- data.frame(term = "necessary_N", value = N_1)
      return(.as_dmar_tbl(Result_Up, conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
    if (!is.null(assurance)) {
      conf_limit_desired_certainty <- ci_R2(R2 = Expected_R2(
        population_R2 = population_R2,
        N = N_1, p = p
      ), alpha_lower = 0, alpha_upper = 1 -
        assurance, N = N_1, p = p, random_predictors = random_predictors)[3, 2]
      N_2 <- N_1
      CI_2 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = conf_limit_desired_certainty,
          N = N_2, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_2, p = p, random_predictors = random_predictors
      )
      w_U <- CI_2[3, 2] - conf_limit_desired_certainty
      Diff <- w_U - width
      while (Diff > tol) {
        N_2 <- N_2 + 1
        CI_2 <- ci_R2(
          R2 = Expected_R2(
            population_R2 = conf_limit_desired_certainty,
            N = N_2, p = p
          ), alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, N = N_2, p = p,
          random_predictors = random_predictors
        )
        w_U <- CI_2[3, 2] - conf_limit_desired_certainty
        Diff <- w_U - width
      }
      w_U_M_2 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_2, p = p, random_predictors = random_predictors
        ),
        alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_2, p = p
      )[3, 2] - Expected_R2(
        population_R2 = population_R2,
        N = N_2, p = p
      )
      Result_Up_2 <- data.frame(term = "necessary_N", value = N_2)
      conf_limit_desired_certainty <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_1, p = p
        ), alpha_lower = 1 - assurance,
        alpha_upper = 0, N = N_1, p = p, random_predictors = random_predictors
      )[1, 2]
      N_3 <- N_1
      CI_2 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = conf_limit_desired_certainty,
          N = N_3, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_3, p = p, random_predictors = random_predictors
      )
      w_U <- CI_2[3, 2] - conf_limit_desired_certainty
      Diff <- w_U - width
      while (Diff > tol) {
        N_3 <- N_3 + 1
        CI_2 <- ci_R2(
          R2 = Expected_R2(
            population_R2 = conf_limit_desired_certainty,
            N = N_3, p = p
          ), alpha_lower = alpha_lower,
          alpha_upper = alpha_upper, N = N_3, p = p,
          random_predictors = random_predictors
        )
        w_U <- CI_2[3, 2] - conf_limit_desired_certainty
        Diff <- w_U - width
      }
      w_U_M_3 <- ci_R2(
        R2 = Expected_R2(
          population_R2 = population_R2,
          N = N_3, p = p
        ), alpha_lower = alpha_lower, alpha_upper = alpha_upper,
        N = N_3, p = p
      )[3, 2] - Expected_R2(
        population_R2 = population_R2,
        N = N_3, p = p, random_predictors = random_predictors
      )
      Result_Up_3 <- data.frame(term = "necessary_N", value = N_3)
      if (N_2 >= N_3) {
        Result_Up <- Result_Up_2
      }
      if (N_3 > N_2) {
        Result_Up <- Result_Up_3
      }
      return(.as_dmar_tbl(Result_Up, conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
  }
  }, warning = function(w) {
    if (inherits(w, "dmar_ncf_clamp")) {
      .clamp_count <<- .clamp_count + 1L
      invokeRestart("muffleWarning")
    }
  })
}
