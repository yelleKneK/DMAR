# This file groups three exported functions for the one-way within-subjects
# design: mauchly_test(), epsilon_corrections(), and anova_within(). They are a
# tight family with shared math, the package's one sanctioned reason to place
# several exported functions in a single source file. The sanctioned
# grouping rationale: mauchly_test(), epsilon_corrections(), and anova_within()
# share the internal within-subjects contrast-covariance helpers
# (.within_to_matrix and .within_contrast_cov), and anova_within() composes the
# other two, so they are intentionally kept in one file. Both public helpers
# and the ANOVA driver reuse the two internal helpers below, .within_to_matrix()
# and .within_contrast_cov(); anova_within() in turn calls mauchly_test() and
# epsilon_corrections() to attach the sphericity diagnostics. Keeping them
# together keeps the contrast-space covariance and the reshape logic in one
# place rather than duplicated across three files.

# --- Internal: turn input (matrix or long-format) into an n x k wide matrix.
.within_to_matrix <- function(x, id = NULL, time = NULL, outcome = NULL) {
  if (is.matrix(x) || is.data.frame(x)) {
    if (is.null(id) && is.null(time) && is.null(outcome)) {
      # Wide format already
      m <- as.matrix(x)
      if (!is.numeric(m)) {
        stop("'x' must be a numeric matrix or data.frame.", call. = FALSE)
      }
      if (anyNA(m)) {
        stop("Missing values are not supported; remove or impute first.",
             call. = FALSE)
      }
      return(m)
    }
    # Long format
    if (is.null(id) || is.null(time) || is.null(outcome)) {
      stop("For long-format input, supply 'id', 'time', and 'outcome' column names.",
           call. = FALSE)
    }
    for (col in c(id, time, outcome)) {
      if (!col %in% names(x)) stop("Column '", col, "' not in 'x'.", call. = FALSE)
    }
    if (anyNA(x[[outcome]]) || anyNA(x[[id]]) || anyNA(x[[time]])) {
      stop("Missing values are not supported; remove or impute first.",
           call. = FALSE)
    }
    m <- stats::reshape(
      data = as.data.frame(x[, c(id, time, outcome)]),
      timevar   = time,
      idvar     = id,
      direction = "wide"
    )
    out <- as.matrix(m[, -1, drop = FALSE])
    rownames(out) <- m[[1]]
    storage.mode(out) <- "double"
    return(out)
  }
  stop("'x' must be a numeric matrix or a data.frame.", call. = FALSE)
}


# --- Internal: residual covariance matrix in orthonormal-contrast space.
# Given an n x k wide matrix, build a (k-1) x k Helmert-like contrast that
# orthonormalizes the (k-1)-dim space orthogonal to the constant vector,
# then return the (k-1) x (k-1) covariance matrix of the contrast scores.
.within_contrast_cov <- function(Y) {
  k <- ncol(Y)
  # Helmert contrasts (k-1 columns), then orthonormalize.
  H <- stats::contr.helmert(k)        # k x (k-1)
  M <- t(qr.Q(qr(H)))                  # (k-1) x k orthonormal rows
  S <- stats::cov(Y %*% t(M))          # (k-1) x (k-1) covariance
  list(S = S, n = nrow(Y), k = k, M = M)
}


#' Mauchly's Test of Sphericity for a One-Way Within-Subjects Design
#'
#' Tests the null hypothesis that the covariance matrix of the orthonormal
#' contrasts among the \eqn{k} repeated measurements is proportional to the
#' identity (the \emph{sphericity} assumption underlying univariate
#' repeated measures \emph{F}-tests).
#'
#' @param x Either an \eqn{n \times k} numeric matrix or
#'   \code{data.frame} (rows = subjects, columns = repeated measurements);
#'   \emph{or} a long-format \code{data.frame} together with \code{id},
#'   \code{time}, and \code{outcome} column names.
#' @param id Column name in \code{x} identifying the subject when \code{x} is in long format (\code{NULL} otherwise).
#' @param time Column name in \code{x} identifying the within-subjects factor level when \code{x} is in long format (\code{NULL} otherwise).
#' @param outcome Column name in \code{x} identifying the dependent variable when \code{x} is in long format (\code{NULL} otherwise).
#'
#' @return A one-row \code{data.frame} with columns \code{W} (Mauchly's
#'   statistic), \code{statistic} (the chi square approximation),
#'   \code{df}, \code{p_value}, \code{n_subjects}, \code{n_levels}, and
#'   \code{method}.
#'
#' @details Sphericity is the assumption that the variances of all pairwise
#'   differences among the \eqn{k} levels are equal, equivalently, that the
#'   covariance matrix \eqn{\Sigma_C} of any orthonormal set of \eqn{k - 1}
#'   contrasts among the levels is proportional to the identity. Mauchly's
#'   (1940) test statistic is
#'   \deqn{W = \frac{\det(\hat\Sigma_C)}{\bigl(\mathrm{tr}(\hat\Sigma_C) / (k - 1)\bigr)^{k - 1}},}
#'   and the chi square approximation
#'   \deqn{X^2 = -\,m \,\log W \quad \mathrm{with}\ m = (n - 1) - \frac{2(k - 1)^2 + (k - 1) + 2}{6\,(k - 1)}}
#'   has approximately \eqn{(k - 1)k/2 - 1} degrees of freedom under
#'   \eqn{H_0}. The reported \emph{p}-value uses Box's (1949) second-order
#'   correction, a weighted combination of the chi square tails on
#'   \eqn{(k - 1)k/2 - 1} and \eqn{(k - 1)k/2 + 3} degrees of freedom, which
#'   improves the first-order approximation in small samples; this matches
#'   \code{\link[stats]{mauchly.test}}.
#'
#'   When sphericity is rejected, the univariate \emph{F} test is liberal;
#'   correct using the Greenhouse-Geisser, Huynh-Feldt, or lower-bound
#'   epsilon adjustments via \code{\link{epsilon_corrections}} or directly
#'   via \code{\link{anova_within}}.
#'
#'   The test is only defined for \eqn{k \ge 3}; with \eqn{k = 2}, sphericity
#'   is trivially true and the function returns \code{W = 1, p = 1}.
#'
#' @references
#' Mauchly, J. W. (1940). Significance test for sphericity of a normal
#'   \eqn{n}-variate distribution. \emph{Annals of Mathematical Statistics,
#'   11}(2), 204--209.
#'
#' Box, G. E. P. (1949). A general distribution theory for a class of
#'   likelihood criteria. \emph{Biometrika, 36}(3/4), 317--346.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 11 for sphericity in within-subjects
#'   designs.)
#'
#' @examples
#' # Wide-format example: simulated within-subjects data with 4 levels.
#' set.seed(113)
#' Y <- matrix(rnorm(20 * 4), nrow = 20)
#' mauchly_test(Y)
#'
#' # Long-format example using built-in nlme::Orthodont (4 ages per subject).
#' mauchly_test(nlme::Orthodont, id = "Subject", time = "age",
#'              outcome = "distance")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{epsilon_corrections}}, \code{\link{anova_within}}
#'
#' @keywords htest
#'
#' @family within-subjects analysis
#' @family hypothesis tests
#'
#' @export
#' @import stats
mauchly_test <- function(x, id = NULL, time = NULL, outcome = NULL) {
  Y <- .within_to_matrix(x, id, time, outcome)
  n <- nrow(Y); k <- ncol(Y)
  if (n < 3L) stop("At least 3 subjects are required.", call. = FALSE)
  if (k < 2L) stop("At least 2 within-subjects levels are required.", call. = FALSE)

  if (k == 2L) {
    # Sphericity is trivially true with two levels.
    return(.as_dmar_tbl(data.frame(
      W = 1, statistic = 0, df = 0, p_value = 1,
      n_subjects = n, n_levels = k,
      method = "Mauchly's test of sphericity",
      stringsAsFactors = FALSE
    )))
  }

  cc  <- .within_contrast_cov(Y)
  S   <- cc$S
  q   <- k - 1
  W   <- as.numeric(det(S) / (sum(diag(S)) / q)^q)
  W   <- max(W, .Machine$double.eps)
  m   <- (n - 1) - (2 * q^2 + q + 2) / (6 * q)
  Xsq <- -m * log(W)
  df  <- q * (q + 1) / 2 - 1
  # Box's (1949) second-order chi square approximation, as used by
  # stats::mauchly.test: the tail probability is a weighted combination of
  # chi square tails on df and df + 4 degrees of freedom, with weight w2. The
  # first-order term alone (Pr1) is biased for small samples; the w2 term is
  # the higher-order correction. Here m = (n - 1) rho is the Bartlett-type
  # multiplier, so R's (n rho q)^2 becomes (m q)^2.
  w2  <- (q + 2) * (q - 1) * (q - 2) * (2 * q^3 + 6 * q^2 + 3 * k + 2) /
    (288 * (m * q)^2)
  Pr1 <- stats::pchisq(Xsq, df = df, lower.tail = FALSE)
  Pr2 <- stats::pchisq(Xsq, df = df + 4, lower.tail = FALSE)
  p   <- Pr1 + w2 * (Pr2 - Pr1)

  out <- data.frame(
    W          = W,
    statistic  = Xsq,
    df         = df,
    p_value    = p,
    n_subjects = n,
    n_levels   = k,
    method     = "Mauchly's test of sphericity",
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
  .as_dmar_tbl(out)
}


#' Greenhouse-Geisser, Huynh-Feldt, and Lower-Bound Epsilon Corrections
#'
#' Computes the three standard sphericity-correction factors for the
#' univariate within-subjects \emph{F} test. When sphericity holds, all three
#' equal 1; departures from sphericity reduce them, deflating the effective
#' degrees of freedom and thus tempering the inflated Type I error rate of
#' the unadjusted univariate test.
#'
#' @param x Either an \eqn{n \times k} numeric matrix or
#'   \code{data.frame} (rows = subjects, columns = repeated measurements);
#'   \emph{or} a long-format \code{data.frame} together with \code{id},
#'   \code{time}, and \code{outcome} column names.
#' @param id Column name in \code{x} identifying the subject when \code{x} is in long format (\code{NULL} otherwise).
#' @param time Column name in \code{x} identifying the within-subjects factor level when \code{x} is in long format (\code{NULL} otherwise).
#' @param outcome Column name in \code{x} identifying the dependent variable when \code{x} is in long format (\code{NULL} otherwise).
#'
#' @return A \code{data.frame} with columns \code{epsilon_method}
#'   (\code{"Greenhouse-Geisser"}, \code{"Huynh-Feldt"}, \code{"lower_bound"})
#'   and \code{epsilon} (the correction factor in \eqn{[1/(k-1), 1]}).
#'
#' @details For an \eqn{(k - 1) \times (k - 1)} covariance matrix
#'   \eqn{\hat\Sigma_C} of orthonormal contrasts among the \eqn{k} repeated
#'   measurements (with eigenvalues \eqn{\lambda_1, \ldots, \lambda_{k-1}}):
#'   \deqn{\hat\varepsilon_{\mathrm{GG}} = \frac{(\sum \lambda_i)^2}{(k - 1)\,\sum \lambda_i^2},}
#'   \deqn{\hat\varepsilon_{\mathrm{HF}} = \min\!\Bigl(1,\ \frac{n(k - 1)\hat\varepsilon_{\mathrm{GG}} - 2}{(k - 1)\bigl(n - 1 - (k - 1)\hat\varepsilon_{\mathrm{GG}}\bigr)}\Bigr),}
#'   \deqn{\hat\varepsilon_{\mathrm{LB}} = \frac{1}{k - 1}.}
#'   The Greenhouse-Geisser \eqn{\hat\varepsilon} tends to be conservative;
#'   the Huynh-Feldt correction adjusts upward to be (approximately)
#'   unbiased; the lower bound is the worst-case adjustment.
#'
#' @references
#' Greenhouse, S. W., & Geisser, S. (1959). On methods in the analysis of
#'   profile data. \emph{Psychometrika, 24}(2), 95--112.
#'
#' Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction for
#'   degrees of freedom from sample data in randomized block and split-plot
#'   designs. \emph{Journal of Educational Statistics, 1}(1), 69--82.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @examples
#' set.seed(113)
#' Y <- matrix(rnorm(20 * 4), nrow = 20)
#' epsilon_corrections(Y)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{mauchly_test}}, \code{\link{anova_within}}
#'
#' @keywords htest
#'
#' @family within-subjects analysis
#'
#' @export
#' @import stats
epsilon_corrections <- function(x, id = NULL, time = NULL, outcome = NULL) {
  Y <- .within_to_matrix(x, id, time, outcome)
  n <- nrow(Y); k <- ncol(Y)
  if (n < 3L) stop("At least 3 subjects are required.", call. = FALSE)
  if (k < 2L) stop("At least 2 within-subjects levels are required.", call. = FALSE)

  cc  <- .within_contrast_cov(Y)
  S   <- cc$S
  lam <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  q   <- k - 1

  if (n - 1 < q) {
    # The contrast covariance matrix has rank at most n - 1 < q, so the
    # sample epsilon is bounded above by (n - 1)/q regardless of the
    # population epsilon: an artifact of rank deficiency, not an
    # estimate. Decline to estimate, as in anova_within_two_way(); the
    # a priori lower bound 1/q needs nothing from S and remains valid
    # (Maxwell, Delaney, & Kelley, 2027, Sections 11.9.1 and 11.9.2).
    warning("With n = ", n, " subjects and ", k, " levels, n - 1 < ",
            q, " so the contrast covariance matrix is singular and the ",
            "Greenhouse-Geisser and Huynh-Feldt epsilons cannot be ",
            "estimated; their rows are NA. The lower bound epsilon = 1/",
            q, " remains valid.", call. = FALSE)
    eps_gg <- NA_real_
    eps_hf <- NA_real_
  } else {
    eps_gg <- (sum(lam))^2 / (q * sum(lam^2))
    hf_num <- n * q * eps_gg - 2
    hf_den <- q * (n - 1 - q * eps_gg)
    eps_hf <- if (hf_den == 0) NA_real_ else min(1, hf_num / hf_den)
  }
  eps_lb <- 1 / q

  out <- data.frame(
    epsilon_method = c("Greenhouse-Geisser", "Huynh-Feldt", "lower_bound"),
    epsilon        = c(eps_gg, eps_hf, eps_lb),
    stringsAsFactors = FALSE
  )
  .as_dmar_tbl(out)
}


#' One Way Within-Subjects ANOVA With Sphericity Diagnostics and Corrections
#'
#' Performs the univariate one-way within-subjects \emph{F} test together
#' with Mauchly's test of sphericity and the three standard
#' \eqn{\varepsilon}-corrected \emph{p}-values (Greenhouse-Geisser,
#' Huynh-Feldt, and lower-bound). Returns everything in a single tidy
#' \code{data.frame} so the user can decide which adjustment to report.
#'
#' @param x Either an \eqn{n \times k} numeric matrix or
#'   \code{data.frame} (rows = subjects, columns = repeated measurements);
#'   \emph{or} a long-format \code{data.frame} together with \code{id},
#'   \code{time}, and \code{outcome} column names.
#' @param id Column name in \code{x} identifying the subject when \code{x} is in long format (\code{NULL} otherwise).
#' @param time Column name in \code{x} identifying the within-subjects factor level when \code{x} is in long format (\code{NULL} otherwise).
#' @param outcome Column name in \code{x} identifying the dependent variable when \code{x} is in long format (\code{NULL} otherwise).
#'
#' @return A \code{data.frame} with one row per reported \emph{F} test:
#'   \code{adjustment} (\code{"none"}, \code{"Greenhouse-Geisser"},
#'   \code{"Huynh-Feldt"}, \code{"lower_bound"}), \code{F_value},
#'   \code{df_1}, \code{df_2}, \code{p_value}, and \code{epsilon} (the
#'   correction factor used; \code{NA} for the unadjusted row).
#'   \code{attr(<output>, "mauchly")} contains the row from
#'   \code{\link{mauchly_test}}, and the partial \eqn{\eta^2} is attached
#'   as \code{attr(<output>, "partial_eta_squared")}.
#'
#' @details The unadjusted within-subjects \emph{F} statistic is the same
#'   regardless of sphericity; corrections shrink the numerator and
#'   denominator degrees of freedom by a factor of
#'   \eqn{\hat\varepsilon \in [1/(k - 1),\, 1]}, and the \emph{p}-value is
#'   recomputed against the adjusted reference \emph{F} distribution. When
#'   Mauchly's test rejects, prefer the Huynh-Feldt-corrected \emph{p}-value
#'   (less conservative than Greenhouse-Geisser).
#'
#'   For multi-factor within-subjects designs or mixed designs, fit the
#'   model with \code{stats::aov(... + Error(id/within))} or with
#'   \code{lme4::lmer()} directly.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 11.)
#'
#' @examples
#' # Simulated within-subjects data with no real effect.
#' set.seed(113)
#' Y <- matrix(rnorm(20 * 4), nrow = 20)
#' anova_within(Y)
#'
#' # Built-in within-subjects example: nlme::Orthodont (distance ~ age).
#' res <- anova_within(nlme::Orthodont,
#'                     id = "Subject", time = "age", outcome = "distance")
#' res
#' attr(res, "mauchly")
#' attr(res, "partial_eta_squared")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{mauchly_test}}, \code{\link{epsilon_corrections}},
#'   \code{\link[stats]{aov}}
#'
#' @keywords htest design
#'
#' @family within-subjects analysis
#' @family hypothesis tests
#'
#' @export
#' @import stats
anova_within <- function(x, id = NULL, time = NULL, outcome = NULL) {
  Y <- .within_to_matrix(x, id, time, outcome)
  n <- nrow(Y); k <- ncol(Y)
  if (n < 2L) stop("At least 2 subjects are required.", call. = FALSE)
  if (k < 2L) stop("At least 2 within-subjects levels are required.", call. = FALSE)

  # Sums of squares decomposition for one-way within-subjects ANOVA.
  grand_mean <- mean(Y)
  cond_means <- colMeans(Y)
  subj_means <- rowMeans(Y)

  ss_total      <- sum((Y - grand_mean)^2)
  ss_subjects   <- k * sum((subj_means - grand_mean)^2)
  ss_condition  <- n * sum((cond_means - grand_mean)^2)
  ss_residual   <- ss_total - ss_subjects - ss_condition

  df_condition  <- k - 1
  df_residual   <- (n - 1) * (k - 1)

  ms_condition  <- ss_condition / df_condition
  ms_residual   <- ss_residual / df_residual
  F_value       <- ms_condition / ms_residual
  partial_eta_sq <- ss_condition / (ss_condition + ss_residual)

  # Sphericity diagnostics and corrections.
  mauchly <- mauchly_test(Y)
  eps_tbl <- epsilon_corrections(Y)
  eps_gg  <- eps_tbl$epsilon[eps_tbl$epsilon_method == "Greenhouse-Geisser"]
  eps_hf  <- eps_tbl$epsilon[eps_tbl$epsilon_method == "Huynh-Feldt"]
  eps_lb  <- eps_tbl$epsilon[eps_tbl$epsilon_method == "lower_bound"]

  build_row <- function(label, eps) {
    if (is.na(eps)) {
      data.frame(adjustment = label, F_value = F_value,
                 df_1 = NA_real_, df_2 = NA_real_,
                 p_value = NA_real_, epsilon = NA_real_,
                 stringsAsFactors = FALSE)
    } else {
      d1 <- df_condition * eps
      d2 <- df_residual  * eps
      p  <- stats::pf(F_value, df1 = d1, df2 = d2, lower.tail = FALSE)
      data.frame(adjustment = label, F_value = F_value,
                 df_1 = d1, df_2 = d2, p_value = p, epsilon = eps,
                 stringsAsFactors = FALSE)
    }
  }

  out <- rbind(
    data.frame(adjustment = "none",
               F_value = F_value,
               df_1 = df_condition, df_2 = df_residual,
               p_value = stats::pf(F_value, df1 = df_condition,
                                   df2 = df_residual,
                                   lower.tail = FALSE),
               epsilon = NA_real_,
               stringsAsFactors = FALSE),
    build_row("Greenhouse-Geisser", eps_gg),
    build_row("Huynh-Feldt",        eps_hf),
    build_row("lower_bound",        eps_lb)
  )
  rownames(out) <- NULL

  attr(out, "mauchly")             <- mauchly
  attr(out, "partial_eta_squared") <- partial_eta_sq
  attr(out, "n_subjects")          <- n
  attr(out, "n_levels")            <- k
  .as_dmar_tbl(out)
}
