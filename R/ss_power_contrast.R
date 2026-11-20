#' Sample Size and Statistical Power for a Contrast in a Fixed-Effects ANOVA
#'
#' Determine the necessary per-group sample size for a contrast in a one-way
#' fixed-effects ANOVA so as to achieve a desired level of statistical power,
#' or, alternatively, compute the achieved power for a given sample size. The
#' contrast is specified by a vector of weights and the population means (or
#' a population contrast value) along with the within-group variance.
#'
#' @param c_weights Vector of contrast weights. Required to satisfy
#'   \code{sum(c_weights) == 0}, \code{sum(c_weights[c_weights > 0]) == 1},
#'   and \code{sum(c_weights[c_weights < 0]) == -1}; that is, write the
#'   contrast in normalized fractional form so that the positive coefficients
#'   sum to 1 and the negative coefficients sum to -1. The contrast estimate
#'   is then directly interpretable as a (weighted) mean difference.
#' @param true_means Vector of population group means, of length \code{length(c_weights)}.
#'   Either \code{true_means} or \code{true_contrast} must be supplied (but not both).
#' @param true_var_error Population within-group variance (\eqn{\sigma^2}); must
#'   be a positive number.
#' @param true_contrast Optional. Directly specify the population contrast value
#'   \eqn{\psi = \sum_j c_j \mu_j} as a single number. Use this when the means
#'   themselves are not of interest, only the contrast value.
#' @param desired_power Target statistical power for the test of the contrast
#'   (default \code{0.85}). Used only when \code{specified_n_per_group} is \code{NULL}.
#' @param alpha_level Type I error rate (default \code{0.05}).
#' @param directional Logical. \code{FALSE} (default) yields a two-sided test.
#'   \code{TRUE} yields a one-sided test in the direction of the population
#'   contrast.
#' @param specified_n_per_group Optional per-group sample size at which to evaluate
#'   power. May be a scalar (the common per-group size, equal across groups)
#'   or a numeric vector of length \code{length(c_weights)} giving each
#'   group's size. When \code{NULL} (default), the function instead solves for
#'   the smallest per-group \emph{n} achieving \code{desired_power}.
#' @param print_iter If \code{TRUE}, print the trial \eqn{n} and the
#'   corresponding power as the iterative search proceeds.
#'
#' @return A \code{data.frame} with two columns, \code{term} and
#'   \code{value}:
#'   \describe{
#'     \item{When solving for \emph{n} (\code{specified_n_per_group = NULL})}{Rows are
#'       \code{necessary_n_per_group}, \code{total_N}, \code{actual_power},
#'       \code{noncentral_t_parm}, and \code{true_f} (Cohen's \eqn{f}
#'       for a one-degree-of-freedom contrast).}
#'     \item{When evaluating power (\code{specified_n_per_group} supplied)}{Rows are
#'       \code{specified_n_per_group} (or \code{NA} for unequal \emph{n}),
#'       \code{total_N}, \code{actual_power}, \code{noncentral_t_parm}, and
#'       \code{true_f}.}
#'   }
#'   The result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} summarize
#'   it in broom convention; the reported size is the per-group \emph{n} (or
#'   \code{NA} when unequal group sizes are supplied).
#'
#' @details
#' Let \eqn{\psi = \sum_j c_j \mu_j} be the population contrast. Under the
#' usual fixed-effects ANOVA model with common within-group variance
#' \eqn{\sigma^2} and per-group sizes \eqn{n_j}, the standard error of the
#' contrast estimate is
#' \deqn{\mathrm{SE}_{\hat\psi} = \sqrt{\,\sigma^2 \sum_j c_j^2 / n_j\,},}
#' the test statistic \eqn{t = \hat\psi / \mathrm{SE}_{\hat\psi}} follows a
#' central \emph{t} distribution with \eqn{N - a} degrees of freedom under
#' \eqn{H_0\!: \psi = 0}, and a noncentral \emph{t} distribution with the
#' same df and noncentrality parameter
#' \eqn{\lambda = \psi / \mathrm{SE}_{\hat\psi}} under the alternative.
#'
#' For a two-sided test at level \eqn{\alpha},
#' \deqn{\text{Power} = \Pr(t > t_{1-\alpha/2,\,df}) + \Pr(t < -t_{1-\alpha/2,\,df}),}
#' computed exactly from the noncentral \eqn{t} distribution; for a one-sided
#' test, only the appropriate tail contributes.
#'
#' Cohen's \eqn{f} for a one-df contrast is reported as
#' \eqn{f = |\lambda| / \sqrt{N}}, equivalently \eqn{f^2 = F_{\text{pop}} / N}
#' (Cohen, 1988, Ch. 8); this matches the value that
#' \code{\link[pwr]{pwr.f2.test}} expects when used with \code{u = 1} and
#' \code{v = N - a}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @examples
#' # Four-group example from Maxwell & Delaney's textbook tradition: contrast
#' # the average of three treatment means with a fourth (control), under
#' # population means (90, 92, 88, 81), within-group variance 144.
#' #
#' # 1. Power achieved at n = 20 per group (total N = 80). Should be ~ .80.
#' ss_power_contrast(
#'   c_weights           = c(1/3, 1/3, 1/3, -1),
#'   true_means          = c(90, 92, 88, 81),
#'   true_var_error      = 144,
#'   specified_n_per_group         = 20
#' )
#'
#' # 2. Per-group sample size needed for power = .90.
#' ss_power_contrast(
#'   c_weights           = c(1/3, 1/3, 1/3, -1),
#'   true_means          = c(90, 92, 88, 81),
#'   true_var_error      = 144,
#'   desired_power       = 0.90
#' )
#'
#' # 3. Same effect size specification using a directly given true_contrast.
#' ss_power_contrast(
#'   c_weights           = c(1/3, 1/3, 1/3, -1),
#'   true_contrast       = 9,
#'   true_var_error      = 144,
#'   specified_n_per_group         = 20
#' )
#'
#' # 4. Unequal per-group sample sizes.
#' ss_power_contrast(
#'   c_weights           = c(0.5, 0.5, -0.5, -0.5),
#'   true_means          = c(90, 92, 88, 81),
#'   true_var_error      = 144,
#'   specified_n_per_group         = c(15, 25, 25, 15)
#' )
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{contrast_test}}, \code{\link{ss_power_reg_coef}},
#'   \code{\link{cv_t}}
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
#' @import stats

ss_power_contrast <- function(
  c_weights,
  true_means          = NULL,
  true_var_error      = NULL,
  true_contrast       = NULL,
  desired_power       = 0.85,
  alpha_level         = 0.05,
  directional         = FALSE,
  specified_n_per_group         = NULL,
  print_iter      = FALSE
) {
  .dmar_check_alpha_level(alpha_level)
  # Validate the contrast weights against the user's required normalization:
  #   sum(c)            = 0,
  #   sum(c[c > 0])     = 1,
  #   sum(c[c < 0])     = -1.
  if (!is.numeric(c_weights)) .dmar_refuse_argument(
    "'c_weights' must be a numeric vector.", arg = "c_weights")
  a <- length(c_weights)
  if (a < 2L) .dmar_refuse_argument(
    "'c_weights' must have length at least 2 (i.e., two or more groups).",
    arg = "c_weights")

  tol <- 1e-8
  pos_sum <- sum(c_weights[c_weights > 0])
  neg_sum <- sum(c_weights[c_weights < 0])
  if (abs(sum(c_weights)) > tol) {
    .dmar_refuse_argument("'c_weights' must sum to zero.", arg = "c_weights")
  }
  .dmar_check_contrast_nonzero(c_weights)
  if (abs(pos_sum - 1) > tol) {
    .dmar_refuse_argument(
      "The positive entries of 'c_weights' must sum to 1 (use fractional, normalized weights).",
      arg = "c_weights")
  }
  if (abs(neg_sum + 1) > tol) {
    .dmar_refuse_argument(
      "The negative entries of 'c_weights' must sum to -1.", arg = "c_weights")
  }

  # Validate true_var_error.
  if (is.null(true_var_error) || !is.numeric(true_var_error) ||
      length(true_var_error) != 1L || true_var_error <= 0) {
    .dmar_refuse_argument(
      "'true_var_error' (within-group variance) must be a single positive number.",
      arg = "true_var_error")
  }

  # Validate alpha_level / desired_power.
  if (alpha_level   <= 0 || alpha_level   >= 1) .dmar_refuse_argument(
    "'alpha_level' must be in (0, 1).", arg = "alpha_level")
  .dmar_check_desired_power(desired_power)

  # Resolve the population contrast value (true_contrast, or from true_means).
  if (is.null(true_contrast)) {
    if (is.null(true_means)) {
      .dmar_refuse_argument(
        "Provide either 'true_contrast' (the population contrast value) or 'true_means' (the population means).",
        arg = c("true_contrast", "true_means"))
    }
    if (length(true_means) != a) {
      .dmar_refuse_argument(
        "'true_means' must have length equal to length(c_weights).",
        arg = "true_means")
    }
    Psi_value <- sum(c_weights * true_means)
  } else {
    if (!is.null(true_means)) {
      .dmar_refuse_argument(
        "Provide either 'true_contrast' or 'true_means', not both.",
        arg = c("true_contrast", "true_means"))
    }
    if (!is.numeric(true_contrast) || length(true_contrast) != 1L) {
      .dmar_refuse_argument(
        "'true_contrast' must be a single number.", arg = "true_contrast")
    }
    Psi_value <- true_contrast
  }

  # Inner helper: exact noncentral t power for a vector of per-group sample
  # sizes. Accepts a scalar (treated as common per-group n) or a vector.
  power_at <- function(n_per) {
    if (length(n_per) == 1L) n_per <- rep(n_per, a)
    SE   <- sqrt(true_var_error * sum(c_weights^2 / n_per))
    df   <- sum(n_per) - a
    ncp  <- Psi_value / SE
    power <- .power_noncentral_t(ncp, df, alpha_level, directional)
    list(power = power, ncp = ncp, df = df, SE = SE,
         total.n = sum(n_per), n_per = n_per)
  }

  if (is.null(specified_n_per_group)) {
    # Solve for the smallest per-group n achieving desired_power.
    if (abs(Psi_value) < tol) {
      .dmar_refuse_argument(
        "Population contrast 'true_contrast' is zero (or numerically zero); statistical power cannot exceed alpha_level.",
        arg = "true_contrast")
    }
    n_i <- 2L
    repeat {
      n_i <- n_i + 1L
      r <- power_at(n_i)
      if (isTRUE(print_iter)) {
        cat(sprintf("specified_n_per_group = %d, total_N = %d, power = %.5f\n",
                    n_i, r$total.n, r$power))
      }
      if (r$power >= desired_power) break
      if (n_i > 100000L) {
        .dmar_refuse_argument(
          "No solution found below specified_n_per_group = 100,000; check effect size and power.",
          arg = NULL)
      }
    }
    f <- abs(r$ncp) / sqrt(r$total.n)
    out <- data.frame(
      term  = c("necessary_n_per_group", "total_N", "actual_power",
                "noncentral_t_parm", "true_f"),
      value = c(n_i, r$total.n, r$power, abs(r$ncp), f),
      stringsAsFactors = FALSE
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  # Otherwise: evaluate power at the user-supplied sample size(s).
  n_per <- specified_n_per_group
  if (length(n_per) == 1L) {
    n_per <- rep(n_per, a)
  } else if (length(n_per) != a) {
    .dmar_refuse_argument(
      "'specified_n_per_group' must be a single number (per-group n) or a numeric vector of length equal to length(c_weights).",
      arg = "specified_n_per_group")
  }
  if (any(n_per < 2)) .dmar_refuse_argument(
    "Each per-group sample size must be at least 2.", arg = NULL)
  if (any(c_weights != 0 & n_per < 2)) .dmar_refuse_argument(
    "Groups with non-zero contrast weight must have n >= 2.", arg = NULL)

  r <- power_at(n_per)
  per_group_value <- if (length(unique(n_per)) == 1L) n_per[1L] else NA_real_
  f <- abs(r$ncp) / sqrt(r$total.n)
  out <- data.frame(
    term  = c("specified_n_per_group", "total_N", "actual_power",
              "noncentral_t_parm", "true_f"),
    value = c(per_group_value, r$total.n, r$power, abs(r$ncp), f),
    stringsAsFactors = FALSE
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
