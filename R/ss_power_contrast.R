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
#' @param mu Vector of population group means, of length \code{length(c_weights)}.
#'   Either \code{mu} or \code{psi} must be supplied (but not both).
#' @param sigma_squared Population within-group variance (\eqn{\sigma^2}); must
#'   be a positive number.
#' @param psi Optional. Directly specify the population contrast value
#'   \eqn{\psi = \sum_j c_j \mu_j} as a single number. Use this when the means
#'   themselves are not of interest, only the contrast value.
#' @param desired_power Target statistical power for the test of the contrast
#'   (default \code{0.85}). Used only when \code{n_per_group} is \code{NULL}.
#' @param alpha_level Type I error rate (default \code{0.05}).
#' @param directional Logical. \code{FALSE} (default) yields a two-sided test.
#'   \code{TRUE} yields a one-sided test in the direction of the population
#'   contrast.
#' @param n_per_group Optional per-group sample size at which to evaluate
#'   power. May be a scalar (the common per-group size, equal across groups)
#'   or a numeric vector of length \code{length(c_weights)} giving each
#'   group's size. When \code{NULL} (default), the function instead solves for
#'   the smallest per-group \emph{n} achieving \code{desired_power}.
#' @param print_progress If \code{TRUE}, print the trial \eqn{n} and the
#'   corresponding power as the iterative search proceeds.
#'
#' @return A \code{data.frame} with two columns, \code{term} and
#'   \code{value}:
#'   \describe{
#'     \item{When solving for \emph{n} (\code{n_per_group = NULL})}{Rows are
#'       \code{necessary_n_per_group}, \code{total_N}, \code{actual_power},
#'       \code{noncentral_t_parm}, and \code{effect_size_f} (Cohen's \eqn{f}
#'       for a one-degree-of-freedom contrast).}
#'     \item{When evaluating power (\code{n_per_group} supplied)}{Rows are
#'       \code{specified_n_per_group} (or \code{NA} for unequal \emph{n}),
#'       \code{total_N}, \code{actual_power}, \code{noncentral_t_parm}, and
#'       \code{effect_size_f}.}
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
#'   c_weights     = c(1/3, 1/3, 1/3, -1),
#'   mu            = c(90, 92, 88, 81),
#'   sigma_squared = 144,
#'   n_per_group   = 20
#' )
#'
#' # 2. Per-group sample size needed for power = .90.
#' ss_power_contrast(
#'   c_weights     = c(1/3, 1/3, 1/3, -1),
#'   mu            = c(90, 92, 88, 81),
#'   sigma_squared = 144,
#'   desired_power = 0.90
#' )
#'
#' # 3. Same effect size specification using a directly given psi.
#' ss_power_contrast(
#'   c_weights     = c(1/3, 1/3, 1/3, -1),
#'   psi           = 9,
#'   sigma_squared = 144,
#'   n_per_group   = 20
#' )
#'
#' # 4. Unequal per-group sample sizes.
#' ss_power_contrast(
#'   c_weights     = c(0.5, 0.5, -0.5, -0.5),
#'   mu            = c(90, 92, 88, 81),
#'   sigma_squared = 144,
#'   n_per_group   = c(15, 25, 25, 15)
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
  mu             = NULL,
  sigma_squared  = NULL,
  psi            = NULL,
  desired_power  = 0.85,
  alpha_level    = 0.05,
  directional    = FALSE,
  n_per_group    = NULL,
  print_progress = FALSE
) {
  # Validate the contrast weights against the user's required normalization:
  #   sum(c)            = 0,
  #   sum(c[c > 0])     = 1,
  #   sum(c[c < 0])     = -1.
  if (!is.numeric(c_weights)) stop("'c_weights' must be a numeric vector.")
  a <- length(c_weights)
  if (a < 2L) stop("'c_weights' must have length at least 2 (i.e., two or more groups).")

  tol <- 1e-8
  pos_sum <- sum(c_weights[c_weights > 0])
  neg_sum <- sum(c_weights[c_weights < 0])
  if (abs(sum(c_weights)) > tol) {
    stop("'c_weights' must sum to zero.")
  }
  if (abs(pos_sum - 1) > tol) {
    stop("The positive entries of 'c_weights' must sum to 1 (use fractional, normalized weights).")
  }
  if (abs(neg_sum + 1) > tol) {
    stop("The negative entries of 'c_weights' must sum to -1.")
  }

  # Validate sigma_squared.
  if (is.null(sigma_squared) || !is.numeric(sigma_squared) ||
      length(sigma_squared) != 1L || sigma_squared <= 0) {
    stop("'sigma_squared' (within-group variance) must be a single positive number.")
  }

  # Validate alpha_level / desired_power.
  if (alpha_level   <= 0 || alpha_level   >= 1) stop("'alpha_level' must be in (0, 1).")
  if (desired_power <= 0 || desired_power >= 1) stop("'desired_power' must be in (0, 1).")

  # Resolve the population contrast value psi.
  if (is.null(psi)) {
    if (is.null(mu)) {
      stop("Provide either 'psi' (the population contrast value) or 'mu' (the population means).")
    }
    if (length(mu) != a) {
      stop("'mu' must have length equal to length(c_weights).")
    }
    Psi_value <- sum(c_weights * mu)
  } else {
    if (!is.null(mu)) {
      stop("Provide either 'psi' or 'mu', not both.")
    }
    if (!is.numeric(psi) || length(psi) != 1L) {
      stop("'psi' must be a single number.")
    }
    Psi_value <- psi
  }

  # Inner helper: exact noncentral t power for a vector of per-group sample
  # sizes. Accepts a scalar (treated as common per-group n) or a vector.
  power_at <- function(n_per) {
    if (length(n_per) == 1L) n_per <- rep(n_per, a)
    SE   <- sqrt(sigma_squared * sum(c_weights^2 / n_per))
    df   <- sum(n_per) - a
    ncp  <- Psi_value / SE
    if (directional) {
      t_crit <- stats::qt(1 - alpha_level, df = df)
      power  <- 1 - stats::pt(t_crit, df = df, ncp = abs(ncp))
    } else {
      t_crit <- stats::qt(1 - alpha_level / 2, df = df)
      power  <- (1 - stats::pt( t_crit, df = df, ncp = abs(ncp))) +
                       stats::pt(-t_crit, df = df, ncp = abs(ncp))
    }
    list(power = power, ncp = ncp, df = df, SE = SE,
         total.n = sum(n_per), n_per = n_per)
  }

  if (is.null(n_per_group)) {
    # Solve for the smallest per-group n achieving desired_power.
    if (abs(Psi_value) < tol) {
      stop("Population contrast 'psi' is zero (or numerically zero); statistical power cannot exceed alpha_level.")
    }
    n_i <- 2L
    repeat {
      n_i <- n_i + 1L
      r <- power_at(n_i)
      if (isTRUE(print_progress)) {
        cat(sprintf("n_per_group = %d, total_N = %d, power = %.5f\n",
                    n_i, r$total.n, r$power))
      }
      if (r$power >= desired_power) break
      if (n_i > 100000L) {
        stop("No solution found below n_per_group = 100,000; check effect size and power.")
      }
    }
    f <- abs(r$ncp) / sqrt(r$total.n)
    out <- data.frame(
      term  = c("necessary_n_per_group", "total_N", "actual_power",
                "noncentral_t_parm", "effect_size_f"),
      value = c(n_i, r$total.n, r$power, abs(r$ncp), f),
      stringsAsFactors = FALSE
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  # Otherwise: evaluate power at the user-supplied sample size(s).
  n_per <- n_per_group
  if (length(n_per) == 1L) {
    n_per <- rep(n_per, a)
  } else if (length(n_per) != a) {
    stop("'n_per_group' must be a single number (per-group n) or a numeric vector of length equal to length(c_weights).")
  }
  if (any(n_per < 2)) stop("Each per-group sample size must be at least 2.")
  if (any(c_weights != 0 & n_per < 2)) stop("Groups with non-zero contrast weight must have n >= 2.")

  r <- power_at(n_per)
  per_group_value <- if (length(unique(n_per)) == 1L) n_per[1L] else NA_real_
  f <- abs(r$ncp) / sqrt(r$total.n)
  out <- data.frame(
    term  = c("specified_n_per_group", "total_N", "actual_power",
              "noncentral_t_parm", "effect_size_f"),
    value = c(per_group_value, r$total.n, r$power, abs(r$ncp), f),
    stringsAsFactors = FALSE
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
