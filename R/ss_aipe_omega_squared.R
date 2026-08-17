# Sample size for AIPE on omega squared (ANOVA effect size).
#' Sample Size for AIPE on Omega Squared (ANOVA Effect Size)
#'
#' Determines the sample size needed for the noncentral \emph{F}
#' confidence interval on the population omega squared (\eqn{\omega^2}) to
#' have a desired width (Accuracy in Parameter Estimation; Kelley, 2008;
#' Steiger, 2004). The function uses the same noncentral \emph{F}
#' machinery as \code{\link{ci_omega_squared}}: at each candidate \eqn{N}
#' it computes the expected CI width by inverting the noncentral \emph{F}
#' distribution and stops at the smallest \eqn{N} that achieves the
#' target width.
#'
#' @param population_omega_squared Anticipated population \eqn{\omega^2}.
#'   Must lie in \eqn{[0, 1)}.
#' @param df_effect Numerator degrees of freedom for the effect (e.g.,
#'   \eqn{a - 1} for an \eqn{a}-group one-way ANOVA, or the appropriate
#'   per-effect numerator df in a factorial design).
#' @param width Desired full width of the CI on \eqn{\omega^2}.
#' @param which_width Whether \code{width} is the \code{"Full"} width of
#'   the interval (default) or a half-width: \code{"Lower"} and
#'   \code{"Upper"} both interpret \code{width} as half the full width,
#'   so they plan for a full width of twice \code{width} and return the
#'   same sample size. Because the interval is not generally symmetric
#'   about the estimate, its realized lower and upper half-widths can
#'   differ from each other and from half the full width; the planner
#'   does not target them separately. A genuinely one-sided width target
#'   is not currently offered.
#' @param conf_level Desired confidence level (default \code{0.95}).
#' @param assurance Optional. Probability that the realized CI is no
#'   wider than \code{width}; when supplied, the sample size is inflated
#'   by the standard chi squared correction (Kelley, 2008).
#'
#' @return A \code{data.frame} with rows for the recommended
#'   \emph{total} sample size \code{necessary_N}, the expected CI width
#'   at that sample size, and the inputs echoed back.
#'
#' @details
#' \strong{Connection to noncentral \emph{F} machinery.} The CI on
#' \eqn{\omega^2} is built by inverting the noncentral \emph{F} sampling
#' distribution of the observed \emph{F} statistic, following Steiger
#' (2004) and Kelley (2007); see \code{\link{ci_omega_squared}}. To plan
#' a sample size, we iterate: for each candidate \eqn{N}, compute the
#' \emph{F} the analyst would observe \emph{at the population effect
#' size}, build its CI on \eqn{\omega^2}, and stop at the smallest \eqn{N}
#' whose CI width is below the target.
#'
#' \strong{Population-effect-to-\emph{F} mapping.} Given a target
#' \eqn{\omega^2}, the expected sample \emph{F} that yields exactly that
#' \eqn{\omega^2} as the point estimate from
#' \eqn{\hat\omega^2 = df_{\text{eff}}(F - 1) / [df_{\text{eff}}(F - 1) + N]}
#' is \eqn{F = 1 + \omega^2 N / [df_{\text{eff}} (1 - \omega^2)]}. This is
#' the \emph{F} value used at each iteration of the search.
#'
#' \strong{Tolerance behavior at small \emph{N}.} For small candidate
#' \emph{N} the noncentral \emph{F} lower limit is often clamped to zero
#' (see \code{?conf_limits_ncf}). The search ignores these clamps in the
#' iteration and reports the final clamp count, if any, as an informational
#' message; this matches the convention in \code{\link{ss_aipe_R2}}.
#'
#' @references
#' Algina, J., Moulder, B. C., & Moser, B. K. (2002). Sample size
#'   requirements for accurate estimation of squared semi-partial
#'   correlation coefficients. \emph{Multivariate Behavioral Research,
#'   37}(1), 37--57. \doi{10.1207/s15327906mbr3701_02}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43}(4),
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Kelley, K., & Preacher, K. J. (2012). On effect size.
#'   \emph{Psychological Methods, 17}, 137--152. \doi{10.1037/a0028086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{\eta^2}, Chapter 7 on
#'   factorial designs, and Chapter 11 on generalized \eqn{\eta^2} for
#'   within-subjects designs.)
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#'   intervals and tests of close fit in the analysis of variance and
#'   contrast analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @seealso \code{\link{ci_omega_squared}}, \code{\link{omega_squared}},
#'   \code{\link{omega_squared_partial}}, \code{\link{ss_aipe_R2}}
#'
#' @examples
#' # 1. Plan total N so the 95% CI on omega^2 has full width <= 0.10
#' #        in a 3-group one-way ANOVA (df_effect = 2), anticipated
#' #        omega^2 = 0.10.
#' ss_aipe_omega_squared(population_omega_squared = 0.10,
#'                       df_effect = 2,
#'                       width = 0.10)
#'
#' # 2. Same problem with 80% assurance:
#' ss_aipe_omega_squared(population_omega_squared = 0.10,
#'                       df_effect = 2,
#'                       width = 0.10,
#'                       assurance = 0.80)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export

ss_aipe_omega_squared <- function(population_omega_squared,
                                  df_effect,
                                  width,
                                  which_width = c("Full", "Lower", "Upper"),
                                  conf_level = 0.95,
                                  assurance = NULL) {
  which_width <- match.arg(which_width)
  if (!is.numeric(population_omega_squared) ||
      length(population_omega_squared) != 1L ||
      population_omega_squared < 0 || population_omega_squared >= 1)
    stop("'population_omega_squared' must be a single value in [0, 1).")
  if (!is.numeric(df_effect) || length(df_effect) != 1L || df_effect < 1)
    stop("'df_effect' must be a single integer >= 1.")
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    stop("'width' must be a single positive number.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  target_width <- if (which_width == "Full") width else 2 * width

  # Map omega^2 + df_effect + N to the F that produces exactly that
  # omega^2 as the point estimate.
  F_from_N <- function(N) {
    1 + population_omega_squared * N /
        (df_effect * (1 - population_omega_squared))
  }

  # Muffle the per-iteration clamp warnings; count and report once.
  .clamp_count <- 0L
  on.exit({
    if (.clamp_count > 0L) {
      message(sprintf(
        "During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in %d intermediate evaluations.",
        .clamp_count))
    }
  }, add = TRUE)

  expected_width_at_N <- function(N) {
    df_error <- N - df_effect - 1
    if (df_error < 3) return(Inf)  # CI machinery is ill-conditioned for tiny df
    F_val <- F_from_N(N)
    res <- tryCatch(
      withCallingHandlers(
        ci_omega_squared(F_value = F_val,
                         df_effect = df_effect,
                         df_error = df_error,
                         N = N,
                         conf_level = conf_level),
        warning = function(w) {
          if (inherits(w, "dmar_ncf_clamp")) {
            .clamp_count <<- .clamp_count + 1L
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) NULL
    )
    if (is.null(res)) return(Inf)
    upper <- res$upper_limit
    lower <- res$lower_limit
    if (length(upper) == 0L || length(lower) == 0L ||
        is.na(upper) || is.na(lower)) return(Inf)
    upper - lower
  }

  # Coarse + fine 1-D search.
  N <- df_effect + 10
  N_max <- 1e6
  width_now <- expected_width_at_N(N)
  iter <- 0L
  while ((is.na(width_now) || width_now > target_width) && N < N_max) {
    N <- ceiling(N * 1.1) + 1
    width_now <- expected_width_at_N(N)
    iter <- iter + 1L
    if (iter > 500L) break
  }
  if (N >= N_max) {
    stop("Required sample size exceeds 1e6; target width may be unattainable.")
  }

  # Refine downward to the smallest N that still satisfies the target.
  for (N_try in seq(N, max(df_effect + 5, N - 100), by = -1L)) {
    if (expected_width_at_N(N_try) > target_width) {
      N <- N_try + 1L
      break
    }
    N <- N_try
  }

  if (!is.null(assurance)) {
    if (assurance <= 0.5 || assurance >= 1)
      stop("'assurance' must be in (0.5, 1).")
    df <- max(1, N - df_effect - 1)
    inflate <- stats::qchisq(assurance, df = df) / df
    N <- ceiling(N * inflate)
  }

  final_width <- expected_width_at_N(N)

  out <- data.frame(
    term  = c("necessary_N", "expected_width",
              "population_omega_squared", "df_effect",
              "width_target", "conf_level"),
    value = c(N, final_width, population_omega_squared, df_effect,
              width, conf_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
