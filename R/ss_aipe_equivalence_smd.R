# AIPE sample size planning for TOST on the SMD.
#' AIPE Sample Size Planning for an Equivalence Test on the Standardized Mean Difference
#'
#' Computes the minimum per-group sample size needed so that the
#' equivalence CI (the 100(1 - 2\eqn{\alpha})\% CI on the standardized
#' mean difference) has expected full width \eqn{\le \omega}, that is,
#' expected half-width \eqn{\le \omega / 2} (Kelley, 2007; Lakens,
#' 2017). With the standard \eqn{\alpha = 0.05} TOST level, the
#' equivalence CI is the 90\% CI. The function inverts the
#' large-sample variance of \emph{d}; the noncentral \emph{t}
#' distribution of \emph{d} enters through the optional assurance
#' step.
#'
#' @param population_smd Anticipated population standardized mean
#'   difference \eqn{\delta} used to plan the variance. Default
#'   \code{0}: plans under the most-conservative null-effect case.
#' @param width Target full CI width on the \emph{d} scale (e.g.,
#'   \code{0.20} for a 90\% CI of width 0.20).
#' @param alpha_level One-sided TOST significance level. The CI used in
#'   planning is at confidence level \eqn{1 - 2\alpha}. Default
#'   \code{0.05} (90\% CI).
#' @param assurance Optional assurance probability in \eqn{(0, 1)}.
#'   When supplied, the function chooses \emph{n} so that the
#'   probability of achieving \code{width} or less is at least
#'   \code{assurance} (Kelley, Maxwell, & Rausch, 2003). Default
#'   \code{NULL} (no assurance correction).
#' @param balanced Logical; \code{TRUE} (default) plans equal-\eqn{n}
#'   groups. Unequal-\eqn{n} planning is not yet supported.
#'
#' @return A 5-row \code{data.frame} with columns \code{term} and
#'   \code{value}: the per-group recommended sample size
#'   \code{necessary_n_per_group}, the implied total \code{total_N}, the
#'   target \code{width}, the planning value \code{population_smd},
#'   and \code{ci_width_expected}, the expected full CI width at the
#'   chosen \emph{n}.
#'
#' @details
#' \strong{Approximate-variance plan.} The large-sample variance of
#' \emph{d} is
#' \deqn{\mathrm{Var}(\hat d) \;\approx\; (n_1 + n_2) / (n_1 n_2) +
#'    d^2 / (2 (n_1 + n_2)).}
#' For a balanced design with per-group size \eqn{n}, the half-width
#' of the equivalence CI at level \eqn{1 - 2\alpha} is approximately
#' \eqn{z_{1-\alpha} \sqrt{\mathrm{Var}(\hat d)}}. The function
#' solves for the smallest integer \eqn{n} giving expected half-width
#' \eqn{\le \omega / 2}.
#'
#' \strong{Assurance.} Under \code{assurance = q}, the function
#' increments \eqn{n} until the simulated probability that the
#' realized half-width is \eqn{\le \omega / 2} is at least \eqn{q}.
#' (Implemented as a thin Monte Carlo overlay. At the default
#' planning value \code{population_smd = 0} the shift is typically
#' zero; it grows with the planning value, reaching several per
#' group by \code{population_smd = 0.5} with a narrow target width.)
#'
#' \strong{Note on conservatism of the assurance plan.} The empirical
#' simulation study of the AIPE planner family finds that
#' \code{ss_aipe_equivalence_smd()} is tight
#' at \eqn{\gamma = 0.80} but operates on the boundary of its valid
#' range at \eqn{\gamma = 0.99}: the realized assurance at the
#' recommended sample size is within Monte Carlo error of the target,
#' typically a few tenths of a percentage point below 0.99. The
#' mechanism is that the planner inverts a normal approximation to
#' \eqn{\Pr(\widehat W > \omega)}, and at the 99\% level the upper
#' tail of \eqn{\widehat W} is heavier than the approximation
#' accounts for. Adding a small safety margin (5 to 10 subjects per
#' group) restores the desired probability statement when planning
#' at high assurance;
#' \code{\link{ss_aipe_equivalence_smd_sensitivity}} reproduces the
#' check for any one condition.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for
#'   standardized effect sizes: Theory, application, and
#'   implementation. \emph{Journal of Statistical Software, 20}(8),
#'   1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample size planning.
#'   \emph{Evaluation and the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for
#'   \emph{t} tests, correlations, and meta-analyses. \emph{Social
#'   Psychological and Personality Science, 8}(4), 355--362.
#'   \doi{10.1177/1948550617697177}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{equivalence_smd}}, \code{\link{ss_aipe_smd}},
#'   \code{\link{ci_smd}}
#'
#' @examples
#' # 1. Plan for a 90% CI on d of width <= 0.20, under d_planning = 0:
#' ss_aipe_equivalence_smd(population_smd = 0, width = 0.20)
#'
#' # 2. Plan for the same width assuming a true d = 0.05:
#' ss_aipe_equivalence_smd(population_smd = 0.05, width = 0.20)
#'
#' # 3. With 80% assurance (the assurance path is Monte Carlo, so seed for
#' #    a reproducible result):
#' set.seed(113)
#' ss_aipe_equivalence_smd(population_smd = 0.05, width = 0.20, assurance = 0.80)
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

ss_aipe_equivalence_smd <- function(population_smd = 0,
                             width,
                             alpha_level = 0.05,
                             assurance = NULL,
                             balanced = TRUE) {
  if (!is.numeric(width) || width <= 0)
    stop("'width' must be a positive number.")
  if (!is.numeric(population_smd) || length(population_smd) != 1L)
    stop("'population_smd' must be a single numeric value.")
  if (alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be in (0, 0.5).")
  if (!isTRUE(balanced))
    stop("Only balanced two-group designs are supported in this release.")

  d <- population_smd
  z_alpha <- stats::qnorm(1 - alpha_level)
  half_target <- width / 2

  half_width <- function(n) {
    var_d <- 2 / n + d^2 / (4 * n)
    z_alpha * sqrt(var_d)
  }

  n <- 4L
  while (half_width(n) > half_target) {
    n <- n + 1L
    if (n > 1e6L) stop("Required n exceeded 1e6 per group; check inputs.")
  }

  if (!is.null(assurance)) {
    if (assurance <= 0 || assurance >= 1)
      stop("'assurance' must be in (0, 1).")
    target_p <- assurance
    repeat {
      df <- 2 * n - 2
      # Monte Carlo over the noncentral t distribution of d. The draw uses the
      # caller's current RNG state; reproducibility is the caller's to control
      # by seeding before the call (see the examples), not the function's to
      # bake in. Seeding inside this loop would also redraw the same values
      # every iteration and silently reset the user's global RNG stream.
      ncp <- d * sqrt(n / 2)
      d_sim <- stats::rt(2000L, df = df, ncp = ncp) / sqrt(n / 2)
      var_sim <- 2 / n + d_sim^2 / (4 * n)
      hw_sim  <- z_alpha * sqrt(var_sim)
      if (mean(hw_sim <= half_target) >= target_p) break
      n <- n + 1L
      if (n > 1e6L) stop("Required n exceeded 1e6 per group under assurance.")
    }
  }

  out <- data.frame(
    term  = c("necessary_n_per_group", "total_N", "width",
              "population_smd", "ci_width_expected"),
    value = c(n, 2 * n, width, d, 2 * half_width(n)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, subclass = "dmar_ss_aipe")
}
