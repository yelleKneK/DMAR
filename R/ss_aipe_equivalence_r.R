# AIPE sample size planning for TOST on the Pearson correlation.
#' AIPE Sample Size Planning for an Equivalence Test on the Pearson Correlation
#'
#' Computes the minimum sample size needed so that the equivalence CI
#' (the 100(1 - 2\eqn{\alpha})\% CI on the Pearson correlation
#' \eqn{\rho}) has expected full width \eqn{\le \omega} (Kelley, 2007;
#' Lakens, 2017). With the standard \eqn{\alpha = 0.05} TOST level,
#' the equivalence CI is the 90\% CI. The interval is the Fisher's
#' \eqn{Z} construction that \code{\link{equivalence_r}} and
#' \code{\link{ci_r}} use, so the plan and the analysis invert the
#' same interval.
#'
#' @param population_r Anticipated population correlation \eqn{\rho}
#'   used to plan the width. Default \code{0}: the confidence interval
#'   for a correlation is widest at \eqn{\rho = 0}, so the default
#'   plans under the widest-interval case and is conservative for any
#'   other value.
#' @param width Target full CI width on the correlation scale (e.g.,
#'   \code{0.20} for a 90\% CI of width 0.20). Must be in
#'   \eqn{(0, 2)}, the width of the correlation scale itself.
#' @param alpha_level One-sided TOST significance level. The CI used in
#'   planning is at confidence level \eqn{1 - 2\alpha}. Default
#'   \code{0.05} (90\% CI).
#' @param assurance Optional assurance probability in \eqn{(0, 1)}.
#'   When supplied, the function chooses \emph{N} so that the
#'   probability of achieving \code{width} or less is at least
#'   \code{assurance} (Kelley, Maxwell, & Rausch, 2003). Default
#'   \code{NULL} (no assurance correction).
#'
#' @return A 4-row \code{data.frame} with columns \code{term} and
#'   \code{value}: the recommended sample size \code{necessary_N}, the
#'   target \code{width}, the planning value \code{population_r}, and
#'   the resulting \code{ci_width_expected} at the chosen \emph{N}.
#'
#' @details
#' \strong{Closed form on the Fisher's \eqn{Z} scale.} The equivalence
#' CI has half-width \eqn{h = z_{1-\alpha} / \sqrt{N - 3}} on the
#' Fisher's \eqn{Z} scale, and its width on the correlation scale is
#' \deqn{w(N) \;=\; \tanh(Z_\rho + h) - \tanh(Z_\rho - h),}
#' where \eqn{Z_\rho = \tanh^{-1}(\rho)}. The function returns the
#' smallest integer \eqn{N \ge 4} with \eqn{w(N) \le \omega}. At
#' \eqn{\rho = 0} this is available in closed form,
#' \eqn{N = \lceil 3 + (z_{1-\alpha} / \tanh^{-1}(\omega / 2))^2
#' \rceil, } and away from zero the back-transform shortens the
#' interval, so the required \emph{N} can only decrease as
#' \eqn{|\rho|} grows.
#'
#' \strong{Choosing the width from equivalence bounds.} To leave room
#' for an equivalence verdict inside bounds \eqn{(-b, b)}, the
#' interval must at minimum fit inside the bounds when centered at the
#' anticipated \eqn{\rho}, so a width somewhat below \eqn{2 b} (for
#' \eqn{\rho} near 0) is the natural target; the Monte Carlo
#' sensitivity sibling \code{\link{ss_aipe_equivalence_r_sensitivity}}
#' reports the realized proportion of equivalence verdicts at the
#' planned \emph{N}.
#'
#' \strong{Assurance.} Under \code{assurance = q}, the function
#' increments \eqn{N} until the Monte Carlo probability that the
#' realized width is \eqn{\le \omega} is at least \eqn{q}, drawing
#' the sampling distribution of \eqn{\widehat Z} as normal with mean
#' \eqn{Z_\rho} and variance \eqn{1 / (N - 3)}.
#'
#' @references
#' Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for
#'   comparing correlation and regression coefficients. \emph{British
#'   Journal of Mathematical and Statistical Psychology, 68}(2),
#'   292--309. \doi{10.1111/bmsp.12045}
#'
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of
#'   association: An equivalence testing approach. \emph{British
#'   Journal of Mathematical and Statistical Psychology, 63}(3),
#'   527--537. \doi{10.1348/000711009X475853}
#'
#' Kelley, K. (2007). Confidence intervals for
#'   standardized effect sizes: Theory, application, and
#'   implementation. \emph{Journal of Statistical Software, 20}(8),
#'   1--24. \doi{10.18637/jss.v020.i08}
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
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{equivalence_r}}, \code{\link{ss_aipe_r}},
#'   \code{\link{ci_r}}, \code{\link{ss_aipe_equivalence_smd}},
#'   \code{\link{ss_aipe_equivalence_r_sensitivity}}
#'
#' @examples
#' # 1. Plan for a 90% CI on the correlation of width <= 0.20, under
#' #    the widest-interval planning value rho = 0:
#' ss_aipe_equivalence_r(population_r = 0, width = 0.20)
#'
#' # 2. The same width assuming a true correlation of 0.30 requires
#' #    fewer participants, since the interval narrows away from zero:
#' ss_aipe_equivalence_r(population_r = 0.30, width = 0.20)
#'
#' # 3. With 80% assurance (the assurance path is Monte Carlo, so seed
#' #    for a reproducible result):
#' set.seed(113)
#' ss_aipe_equivalence_r(population_r = 0.30, width = 0.20, assurance = 0.80)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family AIPE sample size planning
#'
#' @export

ss_aipe_equivalence_r <- function(population_r = 0,
                           width,
                           alpha_level = 0.05,
                           assurance = NULL) {
  if (!is.numeric(width) || length(width) != 1L || width <= 0 || width >= 2)
    stop("'width' must be a single number in (0, 2).", call. = FALSE)
  if (!is.numeric(population_r) || length(population_r) != 1L ||
      abs(population_r) >= 1)
    stop("'population_r' must be a single correlation with |r| < 1.",
         call. = FALSE)
  if (alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be in (0, 0.5).", call. = FALSE)

  Z_rho <- atanh(population_r)
  z_alpha <- stats::qnorm(1 - alpha_level)

  width_r <- function(N) {
    h <- z_alpha / sqrt(N - 3)
    tanh(Z_rho + h) - tanh(Z_rho - h)
  }

  N <- 4L
  while (width_r(N) > width) {
    N <- N + 1L
    if (N > 1e6L) stop("Required N exceeded 1e6; check inputs.", call. = FALSE)
  }

  if (!is.null(assurance)) {
    if (assurance <= 0 || assurance >= 1)
      stop("'assurance' must be in (0, 1).", call. = FALSE)
    repeat {
      # Monte Carlo over the sampling distribution of Fisher's Z. The draw
      # uses the caller's current RNG state; reproducibility is the caller's
      # to control by seeding before the call (see the examples), not the
      # function's to bake in. Seeding inside this loop would also redraw
      # the same values every iteration and silently reset the user's
      # global RNG stream.
      h <- z_alpha / sqrt(N - 3)
      Z_sim <- stats::rnorm(2000L, mean = Z_rho, sd = 1 / sqrt(N - 3))
      w_sim <- tanh(Z_sim + h) - tanh(Z_sim - h)
      if (mean(w_sim <= width) >= assurance) break
      N <- N + 1L
      if (N > 1e6L)
        stop("Required N exceeded 1e6 under assurance.", call. = FALSE)
    }
  }

  out <- data.frame(
    term  = c("necessary_N", "width", "population_r", "ci_width_expected"),
    value = c(N, width, population_r, width_r(N)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, subclass = "dmar_ss_aipe")
}
