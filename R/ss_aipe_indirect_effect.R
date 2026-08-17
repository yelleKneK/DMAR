# Sample size for AIPE on a mediated (indirect) effect ab.
#' Sample Size for AIPE on a Mediated (Indirect) Effect \eqn{ab}
#'
#' Determines the sample size needed for the confidence interval on a
#' mediated effect \eqn{ab} (the product of the \eqn{X \to M} and
#' \eqn{M \to Y} coefficients in a simple three-variable mediation model)
#' to have a desired full width. Two methods are available.
#' \code{"closed_form"} (the default) plans for the symmetric Wald
#' interval built on the delta method standard error of the product
#' (Sobel, 1982) and answers instantly. \code{"monte_carlo"} plans for
#' the Monte Carlo confidence interval (MacKinnon, Lockwood, & Williams,
#' 2004; Tofighi & MacKinnon, 2011), the interval
#' \code{\link{mediation_mbco}} reports under
#' \code{ci_method = "monte_carlo"}, and it plans by a priori Monte
#' Carlo simulation: at each candidate sample size the mediation model
#' is fit to \code{G} simulated data sets and the realized interval
#' widths are recorded. Because the Monte Carlo interval respects the
#' skewness of the sampling distribution of a product, and because the
#' simulation measures the widths that fitted models actually deliver,
#' \code{method = "monte_carlo"} is the recommended way to settle on the
#' final sample size; the closed form is its fast first approximation.
#'
#' @param a Anticipated population coefficient for \eqn{X \to M},
#'   on the standardized scale. Numeric scalar in \eqn{(-1, 1)}.
#' @param b Anticipated population coefficient for \eqn{M \to Y}
#'   controlling for \eqn{X}, standardized. Numeric scalar in
#'   \eqn{(-1, 1)}.
#' @param width Desired full width of the confidence interval on
#'   \eqn{ab}.
#' @param method One of \code{"closed_form"} (default) or
#'   \code{"monte_carlo"}; see Details.
#' @param conf_level Desired confidence level (default \code{0.95}).
#' @param n_max Upper bound on the search; default \code{10000}.
#' @param B Number of Monte Carlo draws forming the interval within each
#'   simulated study when \code{method = "monte_carlo"}; default
#'   \code{5000}. This is the same \code{B} the analysis-stage interval
#'   uses (see \code{\link{mediation_mbco}}).
#' @param G Number of simulated studies per candidate sample size when
#'   \code{method = "monte_carlo"}; default \code{1000}. The mean
#'   simulated width carries a simulation error of about its standard
#'   deviation over \eqn{\sqrt{G}}; raise \code{G} for a sharper answer.
#' @param seed Optional integer seed for the Monte Carlo method, used
#'   locally (the caller's random number generator state is restored on
#'   exit). Default \code{NULL} leaves the random number generator state
#'   alone.
#'
#' @return A \code{data.frame} with rows for the recommended
#'   sample size, the expected CI width at that size, and the inputs
#'   echoed back. Under \code{method = "closed_form"} the expected width
#'   is the delta method width evaluated at the returned sample size;
#'   under \code{method = "monte_carlo"} it is the mean simulated width
#'   there. The method is carried on the returned object as the
#'   \code{ci_method} attribute.
#'
#' @details
#' \strong{The mediation model.} The simple mediator model is
#' \deqn{M = \alpha_1 + a X + \varepsilon_M,}
#' \deqn{Y = \alpha_2 + c' X + b M + \varepsilon_Y,}
#' with the indirect (mediated) effect of \eqn{X} on \eqn{Y} through
#' \eqn{M} equal to \eqn{ab} (MacKinnon, Lockwood, Hoffman, West, &
#' Sheets, 2002). Both methods plan on the standardized scale with no
#' direct effect: the planning population takes \eqn{X}, \eqn{M}, and
#' \eqn{Y} with unit variances and \eqn{c' = 0}, so \code{a} and
#' \code{b} are the standardized paths.
#'
#' \strong{The closed form.} Under the planning population the sampling
#' variance of \eqn{\hat a} is \eqn{(1 - a^2)/(n - 2)}. In the equation
#' for \eqn{Y} the mediator is regressed alongside \eqn{X}, with which
#' it is correlated at \eqn{a}, so the sampling variance of \eqn{\hat b}
#' carries the variance inflation factor \eqn{1/(1 - a^2)}:
#' \deqn{\mathrm{Var}(\hat b) \;=\;
#'         \frac{1 - b^2}{(n - 3)(1 - a^2)}.}
#' The estimators come from two separate equations and are uncorrelated,
#' so the delta method (Sobel, 1982) standard error of the product is
#' \deqn{\mathrm{SE}(\hat a \hat b) \;=\;
#'         \sqrt{\,a^2 \mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a)\,},}
#' and the closed form returns the smallest \eqn{n} at which the Wald
#' width \eqn{2 z_{1 - \alpha/2}\, \mathrm{SE}(\hat a \hat b)} is at or
#' below \code{width}. Two approximations remain. The Wald interval is
#' symmetric while the sampling distribution of a product is skewed, so
#' the Wald interval is not the interval an indirect effect should be
#' reported with (Tofighi & Kelley, 2020). And the closed form evaluates
#' the standard error at the planning values, while a fitted model
#' evaluates it at the estimates, which leaves a discrepancy of a
#' percent or two in realized width at moderate sample sizes. Both are
#' reasons to treat the closed form as the first approximation and to
#' verify the final plan with \code{method = "monte_carlo"}, which
#' measures the realized widths directly.
#'
#' \strong{Planning for the Monte Carlo interval.} With
#' \code{method = "monte_carlo"}, each candidate \eqn{n} is evaluated by
#' a priori Monte Carlo simulation (Muthén & Muthén, 2002; Schoemann,
#' Boulton, & Short, 2017): \code{G} data sets of size \eqn{n} are drawn
#' from the planning population, the two mediation regressions are fit
#' to each, and the Monte Carlo interval is formed by drawing \code{B}
#' pairs \eqn{(\tilde a, \tilde b)} from normal distributions centered
#' at the estimates with the estimated standard errors, multiplying, and
#' reading off the empirical \eqn{(\alpha/2, 1 - \alpha/2)} quantiles
#' (MacKinnon, Lockwood, & Williams, 2004). Since \eqn{\hat a} and
#' \eqn{\hat b} are uncorrelated here, the independent draws realize the
#' joint normal approximation of the estimates, the same construction
#' \code{\link{mediation_mbco}} uses for its Monte Carlo interval. The
#' necessary sample size is the smallest \eqn{n} whose mean simulated
#' width is at or below \code{width}; the search starts from the
#' closed-form answer, brackets the crossing geometrically, and bisects.
#' A planning call fits the mediation model several thousand times and
#' takes a few seconds, which is why the Monte Carlo example below is
#' shown rather than run. The necessary sample size inherits the
#' simulation error of the mean widths; raising \code{G} narrows it, and
#' supplying \code{seed} makes a plan reproducible.
#'
#' \strong{Relation to the MBCO procedure.} The model-based constrained
#' optimization (MBCO) likelihood ratio test of Tofighi and Kelley
#' (2020), implemented in \code{\link{mediation_mbco}}, is the
#' recommended test of a mediation effect, and the intervals that suit
#' an indirect effect are the profile likelihood interval and the Monte
#' Carlo interval, both of which accommodate the skewness of the
#' product. This planner targets the Monte Carlo interval. Planning for
#' the profile likelihood interval would require inverting a pair of
#' constrained optimizations in every simulated study (two constrained
#' \pkg{OpenMx} fits per interval, times \code{G}, times every candidate
#' sample size), while the Monte Carlo interval costs \code{B} products
#' of normal draws per study and is the inexpensive interval that also
#' accommodates the skewness, the one Tofighi and Kelley (2020) report
#' for their memory example. A study planned with
#' \code{method = "monte_carlo"} and analyzed with
#' \code{mediation_mbco(ci_method = "monte_carlo")} is therefore planned
#' and analyzed on the same interval.
#'
#' \strong{Beyond the simple model.} The planning population assumes
#' standardized observed variables, one mediator, no covariates, and no
#' direct effect. With a nonzero direct effect the residual variance of
#' \eqn{Y} is \eqn{1 - b^2 - c'^2 - 2abc'} rather than \eqn{1 - b^2}, so
#' assuming \eqn{c' = 0} errs toward a larger sample whenever
#' \eqn{c'(c' + 2ab) > 0} (consistent mediation) and toward a smaller
#' one otherwise. When the direct effect, covariates, several mediators,
#' or latent variables matter to the design, plan by simulation from the
#' full model with \code{\link{ss_aipe_composite_sem}}, labeling the
#' paths and defining the indirect effect via \code{ab := a*b}; its
#' intervals are the Wald intervals of the fitted model, the same target
#' as the closed form here. \code{\link{ss_aipe_indirect_effect_sensitivity}}
#' quantifies what a plan from this page delivers when the population
#' paths differ from the planning values.
#'
#' @references
#' Fritz, M. S., & MacKinnon, D. P. (2007). Required sample size to
#'   detect the mediated effect. \emph{Psychological Science, 18}(3),
#'   233--239. \doi{10.1111/j.1467-9280.2007.01882.x}
#'
#' Lachowicz, M. J., Preacher, K. J., & Kelley, K. (2018). A novel measure
#'   of effect size for mediation analysis.
#'   \emph{Psychological Methods, 23}, 244--261. \doi{10.1037/met0000165}
#'
#' MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
#'   Sheets, V. (2002). A comparison of methods to test mediation and
#'   other intervening variable effects. \emph{Psychological Methods,
#'   7}(1), 83--104. \doi{10.1037/1082-989X.7.1.83}
#'
#' MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
#'   limits for the indirect effect: Distribution of the product and
#'   resampling methods. \emph{Multivariate Behavioral Research, 39}(1),
#'   99--128. \doi{10.1207/s15327906mbr3901_4}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Muthén, L. K., & Muthén, B. O. (2002). How to use a
#'   Monte Carlo study to decide on sample size and determine power.
#'   \emph{Structural Equation Modeling, 9}(4), 599--620.
#'   \doi{10.1207/S15328007SEM0904_8}
#'
#' Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
#'   models: Quantitative strategies for communicating indirect effects.
#'   \emph{Psychological Methods, 16}(2), 93--115. \doi{10.1037/a0022658}
#'
#' Schoemann, A. M., Boulton, A. J., & Short, S. D. (2017).
#'   Determining power and sample size for simple and complex
#'   mediation models. \emph{Social Psychological and Personality
#'   Science, 8}(4), 379--386. \doi{10.1177/1948550617715068}
#'
#' Sobel, M. E. (1982). Asymptotic confidence intervals for indirect
#'   effects in structural equation models. \emph{Sociological
#'   Methodology, 13}, 290--312.
#'
#' Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
#'   analysis: Introducing the model-based constrained optimization
#'   procedure. \emph{Psychological Methods, 25}, 496--515.
#'   \doi{10.1037/met0000259}
#'
#' Tofighi, D., & MacKinnon, D. P. (2011). RMediation: An R package for
#'   mediation analysis confidence intervals. \emph{Behavior Research
#'   Methods, 43}(3), 692--700. \doi{10.3758/s13428-011-0076-x}
#'
#' @seealso \code{\link{mediation_mbco}} for the analysis the plan
#'   feeds; \code{\link{ss_aipe_composite_sem}} for AIPE planning of an
#'   indirect effect in an arbitrary lavaan model;
#'   \code{\link{ss_aipe_indirect_effect_sensitivity}};
#'   \code{\link{ss_aipe_partial_r}},
#'   \code{\link{ss_aipe_semipartial_r}}, \code{\link{ss_aipe_rc}}
#'
#' @examples
#' # 1. Plan n so the 95% CI on ab has full width <= 0.20, with
#' #        anticipated standardized a = 0.40 and b = 0.40. The closed
#' #        form answers instantly:
#' ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20)
#'
#' # 2. The recommended plan targets the Monte Carlo interval directly:
#' #        every candidate sample size fits the mediation model to G
#' #        simulated data sets and measures the realized widths. The
#' #        call takes a few seconds, so it is shown here rather than
#' #        run. It returns a slightly larger sample size than the
#' #        closed form because the interval it plans for is a little
#' #        wider than the Wald interval:
#' # ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20,
#' #                         method = "monte_carlo", seed = 113)
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

ss_aipe_indirect_effect <- function(a, b, width,
                                    method = c("closed_form",
                                               "monte_carlo"),
                                    conf_level = 0.95,
                                    n_max = 10000L,
                                    B = 5000L,
                                    G = 1000L,
                                    seed = NULL) {
  method <- match.arg(method)
  if (!is.numeric(a) || length(a) != 1L || abs(a) >= 1)
    stop("'a' must be a single value in (-1, 1).")
  if (!is.numeric(b) || length(b) != 1L || abs(b) >= 1)
    stop("'b' must be a single value in (-1, 1).")
  if (!is.numeric(width) || length(width) != 1L || width <= 0)
    stop("'width' must be a single positive number.")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")
  n_min <- 10L
  if (!is.numeric(n_max) || length(n_max) != 1L || is.na(n_max) ||
      n_max != round(n_max) || n_max <= n_min)
    stop("'n_max' must be a single whole number greater than ", n_min, ".")
  n_max <- as.integer(n_max)
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 100 ||
      B != round(B))
    stop("'B' must be a single whole number of at least 100.")
  B <- as.integer(B)
  if (!is.numeric(G) || length(G) != 1L || is.na(G) || G < 10 ||
      G != round(G))
    stop("'G' must be a single whole number of at least 10.")
  G <- as.integer(G)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
      stop("'seed' must be NULL or a single number.")
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(),
                      inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = globalenv(),
                         inherits = FALSE))
        rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  ab_hat  <- a * b

  # Delta method width of the Wald interval at the planning values. The
  # b-hat variance carries the 1/(1 - a^2) variance inflation factor for
  # the correlation between M and X in the Y equation.
  closed_form_width_at_n <- function(n) {
    var_a <- (1 - a^2) / (n - 2)
    var_b <- (1 - b^2) / ((n - 3) * (1 - a^2))
    2 * z_alpha * sqrt(a^2 * var_b + b^2 * var_a)
  }

  # Smallest n at which the closed-form width meets the target; the
  # width is strictly decreasing in n, so bracket by doubling and
  # bisect to adjacent integers.
  closed_form_n <- function() {
    if (closed_form_width_at_n(n_min) <= width) return(n_min)
    lo <- n_min
    hi <- n_min
    repeat {
      hi <- min(n_max, 2L * hi)
      if (closed_form_width_at_n(hi) <= width) break
      if (hi >= n_max)
        stop("Required sample size exceeds n_max = ", n_max,
             "; target width may be unattainable.")
      lo <- hi
    }
    while (hi - lo > 1L) {
      mid <- lo + (hi - lo) %/% 2L
      if (closed_form_width_at_n(mid) <= width) hi <- mid else lo <- mid
    }
    hi
  }

  # Mean simulated width of the Monte Carlo interval at candidate n:
  # G data sets from the planning population, both regressions fit by
  # their normal equations (no lm() allocation in the inner loop), and
  # each study's interval from B products of normal draws.
  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  mc_mean_width_at_n <- function(n) {
    widths <- numeric(G)
    sd_m <- sqrt(1 - a^2)
    sd_y <- sqrt(1 - b^2)
    for (g in seq_len(G)) {
      X <- stats::rnorm(n)
      M <- a * X + stats::rnorm(n, sd = sd_m)
      Y <- b * M + stats::rnorm(n, sd = sd_y)
      x <- X - mean(X)
      m <- M - mean(M)
      y <- Y - mean(Y)
      sxx <- sum(x * x); sxm <- sum(x * m); smm <- sum(m * m)
      smy <- sum(m * y); sxy <- sum(x * y); syy <- sum(y * y)
      a_hat <- sxm / sxx
      se_a  <- sqrt((smm - a_hat * sxm) / (n - 2) / sxx)
      det   <- smm * sxx - sxm * sxm
      b_hat <- (smy * sxx - sxm * sxy) / det
      c_hat <- (smm * sxy - sxm * smy) / det
      se_b  <- sqrt((syy - b_hat * smy - c_hat * sxy) / (n - 3) *
                      sxx / det)
      qs <- stats::quantile(stats::rnorm(B, a_hat, se_a) *
                              stats::rnorm(B, b_hat, se_b),
                            probs = probs, names = FALSE)
      widths[g] <- qs[2L] - qs[1L]
    }
    mean(widths)
  }

  if (method == "closed_form") {
    n <- closed_form_n()
    expected_width <- closed_form_width_at_n(n)
  } else {
    # Start at the closed-form answer, bracket the crossing of the mean
    # simulated width geometrically, and bisect to adjacent integers.
    n0 <- closed_form_n()
    w0 <- mc_mean_width_at_n(n0)
    if (w0 <= width) {
      hi <- n0
      w_hi <- w0
      lo <- n_min - 1L
      while (hi > n_min) {
        cand <- max(n_min, min(as.integer(hi / 1.25), hi - 1L))
        w_c <- mc_mean_width_at_n(cand)
        if (w_c <= width) {
          hi <- cand
          w_hi <- w_c
        } else {
          lo <- cand
          break
        }
      }
    } else {
      lo <- n0
      repeat {
        if (lo >= n_max)
          stop("Required sample size exceeds n_max = ", n_max,
               "; target width may be unattainable.")
        cand <- min(n_max, as.integer(ceiling(lo * 1.25)) + 1L)
        w_c <- mc_mean_width_at_n(cand)
        if (w_c <= width) {
          hi <- cand
          w_hi <- w_c
          break
        }
        lo <- cand
      }
    }
    while (hi - lo > 1L) {
      mid <- lo + (hi - lo) %/% 2L
      w_c <- mc_mean_width_at_n(mid)
      if (w_c <= width) {
        hi <- mid
        w_hi <- w_c
      } else {
        lo <- mid
      }
    }
    n <- hi
    expected_width <- w_hi
  }

  out <- data.frame(
    term  = c("necessary_N", "expected_width", "a", "b", "ab",
              "width_target", "conf_level"),
    value = c(n, expected_width, a, b, ab_hat, width, conf_level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "ci_method") <- method
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
