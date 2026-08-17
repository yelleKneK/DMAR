#' Consequences of a Design: Power, Sign and Magnitude Errors, and Expected Precision
#'
#' Evaluates what a design of a given precision will actually deliver, under
#' both of the package's lenses at once. The \emph{significance lens}: the
#' \code{power} of the two-sided test, the \code{type_s_error} (the
#' probability that a statistically significant estimate has the wrong
#' sign), and the \code{exaggeration_ratio} (Type M: the average factor by
#' which significant estimates overstate the true effect), following the
#' design analysis of Gelman and Carlin (2014). The \emph{precision lens},
#' in the accuracy in parameter estimation (AIPE) tradition: the expected
#' half-width and full width of the \code{conf_level} confidence interval
#' the design will produce, the spread of that realized width, and, when a
#' target width \code{w} is supplied, \code{pct_ci_less_w}, the probability
#' that the realized interval is no wider than the target, the same
#' quantities the \code{ss_aipe_*_sensitivity()} family estimates by Monte
#' Carlo, here in closed form.
#'
#' Together they answer the two questions a chosen design should be
#' interrogated with before data collection: \emph{if I run this study and
#' filter it through a significance test, what will the published record
#' look like?} and \emph{how precisely will I estimate the effect
#' regardless of significance?} An underpowered design fails both: its
#' significant estimates are exaggerated and possibly sign-reversed, and
#' its confidence intervals are too wide to be informative.
#'
#' @param true_effect The assumed true (population) effect, on the scale of
#'   the estimate (a mean difference, a regression coefficient, a
#'   standardized mean difference). May be negative. May be \code{NULL},
#'   in which case the significance-lens rows are returned as \code{NA}
#'   and only the precision lens (which does not involve the true effect)
#'   is informative.
#' @param se The standard error of the estimate the design will produce,
#'   on the same scale as \code{true_effect}. Supply \code{se} directly,
#'   or supply \code{sd} with \code{n_1} (and \code{n_2}) and let the
#'   function derive it.
#' @param sd,n_1,n_2 An alternative to \code{se} for the two most common
#'   cases: with \code{sd} and \code{n_1} only, the one-sample (or
#'   paired-difference) design \code{se = sd / sqrt(n_1)} with
#'   \code{df = n_1 - 1}; with \code{n_2} as well, the two-group design
#'   \code{se = sd * sqrt(1/n_1 + 1/n_2)} with \code{df = n_1 + n_2 - 2}
#'   (\code{sd} is the common within-group standard deviation). A
#'   \code{df} supplied explicitly overrides the derived one.
#' @param alpha_level Two-sided Type I error rate of the significance
#'   test. Defaults to 0.05.
#' @param df Degrees of freedom of the reference \emph{t} distribution.
#'   Defaults to \code{Inf} (the normal case) unless derived from
#'   \code{n_1} / \code{n_2}.
#' @param conf_level Confidence level of the interval evaluated by the
#'   precision lens. Defaults to 0.95.
#' @param w Optional target full width for the confidence interval;
#'   when supplied, \code{pct_ci_less_w} reports the probability that the
#'   realized interval is no wider than \code{w}.
#'
#' @details
#' \strong{Significance lens.} Writing \eqn{\lambda = \theta /
#' \mathrm{se}} and \eqn{c} for the two-sided critical value, the power and
#' the Type S error follow from the two tails of the distribution of the
#' test statistic: the noncentral \emph{t} with noncentrality
#' \eqn{\lambda} when \code{df} is finite (the exact distribution of the
#' \emph{t} statistic when the standard error is estimated from the data,
#' the same sampling model the precision lens uses), and the normal when
#' \code{df = Inf}. The exaggeration ratio is the expected absolute
#' estimate conditional on significance over the absolute true effect,
#' computed exactly: from truncated normal moments when \code{df = Inf},
#' and otherwise by integrating those moments over the chi distribution
#' of the estimated standard error, so no simulation error enters.
#' Gelman and Carlin's \code{retrodesign()} instead evaluates a central
#' \emph{t} shifted by \eqn{\lambda} (and simulates the exaggeration
#' ratio under that model), an approximation that treats the standard
#' error as known; the two agree as \code{df} grows and coincide at
#' \code{df = Inf}, but at small \code{df} they differ: in the
#' underpowered regime the design analysis is aimed at (power below
#' about 0.7), the known-se approximation understates power and
#' overstates the Type S and Type M errors, so their published
#' finite-\code{df} values differ from the exact ones reported here. When \code{true_effect = 0} the power equals
#' \code{alpha_level}, the Type S error is 0.5, and the exaggeration ratio
#' is undefined (\code{NA}).
#'
#' \strong{Precision lens.} The realized interval half-width is
#' \eqn{t_{1-\alpha^*/2,\,\mathit{df}} \cdot \widehat{\mathrm{se}}} with
#' \eqn{\alpha^* = 1 - \mathtt{conf\_level}}, and
#' \eqn{\widehat{\mathrm{se}} = \mathrm{se}\sqrt{W/\mathit{df}}} with
#' \eqn{W \sim \chi^2_{\mathit{df}}}, so the width's mean, median, and
#' standard deviation have closed chi-distribution forms and
#' \deqn{P(\mathrm{width} \le w) \;=\;
#'   P\!\left(W \le \mathit{df}\left[\frac{w}
#'   {2\,t\,\mathrm{se}}\right]^{2}\right).}
#' These are the population versions of the \code{mean_ci_width},
#' \code{median_ci_width}, \code{sd_ci_width}, and \code{pct_ci_less_w}
#' terms that the \code{ss_aipe_*_sensitivity()} functions estimate by
#' Monte Carlo. With \code{df = Inf} the standard error is treated as
#' known, the width is deterministic, and \code{pct_ci_less_w} is a step:
#' 1 when the fixed width is at most \code{w} and 0 otherwise.
#'
#' The function complements the planners rather than replacing them:
#' \code{ss_power_*()} chooses a sample size for detection,
#' \code{ss_aipe_*()} chooses one for precision, and
#' \code{design_consequences()} interrogates whatever design came out (or
#' the design a completed study used). The two lenses are the power and
#' accuracy in parameter estimation approaches to sample size planning
#' reviewed by Maxwell, Kelley, and Rausch (2008).
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with the
#'   significance-lens rows (\code{power}, \code{type_s_error},
#'   \code{exaggeration_ratio}), the precision-lens rows
#'   (\code{expected_half_width}, \code{mean_ci_width},
#'   \code{median_ci_width}, \code{sd_ci_width}, \code{pct_ci_less_w},
#'   \code{target_width}; the last two are \code{NA} when no \code{w} is
#'   supplied), and the design rows (\code{true_effect}, \code{se},
#'   \code{df}, \code{alpha_level}). The confidence level is recorded in
#'   the \code{"conf_level"} attribute. The schema is constant: rows that
#'   do not apply are \code{NA}, never dropped.
#'
#' @references
#' Gelman, A., & Carlin, J. (2014). Beyond power calculations: Assessing
#'   Type S (sign) and Type M (magnitude) errors. \emph{Perspectives on
#'   Psychological Science, 9}(6), 641--651. \doi{10.1177/1745691614551642}
#'   (Their accompanying \code{retrodesign()} function evaluates a
#'   location-shifted central \emph{t}, the known-se approximation, and
#'   obtains the exaggeration ratio by simulation; the finite-df case
#'   here uses the exact noncentral \emph{t} and exact moments instead,
#'   so the two differ at small df. See Details.)
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample size planning.
#'   \emph{Evaluation and the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_smd}} and \code{\link{ss_aipe_smd}} for
#'   choosing the sample size by detection or by precision before this
#'   function interrogates the choice; \code{\link{expected_smd}} for the
#'   unconditional small-sample bias of the standardized mean difference,
#'   a different bias than the significance-filter exaggeration here.
#'
#' @family design utilities
#'
#' @keywords design
#'
#' @examples
#' # ---- Both lenses on one underpowered design --------------------------
#' # True effect 0.1 measured with standard error 0.3 (say, d = .1 with
#' # about 22 per group, 45 in total): power is 6 percent, a significant
#' # result has a 17 percent chance of the wrong sign and overstates the
#' # truth seven-fold, and the 95 percent CI is about 1.2 wide, twelve
#' # times the effect. Bad for detection, bad for precision.
#' design_consequences(true_effect = 0.1, se = 0.3)
#'
#' # ---- The same effect, precisely measured -----------------------------
#' design_consequences(true_effect = 0.1, se = 0.03)
#'
#' # ---- From a planned two-group design (sd and per-group n) ------------
#' # d-type effect 0.4, common sd 1, 60 per group; finite df flows through
#' # both lenses.
#' design_consequences(true_effect = 0.4, sd = 1, n_1 = 60, n_2 = 60)
#'
#' # ---- Will the interval beat a target width? --------------------------
#' # Same design, asking for the probability the realized 95 percent CI is
#' # no wider than 0.7 (the closed-form pct_ci_less_w that the
#' # ss_aipe_smd_sensitivity() simulation estimates).
#' design_consequences(true_effect = 0.4, sd = 1, n_1 = 60, n_2 = 60,
#'                     w = 0.7)
#'
#' # ---- Precision lens alone (no effect assumption needed) --------------
#' design_consequences(true_effect = NULL, sd = 1, n_1 = 60, n_2 = 60,
#'                     w = 0.7)
#'
#' # ---- After an AIPE plan: check the detection side --------------------
#' # Plan for a full width of 0.5 on the SMD at delta = .4, then ask what
#' # that design does under the significance filter.
#' n_plan <- ss_aipe_smd(delta = 0.4, width = 0.5)$value[1]
#' design_consequences(true_effect = 0.4, sd = 1,
#'                     n_1 = n_plan, n_2 = n_plan, w = 0.5)
#'
#' @export
#' @importFrom stats dchisq dnorm integrate pchisq pnorm pt qchisq qnorm qt
design_consequences <- function(true_effect = NULL, se = NULL,
                                sd = NULL, n_1 = NULL, n_2 = NULL,
                                alpha_level = 0.05, df = NULL,
                                conf_level = 0.95, w = NULL) {
  # ---- resolve se and df from the supplied design ----
  if (is.null(se)) {
    if (is.null(sd) || is.null(n_1)) {
      stop("Supply 'se' directly, or 'sd' with 'n_1' (and optionally ",
           "'n_2') to derive it.", call. = FALSE)
    }
    if (!is.numeric(sd) || length(sd) != 1L || is.na(sd) || sd <= 0) {
      stop("'sd' must be a single positive number.", call. = FALSE)
    }
    for (nm in c("n_1", "n_2")) {
      val <- get(nm)
      if (!is.null(val) && (!is.numeric(val) || length(val) != 1L ||
                            is.na(val) || val < 2 || val != round(val))) {
        stop(sprintf("'%s' must be a single integer of at least 2.", nm),
             call. = FALSE)
      }
    }
    if (is.null(n_2)) {
      se <- sd / sqrt(n_1)
      if (is.null(df)) df <- n_1 - 1
    } else {
      se <- sd * sqrt(1 / n_1 + 1 / n_2)
      if (is.null(df)) df <- n_1 + n_2 - 2
    }
  } else if (!is.null(sd) || !is.null(n_1) || !is.null(n_2)) {
    stop("Supply either 'se' or the (sd, n_1, n_2) route, not both.",
         call. = FALSE)
  }
  if (is.null(df)) df <- Inf
  if (!is.numeric(se) || length(se) != 1L || is.na(se) || se <= 0) {
    stop("'se' must be a single positive number.", call. = FALSE)
  }
  if (!is.null(true_effect) &&
      (!is.numeric(true_effect) || length(true_effect) != 1L ||
       is.na(true_effect))) {
    stop("'true_effect' must be a single number, or NULL for the ",
         "precision lens alone.", call. = FALSE)
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop("'df' must be a single positive number (possibly Inf).",
         call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(w) && (!is.numeric(w) || length(w) != 1L || is.na(w) ||
                      w <= 0)) {
    stop("'w' must be a single positive target width.", call. = FALSE)
  }

  # ---- significance lens (Gelman & Carlin's design analysis, exact) ----
  crit <- if (is.finite(df)) qt(1 - alpha_level / 2, df) else
    qnorm(1 - alpha_level / 2)

  if (is.null(true_effect)) {
    power <- type_s <- exaggeration <- theta <- NA_real_
  } else {
    theta  <- true_effect
    lambda <- theta / se
    L      <- abs(lambda)
    if (is.finite(df)) {
      # With an estimated standard error the t statistic is exactly
      # noncentral t(df, ncp = lambda), the sampling model the precision
      # lens below also uses. retrodesign()'s location-shifted central t
      # is the known-se approximation; the two coincide at df = Inf.
      p_same <- pt(crit, df, ncp = L, lower.tail = FALSE)
      p_opp  <- pt(-crit, df, ncp = L)
    } else {
      p_same <- 1 - pnorm(crit - L)
      p_opp  <- pnorm(-crit - L)
    }
    power  <- p_same + p_opp
    type_s <- if (L > 0) p_opp / power else 0.5
    if (lambda == 0) {
      exaggeration <- NA_real_
    } else if (!is.finite(df)) {
      c1 <- crit - L
      c2 <- -crit - L
      e_abs <- (L * (1 - pnorm(c1)) + dnorm(c1)) +
               (dnorm(c2) - L * pnorm(c2))
      exaggeration <- e_abs / (power * L)
    } else {
      # E[|estimate|; significant] in se units. Conditional on the
      # estimated standard error, se_hat = se * u with u = sqrt(W / df)
      # and W ~ chi-square_df (density 2 * df * u * dchisq(df * u^2, df)),
      # the estimate is normal and significance means |estimate| >
      # crit * se * u, so the known-se truncated normal moment with
      # critical value crit * u is integrated over the chi distribution.
      e_abs <- integrate(function(u) {
        c1 <- crit * u - L
        c2 <- -crit * u - L
        ((L * (1 - pnorm(c1)) + dnorm(c1)) +
         (dnorm(c2) - L * pnorm(c2))) *
          2 * df * u * dchisq(df * u^2, df)
      }, lower = 0, upper = Inf, rel.tol = 1e-10)$value
      exaggeration <- e_abs / (power * L)
    }
  }

  # ---- precision lens (closed chi-distribution forms) ----
  t_ci <- if (is.finite(df)) qt(1 - (1 - conf_level) / 2, df) else
    qnorm(1 - (1 - conf_level) / 2)
  fixed_width <- 2 * t_ci * se
  if (is.finite(df)) {
    # E[sqrt(W/df)] for W ~ chi-square_df: the chi mean factor.
    c_df <- sqrt(2 / df) * exp(lgamma((df + 1) / 2) - lgamma(df / 2))
    mean_w   <- fixed_width * c_df
    median_w <- fixed_width * sqrt(qchisq(0.5, df) / df)
    sd_w     <- fixed_width * sqrt(max(0, 1 - c_df^2))
    pct_w <- if (is.null(w)) NA_real_ else
      pchisq(df * (w / fixed_width)^2, df)
  } else {
    mean_w <- median_w <- fixed_width
    sd_w   <- 0
    pct_w  <- if (is.null(w)) NA_real_ else as.numeric(fixed_width <= w)
  }

  .as_dmar_tbl(data.frame(
    term  = c("power", "type_s_error", "exaggeration_ratio",
              "expected_half_width", "mean_ci_width", "median_ci_width",
              "sd_ci_width", "pct_ci_less_w", "target_width",
              "true_effect", "se", "df", "alpha_level"),
    value = c(power, type_s, exaggeration,
              mean_w / 2, mean_w, median_w,
              sd_w, pct_w, if (is.null(w)) NA_real_ else w,
              theta, se, df, alpha_level),
    stringsAsFactors = FALSE
  ), conf_level = conf_level)
}
