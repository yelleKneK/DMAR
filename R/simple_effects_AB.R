# Simple effect F tests for a two-factor between-subjects design.
#' Simple Effect F Tests for a Two-Factor Between-Subjects Design
#'
#' Given a fitted \code{\link[stats]{aov}} or \code{\link[stats]{lm}} object
#' for a two-factor between-subjects design, conventionally written
#' \eqn{Y \sim A \times B}, where \eqn{A} and \eqn{B} are crossed fixed
#' factors, computes the family of \emph{simple main effects}: the
#' effect of \eqn{A} at each level of \eqn{B} and/or the effect of \eqn{B}
#' at each level of \eqn{A}. Each row carries the simple effect \emph{F}
#' test, its (optionally adjusted) \emph{p}-value, and the partial
#' \eqn{\eta^2} with a noncentrality-based confidence interval. The error
#' term can be either the full-model pooled \eqn{\mathit{MS}_W} (the
#' textbook default) or a Welch--Satterthwaite test refitted within each
#' conditioning level (robust to within-level heteroscedasticity).
#'
#' @param object A fitted \code{\link[stats]{aov}} or \code{\link[stats]{lm}}
#'   object whose right-hand side has \emph{exactly two} crossed factors
#'   (e.g.\ \code{iq_gain ~ treatment * grade}). The interaction term is
#'   strongly recommended so that the pooled error is the pure within-cell
#'   \eqn{\mathit{MS}_W}; the function still runs without it but issues a
#'   warning (the additive-model residual includes interaction variance and
#'   inflates the error term used for the pooled simple effect \emph{F}).
#' @param which Which family of simple effects to report:
#'   \describe{
#'     \item{\code{"both"} (default)}{The \eqn{b} tests of \eqn{A} at each
#'       level of \eqn{B} followed by the \eqn{a} tests of \eqn{B} at each
#'       level of \eqn{A} (\eqn{a + b} rows).}
#'     \item{\code{"A_at_B"}}{Only the \eqn{b} tests of the first factor at
#'       each level of the second.}
#'     \item{\code{"B_at_A"}}{Only the \eqn{a} tests of the second factor at
#'       each level of the first.}
#'   }
#'   The first factor on the right-hand side of the model formula is treated
#'   as \eqn{A}; the second as \eqn{B}.
#' @param error_term Error-term strategy for the simple effect \emph{F}:
#'   \describe{
#'     \item{\code{"pooled"} (default)}{Use the full factorial model's
#'       \eqn{\mathit{MS}_W} and its residual \emph{df}. This is Maxwell,
#'       Delaney, and Kelley's preferred default and gives the simple effect
#'       \emph{F} maximum denominator \emph{df}. Validity rests on
#'       homogeneity of variance across all \eqn{a \times b} cells.}
#'     \item{\code{"welch"}}{Refit a Welch--Satterthwaite one-way test on
#'       only the data at the conditioning level (using
#'       \code{\link[stats]{oneway.test}} with \code{var.equal = FALSE}).
#'       Robust to heteroscedasticity within the conditioning level at the
#'       cost of fewer (and fractional) denominator \emph{df}.}
#'   }
#' @param adjust Multiple-comparison adjustment applied to the
#'   \emph{p}-values of the entire family of simple effects returned
#'   (\eqn{a + b} for \code{which = "both"}). One of \code{"none"}
#'   (default), \code{"bonferroni"}, or any sequential method supported by
#'   \code{\link[stats]{p.adjust}}: \code{"holm"}, \code{"hochberg"},
#'   \code{"BH"}, \code{"BY"}.
#' @param conf_level Confidence level for each row's partial \eqn{\eta^2}
#'   interval. Default \code{0.95}.
#'
#' @return A \code{data.frame} with one row per simple effect test and
#'   columns
#'   \code{effect}, \code{focal_factor}, \code{conditioning_factor},
#'   \code{conditioning_level}, \code{F_value}, \code{df_effect},
#'   \code{df_error}, \code{p_value}, \code{p_adjusted},
#'   \code{partial_eta_squared}, \code{lower_limit}, \code{upper_limit},
#'   \code{n_at_level}. Attributes \code{error_term}, \code{adjust},
#'   \code{conf_level}, \code{factor_A}, and \code{factor_B} record the
#'   call options.
#'
#' @details
#' \strong{What a simple effect is.} The simple main effect of \eqn{A} at
#' level \eqn{B = b_j} tests whether the \eqn{a} cell means at that single
#' level of \eqn{B} differ. It is the one-way analysis of variance of
#' \eqn{Y} on \eqn{A} restricted to observations with \eqn{B = b_j}. The
#' counterpart, the simple effect of \eqn{B} at \eqn{A = a_i}, is defined
#' symmetrically.
#'
#' \strong{Test statistic.} For the simple effect of \eqn{A} at \eqn{B = b_j},
#' let \eqn{\mathit{SS}_{A\,|\,b_j} = \sum_i n_{ij}\,(\bar{Y}_{ij\cdot} -
#' \bar{Y}_{\cdot j\cdot})^2} be the between-A sum of squares computed at
#' that level, and let \eqn{\mathit{MS}_{A\,|\,b_j} = \mathit{SS}_{A\,|\,b_j} / (a - 1)}.
#' \itemize{
#'   \item \emph{Pooled} \eqn{\mathit{MS}_W}: \eqn{F = \mathit{MS}_{A\,|\,b_j} /
#'     \mathit{MS}_W} with degrees of freedom \eqn{(a - 1,\, N - ab)},
#'     where \eqn{\mathit{MS}_W} and its \emph{df} come from the fitted full
#'     factorial model.
#'   \item \emph{Welch}: \eqn{F} and its (fractional) denominator
#'     \emph{df} come from \code{\link[stats]{oneway.test}}\code{(y ~ A,
#'     subset = (B == b_j), var.equal = FALSE)}.
#' }
#' Test statistics for \eqn{B} at \eqn{a_i} are computed by interchanging the
#' two factors.
#'
#' \strong{Choosing an error term.} The pooled denominator borrows strength
#' from all \eqn{N} observations and is the textbook default in
#' Maxwell, Delaney, and Kelley's treatment. Its validity rests on
#' homogeneity of variance across \emph{all} \eqn{a \times b} cells, not
#' merely within the conditioning level. When that assumption is doubtful
#', for example, if Levene's or the Brown--Forsythe test flags
#' heteroscedasticity, or if cell variances visibly differ, the Welch
#' option provides a level-conditional test that does not require
#' homogeneity across cells. The trade-off is denominator \emph{df}: pooled
#' carries the full \eqn{N - ab} residual \emph{df}, whereas Welch carries
#' the Welch--Satterthwaite \emph{df} based on the \eqn{a} cell variances
#' at that level only.
#'
#' \strong{Partial \eqn{\eta^2} and its CI.} The point estimate is
#' \deqn{\hat{\eta}^2_p = \frac{df_{\text{effect}}\, F}{df_{\text{effect}}\, F + df_{\text{error}}},}
#' computed from the \emph{F} and the \emph{df} actually used in the test
#' (so it reflects whichever error term was chosen). The confidence
#' interval is built by Steiger's (2004) transformation principle: a CI
#' for the noncentrality parameter \eqn{\lambda} of the \emph{F}
#' distribution is obtained via \code{\link{conf_limits_ncf}} and then
#' mapped through \eqn{\eta^2_p = \lambda / (\lambda + N_{\text{ref}})},
#' with \eqn{N_{\text{ref}}} taken to be the total study \emph{N} for the
#' pooled error term (treating the simple effect as a contrast within the
#' full factorial design) and the level-conditional sample size
#' \eqn{n_{|b_j}} for the Welch error term (since the Welch test uses only
#' those observations). When the lower limit on \eqn{\lambda} is not
#' identified (i.e., the observed \emph{F} is below its one-sided critical
#' value), the lower limit on \eqn{\eta^2_p} is set to 0; when the upper
#' limit is at infinity, the upper limit on \eqn{\eta^2_p} is set to 1.
#'
#' \strong{Multiplicity across the family.} With \code{which = "both"} the
#' family is the \eqn{a + b} simple effects returned in a single call; the
#' adjustment is applied to that entire family. If only one direction is
#' wanted, call the function twice with \code{which = "A_at_B"} and
#' \code{which = "B_at_A"} so each family is adjusted on its own. The
#' \code{"bonferroni"} adjustment is \eqn{p_{\text{adj}} = \min(1, m\, p)}
#' for \eqn{m} rows; the sequential methods (\code{"holm"},
#' \code{"hochberg"}, \code{"BH"}, \code{"BY"}) are computed via
#' \code{\link[stats]{p.adjust}}.
#'
#' \strong{When to perform simple effects.} It is no longer required that
#' a significant omnibus interaction precede simple effect testing
#' (Maxwell, Delaney, & Kelley, 2027). Simple effects are
#' informative whenever the substantive question is conditional on a level
#' of the other factor, and the multiplicity adjustment controls the
#' family-wise error rate independently of any interaction screen.
#'
#' \strong{Scope.} Only fixed-effects between-subjects designs with
#' \emph{exactly two} crossed factors are supported. Within-subjects or
#' mixed designs (\code{aovlist} fits) and three- or higher-way designs
#' are out of scope for this function. Per-cell sample sizes may be
#' unequal.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#' \emph{Designing experiments and analyzing data: A model comparison
#' perspective} (4th ed.). Routledge.
#'
#' Kelley, K. (2007). Confidence intervals for standardized
#' effect sizes: Theory, application, and implementation. \emph{Journal of
#' Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence
#' intervals and tests of close fit in the analysis of variance and
#' contrast analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' Welch, B. L. (1951). On the comparison of several mean values: An
#' alternative approach. \emph{Biometrika, 38}, 330--336.
#'
#' @examples
#' # 2 x 3 factorial: expectancy treatment (A) x grade (B) on the
#' # pygmalion data. Grades 4 through 6 are omitted so the family of
#' # simple effects stays short enough to read at a glance.
#' pyg <- pygmalion[pygmalion$grade <= 3, ]
#' pyg$grade <- factor(pyg$grade)
#' fit <- aov(iq_gain ~ treatment * grade, data = pyg)
#'
#' # Default: pooled MS_W, both families, no adjustment. The expectancy
#' # effect on IQ gain is concentrated in grades 1 and 2; at grade 3 the
#' # F is 0.004, so the lower limit on partial eta squared is clamped to
#' # 0 and the function notes the clamp in a warning.
#' simple_effects_AB(fit)
#'
#' # Only the simple effects of grade within each treatment level, with
#' # a Holm adjustment across that family of two tests.
#' simple_effects_AB(fit, which = "B_at_A", adjust = "holm")
#'
#' # Welch error term: refits a Welch one-way at each conditioning
#' # level. The Welch denominator df fall well below the pooled 157, so
#' # more of the lower limits are clamped to 0.
#' simple_effects_AB(fit, error_term = "welch")
#'
#' # Bonferroni across the full a + b = 5-test family.
#' simple_effects_AB(fit, adjust = "bonferroni")
#'
#' # Simulated 2 x 2 design with a known interaction pattern.
#' set.seed(113)
#' d <- simulate_ancova_factorial_data(
#'   a = 2, b = 2,
#'   mu_y    = c(50, 60, 55, 50),   # crossover at B = 2
#'   mu_x    = matrix(10, nrow = 4, ncol = 1),
#'   sigma_y = 8, sigma_x = 3, rho_y_x = 0,
#'   n       = 30
#' )
#' fit_sim <- aov(y ~ A * B, data = d)
#' simple_effects_AB(fit_sim, conf_level = 0.95)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{contrast_test}} for within-level pairwise or custom
#'   contrasts, \code{\link{eta_squared_partial}} and
#'   \code{\link{ci_eta_squared_partial}} for the omnibus effect size
#'   counterparts, \code{\link{conf_limits_ncf}} for the noncentrality
#'   machinery, \code{\link{ss_power_factorial_anova}} for power
#'   calculations on the omnibus factorial effects.
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#'
#' @export
#' @import stats

simple_effects_AB <- function(
  object,
  which       = "both",
  error_term  = "pooled",
  adjust      = "none",
  conf_level  = 0.95
) {
  which      <- match.arg(which,      c("both", "A_at_B", "B_at_A"))
  error_term <- match.arg(error_term, c("pooled", "welch"))
  adjust     <- match.arg(
    adjust,
    c("none", "bonferroni", "holm", "hochberg", "BH", "BY")
  )

  # Each row calls conf_limits_ncf(), which emits a warning whenever F falls
  # below the alpha_lower critical value of the central F (lower NCP limit
  # clamped to 0). For a family of a + b simple effects this can produce
  # several identical warnings. Muffle them per call and report a single
  # summary warning at the end. Pattern matches ss_aipe_R2().
  .clamp_count <- 0L
  on.exit({
    if (.clamp_count > 0L) {
      warning(sprintf(
        "The conf_limits_ncf() lower-limit clamp fired in %d of the simple effect rows (observed F below the alpha_lower critical value of the central F-distribution); the corresponding lower_limit on partial_eta_squared is clamped to 0. See ?conf_limits_ncf for the meaning of the clamp.",
        .clamp_count
      ), call. = FALSE)
    }
  }, add = TRUE)

  withCallingHandlers({

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!inherits(object, c("aov", "lm"))) {
    stop("'object' must be a fitted aov or lm.", call. = FALSE)
  }
  if (inherits(object, "aovlist")) {
    stop("'simple_effects_AB' supports between-subjects designs only; ",
         "aovlist (multi-stratum) fits are out of scope.", call. = FALSE)
  }

  # Extract the two factors from the model formula.
  tt    <- stats::terms(object)
  preds <- attr(tt, "term.labels")
  main_terms <- preds[!grepl(":", preds, fixed = TRUE)]
  if (length(main_terms) != 2L) {
    stop("'simple_effects_AB' requires a two-factor design; the model has ",
         length(main_terms), " main-effect terms (",
         paste(main_terms, collapse = ", "), ").", call. = FALSE)
  }
  has_interaction <- any(grepl(":", preds, fixed = TRUE))
  if (!has_interaction) {
    warning("Model has no interaction term; the pooled MS_W will include ",
            "interaction variance and the simple effect F tests will be ",
            "conservatively biased. Refit with ", main_terms[1], " * ",
            main_terms[2], " for the textbook pooled error term.",
            call. = FALSE)
  }

  A_name <- main_terms[1]
  B_name <- main_terms[2]

  mf    <- stats::model.frame(object)
  y_vec <- mf[[1]]
  A_vec <- mf[[A_name]]
  B_vec <- mf[[B_name]]
  if (!is.factor(A_vec)) A_vec <- factor(A_vec)
  if (!is.factor(B_vec)) B_vec <- factor(B_vec)
  A_vec  <- droplevels(A_vec)
  B_vec  <- droplevels(B_vec)
  A_levs <- levels(A_vec)
  B_levs <- levels(B_vec)
  if (length(A_levs) < 2L || length(B_levs) < 2L) {
    stop("Both factors must have at least two levels.", call. = FALSE)
  }

  # Pooled error term: MS_W and its df come from the full model's anova.
  tbl <- stats::anova(object)
  resid_row <- base::which(rownames(tbl) == "Residuals")
  if (length(resid_row) == 0L) {
    stop("Could not find a 'Residuals' row in anova(object).", call. = FALSE)
  }
  df_W_pooled <- as.numeric(tbl[resid_row, "Df"])
  MS_W_pooled <- as.numeric(tbl[resid_row, "Mean Sq"])
  N_total     <- stats::nobs(object)

  alpha_lower <- (1 - conf_level) / 2
  alpha_upper <- (1 - conf_level) / 2

  # Build a single simple effect row.
  compute_row <- function(focal_name, cond_name, cond_level) {
    cond_vec  <- mf[[cond_name]]
    focal_vec <- mf[[focal_name]]
    if (!is.factor(cond_vec))  cond_vec  <- factor(cond_vec)
    if (!is.factor(focal_vec)) focal_vec <- factor(focal_vec)

    mask  <- as.character(cond_vec) == cond_level
    y_sub <- y_vec[mask]
    x_sub <- droplevels(focal_vec[mask])
    n_at  <- length(y_sub)
    k     <- length(levels(x_sub))
    if (k < 2L || n_at < k + 1L) {
      return(data.frame(
        effect              = paste0(focal_name, " | ", cond_name, " = ", cond_level),
        focal_factor        = focal_name,
        conditioning_factor = cond_name,
        conditioning_level  = cond_level,
        F_value             = NA_real_,
        df_effect           = NA_real_,
        df_error            = NA_real_,
        p_value             = NA_real_,
        partial_eta_squared = NA_real_,
        lower_limit         = NA_real_,
        upper_limit         = NA_real_,
        n_at_level          = n_at,
        stringsAsFactors    = FALSE
      ))
    }
    df1 <- k - 1L

    if (error_term == "pooled") {
      grand      <- mean(y_sub)
      cell_means <- tapply(y_sub, x_sub, mean)
      cell_n     <- as.numeric(table(x_sub))
      ss_simple  <- sum(cell_n * (cell_means - grand)^2)
      ms_simple  <- ss_simple / df1
      F_val      <- ms_simple / MS_W_pooled
      df2        <- df_W_pooled
      N_for_ci   <- N_total
    } else {
      ow    <- tryCatch(stats::oneway.test(y_sub ~ x_sub, var.equal = FALSE),
                        error = function(e) NULL)
      if (is.null(ow) || is.na(as.numeric(ow$statistic))) {
        F_val <- NA_real_; df2 <- NA_real_; N_for_ci <- n_at
      } else {
        F_val    <- as.numeric(ow$statistic)
        df1      <- as.numeric(ow$parameter[["num df"]])
        df2      <- as.numeric(ow$parameter[["denom df"]])
        N_for_ci <- n_at
      }
    }

    p_unadj <- if (is.na(F_val)) NA_real_ else
      stats::pf(F_val, df1, df2, lower.tail = FALSE)

    eta_p <- if (is.na(F_val)) NA_real_ else
      (df1 * F_val) / (df1 * F_val + df2)

    if (is.na(F_val)) {
      lower_limit <- NA_real_; upper_limit <- NA_real_
    } else {
      ncp_lims <- tryCatch(
        conf_limits_ncf(
          F_value     = F_val,
          df_1        = df1,
          df_2        = df2,
          alpha_lower = alpha_lower,
          alpha_upper = alpha_upper,
          conf_level  = NULL,
          verbose     = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(ncp_lims)) {
        lower_limit <- NA_real_; upper_limit <- NA_real_
      } else {
        lo_ncp <- ncp_lims$value[ncp_lims$term == "lower_limit"]
        up_ncp <- ncp_lims$value[ncp_lims$term == "upper_limit"]
        lower_limit <- if (length(lo_ncp) == 0L || is.na(lo_ncp) || lo_ncp <= 0)
          0 else lo_ncp / (lo_ncp + N_for_ci)
        upper_limit <- if (length(up_ncp) == 0L || is.na(up_ncp))
          NA_real_ else if (is.infinite(up_ncp)) 1 else
          up_ncp / (up_ncp + N_for_ci)
      }
    }

    data.frame(
      effect              = paste0(focal_name, " | ", cond_name, " = ", cond_level),
      focal_factor        = focal_name,
      conditioning_factor = cond_name,
      conditioning_level  = cond_level,
      F_value             = F_val,
      df_effect           = df1,
      df_error            = df2,
      p_value             = p_unadj,
      partial_eta_squared = eta_p,
      lower_limit         = lower_limit,
      upper_limit         = upper_limit,
      n_at_level          = n_at,
      stringsAsFactors    = FALSE
    )
  }

  rows <- list()
  if (which %in% c("both", "A_at_B")) {
    for (bl in B_levs) {
      rows[[length(rows) + 1L]] <- compute_row(A_name, B_name, bl)
    }
  }
  if (which %in% c("both", "B_at_A")) {
    for (al in A_levs) {
      rows[[length(rows) + 1L]] <- compute_row(B_name, A_name, al)
    }
  }
  out <- do.call(rbind, rows)

  # Multiplicity adjustment across the returned family.
  if (adjust == "none") {
    out$p_adjusted <- out$p_value
  } else if (adjust == "bonferroni") {
    m               <- sum(!is.na(out$p_value))
    out$p_adjusted  <- pmin(1, out$p_value * m)
  } else {
    out$p_adjusted  <- stats::p.adjust(out$p_value, method = adjust)
  }

  out <- out[, c("effect", "focal_factor", "conditioning_factor",
                 "conditioning_level", "F_value", "df_effect", "df_error",
                 "p_value", "p_adjusted", "partial_eta_squared",
                 "lower_limit", "upper_limit", "n_at_level")]
  rownames(out) <- NULL

  attr(out, "error_term") <- error_term
  attr(out, "adjust")     <- adjust
  attr(out, "conf_level") <- conf_level
  attr(out, "factor_A")   <- A_name
  attr(out, "factor_B")   <- B_name
  .as_dmar_tbl(out, conf_level = conf_level)
  }, warning = function(w) {
    if (inherits(w, "dmar_ncf_clamp")) {
      .clamp_count <<- .clamp_count + 1L
      invokeRestart("muffleWarning")
    }
  })
}
