# Between-subjects factorial ANOVA for unbalanced designs, with the sum-of-
# squares type chosen explicitly and computed by model comparison (no car).
#' Between-Subjects Factorial ANOVA for Unbalanced Designs
#'
#' Fits a between-subjects factorial ANOVA (two or more crossed factors) and
#' returns each effect's sum of squares, \emph{F} test, and partial
#' \eqn{\eta^2} and partial \eqn{\omega^2} with confidence intervals, using a
#' sum-of-squares type the user chooses. It is built for the case that makes
#' the choice matter, an \emph{unbalanced} design (unequal cell sizes), where
#' the Type I, II, and III sums of squares differ. For a balanced design the
#' three types coincide and the choice is immaterial.
#'
#' Every sum of squares is computed as a model comparison, the increase in
#' error sum of squares when an effect's parameters are removed from a model,
#' following the model comparison development of Maxwell, Delaney, and Kelley
#' (2027, Chapter 7). The computation uses only base R (\code{stats}); it does
#' not depend on the \pkg{car} package, though it agrees with
#' \code{car::Anova()} to numerical precision.
#'
#' @param formula A two-sided \code{\link{formula}} naming the numeric
#'   response and the crossed factors, for example \code{y ~ A * B} or
#'   \code{y ~ A * B * C}. Predictors are coerced to factors. Write the full
#'   crossing you want tested; \code{A * B} expands to \code{A + B + A:B}.
#' @param data A \code{data.frame} containing the response and the factors.
#' @param ss_type The sum-of-squares type: \code{1}, \code{2}, or \code{3}
#'   (equivalently \code{"I"}, \code{"II"}, or \code{"III"}). Default \code{3}
#'   (Type III), the common reporting default. See Details for what each type
#'   tests.
#' @param conf_level Confidence level for the effect size confidence intervals.
#'   Default \code{0.95}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per effect
#'   plus a \code{Residuals} row. Columns are the effect label (\code{effect}),
#'   the sum of squares (\code{SS}), degrees of freedom (\code{df}), the
#'   \emph{F} statistic (\code{F_value}) and its \emph{p}-value
#'   (\code{p_value}), and partial \eqn{\eta^2} and partial \eqn{\omega^2} with
#'   their lower and upper confidence limits. The chosen sum-of-squares type is
#'   recorded on the object (\code{attr(x, "ss_type")}) and printed beneath the
#'   table. The residual row carries only \code{SS} and \code{df}. Stored
#'   values keep full precision; the display rounds (see \code{\link{dmar_tbl}}).
#'
#' @details
#' \strong{The types as model comparisons.} A sum of squares for an effect is
#' the increase in the error sum of squares when the effect's parameters are
#' dropped from the model, \code{SS = E(restricted) - E(full)}. The three
#' conventional types differ only in which other effects the full and
#' restricted models hold in common (Maxwell, Delaney, and Kelley, 2027,
#' Chapter 7; Overall and Spiegel, 1969):
#' \itemize{
#'   \item \strong{Type I (sequential).} Each effect is adjusted only for the
#'     effects listed before it in \code{formula}: \code{SS(A)}, then
#'     \code{SS(B | A)}, then \code{SS(A:B | A, B)}. The parts sum to the model
#'     sum of squares, but the answer depends on the order of the terms.
#'   \item \strong{Type II.} Each effect is adjusted for every other effect that
#'     does not contain it (it respects marginality): a main effect is adjusted
#'     for the other main effects but not for the interactions that contain it.
#'     Type II is the most powerful choice when the interaction is null and does
#'     not depend on how the factors are coded (Overall and Spiegel's Method 2;
#'     Appelbaum and Cramer, 1974).
#'   \item \strong{Type III.} Each effect is adjusted for all other effects,
#'     including the higher order interactions that contain it. This tests each
#'     main effect as a contrast on the unweighted marginal means, so it is
#'     computed here with sum-to-zero contrasts (\code{\link[stats]{contr.sum}}),
#'     which is what makes the Type III main-effect test the intended one. It is
#'     the default reported by many programs.
#' }
#' For a balanced design the effects are orthogonal and the three types are
#' identical; the distinction is a property of unbalanced (nonorthogonal) data.
#'
#' \strong{Reading the main effects when an interaction is present.} When an
#' interaction is real, the marginal main-effect tests, of any type, are usually
#' not the question of interest; examine the interaction and the simple effects
#' instead (Maxwell, Delaney, and Kelley, 2027). See the vignette
#' \emph{Sums of Squares in Nonorthogonal Designs} for a worked comparison.
#'
#' \strong{Effect sizes.} Partial \eqn{\eta^2} and partial \eqn{\omega^2} are
#' formed for each effect from its \emph{F}, its degrees of freedom, and the
#' error degrees of freedom, with noncentral \emph{F} confidence intervals
#' (\code{\link{ci_eta_squared_partial}}, \code{\link{ci_omega_squared}}). The
#' confidence limits are those of the interval for the population proportion of
#' variance the effect accounts for (Kelley, 2007). Partial \eqn{\eta^2} and
#' partial \eqn{\omega^2} are two point estimators of that same population
#' quantity, partial \eqn{\omega^2} correcting the upward bias of partial
#' \eqn{\eta^2}, so the two estimators differ but share the interval.
#'
#' \strong{Estimability.} All cells must be filled. If a factor combination is
#' empty the design is rank deficient and the factorial effects are not all
#' estimable; the function stops with a message rather than return a value that
#' depends on an arbitrary choice.
#'
#' @references
#' Appelbaum, M. I., & Cramer, E. M. (1974). Some problems in the nonorthogonal
#'   analysis of variance. \emph{Psychological Bulletin, 81}(6), 335--343.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective} (4th ed.).
#'   Routledge. (See Chapter 7 on higher order between-subjects designs.)
#'
#' Overall, J. E., & Spiegel, D. K. (1969). Concerning least squares analysis
#'   of experimental data. \emph{Psychological Bulletin, 72}(5), 311--322.
#'
#' @seealso \code{\link{ancova}} for a covariate-adjusted one-way design,
#'   \code{\link{mixed_anova}} for mixed-model (fixed and random) \emph{F}
#'   ratios, \code{\link{manova_split_plot}} for the multivariate mixed design, and
#'   \code{\link{ci_eta_squared_partial}} / \code{\link{ci_omega_squared}} for
#'   the effect size intervals.
#'
#' @examples
#' # An unbalanced two-factor design: IQ gain in the pygmalion expectancy
#' # experiment, by treatment and by lower (grades 1 and 2) versus upper
#' # (grades 3 through 6) grades, where the expectancy effect concentrated
#' # in the lower grades. The cell sizes are unequal, so the Types differ.
#' pyg <- pygmalion
#' pyg$grade_band <- factor(ifelse(pyg$grade <= 2, "lower", "upper"))
#' factorial_anova(iq_gain ~ treatment * grade_band, data = pyg)
#'
#' # The same design under Type II (adjusts each main effect for the other
#' # main effect, but not for the interaction). Here the main effect F
#' # statistics drop, and a warning reports that the affected noncentral F
#' # lower limits are clamped to 0.
#' factorial_anova(iq_gain ~ treatment * grade_band, data = pyg, ss_type = 2)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#'
#' @export
factorial_anova <- function(formula, data, ss_type = 3L, conf_level = 0.95) {
  if (!inherits(formula, "formula"))
    stop("'formula' must be a formula, for example y ~ A * B.", call. = FALSE)
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  ss <- switch(toupper(as.character(ss_type)),
               "1" = 1L, "I" = 1L, "2" = 2L, "II" = 2L, "3" = 3L, "III" = 3L,
               stop("'ss_type' must be 1, 2, or 3 (or \"I\", \"II\", \"III\").",
                    call. = FALSE))

  Terms <- stats::terms(formula)
  if (attr(Terms, "response") != 1L)
    stop("'formula' must have a response, for example y ~ A * B.", call. = FALSE)
  term_labels <- attr(Terms, "term.labels")
  if (length(term_labels) < 1L)
    stop("'formula' must include at least one predictor.", call. = FALSE)

  mf <- stats::model.frame(Terms, data = data)
  y  <- stats::model.response(mf)
  if (!is.numeric(y))
    stop("The response must be numeric.", call. = FALSE)

  # A factorial ANOVA treats every predictor as categorical.
  fvars <- names(mf)[-1L]
  for (v in fvars) if (!is.factor(mf[[v]])) mf[[v]] <- factor(mf[[v]])
  for (v in fvars) if (nlevels(mf[[v]]) < 2L)
    stop(sprintf("Factor '%s' has fewer than two levels.", v), call. = FALSE)

  # Sum-to-zero contrasts make the Type III main-effect tests the intended
  # unweighted-marginal-means tests. Type I and Type II sums of squares are
  # invariant to the within-factor basis, so building all three from this one
  # model matrix is correct for every type.
  contr_list <- stats::setNames(as.list(rep("contr.sum", length(fvars))), fvars)
  X <- stats::model.matrix(Terms, mf, contrasts.arg = contr_list)
  asg <- attr(X, "assign")
  n <- length(y)

  qr_full <- qr(X)
  if (qr_full$rank < ncol(X))
    stop("The model is rank deficient, which usually means some cells are ",
         "empty. The factorial effects are not all estimable; use a design ",
         "with every cell filled or drop the interaction.", call. = FALSE)
  sse_full <- sum(qr.resid(qr_full, y)^2)
  df_resid <- n - qr_full$rank
  mse <- sse_full / df_resid

  # Error sum of squares of the model that keeps the intercept plus the terms
  # indexed by `keep`; the workhorse for every model comparison below.
  sse_of <- function(keep) {
    cols <- which(asg %in% c(0L, keep))
    sum(qr.resid(qr(X[, cols, drop = FALSE]), y)^2)
  }

  # Variable set of each term, for the Type II marginality rule.
  fmat <- attr(Terms, "factors")
  vars_of <- lapply(seq_along(term_labels), function(j) rownames(fmat)[fmat[, j] > 0])

  n_eff <- length(term_labels)
  SS <- df <- Fv <- pv <- numeric(n_eff)
  eta <- eta_lo <- eta_hi <- omega <- omega_lo <- omega_hi <- numeric(n_eff)

  # The per-effect ci_eta_squared_partial() / ci_omega_squared() calls each go
  # through conf_limits_ncf(), whose noncentral F lower limit is often clamped
  # to 0 for a small effect; surfacing that warning once per effect produces
  # many identical messages. Muffle the per-effect warnings, count them, and
  # emit a single summary warning at the end (via on.exit so it fires on any
  # return path). This is the ss_aipe_R2() dedup pattern.
  .clamp_count <- 0L
  on.exit({
    if (.clamp_count > 0L) {
      warning(sprintf(
        "The noncentral F lower-limit clamp in conf_limits_ncf() fired for %d of the effect size confidence intervals; the affected lower limits were clamped to 0. See ?conf_limits_ncf for the meaning of the clamp.",
        .clamp_count
      ), call. = FALSE)
    }
  }, add = TRUE)

  withCallingHandlers(
  for (k in seq_len(n_eff)) {
    df_k <- sum(asg == k)
    ss_k <- switch(ss,
      # Type I: adjust only for the terms entered before k.
      { before <- if (k > 1L) seq_len(k - 1L) else integer(0)
        sse_of(before) - sse_of(c(before, k)) },
      # Type II: adjust for every term that does not contain k.
      { contains_k <- vapply(seq_len(n_eff),
                             function(j) all(vars_of[[k]] %in% vars_of[[j]]),
                             logical(1L))
        m0 <- setdiff(which(!contains_k), k)
        sse_of(m0) - sse_of(c(m0, k)) },
      # Type III: adjust for all other terms, including those containing k.
      { sse_of(setdiff(seq_len(n_eff), k)) - sse_full }
    )
    ss_k <- max(ss_k, 0)
    f_k <- (ss_k / df_k) / mse
    p_k <- stats::pf(f_k, df_k, df_resid, lower.tail = FALSE)
    e <- ci_eta_squared_partial(F_value = f_k, df_effect = df_k,
                                df_error = df_resid, N = n, conf_level = conf_level)
    o <- ci_omega_squared(F_value = f_k, df_effect = df_k,
                          df_error = df_resid, N = n, conf_level = conf_level)
    SS[k] <- ss_k; df[k] <- df_k; Fv[k] <- f_k; pv[k] <- p_k
    eta[k]   <- e$eta_squared_partial; eta_lo[k]   <- e$lower_limit; eta_hi[k]   <- e$upper_limit
    omega[k] <- o$omega_squared;       omega_lo[k] <- o$lower_limit; omega_hi[k] <- o$upper_limit
  },
  warning = function(w) {
    if (inherits(w, "dmar_ncf_clamp")) {
      .clamp_count <<- .clamp_count + 1L
      invokeRestart("muffleWarning")
    }
  })

  out <- data.frame(
    effect = c(term_labels, "Residuals"),
    SS      = c(SS, sse_full),
    df      = c(df, df_resid),
    F_value = c(Fv, NA_real_),
    p_value = c(pv, NA_real_),
    eta_squared_partial       = c(eta, NA_real_),
    eta_squared_partial_lower = c(eta_lo, NA_real_),
    eta_squared_partial_upper = c(eta_hi, NA_real_),
    omega_squared_partial       = c(omega, NA_real_),
    omega_squared_partial_lower = c(omega_lo, NA_real_),
    omega_squared_partial_upper = c(omega_hi, NA_real_),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, ss_type = ss)
}
