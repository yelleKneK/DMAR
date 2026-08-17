# Two-factor within-subjects ANOVA with sphericity adjustments per effect.
#' Two-Factor Within-Subjects ANOVA With Sphericity Adjustments
#'
#' Computes the full two-factor within-subjects ANOVA table, main
#' effects \emph{A}, \emph{B}, and the \eqn{A \times B} interaction
#' with each effect tested against its own residual stratum and with
#' Greenhouse-Geisser, Huynh-Feldt, and lower-bound sphericity
#' adjustments applied per effect. Returns a tidy long-form
#' \code{data.frame} that composes with the rest of DMAR.
#'
#' @param data A \code{data.frame} in long format, with one row per
#'   subject-by-cell observation.
#' @param outcome Character name of the response column in \code{data}.
#' @param factor_A Character name of the first within-subjects factor.
#' @param factor_B Character name of the second within-subjects factor.
#' @param subject Character name of the subject-id column.
#'
#' @return A \code{data.frame} with rows for each of the three
#'   effects (\code{A}, \code{B}, \code{A:B}) crossed with each
#'   sphericity adjustment (\code{none}, \code{Greenhouse-Geisser},
#'   \code{Huynh-Feldt}, \code{lower_bound}). Columns: \code{effect},
#'   \code{adjustment}, \code{F_value}, \code{df_1}, \code{df_2},
#'   \code{p_value}, \code{epsilon}, \code{partial_eta_squared}. When
#'   an effect has too few subjects for its epsilon to be estimable
#'   (\eqn{n - 1} smaller than the effect's numerator degrees of
#'   freedom; see Details), the Greenhouse-Geisser and Huynh-Feldt rows
#'   for that effect carry \code{NA} in \code{epsilon}, \code{df_1},
#'   \code{df_2}, and \code{p_value}, and a single warning names the
#'   condition; the unadjusted and lower-bound rows are unaffected.
#'
#' @details
#' \strong{Design.} The two within-subjects factors \emph{A} (with
#' \eqn{a} levels) and \emph{B} (with \eqn{b} levels) are fully
#' crossed; every subject contributes \eqn{a \cdot b} observations.
#' Each of the three fixed effects is tested against its own
#' subject-by-effect residual stratum:
#' \itemize{
#'   \item \eqn{F_A = MS_A / MS_{A:S}}
#'   \item \eqn{F_B = MS_B / MS_{B:S}}
#'   \item \eqn{F_{A:B} = MS_{A:B} / MS_{A:B:S}}
#' }
#'
#' \strong{Sphericity.} Each effect's univariate \emph{F}-ratio assumes
#' sphericity of its corresponding subject-by-effect residual
#' covariance matrix. Three adjustments are reported per effect:
#' Greenhouse-Geisser (Greenhouse & Geisser, 1959), Huynh-Feldt
#' (Huynh & Feldt, 1976), and the lower bound \eqn{\epsilon = 1 / df},
#' where \eqn{df} is the effect's numerator degrees of freedom; this is
#' the smallest value \eqn{\epsilon} can attain, reached under maximal
#' departure from sphericity.
#'
#' \strong{Subjects needed to estimate epsilon.} The Greenhouse-Geisser
#' epsilon for an effect with \eqn{q} numerator degrees of freedom is
#' estimated from the sample covariance matrix of \eqn{q} orthonormal
#' contrasts among the effect's cell means, a different matrix for each
#' effect (Maxwell, Delaney, & Kelley, 2027, Chapters 11 and 12). That
#' matrix has rank at most \eqn{n - 1}, so when \eqn{n - 1 < q} it is
#' necessarily singular; the same rank deficiency makes the
#' multivariate approach to a within-subjects design mathematically
#' impossible when \eqn{n < a} (Maxwell, Delaney, & Kelley, 2027,
#' Chapter 13). The Greenhouse-Geisser formula still returns a number
#' in that case, but the number is an artifact of the rank deficiency
#' rather than an estimate: it cannot exceed \eqn{(n - 1)/q} no matter
#' what the population epsilon is, even under exact sphericity, where
#' the population value is 1. Rather than report a value the design
#' cannot support, the function reports \code{NA} for the
#' Greenhouse-Geisser and Huynh-Feldt rows of any effect with
#' \eqn{n - 1 < q} and issues a single warning naming the condition
#' (\code{car::Anova} likewise declines to report the corrections for
#' an effect whose error matrix is singular). The unadjusted row and
#' the lower-bound row remain: the lower bound \eqn{1/q} is Geisser and
#' Greenhouse's a priori bound on epsilon, valid no matter how badly
#' sphericity is violated, and it requires no estimate of the
#' covariance matrix (Maxwell, Delaney, & Kelley, 2027, Chapter 11).
#'
#' \strong{Per-effect partial \eqn{\eta^2}.} Computed as
#' \eqn{SS_\mathrm{effect} / (SS_\mathrm{effect} + SS_\mathrm{effect,\, error})}
#' using the appropriate subject-by-effect residual sum of squares.
#'
#' \strong{Balanced data assumed.} The implementation assumes a fully
#' balanced design (every subject observed once in every cell). When
#' the design is unbalanced, the function errors and recommends a
#' mixed-effects fit via \code{\link[lme4]{lmer}}.
#'
#' \strong{Sums of squares are unambiguous here.} Because the design is
#' balanced, the within-subjects factors are orthogonal and the Type I,
#' Type II, and Type III sums of squares for each effect coincide. A
#' Type toggle is therefore not meaningful, and the reported
#' decomposition is unambiguous: the sum of squares attributed to each
#' effect does not depend on the order in which terms enter the model
#' (Maxwell, Delaney, & Kelley, 2027, Chapter 12).
#'
#' @references
#' Greenhouse, S. W., & Geisser, S. (1959). On methods in the analysis
#'   of profile data. \emph{Psychometrika, 24}(2), 95--112.
#'
#' Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction
#'   for degrees of freedom from sample data in randomized block and
#'   split-plot designs. \emph{Journal of Educational Statistics, 1}(1),
#'   69--82.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapters 11--13.)
#'
#' @seealso \code{\link{anova_within}}, \code{\link{mauchly_test}},
#'   \code{\link{epsilon_corrections}}
#'
#' @examples
#' # 1. Balanced 3x2 within-subjects design simulated for illustration.
#' set.seed(113)
#' n_sub <- 12
#' grid  <- expand.grid(subject = factor(1:n_sub),
#'                      A = factor(c("a1", "a2", "a3")),
#'                      B = factor(c("b1", "b2")))
#' grid$y <- with(grid,
#'   3 * (A == "a2") + 1 * (A == "a3") +
#'   2 * (B == "b2") + 1 * ((A == "a3") & (B == "b2")) +
#'   rnorm(nrow(grid), 0, 1) +
#'   rep(rnorm(n_sub, 0, 1.5), times = 6))
#' anova_within_two_way(grid, outcome = "y", factor_A = "A",
#'                      factor_B = "B", subject = "subject")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family within-subjects analysis
#'
#' @export

anova_within_two_way <- function(data, outcome, factor_A, factor_B, subject) {
  for (nm in c("outcome", "factor_A", "factor_B", "subject")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L)
      stop(sprintf("'%s' must be a single column name.", nm))
    if (!v %in% names(data))
      stop(sprintf("Column '%s' not found in 'data'.", v))
  }

  d <- data.frame(
    y   = data[[outcome]],
    A   = factor(data[[factor_A]]),
    B   = factor(data[[factor_B]]),
    sub = factor(data[[subject]]),
    stringsAsFactors = FALSE
  )
  d <- d[stats::complete.cases(d), , drop = FALSE]
  a_lev <- levels(d$A); b_lev <- levels(d$B); s_lev <- levels(d$sub)
  a <- length(a_lev); b <- length(b_lev); n <- length(s_lev)

  if (n < 2L) stop("Need at least 2 subjects.")
  if (a < 2L || b < 2L)
    stop("Both within-subjects factors must have at least 2 levels.")
  if (nrow(d) != n * a * b)
    stop("Design is not balanced: every subject must appear once in every (A, B) cell.",
         "  Use a mixed-effects fit (e.g., lme4::lmer) for unbalanced data.")

  fit <- stats::aov(y ~ A * B + Error(sub / (A * B)), data = d)
  smry <- summary(fit)

  # The aov-with-Error summary returns four strata named like
  # "Error: sub", "Error: sub:A", "Error: sub:B", "Error: sub:A:B".
  # Pull the matching stratum and the requested row out of its table.
  # aov stores row names with trailing whitespace, so trim before matching.
  get_row <- function(s, name) {
    tab <- s[[1L]]
    if (is.null(tab)) return(NULL)
    rn <- trimws(rownames(tab))
    idx <- which(rn == name)
    if (length(idx) == 0L) return(NULL)
    tab[idx, , drop = FALSE]
  }
  strata_names <- names(smry)
  pull <- function(suffix, effect) {
    pat <- paste0("Error: sub:", suffix, "$")
    s_idx <- grep(pat, strata_names)
    if (length(s_idx) != 1L)
      stop(sprintf("Could not locate stratum 'sub:%s' in aov() summary.", suffix))
    get_row(smry[[s_idx]], effect)
  }

  rowA   <- pull("A",   "A")
  errA   <- pull("A",   "Residuals")
  rowB   <- pull("B",   "B")
  errB   <- pull("B",   "Residuals")
  rowAB  <- pull("A:B", "A:B")
  errAB  <- pull("A:B", "Residuals")

  build <- function(rr, ee, effect_label, eps_vec) {
    ss_e   <- rr[["Sum Sq"]]
    ss_err <- ee[["Sum Sq"]]
    df_e   <- rr[["Df"]]
    df_err <- ee[["Df"]]
    ms_e   <- ss_e   / df_e
    ms_err <- ss_err / df_err
    F_v    <- ms_e / ms_err
    p_eta  <- ss_e / (ss_e + ss_err)
    rows <- list()
    rows[["none"]] <- data.frame(
      effect = effect_label, adjustment = "none",
      F_value = F_v, df_1 = df_e, df_2 = df_err,
      p_value = stats::pf(F_v, df_e, df_err, lower.tail = FALSE),
      epsilon = NA_real_, partial_eta_squared = p_eta,
      stringsAsFactors = FALSE)
    for (nm in names(eps_vec)) {
      eps <- eps_vec[[nm]]
      d1 <- df_e * eps; d2 <- df_err * eps
      rows[[nm]] <- data.frame(
        effect = effect_label, adjustment = nm,
        F_value = F_v, df_1 = d1, df_2 = d2,
        p_value = stats::pf(F_v, d1, d2, lower.tail = FALSE),
        epsilon = eps, partial_eta_squared = p_eta,
        stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  }

  # Per-effect epsilon: epsilon needs the subject-by-effect contrast covariance.
  Y_array <- .within_two_way_to_array(d, a_lev, b_lev, s_lev)
  eps_A   <- .epsilon_from_mat_A(Y_array)
  eps_B   <- .epsilon_from_mat_B(Y_array)
  eps_AB  <- .epsilon_from_mat_AB(Y_array)

  # An effect whose numerator degrees of freedom exceed n - 1 has a
  # singular contrast covariance matrix, so its Greenhouse-Geisser and
  # Huynh-Feldt epsilons come back NA from .eps_from_cov(); warn once,
  # naming the condition and every affected effect (see Details).
  eps_na <- c(A = is.na(eps_A[["Greenhouse-Geisser"]]),
              B = is.na(eps_B[["Greenhouse-Geisser"]]),
              `A:B` = is.na(eps_AB[["Greenhouse-Geisser"]]))
  if (any(eps_na)) {
    df_eff <- c(A = a - 1L, B = b - 1L, `A:B` = (a - 1L) * (b - 1L))
    lab <- paste0(names(eps_na)[eps_na], " (df = ", df_eff[eps_na], ")",
                  collapse = ", ")
    warning(sprintf(paste0(
      "With n = %d subjects, n - 1 is smaller than the numerator degrees ",
      "of freedom for %s, so the covariance matrix of the effect's ",
      "contrasts is singular and the Greenhouse-Geisser and Huynh-Feldt ",
      "epsilons cannot be estimated (the sample value cannot exceed ",
      "(n - 1)/df regardless of the population epsilon). Those rows ",
      "report NA; the unadjusted and lower-bound rows do not depend on ",
      "the estimated epsilon and are reported as usual."),
      n, lab))
  }

  out <- rbind(
    build(rowA,  errA,  "A",   eps_A),
    build(rowB,  errB,  "B",   eps_B),
    build(rowAB, errAB, "A:B", eps_AB)
  )
  rownames(out) <- NULL
  out
}

# ---- internals ----

# Build n x a x b array from long-form, in canonical order.
.within_two_way_to_array <- function(d, a_lev, b_lev, s_lev) {
  n <- length(s_lev); a <- length(a_lev); b <- length(b_lev)
  arr <- array(NA_real_, c(n, a, b),
               dimnames = list(s_lev, a_lev, b_lev))
  for (i in seq_len(nrow(d)))
    arr[as.character(d$sub[i]), as.character(d$A[i]), as.character(d$B[i])] <-
      d$y[i]
  if (anyNA(arr)) stop("Missing cells after restructuring; check the design.")
  arr
}

# Greenhouse-Geisser / Huynh-Feldt epsilon from a residual contrast covariance.
.eps_from_cov <- function(S, n) {
  # S is the (q x q) covariance of q = df_effect orthonormal contrasts.
  q <- nrow(S)
  if (q < 1L) return(c(`Greenhouse-Geisser` = 1, `Huynh-Feldt` = 1,
                       `lower_bound` = 1))
  if (q == 1L) return(c(`Greenhouse-Geisser` = 1, `Huynh-Feldt` = 1,
                        `lower_bound` = 1))
  if (n - 1L < q) {
    # S has rank at most n - 1, so it is singular here, and the
    # Greenhouse-Geisser formula returns a value bounded above by
    # (n - 1)/q no matter what the population epsilon is; the number is
    # an artifact of the rank deficiency, not an estimate. Report NA
    # for both estimates (the caller warns once, naming the affected
    # effects). The lower bound 1/q is not estimated from S and
    # remains valid.
    return(c(`Greenhouse-Geisser` = NA_real_, `Huynh-Feldt` = NA_real_,
             `lower_bound` = 1 / q))
  }
  trS  <- sum(diag(S))
  tr2S <- sum(diag(S %*% S))
  eps_gg <- trS^2 / (q * tr2S)
  eps_gg <- min(1, max(1 / q, eps_gg))
  eps_hf <- (n * q * eps_gg - 2) / (q * (n - 1 - q * eps_gg))
  eps_hf <- min(1, max(1 / q, eps_hf))
  eps_lb <- 1 / q
  c(`Greenhouse-Geisser` = eps_gg, `Huynh-Feldt` = eps_hf,
    `lower_bound` = eps_lb)
}

.epsilon_from_mat_A <- function(arr) {
  # Marginal A means per subject (collapse over B):
  M <- apply(arr, c(1, 2), mean)        # n x a
  H <- stats::contr.helmert(ncol(M)); M_orth <- t(qr.Q(qr(H)))
  S <- stats::cov(M %*% t(M_orth))
  .eps_from_cov(S, nrow(M))
}
.epsilon_from_mat_B <- function(arr) {
  M <- apply(arr, c(1, 3), mean)        # n x b
  H <- stats::contr.helmert(ncol(M)); M_orth <- t(qr.Q(qr(H)))
  S <- stats::cov(M %*% t(M_orth))
  .eps_from_cov(S, nrow(M))
}
.epsilon_from_mat_AB <- function(arr) {
  n <- dim(arr)[1]; a <- dim(arr)[2]; b <- dim(arr)[3]
  # Interaction contrast: (a-1)(b-1) interaction scores per subject.
  HA <- t(qr.Q(qr(stats::contr.helmert(a))))   # (a-1) x a
  HB <- t(qr.Q(qr(stats::contr.helmert(b))))   # (b-1) x b
  # matrix(arr, n, a * b) flattens the cells with the A index varying fastest
  # (column j + (k - 1) * a is cell A = j, B = k), so the interaction basis must
  # use kronecker(HB, HA) to match that ordering. kronecker(HA, HB) would order
  # its columns with B fastest, applying each contrast weight to the wrong cell
  # and corrupting the interaction epsilon whenever a and b differ.
  K  <- kronecker(HB, HA)                       # (a-1)(b-1) x ab, A-fastest columns
  Y  <- matrix(arr, n, a * b)
  S  <- stats::cov(Y %*% t(K))
  .eps_from_cov(S, n)
}
