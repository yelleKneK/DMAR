# LRT comparing covariance structures for the same fixed effects.
#' Likelihood-Ratio Comparison of Covariance Structures
#'
#' Fits a long-format within-subjects regression under a menu of
#' variance-covariance structures, from independence through the
#' unstructured form, and returns a comparison table of
#' log-likelihood, AIC, BIC, and pairwise likelihood-ratio tests
#' against the most general structure (UN). Wraps
#' \code{\link[nlme]{gls}}.
#'
#' @param data Long-format \code{data.frame} with one row per
#'   subject-by-condition observation.
#' @param outcome Character name of the response column.
#' @param subject Character name of the subject-id column.
#' @param time Character name of the time / within-subjects factor
#'   column.
#' @param fixed_effects Right-hand-side formula for the fixed effects
#'   (default: \code{~ time}).
#' @param structures Character vector of structures to fit. Any subset
#'   of \code{c("IND", "CS", "CSH", "AR1", "ARH1", "TOEP", "TOEPH",
#'   "UN")} (default: all eight). Matching is case insensitive, so
#'   lowercase aliases such as \code{"cs"}, \code{"ar1"}, \code{"csh"},
#'   \code{"arh1"}, \code{"toep"}, and \code{"un"} are accepted and
#'   normalized to their canonical uppercase labels.
#'
#' @return A \code{data.frame} with one row per structure.
#'   Columns: \code{structure}, \code{log_lik}, \code{AIC}, \code{BIC},
#'   \code{n_par}, \code{LRT_vs_UN_chisq}, \code{LRT_vs_UN_df},
#'   \code{LRT_vs_UN_p}.
#'
#' @details
#' \strong{Structures.} Every structure below is nested in UN, so the
#' likelihood-ratio test against UN is well defined for each.
#' \itemize{
#'   \item \code{IND}: independent observations within subject
#'         (\code{correlation = NULL} in \code{gls}). Provided as a
#'         baseline.
#'   \item \code{CS}: compound symmetry, a constant correlation and a
#'         single variance across time points:
#'         \code{nlme::corCompSymm()}.
#'   \item \code{CSH}: heterogeneous compound symmetry, a constant
#'         correlation with a separate variance at each time point:
#'         \code{nlme::corCompSymm()} with \code{nlme::varIdent()}.
#'   \item \code{AR1}: first-order autoregressive correlation with a
#'         single variance: \code{nlme::corAR1()}.
#'   \item \code{ARH1}: heterogeneous first-order autoregressive
#'         correlation with a separate variance at each time point:
#'         \code{nlme::corAR1()} with \code{nlme::varIdent()}.
#'   \item \code{TOEP}: Toeplitz (banded), a separate correlation at
#'         each lag with a single variance:
#'         \code{nlme::corARMA()} with autoregressive order one less
#'         than the number of time points and no moving-average term.
#'   \item \code{TOEPH}: heterogeneous Toeplitz, the Toeplitz
#'         correlation with a separate variance at each time point:
#'         \code{nlme::corARMA()} with \code{nlme::varIdent()}.
#'   \item \code{UN}: unstructured, every variance and covariance free:
#'         \code{nlme::corSymm()} with \code{nlme::varIdent()}.
#' }
#'
#' \strong{LRT.} Each restricted structure is compared against UN by
#' the likelihood-ratio test. Both fits are re-estimated under ML (not
#' REML) for the LRT, following \code{nlme} convention. The chi square
#' statistic is \eqn{-2 (\log L_{\mathrm{restricted}} - \log
#' L_{\mathrm{UN}})} on degrees of freedom equal to the difference in
#' parameter count.
#'
#' \strong{Caveats.} The likelihood-ratio test against UN is valid
#' because each listed structure is a restriction of UN. Two
#' structures that are not nested in each other (for example CS and
#' AR(1)) should be compared by AIC or BIC rather than by an LRT.
#'
#' @references
#' Littell, R. C., Milliken, G. A., Stroup, W. W., Wolfinger, R. D.,
#'   & Schabenberger, O. (2006). \emph{SAS for mixed models} (2nd
#'   ed.). SAS Institute.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 15.)
#'
#' Pinheiro, J. C., & Bates, D. M. (2000). \emph{Mixed-effects models
#'   in S and S-PLUS}. Springer.
#'
#' @seealso \code{\link[nlme]{gls}}, \code{\link[stats]{logLik}},
#'   \code{\link[stats]{anova}}
#'
#' @examples
#' # Four repeated measures on each of 30 subjects.
#' set.seed(113)
#' n <- 30; k <- 4
#' subj <- factor(rep(1:n, each = k))
#' tm   <- factor(rep(1:k, times = n))
#' y    <- as.vector(t(matrix(rnorm(n * k), n, k) +
#'                      rep(rnorm(n, 0, 1), each = k)))
#' d <- data.frame(y, subj, tm)
#'
#' # All eight structures at once. Read the table by comparing AIC and BIC
#' # across rows, and use the likelihood-ratio test only for the nested
#' # comparison it reports, each structure against UN.
#' compare_cov_structures(d, outcome = "y", subject = "subj",
#'                        time = "tm")
#'
#' # A subset, requested with lowercase aliases (matching is case
#' # insensitive).
#' compare_cov_structures(d, outcome = "y", subject = "subj",
#'                        time = "tm",
#'                        structures = c("cs", "csh", "ar1", "arh1"))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#'
#' @export

compare_cov_structures <- function(data, outcome, subject, time,
                                   fixed_effects = NULL,
                                   structures = c("IND", "CS", "CSH",
                                                  "AR1", "ARH1",
                                                  "TOEP", "TOEPH",
                                                  "UN")) {
  if (!requireNamespace("nlme", quietly = TRUE))
    stop("The 'nlme' package is required. Install it with: ",
         "install.packages('nlme')")
  for (nm in c("outcome", "subject", "time")) {
    v <- get(nm)
    if (!is.character(v) || !v %in% names(data))
      stop(sprintf("'%s' must be a single column name in 'data'.", nm))
  }
  choices <- c("IND", "CS", "CSH", "AR1", "ARH1", "TOEP", "TOEPH", "UN")
  structures <- match.arg(toupper(structures), choices = choices,
                          several.ok = TRUE)

  d <- data
  d[[subject]] <- factor(d[[subject]])
  d[[time]]    <- factor(d[[time]])
  if (is.null(fixed_effects))
    fixed_effects <- stats::as.formula(paste("~", time))
  rhs <- deparse(fixed_effects[[length(fixed_effects)]])
  f   <- stats::as.formula(sprintf("%s ~ %s", outcome, rhs))

  cform <- stats::as.formula(
    sprintf("~ as.integer(%s) | %s", time, subject))
  vform <- stats::as.formula(paste("~ 1 |", time))
  n_time <- nlevels(d[[time]])

  fit_one <- function(struct) {
    corr <- switch(struct,
      IND   = NULL,
      CS    = nlme::corCompSymm(form = cform),
      CSH   = nlme::corCompSymm(form = cform),
      AR1   = nlme::corAR1(form = cform),
      ARH1  = nlme::corAR1(form = cform),
      TOEP  = nlme::corARMA(form = cform, p = n_time - 1, q = 0),
      TOEPH = nlme::corARMA(form = cform, p = n_time - 1, q = 0),
      UN    = nlme::corSymm(form = cform))
    weights <- if (struct %in% c("CSH", "ARH1", "TOEPH", "UN"))
      nlme::varIdent(form = vform) else NULL
    nlme::gls(f, data = d, correlation = corr, weights = weights,
              method = "ML")
  }

  fits <- list()
  for (s in structures) {
    fits[[s]] <- tryCatch(fit_one(s), error = function(e) e)
  }

  # Build the table. UN is the reference for every LRT; when it
  # was not requested (or did not fit) the LRT columns stay NA.
  rows <- list()
  fit_un <- fits[["UN"]]
  un_ok  <- !is.null(fit_un) && !inherits(fit_un, "error")
  ll_un  <- if (un_ok) as.numeric(stats::logLik(fit_un)) else NA_real_
  npar_un <- if (un_ok) attr(stats::logLik(fit_un), "df") else NA_integer_
  for (s in structures) {
    f_s <- fits[[s]]
    if (inherits(f_s, "error")) {
      rows[[s]] <- data.frame(structure = s,
        log_lik = NA_real_, AIC = NA_real_, BIC = NA_real_,
        n_par = NA_integer_,
        LRT_vs_UN_chisq = NA_real_, LRT_vs_UN_df = NA_integer_,
        LRT_vs_UN_p = NA_real_, stringsAsFactors = FALSE)
      next
    }
    ll_s   <- as.numeric(stats::logLik(f_s))
    npar_s <- attr(stats::logLik(f_s), "df")
    aic_s  <- stats::AIC(f_s);  bic_s <- stats::BIC(f_s)
    if (s == "UN" || !un_ok) {
      lrt_chi <- NA_real_; lrt_df <- NA_integer_; lrt_p <- NA_real_
    } else {
      lrt_chi <- 2 * (ll_un - ll_s)
      lrt_df  <- npar_un - npar_s
      lrt_p   <- if (lrt_df > 0)
                   stats::pchisq(lrt_chi, lrt_df, lower.tail = FALSE)
                 else NA_real_
    }
    rows[[s]] <- data.frame(structure = s,
      log_lik = ll_s, AIC = aic_s, BIC = bic_s, n_par = npar_s,
      LRT_vs_UN_chisq = lrt_chi, LRT_vs_UN_df = lrt_df,
      LRT_vs_UN_p = lrt_p, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
