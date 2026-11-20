#' Factor Count Comparison: Model Fit, Parallel Analysis, and MAP
#'
#' Fits the exploratory common-factor model with one factor, two, and so
#' on, and sets beside each count what bears on it: the fit of that model
#' and the test of the previous count against it, its information
#' criteria, the reduced eigenvalue at that rank against its parallel
#' analysis reference, and Velicer's average squared partial correlation
#' after that many components. The table applies no cutoff and chooses
#' nothing: the lines beneath it state what each criterion shows, and the
#' count is a judgment that combines them with theory and the
#' interpretability of the rotated solutions, which
#' \code{\link{fa_exploratory}} gives at a chosen count.
#'
#' @param data Numeric raw observations. Supply exactly one of \code{data}
#'   or \code{S}.
#' @param S Named positive definite covariance matrix, with \code{N}.
#' @param N Sample size for covariance input, or a reference sample size
#'   for the parallel analysis of incomplete raw data (see
#'   \code{\link{parallel_analysis}}).
#' @param items The items to analyze, by name: columns of \code{data}, or
#'   rows and columns of \code{S}, which is then reduced to them, in the
#'   order given; only the reduced matrix need be finite and positive
#'   definite. \code{NULL}, the default, analyzes every column of
#'   \code{data} or of \code{S}.
#'   Names that cannot be read as at least three distinct items are refused
#'   with an error of class \code{dmar_items_invalid} (see
#'   \code{\link{dmar_metadata}}).
#' @param max_factors The largest count fitted. The default is the largest
#'   count whose model keeps at least one degree of freedom, the largest
#'   \eqn{k} with \eqn{(p - k)^2 - (p + k) \ge 2} for \eqn{p} items; a
#'   larger count is accepted up to the one with zero degrees of freedom,
#'   whose model reproduces the covariances exactly.
#' @param estimator,missing,ordered As for \code{\link{fa_exploratory}},
#'   which fits each count with \code{rotation = "none"}, since the fit of
#'   a count does not depend on the rotation. The missing data method is
#'   resolved once, as \code{fa_exploratory()} resolves it, and the same
#'   method serves every fit, the correlations of the eigenvalues and the
#'   MAP, and the parallel analysis. With ordered items and no item value
#'   missing, leaving \code{missing} out fits the items as they stand
#'   (listwise deletion with nothing to delete), exactly as
#'   \code{missing = "listwise"} does, and the \code{method} attribute
#'   records \code{"listwise"}; with a value missing, ordered items require
#'   an explicit \code{"pairwise"} or \code{"listwise"}, and leaving
#'   \code{missing} out is refused before any correlation, fit, or
#'   reference data set is computed.
#' @param criterion The information criterion by which the proper counts
#'   are compared, as \code{\link{compare_sem}} compares a set of models:
#'   \code{"BIC"} (the default), \code{"SABIC"}, \code{"HBIC"}, or
#'   \code{"AIC"}, in either case. IBIC and SPBIC, the information matrix
#'   based criteria of Bollen, Ray, Zavisca, and Harden (2012), are not
#'   offered: lavaan's exploratory fits carry singular information matrices
#'   (their parameterization counts the rotational constraints), which
#'   \code{\link{sem_information_criteria}} refuses for them.
#' @param conf_level Confidence level of the RMSEA intervals, default 0.90,
#'   the conventional level, whose limits decide the tests of close and not
#'   close fit at .05 (see \code{\link{sem_fit_measures}});
#'   \code{conf_level = 0.95} gives 95 percent intervals. It sets nothing
#'   else: the Monte Carlo limits of the parallel analysis reference in the
#'   \code{parallel_analysis} attribute keep the 95 percent default of
#'   \code{\link{parallel_analysis}}.
#' @param B Number of reference datasets for the parallel analysis, a
#'   single whole number, default 1000; \code{0} omits the reference, and
#'   the table then carries the reduced eigenvalues alone. A positive
#'   \code{B} must be at least the number \code{\link{parallel_analysis}}
#'   requires at \code{percentile}: 200 at the default 95th percentile,
#'   1000 at the 99th, 100 at the 90th or at the median (at least 10
#'   reference datasets beyond the quantile and never fewer than 100). A
#'   smaller positive \code{B} is refused before any model is fit (class
#'   \code{dmar_count_invalid}); see \code{\link{dmar_bootstrap}}.
#' @param percentile Reference percentile of the parallel analysis, default
#'   0.95.
#' @param seed Optional seed for the parallel analysis; default \code{NULL}
#'   uses the current random stream. A supplied seed restores the caller's
#'   stream.
#' @param x For \code{print()}, a table returned by \code{factor_count()}.
#' @param digits For \code{print()}, the decimal places of the display
#'   (default \code{getOption("dmar.digits", 3)}; see
#'   \code{\link{dmar_tbl}}).
#' @param \dots For \code{print()}, passed to \code{print.data.frame()}.
#'
#' @details
#' \strong{The fit of each count.} Each count is fitted by
#' \code{\link{fa_exploratory}} with the same estimator and missing data
#' method, and its \code{fit_measures} supply the chi square statistic
#' with its degrees of freedom and p-value, the RMSEA with its interval,
#' the CFI, TLI, and SRMR, under the same version policy (see
#' \code{\link{sem_fit_measures}}; after FIML the indices are the
#' corrected ones of Zhang and Savalei, 2023). The model with one factor
#' fewer is nested in the model with one more, and \code{delta_chi_square},
#' \code{delta_df}, and \code{p_delta} give \code{lavaan::lavTestLRT()}'s
#' test of the previous count against this one (the Satorra-Bentler scaled
#' difference under \code{estimator = "MLR"}). The AIC, BIC, and sample
#' size adjusted BIC are those of \code{\link{sem_information_criteria}};
#' smaller is preferred, and only the differences between rows carry
#' meaning. Under \code{estimator = "ULS"} the statistic has no reference
#' distribution, so the p-values, the RMSEA, CFI, TLI, difference tests,
#' and information criteria are \code{NA} and the SRMR remains (see
#' \code{\link{fa_exploratory}}).
#'
#' \strong{Differences, weights, and evidence.} The values of a criterion
#' mean nothing on their own; what carries the evidence is the difference
#' between two counts. The proper counts are therefore compared by
#' \code{criterion} as \code{\link{compare_sem}} compares a set of models
#' fitted to the same data: \code{criterion_delta} is each count's gap from
#' the lowest value, \code{criterion_weight} is \eqn{\exp(-\Delta/2)}
#' normalized over the proper counts, and, for BIC and HBIC, which
#' approximate \eqn{-2} times the log marginal likelihood (Bollen et al.,
#' 2012), \code{criterion_log_evidence} is \eqn{\Delta/2}, the approximate
#' log Bayes factor of the best count over this one, and the weights
#' approximate posterior probabilities under equal prior odds (Burnham &
#' Anderson, 2004); for AIC the weights are Akaike weights and for SABIC
#' they carry no further interpretation, so the log evidence column is
#' \code{NA} for both. The weights are conditional on the set of proper
#' counts compared, and the line under the table gives the best count, its
#' weight, and how far the nearest other count stands from it. Every count
#' is fitted to the same cases under the same missing data method and mean
#' structure, which \code{compare_sem()} checks. The comparison needs two
#' proper likelihood fits, so under ULS, or with one proper count, these
#' columns are \code{NA}.
#'
#' A count whose model does not converge, or converges to an improper
#' solution (lavaan's post-estimation check fails: a residual variance
#' below zero, or a covariance matrix that is not positive definite), or
#' whose likelihood has no maximum (a residual variance runs off below
#' zero while the likelihood is flat or still rising, and the optimizer
#' stopped somewhere along the way; see \code{\link{fa_exploratory}}),
#' has no fit to report: its \code{solution} column says so, its fit
#' columns are \code{NA}, and a line under the table names it. Such a
#' count is informative in its own right, since the unrestricted model
#' with that many factors drives an item's residual variance to or below
#' zero on these data, but it is not a stopping rule: the counts beyond it
#' are fitted and reported in the same way. \code{\link{fa_exploratory}}
#' refuses an improper solution for the same reason.
#'
#' \strong{Eigenvalues and their reference.} The \code{eigenvalue} column
#' holds the eigenvalues of the reduced correlation matrix, the item
#' correlations with the communalities of a fixed one-factor maximum
#' likelihood solution on the diagonal, as \code{\link{parallel_analysis}}
#' forms them; the \code{reference_quantile} is that function's
#' \code{percentile} of the same eigenvalue in \code{B} independent normal
#' reference datasets, and \code{exceeds_reference} compares the two (Horn,
#' 1965). The correlation matrix is the one \code{parallel_analysis()}
#' uses: the full information maximum likelihood estimate of
#' \code{\link{cov_fiml}} for incomplete continuous items, the polychoric
#' correlations for ordered items, the sample correlations otherwise. The
#' reference percentile is a null threshold, not a confidence interval for
#' an eigenvalue or a confidence statement about the count. The parallel
#' analysis is run when every item is ordered or none is, which is the
#' interface \code{parallel_analysis()} supports.
#'
#' \strong{The minimum average partial.} Velicer's (1976) rule extracts
#' principal components of the correlation matrix one at a time, partials
#' the first \eqn{m} out of it, and averages the squared partial
#' correlations that remain; no component is retained after the average
#' reaches its minimum. The \code{map} column holds the average after as
#' many components as the row has factors; the \code{map} attribute holds
#' it for every \eqn{m} from 0 (the average squared correlation itself) to
#' \eqn{p - 2}, and the line under the table names the \eqn{m} at which it
#' is smallest. The rule concerns components used as a first-stage
#' description of the common variance; it is reported beside the fit of
#' the factor models, not as a substitute for it.
#'
#' \strong{Display.} The printed table, and \code{knit_print()} and
#' \code{\link{as_kable}} of it, show one row per count with its solution
#' status, the chi square statistic, degrees of freedom, and p-value, the
#' RMSEA with its interval (whose level the line beneath the table names),
#' the CFI, TLI, and SRMR, the chosen criterion
#' with its \code{delta} and \code{weight}, the p-value of the difference
#' test against the previous count, the eigenvalue and its reference, and
#' the MAP; columns with nothing to show (the reference when \code{B = 0},
#' the likelihood quantities under ULS, the comparison with one proper
#' count) are left out. The table itself keeps every column. Beneath the
#' table: the counts with the lowest AIC, BIC, and SABIC among the proper
#' solutions, the criterion comparison, how far the eigenvalues exceed
#' their reference, where the MAP is smallest, and the counts without a
#' proper solution.
#'
#' @return A \code{dmar_tbl} data frame (class \code{dmar_factor_count})
#'   with one row per count: \code{term}, \code{factors},
#'   \code{solution} (\code{"proper"}, \code{"improper"},
#'   \code{"no maximum"}, or \code{"not converged"}), \code{chi_square},
#'   \code{df}, \code{p_value}, \code{rmsea}, \code{rmsea_ci_lower},
#'   \code{rmsea_ci_upper}, \code{cfi}, \code{tli}, \code{srmr},
#'   \code{aic}, \code{bic}, \code{sabic}, \code{hbic}, \code{criterion_delta},
#'   \code{criterion_weight}, \code{criterion_log_evidence} (the
#'   comparison by \code{criterion}), \code{delta_chi_square}, \code{delta_df},
#'   \code{p_delta}, \code{eigenvalue}, \code{reference_quantile},
#'   \code{exceeds_reference}, and \code{map}. Its \code{tidy()} keeps
#'   the count as \code{estimate} and adds a \code{conf_level} column, the
#'   level of the RMSEA intervals (see \code{\link{dmar_tidiers}}), which
#'   \code{glance()} also reports. Attributes retain
#'   \code{fits} (the lavaan object of each proper count, \code{NULL}
#'   otherwise), \code{fit_measures} (each count's table, as
#'   \code{fa_exploratory()} returns it), \code{parallel_analysis} (the
#'   \code{\link{parallel_analysis}} table, when \code{B > 0}),
#'   \code{comparison} (the \code{\link{compare_sem}} table of the proper
#'   counts, when two or more were compared), \code{map} (the average squared partial correlation for \eqn{m = 0,
#'   1, \dots, p - 2}), \code{reasons} (why a count has no proper
#'   solution), \code{method}, and \code{analysis_note}.
#' @references
#' Bollen, K. A., Ray, S., Zavisca, J., & Harden, J. J. (2012). A
#'   comparison of Bayes factor approximation methods including two new
#'   methods. \emph{Sociological Methods & Research, 41}(2), 294--324.
#'   \doi{10.1177/0049124112452393}
#'
#' Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
#'   Understanding AIC and BIC in model selection. \emph{Sociological
#'   Methods & Research, 33}(2), 261--304. \doi{10.1177/0049124104268644}
#'
#' Horn, J. L. (1965). A rationale and test for the number of factors in
#'   factor analysis. \emph{Psychometrika, 30}, 179--185.
#'   \doi{10.1007/BF02289447}
#'
#' Velicer, W. F. (1976). Determining the number of components from the
#'   matrix of partial correlations. \emph{Psychometrika, 41}(3),
#'   321--327. \doi{10.1007/BF02293557}
#'
#' Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI
#'   following FIML and TS estimation with missing data.
#'   \emph{Psychological Methods, 28}(2), 263--283.
#'   \doi{10.1037/met0000445}
#' @author Ken Kelley \email{kkelley@@nd.edu}
#' @examples
#' items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
#'            "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
#' # B = 200 is the fewest reference data sets the 95th percentile of the
#' # parallel analysis is reported from; the default is B = 1000.
#' factor_count(holzinger_swineford, items = items, B = 200, seed = 113)
#' @seealso \code{\link{fa_exploratory}} for the solution at a chosen
#'   count, \code{\link{parallel_analysis}} for the eigenvalue reference
#'   with its Monte Carlo limits, \code{\link{compare_sem}} for the
#'   comparison the criterion columns reproduce, and
#'   \code{\link{sem_information_criteria}}.
#' @export
factor_count <- function(data = NULL, S = NULL, N = NULL, items = NULL,
                         max_factors = NULL, estimator = "ML",
                         missing = "fiml", ordered = NULL,
                         criterion = "BIC", conf_level = 0.90, B = 1000L,
                         percentile = 0.95, seed = NULL) {
  # Whether 'missing' was left out, read before match.arg() gives it the
  # value of the formal default: the default is resolved as
  # fa_exploratory() resolves it (.efa_resolve_missing()).
  missing_omitted <- base::missing(missing)
  .fa_require_lavaan()
  .fa_level(conf_level)
  .fa_level(percentile, "percentile")
  .dmar_local_seed(seed)
  if (is.null(data) == is.null(S)) {
    .dmar_refuse_argument(
      "Supply exactly one of 'data' or 'S'.", arg = c("data", "S"))
  }
  estimator <- match.arg(estimator, c("ML", "MLR", "ULS", "WLSMV"))
  missing <- match.arg(missing, c("fiml", "ml", "listwise", "pairwise"))
  criterion <- match.arg(toupper(criterion), c("BIC", "SABIC", "HBIC", "AIC"))
  # B = 0 omits the parallel analysis. Otherwise its reference quantile
  # is read from B reference data sets and B follows the rule of
  # ?dmar_bootstrap: at least 10 of them beyond the quantile and never
  # fewer than 100 (200 at the 95th percentile).
  B <- .dmar_check_B(B, tail = min(percentile, 1 - percentile),
                     zero_ok = TRUE, what = "reference", position = percentile)
  ord <- character()
  if (!is.null(data)) {
    data <- .fa_data(data, items)
    items <- names(data)
    if (isTRUE(ordered)) ord <- items else if (!is.null(ordered)) ord <- ordered
    if (!is.character(ord) || anyNA(ord) || any(!ord %in% items) ||
        anyDuplicated(ord)) .dmar_refuse_argument(
          "'ordered' must name analysis items.", arg = "ordered")
    # The missing data method of every fit, of the correlations of the
    # eigenvalues and the MAP, and of the parallel analysis, resolved once
    # as fa_exploratory() resolves it, before anything is computed: with
    # ordered items and no value missing, the default fits the items as
    # they stand (listwise deletion with nothing to delete), and with a
    # value missing, pairwise or listwise deletion must be chosen. Earlier
    # builds of 1.1.0 passed the formal default "fiml" on explicitly, and
    # complete ordered items were refused.
    missing <- .efa_resolve_missing(missing, missing_omitted, data, ord,
                                    estimator)
    used <- if (missing == "listwise") {
      data[stats::complete.cases(data), , drop = FALSE]
    } else data
    used <- used[rowSums(!is.na(used)) > 0L, , drop = FALSE]
    n_cases <- nrow(used)
    # The correlation matrix of the eigenvalues and the MAP, formed as
    # parallel_analysis() forms it.
    R <- if (length(ord)) {
      if (setequal(ord, items)) .pa_cor(used, ordered = TRUE) else {
        as.matrix(lavaan::lavCor(used, ordered = ord, missing = "pairwise",
                                 output = "cor"))
      }
    } else if (anyNA(used)) {
      stats::cov2cor(cov_fiml(used)$S)
    } else stats::cor(used)
  } else {
    # Named items reduce S to their rows and columns, in the order given,
    # before it is checked, as they select the columns of raw data; the
    # fits, the eigenvalues, and the MAP all read the reduced matrix.
    S <- .fa_covariance(S, N, items)
    items <- colnames(S)
    n_cases <- N
    R <- stats::cov2cor(S)
  }
  p <- length(items)
  df_of <- function(k) ((p - k)^2 - (p + k)) / 2
  with_df <- which(vapply(seq_len(p), function(k) df_of(k) >= 1, logical(1)))
  identified <- which(vapply(seq_len(p), function(k) df_of(k) >= 0, logical(1)))
  if (is.null(max_factors)) {
    if (!length(with_df)) {
      .dmar_refuse_argument(c(
        "No factor count leaves the exploratory model a degree of freedom ",
        "with ", p, " items; the comparison needs at least four."), arg = NULL)
    }
    max_factors <- max(with_df)
  }
  if (!is.numeric(max_factors) || length(max_factors) != 1L ||
      !is.finite(max_factors) || max_factors < 1 || max_factors %% 1 ||
      !max_factors %in% identified) {
    .dmar_refuse_argument(c(
      "'max_factors' must be a whole number of at least one whose ",
      "exploratory model is identified for these ", p, " items (the ",
      "largest is ", if (length(identified)) max(identified) else 0L, ")."),
      arg = "max_factors")
  }
  max_factors <- as.integer(max_factors)
  counts <- seq_len(max_factors)

  # ---- the fit of each count ----------------------------------------
  fits <- measures <- vector("list", max_factors)
  status <- rep("proper", max_factors)
  fit_missing <- NULL
  reasons <- rep(NA_character_, max_factors)
  h1_notes <- character()
  for (k in counts) {
    notices <- new.env(parent = emptyenv())
    notices$qualified <- character()
    res <- withCallingHandlers(
      tryCatch(
        fa_exploratory(data = data, nfactors = k, S = S, N = N, items = items,
                       estimator = estimator, missing = missing,
                       ordered = if (length(ord)) ord else NULL,
                       rotation = "none", rmsea_conf_level = conf_level),
        error = function(e) {
          # A count without a proper solution is a row of the table; any
          # other error concerns the input and stops the comparison.
          if (inherits(e, "dmar_efa_no_solution")) e else stop(e)
        }),
      warning = function(w) {
        if (inherits(w, "dmar_h1_unavailable")) {
          notices$qualified <- c(notices$qualified, conditionMessage(w))
        }
        invokeRestart("muffleWarning")
      })
    if (inherits(res, "error")) {
      status[k] <- res$reason
      # A count whose likelihood has no maximum keeps the short reason; its
      # refusal also prints the lavaan code of the fit, which the call
      # that wants it (fa_exploratory() with that count) shows.
      reasons[k] <- if (inherits(res, "dmar_no_maximum")) {
        .no_max_short(.no_max_fields(res))
      } else {
        conditionMessage(res)
      }
    } else {
      fits[[k]] <- attr(res, "fit", exact = TRUE)
      measures[[k]] <- attr(res, "fit_measures", exact = TRUE)
      if (is.null(fit_missing)) fit_missing <- attr(res, "method", exact = TRUE)$missing
      if (length(notices$qualified)) {
        h1_notes <- c(h1_notes, paste0(.fc_count_word(k), ": ", notices$qualified[1L]))
      }
    }
  }
  pick <- function(term) {
    vapply(measures, function(fm) {
      if (is.null(fm)) return(NA_real_)
      v <- fm$value[match(term, fm$term)]
      if (length(v) != 1L) NA_real_ else v
    }, numeric(1))
  }
  likelihood <- estimator %in% c("ML", "MLR")
  criteria <- lapply(fits, function(f) {
    if (is.null(f) || !likelihood) return(NULL)
    tryCatch(suppressWarnings(suppressMessages(sem_information_criteria(f))),
             error = function(e) NULL)
  })
  ic_value <- function(term) {
    vapply(criteria, function(ic) {
      if (is.null(ic)) NA_real_ else ic$value[match(term, ic$term)]
    }, numeric(1))
  }
  # The comparison of the proper counts by the chosen criterion, as
  # compare_sem() draws it: each count's gap from the lowest value, its
  # relative weight within the set, and (for BIC and HBIC) the approximate
  # log Bayes factor of the best count over it. The fits share the data,
  # the missing data method, and the mean structure, which compare_sem()
  # checks. It needs two proper likelihood fits.
  comparison <- NULL
  crit_delta <- crit_weight <- crit_evidence <- rep(NA_real_, max_factors)
  has_fit <- !vapply(fits, is.null, logical(1))
  if (likelihood && sum(has_fit) >= 2L) {
    proper_fits <- stats::setNames(fits[has_fit], paste0("factors_", counts[has_fit]))
    comparison <- tryCatch(suppressWarnings(suppressMessages(
      compare_sem(proper_fits, criterion = criterion))), error = function(e) NULL)
    if (!is.null(comparison)) {
      idx <- match(paste0("factors_", counts), comparison$term)
      crit_delta <- comparison$delta[idx]
      crit_weight <- comparison$relative_weight[idx]
      crit_evidence <- comparison$log_evidence_best_over_model[idx]
    }
  }
  # The test of the previous count against this one: lavaan's difference
  # test of two nested fits (scaled under MLR). The ULS statistic has no
  # reference distribution, so no difference test.
  delta_chi <- delta_df <- p_delta <- rep(NA_real_, max_factors)
  if (estimator != "ULS") {
    for (k in counts[-1L]) {
      if (is.null(fits[[k - 1L]]) || is.null(fits[[k]])) next
      lrt <- tryCatch(suppressWarnings(suppressMessages(
        lavaan::lavTestLRT(fits[[k - 1L]], fits[[k]]))), error = function(e) NULL)
      if (!is.null(lrt) && nrow(lrt) == 2L) {
        delta_chi[k] <- lrt[2L, "Chisq diff"]
        delta_df[k] <- lrt[2L, "Df diff"]
        p_delta[k] <- lrt[2L, "Pr(>Chisq)"]
      }
    }
  }

  # ---- eigenvalues against their reference, and the MAP ----------------
  eigenvalues <- tryCatch(.pa_reduced_eigen(R, max(3L, n_cases))$values,
                          error = function(e) rep(NA_real_, p))
  reference <- rep(NA_real_, p)
  exceeds <- rep(NA, p)
  pa <- NULL
  pa_notices <- new.env(parent = emptyenv())
  pa_notices$messages <- character()
  pa_skipped <- NULL
  if (B > 0L) {
    if (length(ord) && !setequal(ord, items)) {
      pa_skipped <- paste("Parallel analysis: not run (its ordered interface",
                          "takes every item ordered or none).")
    } else {
      # The missing data method resolved above, which parallel_analysis()
      # takes explicitly for ordered items: "listwise" for complete
      # ordered items with 'missing' left out.
      pa <- withCallingHandlers(
        parallel_analysis(data = data, S = S, N = N, items = items,
                          missing = missing,
                          ordered = if (length(ord)) TRUE else NULL,
                          reference = "normal", B = B,
                          percentile = percentile),
        warning = function(w) {
          # The mark of a run below the required number of reference data
          # sets (options(dmar.bootstrap_exploratory = TRUE)) reaches the
          # caller as its own classed warning; the other notices of the
          # parallel analysis are kept for the notes under the table.
          # It is given under this function's name, the one the user
          # called.
          if (inherits(w, "dmar_bootstrap_exploratory")) {
            w$message <- sub("parallel_analysis(): ", "factor_count(): ",
                             conditionMessage(w), fixed = TRUE)
            warning(w)
            invokeRestart("muffleWarning")
          }
          pa_notices$messages <- c(pa_notices$messages, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      eigenvalues <- pa$estimate
      reference <- pa$reference_quantile
      exceeds <- pa$exceeds_reference
    }
  }
  map <- .fc_map(R)

  # ---- the table -------------------------------------------------------
  out <- data.frame(
    term = paste0("factors_", counts), factors = counts, solution = status,
    chi_square = pick("chi_square"), df = pick("df"),
    p_value = pick("p_chi_square"), rmsea = pick("rmsea"),
    rmsea_ci_lower = pick("rmsea_ci_lower"),
    rmsea_ci_upper = pick("rmsea_ci_upper"),
    cfi = pick("cfi"), tli = pick("tli"), srmr = pick("srmr"),
    aic = ic_value("AIC"), bic = ic_value("BIC"), sabic = ic_value("SABIC"),
    hbic = ic_value("HBIC"), criterion_delta = crit_delta,
    criterion_weight = crit_weight, criterion_log_evidence = crit_evidence,
    delta_chi_square = delta_chi, delta_df = delta_df, p_delta = p_delta,
    eigenvalue = eigenvalues[counts], reference_quantile = reference[counts],
    exceeds_reference = exceeds[counts], map = unname(map[counts + 1L]),
    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  attr(out, "fits") <- fits
  attr(out, "fit_measures") <- measures
  attr(out, "parallel_analysis") <- pa
  # The counts of the reference data sets, in the names every result of
  # the package that ran replications carries (?dmar_bootstrap).
  if (!is.null(pa)) {
    for (nm in c("B", "B_used", "n_dropped", "B_required",
                 "bootstrap_status")) {
      attr(out, nm) <- attr(pa, nm, exact = TRUE)
    }
  }
  attr(out, "comparison") <- comparison
  attr(out, "map") <- map
  attr(out, "reasons") <- stats::setNames(reasons, paste0("factors_", counts))
  attr(out, "method") <- list(
    N = n_cases, p = p, estimator = estimator, criterion = criterion,
    missing = if (!is.null(fit_missing)) fit_missing else if (is.null(data)) {
      "covariance input"
    } else if (missing == "ml") "fiml" else missing,
    ordered = ord, max_factors = max_factors, B = B, percentile = percentile,
    suggested_factors = if (is.null(pa)) NA_integer_ else attr(pa, "suggested_factors", exact = TRUE),
    reference_N = if (is.null(pa)) NA_integer_ else attr(pa, "method", exact = TRUE)$N,
    map_minimum = unname(which.min(map)) - 1L,
    notes = c(pa_skipped, pa_notices$messages, h1_notes))
  attr(out, "analysis_note") <- paste0(
    "No cutoff is applied and no count is chosen: the count is a judgment ",
    "that combines the fit of each model, the difference tests, the ",
    "information criteria, the eigenvalues against their reference, and ",
    "the MAP with theory and the interpretability of the rotated solutions, ",
    "which fa_exploratory() gives at a chosen count. Each model's interval ",
    "and p-values condition on its own count, and the criterion weights are ",
    "conditional on the set of proper counts compared.",
    if (estimator == "ULS") paste0(
      " Under ULS the statistic has no reference distribution, so the ",
      "p-values, RMSEA, CFI, TLI, difference tests, and information ",
      "criteria are NA (see fa_exploratory())."))
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_factor_count")
}

#' @rdname factor_count
#' @export
print.dmar_factor_count <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  if (!.fc_display_ok(x)) return(invisible(NextMethod()))
  .dmar_tbl_write(x, digits = digits, print_args = c(list(right = TRUE), list(...)),
                  view_args = list(digits = digits))
}

# ---- internals ----------------------------------------------------------

# Velicer's (1976) average squared partial correlation after m principal
# components of the correlation matrix R, for m = 0 (the average squared
# correlation itself) to p - 2: the first m components are partialed out
# (C = R - A_m A_m', with A_m the component loadings), the partial
# correlation matrix is C rescaled to unit diagonal, and the average is
# over its off-diagonal elements. NA where a diagonal element of C is not
# positive. The rule retains the m at which the average is smallest.
.fc_map <- function(R) {
  p <- ncol(R)
  e <- eigen(R, symmetric = TRUE)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), p)
  average <- function(M) sum(M[upper.tri(M)]^2) * 2 / (p * (p - 1))
  m_max <- max(0L, p - 2L)
  out <- stats::setNames(rep(NA_real_, m_max + 1L), paste0("m_", 0:m_max))
  out[1L] <- average(R)
  for (m in seq_len(m_max)) {
    C <- R - tcrossprod(A[, seq_len(m), drop = FALSE])
    d <- diag(C)
    if (any(!is.finite(d)) || any(d <= 0)) next
    out[m + 1L] <- average(C / sqrt(outer(d, d)))
  }
  out
}

.fc_count_word <- function(k) paste(k, if (k == 1L) "factor" else "factors")

# Whether a table still carries what the display reads; a subset prints as
# an ordinary DMAR table.
.fc_display_ok <- function(x) {
  m <- attr(x, "method", exact = TRUE)
  is.list(m) && !is.null(m$p) && is.character(m$criterion) &&
    all(c("factors", "solution", "chi_square", "df", "p_value", "rmsea",
          "rmsea_ci_lower", "rmsea_ci_upper", "cfi", "tli", "srmr", "aic",
          "bic", "sabic", "hbic", "criterion_delta", "criterion_weight",
          "p_delta", "eigenvalue", "reference_quantile", "map") %in%
          names(x))
}

# The tidy() method for factor_count(), documented with the other DMAR
# tidiers at ?dmar_tidiers (R/dmar_tidiers.R). The table holds each count's
# RMSEA interval in rmsea_ci_lower and rmsea_ci_upper beside its other fit
# measures. The generic wide tidier does not read those columns as the
# table's interval, and should not, since the row's estimate is the count,
# so the level of the interval, which the table records only as an
# attribute, is added here as a conf_level column, as glance() gives it.

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_factor_count <- function(x, ...) {
  out <- NextMethod()
  if (is.null(out$conf_level)) {
    level <- attr(x, "conf_level", exact = TRUE)
    out$conf_level <- rep(if (is.numeric(level) && length(level) == 1L)
      level else NA_real_, nrow(out))
  }
  out
}

#' @exportS3Method .dmar_tbl_display dmar_factor_count
#' @noRd
.dmar_tbl_display.dmar_factor_count <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  if (!.fc_display_ok(x)) return(NextMethod())
  # Every column but df is measured, so a CFI of 1 or a weight of 1 reads
  # 1.000 beside the others, as a fit index or a weight does in any DMAR
  # table (.dmar_measured_cells() in R/dmar_tbl.R).
  num <- function(v, whole = FALSE) {
    out <- .format_dmar_value(v, digits, whole = whole)
    out[is.na(v)] <- "NA"
    out
  }
  pv <- function(v) ifelse(is.na(v), "NA", .format_dmar_pvalue(v))
  # The criteria to one decimal, since differences of a few units between
  # rows are what they are read for.
  crit <- function(v) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = 1))
  crit_name <- tolower(attr(x, "method", exact = TRUE)$criterion)
  value <- x[[crit_name]]
  disp <- data.frame(
    factors = x$factors, solution = x$solution,
    chi_square = num(x$chi_square), df = num(x$df, whole = TRUE),
    p_value = pv(x$p_value),
    rmsea = .efa_cells(x$rmsea, x$rmsea_ci_lower, x$rmsea_ci_upper, digits),
    cfi = num(x$cfi), tli = num(x$tli), srmr = num(x$srmr),
    stringsAsFactors = FALSE)
  # The chosen criterion under its own name, then its gap from the lowest
  # and the relative weight of the count.
  disp[[crit_name]] <- crit(value)
  disp$delta <- crit(x$criterion_delta)
  disp$weight <- num(x$criterion_weight)
  disp$p_delta <- pv(x$p_delta)
  disp$eigenvalue <- num(x$eigenvalue)
  disp$reference <- num(x$reference_quantile)
  disp$map <- num(x$map)
  # Columns with nothing to show (the reference without a parallel
  # analysis, the likelihood quantities under ULS, the comparison with a
  # single proper count) are left out.
  empty <- vapply(names(x), function(nm) all(is.na(x[[nm]])), logical(1))
  drop <- c(if (empty[["p_value"]]) "p_value", if (empty[["rmsea"]]) "rmsea",
            if (empty[["cfi"]]) "cfi", if (empty[["tli"]]) "tli",
            if (all(is.na(value))) crit_name,
            if (empty[["criterion_delta"]]) c("delta", "weight"),
            if (empty[["p_delta"]]) "p_delta",
            if (empty[["reference_quantile"]]) "reference")
  disp <- disp[setdiff(names(disp), drop)]
  .as_dmar_tbl(disp)
}

#' @exportS3Method .dmar_tbl_footer dmar_factor_count
#' @noRd
.dmar_tbl_footer.dmar_factor_count <- function(x, digits = getOption("dmar.digits", 3L), ...) {
  if (!.fc_display_ok(x)) return(NextMethod())
  m <- attr(x, "method", exact = TRUE)
  level <- .format_conf_pct(attr(x, "conf_level", exact = TRUE))
  lines <- .dmar_footer_add(character(0), c(
    sprintf("Factor counts 1 to %d for %d items, N = %d", max(x$factors), m$p, m$N),
    sprintf("Estimator: %s; missing data: %s; unrotated solutions (the fit of a count does not depend on the rotation)",
            .efa_estimator_label(m), .efa_missing_label(m))), "header")
  has_reference <- !all(is.na(x$reference_quantile))
  value <- x[[tolower(m$criterion)]]
  lines <- .dmar_footer_add(lines, paste0(
    "(each count's model: ",
    if (identical(m$estimator, "ULS")) "ULS statistic (N - 1 times the discrepancy; no reference distribution)" else "chi square test",
    if (!all(is.na(x$rmsea))) sprintf(", RMSEA with %s%% interval", level),
    if (!all(is.na(x$cfi))) ", CFI, TLI", ", SRMR",
    if (!all(is.na(value))) paste0(", ", m$criterion),
    if (!all(is.na(x$criterion_delta))) paste(
      " with delta (its gap from the lowest) and weight (its relative weight",
      "among the proper solutions)"),
    if (!all(is.na(x$p_delta))) "; p_delta: the test of the previous count against this one",
    "; eigenvalue: of the reduced correlation matrix",
    if (has_reference) sprintf(", with its parallel analysis %sth percentile", format(100 * m$percentile)),
    "; map: Velicer's average squared partial correlation after this many components)"),
    "table_note")
  proper <- x$solution == "proper"
  summary <- character()
  if (any(proper & is.finite(x$aic))) {
    lowest <- function(v) x$factors[proper][which.min(v[proper])]
    summary <- c(summary, sprintf(
      "Information criteria over the proper solutions: lowest AIC at %s, lowest BIC at %s, lowest SABIC at %s.",
      .fc_count_word(lowest(x$aic)), .fc_count_word(lowest(x$bic)),
      .fc_count_word(lowest(x$sabic))))
  }
  if (any(is.finite(x$criterion_delta))) {
    ok <- which(is.finite(x$criterion_delta))
    best <- ok[which.min(x$criterion_delta[ok])]
    others <- setdiff(ok, best)
    line <- sprintf("%s comparison over the proper solutions (compare_sem()): lowest at %s, weight %s",
                    m$criterion, .fc_count_word(x$factors[best]),
                    .format_dmar_value(x$criterion_weight[best], digits,
                                       whole = FALSE))
    if (length(others)) {
      runner <- others[which.min(x$criterion_delta[others])]
      # The delta and the log Bayes factor, which is half the delta, to one
      # decimal, as the criterion columns of the table are, so that the two
      # are not shown at two precisions.
      line <- paste0(line, sprintf(
        "; the nearest other count, %s, is %s higher (weight %s%s)",
        .fc_count_word(x$factors[runner]),
        .dmar_fixed_text(x$criterion_delta[runner], 1L),
        .format_dmar_value(x$criterion_weight[runner], digits, whole = FALSE),
        if (is.finite(x$criterion_log_evidence[runner])) sprintf(
          ", log Bayes factor %s for %s over it",
          .dmar_fixed_text(x$criterion_log_evidence[runner], 1L),
          .fc_count_word(x$factors[best])) else ""))
    }
    summary <- c(summary, paste0(line, "."))
  }
  if (has_reference) {
    s <- m$suggested_factors
    summary <- c(summary, sprintf(
      "Parallel analysis (normal reference, B = %d, %sth percentile, reference N = %d): %s",
      m$B, format(100 * m$percentile), m$reference_N,
      if (is.na(s)) "thresholds withheld (see the notes)."
      else if (s == 0L) "no reduced eigenvalue exceeds its reference."
      else sprintf("the reduced eigenvalues exceed their reference through rank %d%s.", s,
                   if (s > max(x$factors)) sprintf(" (the table shows ranks 1 to %d)", max(x$factors)) else "")))
  }
  map <- attr(x, "map", exact = TRUE)
  if (!is.null(map) && any(is.finite(map))) {
    at <- m$map_minimum
    summary <- c(summary, if (at == 0L) sprintf(
      "Velicer's MAP: the average squared partial correlation is smallest before any component is removed (%s), so the rule retains none.",
      .format_dmar_value(map[1L], digits, whole = FALSE)) else sprintf(
      "Velicer's MAP: the average squared partial correlation is smallest after %d %s (%s; before any component is removed it is %s).",
      at, if (at == 1L) "component" else "components",
      .format_dmar_value(map[at + 1L], digits, whole = FALSE),
      .format_dmar_value(map[1L], digits, whole = FALSE)))
  }
  if (any(!proper)) {
    bad <- which(!proper)
    why <- c("not converged" = "the optimizer did not converge",
             improper = paste("improper: a residual variance below zero, or a",
                              "covariance matrix that is not positive definite"),
             "no maximum" = paste("no maximum: a residual variance runs",
                                  "off below zero while the likelihood",
                                  "is flat or still rising"),
             "standard errors unavailable" = "standard errors unavailable")
    summary <- c(summary, paste0(
      "Without a proper solution, so their fit is not reported: ",
      paste(sprintf("%s (%s)", vapply(bad, .fc_count_word, character(1)),
                    unname(why[x$solution[bad]])), collapse = "; "), "."))
  }
  if (length(m$notes)) summary <- c(summary, m$notes)
  lines <- .dmar_footer_add(lines, summary, "summary")
  c(lines, .dmar_footer_add(NextMethod(), attr(x, "analysis_note", exact = TRUE), "analysis"))
}
