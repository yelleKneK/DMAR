#' Parallel Analysis of Common-Factor Eigenvalues
#'
#' Compares factor-analytic eigenvalues with independently generated null
#' references. Each correlation matrix is reduced by replacing its diagonal
#' with the communalities of a fixed one-factor ML solution. The procedure
#' concerns common factors, not principal components. The reference
#' percentile is a null threshold; it is not a sampling confidence interval
#' for an observed eigenvalue or a confidence statement about a factor count.
#'
#' @param data Numeric raw observations, with \code{items} selecting columns.
#'   Supply exactly one of \code{data} or \code{S}.
#' @param S Named positive definite covariance matrix, with \code{N}.
#'   Available only for a normal-reference analysis of continuous items.
#' @param N Reference sample size for covariance input. For raw incomplete
#'   data it may be supplied as an explicit benchmark. With complete data
#'   the default is the number of cases; with FIML it is the complete-case
#'   count, and the output explicitly identifies this benchmark.
#' @param items The items to analyze, by name: columns of \code{data}, or
#'   rows and columns of \code{S}, which is then reduced to them, in the
#'   order given; only the reduced matrix need be finite and positive
#'   definite. \code{NULL}, the default, analyzes every column of
#'   \code{data} or of \code{S}.
#'   Names that cannot be read as at least three distinct items are refused
#'   with an error of class \code{dmar_items_invalid} (see
#'   \code{\link{dmar_metadata}}).
#' @param missing \code{"fiml"} (default), \code{"ml"},
#'   \code{"listwise"}, or \code{"pairwise"}. Continuous FIML obtains
#'   the unrestricted correlation matrix from \code{\link{cov_fiml}}.
#'   Ordered items follow the rule of \code{\link{fa_exploratory}} and
#'   \code{\link{factor_count}}: with no value missing, leaving
#'   \code{missing} out analyzes the items as they stand (listwise
#'   deletion with nothing to delete, recorded as \code{"listwise"}); with
#'   a value missing, the choice of \code{"pairwise"} or
#'   \code{"listwise"} is the user's, and leaving \code{missing} out is
#'   refused; giving \code{"fiml"} or \code{"ml"} for ordered items is
#'   refused whether or not a value is missing. Each refusal has class
#'   \code{dmar_argument_invalid}, as in those functions, and comes before
#'   anything is computed.
#' @param ordered \code{NULL}, \code{TRUE}, or names of ordered items.
#'   The current ordered interface requires all analysis items to be ordered.
#' @param reference \code{"normal"} (default) generates independent normal
#'   responses. With ordered items these are categorized using the observed
#'   univariate category probabilities, then polychoric correlations are
#'   estimated. \code{"resample"} samples each item's observed responses
#'   independently with replacement, preserving its empirical marginal
#'   distribution in expectation; ordered correlations remain polychoric.
#' @param B Number of independent reference datasets, a single whole
#'   number, default 1000 (each one is generated and fit by a one-factor
#'   model). The fewest accepted depends on \code{percentile}: the
#'   reference quantile is reported only when at least 10 reference
#'   datasets lie beyond it and never from fewer than 100, which takes 200
#'   at the default 95th percentile, 1000 at the 99th, and 100 at the 90th
#'   or at the median. A smaller \code{B} is refused before anything is
#'   computed (class \code{dmar_count_invalid}); see \emph{Details} and
#'   \code{\link{dmar_bootstrap}}.
#' @param percentile Null reference percentile, default 0.95.
#' @param conf_level Confidence level for Monte Carlo uncertainty in the
#'   reference percentile, default 0.95. This is not observed-data sampling
#'   uncertainty.
#' @param seed Optional seed, default \code{NULL}. A supplied seed restores
#'   the caller's random stream on exit.
#'
#' @details
#' The fixed communality reduction is computed by \code{stats::factanal()}
#' with its lower uniqueness bound of 0.005. Observed and null matrices
#' use exactly the same reduction. Defined boundary solutions are retained
#' and counted; failed fits are counted separately.
#'
#' \strong{Reference datasets.} The threshold of each rank is the
#' \code{percentile} quantile of the eigenvalues of the \code{B}
#' reference datasets, so it is located by the reference datasets that lie
#' beyond it. The function runs only when \code{B} puts at least 10 of
#' them beyond the quantile, and never with fewer than 100 (200 at the
#' 95th percentile). A reference dataset whose fit fails is dropped and
#' counted, never replaced. The count is held to the number requested.
#' Thresholds and a factor count are reported only when no more than half of
#' the attempted reference datasets failed; otherwise they are \code{NA},
#' the observed eigenvalues are returned,
#' and the warning has class \code{dmar_bootstrap_insufficient}. When
#' some failed and the thresholds were formed the warning has class
#' \code{dmar_bootstrap_replications_dropped}. The result records
#' \code{B}, \code{B_used}, \code{n_dropped}, \code{B_required}, and
#' \code{bootstrap_status} as attributes beside the fuller \code{counts}
#' list. \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every
#' such call gives a warning of class \code{dmar_bootstrap_exploratory}
#' and the line under the table says so. These rules limit the error of
#' drawing a finite number of reference datasets; they are not a coverage
#' guarantee.
#'
#' The suggested factor count is the number of consecutive ranks, starting
#' at one, whose observed eigenvalues exceed their reference thresholds.
#' It is evidence to combine with theory, model fit, and interpretability.
#' Selection itself has no confidence interval here. Monte Carlo limits
#' for the reference quantile use binomial order-statistic bounds; the
#' quantile estimate uses R's type 7 interpolation. When reference fits fail,
#' these quantities describe the null distribution conditional on returning
#' a defined fit. Discarding failures can select a different distribution;
#' the half-return rule does not repair that selection or establish the
#' intended unconditional null quantile.
#'
#' An incomplete-data FIML correlation matrix uses all observed values,
#' whereas the complete-case reference size is only an explicitly labeled
#' benchmark. It is not an effective sample size and is not guaranteed to
#' be conservative. These references do not reproduce a study's missingness
#' mechanism. They do not establish general MAR calibration. If fewer than
#' three complete cases remain, choose a justified reference \code{N}
#' explicitly rather than silently substituting the total sample size.
#'
#' @return A \code{dmar_tbl} with ranks, observed \code{estimate},
#'   \code{reference_quantile}, \code{reference_mean},
#'   \code{reference_quantile_ci_lower}, \code{reference_quantile_ci_upper}
#'   (the Monte Carlo confidence limits for the reference quantile), and
#'   \code{exceeds_reference}. Attributes retain \code{suggested_factors},
#'   \code{reference_draws}, \code{counts}, \code{observed_boundary},
#'   \code{method}, and an explicit \code{analysis_note}, with the counts
#'   every result of the package that ran replications carries: \code{B},
#'   \code{B_used}, \code{n_dropped}, \code{B_required}, and
#'   \code{bootstrap_status} (see \code{\link{dmar_bootstrap}}). Observed sampling
#'   confidence limits are unavailable and are not inferred from the null
#'   distribution.
#' @references
#' Horn, J. L. (1965). A rationale and test for the number of factors in
#'   factor analysis. \emph{Psychometrika, 30}, 179--185.
#'   \doi{10.1007/BF02289447}
#' @author Ken Kelley \email{kkelley@@nd.edu}
#' @examples
#' set.seed(113)
#' items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
#'            "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
#' # B = 200 is the fewest reference data sets the 95th percentile is
#' # reported from and keeps the example quick; the default is B = 1000.
#' parallel_analysis(holzinger_swineford, items = items, B = 200, seed = 113)
#' @export
parallel_analysis <- function(data = NULL, S = NULL, N = NULL, items = NULL,
                              missing = "fiml", ordered = NULL,
                              reference = "normal", B = 1000L,
                              percentile = 0.95, conf_level = 0.95,
                              seed = NULL) {
  # Whether 'missing' was left out, read before match.arg() gives it the
  # value of the formal default: ordered items resolve the default as
  # fa_exploratory() and factor_count() resolve it (.efa_resolve_missing()).
  missing_omitted <- base::missing(missing)
  .fa_level(conf_level)
  .fa_level(percentile, "percentile")
  .dmar_local_seed(seed)
  reference <- match.arg(reference, c("normal", "resample"))
  missing <- match.arg(missing, c("fiml", "ml", "listwise", "pairwise"))
  # The reference quantile is read from the B reference data sets, so B
  # follows the rule of ?dmar_bootstrap: at least 10 of them beyond the
  # quantile and never fewer than 100 (200 at the 95th percentile).
  boot_tail <- min(percentile, 1 - percentile)
  B <- .dmar_check_B(B, tail = boot_tail, what = "reference",
                     position = percentile)
  if (is.null(data) == is.null(S)) .dmar_refuse_argument(
    "Supply exactly one of 'data' or 'S'.", arg = c("data", "S"))
  ord <- character()
  benchmark <- FALSE
  if (!is.null(data)) {
    data <- .fa_data(data, items)
    items <- names(data)
    if (isTRUE(ordered)) ord <- items else if (!is.null(ordered)) ord <- ordered
    if (!is.character(ord) || anyNA(ord) || anyDuplicated(ord) ||
        (length(ord) && !setequal(ord, items))) {
      .dmar_refuse_argument(
        "'ordered' must name all analysis items for this interface.",
        arg = "ordered")
    }
    # Ordered items: the missing data method resolved and checked by the
    # rule of fa_exploratory() and factor_count() (.efa_resolve_missing(),
    # R/fa_exploratory_internals.R), whose polychoric correlations this
    # analysis shares (estimator "WLSMV" names that ordered-data fit).
    # With no value missing, the default analyzes the items as they stand
    # (listwise deletion with nothing to delete); with a value missing,
    # pairwise or listwise deletion must be chosen. Their refusal, the only
    # one the ordered branch of that rule gives, is signaled again with the
    # same class and a sentence for this function, which has no
    # 'estimator' argument. Through the first builds of 1.1.0 complete
    # ordered items with 'missing' left out were refused here, while
    # fa_exploratory() and factor_count() fitted them.
    if (length(ord)) {
      missing <- tryCatch(
        .efa_resolve_missing(missing, missing_omitted, data, ord, "WLSMV"),
        dmar_argument_invalid = function(e) .dmar_refuse_argument(c(
          "Ordered parallel analysis requires, when item values are ",
          "missing, an explicit missing = \"pairwise\" or \"listwise\" ",
          "choice, and does not take missing = \"fiml\" or \"ml\": ",
          "continuous normal-theory FIML is not an ordered-data ",
          "likelihood."), arg = "missing"))
    }
    if (!length(ord) && missing == "pairwise") {
      .dmar_refuse_argument(
        "Continuous parallel analysis uses FIML or listwise deletion.",
        arg = NULL)
    }
    if (missing == "listwise") data <- data[stats::complete.cases(data), , drop = FALSE]
    incomplete <- anyNA(data)
    benchmark <- incomplete || (!is.null(N) && !isTRUE(N == nrow(data)))
    if (is.null(N)) N <- sum(stats::complete.cases(data))
    if (length(ord)) {
      R <- .pa_cor(data, ordered = TRUE)
    } else if (incomplete) {
      R <- stats::cov2cor(cov_fiml(data)$S)
    } else R <- stats::cor(data)
  } else {
    # Named items reduce S to their rows and columns, in the order given,
    # before it is checked, as they select the columns of raw data.
    S <- .fa_covariance(S, N, items)
    if (reference != "normal" || !is.null(ordered)) {
      .dmar_refuse_argument(
        "Covariance input supports the continuous normal reference only.",
        arg = NULL)
    }
    items <- colnames(S)
    R <- stats::cov2cor(S)
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 3 || N %% 1) {
    .dmar_refuse_argument(c(
      "The reference sample size N must be a whole number of at least ",
      "three. With incomplete data, explicitly choose a justified ",
      "benchmark when the complete-case count is smaller."), arg = NULL)
  }
  observed <- .pa_reduced_eigen(R, N)
  p <- length(items)
  draws <- matrix(NA_real_, as.integer(B), p, dimnames = list(NULL, items))
  boundary <- logical(B)
  failures <- character(B)
  draw_warnings <- vector("list", B)
  marginals <- if (!is.null(data)) lapply(data, function(x) x[!is.na(x)]) else NULL
  for (b in seq_len(B)) {
    notices <- new.env(parent = emptyenv())
    notices$warnings <- character()
    run <- withCallingHandlers(tryCatch({
      if (reference == "resample") {
        x <- as.data.frame(lapply(marginals, function(x) sample(x, N, replace = TRUE)))
      } else {
        x <- as.data.frame(matrix(stats::rnorm(N * p), N, p))
        if (length(ord)) {
          for (j in seq_len(p)) {
            values <- sort(unique(marginals[[j]]))
            prob <- as.numeric(table(factor(marginals[[j]], levels = values))) /
              length(marginals[[j]])
            cuts <- stats::qnorm(c(0, cumsum(prob)[-length(prob)], 1))
            x[[j]] <- values[findInterval(x[[j]], cuts, all.inside = TRUE)]
          }
        }
      }
      names(x) <- items
      .pa_reduced_eigen(.pa_cor(x, ordered = length(ord) > 0L), N)
    }, error = function(e) list(error = conditionMessage(e))),
    warning = function(w) {
      notices$warnings <- c(notices$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    draw_warnings[[b]] <- notices$warnings
    if (!is.null(run$error)) failures[b] <- run$error else {
      draws[b, ] <- run$values
      boundary[b] <- length(run$boundary) > 0L
    }
  }
  usable <- stats::complete.cases(draws)
  used <- sum(usable)
  # The verdict of the rule (R/dmar_bootstrap.R): the thresholds are
  # withheld when more than half of the reference data sets failed (the
  # count rule is held to the request, checked above).
  verdict <- .dmar_boot_settle(
    B, used, tail = boot_tail, fun = "parallel_analysis",
    quantity = "the reference eigenvalues",
    dropped_because = "gave no defined reference fit",
    what = "reference", position = percentile,
    replications = "reference data sets", signal = FALSE)
  eligible <- verdict$formed
  quantile <- mean <- lower <- upper <- rep(NA_real_, p)
  if (eligible) {
    d <- draws[usable, , drop = FALSE]
    quantile <- apply(d, 2L, stats::quantile, probs = percentile, names = FALSE, type = 7)
    mean <- colMeans(d)
    a <- (1 - conf_level) / 2
    lo <- stats::qbinom(a, used, percentile)
    hi <- stats::qbinom(1 - a, used, percentile) + 1L
    lower <- apply(d, 2L, function(v) if (lo < 1L) -Inf else sort(v)[lo])
    upper <- apply(d, 2L, function(v) if (hi > used) Inf else sort(v)[hi])
  }
  n_warnings <- sum(lengths(draw_warnings) > 0L)
  # One warning, classed by the verdict: dmar_bootstrap_insufficient when
  # the thresholds are withheld, dmar_bootstrap_replications_dropped when
  # reference data sets failed and the thresholds were formed.
  .dmar_boot_signal_text(verdict, if (used < B || n_warnings > 0L) sprintf(
    "%s of %s parallel-analysis reference datasets returned defined eigenvalues.%s%s",
    used, B,
    if (eligible) {
      if (used < B) .dmar_boot_request_note(B, used, verdict$B_required,
                                            "reference data sets")
      else ""
    } else if (identical(verdict$status, "insufficient_share")) {
      " Thresholds and the factor count are withheld: more than half failed."
    } else {
      sprintf(paste0(
        " Thresholds and the factor count are withheld: fewer than the %d ",
        "a reference quantile at this percentile is run from were ",
        "requested (at least 10 beyond the quantile and never fewer than ",
        "100; see ?dmar_bootstrap), or fewer than two were usable."),
        verdict$B_required)
    },
    if (n_warnings) paste0(" ", n_warnings, " reference datasets issued warnings; see counts$warning_messages.") else ""))
  exceeds <- if (eligible) observed$values > quantile else rep(NA, p)
  count <- if (eligible) sum(cumprod(exceeds)) else NA_integer_
  out <- data.frame(term = paste0("rank_", seq_len(p)), rank = seq_len(p),
    estimate = observed$values, reference_quantile = quantile,
    reference_mean = mean, reference_quantile_ci_lower = lower,
    reference_quantile_ci_upper = upper, exceeds_reference = exceeds)
  note <- paste0("Common-factor parallel analysis: fixed one-factor ML ",
    "communality reduction; ", reference, " reference; N = ", N, "; ", used,
    " of ", B, " reference datasets usable. Reference percentile: ",
    format(100 * percentile), "%. Its Monte Carlo limits have nominal ",
    format(100 * conf_level), "% coverage conditional on defined reference fits; discarded failures can select a different null distribution. They are not sampling CIs for ",
    "the observed eigenvalues or the suggested factor count.",
    if (benchmark) " The incomplete-data correlation matrix uses observed values through the stated estimator; reference N is a benchmark, not an effective N or a guaranteed conservative choice." else "",
    " Defined boundary reference fits: ", sum(boundary), ".",
    if (length(observed$boundary)) paste0(" The observed one-factor reduction reaches the lower uniqueness bound for: ", paste(observed$boundary, collapse = ", "), ". Interpret that reduction with this boundary flag.") else "",
    if (identical(verdict$status, "exploratory")) paste0(
      " Exploratory: ", B, " reference data sets are below the ",
      verdict$B_required, " a reported reference quantile needs.") else "")
  attr(out, "suggested_factors") <- as.integer(count)
  out <- .dmar_boot_attrs(out, verdict)
  attr(out, "reference_draws") <- draws
  attr(out, "counts") <- list(attempted = B, usable = used, failed = B - used,
    boundary = sum(boundary), failure_reasons = table(failures[nzchar(failures)]),
    warned = n_warnings, warning_messages = table(unlist(draw_warnings)))
  attr(out, "observed_boundary") <- observed$boundary
  attr(out, "method") <- list(N = N, reference = reference,
    missing = if (is.null(data)) "covariance input" else missing,
    ordered = ord, percentile = percentile, mc_conf_level = conf_level,
    reduction = "one-factor ML", reference_is_benchmark = benchmark)
  attr(out, "analysis_note") <- note
  .as_dmar_tbl(out, subclass = "dmar_parallel_analysis")
}

.pa_cor <- function(data, ordered) {
  if (!ordered) return(stats::cor(data))
  .fa_require_lavaan()
  as.matrix(lavaan::lavCor(data, ordered = names(data), missing = "pairwise",
                          output = "cor"))
}

.pa_reduced_eigen <- function(R, N) {
  if (any(!is.finite(R)) || min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    .dmar_refuse_argument(
      "The analysis correlation matrix is not positive definite.", arg = NULL)
  }
  fit <- stats::factanal(covmat = list(cov = R, n.obs = N), factors = 1L,
                         rotation = "none")
  reduced <- R
  diag(reduced) <- 1 - fit$uniquenesses
  list(values = eigen(reduced, symmetric = TRUE, only.values = TRUE)$values,
       boundary = names(fit$uniquenesses)[fit$uniquenesses <= 0.005 * (1 + 1e-3)])
}

#' @exportS3Method .dmar_tbl_footer dmar_parallel_analysis
.dmar_tbl_footer.dmar_parallel_analysis <- function(x, ...) {
  .dmar_footer_add(NextMethod(), attr(x, "analysis_note", exact = TRUE), "analysis")
}

#' @rdname dmar_tidiers
#' @exportS3Method generics::glance
glance.dmar_parallel_analysis <- function(x, ...) {
  method <- attr(x, "method")
  counts <- attr(x, "counts")
  data.frame(suggested_factors = attr(x, "suggested_factors"),
    reference_N = method$N, reference = method$reference,
    percentile = method$percentile, mc_conf_level = method$mc_conf_level,
    reference_is_benchmark = method$reference_is_benchmark,
    attempted = counts$attempted, usable = counts$usable,
    failed = counts$failed, boundary = counts$boundary,
    stringsAsFactors = FALSE)
}
