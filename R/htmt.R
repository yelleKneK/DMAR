#' Heterotrait-Monotrait Ratio of Correlations (HTMT)
#'
#' Computes the HTMT discriminant-validity index of Henseler, Ringle, and
#' Sarstedt (2015) for every pair of constructs: the average correlation
#' between items of \emph{different} constructs, divided by the geometric
#' mean of the average correlations among items \emph{within} each
#' construct. Two constructs whose HTMT approaches 1 are empirically
#' indistinguishable however cleanly the model draws them; the customary
#' red flags are 0.85 (strict) or 0.90 (liberal). An optional bootstrap
#' gives a one-sided upper confidence bound, the quantity actually
#' compared against the cutoff in the validity literature.
#'
#' @param data A \code{data.frame} of item responses.
#' @param blocks Named list of character vectors: each element names a
#'   construct and gives its item columns (two or more per construct, two
#'   or more constructs).
#' @param B Number of bootstrap resamples for the upper confidence bound;
#'   \code{0} (default) skips the bootstrap and reports the point
#'   estimates only.
#' @param conf_level Confidence level for the one-sided upper bound.
#'   Defaults to 0.95.
#' @param seed Optional integer seed for the bootstrap, used locally (the
#'   caller's random number generator state is restored on exit).
#'
#' @details
#' All correlations are Pearson, computed on pairwise-complete
#' observations. The statistic uses absolute average heterotrait
#' correlations enter as absolute values, the convention of later
#' implementations (the 2015 proposal used the plain correlations,
#' which can cancel when signs mix; Roemer, Schuberth, & Henseler,
#' 2021, recommend the absolute form for exactly that reason), and
#' values near or above 1 indicating that the two item sets correlate
#' across constructs about as strongly as within them. HTMT is a
#' correlation-based screen, deliberately model-free; the confirmatory
#' companion is the latent correlation between the two factors (see
#' \code{\link{correction_for_attenuation}} and its factor-model
#' discussion).
#'
#' The bootstrap, when requested (\code{B > 0}), resamples the rows of
#' \code{data} with replacement \code{B} times and recomputes every
#' pairwise HTMT on each resample; the reported \code{upper_limit} is
#' the \code{conf_level} empirical quantile of each pair's bootstrap
#' distribution, a one-sided upper percentile bound (Efron &
#' Tibshirani, 1993). That bound is the only interval offered, matching
#' how the validity literature uses HTMT (the question is whether the
#' ratio credibly exceeds the cutoff); no two-sided or bias-corrected
#' and accelerated (BCa) variant is provided. A resample in which some
#' block's average within-construct correlation is not positive leaves
#' HTMT undefined there; such resamples are dropped, a single warning
#' reports how many, and the bound is computed from the resamples that
#' remained (the call stops only when fewer than 100 remain). Bootstrap
#' results vary from run to run; supply \code{seed} for
#' reproducibility.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row
#'   per construct pair: \code{construct_1}, \code{construct_2},
#'   \code{htmt}, and, when \code{B > 0}, \code{upper_limit} (the
#'   one-sided \code{conf_level} bootstrap percentile bound).
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for
#'   assessing discriminant validity in variance-based structural equation
#'   modeling. \emph{Journal of the Academy of Marketing Science, 43}(1),
#'   115--135. \doi{10.1007/s11747-014-0403-8}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{average_variance_extracted}} for the convergent
#'   side of the validity ledger; \code{\link{reliability_omega}} for the
#'   composite reliability of each block;
#'   \code{\link{correction_for_attenuation}} for the latent correlation
#'   route.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Two clean constructs and their six items.
#' set.seed(113)
#' n <- 300
#' f1 <- rnorm(n); f2 <- 0.3 * f1 + sqrt(1 - 0.09) * rnorm(n)
#' d <- data.frame(
#'   a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
#'   a3 = .6 * f1 + rnorm(n, 0, .8),
#'   b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
#'   b3 = .6 * f2 + rnorm(n, 0, .8))
#' h <- htmt(d, blocks = list(A = c("a1", "a2", "a3"),
#'                            B = c("b1", "b2", "b3")))
#' h
#'
#' # The broom verbs: one row per construct pair.
#' generics::tidy(h)
#' generics::glance(h)
#'
#' # The upper confidence bound, which is the quantity the validity
#' # literature compares against 0.85 or 0.90, comes from a bootstrap
#' # that recomputes every pairwise ratio on each of B resamples. It is
#' # shown rather than run; the call is
#' #   htmt(d, blocks = list(A = c("a1", "a2", "a3"),
#' #                         B = c("b1", "b2", "b3")),
#' #        B = 10000, seed = 113)
#' # and a claim about discriminant validity deserves that bound rather
#' # than the point estimate alone.
#'
#' @export
#' @importFrom stats cor quantile
htmt <- function(data, blocks, B = 0, conf_level = 0.95, seed = NULL) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  if (!is.list(blocks) || length(blocks) < 2L ||
      is.null(names(blocks)) || any(names(blocks) == "")) {
    stop("'blocks' must be a named list of two or more constructs.",
         call. = FALSE)
  }
  for (nm in names(blocks)) {
    bl <- blocks[[nm]]
    if (!is.character(bl) || length(bl) < 2L || !all(bl %in% names(data))) {
      stop(sprintf("Block '%s' must name two or more columns of 'data'.",
                   nm), call. = FALSE)
    }
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 0 ||
      B != round(B)) {
    stop("'B' must be a single non-negative integer.", call. = FALSE)
  }
  if (B > 0 && B < 100) {
    stop("'B' must be at least 100 when bootstrapping.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  htmt_from <- function(d) {
    R <- stats::cor(d[, unlist(blocks), drop = FALSE],
                    use = "pairwise.complete.obs")
    pairs <- utils::combn(names(blocks), 2, simplify = FALSE)
    vapply(pairs, function(pr) {
      i <- blocks[[pr[1]]]; j <- blocks[[pr[2]]]
      hetero <- mean(abs(R[i, j]))
      mono_i <- mean(R[i, i][lower.tri(R[i, i])])
      mono_j <- mean(R[j, j][lower.tri(R[j, j])])
      if (mono_i <= 0 || mono_j <= 0) {
        stop("A block's average within-construct correlation is not ",
             "positive; HTMT is undefined for such a block.",
             call. = FALSE)
      }
      hetero / sqrt(mono_i * mono_j)
    }, numeric(1))
  }

  est <- htmt_from(data)
  pairs <- utils::combn(names(blocks), 2)

  out <- data.frame(
    construct_1 = pairs[1, ],
    construct_2 = pairs[2, ],
    htmt        = unname(est),
    stringsAsFactors = FALSE
  )

  if (B > 0) {
    if (!is.null(seed)) {
      has_old <- exists(".Random.seed", envir = globalenv())
      old <- if (has_old) get(".Random.seed", envir = globalenv()) else NULL
      on.exit({
        if (has_old) assign(".Random.seed", old, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv()))
          rm(".Random.seed", envir = globalenv())
      }, add = TRUE)
      set.seed(seed)
    }
    N <- nrow(data)
    # A resample can draw a case mix whose average within-construct
    # correlation is not positive, which leaves HTMT undefined there;
    # drop that resample rather than aborting the whole call (the
    # observed-data estimate above has already passed the same guard).
    boots <- vapply(seq_len(B), function(i) {
      tryCatch(
        htmt_from(data[sample.int(N, N, replace = TRUE), , drop = FALSE]),
        error = function(e) rep(NA_real_, length(est)))
    }, numeric(length(est)))
    boots <- matrix(boots, nrow = length(est))
    ok <- colSums(!is.finite(boots)) == 0L
    n_bad <- sum(!ok)
    if (n_bad > 0L) {
      if (sum(ok) < 100L) {
        stop("Only ", sum(ok), " of ", B, " bootstrap resamples produced ",
             "a defined HTMT for every pair; the upper bound would not ",
             "be trustworthy. Check the constructs whose within-block ",
             "correlations sit near zero.", call. = FALSE)
      }
      warning(n_bad, " of ", B, " bootstrap resamples left HTMT ",
              "undefined for some pair (average within-construct ",
              "correlation not positive) and were dropped; the upper ",
              "bound is computed from the ", sum(ok), " that remained.",
              call. = FALSE)
      boots <- boots[, ok, drop = FALSE]
    }
    out$upper_limit <- apply(boots, 1, stats::quantile,
                             probs = conf_level, names = FALSE)
  }
  res <- .as_dmar_tbl(out)
  if (B > 0) attr(res, "conf_level") <- conf_level
  res
}
