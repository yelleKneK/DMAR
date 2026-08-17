#' Combine Independent P-Values Across Studies
#'
#' Combines the one-tailed \emph{p}-values from \eqn{k} independent tests of
#' a common directional hypothesis into a single test, by any of the four
#' classical methods that Raudenbush (1984) applied to the teacher expectancy
#' experiments: Fisher's (1938) \dQuote{adding logs} chi square, Edgington's
#' (1972) \dQuote{adding \emph{p}s}, the Mosteller and Bush (1954)
#' \dQuote{adding \emph{Z}s} (Stouffer) method, and a weighted adding-Zs
#' variant (weights are typically the studies' degrees of freedom). Combined
#' significance tests answer the narrow question \dQuote{is there an effect
#' in at least some studies?}; they do not estimate its size. Pair them with
#' \code{\link{meta_smd}} or \code{\link{meta_es}} for estimation, which is
#' almost always the more informative summary.
#'
#' @param p Numeric vector of one-tailed \emph{p}-values, each in (0, 1),
#'   oriented so that small values support the common directional
#'   hypothesis.
#' @param method Character vector naming the methods to compute: any of
#'   \code{"fisher"}, \code{"edgington"}, \code{"stouffer"},
#'   \code{"stouffer_weighted"}; the default computes all four (the
#'   \code{"stouffer_weighted"} row appears only when \code{weights} is
#'   supplied).
#' @param weights Optional non-negative weights for
#'   \code{"stouffer_weighted"}, one per study; degrees of freedom are the
#'   conventional choice (Mosteller & Bush, 1954).
#'
#' @details
#' Fisher's statistic is \eqn{-2 \sum \log p_i}, distributed chi square with
#' \eqn{2k} degrees of freedom under the joint null. Edgington's statistic
#' is the plain sum \eqn{\sum p_i}, referred to a normal approximation with
#' mean \eqn{k/2} and variance \eqn{k/12} (accurate for \eqn{k \ge 10}; for
#' smaller \eqn{k} it is conservative in the tails). The Stouffer statistic
#' is \eqn{\sum z_i / \sqrt{k}} with \eqn{z_i = \Phi^{-1}(1 - p_i)}, and the
#' weighted variant is \eqn{\sum w_i z_i / \sqrt{\sum w_i^2}}. All four are
#' reported with one-tailed combined \emph{p}-values, matching the
#' directional inputs.
#'
#' Methods can disagree, and the disagreement is informative: Rosenthal
#' (1978) notes there is no uniformly best test. In the published analysis,
#' Raudenbush (1984) found three of the four rejecting the null at the .05
#' level while the df-weighted variant did not, an early warning that large
#' studies were finding smaller effects. Computed from the study-level
#' \emph{p}-values as tabled, the example below shows two of the four
#' rejecting: Fisher's (\eqn{p = .004}) and Stouffer's (\eqn{p = .014})
#' tests reject, Edgington's sits just above the level (\eqn{p = .051}; the
#' tabled values sum to 7.00 where the paper's Table 2, p. 90, prints a
#' sum of 6.84 with \eqn{p = .04}), and the df-weighted variant is not
#' close (\eqn{p = .192}).
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with, per
#'   requested method, its statistic row(s) and a one-tailed
#'   \code{<method>_p} row, plus a final \code{k} row. The \emph{p} rows
#'   print to fixed decimals via the \code{p_terms} attribute.
#'
#' @references
#' Edgington, E. S. (1972). An additive method for combining probability
#'   values from independent experiments. \emph{The Journal of Psychology,
#'   80}(2), 351--363.
#'
#' Fisher, R. A. (1938). \emph{Statistical methods for research workers}
#'   (7th ed.). Oliver & Boyd.
#'
#' Mosteller, F., & Bush, R. R. (1954). Selected quantitative techniques. In
#'   G. Lindzey (Ed.), \emph{Handbook of social psychology} (Vol. 1).
#'   Addison-Wesley.
#'
#' Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on
#'   pupil IQ as a function of the credibility of expectancy induction: A
#'   synthesis of findings from 18 experiments. \emph{Journal of Educational
#'   Psychology, 76}(1), 85--97.
#'
#' Rosenthal, R. (1978). Combining results of independent studies.
#'   \emph{Psychological Bulletin, 85}(1), 185--193.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{meta_smd}} and \code{\link{meta_es}} for estimating
#'   the pooled effect rather than only testing it;
#'   \code{\link{meta_contrast}} for differences among study effects;
#'   \code{\link{teacher_expectancy}} for the data behind the examples.
#'
#' @family meta-analysis
#'
#' @keywords htest
#'
#' @examples
#' # Raudenbush (1984), Table 2: the four combined tests over the 18
#' # teacher expectancy studies (Pellegrini & Hicks at its study-level
#' # values), weighting the Z method by degrees of freedom.
#' data(teacher_expectancy)
#' study <- teacher_expectancy[-c(4, 5), ]
#' p18  <- append(study$p_one_tailed, .010, after = 3)
#' df18 <- append(study$n_experimental + study$n_control - 2, 42, after = 3)
#' combine_p(p18, weights = df18)
#' # Fisher chi square 62.17 on 36 df; Edgington sum near 7; Stouffer
#' # z near 2.2; and the df-weighted z under 1: the large studies disagree.
#'
#' @export
#' @importFrom stats pchisq pnorm qnorm
combine_p <- function(p,
                      method = c("fisher", "edgington", "stouffer",
                                 "stouffer_weighted"),
                      weights = NULL) {
  if (!is.numeric(p) || length(p) < 2L || anyNA(p) ||
      any(p <= 0) || any(p >= 1)) {
    stop("'p' must be two or more p-values, each strictly inside (0, 1).",
         call. = FALSE)
  }
  method <- match.arg(method, several.ok = TRUE)
  if (is.null(weights)) {
    method <- setdiff(method, "stouffer_weighted")
  } else {
    if (!is.numeric(weights) || length(weights) != length(p) ||
        anyNA(weights) || any(weights < 0) || all(weights == 0)) {
      stop("'weights' must be non-negative, with one weight per p-value ",
           "and at least one positive.", call. = FALSE)
    }
    if (!("stouffer_weighted" %in% method)) {
      method <- c(method, "stouffer_weighted")
    }
  }
  if (length(method) == 0L) {
    stop("No method left to compute: 'stouffer_weighted' needs 'weights'.",
         call. = FALSE)
  }

  k <- length(p)
  term <- character(0); value <- numeric(0); p_terms <- character(0)
  add <- function(t, v) { term <<- c(term, t); value <<- c(value, v) }

  if ("fisher" %in% method) {
    X2 <- -2 * sum(log(p))
    add(c("fisher_chi_square", "fisher_df", "fisher_p"),
        c(X2, 2 * k, pchisq(X2, df = 2 * k, lower.tail = FALSE)))
    p_terms <- c(p_terms, "fisher_p")
  }
  if ("edgington" %in% method) {
    S <- sum(p)
    z <- (S - k / 2) / sqrt(k / 12)
    add(c("edgington_sum_p", "edgington_p"), c(S, pnorm(z)))
    p_terms <- c(p_terms, "edgington_p")
  }
  if ("stouffer" %in% method) {
    z <- sum(qnorm(1 - p)) / sqrt(k)
    add(c("stouffer_z", "stouffer_p"), c(z, pnorm(z, lower.tail = FALSE)))
    p_terms <- c(p_terms, "stouffer_p")
  }
  if ("stouffer_weighted" %in% method) {
    zw <- sum(weights * qnorm(1 - p)) / sqrt(sum(weights^2))
    add(c("stouffer_weighted_z", "stouffer_weighted_p"),
        c(zw, pnorm(zw, lower.tail = FALSE)))
    p_terms <- c(p_terms, "stouffer_weighted_p")
  }
  add("k", k)

  .as_dmar_tbl(data.frame(term = term, value = value,
                          stringsAsFactors = FALSE),
               p_terms = p_terms)
}
