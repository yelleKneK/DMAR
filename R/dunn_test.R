# Dunn's (1964) nonparametric test of all pairwise differences in mean ranks,
# the rank-based follow-up to a Kruskal-Wallis test.
#' Provides Dunn's Rank-Sum Test of All Pairwise Differences Following a Kruskal--Wallis Test
#'
#' @param x Either (a) a numeric vector of the outcome, in which case
#'   \code{group} must also be supplied, or (b) a fitted
#'   \code{\link[stats]{lm}} or \code{\link[stats]{aov}} object with a
#'   single factor predictor, from which the outcome and grouping factor
#'   are taken (the model itself is not used, since the test is
#'   distribution free).
#' @param group When \code{x} is a vector, a factor (or coercible to
#'   factor) giving the group membership of each observation.
#' @param method The multiplicity adjustment applied to the pairwise
#'   \emph{p}-values, passed to \code{\link[stats]{p.adjust}}. Default
#'   \code{"holm"}. Use \code{"none"} to obtain the unadjusted values.
#'
#' @return A \code{data.frame} with one row per pairwise comparison and
#'   columns \code{contrast}, \code{mean_rank_difference}, \code{se},
#'   \code{z_statistic}, \code{p_value}, and \code{p_adjusted}. The table
#'   prints through the \code{\link{dmar_tbl}} display layer.
#'
#' @details
#' \strong{Which Dunn.} Olive Jean Dunn published two different multiple
#' comparison procedures, and both are called \dQuote{Dunn's} in the
#' literature. This function implements the \emph{nonparametric} one of
#' Dunn (1964), which compares mean ranks. It is not the Bonferroni
#' procedure of Dunn (1961), which Maxwell, Delaney, and Kelley (2027,
#' Chapter 5) call Dunn's procedure and which reaches DMAR through the
#' \code{method} arguments of \code{\link{contrast_adjusted}} and
#' \code{\link[stats]{p.adjust}}. The two are unrelated apart from their
#' author, and the \code{method} argument here can apply the 1961 procedure
#' to the 1964 procedure's \emph{p}-values.
#'
#' \strong{What it does.} The Kruskal--Wallis test asks whether \emph{any} of
#' the groups differ; it does not say which. Dunn's test is its pairwise
#' follow-up. All \eqn{N} observations are ranked together, with tied values
#' receiving their average rank. Writing \eqn{\bar R_g} for the mean rank of
#' group \eqn{g}, each pair is compared with
#' \deqn{z = \frac{\bar R_g - \bar R_h}{\sqrt{\left[\frac{N(N+1)}{12} - \frac{\sum_i (t_i^3 - t_i)}{12(N-1)}\right]\left(\frac{1}{n_g} + \frac{1}{n_h}\right)}},}
#' where \eqn{t_i} is the number of observations tied at the \eqn{i}th
#' distinct value; the second term in the brackets is the tie correction and
#' vanishes when there are no ties. The statistic is referred to the standard
#' normal distribution, and the resulting \emph{p}-values are adjusted for
#' multiplicity by \code{method}.
#'
#' The pooled ranking is what makes this the right follow-up. Dunn's test
#' ranks across all groups at once and uses the variance of the ranks implied
#' by the Kruskal--Wallis null, so it is consistent with the omnibus test that
#' preceded it. Running a separate Mann--Whitney test on each pair instead
#' re-ranks the data within every pair, which answers a different question for
#' every comparison and is not coherent with the omnibus result.
#'
#' \strong{When to use it.} Use Dunn's test when the outcome is ordinal, or
#' when it is continuous but the normality or homogeneity assumptions behind
#' \code{\link{ci_tukey_kramer}} and \code{\link{ci_games_howell}} are
#' untenable and a rank-based analysis is preferred to a transformation; the
#' design is between subjects; and a significant Kruskal--Wallis test leaves
#' the question of which groups differ. It is the rank analogue of Tukey's
#' HSD, in the sense of covering all \eqn{a(a-1)/2} pairs.
#'
#' \strong{What it does not tell you.} The test compares mean \emph{ranks},
#' not medians. A significant pair means one group's observations tend to be
#' larger, that is, stochastic dominance; it does not by itself license a
#' statement about medians unless the group distributions have the same shape.
#' It also returns no confidence interval on any quantity in the original
#' units, which is a real cost: where a parametric procedure is defensible,
#' \code{\link{ci_games_howell}} or \code{\link{ci_tukey_kramer}} reports
#' intervals on the mean difference, and Maxwell, Delaney, and Kelley (2027)
#' emphasize interval estimation over test decisions throughout. For a
#' distribution free effect size with a confidence interval, see
#' \code{\link{cliff_delta}}.
#'
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums.
#'   \emph{Technometrics, 6}(3), 241--252.
#'   \doi{10.1080/00401706.1964.10490181}
#'
#' Dunn, O. J. (1961). Multiple comparisons among means. \emph{Journal of
#'   the American Statistical Association, 56}(293), 52--64.
#'   \doi{10.2307/2282330} (The Bonferroni procedure; not what this
#'   function computes.)
#'
#' Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion
#'   variance analysis. \emph{Journal of the American Statistical
#'   Association, 47}(260), 583--621. \doi{10.2307/2280779}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @examples
#' # The omnibus question first: do the three treatment arms of the
#' # drinks_trial data differ at all? The raw drinks_per_week outcome is
#' # markedly right-skewed, which is no obstacle here: ranks are invariant
#' # to a monotone transformation, so the raw and log scales give
#' # identical results.
#' kruskal.test(drinks_per_week ~ treatment, data = drinks_trial)
#'
#' # Then which pairs, on the same pooled ranking the omnibus test used.
#' dunn_test(drinks_trial$drinks_per_week, drinks_trial$treatment)
#'
#' # Without a multiplicity adjustment (rarely what you want).
#' dunn_test(drinks_trial$drinks_per_week, drinks_trial$treatment,
#'           method = "none")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link[stats]{kruskal.test}} for the omnibus test this
#'   follows, \code{\link{ci_games_howell}} and \code{\link{ci_tukey_kramer}}
#'   for the parametric all-pairs procedures, \code{\link{cliff_delta}} for a
#'   distribution free effect size with a confidence interval, and
#'   \code{\link[stats]{p.adjust}} for the \code{method} choices.
#'
#' @keywords htest nonparametric
#'
#' @export
#' @import stats
dunn_test <- function(x, group = NULL, method = "holm") {
  if (!method %in% stats::p.adjust.methods)
    stop("'method' must be one of the stats::p.adjust.methods: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".")

  if (inherits(x, c("lm", "aov"))) {
    mf <- stats::model.frame(x)
    if (ncol(mf) != 2L)
      stop("dunn_test() requires a fitted lm/aov with one factor predictor.")
    y <- mf[[1L]]; g <- mf[[2L]]
    if (!is.factor(g)) g <- factor(g)
  } else if (is.numeric(x) && !is.null(group)) {
    if (length(x) != length(group))
      stop("'x' and 'group' must be the same length.")
    y <- x; g <- factor(group)
  } else {
    stop("Supply either a numeric vector 'x' with a 'group' factor, or a fitted lm/aov.")
  }

  ok <- !is.na(y) & !is.na(g)
  y <- y[ok]; g <- droplevels(g[ok])

  lev <- levels(g)
  a   <- length(lev)
  if (a < 2) stop("At least two groups are required.")
  N <- length(y)

  # One ranking across all groups, as the Kruskal-Wallis test uses; ties take
  # the average rank (Dunn, 1964).
  r         <- rank(y)
  mean_rank <- tapply(r, g, mean)
  ns        <- tapply(r, g, length)

  # Tie correction: sum over distinct values of (t^3 - t), which is 0 with no
  # ties, giving the untied variance N(N+1)/12.
  tie_counts <- as.numeric(table(y))
  tie_term   <- sum(tie_counts^3 - tie_counts)
  sigma2     <- N * (N + 1) / 12 - tie_term / (12 * (N - 1))

  # When every observation is tied at a single value the tie correction cancels
  # the untied variance exactly, leaving sigma2 = 0. The standard error would be
  # zero and every z statistic and p-value an undefined 0/0, so refuse the
  # degenerate input rather than return a table of NaN.
  if (!is.finite(sigma2) || sigma2 <= 0)
    stop("The tie-corrected variance of the ranks is zero because all observations are tied at a single value; Dunn's test is undefined. Provide data with at least two distinct values.")

  rows <- list(); idx <- 0L
  for (i in seq_len(a - 1L)) {
    for (j in (i + 1L):a) {
      idx  <- idx + 1L
      se   <- sqrt(sigma2 * (1 / ns[i] + 1 / ns[j]))
      d_ij <- as.numeric(mean_rank[j] - mean_rank[i])
      z    <- d_ij / as.numeric(se)
      rows[[idx]] <- data.frame(
        contrast             = paste(lev[j], "-", lev[i]),
        mean_rank_difference = d_ij,
        se                   = as.numeric(se),
        z_statistic          = z,
        p_value              = 2 * stats::pnorm(-abs(z)),
        stringsAsFactors     = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$p_adjusted <- stats::p.adjust(out$p_value, method = method)
  rownames(out) <- NULL
  attr(out, "adjustment_method") <- method
  .as_dmar_tbl(out)
}
