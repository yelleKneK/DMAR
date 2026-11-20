#' Bryant--Paulson Simultaneous Confidence Intervals for Contrasts of Adjusted Means in ANCOVA
#'
#' @description
#' Constructs Tukey--Kramer-type \emph{simultaneous} confidence intervals on
#' one or more contrasts of covariate-adjusted means in the analysis of
#' covariance (ANCOVA) when the covariate(s) are \emph{random}, using the
#' Bryant--Paulson generalized studentized range (see
#' \code{\link{bryant_paulson}}). Unlike the per-comparison interval of
#' \code{\link{ci_c_ancova}}, these intervals control the
#' \emph{familywise} error rate over the whole family of comparisons and,
#' through the Bryant--Paulson critical value, correctly account for the
#' extra sampling uncertainty that comes from estimating the covariate
#' adjustment from random covariates. Naively applying Tukey's method to
#' adjusted means ignores that uncertainty and produces intervals that are
#' too narrow (below-nominal coverage). The Bryant--Paulson critical value
#' exceeds the studentized range critical value for the same
#' \eqn{\alpha}, number of means, and error degrees of freedom, so these
#' intervals are wider than Tukey's (Maxwell, Delaney, & Kelley, 2027,
#' Chapter 9). Over the range of Bryant and Paulson's (1976) Table 1, the
#' excess is under 5\% for every number of means and both tabled
#' \eqn{\alpha} levels (.05 and .01) when the error degrees of freedom
#' \eqn{\nu} are at least 14 with one covariate, 30 with two, or 40 with
#' three. It grows as \eqn{\nu} falls and as the numbers of covariates
#' and means rise; at \eqn{\nu = 2} and \eqn{\alpha = .05} it reaches
#' 34\% with one covariate and 84\% with three (92\% at
#' \eqn{\alpha = .01}).
#'
#' @param adj_means A numeric vector of the covariate-adjusted group means
#'   (one per group). The number of groups \eqn{k} is inferred from its
#'   length.
#' @param s_ancova The standard deviation of the errors from the ANCOVA
#'   model, i.e., the square root of the ANCOVA error mean square (use the
#'   standard deviation, \emph{not} the variance from the source table).
#' @param c_weights Optional contrast weights. May be (i) a numeric vector of
#'   length \eqn{k} giving a single contrast, or (ii) a matrix/data.frame with
#'   \eqn{k} columns, each row a contrast. If \code{NULL} (the default), all
#'   \eqn{k(k-1)/2} pairwise comparisons are returned. For each contrast the
#'   weights should sum to zero.
#' @param n_per_group Either a single number giving the common per-group sample size or
#'   a numeric vector of per-group sample sizes. The Bryant--Paulson
#'   distribution is exact for balanced designs; for unequal \eqn{n} a
#'   Tukey--Kramer harmonic adjustment is used (see Details).
#' @param n_covariates The number of random covariates, \eqn{p}.
#'   Default \code{1}.
#' @param df Optional error degrees of freedom \eqn{\nu}. If \code{NULL}
#'   (default) it is computed for a one-way ANCOVA as
#'   \eqn{\nu = \sum n - k - p}. Supply it directly for other designs (e.g.,
#'   a randomized-block ANCOVA, where \eqn{\nu} differs).
#' @param conf_level The simultaneous (familywise) confidence level.
#'   Default \code{0.95}.
#' @param contrast_type One of \code{"pairwise"} (default) or
#'   \code{"allowance"}. \code{"pairwise"} uses the Tukey--Kramer quadratic
#'   standard error, exact for pairwise comparisons. \code{"allowance"} uses
#'   Tukey's allowance, which yields intervals that hold simultaneously over
#'   \emph{all} contrasts, including complex ones (this is the form in
#'   Eq. (2.4) of Bryant and Bruvold, 1980). The width factor each choice
#'   applies is given in Details. The two coincide for pairwise comparisons.
#'
#' @details
#' \strong{The interval.} For a contrast \eqn{\psi = \sum_i c_i \theta_i} of
#' adjusted means, the simultaneous interval is
#' \deqn{\hat\psi \;\pm\; q_{\alpha;\,p,k,\nu}\; \hat\sigma_{y\mid x}\; w(c),}
#' where \eqn{q_{\alpha;\,p,k,\nu}} is the upper-\eqn{\alpha} Bryant--Paulson
#' critical value (\code{\link{qbryant_paulson}}) and the width factor is
#' \eqn{w(c) = \tfrac{1}{\sqrt2}\sqrt{\sum_i c_i^2/n_i}} for
#' \code{contrast_type = "pairwise"} or
#' \eqn{w(c) = \tfrac12 \sum_i |c_i| \sqrt{1/n}} for
#' \code{contrast_type = "allowance"} (balanced \eqn{n}). For a pairwise
#' difference with common \eqn{n} both reduce to
#' \eqn{q_{\alpha;\,p,k,\nu}\,\hat\sigma_{y\mid x}\sqrt{1/n}}, reproducing the
#' critical difference of Bryant and Bruvold (1980).
#'
#' \strong{No per-comparison covariate term.} By design the standard error
#' here does \emph{not} include the
#' \eqn{(\bar X_i - \bar X_j)^2 / SS_{\mathrm{within}(x)}} term that appears
#' in a single-comparison ANCOVA interval (\code{\link{ci_c_ancova}}).
#' In the Bryant--Paulson framework the random-covariate uncertainty is
#' carried by the (larger) critical value, which holds \emph{on average} over
#' the covariate distribution; adding the per-pair term as well would
#' double-count it.
#'
#' \strong{Unequal sample sizes.} For unbalanced designs the
#' \code{"pairwise"} standard error uses \eqn{\sqrt{c_i^2/n_i}} directly
#' (the Tukey--Kramer generalization); coverage is then approximate but
#' typically very close to nominal and slightly conservative.
#'
#' @return
#' A \code{data.frame} (class \code{dmar_tbl}) with one row per
#' contrast and columns \code{contrast} (a label), \code{estimate} (the
#' contrast of adjusted means \eqn{\hat\psi}), \code{ci_lower}, and
#' \code{ci_upper}. The Bryant--Paulson critical value used is stored in
#' the \code{"critical_value"} attribute and the confidence level in the
#' \code{"conf_level"} attribute (printed beneath the table).
#'
#' @references
#' Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method of
#'   multiple comparisons to experimental designs with random concomitant
#'   variables. \emph{Biometrika, 63}, 631--638.
#'
#' Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures in
#'   the analysis of covariance. \emph{Journal of the American Statistical
#'   Association, 75}(372), 874--880. \doi{10.2307/2287175}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{bryant_paulson}} for the underlying distribution;
#' \code{\link{ci_c_ancova}} for the per-comparison interval;
#' \code{\link{ancova}} for an ANCOVA fit.
#'
#' @examples
#' # Bryant and Bruvold (1980) worked example: 6 panels, 1 covariate, nu = 14,
#' # ANCOVA error MS = 0.01326. Here the design is a randomized block with
#' # s = 4 blocks, so the per-group size for the adjusted-mean SE is 4 and the
#' # error df (14) must be supplied directly. The multiplier for these
#' # intervals is a Bryant-Paulson quantile, which has no closed form; it is
#' # obtained by inverting the distribution function with a root search.
#' adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
#' bp <- ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326),
#'                      n_per_group = 4, n_covariates = 1, df = 14)
#' bp
#' # The multiplier is 4.83 and every pairwise critical difference is 0.278,
#' # matching the paper; the 15 intervals hold jointly at the 95 percent level.
#' # The multiplier is kept on the result, so the critical difference can be
#' # rebuilt by hand as q * s_ancova * sqrt(1/n):
#' attr(bp, "critical_value")
#' attr(bp, "critical_value") * sqrt(0.01326) * sqrt(1 / 4)
#'
#' # A complex contrast (say panels 1 and 2 against panels 3 through 6) is
#' # requested by passing its weights to c_weights, together with
#' # contrast_type = "allowance", the all-contrasts form of Eq. (2.4) of
#' # Bryant and Bruvold. That contrast of adjusted means is -0.921, with
#' # simultaneous limits of -1.199 and -0.643.
#' ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326),
#'                c_weights = c(0.5, 0.5, -0.25, -0.25, -0.25, -0.25),
#'                n_per_group = 4, n_covariates = 1, df = 14,
#'                contrast_type = "allowance")
#'
#' @keywords design htest
#'
#' @family confidence intervals for effect sizes
#'
#' @export
ci_c_ancova_bp <- function(adj_means, s_ancova, c_weights = NULL, n_per_group,
                           n_covariates = 1, df = NULL, conf_level = 0.95,
                           contrast_type = c("pairwise", "allowance")) {
  .dmar_check_conf_level(conf_level)
  .dmar_check_numbers("ci_c_ancova_bp", "df")
  contrast_type <- match.arg(contrast_type)
  if (missing(adj_means) || !is.numeric(adj_means))
    .dmar_refuse_argument(
      "'adj_means' must be a numeric vector of adjusted group means.",
      arg = "adj_means")
  # The ANCOVA error standard deviation scales every half-width, so it is
  # required and must be a single positive, finite number
  # (.check_error_sd(), shared with ci_c_ancova()); an NA stopped with
  # "missing value where TRUE/FALSE needed".
  .check_error_sd(if (missing(s_ancova)) NULL else s_ancova, "s_ancova",
                  "ANCOVA")
  k <- length(adj_means)
  if (k < 2L) .dmar_refuse_argument(
    "At least two groups (adjusted means) are required.", arg = NULL)
  if (conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")
  if (length(n_per_group) == 1L) n_per_group <- rep(n_per_group, k)
  if (length(n_per_group) != k)
    .dmar_refuse_argument(
      "'n_per_group' must have length 1 or length(adj_means).",
      arg = "n_per_group")

  # Contrast matrix C (rows = contrasts) and labels.
  if (is.null(c_weights)) {
    pairs <- utils::combn(k, 2)
    C <- matrix(0, nrow = ncol(pairs), ncol = k)
    labels <- character(ncol(pairs))
    for (j in seq_len(ncol(pairs))) {
      C[j, pairs[1, j]] <-  1
      C[j, pairs[2, j]] <- -1
      labels[j] <- sprintf("group_%d - group_%d", pairs[1, j], pairs[2, j])
    }
  } else {
    C <- if (is.matrix(c_weights) || is.data.frame(c_weights))
      as.matrix(c_weights) else matrix(c_weights, nrow = 1)
    if (ncol(C) != k)
      .dmar_refuse_argument(
        "'c_weights' must have length / number of columns equal to the number of groups.",
        arg = "c_weights")
    # A contrast (row) whose weights are all zero defines no contrast; it
    # gave the interval [0, 0].
    .dmar_check_contrast_nonzero(c_weights, margin = 1L)
    if (any(abs(rowSums(C)) > 1e-8))
      .dmar_warn(
        "One or more contrasts do not sum to zero.", "result_qualified")
    labels <- rownames(C)
    if (is.null(labels)) labels <- paste0("contrast_", seq_len(nrow(C)))
  }

  if (is.null(df)) df <- sum(n_per_group) - k - n_covariates
  if (df <= 0) .dmar_refuse_argument(
    "Error degrees of freedom (df) must be positive; supply 'df' explicitly.",
    arg = "df")

  # The simultaneous interval is given by Bryant and Paulson (1976, Equation 8):
  #     sum_i c_i theta_i_hat  +/-  q_{alpha;p,k,nu} * sigma_hat_{y|x} * w(c),
  # where q is the upper-alpha Bryant-Paulson critical value (qbryant_paulson)
  # and the design constant {(K1 - K2) v/v}^{1/2} (sum |c_i|) of their
  # Equation (8) specializes, for a one-way design with per-group n, to the
  # width factors w(c) below. Conditions (i) to (iii) of Bryant and Paulson
  # (1976), the analysis-of-covariance analogs of the balance conditions
  # for Tukey's method, are assumed.
  qcrit <- qbryant_paulson(conf_level, n_covariates = n_covariates,
                           n_groups = k, df = df)

  est <- as.numeric(C %*% adj_means)
  # Width factor per contrast.
  if (contrast_type == "pairwise") {
    # Tukey-Kramer quadratic standard error: w = (1/sqrt2) sqrt(sum c_i^2 / n_i),
    # exact for pairwise differences. For a balanced pairwise contrast this
    # equals sqrt(1/n), so the half-width is q * sigma_hat_{y|x} * sqrt(1/n),
    # the critical difference of Bryant and Paulson (1976, Equation 8) and
    # Bryant and Bruvold (1980).
    w <- sqrt(rowSums(sweep(C^2, 2, n_per_group, "/"))) / sqrt(2)
  } else {                                   # "allowance": all-contrasts form
    # Tukey allowance (1/2) sum |c_i| from Bryant and Paulson (1976,
    # Equation 8); holds simultaneously over all contrasts, including complex
    # ones.
    if (length(unique(n_per_group)) > 1L)
      .dmar_warn(
        "The 'allowance' form assumes a balanced design; using the harmonic mean of 'n_per_group'.",
        "result_qualified")
    n_bar <- length(n_per_group) / sum(1 / n_per_group)
    w <- 0.5 * rowSums(abs(C)) * sqrt(1 / n_bar)
  }
  half <- qcrit * s_ancova * w

  out <- data.frame(
    contrast    = labels,
    estimate    = est,
    ci_lower = est - half,
    ci_upper = est + half,
    stringsAsFactors = FALSE, row.names = NULL
  )
  attr(out, "critical_value") <- qcrit
  if (exists(".as_dmar_tbl", mode = "function")) {
    .as_dmar_tbl(out, conf_level = conf_level)
  } else {
    attr(out, "conf_level") <- conf_level
    out
  }
}
