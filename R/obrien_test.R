# O'Brien's test for homogeneity of variance.
#' O'Brien's Test for Homogeneity of Variance
#'
#' Tests the null hypothesis that two or more groups have equal population
#' variances using O'Brien's (1981) procedure: each observation is transformed
#' into a quantity whose expected value equals the group's variance, and a
#' one-way analysis of variance is then run on those transformed values. The
#' test is generally regarded as more robust to non-normality than Bartlett's
#' test while retaining good power.
#'
#' @param x Either a numeric vector of observations (in which case \code{group}
#'   must also be supplied), or a one-sided formula of the form
#'   \code{y ~ group}, in which case \code{data} is consulted for the variables.
#' @param group A grouping vector or factor of the same length as \code{x}; used
#'   only when \code{x} is a numeric vector.
#' @param data An optional \code{data.frame} containing the variables named in
#'   the formula.
#' @param na_action Function specifying how missing values are handled
#'   (default \code{\link[stats]{na.omit}}).
#'
#' @return A one-row \code{data.frame} with columns \code{statistic} (the
#'   \emph{F}-value from the ANOVA on the transformed scores), \code{df_1},
#'   \code{df_2}, \code{p_value}, \code{n_groups}, \code{n_total}, and
#'   \code{method}.
#'
#' @details Following O'Brien (1981) and the version given in Abdi (2007),
#' each observation \eqn{Y_{ij}} (the \eqn{j}th observation in group \eqn{i},
#' with size \eqn{n_i} and sample variance \eqn{s_i^2}) is transformed to
#' \deqn{r_{ij} = \frac{(n_i - 1.5)\, n_i\, (Y_{ij} - \bar{Y}_i)^2 - 0.5\, s_i^2\, (n_i - 1)}{(n_i - 1)(n_i - 2)}.}
#' The mean of the \eqn{r_{ij}} within group \eqn{i} equals \eqn{s_i^2}, so a
#' one-way ANOVA on the \eqn{r_{ij}} tests
#' \eqn{H_0\!: \sigma_1^2 = \cdots = \sigma_k^2}.
#' Each group must have at least three observations for the transformation
#' to be defined.
#'
#' @references
#' Abdi, H. (2007). O'Brien's test for homogeneity of variance. In N. J.
#' Salkind (Ed.), \emph{Encyclopedia of measurement and statistics}. Sage.
#'
#' O'Brien, R. G. (1981). A simple test for variance effects in experimental
#' designs. \emph{Psychological Bulletin, 89}(3), 570--574.
#'
#' @examples
#' # Hunter's (1964) "one-is-a-bun" peg-word memory experiment, as discussed
#' # by Abdi (2007). Sixty-four participants were assigned to a control group
#' # (no mnemonic instruction) or an experimental group (peg-word mnemonic).
#' # The score is the number of word pairs (out of 10) recalled. Abdi (2007,
#' # Table 6) reports F = 1.29 (df = 1, 62) for the O'Brien test of equal
#' # variances, p = .260 as computed here; the experimental group's apparent
#' # ceiling effect does not produce statistically detectable variance
#' # heterogeneity.
#' hunter_1964 <- data.frame(
#'   group = factor(
#'     c(rep("Control", 32), rep("Experimental", 32)),
#'     levels = c("Control", "Experimental")
#'   ),
#'   recall = c(
#'     # Control group (n = 32):
#'     5, 5, 5, 5, 5,
#'     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
#'     7, 7, 7, 7, 7, 7, 7, 7, 7,
#'     8, 8, 8,
#'     9, 9,
#'     10, 10,
#'     # Experimental group (n = 32):
#'     6,
#'     7, 7,
#'     8, 8, 8, 8,
#'     9, 9, 9, 9, 9, 9, 9, 9, 9,
#'     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
#'   )
#' )
#' obrien_test(recall ~ group, data = hunter_1964)
#'
#' # Comparison against Bartlett's test on the same data.
#' bartlett.test(recall ~ group, data = hunter_1964)
#'
#' # Vector / grouping-variable interface, on DMAR's depression_bdi data.
#' # The wait list variance is about twice the SSRI variance, but with ten
#' # observations per group the test does not reject equal variances.
#' obrien_test(depression_bdi$bdi_post, depression_bdi$condition)
#'
#' @seealso \code{\link[stats]{bartlett.test}}, \code{\link[stats]{var.test}}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @export
#' @import stats

obrien_test <- function(x, group = NULL, data = NULL, na_action = stats::na.omit) {
  # Resolve formula vs. (vector, group) interface to two parallel vectors.
  if (inherits(x, "formula")) {
    if (length(x) != 3L) stop("'x' must be a two-sided formula like 'y ~ group'.")
    mf <- stats::model.frame(x, data = data, na.action = na_action)
    y_vec <- mf[[1]]
    g_vec <- mf[[2]]
  } else {
    if (is.null(group)) stop("When 'x' is a vector, 'group' must be supplied.")
    if (length(x) != length(group)) stop("'x' and 'group' must have the same length.")
    y_vec <- x
    g_vec <- group
    keep <- !is.na(y_vec) & !is.na(g_vec)
    y_vec <- y_vec[keep]
    g_vec <- g_vec[keep]
  }

  if (!is.numeric(y_vec)) stop("Response variable must be numeric.")
  if (!is.factor(g_vec)) g_vec <- factor(g_vec)
  g_vec <- droplevels(g_vec)

  levs <- levels(g_vec)
  k <- length(levs)
  if (k < 2L) stop("At least two groups are required.")

  group_n <- as.integer(table(g_vec))
  if (any(group_n < 3L)) {
    stop("Every group must contain at least 3 observations for O'Brien's transformation.")
  }

  group_means <- as.numeric(tapply(y_vec, g_vec, mean))
  group_vars  <- as.numeric(tapply(y_vec, g_vec, stats::var))

  # O'Brien (1981) / Abdi (2007) transformation, vectorized over observations.
  idx       <- as.integer(g_vec)
  n_i       <- group_n[idx]
  ybar_i    <- group_means[idx]
  s2_i      <- group_vars[idx]
  numer     <- (n_i - 1.5) * n_i * (y_vec - ybar_i)^2 - 0.5 * s2_i * (n_i - 1)
  denom     <- (n_i - 1) * (n_i - 2)
  r_ij      <- numer / denom

  # One way ANOVA on the transformed scores.
  fit <- stats::aov(r_ij ~ g_vec)
  tbl <- stats::anova(fit)
  F_stat  <- tbl[1, "F value"]
  df_1    <- tbl[1, "Df"]
  df_2    <- tbl["Residuals", "Df"]
  p_value <- tbl[1, "Pr(>F)"]

  out <- data.frame(
    statistic = F_stat,
    df_1      = df_1,
    df_2      = df_2,
    p_value   = p_value,
    n_groups  = k,
    n_total   = length(y_vec),
    method    = "O'Brien's test for homogeneity of variance",
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}
