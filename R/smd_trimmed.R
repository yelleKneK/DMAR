# Robust standardized mean difference with trimmed means and Winsorized SD.
#' Robust Standardized Mean Difference (Algina-Keselman-Penfield)
#'
#' Computes the Algina, Keselman, and Penfield (2005) robust
#' standardized mean difference, which replaces the sample means and
#' pooled SD in Cohen's \emph{d} with their trimmed-mean and
#' Winsorized-SD counterparts:
#' \deqn{d_{R} \;=\; 0.642 \cdot \frac{\bar X_{t,\, 1} - \bar X_{t,\, 2}}
#'                                     {s_{W,\, p}},}
#' where \eqn{\bar X_{t,\, j}} is the trimmed mean of group \eqn{j},
#' \eqn{s_{W,\, p}} is the pooled Winsorized standard deviation, and
#' \eqn{0.642} is the Algina-Keselman-Penfield (2005) constant chosen
#' so that \eqn{d_R} equals Cohen's \eqn{\delta} when the data are
#' normal. Returns the point estimate, a noncentral \emph{t}
#' confidence interval, and the trimmed / Winsorized summary
#' statistics.
#'
#' @param x,y Numeric vectors of observations from the two groups.
#' @param trim Proportion to trim from each tail (and Winsorize from
#'   each tail). Must be in \eqn{[0, 0.5)}. Default \code{0.20}
#'   (Wilcox's 2017 recommended setting).
#' @param conf_level Confidence level for the CI. Default
#'   \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the robust \emph{d}
#'   estimate, the lower/upper CI bounds, the per-group trimmed
#'   means, the per-group Winsorized SDs, the pooled Winsorized SD,
#'   and the effective sample sizes (after trimming).
#'
#' @details
#' \strong{Why robust.} Under heavy-tailed or skewed marginal
#' distributions, the conventional Cohen's \emph{d} has very large
#' standard error and biased coverage. Kelley (2005) documents the
#' coverage distortion of parametric confidence intervals for the
#' standardized mean difference under nonnormal distributions.
#' Replacing means by 20\%-
#' trimmed means and SD by 20\%-Winsorized SD yields an estimator
#' whose efficiency under normality is roughly 96\% (Wilcox, 2017,
#' ch. 5) and whose efficiency under heavy-tailed contamination is
#' substantially higher than Cohen's \emph{d}.
#'
#' \strong{The 0.642 constant.} \eqn{0.642 = \mathrm{SD}(X_W) /
#' \mathrm{SD}(X) = \sqrt{\mathrm{Var}(X_W) / \mathrm{Var}(X)}} when
#' \eqn{X \sim N(0, 1)} and \eqn{X_W} is the 20\%-Winsorized version.
#' Choosing this constant makes \eqn{d_R = \delta} when the data are
#' normal, so the new estimator is on the same scale as Cohen's
#' \emph{d}.
#'
#' \strong{CI.} The CI follows the construction of Keselman, Algina,
#' Lix, Wilcox, and Deering (2008): Yuen's (1974) \emph{t}-statistic
#' on the trimmed-mean difference (their Equation 8) is referred to a
#' noncentral \emph{t} distribution with the Yuen-Welch approximate
#' degrees of freedom (their Equation 9), the noncentrality
#' parameters whose tail probabilities bracket the observed statistic
#' are located with \code{\link{conf_limits_nct}}, and those limits
#' are rescaled to the \eqn{d_R} metric. The degrees of freedom are
#' reported in the \code{df_yuen} row of the returned table. At
#' \code{trim = 0} the construction reduces to the Welch approximate
#' degrees of freedom interval; for the exact equal-variance interval
#' on the untrimmed standardized mean difference use
#' \code{\link{ci_smd}}.
#'
#' @references
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An
#'   alternative to Cohen's standardized mean difference effect
#'   size: A robust parameter and confidence interval in the two
#'   independent groups case. \emph{Psychological Methods, 10}(3),
#'   317--328. \doi{10.1037/1082-989X.10.3.317}
#'
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence
#'   intervals around the standardized mean difference: Bootstrap and
#'   parametric confidence intervals.
#'   \emph{Educational and Psychological Measurement, 65}(1), 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Keselman, H. J., Algina, J., Lix, L. M., Wilcox, R. R., & Deering,
#'   K. N. (2008). A generally robust approach for testing hypotheses
#'   and setting confidence intervals for effect sizes.
#'   \emph{Psychological Methods, 13}(2), 110--129.
#'   \doi{10.1037/1082-989X.13.2.110}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Wilcox, R. R. (2017). \emph{Introduction to robust estimation and
#'   hypothesis testing} (4th ed.). Academic Press.
#'
#' Yuen, K. K. (1974). The two-sample trimmed \emph{t} for unequal
#'   population variances. \emph{Biometrika, 61}(1), 165--170.
#'
#' @seealso \code{\link{smd}}, \code{\link{var_smd_trimmed}},
#'   \code{\link{ci_smd}}, \code{\link{conf_limits_nct}}
#'
#' @examples
#' # 1. Two normal groups: robust d agrees closely with Cohen's d.
#' set.seed(113)
#' x <- rnorm(40, 0,   1); y <- rnorm(40, 0.5, 1)
#' smd_trimmed(x, y)
#'
#' # 2. Contaminated y: a few outliers; robust d shifts much less
#' #        than Cohen's d.
#' set.seed(113)
#' x <- rnorm(40, 0, 1)
#' y <- c(rnorm(38, 0.5, 1), 30, -25)
#' smd_trimmed(x, y)
#' smd(x, y)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family effect size estimates
#'
#' @export

smd_trimmed <- function(x, y, trim = 0.20, conf_level = 0.95) {
  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must be numeric vectors.")
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 0.5)
    stop("'trim' must be a single value in [0, 0.5).")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  n_1 <- length(x); n_2 <- length(y)
  if (n_1 < 4L || n_2 < 4L)
    stop("Each group needs at least 4 non-missing observations.")

  m_t1 <- mean(x, trim = trim); m_t2 <- mean(y, trim = trim)
  s_W1 <- .winsorized_sd(x, trim); s_W2 <- .winsorized_sd(y, trim)

  h_1 <- n_1 - 2L * floor(trim * n_1)
  h_2 <- n_2 - 2L * floor(trim * n_2)

  s_W_pool <- sqrt(((h_1 - 1) * s_W1^2 + (h_2 - 1) * s_W2^2) /
                   (h_1 + h_2 - 2))

  c_const <- .akp_constant(trim)
  d_R <- c_const * (m_t1 - m_t2) / s_W_pool

  # CI (Keselman, Algina, Lix, Wilcox, & Deering, 2008): Yuen's
  # (1974) t statistic on the trimmed-mean difference (their Equation
  # 8, whose denominator carries the scaled Winsorized variances of
  # their Equation 7) is referred to a noncentral t distribution with
  # the Yuen-Welch approximate degrees of freedom (their Equation 9),
  # and the noncentrality limits are rescaled to the d_R metric
  # through the identity d_R = scale_to_d * t_yuen.
  v_1 <- s_W1^2 * (n_1 - 1) / (h_1 * (h_1 - 1))
  v_2 <- s_W2^2 * (n_2 - 1) / (h_2 * (h_2 - 1))
  t_yuen <- (m_t1 - m_t2) / sqrt(v_1 + v_2)
  df_y   <- (v_1 + v_2)^2 / (v_1^2 / (h_1 - 1) + v_2^2 / (h_2 - 1))
  scale_to_d <- c_const * sqrt(v_1 + v_2) / s_W_pool
  lims <- conf_limits_nct(ncp = t_yuen, df = df_y,
                          conf_level = conf_level, verbose = FALSE)
  lo <- scale_to_d * lims$value[lims$term == "lower_limit"]
  hi <- scale_to_d * lims$value[lims$term == "upper_limit"]

  out <- data.frame(
    term  = c("smd_trimmed", "lower_limit", "upper_limit",
              "trimmed_mean_x", "trimmed_mean_y",
              "winsorized_sd_x", "winsorized_sd_y",
              "winsorized_sd_pooled",
              "h_x", "h_y", "trim", "df_yuen"),
    value = c(d_R, lo, hi, m_t1, m_t2, s_W1, s_W2, s_W_pool,
              h_1, h_2, trim, df_y),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}

# ---- internals ----

.winsorized_sd <- function(z, trim) {
  z <- sort(z); n <- length(z); g <- floor(trim * n)
  if (g > 0) {
    z[seq_len(g)]          <- z[g + 1L]
    z[(n - g + 1L):n]      <- z[n - g]
  }
  stats::sd(z)
}

# Algina-Keselman-Penfield constant: SD(X_W) / SD(X), that is
# sqrt(Var(X_W) / Var(X)), for X ~ N(0, 1) Winsorized at 'trim' in
# each tail; 0.642 at trim = 0.20.
.akp_constant <- function(trim) {
  if (trim == 0) return(1)
  z <- stats::qnorm(1 - trim)
  win_var <- 1 - 2 * trim - 2 * z * stats::dnorm(z) +
             2 * trim * z^2
  sqrt(win_var)
}
