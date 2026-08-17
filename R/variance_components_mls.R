# Modified-large-sample CIs on variance components.
#' Modified-Large-Sample Confidence Intervals on Variance Components
#'
#' Computes modified-large-sample (MLS) confidence intervals on the
#' between-group and within-group variance components of a balanced
#' one-way random-effects ANOVA, following Burdick & Graybill (1992).
#' MLS intervals have substantially better coverage than the
#' Satterthwaite or simple-Wald intervals when the components are far
#' from zero, and they are the standard interval method in
#' generalizability theory (Brennan, 2001).
#'
#' @param ms_between Mean square between groups (numerator of the
#'   ANOVA \emph{F}).
#' @param ms_within Mean square within groups (denominator of the
#'   ANOVA \emph{F}).
#' @param df_between Degrees of freedom for the between-group MS
#'   (typically \eqn{a - 1} for \eqn{a} groups).
#' @param df_within Degrees of freedom for the within-group MS
#'   (typically \eqn{a (n - 1)} for \eqn{a} groups of size \eqn{n}).
#' @param n Number of observations per group (assumed balanced).
#' @param conf_level Confidence level for the CIs. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the point estimates
#'   and MLS lower / upper CIs of the between-group variance component
#'   (\eqn{\sigma^2_b}), the within-group component (\eqn{\sigma^2_w}),
#'   and the implied intraclass correlation
#'   (\eqn{\rho = \sigma^2_b / (\sigma^2_b + \sigma^2_w)}).
#'
#' @details
#' \strong{Point estimates.} For a one-way random-effects ANOVA on
#' \eqn{a} groups of size \eqn{n}, the method-of-moments estimators are
#' \deqn{\hat\sigma^2_b \;=\; \max(0,\, (\mathit{MS}_b - \mathit{MS}_w)/n),
#'   \qquad \hat\sigma^2_w \;=\; \mathit{MS}_w.}
#'
#' \strong{Modified-large-sample CIs (Burdick-Graybill 1992).} The MLS
#' interval for \eqn{\sigma^2_b} is
#' \deqn{\left[\frac{\mathit{MS}_b - \mathit{MS}_w - \sqrt{V_L}}{n},\;\;
#'   \frac{\mathit{MS}_b - \mathit{MS}_w + \sqrt{V_U}}{n}\right],}
#' with
#' \deqn{V_L \;=\; G_1^2 \mathit{MS}_b^2 + H_2^2 \mathit{MS}_w^2
#'     + G_{12} \mathit{MS}_b \mathit{MS}_w, \qquad
#'   V_U \;=\; H_1^2 \mathit{MS}_b^2 + G_2^2 \mathit{MS}_w^2
#'     + H_{12} \mathit{MS}_b \mathit{MS}_w,}
#' where the constants \eqn{G_1}, \eqn{G_2}, \eqn{H_1}, \eqn{H_2} and
#' the cross-term constants \eqn{G_{12}}, \eqn{H_{12}} depend on the
#' degrees of freedom and on \eqn{\chi^2} and \emph{F} quantiles at the
#' chosen confidence level (Burdick & Graybill, 1992, equations
#' 2.4.1--2.4.5 give the explicit formulas). The lower limit is
#' truncated at zero. For the within-group component, the standard
#' \eqn{\chi^2}-based CI on \eqn{\mathit{MS}_w} (Searle, Casella, &
#' McCulloch, 1992) is used.
#'
#' \strong{Caveats.} MLS intervals assume balanced data and homogeneous
#' variances within groups. For unbalanced data the appropriate analog
#' is the Burdick-Graybill MLS extension to unequal sample sizes
#' (Burdick & Graybill, 1992, Section 2.5), which is not implemented
#' here.
#'
#' @references
#' Brennan, R. L. (2001). \emph{Generalizability theory}. Springer.
#'
#' Burdick, R. K., & Graybill, F. A. (1992). \emph{Confidence intervals
#'   on variance components}. Marcel Dekker.
#'
#' Searle, S. R., Casella, G., & McCulloch, C. E. (1992).
#'   \emph{Variance components}. Wiley.
#'
#' @seealso \code{\link{icc}}, \code{\link{var_icc}}, \code{\link{ss_aipe_icc}}
#'
#' @examples
#' # 1. Balanced one-way random-effects ANOVA: a = 10 groups, n = 5.
#' #        Hypothetical MS_b = 6.0, MS_w = 1.5.
#' variance_components_mls(ms_between = 6.0, ms_within = 1.5,
#'                         df_between = 9, df_within = 40, n = 5)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family agreement and measurement
#'
#' @export

variance_components_mls <- function(ms_between, ms_within,
                                    df_between, df_within, n,
                                    conf_level = 0.95) {
  if (!is.numeric(ms_between) || ms_between < 0)
    stop("'ms_between' must be a non-negative number.")
  if (!is.numeric(ms_within) || ms_within < 0)
    stop("'ms_within' must be a non-negative number.")
  if (!is.numeric(df_between) || df_between < 1)
    stop("'df_between' must be >= 1.")
  if (!is.numeric(df_within) || df_within < 1)
    stop("'df_within' must be >= 1.")
  if (!is.numeric(n) || n < 2)
    stop("'n' (per-group sample size) must be >= 2.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  alpha <- 1 - conf_level

  # Point estimates:
  sigma2_b <- max(0, (ms_between - ms_within) / n)
  sigma2_w <- ms_within
  rho_hat  <- if (sigma2_b + sigma2_w > 0)
                sigma2_b / (sigma2_b + sigma2_w) else 0

  # Burdick-Graybill (1992) MLS constants.
  G1 <- 1 - 1 / stats::qf(1 - alpha / 2, df_between, Inf)
  G2 <- 1 - 1 / stats::qf(1 - alpha / 2, df_within,  Inf)
  H1 <- 1 / stats::qf(alpha / 2, df_between, Inf) - 1
  H2 <- 1 / stats::qf(alpha / 2, df_within,  Inf) - 1
  G12 <- ((stats::qf(1 - alpha / 2, df_between, df_within) - 1)^2 -
          (G1 * stats::qf(1 - alpha / 2, df_between, df_within))^2 - H2^2) /
         stats::qf(1 - alpha / 2, df_between, df_within)
  H12 <- ((1 - stats::qf(alpha / 2, df_between, df_within))^2 -
          (H1 * stats::qf(alpha / 2, df_between, df_within))^2 - G2^2) /
         stats::qf(alpha / 2, df_between, df_within)

  V_L <- G1^2 * ms_between^2 + H2^2 * ms_within^2 + G12 * ms_between * ms_within
  V_U <- H1^2 * ms_between^2 + G2^2 * ms_within^2 + H12 * ms_between * ms_within

  lo_b <- (ms_between - ms_within - sqrt(V_L)) / n
  hi_b <- (ms_between - ms_within + sqrt(V_U)) / n
  lo_b <- max(0, lo_b)
  hi_b <- max(lo_b, hi_b)

  # Within-component CI: chi squared on MS_w * df_within / sigma^2_w.
  lo_w <- df_within * ms_within / stats::qchisq(1 - alpha / 2, df_within)
  hi_w <- df_within * ms_within / stats::qchisq(alpha / 2,     df_within)

  # ICC implied bounds (transform-based, conservative):
  rho_lo <- if (sigma2_w > 0)
              max(0, lo_b / (lo_b + hi_w)) else NA_real_
  rho_hi <- if (sigma2_w > 0)
              min(1, hi_b / (hi_b + lo_w)) else NA_real_

  out <- data.frame(
    term  = c("sigma2_between", "sigma2_b_lower", "sigma2_b_upper",
              "sigma2_within",  "sigma2_w_lower", "sigma2_w_upper",
              "icc", "icc_lower", "icc_upper"),
    value = c(sigma2_b, lo_b, hi_b,
              sigma2_w, lo_w, hi_w,
              rho_hat, rho_lo, rho_hi),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
