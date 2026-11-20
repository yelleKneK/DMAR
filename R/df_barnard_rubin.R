#' Barnard and Rubin's Small-Sample Degrees of Freedom for Multiple
#' Imputation
#'
#' Computes the degrees of freedom of the \emph{t} reference distribution
#' for an estimate pooled over \emph{m} imputed data sets by Rubin's (1987)
#' rules, using Barnard and Rubin's (1999) small-sample formula, together
#' with the quantities it is built from: the total variance, the fraction
#' of missing information the imputations estimate, the relative increase
#' in variance due to the missing values, and the two components of the
#' degrees of freedom. The function is vectorized over parameters.
#'
#' \strong{Rubin's rules.} Each of the \emph{m} completed data sets gives
#' an estimate \eqn{\hat{\theta}_d}{theta_d} and a squared standard error
#' \eqn{U_d}. The pooled estimate \eqn{\bar{\theta}}{theta_bar} is their
#' mean, the within-imputation variance is
#' \eqn{\bar{U} = m^{-1} \sum_d U_d}{U_bar = sum(U_d) / m}, the
#' between-imputation variance is
#' \eqn{B = (m - 1)^{-1} \sum_d (\hat{\theta}_d - \bar{\theta})^2}{B = sum((theta_d - theta_bar)^2) / (m - 1)},
#' and the total variance of the pooled estimate is
#' \deqn{T = \bar{U} + (1 + 1/m) B.}{T = U_bar + (1 + 1/m) B.}
#' The share of the total variance due to the missing values,
#' \deqn{\gamma = (1 + 1/m) B / T,}
#' is the multiple imputation estimate of the fraction of missing
#' information, which \pkg{mice} reports as \code{lambda}, and
#' \eqn{r = (1 + 1/m) B / \bar{U}}{r = (1 + 1/m) B / U_bar} is the relative
#' increase in variance due to nonresponse (\pkg{mice}'s \code{riv}).
#' \pkg{mice}'s \code{fmi} column is a different quantity,
#' \eqn{(r + 2/(\nu + 3))/(r + 1)}, which adjusts \eqn{\gamma} for the
#' degrees of freedom. Von Hippel's (2016) Eq. 25 prints the
#' within-imputation variance in place of \eqn{B} in the numerator of
#' \eqn{\gamma}, a typographical error: that ratio would exceed 1 with no
#' missing values at all.
#'
#' \strong{The degrees of freedom.} Rubin's (1987) degrees of freedom,
#' \deqn{\nu_m = (m - 1) / \gamma^2,}
#' reflect only the finite number of imputations: they treat each
#' completed data analysis as if it had infinitely many degrees of
#' freedom, and when \eqn{\gamma} is small they can exceed the degrees of
#' freedom the analysis would have with no missing values at all, which
#' cannot be right in a small sample. Barnard and Rubin (1999) combine
#' them with observed data degrees of freedom,
#' \deqn{\nu_{obs} = \nu_{com} (1 - \gamma) \frac{\nu_{com} + 1}{\nu_{com} + 3},}{nu_obs = nu_com (1 - gamma) (nu_com + 1) / (nu_com + 3),}
#' where \eqn{\nu_{com}}, given as \code{df_complete}, is the complete
#' data degrees of freedom (for a regression coefficient, \eqn{n - k},
#' with \eqn{k} the number of regression parameters), as
#' \deqn{\nu_{BR} = \left(\frac{1}{\nu_m} + \frac{1}{\nu_{obs}}\right)^{-1},}{nu_BR = 1 / (1/nu_m + 1/nu_obs),}
#' which never exceeds \eqn{\nu_{obs}} and hence never exceeds
#' \eqn{\nu_{com}}. The formulas are as von Hippel (2016, Eqs. 23 and 24)
#' writes them, and they agree with \pkg{mice}'s \code{barnard.rubin()}.
#' With \code{df_complete = Inf}, the default, \eqn{\nu_{obs}} is infinite and
#' \eqn{\nu_{BR}} is Rubin's \eqn{\nu_m}. The observed data degrees of
#' freedom are \code{\link{df_missing}}'s \code{"complete_df"} rule, and
#' von Hippel's degrees of freedom for maximum likelihood are the limit of
#' \eqn{\nu_{BR}} as the number of imputations grows, where
#' \eqn{\nu_m} becomes infinite.
#'
#' \strong{The floor.} An estimate of \eqn{\gamma} near 1 drives the degrees
#' of freedom toward 0, where the \emph{t} quantile grows without bound. In
#' von Hippel's (2016, p. 431) simulation of 25 cases with half the values
#' of one variable missing at random, 24\% of the Barnard and Rubin
#' degrees of freedom were below 3, 11\% below 2, and 2\% below 1, and some
#' intervals were orders of magnitude longer than coverage required. He
#' bounds the degrees of freedom below at 3 (his Eq. 27), which is the
#' default \code{floor = 3} here; \code{floor = 0} reproduces Barnard and
#' Rubin's formula, and \pkg{mice}, exactly. The unbounded value is kept in
#' \code{df_unbounded}, and \code{at_floor} marks the rows where the bound
#' applied.
#'
#' \pkg{DMAR} does not impute; this function serves analyses pooled
#' elsewhere and the comparison of maximum likelihood with multiple
#' imputation. For a FIML fit, \code{\link{df_missing}} gives the
#' corresponding degrees of freedom from the fit's own fraction of missing
#' information, without imputations.
#'
#' @param between The between-imputation variance \eqn{B} of each
#'   parameter: the variance, with divisor \eqn{m - 1}, of its \emph{m}
#'   estimates. A numeric vector, one value per parameter, whose names
#'   label the rows. (It is not called \code{B} because \pkg{DMAR} keeps
#'   \code{B} for the number of bootstrap replications.)
#' @param within The within-imputation variance \eqn{\bar{U}}{U_bar} of each
#'   parameter: the mean of its \emph{m} squared standard errors. A numeric
#'   vector of the same length as \code{between}.
#' @param m The number of imputations, a single whole number of 2 or more.
#' @param df_complete The complete data degrees of freedom, of length 1 or the
#'   length of \code{between}. Defaults to \code{Inf}, which gives Rubin's
#'   (1987) large-sample degrees of freedom; a small-sample analysis should
#'   supply its complete data error degrees of freedom.
#' @param floor The lower bound on the degrees of freedom: a single number
#'   of 0 or more. Defaults to 3, von Hippel's (2016) recommendation;
#'   \code{0} reproduces Barnard and Rubin (1999).
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   parameter and the columns \code{term}, \code{df} (bounded below by
#'   \code{floor}), \code{df_unbounded} (\eqn{\nu_{BR}}), \code{gamma},
#'   \code{relative_increase} (\eqn{r}), \code{total_variance}
#'   (\eqn{T}), \code{nu_m}, \code{nu_obs}, and \code{at_floor}. The
#'   attributes \code{m}, \code{df_complete}, and \code{floor} record the call,
#'   and the printed table, at the console or knitted, is followed by a
#'   line naming the rule, \emph{m}, \code{df_complete}, and the floor. A
#'   missing input value gives missing values in its row.
#'
#' @references
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#'   multiple imputation. \emph{Biometrika, 86}(4), 948--955.
#'   \doi{10.1093/biomet/86.4.948}
#'
#' Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in
#'   surveys}. Wiley.
#'
#' van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
#'   imputation by chained equations in R. \emph{Journal of Statistical
#'   Software, 45}(3), 1--67. \doi{10.18637/jss.v045.i03}
#'
#' von Hippel, P. T. (2016). New confidence intervals and bias comparisons
#'   show that maximum likelihood can beat multiple imputation in small
#'   samples. \emph{Structural Equation Modeling, 23}(3), 422--437.
#'   \doi{10.1080/10705511.2015.1047931}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{df_missing}} for the degrees of freedom of a FIML
#'   estimate.
#'
#' @family missing data
#'
#' @keywords htest
#'
#' @examples
#' # One coefficient from m = 5 imputations, with complete data error
#' # degrees of freedom of 23: B = 0.04 and a within-imputation variance of
#' # 0.10. floor = 0 gives Barnard and Rubin's formula itself.
#' df_barnard_rubin(between = 0.04, within = 0.10, m = 5, df_complete = 23,
#'                  floor = 0)
#'
#' # Rubin's (1987) large-sample degrees of freedom are 38, well above the
#' # 23 a complete data analysis of these cases would have.
#' df_barnard_rubin(between = 0.04, within = 0.10, m = 5)
#'
#' # Pooling by hand: five estimates of each of two coefficients and their
#' # squared standard errors (illustrative values), one column per
#' # coefficient.
#' est <- cbind(slope = c(0.42, 0.51, 0.38, 0.47, 0.55),
#'              intercept = c(1.90, 1.94, 1.88, 1.93, 1.91))
#' u <- cbind(slope = c(0.010, 0.011, 0.009, 0.010, 0.012),
#'            intercept = c(0.020, 0.021, 0.019, 0.020, 0.022))
#' df_barnard_rubin(between = apply(est, 2, var), within = colMeans(u),
#'                  m = 5, df_complete = 23)
#'
#' @export
df_barnard_rubin <- function(between, within, m, df_complete = Inf, floor = 3) {
  .df_check_floor(floor)
  if (!is.numeric(between) || !length(between)) {
    .dmar_refuse_argument(c(
      "'between' must be a numeric vector of between-imputation ",
      "variances."), arg = "between")
  }
  if (!is.numeric(within) || length(within) != length(between)) {
    .dmar_refuse_argument(c(
      "'within' must be a numeric vector of within-imputation ",
      "variances, one per value of 'between'."), arg = c("within", "between"))
  }
  if (any(is.infinite(between)) || any(is.infinite(within))) {
    .dmar_refuse_argument(
      "'between' and 'within' must be finite.", arg = c("between", "within"))
  }
  if (any(between < 0, na.rm = TRUE) || any(within < 0, na.rm = TRUE)) {
    .dmar_refuse_argument(
      "'between' and 'within' are variances and cannot be negative.",
      arg = c("between", "within"))
  }
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m < 2 || m != round(m)) {
    .dmar_refuse_argument(c(
      "'m' must be a single whole number of 2 or more: the ",
      "between-imputation variance needs at least two imputations."), arg = "m")
  }
  n <- length(between)
  if (!is.numeric(df_complete) || anyNA(df_complete) || any(df_complete <= 0)) {
    .dmar_refuse_argument(c(
      "'df_complete' must be positive (Inf, the default, gives Rubin's ",
      "large-sample degrees of freedom)."), arg = "df_complete")
  }
  df_complete <- .df_recycle(as.numeric(df_complete), n, "df_complete")

  total <- within + (1 + 1 / m) * between
  if (any(total <= 0, na.rm = TRUE)) {
    .dmar_refuse_argument(c(
      "The total variance is 0 for some parameter (both 'between' and ",
      "'within' are 0); its degrees of freedom are undefined."),
      arg = c("between", "within"))
  }
  gamma <- (1 + 1 / m) * between / total
  relative_increase <- (1 + 1 / m) * between / within
  # Rubin's (m - 1) / gamma^2 is infinite when the imputations agree
  # exactly (between = 0); 1 / Inf = 0 then leaves the observed data
  # degrees of freedom alone in the harmonic combination.
  nu_m <- (m - 1) / gamma^2
  nu_obs <- .df_complete(gamma, df_complete)
  df_raw <- 1 / (1 / nu_m + 1 / nu_obs)

  bounded <- .df_bound(df_raw, floor)
  term <- if (!is.null(names(between))) {
    .df_terms(between, n)
  } else {
    .df_terms(within, n)
  }
  out <- data.frame(term              = term,
                    df                = bounded$df,
                    df_unbounded      = df_raw,
                    gamma             = gamma,
                    relative_increase = relative_increase,
                    total_variance    = total,
                    nu_m              = nu_m,
                    nu_obs            = nu_obs,
                    at_floor          = bounded$at_floor,
                    stringsAsFactors = FALSE)
  attr(out, "m") <- m
  attr(out, "df_complete") <- df_complete
  attr(out, "floor") <- floor
  .as_dmar_tbl(out, subclass = "dmar_df_barnard_rubin")
}

# The rule behind the df column and its floor, in one paragraph under the
# table, at the console and in a knitted report alike (see
# .dmar_tbl_footer() in R/dmar_output_helpers.R). A column subset, which
# loses the m, df_complete, and floor attributes, prints no line.
#' @exportS3Method .dmar_tbl_footer dmar_df_barnard_rubin
#' @noRd
.dmar_tbl_footer.dmar_df_barnard_rubin <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .df_barnard_rubin_footer(x), "df_rule")
}

# The line itself: with every df_complete infinite, Rubin's (1987) degrees of
# freedom, which nu_BR then equals; otherwise Barnard and Rubin's (1999)
# with the df_complete values used. Broken for the console (a knitted report
# reads the break as a space), with the floor.
.df_barnard_rubin_footer <- function(x) {
  m <- attr(x, "m", exact = TRUE)
  df_complete <- attr(x, "df_complete", exact = TRUE)
  floor <- attr(x, "floor", exact = TRUE)
  if (is.null(m) || !length(df_complete) || is.null(floor)) return(NULL)
  fmt <- function(v) format(v, digits = 6, drop0trailing = TRUE, trim = TRUE)
  nu <- vapply(unique(df_complete), fmt, character(1L))
  rule <- if (all(is.infinite(df_complete))) {
    sprintf(paste0("Rubin's (1987) large-sample rule for m = %s ",
                   "imputations,\n(m - 1) / gamma^2 (df_complete = Inf)"),
            fmt(m))
  } else {
    k <- length(nu)
    if (k > 1L) {
      nu <- paste0(paste(nu[-k], collapse = ", "),
                   if (k > 2L) ", and " else " and ", nu[k])
    }
    sprintf(paste0("Barnard and Rubin's (1999) rule for m = %s ",
                   "imputations,\n1 / (1 / nu_m + 1 / nu_obs) with ",
                   "df_complete = %s"), fmt(m), nu)
  }
  sprintf("Degrees of freedom: %s, bounded below at %s", rule, fmt(floor))
}


# The tidy() and glance() methods for df_barnard_rubin(), documented at
# ?dmar_tidiers; see the note above tidy.dmar_df_missing().

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_df_barnard_rubin <- function(x, ...) {
  .df_tidy(x)
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_df_barnard_rubin <- function(x, ...) {
  data.frame(
    n_terms    = nrow(x),
    n_at_floor = sum(x$at_floor, na.rm = TRUE),
    m          = attr(x, "m") %||% NA_real_,
    floor      = attr(x, "floor") %||% NA_real_,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}
