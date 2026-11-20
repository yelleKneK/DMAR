#' Small-Sample Degrees of Freedom for Maximum Likelihood Intervals With
#' Missing Data
#'
#' Computes the degrees of freedom of a \emph{t} reference distribution for
#' a Wald confidence interval from a full information maximum likelihood
#' (FIML) fit, from each parameter's fraction of missing information (FMI),
#' by either of von Hippel's (2016) two rules. The interval keeps the FIML
#' estimate and its standard error and replaces the normal quantile with
#' the \emph{t} quantile on these degrees of freedom, which lengthens it by
#' an amount that grows as the sample shrinks and as the share of the
#' information lost to missing values grows.
#'
#' \strong{Why a \emph{t} reference distribution.} A maximum likelihood
#' Wald interval is the estimate plus or minus a normal quantile times the
#' standard error; the normal quantile comes from the large-sample
#' distribution of the estimate. With incomplete data in small samples that
#' interval is too short. In
#' von Hippel's (2016) simulations of bivariate normal data with half the
#' values of one variable missing, the normal FIML intervals at
#' \emph{n} = 25 covered 5 to 6 points below the nominal 95\% for the
#' intercept and slope of a regression and 8 to 14 points below it for the
#' residual variance, and at \emph{n} = 100 they were 0 to 3 points below
#' (pp. 432--433). What fails is the reference distribution, not the
#' estimate or its standard error: the normal quantile treats the standard
#' error as known, and with few cases, or few effective cases once the
#' missing information is subtracted, the standard error is itself
#' estimated with error, as it is in a complete data regression whose
#' coefficients are referred to \emph{t} on \eqn{N - k} degrees of freedom.
#' The \emph{t} distribution is an approximation here, since the
#' small-sample distribution of a FIML estimate is not exactly \emph{t}
#' (von Hippel, 2016, p. 431). Von Hippel's \emph{t} intervals brought the
#' coverage of most parameters to within 2 points of nominal (p. 434). At
#' \emph{n} = 25 it stayed below nominal for the parameters whose point
#' estimates are biased, because a longer interval covers the variability
#' of an estimate, not its bias.
#'
#' \strong{The two rules.} With \eqn{\gamma} the parameter's fraction of
#' missing information, \code{method = "effective_n"} (von Hippel's
#' \eqn{\nu^*_{ML}}, his Eq. 31) is
#' \deqn{\nu^* = N(1 - \gamma) - k,}
#' the complete data regression degrees of freedom \eqn{N - k} with
#' \eqn{N} replaced by the effective sample size \eqn{N(1 - \gamma)}, the
#' number of complete cases that would carry the same information about
#' the parameter (Savalei & Rhemtulla, 2012, p. 488). Here \eqn{k} is the
#' number of regression parameters, slopes and intercept. With no missing
#' information it is \eqn{N - k} exactly. \code{method = "complete_df"}
#' (von Hippel's \eqn{\nu_{ML}}, his Eq. 32) is
#' \deqn{\nu = \nu_{com} (1 - \gamma) \frac{\nu_{com} + 1}{\nu_{com} + 3},}{nu = nu_com (1 - gamma) (nu_com + 1) / (nu_com + 3),}
#' where \eqn{\nu_{com}}, given as \code{df_complete}, is the degrees of
#' freedom the interval would use if no value were missing: the error
#' degrees of freedom, not the model degrees of freedom of a structural
#' equation model (p. 423). It is the
#' limit of Barnard and Rubin's (1999) degrees of freedom for multiple
#' imputation as the number of imputations grows, which von Hippel derives
#' from the large-sample equivalence of maximum likelihood and multiple
#' imputation with infinitely many imputations, and it is the observed
#' data part of that formula (see \code{\link{df_barnard_rubin}}). With no
#' missing information it gives \eqn{\nu_{com}(\nu_{com} + 1)/(\nu_{com} +
#' 3)} rather than \eqn{\nu_{com}} (21.2 rather than 23 when
#' \eqn{\nu_{com}} = 23); the factor comes from the observed data degrees
#' of freedom and does not vanish with the missing information. The two
#' rules agree closely: for a simple regression with \emph{n} = 25,
#' \emph{k} = 2, and \eqn{\gamma} = .5 they give 10.6 and 10.5 (p. 432).
#' Von Hippel's Table 5 used the effective sample size rule, and he reports
#' very similar results with the general one (p. 433).
#'
#' \strong{The floor.} Either rule can return a small or negative value
#' when the FMI is estimated near 1 or, for \code{"effective_n"}, when the
#' effective sample size approaches \emph{k}. A \emph{t} distribution with
#' 2 or fewer degrees of freedom has no finite variance, and its quantiles
#' grow without bound as the degrees of freedom approach 0 (a 95\%
#' interval on 1 degree of freedom uses 12.7), so von Hippel bounds the
#' degrees of freedom below at 3 (his Eq. 34). \code{floor = 3} is that
#' bound. The unbounded value is kept in \code{df_unbounded}, and
#' \code{at_floor} marks the rows where the bound applied. \code{floor = 0}
#' returns the unbounded \code{"complete_df"} value, which is never
#' negative; a negative \code{"effective_n"} value becomes 0 there, which
#' is not a usable reference distribution.
#'
#' \strong{Where the FMI comes from.} The FMI of a parameter is the share
#' of its sampling variance due to the missing values,
#' \eqn{\gamma = 1 - SE_C^2 / SE_O^2}, with \eqn{SE_C} its standard error
#' had the data been complete and \eqn{SE_O} its standard error from the
#' incomplete data (Savalei & Rhemtulla, 2012, Eq. 4). It differs from
#' parameter to parameter, and so do the degrees of freedom: a parameter
#' estimated mostly from fully observed variables keeps nearly all of its
#' degrees of freedom however many values are missing elsewhere.
#' \code{\link{fraction_missing_information}} estimates it from a FIML fit, and
#' its table can be passed as \code{fmi} directly; \pkg{lavaan}'s
#' \code{parameterEstimates(fit, fmi = TRUE)} gives another estimate. Some
#' estimators return negative FMIs in small samples (Chen & Savalei, 2021);
#' a negative FMI has no degrees of freedom interpretation, so a value
#' below \eqn{-10^{-8}} is refused, and a value between that and 0, which
#' is rounding error around a parameter with no missing information, is
#' set to 0.
#'
#' \strong{Which parameters the rules are supported for.} Von Hippel's
#' simulations used bivariate normal data with correlation .5, \emph{n} =
#' 25 and 100, and half the values of one variable missing, completely at
#' random or whenever the other variable was negative (missing at random).
#' The intervals evaluated were for the intercepts and slopes of the
#' regression of each variable on the other, the mean of the incomplete
#' variable, the covariance, and the variances (on a cube root scale). The
#' rules therefore have published support for regression coefficients and
#' for moments: means, variances, and covariances. They have none for
#' factor loadings or for reliability coefficients such as omega, whose
#' complete data intervals are themselves asymptotic, so that no complete
#' data degrees of freedom are defined for them. For those, DMAR's own
#' calibration study (in the package's source repository,
#' \code{methodology/missing_data}; \code{\link{fa_confirmatory}} describes its
#' design) simulated one factor models fit by FIML and found that
#' \code{method = "effective_n"} with \code{k = 1} brought the coverage
#' of the loadings and of omega closer to .95 than the normal reference
#' did: no cell's coverage, averaged over its loadings, exceeded .9502,
#' and omega's .9550 (single loadings reached .9595 in the cells with
#' missing values and .963 in the complete ones, where the rule is
#' \emph{t} on \eqn{N - 1} and the normal reference reached .962); it is
#' the rule \code{\link{fa_confirmatory}} and
#' \code{\link{reliability_omega}} apply under \code{wald_df =
#' "missing_information"}. For coefficient alpha no rule passed the
#' study's decision rule, fixed before the run, which a coverage of .966
#' in one cell decided (see \code{\link{reliability_alpha}}).
#'
#' These numerical results are conditional on the H1 eligibility rule the
#' study applied. The affine-collapse diagnostic sets aside additional
#' draws, and the study has not established coverage under that selection.
#' A raw-data eligibility replay alone cannot validate the FMI refits and
#' information diagnostics. The keyword is an approximate, opt-in rule, not
#' a coverage guarantee.
#'
#' The function has no default for \code{k} or \code{df_complete}, because
#' the right value depends on the parameter.
#'
#' @param fmi The fraction of missing information of each parameter: a
#'   numeric vector of values below 1, or the table returned by
#'   \code{\link{fraction_missing_information}}, whose \code{fmi} and \code{N}
#'   columns are used and whose \code{term} column labels the rows. The
#'   names of a named vector label the rows. A missing value gives missing
#'   degrees of freedom for that row.
#' @param N The number of cases the fit used (with \pkg{lavaan}, the rows
#'   with at least one observed variable, \code{lavInspect(fit, "nobs")}),
#'   of length 1 or the length of \code{fmi}. Required by
#'   \code{method = "effective_n"} unless \code{fmi} is a table with an
#'   \code{N} column; not used by \code{"complete_df"}.
#' @param k The number of parameters subtracted from the effective sample
#'   size, of length 1 or the length of \code{fmi}: for a regression
#'   coefficient, the number of regression parameters (slopes and
#'   intercept), so that the complete data degrees of freedom are
#'   \eqn{N - k}; for a factor loading or omega, 1, the value DMAR's
#'   calibration chose (see \emph{Details}). Required by \code{method =
#'   "effective_n"}, with no default.
#' @param df_complete The complete data degrees of freedom, of length 1 or the
#'   length of \code{fmi}: the error degrees of freedom the interval would
#'   use if no value were missing (\eqn{n - k} for a regression
#'   coefficient, \eqn{n - 1} for a mean), not the model degrees of
#'   freedom of a structural equation model. \code{Inf} is allowed and
#'   gives the normal reference. Required by \code{method = "complete_df"},
#'   with no default.
#' @param method \code{"effective_n"} (the default; von Hippel's Eq. 31,
#'   the rule his Table 5 used) or \code{"complete_df"} (his Eq. 32).
#' @param floor The lower bound on the degrees of freedom: a single number
#'   of 0 or more. Defaults to 3, von Hippel's (2016) recommendation;
#'   \code{0} removes it.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   value of \code{fmi}, in input order, and the columns \code{term},
#'   \code{df} (the degrees of freedom, bounded below by \code{floor}),
#'   \code{df_unbounded}, \code{fmi} (the value used), \code{N},
#'   \code{k}, \code{df_complete}, and \code{at_floor} (\code{TRUE} where the
#'   bound applied). A column the method does not use is \code{NA}:
#'   \code{df_complete} under \code{"effective_n"}, \code{k} under
#'   \code{"complete_df"}, and \code{N} under \code{"complete_df"} unless a
#'   table supplied it. The attributes \code{method} and \code{floor}
#'   record the call, and the printed table, at the console or knitted,
#'   is followed by a line naming the rule and the floor.
#'
#' @references
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#'   multiple imputation. \emph{Biometrika, 86}(4), 948--955.
#'   \doi{10.1093/biomet/86.4.948}
#'
#' Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
#'   missing information from full information maximum likelihood.
#'   \emph{Frontiers in Psychology, 12}, 667802.
#'   \doi{10.3389/fpsyg.2021.667802}
#'
#' Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
#'   fraction of missing information from full information maximum
#'   likelihood. \emph{Structural Equation Modeling, 19}(3), 477--494.
#'   \doi{10.1080/10705511.2012.687669}
#'
#' von Hippel, P. T. (2016). New confidence intervals and bias comparisons
#'   show that maximum likelihood can beat multiple imputation in small
#'   samples. \emph{Structural Equation Modeling, 23}(3), 422--437.
#'   \doi{10.1080/10705511.2015.1047931}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{df_barnard_rubin}} for the degrees of freedom of a
#'   multiple imputation estimate; \code{\link{mlmr}} for FIML regression.
#'
#' @family missing data
#'
#' @keywords htest
#'
#' @examples
#' # von Hippel's (2016, p. 432) worked example: a simple regression with
#' # n = 25 cases, an intercept and a slope (k = 2), and half of the
#' # information about the slope missing.
#' df_missing(0.5, N = 25, k = 2)
#' df_missing(0.5, df_complete = 23, method = "complete_df")
#'
#' # Each parameter gets the degrees of freedom of its own FMI. With no
#' # missing information the effective sample size rule returns the
#' # complete data degrees of freedom, N - k = 23; near an FMI of 1 the
#' # floor of 3 takes over.
#' df_missing(c(slope = 0.5, intercept = 0.1, complete = 0, sparse = 0.95),
#'            N = 25, k = 2)
#'
#' # A regression whose outcome is observed for only the 145 Grant-White
#' # students of 301: t26_flags is missing for every Pasteur student. mlmr()
#' # fits it by full information maximum likelihood (with the Wald interval,
#' # whose normal quantile the t reference replaces, and without the effect
#' # size refits, which this example does not use), and
#' # fraction_missing_information() gives each slope's FMI with the N behind
#' # it, so its rows go to df_missing() as they are; k = 3 counts the two
#' # slopes and the intercept.
#' fit <- mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
#'             data = holzinger_swineford, ci_method = "wald",
#'             effect_sizes = FALSE)
#' fmi <- fraction_missing_information(fit)
#' slopes <- fmi[fmi$op == "~", ]
#' dfs <- df_missing(slopes, k = 3)
#' dfs
#'
#' # The t interval keeps the FIML estimate and standard error and replaces
#' # the normal quantile. With more than a hundred effective cases the
#' # change is small; it matters when the effective sample size is small.
#' est <- coef(fit)[-1]
#' se <- sqrt(diag(vcov(fit)))[-1]
#' crit <- qt(0.975, dfs$df)
#' data.frame(term = names(est), df = dfs$df,
#'            t_ci_lower = est - crit * se, t_ci_upper = est + crit * se,
#'            normal_ci_lower = est - qnorm(0.975) * se,
#'            normal_ci_upper = est + qnorm(0.975) * se, row.names = NULL)
#'
#' @export
df_missing <- function(fmi, N = NULL, k = NULL, df_complete = NULL,
                       method = c("effective_n", "complete_df"),
                       floor = 3) {
  method <- match.arg(method)
  .df_check_floor(floor)

  # The table from fraction_missing_information(), or any data frame with
  # an fmi column, supplies the FMIs, their row labels, and the N basis of
  # each row (a group's N for a group-specific parameter, the summed N of
  # the groups an equality constraint spans).
  N_from_table <- NULL
  if (is.data.frame(fmi)) {
    if (!"fmi" %in% names(fmi)) {
      .dmar_refuse_argument(c(
        "A data frame passed as 'fmi' must have an fmi column, as the ",
        "table fraction_missing_information() returns does."), arg = "fmi")
    }
    term <- if ("term" %in% names(fmi)) {
      as.character(fmi$term)
    } else {
      paste0("parameter_", seq_len(nrow(fmi)))
    }
    if ("N" %in% names(fmi)) N_from_table <- as.numeric(fmi$N)
    gamma <- .df_clean_fmi(fmi$fmi)
  } else {
    gamma <- .df_clean_fmi(fmi)
    term <- .df_terms(fmi, length(gamma))
  }
  n <- length(gamma)

  if (method == "effective_n") {
    if (!is.null(df_complete)) {
      .dmar_refuse_argument(c(
        "'df_complete' is used by method = \"complete_df\"; ",
        "method = \"effective_n\" takes 'N' and 'k'."),
        arg = c("df_complete", "N", "k"))
    }
    if (!is.null(N) && !is.null(N_from_table)) {
      .dmar_refuse_argument(c(
        "'fmi' is a table with an N column, which gives the N of each ",
        "row; leave 'N' unset, or pass the fmi column alone with 'N'."),
        arg = c("fmi", "N"))
    }
    if (is.null(N)) N <- N_from_table
    if (is.null(N)) {
      .dmar_refuse_argument(c(
        "method = \"effective_n\" needs 'N', the number of cases the ",
        "fit used."), arg = "N")
    }
    if (is.null(k)) {
      .dmar_refuse_argument(c(
        "method = \"effective_n\" needs 'k', the number of parameters ",
        "subtracted from the effective sample size (for a regression ",
        "coefficient, the number of regression parameters, slopes and ",
        "intercept; for a factor loading or omega, 1). It has no ",
        "default, because the right value depends on the parameter."),
        arg = "k")
    }
    if (!is.numeric(N) || anyNA(N) || any(!is.finite(N)) || any(N <= 0)) {
      .dmar_refuse_argument("'N' must be finite and positive.", arg = "N")
    }
    if (!is.numeric(k) || anyNA(k) || any(!is.finite(k)) || any(k < 0) ||
        any(k != round(k))) {
      .dmar_refuse_argument(
        "'k' must be a whole number of 0 or more.", arg = "k")
    }
    N <- .df_recycle(as.numeric(N), n, "N")
    k <- .df_recycle(as.numeric(k), n, "k")
    if (any(N <= k)) {
      .dmar_refuse_argument(c(
        "'N' must exceed 'k': with N <= k even complete data leave no ",
        "degrees of freedom."), arg = c("N", "k"))
    }
    df_raw <- .df_effective_n(gamma, N, k)
    nu_out <- rep(NA_real_, n)
    k_out  <- k
    N_out  <- N
  } else {
    if (!is.null(k)) {
      .dmar_refuse_argument(c(
        "'k' is used by method = \"effective_n\"; ",
        "method = \"complete_df\" takes 'df_complete'."), arg = c("k", "df_complete"))
    }
    if (!is.null(N)) {
      .dmar_refuse_argument(c(
        "'N' is used by method = \"effective_n\"; ",
        "method = \"complete_df\" takes 'df_complete'."), arg = c("N", "df_complete"))
    }
    if (is.null(df_complete)) {
      .dmar_refuse_argument(c(
        "method = \"complete_df\" needs 'df_complete', the degrees of ",
        "freedom the interval would use with complete data (the error ",
        "degrees of freedom, not the model degrees of freedom of a ",
        "structural equation model). It has no default, because the ",
        "right value depends on the parameter: n - k for a regression ",
        "coefficient, n - 1 for a mean."), arg = "df_complete")
    }
    if (!is.numeric(df_complete) || anyNA(df_complete) || any(df_complete <= 0)) {
      .dmar_refuse_argument(
        "'df_complete' must be positive (Inf is allowed).", arg = "df_complete")
    }
    df_complete <- .df_recycle(as.numeric(df_complete), n, "df_complete")
    df_raw <- .df_complete(gamma, df_complete)
    nu_out <- df_complete
    k_out  <- rep(NA_real_, n)
    N_out  <- if (is.null(N_from_table)) rep(NA_real_, n) else N_from_table
  }

  bounded <- .df_bound(df_raw, floor)
  out <- data.frame(term         = term,
                    df           = bounded$df,
                    df_unbounded = df_raw,
                    fmi          = gamma,
                    N            = N_out,
                    k            = k_out,
                    df_complete  = nu_out,
                    at_floor     = bounded$at_floor,
                    stringsAsFactors = FALSE)
  attr(out, "method") <- method
  attr(out, "floor") <- floor
  .as_dmar_tbl(out, subclass = "dmar_df_missing")
}

# The rule behind the df column and its floor, in one paragraph under the
# table, at the console and in a knitted report alike (see
# .dmar_tbl_footer() in R/dmar_output_helpers.R). A column subset, which
# loses the method and floor attributes, prints no line.
#' @exportS3Method .dmar_tbl_footer dmar_df_missing
#' @noRd
.dmar_tbl_footer.dmar_df_missing <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .df_missing_footer(x), "df_rule")
}

# The line itself: the rule by name and formula, broken for the console
# (a knitted report reads the break as a space), and the floor.
.df_missing_footer <- function(x) {
  method <- attr(x, "method", exact = TRUE)
  floor <- attr(x, "floor", exact = TRUE)
  if (!is.character(method) || length(method) != 1L || is.null(floor)) {
    return(NULL)
  }
  rule <- switch(method,
    effective_n = paste0("von Hippel's (2016, Eq. 31) effective sample ",
                         "size rule,\nN(1 - fmi) - k"),
    complete_df = paste0("von Hippel's (2016, Eq. 32) complete data ",
                         "rule,\ndf_complete (1 - fmi) (df_complete + 1) / ",
                         "(df_complete + 3)"),
    NULL)
  if (is.null(rule)) return(NULL)
  sprintf("Degrees of freedom: %s, bounded below at %s", rule,
          format(floor, digits = 6, drop0trailing = TRUE, trim = TRUE))
}


# The tidy() and glance() methods for df_missing(). They are documented
# with the other DMAR tidiers at ?dmar_tidiers (R/dmar_tidiers.R). The
# bespoke pair exists because the logical at_floor column would otherwise
# be read as a row label by the generic wide tidier.

#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_df_missing <- function(x, ...) {
  .df_tidy(x)
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_df_missing <- function(x, ...) {
  data.frame(
    n_terms    = nrow(x),
    n_at_floor = sum(x$at_floor, na.rm = TRUE),
    method     = attr(x, "method") %||% NA_character_,
    floor      = attr(x, "floor") %||% NA_real_,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}
