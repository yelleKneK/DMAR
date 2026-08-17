# ANCOVA: adjusted means, covariate-adjusted omnibus F, effect size CIs.
#' Analysis of Covariance (ANCOVA)
#'
#' Fits a one-way analysis of covariance so the covariate-adjusted group
#' comparison is available from a single call, without assembling the
#' adjusted means, the omnibus test, and the effect sizes by hand. It
#' returns the adjusted (covariate-corrected) group means with standard
#' errors, the covariate-adjusted omnibus \emph{F} for the group effect
#' (Type III sums of squares), partial \eqn{\eta^2} and partial
#' \eqn{\omega^2} with noncentral \emph{F} confidence intervals, and a
#' homogeneity-of-regression check, returned in one \code{data.frame}.
#'
#' @param data A \code{data.frame} containing the response, the
#'   treatment factor, and the covariate(s).
#' @param outcome Character name of the response column.
#' @param treatment Character name of the grouping factor column (the
#'   groups being compared, for example treatment arms); a factor or
#'   character column.
#' @param covariates Character vector of one or more covariate
#'   column names.
#' @param conf_level Confidence level for the effect size CIs.
#'   Default \code{0.95}.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows for the
#'   omnibus test (\code{F_value}, \code{df_1}, \code{df_2},
#'   \code{p_value}), the sum-of-squares type used
#'   (\code{sum_of_squares_type}; 3 for Type III), the point estimates and
#'   confidence intervals of partial \eqn{\eta^2} and partial
#'   \eqn{\omega^2}, the adjusted group means and their standard errors
#'   (one row per level of \code{treatment}), and the
#'   homogeneity-of-regression \emph{F}-test. The result carries the
#'   \code{dmar_tbl} class, so it
#'   prints to 3 significant figures with whole numbers (such as the
#'   degrees of freedom) shown without a decimal part and \emph{p}-values
#'   to 4 decimal places (a \emph{p}-value below 0.0001 prints as
#'   \dQuote{< 0.0001}); the stored values keep full precision. Control
#'   the display with \code{print(x, digits = )} or globally with
#'   \code{options(dmar.digits = )} (see \code{\link{dmar_tbl}}).
#'
#' @details
#' \strong{Covariate-adjusted test (Type III sums of squares).} The omnibus
#' \emph{F} tests the group effect after adjusting for the covariate(s),
#' that is, the Type III sum of squares for the grouping factor. For a
#' one-way ANCOVA (one grouping factor, with the covariate slopes held
#' constant) the Type II and Type III sums of squares for the group effect
#' coincide, and both equal the sequential sum of squares obtained with the
#' covariate(s) entered first and the grouping factor last, which is how it
#' is computed here; the value matches \code{car::Anova(fit, type = 3)}. The
#' choice of sum-of-squares type changes the result only in designs with
#' more than one factor or with interactions among factors (Maxwell,
#' Delaney, and Kelley, 2027, Chapter 7); for those, the two-way and mixed
#' analyses report their sum-of-squares type and allow Type I, II, or III.
#'
#' \strong{Adjusted means.} The adjusted mean for treatment level
#' \eqn{j} is the model-predicted response at \eqn{X = \bar X} (the
#' covariate grand mean):
#' \deqn{\hat \mu_j^{\mathrm{adj}} \;=\; \hat\mu_j -
#'   \sum_k \hat\beta_k (\bar X_{kj} - \bar X_k),}
#' where \eqn{\hat\beta_k} is the within-cell slope on covariate \eqn{k}
#' and \eqn{\bar X_{kj}, \bar X_k} are the per-cell and grand means of
#' covariate \eqn{k}.
#'
#' \strong{Homogeneity of regression.} The model fit here holds the
#' within-group covariate slopes \eqn{\beta_k} constant across groups. This
#' is a property of the particular model being fit, not an assumption of
#' analysis of covariance in general: it is a testable claim. Adding all
#' group-by-covariate interactions gives an expanded model, and a model
#' comparison \emph{F}-test of the additive model against the expanded one
#' (\code{stats::anova}) assesses whether the slopes differ across groups
#' (Maxwell, Delaney, and Kelley, 2027, Chapter 9). A large \emph{F}
#' indicates the slopes are not constant, in which case the single adjusted
#' comparison is not the whole story and the interaction model should be
#' entertained directly. The check is reported in the
#' \code{F_homogeneity_of_regression} rows.
#'
#' \strong{Effect size CIs.} Partial \eqn{\eta^2} and partial
#' \eqn{\omega^2} use the noncentral \emph{F} framework
#' (\code{\link{ci_eta_squared_partial}},
#' \code{\link{ci_omega_squared}}).
#'
#' @references
#' Huitema, B. E. (2011). \emph{The analysis of covariance and
#'   alternatives} (2nd ed.). Wiley.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 9 on analysis of
#'   covariance and Chapter 7 on higher-order designs.)
#'
#' @seealso \code{\link{ci_c_ancova}}, \code{\link{ci_sc_ancova}},
#'   \code{\link{ss_aipe_sc_ancova}}, \code{\link{omega_squared_partial}}
#'
#' @examples
#' # 1. Compare the two groups in the Pygmalion data on eighth-grade IQ,
#' #    adjusting for pre-test IQ.
#' ancova(outcome = "iq_8", treatment = "treatment",
#'        covariates = "iq_pre", data = pygmalion)
#'
#' # 2. The same comparison adjusting for two covariates (pre-test IQ and
#' #    grade).
#' ancova(outcome = "iq_8", treatment = "treatment",
#'        covariates = c("iq_pre", "grade"), data = pygmalion)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest design
#'
#' @family hypothesis tests
#'
#' @export

ancova <- function(data, outcome, treatment, covariates,
                   conf_level = 0.95) {
  if (!is.character(outcome) || !outcome %in% names(data))
    stop("'outcome' must be a column name in 'data'.")
  if (!is.character(treatment) || !treatment %in% names(data))
    stop("'treatment' must be a column name in 'data'.")
  if (!is.character(covariates) || length(covariates) < 1L ||
      !all(covariates %in% names(data)))
    stop("'covariates' must be one or more column names in 'data'.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  d <- data[, c(outcome, treatment, covariates), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d[[treatment]] <- factor(d[[treatment]])
  for (cv in covariates)
    if (!is.numeric(d[[cv]]))
      stop(sprintf("Covariate '%s' must be numeric.", cv))

  # Fit ANCOVA and the full model (covariate x treatment) for HoR.
  rhs_cov <- paste(covariates, collapse = " + ")
  f_ancova <- stats::as.formula(sprintf("%s ~ %s + %s",
                                        outcome, treatment, rhs_cov))
  f_full   <- stats::as.formula(sprintf("%s ~ %s * (%s)",
                                        outcome, treatment, rhs_cov))
  fit_ancova <- stats::lm(f_ancova, data = d)
  fit_full   <- stats::lm(f_full,   data = d)

  # Omnibus F for the covariate-adjusted treatment effect. With the
  # covariate(s) entered first and the grouping factor last, the factor's
  # sequential sum of squares is its covariate-adjusted sum of squares,
  # which for this one-way ANCOVA equals the Type III (and Type II) sum of
  # squares for the factor and matches car::Anova(fit, type = 3). fit_ancova
  # is retained for the adjusted means and the homogeneity-of-regression
  # comparison.
  fit_for_F <- stats::lm(
    stats::as.formula(sprintf("%s ~ %s + %s", outcome, rhs_cov, treatment)),
    data = d)
  # A covariate collinear with the treatment factor (e.g., constant within
  # groups) makes the design rank deficient; lm() then silently aliases a
  # treatment dummy (its coefficient is NA) and the sequential table reports a
  # reduced-df remainder as the "treatment" test. Refuse rather than return a
  # mislabeled effect.
  if (anyNA(stats::coef(fit_for_F))) {
    stop("The covariate(s) are collinear with the treatment factor, so the ",
         "covariate-adjusted treatment effect is not estimable (the ANCOVA ",
         "design is rank deficient). Check for a covariate that is constant ",
         "within groups or otherwise aliased with the treatment.", call. = FALSE)
  }
  tab <- stats::anova(fit_for_F)
  rn  <- trimws(rownames(tab))
  ti  <- which(rn == treatment)
  ri  <- which(rn == "Residuals")
  F_treat <- tab[ti, "F value"]
  df_1    <- tab[ti, "Df"]
  df_2    <- tab[ri, "Df"]
  p_treat <- tab[ti, "Pr(>F)"]

  N <- stats::nobs(fit_ancova)

  # Effect size CIs:
  eta_p   <- ci_eta_squared_partial(F_value = F_treat, df_effect = df_1,
                                    df_error = df_2, N = N,
                                    conf_level = conf_level)
  omega_p <- ci_omega_squared(F_value = F_treat, df_effect = df_1,
                              df_error = df_2, N = N,
                              conf_level = conf_level)

  # Adjusted means: predict response at covariate grand means.
  newdata <- expand.grid(stats::setNames(
    c(list(levels(d[[treatment]])),
      lapply(covariates, function(cv) mean(d[[cv]]))),
    c(treatment, covariates)
  ))
  pred <- stats::predict(fit_ancova, newdata = newdata, se.fit = TRUE)
  adj_mean_df <- data.frame(
    term  = paste0("adjusted_mean[", as.character(newdata[[treatment]]), "]"),
    value = pred$fit,
    stringsAsFactors = FALSE, row.names = NULL
  )
  adj_se_df <- data.frame(
    term  = paste0("se_adjusted_mean[", as.character(newdata[[treatment]]), "]"),
    value = pred$se.fit,
    stringsAsFactors = FALSE, row.names = NULL
  )

  # Homogeneity-of-regression test:
  hor <- stats::anova(fit_ancova, fit_full)
  F_hor  <- hor[["F"]][2]
  df_hor <- hor[["Df"]][2]
  p_hor  <- hor[["Pr(>F)"]][2]

  out <- rbind(
    data.frame(
      term  = c("F_value", "df_1", "df_2", "p_value",
                "sum_of_squares_type",
                "eta_squared_partial",
                "eta_squared_partial_lower",
                "eta_squared_partial_upper",
                "omega_squared_partial",
                "omega_squared_partial_lower",
                "omega_squared_partial_upper"),
      value = c(F_treat, df_1, df_2, p_treat,
                3,
                eta_p$eta_squared_partial,
                eta_p$lower_limit,
                eta_p$upper_limit,
                omega_p$omega_squared,
                omega_p$lower_limit,
                omega_p$upper_limit),
      stringsAsFactors = FALSE, row.names = NULL
    ),
    adj_mean_df,
    adj_se_df,
    data.frame(
      term  = c("F_homogeneity_of_regression", "df_homogeneity_of_regression",
                "p_homogeneity_of_regression"),
      value = c(F_hor, df_hor, p_hor),
      stringsAsFactors = FALSE, row.names = NULL
    )
  )

  # The effect size rows (eta_squared_partial / omega_squared_partial and
  # their limits) are reported at conf_level; recording it lets the
  # dmar_tbl print method note the confidence level beneath the table.
  attr(out, "conf_level") <- conf_level
  # Name the p-value rows so dmar_tbl prints them to fixed decimal places
  # (4 by default) rather than significant figures; this is a long-format
  # table, so the rows are flagged by term name, not a dedicated column.
  attr(out, "p_terms") <- c("p_value", "p_homogeneity_of_regression")
  class(out) <- c("dmar_tbl", "data.frame")
  out
}
