#' Structural Relations Estimated After the Measurement Models
#'
#' Fits local structural-after-measurement (SAM) models using
#' \code{lavaan::sam()}. Each measurement block is estimated separately,
#' and the structural relations use a covariance matrix corrected for
#' measurement error. Two-step standard errors account for estimating the
#' measurement parameters. This offers an alternative to treating estimated
#' factor scores as perfectly measured observations. A latent SEM remains
#' available when this wrapper's narrower scope does not fit the study.
#'
#' @param model A single-group lavaan model string with at least two
#'   simple-structure factors, at least three distinct continuous indicators
#'   per factor, and structural regressions among those factors. Indicators
#'   cannot cross-load. Derived parameters, explicit equality constraints,
#'   observed covariates, residual covariances, feedback paths, and latent
#'   interactions are outside this interface's validated scope.
#' @param data Raw numeric observations on the model's indicators.
#'   Incomplete observations are refused; the function never drops cases
#'   automatically. Use a joint latent SEM with FIML for incomplete data.
#'   An indicator of \code{model} that \code{data} lacks is refused with an
#'   error of class \code{dmar_items_invalid} (see
#'   \code{\link{dmar_metadata}}), as are columns of \code{data} that
#'   share a name.
#' @param estimator \code{"ML"} (default) with normal-theory
#'   \code{"twostep"} errors, or \code{"MLR"} with
#'   \code{"twostep.robust"} errors. Both use local ML mapping and the
#'   ML measurement estimators; the covariance assumptions for inference
#'   differ. Numerical optimization versus an analytic solution does not
#'   determine which SE is appropriate.
#' @param conf_level Confidence level, default 0.95.
#' @param gamma_unbiased Logical, default \code{TRUE}. For robust inference,
#'   use lavaan's finite-sample unbiased moment-covariance correction.
#'   \code{FALSE} uses the empirical moment covariance (HC0 convention).
#'
#' @details
#' This interface currently supports unweighted, single-level, complete
#' continuous data. It verifies the actual local mapping and two-step SE
#' method after fitting and refuses a substitution. These restrictions
#' distinguish numerical agreement with independent complete-data oracles
#' from evidence for other models or missing-data uncertainty calculations.
#' Small or moderate samples can benefit from separating measurement and
#' structural estimation (Rosseel & Loh, 2024), but SAM does not repair a
#' misspecified or unidentified measurement model.
#'
#' Returned intervals concern individual, original-unit parameters.
#' The underlying fit's joint covariance matrix does not populate every
#' cross-stage covariance. Consequently this function refuses user-defined
#' parameters and does not report standardized structural intervals or
#' intervals for effects combining measurement and structural estimates.
#' The attached fit is retained for inspection with that limitation, not as
#' a complete covariance matrix for arbitrary delta-method calculations.
#'
#' @return A \code{dmar_tbl} with parameter labels, \code{estimate},
#'   \code{se}, \code{statistic}, \code{p_value}, confidence limits,
#'   \code{conf_level}, and \code{stage}. Attributes retain the
#'   \code{fit}, requested and actual \code{method}, and the covariance
#'   limitation. \code{tidy()} retains the rows and \code{glance()} gives
#'   the method and sample size.
#' @references
#' Rosseel, Y., & Loh, W. W. (2024). A structural after measurement (SAM)
#'   approach to structural equation modeling. \emph{Psychological
#'   Methods, 29}(3), 561--588. \doi{10.1037/met0000503}
#' @author Ken Kelley \email{kkelley@@nd.edu}
#' @examples
#' model <- "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
#'           verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
#'           verbal ~ visual"
#' sem_sam(model, holzinger_swineford)
#' @export
sem_sam <- function(model, data, estimator = "ML", conf_level = 0.95,
                    gamma_unbiased = TRUE) {
  .fa_require_lavaan()
  .fa_level(conf_level)
  estimator <- match.arg(estimator, c("ML", "MLR"))
  if (!is.logical(gamma_unbiased) || length(gamma_unbiased) != 1L ||
      is.na(gamma_unbiased)) .dmar_refuse_argument(
        "'gamma_unbiased' must be TRUE or FALSE.", arg = "gamma_unbiased")
  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(model)) {
    .dmar_refuse_argument(
      "'model' must be one lavaan model string.", arg = "model")
  }
  pt <- lavaan::lavaanify(model)
  bad <- !pt$op %in% c("=~", "~", "~~", "~1") |
    grepl(":|\\|", pt$lhs) | grepl(":|\\|", pt$rhs)
  if (any(bad) || any(pt$group > 1L) || any(pt$block > 1L)) {
    .dmar_refuse_argument(c(
      "This local SAM interface supports one continuous group without ",
      "constraints, derived parameters, or interactions."), arg = NULL)
  }
  mm <- pt[pt$op == "=~", , drop = FALSE]
  factors <- unique(mm$lhs)
  indicators <- mm$rhs
  if (length(factors) < 2L || any(table(mm$lhs) < 3L) ||
      anyDuplicated(indicators) || any(indicators %in% factors)) {
    .dmar_refuse_argument(c(
      "Local SAM requires at least two simple-structure factors with ",
      "at least three distinct indicators each and no cross-loadings."),
      arg = NULL)
  }
  paths <- pt[pt$op == "~", , drop = FALSE]
  if (!nrow(paths) || any(!paths$lhs %in% factors | !paths$rhs %in% factors) ||
      any(pt$op == "~~" & pt$lhs != pt$rhs &
            (!pt$lhs %in% factors | !pt$rhs %in% factors))) {
    .dmar_refuse_argument(c(
      "Structural regressions must link the latent factors only. ",
      "Observed covariates and residual covariances are outside this ",
      "local SAM interface's scope."), arg = NULL)
  }
  # A transitive closure detects feedback without depending on factor order.
  reach <- matrix(FALSE, length(factors), length(factors))
  reach[cbind(match(paths$rhs, factors), match(paths$lhs, factors))] <- TRUE
  for (k in seq_along(factors)) reach <- reach | outer(reach[, k], reach[k, ], "&")
  if (any(diag(reach))) .dmar_refuse_argument(
    "Feedback paths are outside this local SAM interface's scope.", arg = NULL)
  labelled <- pt$label[nzchar(pt$label)]
  if (anyDuplicated(labelled)) {
    .dmar_refuse_argument(c(
      "Repeated parameter labels impose constraints outside this local ",
      "SAM interface's scope."), arg = NULL)
  }
  # An indicator the data lack is named here, before .fa_data(), whose
  # message speaks of an 'items' argument this function does not have;
  # the class and reason are the ones .fa_data() gives it
  # (dmar_items_invalid, "unknown"; .fa_items_absent()). The names are
  # read as .fa_data() reads them, from the data frame the input becomes.
  if (is.data.frame(data) || is.matrix(data)) {
    absent <- setdiff(indicators, names(as.data.frame(data)))
    if (length(absent)) {
      .fa_items_absent(c("Indicators not found in 'data': ",
                         paste(absent, collapse = ", "), "."), absent)
    }
  }
  data <- .fa_data(data, indicators)
  if (anyNA(data)) {
    .dmar_refuse_argument(c(
      "Local SAM currently requires complete indicator data. No cases ",
      "were deleted. Use a joint latent SEM with FIML for incomplete ",
      "data; the missing-data SAM uncertainty route needs independent ",
      "validation before this interface can report it."), arg = NULL)
  }
  se <- if (estimator == "ML") "twostep" else "twostep.robust"
  fit <- lavaan::sam(model, data = data, estimator = estimator,
    missing = "listwise", sam_method = "local", se = se,
    gamma.unbiased = gamma_unbiased,
    mm_args = list(estimator = estimator, bounds = "none"),
    struc_args = list(estimator = "ML"),
    local_options = list(m_method = "ML", lambda_correction = TRUE,
      lambda1_floor = "debias", lambda2_floor = "debias", alpha_correction = 0L))
  actual <- fit@internal
  if (!identical(actual$sam.method, "local") ||
      !identical(actual$sam.local.options$m_method, "ML") ||
      !identical(fit@Options$se, se) || !identical(fit@vcov$se, se) ||
      !identical(actual$sam.lavoptions$estimator.orig, estimator)) {
    .dmar_refuse_argument(c(
      "lavaan changed the requested local SAM mapping or uncertainty ",
      "method. No substituted analysis is returned."), arg = NULL)
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged")) ||
      !isTRUE(lavaan::lavInspect(fit, "post.check")) ||
      any(unlist(actual$sam_trunc) != 0)) {
    .dmar_refuse_argument(c(
      "The SAM model did not yield an untruncated, proper solution. ",
      "Review the measurement blocks before interpreting structural ",
      "relations."), arg = NULL)
  }
  pe <- lavaan::parameterEstimates(fit, level = conf_level, remove_step1 = FALSE)
  out <- .fa_parameter_table(pe, conf_level)
  out$stage <- ifelse(pe$op == "=~" | pe$lhs %in% indicators,
                      "measurement", "structural")
  free <- !is.na(pe$se) & pe$se > 0
  if (!any(free) || any(!is.finite(pe$se))) {
    .dmar_refuse_argument(
      "Two-step parameter uncertainty is unavailable for this SAM fit.",
      arg = NULL)
  }
  limitation <- paste0("Local SAM: ", estimator, "; SE: ", se,
    "; complete continuous data; gamma_unbiased = ", gamma_unbiased,
    ". Intervals concern individual original-unit parameters. Cross-stage ",
    "covariance cells are incomplete; do not use the attached fit for ",
    "arbitrary derived or standardized structural intervals.")
  attr(out, "fit") <- fit
  attr(out, "method") <- list(N = nrow(data), estimator = estimator,
    se = se, missing = "complete data", sam_method = actual$sam.method,
    mapping = actual$sam.local.options$m_method, gamma_unbiased = gamma_unbiased,
    measurement_args = actual$sam.mm.args, structural_args = actual$sam.struc.args)
  attr(out, "analysis_note") <- limitation
  .as_dmar_tbl(out, conf_level = conf_level, p_terms = "p_value",
                subclass = "dmar_factor_analysis")
}
