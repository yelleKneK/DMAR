#' Reliability Coefficient With a Confidence Interval (General Dispatch)
#'
#' @description
#' General-purpose entry point for the reliability family. Dispatches to
#' \code{\link{reliability_alpha}},
#' \code{\link{reliability_kr20}},
#' \code{\link{reliability_omega}},
#' \code{\link{reliability_omega_categorical}}, or
#' \code{\link{reliability_H}} according to the requested
#' \code{type}. When \code{type} is not specified, the function picks a
#' reasonable default from the supplied input following the
#' recommendations of Kelley and Pornprasertmanit (2016).
#'
#' @details
#' Auto-detection rules (used only when \code{type = NULL}):
#' \itemize{
#'   \item If raw items are integer-valued and every item has at most
#'         10 distinct values, \code{type = "omega_categorical"} (categorical
#'         omega; appropriate when items are ordered categorical and
#'         the relationship between the underlying factor and the
#'         observed items is non-linear).
#'   \item Otherwise, \code{type = "omega"} (McDonald's coefficient
#'         \eqn{\omega} from a single-factor CFA). This includes the
#'         case where only a covariance matrix \code{S} and sample
#'         size \code{N} are supplied; \pkg{lavaan} fits the CFA on
#'         the covariance matrix.
#' }
#' The rule reads the items only. The columns of \code{data} named in
#' \code{aux} are auxiliary variables, which never enter the composite,
#' so their scale says nothing about the items' response format and they
#' take no part.
#' Auto-detection emits a single \code{message()} indicating which
#' \code{type} was chosen, so the choice is stated rather than made
#' silently. It
#' never chooses \code{type = "H"}, which describes a different composite
#' (the optimally weighted one) rather than a different set of
#' assumptions about the summed one, and is computed only when asked for.
#'
#' The selected family function determines which \code{ci_method} values
#' are accepted; see the help page for the chosen function for the full
#' list. When \code{ci_method} is left at its default (\code{NULL}), the
#' family function's own default is used:
#' \itemize{
#'   \item \code{reliability_alpha}: \code{"bonett"}, or \code{"ml"} when
#'         incomplete items are analyzed by full information maximum
#'         likelihood (the default of \code{missing} for them).
#'   \item \code{reliability_alpha(estimator = "model_implied")}: \code{"mlr"} with raw
#'         data; \code{"ml"} with covariance input.
#'   \item \code{reliability_kr20}: \code{"feldt"}.
#'   \item \code{reliability_omega}: for
#'         \code{denominator = "observed"} (robust omega, the default),
#'         the point estimate with no interval, since its interval is
#'         bootstrap based and no bootstrap runs unless requested; for
#'         \code{denominator = "model_implied"}, \code{"mlr_logistic"}
#'         with raw data and \code{"ml_logistic"} with covariance input,
#'         the limits formed on the logit scale (\code{"mlr"} and
#'         \code{"ml"} under \code{wald_df = "missing_information"}).
#'   \item \code{reliability_omega_categorical}: the point estimate with
#'         no interval, for the same reason; request \code{"bca"}.
#'   \item \code{reliability_H}: \code{"mlr_logistic"} with raw data;
#'         \code{"ml_logistic"} with covariance input.
#' }
#' A bootstrap is never run by default anywhere in the family. When a
#' bootstrap method is requested, \code{B = 10000} replications is the
#' default, except for categorical omega under \code{missing = "mi"},
#' whose bootstrap then impute replications each impute the data and
#' default to \code{B = 1000}.
#'
#' Several reliability coefficients exist because their assumptions
#' differ. Coefficient \eqn{\alpha} (and its dichotomous specialization
#' KR-20) equals the population reliability under essential
#' \eqn{\tau}-equivalence (equal loadings); McDonald's \eqn{\omega}
#' relaxes that assumption to a congeneric single-factor model; and
#' \code{reliability_omega(denominator = "observed")} further relaxes
#' the requirement that the single-factor model fit perfectly by using
#' the observed composite variance in the denominator (see the
#' \code{\link{reliability_omega}} help page for the properties of that
#' choice). For well-behaved homogeneous measurement
#' instruments these coefficients typically yield very similar values.
#' For ordered categorical items the relationship between the latent
#' factor and the observed responses is non-linear, and
#' \code{reliability_omega_categorical} handles that case explicitly via a
#' probit-link single-factor model.
#'
#' @param data A numeric matrix or data frame of item scores, or
#'   \code{NULL}.
#' @param S A symmetric covariance matrix among the items, or
#'   \code{NULL}. If supplied, \code{N} must also be supplied; methods
#'   that require raw data are then unavailable.
#' @param N Total sample size; required when \code{S} is supplied.
#' @param type Character; one of \code{"alpha"},
#'   \code{"kr20"},
#'   \code{"omega"}, \code{"omega_categorical"} (\code{"omega_c"} is
#'   accepted as a shorthand), \code{"H"} (coefficient \emph{H}, the
#'   reliability of the optimally weighted composite, from raw data or a
#'   covariance matrix; \code{\link{reliability_H}} also takes a fitted
#'   model or standardized loadings), or \code{NULL}
#'   for auto-detection. See \emph{Details}.
#' @param estimator For \code{type = "alpha"} only: how coefficient
#'   \eqn{\alpha} is estimated, \code{"analytic"} (default; the
#'   closed-form equation applied to the observed covariance matrix) or
#'   \code{"model_implied"} (the reliability implied by the
#'   \eqn{\tau}-equivalent single-factor model fit by maximum
#'   likelihood); forwarded to \code{\link{reliability_alpha}}, whose
#'   help page discusses the choice and which interval methods each
#'   estimator supports. Supplying it with any other \code{type} is an
#'   error, as supplying \code{denominator} is for any type but
#'   \code{"omega"} and \code{"omega_categorical"}.
#' @param denominator For \code{type = "omega"} and
#'   \code{type = "omega_categorical"} only: how the total variance in
#'   the denominator is estimated, \code{"observed"} (default; robust
#'   omega, and for categorical omega the unrestricted polychoric
#'   correlations) or \code{"model_implied"} (for categorical omega, the
#'   coefficient as Green and Yang, 2009, defined it); forwarded to
#'   \code{\link{reliability_omega}} or
#'   \code{\link{reliability_omega_categorical}}, whose help pages
#'   discuss the choice.
#' @param missing How incomplete rows of \code{data} are handled,
#'   forwarded to the family function when supplied, which checks it
#'   against its own choices and whose help page discusses them: for
#'   \code{type = "alpha"} and \code{type = "omega"}, \code{"listwise"},
#'   \code{"fiml"} (full information maximum likelihood), or
#'   \code{"two_stage"} (two-stage maximum likelihood); for
#'   \code{type = "H"}, \code{"listwise"} or \code{"fiml"}; for
#'   \code{type = "omega_categorical"}, \code{"listwise"},
#'   \code{"pairwise"}, or \code{"mi"} (multiple imputation). The
#'   default, \code{NULL}, leaves the family function's default in place:
#'   for alpha, omega, and H, full information maximum likelihood when an
#'   item value is missing and the complete data analysis when none is,
#'   since the complete-data intervals of alpha and omega need that
#'   analysis (\code{\link{fa_confirmatory}} fits complete data by FIML
#'   too; a row with every item missing does not count, since both
#'   treatments drop it), and \code{"listwise"} for categorical omega. Supplying it with
#'   \code{type = "kr20"} is an error.
#' @param aux Optional character vector naming auxiliary variable
#'   columns of \code{data}, forwarded to the family function: saturated
#'   correlates under full information maximum likelihood, and variables
#'   entered into the first stage under two-stage estimation, for
#'   \code{type = "alpha"}, \code{type = "omega"}, and \code{type =
#'   "H"} (supplying \code{aux} without \code{missing} implies
#'   \code{missing = "fiml"}); predictors in the imputation model for
#'   \code{type = "omega_categorical"} (it implies \code{missing =
#'   "mi"}). Supplying it with \code{type = "kr20"} is an error.
#' @param ci_method Method for constructing the confidence interval, or
#'   \code{NULL} to use the chosen family function's default. See the
#'   help page for the chosen \code{reliability_*} function for the full
#'   list of accepted values.
#' @param conf_level Confidence level. Defaults to \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected. The default, \code{NULL}, uses the family function's
#'   default: 10000, or 1000 for categorical omega under
#'   \code{missing = "mi"}. The family function checks \code{B} against the number
#'   its interval method and level require (400 for \code{"percentile"},
#'   2000 for \code{"bca"}, 100 for the two bootstrap standard error
#'   methods, at the default \code{conf_level}) and refuses a smaller one
#'   (class \code{dmar_count_invalid}); see
#'   \code{\link{dmar_bootstrap}}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, under which the resamples are drawn from,
#'   and advance, the current random number stream; supply an integer
#'   for a reproducible interval, which is set for the call only, the
#'   caller's random number state being restored when it returns.
#' @inheritParams reliability_omega
#' @param ncpus The number of processes over which the bootstrap
#'   replications are spread; forwarded to every family function. The
#'   default \code{1} runs them in this R session, and a given
#'   \code{seed} returns the same interval whatever \code{ncpus} is (see
#'   \code{\link{reliability_omega}}).
#' @param improper For \code{type = "alpha"} and \code{type = "omega"}
#'   only: how a bootstrap replication whose single-factor solution is
#'   improper but still yields a coefficient is treated, \code{"keep"}
#'   (the default) or \code{"drop"}; forwarded to
#'   \code{\link{reliability_alpha}} or \code{\link{reliability_omega}},
#'   whose help pages describe it. Categorical omega is not defined for
#'   an improper solution, and KR-20 fits no model, so supplying it with
#'   another \code{type} is an error.
#' @param impute,imputations For \code{type = "omega_categorical"} with
#'   \code{missing = "mi"} only: the imputation model and the number of
#'   imputations of each bootstrap replication, forwarded to
#'   \code{\link{reliability_omega_categorical}}, whose help page
#'   describes them. The default, \code{NULL}, leaves that function's
#'   defaults (\code{"pmm"} and 5) in place; supplying either with
#'   another \code{type} is an error.
#' @param wald_df For \code{type = "alpha"} and \code{type = "omega"}
#'   only: the reference distribution of a Wald interval, \code{Inf} (the
#'   default, the normal reference), a single positive number of degrees
#'   of freedom of a \emph{t} reference, or, for \code{type = "omega"}
#'   with \code{denominator = "model_implied"} only,
#'   \code{"missing_information"}, the degrees of freedom of a \emph{t}
#'   reference from omega's fraction of missing information; forwarded to
#'   \code{\link{reliability_alpha}} or \code{\link{reliability_omega}},
#'   whose help pages describe it (the keyword is an error for alpha,
#'   for which no rule passed the calibration study's decision rule). A
#'   finite value or the keyword with another \code{type} is an error.
#'
#' @return The \code{data.frame} returned by the dispatched
#' \code{reliability_*} function (rows: \code{estimate}, \code{se},
#' \code{ci_lower}, \code{ci_upper}, \code{conf_level}, \code{N},
#' \code{N_complete}, \code{n_items}). The \code{se} row is on the coefficient
#' scale; the transformation-based intervals (\code{"fisher"},
#' \code{"bonett"}, \code{"hakstian_whalen"}) add an
#' \code{se_transformed} row carrying the transformation-scale standard
#' error, with the scale named in the \code{se_transform_scale}
#' attribute. The \code{coefficient} attribute
#' identifies which coefficient was computed, coefficients alpha and omega
#' record \code{wald_df}, and a bootstrap call also
#' carries the \code{bootstrap} attribute with its replication counts
#' (see \code{\link{reliability_alpha}}).
#'
#' @references
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores
#'   using structural equation modeling: An alternative to coefficient
#'   alpha. \emph{Psychometrika, 74}, 155--167.
#'   \doi{10.1007/s11336-008-9099-3}
#'
#' Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
#'   formation for reliability coefficients of homogeneous measurement
#'   instruments. \emph{Methodology, 8}, 39--50.
#'   \doi{10.1027/1614-2241/a000036}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{reliability_alpha}}, \code{\link{reliability_kr20}},
#' \code{\link{reliability_omega}},
#' \code{\link{reliability_omega_categorical}},
#' \code{\link{ss_aipe_reliability}}, \code{\link{cfa_1}}.
#'
#' @examples
#' set.seed(113)
#' J <- 6
#' loadings <- seq(0.4, 0.8, length.out = J)
#' eta <- rnorm(200)
#' errors <- matrix(rnorm(200 * J), 200, J) %*% diag(sqrt(1 - loadings^2))
#' items <- sweep(matrix(rep(eta, J), 200, J), 2, loadings, `*`) + errors
#' colnames(items) <- paste0("y", seq_len(J))
#'
#' # Auto-detection picks coefficient omega for continuous data.
#' reliability(data = items)
#'
#' # Explicit type.
#' reliability(data = items, type = "alpha")
#'
#' # A covariance matrix and its sample size stand in for raw data, and
#' # auto-detection again picks coefficient omega. The point estimate is
#' # the same as from the raw data; the message explains why no interval
#' # accompanies it, since the bootstrap behind robust omega resamples
#' # rows that a covariance matrix does not carry.
#' reliability(S = cov(items), N = 200)
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability <- function(data = NULL, S = NULL, N = NULL,
                        type = NULL,
                        estimator = c("analytic", "model_implied"),
                        denominator = c("observed", "model_implied"),
                        missing = NULL,
                        aux = NULL,
                        ci_method = NULL,
                        conf_level = 0.95,
                        B = NULL,
                        seed = NULL,
                        ncpus = 1L,
                        improper = c("keep", "drop"),
                        impute = NULL,
                        imputations = NULL,
                        wald_df = Inf,
                        progress = NULL) {
  .dmar_check_conf_level(conf_level)

  # SPSS user-defined missing values are missing, both in the choice of
  # coefficient and in the function called (.relia_declared_missing()).
  if (!is.null(data)) data <- .relia_declared_missing(data)
  type <- .resolve_relia_type(type, data, S, aux)
  est_supplied <- !missing(estimator)
  estimator <- match.arg(estimator)
  if (est_supplied && type != "alpha") {
    .dmar_refuse_argument(
      "'estimator' applies only to type = \"alpha\".", arg = "estimator")
  }
  denom_supplied <- !missing(denominator)
  denominator <- match.arg(denominator)
  if (denom_supplied && !type %in% c("omega", "omega_categorical")) {
    .dmar_refuse_argument(c(
      "'denominator' applies only to type = \"omega\" and ",
      "type = \"omega_categorical\"."), arg = "denominator")
  }
  # `missing` is forwarded as given, and each family function checks it
  # against its own choices (alpha and omega: listwise, fiml, two_stage;
  # H: listwise, fiml; categorical omega: listwise, pairwise, mi). KR-20
  # has no missing data treatment beyond listwise deletion.
  missing_supplied <- !is.null(missing)
  if ((missing_supplied || !is.null(aux)) && type == "kr20") {
    .dmar_refuse_argument(c(
      "'missing' and 'aux' apply only to type = \"alpha\", ",
      "\"omega\", \"omega_categorical\", and \"H\"; KR-20 analyzes the ",
      "complete rows only."), arg = c("missing", "aux"))
  }
  if ((!is.null(impute) || !is.null(imputations)) &&
      type != "omega_categorical") {
    .dmar_refuse_argument(c(
      "'impute' and 'imputations' apply only to type = ",
      "\"omega_categorical\" with missing = \"mi\"."),
      arg = c("impute", "imputations"))
  }

  # When the user did not specify a CI method, the family function's own
  # default applies, with one exception: the model implied alpha estimator
  # defaults to the robust ML CI ("mlr"), which requires raw data. If the
  # wrapper is called with only a covariance matrix, fall back to the
  # maximum likelihood CI ("ml") so the call succeeds without forcing the
  # user to know which methods need raw data. reliability_omega() makes
  # the corresponding choice itself since DMAR 1.1.0: with a covariance
  # matrix its model implied default is "ml_logistic" ("ml" under wald_df =
  # "missing_information"), the counterpart of its raw data default
  # "mlr_logistic" ("ml" through 1.0.x on this path, which the wrapper
  # chose then). (With denominator = "observed", reliability_omega
  # resolves its own default to no interval with a message.)
  if (is.null(ci_method) && is.null(data) &&
      type == "alpha" && estimator == "model_implied") {
    ci_method <- "ml"
  }

  # Build the argument list for the chosen family function. ci_method is
  # forwarded only when the caller specified it (or the wrapper has
  # picked a sensible override above); otherwise the family function's
  # own default applies. `missing` is forwarded only when the caller
  # specified it, so the family function still sees the aux-implies-fiml
  # shorthand (an unsupplied `missing` plus `aux`) exactly as a direct
  # call would.
  # `improper` is forwarded only when the caller supplied it, and only
  # to the two functions that take it: categorical omega and H set every
  # improper replication aside (neither coefficient is defined for one)
  # and KR-20 fits no model.
  improper_supplied <- !missing(improper)
  improper <- match.arg(improper)
  if (improper_supplied && !type %in% c("alpha", "omega")) {
    .dmar_refuse_argument(c(
      "'improper' applies only to type = \"alpha\" and type = ",
      "\"omega\": ", if (type == "kr20") {
        "KR-20 fits no model, so no solution can be improper."
      } else {
        paste0(if (type == "H") "coefficient H" else "categorical omega",
               " is not defined for an improper ",
               "solution, so such a bootstrap replication is always ",
               "set aside.")
      }), arg = "improper")
  }
  # `wald_df` is forwarded only when supplied, and a t reference only to
  # the two functions whose Wald intervals take one. KR-20, categorical
  # omega, and H take no wald_df, and their intervals are their own: by
  # default Feldt's F interval, a BCa bootstrap interval, and a Wald
  # interval on the logit scale with the normal reference.
  wald_df_supplied <- !missing(wald_df)
  own_interval <- switch(type,
    kr20 = "Feldt's F interval by default",
    H = paste0("by default a Wald interval on the logit scale with the ",
               "normal reference"),
    omega_categorical = paste0("by default a bias-corrected and ",
                               "accelerated bootstrap interval"),
    NULL)
  no_wald_df <- paste0("type = \"", type, "\" takes no wald_df, and its ",
                       "interval is its own (", own_interval, ").")
  if (wald_df_supplied && identical(wald_df, "missing_information") &&
      !type %in% c("alpha", "omega")) {
    .dmar_refuse_argument(c(
      "wald_df = \"missing_information\" applies only to type = ",
      "\"omega\" (with denominator = \"model_implied\"), the ",
      "coefficient the calibration study behind it covers; ", no_wald_df),
      arg = NULL)
  }
  if (wald_df_supplied && !.wald_df_is_normal(wald_df) &&
      !type %in% c("alpha", "omega")) {
    .dmar_refuse_argument(c(
      "'wald_df' applies only to type = \"alpha\" and type = ",
      "\"omega\", whose Wald intervals take a t reference ",
      "distribution; ", no_wald_df), arg = "wald_df")
  }
  # B is forwarded only when supplied, so each family function's default
  # applies otherwise (10000, or 1000 for categorical omega under
  # multiple imputation).
  common <- list(conf_level = conf_level, seed = seed, ncpus = ncpus)
  if (!is.null(progress)) common$progress <- .check_progress(progress)
  if (!is.null(B)) common$B <- B
  if (improper_supplied) common$improper <- improper
  if (wald_df_supplied && type %in% c("alpha", "omega")) {
    common$wald_df <- wald_df
  }
  if (!is.null(ci_method)) common$ci_method <- ci_method
  if (missing_supplied) common$missing <- missing
  if (!is.null(aux)) common$aux <- aux

  if (type == "alpha") {
    args <- c(list(data = data, S = S, N = N, estimator = estimator), common)
    return(do.call(reliability_alpha, args))
  }
  if (type == "kr20") {
    if (!is.null(S)) {
      .dmar_refuse_argument(c(
        "type = \"kr20\" requires raw 'data'; a covariance matrix is ",
        "not sufficient."), arg = "data")
    }
    args <- c(list(data = data), common)
    return(do.call(reliability_kr20, args))
  }
  if (type == "omega") {
    args <- c(list(data = data, S = S, N = N,
                   denominator = denominator), common)
    return(do.call(reliability_omega, args))
  }
  if (type == "omega_categorical") {
    if (!is.null(S)) {
      .dmar_refuse_argument(c(
        "type = \"omega_categorical\" requires raw 'data'; a covariance matrix ",
        "is not sufficient."), arg = "data")
    }
    args <- c(list(data = data, denominator = denominator), common)
    if (!is.null(impute)) args$impute <- impute
    if (!is.null(imputations)) args$imputations <- imputations
    return(do.call(reliability_omega_categorical, args))
  }
  if (type == "H") {
    # Coefficient H from raw data or a covariance matrix; its own defaults
    # (the robust ML interval on the logit scale with data, the ML one
    # with a covariance matrix) need no override here.
    args <- c(list(data = data, S = S, N = N), common)
    return(do.call(reliability_H, args))
  }
}

# .resolve_relia_type(type, data, S, aux = NULL)
#
# The coefficient reliability() computes: `type` when it is given, and
# otherwise the one the auto-detection rule of the help page picks. The
# rule reads the items only, the columns of `data` not named in `aux`.
# Through DMAR 1.0.1 it read every column, so a continuous auxiliary
# variable (or an integer one with more than 10 values) beside ordered
# items sent them to continuous omega under FIML rather than categorical
# omega, and a factor auxiliary stopped the call with a base R error.
.resolve_relia_type <- function(type, data, S, aux = NULL) {
  if (!is.null(type)) {
    type <- match.arg(type, c("alpha", "kr20", "omega",
                              "omega_categorical", "omega_c", "H"))
    # "omega_c" is accepted as a shorthand for the full name.
    if (type == "omega_c") type <- "omega_categorical"
    return(type)
  }

  if (is.null(data) && is.null(S)) {
    .dmar_refuse_argument(c(
      "Either 'data' or 'S' must be supplied so reliability() can ",
      "select an appropriate coefficient."), arg = c("data", "S"))
  }

  # When raw items are available, check whether they look ordered-
  # categorical (integer-valued and few distinct values per item). Only
  # in that case do we override the default of coefficient omega; the
  # categorical-omega estimator treats the items via a probit-link
  # single-factor model, which the continuous-item omega cannot do.
  if (!is.null(data)) {
    dat <- as.data.frame(data)
    if (is.character(aux)) {
      dat <- dat[, setdiff(colnames(dat), aux), drop = FALSE]
    }
    dat <- as.matrix(dat)
    is_integer_valued <- ncol(dat) > 0L &&
      isTRUE(all(dat == round(dat), na.rm = TRUE))
    max_unique <- if (ncol(dat) > 0L) {
      max(apply(dat, 2, function(x) length(unique(stats::na.omit(x)))))
    } else {
      Inf
    }
    if (is_integer_valued && max_unique <= 10) {
      message("Auto-detected type = \"omega_categorical\" (integer items ",
              "with <= 10 distinct values; treating items as ordered ",
              "categorical).")
      return("omega_categorical")
    }
  }

  message("Auto-detected type = \"omega\" (McDonald's coefficient omega ",
          "from a single-factor CFA).")
  "omega"
}


#' A Reliability Coefficient Estimate
#'
#' Returns a one-row \code{data.frame} in the column convention
#' used by the \pkg{broom} ecosystem (\code{term},
#' \code{estimate}, \code{se}, \code{ci_lower},
#' \code{ci_upper}, \code{conf_level}, \code{wald_df}). The \code{term} is the
#' coefficient name (\code{"alpha"}, \code{"omega"}, etc.),
#' \code{conf_level} is the level of the interval, so that
#' \code{\link{results_sentence}} on the row states it, and
#' \code{wald_df} the degrees of freedom of the \emph{t} reference of a
#' Wald interval (the \code{wald_df} argument of
#' \code{\link{reliability_alpha}}, \code{\link{reliability_omega}}, and
#' \code{\link{reliability}}), \code{Inf} for the normal reference or an
#' interval that is not a Wald interval, so that rows on the two
#' references can be told apart after they are stacked.
#'
#' @param x A \code{dmar_reliability} object returned by any of
#'   \code{\link{reliability_alpha}}, \code{\link{reliability_kr20}},
#'   \code{\link{reliability_omega}},
#'   \code{\link{reliability_omega_categorical}}, or \code{\link{reliability}}.
#' @param \dots Unused.
#'
#' @return A one-row \code{data.frame}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Coefficient alpha for the three verbal tests of the Holzinger and
#' # Swineford battery, from their covariance matrix.
#' S <- cov(holzinger_swineford[, c("t6_paragraph_comprehension",
#'                                  "t7_sentence", "t9_word_meaning")])
#' res <- reliability_alpha(S = S, N = 301, ci_method = "feldt")
#' generics::tidy(res)
#' generics::glance(res)
#'
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_reliability <- function(x, ...) {
  pick <- function(term) {
    val <- x$value[x$term == term]
    if (length(val) == 0L) NA_real_ else val[1L]
  }
  # The level of the interval travels with it, as glance() reports it, so
  # that results_sentence() on the tidy row states the level of its limits.
  data.frame(
    term      = attr(x, "coefficient") %||% "reliability",
    estimate  = pick("estimate"),
    se = pick("se"),
    ci_lower  = pick("ci_lower"),
    ci_upper = pick("ci_upper"),
    conf_level = pick("conf_level"),
    wald_df = .wald_df_column(attr(x, "wald_df", exact = TRUE), 1L),
    h1_diagnostic = attr(x, "h1_diagnostic", exact = TRUE) %||% NA_character_,
    wald_calibration_note = attr(x, "wald_calibration_note", exact = TRUE) %||% NA_character_,
    stringsAsFactors = FALSE
  )
}


#' Glance at a Reliability Coefficient Estimate
#'
#' Returns a one-row \code{data.frame} of model-level summaries in
#' the column convention used by the \pkg{broom} ecosystem
#' (\code{coefficient}, \code{estimate}, \code{se},
#' \code{ci_lower}, \code{ci_upper}, \code{conf_level},
#' \code{nobs}, \code{n_items}, \code{ci_method}, \code{wald_df}, the
#' degrees of freedom of the \emph{t} reference of a Wald interval,
#' \code{Inf} for the normal reference or another interval).
#'
#' @param x A \code{dmar_reliability} object.
#' @param \dots Unused.
#'
#' @return A one-row \code{data.frame}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_reliability <- function(x, ...) {
  pick <- function(term) {
    val <- x$value[x$term == term]
    if (length(val) == 0L) NA_real_ else val[1L]
  }
  data.frame(
    coefficient = attr(x, "coefficient") %||% NA_character_,
    estimate    = pick("estimate"),
    se   = pick("se"),
    ci_lower    = pick("ci_lower"),
    ci_upper   = pick("ci_upper"),
    conf_level  = pick("conf_level"),
    nobs        = pick("N"),
    n_items     = pick("n_items"),
    ci_method   = attr(x, "ci_method") %||% NA_character_,
    wald_df     = .wald_df_column(attr(x, "wald_df", exact = TRUE), 1L),
    h1_diagnostic = attr(x, "h1_diagnostic", exact = TRUE) %||% NA_character_,
    wald_calibration_note = attr(x, "wald_calibration_note", exact = TRUE) %||% NA_character_,
    stringsAsFactors = FALSE
  )
}


# `%||%` from rlang; reimplemented locally so we do not need to
# import rlang for this single use.
`%||%` <- function(a, b) if (is.null(a)) b else a
