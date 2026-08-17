#' Reliability Coefficient With a Confidence Interval (General Dispatch)
#'
#' @description
#' General-purpose entry point for the reliability family. Dispatches to
#' \code{\link{reliability_alpha}},
#' \code{\link{reliability_kr20}},
#' \code{\link{reliability_omega}}, or
#' \code{\link{reliability_omega_categorical}} according to the requested
#' \code{type}. When \code{type} is not specified, the function picks a
#' reasonable default from the supplied input following the
#' recommendations of Kelley and Pornprasertmanit (2016).
#'
#' @details
#' Auto-detection rules (used only when \code{type = NULL}):
#' \itemize{
#'   \item If raw items are integer-valued and every column has at most
#'         10 distinct values, \code{type = "omega_categorical"} (categorical
#'         omega; appropriate when items are ordered-categorical and
#'         the relationship between the underlying factor and the
#'         observed items is non-linear).
#'   \item Otherwise, \code{type = "omega"} (McDonald's coefficient
#'         \eqn{\omega} from a single-factor CFA). This includes the
#'         case where only a covariance matrix \code{S} and sample
#'         size \code{N} are supplied; \pkg{lavaan} fits the CFA on
#'         the covariance matrix.
#' }
#' Auto-detection emits a single \code{message()} indicating which
#' \code{type} was chosen, so it never surprises the user silently.
#'
#' The selected family function determines which \code{ci_method} values
#' are accepted; see the help page for the chosen function for the full
#' list. When \code{ci_method} is left at its default (\code{NULL}), the
#' family function's own default is used:
#' \itemize{
#'   \item \code{reliability_alpha}: \code{"bonett"}.
#'   \item \code{reliability_alpha(estimator = "model_implied")}: \code{"mlr"} with raw
#'         data; \code{"ml"} with covariance input.
#'   \item \code{reliability_kr20}: \code{"feldt"}.
#'   \item \code{reliability_omega}: for
#'         \code{denominator = "observed"} (robust omega, the default),
#'         the point estimate with no interval, since its interval is
#'         bootstrap based and no bootstrap runs unless requested; for
#'         \code{denominator = "model_implied"}, \code{"mlr"}.
#'   \item \code{reliability_omega_categorical}: the point estimate with
#'         no interval, for the same reason; request \code{"bca"}.
#' }
#' A bootstrap is never run by default anywhere in the family. When a
#' bootstrap method is requested, \code{B = 10000} replications is the
#' default.
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
#' For ordered-categorical items the relationship between the latent
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
#'   accepted as a shorthand), or \code{NULL}
#'   for auto-detection. See \emph{Details}.
#' @param estimator For \code{type = "alpha"} only: how coefficient
#'   \eqn{\alpha} is estimated, \code{"analytic"} (default; the
#'   closed-form equation applied to the observed covariance matrix) or
#'   \code{"model_implied"} (the reliability implied by the
#'   \eqn{\tau}-equivalent single-factor model fit by maximum
#'   likelihood); forwarded to \code{\link{reliability_alpha}}, whose
#'   help page discusses the choice and which interval methods each
#'   estimator supports. Supplying it with any other \code{type} is an
#'   error, as \code{denominator} is for any type but \code{"omega"}.
#' @param denominator For \code{type = "omega"} only: how the total
#'   variance in the denominator of \eqn{\omega} is estimated,
#'   \code{"observed"} (default; robust omega) or
#'   \code{"model_implied"}; forwarded to
#'   \code{\link{reliability_omega}}, whose help page discusses the
#'   choice.
#' @param missing For \code{type = "alpha"} and \code{type = "omega"}
#'   only: how incomplete rows of \code{data} are handled,
#'   \code{"listwise"} (the default) or \code{"fiml"} (full information
#'   maximum likelihood); forwarded to the family function, whose help
#'   page discusses the choice. Supplying it with any other \code{type}
#'   is an error.
#' @param aux For \code{type = "alpha"} and \code{type = "omega"} only:
#'   optional character vector naming auxiliary variable columns of
#'   \code{data}, entered as saturated correlates under full
#'   information maximum likelihood; forwarded to the family function.
#'   Supplying \code{aux} implies \code{missing = "fiml"}.
#' @param ci_method Method for constructing the confidence interval, or
#'   \code{NULL} to use the chosen family function's default. See the
#'   help page for the chosen \code{reliability_*} function for the full
#'   list of accepted values.
#' @param conf_level Confidence level. Defaults to \code{0.95}.
#' @param B Number of bootstrap replications when a bootstrap method is
#'   selected. Defaults to \code{10000}.
#' @param seed Random number seed used for bootstrap reproducibility.
#'   Defaults to \code{NULL}, which leaves the user's current RNG state intact; supply an integer for reproducibility.
#'
#' @return The \code{data.frame} returned by the dispatched
#' \code{reliability_*} function (rows: \code{estimate}, \code{se},
#' \code{lower_limit}, \code{upper_limit}, \code{conf_level}, \code{N},
#' \code{N_complete}, \code{J}). The \code{se} row is on the coefficient
#' scale; the transformation-based intervals (\code{"fisher"},
#' \code{"bonett"}, \code{"hakstian_whalen"}) add an
#' \code{se_transformed} row carrying the transformation-scale standard
#' error, with the scale named in the \code{se_transform_scale}
#' attribute. The \code{coefficient} attribute
#' identifies which coefficient was computed.
#'
#' @references
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
#' # auto-detection again picks coefficient omega. The call fits the
#' # single factor model a second time, so it is shown rather than run:
#' #   reliability(S = cov(items), N = 200)
#'
#' @keywords htest multivariate
#' @family reliability
#'
#' @export

reliability <- function(data = NULL, S = NULL, N = NULL,
                        type = NULL,
                        estimator = c("analytic", "model_implied"),
                        denominator = c("observed", "model_implied"),
                        missing = c("listwise", "fiml"),
                        aux = NULL,
                        ci_method = NULL,
                        conf_level = 0.95,
                        B = 10000,
                        seed = NULL) {

  type <- .resolve_relia_type(type, data, S)
  est_supplied <- !missing(estimator)
  estimator <- match.arg(estimator)
  if (est_supplied && type != "alpha") {
    stop("'estimator' applies only to type = \"alpha\".", call. = FALSE)
  }
  denom_supplied <- !missing(denominator)
  denominator <- match.arg(denominator)
  if (denom_supplied && type != "omega") {
    stop("'denominator' applies only to type = \"omega\".", call. = FALSE)
  }
  missing_supplied <- !missing(missing)
  missing <- match.arg(missing)
  if ((missing_supplied || !is.null(aux)) &&
      !type %in% c("alpha", "omega")) {
    stop("'missing' and 'aux' apply only to type = \"alpha\" and ",
         "type = \"omega\"; the FIML treatment has not been extended ",
         "to type = \"", type, "\".", call. = FALSE)
  }

  # When the user did not specify a CI method, the family function's own
  # default applies, with one exception: reliability_omega with the
  # model implied denominator defaults to the robust ML CI ("mlr"), which
  # requires raw data, and the model implied alpha estimator shares that
  # default. If the wrapper is called with only a covariance matrix,
  # fall back to the maximum likelihood CI ("ml") so the call succeeds
  # without forcing the user to know which methods need raw data.
  # (With denominator = "observed", reliability_omega already resolves
  # its own default to no interval with a message.)
  if (is.null(ci_method) && is.null(data) &&
      ((type == "omega" && denominator == "model_implied") ||
       (type == "alpha" && estimator == "model_implied"))) {
    ci_method <- "ml"
  }

  # Build the argument list for the chosen family function. ci_method is
  # forwarded only when the caller specified it (or the wrapper has
  # picked a sensible override above); otherwise the family function's
  # own default applies. `missing` is forwarded only when the caller
  # specified it, so the family function still sees the aux-implies-fiml
  # shorthand (an unsupplied `missing` plus `aux`) exactly as a direct
  # call would.
  common <- list(conf_level = conf_level, B = B, seed = seed)
  if (!is.null(ci_method)) common$ci_method <- ci_method
  if (missing_supplied) common$missing <- missing
  if (!is.null(aux)) common$aux <- aux

  if (type == "alpha") {
    args <- c(list(data = data, S = S, N = N, estimator = estimator), common)
    return(do.call(reliability_alpha, args))
  }
  if (type == "kr20") {
    if (!is.null(S)) {
      stop("type = \"kr20\" requires raw 'data'; a covariance matrix is ",
           "not sufficient.", call. = FALSE)
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
      stop("type = \"omega_categorical\" requires raw 'data'; a covariance matrix ",
           "is not sufficient.", call. = FALSE)
    }
    args <- c(list(data = data), common)
    return(do.call(reliability_omega_categorical, args))
  }
}

.resolve_relia_type <- function(type, data, S) {
  if (!is.null(type)) {
    type <- match.arg(type, c("alpha", "kr20", "omega",
                              "omega_categorical", "omega_c"))
    # "omega_c" is accepted as a shorthand for the full name.
    if (type == "omega_c") type <- "omega_categorical"
    return(type)
  }

  if (is.null(data) && is.null(S)) {
    stop("Either 'data' or 'S' must be supplied so reliability() can ",
         "select an appropriate coefficient.", call. = FALSE)
  }

  # When raw items are available, check whether they look ordered-
  # categorical (integer-valued and few distinct values per item). Only
  # in that case do we override the default of coefficient omega; the
  # categorical-omega estimator treats the items via a probit-link
  # single-factor model, which the continuous-item omega cannot do.
  if (!is.null(data)) {
    dat <- as.matrix(data)
    is_integer_valued <- isTRUE(all(dat == round(dat), na.rm = TRUE))
    max_unique <- max(apply(dat, 2, function(x)
                            length(unique(stats::na.omit(x)))))
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
#' \code{ci_upper}). The \code{term} is the coefficient name
#' (\code{"alpha"}, \code{"omega"}, etc.).
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
  data.frame(
    term      = attr(x, "coefficient") %||% "reliability",
    estimate  = pick("estimate"),
    se = pick("se"),
    ci_lower  = pick("lower_limit"),
    ci_upper = pick("upper_limit"),
    stringsAsFactors = FALSE
  )
}


#' Glance at a Reliability Coefficient Estimate
#'
#' Returns a one-row \code{data.frame} of model-level summaries in
#' the column convention used by the \pkg{broom} ecosystem
#' (\code{coefficient}, \code{estimate}, \code{se},
#' \code{ci_lower}, \code{ci_upper}, \code{conf_level},
#' \code{nobs}, \code{n_items}, \code{ci_method}).
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
    ci_lower    = pick("lower_limit"),
    ci_upper   = pick("upper_limit"),
    conf_level  = pick("conf_level"),
    nobs        = pick("N"),
    n_items     = pick("J"),
    ci_method   = attr(x, "ci_method") %||% NA_character_,
    stringsAsFactors = FALSE
  )
}


# `%||%` from rlang; reimplemented locally so we do not need to
# import rlang for this single use.
`%||%` <- function(a, b) if (is.null(a)) b else a
