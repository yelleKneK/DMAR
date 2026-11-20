#' Compare Structural Equation Models by an Information Criterion
#'
#' Ranks a set of structural equation models fitted to the same data by one
#' information criterion of \code{\link{sem_information_criteria}} (BIC by
#' default) and reports, for each model, its difference from the best, its
#' relative weight within the set, and, for the criteria that approximate
#' a log marginal likelihood, the approximate log Bayes factor of the best
#' model over it. The models may be nested or not. A criterion compares
#' models only through their log likelihoods, so before anything is
#' computed the function checks that the log likelihoods are comparable:
#' the same cases and data values, the same observed variables and groups,
#' the same missing data method and mean structure, and the same exogenous
#' observed variables conditioned on. Requires \pkg{lavaan} (in Suggests).
#'
#' @param fits A named list of two or more fitted \pkg{lavaan} models:
#'   maximum likelihood fits of continuous variables (estimator ML or its
#'   robust variants MLM, MLMV, and MLR, which maximize the same likelihood)
#'   fitted to the same data. That is, the same cases, data values,
#'   observed variables (in any order), and groups; the same missing data
#'   method (\code{missing}) and mean structure (\code{meanstructure}); and,
#'   when \code{fixed.x} or \code{conditional.x} is \code{TRUE}, the same
#'   exogenous observed variables. Models fitted to summary statistics
#'   (\code{sample.cov}) are compared when they share the sample moments and
#'   \emph{N}. Two-level models and models fitted with sampling weights are
#'   refused, and the criterion can impose further conditions (see
#'   \code{\link{sem_information_criteria}}).
#' @param criterion The criterion: \code{"BIC"} (the default),
#'   \code{"SABIC"}, \code{"HBIC"}, \code{"IBIC"}, \code{"SPBIC"}, or
#'   \code{"AIC"}, in either case.
#' @param prior_means For SPBIC, an optional list of prior mean vectors, one
#'   per model in the order of \code{fits} (and with its names, if named),
#'   each in the order of \code{lavaan::coef()} of its model. \code{NULL}
#'   applies the default rule of \code{\link{sem_information_criteria}}:
#'   free loadings and regression coefficients centered at zero and the
#'   other free parameters at their estimates. A saturated model needs an
#'   explicit prior mean.
#'
#' @details
#' \subsection{Differences, Weights, and Evidence}{
#' For the criterion values \eqn{C_k} of the models, \eqn{\Delta_k = C_k -
#' \min_j C_j} (column \code{delta}), and the relative weight of model
#' \emph{k} is \eqn{\exp(-\Delta_k/2)/\sum_j \exp(-\Delta_j/2)}. BIC, HBIC,
#' IBIC, and SPBIC each approximate \eqn{-2} times the log marginal
#' likelihood of a model (Bollen, Ray, Zavisca, & Harden, 2012), so
#' \eqn{\exp(-\Delta_k/2)} approximates the Bayes factor of model \emph{k}
#' against the best model, the weights approximate the posterior
#' probabilities of the models when every model in the set has the same
#' prior probability (Burnham & Anderson, 2004), and \eqn{\Delta_k/2}
#' (column \code{log_evidence_best_over_model}) approximates the log Bayes
#' factor of the best model over model \emph{k}. For AIC the weights are
#' Akaike weights, the weight of evidence that each model is the one of the
#' set closest to the data generating process in Kullback-Leibler
#' discrepancy (Burnham & Anderson, 2004); they are not posterior
#' probabilities under equal prior odds, and the log evidence column is
#' \code{NA}. SABIC rests on the minimum description length principle
#' rather than on a marginal likelihood, so its weights carry no further
#' interpretation and its log evidence column is \code{NA} too. So it is
#' for every criterion of models fitted by two-stage maximum likelihood
#' (\code{missing = "two.stage"} or \code{"robust.two.stage"}), whose
#' criteria are those of Lai (2021): built from the complete data normal
#' log likelihood at the EM estimates of the means and covariances, they
#' select in large samples the model complete data would select, but they
#' are not a likelihood of the observed data and approximate no marginal
#' likelihood (see \code{\link{sem_information_criteria}}).
#'
#' Every weight is conditional on the set of models compared: adding or
#' removing a model changes all of them, and none of them says how well a
#' model fits in absolute terms (\code{\link{sem_fit_measures}} speaks to
#' that). Nor do the weights carry the sampling variability of the
#' criteria; the uncertainty of the selection itself is not estimated.
#' }
#'
#' \subsection{What Is Checked}{
#' Each model must have converged to an admissible solution and have a log
#' likelihood the criteria are defined on (see
#' \code{\link{sem_information_criteria}} for what is refused and why).
#' Then every model is compared with the first, and the first difference
#' found is reported with the two models it concerns and what to do.
#' \describe{
#' \item{Data}{The same groups, the same observed variables, the same
#'   number of cases, the same cases, and the same data values (or, for
#'   models fitted to summary statistics, the same sample moments). The
#'   order of the variables does not matter; a set of different variables
#'   does, since a log likelihood is a density of the variables it models.}
#' \item{Missing data}{The same missing data method. Listwise deletion,
#'   FIML, and the two-stage method give different likelihoods, and two
#'   listwise fits of different variables keep different cases. The two
#'   two-stage settings, \code{missing = "two.stage"} and
#'   \code{"robust.two.stage"}, count as one method: they fit the same EM
#'   moments by the same likelihood and differ only in their standard errors
#'   and test statistic, which no criterion uses. Two-stage fits must also
#'   share the auxiliary variables of their first stage (lavaan's
#'   \code{aux}), since those change the EM estimates their log likelihoods
#'   are evaluated at. A two-stage fit whose first stage read cases
#'   observed on an auxiliary variable but on none of the model's
#'   variables, which lavaan leaves out of its \emph{N}, is compared by
#'   the criterion \code{\link{sem_information_criteria}} computes for it
#'   with every case the first stage read, its log likelihood and the
#'   \emph{N} of the penalty included (in several groups held together by
#'   a constraint, at the estimates with every case read, which move with
#'   the groups' counts), and a line under the table says so; such a fit and
#'   the fit \code{\link{fa_confirmatory}(missing = "two_stage")} or
#'   \code{\link{measurement_invariance}} returns for the same data count
#'   the same cases.}
#' \item{Mean structure}{The same \code{meanstructure} setting: without a
#'   mean structure lavaan leaves the saturated means out of the parameter
#'   count, which moves every criterion by an amount that depends on the
#'   setting.}
#' \item{Exogenous variables}{With \code{fixed.x = TRUE} (the default of
#'   \code{lavaan::sem()}) or \code{conditional.x = TRUE}, lavaan's log
#'   likelihood is the density conditional on the exogenous observed
#'   variables, the joint density less the marginal density of those
#'   variables, and its parameter count leaves out their moments. Two
#'   models that treat different variables as exogenous therefore have log
#'   likelihoods of different conditional distributions, whose difference
#'   can be hundreds of points for two saturated models of the same three
#'   variables. Models with the same exogenous variables are compared, and
#'   the others are refused with the advice to refit every model with
#'   \code{fixed.x = FALSE}.}
#' \item{Identification}{For IBIC, and for SPBIC under its default prior,
#'   every latent variable of every model must be identified the same way
#'   (each by a marker indicator, or each by a fixed variance), since both
#'   criteria change with the identification even though the likelihood
#'   does not (see \code{\link{sem_information_criteria}}). A set that uses
#'   both is refused, also when each model mixes them in the same way: one
#'   identification throughout is a convention that can be stated and
#'   reported, and it catches the models that differ, which may share no
#'   latent variable to compare. Within one data set the units of the
#'   observed variables are shared, but an IBIC ranking still depends on
#'   them.}
#' }
#' }
#'
#' @return A \code{dmar_tbl} (class \code{dmar_compare_sem}) with one row
#'   per model, sorted by the criterion (ties keep the order of
#'   \code{fits}), and the columns \code{term} (the model's name),
#'   \code{value} (its criterion), \code{delta}, \code{relative_weight},
#'   \code{log_evidence_best_over_model}, and \code{npar} (its number of
#'   free parameters). The \code{criterion} attribute names the criterion,
#'   the \code{comparison} attribute records what was checked (\emph{N}, the
#'   number of groups, whether the models were fitted to raw data, the
#'   \code{missing} settings used, the mean structure, the exogenous observed
#'   variables conditioned on, and for IBIC and SPBIC the identification),
#'   and the \code{sem_information} attribute holds the details of each
#'   model's \code{\link{sem_information_criteria}} result, named by model.
#'   \code{tidy()} returns every column with the criterion beside it, and
#'   \code{glance()} one row with the criterion, the number of models,
#'   \emph{N}, the best model, its weight, and the difference of the second
#'   best.
#'
#' @references
#' Bollen, K. A., Harden, J. J., Ray, S., & Zavisca, J. (2014). BIC and
#'   alternative Bayesian information criteria in the selection of
#'   structural equation models. \emph{Structural Equation Modeling,
#'   21}(1), 1--19. \doi{10.1080/10705511.2014.856691}
#'
#' Bollen, K. A., Ray, S., Zavisca, J., & Harden, J. J. (2012). A
#'   comparison of Bayes factor approximation methods including two new
#'   methods. \emph{Sociological Methods & Research, 41}(2), 294--324.
#'   \doi{10.1177/0049124112452393}
#'
#' Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
#'   Understanding AIC and BIC in model selection. \emph{Sociological
#'   Methods & Research, 33}(2), 261--304. \doi{10.1177/0049124104268644}
#'
#' Lai, K. (2021). Using information criteria under missing data: Full
#'   information maximum likelihood versus two-stage estimation.
#'   \emph{Structural Equation Modeling, 28}(2), 278--291.
#'   \doi{10.1080/10705511.2020.1780925}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{sem_information_criteria}} for the criteria of one
#'   model and what each assumes, \code{\link{sem_fit_measures}} for the fit
#'   of each model in absolute terms, and the vignette
#'   \code{vignette("sem_fit_measures", package = "DMAR")}.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # One factor against two correlated factors for six of the
#' # holzinger_swineford tests.
#' one <- lavaan::cfa(
#'   "g =~ t1_visual_perception + t2_cubes + t4_lozenges +
#'         t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
#'   data = holzinger_swineford)
#' two <- lavaan::cfa(
#'   "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
#'    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
#'   data = holzinger_swineford)
#' x <- compare_sem(list(one_factor = one, two_factors = two))
#' x
#' generics::glance(x)
#'
#' # The same comparison by SPBIC, under its default prior.
#' compare_sem(list(one_factor = one, two_factors = two), criterion = "SPBIC")
#'
#' @export
compare_sem <- function(fits, criterion = c("BIC", "SABIC", "HBIC", "IBIC", "SPBIC", "AIC"),
                        prior_means = NULL) {
  criterion <- if (identical(criterion, c("BIC", "SABIC", "HBIC", "IBIC", "SPBIC", "AIC"))) {
    "BIC"
  } else .sem_ic_match(criterion, "criterion", several = FALSE)
  if (!is.list(fits) || inherits(fits, "lavaan") || length(fits) < 2L ||
      is.null(names(fits)) || anyNA(names(fits)) || any(!nzchar(names(fits))) ||
      anyDuplicated(names(fits))) {
    .dmar_refuse_argument(
      "'fits' must be a uniquely named list of at least two models.",
      arg = "fits")
  }
  if (!is.null(prior_means)) {
    if (criterion != "SPBIC") {
      .dmar_refuse_argument(
        "'prior_means' is used only with criterion = \"SPBIC\".",
        arg = "prior_means")
    }
    if (!is.list(prior_means) || length(prior_means) != length(fits)) {
      .dmar_refuse_argument(c(
        "'prior_means' must be a list with one prior mean vector per model, ",
        "in the order of 'fits'."), arg = c("prior_means", "fits"))
    }
    if (!is.null(names(prior_means)) && !identical(names(prior_means), names(fits))) {
      .dmar_refuse_argument(
        "Named 'prior_means' must match the names and order of 'fits'.",
        arg = c("prior_means", "fits"))
    }
  }
  labels <- paste0("fits$", names(fits))
  liks <- lapply(seq_along(fits), function(i) {
    options <- .sem_fit_check(fits[[i]], labels[i])
    lik <- .sem_likelihood_check(fits[[i]], options, labels[i], criterion)
    if (lavaan::lavInspect(fits[[i]], "nlevels") > 1L) {
      .dmar_refuse_argument(c(
        "The model in ", labels[i], " is a two-level model, which ",
        "compare_sem() does not compare; sem_information_criteria() gives ",
        "its AIC."), arg = NULL)
    }
    lik
  })
  signature <- stats::setNames(
    lapply(seq_along(fits), function(i) .sem_comparison_signature(fits[[i]], liks[[i]])),
    names(fits))
  .sem_comparison_check(signature)
  identification <- NULL
  if (criterion == "IBIC" || (criterion == "SPBIC" && is.null(prior_means))) {
    identification <- lapply(seq_along(fits), function(i) {
      .sem_information_fit_check(fits[[i]], lavaan::lavInspect(fits[[i]], "options"),
                                 criterion, labels[i])
      .sem_identification(fits[[i]])
    })
    names(identification) <- names(fits)
    kinds <- lapply(identification, unique)
    used <- unique(unlist(kinds))
    # The rule is one identification for every latent variable of every
    # model, so a set is refused when it uses more than one, whether the
    # models differ or each model mixes the two in the same way.
    if (length(used) > 1L) {
      describe <- vapply(seq_along(kinds), function(i) {
        paste0(labels[i], ": ", if (length(kinds[[i]])) paste(kinds[[i]], collapse = " and ")
               else "no latent variables")
      }, character(1))
      .dmar_refuse_argument(c(
        criterion, " depends on how the latent variables are identified, so ",
        "compare_sem() requires every latent variable of every model to be ",
        "identified the same way, and these models use more than one ",
        "identification (", paste(describe, collapse = "; "),
        "). Refit every model with one identification: std.lv = TRUE for all ",
        "(each latent variance fixed at 1) or the marker indicator for all ",
        "(lavaan's default)."), arg = NULL)
    }
  }
  tables <- lapply(seq_along(fits), function(i) {
    withCallingHandlers(
      sem_information_criteria(fit = fits[[i]], measures = criterion,
                               prior_mean = if (is.null(prior_means)) NULL else prior_means[[i]]),
      dmar_two_stage_aux_only = function(w) invokeRestart("muffleWarning"))
  })
  # A two-stage fit whose first stage read cases its N leaves out (cases
  # observed on an auxiliary variable alone) has its log likelihood and N
  # computed with every case the first stage read by
  # sem_information_criteria() (.ts_all_rows() in R/two_stage_internals.R).
  # Only when that could not be done does the fit get a warning, which
  # names it, in place of the one its sem_information_criteria() call above
  # gave unnamed and was muffled. It is given once every criterion has been
  # computed, so a comparison that is refused gives no warning before its
  # error. Through the first versions of 1.1.0 every such fit was warned
  # about and its criterion was lavaan's.
  ts_rows <- lapply(tables, function(x) attr(x, "sem_information")$two_stage_rows)
  for (i in seq_along(fits)) {
    if (!identical(ts_rows[[i]]$corrected, FALSE)) next
    .ts_aux_only_warning(
      fits[[i]],
      rest = paste0("its log likelihood, which lavaan forms as N times the ",
                    "complete data log likelihood per case, ",
                    if (criterion == "AIC") "rests" else
                      paste0("and the sample size in the ", criterion,
                             " penalty rest")),
      consequence = paste0("its ", criterion, " weighs fit against ",
                           "complexity at a sample size other than the one ",
                           "the moments describe (DMAR could not compute ",
                           "the fit again with every case)"),
      where = paste0("the two-stage fit ", labels[i]))
  }
  vals <- vapply(tables, function(x) x$value[x$term == criterion], numeric(1))
  npar <- vapply(tables, function(x) x$value[x$term == "npar"], numeric(1))
  delta <- vals - min(vals)
  weights <- exp(-delta / 2)
  # The two-stage criteria of Lai (2021) are built from the complete data
  # log likelihood at the EM estimates, not from a likelihood of the
  # observed data, so they approximate no marginal likelihood and no
  # Bayes factor is reported for them.
  two_stage <- signature[[1L]]$missing %in% c("two.stage", "robust.two.stage")
  evidence <- criterion %in% c("BIC", "HBIC", "IBIC", "SPBIC") && !two_stage
  out <- data.frame(term = names(fits), value = vals, delta = delta,
                    relative_weight = weights / sum(weights),
                    log_evidence_best_over_model = if (evidence) delta / 2 else NA_real_,
                    npar = npar, row.names = NULL, stringsAsFactors = FALSE)
  out <- out[order(out$value), , drop = FALSE]
  rownames(out) <- NULL
  first <- signature[[1L]]
  # Both two-stage settings may appear in one comparison (see
  # .sem_missing_method()); the record names each setting used.
  missing <- unique(vapply(signature, function(s) s$missing, character(1)))
  # The N of the criteria: with cases observed on an auxiliary variable
  # alone, every case the first stage read (the same for every fit, since
  # the fits share their data and auxiliary variables).
  comparison <- list(N = tables[[1L]]$value[tables[[1L]]$term == "N"],
                     groups = first$groups, raw_data = first$raw,
                     missing = missing, aux = first$aux,
                     meanstructure = first$meanstructure,
                     conditioning = first$conditioning,
                     likelihood = liks[[1L]]$likelihood)
  aux_only <- vapply(ts_rows, function(r) {
    if (isTRUE(r$corrected)) as.numeric(r$extra) else 0
  }, numeric(1L))
  if (any(aux_only > 0)) {
    comparison$two_stage_extra <- max(aux_only)
    # A fit whose groups a constraint holds together was estimated again
    # with every case read (sem_information_criteria()).
    comparison$two_stage_reestimated <- any(vapply(ts_rows, function(r) {
      isTRUE(r$corrected) && isTRUE(r$reestimated)
    }, logical(1L)))
  }
  if (!is.null(identification)) comparison$identification <- identification
  comparison$footer <- .compare_sem_footer(criterion, comparison, nrow(out))
  # Every row of the criterion column prints to fixed decimals, where a
  # difference of a few points between models would otherwise round away.
  out <- .as_dmar_tbl(out, fixed_terms = out$term, subclass = "dmar_compare_sem")
  attr(out, "criterion") <- criterion
  attr(out, "comparison") <- comparison
  attr(out, "sem_information") <- stats::setNames(lapply(tables, attr, "sem_information"),
                                                  names(fits))
  out
}

# The lines written under a compare_sem() table: the criterion and sample,
# what was checked, and what the weights and the evidence column mean.
.compare_sem_footer <- function(criterion, comparison, n_models) {
  checked <- paste0(
    "Checked: the same ", if (comparison$raw_data) "cases and data values" else
      "sample moments", ", observed variables, ",
    if (comparison$groups > 1L) paste0(comparison$groups, " groups, ") else "",
    "missing data method (missing = ",
    paste0("\"", comparison$missing, "\"", collapse = " and "),
    if (length(comparison$aux)) paste0(", with the auxiliary variables ",
      paste(comparison$aux, collapse = ", ")) else "",
    "), and mean structure (meanstructure = ", comparison$meanstructure, "); ",
    if (length(comparison$conditioning)) paste0("conditional on the exogenous ",
      "variables ", paste(comparison$conditioning, collapse = ", "))
    else "no exogenous variables conditioned on",
    if (!is.null(comparison$identification)) paste0("; identification: ",
      paste(unique(unlist(comparison$identification)), collapse = ", ")) else "",
    ".")
  two_stage <- any(comparison$missing %in% c("two.stage", "robust.two.stage"))
  weights <- if (two_stage && criterion != "AIC") {
    paste0("relative_weight: exp(-delta/2) normalized within this set; the ",
           "two-stage criteria of Lai (2021) rest on the complete data log ",
           "likelihood at the EM estimates rather than on a likelihood of the ",
           "observed data, so no evidence is reported.")
  } else switch(criterion,
    AIC = paste0("relative_weight: the Akaike weight of each model within this ",
                 "set, not a posterior probability under equal prior odds."),
    SABIC = paste0("relative_weight: exp(-delta/2) normalized within this set; ",
                   "SABIC approximates no marginal likelihood, so no evidence ",
                   "is reported."),
    paste0("relative_weight: approximate posterior probability of each model ",
           "within this set under equal prior odds; log_evidence_best_over_model: ",
           "approximate log Bayes factor of the best model over each."))
  c(sprintf("Criterion: %s%s (smaller is preferred); N = %s; %d models.",
            criterion, if (two_stage) ", the two-stage criterion of Lai (2021)" else "",
            format(comparison$N), n_models),
    checked, weights,
    if (!is.null(comparison$two_stage_extra)) {
      paste0("N and the log likelihoods count the ", comparison$two_stage_extra,
             " cases the first stage read on an auxiliary variable alone, ",
             "which lavaan's N leaves out",
             if (isTRUE(comparison$two_stage_reestimated)) {
               paste0(", at the estimates with every case read of a fit ",
                      "whose groups a constraint holds together")
             },
             " (see sem_information_criteria()).")
    })
}

# Every display of a compare_sem() table (the console, knit_print(),
# as_kable()) writes these lines with it; see .dmar_tbl_footer() in
# R/dmar_output_helpers.R.
#' @exportS3Method .dmar_tbl_footer dmar_compare_sem
#' @noRd
.dmar_tbl_footer.dmar_compare_sem <- function(x, ...) {
  .dmar_footer_add(NextMethod(), attr(x, "comparison")$footer, "compare_sem")
}

#' @rdname compare_sem
#' @param x A \code{dmar_compare_sem} table returned by \code{compare_sem()}.
#' @param ... Unused.
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_compare_sem <- function(x, ...) {
  data.frame(term = as.character(x$term), estimate = as.numeric(x$value),
             delta = x$delta, relative_weight = x$relative_weight,
             log_evidence_best_over_model = x$log_evidence_best_over_model,
             npar = x$npar, criterion = attr(x, "criterion") %||% NA_character_,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' @rdname compare_sem
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_compare_sem <- function(x, ...) {
  data.frame(criterion = attr(x, "criterion") %||% NA_character_,
             n_models = nrow(x),
             N = attr(x, "comparison")$N %||% NA_real_,
             best_model = as.character(x$term[1L]),
             best_weight = x$relative_weight[1L],
             delta_second = if (nrow(x) > 1L) x$delta[2L] else NA_real_,
             stringsAsFactors = FALSE)
}
