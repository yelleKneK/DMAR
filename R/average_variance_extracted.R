#' Average Variance Extracted (AVE)
#'
#' The average variance extracted (AVE) is the mean proportion of indicator
#' variance a factor accounts for, a standard convergent validity summary
#' for a reflective measurement block in confirmatory factor analysis and
#' structural equation modeling. With standardized loadings \eqn{\ell_j},
#' \deqn{\mathrm{AVE} = \frac{1}{J} \sum_j \ell_j^2.}
#' The quantity itself is elementary. In a model with cross-loadings, an
#' item contributes its loading to the AVE of every factor it loads on;
#' the same shared variance then counts toward each factor's summary, so
#' compare AVE values across factors of such a model with that overlap in
#' mind. A loading the model holds at zero, fixed (\code{0*item}) or
#' constrained (\code{a*item} with \code{a == 0}), is not a loading: the
#' model written with it is the model written without it, and the item is
#' not one of the factor's indicators. Fornell and Larcker (1981) are
#' credited for establishing it as a validity criterion: a construct shows
#' convergent validity when its AVE reaches the conventional 0.50 (the
#' construct explains at least half its indicators' variance), and the
#' Fornell-Larcker discriminant criterion compares each construct's AVE with
#' its squared correlations with the other constructs. AVE is closely
#' related to composite reliability (omega); for discriminant validity,
#' \code{\link{htmt}} gives the heterotrait-monotrait ratio of
#' correlations, a criterion that does not use the AVE.
#'
#' @param fit Optional \pkg{lavaan} fit of one group and one level (for
#'   example from \code{\link{cfa_1}} or \code{lavaan::cfa}); the
#'   standardized loadings of every factor are extracted and one AVE is
#'   reported per factor. A fit whose optimizer did not converge, and a
#'   multiple-group or multilevel fit, are refused (see \emph{The fit}
#'   under Details).
#' @param loadings Optional numeric vector of standardized loadings for a
#'   single block, as an alternative to \code{fit}. Supply exactly one of
#'   \code{fit} and \code{loadings}. No interval can be constructed
#'   from loadings alone (their sampling variability is not carried by
#'   the numbers), so \code{ci_method = "percentile"} requires
#'   \code{fit}.
#' @param conf_level Confidence level for the bootstrap interval
#'   (default \code{0.95}); used when \code{ci_method = "percentile"}.
#' @param ci_method Interval method: \code{"none"} (the default) or
#'   \code{"percentile"}. No closed-form interval for the AVE is in
#'   common use; the percentile bootstrap is the standard route in the
#'   validity literature. The cases in the fitted data are resampled
#'   with replacement \code{B} times, the model is refit to each
#'   resample, and each factor's interval is the pair of empirical
#'   quantiles of its \code{B} AVE values (Efron & Tibshirani, 1993).
#'   Replications whose refit fails or does not converge are set aside,
#'   and the interval is computed from those that return a value; a
#'   single warning reports how many were set aside, and gives the
#'   caution described under \emph{Details} when it holds. A replication
#'   whose solution is improper for a factor is kept and counted (see
#'   \emph{Details}).
#'   An interval requires that no more than half of the attempted
#'   replications were set aside and that at least the number the level
#'   requires returned an estimate (400 at the 95 percent level; see
#'   \code{B}). Otherwise the point estimates and counts are retained,
#'   with \code{NA} limits and a warning of class
#'   \code{dmar_bootstrap_insufficient}. These floors limit the error of
#'   drawing a finite number of resamples; they are not coverage
#'   guarantees.
#' @param B Number of bootstrap replications when \code{ci_method =
#'   "percentile"}, a single whole number. Default \code{1000}, smaller than
#'   the package's usual \code{10000} because every replication refits the
#'   model; raise it for a reported analysis when time allows. The fewest
#'   accepted is 400 at the default \code{conf_level} (200 at 0.90, and 2000
#'   at 0.99, which the default does not meet); a smaller or fractional
#'   \code{B} is refused before anything is computed (class
#'   \code{dmar_count_invalid}). See \emph{Details} and
#'   \code{\link{dmar_bootstrap}}.
#' @param seed Optional integer seed for the bootstrap. The default
#'   \code{NULL} uses the current state of the random number generator;
#'   a supplied seed is set internally and the prior state restored on
#'   exit.
#'
#' @details
#' \strong{The fit.} The AVE is read from a solution, so a fit whose
#' optimizer did not converge (\code{lavaan::lavInspect(fit,
#' "converged")} is \code{FALSE}, as when it reached its iteration limit)
#' is refused with an error of class \code{dmar_not_converged}: its
#' estimates mark where the optimizer stopped. A fit that converged to an
#' improper solution is not refused; it is described under
#' \emph{Improper solutions} below.
#'
#' A factor's standardized loadings are those of one group at one level.
#' A factor of a multiple-group fit therefore has an AVE in each group,
#' and a factor of a two-level fit has one at each level, even when it
#' has the same name at both. Neither the unweighted mean of those values
#' nor a mean weighted by group size is the AVE of any group, so the
#' function describes a fit of one group and one level and refuses a
#' multiple-group or multilevel fit with an error of class
#' \code{dmar_argument_invalid} that names the groups (as \code{groups})
#' or the levels (as \code{levels}). The AVE of one group or level of
#' such a model is \code{average_variance_extracted(loadings = )} given
#' the factor's standardized loadings in that group or at that level:
#' the factor's column of that group's or level's \code{lambda} matrix
#' in \code{lavaan::lavInspect(fit, "std")}, without the zeros of the
#' items that do not load on it. For a group, a fit of the model to that
#' group's data alone also gives it, without any equality constraints
#' across groups.
#'
#' \strong{Improper solutions.} When the error variance of an indicator
#' is estimated below zero (a Heywood case; for an ordered item, the
#' residual variance of its latent response), its communality, one minus
#' its standardized error variance, exceeds one, and so does its squared
#' standardized loading in a model without cross-loadings. The AVE
#' averages those communalities, so for a factor with such an indicator
#' it is not a proportion of variance: its \code{ave} (and its interval)
#' is \code{NA}, with a warning naming the indicators and their
#' communalities, as in the measurement table of \code{\link{fa_confirmatory}}.
#' On \code{holzinger_swineford}, a factor measured by visual perception
#' and lozenges alone, beside one measured by general information,
#' paragraph comprehension, and sentence completion, gives visual
#' perception an error variance below zero and a communality of 1.30,
#' so the first factor's AVE is \code{NA} (the mean of its squared
#' standardized loadings would be .726), while the verbal factor's, .705,
#' is reported. A larger sample, or reconsidering the model for those
#' indicators, are the usual remedies.
#'
#' No interval is formed for such a factor, so when every factor of the
#' fit is improper the bootstrap is not run, and the table records no
#' confidence level.
#'
#' The bootstrap meets the same problem one resample at a time. A
#' replication whose solution is improper for a factor that is proper in
#' the observed fit is kept, as the reliability family keeps one under its
#' default \code{improper = "keep"} (see \code{\link{reliability_omega}}),
#' and counted. Its AVE for that factor averages a communality above one
#' and can itself exceed one, so such replications gather in the upper
#' tail of the factor's bootstrap distribution. The line printed under the
#' table says how many replications were attempted, how many returned a
#' value, and how many of those came from an improper solution, with, for
#' two or more factors, the count for each; when the replications set
#' aside and the improper ones kept for a factor reach the tail area beyond
#' a limit, \eqn{(1 - C)/2} of \code{B} for confidence level \eqn{C}, it
#' adds a caution, since the limits then depend on how those replications
#' were treated; with two or more factors the caution gives the share of
#' each factor it holds for, counted for that factor alone. The warning
#' the call gives says the same, as the
#' reliability family's does, naming each factor the caution holds for
#' with its count. On 80 pupils drawn from \code{holzinger_swineford}
#' (\code{set.seed(113)}), a factor of visual perception, cubes, and
#' lozenges beside a verbal factor is proper in the sample, with an AVE of
#' .449, but 101 of the 198 replications that returned a value were
#' improper for it (\code{B = 200}, \code{seed = 113}), and its upper
#' limit, 1.19, is not a proportion of variance; the verbal factor had
#' none. With the 2 replications whose refit failed or did not converge,
#' 103 of the 200 are in question for the spatial factor, against the 5
#' that lie beyond each limit of a .95 interval, and the warning names it.
#' The counts are those the console and a knitted report print, and they
#' are kept in the \code{bootstrap} attribute.
#'
#' \strong{Bootstrap replications.} The percentile interval reads its two
#' limits from the quantiles of the \code{B} replicated AVE values of each
#' factor, so it is run only when \code{B} puts at least 10 replications
#' beyond each limit (400 at the 95 percent level). A replication in which
#' an AVE is undefined (the refit failed or did not converge) is dropped
#' and counted, never replaced, and the interval is computed from the rest.
#' The result records \code{B} (requested), \code{B_used} (the number the
#' interval rests on), \code{n_dropped}, and \code{B_required} (the fewest
#' replications a request may ask for; that count is held to \code{B}, not
#' to \code{B_used}) as attributes, beside
#' the fuller \code{bootstrap} list, with \code{bootstrap_status}
#' (\code{"complete"}, \code{"dropped"}, \code{"exploratory"},
#' \code{"insufficient_count"}, or \code{"insufficient_share"}). One
#' warning of class \code{dmar_bootstrap_replications_dropped} states the
#' counts when any replication was dropped, with a caution when the dropped
#' replications are at least the 2.5 percent of \code{B} that lies beyond
#' one limit at the 95 percent level, because they could then decide by
#' themselves where that limit falls. When more than half were dropped, the
#' limits of every factor are
#' \code{NA}, the estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a smaller
#' \code{B} (2 or more) run for a first look or a demonstration; every such
#' call gives a warning of class \code{dmar_bootstrap_exploratory} and the
#' result is not one to report. See \code{\link{dmar_bootstrap}} for the
#' rule and its reasons.
#'
#' The percentile bootstrap interval resamples the cases behind
#' \code{fit} and refits the model once per replication, so its cost is
#' \code{B} model fits. That refitting is why the examples below stop at
#' the point estimates: the fewest replications a 95 percent interval is
#' reported from, \code{B = 400}, run well beyond the time an example is
#' given on the two-factor model there. To obtain the
#' interval, pass \code{ci_method = "percentile"} together with a
#' \code{seed}, as in \code{average_variance_extracted(fit, ci_method =
#' "percentile", seed = 113)}, and keep the default \code{B = 1000} or
#' more for a reported analysis.
#'
#' A supplied two-stage fit requires eligible unrestricted first-stage
#' moments. Bootstrap intervals currently refuse separately stored
#' two-stage auxiliary variables because the refits do not carry their
#' resampled values. Resample and refit items and auxiliaries together
#' explicitly when that analysis is needed.
#'
#' @return A \code{data.frame} (classes \code{dmar_ave} and
#'   \code{dmar_tbl}) with one row
#'   per factor: \code{factor} (label), \code{ave}, and
#'   \code{ci_lower} / \code{ci_upper} (the percentile bootstrap
#'   limits; \code{NA} when \code{ci_method = "none"}). All three are
#'   \code{NA} for a factor whose solution is improper (see Details).
#'   With \code{ci_method = "percentile"} and at least one factor whose
#'   solution is proper, so that an interval is formed, the table carries
#'   the attributes \code{conf_level}, \code{B} (requested),
#'   \code{B_used} (the replications the intervals rest on),
#'   \code{n_dropped}, \code{B_required}, and \code{bootstrap_status}
#'   (see \emph{Bootstrap replications} under Details), and
#'   \code{bootstrap}, the counts in the form the
#'   reliability family records them (see \code{\link{reliability_omega}}):
#'   \code{attempted} (\code{B}), \code{returned}, \code{failed} (the
#'   replications whose refit failed or did not converge, which
#'   \pkg{lavaan}'s bootstrap does not tell apart, counted in
#'   \code{reasons} as non-finite), \code{reasons}, \code{error_messages},
#'   \code{undetermined}, \code{improper}, and \code{no_solution} (all
#'   zero here), \code{improper_kept} (the replications used whose
#'   solution was improper for at least one factor, including when its interval is withheld),
#'   \code{not_converged} (zero), \code{refit_failed} (the same count as
#'   \code{failed}, under the reason that the line under the table and
#'   the warning give for both causes, that the refit failed or did not
#'   converge), \code{improper_rule} (\code{"keep"}), \code{caution}, and
#'   \code{improper_kept_by_factor}, the improper replications kept for
#'   each factor, with \code{factors_with_intervals} recording which intervals were eligible.
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#'   bootstrap}. New York, NY: Chapman & Hall/CRC.
#'
#' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation
#'   models with unobservable variables and measurement error.
#'   \emph{Journal of Marketing Research, 18}(1), 39--50.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{htmt}} for discriminant validity;
#'   \code{\link{reliability_omega}} for the composite reliability of the
#'   same block (coefficient omega is what the composite-reliability
#'   literature computes); \code{\link{cfa_1}} to obtain the fit.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' # Directly from the standardized loadings a paper reports:
#' average_variance_extracted(loadings = c(.8, .7, .6))
#'
#' # From a fitted model, one AVE per factor (requires lavaan).
#' data(holzinger_swineford)
#' fit <- lavaan::cfa(
#'   "verbal    =~ t6_paragraph_comprehension + t7_sentence +
#'                 t9_word_meaning
#'    deduction =~ t20_deduction + t22_problem_reasoning +
#'                 t23_series_completion",
#'   data = holzinger_swineford)
#' ave_tbl <- average_variance_extracted(fit)
#' ave_tbl
#'
#' # The Fornell and Larcker (1981) discriminant criterion compares each
#' # factor's AVE with the squared correlation between the factors: a
#' # factor should account for more of its own indicators' variance than
#' # it shares with the other factor. Here the comparison favors verbal
#' # and goes against deduction, whose AVE falls below the shared
#' # variance. Fitting with fa_confirmatory(..., output = "measurement") puts the
#' # AVE values and the latent correlations in one table.
#' lavaan::lavInspect(fit, "cor.lv")["verbal", "deduction"]^2
#'
#' # The broom verbs: one row per factor.
#' generics::tidy(ave_tbl)
#' generics::glance(ave_tbl)
#'
#' @export
average_variance_extracted <- function(fit = NULL, loadings = NULL,
                                       conf_level = 0.95,
                                       ci_method = c("none", "percentile"),
                                       B = 1000L, seed = NULL) {
  .dmar_check_conf_level(conf_level)
  ci_method <- match.arg(ci_method)
  if (is.null(fit) == is.null(loadings)) {
    .dmar_refuse_argument(
      "Supply exactly one of 'fit' or 'loadings'.", arg = c("fit", "loadings"))
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  if (!is.null(loadings)) {
    if (ci_method != "none") {
      .dmar_refuse_argument(c(
        "No interval can be constructed from 'loadings' alone; their ",
        "sampling variability is not carried by the numbers. Supply ",
        "the fitted model as 'fit' to bootstrap."), arg = c("loadings", "fit"))
    }
    if (!is.numeric(loadings) || length(loadings) < 2L ||
        anyNA(loadings) || any(abs(loadings) > 1)) {
      .dmar_refuse_argument(c(
        "'loadings' must be two or more standardized loadings in ",
        "[-1, 1]."), arg = "loadings")
    }
    out <- data.frame(factor = "f", ave = mean(loadings^2),
                      ci_lower = NA_real_, ci_upper = NA_real_,
                      stringsAsFactors = FALSE)
    return(.as_dmar_tbl(out, subclass = "dmar_ave"))
  }
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "Package 'lavaan' is required when supplying a fit. Install it ",
      "with install.packages(\"lavaan\")."),
      "package_required", package = "lavaan")
  }
  if (!inherits(fit, "lavaan")) {
    .dmar_refuse_argument("'fit' must be a lavaan fit object.", arg = "fit")
  }
  # A fit whose optimizer did not converge holds the point where it
  # stopped, not a solution (.lavaan_check_converged(),
  # R/no_maximum_internals.R). Through DMAR 1.0.1 and the first builds of
  # 1.1.0 its AVE was reported as any other: .412 for a one-factor fit
  # stopped after one iteration.
  .lavaan_check_converged(fit, "average_variance_extracted()")
  # One group and one level (.ave_check_one_block()). Through DMAR 1.0.1
  # and the first builds of 1.1.0 the loadings of every group, or of both
  # levels, of a factor were averaged into one unlabeled row: .450 for a
  # factor whose AVE is .810 in one group and .090 in the other.
  .ave_check_one_block(fit)
  ave_from <- function(f) {
    # Only est.std is used, so the delta method standard errors, z
    # statistics, p values, and intervals are switched off; that roughly
    # halves the cost of each bootstrap refit.
    .fiml_require_h1(f, "Two-stage AVE")
    std <- lavaan::standardizedSolution(f, se = FALSE, zstat = FALSE,
                                        pvalue = FALSE, ci = FALSE)
    lam <- std[std$op == "=~" & !.ave_zero_loadings(f, std),
               c("lhs", "est.std")]
    if (nrow(lam) == 0L) return(NULL)
    tapply(lam$est.std^2, lam$lhs, mean)
  }
  ave <- ave_from(fit)
  if (is.null(ave)) {
    .dmar_refuse_argument(
      "The fit contains no measurement (=~) loadings.", arg = NULL)
  }
  out <- data.frame(factor = names(ave), ave = as.numeric(ave),
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    stringsAsFactors = FALSE, row.names = NULL)
  # A factor with an indicator whose error variance is negative has a
  # communality above one among those it averages, so its AVE is not a
  # proportion of variance; it is NA, with the warning fa_confirmatory() gives for
  # its ave row. Through 1.0.1 it was reported (.726 for visual
  # perception and lozenges beside a verbal factor of
  # holzinger_swineford, whose measurement table in fa_confirmatory() says NA).
  improper <- .ave_improper(fit)
  improper <- improper[intersect(names(improper), out$factor)]
  for (f in names(improper)) {
    out$ave[out$factor == f] <- NA_real_
    .dmar_warn(.ave_improper_message(f, improper[[f]]), "heywood_warning")
  }

  if (ci_method == "percentile") {
    .fiml_check_boot_aux(fit, "AVE bootstrap")
    # The rule of ?dmar_bootstrap: 10 replications beyond each limit of
    # the percentile interval (400 at the 95 percent level).
    B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
    # A multiple-group or multilevel fit was refused above, before the
    # point estimate, so the bootstrap meets only fits of one group and
    # one level, whose rows lavaan's bootstrap resamples as independent
    # cases.
  }
  # No interval is formed for a factor that is improper in the observed
  # fit, so when no factor is proper the bootstrap is not run and no
  # confidence level is recorded, as krippendorff_alpha() forms no
  # interval for an undefined coefficient. Through 1.0.1 the call refit
  # the model B times and then reported the AVE of every factor and its
  # interval regardless; the first versions of 1.1.0 refit it B times for
  # nothing and printed "Confidence level: 95%" and a bootstrap line
  # beneath a table whose estimates and limits were all NA.
  formed <- ci_method == "percentile" &&
    length(setdiff(names(ave), names(improper))) > 0L
  if (formed) {
    .dmar_local_seed(seed)
    fac <- names(ave)
    k <- length(fac)
    # Each replication returns the AVE of every factor and then, for each
    # factor, whether the resample's solution is improper for it (an
    # indicator with a negative error variance; .ave_improper()). Such a
    # replication is kept, as the reliability family's improper = "keep"
    # keeps one, and counted. The check cannot fail a replication that
    # returned its AVE values (an error leaves its flags NA, counted as
    # proper), and it draws no random numbers, so the resamples and the
    # interval are those of 1.0.x.
    boot_fun <- function(f) {
      a <- tryCatch(ave_from(f), error = function(e) NULL)
      if (is.null(a)) return(rep(NA_real_, 2L * k))
      flags <- tryCatch(as.numeric(fac %in% names(.ave_improper(f))),
                        error = function(e) rep(NA_real_, k))
      c(as.numeric(a[fac]), flags)
    }
    boots <- try(suppressWarnings(
      lavaan::lavBootstrap(fit, r = as.integer(B), fun = boot_fun)),
      silent = TRUE)
    if (inherits(boots, "try-error")) {
      .dmar_refuse(c(
        "The bootstrap failed: ",
        conditionMessage(attr(boots, "condition")),
        " (the fit must carry its data; refit with the raw data ",
        "available)."), "computation_failed")
    }
    boots <- matrix(as.numeric(boots), ncol = 2L * k)
    flags <- boots[, k + seq_len(k), drop = FALSE]
    flags <- !is.na(flags) & flags == 1
    boots <- boots[, seq_len(k), drop = FALSE]
    ok <- rowSums(!is.finite(boots)) == 0L
    n_bad <- as.integer(B) - sum(ok)
    # The verdict of the rule (R/dmar_bootstrap.R): no interval when more
    # than half of the replications were set aside. The number the level
    # requires is held to the request, checked above, and not to the
    # replications that returned a value. Through the first versions of
    # 1.1.0 the floor on those returned was max(100, B / 2).
    verdict <- .dmar_boot_settle(
      B, sum(ok), tail = (1 - conf_level) / 2, conf_level = conf_level,
      fun = "average_variance_extracted", quantity = "AVE for some factor",
      dropped_because = "gave a refit that failed or did not converge",
      signal = FALSE)
    enough <- verdict$formed
    boots <- boots[ok, , drop = FALSE]
    if (enough) {
      alpha_2 <- (1 - conf_level) / 2
      out$ci_lower <- apply(boots, 2, stats::quantile, probs = alpha_2, names = FALSE)
      out$ci_upper <- apply(boots, 2, stats::quantile, probs = 1 - alpha_2, names = FALSE)
    } else {
      out$ci_lower <- out$ci_upper <- NA_real_
    }
    # No interval for an AVE that is not a proportion of variance.
    out[out$factor %in% names(improper), c("ci_lower", "ci_upper")] <-
      NA_real_
    bootstrap <- .ave_bootstrap_counts(flags[ok, , drop = FALSE], fac,
                                       if (enough) setdiff(fac, names(improper)) else character(0),
                                       B, conf_level)
    bootstrap$require_half <- TRUE
    bootstrap$minimum_returned <- as.integer(max(2L, ceiling(B / 2)))
    bootstrap$interval_withheld <- !enough
    bootstrap$required <- verdict$B_required
    bootstrap$bootstrap_status <- verdict$status
    bootstrap$B_used <- verdict$B_used
    # One warning for the bootstrap, as the reliability family gives: the
    # replications set aside, in the words of the line under the table,
    # and the caution, naming each factor it holds for
    # (.ave_bootstrap_warning()). Through 1.0.1 the warning said only how
    # many replications "failed to refit (or did not converge) and were
    # dropped", and the improper replications kept were neither counted
    # nor reported; the first versions of 1.1.0 printed the caution under
    # the table only.
    # Its class follows the verdict: dmar_bootstrap_insufficient when the
    # interval is withheld, dmar_bootstrap_replications_dropped when
    # replications were set aside and it was formed.
    note <- .ave_bootstrap_warning(bootstrap, conf_level)
    .dmar_boot_signal_text(verdict, note)
  }
  out <- .as_dmar_tbl(out, subclass = "dmar_ave")
  if (formed) {
    attr(out, "conf_level") <- conf_level
    out <- .dmar_boot_attrs(out, verdict)
    attr(out, "bootstrap") <- bootstrap
  }
  out
}

# The `bootstrap` attribute of a percentile interval: the counts of the
# reliability family's attribute (.bootstrap_ci() in
# R/reliability_internals.R), so that its footer line
# (.bootstrap_footer_line()) reports them, and the improper replications
# kept for each factor.
#
# flags      : the replications that returned a value (rows) by factor
#              (columns, in the order of `fac`), TRUE where the resample's
#              solution was improper for the factor.
# fac        : every factor.
# with_ci    : the factors whose interval was formed, those proper in the
#              observed fit. All defined improper fits remain counted;
#              only the interval caution is restricted to these factors.
# B          : the replications attempted.
# conf_level : the confidence level, for the caution.
#
# lavaan's bootstrap returns NA for a replication whose refit fails or
# does not converge, without saying which, so the replications set aside
# are counted under reasons[["non_finite"]], and none under `error` or
# `not_converged`; `refit_failed`, the same count, names them for both
# causes in the line under the table and in the warning
# (.bootstrap_reason_phrase()), where the line once said they "returned a
# non-finite value". No replication is set aside for being improper
# (`improper` is 0 and `improper_rule` "keep"): `improper_kept` counts
# those used whose solution was improper for at least one factor, and `improper_kept_by_factor` (named by factor, one element
# per factor) counts them for each; `factors_with_intervals` records eligibility. The caution of the
# family, that the replications in question reach the tail area beyond a
# limit (.boot_caution()), is TRUE when it holds for the interval of some
# factor: those set aside and the improper ones kept for that factor.
.ave_bootstrap_counts <- function(flags, fac, with_ci, B, conf_level) {
  B <- as.integer(B)
  returned <- nrow(flags)
  failed <- B - returned
  cols <- match(with_ci, fac)
  by_factor <- stats::setNames(as.integer(colSums(flags)), fac)
  kept <- as.integer(sum(rowSums(flags) > 0))
  caution <- any(vapply(by_factor[cols], function(n)
    .boot_caution(failed + n, B, conf_level), logical(1L)))
  list(attempted = B,
       returned = as.integer(returned),
       failed = as.integer(failed),
       reasons = c(error = 0L, non_finite = as.integer(failed)),
       error_messages = integer(),
       undetermined = 0L,
       improper = 0L,
       improper_kept = kept,
       no_solution = 0L,
       not_converged = 0L,
       refit_failed = as.integer(failed),
       improper_rule = "keep",
       caution = caution,
       improper_kept_by_factor = by_factor,
       factors_with_intervals = with_ci)
}

# The warning a percentile interval gives (NULL when there is nothing to
# say), as .bootstrap_ci() in R/reliability_internals.R gives one for the
# reliability family: one warning per call, its notes joined.
#
#   - The replications set aside, in the words of the family's warning,
#     with the reason the line under the table gives
#     (.bootstrap_reason_phrase(), whose `refit_failed` names both causes
#     lavaan's bootstrap cannot tell apart), and the improper ones kept.
#   - The caution, for each factor with an interval for which the
#     replications set aside and the improper ones kept for it reach the
#     tail area beyond a limit (.boot_caution()), naming the factor and
#     its count: its limits then depend on how those replications were
#     treated.
#
# bs         : the `bootstrap` attribute (.ave_bootstrap_counts()).
# conf_level : the confidence level.
.ave_bootstrap_warning <- function(bs, conf_level) {
  by_factor <- bs$improper_kept_by_factor
  if (!is.null(bs$factors_with_intervals)) {
    by_factor <- by_factor[intersect(names(by_factor), bs$factors_with_intervals)]
  }
  n_ci <- length(by_factor)
  notes <- character()
  if (isTRUE(bs$interval_withheld)) {
    notes <- c(notes, sprintf(paste0(
      "Only %d of %d replications returned estimates; at least %d are ",
      "required (half of those requested; see ?dmar_bootstrap), so the ",
      "interval is withheld"),
      bs$returned, bs$attempted, bs$minimum_returned))
  }
  if (bs$failed > 0L) {
    notes <- c(notes, sprintf(
      "%d of %d bootstrap replications %s set aside (%s)%s",
      bs$failed, bs$attempted, if (bs$failed == 1L) "was" else "were",
      .bootstrap_reason_phrase(bs$reasons, bs$error_messages,
                               refit_failed = bs$refit_failed),
      if (n_ci == 0L) {
        ""
      } else {
        sprintf("; the %s computed from the %d that returned a value%s%s",
                if (n_ci == 1L) "interval is" else "intervals are",
                bs$returned,
                if (bs$improper_kept > 0L) {
                  sprintf(paste0(", %d of them from an improper (Heywood) ",
                                 "solution and kept"), bs$improper_kept)
                } else {
                  ""
                },
                # The count is held to the request; the note says so and
                # says to request more where replications can fail.
                sub("\\.$", "", sub(
                  "^ The", "; the",
                  .dmar_boot_request_note(bs$attempted, bs$returned,
                                          bs$required, "replications"))))
      }))
  }
  in_question <- bs$failed + by_factor
  cautious <- names(by_factor)[vapply(in_question, function(n)
    .boot_caution(n, bs$attempted, conf_level), logical(1L))]
  if (length(cautious)) {
    n <- in_question[cautious]
    # Each share is formatted on its own: .pct() on the vector would give
    # them a common number of decimals ("10.0%" beside "7.5%").
    each <- sprintf("for %s, %d of %d%s (%s)", cautious, n, bs$attempted,
                    c(" replications", rep("", length(n) - 1L)),
                    vapply(n / bs$attempted, .pct, character(1L)))
    one <- length(each) == 1L
    notes <- c(notes, sprintf(paste0(
      "Caution: %s came from a solution improper for the factor or were ",
      "set aside, at least the %s of the bootstrap distribution that lies ",
      "beyond each limit of a %s%% interval, so the limits of %s depend on ",
      "how those replications were treated; report %s with these counts"),
      if (one) {
        each
      } else {
        paste0(paste(each[-length(each)], collapse = ", "),
               if (length(each) > 2L) "," else "", " and ",
               each[length(each)])
      },
      .pct((1 - conf_level) / 2), .format_conf_pct(conf_level),
      if (one) "its interval" else "their intervals",
      if (one) "the interval" else "the intervals"))
  }
  if (!length(notes)) return(NULL)
  paste0(paste(notes, collapse = ". "), ".")
}

# The line a table of more than one factor with an interval adds under the
# family's bootstrap line (.bootstrap_footer_line()), when some replication
# was improper for one of them: which factor's interval the improper
# replications kept belong to. NULL otherwise.
.ave_bootstrap_footer <- function(x) {
  bs <- attr(x, "bootstrap", exact = TRUE)
  by_factor <- if (is.list(bs)) bs$improper_kept_by_factor
  if (length(by_factor) < 2L || !any(by_factor > 0L)) return(NULL)
  paste0("Improper replications kept, by factor (a resample in which an ",
         "indicator of the factor had a negative error variance, so that ",
         "its AVE averages a communality above one): ",
         paste(names(by_factor), by_factor, collapse = ", "))
}

# The lines under the table: those every table carries, the bootstrap
# line among them, and then the improper replications kept for each
# factor (.ave_bootstrap_footer()) in the same block.
#' @exportS3Method .dmar_tbl_footer dmar_ave
#' @noRd
.dmar_tbl_footer.dmar_ave <- function(x, ...) {
  .dmar_footer_add(NextMethod(), .ave_bootstrap_footer(x), "bootstrap")
}

# The refusal of a fit of more than one group or more than one level. A
# factor's standardized loadings are those of one group at one level:
# each group has its own, standardized by that group's variances, and in
# a two-level model the within and between parts of the items have their
# own, whatever the factors are named. A factor therefore has an AVE in
# each group (or at each level). Their unweighted mean is not the AVE of
# any group, a mean weighted by group size would change with the group
# sizes when the groups did not, and a table with a row per group or
# level would take a column that a fit of one group does not have. The
# function describes a fit of one group and one level, and refuses the
# rest with class dmar_argument_invalid, naming the groups (as `groups`)
# or the levels (as `levels`) and the two routes to the AVE of one of
# them: that group's or level's standardized loadings given as
# `loadings`, or a fit to that group's data alone. Through DMAR 1.0.1 and
# the first builds of 1.1.0 the loadings of all groups, or of both
# levels, were averaged into one unlabeled row (.450 for a factor whose
# AVE is .810 in one group and .090 in the other; .591 for a factor whose
# AVE is .356 within clusters and .826 between them).
.ave_check_one_block <- function(fit) {
  n_groups <- lavaan::lavInspect(fit, "ngroups")
  if (n_groups > 1L) {
    labels <- as.character(lavaan::lavInspect(fit, "group.label"))
    .dmar_refuse_argument(c(
      "'fit' is a multiple-group fit (", n_groups, " groups: ",
      paste(labels, collapse = ", "), "). Each group has its own ",
      "standardized loadings, so a factor has an AVE in each group, and ",
      "neither their unweighted mean nor a mean weighted by group size ",
      "is the AVE of any group; average_variance_extracted() describes a ",
      "fit of one group. For one group of this model, give the ",
      "standardized loadings of the factor in that group (the est.std ",
      "of its =~ rows for the group in lavaan::standardizedSolution(fit)) ",
      "as 'loadings', or fit the model to that group's data alone (which ",
      "drops any equality constraints across groups)."),
      arg = "fit", groups = labels)
  }
  n_levels <- lavaan::lavInspect(fit, "nlevels")
  if (n_levels > 1L) {
    levels <- tryCatch(as.character(lavaan::lavInspect(fit, "level.label")),
                       error = function(e) as.character(seq_len(n_levels)))
    .dmar_refuse_argument(c(
      "'fit' is a multilevel fit (", n_levels, " levels: ",
      paste(levels, collapse = ", "), "). The items' variation at each ",
      "level has its own standardized loadings, even for a factor of the ",
      "same name at both, so a factor has an AVE at each level, and their ",
      "mean is the AVE of neither; average_variance_extracted() describes ",
      "a fit of one level. For one level of this model, give the ",
      "standardized loadings of the factor at that level as 'loadings': ",
      "the factor's column of that level's lambda matrix in ",
      "lavaan::lavInspect(fit, \"std\"), without the zeros of items that ",
      "do not load on it."),
      arg = "fit", levels = levels)
  }
  invisible(NULL)
}

# Which rows of `tab`, a parameterEstimates() or standardizedSolution()
# table of `fit`, are loadings the model holds at zero (.pt_zero_rows()
# in R/reliability_internals.R: fixed at zero, 0*x, or constrained to
# zero, a*x with a == 0). Such a row is not a loading: the model written
# with it is the model written without it, as factor_scores(),
# bifactor_indices(), and reliability_H() read it, so its item is not an
# indicator of the factor, in the AVE or in the check for an improper
# solution. Through 1.0.1 the AVE averaged it in as a loading of zero.
.ave_zero_loadings <- function(fit, tab) {
  pt <- lavaan::parTable(fit)
  group_of <- function(x) if (is.null(x$group)) rep(1L, nrow(x)) else x$group
  zero <- pt$op == "=~" & .pt_zero_rows(pt)
  paste(tab$lhs, tab$op, tab$rhs, group_of(tab)) %in%
    paste(pt$lhs, pt$op, pt$rhs, group_of(pt))[zero]
}

# The indicators of each factor whose error variance is estimated below
# zero in a lavaan fit (for an ordered item, the residual variance of its
# latent response, which the delta parameterization reports as a derived
# value and the theta parameterization estimates; for a factor of a
# higher order model, the disturbance variance of a first order factor),
# with their communalities, one minus the standardized error variance,
# which a negative error variance puts above one. The fits it meets have
# one group and one level (.ave_check_one_block() refuses the others
# before it is called); the group key of the matching below is kept so
# that a row is never matched with another group's. Returns a list named
# by the factors that have such indicators, each element a numeric vector
# of the communalities named by indicator; an empty list for a proper
# solution. The counterpart of .cfa_k_improper() in R/fa_confirmatory.R for any
# fit, whatever its identification.
.ave_improper <- function(fit) {
  pe <- lavaan::parameterEstimates(fit, se = FALSE, zstat = FALSE,
                                   pvalue = FALSE, ci = FALSE)
  std <- lavaan::standardizedSolution(fit, se = FALSE, zstat = FALSE,
                                      pvalue = FALSE, ci = FALSE)
  group_of <- function(x) if (is.null(x$group)) rep(1L, nrow(x)) else x$group
  pe$group <- group_of(pe)
  std$group <- group_of(std)
  load <- pe[pe$op == "=~" & !.ave_zero_loadings(fit, pe),
             c("lhs", "rhs", "group"), drop = FALSE]
  vars <- pe[pe$op == "~~" & pe$lhs == pe$rhs, c("lhs", "group", "est"),
             drop = FALSE]
  vars_std <- std[std$op == "~~" & std$lhs == std$rhs,
                  c("lhs", "group", "est.std"), drop = FALSE]
  out <- list()
  for (f in unique(load$lhs)) {
    rows <- load[load$lhs == f, , drop = FALSE]
    key <- paste(rows$rhs, rows$group)
    psi <- vars$est[match(key, paste(vars$lhs, vars$group))]
    bad <- !is.na(psi) & psi < 0
    if (!any(bad)) next
    resid_std <- vars_std$est.std[match(key[bad], paste(vars_std$lhs,
                                                        vars_std$group))]
    communality <- tapply(1 - resid_std, rows$rhs[bad], max)
    communality <- communality[unique(rows$rhs[bad])]
    out[[f]] <- structure(as.numeric(communality),
                          names = names(communality),
                          kind = .ave_indicator_kind(fit,
                                                     names(communality)))
  }
  out
}

# What each indicator is, for the wording of the warning: "ordered" (an
# ordered item, whose residual variance is that of its latent response),
# "observed" (a continuous item, with an error variance), or "latent" (a
# factor indicating a higher order factor, with a residual variance).
.ave_indicator_kind <- function(fit, indicators) {
  ord <- lavaan::lavNames(fit, "ov.ord")
  ov <- lavaan::lavNames(fit, "ov")
  ifelse(indicators %in% ord, "ordered",
         ifelse(indicators %in% ov, "observed", "latent"))
}

# The warning average_variance_extracted() gives for a factor whose
# solution is improper (see .ave_improper()): the factor, each indicator
# with its communality, and why the AVE is undefined. Worded as fa_confirmatory()
# words its warning (.kfit_improper_message() in
# R/reliability_internals.R); the communalities are formatted for the
# message only.
.ave_improper_message <- function(factor, communality) {
  items <- names(communality)
  vals <- formatC(as.numeric(communality), format = "f", digits = 3)
  kind <- attr(communality, "kind") %||% rep("observed", length(items))
  one <- length(items) == 1L
  which <- if (all(kind == "ordered") && one) {
    paste0("the latent response behind ", items, " has a negative ",
           "residual variance (its communality is ", vals, ")")
  } else if (all(kind == "ordered")) {
    paste0("the latent responses behind ", paste(items, collapse = " and "),
           " have negative residual variances (their communalities are ",
           paste(vals, collapse = " and "), ")")
  } else if (one) {
    paste0(items, " has a negative ",
           if (kind == "observed") "error" else "residual",
           " variance (its communality is ", vals, ")")
  } else {
    paste0(paste(items, collapse = " and "), " have negative ",
           if (all(kind == "observed")) "error" else "error or residual",
           " variances (their communalities are ",
           paste(vals, collapse = " and "), ")")
  }
  paste0("The AVE of ", factor, " is NA: the solution is improper (a ",
         "Heywood case): ", which, ". The average variance extracted ",
         "averages the communalities, so it is not a proportion of ",
         "variance here. A larger sample, or reconsidering the model for ",
         "these items, are the usual remedies.")
}
