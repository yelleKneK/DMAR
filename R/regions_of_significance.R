# Regions of significance for a covariate-by-group interaction.
#
# When the within-group regression slopes differ, the group difference is a
# function of the covariate rather than a single number. Setting the squared
# difference equal to the squared critical value times its sampling variance
# gives a quadratic in the covariate whose real roots bound the covariate
# values at which the groups differ significantly.

# Internal: solve A x^2 + B x + C = 0 for the boundaries of the region of
# significance, where f(x) = D(x)^2 - tcrit^2 * Var(D(x)) and the groups
# differ significantly wherever f(x) > 0. Returns the sorted boundaries
# (NA when a boundary does not exist), the number of boundaries, and a
# numeric region code, documented in regions_of_significance().
#
#   1 = two boundaries, significant outside  [lower_bound, upper_bound]
#   2 = two boundaries, significant between  (lower_bound, upper_bound)
#   3 = no boundary, significant at every covariate value
#   4 = no boundary, significant at no covariate value
#   5 = one boundary, significant above it   (degenerate, linear case)
#   6 = one boundary, significant below it   (degenerate, linear case)
#
# The sign of A is the whole geometry: A = d1^2 - tcrit^2 * Var(d1) is
# positive exactly when the slope difference itself clears the critical
# value, and a positive leading coefficient makes the parabola open upward,
# so significance lies outside the two boundaries rather than between them.
.regions_solve <- function(d0, d1, v00, v01, v11, tcrit) {
  A <- d1^2 - tcrit^2 * v11
  B <- 2 * (d0 * d1 - tcrit^2 * v01)
  C <- d0^2 - tcrit^2 * v00

  # Degenerate case: the quadratic collapses to a line (or to a constant)
  # because the slope difference sits exactly on the critical value.
  if (A == 0) {
    if (B == 0) {
      return(list(lower = NA_real_, upper = NA_real_, n = 0L,
                  code = if (C > 0) 3L else 4L))
    }
    root <- -C / B
    return(list(lower = root, upper = NA_real_, n = 1L,
                code = if (B > 0) 5L else 6L))
  }

  disc <- B^2 - 4 * A * C
  if (disc <= 0) {
    # No crossing: f(x) keeps the sign of A over the whole real line.
    return(list(lower = NA_real_, upper = NA_real_, n = 0L,
                code = if (A > 0) 3L else 4L))
  }

  root_span <- sqrt(disc)
  roots <- sort(c((-B - root_span) / (2 * A), (-B + root_span) / (2 * A)))
  list(lower = roots[1L], upper = roots[2L], n = 2L,
       code = if (A > 0) 1L else 2L)
}

# Internal: a compact number for the human-readable region string.
.regions_fmt <- function(v) format(signif(v, 4L), trim = TRUE)

# Internal: the plain-language reading of a region code.
.regions_describe <- function(code, lower, upper, covariate) {
  switch(as.character(code),
    "1" = sprintf("The groups differ for %s < %s or %s > %s",
                  covariate, .regions_fmt(lower),
                  covariate, .regions_fmt(upper)),
    "2" = sprintf("The groups differ for %s < %s < %s",
                  .regions_fmt(lower), covariate, .regions_fmt(upper)),
    "3" = sprintf("The groups differ at every value of %s", covariate),
    "4" = sprintf("The groups do not differ at any value of %s", covariate),
    "5" = sprintf("The groups differ for %s > %s",
                  covariate, .regions_fmt(lower)),
    "6" = sprintf("The groups differ for %s < %s",
                  covariate, .regions_fmt(lower)),
    NA_character_
  )
}

#' Regions of Significance for a Covariate by Group Interaction
#'
#' Finds the values of a covariate at which two groups differ
#' significantly when the within-group regression slopes are \emph{not}
#' equal, that is, when there is a covariate-by-group interaction
#' (heterogeneity of regression). With more than two groups the
#' calculation is carried out for every pair of groups.
#'
#' @param object Either a fitted \code{lm} or \code{aov} of the form
#'   \code{y ~ x * group}, or a model formula of that form. When a
#'   formula is supplied, \code{data} must be supplied too and the model
#'   is fit with \code{\link[stats]{lm}}.
#' @param data A \code{data.frame} holding the variables in the formula.
#'   Used, and required, only when \code{object} is a formula; ignored
#'   with a warning when \code{object} is already a fitted model.
#' @param conf_level Confidence level for the boundaries. Default
#'   \code{0.95}.
#' @param method Character string naming the critical value.
#'   \code{"simultaneous"} (the default) uses Potthoff's (1964)
#'   simultaneous critical value \eqn{\sqrt{2 F_{1 - \alpha; 2,
#'   \nu}}}, which holds the error rate over the whole covariate
#'   continuum at once. \code{"pointwise"} uses the classic critical
#'   value \eqn{t_{1 - \alpha/2; \nu}}, which holds it at one
#'   covariate value chosen in advance.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with one row per
#'   reported quantity per pair of groups and columns \code{pair},
#'   \code{term}, and a numeric \code{value}. The terms, for each pair,
#'   are
#'   \describe{
#'     \item{\code{lower_bound}, \code{upper_bound}}{The boundaries of
#'       the region, sorted, on the scale of the covariate.
#'       \code{NA} when a boundary does not exist; when there is a
#'       single boundary it is reported as \code{lower_bound} and
#'       \code{upper_bound} is \code{NA}. Read them with
#'       \code{region_code}: the boundaries are the covariate values at
#'       which the group difference sits exactly on the critical value,
#'       and it is \code{region_code} that says on which side of them
#'       the groups differ.}
#'     \item{\code{region_code}}{How to read the boundaries. \code{1}:
#'       two boundaries, the groups differ \emph{outside} them.
#'       \code{2}: two boundaries, the groups differ \emph{between}
#'       them. \code{3}: no boundary, the groups differ at every
#'       covariate value. \code{4}: no boundary, the groups differ at no
#'       covariate value. \code{5} and \code{6}: one boundary, the
#'       groups differ above it (\code{5}) or below it (\code{6}).}
#'     \item{\code{n_boundaries}}{How many real boundaries exist: 0, 1,
#'       or 2.}
#'     \item{\code{lower_bound_in_range}, \code{upper_bound_in_range}}{\code{1}
#'       when the boundary falls inside the covariate values actually
#'       observed in the two groups, \code{0} when it falls outside them,
#'       \code{NA} when the boundary does not exist. A boundary outside
#'       the observed range is an extrapolation of the fitted lines and
#'       should not be interpreted as a covariate value at which anything
#'       was or could be seen.}
#'     \item{\code{difference_intercept}, \code{difference_slope}}{The
#'       intercept \eqn{d_0} and slope \eqn{d_1} of the group difference
#'       as a function of the covariate, so that the estimated difference
#'       at covariate value \eqn{x} is \eqn{d_0 + d_1 x}.}
#'     \item{\code{critical_value}}{The critical value used, \eqn{t_{crit}}.}
#'     \item{\code{df_error}}{Error degrees of freedom of the fitted model.}
#'     \item{\code{conf_level}}{The confidence level.}
#'   }
#'   Everything that is not a number is carried on attributes rather
#'   than forced into \code{value}: \code{method}, \code{outcome},
#'   \code{covariate}, \code{group}, the overall observed
#'   \code{covariate_range}, a \code{pairs} \code{data.frame} with the
#'   group labels, the per-pair observed covariate range, and a
#'   plain-language \code{region} string for each pair, and a
#'   \code{geometry} \code{data.frame} holding \eqn{d_0}, \eqn{d_1}, and
#'   the three elements of their covariance matrix, which is what
#'   \code{\link{plot_regions_of_significance}} draws.
#'
#' @details
#' \strong{Why the procedure exists.} An analysis of covariance that
#' assumes a common within-group slope reports one adjusted mean
#' difference, and that single number is a complete summary of the group
#' comparison only if the slopes really are common. When the covariate
#' interacts with the group factor the slopes are not common, the two
#' fitted lines converge or cross, and there is no such thing as
#' \dQuote{the} treatment effect: the difference between the groups
#' depends on where along the covariate you look. Reporting the adjusted
#' mean difference anyway reports the difference at one covariate value,
#' the covariate grand mean, and says nothing about the rest of the
#' range. The question worth answering is instead \emph{where} on the
#' covariate the groups differ, and that is what this function answers.
#'
#' \strong{The calculation.} For two groups, write the estimated
#' difference at covariate value \eqn{x} as the line
#' \deqn{\hat D(x) = \hat d_0 + \hat d_1 x,}
#' where \eqn{\hat d_0} is the difference in intercepts and \eqn{\hat
#' d_1} the difference in slopes. Because \eqn{\hat D(x)} is a linear
#' combination of the regression coefficients, its sampling variance
#' follows from their covariance matrix,
#' \deqn{\mathrm{Var}[\hat D(x)] = \mathrm{Var}(\hat d_0) +
#'   2 x \, \mathrm{Cov}(\hat d_0, \hat d_1) + x^2 \mathrm{Var}(\hat d_1).}
#' The groups differ significantly at \eqn{x} exactly when \eqn{\hat
#' D(x)^2 > t_{crit}^2 \, \mathrm{Var}[\hat D(x)]}. Setting the two sides
#' equal gives a quadratic in \eqn{x},
#' \deqn{(\hat d_1^2 - t_{crit}^2 \mathrm{Var}(\hat d_1)) x^2 +
#'   2(\hat d_0 \hat d_1 - t_{crit}^2 \mathrm{Cov}(\hat d_0, \hat d_1)) x +
#'   (\hat d_0^2 - t_{crit}^2 \mathrm{Var}(\hat d_0)) = 0,}
#' whose real roots are the boundaries of the region of significance.
#'
#' \strong{Every geometry is possible, and the leading coefficient
#' decides which.} The coefficient on \eqn{x^2} is positive exactly when
#' the slope difference itself clears the critical value. When it is
#' positive the parabola opens upward and the groups differ
#' \emph{outside} the two boundaries, the familiar picture of two lines
#' that cross somewhere in the middle of the covariate and separate at
#' both ends. When it is negative the parabola opens downward and the
#' groups differ \emph{between} the boundaries, a middle band of
#' covariate values where the two lines are far enough apart relative to
#' the precision available there. When there are no real roots the sign
#' never changes, so the groups differ either everywhere or nowhere. All
#' of these are reported through \code{region_code} rather than being
#' treated as failures.
#'
#' Of those, \dQuote{everywhere} is a case the code enumerates but the
#' mathematics rules out whenever the slopes genuinely differ. The
#' estimated difference \eqn{\hat D(x)} is then a nonconstant line in
#' \eqn{x}, so it crosses zero at some covariate value, and at that value
#' the squared difference is zero while the critical bound is positive.
#' The quadratic therefore always has two real roots, and the significant
#' set is the pair of tails outside them or the band between them, never
#' the whole line. Note this is a statement about the covariate axis
#' extended without limit, not about the observed data: within the range
#' actually observed the groups may well differ everywhere, which is why
#' the next paragraph matters.
#'
#' \strong{Boundaries outside the observed data.} A boundary is a root of
#' an equation, and the equation is happy to place it far outside the
#' covariate values that were actually observed. Such a boundary is an
#' extrapolation of two fitted lines into a region where there are no
#' data to support them, and it should not be read as a covariate value
#' at which anything can be claimed. The boundary is still reported,
#' since suppressing it would hide the shape of the result, but it is
#' flagged in \code{lower_bound_in_range} and
#' \code{upper_bound_in_range} and noted in the \code{region} string.
#'
#' \strong{Simultaneous versus pointwise.} The classic critical value is
#' \eqn{t_{1 - \alpha/2}}, which controls the Type I error rate at a
#' \emph{single} covariate value fixed in advance. That is not what
#' anyone actually does: the whole point of the procedure is to scan the
#' covariate continuum and read off where the groups differ, which is a
#' search over infinitely many tests. Potthoff (1964) gave the
#' simultaneous critical value \eqn{\sqrt{2 F_{1 - \alpha; 2, \nu}}},
#' which holds the error rate over the entire covariate range at once
#' and is therefore the default here, and the value used in Maxwell,
#' Delaney, and Kelley (2027, Chapter 9). The simultaneous critical
#' value is always the larger of the two, so its region of significance
#' is always the more conservative one. With more than two groups the
#' simultaneous guarantee applies to each pair over the covariate range,
#' not to the family of pairs. For a Bonferroni protection across the
#' pairs, divide the Type I error rate by the number of pairs before
#' choosing \code{conf_level}: with three groups (three pairs) and a
#' familywise rate of .05, pass \code{conf_level = 1 - .05 / 3}.
#'
#' \strong{What the model may contain.} The model must contain exactly
#' one interaction between a numeric covariate and a grouping factor,
#' and the grouping factor may not appear in any other term. Additional
#' predictors that do not interact with the group are allowed and drop
#' out of the group difference, so \code{y ~ x * group + block} is fine
#' while \code{y ~ x * group * block} is not. A grouping variable stored
#' as a number (0/1, say) is a numeric predictor to R, not a factor, so
#' convert it with \code{\link[base]{factor}} first.
#'
#' @references
#' Johnson, P. O., & Neyman, J. (1936). Tests of certain linear
#'   hypotheses and their application to some educational problems.
#'   \emph{Statistical Research Memoirs, 1}, 57--93.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 9 and its extension
#'   on heterogeneity of regression.)
#'
#' Potthoff, R. F. (1964). On the Johnson-Neyman technique and some
#'   extensions thereof. \emph{Psychometrika, 29}(3), 241--256.
#'   \doi{10.1007/BF02289721}
#'
#' Rogosa, D. (1980). Comparing nonparallel regression lines.
#'   \emph{Psychological Bulletin, 88}(2), 307--321.
#'   \doi{10.1037/0033-2909.88.2.307}
#'
#' @seealso \code{\link{plot_regions_of_significance}} to see the group
#'   difference and its confidence band across the covariate;
#'   \code{\link{ancova}} for the common-slope analysis and its
#'   homogeneity-of-regression test; \code{\link{pygmalion}} for the
#'   data used below.
#'
#' @examples
#' # ---- The Pygmalion teacher-expectancy data ----
#' # Post-test IQ, averaged over the two follow-up assessments, on
#' # pretest IQ, separately by condition. The slopes differ, so the
#' # expectancy effect depends on where the child started.
#' data(pygmalion)
#' pygmalion$iq_post <- (pygmalion$iq_4 + pygmalion$iq_8) / 2
#' fit <- lm(iq_post ~ iq_pre * treatment, data = pygmalion)
#'
#' regions_of_significance(fit)
#'
#' # The plain-language reading of each pair is on an attribute.
#' attr(regions_of_significance(fit), "pairs")$region
#'
#' # The pointwise (classic) critical value gives a wider region,
#' # because it does not pay for scanning the whole covariate.
#' regions_of_significance(fit, method = "pointwise")
#'
#' # ---- Formula interface, and more than two groups ----
#' set.seed(113)
#' n <- 150
#' g <- factor(rep(c("control", "low", "high"), each = n / 3))
#' x <- rnorm(n, 50, 10)
#' y <- 2 + 0.5 * x + (g == "high") * (0.4 * x - 15) + rnorm(n, 0, 5)
#' d <- data.frame(y, x, g)
#' regions_of_significance(y ~ x * g, data = d)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest models
#'
#' @family hypothesis tests
#'
#' @importFrom stats qf qt
#'
#' @export

regions_of_significance <- function(object, data = NULL, conf_level = 0.95,
                                    method = c("simultaneous", "pointwise")) {

  method <- match.arg(method)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be a single number strictly between 0 and 1.",
         call. = FALSE)

  # ---------- Fitted model, or formula plus data ----------
  if (inherits(object, c("lm", "aov"))) {
    if (!is.null(data))
      warning("'object' is already a fitted model, so 'data' is ignored.",
              call. = FALSE)
    fit <- object
  } else if (inherits(object, "formula")) {
    if (is.null(data) || !is.data.frame(data))
      stop("When 'object' is a formula, 'data' must be a data.frame ",
           "containing the variables in it.", call. = FALSE)
    fit <- stats::lm(object, data = data)
  } else {
    stop("'object' must be a fitted lm or aov of the form y ~ x * group, ",
         "or a formula of that form with 'data' supplied.", call. = FALSE)
  }

  beta <- stats::coef(fit)
  if (anyNA(beta))
    stop("The fitted model is rank deficient (some coefficients are NA), ",
         "so the group difference is not estimable. Check for empty ",
         "group-by-covariate combinations or collinear predictors.",
         call. = FALSE)
  vcov_beta <- stats::vcov(fit)
  df_error  <- stats::df.residual(fit)
  if (!is.finite(df_error) || df_error < 1)
    stop("The fitted model has no error degrees of freedom.", call. = FALSE)

  # ---------- Locate the covariate-by-group interaction ----------
  model_terms  <- stats::terms(fit)
  term_labels  <- attr(model_terms, "term.labels")
  term_factors <- attr(model_terms, "factors")
  term_order   <- attr(model_terms, "order")
  mf <- stats::model.frame(fit)

  is_grouping <- function(v) {
    is.factor(mf[[v]]) || is.character(mf[[v]]) || is.logical(mf[[v]])
  }
  vars_in <- function(lab) {
    rownames(term_factors)[term_factors[, lab] > 0]
  }

  two_way <- term_labels[term_order == 2L]
  candidates <- character(0)
  cov_of <- grp_of <- character(0)
  for (lab in two_way) {
    vs <- vars_in(lab)
    if (length(vs) != 2L) next
    grouping <- vapply(vs, is_grouping, logical(1L))
    if (sum(grouping) == 1L && is.numeric(mf[[vs[!grouping]]])) {
      candidates <- c(candidates, lab)
      cov_of  <- c(cov_of,  vs[!grouping])
      grp_of  <- c(grp_of,  vs[grouping])
    }
  }

  if (length(candidates) == 0L)
    stop("No covariate-by-group interaction was found in the model. ",
         "Regions of significance describe how a group difference changes ",
         "with a covariate, so the model must be of the form y ~ x * group ",
         "with a numeric x and a factor group. A grouping variable stored ",
         "as a number is a numeric predictor to R; convert it with ",
         "factor() first.", call. = FALSE)
  if (length(candidates) > 1L)
    stop("The model contains more than one numeric-by-factor interaction (",
         paste(candidates, collapse = ", "), "). Fit a model with a single ",
         "covariate-by-group interaction.", call. = FALSE)

  interaction_label <- candidates
  cov_name <- cov_of
  grp_name <- grp_of

  # Every other appearance of the grouping factor would leave a term in the
  # group difference that the two-point construction below cannot see, so
  # refuse rather than silently report the wrong line.
  grp_terms <- term_labels[term_factors[grp_name, term_labels] > 0]
  extra <- setdiff(grp_terms, c(grp_name, interaction_label))
  if (length(extra))
    stop("The grouping factor '", grp_name, "' also appears in the term(s) ",
         paste0("'", extra, "'", collapse = ", "), ". Regions of ",
         "significance are defined for a single covariate-by-group ",
         "interaction; other predictors may enter the model but may not ",
         "interact with the group.", call. = FALSE)

  group_column <- mf[[grp_name]]
  levs <- fit$xlevels[[grp_name]]
  if (is.null(levs)) {
    levs <- if (is.factor(group_column)) levels(group_column) else
      sort(unique(as.character(group_column)))
  }
  n_levels <- length(levs)
  if (n_levels < 2L)
    stop("The grouping factor '", grp_name, "' has fewer than two levels.",
         call. = FALSE)

  # ---------- Difference lines, one per pair of groups ----------
  # Evaluate the model matrix at covariate 0 and covariate 1 for each level,
  # holding every other predictor fixed. The row at 0 gives the coefficient
  # weights of a group's intercept and the difference of the two rows gives
  # the weights of its slope, so any pair of groups yields the exact linear
  # combinations behind d_0 and d_1. Predictors that do not interact with the
  # group take the same value in both groups and cancel.
  grid <- mf[rep(1L, 2L * n_levels), , drop = FALSE]
  grid[[cov_name]] <- rep(c(0, 1), times = n_levels)
  level_seq <- rep(levs, each = 2L)
  grid[[grp_name]] <-
    if (is.factor(group_column)) factor(level_seq, levels = levs)
    else if (is.logical(group_column)) as.logical(level_seq)
    else level_seq
  rownames(grid) <- NULL

  mm <- stats::model.matrix(model_terms, grid, contrasts.arg = fit$contrasts)
  if (!identical(colnames(mm), names(beta)))
    stop("The model matrix could not be reconstructed at chosen covariate ",
         "values, so the group difference cannot be formed. Refit the model ",
         "with plain variables (no transformations applied in the formula) ",
         "and try again.", call. = FALSE)

  at_zero  <- mm[seq(1L, 2L * n_levels, by = 2L), , drop = FALSE]
  per_unit <- mm[seq(2L, 2L * n_levels, by = 2L), , drop = FALSE] - at_zero

  tcrit <- if (method == "simultaneous") {
    sqrt(2 * stats::qf(conf_level, 2, df_error))
  } else {
    stats::qt(1 - (1 - conf_level) / 2, df_error)
  }

  covariate_values <- mf[[cov_name]]
  group_as_char <- as.character(group_column)
  covariate_range <- range(covariate_values, na.rm = TRUE)

  rows <- list()
  pair_info <- list()
  geometry  <- list()

  for (i in seq_len(n_levels - 1L)) {
    for (j in seq(i + 1L, n_levels)) {
      # Later level minus earlier level, the TukeyHSD reading order.
      contrast_0 <- at_zero[j, ]  - at_zero[i, ]
      contrast_1 <- per_unit[j, ] - per_unit[i, ]
      L <- rbind(contrast_0, contrast_1)
      d  <- as.numeric(L %*% beta)
      Vd <- L %*% vcov_beta %*% t(L)
      d0 <- d[1L]; d1 <- d[2L]
      v00 <- Vd[1L, 1L]; v01 <- Vd[1L, 2L]; v11 <- Vd[2L, 2L]

      sol <- .regions_solve(d0, d1, v00, v01, v11, tcrit)

      in_pair <- group_as_char %in% c(levs[i], levs[j])
      pair_range <- range(covariate_values[in_pair], na.rm = TRUE)
      in_range <- function(b) {
        if (is.na(b)) NA_real_ else
          as.numeric(b >= pair_range[1L] && b <= pair_range[2L])
      }
      lower_in <- in_range(sol$lower)
      upper_in <- in_range(sol$upper)

      pair_label <- paste(levs[j], "-", levs[i])
      region <- .regions_describe(sol$code, sol$lower, sol$upper, cov_name)
      if (isTRUE(lower_in == 0) || isTRUE(upper_in == 0))
        region <- paste0(region, "; a boundary lies outside the observed ",
                         "range of ", cov_name, " [",
                         .regions_fmt(pair_range[1L]), ", ",
                         .regions_fmt(pair_range[2L]),
                         "] and is an extrapolation")

      rows[[pair_label]] <- data.frame(
        pair = pair_label,
        term = c("lower_bound", "upper_bound", "region_code",
                 "n_boundaries", "lower_bound_in_range",
                 "upper_bound_in_range", "difference_intercept",
                 "difference_slope", "critical_value", "df_error",
                 "conf_level"),
        value = c(sol$lower, sol$upper, as.numeric(sol$code),
                  as.numeric(sol$n), lower_in, upper_in, d0, d1,
                  tcrit, df_error, conf_level),
        stringsAsFactors = FALSE
      )

      pair_info[[pair_label]] <- data.frame(
        pair = pair_label, group_1 = levs[j], group_2 = levs[i],
        region = region, covariate_min = pair_range[1L],
        covariate_max = pair_range[2L], stringsAsFactors = FALSE
      )

      geometry[[pair_label]] <- data.frame(
        pair = pair_label,
        difference_intercept = d0, difference_slope = d1,
        var_intercept = v00, cov_intercept_slope = v01, var_slope = v11,
        lower_bound = sol$lower, upper_bound = sol$upper,
        region_code = as.numeric(sol$code),
        critical_value = tcrit,
        covariate_min = pair_range[1L], covariate_max = pair_range[2L],
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  attr(out, "method")          <- method
  attr(out, "outcome")         <- as.character(attr(model_terms,
                                                    "variables"))[2L]
  attr(out, "covariate")       <- cov_name
  attr(out, "group")           <- grp_name
  attr(out, "covariate_range") <- covariate_range
  row_bind <- function(x) do.call(rbind, c(x, list(make.row.names = FALSE)))
  attr(out, "pairs")           <- row_bind(pair_info)
  attr(out, "geometry")        <- row_bind(geometry)

  out <- .as_dmar_tbl(out, conf_level = conf_level)
  class(out) <- c("dmar_regions_of_significance", class(out))
  out
}
