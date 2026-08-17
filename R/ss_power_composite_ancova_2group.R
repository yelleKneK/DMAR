#' Sample Size or Composite Power for a Two-Group ANCOVA With a Covariate
#'
#' Determine the necessary per-group sample size to achieve a desired level of
#' composite statistical power in a two-group analysis of covariance, or, given a
#' per-group sample size, return the realized composite power. Composite power is
#' the probability that every effect named in \code{composite_terms} is
#' statistically significant in the same study, which is the quantity a design
#' has to be planned against when its conclusion requires more than one result to
#' hold at once.
#'
#' This is the two-group special case, kept under its own name for the simple
#' \code{smd}/\code{rho} interface it allows. For a one-way design with more than
#' two groups, or a factorial design, use the general
#' \code{\link{ss_power_composite_ancova}}.
#'
#' The model is
#' \deqn{Y = b_0 + b_{group} G + b_{cov} X + b_{group \times cov} G X + e.}
#' The group effect, the covariate effect, and the group by covariate interaction
#' can each be named in the composite, alone or in any combination.
#'
#' Calling \code{\link[=plot.dmar_composite_power]{plot}} on the result draws the
#' population effects the planning values describe. The figure needs only the
#' effect sizes, and nothing is simulated to draw it.
#'
#' @param smd Supposed standardized mean difference (Cohen's \emph{d}) between
#'   the two groups at the mean of the covariate, standardized by \code{sigma}: a
#'   value the researcher posits for the population, either a minimally important
#'   effect or a value believed to be true, never a sample estimate. Defaults to
#'   0, which is no group effect and therefore power equal to \code{alpha_level}
#'   for that test.
#' @param rho Supposed within-group population correlation between the covariate
#'   and the outcome. Length 1 for the same correlation in both groups, in which
#'   case the two slopes are equal and the interaction is zero, or length 2 for
#'   one correlation per group, in which case the slopes differ and the
#'   interaction carries the difference. Each element must lie in (-1, 1).
#'   Defaults to 0.
#' @param sigma Within-group population standard deviation of the outcome: the
#'   same \eqn{\sigma} as in a one-way ANOVA on the outcome, the same in both
#'   groups, and the standardizer of \code{smd}. Defaults to 1, which puts
#'   \code{smd} and the coefficients on a standardized scale.
#' @param sd_cov Population standard deviation of the covariate. Defaults to 1.
#'   A correlation is scale free, so \code{sd_cov} does not affect power; it sets
#'   the units the slopes and the figure are expressed in.
#' @param composite_terms Character vector naming the effects that must all be
#'   statistically significant. Any subset of \code{"group"}, \code{"covariate"},
#'   and \code{"group_by_covariate"}. A single term returns that test's ordinary
#'   power. Defaults to all three.
#' @param include_interaction Logical: whether the fitted model contains the
#'   group by covariate interaction. Defaults to \code{TRUE}. Dropping it returns
#'   one residual degree of freedom.
#' @param desired_power Desired composite statistical power (default 0.85). Used
#'   only when \code{n} is \code{NULL}.
#' @param alpha_level Type I error rate for each individual test (default 0.05).
#'   This is the per-test rate, not a rate for the composite event.
#' @param n Per-group sample size (assumed balanced); if specified, the realized
#'   power is returned rather than a sample size planned.
#' @param directional Logical: \code{TRUE} for one-sided tests, each in the
#'   direction of its own supposed effect, \code{FALSE} (default) for two-sided
#'   tests.
#'
#' @details
#' Coding the group factor as -1/2 and +1/2 and centering the covariate makes
#' each coefficient read directly: \eqn{b_{group}} is the difference between the
#' group means at the covariate mean, \eqn{b_{cov}} is the average of the two
#' within-group slopes, and \eqn{b_{group \times cov}} is the difference between
#' them. With balanced groups and a covariate whose distribution does not differ
#' across them, the columns of the design are mutually orthogonal, so the three
#' tests are orthogonal and the coefficient estimates are uncorrelated.
#'
#' Orthogonal effects do not give independent tests. Every test divides by the
#' same estimated error standard deviation, so an error estimate that lands low
#' inflates all of the test statistics together. The tests are positively
#' dependent, and composite power is strictly larger than the product of the
#' marginal powers. Multiplying the marginal powers understates the composite;
#' the gap closes as the residual degrees of freedom grow and the error estimate
#' stabilizes. Composite power can never exceed the least powerful test in the
#' set, so the weakest effect governs the design.
#'
#' Conditional on the error estimate the tests are independent, which reduces the
#' composite to a one-dimensional integral over the chi square distribution of
#' that estimate. Adaptive quadrature evaluates the integral, so no data are
#' simulated and the result is deterministic to quadrature precision.
#'
#' Within group \emph{g} the slope is \eqn{\rho_g \sigma / \sigma_X} and the
#' residual variance is \eqn{\sigma^2 (1 - \rho_g^2)}, so correlations that
#' differ across the groups give slopes that differ. The error variance the
#' ANCOVA pools and estimates is the average of the two,
#' \eqn{\sigma^2_{adj} = \sigma^2 (1 - \bar{\rho^2})} with \eqn{\bar{\rho^2}} the
#' mean of \eqn{\rho_1^2} and \eqn{\rho_2^2}, which is the familiar
#' \eqn{\sigma \sqrt{1 - \rho^2}} whenever the two correlations are equal in
#' absolute value. Correlations that are not make the residual variance differ
#' across the groups, and that has consequences the section on unequal residual
#' variances below spells out.
#'
#' The noncentralities are formed as \eqn{\sqrt{N} f}, the convention the rest of
#' the \code{ss_power_*} family uses (see \code{\link{ss_power_reg_coef}}), and
#' the residual degrees of freedom are \eqn{N - p - 1}. When the correlations are
#' equal and the interaction is dropped, the group test reproduces
#' \code{\link{ss_power_c_ancova}} with contrast weights \code{c(1, -1)}.
#'
#' Two approximations are worth separating, because they have different causes.
#' The first is the conditioning on the covariate, which is present at every set
#' of planning values and is described next. The second appears only when the
#' correlations differ in absolute value, and has its own section below.
#'
#' All three noncentralities condition on the covariate, substituting the
#' expected cross-product matrix for the expectation of its inverse. Inversion is
#' convex, so the assumed sampling variance is too small and every power in the
#' table is overstated, by an amount of order \eqn{1/N}. The group test is not
#' exempt: \eqn{b_{group}} is the difference at the covariate mean, and the
#' realized group covariate means are not exactly equal, so the group coefficient
#' inherits their sampling variability even under balanced groups and a covariate
#' independent of them. What the group test's noncentrality does not involve is
#' \code{sd_cov}, because a correlation is scale free; that is not the same as
#' exactness. The Type I error rate is unaffected, so it is specifically power
#' that is overstated, not the calibration of the tests. This is the same
#' approximation, and the same convention, that \code{\link{ss_power_c_ancova}}
#' and \code{\link{ss_power_reg_coef}} already use, and the correction
#' \code{ss_power_c_ancova} names in its own documentation is one instance of it.
#'
#' The size of that overstatement was measured against simulation for a two-term
#' composite with equal correlations, which isolates the conditioning from
#' everything else, at 200,000 replications per cell: about 2 points of composite
#' power at \emph{n} = 25 per group, under half a point at \emph{n} = 50, and
#' within Monte Carlo error from \emph{n} = 100 up. The constant grows with the
#' number of tests in the composite, because each contributes the same
#' approximation, so a three-term composite is worse than this at any given
#' \emph{n}. Treat these as the scale of the effect rather than a correction to
#' apply: for planned samples of a few dozen per group the reported power is
#' optimistic by a point or two, and for the sample sizes a three-term composite
#' usually requires the effect has died away.
#'
#' Because every test divides by the same error estimate, composite power is not
#' monotone in \emph{n} at the smallest residual degrees of freedom: an error
#' estimate that lands low at one or two degrees of freedom inflates all of the
#' statistics at once, so the composite there can exceed its value at slightly
#' larger \emph{n}. \code{necessary_n_per_group} (or
#' \code{approximate_n_per_group}, see the section on unequal residual variances)
#' is the smallest \emph{n} attaining \code{desired_power}, which for a target
#' below \code{alpha_level} need not be a size that every larger \emph{n} also
#' attains.
#'
#' @section Unequal Residual Variances When the Correlations Differ:
#' Within group \emph{g} the residual variance is
#' \eqn{\tau_g^2 = \sigma^2 (1 - \rho_g^2)}. Correlations that differ in
#' absolute value therefore leave the two groups with different residual
#' variances, and two things the composite integral assumes stop being true.
#'
#' The pooled error is no longer one scaled chi square. The residual sum of
#' squares is \eqn{\tau_1^2 Q_1 + \tau_2^2 Q_2} with \eqn{Q_1} and \eqn{Q_2}
#' independent chi square variables on \eqn{n - 2} degrees of freedom each, a
#' mixture of two scaled chi squares. Averaging the squared correlations gets
#' its mean exactly right and its spread wrong: the mixture's variance exceeds
#' that of the single scaled chi square the integral uses by
#' \eqn{(n - 2)(\tau_1^2 - \tau_2^2)^2}, which is zero exactly when the two
#' residual variances agree.
#'
#' The numerators stop being uncorrelated. The covariate coefficient is the
#' average of the two within-group slopes and the interaction coefficient is
#' their difference, and slopes estimated with different residual variances
#' leave those two estimators correlated:
#' \eqn{\mathrm{Cov} (\hat{b}_{cov}, \hat{b}_{group \times cov}) =
#' (\tau_2^2 - \tau_1^2) / (2 n \sigma_X^2)}, again zero exactly when the
#' residual variances agree. The tests are then not independent even given the
#' error estimate, which is the step that reduced the composite to a
#' one-dimensional integral in the first place.
#'
#' What the function reports in that case is therefore the composite power of a
#' design whose pooled error is a single scaled chi square with the right mean
#' and whose tests are conditionally independent, which is a near neighbor of
#' the design described but not that design. It is an approximation, and the
#' output says so: every row carrying a power is renamed with an
#' \code{approximate_} prefix, and a planned sample size is reported as
#' \code{approximate_n_per_group} and \code{approximate_N} rather than as
#' \code{necessary_n_per_group} and \code{necessary_N}, because it is the
#' smallest \emph{n} at which the approximation reaches \code{desired_power},
#' not an \emph{n} known to attain it. The numbers are the same numbers; only
#' the names change, and only in this case.
#'
#' The sign of the error is not guaranteed. Two fixed-covariate simulations, run
#' with the covariate values held to the population moments so that the
#' conditioning approximation above plays no part, bracket it. With
#' \code{smd = 0.30}, \code{rho = c(0, 0.9)} and \emph{n} = 8 per group, the
#' composite of the group and covariate effects is reported as 0.0752 against a
#' simulated 0.0865 (Monte Carlo standard error 0.0002), so the report is
#' conservative. With \code{smd = 2.50}, \code{rho = c(0, 0.95)} and \emph{n} = 6
#' per group, the group effect alone is reported as 0.9990 against a simulated
#' 0.9979 (Monte Carlo standard error 0.0001), so the report is optimistic.
#' Both are deliberately severe: a correlation gap of 0.9 and a handful of cases
#' per group. At a gap a covariate plausibly shows, \code{smd = 0.50} with
#' \code{rho = c(0.1, 0.5)}
#' and \emph{n} = 25 per group, the composite of the group effect and the
#' interaction is reported as 0.1515 against a simulated 0.1508 (Monte Carlo
#' standard error 0.0006), a difference inside simulation error.
#'
#' What to do with the number, then. Read it as an approximation whose accuracy
#' degrades with the gap between the correlations and not with \emph{N}, and
#' confirm a design you intend to run by simulating it: draw each group's errors
#' with its own residual standard deviation \eqn{\sigma \sqrt{1 - \rho_g^2}},
#' fit the same model, and count the replications in which every test in the
#' composite rejects. Equal absolute correlations need none of this. Two
#' correlations of the same magnitude and opposite sign, \code{rho = c(0.5,
#' -0.5)}, give a large interaction and still equal residual variances, so that
#' design is exact and its rows keep the ordinary names.
#'
#' @section Planning Without the Composite:
#' Composite power is the right quantity only when the conclusion needs several
#' results at once. When one effect carries the argument, DMAR already plans for
#' it and this function is unnecessary.
#'
#' For the group effect, \code{\link{ss_power_c_ancova}} is the planner: a
#' two-group comparison is the contrast \code{c(1, -1)} on the adjusted means.
#' Naming one term here reproduces it exactly, which is the check the tests
#' assert:
#'
#' \preformatted{
#' ss_power_c_ancova(psi = 0.5, c_weights = c(1, -1), sigma = 1, rho = 0.3,
#'                   n = 30)
#' ss_power_composite_ancova_2group(smd = 0.5, rho = 0.3, n = 30,
#'                           composite_terms = "group",
#'                           include_interaction = FALSE)
#' }
#'
#' Both return 0.5143. The interaction is dropped in the second call because
#' \code{ss_power_c_ancova} plans for the model without it, and carrying a term
#' the other function does not have would spend a residual degree of freedom on
#' nothing. With more than two groups, or a contrast other than a simple
#' difference, \code{ss_power_c_ancova} is the only one of the two that applies.
#'
#' For the design with no covariate at all, \code{\link{ss_power_smd}} is the
#' planner, and comparing the two is the cleanest way to see what a covariate
#' buys:
#'
#' \preformatted{
#' ss_power_smd(smd = 0.5, n_1 = 30)                    # no covariate
#' ss_power_c_ancova(psi = 0.5, c_weights = c(1, -1),
#'                   sigma = 1, rho = 0.5, n = 30)      # covariate, rho = .5
#' }
#'
#' Power rises from 0.4779 to 0.5942 because the covariate removes
#' \eqn{\rho^2} of the error variance, at the cost of one degree of freedom.
#' That is the ANCOVA bargain, and it is worth making before reaching for a
#' composite.
#'
#' @return
#' A \code{data.frame} with \code{term} and \code{value} columns. The design
#' result comes first, then the marginal power and noncentrality of each test in
#' the composite, then rows echoing the planning values, so the assumptions the
#' power was evaluated under travel with the result. The \code{tails} row is 2 for
#' nondirectional tests and 1 for directional tests. The names of the composite
#' terms, the implied coefficients, the per-group slopes, and \eqn{\sigma_{adj}}
#' are carried as attributes rather than rows, keeping the \code{value} column
#' numeric. Call \code{\link[=plot.dmar_composite_power]{plot}} on the result to
#' draw the population effects.
#'
#' When the two correlations differ in absolute value the powers are
#' approximations, and the row names say so: \code{composite_power} is reported
#' as \code{approximate_composite_power}, each \code{power_<term>} as
#' \code{approximate_power_<term>}, and a planned size as
#' \code{approximate_n_per_group} and \code{approximate_N}. The section on
#' unequal residual variances explains why. An \code{approximate} attribute
#' carries the same flag for a program to test without parsing row names.
#' \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} read both
#' sets of names, so a relabeled table still summarizes to its sample size and
#' its power, and \code{glance()} keeps the \code{approximate_} names on the
#' columns it carries through. Nothing changes when the correlations are equal
#' in absolute value, which is the exact case.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective} (4th ed.).
#'   Routledge. (See Chapter 9 on analysis of covariance, and Chapter 3 on
#'   statistical power and the noncentral distributions the tests follow.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_composite_ancova}} for the general one-way or
#'   factorial ANCOVA composite (of which this is the two-group case);
#'   \code{\link{ss_power_c_ancova}} for a single contrast on the adjusted
#'   means, which is the non-composite planner for this design;
#'   \code{\link{ss_power_smd}} for the two-group design with no covariate;
#'   \code{\link{ss_power_c}} for a contrast with no covariate;
#'   \code{\link{ss_power_reg_coef}}; \code{\link{ci_c_ancova}} and
#'   \code{\link{ci_sc_ancova}} for intervals on the adjusted means;
#'   \code{\link{ancova}} to fit the model the plan is for
#'
#' @examples
#' # A covariate correlating 0.10 with the outcome in one group and 0.40 in the
#' # other. The correlations differ, so the slopes differ and there is an
#' # interaction to detect; nothing had to be assumed about that in advance.
#' # They also differ in absolute value, so the groups' residual variances
#' # differ, the powers are approximations, and the rows are named to say so.
#' ss_power_composite_ancova_2group(smd = 0.20, rho = c(0.10, 0.40), n = 100)
#'
#' # The covariate is nearly certain to be detected and the interaction is
#' # better than even, but the group effect is weak, and the composite of all
#' # three is far below any of them. No single marginal power reveals that.
#'
#' # Equal correlations mean equal slopes, so the interaction is zero and its
#' # test rejects only at the Type I error rate. Asking for a composite that
#' # includes it therefore asks for something that cannot happen often.
#' ss_power_composite_ancova_2group(smd = 0.20, rho = 0.10, n = 100)
#'
#' # The composite of the two effects that are actually present.
#' ss_power_composite_ancova_2group(smd = 0.20, rho = 0.10, n = 100,
#'                           composite_terms = c("group", "covariate"))
#'
#' # Per-group sample size for composite power of 0.80 on the group effect and
#' # the interaction together.
#' ss_power_composite_ancova_2group(smd = 0.50, rho = c(0.10, 0.50),
#'                           composite_terms = c("group", "group_by_covariate"),
#'                           desired_power = 0.80)
#'
#' # Composite power is not the product of the marginal powers. The tests share
#' # one error estimate, so they are positively dependent and the composite is
#' # larger than the product. Multiplying would understate the design. These
#' # correlations differ in absolute value, so the rows carry the approximate_
#' # prefix; see the section on unequal residual variances for what that means.
#' plan <- ss_power_composite_ancova_2group(smd = 0.50, rho = c(0.10, 0.50), n = 95,
#'                                   composite_terms = c("group",
#'                                                       "group_by_covariate"))
#' plan$value[plan$term == "approximate_composite_power"]
#' prod(plan$value[plan$term %in% c("approximate_power_group",
#'                                  "approximate_power_group_by_covariate")])
#'
#' # Correlations of equal magnitude and opposite sign leave the residual
#' # variances equal, so this design is exact and keeps the ordinary row names,
#' # even though the two slopes could hardly differ more.
#' exact <- ss_power_composite_ancova_2group(smd = 0.50, rho = c(0.40, -0.40),
#'                                    n = 95,
#'                                    composite_terms = c("group",
#'                                                        "group_by_covariate"))
#' exact$value[exact$term == "composite_power"]
#'
#' # Draw the population effects a result was planned on. The figure needs only
#' # the effect sizes, and nothing is simulated to draw it.
#' plot(ss_power_composite_ancova_2group(smd = 0.20, rho = c(0.10, 0.40), n = 100))
#'
#' # Planning without the composite: one term returns that test's ordinary
#' # power, and reproduces ss_power_c_ancova once the interaction is dropped.
#' ss_power_composite_ancova_2group(smd = 0.50, rho = 0.30, n = 30,
#'                           composite_terms = "group",
#'                           include_interaction = FALSE)
#' ss_power_c_ancova(psi = 0.50, c_weights = c(1, -1), sigma = 1,
#'                   rho = 0.30, n = 30)
#'
#' # What the covariate buys, against the same design with no covariate.
#' ss_power_smd(smd = 0.50, n_1 = 30)
#'
#' # The broom verbs summarize the plan in one row.
#' generics::tidy(ss_power_composite_ancova_2group(smd = 0.50, rho = c(0.10, 0.50),
#'                                          composite_terms = c("group",
#'                                                          "group_by_covariate"),
#'                                          desired_power = 0.80))
#'
#' @keywords design htest
#'
#' @family sample size for power
#'
#' @family composite power
#'
#' @export
ss_power_composite_ancova_2group <- function(smd = 0,
                                      rho = 0,
                                      sigma = 1,
                                      sd_cov = 1,
                                      composite_terms = c("group", "covariate",
                                                          "group_by_covariate"),
                                      include_interaction = TRUE,
                                      desired_power = 0.85,
                                      alpha_level = 0.05,
                                      n = NULL,
                                      directional = FALSE) {
  rho <- .check_ancova_args(smd, rho, sigma, sd_cov, alpha_level)

  allowed <- c("group", "covariate", "group_by_covariate")
  if (!is.character(composite_terms) || !length(composite_terms)) {
    stop("'composite_terms' must be a character vector naming at least one of: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
  bad <- setdiff(composite_terms, allowed)
  if (length(bad)) {
    stop("'composite_terms' contains terms that are not in the model: ",
         paste(bad, collapse = ", "), ". Allowed terms are: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
  if (anyDuplicated(composite_terms)) {
    stop("'composite_terms' contains duplicate terms.", call. = FALSE)
  }
  if (!is.logical(include_interaction) || length(include_interaction) != 1L ||
      is.na(include_interaction)) {
    stop("'include_interaction' must be TRUE or FALSE.", call. = FALSE)
  }
  if ("group_by_covariate" %in% composite_terms && !include_interaction) {
    stop("'group_by_covariate' is named in 'composite_terms' but ",
         "'include_interaction' is FALSE; an effect cannot be tested unless it ",
         "is in the model.", call. = FALSE)
  }
  # Correlations that differ put a nonzero interaction in the population. A
  # model that omits it does not make it go away: the omitted term is orthogonal
  # to the rest, so it biases nothing, but its variance moves into the error the
  # fitted model estimates. Planning against that model with the error variance
  # of the model that keeps the term would overstate every power in the table.
  # Rather than quietly plan for a misspecified model, refuse the combination.
  if (!isTRUE(all.equal(rho[1L], rho[2L])) && !include_interaction) {
    stop("'rho' differs across the groups, which puts a nonzero group by ",
         "covariate interaction in the population, but 'include_interaction' ",
         "is FALSE. A fitted model that omits a truly nonzero term absorbs it ",
         "into the error, so power planned that way would be overstated. ",
         "Either keep the interaction in the model, or supply a single 'rho' ",
         "if the slopes are taken to be equal.", call. = FALSE)
  }
  # Unequal absolute correlations make the two groups' residual variances differ,
  # so the pooled error the tests share is a mixture of two scaled chi squares
  # rather than the single one the composite integral assumes, and the average
  # slope and slope difference estimators become correlated. Every power in the
  # table is then an approximation, and a resolved sample size is the smallest
  # one at which the approximation reaches the target rather than a size known
  # to attain it. Warn, and carry the same statement into the row names so it
  # survives being read out of the table (see the "Unequal residual variances"
  # section of the documentation).
  approximate <- !isTRUE(all.equal(abs(rho[1L]), abs(rho[2L])))
  if (approximate) {
    warning("The group correlations differ in absolute value, so the groups' ",
            "residual variances differ and the shared-error composite power is ",
            "an approximation whose error grows with the gap between the ",
            "correlations. The powers are reported as approximate_* rows, and ",
            "a planned sample size as approximate_n_per_group, because neither ",
            "is exact; confirm a final design by simulation.", call. = FALSE)
  }

  # Keep the reported order fixed rather than dependent on the order the terms
  # were typed in.
  composite_terms <- allowed[allowed %in% composite_terms]

  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 2 ||
        n != round(n)) {
      stop("'n' must be a single whole number of at least 2.", call. = FALSE)
    }
    n <- as.integer(n)
  }

  evaluate_at <- function(n) {
    d <- .ancova_design(smd, rho, sigma, sd_cov, n, include_interaction)
    res <- .composite_power_shared_sigma(d$delta[composite_terms], d$df,
                                         alpha_level, directional)
    list(composite = res$composite,
         marginal = stats::setNames(res$marginal, composite_terms),
         design = d)
  }

  if (!is.null(n)) {
    res <- evaluate_at(n)
    if (res$design$df < 1) {
      stop("'n' is too small to leave residual degrees of freedom for this ",
           "model; increase 'n'.", call. = FALSE)
    }
    out <- data.frame(
      term = c("specified_n_per_group", "specified_N", "composite_power",
               "residual_df",
               paste0("power_", composite_terms),
               paste0("noncentral_t_parm_", composite_terms)),
      value = c(n, res$design$N, res$composite, res$design$df,
                unname(res$marginal),
                unname(res$design$delta[composite_terms])))
    return(.finish_composite_ancova(out, res$design, smd, rho, sigma, sd_cov,
                                    alpha_level, directional, composite_terms,
                                    include_interaction, approximate))
  }

  if (!is.numeric(desired_power) || length(desired_power) != 1L ||
      desired_power <= 0 || desired_power >= 1) {
    stop("'desired_power' must be a single numeric value in (0, 1).",
         call. = FALSE)
  }
  # A term whose supposed effect is zero rejects at the Type I error rate at
  # every sample size, so the composite has a ceiling below any desired power
  # worth planning for and the search would not terminate.
  zero_terms <- composite_terms[
    abs(.ancova_design(smd, rho, sigma, sd_cov, 100L,
                       include_interaction)$delta[composite_terms]) < 1e-12]
  if (length(zero_terms)) {
    stop("These terms in 'composite_terms' have a supposed effect of zero: ",
         paste(zero_terms, collapse = ", "),
         ". A null effect is significant only at the Type I error rate, at ",
         "every sample size, so no sample size attains 'desired_power'. Supply ",
         "a nonzero effect for these terms or drop them from ",
         "'composite_terms'.", call. = FALSE)
  }

  n_i <- 2L
  res <- evaluate_at(n_i)
  while (is.na(res$composite) || res$composite < desired_power) {
    n_i <- n_i + 1L
    res <- evaluate_at(n_i)
    if (n_i > 1e6) {
      stop("Failed to converge within a reasonable sample size.", call. = FALSE)
    }
  }
  out <- data.frame(
    term = c("necessary_n_per_group", "necessary_N", "composite_power",
             "residual_df",
             paste0("power_", composite_terms),
             paste0("noncentral_t_parm_", composite_terms)),
    value = c(n_i, res$design$N, res$composite, res$design$df,
              unname(res$marginal),
              unname(res$design$delta[composite_terms])))
  out <- rbind(out, data.frame(term = "desired_power", value = desired_power))
  .finish_composite_ancova(out, res$design, smd, rho, sigma, sd_cov,
                           alpha_level, directional, composite_terms,
                           include_interaction, approximate)
}


# Append the rows echoing the planning values, attach the metadata that must
# stay off the numeric value column, and route through the package display
# layer. The leading class is set before .as_dmar_tbl() so the helper inserts
# dmar_tbl ahead of data.frame and leaves the dispatch classes in front.
#
# When the correlations leave the groups with unequal residual variances the
# rows carrying a power, and a sample size resolved against one, are renamed
# with an approximate_ prefix. The numbers are the same either way; the names
# are what stop the table from claiming more than the method delivers.
.finish_composite_ancova <- function(out, design, smd, rho, sigma, sd_cov,
                                     alpha_level, directional, composite_terms,
                                     include_interaction,
                                     approximate = FALSE) {
  out <- rbind(out, data.frame(
    term = c("supposed_smd", "supposed_rho_group_1", "supposed_rho_group_2",
             "sigma", "sd_cov", "alpha_level", "tails"),
    value = c(smd, rho[1L], rho[2L], sigma, sd_cov, alpha_level,
              if (directional) 1 else 2)))
  if (approximate) out$term <- .composite_approximate_terms(out$term)
  # dmar_composite_power leads so plot() dispatches without claiming every
  # ss_power_* return; dmar_ss_power brings the shared tidy() and glance()
  # methods, which read the sample size and the power off named rows. Both the
  # exact and the approximate names are listed in .SS_POWER_SIZE_TERMS and
  # .SS_POWER_POWER_TERMS (R/dmar_tidiers.R), so the relabeling above does not
  # cost those verbs their rows.
  class(out) <- c("dmar_composite_power", "dmar_ss_power", "data.frame")
  attr(out, "composite_terms")     <- composite_terms
  attr(out, "include_interaction") <- include_interaction
  attr(out, "directional")         <- directional
  attr(out, "approximate")         <- approximate
  attr(out, "sigma_adj")           <- design$sigma_adj
  attr(out, "slopes")              <- c(group_1 = design$slope_1,
                                        group_2 = design$slope_2)
  attr(out, "coefficients")        <- c(group = design$b_group,
                                        covariate = design$b_cov,
                                        group_by_covariate =
                                          design$b_group_by_cov)
  .as_dmar_tbl(out)
}


#' @describeIn ss_power_composite_ancova_2group Draw the population effects a
#'   result was planned on, reached from a result already in hand. With
#'   \code{show_power = TRUE} (the default) the figure is annotated with the
#'   power the plan delivers.
#'
#' @param x An object returned by \code{ss_power_composite_ancova_2group}.
#' @param ... Further arguments passed to the figure: \code{cov_range},
#'   \code{palette}, \code{group_labels}, and \code{show_power}.
#'
#' @export
plot.dmar_composite_power <- function(x, ...) {
  dots <- list(...)
  arg <- function(nm, default) if (is.null(dots[[nm]])) default else dots[[nm]]

  # A table planned under unequal residual variances carries approximate_ row
  # names. Restoring the exact names here lets the lookups below name one thing.
  x <- .composite_as_exact(x)
  val <- function(term) {
    v <- x$value[x$term == term]
    if (!length(v)) NA_real_ else v[1L]
  }
  n <- val("specified_n_per_group")
  if (is.na(n)) n <- val("necessary_n_per_group")

  terms <- attr(x, "composite_terms")
  design <- .ancova_design(val("supposed_smd"),
                           c(val("supposed_rho_group_1"),
                             val("supposed_rho_group_2")),
                           val("sigma"), val("sd_cov"), n,
                           attr(x, "include_interaction"))

  power_info <- NULL
  if (arg("show_power", TRUE)) {
    power_info <- list(
      marginal = stats::setNames(
        vapply(terms, function(tm) val(paste0("power_", tm)), numeric(1)),
        terms),
      composite = val("composite_power"), n = n,
      alpha_level = val("alpha_level"), tails = val("tails"))
  }
  .plot_ancova_effects(design, val("supposed_smd"), val("sigma"),
                       val("sd_cov"), arg("cov_range", 2),
                       arg("palette", "okabe_ito"),
                       arg("group_labels", c("Group 1", "Group 2")),
                       power_info)
}
