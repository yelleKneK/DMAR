#' Mediation Analysis via Model-Based Constrained Optimization
#'
#' Tests hypotheses about mediation effects with the model-based
#' constrained optimization (MBCO) procedure of Tofighi and Kelley (2020):
#' a likelihood ratio test that compares the full mediation model to a
#' null model fit subject to the nonlinear constraint that the effect of
#' interest (an indirect effect, a total effect, or any smooth function of
#' path coefficients) equals zero. The model is specified in
#' \pkg{lavaan} syntax but is fit with \pkg{OpenMx}, whose optimizers
#' support the nonlinear equality constraints the procedure requires (see
#' Details). The function enumerates the mediation effects implied by the
#' model (the total effect, the direct effect, the total indirect effect,
#' and every specific indirect pathway from \code{x} to \code{y}),
#' reports each with a confidence interval, and tests each with the MBCO
#' likelihood ratio statistic and its \emph{p}-value. Models may include
#' observed or latent variables, one or many mediators, in parallel or in
#' sequence, and raw data or summary statistics may be supplied. With a
#' grouping variable the model becomes a multiple-group SEM: every
#' effect is estimated within each group and the between-group
#' difference of each effect is tested, which is moderated mediation
#' with a categorical moderator. With a numeric \code{moderator},
#' interactions stated in the model are probed: conditional effects at
#' chosen moderator values, the index of moderated mediation, and a
#' joint constancy test, all with MBCO likelihood ratio inference (see
#' Details).
#'
#' @param model A character string with the model in \pkg{lavaan} syntax
#'   (regressions \code{~}, latent variable definitions \code{=~},
#'   residual covariances \code{~~}, optionally labeled coefficients and
#'   defined quantities via \code{:=}), or a fitted \pkg{lavaan} object
#'   whose parameter table and data are then reused. Parameter labels
#'   must be valid R names (e.g., \code{b1}, not \code{1b}).
#' @param data A \code{data.frame} of raw data containing the observed
#'   variables. Rows with missing values are retained and handled by full
#'   information maximum likelihood. Supply either \code{data} or
#'   \code{S} (with \code{N}), not both.
#' @param S A covariance matrix of the observed variables (with dimnames
#'   naming the variables), used with \code{N} and optionally \code{M}
#'   when raw data are not available. Complete-data maximum likelihood
#'   depends on the data only through the sample means and covariance
#'   matrix, so an analysis from summary statistics reproduces the
#'   raw-data analysis exactly (see Details). For a multiple-group
#'   analysis, a named list of such matrices, one per group.
#' @param M Optional named numeric vector of variable means accompanying
#'   \code{S} (for multiple groups, a named list of such vectors). When
#'   omitted, means of zero are used; intercepts are then uninteresting
#'   but every effect, test, and confidence interval is unaffected.
#' @param N Total sample size. Required with \code{S}; ignored for raw
#'   data. For a multiple-group analysis, a vector with one sample size
#'   per group in \code{S} (matched by name when named).
#' @param group Optional name (single character string) of the grouping
#'   variable in \code{data} for a multiple-group analysis. With
#'   summary input, supply \code{S}, \code{M}, and \code{N} as lists
#'   instead and leave \code{group} alone (it then merely names the
#'   constructed grouping column). See the Details section on multiple
#'   groups, including how lavaan's labeling rules decide which
#'   parameters differ by group.
#' @param moderator Optional name (single character string) of a
#'   numeric moderator variable, turning on probing: every pathway
#'   effect that involves an interaction with the moderator is
#'   estimated at each probe value, the change in the effect per unit
#'   moderator is reported (for a pathway moderated in one place, the
#'   index of moderated mediation), and a constancy test asks whether
#'   the pathway effect depends on the moderator at all. State the
#'   interaction in the model syntax with a \code{:} term (e.g.,
#'   \code{m ~ x + w + x:w}) or include a precomputed product column
#'   among the predictors; both are recognized. A categorical
#'   moderator goes through \code{group} instead.
#' @param probe_values Numeric vector of moderator values at which
#'   moderated effects are estimated; names, when given, label the
#'   rows. Defaults to the moderator's mean and one standard deviation
#'   either side (labeled \code{low}, \code{mean}, \code{high}), or to
#'   the two observed values of a two-valued moderator.
#' @param x,y Names (single character strings) of the focal predictor and
#'   the outcome between which effects are traced. When both are omitted
#'   the function looks for a unique variable with no incoming regression
#'   path (the source) and a unique variable with no outgoing regression
#'   path (the outcome); if either is ambiguous an error asks for
#'   \code{x} and \code{y} explicitly. When the model defines quantities
#'   via \code{:=} or \code{hypotheses} are supplied, \code{x} and
#'   \code{y} may be omitted entirely and only those quantities are
#'   estimated and tested.
#' @param hypotheses Optional named character vector of additional
#'   quantities to estimate and test against zero, written as R
#'   expressions in the parameter labels and effect names (e.g.,
#'   \code{c(difference = "indirect_via_imagery -
#'   indirect_via_repetition")} tests the equality of two specific
#'   indirect effects). A named \emph{list} may mix such scalar
#'   expressions with character vectors of several expressions; a
#'   multi-expression element is tested jointly (all constrained to
#'   zero at once, with as many degrees of freedom as expressions),
#'   and its row reports the test alone, with no scalar estimate. Each
#'   expression of a joint test must restrict the parameters beyond what
#'   the others do. An expression given twice, a multiple of another
#'   (\code{"a"} and \code{"2*a"}), a combination of others (\code{"a"},
#'   \code{"b"}, and \code{"a + b"}; \code{"a - b"} and \code{"b - a"}),
#'   or an expression that is zero wherever the others are (\code{"a*b"}
#'   with \code{"b"}; \code{"a*cp"} with \code{"a"}; \code{"a*b*cp"} with
#'   \code{"a*b"}) would be counted as a degree of freedom while
#'   restricting nothing, so such a set is refused before any model is
#'   fit, with a condition of class \code{dmar_argument_invalid} that names
#'   the redundant expressions (as \code{redundant}). Equivalent sets of
#'   expressions that each restrict, such as \code{c("a", "b")} and
#'   \code{c("a + b", "a - b")}, give the same test. See \emph{Constraints
#'   that add no restriction} under Details.
#' @param ci_method Confidence interval method for every reported
#'   effect: \code{"profile_likelihood"} (default; Neale & Miller, 1997,
#'   computed by \pkg{OpenMx}), \code{"monte_carlo"} (simulate the
#'   coefficients from their joint normal approximation and form
#'   percentile intervals of the effect; MacKinnon, Lockwood, &
#'   Williams, 2004), or \code{"wald"} (delta method symmetric
#'   interval; reported for comparison, not recommended for indirect
#'   effects, whose sampling distributions are skewed).
#' @param conf_level Confidence level. Defaults to 0.95.
#' @param B Number of Monte Carlo draws when \code{ci_method =
#'   "monte_carlo"}, a single whole number (not read, and not checked, for
#'   the other intervals). Defaults to 10000, because a draw costs one
#'   evaluation of each effect and no refit. The fewest accepted is 400 at
#'   the default \code{conf_level} (200 at 0.90, 2000 at 0.99); a smaller
#'   \code{B} is refused before the model is fit (class
#'   \code{dmar_count_invalid}). See \emph{Monte Carlo draws} under
#'   \emph{Details} and \code{\link{dmar_bootstrap}}.
#' @param optimizer Optimizer used by \pkg{OpenMx} for all fits:
#'   \code{"SLSQP"} (default), \code{"CSOLNP"}, or \code{"NPSOL"}.
#'   The CRAN build of \pkg{OpenMx} ships SLSQP and CSOLNP; NPSOL, the
#'   optimizer used in Tofighi and Kelley (2020), is available only in
#'   the build distributed from the OpenMx website.
#' @param seed Optional integer seed for the Monte Carlo interval, used
#'   locally (the caller's random number generator state is restored on
#'   exit). Default \code{NULL} draws the Monte Carlo replications from the
#'   current random number stream, which advances that stream.
#'
#' @details
#' \strong{The MBCO procedure.} A hypothesis about mediation is recast as
#' a model comparison (Maxwell, Delaney, & Kelley, 2027). The full model
#' estimates all path coefficients freely. The null model is the same
#' model estimated subject to the constraint that the effect of interest
#' is zero, for example \eqn{H_0\!: \beta_1 \beta_2 = 0} for the indirect
#' effect of \eqn{X} on \eqn{Y} through one mediator. Writing
#' \eqn{D = -2\,\mathrm{log\,likelihood}} for the deviance of each fitted
#' model, the test statistic is
#' \deqn{\mathrm{LRT}_{\mathrm{MBCO}} = D_{\mathrm{null}} -
#'   D_{\mathrm{full}},}
#' and the \emph{p}-value refers it to the chi square distribution with
#' degrees of freedom equal to the number of constraints (here 1). That
#' distribution is the large sample distribution of the statistic when the
#' null hypothesis is true (Wilks, 1938) and regular: the constraints are
#' smooth and functionally independent, and the population parameter is
#' identified, lies in the interior of the parameter space, and is a point
#' near which the null set is a smooth surface with one dimension fewer
#' per constraint than the full model. The null hypothesis of no indirect
#' effect is regular at every point of its null set except those where both
#' of its coefficients are zero (whatever the values of the other
#' parameters), where the chi square reference does not hold (see
#' \emph{The reference distribution where the branches cross} below). The
#' \emph{p}-value is a
#' continuous measure of the compatibility of the data with the null
#' model, not merely a reject or fail-to-reject decision. In the
#' simulations of Tofighi and Kelley (2020, p. 508), with one coefficient
#' of the indirect effect zero, the Type I error rate of the test stayed
#' close to the nominal .05 for every sample size (50 to 200) and size of
#' the other coefficient, whereas the percentile bootstrap, Monte Carlo,
#' and profile likelihood intervals and the joint significance test
#' rejected too rarely when the other coefficient was small (0.14) and the
#' sample size was 100 or less.
#'
#' \strong{Why the constraint is nonlinear.} Null hypotheses about single
#' parameters are linear constraints: fixing \eqn{\beta_3 = 0} restricts
#' the parameter space to a flat hyperplane, something every structural
#' equation modeling program does by fixing the parameter and refitting.
#' The null hypothesis of no mediation is different in kind. The
#' constraint \eqn{\beta_1 \beta_2 = 0} involves a product of parameters,
#' so it is a nonlinear function of the parameter vector, and its
#' solution set is not a hyperplane but the union of two hyperplanes: the
#' set where \eqn{\beta_1 = 0} (with \eqn{\beta_2} free) and the set
#' where \eqn{\beta_2 = 0} (with \eqn{\beta_1} free). Mediation is absent
#' in infinitely many ways, and no single fixed parameter expresses them
#' all: fixing \eqn{\beta_1 = 0} alone tests a more restrictive
#' hypothesis than \eqn{H_0\!: \beta_1 \beta_2 = 0}. Maximizing the
#' likelihood over such a constraint set requires an optimizer that
#' handles general nonlinear equality constraints (sequential quadratic
#' programming methods such as SLSQP and NPSOL, or CSOLNP); among R
#' packages for structural equation modeling, \pkg{OpenMx} provides
#' them, which is why the model is fit there even though it is specified
#' in \pkg{lavaan} syntax. The same geometry explains why Wald-type
#' (delta method) tests of a product behave badly near the null: at
#' \eqn{\beta_1 = \beta_2 = 0} the two hyperplanes intersect, the
#' gradient of \eqn{\beta_1 \beta_2} vanishes, and the usual normal
#' approximation for the product breaks down. The MBCO procedure does not
#' use that approximation, since it compares maximized likelihoods
#' directly; the chi square reference of its statistic is itself a large
#' sample approximation, though, and it too fails at that point, as the
#' next paragraph explains.
#'
#' \strong{The reference distribution where the branches cross.} Where
#' exactly one of \eqn{\beta_1} and \eqn{\beta_2} is zero, the null set
#' near the population values is the single hyperplane of that
#' coefficient, the null hypothesis is regular there, and the chi square
#' distribution on 1 degree of freedom is the large sample distribution of
#' the statistic. At \eqn{\beta_1 = \beta_2 = 0} the two hyperplanes cross
#' and the null hypothesis is not regular. The best-fitting null model
#' lies on one branch or the other, so the statistic is the smaller of the
#' two branch statistics, each the likelihood ratio statistic for fixing
#' one coefficient at zero. When the estimates of the two coefficients are
#' independent, as in the single mediator model with observed variables,
#' the branch statistics are asymptotically the squares \eqn{Z_1^2} and
#' \eqn{Z_2^2} of two independent standard normal variables, and the
#' statistic is asymptotically \eqn{\min(Z_1^2, Z_2^2)} rather than a chi
#' square variable on 1 degree of freedom. It exceeds the .95 quantile of
#' that chi square distribution with probability \eqn{.05^2 = .0025}. At
#' \eqn{\beta_1 = \beta_2 = 0}, then, the test at the .05 level rejects
#' with probability near .0025 in large samples, not .05. When the two
#' estimates are correlated, as they can be in models with latent variables
#' and in some multiple-group or moderated designs, \eqn{Z_1} and
#' \eqn{Z_2} are correlated standard normal variables, and the probability
#' lies between .0025 and .05, nearer .05 the stronger the correlation
#' (.0049, .0125, and .0296 for correlations of .3, .6, and .9). The test
#' is conservative at the crossing in every case: its \emph{p}-value is
#' larger than one computed from the distribution that holds there, and
#' its power against an indirect effect whose coefficients are both near
#' zero is lower accordingly. For a chain of three or more coefficients
#' the same holds wherever two or more of them are zero, and the test is
#' more conservative still where more of them are zero (with independent
#' estimates, \eqn{.05^3} when three are zero). The function reports the
#' chi square \emph{p}-value throughout and does not adjust it, because
#' which point
#' of the null set holds is not known, and a critical value taken from the
#' distribution at the crossing would make the test reject too often where
#' only one coefficient is zero. The same qualification applies to
#' every test of a product of coefficients the function reports (specific
#' indirect effects, conditional effects at a probe value, and moderation
#' coefficients that are products), and to a joint test whose constraints
#' are satisfied together at a point where their gradients are linearly
#' dependent, such as the crossing of two branches.
#'
#' \strong{The rate at the crossing in the simulations of Tofighi and
#' Kelley (2020).} Their simulations report a different rate at
#' \eqn{\beta_1 = \beta_2 = 0}. For the single mediator model with
#' \eqn{X} binary (probability .5) and standard normal residuals, 5,000
#' replications at each of \emph{N} = 50, 75, 100, and 200, the .05
#' level, and their test fitted by the NPSOL optimizer of \pkg{OpenMx}
#' 2.9.6 (p. 507 and its footnote 4), the empirical Type I error rate of
#' \eqn{\mathrm{LRT}_{\mathrm{MBCO}}} with both coefficients zero fell
#' within the interval from .025 to .075 that they took as close enough
#' to .05 for practical purposes, at every sample size (p. 508; about .03
#' to .04 in their Figure 3, p. 509), and the other five methods they
#' studied rejected at rates below .01. The two numbers describe
#' different things. The .0025 above is a limit: the
#' rate that the statistic defined by the best-fitting null model
#' approaches as \emph{N} grows at that point. The published rates are
#' estimates from a simulation of one design at four sample sizes.
#' Neither is the rate at every sample size. For the statistic this
#' function computes, the rate in that design is already close to the
#' limit at these sample sizes: in 400,000 simulated data sets at each
#' \emph{N}, the smaller of the two branch statistics (the statistic
#' \code{mediation_mbco()} returns for this model, to within
#' \eqn{10^{-4}} on 630 further data sets of the design at \emph{N} = 50
#' to 200) exceeded the .95 chi square quantile at rates of
#' .0032, .0029, .0030, and .0028 for \emph{N} = 50, 75, 100, and 200
#' (Monte Carlo standard error .00009), so sample size in this range does
#' not account for the difference. The statistic of the worse-fitting
#' branch, the one a search reports when it stops on that branch (see
#' \emph{Local solutions and the search strategy} below), rejected at
#' .100 to .111 in the same simulations, and the published rates lie
#' between the two.
#'
#' \strong{Local solutions and the search strategy.} Because the null
#' set is a union of surfaces, the constrained deviance surface generally
#' has one local minimum per branch (one where the mediator does not
#' respond to \eqn{X}, one where the outcome does not respond to the
#' mediator, and so on along longer chains). A constrained optimizer
#' started from the full-model estimates can converge to whichever branch
#' is nearest rather than to the branch that fits best. The likelihood
#' ratio test is defined by the globally best-fitting null model, so
#' \code{mediation_mbco()} refits each null model from several starting
#' configurations (the full-model estimates, and the estimates with each
#' coefficient entering the constrained effect set to zero in turn),
#' verifies that each candidate solution actually satisfies the
#' constraint, and keeps the feasible solution with the smallest
#' deviance. In the memory example below this matters: for the single
#' mediator model the best-fitting null model sets the imagery-to-recall
#' path to zero (\eqn{\mathrm{LRT}_{\mathrm{MBCO}} = 71.31}, the
#' statistic the example reports), while the branch that sets the
#' instruction-to-imagery path to zero fits worse
#' (\eqn{\mathrm{LRT}_{\mathrm{MBCO}} = 179.02}). Tofighi and Kelley
#' (2020) report 175.77 for this example, a value from that
#' worse-fitting branch, on which their optimizer stopped; their
#' statistic also differs from the 179.02 here because the example runs
#' from the published (rounded) summary statistics of their Table 1
#' rather than the full-precision moments (which give 72.54 and 175.77
#' for the two branches). Either way the null model is overwhelmingly
#' incompatible with the data, so the substantive conclusion is the
#' same.
#'
#' The null model of a product is fit branch by branch: a product is zero
#' exactly where one of its factors is, so for each way of choosing one
#' factor of each constraint, a factor that is a coefficient is fixed at
#' zero and a factor that is a sum (the conditional coefficient
#' \eqn{a_0 + a_1 w} of a moderated path at a probe value \eqn{w}, for
#' example) is kept as the constraint of that fit. A label common to every
#' term of a sum is taken out as a factor
#' (\eqn{a_0 b + a_1 b w = b (a_0 + a_1 w)}). The null set is the union of
#' the branches, so when every branch fit is a valid solution the smallest
#' of them is the constrained minimum, and the constrained searches from
#' the starting configurations are run only when a branch fit is not, or
#' when a constraint has no branches (a total effect, a difference between
#' groups). Such a constraint is also searched from the full-model
#' estimates moved onto its null set by Gauss-Newton steps.
#'
#' \strong{Badly scaled variables.} When the observed variables differ in
#' scale by orders of magnitude (the ratio of the largest standard
#' deviation to the smallest is at least 1000, or the largest is at least
#' 1000, or the smallest at most 0.001, the rule \pkg{lavaan} uses for its
#' \code{rescale_data} option), the free parameters do too, and the
#' \pkg{OpenMx} optimizers can stop far from the maximum. The \pkg{OpenMx}
#' models are then fit to the variables each divided by a power of 2 near
#' its standard deviation (\eqn{2^{\mathrm{round}(\log_2 s)}} for standard
#' deviation \eqn{s}, over all groups), with every parameter in the
#' matching metric (a latent variable takes the metric of the indicator
#' whose loading is fixed), and every effect, constraint, and Monte Carlo
#' draw is written in the parameters of the original metric, so the
#' estimates, standard errors, intervals, and likelihood ratio statistics
#' are those of the variables as given, and the deviance, AIC, and BIC are
#' returned to the original metric as well. A likelihood ratio statistic
#' does not depend on the units of the variables, and the rescaling changes
#' it only within the precision of the optimizer. A profile likelihood
#' interval is profiled for the effect in the metric of the fit and its
#' limits returned to the original metric. The constraints of a
#' joint test are then handed to the optimizer as linear combinations
#' whose gradients are orthonormal in the metric of the fit at the
#' full-model estimates, which leaves their null set unchanged. Nothing is
#' rescaled when a label shared by parameters in different metrics (an
#' equality constraint across them) would need two multipliers. The
#' divisors and the multipliers are returned in the \code{"rescale"}
#' attribute. Scales that differ by less than that rule can still slow the
#' optimizers; the branch fits and the start on the null set described
#' above are the remedy there, and a variable can be put on a scale near
#' that of the others before the call (the statistics do not depend on its
#' units, and the estimates change by the factor applied).
#'
#' \strong{Effects estimated and tested.} With \code{x} and \code{y}
#' resolved, the function enumerates every directed pathway from
#' \code{x} to \code{y} along regression (\code{~}) paths. Each pathway
#' through at least one intermediate variable contributes a specific
#' indirect effect, the product of its path coefficients, reported as
#' \code{indirect_via_} followed by the intermediate variable names. The
#' direct effect is the \code{x} to \code{y} coefficient when that path
#' is in the model. The total indirect effect is the sum of the specific
#' indirect effects (reported when there are two or more), and the total
#' effect is the direct effect plus the total indirect effect (reported
#' when the direct path is in the model). Quantities defined in the
#' model syntax via \code{:=} and any \code{hypotheses} are estimated
#' and tested the same way, so contrasts of indirect effects, proportions
#' mediated, or any other smooth function of parameters can be examined;
#' each is tested against zero, so write an equality of two effects as
#' their difference. Every reported row carries its estimate, a delta
#' method standard error (via \code{\link[OpenMx]{mxSE}}), the
#' \code{ci_method} confidence interval, and the MBCO likelihood ratio
#' test with its degrees of freedom and \emph{p}-value.
#'
#' \strong{The parallel two-mediator model.} Tofighi and Kelley (2020)
#' continue the memory example with imagery and repetition as parallel
#' mediators and a residual covariance between them, asking whether the
#' indirect effect through repetition is zero and whether the two
#' specific indirect effects differ (their Research Questions 2 and 3).
#' That analysis is the single-mediator call of the examples with the
#' parallel model in place of the single one and the contrast supplied
#' through \code{hypotheses}:
#' \preformatted{parallel <- "
#'   imagery    ~ b1*instruction
#'   repetition ~ b3*instruction
#'   recall     ~ b2*imagery + b4*repetition + b5*instruction
#'   imagery ~~ repetition
#' "
#' mediation_mbco(parallel, S = S_tk, M = M_tk, N = 369,
#'                x = "instruction", y = "recall",
#'                hypotheses = c(imagery_minus_repetition =
#'                  "indirect_via_imagery - indirect_via_repetition"),
#'                ci_method = "monte_carlo", seed = 113)}
#' Six effects are reported, so six constrained null models are fit,
#' about twice the cost of the single-mediator analysis. The indirect
#' effect through repetition is near zero (the paper reports
#' \eqn{\mathrm{LRT}_{\mathrm{MBCO}} = 0.083}, \emph{p} = .773), while
#' the contrast shows the imagery pathway is larger (the paper reports
#' \eqn{\mathrm{LRT}_{\mathrm{MBCO}} = 25.828}, difference = 2.222,
#' SE = 0.445).
#'
#' \strong{Choosing the confidence interval.} The default profile
#' likelihood interval inverts the likelihood ratio test for the effect
#' itself (Neale & Miller, 1997), so its limits are free to sit
#' asymmetrically about the estimate, as the skewed sampling
#' distribution of a product of coefficients calls for, and the interval
#' and the MBCO test tell the same story. Each bound is a constrained
#' search of its own, which is what makes it the expensive choice. The
#' Monte Carlo interval is the inexpensive alternative that also
#' accommodates the skewness: \code{B} coefficient vectors are drawn
#' from the joint normal approximation of the estimates, the effect is
#' evaluated in each draw, and the limits are the empirical
#' \eqn{(\alpha/2, 1 - \alpha/2)} quantiles of those \code{B} values
#' (MacKinnon, Lockwood, & Williams, 2004). Ask for it when the profile
#' searches are slow, when a profile bound fails to converge, or when
#' comparing with a published analysis that reports one, as Tofighi and
#' Kelley (2020) do for the memory data in the examples.
#'
#' \strong{Monte Carlo draws.} The limits of each effect are quantiles
#' of its \code{B} drawn values, so the Monte Carlo interval is run only
#' when \code{B} puts at least 10 draws beyond each limit (400 at the 95
#' percent level). A draw in which an effect is undefined (a
#' division by a draw near zero, or a function outside its domain) is
#' dropped and counted, never replaced, and the interval is computed from
#' the rest. The result records \code{B} (requested), \code{B_used} (the
#' number the interval rests on), \code{n_dropped}, and \code{B_required}
#' (the fewest replications a request may ask for; that count is held to
#' \code{B}, not to \code{B_used}) as
#' attributes, and the line under the table states them, with
#' \code{bootstrap_status} (\code{"complete"}, \code{"dropped"},
#' \code{"exploratory"}, \code{"insufficient_count"}, or
#' \code{"insufficient_share"}). A call gives at most one such warning,
#' whatever the number of effects. A warning of class
#' \code{dmar_bootstrap_replications_dropped} states the counts when any
#' replication was dropped, with a caution when the dropped replications
#' are at least the 2.5 percent of \code{B} that lies beyond one limit at
#' the 95 percent level, because they could then decide by themselves where
#' that limit falls. When more than half were dropped, the limits of that
#' effect are \code{NA}, the
#' estimate is returned, and the warning has class
#' \code{dmar_bootstrap_insufficient}. The draws are counted effect by
#' effect: an effect whose interval is withheld does not withhold the
#' others, \code{B_used} is the fewest draws any effect's limits rest on,
#' and \code{bootstrap_status} is the least favorable status among the
#' effects. When several effects are affected the one warning lists them
#' (draws dropped, draws used, and whether the interval is reported),
#' gives in full the account of the effect with the fewest usable draws,
#' has the least favorable class among them, and carries the list as the
#' data frame \code{effects}.
#' \code{options(dmar.bootstrap_exploratory = TRUE)} lets a
#' smaller \code{B} (2 or more) run for a first look or a demonstration;
#' every such call gives a warning of class
#' \code{dmar_bootstrap_exploratory} and the result is not one to report.
#' See \code{\link{dmar_bootstrap}} for the rule and its reasons.
#'
#' The Wald
#' interval is the symmetric delta method interval and is reported only
#' for comparison.
#'
#' \strong{Multiple groups (moderated mediation across groups).} With
#' \code{group} (or list-form \code{S}, \code{M}, \code{N}), the model
#' is fit as a multiple-group SEM: the same structure in every group,
#' with group-specific parameters. Every enumerated effect is then
#' estimated within each group (its term carries the group label as a
#' suffix), and for each effect the difference from the reference group
#' (the first group label) is estimated and tested with the same
#' constrained-optimization machinery, since a difference of two
#' products is itself a nonlinear function of the parameters. That
#' difference test is moderated mediation with a categorical
#' moderator: "does the indirect effect differ across groups?" is
#' exactly the between-group contrast of the conditional indirect
#' effects. lavaan's labeling rules decide what varies: a single label
#' such as \code{b1} on a path applies to every group and therefore
#' imposes cross-group equality (the corresponding difference is
#' identically zero and its row is dropped); to let a path differ by
#' group, leave it unlabeled or give per-group labels with the vector
#' form \code{c("b1_f", "b1_m")*x}. With more than two groups, each
#' non-reference group is compared with the reference group. The call
#' below fits the simple model in both groups of a data frame
#' \code{d} whose \code{condition} column takes two values, and reports
#' every effect per group beside its between-group difference; leaving
#' the paths unlabeled lets them differ by group.
#' \preformatted{mediation_mbco("m ~ x
#'                 y ~ m + x",
#'                data = d, group = "condition", x = "x", y = "y",
#'                ci_method = "wald")}
#' Each reported row costs its own constrained null model fit, so this
#' two-group analysis fits nine null models where the single-group
#' analysis fits three.
#'
#' \strong{Moderated mediation, probed.} With \code{moderator}, every
#' regression coefficient along a pathway becomes a linear function of
#' the moderator wherever the model contains the matching interaction
#' (stated as \code{x:w} in the syntax, or as a product column among
#' the predictors). A pathway effect, the product of its edge
#' coefficients, is then a polynomial in the moderator, and the
#' function derives that polynomial symbolically. Three kinds of rows
#' follow for each moderated effect. First, the conditional effect at
#' each probe value, each with its own confidence interval and MBCO
#' test. Second, the polynomial's moderator coefficients: for a
#' pathway moderated in one place the single such coefficient is
#' exactly the index of moderated mediation (Hayes, 2015), here tested
#' by likelihood ratio rather than bootstrap; a pathway moderated in
#' several places gets one row per power of the moderator. Third, for
#' a pathway moderated in several places (so the conditional effect is
#' curved in the moderator), a joint constancy test of all moderator
#' coefficients at once, with as many degrees of freedom as
#' constraints, asking whether the pathway effect depends on the
#' moderator at all; as a joint test its row reports no scalar
#' estimate. Unmoderated effects in the same model keep their single
#' rows. The null set of a conditional-effect constraint is again a
#' union of branches (one edge's conditional coefficient or another's
#' must vanish at the probed value), and the null fits are started on
#' each branch, and at the model with all interactions removed for the
#' moderation tests, so the reported statistics reflect the
#' best-fitting null models. For example:
#' \preformatted{mediation_mbco("m ~ x + w + x:w
#'                 y ~ m + x + w",
#'                data = d, x = "x", y = "y", moderator = "w")}
#' \code{\link{plot_mediation_mbco}} draws the same conditional
#' effects as curves over the moderator's whole range with a
#' confidence band, the visual companion to the probe rows. Bespoke
#' conditional quantities can still be written directly with
#' \code{:=} definitions or \code{hypotheses} when the built-in
#' probing does not cover them.
#'
#' \strong{Moderated mediation and mediated moderation.} Moderated
#' mediation asks whether an indirect effect depends on a moderator
#' \eqn{W}: the conditional indirect effect varies with \eqn{w}
#' (Muller, Judd, & Yzerbyt, 2005; Preacher, Rucker, & Hayes, 2007).
#' Mediated moderation asks whether an observed \eqn{X \times W}
#' interaction on \eqn{Y} is transmitted through the mediator. The two
#' share their algebra: with first-stage moderation
#' \eqn{M = a_1 X + a_2 W + a_3 XW + e_M} and \eqn{Y = b M + \ldots},
#' the quantity \eqn{a_3 b} is at once the index of moderated mediation
#' (Hayes, 2015) and the indirect effect of the product term through
#' \eqn{M}; what differs is the question and the reporting emphasis.
#' This function reports the conditional-indirect-effect framing: the
#' probe rows and the moderation rows above. For a categorical
#' moderator, use \code{group}.
#'
#' \strong{Refusals, warnings, and judgment calls.} The function stops
#' where a computed number would be mislabeled: a nonrecursive
#' (feedback) regression structure, categorical-endogenous syntax
#' (thresholds), explicit \code{==} constraints, a model with no
#' indirect pathway between \code{x} and \code{y}. It stops as well, with
#' an error of class \code{dmar_no_maximum}, when the \pkg{lavaan} fit
#' that the \pkg{OpenMx} fits start from has no maximum, an error variance
#' running off below zero while the likelihood is flat or still rising
#' (see \code{\link{fa_bifactor}}), printing the \pkg{lavaan} code of that
#' fit, or none for a fit given as \code{model}. It warns where the
#' data make trouble detectable: an endogenous variable with only two
#' distinct values (a binary mediator or outcome, where the product of
#' coefficients is not the causal indirect effect), an interaction
#' term among the predictors that no declared \code{moderator}
#' accounts for (declaring the moderator resolves the warning by
#' probing the moderation), an interaction included without its
#' matching main effect (the principle of marginality), null models
#' that converge imperfectly, profile bounds that fail. What it cannot
#' check is left to the analyst, and stated rather than assumed: the
#' no omitted confounder assumption, linearity, and the causal
#' direction of the arrows come from the design, not from the fit.
#'
#' \strong{Information criteria.} The columns \code{delta_aic} and
#' \code{delta_bic} report AIC and BIC for the null model minus the same
#' criterion for the full model; positive values favor the full model
#' (the effect improves fit by more than the parsimony penalty). A
#' scalar equality constraint reduces the effective number of free
#' parameters by one, so the differences equal
#' \eqn{\mathrm{LRT}_{\mathrm{MBCO}} - 2} for the AIC and
#' \eqn{\mathrm{LRT}_{\mathrm{MBCO}} - \log N} for the BIC; for a joint
#' test the 2 and the \eqn{\log N} are multiplied by its degrees of
#' freedom. (The
#' \pkg{OpenMx} summary counts the same number of estimated parameters
#' in both models, so its printed AIC difference equals the likelihood
#' ratio statistic; DMAR counts the constraint against the null model,
#' consistent with the degrees of freedom of the test.)
#'
#' \strong{Constraints that add no restriction.} The degrees of freedom of
#' a test are the number of its constraints, which is the number of
#' restrictions only when each constraint restricts the parameters beyond
#' what the others do. Before any \pkg{OpenMx} model is fit, the
#' constraints of every test, written in the parameter labels, are checked
#' twice. First, the function forms the Jacobian of the constraints (the
#' matrix of their derivatives with respect to the free parameters) at the
#' unrestricted estimates of the \pkg{lavaan} fit, which lie off the null
#' set, and reads its rank after scaling every column by its largest entry
#' and every row to unit length, so that the rank does not depend on the
#' units of the parameters or the scale of the constraints. A singular
#' value below \eqn{10^{-7}} times the largest counts as zero (the same
#' number as the default \code{tol} of \code{\link[base]{qr}}). The rank is
#' read again at a fixed point near the estimates, so that an estimate of
#' exactly zero cannot make independent constraints look dependent, and
#' the larger rank is used. A test whose constraints have a rank below
#' their number (a constraint that is a function of the others) is refused
#' (see \code{hypotheses}); so is a single constraint whose gradient is
#' zero, which restricts nothing. The rank is not read on the null set: at
#' the crossing of the branches of a product constraint its gradient
#' vanishes, and a rank read there would call every set that contains it
#' dependent.
#'
#' Second, constraints whose Jacobian has full rank can still restrict the
#' parameters in fewer directions than their number, when one of them is
#' zero wherever the others are: \code{c("a*b", "b")} has the null set
#' \eqn{b = 0}, \code{c("a", "a*cp")} the null set \eqn{a = 0}, and
#' \code{c("a*b", "a*b*cp")} the null set \eqn{ab = 0}, one restriction
#' each. For every constraint of a joint test the function finds points of
#' the null set of the others, by Gauss-Newton steps from fixed starting
#' points that cover both signs and a range of magnitudes of the parameters
#' (in coordinates scaled by the estimates, and again in coordinates scaled
#' by the derivatives of the constraints, so that neither the units of the
#' parameters nor the starting points of one scaling decide the result),
#' and refuses the test when the constraint is zero at all of them, to
#' within \eqn{10^{-6}} of its typical size, under both scalings. No random
#' numbers are drawn. A set such as \code{c("a*b", "a*cp")}, whose null set
#' is the union of \eqn{a = 0} (one restriction) and \eqn{b = cp = 0}
#' (two), is accepted with 2 degrees of freedom: the chi square reference
#' on 2 degrees of freedom is right where the second branch holds, and the
#' test is conservative where only \eqn{a = 0} holds, since the statistic
#' there is asymptotically chi square on 1 degree of freedom, which exceeds
#' the .95 quantile of the chi square distribution on 2 degrees of freedom
#' with probability .0144. The checks do not detect a set none of whose
#' constraints is zero wherever the others are but whose null set still
#' imposes fewer restrictions than there are constraints, such as
#' \code{c("a*b", "a*cp", "b*cp")} (any two of the three coefficients zero,
#' two restrictions); its test is conservative. Neither check establishes
#' that the null hypothesis is regular (see \emph{The reference
#' distribution where the branches cross} above).
#'
#' \strong{Change in explained variance.} Following the reporting
#' recommendation of Tofighi and Kelley (2020), the function computes
#' \eqn{R^2} for every endogenous variable under the full model (the
#' \code{"R2"} attribute) and the drop in each \eqn{R^2} under every
#' null model (the \code{"delta_R2"} attribute, variables by tested
#' effects), so the fit cost of removing an effect can be read as a
#' change in effect size and not only as a test statistic.
#'
#' \strong{Summary statistics input.} With complete data the multivariate
#' normal log likelihood depends on the data only through the sample
#' means and the sample covariance matrix. When \code{S} (and optionally
#' \code{M}) is supplied, the function therefore constructs an internal
#' data set with exactly those moments (via \code{\link[MASS]{mvrnorm}}
#' with \code{empirical = TRUE}) and proceeds as with raw data; the
#' estimates, likelihood ratio tests, and confidence intervals are
#' identical to what the raw data would give, whatever internal data set
#' realizes the moments. \code{S} is treated as the unbiased
#' (divisor \eqn{N - 1}) covariance matrix, the form reported in
#' articles. This is how a published mediation analysis can be
#' reproduced, and its hypotheses re-tested with the MBCO procedure,
#' from a table of descriptive statistics alone.
#'
#' \strong{Assumptions.} The causal reading of any mediation analysis
#' rests on the no omitted confounder assumption for the predictor to
#' mediator and mediator to outcome relations, in addition to the usual
#' distributional assumptions (residuals multivariate normal, linear
#' relations, no treatment by mediator interaction); randomizing
#' \eqn{X} supports the first link but not the second (MacKinnon, 2008;
#' Tofighi & Kelley, 2020). With a binary randomized \eqn{X} the
#' variable enters the model as numeric 0/1, its exogenous variance
#' freely estimated; the normality assumption applies to the residuals
#' of the endogenous variables.
#'
#' @return A \code{data.frame} (classes \code{dmar_mediation_mbco},
#'   \code{dmar_tbl}) with one row per effect and columns
#'   \code{pathway} (the traced pathway or defining expression),
#'   \code{term}, \code{estimate}, \code{se} (delta method),
#'   \code{ci_lower}, \code{ci_upper}, \code{lrt} (the MBCO likelihood
#'   ratio statistic), \code{df}, \code{p_value}, \code{delta_aic}, and
#'   \code{delta_bic}. A joint test row (a several-constraint
#'   moderation or hypothesis test) reports the test columns with
#'   \code{df} equal to the number of constraints, none of which is
#'   redundant (see \code{hypotheses}); its \code{estimate},
#'   \code{se}, and interval are \code{NA} because no single number
#'   summarizes several constraints. Attributes: \code{"conf_level"},
#'   \code{"ci_method"}, with the Monte Carlo interval \code{"B"},
#'   \code{"B_used"}, \code{"n_dropped"}, \code{"B_required"}, and
#'   \code{"bootstrap_status"} (see \emph{Monte Carlo draws} under
#'   Details), \code{"optimizer"}, \code{"x"}, \code{"y"},
#'   \code{"N"}, \code{"deviance"}, \code{"aic"}, \code{"bic"},
#'   \code{"n_par"} (full model fit information), \code{"groups"} (the
#'   group labels, reference first; multiple-group fits only),
#'   \code{"R2"} (named vector, endogenous variables under the full
#'   model, per group when grouped), \code{"delta_R2"} (matrix of full
#'   minus null \eqn{R^2}, variables by tested effects),
#'   \code{"moderation"} (for a moderated analysis: the moderator
#'   name, the probe values, and each moderated effect's moderator
#'   polynomial, which is what \code{\link{plot_mediation_mbco}}
#'   draws), \code{"rescale"} (only when badly scaled variables were
#'   rescaled for the fits, see \emph{Badly scaled variables} under Details:
#'   \code{divisors}, the power of 2 each observed variable was
#'   divided by, and \code{multipliers}, the factor that turns each
#'   parameter of the rescaled fit into the parameter of the original
#'   metric), and \code{"mx_model"} (the fitted \pkg{OpenMx} full model, in
#'   the rescaled metric when there is a \code{"rescale"} attribute, an
#'   escape hatch for further OpenMx work). Use
#'   \code{tidy()} and \code{glance()} for broom-style views.
#'
#' @references
#' Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
#'   analysis: Introducing the model-based constrained optimization
#'   procedure. \emph{Psychological Methods, 25}(4), 496--515.
#'   \doi{10.1037/met0000259}
#'
#' Tofighi, D., & Kelley, K. (2020). Indirect effects in sequential
#'   mediation models: Evaluating methods for hypothesis testing and
#'   confidence interval formation. \emph{Multivariate Behavioral
#'   Research, 55}(2), 188--210. \doi{10.1080/00273171.2019.1618545}
#'
#' Hayes, A. F. (2015). An index and test of linear moderated
#'   mediation. \emph{Multivariate Behavioral Research, 50}(1), 1--22.
#'   \doi{10.1080/00273171.2014.962683}
#'
#' MacKinnon, D. P. (2008). \emph{Introduction to statistical mediation
#'   analysis}. Erlbaum.
#'
#' MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
#'   limits for the indirect effect: Distribution of the product and
#'   resampling methods. \emph{Multivariate Behavioral Research, 39}(1),
#'   99--128. \doi{10.1207/s15327906mbr3901_4}
#'
#' MacKinnon, D. P., Valente, M. J., & Wurpts, I. C. (2018). Benchmark
#'   validation of statistical models: Application to mediation analysis
#'   of imagery and memory. \emph{Psychological Methods, 23}(4),
#'   654--671. \doi{10.1037/met0000174}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Muller, D., Judd, C. M., & Yzerbyt, V. Y. (2005). When moderation is
#'   mediated and mediation is moderated. \emph{Journal of Personality
#'   and Social Psychology, 89}(6), 852--863.
#'   \doi{10.1037/0022-3514.89.6.852}
#'
#' Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick,
#'   T. R., Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes,
#'   H. H., & Boker, S. M. (2016). OpenMx 2.0: Extended structural
#'   equation and statistical modeling. \emph{Psychometrika, 81}(2),
#'   535--549. \doi{10.1007/s11336-014-9435-8}
#'
#' Neale, M. C., & Miller, M. B. (1997). The use of likelihood-based
#'   confidence intervals in genetic models. \emph{Behavior Genetics,
#'   27}(2), 113--120. \doi{10.1023/A:1025681223921}
#'
#' Preacher, K. J., Rucker, D. D., & Hayes, A. F. (2007). Addressing
#'   moderated mediation hypotheses: Theory, methods, and
#'   prescriptions. \emph{Multivariate Behavioral Research, 42}(1),
#'   185--227. \doi{10.1080/00273170701341316}
#'
#' Wilks, S. S. (1938). The large-sample distribution of the likelihood
#'   ratio for testing composite hypotheses. \emph{Annals of
#'   Mathematical Statistics, 9}(1), 60--62.
#'   \doi{10.1214/aoms/1177732360}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{plot_mediation_mbco}} for the
#'   conditional-effect display of a moderated analysis;
#'   \code{\link{mediate}} for the regression-based simple mediation
#'   model with bootstrap intervals;
#'   \code{\link{ss_aipe_indirect_effect}} and
#'   \code{\link{ss_power_indirect_effect}} for planning;
#'   \code{\link{var_indirect_effect}} for the delta method variance.
#'
#' @family mediation
#'
#' @keywords models
#'
#' @examples
#' # Replicate the memory experiment analyses of Tofighi and Kelley
#' # (2020) from the summary statistics in their Table 1 (data from
#' # MacKinnon, Valente, & Wurpts, 2018; N = 369). Instruction is 1 for
#' # imagery rehearsal instructions and 0 for repetition instructions.
#' # Requires the OpenMx and lavaan packages to be installed.
#' vars <- c("instruction", "imagery", "repetition", "recall")
#' sds  <- c(0.50, 2.96, 2.84, 3.40)
#' R_tk <- matrix(c(1.00,  .62, -.67,  .32,
#'                   .62, 1.00, -.56,  .51,
#'                  -.67, -.56, 1.00, -.28,
#'                   .32,  .51, -.28, 1.00), 4, 4,
#'                dimnames = list(vars, vars))
#' S_tk <- outer(sds, sds) * R_tk
#' M_tk <- c(instruction = 0.51, imagery = 5.66, repetition = 6.08,
#'           recall = 12.07)
#'
#' # Single-mediator model: instruction -> imagery -> recall. Every
#' # reported effect costs its own constrained null model fit in OpenMx,
#' # which is where the run time goes; the Monte Carlo interval itself
#' # is inexpensive at any B. The indirect effect is about 2.1 words
#' # (the paper reports 2.121, SE = 0.276, 95% Monte Carlo CI
#' # [1.600, 2.682]), and its likelihood ratio statistic of 71.31 is the
#' # value discussed under Details.
#' single <- "
#'   imagery ~ b1*instruction
#'   recall  ~ b2*imagery + b3*instruction
#' "
#' mediation_mbco(single, S = S_tk, M = M_tk, N = 369,
#'                x = "instruction", y = "recall",
#'                ci_method = "monte_carlo", seed = 113)
#'
#' # The parallel two-mediator model of the same paper, with the contrast
#' # of its two specific indirect effects, is described under Details;
#' # it is fit the same way from the same summary statistics, adding the
#' # 'hypotheses' argument. Raw data go in through 'data' rather than
#' # 'S', 'M', and 'N'; adding 'group' fits the model in every group and
#' # tests the between-group difference of each effect, which is
#' # moderated mediation with a categorical moderator (see Details). A
#' # continuous moderator goes in through 'moderator', and the help page
#' # for plot_mediation_mbco fits that analysis from a data frame and
#' # draws the conditional effects it estimates.
#'
#' @export
#' @importFrom stats coef vcov qnorm pchisq quantile setNames
#' @importFrom utils capture.output combn
mediation_mbco <- function(model, data = NULL, S = NULL, M = NULL,
                           N = NULL, group = NULL, x = NULL, y = NULL,
                           moderator = NULL, probe_values = NULL,
                           hypotheses = NULL,
                           ci_method = c("profile_likelihood",
                                         "monte_carlo", "wald"),
                           conf_level = 0.95, B = 10000,
                           optimizer = c("SLSQP", "CSOLNP", "NPSOL"),
                           seed = NULL) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    .dmar_refuse(c(
      "The package 'lavaan' is needed; please install the package ",
      "and try again."), "package_required", package = "lavaan")
  }
  if (!requireNamespace("OpenMx", quietly = TRUE)) {
    .dmar_refuse(c(
      "The package 'OpenMx' is needed; please install the package ",
      "and try again."), "package_required", package = "OpenMx")
  }
  ci_method <- match.arg(ci_method)
  optimizer <- match.arg(optimizer)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .dmar_refuse_argument(
      "'conf_level' must be a single number in (0, 1).", arg = "conf_level")
  }
  # B is read only by the Monte Carlo interval, and is checked only when
  # that interval is requested, by the rule of ?dmar_bootstrap: 10 draws
  # beyond each limit (400 at the 95 percent level). Through DMAR 1.0.x a
  # B of at least 100 was required whatever the interval.
  if (ci_method == "monte_carlo") {
    B <- .dmar_check_B(B, tail = (1 - conf_level) / 2,
                       conf_level = conf_level)
  }
  # A supplied seed is local to this call (covering the Monte Carlo
  # interval and any optimizer restarts); the caller's random number
  # generator state is restored on exit.
  .dmar_local_seed(seed)
  if (!optimizer %in% OpenMx::mxAvailableOptimizers()) {
    .dmar_refuse_argument(c(
      "The '", optimizer, "' optimizer is not available in this ",
      "build of OpenMx. The CRAN build ships 'SLSQP' and 'CSOLNP'; ",
      "'NPSOL' requires the build distributed from the OpenMx ",
      "website."), arg = NULL)
  }

  # ---- resolve the data --------------------------------------------------
  lav_input <- inherits(model, "lavaan")
  # Whether the data are built from summary statistics, for the lavaan code
  # of a refusal (.mbco_no_max_code()).
  from_S <- is.null(data) && !is.null(S)
  products <- character(0)
  if (lav_input && is.null(data) && is.null(S)) {
    raw <- try(lavaan::lavInspect(model, "data"), silent = TRUE)
    if (!inherits(raw, "try-error") &&
        is.list(raw) && !is.data.frame(raw)) {
      # A fitted multiple-group lavaan object returns one matrix per
      # group; rebuild the stacked data and recover the group variable.
      group <- lavaan::lavInspect(model, "group")
      gl <- as.character(lavaan::lavInspect(model, "group.label"))
      pieces <- lapply(seq_along(raw), function(g) {
        dg <- as.data.frame(raw[[g]])
        dg[[group]] <- gl[g]
        dg
      })
      data <- do.call(rbind, pieces)
      data[[group]] <- factor(data[[group]], levels = gl)
    } else if (!inherits(raw, "try-error") && !is.null(dim(raw))) {
      data <- as.data.frame(raw)
    } else {
      .dmar_refuse_argument(c(
        "The raw data could not be recovered from the fitted lavaan ",
        "object; supply 'data' (or 'S' and 'N') directly."),
        arg = c("data", "S", "N"))
    }
  }
  if (!is.null(data) && !is.null(S)) {
    .dmar_refuse_argument(
      "Supply either 'data' or 'S' (with 'N'), not both.",
      arg = c("data", "S", "N"))
  }
  if (is.null(data)) {
    if (is.null(S)) {
      .dmar_refuse_argument(c(
        "Supply raw data via 'data' or a covariance matrix via 'S' ",
        "(with 'N')."), arg = c("data", "S", "N"))
    }
    if (is.list(S) && !is.data.frame(S)) {
      # Multiple groups from summary statistics: one covariance matrix
      # (and optionally one mean vector) per group.
      gl <- names(S)
      if (is.null(gl) || any(!nzchar(gl))) {
        .dmar_refuse_argument(
          "'S' as a list must have one named element per group.", arg = "S")
      }
      if (is.null(N) || length(N) != length(S)) {
        .dmar_refuse_argument(
          "'N' must give one sample size per group in 'S'.", arg = c("N", "S"))
      }
      if (!is.null(names(N))) {
        if (!all(gl %in% names(N))) {
          .dmar_refuse_argument(
            "The names of 'N' must match the names of 'S'.", arg = c("N", "S"))
        }
        N <- N[gl]
      }
      if (!is.null(M) && (!is.list(M) || !all(gl %in% names(M)))) {
        .dmar_refuse_argument(c(
          "'M' must be a named list with one mean vector per group ",
          "in 'S'."), arg = c("M", "S"))
      }
      if (is.null(group)) group <- "group"
      pieces <- lapply(seq_along(gl), function(g) {
        dg <- .mbco_surrogate(S[[g]],
                              if (is.null(M)) NULL else M[[gl[g]]],
                              N[[g]])
        dg[[group]] <- gl[g]
        dg
      })
      data <- do.call(rbind, pieces)
      data[[group]] <- factor(data[[group]], levels = gl)
    } else {
      data <- .mbco_surrogate(S, M, N)
    }
  }
  data <- as.data.frame(data)
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L ||
        !group %in% names(data)) {
      .dmar_refuse_argument(
        "'group' must name one column of the data.", arg = "group")
    }
  }
  if (!is.null(moderator) && !is.null(group)) {
    .dmar_refuse_argument(c(
      "'moderator' and 'group' cannot yet be combined; probe one ",
      "moderator at a time."), arg = c("moderator", "group"))
  }
  if (!is.null(probe_values) && is.null(moderator)) {
    .dmar_refuse_argument(
      "'probe_values' requires 'moderator'.",
      arg = c("probe_values", "moderator"))
  }

  # Interaction syntax: rewrite each 'a:b' term into a product column
  # so the model can state the moderation directly.
  if (is.character(model)) {
    prep <- .mbco_product_terms(model, data)
    model <- prep$model
    data <- prep$data
    products <- prep$products
  }
  if (!is.null(moderator)) {
    if (!is.character(moderator) || length(moderator) != 1L ||
        !moderator %in% names(data)) {
      .dmar_refuse_argument(
        "'moderator' must name one column of the data.", arg = "moderator")
    }
    if (!is.numeric(data[[moderator]])) {
      .dmar_refuse_argument(c(
        "The moderator must be numeric; code a categorical ",
        "moderator through 'group' instead."), arg = "group")
    }
  }

  # ---- lavaan: parse (and pre-fit) the model -----------------------------
  if (lav_input) {
    lav <- model
    if (is.null(group) && lavaan::lavInspect(lav, "ngroups") > 1L) {
      group <- lavaan::lavInspect(lav, "group")
    }
  } else {
    if (!is.character(model) || length(model) != 1L) {
      .dmar_refuse_argument(c(
        "'model' must be a lavaan model syntax string or a fitted ",
        "lavaan object."), arg = "model")
    }
    miss <- if (anyNA(data)) "ml" else "listwise"
    lav <- try(lavaan::sem(model, data = data, meanstructure = TRUE,
                           fixed.x = FALSE, missing = miss, se = "none",
                           group = group, warn = FALSE),
               silent = TRUE)
    if (inherits(lav, "try-error")) {
      .dmar_refuse_argument(c(
        "lavaan could not process 'model': ",
        attr(lav, "condition")$message), arg = "model")
    }
  }
  ngroups <- lavaan::lavInspect(lav, "ngroups")
  glabels <- if (ngroups > 1L) {
    as.character(lavaan::lavInspect(lav, "group.label"))
  } else {
    character(0)
  }
  pt <- lavaan::parTable(lav)
  bad_ops <- setdiff(unique(pt$op), c("=~", "~", "~~", "~1", ":=", "=="))
  if (length(bad_ops) > 0L) {
    .dmar_refuse_argument(c(
      "The model uses operator(s) not supported by mediation_mbco(): ",
      paste(sQuote(bad_ops), collapse = ", "), "."), arg = NULL)
  }

  # One label per parameter; user labels are kept (shared labels keep
  # their meaning as equality constraints), unlabeled parameters get
  # generated names so effects can reference them in algebra.
  labels <- pt$label
  auto <- !nzchar(labels) & pt$op != ":="
  labels[auto] <- paste0("prm", pt$id[auto])
  user_labels <- unique(labels[!auto & pt$op != ":="])
  not_names <- user_labels[user_labels != make.names(user_labels)]
  if (length(not_names) > 0L) {
    .dmar_refuse_argument(c(
      "Parameter labels must be valid R names; rename: ",
      paste(sQuote(not_names), collapse = ", "), "."), arg = NULL)
  }
  pt$mx_label <- labels

  # lavaan records a shared label as an explicit "==" row pairing the
  # two parameters. Sharing the label in OpenMx enforces the equality
  # already, so such rows are dropped; any other "==" constraint has no
  # shared-label translation and is refused.
  if ("==" %in% pt$op) {
    plab_map <- stats::setNames(
      pt$mx_label[pt$op != "=="], pt$plabel[pt$op != "=="])
    eqs <- pt[pt$op == "==", , drop = FALSE]
    redundant <- eqs$lhs %in% names(plab_map) &
      eqs$rhs %in% names(plab_map) &
      plab_map[eqs$lhs] == plab_map[eqs$rhs]
    if (!all(redundant)) {
      .dmar_refuse_argument(c(
        "Explicit '==' constraints are not supported; equate ",
        "parameters by giving them the same label."), arg = NULL)
    }
    pt <- pt[pt$op != "==", , drop = FALSE]
  }

  # Rows sharing a label are one parameter (an equality constraint);
  # OpenMx requires them to carry bitwise-identical starting values,
  # while lavaan's constrained estimates agree only to solver
  # precision. Harmonize on the first row's value.
  par_rows <- which(pt$op != ":=")
  dup_labs <- unique(pt$mx_label[par_rows][
    duplicated(pt$mx_label[par_rows])])
  for (lb in dup_labs) {
    rows <- par_rows[pt$mx_label[par_rows] == lb]
    pt$est[rows] <- pt$est[rows[1L]]
  }

  latents  <- unique(pt$lhs[pt$op == "=~"])
  manifest <- setdiff(unique(c(pt$lhs[pt$op %in% c("~", "~~", "~1", "=~")],
                               pt$rhs[pt$op %in% c("~", "~~", "=~")])),
                      c(latents, ""))
  missing_vars <- setdiff(manifest, names(data))
  if (length(missing_vars) > 0L) {
    .dmar_refuse_argument(c(
      "Variable(s) not found in the data: ",
      paste(sQuote(missing_vars), collapse = ", "), "."), arg = NULL)
  }
  if (!all(vapply(data[manifest], is.numeric, logical(1L)))) {
    .dmar_refuse_argument(c(
      "All modeled variables must be numeric; convert factors ",
      "before calling mediation_mbco()."), arg = NULL)
  }

  # ---- guardrails on the data ------------------------------------------
  # An endogenous variable with two distinct values is fit as if
  # continuous; the product-of-coefficients estimand is then not the
  # causal indirect effect (see the paper's limitations).
  endo_obs <- intersect(unique(pt$lhs[pt$op == "~"]), manifest)
  two_level <- endo_obs[vapply(endo_obs, function(v) {
    vals <- data[[v]]
    length(unique(vals[is.finite(vals)])) <= 2L
  }, logical(1L))]
  if (length(two_level) > 0L) {
    .dmar_warn(c(
      "Endogenous variable(s) ",
      paste(sQuote(two_level), collapse = ", "),
      " take only two distinct values. The linear normal model ",
      "treats them as continuous, and for a binary mediator or ",
      "outcome the product of coefficients is not the causal ",
      "indirect effect; potential outcome methods are needed for ",
      "that reading."), "result_qualified")
  }

  # ---- effects: enumerate pathways from x to y ---------------------------
  reg <- pt[pt$op == "~", , drop = FALSE]
  if (nrow(reg) == 0L) {
    .dmar_refuse_argument(c(
      "The model contains no regression ('~') paths; there is ",
      "nothing to mediate."), arg = NULL)
  }
  mod <- NULL
  if (!is.null(moderator)) {
    mod <- .mbco_moderation(reg, data, manifest, moderator,
                            probe_values)
    if (is.null(mod)) moderator <- NULL
  }
  defined  <- pt[pt$op == ":=", , drop = FALSE]
  eff <- .mbco_effects(reg, x = x, y = y,
                       have_extra = nrow(defined) > 0L ||
                         length(hypotheses) > 0L,
                       glabels = glabels, mod = mod)
  x <- eff$x
  y <- eff$y

  # A predictor that is the elementwise product of two other modeled
  # variables is an interaction term: the enumerated pathway effects
  # ignore that moderation, so point to the conditional-effect recipe.
  if (!is.null(x)) {
    skip_v <- if (is.null(mod)) character(0) else mod$moderator
    skip_p <- if (is.null(mod)) character(0) else mod$prods
    pred_obs <- setdiff(intersect(unique(reg$rhs), manifest), skip_p)
    hits <- character(0)
    pairs <- if (length(manifest) >= 3L) {
      utils::combn(manifest, 2L, simplify = FALSE)
    } else {
      list()
    }
    for (v in pred_obs) {
      vv <- data[[v]]
      for (pr in pairs) {
        if (v %in% pr || any(pr %in% skip_v)) next
        a <- data[[pr[1L]]]
        b <- data[[pr[2L]]]
        ok <- is.finite(a) & is.finite(b) & is.finite(vv)
        if (sum(ok) < 3L) next
        if (max(abs(vv[ok] - a[ok] * b[ok])) <
            1e-8 * max(1, max(abs(vv[ok])))) {
          hits <- c(hits, paste0(sQuote(v), " = ", sQuote(pr[1L]),
                                 " * ", sQuote(pr[2L])))
        }
      }
    }
    if (length(hits) > 0L) {
      .dmar_warn(c(
        "Interaction term(s) detected among the predictors (",
        paste(unique(hits), collapse = "; "),
        "). The automatically enumerated effects average over ",
        "that moderation; declare the moderator through the ",
        "'moderator' argument to probe it, or define conditional ",
        "effects with ':=' or 'hypotheses' (see ",
        "?mediation_mbco)."), "result_qualified")
    }
  }

  # User-defined quantities and extra hypotheses join the effect list.
  # Each row's 'test' column names the algebra(s) its null model
  # constrains to zero; 'df' is the number of constraints. A row with
  # several constraints (a joint test) has no scalar estimate.
  hidden <- data.frame(name = character(0), expr = character(0),
                       stringsAsFactors = FALSE)
  user_terms <- character(0)  # the rows that come from 'hypotheses'
  if (nrow(defined) > 0L) {
    eff$table <- rbind(eff$table, data.frame(
      term = defined$label, pathway = defined$rhs, expr = defined$rhs,
      test = defined$label, df = 1L, stringsAsFactors = FALSE))
  }
  if (length(hypotheses) > 0L) {
    if (is.character(hypotheses)) hypotheses <- as.list(hypotheses)
    if (!is.list(hypotheses) ||
        !all(vapply(hypotheses, is.character, logical(1L))) ||
        !all(lengths(hypotheses) >= 1L)) {
      .dmar_refuse_argument(c(
        "'hypotheses' must be a named character vector, or a named ",
        "list of character vectors for joint tests."), arg = "hypotheses")
    }
    hnames <- names(hypotheses)
    if (is.null(hnames)) hnames <- rep("", length(hypotheses))
    hnames[!nzchar(hnames)] <-
      paste0("hypothesis_", seq_len(sum(!nzchar(hnames))))
    user_terms <- hnames
    for (h in seq_along(hypotheses)) {
      exprs <- hypotheses[[h]]
      if (length(exprs) == 1L) {
        eff$table <- rbind(eff$table, data.frame(
          term = hnames[h], pathway = unname(exprs),
          expr = unname(exprs), test = hnames[h], df = 1L,
          stringsAsFactors = FALSE))
      } else {
        members <- paste0(hnames[h], "_h", seq_along(exprs))
        hidden <- rbind(hidden, data.frame(
          name = members, expr = unname(exprs),
          stringsAsFactors = FALSE))
        eff$table <- rbind(eff$table, data.frame(
          term = hnames[h],
          pathway = paste0("joint: ",
                           paste(exprs, "= 0", collapse = ", ")),
          expr = NA_character_,
          test = paste(members, collapse = ";"),
          df = length(exprs), stringsAsFactors = FALSE))
      }
    }
  }
  effects <- eff$table
  if (anyDuplicated(effects$term)) {
    .dmar_refuse_argument(c(
      "Effect names collide: ",
      paste(sQuote(effects$term[duplicated(effects$term)]),
            collapse = ", "),
      ". Rename the clashing ':=' definition or hypothesis."), arg = NULL)
  }
  alg_defs <- rbind(
    data.frame(name = effects$term[!is.na(effects$expr)],
               expr = effects$expr[!is.na(effects$expr)],
               stringsAsFactors = FALSE),
    hidden)
  clash <- intersect(alg_defs$name, pt$mx_label[pt$op != ":="])
  if (length(clash) > 0L) {
    .dmar_refuse_argument(c(
      "Name(s) used both as a parameter label and as an effect: ",
      paste(sQuote(clash), collapse = ", "), "."), arg = NULL)
  }

  # Expand every expression down to parameter labels (":=" definitions
  # and effect names may reference each other), for the Monte Carlo
  # interval and for constraint-aware starting values.
  defs <- stats::setNames(
    lapply(alg_defs$expr, function(e) parse(text = e)[[1L]]),
    alg_defs$name)
  free_labels <- unique(pt$mx_label[pt$free > 0])
  expanded <- vapply(alg_defs$name, function(nm) {
    deparse1(.mbco_expand(defs[[nm]], defs))
  }, character(1L))
  for (i in seq_along(expanded)) {
    unknown <- setdiff(all.vars(parse(text = expanded[i])[[1L]]),
                       free_labels)
    if (length(unknown) > 0L) {
      .dmar_refuse_argument(c(
        "Effect ", sQuote(alg_defs$name[i]),
        " references unknown parameter(s): ",
        paste(sQuote(unknown), collapse = ", "), "."), arg = NULL)
    }
  }

  # ---- a likelihood with no maximum --------------------------------------
  # The OpenMx fits below start from the estimates of the lavaan fit, of
  # the same model and data. When that fit has no maximum, an error
  # variance running off below zero along a ridge with the likelihood flat
  # or still rising (.lavaan_no_maximum(), R/no_maximum_internals.R), the
  # likelihood has none for OpenMx either, and every effect, test, and
  # interval would be read from wherever an optimizer stopped; the call is
  # refused before any OpenMx fit, with the lavaan code of DMAR's fit
  # (.mbco_no_max_code()). A lavaan fit given as 'model' is checked as it
  # stands and refused without code, since there is no DMAR syntax to
  # print, as ci_bifactor_indices() refuses one. The OpenMx model carries
  # no bounds, so lavaan's bounds = "pos.var" is not offered as a choice
  # here. A finite negative estimate, at a maximum, goes on as before.
  no_max <- .lavaan_no_maximum(lav)
  if (!is.null(no_max)) {
    code <- if (!lav_input) {
      .mbco_no_max_code(model, data, S, M, N, group, miss, from_S,
                        products)
    }
    .no_max_refuse(
      "mediation_mbco", no_max, syntax = if (!lav_input) model,
      lavaan_call = if (!lav_input) {
        .no_max_call_text("lavaan::sem", code$args, select = code$select)
      },
      where = if (lav_input) {
        "the lavaan fit given as model"
      } else {
        "the lavaan fit the OpenMx fits start from"
      },
      placeholder = code$placeholder,
      engine = paste0(
        "DMAR fits the model with OpenMx, the estimation engine, starting ",
        "from the estimates of this lavaan fit of the same model and data; ",
        "see lavaan for more details. The OpenMx model has no bounds on its ",
        "variances, such as bounds = \"pos.var\" sets in lavaan, so the ",
        "choice is a different model."),
      defined = "keep")
  }

  # ---- every constraint of a test restricts the parameters --------------
  # A test refers its statistic to a chi square distribution with as many
  # degrees of freedom as constraints, which is right only when no
  # constraint is a function of the others. Through the first builds of
  # 1.1.0, hypotheses = list(h = c("a", "a")) or c("a", "2*a") gave the
  # likelihood ratio statistic of "a" alone with 2 degrees of freedom
  # (p = .083 where "a" alone gave .026), and c("a - b", "b - a") and
  # c("a", "b", "a + b") counted a constraint that restricts nothing. The
  # check reads the Jacobian of each set of constraints with respect to the
  # free parameters at the unrestricted (lavaan) estimate, a point off the
  # null set, and refuses a set whose rank is below its number of
  # constraints, before any OpenMx fit; see .mbco_redundant().
  theta_hat <- .mbco_free_values(pt)
  shown <- c(stats::setNames(hidden$expr, hidden$name),
             stats::setNames(effects$expr[!is.na(effects$expr)],
                             effects$term[!is.na(effects$expr)]))
  for (i in seq_len(nrow(effects))) {
    members <- strsplit(effects$test[i], ";", fixed = TRUE)[[1L]]
    red <- .mbco_redundant(expanded[members], theta_hat)
    if (length(red$redundant) == 0L && length(members) > 1L) {
      # Independent off the null set; is any constraint 0 wherever the
      # others are? See .mbco_vanishing().
      rank <- red$rank
      red <- .mbco_vanishing(expanded[members], theta_hat)
      red$rank <- rank
    }
    if (length(red$redundant) == 0L) next
    .mbco_refuse_redundant(effects$term[i], unname(shown[members]),
                           unname(expanded[members]), red,
                           user = effects$term[i] %in% user_terms)
  }

  # ---- the scale of the variables ---------------------------------------
  # When the observed variables differ in scale, the free parameters do too
  # (a coefficient of 2e-7 beside a variance of 1e6), and the OpenMx
  # optimizers stop short: through the first builds of 1.1.0, with x in
  # thousands and m in thousandths, the statistic for a = 0 was 1316 where
  # least squares gives 4.98, a joint test 22.0 where it gives 11.9, and
  # the indirect effect 7.56 where it gives 4.98. Each observed variable is
  # divided by a power of 2 near its standard deviation for the OpenMx
  # fits, every parameter is fit in the matching metric, and every effect
  # and constraint is written in the parameters of the original metric, so
  # estimates, standard errors, intervals, and likelihood ratio statistics
  # are those of the original variables. See .mbco_rescale_prep().
  resc <- .mbco_rescale_prep(pt, data, manifest)
  pt_fit <- pt
  data_fit <- data
  alg_fit <- alg_defs
  mult <- numeric(0)  # original parameter = mult * fitted parameter
  dev_shift <- 0      # original deviance = fitted deviance + dev_shift
  if (!is.null(resc)) {
    mult <- resc$mult
    par_rows <- pt$op != ":="
    pt_fit$est[par_rows] <- pt$est[par_rows] / resc$mono[par_rows]
    pt_fit$start[par_rows] <- pt$start[par_rows] / resc$mono[par_rows]
    for (v in manifest) data_fit[[v]] <- data[[v]] / resc$cj[[v]]
    alg_fit$expr <- vapply(alg_defs$expr, .mbco_rescale_expr, character(1L),
                           mult = mult, USE.NAMES = FALSE)
    dev_shift <- 2 * sum(vapply(manifest, function(v) {
      sum(is.finite(data[[v]])) * log(resc$cj[[v]])
    }, numeric(1L)))
  }
  # The size of each constraint in the metric of the fit, for the
  # feasibility tolerance of the null fits (1 without rescaling).
  units <- vapply(expanded, .mbco_unit, numeric(1L), mult = mult)

  # ---- build and fit the full model in OpenMx ----------------------------
  old_opt <- OpenMx::mxOption(NULL, "Default optimizer")
  on.exit(OpenMx::mxOption(NULL, "Default optimizer", old_opt),
          add = TRUE)
  OpenMx::mxOption(NULL, "Default optimizer", optimizer)

  full <- .mbco_build_mx(pt_fit, data_fit, manifest, latents, alg_fit,
                         group = group, glabels = glabels)
  # The profile likelihood interval of each effect. When the variables were
  # rescaled, the interval is profiled for the effect divided by its unit
  # (an algebra in the metric of the fit), so that the interval search
  # works on a value near the size of the parameters; its limits are then
  # multiplied by the unit.
  ci_term <- effects$term[!is.na(effects$expr)]
  ci_alg <- stats::setNames(ci_term, ci_term)
  if (ci_method == "profile_likelihood") {
    if (length(mult) > 0L) {
      ci_alg[] <- paste0("mbco_ci_", seq_along(ci_term))
      full <- do.call(OpenMx::mxModel, c(list(full), lapply(
        seq_along(ci_term), function(j) {
          OpenMx::mxAlgebraFromString(
            paste0(ci_term[j], " / ",
                   sprintf("%.17g", units[[ci_term[j]]])),
            name = ci_alg[[j]])
        })))
    }
    full <- OpenMx::mxModel(full,
      OpenMx::mxCI(unname(ci_alg), interval = conf_level))
  }
  fit_full <- suppressWarnings(OpenMx::mxRun(
    full, silent = TRUE, intervals = ci_method == "profile_likelihood"))
  if (!fit_full@output$status$code %in% c(0L, 1L)) {
    fit_full <- .mbco_try_hard(
      full, extraTries = 15L,
      intervals = ci_method == "profile_likelihood")
  }
  D_full <- fit_full@output$Minus2LogLikelihood
  # Cross-check against lavaan, which fits the same likelihood; if
  # lavaan found a better optimum, push OpenMx to at least match it.
  D_lav <- try(-2 * as.numeric(lavaan::fitMeasures(lav, "logl")),
               silent = TRUE)
  if (!inherits(D_lav, "try-error") && length(D_lav) == 1L &&
      is.finite(D_lav) && !lav_input && D_full + dev_shift > D_lav + 0.01) {
    fit_full <- .mbco_try_hard(
      full, extraTries = 25L,
      intervals = ci_method == "profile_likelihood")
    D_full <- fit_full@output$Minus2LogLikelihood
    if (D_full + dev_shift > D_lav + 0.01) {
      .dmar_warn(c(
        "The OpenMx full-model deviance (", round(D_full + dev_shift, 3),
        ") did not reach the lavaan value (", round(D_lav, 3),
        "); results may reflect a local solution."), "result_qualified")
    }
  }
  if (!fit_full@output$status$code %in% c(0L, 1L)) {
    .dmar_warn(c(
      "The full model finished with OpenMx status code ",
      fit_full@output$status$code,
      "; interpret the results with care."), "result_qualified")
  }
  n_obs <- .mbco_numobs(fit_full)
  n_par <- length(coef(fit_full))
  aic_full <- D_full + dev_shift + 2 * n_par
  bic_full <- D_full + dev_shift + log(n_obs) * n_par

  # ---- point estimates and standard errors -------------------------------
  # Values of every defined algebra at the full-model solution (also
  # the scale reference for the constraint feasibility tolerance).
  alg_est <- vapply(alg_defs$name, function(nm) {
    as.numeric(OpenMx::mxEvalByName(nm, fit_full, compute = TRUE)[1L])
  }, numeric(1L))
  k <- nrow(effects)
  estimate <- ifelse(is.na(effects$expr), NA_real_,
                     alg_est[effects$term])
  se <- vapply(seq_len(k), function(i) {
    if (is.na(effects$expr[i])) return(NA_real_)
    out <- try(suppressWarnings(
      as.numeric(OpenMx::mxSE(effects$term[i], fit_full,
                              forceName = TRUE,
                              silent = TRUE)[1L])), silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else out
  }, numeric(1L))

  # ---- confidence intervals ----------------------------------------------
  alpha <- 1 - conf_level
  z <- qnorm(1 - alpha / 2)
  ci_lower <- estimate - z * se
  ci_upper <- estimate + z * se
  mc_summary <- NULL  # the counts of a Monte Carlo interval, when one ran
  if (ci_method == "profile_likelihood") {
    ints <- fit_full@output$confidenceIntervals
    for (i in seq_len(k)) {
      if (is.na(effects$expr[i])) next
      nm <- ci_alg[[effects$term[i]]]
      row <- grep(paste0("\\.", nm, "\\["), rownames(ints))
      if (length(row) == 0L) {
        row <- match(nm, rownames(ints))
      }
      u_i <- if (nm == effects$term[i]) 1 else units[[effects$term[i]]]
      lo <- NA_real_
      hi <- NA_real_
      if (length(row) >= 1L && !is.na(row[1L])) {
        lo <- ints[row[1L], "lbound"] * u_i
        hi <- ints[row[1L], "ubound"] * u_i
      }
      if (is.na(lo) || is.na(hi)) {
        .dmar_warn(c(
          "The profile likelihood interval for ",
          sQuote(effects$term[i]),
          " did not fully converge; the Wald limit is reported ",
          "for the missing bound."), "result_qualified")
      }
      if (!is.na(lo)) ci_lower[i] <- lo
      if (!is.na(hi)) ci_upper[i] <- hi
    }
  } else if (ci_method == "monte_carlo") {
    V <- try(suppressWarnings(vcov(fit_full)), silent = TRUE)
    if (inherits(V, "try-error")) {
      .dmar_warn(c(
        "The parameter covariance matrix is unavailable; Wald ",
        "intervals are reported instead of Monte Carlo."), "result_qualified")
    } else {
      draws <- as.data.frame(
        MASS::mvrnorm(B, mu = coef(fit_full), Sigma = V))
      # Draws of the fitted parameters, in the original metric.
      for (lb in intersect(names(mult), names(draws))) {
        draws[[lb]] <- draws[[lb]] * mult[[lb]]
      }
      # A draw at which an effect is not a finite number (a division by a
      # draw near zero, a function outside its domain) has no value for
      # that effect: it is dropped and counted, effect by effect, and the
      # rule of R/dmar_bootstrap.R decides whether the effect's interval
      # is formed. Through DMAR 1.0.x such draws were removed without a
      # count or a warning.
      mc_verdicts <- list()
      for (i in seq_len(k)) {
        if (is.na(effects$expr[i])) next
        vals <- eval(parse(text = expanded[[effects$term[i]]])[[1L]],
                     envir = draws)
        vals <- vals[is.finite(vals)]
        v <- .dmar_boot_settle(
          B, length(vals), tail = alpha / 2, conf_level = conf_level,
          fun = "mediation_mbco",
          quantity = paste0("the effect '", effects$term[i], "'"),
          dropped_because = paste0(
            "gave a value of '", effects$term[i], "' that is not a finite ",
            "number (a division by a draw near zero, or a function ",
            "outside its domain)"),
          replications = "Monte Carlo draws", signal = FALSE)
        mc_verdicts[[effects$term[i]]] <- v
        if (v$formed) {
          lims <- quantile(vals, c(alpha / 2, 1 - alpha / 2), names = FALSE)
          ci_lower[i] <- lims[1L]
          ci_upper[i] <- lims[2L]
        } else {
          ci_lower[i] <- ci_upper[i] <- NA_real_
        }
      }
      mc_summary <- .mbco_mc_settle(mc_verdicts)
    }
  }

  # ---- MBCO likelihood ratio test per effect -----------------------------
  lrt <- rep(NA_real_, k)
  delta_r2 <- NULL
  r2_full <- .mbco_r2(fit_full, glabels)
  for (i in seq_len(k)) {
    tests <- strsplit(effects$test[i], ";", fixed = TRUE)[[1L]]
    nul <- .mbco_fit_null(fit_full, tests, expanded[tests],
                          alg_est[tests], D_full,
                          label = effects$term[i],
                          extra_starts = eff$starts[[effects$term[i]]],
                          mult = mult, units = units[tests])
    lrt[i] <- nul$lrt
    dr2 <- if (is.null(nul$fit)) {
      rep(NA_real_, length(r2_full))
    } else {
      # Align by name: a branch fit that fixes a variable's only
      # incoming path to zero leaves that variable exogenous in the
      # null model, so its entry is absent there; its null R2 is 0
      # (the null model explains none of its variance), and the drop
      # is the full-model R2 itself.
      r2_null <- .mbco_r2(nul$fit, glabels)
      r2_full - ifelse(is.na(r2_null[names(r2_full)]), 0,
                       r2_null[names(r2_full)])
    }
    delta_r2 <- cbind(delta_r2, dr2)
  }
  colnames(delta_r2) <- effects$term
  rownames(delta_r2) <- names(r2_full)
  df <- as.integer(effects$df)
  p_value <- pchisq(lrt, df = df, lower.tail = FALSE)

  out <- data.frame(
    pathway   = effects$pathway,
    term      = effects$term,
    estimate  = unname(estimate),
    se        = unname(se),
    ci_lower  = unname(ci_lower),
    ci_upper  = unname(ci_upper),
    lrt       = lrt,
    df        = df,
    p_value   = p_value,
    delta_aic = lrt - 2 * df,
    delta_bic = lrt - df * log(n_obs),
    stringsAsFactors = FALSE
  )
  class(out) <- c("dmar_mediation_mbco", "data.frame")
  out <- .as_dmar_tbl(out, conf_level = conf_level)
  attr(out, "ci_method") <- ci_method
  if (!is.null(mc_summary)) out <- .dmar_boot_attrs(out, mc_summary)
  attr(out, "optimizer") <- optimizer
  attr(out, "x") <- x
  attr(out, "y") <- y
  attr(out, "N") <- n_obs
  if (ngroups > 1L) attr(out, "groups") <- glabels
  if (!is.null(mod)) {
    attr(out, "moderation") <- list(
      moderator = mod$moderator,
      values = stats::setNames(mod$values, mod$labels),
      effects = eff$polys)
  }
  attr(out, "deviance") <- D_full + dev_shift
  attr(out, "aic") <- aic_full
  attr(out, "bic") <- bic_full
  attr(out, "n_par") <- n_par
  attr(out, "R2") <- r2_full
  attr(out, "delta_R2") <- delta_r2
  attr(out, "mx_model") <- fit_full
  if (!is.null(resc)) {
    attr(out, "rescale") <- list(divisors = resc$cj, multipliers = mult)
  }
  out
}

# Enumerate the mediation effects implied by the regression structure.
# Returns list(x, y, table) where table has columns term / pathway /
# expr / test / df. With multiple groups, every effect appears once per
# group plus a between-group difference against the first (reference)
# group. With a moderator (mod), each edge coefficient is a linear
# function of the moderator, each pathway effect a polynomial in it;
# effects are then probed at the moderator values and the polynomial's
# nonconstant coefficients form the moderation (constancy) test.
.mbco_effects <- function(reg, x, y, have_extra, glabels = character(0),
                          mod = NULL) {
  ngroups <- max(1L, length(glabels))
  reg1 <- reg[reg$group <= 1L, , drop = FALSE]
  from <- reg1$rhs
  to   <- reg1$lhs
  vars <- unique(c(from, to))
  empty <- data.frame(term = character(0), pathway = character(0),
                      expr = character(0), test = character(0),
                      df = integer(0), stringsAsFactors = FALSE)
  if (is.null(x) != is.null(y)) {
    .dmar_refuse_argument(
      "Supply both 'x' and 'y', or neither.", arg = c("x", "y"))
  }
  if (is.null(x)) {
    sources <- setdiff(unique(from), unique(to))
    if (!is.null(mod)) {
      sources <- setdiff(sources, c(mod$moderator, mod$prods))
    }
    sinks <- setdiff(unique(to), unique(from))
    if (length(sources) == 1L && length(sinks) == 1L) {
      x <- sources
      y <- sinks
    } else if (have_extra) {
      return(list(x = NULL, y = NULL, table = empty, starts = list(),
                  polys = list()))
    } else {
      .dmar_refuse_argument(c(
        "The focal predictor and outcome are ambiguous (sources: ",
        paste(sQuote(sources), collapse = ", "), "; outcomes: ",
        paste(sQuote(sinks), collapse = ", "),
        "); supply 'x' and 'y'."), arg = c("x", "y"))
    }
  }
  for (nm in c(x, y)) {
    if (!nm %in% vars) {
      .dmar_refuse_argument(c(
        sQuote(nm), " has no regression path in the model."), arg = NULL)
    }
  }
  # A feedback loop makes "the pathways from x to y" ill defined (total
  # effects would be infinite series, not sums of simple paths), so
  # automatic enumeration refuses rather than mislabel.
  # The depth-first search keeps its vertex colors and the cycle flag in
  # one state object that every recursive call reads and updates.
  dfs <- new.env(parent = emptyenv())
  dfs$color <- stats::setNames(rep(0L, length(vars)), vars)
  dfs$cyc <- FALSE
  visit <- function(v) {
    dfs$color[[v]] <- 1L
    for (w in to[from == v]) {
      if (dfs$cyc) return(invisible(NULL))
      if (dfs$color[[w]] == 1L) {
        dfs$cyc <- TRUE
        return(invisible(NULL))
      }
      if (dfs$color[[w]] == 0L) visit(w)
    }
    dfs$color[[v]] <- 2L
  }
  for (v in vars) if (dfs$color[[v]] == 0L) visit(v)
  if (dfs$cyc) {
    .dmar_refuse_argument(c(
      "The regression structure is nonrecursive (a feedback loop ",
      "connects the equations), so enumerating pathways from ",
      sQuote(x), " to ", sQuote(y), " is not meaningful. State the ",
      "quantities to test through ':=' definitions or 'hypotheses' ",
      "and omit 'x' and 'y'."), arg = c("hypotheses", "x", "y"))
  }
  # Depth-first enumeration of simple directed pathways x -> ... -> y;
  # each call returns the pathways found below its node, in edge order.
  walk <- function(node, seen) {
    found <- list()
    for (i in which(from == node)) {
      nxt <- to[i]
      if (nxt %in% seen) next
      if (nxt == y) {
        found[[length(found) + 1L]] <- c(seen, nxt)
      } else {
        found <- c(found, walk(nxt, c(seen, nxt)))
      }
    }
    found
  }
  paths <- walk(x, x)
  if (length(paths) == 0L) {
    .dmar_refuse_argument(c(
      "No directed pathway from ", sQuote(x), " to ", sQuote(y),
      " exists in the model."), arg = NULL)
  }
  n_edges  <- lengths(paths) - 1L
  direct   <- paths[n_edges == 1L]
  indirect <- paths[n_edges > 1L]
  if (length(indirect) == 0L && !have_extra) {
    .dmar_refuse_argument(c(
      "Only the direct path connects ", sQuote(x), " and ",
      sQuote(y), "; there is no indirect pathway to test."), arg = NULL)
  }
  san <- function(s) gsub("[^A-Za-z0-9]", "_", s)
  arrow <- function(nodes) paste(nodes, collapse = " -> ")
  label_edge <- function(f, t, g) {
    hit <- reg$rhs == f & reg$lhs == t
    if (ngroups > 1L) hit <- hit & reg$group == g
    reg$mx_label[hit][1L]
  }
  # An edge coefficient as a polynomial in the moderator: the base
  # label, plus the interaction label when the edge is moderated.
  edge_poly <- function(f, t, g) {
    base <- label_edge(f, t, g)
    inter <- if (is.null(mod)) NULL else mod$inter[[paste(t, f,
                                                          sep = "\r")]]
    c(base, inter)
  }
  path_edges <- function(nodes, g) {
    lapply(seq_len(length(nodes) - 1L), function(i) {
      edge_poly(nodes[i], nodes[i + 1L], g)
    })
  }
  poly_of <- function(edges) Reduce(.mbco_poly_mul, edges)
  poly_sum <- function(polys) {
    vapply(seq_len(max(lengths(polys))), function(j) {
      parts <- unlist(lapply(polys, function(p) {
        if (j <= length(p)) paste0("(", p[j], ")") else NULL
      }))
      paste(parts, collapse = " + ")
    }, character(1L))
  }
  inters_of <- function(edges) {
    unique(unlist(lapply(edges, function(e) {
      if (length(e) > 1L) e[2L] else NULL
    })))
  }
  # Branch-aware starting values for the constrained null fits, keyed
  # by row term: symbolic assignments evaluated at the full-model
  # estimates (see .mbco_fit_null).
  starts_env <- new.env(parent = emptyenv())
  # The moderator polynomial of every moderated base effect, for the
  # conditional-effect display (see plot_mediation_mbco).
  polys_env <- new.env(parent = emptyenv())
  polys_env$items <- list()
  # Rows for one base effect. A constant effect is one row; a
  # moderated effect gets one row per probe value, one row per
  # nonconstant polynomial coefficient (the linear one is the index of
  # moderated mediation), and, when there are several coefficients, a
  # joint constancy test. 'edges' (a single pathway's factor
  # structure) seeds one start per branch of a probe constraint;
  # 'inters' seeds the moderation tests with all interactions zeroed.
  rows_for <- function(term, pathway, poly, edges = NULL,
                       inters = character(0)) {
    if (length(poly) == 1L) {
      if (!is.null(edges)) {
        assign(term, lapply(edges, function(e) {
          stats::setNames("0", e[1L])
        }), envir = starts_env)
      }
      return(data.frame(term = term, pathway = pathway, expr = poly,
                        test = term, df = 1L, stringsAsFactors = FALSE))
    }
    w <- mod$moderator
    polys_env$items[[term]] <- list(term = term, pathway = pathway,
                                    coefs = poly)
    out <- empty
    for (j in seq_along(mod$values)) {
      v <- mod$values[j]
      at_expr <- paste0("(", poly[1L], ")")
      for (dgr in seq_len(length(poly) - 1L)) {
        at_expr <- paste0(at_expr, " + (", poly[dgr + 1L], ")*",
                          sprintf("(%.15g)", v^dgr))
      }
      tm <- paste0(term, "_at_", mod$labels[j])
      out <- rbind(out, data.frame(
        term = tm,
        pathway = paste0(pathway, " at ", w, " = ", signif(v, 4L)),
        expr = at_expr, test = tm, df = 1L, stringsAsFactors = FALSE))
      if (!is.null(edges)) {
        # One start per branch: that edge's conditional coefficient
        # (base + inter w) is zero at the probed value.
        assign(tm, lapply(edges, function(e) {
          if (length(e) > 1L) {
            stats::setNames(sprintf("-(%s)*(%.15g)", e[2L], v), e[1L])
          } else {
            stats::setNames("0", e[1L])
          }
        }), envir = starts_env)
      }
    }
    inter_zero <- if (length(inters) > 0L) {
      list(stats::setNames(rep("0", length(inters)), inters))
    } else {
      list()
    }
    if (length(poly) == 2L) {
      tm <- paste0(term, "_moderation")
      out <- rbind(out, data.frame(
        term = tm,
        pathway = paste0(pathway, ", change per unit ", w),
        expr = poly[2L], test = tm, df = 1L, stringsAsFactors = FALSE))
      assign(tm, inter_zero, envir = starts_env)
    } else {
      coef_terms <- paste0(term, "_moderation_w",
                           c("", seq_len(length(poly) - 2L) + 1L))
      for (dgr in seq_len(length(poly) - 1L)) {
        out <- rbind(out, data.frame(
          term = coef_terms[dgr],
          pathway = paste0(pathway, ", coefficient of ", w,
                           if (dgr > 1L) paste0("^", dgr) else ""),
          expr = poly[dgr + 1L], test = coef_terms[dgr], df = 1L,
          stringsAsFactors = FALSE))
        assign(coef_terms[dgr], inter_zero, envir = starts_env)
      }
      tm <- paste0(term, "_moderation")
      out <- rbind(out, data.frame(
        term = tm,
        pathway = paste0(pathway, ", constant over ", w,
                         " (joint test)"),
        expr = NA_character_,
        test = paste(coef_terms, collapse = ";"),
        df = length(poly) - 1L, stringsAsFactors = FALSE))
      assign(tm, inter_zero, envir = starts_env)
    }
    out
  }
  # Base effect list for one group.
  base_tab <- function(g) {
    tab <- empty
    spec_terms <- vapply(indirect, function(nodes) {
      mid <- nodes[-c(1L, length(nodes))]
      paste0("indirect_via_", paste(san(mid), collapse = "_"))
    }, character(1L))
    spec_edges <- lapply(indirect, path_edges, g = g)
    spec_polys <- lapply(spec_edges, poly_of)
    all_edges <- do.call(c, spec_edges)
    if (length(direct) == 1L) {
      d_edges <- path_edges(direct[[1L]], g)
      d_poly <- poly_of(d_edges)
      all_edges <- c(all_edges, d_edges)
      if (length(indirect) > 0L) {
        tab <- rbind(tab, rows_for(
          "total_effect", paste0(arrow(c(x, y)), " (all pathways)"),
          poly_sum(c(list(d_poly), spec_polys)),
          inters = inters_of(all_edges)))
      }
      tab <- rbind(tab, rows_for("direct_effect", arrow(c(x, y)),
                                 d_poly, edges = d_edges,
                                 inters = inters_of(d_edges)))
    }
    if (length(indirect) > 1L) {
      tab <- rbind(tab, rows_for(
        "total_indirect",
        paste0(arrow(c(x, y)), " (all indirect pathways)"),
        poly_sum(spec_polys),
        inters = inters_of(do.call(c, spec_edges))))
    }
    for (s in seq_along(indirect)) {
      tab <- rbind(tab, rows_for(spec_terms[s],
                                 arrow(indirect[[s]]),
                                 spec_polys[[s]],
                                 edges = spec_edges[[s]],
                                 inters = inters_of(spec_edges[[s]])))
    }
    tab
  }
  if (ngroups == 1L) {
    tab1 <- base_tab(1L)
    return(list(x = x, y = y, table = tab1,
                starts = as.list(starts_env),
                polys = polys_env$items))
  }
  # Multiple groups: per-group rows, then differences against the
  # reference (first) group. A difference whose two expressions are
  # identical (the path labels are shared across groups, an equality
  # constraint) is degenerate and skipped.
  gsuf <- make.unique(vapply(glabels, san, character(1L)), sep = "_")
  per_group <- lapply(seq_len(ngroups), base_tab)
  tab <- empty
  for (r in seq_len(nrow(per_group[[1L]]))) {
    base <- per_group[[1L]]$term[r]
    for (g in seq_len(ngroups)) {
      row_g <- per_group[[g]][r, ]
      tab <- rbind(tab, data.frame(
        term = paste0(base, "_", gsuf[g]),
        pathway = paste0(row_g$pathway, " [", glabels[g], "]"),
        expr = row_g$expr, test = paste0(base, "_", gsuf[g]), df = 1L,
        stringsAsFactors = FALSE))
    }
    for (g in seq_len(ngroups)[-1L]) {
      if (per_group[[g]]$expr[r] == per_group[[1L]]$expr[r]) next
      tm <- paste0(base, "_", gsuf[g], "_minus_", gsuf[1L])
      tab <- rbind(tab, data.frame(
        term = tm,
        pathway = paste0(per_group[[1L]]$pathway[r], " [",
                         glabels[g], " - ", glabels[1L], "]"),
        expr = paste0("(", per_group[[g]]$expr[r], ") - (",
                      per_group[[1L]]$expr[r], ")"),
        test = tm, df = 1L, stringsAsFactors = FALSE))
    }
  }
  list(x = x, y = y, table = tab, starts = list(), polys = list())
}

# Multiply a coefficient polynomial in the moderator by one edge's
# polynomial (both held as expression strings, constant term first).
.mbco_poly_mul <- function(p, e) {
  deg <- length(p) + length(e) - 1L
  out <- character(deg)
  for (j in seq_len(deg)) {
    parts <- character(0)
    for (a in seq_along(e)) {
      i <- j - a + 1L
      if (i >= 1L && i <= length(p)) {
        parts <- c(parts, paste0("(", p[i], ")*(", e[a], ")"))
      }
    }
    out[j] <- paste(parts, collapse = " + ")
  }
  out
}

# Rewrite 'a:b' interaction terms in lavaan model syntax into product
# columns named a_by_b, so the model can state moderation directly. Also
# returns the columns added to the data, as "a * b" named by the column
# (`products`), which the lavaan code of a refusal describes.
.mbco_product_terms <- function(model, data) {
  pat <- "([A-Za-z][A-Za-z0-9._]*)[ \t]*:[ \t]*([A-Za-z][A-Za-z0-9._]*)"
  hits <- regmatches(model, gregexpr(pat, model))[[1L]]
  products <- character(0)
  for (h in unique(hits)) {
    parts <- strsplit(gsub("[ \t]", "", h), ":", fixed = TRUE)[[1L]]
    if (!all(parts %in% names(data))) next
    nm <- paste0(parts[1L], "_by_", parts[2L])
    prod <- data[[parts[1L]]] * data[[parts[2L]]]
    if (nm %in% names(data)) {
      if (!isTRUE(all.equal(data[[nm]], prod, tolerance = 1e-10))) {
        .dmar_refuse_argument(c(
          "Cannot create the interaction column ", sQuote(nm),
          " for the term ", sQuote(h), "; a different column with ",
          "that name already exists."), arg = NULL)
      }
    } else {
      data[[nm]] <- prod
      products[nm] <- paste(parts[1L], "*", parts[2L])
    }
    model <- gsub(h, nm, model, fixed = TRUE)
  }
  list(model = model, data = data, products = products)
}

# .mbco_no_max_code(model, data, S, M, N, group, miss, from_S, products)
#
# The lavaan code of mediation_mbco()'s own lavaan fit, for the refusal of
# a fit whose likelihood has no maximum (.no_max_refuse()): the arguments
# of the call as DMAR passed them (`args`), the variables to select from a
# covariance matrix and its means (`select`), and the words that say what
# the placeholders stand for (`placeholder`, NULL for lavaan's own "your
# data in place of my_data" or "your covariance matrix in place of my_S
# and its means in place of my_M"). DMAR fits raw data, and from summary
# statistics the data .mbco_surrogate() builds with exactly their moments.
# Complete data maximum likelihood depends on the data only through the
# means and the covariance matrix, and lavaan rescales a covariance matrix
# it is given to the divisor N it uses for data, so the code for the
# summary statistics of one group fits them directly, the same fit. Means
# are paired with the rows of the covariance matrix by position, so means
# named in another order are selected by name. Without M the surrogate
# has zero means, and so has the code. A covariance matrix given as a data
# frame is one group's, as mediation_mbco() reads it, and the placeholder
# says to pass it as a matrix, since lavaan reads a list (a data frame is
# one) given as sample.cov as one matrix per group. Summary statistics of
# several groups, or a model whose product terms DMAR rewrote into columns
# of the surrogate, are printed with the data DMAR built, which the
# placeholder describes as .mbco_surrogate() builds it, seed included: a
# product column's moments depend on the rows drawn, which the seed of 1
# fixes, so the fit is reproduced only from those rows. `from_S` says the
# data came from summary statistics, and `products` names the columns
# added (.mbco_product_terms()).
.mbco_no_max_code <- function(model, data, S, M, N, group, miss, from_S,
                              products) {
  multi <- is.list(S) && !is.data.frame(S)
  added <- if (length(products)) {
    paste0(", once ", .no_max_and(paste(names(products), "=", products)),
           if (length(products) == 1L) " is" else " are", " added to it")
  }
  if (from_S && !multi && !length(products)) {
    frame <- is.data.frame(S)
    S <- as.matrix(S)
    vars <- if (!is.null(colnames(S))) colnames(S) else rownames(S)
    mu <- if (is.null(M)) stats::setNames(rep(0, length(vars)), vars) else M
    return(list(
      args = list(model = model, sample.cov = S, sample.mean = mu,
                  sample.nobs = N, meanstructure = TRUE, fixed.x = FALSE,
                  se = "none", warn = FALSE),
      select = if (!is.null(M) && !identical(names(M), vars)) {
        list(sample.cov = vars, sample.mean = vars)
      },
      placeholder = if (frame || is.null(M)) {
        paste0(if (frame) "as.matrix() of ", "your covariance matrix in ",
               "place of my_S and ",
               if (is.null(M)) "a zero mean for each of its rows" else
                 "its means", " in place of my_M")
      }))
  }
  draw <- paste0("set.seed(1) and then MASS::mvrnorm(N, mu = ",
                 "M[colnames(S)], Sigma = S, empirical = TRUE)",
                 if (multi) " on that group's N, M, and S",
                 ", with mu = 0 without M")
  list(
    args = list(model = model, data = data, meanstructure = TRUE,
                fixed.x = FALSE, missing = miss, se = "none", group = group,
                warn = FALSE),
    select = list(),
    placeholder = if (from_S) {
      paste0("the data mediation_mbco() builds from the summary ",
             "statistics in place of my_data (",
             if (multi) {
               paste0("for each group, ", draw, ", the groups stacked in ",
                      "the order of S with the group's name in the column ",
                      group)
             } else {
               draw
             }, ")", added)
    } else if (length(products)) {
      paste0("your data in place of my_data", added)
    })
}

# Identify which regression edges the moderator moderates: an edge
# A -> child is moderated when the child's equation also contains a
# predictor that is elementwise A times the moderator. Resolves the
# probe values (mean and one standard deviation either side by
# default; the two observed values for a two-valued moderator).
.mbco_moderation <- function(reg, data, manifest, moderator,
                             probe_values) {
  reg1 <- reg[reg$group <= 1L, , drop = FALSE]
  if (moderator %in% reg1$lhs) {
    .dmar_refuse_argument(c(
      "The moderator ", sQuote(moderator), " is endogenous (it has ",
      "an incoming regression path); probing requires an exogenous ",
      "moderator."), arg = NULL)
  }
  w <- data[[moderator]]
  inter <- list()
  prods <- character(0)
  no_main <- character(0)
  for (i in seq_len(nrow(reg1))) {
    child <- reg1$lhs[i]
    P <- reg1$rhs[i]
    if (!P %in% names(data)) next
    vv <- data[[P]]
    for (A in setdiff(manifest, c(P, moderator))) {
      a <- data[[A]]
      ok <- is.finite(a) & is.finite(w) & is.finite(vv)
      if (sum(ok) < 3L) next
      if (max(abs(vv[ok] - a[ok] * w[ok])) <
          1e-8 * max(1, max(abs(vv[ok])))) {
        key <- paste(child, A, sep = "\r")
        if (is.null(inter[[key]])) inter[[key]] <- reg1$mx_label[i]
        prods <- c(prods, P)
        if (!any(reg1$lhs == child & reg1$rhs == A)) {
          no_main <- c(no_main, paste0(sQuote(A), " in the equation ",
                                       "for ", sQuote(child)))
        }
        break
      }
    }
  }
  if (length(inter) == 0L) {
    .dmar_warn(c(
      "No interaction between ", sQuote(moderator), " and a ",
      "modeled variable was found in the model; the moderator ",
      "is ignored. State the interaction as, e.g., 'm ~ x + w ",
      "+ x:w' or include the product column as a predictor."),
      "argument_ignored")
    return(NULL)
  }
  if (length(no_main) > 0L) {
    .dmar_warn(c(
      "Interaction without the matching main effect of ",
      paste(unique(no_main), collapse = ", "),
      " (the principle of marginality); conditional effects ",
      "along that path are not probed."), "result_qualified")
  }
  if (is.null(probe_values)) {
    uw <- sort(unique(w[is.finite(w)]))
    if (length(uw) == 2L) {
      values <- uw
      labels <- vapply(uw, .mbco_value_label, character(1L))
    } else {
      mw <- mean(w, na.rm = TRUE)
      sw <- stats::sd(w, na.rm = TRUE)
      values <- c(mw - sw, mw, mw + sw)
      labels <- c("low", "mean", "high")
    }
  } else {
    if (!is.numeric(probe_values) || length(probe_values) < 1L ||
        anyNA(probe_values)) {
      .dmar_refuse_argument(c(
        "'probe_values' must be a numeric vector without missing ",
        "values."), arg = "probe_values")
    }
    values <- as.numeric(probe_values)
    labels <- names(probe_values)
    if (is.null(labels)) labels <- rep("", length(values))
    labels[!nzchar(labels)] <- vapply(values[!nzchar(labels)],
                                      .mbco_value_label, character(1L))
    labels <- make.unique(gsub("[^A-Za-z0-9]", "_", labels),
                          sep = "_")
  }
  list(moderator = moderator, inter = inter, prods = unique(prods),
       values = values, labels = labels)
}

# A short, name-safe rendering of a numeric probe value ("m" for the
# minus sign): -1.5 becomes "m1_5".
.mbco_value_label <- function(v) {
  s <- formatC(v, format = "g", digits = 6L)
  s <- sub("^-", "m", s)
  gsub("[^A-Za-z0-9]", "_", s)
}

# Substitute effect and ":=" names in an expression by their definitions
# until only parameter labels remain.
.mbco_expand <- function(e, defs) {
  if (is.name(e)) {
    nm <- as.character(e)
    if (nm %in% names(defs)) {
      return(call("(", .mbco_expand(defs[[nm]], defs)))
    }
    return(e)
  }
  if (is.call(e)) {
    for (i in seq_along(e)[-1L]) e[[i]] <- .mbco_expand(e[[i]], defs)
  }
  e
}

# The free parameters at the estimates of the lavaan fit, named by their
# OpenMx labels (a label shared by several rows is one parameter, whose
# rows carry one value). These are the unrestricted maximum likelihood
# estimates, from the fit whose likelihood the OpenMx full model repeats.
.mbco_free_values <- function(pt) {
  rows <- which(pt$op != ":=" & pt$free > 0)
  rows <- rows[!duplicated(pt$mx_label[rows])]
  stats::setNames(as.numeric(pt$est[rows]), pt$mx_label[rows])
}

# The Jacobian of constraints (character expressions in the parameter
# labels, as .mbco_expand() leaves them) with respect to the parameters
# they involve, at `theta`. The derivative of each expression is symbolic
# (stats::deriv()), so a product, ratio, or polynomial of coefficients is
# differentiated exactly; an expression with a function outside R's
# derivative table falls back to a central difference. Columns are the
# parameters that appear in any of the expressions. A point outside the
# domain of an expression (an evaluation that warns, such as log() of a
# negative number, or fails) leaves NA entries, and the point is not used.
.mbco_jacobian <- function(exprs, theta) {
  calls <- lapply(exprs, function(e) parse(text = e)[[1L]])
  vars <- intersect(unique(unlist(lapply(calls, all.vars))), names(theta))
  J <- matrix(0, length(calls), length(vars), dimnames = list(NULL, vars))
  if (length(vars) == 0L) return(J)
  env <- as.list(theta[vars])
  for (i in seq_along(calls)) {
    g <- tryCatch({
      val <- eval(stats::deriv(calls[[i]], vars), env)
      as.numeric(attr(val, "gradient"))
    }, warning = function(w) NULL, error = function(e) NULL)
    if (length(g) != length(vars)) {
      g <- vapply(vars, function(v) {
        h <- .Machine$double.eps^(1 / 3) * max(abs(env[[v]]), 1)
        up <- env
        dn <- env
        up[[v]] <- up[[v]] + h
        dn[[v]] <- dn[[v]] - h
        out <- tryCatch(
          (eval(calls[[i]], up) - eval(calls[[i]], dn)) / (2 * h),
          warning = function(w) NA_real_, error = function(e) NA_real_)
        if (length(out) == 1L) as.numeric(out) else NA_real_
      }, numeric(1L))
    }
    J[i, ] <- g
  }
  J
}

# The numerical rank of a Jacobian, read so that it does not depend on the
# units of the parameters or on the scale of the constraints: every
# column is divided by its largest absolute entry (a parameter in other
# units multiplies its column by a constant) and every row by its length
# (a constraint multiplied by a constant multiplies its row), and the rank
# is the number of singular values above 1e-7 times the largest, the
# relative tolerance of base R's qr(). Exact dependence (a repeated
# constraint, a multiple, a sum of others) leaves singular values at
# rounding level, far below that. NA when an entry is not finite.
.mbco_rank <- function(J, tol = 1e-7) {
  if (any(!is.finite(J))) return(NA_integer_)
  J <- J[, colSums(abs(J)) > 0, drop = FALSE]
  if (ncol(J) == 0L || nrow(J) == 0L) return(0L)
  J <- sweep(J, 2L, apply(abs(J), 2L, max), "/")
  len <- sqrt(rowSums(J^2))
  J[len > 0, ] <- J[len > 0, , drop = FALSE] / len[len > 0]
  d <- svd(J, nu = 0L, nv = 0L)$d
  if (max(d) == 0) return(0L)
  as.integer(sum(d > tol * max(d)))
}

# Which of a test's constraints add no restriction to the ones before
# them. Constraints are functionally independent when their Jacobian with
# respect to the free parameters has full row rank at a generic point;
# the rank is read at the unrestricted estimate, which lies off the null
# set, and is never read on the null set, where a product constraint has
# a vanishing gradient at the crossing of its branches whatever the
# hypothesis. To keep an estimate that falls exactly on a special value
# (a coefficient of exactly 0 from rounded summary statistics) from
# making independent constraints look dependent, the rank is also read at
# a fixed point near the estimate, and the larger of the two ranks is
# used: the Jacobian of analytic constraints (sums, products, ratios, and
# the functions of R's derivative table) has its largest rank at all but
# special points. Constraints are then taken in the order given, and one
# that does not raise the rank of those before it is redundant. Returns
# list(redundant, kept, rank, n), with no redundant constraints when the
# Jacobian cannot be evaluated at either point (the test then runs as
# before).
.mbco_redundant <- function(exprs, theta) {
  n <- length(exprs)
  none <- list(redundant = integer(0), kept = seq_len(n), rank = n, n = n)
  if (n == 0L) return(none)
  u <- (seq_along(theta) * (sqrt(5) - 1) / 2) %% 1
  near <- theta + (0.1 + 0.2 * u) * pmax(abs(theta), 1)
  best <- NULL
  for (point in list(theta, near)) {
    J <- .mbco_jacobian(exprs, point)
    r <- .mbco_rank(J)
    if (is.na(r)) next
    if (is.null(best) || r > best$rank) best <- list(J = J, rank = r)
    if (r == n) return(none)
  }
  if (is.null(best)) return(none)
  kept <- integer(0)
  redundant <- integer(0)
  for (i in seq_len(n)) {
    if (.mbco_rank(best$J[c(kept, i), , drop = FALSE]) > length(kept)) {
      kept <- c(kept, i)
    } else {
      redundant <- c(redundant, i)
    }
  }
  list(redundant = redundant, kept = kept, rank = best$rank, n = n,
       found_by = "jacobian")
}

# Which of a test's constraints are 0 wherever the others are. A set can
# be functionally independent (its Jacobian has full rank off the null
# set, so .mbco_redundant() passes it) and still have a null set with
# fewer restrictions than constraints: c("a*b", "b") has the null set
# b = 0, c("a", "a*cp") the null set a = 0, and c("a*b", "a*b*cp") the
# null set ab = 0, each one restriction where 2 degrees of freedom would
# be counted. A constraint adds no restriction to the others when it is 0
# at every point of their null set (it then vanishes on the joint null
# set, which is the null set of the others alone).
#
# The null set of the others is the union of branches, one for each way of
# choosing one factor of each of their constraints (a product is 0 exactly
# where one of its factors is; a constant multiple, a sign, a positive
# power, and a denominator do not change where a constraint is 0), the fact
# the branch fits of .mbco_fit_null() use. A factor that is a parameter
# label is set to 0 exactly; a factor that is any other expression (a sum
# such as a3*b1 + a1*b3, a difference, exp(a) - 1) is reached by
# Gauss-Newton projection (minimum norm steps in scaled coordinates,
# halved until the largest relative factor value falls). Points of a
# branch are taken from fixed starting points that cover both signs of
# every parameter and scaled magnitudes from 1/4 to 2 on a log scale, from
# a Kronecker sequence (fractional parts of multiples of square roots of
# primes; no random number generator is used, so the caller's stream is
# untouched), so they are generic points of the branch and not special
# points where more than its own factors vanish. A projected point counts
# when its largest factor value is below 1e-10 of that factor's typical
# size (its median absolute value over the starting points). The candidate
# is read at each point relative to its own typical size, and it vanishes
# on the branch when it is below 1e-6 at every point (a constraint that is
# not 0 on a branch is of the order of its typical size at a generic point
# of it). It adds no restriction when it vanishes on every branch; a
# branch with general factors that no projection reaches leaves the answer
# open, and the constraint is then kept, as is every constraint of a set
# with more than 256 branches.
#
# Two scalings are used, and a constraint is redundant only when it
# vanishes under both: the magnitude of each parameter's unrestricted
# estimate, and the reciprocal of the largest derivative of the
# constraints with respect to the parameter (at the estimate and at the
# fixed point near it of .mbco_redundant()), which balances the parameters
# as the constraints see them. Each scaling moves with the units of the
# parameters, so the result does not depend on them; the second keeps two
# independent constraints from looking alike where the parameters differ
# in size by orders of magnitude (c("a + b", "a - b") with a near 1e-7 and
# b near 250 are both close to b = 0 at points scaled by the estimates).
# Constraints are taken in the order given, and one found redundant is set
# aside before the next is read. A point outside the domain of an
# expression (an evaluation that warns or fails) is not used. The check
# cannot see a set none of whose constraints vanishes on the null set of
# the others but whose null set has fewer restrictions than constraints,
# such as c("a*b", "a*cp", "b*cp") (null set: two of a, b, cp zero, two
# restrictions); see the help page.
.mbco_vanishing <- function(exprs, theta, n_starts = 32L) {
  n <- length(exprs)
  none <- list(redundant = integer(0), kept = seq_len(n), n = n,
               found_by = "null_set")
  if (n < 2L) return(none)
  calls <- lapply(exprs, function(e) parse(text = e)[[1L]])
  vars <- intersect(unique(unlist(lapply(calls, all.vars))), names(theta))
  p <- length(vars)
  if (p == 0L) return(none)
  # The factors of each constraint, as deparsed expressions.
  factors_of <- function(e) {
    if (is.call(e)) {
      op <- as.character(e[[1L]])
      if (op == "(") return(factors_of(e[[2L]]))
      if (op %in% c("-", "+") && length(e) == 2L) return(factors_of(e[[2L]]))
      if (op == "*") return(unique(c(factors_of(e[[2L]]),
                                     factors_of(e[[3L]]))))
      if (op == "/") return(factors_of(e[[2L]]))
      if (op == "^" && is.numeric(e[[3L]]) && length(e[[3L]]) == 1L &&
          e[[3L]] > 0) {
        return(factors_of(e[[2L]]))
      }
    }
    if (is.numeric(e)) return(character(0))
    deparse1(e)
  }
  fac <- lapply(calls, factors_of)
  fac_txt <- unique(unlist(fac))
  fac_calls <- lapply(fac_txt, function(f) parse(text = f)[[1L]])
  # The two scalings.
  s_est <- abs(as.numeric(theta[vars]))
  s_est[!is.finite(s_est) | s_est == 0] <- 1
  u <- (seq_along(theta) * (sqrt(5) - 1) / 2) %% 1
  near <- theta + (0.1 + 0.2 * u) * pmax(abs(theta), 1)
  colmax <- rep(0, p)
  for (point in list(theta, near)) {
    J0 <- abs(.mbco_jacobian(exprs, point)[, vars, drop = FALSE])
    J0[!is.finite(J0)] <- 0
    colmax <- pmax(colmax, apply(J0, 2L, max))
  }
  s_jac <- ifelse(colmax > 0, 1 / colmax, s_est)
  # Values (and gradients) of a list of expressions at a point. Symbolic
  # gradients are compiled once; a function outside R's derivative table
  # falls back to a central difference.
  quiet_eval <- function(e, env) {
    tryCatch(eval(e, env), warning = function(w) NULL,
             error = function(e) NULL)
  }
  compile <- function(cls) {
    lapply(cls, function(cl) {
      tryCatch(stats::deriv(cl, vars), error = function(e) NULL)
    })
  }
  d_cons <- compile(calls)
  d_fac <- compile(fac_calls)
  ev <- function(th, cls, dcls, s, grad = TRUE) {
    env <- stats::setNames(as.list(th), vars)
    g <- numeric(length(cls))
    J <- matrix(0, length(cls), p)
    for (r in seq_along(cls)) {
      val <- if (grad && !is.null(dcls[[r]])) quiet_eval(dcls[[r]], env)
      if (!is.null(val) && length(val) == 1L &&
          length(attr(val, "gradient")) == p) {
        g[r] <- as.numeric(val)
        J[r, ] <- as.numeric(attr(val, "gradient"))
        next
      }
      val <- quiet_eval(cls[[r]], env)
      if (length(val) != 1L) return(NULL)
      g[r] <- as.numeric(val)
      if (!grad) next
      for (j in seq_len(p)) {
        h <- .Machine$double.eps^(1 / 3) * max(abs(th[j]), s[j])
        up <- env
        dn <- env
        up[[j]] <- up[[j]] + h
        dn[[j]] <- dn[[j]] - h
        fu <- quiet_eval(cls[[r]], up)
        fd <- quiet_eval(cls[[r]], dn)
        if (length(fu) != 1L || length(fd) != 1L) return(NULL)
        J[r, j] <- (as.numeric(fu) - as.numeric(fd)) / (2 * h)
      }
    }
    if (any(!is.finite(g)) || any(!is.finite(J))) return(NULL)
    list(g = g, J = J)
  }
  # The fixed starting points, as scaled coordinates.
  primes <- function(k) {
    out <- integer(0)
    cand <- 2L
    while (length(out) < k) {
      if (all(cand %% out[out <= sqrt(cand)] != 0L)) out <- c(out, cand)
      cand <- cand + 1L
    }
    out
  }
  pr <- sqrt(primes(2L * p))
  al <- pr[seq_len(p)] %% 1
  be <- pr[p + seq_len(p)] %% 1
  unit_starts <- t(vapply(seq_len(n_starts), function(i) {
    mag <- 2^(-2 + 3 * ((i * al) %% 1))
    sgn <- ifelse((i * be) %% 1 < 0.5, -1, 1)
    sgn * mag
  }, numeric(p)))
  typical <- function(starts, cls, dcls, s) {
    sizes <- vapply(seq_len(n_starts), function(i) {
      z <- ev(starts[i, ], cls, dcls, s, grad = FALSE)
      if (is.null(z)) rep(NA_real_, length(cls)) else abs(z$g)
    }, numeric(length(cls)))
    G <- apply(matrix(sizes, nrow = length(cls)), 1L, stats::median,
               na.rm = TRUE)
    G[!is.finite(G) | G == 0] <- 1
    G
  }
  # Gauss-Newton projection onto the common zero set of expressions.
  project <- function(th, cls, dcls, G, s) {
    cur <- ev(th, cls, dcls, s)
    if (is.null(cur)) return(NULL)
    rel <- function(z) max(abs(z$g) / G)
    r <- rel(cur)
    slow <- 0L
    for (it in seq_len(200L)) {
      if (r <= 1e-30) break
      sv <- svd(cur$J * rep(s, each = nrow(cur$J)))
      if (max(sv$d) == 0) break
      keep <- sv$d > 1e-12 * max(sv$d)
      step <- -s * as.numeric(
        sv$v[, keep, drop = FALSE] %*%
          (crossprod(sv$u[, keep, drop = FALSE], cur$g) / sv$d[keep]))
      lam <- 1
      new <- NULL
      while (lam >= 2^-30) {
        cand <- ev(th + lam * step, cls, dcls, s)
        if (!is.null(cand) && rel(cand) < r) {
          new <- cand
          break
        }
        lam <- lam / 2
      }
      if (is.null(new)) break
      th <- th + lam * step
      r_new <- rel(new)
      slow <- if (r_new > 0.5 * r) slow + 1L else 0L
      cur <- new
      r <- r_new
      if (slow >= 5L) break
    }
    list(theta = th, r = r)
  }
  # Does constraint k vanish on the null set of `others`, read with the
  # scaling s? TRUE, FALSE, or NA (open).
  vanishes_at <- function(k, others, s) {
    combos <- prod(lengths(fac[others]))
    if (combos == 0L || combos > 256L) return(NA)
    grid <- expand.grid(fac[others], stringsAsFactors = FALSE)
    branches <- unique(lapply(seq_len(nrow(grid)), function(i) {
      sort(unique(as.character(unlist(grid[i, ], use.names = FALSE))))
    }))
    starts <- unit_starts * rep(s, each = n_starts)
    G_k <- typical(starts, calls[k], d_cons[k], s)
    for (br in branches) {
      zero <- intersect(br, vars)
      general <- setdiff(br, vars)
      idx <- match(general, fac_txt)
      G_b <- if (length(idx)) typical(starts, fac_calls[idx], d_fac[idx], s)
      counted <- 0L
      for (i in seq_len(if (length(general)) n_starts else 4L)) {
        th <- starts[i, ]
        th[match(zero, vars)] <- 0
        if (length(general)) {
          # The labels of the branch are kept at 0 by including them as
          # factors of the projection.
          cls <- c(fac_calls[idx], lapply(zero, as.name))
          dcls <- c(d_fac[idx], compile(lapply(zero, as.name)))
          pr_i <- project(th, cls, dcls, c(G_b, s[match(zero, vars)]), s)
          if (is.null(pr_i) || pr_i$r > 1e-10) next
          th <- pr_i$theta
        }
        z <- ev(th, calls[k], d_cons[k], s, grad = FALSE)
        if (is.null(z)) next
        counted <- counted + 1L
        if (abs(z$g) / G_k > 1e-6) return(FALSE)
      }
      if (counted == 0L) return(NA)
    }
    TRUE
  }
  active <- seq_len(n)
  redundant <- integer(0)
  for (k in seq_len(n)) {
    others <- setdiff(active, k)
    if (length(others) == 0L) next
    if (isTRUE(vanishes_at(k, others, s_est)) &&
        isTRUE(vanishes_at(k, others, s_jac))) {
      redundant <- c(redundant, k)
      active <- others
    }
  }
  list(redundant = redundant, kept = active, n = n, found_by = "null_set")
}

# The refusal of a test whose constraints are not functionally
# independent (see .mbco_redundant()), or one of whose constraints is 0
# wherever the others are (see .mbco_vanishing()). `shown` holds the
# constraints as the user wrote them (or as the model's moderation
# produced them), `expanded` the same in parameter labels, `user` whether
# the test comes from 'hypotheses'.
.mbco_refuse_redundant <- function(term, shown, expanded, red, user) {
  q <- function(x) paste0("'", x, "'", collapse = ", ")
  shown <- trimws(shown)
  n <- red$n
  red_txt <- shown[red$redundant]
  if (identical(red$found_by, "null_set")) {
    kept_txt <- shown[red$kept]
    msg <- c(
      if (user) "The joint hypothesis '" else "The joint test '", term,
      "' constrains ", n, " expressions to 0 (", q(shown), "), but ",
      q(red_txt), if (length(red_txt) == 1L) " adds" else " add",
      " no restriction to ", q(kept_txt), ": ",
      if (length(red_txt) == 1L) "it is" else "each is",
      " 0 at every point where ", q(kept_txt),
      if (length(kept_txt) == 1L) " is 0" else " are all 0",
      ", so the null hypothesis is the same without ",
      if (length(red_txt) == 1L) "it" else "them",
      ", and it restricts the parameters in fewer than ", n,
      " directions. The likelihood ratio statistic would be referred to a ",
      "chi square distribution with ", n, " degrees of freedom, which ",
      "counts every expression as a restriction, and its p-value, ",
      "delta_aic, and delta_bic would be wrong. ")
  } else if (n == 1L) {
    msg <- c(
      "The quantity '", term, "' (", shown, ") does not depend on the ",
      "free parameters of the model: its gradient with respect to them is ",
      "0 at the unrestricted estimate. Setting it to 0 restricts nothing, ",
      "so a likelihood ratio test of it with 1 degree of freedom would be ",
      "wrong. ")
  } else {
    is_repeat <- vapply(red$redundant, function(i) {
      any(expanded[red$kept] == expanded[i])
    }, logical(1L))
    rep_txt <- shown[red$redundant[is_repeat]]
    dep_txt <- shown[red$redundant[!is_repeat]]
    msg <- c(
      if (user) "The joint hypothesis '" else "The joint test '", term,
      "' constrains ", n, " expressions to 0 (", q(shown), "), but they ",
      "are not functionally independent: they restrict the parameters in ",
      red$rank, " independent direction", if (red$rank != 1L) "s",
      ", not ", n, ". ",
      if (length(rep_txt)) c(
        q(rep_txt), if (length(rep_txt) == 1L) " repeats an expression" else
          " repeat expressions", " already given",
        if (length(dep_txt)) "; " else ". "),
      if (length(dep_txt)) c(
        q(dep_txt), if (length(dep_txt) == 1L) " adds" else " add",
        " no restriction to ",
        if (length(red$kept)) q(shown[red$kept]) else "the others",
        ": at the unrestricted estimate, the gradient with respect to the ",
        "free parameters of each such expression is a linear combination ",
        "of the gradients of the others. "),
      "The likelihood ratio statistic would be referred to a chi square ",
      "distribution with ", n, " degrees of freedom, which counts every ",
      "expression as a restriction, and its p-value, delta_aic, and ",
      "delta_bic would be wrong. ")
  }
  msg <- c(msg, if (user && n == 1L) {
    c("Remove it from 'hypotheses'. ")
  } else if (user) {
    c("Remove ", q(red_txt), " from 'hypotheses' and test the remaining ",
      "expressions. ")
  } else if (n == 1L) {
    c("Check its ':=' definition and the labels of the model. ")
  } else {
    c("The labels of the model make these constraints dependent; give ",
      "the paths that enter them labels of their own. ")
  }, "Nothing was fitted.")
  .dmar_refuse_argument(msg, arg = if (user) "hypotheses" else NULL,
                        term = term, constraints = shown,
                        redundant = red_txt, kept = shown[red$kept],
                        n_constraints = n, rank = red$rank,
                        found_by = red$found_by)
}

# The factors of a sum of products: the parameter labels common to every
# term, and the sum that remains when one of each is taken out of every
# term. A conditional indirect effect at a probe value is written as
# (a0*b) + (a1*b)*(v); its factors are b and a0 + a1*v, so its null set is
# b = 0 or a0 + a1*v = 0. Each term must be a product of labels and
# numbers (with signs and parentheses); otherwise, or when no label is
# common to every term, the sum is one factor, returned as its text.
.mbco_common_factors <- function(e) {
  # The terms of the sum, each with its sign.
  collect <- function(x, sign) {
    if (is.call(x) && identical(x[[1L]], as.name("("))) {
      return(collect(x[[2L]], sign))
    }
    if (is.call(x) && as.character(x[[1L]]) %in% c("+", "-")) {
      flip <- if (identical(x[[1L]], as.name("-"))) -sign else sign
      if (length(x) == 2L) return(collect(x[[2L]], flip))
      return(c(collect(x[[2L]], sign), collect(x[[3L]], flip)))
    }
    list(list(sign = sign, expr = x))
  }
  terms <- collect(e, 1)
  mono <- function(x) {
    if (is.name(x)) return(list(labels = as.character(x), const = 1))
    if (is.numeric(x) && length(x) == 1L) {
      return(list(labels = character(0), const = as.numeric(x)))
    }
    if (is.call(x) && identical(x[[1L]], as.name("("))) return(mono(x[[2L]]))
    if (is.call(x) && identical(x[[1L]], as.name("-")) && length(x) == 2L) {
      m <- mono(x[[2L]])
      if (is.null(m)) return(NULL)
      m$const <- -m$const
      return(m)
    }
    if (is.call(x) && identical(x[[1L]], as.name("*"))) {
      a <- mono(x[[2L]])
      b <- mono(x[[3L]])
      if (is.null(a) || is.null(b)) return(NULL)
      return(list(labels = c(a$labels, b$labels), const = a$const * b$const))
    }
    NULL
  }
  monos <- lapply(terms, function(t) {
    m <- mono(t$expr)
    if (!is.null(m)) m$const <- m$const * t$sign
    m
  })
  whole <- deparse1(e)
  if (length(monos) < 2L || any(vapply(monos, is.null, logical(1L)))) {
    return(whole)
  }
  common <- Reduce(intersect, lapply(monos, `[[`, "labels"))
  if (length(common) == 0L) return(whole)
  rest <- vapply(monos, function(m) {
    labs <- m$labels
    for (cl in common) labs <- labs[-match(cl, labs)]
    paste(c(sprintf("(%.17g)", m$const), labs), collapse = "*")
  }, character(1L))
  c(common, paste(rest, collapse = " + "))
}

# A starting point on the null set for a constrained fit: the parameters
# `est` (named, in the metric of the fit) moved onto the common zero set of
# the constraints `exprs` (in the parameter labels of the original metric)
# by Gauss-Newton steps in the metric of the fit (minimum norm steps,
# halved until the largest weighted constraint value falls; at most 50),
# with the constraints weighted by the matrix W (the diagonal of 1/units,
# or the T of .mbco_fit_null() for a joint test). Parameters in `fixed` do
# not move. The minimum norm step moves the estimates the least distance
# that meets the linearized constraints, so the start is near the
# estimates. Returns `est` with the moved values, or NULL when the
# constraints cannot be differentiated symbolically, an evaluation warns
# or fails, or the weighted constraint values do not fall below 1e-10
# times max(1, the largest moved parameter).
.mbco_project_start <- function(exprs, est, mult, W, fixed = character(0)) {
  calls <- lapply(exprs, function(e) parse(text = e)[[1L]])
  vars <- setdiff(intersect(unique(unlist(lapply(calls, all.vars))),
                            names(est)), fixed)
  if (length(vars) == 0L) return(NULL)
  allv <- intersect(unique(unlist(lapply(calls, all.vars))), names(est))
  d <- lapply(calls, function(cl) {
    tryCatch(stats::deriv(cl, vars), error = function(e) NULL)
  })
  if (any(vapply(d, is.null, logical(1L)))) return(NULL)
  m_of <- function(lb) {
    out <- rep(1, length(lb))
    hit <- lb %in% names(mult)
    out[hit] <- mult[lb[hit]]
    out
  }
  mv <- m_of(vars)
  base <- as.list(est[allv] * m_of(allv))
  h_of <- function(th) {
    env <- base
    env[vars] <- as.list(th * mv)
    vals <- lapply(d, function(dd) {
      tryCatch(eval(dd, env), warning = function(w) NULL,
               error = function(e) NULL)
    })
    if (any(vapply(vals, is.null, logical(1L)))) return(NULL)
    g <- vapply(vals, function(v) as.numeric(v)[1L], numeric(1L))
    J <- do.call(rbind, lapply(vals, function(v) {
      as.numeric(attr(v, "gradient"))
    })) * rep(mv, each = length(vals))
    h <- as.numeric(W %*% g)
    if (any(!is.finite(h)) || any(!is.finite(J))) return(NULL)
    list(h = h, J = W %*% J)
  }
  th <- as.numeric(est[vars])
  cur <- h_of(th)
  if (is.null(cur)) return(NULL)
  r <- max(abs(cur$h))
  for (it in seq_len(50L)) {
    if (r <= 1e-12 * max(1, abs(th))) break
    sv <- svd(cur$J)
    if (max(sv$d) == 0) return(NULL)
    keep <- sv$d > 1e-10 * max(sv$d)
    step <- -as.numeric(sv$v[, keep, drop = FALSE] %*%
      (crossprod(sv$u[, keep, drop = FALSE], cur$h) / sv$d[keep]))
    lam <- 1
    new <- NULL
    while (lam >= 2^-20) {
      cand <- h_of(th + lam * step)
      if (!is.null(cand) && max(abs(cand$h)) < r) {
        new <- cand
        break
      }
      lam <- lam / 2
    }
    if (is.null(new)) break
    th <- th + lam * step
    cur <- new
    r <- max(abs(cur$h))
  }
  if (!is.finite(r) || r > 1e-10 * max(1, abs(th))) return(NULL)
  est[vars] <- th
  est
}

# The rescaling of badly scaled observed variables for the OpenMx fits, as
# lavaan's rescale_data does (lav_rescale_prep() and lav_rescale_mono() of
# lavaan 0.7). The variables are badly scaled when the ratio of the largest
# to the smallest standard deviation (over all rows, all groups) is at
# least 1000, or the largest is at least 1000, or the smallest at most
# .001, lavaan's rule. Each observed variable v is then divided by
# c_v = 2^round(log2(sd_v)), a power of 2, so the rescaling and its
# inverse are exact in floating point. Rescaling whenever some c_v is not 1
# was tried and set aside: on the memory study moments of Tofighi and
# Kelley (2020; standard deviations .5 to 3.4) the null fits in the
# rescaled metric ended with status code 6 more often, and the feasibility
# tolerance of .mbco_fit_null(), read in the units of the fit, admitted a
# total effect statistic of 39.692 where the exact value (the constraint
# is cov(x, y) = 0 in that saturated model) is 39.703. Scales that differ
# by a factor of 10 to 300 are left to the search of .mbco_fit_null()
# (branch fits and a start on the null set).
# A latent variable takes the scale of
# the indicator whose loading is fixed (the first such indicator; 1 when
# none is, as for a latent variable whose variance is fixed). A parameter
# of the original metric is the parameter of the rescaled metric times
# its multiplier: c(lhs) / c(rhs) for a regression lhs ~ rhs,
# c(indicator) / c(latent) for a loading, c(lhs) c(rhs) for a variance or
# covariance, and c(v) for an intercept. Fixed values are divided by their
# multipliers too, so the rescaled model is the same model. Returns NULL
# (no rescaling) when the variables are not badly scaled or every c_v is
# 1, when a standard deviation is not a positive number, or when a label
# shared by several
# rows (an equality constraint) would need different multipliers, which
# no common parameter of the rescaled metric can carry; otherwise
# list(cj, mono, mult): the divisors, the multiplier of every row of the
# parameter table, and the multiplier of every free label.
.mbco_rescale_prep <- function(pt, data, manifest) {
  sds <- vapply(manifest, function(v) stats::sd(data[[v]], na.rm = TRUE),
                numeric(1L))
  if (any(!is.finite(sds)) || any(sds <= 0)) return(NULL)
  if (max(sds) / min(sds) < 1000 && max(sds) < 1000 && min(sds) > 0.001) {
    return(NULL)
  }
  cj <- 2^round(log2(sds))
  if (all(cj == 1)) return(NULL)
  par_rows <- which(pt$op != ":=")
  grp <- if (is.null(pt$group)) rep(1L, nrow(pt)) else pt$group
  mono <- rep(1, nrow(pt))
  for (g in unique(grp[par_rows])) {
    rows <- par_rows[grp[par_rows] == g]
    lv <- unique(pt$lhs[rows][pt$op[rows] == "=~"])
    s_lv <- stats::setNames(rep(NA_real_, length(lv)), lv)
    scale_of <- function(v) {
      if (v %in% names(cj)) return(cj[[v]])
      if (v %in% lv && !is.na(s_lv[[v]])) return(s_lv[[v]])
      1
    }
    for (pass in 1:5) {
      for (l in lv) {
        if (!is.na(s_lv[[l]])) next
        mk <- rows[pt$op[rows] == "=~" & pt$lhs[rows] == l &
                     pt$free[rows] == 0 & is.finite(pt$est[rows]) &
                     pt$est[rows] != 0]
        if (length(mk) == 0L) {
          s_lv[[l]] <- 1
          next
        }
        ind <- pt$rhs[mk[1L]]
        if (ind %in% lv && is.na(s_lv[[ind]])) next
        s_lv[[l]] <- scale_of(ind)
      }
      if (!anyNA(s_lv)) break
    }
    s_lv[is.na(s_lv)] <- 1
    for (i in rows) {
      mono[i] <- switch(pt$op[i],
        "=~" = scale_of(pt$rhs[i]) / scale_of(pt$lhs[i]),
        "~"  = scale_of(pt$lhs[i]) / scale_of(pt$rhs[i]),
        "~~" = scale_of(pt$lhs[i]) * scale_of(pt$rhs[i]),
        "~1" = scale_of(pt$lhs[i]),
        1)
    }
  }
  labs <- pt$mx_label[par_rows]
  by_label <- split(mono[par_rows], labs)
  if (any(vapply(by_label, function(v) any(v != v[1L]), logical(1L)))) {
    return(NULL)
  }
  mult <- vapply(by_label, `[`, numeric(1L), 1L)
  free_labels <- unique(pt$mx_label[par_rows][pt$free[par_rows] > 0])
  mult <- mult[names(mult) %in% free_labels & mult != 1]
  list(cj = cj, mono = mono, mult = mult)
}

# An algebra of the effects or constraints, rewritten for the rescaled fit:
# every free label L with a multiplier m becomes (m*L), the parameter of
# the original metric, so the algebra keeps its value in the original
# metric. The multipliers are powers of 2, written with 17 significant
# digits, so they are read back exactly.
.mbco_rescale_expr <- function(e, mult) {
  if (length(mult) == 0L) return(e)
  sub <- function(x) {
    if (is.name(x)) {
      nm <- as.character(x)
      if (nm %in% names(mult)) {
        return(call("(", call("*", unname(mult[[nm]]), x)))
      }
      return(x)
    }
    if (is.call(x)) {
      for (i in seq_along(x)[-1L]) x[[i]] <- sub(x[[i]])
    }
    x
  }
  paste(deparse(sub(parse(text = e)[[1L]]),
                control = c("keepNA", "keepInteger", "niceNames",
                            "showAttributes", "digits17"),
                width.cutoff = 500L), collapse = " ")
}

# The size of a constraint's unit in the metric of the fit: the factor
# that turns a value of the constraint in the rescaled metric into its
# value in the original metric. A label contributes its multiplier, a
# product or ratio the product or ratio of its factors' units, a power
# the power of its base's unit, and a numeric constant nothing; a sum of
# terms with one unit has that unit, and a sum of terms with different
# units (a constraint that compares parameters of different units) the
# geometric mean of theirs. 1 for every constraint when nothing was
# rescaled.
.mbco_unit <- function(e, mult) {
  if (length(mult) == 0L) return(1)
  u <- function(x) {
    if (is.name(x)) {
      nm <- as.character(x)
      return(if (nm %in% names(mult)) unname(mult[[nm]]) else 1)
    }
    if (is.numeric(x)) return(NA_real_)
    if (!is.call(x)) return(1)
    op <- as.character(x[[1L]])
    if (op == "(") return(u(x[[2L]]))
    if (op %in% c("+", "-") && length(x) == 2L) return(u(x[[2L]]))
    if (op %in% c("+", "-")) {
      v <- c(u(x[[2L]]), u(x[[3L]]))
      v <- v[!is.na(v)]
      if (length(v) == 0L) return(NA_real_)
      if (length(v) == 1L || v[1L] == v[2L]) return(v[1L])
      return(sqrt(v[1L] * v[2L]))
    }
    if (op %in% c("*", "/")) {
      a <- u(x[[2L]])
      b <- u(x[[3L]])
      if (is.na(a) && is.na(b)) return(NA_real_)
      if (is.na(a)) a <- 1
      if (is.na(b)) b <- 1
      return(if (op == "*") a * b else a / b)
    }
    if (op == "^" && is.numeric(x[[3L]])) return(u(x[[2L]])^x[[3L]])
    args <- vapply(as.list(x)[-1L], u, numeric(1L))
    if (length(args) && all(is.na(args))) NA_real_ else 1
  }
  out <- u(parse(text = e)[[1L]])
  if (is.na(out) || !is.finite(out) || out <= 0) 1 else out
}

# One mxPath per parameter-table row (":=" rows are handled as
# algebras, not paths).
.mbco_mx_paths <- function(pt) {
  paths <- list()
  for (i in seq_len(nrow(pt))) {
    op <- pt$op[i]
    if (op == ":=") next
    from <- switch(op, "=~" = pt$lhs[i], "~" = pt$rhs[i],
                   "~~" = pt$lhs[i], "~1" = "one")
    to   <- switch(op, "=~" = pt$rhs[i], "~" = pt$lhs[i],
                   "~~" = pt$rhs[i], "~1" = pt$lhs[i])
    val <- pt$est[i]
    if (!is.finite(val)) val <- pt$start[i]
    if (!is.finite(val)) val <- 0
    paths[[length(paths) + 1L]] <- OpenMx::mxPath(
      from = from, to = to, arrows = if (op == "~~") 2L else 1L,
      free = pt$free[i] > 0, values = val, labels = pt$mx_label[i])
  }
  paths
}

# Translate a lavaan parameter table into an OpenMx RAM model (one
# submodel per group when there are several) with one mxAlgebra per
# defined quantity at the top level.
.mbco_build_mx <- function(pt, data, manifest, latents, alg_defs,
                           group = NULL, glabels = character(0)) {
  algebras <- lapply(seq_len(nrow(alg_defs)), function(i) {
    OpenMx::mxAlgebraFromString(alg_defs$expr[i],
                                name = alg_defs$name[i])
  })
  if (length(glabels) == 0L) {
    return(do.call(OpenMx::mxModel, c(
      list("mbco_full", type = "RAM", manifestVars = manifest,
           latentVars = latents),
      .mbco_mx_paths(pt),
      algebras,
      list(OpenMx::mxData(observed = data[, manifest, drop = FALSE],
                          type = "raw")))))
  }
  subs <- lapply(seq_along(glabels), function(g) {
    ptg <- pt[pt$group == g, , drop = FALSE]
    dg <- data[data[[group]] == glabels[g], manifest, drop = FALSE]
    do.call(OpenMx::mxModel, c(
      list(paste0("mbcog", g), type = "RAM", manifestVars = manifest,
           latentVars = latents),
      .mbco_mx_paths(ptg),
      list(OpenMx::mxData(observed = dg, type = "raw"))))
  })
  do.call(OpenMx::mxModel, c(
    list("mbco_full"), subs, algebras,
    list(OpenMx::mxFitFunctionMultigroup(
      paste0("mbcog", seq_along(glabels))))))
}

# Fit the null model for one test (one or several algebras jointly
# constrained to zero) from several starting configurations, verify
# feasibility, and keep the best solution.
.mbco_fit_null <- function(fit_full, terms, expanded, est_full, D_full,
                           label = terms[1L], extra_starts = NULL,
                           mult = numeric(0), units = rep(1, length(terms))) {
  # A constraint in units other than those of the fit (see
  # .mbco_rescale_prep()) is handed to the optimizer divided by its unit,
  # so that the optimizer's tolerance on it is a tolerance in the metric of
  # the fit; it is the same constraint.
  # Several constraints on parameters of very different metrics can be
  # nearly parallel in the metric of the fit (a + b and a - b with a and b
  # differing by orders of magnitude), and an optimizer then meets them
  # only approximately. When the variables were rescaled, a joint test's
  # constraints g are handed to the optimizer, and read for feasibility,
  # as h = T g with T = D^-1 U' from the singular value decomposition
  # U D V' of their Jacobian with respect to the fitted parameters at the
  # full-model estimates: the rows of the Jacobian of h are orthonormal
  # there, T is a constant invertible matrix, so h = 0 is the same null
  # set, and h is a distance in the metric of the fit (exactly so for
  # linear constraints). T is not used when that Jacobian cannot be formed
  # or, with every column divided by its largest entry (so that the units
  # of the parameters, which T is there to undo, do not count), is not of
  # full rank to 1e-8, as near the crossing of the branches of a product.
  Tm <- NULL
  if (length(terms) > 1L && length(mult) > 0L) {
    m_lab <- function(lb) {
      out <- rep(1, length(lb))
      hit <- lb %in% names(mult)
      out[hit] <- mult[lb[hit]]
      out
    }
    est0 <- coef(fit_full)
    J <- .mbco_jacobian(expanded, est0 * m_lab(names(est0)))
    if (all(is.finite(J)) && ncol(J) >= nrow(J)) {
      J <- J * rep(m_lab(colnames(J)), each = nrow(J))
      cm <- apply(abs(J), 2L, max)
      Jn <- J[, cm > 0, drop = FALSE] /
        rep(cm[cm > 0], each = nrow(J))
      dn <- if (ncol(Jn) >= nrow(Jn)) svd(Jn, nu = 0L, nv = 0L)$d else 0
      if (min(dn) > 1e-8 * max(dn)) {
        sv <- svd(J)
        Tm <- diag(1 / sv$d, nrow = length(sv$d)) %*% t(sv$u)
      }
    }
  }
  cons <- lapply(seq_along(terms), function(j) {
    lhs <- if (!is.null(Tm)) {
      Reduce(function(acc, k) {
        call("+", acc, call("*", Tm[j, k], as.name(terms[k])))
      }, seq_along(terms)[-1L],
      call("*", Tm[j, 1L], as.name(terms[1L])))
    } else if (units[j] == 1) {
      as.name(terms[j])
    } else {
      call("/", as.name(terms[j]), units[j])
    }
    eval(bquote(OpenMx::mxConstraint(.(lhs) == 0,
                                     name = .(paste0("mbco_h0_", j)))))
  })
  null_model <- do.call(OpenMx::mxModel, c(
    list(fit_full,
         name = paste0("mbco_null_", gsub("[^A-Za-z0-9_]", "_",
                                          label))),
    cons))
  est <- coef(fit_full)
  used <- unique(unlist(lapply(expanded, function(e) {
    all.vars(parse(text = e)[[1L]])
  })))
  used <- intersect(used, names(est))
  starts <- c(list(est), lapply(used, function(lb) {
    s <- est
    s[lb] <- 0
    s
  }))
  # Branch-aware starts: symbolic assignments (expressions in the
  # parameter labels) evaluated at the full-model estimates, one start
  # per branch of the constraint's null set. The expressions are in the
  # original metric (see .mbco_rescale_prep()); a value is divided by the
  # multiplier of the parameter it is assigned to.
  m_of <- function(lb) {
    out <- rep(1, length(lb))
    hit <- lb %in% names(mult)
    out[hit] <- mult[lb[hit]]
    out
  }
  est_orig <- est * m_of(names(est))
  for (a in extra_starts) {
    s <- est
    keep <- intersect(names(a), names(est))
    if (length(keep) == 0L) next
    vals <- vapply(a[keep], function(e) {
      out <- try(eval(parse(text = e)[[1L]], envir = as.list(est_orig)),
                 silent = TRUE)
      if (inherits(out, "try-error") || !is.finite(out)) NA_real_
      else as.numeric(out)
    }, numeric(1L))
    if (anyNA(vals)) next
    s[keep] <- vals / m_of(keep)
    starts[[length(starts) + 1L]] <- s
  }
  # The factors of each constraint, for the branch fits below: a product
  # is split into its factors, and a sum into the labels common to all its
  # terms and the rest (.mbco_common_factors()).
  fac_of <- function(e) {
    if (is.call(e)) {
      op <- as.character(e[[1L]])
      if (op == "(") return(fac_of(e[[2L]]))
      if (op %in% c("-", "+") && length(e) == 2L) return(fac_of(e[[2L]]))
      if (op %in% c("-", "+")) return(.mbco_common_factors(e))
      if (op == "*") return(c(fac_of(e[[2L]]), fac_of(e[[3L]])))
      if (op == "/") return(fac_of(e[[2L]]))
      if (op == "^" && is.numeric(e[[3L]]) && length(e[[3L]]) == 1L &&
          e[[3L]] > 0) {
        return(fac_of(e[[2L]]))
      }
    }
    if (is.numeric(e)) return(character(0))
    deparse1(e)
  }
  term_factors <- lapply(expanded, function(e) {
    unique(fac_of(parse(text = e)[[1L]]))
  })
  # No branches when every constraint is one expression that is not a
  # label (a total effect): its one branch is the constrained fit itself.
  has_branches <- all(lengths(term_factors) >= 1L) &&
    (any(unlist(term_factors) %in% names(est)) ||
       any(lengths(term_factors) > 1L))
  # A start on the null set: the full-model estimates moved onto it by
  # Gauss-Newton steps (see .mbco_project_start()). A sum or difference of
  # effects (a total effect, a difference between groups) has no branch
  # fits below, and through the first builds of 1.1.0 its constrained runs,
  # all started off the null set, could each stop off it, leaving the test
  # NA (10 of 100 total effects at a = b = cp = 0, N = 150). The start is
  # not needed for a constraint that has branch fits: a branch that fixes
  # labels only is an exact fit, and a branch with another factor starts on
  # that factor's null set.
  if (!has_branches) {
    W <- if (is.null(Tm)) diag(1 / units, nrow = length(terms)) else Tm
    proj <- .mbco_project_start(expanded, est, mult, W)
    if (!is.null(proj)) starts[[length(starts) + 1L]] <- proj
  }
  # Feasibility tolerance for each constraint at the solution: strict
  # enough to reject convergence failures that drift off the null set
  # (whose deviance can undercut the true constrained minimum), loose
  # enough to accept the small residual a converged sequential
  # quadratic programming step leaves. When the variables were rescaled
  # (see .mbco_rescale_prep()), a constraint must also be met in the metric
  # of the fit (its value divided by `units`, or T g for a joint test), so
  # that a constraint in small units is not counted as met at the
  # unrestricted estimate; the original rule is kept as well, so the
  # tolerance is never looser than without rescaling. Each candidate's
  # `g` is the largest ratio of a constraint value to its tolerance, and a
  # candidate is feasible when it is at most 1.
  tol_orig <- 1e-4 * pmax(1, abs(est_full))
  tol_fit <- if (is.null(Tm)) {
    1e-4 * pmax(1, abs(est_full) / units)
  } else {
    1e-4 * pmax(1, abs(as.numeric(Tm %*% est_full)))
  }
  tol <- 1
  cand_of <- function(f) {
    if (inherits(f, "try-error") || is.null(f)) return(NULL)
    g <- rep(NA_real_, length(terms))
    for (j in seq_along(terms)) {
      gj <- try(as.numeric(
        OpenMx::mxEvalByName(terms[j], f, compute = TRUE)[1L]),
        silent = TRUE)
      if (inherits(gj, "try-error") || !is.finite(gj)) return(NULL)
      g[j] <- gj
    }
    g_fit <- if (is.null(Tm)) abs(g) / units else abs(as.numeric(Tm %*% g))
    g <- max(abs(g) / tol_orig, g_fit / tol_fit)
    list(fit = f, dev = f@output$Minus2LogLikelihood, g = g,
         code = f@output$status$code)
  }
  # Each start is run and then polished by a warm restart from its own
  # solution; the restart routinely tightens feasibility and can step
  # from a stalled point down to the constrained optimum.
  run_one <- function(vals, hard = FALSE) {
    m <- OpenMx::omxSetParameters(null_model, labels = names(vals),
                                  values = vals)
    f <- if (hard) {
      try(.mbco_try_hard(m, extraTries = 12L), silent = TRUE)
    } else {
      try(suppressWarnings(OpenMx::mxRun(m, silent = TRUE)),
          silent = TRUE)
    }
    out <- list(cand_of(f))
    if (!inherits(f, "try-error")) {
      f2 <- try(suppressWarnings(OpenMx::mxRun(f, silent = TRUE)),
                silent = TRUE)
      out <- c(out, list(cand_of(f2)))
    }
    Filter(Negate(is.null), out)
  }
  valid <- function(z, codes) {
    all(z$g <= tol) && z$code %in% codes && z$dev >= D_full - 1e-4
  }
  cands <- list()
  # Deterministic branch fits. A product is 0 exactly where one of its
  # factors is, so the null set of constraints that are products is the
  # union of branches, one for each way of choosing one factor of each
  # constraint (a constant multiple, a sign, a positive power, and a
  # denominator do not change where a constraint is 0). A factor that is a
  # parameter label is fixed at 0; a factor that is any other expression
  # (the conditional coefficient a0 + a1*w of a moderated path, a sum) is
  # kept as a constraint of the branch fit, from a start on it. When every
  # factor is a label, each branch is an ordinary unconstrained fit: no
  # sequential quadratic programming step that a different optimizer build
  # might carry into the wrong basin. The reported likelihood ratio is
  # then defined by the best-fitting branch on every platform and OpenMx
  # version, not by whichever basin the constrained search happens to
  # reach. Through the first builds of 1.1.0 only products of labels had
  # branch fits; a conditional indirect effect (a0 + a1*w)*b at a probe
  # value was searched by the constrained runs alone, which could miss the
  # branch b = 0 (first-stage moderation, seed 3812, w = 1: 37.294 where
  # the smaller branch statistic, b = 0, is 33.631). Constraints with no
  # product (a total effect) skip this path; they keep the multi-start
  # constrained machinery.
  branches_done <- FALSE
  if (has_branches && prod(lengths(term_factors)) <= 64L) {
    branches_done <- TRUE
    combos <- expand.grid(term_factors, stringsAsFactors = FALSE)
    branch_sets <- unique(lapply(seq_len(nrow(combos)), function(i) {
      sort(unique(as.character(unlist(combos[i, ], use.names = FALSE))))
    }))
    for (br in branch_sets) {
      labs <- intersect(br, names(est))
      general <- setdiff(br, labs)
      m <- fit_full
      if (length(labs)) {
        m <- OpenMx::omxSetParameters(m, labels = labs, values = 0,
                                      free = FALSE)
      }
      if (length(general)) {
        u_g <- vapply(general, .mbco_unit, numeric(1L), mult = mult)
        parts <- unlist(lapply(seq_along(general), function(k) {
          nm <- paste0("mbco_bf_", k)
          lhs <- if (u_g[k] == 1) {
            as.name(nm)
          } else {
            call("/", as.name(nm), u_g[k])
          }
          list(OpenMx::mxAlgebraFromString(
                 .mbco_rescale_expr(general[k], mult), name = nm),
               eval(bquote(OpenMx::mxConstraint(
                 .(lhs) == 0, name = .(paste0("mbco_bfc_", k))))))
        }), recursive = FALSE)
        m <- do.call(OpenMx::mxModel, c(list(m), parts))
        s <- est
        s[labs] <- 0
        s_proj <- .mbco_project_start(general, s, mult,
                                      diag(1 / u_g, nrow = length(u_g)),
                                      fixed = labs)
        if (!is.null(s_proj)) s <- s_proj
        free_now <- setdiff(names(s), labs)
        m <- OpenMx::omxSetParameters(m, labels = free_now,
                                      values = s[free_now])
      }
      m <- OpenMx::mxModel(m, name = paste0(
        "mbco_branch_", gsub("[^A-Za-z0-9_]", "_", label)))
      f <- try(suppressWarnings(OpenMx::mxRun(m, silent = TRUE)),
               silent = TRUE)
      z <- cand_of(f)
      if (!is.null(z)) cands[[length(cands) + 1L]] <- z
      if (is.null(z) || !valid(z, c(0L, 1L))) branches_done <- FALSE
    }
  }
  # The null set is the union of the branches, so when every branch fit is
  # a valid solution the smallest of them is the constrained minimum, and
  # the constrained runs from the starting configurations could add only
  # the same minima or a point just off the null set that the feasibility
  # tolerance admits with a lower deviance. They run when some branch fit
  # is not valid, and for constraints without branches.
  if (!branches_done) cands <- c(cands, do.call(c, lapply(starts, run_one)))
  ok <- Filter(function(z) valid(z, c(0L, 1L)), cands)
  if (length(ok) == 0L) {
    cands <- c(cands,
               do.call(c, lapply(starts, run_one, hard = TRUE)))
    ok <- Filter(function(z) valid(z, c(0L, 1L)), cands)
  }
  if (length(ok) == 0L) {
    # Accept a feasible solution with an imperfect status code, with a
    # warning, rather than fail silently.
    ok <- Filter(function(z) valid(z, 0L:10L), cands)
    if (length(ok) > 0L) {
      .dmar_warn(c(
        "The null model for ", sQuote(label),
        " converged with OpenMx status code ",
        ok[[which.min(vapply(ok, `[[`, numeric(1L), "dev"))]]$code,
        "; interpret its test with care."), "result_qualified")
    }
  }
  if (length(ok) == 0L) {
    .dmar_warn(c(
      "The null model for ", sQuote(label),
      " could not be fit subject to the constraint; its test is ",
      "reported as NA."), "result_qualified")
    return(list(lrt = NA_real_, fit = NULL))
  }
  best <- ok[[which.min(vapply(ok, `[[`, numeric(1L), "dev"))]]
  list(lrt = max(0, best$dev - D_full), fit = best$fit)
}

# R-squared for every endogenous variable (observed or latent) of a
# fitted RAM model: 1 - residual variance / model implied variance.
# For a multiple-group fit, one entry per variable and group.
.mbco_r2 <- function(fit, glabels = character(0)) {
  if (length(fit@submodels) > 0L) {
    return(unlist(lapply(seq_along(fit@submodels), function(g) {
      r2 <- .mbco_r2_one(fit@submodels[[g]])
      names(r2) <- paste0(names(r2), " (", glabels[g], ")")
      r2
    })))
  }
  .mbco_r2_one(fit)
}

.mbco_r2_one <- function(fit) {
  A <- fit$A$values
  Ssym <- fit$S$values
  eye <- diag(nrow(A))
  inv <- solve(eye - A)
  C <- inv %*% Ssym %*% t(inv)
  endo <- rownames(A)[rowSums(A != 0 | fit$A$free) > 0]
  r2 <- 1 - diag(Ssym)[endo] / diag(C)[endo]
  names(r2) <- endo
  r2
}

# Total observations across a fit that may hold its data in submodels.
.mbco_numobs <- function(fit) {
  if (length(fit@submodels) > 0L) {
    return(sum(vapply(fit@submodels, function(m) m@data@numObs,
                      numeric(1L))))
  }
  fit@data@numObs
}

# Build a data set whose sample moments equal the supplied summary
# statistics exactly; complete-data maximum likelihood then reproduces
# the raw-data analysis. Any seed gives the same fit (the moments are
# matched exactly), so a fixed local seed is used and the caller's RNG
# state is restored.
.mbco_surrogate <- function(S, M, N) {
  S <- as.matrix(S)
  if (is.null(dimnames(S)[[1L]]) && is.null(dimnames(S)[[2L]])) {
    .dmar_refuse_argument(
      "'S' must have dimnames naming the variables.", arg = "S")
  }
  if (is.null(dimnames(S)[[1L]])) rownames(S) <- colnames(S)
  if (is.null(dimnames(S)[[2L]])) colnames(S) <- rownames(S)
  if (!isSymmetric(unname(S), tol = 1e-8)) {
    .dmar_refuse_argument(
      "'S' must be a symmetric covariance matrix.", arg = "S")
  }
  if (is.null(N) || !is.numeric(N) || length(N) != 1L || is.na(N) ||
      N <= ncol(S) || N != round(N)) {
    .dmar_refuse_argument(c(
      "'N' must be a single whole number larger than the number ",
      "of variables when 'S' is supplied."), arg = c("N", "S"))
  }
  mu <- rep(0, ncol(S))
  names(mu) <- colnames(S)
  if (!is.null(M)) {
    if (is.null(names(M)) || !all(colnames(S) %in% names(M))) {
      .dmar_refuse_argument(
        "'M' must be a named vector covering the variables in 'S'.",
        arg = c("M", "S"))
    }
    mu <- M[colnames(S)]
  }
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    .dmar_refuse_argument(c(
      "'S' is not positive definite; check the covariance ",
      "(or correlation) matrix."), arg = "S")
  }
  .dmar_local_seed(1L)
  as.data.frame(MASS::mvrnorm(N, mu = mu, Sigma = S, empirical = TRUE))
}

# mxTryHard writes progress to the console even when asked to be
# quiet; capture it so iterative refits do not flood the user.
.mbco_try_hard <- function(model, extraTries, intervals = FALSE) {
  fit <- NULL
  utils::capture.output(
    fit <- suppressWarnings(suppressMessages(OpenMx::mxTryHard(
      model, extraTries = extraTries, silent = TRUE, verbose = FALSE,
      intervals = intervals))),
    file = NULL, type = "output")
  fit
}

# The lines that introduce an MBCO table (the test, the effects and N,
# the groups, the full model's fit and R2, how to read each p_value) are
# the "header" block of .dmar_tbl_footer() (R/dmar_output_helpers.R),
# which print.dmar_tbl() prints above the table, knit_print() sets as
# paragraphs before it, and as_kable() puts first in the caption, so the
# console and a knitted report say the same thing. A column subset of
# the table loses the fit attributes, and only what exists is written.
#' @exportS3Method .dmar_tbl_footer dmar_mediation_mbco
#' @noRd
.dmar_tbl_footer.dmar_mediation_mbco <- function(x, ...) {
  header <- "Mediation tests via model-based constrained optimization (MBCO)"
  xv <- attr(x, "x", exact = TRUE)
  yv <- attr(x, "y", exact = TRUE)
  if (!is.null(attr(x, "N", exact = TRUE))) {
    header <- c(header, paste0(
      if (!is.null(xv) && !is.null(yv)) {
        paste0("  Effects of ", xv, " on ", yv)
      } else {
        "  User-defined effects"
      },
      "; N = ", attr(x, "N", exact = TRUE), "; optimizer: ",
      attr(x, "optimizer", exact = TRUE)))
  }
  groups <- attr(x, "groups", exact = TRUE)
  if (!is.null(groups)) {
    header <- c(header, paste0("  Groups (reference first): ",
                               paste(groups, collapse = ", ")))
  }
  if (!is.null(attr(x, "deviance", exact = TRUE))) {
    header <- c(header, paste0(
      "  Full model: deviance = ",
      format(attr(x, "deviance", exact = TRUE), nsmall = 3, digits = 10),
      ", AIC = ", format(attr(x, "aic", exact = TRUE), nsmall = 3,
                         digits = 10),
      ", BIC = ", format(attr(x, "bic", exact = TRUE), nsmall = 3,
                         digits = 10)))
  }
  r2 <- attr(x, "R2", exact = TRUE)
  if (length(r2) > 0L) {
    header <- c(header, paste0(
      "  R2 (full model): ",
      paste(names(r2), "=", formatC(r2, digits = 3, format = "f"),
            collapse = ", ")))
  }
  # One paragraph, broken for the console; a knitted report reads the
  # breaks as spaces.
  header <- c(header, paste0(
    "  Each p_value is from the MBCO likelihood ratio test of the\n",
    "  null model constraining that effect to zero. Causal readings\n",
    "  rest on the no omitted confounder assumption."))
  lines <- c(.dmar_footer_add(character(0), header, "header"), NextMethod())
  .dmar_footer_add(lines, .dmar_boot_footer(x, "Monte Carlo interval"),
                   "bootstrap")
}

# The condition and the summary of the Monte Carlo intervals of a call,
# from the per-effect verdicts of .dmar_boot_settle(). A call gives at
# most one warning, whatever the number of effects (through the first
# 1.1.0 builds each affected effect gave its own). Its class is the least
# favorable that applies: dmar_bootstrap_insufficient when the interval
# of any effect is withheld, dmar_bootstrap_replications_dropped when
# every interval is formed and some effect had dropped draws, and
# dmar_bootstrap_exploratory when nothing was dropped and the run is
# below the required number of draws under
# options(dmar.bootstrap_exploratory = TRUE). With one affected effect
# the message is that effect's own; with several it opens with one line
# per effect (draws dropped, draws used, whether the interval is
# reported) and then gives in full the account of the effect with the
# fewest usable draws. The condition carries the counts of that effect,
# the status of the call, and `effects`, a data frame with one row per
# affected effect (term, B_used, n_dropped, formed, status). Returns the
# verdict of the effect with the fewest usable draws, with the status of
# the call (the least favorable among the effects), which is what the
# attributes record; NULL without effects.
.mbco_mc_settle <- function(verdicts) {
  if (!length(verdicts)) return(NULL)
  used <- vapply(verdicts, function(v) v$B_used, integer(1))
  least <- verdicts[[which.min(used)]]
  order <- c("insufficient_share", "insufficient_count", "exploratory",
             "dropped", "complete")
  status <- order[min(match(vapply(verdicts, function(v) v$status,
                                   character(1)), order))]
  least$status <- status

  formed  <- vapply(verdicts, function(v) isTRUE(v$formed), logical(1))
  dropped <- vapply(verdicts, function(v) v$n_dropped > 0L, logical(1))
  explore <- formed & vapply(verdicts, function(v)
    identical(v$status, "exploratory"), logical(1))
  affected <- !formed | dropped
  if (!any(affected) && !any(explore)) return(least)

  bare <- function(text) sub("^mediation_mbco\\(\\): ", "", text)
  # A stitched account starts a sentence, so its first letter is a capital
  # ("No interval is reported.").
  sentence <- function(text) {
    paste0(toupper(substring(text, 1L, 1L)), substring(text, 2L))
  }
  own  <- function(v) {
    .dmar_boot_message(v, if (isTRUE(v$formed)) "dropped" else v$status)
  }
  text <- character(0)
  if (sum(affected) == 1L) {
    text <- own(verdicts[[which(affected)]])
  } else if (sum(affected) > 1L) {
    whole <- function(x) format(x, scientific = FALSE, trim = TRUE)
    lines <- vapply(which(affected), function(i) {
      v <- verdicts[[i]]
      paste0("'", names(verdicts)[i], "': ", whole(v$n_dropped), " of ",
             whole(v$B), " dropped, ", whole(v$B_used), " used, ",
             if (isTRUE(v$formed)) "interval reported" else
               "no interval reported")
    }, character(1))
    text <- paste0(
      "mediation_mbco(): ", sum(affected), " of the ", length(verdicts),
      " effects with a Monte Carlo interval had draws without a finite ",
      "value, or too few usable draws. By effect: ",
      paste(lines, collapse = "; "), ". The effects are counted ",
      "separately, and an effect whose interval is withheld does not ",
      "withhold the others. The account of the effect with the fewest ",
      "usable draws follows. ", sentence(bare(own(least))))
  }
  if (any(explore)) {
    marked <- verdicts[[which(explore)[1L]]]
    marked$context$quantity <- "an effect"
    note <- .dmar_boot_message(marked, "exploratory")
    text <- if (length(text)) paste(text, sentence(bare(note))) else note
  }
  code <- if (any(!formed)) "bootstrap_insufficient" else
    if (any(dropped)) "bootstrap_replications_dropped" else
      "bootstrap_exploratory"
  .dmar_warn(text, code, B = least$B, B_used = least$B_used,
             n_dropped = least$n_dropped, B_required = least$B_required,
             share_dropped = least$share_dropped, caution = least$caution,
             status = status,
             effects = data.frame(
               term = names(verdicts)[affected],
               B_used = used[affected],
               n_dropped = vapply(verdicts[affected], function(v)
                 v$n_dropped, integer(1)),
               formed = formed[affected],
               status = vapply(verdicts[affected], function(v) v$status,
                               character(1)),
               stringsAsFactors = FALSE, row.names = NULL))
  least
}

#' Tidy an MBCO Mediation Table
#'
#' Returns the effect rows of a \code{\link{mediation_mbco}} table in
#' the column convention used by the \pkg{broom} ecosystem. The
#' \code{statistic} column is the MBCO likelihood ratio statistic.
#'
#' @param x A \code{dmar_mediation_mbco} object returned by
#'   \code{\link{mediation_mbco}}.
#' @param \dots Unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'   \code{estimate}, \code{se}, \code{statistic},
#'   \code{p_value}, \code{ci_lower}, \code{ci_upper}, and
#'   \code{conf_level}, the level of the intervals.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_mediation_mbco <- function(x, ...) {
  data.frame(
    term      = x$term,
    estimate  = x$estimate,
    se = x$se,
    statistic = x$lrt,
    p_value   = x$p_value,
    ci_lower  = x$ci_lower,
    ci_upper = x$ci_upper,
    # The level of the intervals, so that a tidy row states the level of
    # its limits.
    conf_level = rep(attr(x, "conf_level", exact = TRUE) %||% NA_real_,
                     nrow(x)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Glance at an MBCO Mediation Fit
#'
#' Returns a one-row \code{data.frame} of full-model summaries from a
#' \code{\link{mediation_mbco}} table, in the column convention used by
#' the \pkg{broom} ecosystem.
#'
#' @param x A \code{dmar_mediation_mbco} object returned by
#'   \code{\link{mediation_mbco}}.
#' @param \dots Unused.
#'
#' @return A one-row \code{data.frame} with columns \code{nobs},
#'   \code{npar}, \code{deviance}, \code{AIC}, \code{BIC},
#'   \code{logLik}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_mediation_mbco <- function(x, ...) {
  dev <- attr(x, "deviance")
  data.frame(
    nobs     = attr(x, "N"),
    npar     = attr(x, "n_par"),
    deviance = dev,
    AIC      = attr(x, "aic"),
    BIC      = attr(x, "bic"),
    logLik   = -dev / 2,
    stringsAsFactors = FALSE
  )
}
