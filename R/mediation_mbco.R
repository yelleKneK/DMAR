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
#'   and its row reports the test alone, with no scalar estimate.
#' @param ci_method Confidence interval method for every reported
#'   effect: \code{"profile_likelihood"} (default; Neale & Miller, 1997,
#'   computed by \pkg{OpenMx}), \code{"monte_carlo"} (simulate the
#'   coefficients from their joint normal approximation and form
#'   percentile intervals of the effect; MacKinnon, Lockwood, &
#'   Williams, 2004), or \code{"wald"} (delta method symmetric
#'   interval; reported for comparison, not recommended for indirect
#'   effects, whose sampling distributions are skewed).
#' @param conf_level Confidence level. Defaults to 0.95.
#' @param B Number of Monte Carlo replications when
#'   \code{ci_method = "monte_carlo"}. Defaults to 10000.
#' @param optimizer Optimizer used by \pkg{OpenMx} for all fits:
#'   \code{"SLSQP"} (default), \code{"CSOLNP"}, or \code{"NPSOL"}.
#'   The CRAN build of \pkg{OpenMx} ships SLSQP and CSOLNP; NPSOL, the
#'   optimizer used in Tofighi and Kelley (2020), is available only in
#'   the build distributed from the OpenMx website.
#' @param seed Optional integer seed for the Monte Carlo interval, used
#'   locally (the caller's random number generator state is restored on
#'   exit). Default \code{NULL} leaves the random number generator state
#'   alone.
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
#' which has a large sample chi square distribution with degrees of
#' freedom equal to the number of constraints (here 1) when the null
#' hypothesis is true (Wilks, 1938). The resulting \emph{p}-value is a
#' continuous measure of the compatibility of the data with the null
#' model, not merely a reject or fail-to-reject decision, and in the
#' simulations of Tofighi and Kelley (2020) the test controls the Type I
#' error rate more robustly than the tests based on confidence intervals in
#' common use, especially when the indirect effect is truly zero.
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
#' approximation for the product breaks down. The MBCO procedure sidesteps
#' that approximation by comparing maximized likelihoods directly.
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
#' Kelley (2020) do for the memory data in the examples. The Wald
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
#' non-reference group is compared with the reference group.
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
#' indirect pathway between \code{x} and \code{y}. It warns where the
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
#' \eqn{\mathrm{LRT}_{\mathrm{MBCO}} - \log N} for the BIC. (The
#' \pkg{OpenMx} summary counts the same number of estimated parameters
#' in both models, so its printed AIC difference equals the likelihood
#' ratio statistic; DMAR counts the constraint against the null model,
#' consistent with the degrees of freedom of the test.)
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
#'   \code{df} equal to the number of constraints; its \code{estimate},
#'   \code{se}, and interval are \code{NA} because no single number
#'   summarizes several constraints. Attributes: \code{"conf_level"},
#'   \code{"ci_method"}, \code{"optimizer"}, \code{"x"}, \code{"y"},
#'   \code{"N"}, \code{"deviance"}, \code{"aic"}, \code{"bic"},
#'   \code{"n_par"} (full model fit information), \code{"groups"} (the
#'   group labels, reference first; multiple-group fits only),
#'   \code{"R2"} (named vector, endogenous variables under the full
#'   model, per group when grouped), \code{"delta_R2"} (matrix of full
#'   minus null \eqn{R^2}, variables by tested effects),
#'   \code{"moderation"} (for a moderated analysis: the moderator
#'   name, the probe values, and each moderated effect's moderator
#'   polynomial, which is what \code{\link{plot_mediation_mbco}}
#'   draws), and \code{"mx_model"} (the fitted \pkg{OpenMx} full
#'   model, an escape hatch for further OpenMx work). Use
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
#' # Single-mediator model: instruction -> imagery -> recall. The fit is
#' # not run here because every reported effect costs its own constrained
#' # null model fit in OpenMx, on top of the B Monte Carlo replications
#' # the interval draws. The call is:
#' # single <- "
#' #   imagery ~ b1*instruction
#' #   recall  ~ b2*imagery + b3*instruction
#' # "
#' # mediation_mbco(single, S = S_tk, M = M_tk, N = 369,
#' #                x = "instruction", y = "recall",
#' #                ci_method = "monte_carlo", seed = 113)
#' # The indirect effect is about 2.1 words (the paper reports 2.121,
#' # SE = 0.276, 95% Monte Carlo CI [1.600, 2.682]).
#'
#' # Parallel two-mediator model, with the contrast of the two specific
#' # indirect effects (the paper's Research Questions 2 and 3). Not run
#' # here for the same reason; the call is:
#' # parallel <- "
#' #   imagery    ~ b1*instruction
#' #   repetition ~ b3*instruction
#' #   recall     ~ b2*imagery + b4*repetition + b5*instruction
#' #   imagery ~~ repetition
#' # "
#' # mediation_mbco(parallel, S = S_tk, M = M_tk, N = 369,
#' #                x = "instruction", y = "recall",
#' #                hypotheses = c(imagery_minus_repetition =
#' #                  "indirect_via_imagery - indirect_via_repetition"),
#' #                ci_method = "monte_carlo", seed = 113)
#' # The indirect effect through repetition is near zero (the paper
#' # reports LRT = 0.083, p = .773), while the contrast shows the
#' # imagery pathway is larger (the paper reports LRT = 25.828,
#' # difference = 2.222, SE = 0.445).
#'
#' # Raw data go in through 'data' rather than 'S', 'M', and 'N'. Adding
#' # 'group' fits the model in every group and tests the between-group
#' # difference of each effect, which is moderated mediation with a
#' # categorical moderator. Leaving the paths unlabeled lets them differ
#' # by group. It is not run here because each reported difference costs
#' # its own constrained null model fit:
#' #   mediation_mbco("m ~ x \n y ~ m + x", data = two_groups,
#' #                  group = "condition", x = "x", y = "y",
#' #                  ci_method = "wald")
#' #
#' # A continuous moderator goes in through 'moderator'. That analysis,
#' # fit from a data frame with the conditional effects it estimates
#' # drawn as curves over the moderator's range, is shown at
#' # ?plot_mediation_mbco.
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
    stop("The package 'lavaan' is needed; please install the package ",
         "and try again.", call. = FALSE)
  }
  if (!requireNamespace("OpenMx", quietly = TRUE)) {
    stop("The package 'OpenMx' is needed; please install the package ",
         "and try again.", call. = FALSE)
  }
  ci_method <- match.arg(ci_method)
  optimizer <- match.arg(optimizer)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 100 ||
      B != round(B)) {
    stop("'B' must be a single integer of at least 100.", call. = FALSE)
  }
  # Local RNG: seed if asked (covering the Monte Carlo interval and any
  # optimizer restarts), and restore the caller's state on exit.
  if (!is.null(seed)) {
    has_old_seed <- exists(".Random.seed", envir = globalenv())
    old_seed <- if (has_old_seed) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit({
      if (has_old_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  if (!optimizer %in% OpenMx::mxAvailableOptimizers()) {
    stop("The '", optimizer, "' optimizer is not available in this ",
         "build of OpenMx. The CRAN build ships 'SLSQP' and 'CSOLNP'; ",
         "'NPSOL' requires the build distributed from the OpenMx ",
         "website.", call. = FALSE)
  }

  # ---- resolve the data --------------------------------------------------
  lav_input <- inherits(model, "lavaan")
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
      stop("The raw data could not be recovered from the fitted lavaan ",
           "object; supply 'data' (or 'S' and 'N') directly.",
           call. = FALSE)
    }
  }
  if (!is.null(data) && !is.null(S)) {
    stop("Supply either 'data' or 'S' (with 'N'), not both.",
         call. = FALSE)
  }
  if (is.null(data)) {
    if (is.null(S)) {
      stop("Supply raw data via 'data' or a covariance matrix via 'S' ",
           "(with 'N').", call. = FALSE)
    }
    if (is.list(S) && !is.data.frame(S)) {
      # Multiple groups from summary statistics: one covariance matrix
      # (and optionally one mean vector) per group.
      gl <- names(S)
      if (is.null(gl) || any(!nzchar(gl))) {
        stop("'S' as a list must have one named element per group.",
             call. = FALSE)
      }
      if (is.null(N) || length(N) != length(S)) {
        stop("'N' must give one sample size per group in 'S'.",
             call. = FALSE)
      }
      if (!is.null(names(N))) {
        if (!all(gl %in% names(N))) {
          stop("The names of 'N' must match the names of 'S'.",
               call. = FALSE)
        }
        N <- N[gl]
      }
      if (!is.null(M) && (!is.list(M) || !all(gl %in% names(M)))) {
        stop("'M' must be a named list with one mean vector per group ",
             "in 'S'.", call. = FALSE)
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
      stop("'group' must name one column of the data.", call. = FALSE)
    }
  }
  if (!is.null(moderator) && !is.null(group)) {
    stop("'moderator' and 'group' cannot yet be combined; probe one ",
         "moderator at a time.", call. = FALSE)
  }
  if (!is.null(probe_values) && is.null(moderator)) {
    stop("'probe_values' requires 'moderator'.", call. = FALSE)
  }

  # Interaction syntax: rewrite each 'a:b' term into a product column
  # so the model can state the moderation directly.
  if (is.character(model)) {
    prep <- .mbco_product_terms(model, data)
    model <- prep$model
    data <- prep$data
  }
  if (!is.null(moderator)) {
    if (!is.character(moderator) || length(moderator) != 1L ||
        !moderator %in% names(data)) {
      stop("'moderator' must name one column of the data.",
           call. = FALSE)
    }
    if (!is.numeric(data[[moderator]])) {
      stop("The moderator must be numeric; code a categorical ",
           "moderator through 'group' instead.", call. = FALSE)
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
      stop("'model' must be a lavaan model syntax string or a fitted ",
           "lavaan object.", call. = FALSE)
    }
    miss <- if (anyNA(data)) "ml" else "listwise"
    lav <- try(lavaan::sem(model, data = data, meanstructure = TRUE,
                           fixed.x = FALSE, missing = miss, se = "none",
                           group = group, warn = FALSE),
               silent = TRUE)
    if (inherits(lav, "try-error")) {
      stop("lavaan could not process 'model': ",
           attr(lav, "condition")$message, call. = FALSE)
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
    stop("The model uses operator(s) not supported by mediation_mbco(): ",
         paste(sQuote(bad_ops), collapse = ", "), ".", call. = FALSE)
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
    stop("Parameter labels must be valid R names; rename: ",
         paste(sQuote(not_names), collapse = ", "), ".", call. = FALSE)
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
      stop("Explicit '==' constraints are not supported; equate ",
           "parameters by giving them the same label.", call. = FALSE)
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
    stop("Variable(s) not found in the data: ",
         paste(sQuote(missing_vars), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!all(vapply(data[manifest], is.numeric, logical(1L)))) {
    stop("All modeled variables must be numeric; convert factors ",
         "before calling mediation_mbco().", call. = FALSE)
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
    warning("Endogenous variable(s) ",
            paste(sQuote(two_level), collapse = ", "),
            " take only two distinct values. The linear normal model ",
            "treats them as continuous, and for a binary mediator or ",
            "outcome the product of coefficients is not the causal ",
            "indirect effect; potential outcome methods are needed for ",
            "that reading.", call. = FALSE)
  }

  # ---- effects: enumerate pathways from x to y ---------------------------
  reg <- pt[pt$op == "~", , drop = FALSE]
  if (nrow(reg) == 0L) {
    stop("The model contains no regression ('~') paths; there is ",
         "nothing to mediate.", call. = FALSE)
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
      warning("Interaction term(s) detected among the predictors (",
              paste(unique(hits), collapse = "; "),
              "). The automatically enumerated effects average over ",
              "that moderation; declare the moderator through the ",
              "'moderator' argument to probe it, or define conditional ",
              "effects with ':=' or 'hypotheses' (see ",
              "?mediation_mbco).", call. = FALSE)
    }
  }

  # User-defined quantities and extra hypotheses join the effect list.
  # Each row's 'test' column names the algebra(s) its null model
  # constrains to zero; 'df' is the number of constraints. A row with
  # several constraints (a joint test) has no scalar estimate.
  hidden <- data.frame(name = character(0), expr = character(0),
                       stringsAsFactors = FALSE)
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
      stop("'hypotheses' must be a named character vector, or a named ",
           "list of character vectors for joint tests.", call. = FALSE)
    }
    hnames <- names(hypotheses)
    if (is.null(hnames)) hnames <- rep("", length(hypotheses))
    hnames[!nzchar(hnames)] <-
      paste0("hypothesis_", seq_len(sum(!nzchar(hnames))))
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
    stop("Effect names collide: ",
         paste(sQuote(effects$term[duplicated(effects$term)]),
               collapse = ", "),
         ". Rename the clashing ':=' definition or hypothesis.",
         call. = FALSE)
  }
  alg_defs <- rbind(
    data.frame(name = effects$term[!is.na(effects$expr)],
               expr = effects$expr[!is.na(effects$expr)],
               stringsAsFactors = FALSE),
    hidden)
  clash <- intersect(alg_defs$name, pt$mx_label[pt$op != ":="])
  if (length(clash) > 0L) {
    stop("Name(s) used both as a parameter label and as an effect: ",
         paste(sQuote(clash), collapse = ", "), ".", call. = FALSE)
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
      stop("Effect ", sQuote(alg_defs$name[i]),
           " references unknown parameter(s): ",
           paste(sQuote(unknown), collapse = ", "), ".", call. = FALSE)
    }
  }

  # ---- build and fit the full model in OpenMx ----------------------------
  old_opt <- OpenMx::mxOption(NULL, "Default optimizer")
  on.exit(OpenMx::mxOption(NULL, "Default optimizer", old_opt),
          add = TRUE)
  OpenMx::mxOption(NULL, "Default optimizer", optimizer)

  full <- .mbco_build_mx(pt, data, manifest, latents, alg_defs,
                         group = group, glabels = glabels)
  if (ci_method == "profile_likelihood") {
    full <- OpenMx::mxModel(full,
      OpenMx::mxCI(effects$term[!is.na(effects$expr)],
                   interval = conf_level))
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
      is.finite(D_lav) && !lav_input && D_full > D_lav + 0.01) {
    fit_full <- .mbco_try_hard(
      full, extraTries = 25L,
      intervals = ci_method == "profile_likelihood")
    D_full <- fit_full@output$Minus2LogLikelihood
    if (D_full > D_lav + 0.01) {
      warning("The OpenMx full-model deviance (", round(D_full, 3),
              ") did not reach the lavaan value (", round(D_lav, 3),
              "); results may reflect a local solution.", call. = FALSE)
    }
  }
  if (!fit_full@output$status$code %in% c(0L, 1L)) {
    warning("The full model finished with OpenMx status code ",
            fit_full@output$status$code,
            "; interpret the results with care.", call. = FALSE)
  }
  n_obs <- .mbco_numobs(fit_full)
  n_par <- length(coef(fit_full))
  aic_full <- D_full + 2 * n_par
  bic_full <- D_full + log(n_obs) * n_par

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
  if (ci_method == "profile_likelihood") {
    ints <- fit_full@output$confidenceIntervals
    for (i in seq_len(k)) {
      if (is.na(effects$expr[i])) next
      row <- grep(paste0("\\.", effects$term[i], "\\["), rownames(ints))
      if (length(row) == 0L) {
        row <- match(effects$term[i], rownames(ints))
      }
      lo <- if (length(row) >= 1L) ints[row[1L], "lbound"] else NA_real_
      hi <- if (length(row) >= 1L) ints[row[1L], "ubound"] else NA_real_
      if (is.na(lo) || is.na(hi)) {
        warning("The profile likelihood interval for ",
                sQuote(effects$term[i]),
                " did not fully converge; the Wald limit is reported ",
                "for the missing bound.", call. = FALSE)
      }
      if (!is.na(lo)) ci_lower[i] <- lo
      if (!is.na(hi)) ci_upper[i] <- hi
    }
  } else if (ci_method == "monte_carlo") {
    V <- try(suppressWarnings(vcov(fit_full)), silent = TRUE)
    if (inherits(V, "try-error")) {
      warning("The parameter covariance matrix is unavailable; Wald ",
              "intervals are reported instead of Monte Carlo.",
              call. = FALSE)
    } else {
      draws <- as.data.frame(
        MASS::mvrnorm(B, mu = coef(fit_full), Sigma = V))
      for (i in seq_len(k)) {
        if (is.na(effects$expr[i])) next
        vals <- eval(parse(text = expanded[[effects$term[i]]])[[1L]],
                     envir = draws)
        lims <- quantile(vals, c(alpha / 2, 1 - alpha / 2),
                         names = FALSE, na.rm = TRUE)
        ci_lower[i] <- lims[1L]
        ci_upper[i] <- lims[2L]
      }
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
                          extra_starts = eff$starts[[effects$term[i]]])
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
  attr(out, "deviance") <- D_full
  attr(out, "aic") <- aic_full
  attr(out, "bic") <- bic_full
  attr(out, "n_par") <- n_par
  attr(out, "R2") <- r2_full
  attr(out, "delta_R2") <- delta_r2
  attr(out, "mx_model") <- fit_full
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
    stop("Supply both 'x' and 'y', or neither.", call. = FALSE)
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
      stop("The focal predictor and outcome are ambiguous (sources: ",
           paste(sQuote(sources), collapse = ", "), "; outcomes: ",
           paste(sQuote(sinks), collapse = ", "),
           "); supply 'x' and 'y'.", call. = FALSE)
    }
  }
  for (nm in c(x, y)) {
    if (!nm %in% vars) {
      stop(sQuote(nm), " has no regression path in the model.",
           call. = FALSE)
    }
  }
  # A feedback loop makes "the pathways from x to y" ill defined (total
  # effects would be infinite series, not sums of simple paths), so
  # automatic enumeration refuses rather than mislabel.
  color <- stats::setNames(rep(0L, length(vars)), vars)
  cyc <- FALSE
  visit <- function(v) {
    color[[v]] <<- 1L
    for (w in to[from == v]) {
      if (cyc) return(invisible(NULL))
      if (color[[w]] == 1L) {
        cyc <<- TRUE
        return(invisible(NULL))
      }
      if (color[[w]] == 0L) visit(w)
    }
    color[[v]] <<- 2L
  }
  for (v in vars) if (color[[v]] == 0L) visit(v)
  if (cyc) {
    stop("The regression structure is nonrecursive (a feedback loop ",
         "connects the equations), so enumerating pathways from ",
         sQuote(x), " to ", sQuote(y), " is not meaningful. State the ",
         "quantities to test through ':=' definitions or 'hypotheses' ",
         "and omit 'x' and 'y'.", call. = FALSE)
  }
  # Depth-first enumeration of simple directed pathways x -> ... -> y.
  paths <- list()
  walk <- function(node, seen) {
    for (i in which(from == node)) {
      nxt <- to[i]
      if (nxt %in% seen) next
      if (nxt == y) {
        paths[[length(paths) + 1L]] <<- c(seen, nxt)
      } else {
        walk(nxt, c(seen, nxt))
      }
    }
  }
  walk(x, x)
  if (length(paths) == 0L) {
    stop("No directed pathway from ", sQuote(x), " to ", sQuote(y),
         " exists in the model.", call. = FALSE)
  }
  n_edges  <- lengths(paths) - 1L
  direct   <- paths[n_edges == 1L]
  indirect <- paths[n_edges > 1L]
  if (length(indirect) == 0L && !have_extra) {
    stop("Only the direct path connects ", sQuote(x), " and ",
         sQuote(y), "; there is no indirect pathway to test.",
         call. = FALSE)
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
# columns named a_by_b, so the model can state moderation directly.
.mbco_product_terms <- function(model, data) {
  pat <- "([A-Za-z][A-Za-z0-9._]*)[ \t]*:[ \t]*([A-Za-z][A-Za-z0-9._]*)"
  hits <- regmatches(model, gregexpr(pat, model))[[1L]]
  for (h in unique(hits)) {
    parts <- strsplit(gsub("[ \t]", "", h), ":", fixed = TRUE)[[1L]]
    if (!all(parts %in% names(data))) next
    nm <- paste0(parts[1L], "_by_", parts[2L])
    prod <- data[[parts[1L]]] * data[[parts[2L]]]
    if (nm %in% names(data)) {
      if (!isTRUE(all.equal(data[[nm]], prod, tolerance = 1e-10))) {
        stop("Cannot create the interaction column ", sQuote(nm),
             " for the term ", sQuote(h), "; a different column with ",
             "that name already exists.", call. = FALSE)
      }
    } else {
      data[[nm]] <- prod
    }
    model <- gsub(h, nm, model, fixed = TRUE)
  }
  list(model = model, data = data)
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
    stop("The moderator ", sQuote(moderator), " is endogenous (it has ",
         "an incoming regression path); probing requires an exogenous ",
         "moderator.", call. = FALSE)
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
    warning("No interaction between ", sQuote(moderator), " and a ",
            "modeled variable was found in the model; the moderator ",
            "is ignored. State the interaction as, e.g., 'm ~ x + w ",
            "+ x:w' or include the product column as a predictor.",
            call. = FALSE)
    return(NULL)
  }
  if (length(no_main) > 0L) {
    warning("Interaction without the matching main effect of ",
            paste(unique(no_main), collapse = ", "),
            " (the principle of marginality); conditional effects ",
            "along that path are not probed.", call. = FALSE)
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
      stop("'probe_values' must be a numeric vector without missing ",
           "values.", call. = FALSE)
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
                           label = terms[1L], extra_starts = NULL) {
  cons <- lapply(seq_along(terms), function(j) {
    eval(bquote(OpenMx::mxConstraint(.(as.name(terms[j])) == 0,
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
  # per branch of the constraint's null set.
  for (a in extra_starts) {
    s <- est
    keep <- intersect(names(a), names(est))
    if (length(keep) == 0L) next
    vals <- vapply(a[keep], function(e) {
      out <- try(eval(parse(text = e)[[1L]], envir = as.list(est)),
                 silent = TRUE)
      if (inherits(out, "try-error") || !is.finite(out)) NA_real_
      else as.numeric(out)
    }, numeric(1L))
    if (anyNA(vals)) next
    s[keep] <- vals
    starts[[length(starts) + 1L]] <- s
  }
  # Feasibility tolerance for each constraint at the solution: strict
  # enough to reject convergence failures that drift off the null set
  # (whose deviance can undercut the true constrained minimum), loose
  # enough to accept the small residual a converged sequential
  # quadratic programming step leaves.
  tol <- 1e-4 * pmax(1, abs(est_full))
  cand_of <- function(f) {
    if (inherits(f, "try-error") || is.null(f)) return(NULL)
    g <- rep(NA_real_, length(terms))
    for (j in seq_along(terms)) {
      gj <- try(abs(as.numeric(
        OpenMx::mxEvalByName(terms[j], f, compute = TRUE)[1L])),
        silent = TRUE)
      if (inherits(gj, "try-error") || !is.finite(gj)) return(NULL)
      g[j] <- gj
    }
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
  cands <- do.call(c, lapply(starts, run_one))
  # Deterministic branch fits. When every constrained algebra is a pure
  # product of free parameters, the null set is exactly the union of
  # branches obtained by fixing one factor of each product to zero, and
  # each branch is an ordinary unconstrained fit: no sequential
  # quadratic programming step that a different optimizer build might
  # carry into the wrong basin. The branch fits are added as candidates
  # beside the constrained runs, so the reported likelihood ratio is
  # defined by the best-fitting branch on every platform and OpenMx
  # version, not by whichever basin the constrained search happens to
  # reach. Non-product constraints (a sum of indirect effects, for
  # example) skip this path, because zeroing one factor does not cover
  # their null set; they keep the multi-start constrained machinery.
  factors_of <- function(e) {
    if (is.name(e)) return(as.character(e))
    if (is.call(e) && identical(e[[1L]], as.name("("))) {
      return(factors_of(e[[2L]]))
    }
    if (is.call(e) && identical(e[[1L]], as.name("*"))) {
      lhs <- factors_of(e[[2L]])
      rhs <- factors_of(e[[3L]])
      if (is.null(lhs) || is.null(rhs)) return(NULL)
      return(c(lhs, rhs))
    }
    NULL
  }
  term_factors <- lapply(expanded, function(e) {
    f <- factors_of(parse(text = e)[[1L]])
    if (is.null(f) || !all(f %in% names(est))) NULL else unique(f)
  })
  if (!any(vapply(term_factors, is.null, logical(1L))) &&
      prod(lengths(term_factors)) <= 64L) {
    combos <- expand.grid(term_factors, stringsAsFactors = FALSE)
    branch_sets <- unique(lapply(seq_len(nrow(combos)), function(i) {
      sort(unique(as.character(unlist(combos[i, ], use.names = FALSE))))
    }))
    for (labs in branch_sets) {
      m <- OpenMx::omxSetParameters(fit_full, labels = labs,
                                    values = 0, free = FALSE)
      m <- OpenMx::mxModel(m, name = paste0(
        "mbco_branch_", gsub("[^A-Za-z0-9_]", "_", label)))
      f <- try(suppressWarnings(OpenMx::mxRun(m, silent = TRUE)),
               silent = TRUE)
      z <- cand_of(f)
      if (!is.null(z)) cands[[length(cands) + 1L]] <- z
    }
  }
  valid <- function(z, codes) {
    all(z$g <= tol) && z$code %in% codes && z$dev >= D_full - 1e-4
  }
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
      warning("The null model for ", sQuote(label),
              " converged with OpenMx status code ",
              ok[[which.min(vapply(ok, `[[`, numeric(1L), "dev"))]]$code,
              "; interpret its test with care.", call. = FALSE)
    }
  }
  if (length(ok) == 0L) {
    warning("The null model for ", sQuote(label),
            " could not be fit subject to the constraint; its test is ",
            "reported as NA.", call. = FALSE)
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
    stop("'S' must have dimnames naming the variables.", call. = FALSE)
  }
  if (is.null(dimnames(S)[[1L]])) rownames(S) <- colnames(S)
  if (is.null(dimnames(S)[[2L]])) colnames(S) <- rownames(S)
  if (!isSymmetric(unname(S), tol = 1e-8)) {
    stop("'S' must be a symmetric covariance matrix.", call. = FALSE)
  }
  if (is.null(N) || !is.numeric(N) || length(N) != 1L || is.na(N) ||
      N <= ncol(S) || N != round(N)) {
    stop("'N' must be a single whole number larger than the number ",
         "of variables when 'S' is supplied.", call. = FALSE)
  }
  mu <- rep(0, ncol(S))
  names(mu) <- colnames(S)
  if (!is.null(M)) {
    if (is.null(names(M)) || !all(colnames(S) %in% names(M))) {
      stop("'M' must be a named vector covering the variables in 'S'.",
           call. = FALSE)
    }
    mu <- M[colnames(S)]
  }
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    stop("'S' is not positive definite; check the covariance ",
         "(or correlation) matrix.", call. = FALSE)
  }
  .mbco_local_rng(1L, {
    as.data.frame(MASS::mvrnorm(N, mu = mu, Sigma = S,
                                empirical = TRUE))
  })
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

# Evaluate an expression under a locally seeded RNG, restoring the
# caller's state on exit. With seed = NULL the state is still saved and
# restored, so internal draws leave the caller's stream untouched.
.mbco_local_rng <- function(seed, expr) {
  has_old <- exists(".Random.seed", envir = globalenv())
  old <- if (has_old) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (has_old) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv())) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  expr
}

#' @export
print.dmar_mediation_mbco <- function(x, ...) {
  cat("Mediation tests via model-based constrained optimization",
      "(MBCO)\n")
  # A column subset loses the fit attributes; print only what exists.
  xv <- attr(x, "x")
  yv <- attr(x, "y")
  if (!is.null(attr(x, "N"))) {
    if (!is.null(xv) && !is.null(yv)) {
      cat("  Effects of ", xv, " on ", yv, sep = "")
    } else {
      cat("  User-defined effects")
    }
    cat("; N = ", attr(x, "N"), "; optimizer: ", attr(x, "optimizer"),
        "\n", sep = "")
  }
  if (!is.null(attr(x, "groups"))) {
    cat("  Groups (reference first): ",
        paste(attr(x, "groups"), collapse = ", "), "\n", sep = "")
  }
  if (!is.null(attr(x, "deviance"))) {
    cat("  Full model: deviance = ", format(attr(x, "deviance"),
                                            nsmall = 3, digits = 10),
        ", AIC = ", format(attr(x, "aic"), nsmall = 3, digits = 10),
        ", BIC = ", format(attr(x, "bic"), nsmall = 3, digits = 10),
        "\n", sep = "")
  }
  r2 <- attr(x, "R2")
  if (length(r2) > 0L) {
    cat("  R2 (full model): ",
        paste(names(r2), "=", formatC(r2, digits = 3, format = "f"),
              collapse = ", "), "\n", sep = "")
  }
  cat("  Each p_value is from the MBCO likelihood ratio test of the\n",
      "  null model constraining that effect to zero. Causal readings\n",
      "  rest on the no omitted confounder assumption.\n\n", sep = "")
  NextMethod()
  invisible(x)
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
#'   \code{p_value}, \code{ci_lower}, \code{ci_upper}.
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
